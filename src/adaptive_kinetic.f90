module adaptive_kinetic
  !! Runtime harmonic ranges around the shared, unchanged susceptibility algebra.
  use vars
  use kinetic_core, only: kinetic_normalization
  use harmonic_selection, only: harmonic_tail_bounds
  use bessel_cache
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private
  public :: dynamic_harmonics, calc_chi_adaptive, disp_adaptive, validate_bessel_input
  public :: physical_bessel_argument, finish_bessel_diagnostics
  abstract interface
     real function bessel_function(n,x)
       integer :: n
       real :: x
     end function
     complex function zeta_function(k,z)
       real :: k
       complex :: z
     end function
  end interface
contains
  logical function dynamic_harmonics()
    ! Fixed +/-128 preserves the original compiler context and Bessel evaluator.
    dynamic_harmonics=bessel_mode/='fixed' .or. bessel_nfixed/=128
  end function

  subroutine validate_bessel_input()
    if (bessel_mode=='fixed128') then
       if (bessel_nfixed/=128) error stop 'ERROR: fixed128 requires bessel_nfixed=128.'
       bessel_mode='fixed'
    endif
    if (bessel_mode/='adaptive' .and. bessel_mode/='fixed') &
         error stop 'ERROR: bessel_mode must be adaptive or fixed.'
    if (.not.ieee_is_finite(bessel_tol)) error stop 'ERROR: bessel_tol must be finite.'
    if (bessel_tol<=0. .or. bessel_tol>=1.) error stop 'ERROR: require 0 < bessel_tol < 1.'
    ! Integer n*n occurs in the historical equations; reject integer overflow.
    if (bessel_nmax<1 .or. bessel_nmax>int(sqrt(real(huge(1))))-2) &
         error stop 'ERROR: bessel_nmax outside representable harmonic range (1..46338 on int32).'
    if (bessel_nfixed<1 .or. bessel_nfixed>int(sqrt(real(huge(1))))-2) &
         error stop 'ERROR: bessel_nfixed outside representable harmonic range (1..46338 on int32).'
  end subroutine

  real function physical_bessel_argument(species,kperp_input,alpha_reference,inertial) result(x)
    type(specie), intent(in) :: species
    real, intent(in) :: kperp_input,alpha_reference
    logical, intent(in) :: inertial
    real :: perpendicular_energy
    logical :: halting(3),flags(3),raised(3)
    type(ieee_flag_type), parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_get_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    ! One physical definition: k_perp**2 * v_tperp,s**2 / (2 Omega_s**2).
    ! rho_ref supplies inverse temperature and the perpendicular reference scale.
    if (inertial) then
       perpendicular_energy=species%alph_s*species%tau_s
    else
       perpendicular_energy=species%alph_s/(species%tau_s*alpha_reference)
    endif
    x=(kperp_input**2/2.)*species%Q_s**2*perpendicular_energy/species%mu_s
    call ieee_get_flag(exceptions,raised)
    if (any(raised)) x=ieee_value(0.,ieee_quiet_nan)
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end function

  subroutine calc_chi_adaptive(species,omega,state,argument,slot,count,collisions,compute_low, &
       chi,chi_low,ierr,bessel,zet_in)
    type(specie), intent(in) :: species
    complex, intent(in) :: omega
    type(kinetic_normalization), intent(in) :: state
    real, intent(in) :: argument
    integer, intent(in) :: slot,count,collisions
    logical, intent(in) :: compute_low
    complex, intent(out) :: chi(3,3),chi_low(3,3,-1:1)
    integer, intent(out) :: ierr
    procedure(bessel_function) :: bessel
    procedure(zeta_function) :: zet_in
    type(bessel_entry), pointer :: entry
    logical :: halting(3),flags(3),raised(3)
    type(ieee_flag_type), parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
    real :: nan
    nan=ieee_value(0.,ieee_quiet_nan)
    chi=cmplx(nan,nan); chi_low=chi(1,1); ierr=1
    if (.not.all(ieee_is_finite([real(omega),aimag(omega),argument]))) return
    if (argument<=0. .or. state%kpar==0. .or. omega==cmplx(0.,0.)) return
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_get_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    call fetch_bessel(slot,count,argument,bessel_mode,bessel_tol,bessel_nmax,bessel_nfixed,entry,ierr)
    if (ierr==0) call evaluate(species,omega,state,entry,collisions,compute_low,chi,chi_low,bessel,zet_in)
    call ieee_get_flag(exceptions,raised)
    if (any(raised) .or. .not.all(ieee_is_finite(real(chi))) .or. &
         .not.all(ieee_is_finite(aimag(chi)))) then
       ierr=5; chi=cmplx(nan,nan); chi_low=chi(1,1)
    endif
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine

  subroutine evaluate(species,om,state,entry,collisions,low_n,chi,chi_low,bessel,zet_in)
    type(specie), intent(in) :: species
    complex, intent(in) :: om
    type(kinetic_normalization), intent(in) :: state
    type(bessel_entry), target, intent(inout) :: entry
    integer, intent(in) :: collisions
    logical, intent(in) :: low_n
    complex, intent(out) :: chi(3,3),chi_low(3,3,-1:1)
    procedure(bessel_function) :: bessel
    procedure(zeta_function) :: zet_in
    real :: betap,vtp,kpar,kperp,alphp,lambdap,lambdas,Vdrifts
    real :: disp_tau,disp_mu,disp_alph,disp_Q,disp_D,disp_v,disp_nuns
    real, pointer :: jn(:,:),jpn(:,:)
    complex, pointer :: tsi(:),zz(:)
    complex :: norm(6),eps_xx_t,eps_xy_t,eps_xz_t,eps_yy_t,eps_yz_t,eps_zz_t
    complex :: susc(1,3,3),susc_low(1,3,3,-1:1)
    complex, parameter :: c_i=(0.,1.)
    logical :: reuse_bessel(1)
    integer :: is,n,nbrack
    betap=state%beta; vtp=state%velocity; kpar=state%kpar; kperp=state%kperp; alphp=state%alpha
    lambdap=kperp**2./2.
    disp_tau=species%tau_s; disp_mu=species%mu_s; disp_alph=species%alph_s
    disp_Q=species%Q_s; disp_D=species%D_s; disp_v=species%vv_s; disp_nuns=species%nu_ns
    is=1; nbrack=entry%active; reuse_bessel=.true.; susc_low=0.
    jn=>entry%jn; jpn=>entry%jpn; tsi=>entry%tsi; zz=>entry%zz
    if (collisions==0) then
       include 'kinetic_species.inc'
    else
       include 'kinetic_collisional_species.inc'
    endif
    chi=susc(1,:,:); chi_low=susc_low(1,:,:,:)
  end subroutine

  complex function disp_adaptive(omega,bessel,zet_in) result(value)
    complex, intent(in) :: omega
    procedure(bessel_function) :: bessel
    procedure(zeta_function) :: zet_in
    type(kinetic_normalization) :: state
    complex :: response(3,3),low(3,3,-1:1),a(3,3),enx2,enz2,enxnz
    real :: argument,nan
    integer :: j
    state=kinetic_normalization(solver_betap,solver_vtp,solver_kpar,solver_kperp,solver_spec(1)%alph_s)
    dispersion_status=0
    do j=1,nspec
       if (spec(j)%tau_s==0.) cycle
       argument=physical_bessel_argument(spec(j),kperp,spec(1)%alph_s,normalization=='d_ref')
       call calc_chi_adaptive(solver_spec(j),omega,state,argument,j,nspec,collision_type,low_n, &
            response,low,dispersion_status,bessel,zet_in)
       if (dispersion_status/=0) then
          nan=ieee_value(0.,ieee_quiet_nan)
          value=cmplx(nan,nan); susc=value; lam=value
          if (low_n) susc_low=value
          return
       endif
       susc(j,:,:)=response
       if (low_n) susc_low(j,:,:,:)=low
    enddo
    a=sum(susc,dim=1)
    enx2=((state%kperp/(omega*state%velocity))**2.)/state%alpha
    enz2=((state%kpar/(omega*state%velocity))**2.)/state%alpha
    enxnz=(state%kperp*state%kpar/(omega*state%velocity)**2.)/state%alpha
    a(1,1)=1.+a(1,1)-enz2
    a(2,2)=1.+a(2,2)-enx2-enz2
    a(3,3)=1.+a(3,3)-enx2
    a(1,3)=a(1,3)+enxnz; a(3,1)=a(1,3)
    lam=a
    value=a(1,1)*(a(2,2)*a(3,3)+a(2,3)**2.)+ &
         2.*a(1,2)*a(2,3)*a(1,3)-a(1,3)**2.*a(2,2)+a(1,2)**2.*a(3,3)
  end function

  subroutine finish_bessel_diagnostics()
    integer :: unit,status,j
    logical :: unresolved
    real :: argument,e0,e2
    unresolved=allocated(cap_events)
    if (allocated(cache)) then
       do j=1,size(cache)
          if (cache(j)%valid .and. .not.cache(j)%converged) unresolved=.true.
       enddo
    endif
    if (.not.dynamic_harmonics()) then
       do j=1,nspec
          if (spec(j)%tau_s==0.) cycle
          argument=physical_bessel_argument(spec(j),kperp,spec(1)%alph_s,normalization=='d_ref')
          if (.not.ieee_is_finite(argument)) then
             unresolved=.true.
             cycle
          endif
          call harmonic_tail_bounds(argument,128,e0,e2)
          if (max(e0,e2)>=bessel_tol) unresolved=.true.
       enddo
    endif
    ! Always leave an auditable mode/tolerance record; full statistics optional
    ! except that unresolved sums must be visible even without diagnostics.
    open(newunit=unit,file='data/'//trim(dataName)//'/'//trim(outputName)//'.bessel', &
         status='replace',action='write',iostat=status)
    if (status/=0) then
       write(*,'(a)') 'WARNING: unable to write Bessel metadata.'
       return
    endif
    write(unit,'(a,a)') 'bessel_mode=',trim(bessel_mode)
    write(unit,'(a,es24.15)') 'bessel_tol=',bessel_tol
    write(unit,'(a,i0)') 'bessel_nmax=',bessel_nmax
    write(unit,'(a,i0)') 'bessel_nfixed=',bessel_nfixed
    write(unit,'(a,l1)') 'unresolved=',unresolved
    if (.not.dynamic_harmonics()) write(unit,'(a)') 'historical fixed-128 evaluator; cache statistics unavailable'
    if (bessel_diagnostics .or. unresolved) call write_bessel_stats(unit)
    write(unit,'(2a)') '# z_backend=',trim(z_backend)
    close(unit)
  end subroutine
end module adaptive_kinetic
