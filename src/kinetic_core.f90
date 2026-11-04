module kinetic_core
  !! Explicit ordinary kinetic susceptibility; shares the legacy expressions.
  !! No physical input globals, reference-species lookup, or persistent cache.
  use nrtype, only: SP, SPC
  use vars, only: specie
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private
  public :: kinetic_normalization, calc_chi_kinetic, guarded_parallel
  include 'kinetic_harmonics.inc'
  type :: kinetic_normalization
     real(SP) :: beta, velocity, kpar, kperp, alpha
  end type
  abstract interface
     real function bessel_function(n,x)
       integer :: n
       real :: x
     end function
     complex function zeta_function(kpar,z)
       real :: kpar
       complex :: z
     end function
  end interface
contains
  subroutine calc_chi_kinetic(species,omega,state,chi,ierr,bessel,zet_in)
    type(specie), intent(in) :: species
    complex(SPC), intent(in) :: omega
    type(kinetic_normalization), intent(in) :: state
    complex(SPC), intent(out) :: chi(3,3)
    integer, intent(out) :: ierr !! 0 valid, 1 invalid input, 5 arithmetic failure.
    procedure(bessel_function) :: bessel
    procedure(zeta_function) :: zet_in
    logical :: halting(3), flags(3), raised(3)
    type(ieee_flag_type), parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
    real(SP) :: nan
    nan=ieee_value(0._SP,ieee_quiet_nan)
    chi=cmplx(nan,nan,SPC)
    ierr=1
    if (.not. all(ieee_is_finite([real(omega),aimag(omega),state%beta,state%velocity, &
         state%kpar,state%kperp,state%alpha,species%tau_s,species%mu_s, &
         species%alph_s,species%Q_s,species%D_s,species%vv_s]))) return
    if (state%beta<=0. .or. state%velocity<=0. .or. state%alpha<=0. .or. &
         species%tau_s<=0. .or. species%mu_s<=0. .or. species%alph_s<=0. .or. &
         species%Q_s==0. .or. species%D_s<0. .or. state%kpar==0. .or. state%kperp==0. .or. &
         omega==cmplx(0.,0.)) return
    ! Return failures even with caller traps enabled. Restore its exception state
    ! after the isolated calculation so a failed trial cannot poison the next one.
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_get_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    call evaluate(species,omega,state,chi,bessel,zet_in)
    call ieee_get_flag(exceptions,raised)
    ierr=0
    if (any(raised) .or. .not. all(ieee_is_finite(real(chi))) .or. &
         .not. all(ieee_is_finite(aimag(chi)))) then
       ierr=5
       chi=cmplx(nan,nan,SPC)
    endif
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine

  subroutine guarded_parallel(omega,value,ierr,evaluation)
    !! Protect the existing native parallel limit only for the new pathway.
    !! Its equations and the warm-reference caller are unchanged.
    complex, intent(in) :: omega
    complex, intent(out) :: value
    integer, intent(out) :: ierr
    interface
       complex function evaluation(w)
         complex, intent(in) :: w
       end function
    end interface
    logical :: halting(3),flags(3),raised(3)
    type(ieee_flag_type), parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_get_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    value=evaluation(omega)
    call ieee_get_flag(exceptions,raised)
    ierr=0
    if (any(raised) .or. .not. ieee_is_finite(real(value)) .or. .not. ieee_is_finite(aimag(value))) ierr=5
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine guarded_parallel

  subroutine evaluate(species,om,state,chi,bessel,zet_in)
    type(specie), intent(in) :: species
    complex, intent(in) :: om
    type(kinetic_normalization), intent(in) :: state
    complex, intent(out) :: chi(3,3)
    procedure(bessel_function) :: bessel
    procedure(zeta_function) :: zet_in
    real :: betap,vtp,kpar,kperp,alphp,lambdap,lambdas,Vdrifts
    real :: disp_tau,disp_mu,disp_alph,disp_Q,disp_D,disp_v
    real :: jn(-nbrack-1:nbrack+1,1),jpn(-nbrack:nbrack,1)
    complex :: tsi(-nbrack:nbrack),zz(-nbrack:nbrack),norm(6)
    complex :: eps_xx_t,eps_xy_t,eps_xz_t,eps_yy_t,eps_yz_t,eps_zz_t
    complex :: susc(1,3,3),susc_low(1,3,3,-1:1)
    complex, parameter :: c_i=(0.,1.)
    logical, parameter :: low_n=.false.
    logical :: reuse_bessel(1)
    integer :: is,n
    betap=state%beta; vtp=state%velocity; kpar=state%kpar; kperp=state%kperp; alphp=state%alpha
    lambdap=kperp**2./2.
    disp_tau=species%tau_s; disp_mu=species%mu_s; disp_alph=species%alph_s
    disp_Q=species%Q_s; disp_D=species%D_s; disp_v=species%vv_s
    is=1
    reuse_bessel=.false.
    include 'kinetic_species.inc'
    chi=susc(1,:,:)
  end subroutine
end module kinetic_core
