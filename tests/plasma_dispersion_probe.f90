module plasma_dispersion_probe
  !! Audit-only C ABI. Linked to unchanged PLUME modules in a temporary build.
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  use vars, only: specie,z_backend
  use disprels, only: production_z=>zet_in,bessel
  use audit_zet_variant
  implicit none
  integer :: z_method,call_index
  procedure(reference_interface), pointer :: reference_z
  complex, allocatable :: captured_argument(:),captured_value(:)
  abstract interface
     subroutine reference_interface(k,x,y,vr,vi) bind(C)
       import c_double
       real(c_double),value :: k,x,y
       real(c_double),intent(out) :: vr,vi
     end subroutine
  end interface
contains
  subroutine audit_z(k,x,y,tolerance,method,repeats,value,info) bind(C)
    real(c_double),value :: k,x,y,tolerance
    integer(c_int),value :: method,repeats
    real(c_double),intent(out) :: value(3)
    integer(c_int),intent(out) :: info(8)
    complex :: z,v,diagnostic
    real :: started,finished
    logical :: flags(4),oldflags(4),halting(4)
    type(ieee_flag_type),parameter :: exceptions(4)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow,ieee_underflow]
    integer :: j
    z_backend='legacy'
    z=cmplx(x,y);audit_tolerance=tolerance
    call ieee_get_flag(exceptions,oldflags)
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_set_halting_mode(exceptions,.false.)
    diagnostic=variant_in(k,z)
    info(1:4)=[audit_branch,audit_iterations,audit_stop,0]
    call ieee_set_flag(exceptions,.false.)
    call cpu_time(started)
    do j=1,repeats
       if (method==0) then
          v=production_z(k,z)
       else
          v=variant_in(k,z)
       endif
    enddo
    call cpu_time(finished)
    value=[real(v),aimag(v),finished-started]
    call ieee_get_flag(exceptions,flags)
    info(5:8)=merge(1,0,flags)
    if (.not.ieee_is_finite(real(v)) .or. .not.ieee_is_finite(aimag(v))) info(4)=1
    call ieee_set_flag(exceptions,oldflags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine

  subroutine audit_tensor(species,parameters,n,collisions,method,tolerance,gamma,callback,repeats, &
       output,arguments,values,info,seconds,low_output) bind(C)
    real(c_double),intent(in) :: species(7),parameters(7),gamma(0:n+1)
    integer(c_int),value :: n,collisions,method,repeats
    real(c_double),value :: tolerance
    type(c_funptr),value :: callback
    real(c_double),intent(out) :: output(18),arguments(2,2*n+1),values(2,2*n+1),seconds
    real(c_double),intent(out) :: low_output(54)
    integer(c_int),intent(out) :: info(6)
    type(specie) :: s
    complex :: chi(3,3),low(3,3,-1:1)
    real :: started,finished
    logical :: flags(4),oldflags(4),halting(4)
    type(ieee_flag_type),parameter :: exceptions(4)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow,ieee_underflow]
    integer :: iteration,i,j,index,harmonic
    z_backend='legacy'
    s=specie(species(1),species(2),species(3),species(4),species(5),species(6),species(7))
    z_method=method;audit_tolerance=tolerance
    if (method==2) call c_f_procpointer(callback,reference_z)
    if (allocated(captured_argument)) deallocate(captured_argument,captured_value)
    allocate(captured_argument(2*n+1),captured_value(2*n+1))
    call ieee_get_flag(exceptions,oldflags)
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    info=0
    call cpu_time(started)
    do iteration=1,repeats
       call_index=0
       call shared_tensor(s,parameters,n,collisions,gamma,chi,low)
    enddo
    index=0
    do harmonic=-1,1
       do i=1,3
          do j=1,3
             low_output(index+1:index+2)=[real(low(i,j,harmonic)),aimag(low(i,j,harmonic))]
             index=index+2
          enddo
       enddo
    enddo
    call cpu_time(finished)
    seconds=finished-started
    index=0
    do i=1,3
       do j=1,3
          output(index+1:index+2)=[real(chi(i,j)),aimag(chi(i,j))]
          index=index+2
       enddo
    enddo
    arguments(1,:)=real(captured_argument);arguments(2,:)=aimag(captured_argument)
    values(1,:)=real(captured_value);values(2,:)=aimag(captured_value)
    call ieee_get_flag(exceptions,flags)
    info(1:4)=merge(1,0,flags)
    info(5)=call_index
    if (.not.all(ieee_is_finite(output))) info(6)=1
    call ieee_set_flag(exceptions,oldflags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine

  complex function audit_zet(k,z) result(v)
    real :: k
    complex :: z
    real(c_double) :: vr,vi
    if (z_method==0) then
       v=production_z(k,z)
    else if (z_method==1) then
       v=variant_in(k,z)
    else
       call reference_z(k,real(z),aimag(z),vr,vi)
       v=cmplx(vr,vi)
    endif
    call_index=call_index+1
    captured_argument(call_index)=z
    captured_value(call_index)=v
  end function

  subroutine shared_tensor(species,parameters,nbrack,collisions,coefficients,chi,low)
    !! Holds all Bessel coefficients and N fixed. Includes production equations
    !! verbatim; only the Z callback is selected by the diagnostic caller.
    type(specie),intent(in) :: species
    real,intent(in) :: parameters(7),coefficients(0:nbrack+1)
    integer,intent(in) :: nbrack,collisions
    complex,intent(out) :: chi(3,3),low(3,3,-1:1)
    real :: betap,vtp,kpar,kperp,alphp,lambdap,lambdas,Vdrifts
    real :: disp_tau,disp_mu,disp_alph,disp_Q,disp_D,disp_v,disp_nuns
    real :: jn(-nbrack-1:nbrack+1,1),jpn(-nbrack:nbrack,1)
    complex :: tsi(-nbrack:nbrack),zz(-nbrack:nbrack),norm(6),om
    complex :: eps_xx_t,eps_xy_t,eps_xz_t,eps_yy_t,eps_yz_t,eps_zz_t
    complex :: susc(1,3,3),susc_low(1,3,3,-1:1)
    complex,parameter :: c_i=(0.,1.)
    logical,parameter :: low_n=.true.
    logical :: reuse_bessel(1)
    integer :: is,n
    om=cmplx(parameters(1),parameters(2))
    betap=parameters(6);alphp=parameters(7)
    vtp=parameters(5)*sqrt(betap)
    kpar=parameters(3)*sqrt(betap*alphp)
    kperp=parameters(4)*sqrt(betap*alphp)
    lambdap=kperp**2./2.
    disp_tau=betap/species%tau_s;disp_mu=species%mu_s;disp_alph=species%alph_s
    disp_Q=species%Q_s;disp_D=species%D_s;disp_v=species%vv_s;disp_nuns=species%nu_ns
    is=1;reuse_bessel=.true.
    ! Poison every slot to detect missing/overwritten harmonic assignments.
    susc_low=cmplx(ieee_value(0.,ieee_quiet_nan),ieee_value(0.,ieee_quiet_nan))
    do n=-nbrack-1,nbrack+1
       jn(n,1)=coefficients(abs(n))
    enddo
    do n=-nbrack,nbrack
       jpn(n,1)=.5*(jn(n-1,1)+jn(n+1,1))
    enddo
    if (collisions==0) then
       call collisionless()
    else
       call collisional()
    endif
    chi=susc(1,:,:)
    low=susc_low(1,:,:,:)
  contains
    subroutine collisionless()
      include 'kinetic_species.inc'
    end subroutine
    subroutine collisional()
      include 'kinetic_collisional_species.inc'
    end subroutine
    complex function zet_in(k,z)
      real :: k
      complex :: z
      zet_in=audit_zet(k,z)
    end function
  end subroutine
end module plasma_dispersion_probe
