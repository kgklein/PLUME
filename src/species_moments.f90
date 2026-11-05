module species_moments
  !! Species moments in the historical units v=deltaU/(c*deltaEx/B0),
  !! density=(delta n/n0)/(deltaEx/B0). All inputs use physical d_ref units.
  use vars, only: specie
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private
  public :: calc_cold_moments, calc_kinetic_moments, calc_species_moments
  integer, parameter, public :: moment_ok=0, moment_invalid=1, moment_zero_frequency=2, &
       moment_zero_doppler=3, moment_cyclotron_pole=4, moment_numerical=5, moment_absent=6
  type(ieee_flag_type),parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
  type :: arithmetic_guard
     logical :: flags(3),halting(3)
  end type
contains
  subroutine calc_cold_moments(species,w,kpar_d,kperp_d,a,e,b,velocity,density,ierr)
    !! Direct cold-fluid momentum/continuity solve, independent of chi.
    !! Cold temperature and anisotropy fields are never read.
    type(specie),intent(in) :: species
    complex,intent(in) :: w,e(3),b(3)
    real,intent(in) :: kpar_d,kperp_d,a
    complex,intent(out) :: velocity(3),density
    integer,intent(out) :: ierr
    complex :: delta,denominator,force(3)
    real :: gyro
    type(arithmetic_guard) :: guard
    complex,parameter :: imaginary=(0.,1.)
    velocity=invalid();density=invalid()
    ierr=validate(species,w,kpar_d,kperp_d,a,e)
    if (ierr/=moment_ok) return
    if (.not.all(finite(b))) then
       ierr=moment_invalid;return
    endif
    call begin_arithmetic(guard)
    gyro=species%mu_s/species%Q_s
    delta=w-kpar_d*species%vv_s
    ierr=moment_zero_doppler
    if (delta==cmplx(0.,0.)) goto 900
    ierr=moment_cyclotron_pole
    if (delta==gyro .or. delta==-gyro) goto 900
    force=e+a*species%vv_s*[-b(2),b(1),cmplx(0.,0.)]
    denominator=(delta-gyro)*(delta+gyro)
    velocity(1)=gyro*(imaginary*delta*force(1)-gyro*force(2))/denominator
    velocity(2)=gyro*(gyro*force(1)+imaginary*delta*force(2))/denominator
    velocity(3)=imaginary*gyro*force(3)/delta
    density=(kperp_d*velocity(1)+kpar_d*velocity(3))/(a*delta)
    ierr=moment_ok
900 call end_arithmetic(guard,velocity,density,ierr)
  end subroutine

  subroutine calc_kinetic_moments(species,w,kpar_d,kperp_d,a,e,chi,velocity,density,ierr)
    !! Existing calc_eigen current/continuity method in physical coordinates.
    !! F=v+a*u*density*z; removing the density current is essential for drift.
    !! This is algebraically its original longitudinal elimination, without
    !! intermediate removable Doppler divisions. No kinetic closure is changed.
    type(specie),intent(in) :: species
    complex,intent(in) :: w,e(3),chi(3,3)
    real,intent(in) :: kpar_d,kperp_d,a
    complex,intent(out) :: velocity(3),density
    integer,intent(out) :: ierr
    complex :: current(3)
    type(arithmetic_guard) :: guard
    velocity=invalid();density=invalid()
    ierr=validate(species,w,kpar_d,kperp_d,a,e)
    if (ierr/=moment_ok) return
    if (.not.all(finite(chi))) then
       ierr=moment_invalid;return
    endif
    call begin_arithmetic(guard)
    current=-cmplx(0.,1.)*w*(species%Q_s/species%D_s)*a**2*matmul(chi,e)
    density=(kperp_d*current(1)+kpar_d*current(3))/(a*w)
    velocity=current
    velocity(3)=current(3)-a*species%vv_s*density
    ierr=moment_ok
    call end_arithmetic(guard,velocity,density,ierr)
  end subroutine

  subroutine calc_species_moments(species,w,kpar_d,kperp_d,a,e,b,chi,velocity,density,status)
    !! Explicit state; original species ordering and per-species failure isolation.
    type(specie),intent(in) :: species(:)
    complex,intent(in) :: w,e(3),b(3),chi(:,:,:)
    real,intent(in) :: kpar_d,kperp_d,a
    complex,intent(out) :: velocity(:,:),density(:)
    integer,intent(out) :: status(:)
    integer :: s,n
    velocity=invalid();density=invalid();status=moment_invalid
    n=size(species)
    if (size(velocity,1)/=3 .or. size(velocity,2)/=n .or. size(density)/=n .or. size(status)/=n) return
    if (size(chi,1)/=n .or. size(chi,2)/=3 .or. size(chi,3)/=3) return
    if (.not.all(finite(e)) .or. .not.all(finite(b))) return
    do s=1,n
       if (.not.ieee_is_finite(species(s)%tau_s)) cycle
       if (species(s)%tau_s==0.) then
          call calc_cold_moments(species(s),w,kpar_d,kperp_d,a,e,b,velocity(:,s),density(s),status(s))
       else if (species(s)%tau_s>0.) then
          call calc_kinetic_moments(species(s),w,kpar_d,kperp_d,a,e,chi(s,:,:),velocity(:,s),density(s),status(s))
       endif
    enddo
  end subroutine

  integer function validate(species,w,kpar,kperp,a,e) result(ierr)
    type(specie),intent(in) :: species
    complex,intent(in) :: w,e(3)
    real,intent(in) :: kpar,kperp,a
    ierr=moment_invalid
    if (.not.all(ieee_is_finite([kpar,kperp,a,species%mu_s,species%Q_s,species%D_s,species%vv_s]))) return
    if (.not.finite(w) .or. .not.all(finite(e))) return
    if (a<=0. .or. species%mu_s<=0. .or. species%Q_s==0. .or. species%D_s<0.) return
    ierr=moment_absent
    if (species%D_s==0.) return
    ierr=moment_zero_frequency
    if (w==cmplx(0.,0.)) return
    ierr=moment_ok
  end function

  subroutine begin_arithmetic(guard)
    type(arithmetic_guard),intent(out) :: guard
    call ieee_get_flag(exceptions,guard%flags)
    call ieee_get_halting_mode(exceptions,guard%halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
  end subroutine

  subroutine end_arithmetic(guard,velocity,density,ierr)
    type(arithmetic_guard),intent(in) :: guard
    complex,intent(inout) :: velocity(3),density
    integer,intent(inout) :: ierr
    logical :: raised(3)
    call ieee_get_flag(exceptions,raised)
    if (any(raised)) ierr=moment_numerical
    if (ierr==moment_ok) then
       if (.not.all(finite(velocity)) .or. .not.finite(density)) ierr=moment_numerical
    endif
    if (ierr/=moment_ok) then
       velocity=invalid();density=invalid()
    endif
    call ieee_set_flag(exceptions,guard%flags)
    call ieee_set_halting_mode(exceptions,guard%halting)
  end subroutine

  elemental logical function finite(z)
    complex,intent(in) :: z
    finite=ieee_is_finite(real(z)).and.ieee_is_finite(aimag(z))
  end function

  function invalid() result(z)
    complex :: z
    real :: nan
    nan=ieee_value(0.,ieee_quiet_nan);z=cmplx(nan,nan)
  end function
end module species_moments
