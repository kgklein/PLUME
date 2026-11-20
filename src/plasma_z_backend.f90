module plasma_z_backend
  !! Backend adapter and accurate-path evaluation guards. Legacy math is untouched.
  use vars, only: z_backend,dispersion_status,susc,lam,susc_low
  use plasma_z_legacy, only: zetout_legacy
  use plasma_z_faddeeva
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private
  public :: zetout
  public :: z_evaluation_guard,begin_z_evaluation,end_z_evaluation
  type :: z_evaluation_guard
     logical :: active=.false.
     logical :: flags(3),halting(3)
  end type
  type(ieee_flag_type),parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
contains
  pure complex function zetout(z) result(value)
    complex,intent(in) :: z
    if (z_backend=='legacy') then
       value=zetout_legacy(z)
    else
       value=z_faddeeva_value(1.,z)
    endif
  end function
  subroutine begin_z_evaluation(omega,guard,value,reject)
    complex,intent(in) :: omega
    type(z_evaluation_guard),intent(out) :: guard
    complex,intent(inout) :: value
    logical,intent(out) :: reject
    guard%active=z_backend=='faddeeva';reject=.false.
    if (.not.guard%active) return
    call ieee_get_flag(exceptions,guard%flags)
    call ieee_get_halting_mode(exceptions,guard%halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    dispersion_status=0
    if (.not.ieee_is_finite(real(omega)) .or. .not.ieee_is_finite(aimag(omega))) then
       dispersion_status=1
       reject=.true.
    else if (omega==cmplx(0.,0.)) then
       dispersion_status=2
       reject=.true.
    endif
    if (reject) then
       value=invalid_value()
       call end_z_evaluation(guard,value)
    endif
  end subroutine

  subroutine end_z_evaluation(guard,value)
    type(z_evaluation_guard),intent(inout) :: guard
    complex,intent(inout) :: value
    logical :: raised(3)
    if (.not.guard%active) return
    call ieee_get_flag(exceptions,raised)
    if (dispersion_status==0) then
       if (any(raised) .or. .not.ieee_is_finite(real(value)) .or. &
            .not.ieee_is_finite(aimag(value))) dispersion_status=5
    endif
    if (dispersion_status/=0) then
       value=invalid_value()
       if (allocated(susc)) susc=value
       lam=value
       if (allocated(susc_low)) susc_low=value
    endif
    call ieee_set_flag(exceptions,guard%flags)
    call ieee_set_halting_mode(exceptions,guard%halting)
    guard%active=.false.
  end subroutine

  complex function invalid_value() result(value)
    real :: nan
    nan=ieee_value(0.,ieee_quiet_nan)
    value=cmplx(nan,nan)
  end function
end module plasma_z_backend
