module harmonic_selection
  !! Positive Chernoff bounds for the ordinary and n**2-weighted Skellam tails.
  !! No subtraction from unity, Bessel approximation, or thermal normalization.
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: select_harmonics, harmonic_tail_bounds
contains
  subroutine harmonic_tail_bounds(lambda,n,e0,e2)
    real, intent(in) :: lambda
    integer, intent(in) :: n
    real, intent(out) :: e0,e2
    real :: m,h,t,ratio,exponent0,exponent2
    if (lambda==0.) then
       e0=0.; e2=0.
       return
    endif
    m=real(n)+1.
    ! hypot(lambda,m) without squaring large lambda; asinh(m/lambda)
    ! without overflowing the ratio for a subnormal positive lambda.
    h=max(lambda,m)*sqrt(1.+(min(lambda,m)/max(lambda,m))**2)
    ratio=log(m)-log(lambda)
    if (ratio>20.) then
       t=log(2.)+ratio
    else
       t=asinh(m/lambda)
    endif
    exponent0=log(2.)+(m/h)*m/(1.+lambda/h)-m*t
    exponent2=exponent0+log(h+m*m)-log(lambda)
    ! Bounds can exceed one. Underflow is harmless and not a convergence floor.
    e0=exp(min(0.,exponent0)); e2=exp(min(0.,exponent2))
  end subroutine

  subroutine select_harmonics(lambda,tolerance,maximum,n,e0,e2,converged,ierr)
    real, intent(in) :: lambda,tolerance
    integer, intent(in) :: maximum
    integer, intent(out) :: n,ierr
    real, intent(out) :: e0,e2
    logical, intent(out) :: converged
    integer :: lower,upper,middle
    ierr=1; n=0; e0=1.; e2=1.; converged=.false.
    if (.not. ieee_is_finite(lambda) .or. .not. ieee_is_finite(tolerance)) return
    if (lambda<0. .or. tolerance<=0. .or. tolerance>=1. .or. maximum<0) return
    ierr=0
    call harmonic_tail_bounds(lambda,maximum,e0,e2)
    n=maximum
    if (max(e0,e2)>=tolerance) return
    ! The bound is decreasing in the convergent region. Binary search avoids
    ! work proportional to the user cap. N=0 is exact only at lambda=0.
    lower=-1; upper=maximum
    do while (upper-lower>1)
       middle=lower+(upper-lower)/2
       call harmonic_tail_bounds(lambda,middle,e0,e2)
       if (max(e0,e2)<tolerance) then
          upper=middle
       else
          lower=middle
       endif
    enddo
    n=upper
    call harmonic_tail_bounds(lambda,n,e0,e2)
    converged=.true.
  end subroutine
end module harmonic_selection
