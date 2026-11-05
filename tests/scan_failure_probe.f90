program scan_failure_probe
  use vars, only: nspec,susc,lam,kperp,kpar
  use scan_support
  use, intrinsic :: ieee_arithmetic
  implicit none
  complex :: seed,w,x1,x2
  integer :: flag,attempt
  nspec=1;allocate(susc(1,3,3));kpar=.1;kperp=.1
  seed=(1.,.1)
  do attempt=1,5
     x1=seed*(1.-1.e-7);x2=seed*(1.+1.e-7)
     call refine_scan_point(x1,x2,1.e-13,w,flag,controlled_solver,controlled_dispersion)
     if (flag>0) seed=w
     if (attempt>=2 .and. attempt<=4) then
        if (flag>=0 .or. ieee_is_finite(real(w))) error stop 'Failed point accepted'
        if (seed/=cmplx(1.,.1)) error stop 'Failed point contaminated seed'
     else
        if (flag<=0 .or. w/=cmplx(1.,.1)) error stop 'Valid point did not recover'
     endif
     write(*,'(a,i0,a,i0)') 'ATTEMPT ',attempt,' FLAG ',flag
  enddo
contains
  complex function controlled_solver(func,a,b,tol,iflag) result(w)
    complex :: a,b
    real :: tol
    integer :: iflag
    interface
       complex function func(w)
         complex :: w
       end function
    end interface
    if (attempt==2) then
       w=cmplx(ieee_value(0.,ieee_quiet_nan),0.);iflag=-1
    else
       w=(1.,.1);iflag=3
       if (attempt==3) w=(1.1,.1) ! Finite, but not a root.
       if (attempt==4) iflag=75 ! Legacy iteration exhaustion must not converge.
    endif
  end function
  complex function controlled_dispersion(w) result(d)
    use vars, only: dispersion_status
    complex :: w
    dispersion_status=0;susc=0.;lam=0.
    lam(1,1)=1.;lam(2,2)=1.;lam(3,3)=w-cmplx(1.,.1)
    d=lam(3,3)
  end function
end program
