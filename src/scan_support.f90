module scan_support
  !! Stage 4A geometry/validation and failure guards, not a scan or root solver.
  use vars, only: scanner, scan, nscan, option, kperp, kpar, has_cold_species, &
       dispersion_status, susc, susc_low, lam
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  use vars, only: heating_convention,collision_type
  implicit none
  private
  public :: scan_frame, wavevector_scan, validate_scan_inputs, prepare_scan_start
  public :: initialize_scan_frame, scan_coordinates, check_scan_kpar, refine_scan_point
  public :: warn_unsupported_parallel
  logical, save :: warned_parallel=.false.
  type :: scan_frame
     real :: perpendicular, parallel, angle, magnitude
  end type
  abstract interface
     complex function dispersion_function(w)
       complex :: w
     end function
     complex function secant_function(func,x1,x2,xacc,iflag)
       import dispersion_function
       procedure(dispersion_function) :: func
       complex :: x1,x2
       real :: xacc
       integer :: iflag
     end function
  end interface
contains
  subroutine warn_unsupported_parallel(kp)
    real,intent(in) :: kp
    if (kp/=0. .or. has_cold_species() .or. warned_parallel) return
    write(*,'(a)') 'WARNING: exact-parallel all-kinetic calculations (kperp=0) are unsupported; '// &
         'the historical kinetic evaluator may return invalid results. No wavevector floor is applied.'
    warned_parallel=.true.
  end subroutine

  logical function wavevector_scan(s)
    type(scanner),intent(in) :: s
    wavevector_scan=(s%style_s==0 .and. s%type_s>=0 .and. s%type_s<=1) .or. &
         (s%style_s==-1 .and. s%type_s>=0 .and. s%type_s<=2)
  end function

  subroutine check_scan_kpar()
    call warn_unsupported_parallel(kperp)
    if (.not.ieee_is_finite(kpar) .or. .not.ieee_is_finite(kperp)) &
         error stop 'ERROR: wavevector scan contains nonfinite coordinates.'
    if (kpar==0.) error stop 'ERROR: k_parallel=0 scan points: exact perpendicular scans are outside Stage 4A.'
  end subroutine

  subroutine initialize_scan_frame(frame,kp,kz)
    type(scan_frame),intent(out) :: frame
    real,intent(in) :: kp,kz
    if (kz==0.) error stop 'ERROR: k_parallel=0: exact perpendicular scans are outside Stage 4A.'
    frame%perpendicular=kp; frame%parallel=kz
    frame%angle=atan(kp/kz)
    frame%magnitude=kz/cos(frame%angle) ! Historical signed magnitude coordinate.
  end subroutine

  real function interpolate(first,last,index,total,logarithmic,exclude_zero) result(value)
    real,intent(in) :: first,last
    integer,intent(in) :: index,total
    logical,intent(in) :: logarithmic
    logical,intent(in),optional :: exclude_zero
    if (logarithmic .and. (first<=0. .or. last<=0.)) &
         error stop 'ERROR: logarithmic wavevector scans require positive scanned coordinates.'
    if (present(exclude_zero)) then
       if (exclude_zero .and. .not.logarithmic) then
          ! Check the linear numerator before division: e.g. -.1 -> .2 in
          ! three steps can otherwise round its zero point to 1.e-17.
          if (first*real(total-index)+last*real(index)==0.) &
               error stop 'ERROR: k_parallel=0 scan points: exact perpendicular scans are outside Stage 4A.'
       endif
    endif
    if (index==0) then
       value=first
    else if (index==total) then
       value=last
    else if (logarithmic) then
       value=10.**(log10(first)+(log10(last)-log10(first))*real(index)/real(total))
    else
       value=first+(last-first)*real(index)/real(total)
    endif
  end function

  subroutine apply_coordinate(frame,s,index,kp,kz,angle,magnitude,polar)
    type(scan_frame),intent(in) :: frame
    type(scanner),intent(in) :: s
    integer,intent(in) :: index
    real,intent(inout) :: kp,kz,angle,magnitude
    logical,intent(inout) :: polar
    real :: degrees,pi
    integer :: total
    pi=4.*atan(1.)
    total=s%n_scan*s%n_res
    if (s%style_s==0) then
       if (s%type_s==0) kp=interpolate(s%range_i,s%range_f,index,total,s%log_scan)
       if (s%type_s==1) kz=interpolate(s%range_i,s%range_f,index,total,s%log_scan,.true.)
    else
       select case(s%type_s)
       case(0)
          kp=interpolate(frame%perpendicular,s%range_i,index,total,s%log_scan)
          kz=interpolate(frame%parallel,s%range_f,index,total,s%log_scan,.true.)
       case(1)
          degrees=interpolate(frame%angle*180./pi,s%range_i,index,total,s%log_scan)
          ! Exact degree-space check: cos(pi/2) need not round to floating zero.
          if (mod(abs(degrees),180.)==90.) &
               error stop 'ERROR: k_parallel=0 angle: exact perpendicular scans are outside Stage 4A.'
          angle=degrees*pi/180.;polar=.true.
       case(2)
          magnitude=interpolate(frame%magnitude,s%range_f,index,total,s%log_scan,.true.)
          polar=.true.
       end select
    endif
  end subroutine

  subroutine scan_coordinates(frame,first,i,kp,kz,second,j)
    type(scan_frame),intent(in) :: frame
    type(scanner),intent(in) :: first
    integer,intent(in) :: i
    real,intent(out) :: kp,kz
    type(scanner),intent(in),optional :: second
    integer,intent(in),optional :: j
    real :: angle,magnitude
    logical :: polar
    kp=frame%perpendicular;kz=frame%parallel
    angle=frame%angle;magnitude=frame%magnitude;polar=.false.
    call apply_coordinate(frame,first,i,kp,kz,angle,magnitude,polar)
    if (present(second)) call apply_coordinate(frame,second,j,kp,kz,angle,magnitude,polar)
    if (polar) then
       kp=magnitude*sin(angle);kz=magnitude*cos(angle)
    endif
    if (kz==0.) error stop 'ERROR: k_parallel=0 scan points: exact perpendicular scans are outside Stage 4A.'
    if (.not.ieee_is_finite(kp) .or. .not.ieee_is_finite(kz)) &
         error stop 'ERROR: wavevector scan produced nonfinite coordinates.'
    call warn_unsupported_parallel(kp)
  end subroutine

  subroutine validate_scan_inputs()
    type(scan_frame) :: frame
    real :: kp,kz
    integer :: s,i,j,steps
    logical :: both_wave,cartesian,polar
    kp=kperp;kz=kpar
    do s=1,nscan
       if (option==4) then
          if (scan(s)%eigen_s .or. scan(s)%tensor_s) write(*,'(a,i0,a)') &
               'WARNING: frequency-map scan ',s,': eigen/tensor requests are ignored; map files contain no such diagnostics.'
          if (scan(s)%heat_s .and. heating_convention=='legacy') write(*,'(a,i0,a)') &
               'WARNING: frequency-map scan ',s,': legacy heating request is ignored; map files contain no heating.'
       endif
       if (has_cold_species()) then
          if (.not.wavevector_scan(scan(s))) &
               error stop 'ERROR: Stage 4A cold-containing scans support wavevectors only; plasma-parameter scans are deferred.'
          if (scan(s)%heat_s .and. heating_convention=='legacy') &
               error stop 'ERROR: legacy cold-containing scan heating is unavailable; Stage 4B.3 corrected mode is required.'
       endif
       if (scan(s)%heat_s .and. heating_convention=='corrected') then
          if (collision_type/=0) &
               error stop 'ERROR: corrected heating requires collision_type=0; use legacy for collisional diagnostics.'
          if (option==4) error stop 'ERROR: corrected heating is not supported in frequency maps; map formats remain unchanged.'
          if (.not.wavevector_scan(scan(s))) error stop 'ERROR: corrected plasma-parameter heating scans are deferred.'
       endif
       if (.not.wavevector_scan(scan(s))) cycle
       if (.not.ieee_is_finite(scan(s)%range_i) .or. .not.ieee_is_finite(scan(s)%range_f)) &
            error stop 'ERROR: wavevector scan bounds must be finite.'
       if (scan(s)%n_scan<=0 .or. scan(s)%n_res<=0) &
            error stop 'ERROR: wavevector scans require positive ns and nres.'
       if (scan(s)%n_scan>huge(steps)/scan(s)%n_res) &
            error stop 'ERROR: wavevector scan step count exceeds integer range.'
    enddo
    if (option==2) then
       if (nscan/=2) error stop 'ERROR: two-dimensional scans require nscan=2.'
       both_wave=wavevector_scan(scan(1)).and.wavevector_scan(scan(2))
       if (both_wave) then
          cartesian=all(scan%style_s==0) .and. scan(1)%type_s/=scan(2)%type_s
          polar=all(scan%style_s==-1) .and. scan(1)%type_s==1 .and. scan(2)%type_s==2
          if (.not.cartesian .and. .not.polar) &
               error stop 'ERROR: wavevector grids require distinct kperp/kpar axes or ordered (theta,magnitude) axes.'
          call initialize_scan_frame(frame,kp,kz)
          ! The supported axes are independent; no O(N*M) preflight is needed.
          do i=0,scan(1)%n_scan*scan(1)%n_res
             call scan_coordinates(frame,scan(1),i,kp,kz,scan(2),0)
          enddo
          do j=0,scan(2)%n_scan*scan(2)%n_res
             call scan_coordinates(frame,scan(1),0,kp,kz,scan(2),j)
          enddo
          return
       endif
    endif
    do s=1,nscan
       if (.not.wavevector_scan(scan(s))) cycle
       call initialize_scan_frame(frame,kp,kz)
       do i=0,scan(s)%n_scan*scan(s)%n_res
          call scan_coordinates(frame,scan(s),i,kp,kz)
       enddo
    enddo
  end subroutine

  subroutine prepare_scan_start()
    type(scan_frame) :: frame
    ! Only the newly enabled compositions change initial-state setup. Keep the
    ! historical all-kinetic initialization and branch-continuation path intact.
    if (.not.has_cold_species() .or. nscan==0) return
    call initialize_scan_frame(frame,kperp,kpar)
    if (option==2) then
       call scan_coordinates(frame,scan(1),0,kperp,kpar,scan(2),0)
    else
       call scan_coordinates(frame,scan(1),0,kperp,kpar)
    endif
  end subroutine

  subroutine refine_scan_point(x1,x2,tolerance,omega,iflag,solver,evaluate)
    complex :: x1,x2,omega
    real :: tolerance
    integer,intent(out) :: iflag
    procedure(secant_function) :: solver
    procedure(dispersion_function) :: evaluate
    complex :: determinant,nan_value
    real :: scale,nan
    logical :: flags(3),halting(3),raised(3),valid
    type(ieee_flag_type),parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
    call ieee_get_flag(exceptions,flags)
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    omega=solver(evaluate,x1,x2,tolerance,iflag)
    valid=iflag>0 .and. iflag<75 .and. ieee_is_finite(real(omega)) .and. ieee_is_finite(aimag(omega))
    if (valid) then
       call ieee_set_flag(exceptions,.false.)
       determinant=evaluate(omega)
       call ieee_get_flag(exceptions,raised)
       valid=dispersion_status==0 .and. .not.any(raised)
       valid=valid .and. ieee_is_finite(real(determinant)) .and. ieee_is_finite(aimag(determinant))
       if (valid) valid=all(ieee_is_finite(real(lam))) .and. all(ieee_is_finite(aimag(lam)))
       if (valid) then
          scale=maxval(abs(lam))
          valid=ieee_is_finite(scale) .and. scale>0.
          if (valid) valid=abs(((determinant/scale)/scale)/scale)<1.e-8
       endif
       if (.not.valid) iflag=-3
    endif
    if (.not.valid) then
       if (iflag>=0) iflag=-2
       nan=ieee_value(0.,ieee_quiet_nan);nan_value=cmplx(nan,nan)
       omega=nan_value;susc=nan_value;lam=nan_value
       if (allocated(susc_low)) susc_low=nan_value
       if (dispersion_status==0) dispersion_status=5
       write(*,'(a,2es15.6,a,i0,a)') 'WARNING: scan root failed at (kperp,kpar)=', &
            kperp,kpar,'; flag=',iflag,'; output is NaN, retaining last valid continuation seed.'
    endif
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine
end module scan_support
