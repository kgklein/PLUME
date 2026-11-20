module scan_support
  use vars, only: emit_warning
  !! Native scan-state application, geometry validation and guarded root refinement.
  use vars, only: scanner, scan, nscan, option, kperp, kpar, has_cold_species, &
       dispersion_status, susc, susc_low, lam
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  use vars, only: heating_convention,collision_type,normalization,betap,vtp,spec,nspec,specie,sync_solver_state
  implicit none
  private
  public :: scan_frame, wavevector_scan, validate_scan_inputs, prepare_scan_start
  public :: initialize_scan_frame, scan_coordinates, check_scan_kpar, refine_scan_point
  public :: warn_unsupported_parallel
  public :: parameter_scan_workflow, parameter_scan_state, capture_scan_state, apply_parameter_point
  logical, save :: warned_parallel=.false.
  logical, save :: warned_scan_charge=.false.,warned_scan_current=.false.,warned_scan_frame=.false.
  type :: scan_frame
     real :: perpendicular, parallel, angle, magnitude
  end type
  type :: parameter_scan_state
     real :: perpendicular,parallel,beta,velocity
     type(specie),allocatable :: populations(:)
     type(scan_frame) :: geometry
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
    block
       character(2048) :: warning_text
       write(warning_text,'(a)') 'WARNING: exact-parallel all-kinetic calculations (kperp=0) are unsupported; '// &
         'the historical kinetic evaluator may return invalid results. No wavevector floor is applied.'
       call emit_warning(trim(warning_text))
    end block
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

  logical function parameter_scan_workflow() result(enabled)
    !! New support is deliberately separate from historical/other-normalization paths.
    integer :: s
    enabled=.false.
    if (normalization/='rho_ref' .or. heating_convention/='corrected' .or. collision_type/=0) return
    if (has_cold_species() .or. (option/=1 .and. option/=2)) return
    if (.not.allocated(scan)) return
    do s=1,nscan
       if (.not.wavevector_scan(scan(s))) enabled=.true.
    enddo
  end function

  subroutine capture_scan_state(state)
    type(parameter_scan_state),intent(out) :: state
    call validate_native_parameters()
    state%perpendicular=kperp;state%parallel=kpar;state%beta=betap;state%velocity=vtp
    state%populations=spec
    call initialize_scan_frame(state%geometry,kperp,kpar)
  end subroutine

  subroutine set_scalar_axis(s,index)
    type(scanner),intent(in) :: s
    integer,intent(in) :: index
    real :: value
    value=interpolate(s%range_i,s%range_f,index,s%n_scan*s%n_res,s%log_scan)
    if (s%style_s==0) then
       select case(s%type_s)
       case(2)
          betap=value
       case(3)
          vtp=value
       end select
    else
       select case(s%type_s)
       case(0)
          spec(s%style_s)%tau_s=value
       case(1)
          spec(s%style_s)%mu_s=value
       case(2)
          spec(s%style_s)%alph_s=value
       case(3)
          spec(s%style_s)%Q_s=value
       case(4)
          spec(s%style_s)%D_s=value
       case(5)
          spec(s%style_s)%vv_s=value
       end select
    endif
  end subroutine

  subroutine apply_parameter_point(state,first,i,second,j,preflight)
    !! Rebuild native inputs, then refresh derived normalization. Diagnostic
    !! frequency derivatives subsequently hold this entire plasma fixed.
    type(parameter_scan_state),intent(in) :: state
    type(scanner),intent(in) :: first
    integer,intent(in) :: i
    type(scanner),intent(in),optional :: second
    integer,intent(in),optional :: j
    logical,intent(in),optional :: preflight
    real :: kp,kz,angle,magnitude
    logical :: polar,checking
    checking=.false.
    if (present(preflight)) checking=preflight
    betap=state%beta;vtp=state%velocity;spec(:)=state%populations
    kp=state%perpendicular;kz=state%parallel
    angle=state%geometry%angle;magnitude=state%geometry%magnitude;polar=.false.
    if (wavevector_scan(first)) then
       call apply_coordinate(state%geometry,first,i,kp,kz,angle,magnitude,polar)
    else
       call set_scalar_axis(first,i)
    endif
    if (present(second)) then
       if (wavevector_scan(second)) then
          call apply_coordinate(state%geometry,second,j,kp,kz,angle,magnitude,polar)
       else
          call set_scalar_axis(second,j)
       endif
    endif
    if (polar) then
       kp=magnitude*sin(angle);kz=magnitude*cos(angle)
    endif
    kperp=kp;kpar=kz
    call check_scan_kpar()
    call validate_native_parameters()
    call sync_solver_state()
    if (.not.checking) call warn_parameter_equilibrium()
  end subroutine

  subroutine validate_native_parameters()
    integer :: s
    if (.not.all(ieee_is_finite([betap,vtp,kperp,kpar]))) &
         error stop 'ERROR: parameter scan plasma inputs must be finite.'
    if (betap<=0. .or. vtp<=0.) error stop 'ERROR: parameter scan beta and vtp must be positive.'
    do s=1,nspec
       if (.not.all(ieee_is_finite([spec(s)%tau_s,spec(s)%mu_s,spec(s)%alph_s, &
            spec(s)%Q_s,spec(s)%D_s,spec(s)%vv_s]))) &
            error stop 'ERROR: parameter scan species inputs must be finite.'
       if (min(spec(s)%tau_s,spec(s)%mu_s,spec(s)%alph_s)<=0.) &
            error stop 'ERROR: rho_ref parameter scans require warm species with positive mass and anisotropy.'
       if (spec(s)%Q_s==0. .or. spec(s)%D_s<0.) &
            error stop 'ERROR: parameter scan charge must be nonzero and density nonnegative.'
    enddo
  end subroutine

  subroutine warn_parameter_equilibrium()
    !! A changing charge/density/drift does not imply a compensating change to
    !! another population. Warn once for each condition, without modifying inputs.
    real :: charge,current
    character(200) :: message
    charge=sum(spec%D_s/spec%Q_s);current=sum(spec%D_s*spec%vv_s/spec%Q_s)
    if (abs(charge)>=1.e-6 .and. .not.warned_scan_charge) then
       write(message,'(a,es24.16)') 'WARNING: parameter scan is not charge neutral; Sum [D_s/Q_s] = ',charge
       call emit_warning(trim(message));warned_scan_charge=.true.
    endif
    if (abs(current)>=1.e-6 .and. .not.warned_scan_current) then
       write(message,'(a,es24.16)') 'WARNING: parameter scan has equilibrium current; Sum [D_s*vv_s/Q_s] = ',current
       call emit_warning(trim(message));warned_scan_current=.true.
    endif
    if (spec(1)%vv_s/=0. .and. .not.warned_scan_frame) then
       call emit_warning('WARNING: parameter scan leaves the reference-species rest frame; no frame change is applied.')
       warned_scan_frame=.true.
    endif
  end subroutine

  subroutine validate_parameter_workflow()
    type(parameter_scan_state) :: initial,segment
    integer :: s,i,steps
    logical :: positive,nonnegative,charge
    call capture_scan_state(initial)
    segment=initial
    do s=1,nscan
       if (scan(s)%style_s < -1 .or. scan(s)%style_s > nspec) &
            error stop 'ERROR: invalid parameter scan species/style index.'
       if (scan(s)%type_s<0) error stop 'ERROR: invalid parameter scan type.'
       if (scan(s)%style_s==-1 .and. scan(s)%type_s>2) error stop 'ERROR: invalid geometry scan type.'
       if (scan(s)%style_s==0 .and. scan(s)%type_s>3) error stop 'ERROR: invalid global scan type.'
       if (scan(s)%style_s>0 .and. scan(s)%type_s>5) error stop 'ERROR: invalid species scan type.'
       if (.not.all(ieee_is_finite([scan(s)%range_i,scan(s)%range_f]))) &
            error stop 'ERROR: parameter scan bounds must be finite.'
       if (scan(s)%n_scan<=0 .or. scan(s)%n_res<=0) &
            error stop 'ERROR: parameter scans require positive ns and nres.'
       if (scan(s)%n_scan>(huge(steps)-1)/scan(s)%n_res) &
            error stop 'ERROR: parameter scan step count exceeds integer range.'
       if (wavevector_scan(scan(s))) cycle
       positive=scan(s)%style_s==0
       nonnegative=.false.;charge=.false.
       if (scan(s)%style_s>0) then
          positive=scan(s)%type_s<=2
          charge=scan(s)%type_s==3
          nonnegative=scan(s)%type_s==4
       endif
       if (positive .and. min(scan(s)%range_i,scan(s)%range_f)<=0.) &
            error stop 'ERROR: beta, vtp, temperature, mass and anisotropy scan bounds must be positive.'
       if (nonnegative .and. min(scan(s)%range_i,scan(s)%range_f)<0.) &
            error stop 'ERROR: density scan bounds must be nonnegative.'
       if (charge) then
          if (scan(s)%range_i==0. .or. scan(s)%range_f==0.) error stop 'ERROR: charge scans cannot reach zero.'
          if ((scan(s)%range_i>0.) .neqv. (scan(s)%range_f>0.)) &
               error stop 'ERROR: charge scans cannot cross zero.'
       endif
       if (scan(s)%log_scan .and. min(scan(s)%range_i,scan(s)%range_f)<=0.) &
            error stop 'ERROR: logarithmic parameter scan bounds must be positive.'
    enddo
    if (option==2) then
       if (scan(1)%style_s==scan(2)%style_s .and. scan(1)%type_s==scan(2)%type_s) &
            error stop 'ERROR: two-dimensional scan axes target the same parameter.'
       if (wavevector_scan(scan(1)).and.wavevector_scan(scan(2))) &
            error stop 'ERROR: overlapping wavevector axes in parameter workflow.'
       ! Independent axes need O(N+M) preflight, not a full Cartesian traversal.
       do i=0,scan(1)%n_scan*scan(1)%n_res
          call apply_parameter_point(initial,scan(1),i,scan(2),0,preflight=.true.)
       enddo
       do i=0,scan(2)%n_scan*scan(2)%n_res
          call apply_parameter_point(initial,scan(1),0,scan(2),i,preflight=.true.)
       enddo
    else
       do s=1,nscan
          do i=0,scan(s)%n_scan*scan(s)%n_res
             call apply_parameter_point(segment,scan(s),i,preflight=.true.)
          enddo
          call capture_scan_state(segment)
       enddo
    endif
    ! Validation must not move the input plasma before initial root selection.
    kperp=initial%perpendicular;kpar=initial%parallel;betap=initial%beta;vtp=initial%velocity
    spec(:)=initial%populations
    call sync_solver_state()
  end subroutine

  subroutine validate_scan_inputs()
    type(scan_frame) :: frame
    real :: kp,kz
    integer :: s,i,j,steps
    logical :: both_wave,cartesian,polar,requested_heat
    if (option==2 .and. nscan/=2) error stop 'ERROR: two-dimensional scans require nscan=2.'
    if (parameter_scan_workflow()) then
       call validate_parameter_workflow()
       return
    endif
    kp=kperp;kz=kpar
    do s=1,nscan
       requested_heat=scan(s)%heat_s
       if (option==2) requested_heat=scan(2)%heat_s
       if (option==4) then
          if (scan(s)%eigen_s .or. scan(s)%tensor_s) then
             block
                character(2048) :: warning_text
                write(warning_text,'(a,i0,a)') &
               'WARNING: frequency-map scan ',s,': eigen/tensor requests are ignored; map files contain no such diagnostics.'
                call emit_warning(trim(warning_text))
             end block
          end if
          if (scan(s)%heat_s .and. heating_convention=='legacy') then
             block
                character(2048) :: warning_text
                write(warning_text,'(a,i0,a)') &
               'WARNING: frequency-map scan ',s,': legacy heating request is ignored; map files contain no heating.'
                call emit_warning(trim(warning_text))
             end block
          end if
       endif
       if (has_cold_species()) then
          if (.not.wavevector_scan(scan(s))) &
               error stop 'ERROR: Stage 4A cold-containing scans support wavevectors only; plasma-parameter scans are deferred.'
          if (scan(s)%heat_s .and. heating_convention=='legacy') &
               error stop 'ERROR: legacy cold-containing scan heating is unavailable; Stage 4B.3 corrected mode is required.'
       endif
       if (requested_heat .and. heating_convention=='corrected') then
          if (collision_type/=0) &
               error stop 'ERROR: corrected heating requires collision_type=0; use legacy for collisional diagnostics.'
          if (option==4) error stop 'ERROR: corrected heating is not supported in frequency maps; map formats remain unchanged.'
          if (.not.wavevector_scan(scan(s))) &
               error stop 'ERROR: corrected plasma-parameter heating requires collisionless all-kinetic rho_ref scans.'
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
    type(parameter_scan_state) :: state
    if (parameter_scan_workflow()) then
       call capture_scan_state(state)
       if (option==2) then
          call apply_parameter_point(state,scan(1),0,scan(2),0)
       else
          call apply_parameter_point(state,scan(1),0)
       endif
       return
    endif
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
       block
          character(2048) :: warning_text
          write(warning_text,'(a,2es15.6,a,i0,a)') 'WARNING: scan root failed at (kperp,kpar)=', &
            kperp,kpar,'; flag=',iflag,'; output is NaN, retaining last valid continuation seed.'
          call emit_warning(trim(warning_text))
       end block
    endif
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine
end module scan_support
