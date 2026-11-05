module electromagnetic_fields
  !! Existing PLUME cofactor polarization and Faraday fields, isolated from
  !! kinetic moments/heating and from the reference thermal normalization.
  use vars, only: lam_state=>lam,spec,kpar_state=>kpar,kperp_state=>kperp, &
       vtp,normalization,dispersion_status
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private
  public :: calc_em_fields,calc_em_eigen,calc_mixed_eigen
  public :: em_ok,em_invalid,em_degenerate,em_gauge,em_not_mode,em_numerical
  integer,parameter :: em_ok=0,em_invalid=1,em_degenerate=2,em_gauge=3,em_not_mode=4,em_numerical=5
  type(ieee_flag_type),parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
  abstract interface
     complex function dispersion_function(w)
       complex :: w
     end function
  end interface
contains
  subroutine calc_em_fields(wave,omega,kpar,kperp,velocity,mode,alpha_ref,electric,magnetic,ierr)
    complex,intent(in) :: wave(3,3),omega
    real,intent(in) :: kpar,kperp,velocity,alpha_ref
    character(*),intent(in) :: mode
    complex,intent(out) :: electric(3),magnetic(3)
    integer,intent(out) :: ierr
    complex :: lam(3,3),cofactors(3,3)
    real :: scale,row_scale,speed,rank_tol,largest,residual
    integer :: row,best
    logical :: flags(3),halting(3),raised(3)
    electric=invalid_value();magnetic=electric;ierr=em_invalid
    if (.not.finite(omega)) return
    if (omega==cmplx(0.,0.)) return
    if (.not.all(ieee_is_finite(real(wave))) .or. .not.all(ieee_is_finite(aimag(wave)))) return
    if (.not.ieee_is_finite(kpar) .or. .not.ieee_is_finite(kperp)) return
    if (.not.ieee_is_finite(velocity)) return
    if (velocity<=0.) return
    if (mode/='rho_ref' .and. mode/='d_ref') return
    if (mode=='rho_ref') then
       if (.not.ieee_is_finite(alpha_ref)) return
       if (alpha_ref<=0.) return
    endif
    call ieee_get_flag(exceptions,flags)
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    ierr=em_numerical
    speed=velocity
    if (mode=='rho_ref') speed=velocity*sqrt(alpha_ref)
    if (.not.ieee_is_finite(speed) .or. speed<=0.) goto 900
    scale=maxval(abs(wave))
    if (.not.ieee_is_finite(scale)) goto 900
    ierr=em_degenerate
    if (scale==0.) goto 900
    lam=wave/scale
    ! Row equilibration preserves the nullspace and avoids a large longitudinal
    ! element masking a poor transverse solution. No tensor is modified.
    do row=1,3
       row_scale=maxval(abs(lam(row,:)))
       if (row_scale>0.) lam(row,:)=lam(row,:)/row_scale
    enddo
    cofactors(:,1)=cross(lam(2,:),lam(3,:))
    cofactors(:,2)=cross(lam(3,:),lam(1,:))
    cofactors(:,3)=cross(lam(1,:),lam(2,:))
    rank_tol=64.*epsilon(1.)
    largest=maxval(abs(cofactors))
    if (largest<=rank_tol) goto 900
    ierr=em_gauge
    if (maxval(abs(cofactors(1,:)))<=rank_tol*largest) goto 900
    best=maxloc(abs(cofactors(1,:)),dim=1)
    residual=huge(1.)
    if (abs(cofactors(1,1))>rank_tol*largest .and. abs(lam(3,2))>rank_tol) then
       include 'electromagnetic_electric.inc'
       if (all(ieee_is_finite(real(electric))) .and. all(ieee_is_finite(aimag(electric)))) &
            residual=row_residual(lam,electric)
    endif
    if (residual>1.e-8) then
       ! Same 2x2-minor method, with an equivalent row pair/divisor. In
       ! particular the parallel limit no longer divides by lam(3,2)=0.
       electric=cofactors(:,best)/cofactors(1,best)
       residual=row_residual(lam,electric)
    endif
    ierr=em_not_mode
    if (.not.ieee_is_finite(residual) .or. residual>1.e-8) goto 900
    magnetic(1)=-kpar*electric(2)/(omega*speed)
    magnetic(2)=-(kperp*electric(3)-kpar*electric(1))/(omega*speed)
    magnetic(3)=kperp*electric(2)/(omega*speed)
    ierr=em_ok
900 continue
    call ieee_get_flag(exceptions,raised)
    if (any(raised)) ierr=em_numerical
    if (.not.all(ieee_is_finite(real(electric))) .or. .not.all(ieee_is_finite(aimag(electric))) .or. &
         .not.all(ieee_is_finite(real(magnetic))) .or. .not.all(ieee_is_finite(aimag(magnetic)))) then
       if (ierr==em_ok) ierr=em_numerical
    endif
    if (ierr/=em_ok) then
       electric=invalid_value();magnetic=electric
    endif
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine

  subroutine calc_em_eigen(omega,electric,magnetic,ierr,evaluate)
    !! Fresh complete production matrix; leaves species susceptibilities at omega.
    complex,intent(in) :: omega
    complex,intent(out) :: electric(3),magnetic(3)
    integer,intent(out) :: ierr
    procedure(dispersion_function) :: evaluate
    complex :: determinant,w
    logical :: flags(3),halting(3),raised(3)
    electric=invalid_value();magnetic=electric;ierr=em_invalid
    if (.not.allocated(spec)) return
    if (size(spec)==0 .or. .not.finite(omega)) return
    if (omega==cmplx(0.,0.)) return
    call ieee_get_flag(exceptions,flags)
    call ieee_get_halting_mode(exceptions,halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
    w=omega;determinant=evaluate(w)
    call ieee_get_flag(exceptions,raised)
    if (dispersion_status==0 .and. finite(determinant) .and. .not.any(raised)) then
       call calc_em_fields(lam_state,omega,kpar_state,kperp_state,vtp,normalization, &
            spec(1)%alph_s,electric,magnetic,ierr)
    endif
    call ieee_set_flag(exceptions,flags)
    call ieee_set_halting_mode(exceptions,halting)
  end subroutine

  subroutine calc_mixed_eigen(omega,electric,magnetic,velocity,density,heating,split,heat_requested,ierr,evaluate)
    complex,intent(in) :: omega
    complex,intent(out) :: electric(3),magnetic(3),velocity(:,:),density(:)
    real,intent(out) :: heating(:),split(:,:)
    logical,intent(in) :: heat_requested
    integer,intent(out) :: ierr
    procedure(dispersion_function) :: evaluate
    if (heat_requested) error stop 'ERROR: mixed heating is unavailable until Stage 4B.3.'
    velocity=invalid_value();density=invalid_value()
    heating=real(invalid_value());split=real(invalid_value())
    call calc_em_eigen(omega,electric,magnetic,ierr,evaluate)
    if (ierr/=em_ok) write(*,'(a,i0,a,2es15.6)') &
         'WARNING: electromagnetic eigenfunction unavailable; status=',ierr,'; kperp,kpar=',kperp_state,kpar_state
  end subroutine

  pure function cross(a,b) result(c)
    complex,intent(in) :: a(3),b(3)
    complex :: c(3)
    c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function

  real function row_residual(matrix,e) result(residual)
    complex,intent(in) :: matrix(3,3),e(3)
    real :: row_norm,vector_norm
    integer :: i
    residual=0.;vector_norm=sqrt(sum(abs(e)**2))
    do i=1,3
       row_norm=sqrt(sum(abs(matrix(i,:))**2))
       if (row_norm==0.) cycle
       residual=max(residual,abs(sum(matrix(i,:)*e))/(row_norm*vector_norm))
    enddo
  end function

  logical function finite(z)
    complex,intent(in) :: z
    finite=ieee_is_finite(real(z)).and.ieee_is_finite(aimag(z))
  end function

  complex function invalid_value() result(z)
    real :: nan
    nan=ieee_value(0.,ieee_quiet_nan);z=cmplx(nan,nan)
  end function
end module electromagnetic_fields
