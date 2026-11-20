module heating_response
  !! Explicit collisionless diagnostic snapshots. Never writes production
  !! susc/lam/solver_* state; shared kinetic equations and caches are reused.
  use vars, only: specie
  use cold_reference, only: calc_warm_auxiliary
  use cold_susceptibility, only: calc_chi_cold
  use species_moments, only: calc_cold_moments
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
  implicit none
  private
  public :: response_snapshot,cold_frequency_derivative,hermitian,absorptive,quadratic,diagnostic_nan
  public :: bessel_callback,z_callback,diagnostic_guard,start_guard,finish_guard
  type(ieee_flag_type),parameter :: exceptions(3)=[ieee_invalid,ieee_divide_by_zero,ieee_overflow]
  type :: diagnostic_guard
     logical :: flags(3),halting(3)
  end type
  abstract interface
     real function bessel_callback(n,x)
       integer :: n
       real :: x
     end function
     complex function z_callback(k,z)
       real :: k
       complex :: z
     end function
  end interface
contains
  subroutine response_snapshot(spec,om,kpar,kperp,vtp,susc,low,ierr,bessel,zet_in,warm_only,auxiliary_beta)
    type(specie),intent(in) :: spec(:)
    complex,intent(in) :: om
    real,intent(in) :: kpar,kperp,vtp
    complex,intent(out) :: susc(:,:,:),low(:,:,:,-1:)
    integer,intent(out) :: ierr
    procedure(bessel_callback) :: bessel
    procedure(z_callback) :: zet_in
    logical,intent(in),optional :: warm_only
    real,intent(in),optional :: auxiliary_beta
    type(diagnostic_guard) :: guard
    real :: thermal,gyro,plasma
    complex :: delta,zeta,response,term
    complex,parameter :: imaginary=(0.,1.)
    integer :: is,n
    logical :: skip_cold
    susc=cmplx(0.,0.);low=diagnostic_nan();ierr=1
    if (size(susc,1)/=size(spec) .or. size(susc,2)/=3 .or. size(susc,3)/=3) return
    if (size(low,1)/=size(spec) .or. size(low,2)/=3 .or. size(low,3)/=3 .or. size(low,4)/=3) return
    skip_cold=.false.
    if (present(warm_only)) skip_cold=warm_only
    call start_guard(guard)
    ierr=0
    do is=1,size(spec)
       if (spec(is)%tau_s/=0. .or. skip_cold) cycle
       call calc_chi_cold(spec(is),om,kpar,kperp,vtp,susc(is,:,:),ierr)
       if (ierr/=0) goto 900
    enddo
    if (.not.any(spec%tau_s>0.)) goto 900
    if (kpar==0.) then
       ierr=1;goto 900
    endif
    if (kperp/=0.) then
       call calc_warm_auxiliary(spec,om,kpar,kperp,vtp,susc,ierr,bessel,zet_in,low,auxiliary_beta)
    else
       ! Exact existing parallel scalar responses, with initialized signed
       ! harmonics. The same scalar statements are included by the solver.
       do is=1,size(spec)
          if (spec(is)%tau_s==0.) cycle
          low(is,:,:,:)=cmplx(0.,0.);susc(is,:,:)=cmplx(0.,0.)
          if (spec(is)%D_s==0.) cycle
          thermal=sqrt(spec(is)%mu_s*spec(is)%tau_s)
          gyro=spec(is)%mu_s/spec(is)%Q_s
          plasma=spec(is)%D_s*spec(is)%mu_s/(spec(is)%Q_s**2*vtp**2)
          delta=om-kpar*spec(is)%vv_s
          do n=-1,1,2
             include 'kinetic_parallel_response.inc'
             term=plasma*response/(2.*om**2)
             low(is,1,1,n)=term;low(is,2,2,n)=term
             low(is,1,2,n)=imaginary*n*term;low(is,2,1,n)=-low(is,1,2,n)
             susc(is,1,1)=susc(is,1,1)+term
             susc(is,1,2)=susc(is,1,2)+imaginary*n*term
          enddo
          susc(is,2,2)=susc(is,1,1);susc(is,2,1)=-susc(is,1,2)
          include 'kinetic_parallel_longitudinal.inc'
          low(is,3,3,0)=susc(is,3,3)
       enddo
    endif
900 if (ierr==0) then
       if (.not.all(ieee_is_finite(real(susc))) .or. .not.all(ieee_is_finite(aimag(susc)))) ierr=5
    endif
    call finish_guard(guard,ierr)
    if (ierr/=0) then
       susc=diagnostic_nan();low=diagnostic_nan()
    endif
  end subroutine

  subroutine cold_frequency_derivative(s,w,kpar,kperp,a,derivative,ierr)
    !! Differentiate the independently validated fluid momentum/continuity
    !! solution: M v'=g F'+i v; eta'=k.v'/(a delta)-eta/delta.
    !! Then d(w chi E)/dw=i D/(Q a^2)*(v'+a*u*eta' z).
    type(specie),intent(in) :: s
    real,intent(in) :: w,kpar,kperp,a
    complex,intent(out) :: derivative(3,3)
    integer,intent(out) :: ierr
    complex :: chi(3,3),e(3),b(3),v(3),dv(3),eta,deta,force_prime(3),effective(3),zero(3)
    real :: gyro,delta
    integer :: col
    type(diagnostic_guard) :: guard
    derivative=diagnostic_nan();zero=cmplx(0.,0.)
    call start_guard(guard)
    call calc_chi_cold(s,cmplx(w,0.),kpar,kperp,a,chi,ierr)
    if (ierr/=0) goto 900
    if (s%D_s==0.) then
       derivative=0.;goto 900
    endif
    gyro=s%mu_s/s%Q_s;delta=w-kpar*s%vv_s
    do col=1,3
       e=zero;e(col)=1.
       b=[-kpar*e(2),kpar*e(1)-kperp*e(3),kperp*e(2)]/(a*w)
       call calc_cold_moments(s,cmplx(w,0.),kpar,kperp,a,e,b,v,eta,ierr)
       if (ierr/=0) goto 900
       force_prime=-a*s%vv_s*[-b(2),b(1),cmplx(0.,0.)]/w
       effective=force_prime+cmplx(0.,1.)*v/gyro
       call calc_cold_moments(s,cmplx(w,0.),kpar,kperp,a,effective,zero,dv,deta,ierr)
       if (ierr/=0) goto 900
       deta=deta-eta/delta
       derivative(:,col)=cmplx(0.,1.)*(s%D_s/(s%Q_s*a**2))* &
            (dv+[cmplx(0.,0.),cmplx(0.,0.),a*s%vv_s*deta])
    enddo
    derivative=hermitian(derivative)
900 if (ierr==0) then
       if (.not.all(ieee_is_finite(real(derivative))) .or. .not.all(ieee_is_finite(aimag(derivative)))) ierr=5
    endif
    call finish_guard(guard,ierr)
    if (ierr/=0) derivative=diagnostic_nan()
  end subroutine

  pure function hermitian(matrix) result(h)
    complex,intent(in) :: matrix(3,3)
    complex :: h(3,3)
    h=(matrix+transpose(conjg(matrix)))*.5
  end function

  real function quadratic(e,matrix) result(value)
    complex,intent(in) :: e(3),matrix(3,3)
    value=real(dot_product(e,matmul(matrix,e)))
  end function

  real function absorptive(e,matrix) result(value)
    complex,intent(in) :: e(3),matrix(3,3)
    value=aimag(dot_product(e,matmul(matrix,e)))
  end function

  complex function diagnostic_nan() result(z)
    real :: nan
    nan=ieee_value(0.,ieee_quiet_nan);z=cmplx(nan,nan)
  end function

  subroutine start_guard(guard)
    type(diagnostic_guard),intent(out) :: guard
    call ieee_get_flag(exceptions,guard%flags)
    call ieee_get_halting_mode(exceptions,guard%halting)
    call ieee_set_halting_mode(exceptions,.false.)
    call ieee_set_flag(exceptions,.false.)
  end subroutine

  subroutine finish_guard(guard,ierr)
    type(diagnostic_guard),intent(in) :: guard
    integer,intent(inout) :: ierr
    logical :: raised(3)
    call ieee_get_flag(exceptions,raised)
    if (any(raised)) ierr=5
    call ieee_set_flag(exceptions,guard%flags)
    call ieee_set_halting_mode(exceptions,guard%halting)
  end subroutine
end module heating_response
