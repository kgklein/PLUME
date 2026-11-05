module heating_energy
  !! Corrected collisionless absorption and formal real-frequency wave energy.
  !! No production dispersion/root state is read or written by this module.
  use vars, only: specie
  use heating_response
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: heating_result,calculate_heating,reset_heating,finish_energy,heating_channels
  public :: heat_ok,heat_invalid,heat_no_field,heat_zero_real,heat_response,heat_derivative
  public :: heat_unresolved,heat_collision,heat_root_failed,heat_arithmetic
  public :: species_heat_ok,species_heat_cold,species_heat_absent,species_heat_failed
  public :: species_heat_derivative,species_heat_energy,species_heat_applicability
  integer,parameter :: heat_ok=0,heat_invalid=1,heat_no_field=2,heat_zero_real=3,heat_response=4
  integer,parameter :: heat_derivative=5,heat_unresolved=6,heat_collision=7,heat_root_failed=8,heat_arithmetic=9
  integer,parameter :: species_heat_ok=0,species_heat_cold=1,species_heat_absent=2,species_heat_failed=3
  integer,parameter :: species_heat_derivative=4,species_heat_energy=5,species_heat_applicability=6
  real,parameter :: derivative_tolerance=1.e-6
  integer,parameter :: derivative_levels=13
  type :: heating_result
     integer :: status=heat_invalid,field_status=-1,applicability=0,response_status=0
     real :: total,electric,magnetic,species_total,f_e,f_b,f_species,f_e_abs,uncertainty
     real :: total_power,gamma_est,ratio
     real,allocatable :: species(:),fractions(:),q(:),power(:),channels(:,:),raw_channels(:,:)
     real,allocatable :: matrix_error(:),projection_error(:),step(:)
     integer,allocatable :: species_status(:)
  end type
contains
  subroutine reset_heating(r,count,status)
    type(heating_result),intent(out) :: r
    integer,intent(in) :: count,status
    real :: nan
    nan=real(diagnostic_nan())
    r%status=status;r%field_status=-1;r%applicability=0;r%response_status=0
    r%total=nan;r%electric=nan;r%magnetic=nan;r%species_total=nan
    r%f_e=nan;r%f_b=nan;r%f_species=nan;r%f_e_abs=nan;r%uncertainty=nan
    r%total_power=nan;r%gamma_est=nan;r%ratio=nan
    allocate(r%species(count),r%fractions(count),r%q(count),r%power(count),r%channels(7,count), &
         r%raw_channels(7,count),r%matrix_error(count),r%projection_error(count),r%step(count),r%species_status(count))
    r%species=nan;r%fractions=nan;r%q=nan;r%power=nan;r%channels=nan;r%raw_channels=nan
    r%matrix_error=nan;r%projection_error=nan;r%step=nan;r%species_status=species_heat_failed
  end subroutine

  subroutine calculate_heating(spec,w,kpar,kperp,a,e,r,bessel,zet_in,auxiliary_beta)
    type(specie),intent(in) :: spec(:)
    complex,intent(in) :: w,e(3)
    real,intent(in) :: kpar,kperp,a
    type(heating_result),intent(out) :: r
    procedure(bessel_callback) :: bessel
    procedure(z_callback) :: zet_in
    real,intent(in),optional :: auxiliary_beta
    complex :: base(size(spec),3,3),low(size(spec),3,3,-1:1)
    complex :: plus(size(spec),3,3),minus(size(spec),3,3),work_low(size(spec),3,3,-1:1)
    complex :: center(3,3),prev_center(size(spec),3,3),richardson(3,3),prev_r(size(spec),3,3)
    complex :: best(size(spec),3,3),older_r(size(spec),3,3),difference(3,3),b(3)
    real :: previous_score(size(spec)),best_score(size(spec)),absolute_error(size(spec))
    real :: wr,h,relative,score,pair_score,merr,perr,value,nan,enorm,roundoff
    logical :: resolved(size(spec))
    type(diagnostic_guard) :: guard
    integer :: n,s,level,ierr,arithmetic_status
    n=size(spec);call reset_heating(r,n,heat_invalid);nan=real(diagnostic_nan())
    do s=1,n
       if (.not.ieee_is_finite(spec(s)%tau_s)) return
       if (spec(s)%tau_s==0.) r%species_status(s)=species_heat_cold
       if (spec(s)%D_s==0.) r%species_status(s)=species_heat_absent
    enddo
    if (n==0) return
    if (.not.all(ieee_is_finite([real(w),aimag(w),kpar,kperp,a,real(e),aimag(e)]))) return
    if (a<=0. .or. e(1)==cmplx(0.,0.)) return
    do s=1,n
       if (.not.all(ieee_is_finite([spec(s)%mu_s,spec(s)%Q_s,spec(s)%D_s,spec(s)%vv_s]))) return
       if (spec(s)%tau_s<0. .or. spec(s)%mu_s<=0. .or. spec(s)%Q_s==0. .or. spec(s)%D_s<0.) return
       if (spec(s)%tau_s>0.) then
          if (.not.ieee_is_finite(spec(s)%alph_s)) return
          if (spec(s)%alph_s<=0. .or. kpar==0.) return
       endif
    enddo
    r%field_status=0
    wr=real(w)
    r%status=heat_zero_real
    if (wr==0.) return
    call start_guard(guard)
    r%field_status=0;r%ratio=abs(aimag(w)/wr)
    r%electric=sum(abs(e)**2)
    b=[-kpar*e(2),kpar*e(1)-kperp*e(3),kperp*e(2)]/(a*wr)
    r%magnetic=sum(abs(b)**2);enorm=sum(abs(e)**2)
    call response_snapshot(spec,cmplx(wr,0.),kpar,kperp,a,base,low,ierr,bessel,zet_in, &
         auxiliary_beta=auxiliary_beta)
    r%status=heat_response
    if (ierr/=0) then
       r%response_status=ierr;goto 900
    endif
    resolved=.false.;best=diagnostic_nan();absolute_error=0.
    previous_score=huge(1.);best_score=huge(1.)
    do s=1,n
       if (spec(s)%D_s==0.) then
          best(s,:,:)=0.;r%species(s)=0.;resolved(s)=.true.
          r%matrix_error(s)=0.;r%projection_error(s)=0.;r%step(s)=0.
       else if (spec(s)%tau_s==0.) then
          call cold_frequency_derivative(spec(s),wr,kpar,kperp,a,best(s,:,:),ierr)
          if (ierr/=0) then
             r%response_status=ierr;goto 900
          endif
          r%species(s)=quadratic(e,best(s,:,:));resolved(s)=.true.
          r%matrix_error(s)=0.;r%projection_error(s)=0.;r%step(s)=0.
          absolute_error(s)=128.*epsilon(1.)*sqrt(sum(abs(best(s,:,:))**2))*enorm
       else
          r%q(s)=absorptive(e,base(s,:,:))
          call heating_channels(e,low(s,:,:,:),r%raw_channels(:,s))
       endif
    enddo
    if (.not.all(resolved)) then
       do level=1,derivative_levels
          relative=.064/2.**(level-1);h=wr*relative
          call response_snapshot(spec,cmplx(wr+h,0.),kpar,kperp,a,plus,work_low,ierr,bessel,zet_in, &
               warm_only=.true.,auxiliary_beta=auxiliary_beta)
          if (ierr/=0) then
             r%status=heat_response;r%response_status=ierr;goto 900
          endif
          call response_snapshot(spec,cmplx(wr-h,0.),kpar,kperp,a,minus,work_low,ierr,bessel,zet_in, &
               warm_only=.true.,auxiliary_beta=auxiliary_beta)
          if (ierr/=0) then
             r%status=heat_response;r%response_status=ierr;goto 900
          endif
          do s=1,n
             if (spec(s)%tau_s==0. .or. spec(s)%D_s==0.) cycle
             center=hermitian(((wr+h)*plus(s,:,:)-(wr-h)*minus(s,:,:))/(2.*h))
             if (level>1) then
                richardson=(4.*center-prev_center(s,:,:))/3.
                if (level>2) then
                   difference=richardson-prev_r(s,:,:)
                   merr=sqrt(sum(abs(difference)**2))/max(1.,sqrt(sum(abs(prev_r(s,:,:))**2)))
                   perr=abs(quadratic(e,difference))/max(1.,abs(quadratic(e,prev_r(s,:,:))))
                   score=max(merr,perr)/derivative_tolerance
                   pair_score=max(score,previous_score(s))
                   if (level>3 .and. pair_score<=1. .and. pair_score<best_score(s)) then
                      best_score(s)=pair_score;best(s,:,:)=prev_r(s,:,:);resolved(s)=.true.
                      r%species(s)=quadratic(e,best(s,:,:));r%step(s)=2.*relative
                      r%matrix_error(s)=max(merr,previous_score(s)*derivative_tolerance)
                      r%projection_error(s)=max(perr,previous_score(s)*derivative_tolerance)
                      absolute_error(s)=max(abs(quadratic(e,difference)), &
                           abs(quadratic(e,prev_r(s,:,:)-older_r(s,:,:))))
                   endif
                   previous_score(s)=score
                endif
                if (level>2) older_r(s,:,:)=prev_r(s,:,:)
                prev_r(s,:,:)=richardson
             endif
             prev_center(s,:,:)=center
          enddo
       enddo
    endif
    r%status=heat_derivative
    do s=1,n
       if (.not.resolved(s)) r%species_status(s)=species_heat_derivative
    enddo
    if (.not.all(resolved)) goto 900
    roundoff=128.*epsilon(1.)*(abs(r%electric)+abs(r%magnetic)+sum(abs(r%species)))
    r%uncertainty=sum(absolute_error)+roundoff
    call finish_energy(r)
    if (r%status/=heat_ok) goto 900
    if (r%ratio>.05) r%applicability=ior(r%applicability,1)
    value=sum(pack(r%q,spec%tau_s>0. .and. spec%D_s>0.))
    r%gamma_est=-wr*value/r%total
    ! A consistency warning, never a fit: no observed gamma enters H or Q.
    roundoff=64.*epsilon(1.)*(1.+abs(w))+ &
         abs(wr)*64.*epsilon(1.)*sum(abs(base))*enorm/abs(r%total)+ &
         abs(r%gamma_est)*r%uncertainty/abs(r%total)
    if (abs(r%gamma_est-aimag(w))>max(roundoff,.1*(abs(r%gamma_est)+abs(aimag(w))))) &
         r%applicability=ior(r%applicability,2)
    do s=1,n
       if (spec(s)%tau_s==0. .or. spec(s)%D_s==0.) cycle
       if (r%applicability/=0) then
          r%species_status(s)=species_heat_applicability
       else
          r%power(s)=r%q(s)/r%total;r%channels(:,s)=r%raw_channels(:,s)/r%total
          r%species_status(s)=species_heat_ok
       endif
    enddo
    if (any(spec%tau_s>0. .and. spec%D_s>0.) .and. r%applicability==0) &
         r%total_power=sum(pack(r%power,spec%tau_s>0. .and. spec%D_s>0.))
900 arithmetic_status=0
    call finish_guard(guard,arithmetic_status)
    if (arithmetic_status/=0) r%status=heat_arithmetic
    where (.not.ieee_is_finite(r%species)) r%species=nan
    where (.not.ieee_is_finite(r%q)) r%q=nan
    where (.not.ieee_is_finite(r%raw_channels)) r%raw_channels=nan
    if (.not.ieee_is_finite(r%total)) r%total=nan
    if (.not.ieee_is_finite(r%electric)) r%electric=nan
    if (.not.ieee_is_finite(r%magnetic)) r%magnetic=nan
    if (.not.ieee_is_finite(r%species_total)) r%species_total=nan
    if (r%status/=heat_ok) then
       r%power=nan;r%channels=nan;r%total_power=nan
       do s=1,n
          if (spec(s)%tau_s>0. .and. spec(s)%D_s>0. .and. r%species_status(s)/=species_heat_derivative) &
               r%species_status(s)=species_heat_energy
       enddo
    endif
  end subroutine

  subroutine finish_energy(r)
    !! Explicit component/fraction API, also useful for cancellation tests.
    type(heating_result),intent(inout) :: r
    real :: scale,nan
    integer :: ierr
    type(diagnostic_guard) :: guard
    nan=real(diagnostic_nan());r%status=heat_unresolved
    r%f_e=nan;r%f_b=nan;r%f_species=nan;r%fractions=nan;r%f_e_abs=nan
    if (.not.all(ieee_is_finite([r%electric,r%magnetic,r%species,r%uncertainty]))) return
    if (r%uncertainty<0.) return
    call start_guard(guard)
    r%species_total=sum(r%species);r%total=r%electric+r%magnetic+r%species_total
    scale=abs(r%electric)+abs(r%magnetic)+sum(abs(r%species))
    if (.not.ieee_is_finite(scale) .or. .not.ieee_is_finite(r%total)) goto 800
    r%uncertainty=max(r%uncertainty,128.*epsilon(1.)*scale)
    if (scale>0.) r%f_e_abs=r%electric/scale
    if (abs(r%total)<=10.*r%uncertainty) goto 800
    r%f_e=r%electric/r%total;r%f_b=r%magnetic/r%total
    r%f_species=r%species_total/r%total;r%fractions=r%species/r%total
    r%status=heat_ok
800 ierr=0
    call finish_guard(guard,ierr)
    if (ierr/=0) then
       r%status=heat_arithmetic;r%f_e=nan;r%f_b=nan;r%f_species=nan;r%fractions=nan;r%f_e_abs=nan
    endif
  end subroutine

  subroutine heating_channels(e,low,channels)
    complex,intent(in) :: e(3),low(3,3,-1:1)
    real,intent(out) :: channels(7)
    complex :: transverse(3)
    transverse=[e(1),e(2),cmplx(0.,0.)]
    channels(1)=abs(e(2))**2*aimag(low(2,2,0))
    channels(2)=aimag(conjg(e(2))*e(3)*low(2,3,0))
    channels(3)=aimag(conjg(e(3))*e(2)*low(3,2,0))
    channels(4)=abs(e(3))**2*aimag(low(3,3,0))
    channels(5)=absorptive(e,low(:,:,0))
    channels(6)=absorptive(transverse,low(:,:,1));channels(7)=absorptive(transverse,low(:,:,-1))
  end subroutine
end module heating_energy
