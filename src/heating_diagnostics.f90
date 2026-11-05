module heating_diagnostics
  !! Production corrected workflow and versioned companion IO, kept separate
  !! from the historical calc_eigen/FPC and the root/dispersion implementation.
  use vars, only: spec,susc,kpar,kperp,vtp,betap,normalization,nspec,collision_type, &
       heating_convention,dataName,outputName,option,z_backend,bessel_mode
  use electromagnetic_fields, only: calc_em_eigen
  use species_moments, only: calc_species_moments
  use heating_response
  use heating_energy
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: calc_corrected_eigen,native_heating_state,reset_energy_output,write_energy_record,failed_heating_result
  integer,save :: record_count=0
  logical,save :: output_started=.false.
  abstract interface
     complex function dispersion_callback(w)
       complex :: w
     end function
  end interface
contains
  subroutine native_heating_state(input,mode,beta,velocity,kz,kp,physical,a,kpar_d,kperp_d,ierr)
    use vars, only: specie
    type(specie),intent(in) :: input(:)
    character(*),intent(in) :: mode
    real,intent(in) :: beta,velocity,kz,kp
    type(specie),intent(out) :: physical(:)
    real,intent(out) :: a,kpar_d,kperp_d
    integer,intent(out) :: ierr
    real :: length
    integer :: j
    type(diagnostic_guard) :: guard
    ierr=1;a=real(diagnostic_nan());kpar_d=a;kperp_d=a
    if (size(input)/=size(physical) .or. size(input)==0) return
    if (.not.all(ieee_is_finite([velocity,kz,kp]))) return
    if (velocity<=0.) return
    physical=input
    if (mode=='d_ref') then
       a=velocity;kpar_d=kz;kperp_d=kp;ierr=0;return
    endif
    if (mode/='rho_ref') return
    if (.not.all(ieee_is_finite([beta,input(1)%alph_s,input%tau_s]))) return
    if (beta<=0. .or. input(1)%alph_s<=0. .or. any(input%tau_s<=0.)) return
    call start_guard(guard)
    length=sqrt(beta*input(1)%alph_s)
    a=velocity/sqrt(beta);kpar_d=kz/length;kperp_d=kp/length
    do j=1,size(input)
       physical(j)%tau_s=beta/input(j)%tau_s
    enddo
    ierr=0
    if (.not.all(ieee_is_finite([a,kpar_d,kperp_d,physical%tau_s]))) ierr=1
    call finish_guard(guard,ierr)
  end subroutine

  subroutine calc_corrected_eigen(w,e,b,velocity,density,power,split,r,evaluate,bessel,zet_in)
    use vars, only: specie
    complex,intent(in) :: w
    complex,intent(out) :: e(3),b(3),velocity(:,:),density(:)
    real,intent(out) :: power(:),split(:,:)
    type(heating_result),intent(out) :: r
    procedure(dispersion_callback) :: evaluate
    procedure(bessel_callback) :: bessel
    procedure(z_callback) :: zet_in
    type(specie) :: physical(size(spec))
    real :: a,kz,kp
    integer :: status,moment_status(size(spec)),j
    e=diagnostic_nan();b=e;velocity=diagnostic_nan();density=diagnostic_nan()
    power=real(diagnostic_nan());split=real(diagnostic_nan())
    call failed_heating_result(r,heat_invalid)
    call native_heating_state(spec,normalization,betap,vtp,kpar,kperp,physical,a,kz,kp,status)
    if (status/=0) return
    call calc_em_eigen(w,e,b,status,evaluate)
    if (status/=0) then
       r%status=heat_no_field;r%field_status=status;return
    endif
    call calc_species_moments(physical,w,kz,kp,a,e,b,susc,velocity,density,moment_status)
    do j=1,size(spec)
       if (moment_status(j)/=0) write(*,'(a,i0,a,i0)') &
            'WARNING: species moment unavailable; species=',j,'; status=',moment_status(j)
    enddo
    if (collision_type/=0) then
       r%status=heat_collision;r%field_status=0;return
    endif
    call calculate_heating(physical,w,kz,kp,a,e,r,bessel,zet_in)
    r%field_status=0;power=r%power;split=r%channels
  end subroutine

  subroutine failed_heating_result(r,status)
    type(heating_result),intent(out) :: r
    integer,intent(in) :: status
    integer :: j
    call reset_heating(r,nspec,status)
    do j=1,nspec
       if (spec(j)%tau_s==0.) r%species_status(j)=species_heat_cold
       if (spec(j)%D_s==0.) r%species_status(j)=species_heat_absent
    enddo
  end subroutine

  subroutine reset_energy_output()
    !! Prevent a stale companion from a previous run being mistaken for new data.
    integer :: u
    logical :: exists
    character(300) :: path
    record_count=0;output_started=.false.
    if (option/=-1 .and. option/=1 .and. option/=2) return
    path=energy_path()
    inquire(file=trim(path),exist=exists)
    if (exists) then
       open(newunit=u,file=trim(path),status='old');close(u,status='delete')
    endif
  end subroutine

  function energy_path() result(path)
    character(300) :: path
    path='data/'//trim(dataName)//'/'//trim(outputName)//'.energy'
  end function

  subroutine write_energy_record(r,w,workflow,scan_index,outer,point,mode,root_valid,source)
    type(heating_result),intent(in) :: r
    complex,intent(in) :: w
    integer,intent(in) :: workflow,scan_index,outer,point,mode
    logical,intent(in) :: root_valid
    character(*),intent(in) :: source
    character(300) :: path,mode_file
    character(30) :: index_text
    integer :: u,j,c
    if (heating_convention/='corrected') return
    path=energy_path()
    if (.not.output_started) then
       open(newunit=u,file=trim(path),status='replace',action='write')
       write(u,'(a)') '# schema=plume-energy-v1'
       write(u,'(a)') '# heating_convention=corrected; normalization='//trim(normalization)
       write(u,'(a)') '# z_backend='//trim(z_backend)//'; bessel_mode='//trim(bessel_mode)
       write(u,'(a,i0)') '# nspec=',nspec
       write(u,'(a)') '# H = 16*pi*W/abs(deltaEx)^2; formal real-frequency weak-mode expression.'
       write(u,'(a)') '# Signed fractions are not clamped; f_E_abs is a noncancelling contribution measure, not energy fraction.'
       write(u,'(a)') '# energy_status: 0=resolved,1=input,2=fields,3=wr_zero,4=response,5=derivative,'// &
            '6=unresolved_H,7=collisions,8=root,9=arithmetic'
       write(u,'(a)') '# applicability bits: 1=abs(gamma/wr)>0.05,2=weak-balance/near-reactive warning.'
       write(u,'(a)') '# species_status: 0=kinetic,1=cold_NA,2=absent,3=input/response,4=derivative,'// &
            '5=energy,6=applicability'
       write(u,'(a)') '# response_code uses existing kernel statuses (including exact Doppler/cyclotron poles).'
       write(u,'(a)',advance='no') '# columns: record workflow scan outer point mode root_valid field_status '// &
            'energy_status applicability response_code kperp kpar omega_r gamma H_total H_E H_B H_species_total '// &
            'f_E f_B f_species f_E_abs uncertainty gamma_ratio gamma_est total_Ps '
       do j=1,nspec
          write(index_text,'(i0)') j
          write(u,'(a)',advance='no') 'H_s'//trim(index_text)//' f_s'//trim(index_text)//' Q_s'//trim(index_text)// &
               ' Ps_s'//trim(index_text)//' status_s'//trim(index_text)//' dmatrix_s'//trim(index_text)// &
               ' dproject_s'//trim(index_text)//' dstep_s'//trim(index_text)//' '
          do c=1,7
             write(u,'(a,i0,a)',advance='no') 'C',c,'_s'//trim(index_text)//' '
          enddo
       enddo
       write(u,'(a)') 'source_file'
       do j=1,nspec
          write(index_text,'(i0)') j
          if (spec(j)%tau_s==0.) then
             write(u,'(a)') '# species '//trim(index_text)//' model=cold-fluid (kinetic heating unavailable)'
          else
             write(u,'(a)') '# species '//trim(index_text)//' model=kinetic'
          endif
       enddo
       output_started=.true.
    else
       open(newunit=u,file=trim(path),status='old',position='append',action='write')
    endif
    record_count=record_count+1
    write(u,'(11(i0,1x),*(es25.16e3,1x))',advance='no') record_count,workflow,scan_index,outer,point,mode, &
         merge(1,0,root_valid),r%field_status,r%status,r%applicability,r%response_status, &
         kperp,kpar,real(w),aimag(w),r%total,r%electric,r%magnetic,r%species_total,r%f_e,r%f_b, &
         r%f_species,r%f_e_abs,r%uncertainty,r%ratio,r%gamma_est,r%total_power, &
         (r%species(j),r%fractions(j),r%q(j),r%power(j),real(r%species_status(j)),r%matrix_error(j), &
         r%projection_error(j),r%step(j),r%channels(:,j),j=1,nspec)
    mode_file=source
    if (workflow/=-1) then
       write(index_text,'(i0)') mode
       mode_file=trim(source)//'.mode'//trim(index_text)
    endif
    write(u,'(a)') '"'//trim(mode_file)//'"'
    close(u)
    if (r%status/=0 .or. r%applicability/=0) write(*,'(a,4i5)') &
         'WARNING: corrected heating/energy status, applicability, point, mode:',r%status,r%applicability,point,mode
  end subroutine
end module heating_diagnostics
