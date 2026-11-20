program corrected_heating_probe
  use vars
  use functions, only: read_in_params
  use disprels, only: disp,rtsec,bessel,zet_in
  use heating_energy
  use heating_diagnostics
  use electromagnetic_fields, only: calc_em_eigen
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(heating_result) :: r
  type(specie),allocatable :: physical(:)
  complex :: w,value,e(3),b(3),old_lam(3,3)
  complex,allocatable :: velocity(:,:),density(:),saved(:,:,:),saved_low(:,:,:,:)
  real,allocatable :: power(:),split(:,:)
  real :: wr,wi,a,kz,kp,bstar
  integer :: mode,flag,ios,i,s,status,old_status
  call read_in_params
  allocate(physical(nspec),velocity(3,nspec),density(nspec),saved(nspec,3,3),power(nspec),split(7,nspec))
  if (allocated(susc_low)) allocate(saved_low(nspec,3,3,-1:1))
  do
     read(*,*,iostat=ios) wr,wi,mode
     if (ios/=0) exit
     w=cmplx(wr,wi)
     if (mode==1) then
        w=rtsec(disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,flag)
        if (flag<0) error stop 'Corrected probe root failed'
     endif
     value=disp(w);saved=susc;old_lam=lam;old_status=dispersion_status
     if (allocated(saved_low)) saved_low=susc_low
     if (mode==2) then
        read(*,*) e
        call native_heating_state(spec,normalization,betap,vtp,kpar,kperp,physical,a,kz,kp,status)
        if (status/=0) error stop 'Bad native state'
        call calculate_heating(physical,w,kz,kp,a,e,r,bessel,zet_in)
        b=0.;velocity=0.;density=0.
     else
        call calc_corrected_eigen(w,e,b,velocity,density,power,split,r,disp,bessel,zet_in)
        if (mode==3 .or. mode==4) then
           bstar=1.
           if (mode==3) read(*,*) bstar
           call native_heating_state(spec,normalization,betap,vtp,kpar,kperp,physical,a,kz,kp,status)
           if (mode==3) then
              call calculate_heating(physical,w,kz,kp,a,e,r,bessel,zet_in,bstar)
           else
              call calculate_heating(physical,w,kz,kp,a,e,r,bessel,unresolved_z)
           endif
        endif
     endif
     if (any(transfer(saved,[0_int64],18*nspec)/=transfer(susc,[0_int64],18*nspec))) &
          error stop 'Corrected heating changed species tensors'
     if (any(transfer(old_lam,[0_int64],18)/=transfer(lam,[0_int64],18))) &
          error stop 'Corrected heating changed wave matrix'
     if (allocated(saved_low)) then
        if (any(transfer(saved_low,[0_int64],54*nspec)/=transfer(susc_low,[0_int64],54*nspec))) &
             error stop 'Corrected heating changed mode harmonics'
     endif
     if (dispersion_status/=old_status) error stop 'Corrected heating changed solver status'
     write(*,'(a,2es26.17e3)') 'MODE ',w
     write(*,'(a,4i5)') 'STATUS ',r%status,r%applicability,r%field_status,r%response_status
     write(*,'(a,12es26.17e3)') 'ENERGY ',r%total,r%electric,r%magnetic,r%species_total,r%f_e,r%f_b, &
          r%f_species,r%f_e_abs,r%uncertainty,r%total_power,r%gamma_est,r%ratio
     do i=1,3
        write(*,'(a,i4,4es26.17e3)') 'FIELD ',i,e(i),b(i)
     enddo
     do s=1,nspec
        write(*,'(a,2i4,7es26.17e3)') 'SPEC ',s,r%species_status(s),r%species(s),r%fractions(s),r%q(s),r%power(s), &
             r%matrix_error(s),r%projection_error(s),r%step(s)
        write(*,'(a,i4,7es26.17e3)') 'RAW_CHANNEL ',s,r%raw_channels(:,s)
        write(*,'(a,i4,7es26.17e3)') 'CHANNEL ',s,r%channels(:,s)
        write(*,'(a,i4,8es26.17e3)') 'MOM ',s,velocity(:,s),density(s)
     enddo
  enddo
contains
  complex function unresolved_z(k,z) result(value)
    ! Deterministic test-only frequency noise: exercise derivative failure,
    ! never a physical reference or production numerical correction.
    real :: k
    complex :: z
    value=zet_in(k,z)+cmplx(.01*sin(1.e8*real(z)),0.)
  end function
end program
