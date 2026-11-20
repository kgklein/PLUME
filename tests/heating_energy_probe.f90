program heating_energy_probe
  !! Audit only: linked production calls, snapshots at explicitly requested w.
  use vars
  use functions, only: read_in_params
  use disprels, only: disp,rtsec,calc_eigen,bessel,zet_in
  use electromagnetic_fields, only: calc_em_eigen
  use auxiliary_normalization
  use kinetic_core, only: kinetic_normalization
  use adaptive_kinetic, only: calc_chi_adaptive,physical_bessel_argument
  use, intrinsic :: ieee_arithmetic
  implicit none
  complex :: w,d,e(3),b(3),old_e(3),old_b(3),nan
  complex,allocatable :: v(:,:),density(:),saved(:,:,:),saved_low(:,:,:,:),low(:,:,:,:),response(:,:,:)
  real,allocatable :: power(:),split(:,:)
  type(auxiliary_scale) :: auxiliary
  type(kinetic_normalization) :: state
  type(specie) :: kinetic
  real :: wr,wi,argument,bstar
  integer :: mode,ios,flag,status,s,i,j,n
  logical :: finite_low
  call read_in_params
  allocate(v(3,nspec),density(nspec),saved(nspec,3,3),saved_low(nspec,3,3,-1:1), &
       low(nspec,3,3,-1:1),response(nspec,3,3),power(nspec),split(7,nspec))
  nan=cmplx(ieee_value(0.,ieee_quiet_nan),ieee_value(0.,ieee_quiet_nan))
  do
     read(*,*,iostat=ios) wr,wi,mode
     if (ios/=0) exit
     w=cmplx(wr,wi)
     ! Deliberately control scratch initialization in an audit, not production.
     if (low_n) susc_low=cmplx(0.,0.)
     if (mode==2 .and. low_n) susc_low=cmplx(0.,2.)
     if (mode==1) then
        w=rtsec(disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,flag)
        if (flag<0) error stop 'Audit root failed'
     endif
     d=disp(w);saved=susc
     if (low_n) saved_low=susc_low
     if (mode==2 .or. mode==4) then
        status=0;e=0.;b=0.;power=0.;split=0.
     else if (has_cold_species()) then
        call calc_em_eigen(w,e,b,status,disp)
        power=real(nan);split=real(nan)
     else
        call calc_eigen(w,e,b,v,density,power,split,.false.,.true.)
        status=0
     endif
     write(*,'(a,i3,2es26.17e3)') 'POINT ',status,w
     do i=1,3
        write(*,'(a,i3,4es26.17e3)') 'FIELD ',i,e(i),b(i)
     enddo
     do s=1,nspec
        write(*,'(a,i3,8es26.17e3)') 'HEAT ',s,power(s),split(:,s)
        do i=1,3
           do j=1,3
              write(*,'(a,3i3,2es26.17e3)') 'CHI ',s,i,j,saved(s,i,j)
              write(*,'(a,3i3,2es26.17e3)') 'AFTER ',s,i,j,susc(s,i,j)
              if (low_n) then
                 do n=-1,1
                    write(*,'(a,4i3,4es26.17e3)') 'LOW ',s,n,i,j,saved_low(s,i,j,n),susc_low(s,i,j,n)
                 enddo
              endif
           enddo
        enddo
     enddo
     if (mode/=3) cycle
     ! Isolated recoverability check only: existing shared kinetic API with
     ! compute_low=true, never exported to production dispersion/heating state.
     bstar=0.
     do s=1,nspec
        if (spec(s)%tau_s>0.) then
           bstar=spec(s)%tau_s;exit
        endif
     enddo
     if (bstar<=0.) cycle
     call make_auxiliary_scale(kpar,kperp,vtp,bstar,auxiliary,status)
     if (status/=0) error stop 'Audit auxiliary scale failed'
     call auxiliary%export(state%beta,state%velocity,state%kpar,state%kperp,state%alpha)
     response=nan;low=nan
     do s=1,nspec
        if (spec(s)%tau_s==0.) cycle
        call auxiliary_species(spec(s),auxiliary,kinetic,status)
        if (status/=0) error stop 'Audit auxiliary species failed'
        argument=physical_bessel_argument(spec(s),kperp,1.,.true.)
        call calc_chi_adaptive(kinetic,w,state,argument,s,nspec,0,.true.,response(s,:,:),low(s,:,:,:), &
             status,bessel,zet_in)
        if (status/=0) error stop 'Audit isolated harmonics failed'
        do i=1,3
           do j=1,3
              do n=-1,1
                 write(*,'(a,4i3,2es26.17e3)') 'RECOVER ',s,n,i,j,low(s,i,j,n)
              enddo
           enddo
        enddo
     enddo
  enddo
end program
