program moments_eigen_probe
  use vars
  use functions, only: read_in_params
  use disprels, only: disp,rtsec,calc_eigen
  use electromagnetic_fields
  use species_moments
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  complex :: w,determinant,e(3),b(3),old_e(3),old_b(3)
  complex,allocatable :: velocity(:,:),density(:),old_v(:,:),old_n(:),saved(:,:,:)
  real,allocatable :: heating(:),split(:,:)
  integer,allocatable :: status(:)
  real :: wr,wi,kz,kp,a,length
  integer :: request,flag,ios,i,j,s,em_status
  call read_in_params
  allocate(velocity(3,nspec),density(nspec),old_v(3,nspec),old_n(nspec),saved(nspec,3,3), &
       heating(nspec),split(7,nspec),status(nspec))
  do
     read(*,*,iostat=ios) wr,wi,request
     if (ios/=0) exit
     w=cmplx(wr,wi)
     if (request==1) then
        w=rtsec(disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,flag)
        if (flag<0) error stop 'Probe root failed'
     endif
     determinant=disp(w);saved=susc
     a=vtp;kz=kpar;kp=kperp
     if (normalization=='rho_ref') then
        a=vtp/sqrt(betap);length=sqrt(betap*spec(1)%alph_s)
        kz=kpar/length;kp=kperp/length
     endif
     if (request==3) then
        read(*,*) e,b
        em_status=0
        call calc_species_moments(spec,w,kz,kp,a,e,b,susc,velocity,density,status)
     else if (has_cold_species()) then
        call calc_mixed_eigen(w,e,b,velocity,density,heating,split,.false.,em_status,disp,status)
     else
        call calc_em_eigen(w,e,b,em_status,disp)
        call calc_species_moments(spec,w,kz,kp,a,e,b,susc,velocity,density,status)
     endif
     if (em_status==0) then
        if (any(transfer(saved,[0_int64],size=18*nspec)/=transfer(susc,[0_int64],size=18*nspec))) &
             error stop 'Moment diagnostics changed species susceptibility'
     endif
     write(*,'(a,i0,2es26.17e3)') 'EM ',em_status,w
     do i=1,3
        write(*,'(a,i3,4es26.17e3)') 'FIELD ',i,e(i),b(i)
        do j=1,3
           write(*,'(a,2i3,2es26.17e3)') 'LAM ',i,j,lam(i,j)
           do s=1,nspec
              write(*,'(a,3i3,2es26.17e3)') 'CHI ',s,i,j,susc(s,i,j)
           enddo
        enddo
     enddo
     do s=1,nspec
        write(*,'(a,2i4,8es26.17e3)') 'MOM ',s,status(s),velocity(:,s),density(s)
     enddo
     if (request==2 .and. .not.has_cold_species()) then
        call calc_eigen(w,old_e,old_b,old_v,old_n,heating,split,.true.,.false.)
        do s=1,nspec
           write(*,'(a,i4,8es26.17e3)') 'OLD ',s,old_v(:,s),old_n(s)
        enddo
     endif
     if (has_cold_reference() .and. allocated(solver_spec)) error stop 'Cold reference constructed thermal state'
  enddo
end program
