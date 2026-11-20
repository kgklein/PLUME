program em_fields_probe
  use vars
  use functions, only: read_in_params
  use disprels, only: disp,rtsec,calc_eigen
  use electromagnetic_fields
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  complex :: omega,determinant,ef(3),bf(3),old_e(3),old_b(3)
  complex,allocatable :: saved(:,:,:),density(:),velocity(:,:)
  real,allocatable :: heating(:),split(:,:)
  real :: wr,wi
  integer :: request,ios,flag,i,j,k,status
  call read_in_params
  allocate(saved(nspec,3,3),density(nspec),velocity(3,nspec),heating(nspec),split(7,nspec))
  do
     read(*,*,iostat=ios) wr,wi,request
     if (ios/=0) exit
     omega=cmplx(wr,wi)
     if (request==1) then
        omega=rtsec(disp,omega*(1.-1.e-7),omega*(1.+1.e-7),1.e-13,flag)
        if (flag<0) error stop 'Probe root failed'
     endif
     determinant=disp(omega);saved=susc
     call calc_em_eigen(omega,ef,bf,status,disp)
     if (status==0) then
        if (any(transfer(saved,[0_int64],size=18*nspec)/=transfer(susc,[0_int64],size=18*nspec))) &
             error stop 'EM calculation changed species responses'
     endif
     write(*,'(a,i0,2es26.17e3)') 'EM ',status,omega
     do i=1,3
        write(*,'(a,i3,4es26.17e3)') 'FIELD ',i,ef(i),bf(i)
        do j=1,3
           write(*,'(a,2i3,2es26.17e3)') 'LAM ',i,j,lam(i,j)
           do k=1,nspec
              write(*,'(a,3i3,2es26.17e3)') 'CHI ',k,i,j,susc(k,i,j)
           enddo
        enddo
     enddo
     if (request==2 .and. .not.has_cold_species()) then
        call calc_eigen(omega,old_e,old_b,velocity,density,heating,split,.false.,.false.)
        do i=1,3
           write(*,'(a,i3,4es26.17e3)') 'LEGACY_FIELD ',i,old_e(i),old_b(i)
        enddo
     endif
     if (has_cold_reference() .and. allocated(solver_spec)) error stop 'EM constructed cold thermal scale'
  enddo
end program
