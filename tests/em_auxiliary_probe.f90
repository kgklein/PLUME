program em_auxiliary_probe
  !! Same physical species and frequency; explicit auxiliary scales are test-only.
  use vars
  use functions, only: read_in_params
  use disprels, only: bessel,zet_in
  use auxiliary_normalization
  use kinetic_core
  use cold_susceptibility, only: calc_chi_cold
  use electromagnetic_fields
  use species_moments
  implicit none
  type(auxiliary_scale) :: scale
  type(kinetic_normalization) :: state
  type(specie) :: warm
  complex :: w,ef(3),bf(3),nx,nz
  real :: wr,wi,bstar
  integer :: ios,i,j,status
  complex,allocatable :: velocity(:,:),density(:)
  integer,allocatable :: moment_status(:)
  call read_in_params
  allocate(velocity(3,nspec),density(nspec),moment_status(nspec))
  do
     read(*,*,iostat=ios) wr,wi,bstar
     if (ios/=0) exit
     w=cmplx(wr,wi)
     call make_auxiliary_scale(kpar,kperp,vtp,bstar,scale,status)
     if (status/=0) error stop 'Bad auxiliary scale'
     call scale%export(state%beta,state%velocity,state%kpar,state%kperp,state%alpha)
     do i=1,nspec
        if (spec(i)%tau_s==0.) then
           call calc_chi_cold(spec(i),w,kpar,kperp,vtp,susc(i,:,:),status)
        else
           call auxiliary_species(spec(i),scale,warm,status)
           call calc_chi_kinetic(warm,w,state,susc(i,:,:),status,bessel,zet_in)
        endif
        if (status/=0) error stop 'Bad species response'
     enddo
     lam=sum(susc,dim=1)
     do j=1,3
        lam(j,j)=lam(j,j)+1.
     enddo
     nx=kperp/(vtp*w);nz=kpar/(vtp*w)
     lam(1,1)=lam(1,1)-nz*nz;lam(2,2)=lam(2,2)-nx*nx-nz*nz;lam(3,3)=lam(3,3)-nx*nx
     lam(1,3)=lam(1,3)+nx*nz;lam(3,1)=lam(3,1)+nx*nz
     call calc_em_fields(lam,w,kpar,kperp,vtp,'d_ref',spec(1)%alph_s,ef,bf,status)
     write(*,'(a,i3,es26.17e3,12es26.17e3)') 'AUX ',status,bstar,ef,bf
     call calc_species_moments(spec,w,kpar,kperp,vtp,ef,bf,susc,velocity,density,moment_status)
     do i=1,nspec
        write(*,'(a,es26.17e3,2i4,8es26.17e3)') 'AUX_MOM ',bstar,i,moment_status(i),velocity(:,i),density(i)
     enddo
  enddo
end program
