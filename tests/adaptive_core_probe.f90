program adaptive_core_probe
  use vars, only: z_backend
  use vars, only: specie,bessel_mode,bessel_tol
  use kinetic_core, only: kinetic_normalization
  use adaptive_kinetic, only: calc_chi_adaptive,physical_bessel_argument
  use auxiliary_normalization
  use bessel_cache, only: cache,write_bessel_stats
  use disprels, only: bessel,zet_in
  implicit none
  type(specie) :: physical,converted
  type(auxiliary_scale) :: scale
  type(kinetic_normalization) :: state
  real :: kp,kt,a,b,wr,wi,x
  complex :: chi(3,3),low(3,3,-1:1)
  integer :: slot,count,status,ios,j,k
  z_backend='legacy'
  bessel_mode='adaptive'
  do
     read(*,*,iostat=ios) slot,count,kp,kt,a,b,physical%tau_s,physical%mu_s, &
          physical%alph_s,physical%Q_s,physical%D_s,physical%vv_s,wr,wi,bessel_tol
     if (ios/=0) exit
     physical%nu_ns=0.
     call make_auxiliary_scale(kp,kt,a,b,scale,status)
     if (status==0) call auxiliary_species(physical,scale,converted,status)
     if (status==0) then
        call scale%export(state%beta,state%velocity,state%kpar,state%kperp,state%alpha)
        x=physical_bessel_argument(physical,kt,1.,.true.)
        call calc_chi_adaptive(converted,cmplx(wr,wi),state,x,slot,count,0,.false.,chi,low,status,bessel,zet_in)
     endif
     write(*,'(a,i4)') 'STATUS ',status
     if (status/=0) cycle
     do j=1,3
        do k=1,3
           write(*,'(a,2i4,2es26.17e3)') 'CHI ',j,k,real(chi(j,k)),aimag(chi(j,k))
        enddo
     enddo
     write(*,'(a,3i12)') 'COUNTERS ',cache(slot)%hits,cache(slot)%misses,cache(slot)%evaluations
     call write_bessel_stats(6)
  enddo
end program
