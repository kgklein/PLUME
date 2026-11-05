program probe
  use harmonic_selection
  implicit none
  real :: x,t,e0,e2
  integer :: cap,n,status,ios
  logical :: converged
  do
     read(*,*,iostat=ios) x,t,cap
     if (ios/=0) exit
     call select_harmonics(x,t,cap,n,e0,e2,converged,status)
     write(*,'(3i10,2es27.17e3)') status,n,merge(1,0,converged),e0,e2
  enddo
end program
