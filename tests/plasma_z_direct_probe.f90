program direct_z_probe
  !! Links the current production objects without recompiling their arithmetic.
  use disprels,only: zet_in
  use vars,only: z_backend
  implicit none
  real :: k,x,y
  complex :: value
  integer :: ios
  z_backend='legacy'
  do
     read(*,*,iostat=ios) k,x,y
     if (ios/=0) exit
     value=zet_in(k,cmplx(x,y))
     write(*,'(2es27.17e3)') real(value),aimag(value)
  enddo
end program
