program selected_z_probe
  use vars,only: z_backend
  use disprels,only: zet_in
  implicit none
  real :: k,x,y
  complex :: value
  integer :: ios
  do
     read(*,*,iostat=ios) z_backend,k,x,y
     if (ios/=0) exit
     value=zet_in(k,cmplx(x,y))
     write(*,'(2es27.17e3)') real(value),aimag(value)
  enddo
end program
