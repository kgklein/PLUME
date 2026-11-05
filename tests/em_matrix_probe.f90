program em_matrix_probe
  use electromagnetic_fields
  implicit none
  complex :: matrix(3,3),omega,e(3),b(3)
  real :: kpar,kperp,velocity,alpha
  character(16) :: mode
  integer :: i,j,ios,status
  do
     read(*,*,iostat=ios) omega,kpar,kperp,velocity,alpha,mode
     if (ios/=0) exit
     do i=1,3
        read(*,*) matrix(i,:)
     enddo
     call calc_em_fields(matrix,omega,kpar,kperp,velocity,mode,alpha,e,b,status)
     write(*,'(i3,12es26.17e3)') status,e,b
  enddo
end program
