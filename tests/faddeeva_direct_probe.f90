program faddeeva_direct_probe
  use plasma_z_faddeeva
  use, intrinsic :: iso_c_binding
  implicit none
  real(c_double) :: k,x,y
  complex(c_double_complex) :: value
  integer :: ios,status
  integer(c_int64_t) :: bits(2)
  do
     read(*,*,iostat=ios) k,x,y
     if (ios/=0) exit
     call calc_z_faddeeva_signed(k,cmplx(x,y,c_double_complex),value,status)
     bits=transfer(value,bits)
     write(*,'(i4,2es27.17e3,2(1x,z16.16))') status,real(value,c_double),aimag(value),bits
  enddo
end program
