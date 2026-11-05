program probe
  use bessel_cache
  implicit none
  type(bessel_entry), pointer :: entry
  real :: x,t
  integer :: s,count,cap,fixed,status,ios,n
  character(16) :: mode
  do
     read(*,*,iostat=ios) s,count,x,mode,t,cap,fixed
     if (ios/=0) exit
     call fetch_bessel(s,count,x,mode,t,cap,fixed,entry,status)
     write(*,'(a,i0)') 'STATUS ',status
     if (status/=0) cycle
     call write_bessel_stats(6)
     write(*,'(a,*(es27.17e3,1x))') 'VALUES ',entry%jn(0:entry%active+1,1)
  enddo
end program
