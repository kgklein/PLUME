program harmonic_probe
  use vars, only: nspec, susc, lam, dispersion_status
  use functions, only: read_in_params
  use disprels, only: disp, rtsec
  use bessel_cache, only: write_bessel_stats,cache
  use adaptive_kinetic, only: finish_bessel_diagnostics
  implicit none
  complex :: w, determinant, root
  real :: wr, wi,started,finished
  integer :: is, row, col, ios, root_request, iflag,repeats,iteration
  call read_in_params
  do
     read(*,*,iostat=ios) wr,wi,root_request,repeats
     if (ios < 0) exit
     if (ios /= 0) stop 1
     call cpu_time(started)
     do iteration=1,abs(repeats)
     w=cmplx(wr,wi)
     if (root_request == 1) then
        if (repeats<0) then
           root=rtsec(uncached_disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,iflag)
        else
           root=rtsec(disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,iflag)
        endif
        if (iteration==abs(repeats)) write(*,'(a,i4,2es26.17e3)') 'ROOT ',iflag,real(root),aimag(root)
        if (iflag < 0) cycle
        w=root
     endif
     determinant=disp(w)
     enddo
     call cpu_time(finished)
     write(*,'(a,es24.15)') 'CPU_SECONDS ',finished-started
     write(*,'(a,i4,2es26.17e3)') 'DET ',dispersion_status,real(determinant),aimag(determinant)
     do row=1,3
        do col=1,3
           write(*,'(a,2i4,2es26.17e3)') 'LAM ',row,col,real(lam(row,col)),aimag(lam(row,col))
           do is=1,nspec
              write(*,'(a,3i4,2es26.17e3)') 'CHI ',is,row,col,real(susc(is,row,col)),aimag(susc(is,row,col))
           enddo
        enddo
     enddo
     call write_bessel_stats(6)
  enddo
  call finish_bessel_diagnostics
contains
  complex function uncached_disp(w)
    complex :: w
    if (allocated(cache)) cache%valid=.false.
    uncached_disp=disp(w)
  end function
end program harmonic_probe
