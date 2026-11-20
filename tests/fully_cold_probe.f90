program fully_cold_probe
  use vars, only: nspec, susc, lam, dispersion_status
  use vars, only: solver_spec, solver_vtp, solver_betap, solver_kpar, solver_kperp
  use cold_susceptibility, only: invalid_complex
  use functions, only: read_in_params
  use disprels, only: disp, rtsec
  implicit none
  complex :: w, determinant, root
  real :: wr, wi
  integer :: is, row, col, ios, root_request, iflag
  solver_vtp=real(invalid_complex())
  solver_betap=solver_vtp
  solver_kpar=solver_vtp
  solver_kperp=solver_vtp
  call read_in_params
  if (allocated(solver_spec)) error stop "All-cold input entered thermal adapter"
  do
     read(*,*,iostat=ios) wr,wi,root_request
     if (ios < 0) exit
     if (ios /= 0) stop 1
     w=cmplx(wr,wi)
     if (root_request == 1) then
        root=rtsec(disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,iflag)
        write(*,'(a,i4,2es26.17e3)') 'ROOT ',iflag,real(root),aimag(root)
        if (iflag < 0) cycle
        w=root
     endif
     determinant=disp(w)
     if (allocated(solver_spec)) error stop "All-cold evaluation entered thermal adapter"
     write(*,'(a,i4,2es26.17e3)') 'DET ',dispersion_status,real(determinant),aimag(determinant)
     do row=1,3
        do col=1,3
           write(*,'(a,2i4,2es26.17e3)') 'LAM ',row,col,real(lam(row,col)),aimag(lam(row,col))
           do is=1,nspec
              write(*,'(a,3i4,2es26.17e3)') 'CHI ',is,row,col,real(susc(is,row,col)),aimag(susc(is,row,col))
           enddo
        enddo
     enddo
  enddo
end program fully_cold_probe
