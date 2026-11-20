program faddeeva_production_probe
  use, intrinsic :: iso_c_binding, only: c_int
  use vars, only: nspec, susc, lam, dispersion_status,spec,kperp,normalization
  use adaptive_kinetic, only: physical_bessel_argument
  use bessel_cache, only: cache
  use functions, only: read_in_params
  use disprels, only: disp, rtsec,bessel
  implicit none
  interface
     function count_z_calls() bind(C,name='plume_test_z_calls') result(count)
       import c_int
       integer(c_int) :: count
     end function
  end interface
  complex :: w, determinant, root
  real :: argument,coefficient
  integer :: harmonic,cutoff
  real :: wr, wi
  real :: started,finished
  integer :: repeats,iteration
  integer :: is, row, col, ios, root_request, iflag
  call read_in_params
  do
     read(*,*,iostat=ios) wr,wi,root_request,repeats
     if (ios < 0) exit
     if (ios /= 0) stop 1
     call cpu_time(started)
     do iteration=1,repeats
     w=cmplx(wr,wi)
     if (root_request == 1) then
        root=rtsec(disp,w*(1.-1.e-7),w*(1.+1.e-7),1.e-13,iflag)
        write(*,'(a,i4,2es26.17e3)') 'ROOT ',iflag,real(root),aimag(root)
        if (iflag < 0) cycle
        w=root
     endif
     determinant=disp(w)
     enddo
     write(*,'(a,i0)') 'Z_CALLS ',count_z_calls()
     call cpu_time(finished)
     write(*,'(a,es26.17e3)') 'CPU_SECONDS ',finished-started
     write(*,'(a,i4,2es26.17e3)') 'DET ',dispersion_status,real(determinant),aimag(determinant)
     do row=1,3
        do col=1,3
           write(*,'(a,2i4,2es26.17e3)') 'LAM ',row,col,real(lam(row,col)),aimag(lam(row,col))
           do is=1,nspec
              write(*,'(a,3i4,2es26.17e3)') 'CHI ',is,row,col,real(susc(is,row,col)),aimag(susc(is,row,col))
           enddo
        enddo
     enddo
     if (dispersion_status==0) then
        do is=1,nspec
           if (spec(is)%tau_s==0.) cycle
           argument=physical_bessel_argument(spec(is),kperp,spec(1)%alph_s,normalization=='d_ref')
           cutoff=128
           if (allocated(cache)) then
              if (cache(is)%valid) cutoff=cache(is)%active
           endif
           do harmonic=0,cutoff+1
              if (allocated(cache)) then
                 if (cache(is)%valid) then
                    coefficient=cache(is)%jn(harmonic,1)
                 else
                    coefficient=bessel(harmonic,argument)
                 endif
              else
                 coefficient=bessel(harmonic,argument)
              endif
              write(*,'(a,2i6,es26.17e3)') 'GAMMA ',is,harmonic,coefficient
           enddo
        enddo
     endif
  enddo
end program faddeeva_production_probe
