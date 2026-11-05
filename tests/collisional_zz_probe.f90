program collisional_zz_probe
  !! Production dispatch, with test-only per-frequency collision-rate injection.
  use vars, only: nspec,spec,susc,susc_low,lam,dispersion_status,solver_spec, &
       solver_kperp,bessel_mode,bessel_nfixed,low_n
  use functions, only: read_in_params
  use disprels, only: disp,bessel
  use bessel_cache, only: cache
  use, intrinsic :: ieee_arithmetic
  implicit none
  real :: wr,wi,nu,argument,nan
  complex :: value
  integer :: ios,s,i,j,h,n,cutoff
  call read_in_params
  if (.not.low_n) error stop 'Probe requires low_n'
  nan=ieee_value(0.,ieee_quiet_nan)
  do
     read(*,*,iostat=ios) wr,wi,nu
     if (ios/=0) exit
     spec(:)%nu_ns=nu
     susc_low=cmplx(nan,nan)
     value=disp(cmplx(wr,wi))
     write(*,'(a,i4,2es26.17e3)') 'DET ',dispersion_status,real(value),aimag(value)
     do s=1,nspec
        do i=1,3
           do j=1,3
              write(*,'(a,3i4,2es26.17e3)') 'CHI ',s,i,j,real(susc(s,i,j)),aimag(susc(s,i,j))
              do h=-1,1
                 write(*,'(a,4i4,2es26.17e3)') 'LOW ',s,i,j,h, &
                      real(susc_low(s,i,j,h)),aimag(susc_low(s,i,j,h))
              enddo
           enddo
        enddo
        if (bessel_mode=='fixed' .and. bessel_nfixed==128) then
           cutoff=128
           argument=(solver_kperp**2/2.)*solver_spec(s)%Q_s**2*solver_spec(s)%alph_s/ &
                (solver_spec(s)%mu_s*solver_spec(s)%tau_s*solver_spec(1)%alph_s)
        else
           cutoff=cache(s)%active
        endif
        write(*,'(a,2i6)') 'HARMONICS ',s,cutoff
        do n=0,cutoff+1
           if (bessel_mode=='fixed' .and. bessel_nfixed==128) then
              write(*,'(a,2i6,es26.17e3)') 'GAMMA ',s,n,bessel(n,argument)
           else
              write(*,'(a,2i6,es26.17e3)') 'GAMMA ',s,n,cache(s)%jn(n,1)
           endif
        enddo
     enddo
  enddo
end program
