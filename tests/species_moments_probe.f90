program species_moments_probe
  use vars
  use species_moments
  use cold_susceptibility, only: calc_chi_cold
  implicit none
  type(specie) :: s
  complex :: w,e(3),b(3),v(3),density,chi(3,3)
  real :: kp,kz,a
  integer :: mode,ios,status,chi_status,i,j
  betap=0.;vtp=0.;normalization='d_ref'
  do
     read(*,*,iostat=ios) mode,w,kz,kp,a,s%tau_s,s%mu_s,s%alph_s,s%Q_s,s%D_s,s%vv_s
     if (ios/=0) exit
     read(*,*) e,b
     chi=cmplx(0.,0.)
     if (mode==0) then
        call calc_cold_moments(s,w,kz,kp,a,e,b,v,density,status)
        if (status==0) call calc_chi_cold(s,w,kz,kp,a,chi,chi_status)
     else
        do i=1,3
           read(*,*) chi(i,:)
        enddo
        call calc_kinetic_moments(s,w,kz,kp,a,e,chi,v,density,status)
     endif
     write(*,'(a,i0,8es26.17e3)') 'MOM ',status,v,density
     do i=1,3
        do j=1,3
           write(*,'(a,2i3,2es26.17e3)') 'CHI ',i,j,chi(i,j)
        enddo
     enddo
     if (allocated(spec) .or. allocated(solver_spec)) error stop 'Moment kernel constructed thermal state'
  enddo
end program
