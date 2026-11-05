module cold_reference
  !! Auxiliary warm responses only. Cold preflight, wave assembly and root
  !! handling remain in disprels. No authoritative input or solver_* state changes.
  use vars, only: specie
  use auxiliary_normalization
  use adaptive_kinetic, only: dynamic_harmonics, calc_chi_adaptive, physical_bessel_argument
  use kinetic_core
  implicit none
  private
  public :: calc_warm_auxiliary
contains
  subroutine calc_warm_auxiliary(spec,omega,kpar,kperp,va,chi,ierr,bessel,zet_in)
    type(specie), intent(in) :: spec(:)
    complex, intent(in) :: omega
    real, intent(in) :: kpar,kperp,va
    complex, intent(inout) :: chi(:,:,:)
    integer, intent(out) :: ierr
    interface
       real function bessel(n,x)
         integer :: n
         real :: x
       end function
       complex function zet_in(k,z)
         real :: k
         complex :: z
       end function
    end interface
    type(auxiliary_scale) :: scale
    type(kinetic_normalization) :: state
    type(specie) :: kinetic
    complex :: response(3,3),low(3,3,-1:1)
    real :: argument
    integer :: j,first_warm
    first_warm=0
    do j=1,size(spec)
       if (spec(j)%tau_s>0.) then
          first_warm=j
          exit
       endif
    enddo
    ierr=1
    if (first_warm==0) return
    call make_auxiliary_scale(kpar,kperp,va,spec(first_warm)%tau_s,scale,ierr)
    if (ierr/=0) then
       if (ierr==auxiliary_range) ierr=5
       return
    endif
    call scale%export(state%beta,state%velocity,state%kpar,state%kperp,state%alpha)
    do j=1,size(spec)
       if (spec(j)%tau_s==0.) cycle
       call auxiliary_species(spec(j),scale,kinetic,ierr)
       if (ierr/=0) then
          if (ierr==auxiliary_range) ierr=5
          return
       endif
       if (dynamic_harmonics()) then
          argument=physical_bessel_argument(spec(j),kperp,1.,.true.)
          call calc_chi_adaptive(kinetic,omega,state,argument,j,size(spec),0,.false., &
               response,low,ierr,bessel,zet_in)
       else
          call calc_chi_kinetic(kinetic,omega,state,response,ierr,bessel,zet_in)
       endif
       if (ierr/=0) return
       chi(j,:,:)=response
    enddo
  end subroutine
end module cold_reference
