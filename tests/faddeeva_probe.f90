module faddeeva_probe
  use, intrinsic :: iso_c_binding
  use plasma_z_faddeeva
  use plasma_dispersion_probe, only: audit_tensor
  implicit none
contains
  subroutine probe_faddeeva(k,x,y,repeats,value,info) bind(C)
    real(c_double),value :: k,x,y
    integer(c_int),value :: repeats
    real(c_double),intent(out) :: value(3)
    integer(c_int),intent(out) :: info(3)
    complex(c_double_complex) :: z,result
    real :: started,finished
    integer :: j,status
    z=cmplx(x,y,c_double_complex)
    call cpu_time(started)
    do j=1,repeats
       call calc_z_faddeeva_signed(k,z,result,status)
    enddo
    call cpu_time(finished)
    value=[real(result,c_double),aimag(result),finished-started]
    info=[status,storage_size(result),storage_size(cmplx(0.,0.))]
  end subroutine

  subroutine reference_native(k,x,y,vr,vi) bind(C)
    real(c_double),value :: k,x,y
    real(c_double),intent(out) :: vr,vi
    complex(c_double_complex) :: result
    integer :: status
    call calc_z_faddeeva_signed(k,cmplx(x,y,c_double_complex),result,status)
    vr=real(result,c_double);vi=aimag(result)
  end subroutine

  subroutine native_tensor(species,parameters,n,collisions,gamma,repeats, &
       output,arguments,values,info,seconds,low) bind(C)
    integer(c_int),value :: n,collisions,repeats
    real(c_double),intent(in) :: species(7),parameters(7),gamma(0:n+1)
    real(c_double),intent(out) :: output(18),arguments(2,2*n+1),values(2,2*n+1),seconds,low(54)
    integer(c_int),intent(out) :: info(6)
    call audit_tensor(species,parameters,n,collisions,2,1.e-7_c_double,gamma, &
         c_funloc(reference_native),repeats,output,arguments,values,info,seconds,low)
  end subroutine
end module
