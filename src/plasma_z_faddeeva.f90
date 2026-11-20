module plasma_z_faddeeva
  !! Standalone binary64 Landau Z=i*sqrt(pi)*w(z), with signed k convention.
  !! No plasma globals. The C bridge preserves the caller's floating-point state,
  !! including FTZ/DAZ and traps, and returns assigned quiet NaNs on failure.
  use, intrinsic :: iso_c_binding, only: c_double,c_double_complex,c_int
  implicit none
  private
  public :: calc_z_faddeeva,calc_z_faddeeva_signed,z_faddeeva_value
  public :: z_ok,z_invalid,z_zero_k,z_range,z_numerical
  integer,parameter :: z_ok=0,z_invalid=1,z_zero_k=2,z_range=3,z_numerical=4
  interface
     pure function z_faddeeva_value(kpar,z) bind(C,name='plume_faddeeva_z_value') result(value)
       import c_double,c_double_complex
       real(c_double),value :: kpar
       complex(c_double_complex),value :: z
       complex(c_double_complex) :: value
     end function
     function bridge(kpar,z,status) bind(C,name='plume_faddeeva_z') result(value)
       import c_double,c_double_complex,c_int
       real(c_double),value :: kpar
       complex(c_double_complex),value :: z
       integer(c_int),intent(out) :: status
       complex(c_double_complex) :: value
     end function
  end interface
contains
  subroutine calc_z_faddeeva(z,value,status)
    complex(c_double_complex),intent(in) :: z
    complex(c_double_complex),intent(out) :: value
    integer,intent(out) :: status
    call calc_z_faddeeva_signed(1._c_double,z,value,status)
  end subroutine

  subroutine calc_z_faddeeva_signed(kpar,z,value,status)
    real(c_double),intent(in) :: kpar
    complex(c_double_complex),intent(in) :: z
    complex(c_double_complex),intent(out) :: value
    integer,intent(out) :: status
    integer(c_int) :: c_status
    value=bridge(kpar,z,c_status)
    status=int(c_status)
  end subroutine
end module plasma_z_faddeeva
