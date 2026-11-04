module cold_susceptibility
  !! Standalone, single-species cold-fluid susceptibility in inertial units.
  !! Independently implements Verscharen & Chandran, ApJ 764, 88 (2013),
  !! equations (5)-(9). See tests/COLD_SUSCEPTIBILITY_DESIGN.md for conventions.
  !! Called before thermal evaluation for cold components in Stage 2.
  use nrtype, only: SP, SPC
  use vars, only: specie
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: calc_chi_cold
  public :: finite_complex, invalid_complex
  integer, parameter, public :: cold_ok = 0, cold_invalid_input = 1
  integer, parameter, public :: cold_zero_frequency = 2, cold_zero_doppler = 3
  integer, parameter, public :: cold_cyclotron_pole = 4

contains

  logical function finite_complex(value)
    !! Keep IEEE checks outside fast-math compilation units.
    complex(SPC), value, intent(in) :: value
    finite_complex = ieee_is_finite(real(value,SP)) .and. ieee_is_finite(aimag(value))
  end function finite_complex

  complex(SPC) function invalid_complex()
    real(SP) :: nan
    nan = ieee_value(0.0_SP, ieee_quiet_nan)
    invalid_complex = cmplx(nan,nan,kind=SPC)
  end function invalid_complex

  subroutine calc_chi_cold(species, omega, kpar_d, kperp_d, va_ref_c, chi, ierr)
    !! Return ordinary chi_s, without identity or wave-matrix terms.
    !! Used inputs must be finite; mu_s and va_ref_c positive, Q_s nonzero,
    !! D_s nonnegative. Signed Q_s determines the signed gyrofrequency.
    !! Zero density returns zero after input validation. For populated species,
    !! exact frequency, Doppler and cyclotron poles return distinct statuses.
    !! Every error returns quiet NaNs; the caller must inspect ierr.
    !! No temperature, anisotropy, collision rate or global state is read.
    type(specie), intent(in) :: species
    complex(SPC), intent(in) :: omega !! omega/Omega_R, exp(-i*omega*t).
    real(SP), intent(in) :: kpar_d, kperp_d !! Signed components of k*d_R.
    real(SP), intent(in) :: va_ref_c !! v_A,R/c, independent of temperature.
    complex(SPC), intent(out) :: chi(3,3)
    integer, intent(out) :: ierr
    real(SP) :: gyro, plasma, perpendicular_drift, nan
    complex(SPC) :: delta, denominator, factor, transverse, hall, coupling, magnetic
    complex(SPC), parameter :: imaginary = (0.0_SP, 1.0_SP)

    nan = ieee_value(0.0_SP, ieee_quiet_nan)
    chi = cmplx(nan, nan, kind=SPC)
    ierr = cold_invalid_input
    if (.not. all(ieee_is_finite([real(omega,SP), aimag(omega), kpar_d, kperp_d, &
         va_ref_c, species%mu_s, species%Q_s, species%D_s, species%vv_s]))) return
    if (species%mu_s <= 0.0_SP .or. va_ref_c <= 0.0_SP .or. &
         species%Q_s == 0.0_SP .or. species%D_s < 0.0_SP) return

    if (species%D_s == 0.0_SP) then
       chi = cmplx(0.0_SP, 0.0_SP, kind=SPC)
       ierr = cold_ok
       return
    endif
    ierr = cold_zero_frequency
    if (omega == cmplx(0.0_SP, 0.0_SP, kind=SPC)) return
    delta = omega - kpar_d*species%vv_s
    ierr = cold_zero_doppler
    if (delta == cmplx(0.0_SP, 0.0_SP, kind=SPC)) return
    gyro = species%mu_s/species%Q_s
    ierr = cold_cyclotron_pole
    if (delta == gyro .or. delta == -gyro) return

    plasma = species%D_s*species%mu_s/(species%Q_s**2*va_ref_c**2)
    perpendicular_drift = kperp_d*species%vv_s
    denominator = (delta-gyro)*(delta+gyro)
    factor = (plasma/omega**2)/denominator
    ! Combined R/L expressions avoid cancellation in S=(R+L)/2, H=(R-L)/2.
    transverse = -factor*delta**2
    hall = factor*delta*gyro
    coupling = -factor*perpendicular_drift*delta
    magnetic = imaginary*factor*perpendicular_drift*gyro

    chi(1,1) = transverse
    chi(1,2) = -imaginary*hall
    chi(1,3) = coupling
    chi(2,1) = -chi(1,2)
    chi(2,2) = transverse
    chi(2,3) = magnetic
    chi(3,1) = coupling
    chi(3,2) = -magnetic
    ! Cancel the redundant omega**2 in the longitudinal electrostatic term.
    chi(3,3) = -plasma/delta**2 - factor*perpendicular_drift**2
    ierr = cold_ok
  end subroutine calc_chi_cold

end module cold_susceptibility
