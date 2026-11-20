module auxiliary_normalization
  !! Auxiliary thermal adapter, validated in Stage 3B.1 and integrated in 3B.2.
  !! Native d_ref inputs retain the physical reference mass/charge/density.
  !! No plasma globals, reference temperature, or reference anisotropy are read.
  use nrtype, only: SP
  use vars, only: specie
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: auxiliary_scale, make_auxiliary_scale, auxiliary_species
  integer, parameter, public :: auxiliary_ok=0, auxiliary_invalid=1
  integer, parameter, public :: auxiliary_not_warm=2, auxiliary_range=3

  type :: auxiliary_scale
     private
     real(SP) :: beta, velocity, parallel, perpendicular
   contains
     procedure :: export => export_scale
  end type auxiliary_scale

contains

  subroutine make_auxiliary_scale(kpar_d, kperp_d, va_ref_c, b_star, scale, ierr)
    !! The isotropic auxiliary scale has alpha_*=1. The output is unusable on
    !! failure. This conversion does not assess the range of all subsequent
    !! legacy kinetic intermediate products; see the conditioning report.
    real(SP), intent(in) :: kpar_d, kperp_d, va_ref_c, b_star
    type(auxiliary_scale), intent(out) :: scale
    integer, intent(out) :: ierr
    real(SP) :: root_b, nan
    nan=ieee_value(0._SP,ieee_quiet_nan)
    scale=auxiliary_scale(nan,nan,nan,nan)
    ierr=auxiliary_invalid
    if (.not. all(ieee_is_finite([kpar_d,kperp_d,va_ref_c,b_star]))) return
    if (b_star<=0._SP .or. va_ref_c<=0._SP) return
    root_b=sqrt(b_star)
    ierr=auxiliary_range
    ! Prevent overflow before multiplying, including in a checked build.
    if (root_b>1._SP) then
       if (max(abs(kpar_d),abs(kperp_d),va_ref_c)>huge(root_b)/root_b) return
    endif
    scale%beta=b_star
    scale%velocity=va_ref_c*root_b
    scale%parallel=kpar_d*root_b
    scale%perpendicular=kperp_d*root_b
    if (scale%velocity==0._SP) return
    if (kpar_d/=0._SP .and. scale%parallel==0._SP) return
    if (kperp_d/=0._SP .and. scale%perpendicular==0._SP) return
    ierr=auxiliary_ok
  end subroutine make_auxiliary_scale

  subroutine export_scale(scale, beta, velocity, kpar, kperp, normalization_alpha)
    class(auxiliary_scale), intent(in) :: scale
    real(SP), intent(out) :: beta, velocity, kpar, kperp, normalization_alpha
    beta=scale%beta
    velocity=scale%velocity
    kpar=scale%parallel
    kperp=scale%perpendicular
    normalization_alpha=1._SP
  end subroutine export_scale

  subroutine auxiliary_species(physical, scale, kinetic, ierr)
    !! Convert one warm species, leaving its physical anisotropy and all other
    !! species fields intact. A cold species is rejected BEFORE division.
    type(specie), intent(in) :: physical
    type(auxiliary_scale), intent(in) :: scale
    type(specie), intent(out) :: kinetic
    integer, intent(out) :: ierr
    real(SP) :: nan
    nan=ieee_value(0._SP,ieee_quiet_nan)
    kinetic=specie(nan,nan,nan,nan,nan,nan,nan)
    ierr=auxiliary_invalid
    if (.not. ieee_is_finite(physical%tau_s)) return
    if (physical%tau_s<0._SP) return
    if (physical%tau_s==0._SP) then
       ierr=auxiliary_not_warm
       return
    endif
    if (.not. all(ieee_is_finite([scale%beta,scale%velocity,scale%parallel,scale%perpendicular, &
         physical%mu_s,physical%alph_s,physical%Q_s,physical%D_s,physical%vv_s]))) return
    if (scale%beta<=0._SP .or. scale%velocity<=0._SP) return
    if (physical%mu_s<=0._SP .or. physical%alph_s<=0._SP .or. &
         physical%Q_s==0._SP .or. physical%D_s<0._SP) return
    ierr=auxiliary_range
    if (physical%tau_s<1._SP) then
       if (scale%beta>huge(nan)*physical%tau_s) return
    endif
    kinetic=physical
    kinetic%tau_s=scale%beta/physical%tau_s
    if (kinetic%tau_s==0._SP) return
    ierr=auxiliary_ok
  end subroutine auxiliary_species
end module auxiliary_normalization
