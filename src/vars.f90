!=============================================================================!
!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kris.klein@gmail.com                                                       !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!!*VARIABLES                                                                *!!
!=============================================================================!
!=============================================================================!
module vars
  implicit none
  private

  !INPUT PARAMETERS
  !      The following parameters will be read in from 'scripts/*.in'
  !Global Input Parameters
  real, target :: betap     !beta_p=8 pi n_ref T_p/B^2
  real, target :: kperp     !k_perp rho_ref
  real, target :: kpar      !k_parallel rho_ref
  real, target :: vtp       !v_t,ref||/c

  !Numer of species
  integer :: nspec=3
  
  !Number of parameter scans
  integer :: nscan = 0

  !Species Input Parameters
  public :: specie
     type :: specie
        real :: tau_s     !T_ref/T_s|_parallel
        real :: mu_s      !m_ref/m_s
        real :: alph_s    !T_perp/T_parallel_s
        real :: Q_s       !q_ref/q_s
        real :: D_s       !n_s/n_ref
        real :: vv_s      !v_drift/c
     end type specie
  type (specie), dimension (:), allocatable, target :: spec 
  type (specie), dimension (:,:), allocatable, target :: rad_spec 

  !Scan Parameters
  public :: scanner
     type :: scanner
        real :: range_i       !initial value
        real :: range_f       !final value
        logical :: log_scan   !T-> log, F-> linear scan
        logical :: heat_s     !T-> heating calc; F-> no heating
        logical :: eigen_s    !T-> eigen calc;   F-> no eigen
        logical :: tensor_s   !T-> output tensor; F-> keep silent
        integer :: type_s     !Type of parameter scan
        integer :: style_s    !Global (0) or species (1 to nspec) parameter scan
        integer :: n_scan     !Number of steps
        integer :: n_res      !scan resolution
        real :: diff          !step size for scanned parameter
     end type scanner
     !-=-=-=-=-=-=-=-=-
     !Defines nature of parameter scans:
     !Style: -1- Global two component Scan:
     !     Type: 0 k_0-> k_1           
     !           1 theta_0 -> theta_1
     !           2 k_fixed angle
     !Style: 0- Global Scan:
     !     Type: 0 kperp
     !           1 kpar
     !           2 betap
     !           3 vtp
     !Style: 1-nspec-> Species Scan:
     !     Type: 0 tau_s
     !           1 mu_s
     !           2 alph_s
     !           3 Q_s
     !           4 D_s
     !           5 vv_s
     !-=-=-=-=-=-=-=-=-    
  type (scanner), dimension (:), allocatable :: scan
  type (scanner), dimension (:), allocatable :: rad_scan

  real, pointer :: sw,sw2,sw3,sw4    !Sweep parameter value

 !Susceptibility and elements in tensor form of wave equation
  complex, dimension(:,:,:), allocatable:: susc   !Susceptibility tensor
  complex, dimension(:,:,:,:), allocatable:: susc_low   !low-n Susceptibility tensor
  logical :: low_n=.true. !toggle on low-n susceptibility
  !>>>GGH: 1/18/23
  logical :: new_low_n=.true. !Flag to use Revised low_n for LD/TTD separation
  !<<<GGH: 1/18/23
  complex, dimension(3,3) :: lam                  !Matrix in Wave equation

  !Option
  ! Selection for the type of dispersion calculation to be undertaken
  ! Values for option are listed in README.cases
  integer :: option
  logical :: writeOut
  character(100) :: dataName,outputName,print_Name !Output file names

  !T-> use mapping routine
  !F-> read nroot_max initial guesses for (om,gam) solutions
  logical :: use_map
  !Set log or linear spacing for dispersion root search
  logical :: loggridw, loggridg
  !Bounds on read (om) and imaginary (gam) dispersion root search
  real :: omi,omf,gami,gamf
  logical :: positive_roots=.false.
  !Number of grid point in complex frequency space
  integer, parameter :: nr=128     !Number points along Re axis: 128 is a good value
  integer, parameter :: ni=128     !Number points along Im axis: 128 is a good value
  integer, parameter :: numroots=500  !Number of minima to keep
  !numerical calculation constants
  integer, parameter :: nbesmax=15 !maximum bessel sum counter 1 (from -nbesmax, to nbesmax) used when calcuating fs1 in JET-PLUME
                                   !This has no impact on how the dispersion relation is calculated!!!
                                   !rule of thumb: As j_n(b) is small for b<n/2, nbesmax/2 should be greater than or equal to b_s,max = |(kperp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)| for all species

  !variables for radial scan of solar wind models
  !rad_spec (species radial parameters) defined above
  integer :: nRad !number of radial points
  character(100) :: modelName !Input file name for radial model
  !Radial distance from the Sun, in Rs
  real, dimension (:), allocatable :: radius
  !Model for beta_||p and vtp, as a function of radial distance (Rs)shin
  real, dimension (:), allocatable :: beta_rad, vtp_rad
  !turns on/off heating and eigenfunction diagnostics for
  !radial scan
  logical :: radial_heating, radial_eigen
  !determines k values to include in radial scan
  integer :: k_scan 
  !0: single point in kperp, kpar space
  !1: fixed kperp, scan over kpar
  !2: fixed kpar,  scan over kperp
  !3: fixed theta, scan over k
  !4: fixed k, scan over theta
  !5: plane scan over (kperp, kpar)
  !6: plane scan over (k, theta)

  real :: pi

  !variables for fpc
  real    :: vperpmin,vperpmax,vparmin,vparmax   !upper and lowerbounds of normalized velocity (v/vts) space samples (gyro coords)
  real    :: vxmin,vxmax,vymin,vymax,vzmin,vzmax !upper and lowerbounds of normalized velocity (v/vts) space samples (cart coords)
  real    :: delv                                !normalized velocity (delv/vts) space grid spacing
  real    :: elecdircontribution                 !Sets components of Electric field (0 (DEFAULT) (or any other value) = Do not modify, 1=Keep only Ex(i.e.Eperp1), 2=Keep only Ey(i.e.Eperp2), 3=Keep only Ez(i.e.Epar))


  integer :: nroots                              !Number of roots found
  real, dimension(1:2,1:numroots) :: wroots      !Omega space of roots
  integer :: nroot_max                           !input specified nroots

  public :: betap,kperp,kpar,vtp,nspec,spec,susc,lam,option,writeOut
  public :: loggridw,loggridg,omi,omf,gami,gamf,nr,ni,numroots,nroot_max,nbesmax
  public :: nroots, wroots,dataName,outputName,use_map,print_Name
  public :: nscan,scan,sw,sw2,sw3,sw4, k_scan, rad_scan, positive_roots
  public :: nRad,modelName,rad_spec,radius, beta_rad, vtp_rad
  public :: radial_heating, radial_eigen, pi
  public :: low_n, susc_low
  public :: vperpmin,vperpmax,vparmin,vparmax,delv
  public :: vxmin,vxmax,vymin,vymax,vzmin,vzmax
  public :: elecdircontribution
  public :: new_low_n

end module vars
