!=============================================================================!
!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kgklein@arizona.edu                                                        !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!!*VARIABLES                                                                *!!
!=============================================================================!
!=============================================================================!
module vars
  !! Contains all global variables.
  implicit none
  private

  !INPUT PARAMETERS
  !The following parameters will be read in from 'scripts/*.in'
  !Global Input Parameters

  real, target :: betap
  !!Reference Species Parallel Thermal-to-Magnetic Pressure Ratio.
  !!8 pi n_ref T_ref,par/B^2
  
  real, target :: kperp
  !!Wavenumber perpendicular to B0, normalized to reference gyroradius.
  !!k_perp rho_ref = k_perp w_perp,ref/Omega_ref
  
  real, target :: kpar      
  !!Wavenumber parallel to B0, normalized to reference gyroradius.
  !!k_perp rho_ref = k_perp w_t,ref,perp/Omega_ref
  
  real, target :: vtp
  !!Parallel reference thermal velocity normalized to the speed of light.
  !!v_t,ref,par/c

  integer :: nspec=3
  !!Number of species/components to be included in calculation.
  
  integer :: nscan = 0
  !!Number of parameter scans.

  public :: specie
  type :: specie
     !!Species Input Parameters
     
     real :: tau_s
     !!Relative Temperature ratio.
     !!\(T_{ref}/T_{s}|_{\parallel}\)

     real :: mu_s
     !!Relative Mass ratio.
     !!\(m_{ref}/m_{s}\)

     real :: alph_s
     !!Temperature Anisotropy.
     !!\(T_{\perp}/T_{\parallel}_s\)

     real :: Q_s
     !!Relative charge ratio.
     !!\(q_{ref}/q_{s}\)

     real :: D_s
     !!Density Ratio.
     !!\(n_{s}/n_{ref}\)

     real :: vv_s
     !!Relative Drift, normalized to reference Alfven velocity
     !!\(v_{drift}/v_{A,ref}\)
     !! with \(v_{A,ref} = B/\sqrt{4 \pi n_{ref} m_{ref}}\).

  end type specie

  type (specie), dimension (:), allocatable, target :: spec
  !! Dimensionless Species/Component Parameters.

  type (specie), dimension (:,:), allocatable, target :: rad_spec
  !! Array for Varying Species/Component Parameters (under development).

  public :: scanner
  type :: scanner
     !!Parameters to control parameter scans.

     real :: range_i
     !!Initial value of scanned parameter.

     real :: range_f
     !!Final value of scanned parameter.

     logical :: log_scan
     !! Linear or Logarithmic scan.
     !!T-> log, F-> linear scan

     logical :: heat_s
     !! Controls supplementary heating calculation.
     !!T-> Heating Calculation; F-> No heating calculation.

     logical :: Eugene_s
     !! Controls supplementary eigenfunction calculation.
     !!T-> Eigenfunction calculation;   F-> No eigenfunction Calculation.

     logical :: tensor_s
     !! Controls supplementary output of susceptibility tensor.
     !!T-> Output tensor; F-> Supress output.

     integer :: type_s

     !!Defines nature of parameter scans.
     !! Style: -1- Global two component Scan:
     !!     Type: 0 k_0-> k_1           
     !!           1 theta_0 -> theta_1
     !!           2 k_fixed angle
     !!Style: 0- Global Scan:
     !!     Type: 0 kperp
     !!           1 kpar
     !!           2 betap
     !!           3 vtp
     !!Style: 1-nspec-> Species Scan:
     !!     Type: 0 tau_s
     !!           1 mu_s
     !!           2 alph_s
     !!           3 Q_s
     !!           4 D_s
     !!           5 vv_s

     integer :: style_s
     !! Defines nature of parameter scan.
     !!-1: Global two-component scan
     !! 0: Global one-component
     !! 1 to nspec: species specific parameter scan
     
     integer :: n_scan
     !!Number of output steps.
     !!n_scan*n_res Total steps taken.
     
     integer :: n_res
     !!Scan resolution between output steps.    
     !!n_scan*n_res Total steps taken.
     
     real :: diff
     !!Step size for scanned parameter.
     !! Either (swf-swi)/(n_scan*n_res)
     !! or
     !! (log10(swf)-log10(swi))/(n_scan*n_res).
     
  end type scanner

     !-=-=-=-=-=-=-=-=-    
  type (scanner), dimension (:), allocatable :: scan
  !! Array of scan parameters for all scans to be calculated.
  
  type (scanner), dimension (:), allocatable :: rad_scan
  !! Array of scan parameters for extended parameter scans (under development).

  real, pointer :: sw,sw2,sw3,sw4
  !!Parameter Sweep parameter values.

  !Susceptibility and elements in tensor form of wave equation

  complex, dimension(:,:,:), allocatable:: susc
  !! Susceptibility tensor.
  !! (1:nspec,1:3,1:3) with the 3x3 subarray arranged as:
  !! (1,1) xx; (1,2) xy; (1,3) xz;
  !! (2,1) yx; (2,2) yy; (2,3) yz;
  !! (3,1) zx; (3,2) zy; (3,3) zz;
  
  complex, dimension(:,:,:,:), allocatable:: susc_low
  !! low-n components of the susceptibility tensor.
  !! (1:nspec,1:3,1:3,0:1) with the 3x3 subarray arranged as susc.
  !! The final index contains the n=0 and n=\pm 1 contributions.
  
  logical :: low_n=.false.
  !!Toggle on low-n susceptibility suplementary calculation.
  
  logical :: new_low_n=.false. 
  !!GGH: 1/18/23
  !!Flag to use Revised low_n for LD/TTD separation.
  
  complex, dimension(3,3) :: lam
  !!Matrix in Wave equation.

  integer :: option
  !! Selection for the type of dispersion calculation to be undertaken.
  !!-1: Calculate disp(om) at a single (omega, gamma)
  !! 0: Calculate Roots for input plasma parameters.
  !! 1: Calculate Roots for input plasma parameters or Reads in root value
  !!    and then scan over plasma parameters, with range and type specified in *.in file.
  !! 2: Calculate Roots for input plasma parameters or Reads in root value
  !!    and then scan over two-dimensional plasma parameters space,
  !!    with range and type specified in *.in file.
  !! 3: deprecated
  !! 4: Make multiple maps of complex frequency space.
  !! 5: Find roots for parameters along a prescribed path
  !!    Path is set by solar wind models, with values calculated and
  !!    output by helper function (in development, the radial scan function.).
  
  logical :: writeOut
  !! Enables or suppressed output to screen.
  
  character(100) :: dataName
  !! Data Subdirectory where output is stored.

  character(100) :: outputName
  !! Common name string for output files.

  character(100) :: print_Name
  !! Additional string for output files.

  logical :: use_map
  !!Determines method for selecting initial solutions.
  !!T-> Use map routine to determine solutions in defined regions of complex frequency space.
  !!    Map parameters determined by &maps list in *.in file.
  !!F-> Read nroot_max initial guesses for complex frequency solutions.
  !!    Guesses determined by &guess_N lists in *.in file.
  
  logical :: loggridw
  !!Set log or linear spacing for real frequency axis of the map search.

  logical :: loggridg
  !!Set log or linear spacing for imaginary frequency axis of the map search.
  
  real :: omi
  !!Lower bound on real frequency map search axis.

  real :: omf
  !!Upper bound on real frequency map search axis.

  real :: gami
  !!Lower bound on imaginary frequency map search axis.

  real :: gamf
  !!Upper bound on imaginary frequency map search axis.
  
  logical :: positive_roots=.false.
  !! Consider all solutions (false) or only solutions with positive real frequencies (true).

  integer, parameter :: nr=128
  !!Number of grid points along real frequency axis
  
  integer, parameter :: ni=128
  !!Number of grid points along imaginary frequency axis
  
  integer, parameter :: numroots=500
  !!Maximum number of minima to keep for a further refinement from a map search.

  !Variables for radial scan of solar wind models: See Option 5.

  integer :: nRad
  !!Number of points to scan in radial models.
  
  character(100) :: modelName
  !!Input file name for radial model
  
  real, dimension (:), allocatable :: radius
  !!Radial distance from the Sun, in Rs.

  real, dimension (:), allocatable :: beta_rad, vtp_rad
  !Model for beta_||p and vtp, as a function of radial distance (Rs)shin
  
  !turns on/off heating and eigenfunction diagnostics for
  !radial scan
  logical :: radial_heating, radial_eigen

  integer :: k_scan
  !!Determines wavevector values to include in radial scan
  !!0: single point in kperp, kpar space
  !!1: fixed kperp, scan over kpar
  !!2: fixed kpar,  scan over kperp
  !!3: fixed theta, scan over k
  !!4: fixed k, scan over theta
  !!5: plane scan over (kperp, kpar)
  !!6: plane scan over (k, theta)

  real :: pi
  !! 4.*atan(1)

  integer :: nroots
  !!Number of roots found.  

  real, dimension(1:2,1:numroots) :: wroots
  !!Real and Imaginary components of each solution.

  integer :: nroot_max
  !!Input specified nroots to follow.

  public :: betap,kperp,kpar,vtp,nspec,spec,susc,lam,option,writeOut
  public :: loggridw,loggridg,omi,omf,gami,gamf,nr,ni,numroots,nroot_max
  public :: nroots, wroots,dataName,outputName,use_map,print_Name
  public :: nscan,scan,sw,sw2,sw3,sw4, k_scan, rad_scan, positive_roots
  public :: nRad,modelName,rad_spec,radius, beta_rad, vtp_rad
  public :: radial_heating, radial_eigen, pi
  public :: low_n, susc_low
  public :: new_low_n   !GGH: 1/18/23


end module vars
