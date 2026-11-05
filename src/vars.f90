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
  !!Wavenumber perpendicular to B0, in the selected input length units.
  !!rho_ref uses v_tperp,ref/|Omega_ref|; d_ref uses c/omega_p,ref.

   real, target :: kpar
  !!Wavenumber parallel to B0, in the selected input length units.

   real, target :: vtp
  !!rho_ref: v_t,ref,par/c; d_ref: v_A,ref/c.

   character(16) :: z_backend = 'faddeeva'
   !! Independent Z evaluator: faddeeva (default) or legacy.
   public :: z_backend

   character(16) :: bessel_mode = 'adaptive'
   real :: bessel_tol = 1.e-10
   integer :: bessel_nmax = 512, bessel_nfixed = 128
   logical :: bessel_diagnostics = .false.
   public :: bessel_mode,bessel_tol,bessel_nmax,bessel_nfixed,bessel_diagnostics

   character(16) :: normalization = 'rho_ref'
   !! Input convention: rho_ref (legacy) or d_ref (inertial length).

   ! Input values above remain in their selected convention. The kinetic
   ! calculation uses these derived rho_ref values; fully cold runs bypass them.
   real :: solver_betap, solver_kperp, solver_kpar, solver_vtp
   integer :: reference_scan_source = 0
   !! 0: validate reference beta; 1: beta scan; 2: reference tau scan.
   integer :: dispersion_status = 0
   !! Last cold-containing dispersion evaluation: cold kernel codes 1-4, 5=nonfinite arithmetic.

   integer :: nspec
  !!Number of species/components to be included in calculation.

   integer :: nscan = 0
   !!Number of parameter scans.

   integer :: collision_type = 0
   !! Type of collisions for system.
   !! 0: no collisions; standard PLUME
   !! 1: non-conservative krook operator,
   !!    neutral-charged collisions only.

   real, target :: Kn
   !! Knudson Number for neutral particles over length $\rho_{ref}$
   !! Kn^{-1} = n_0 \sigma \rho_{ref}
   
   public :: specie
   type :: specie
     !!Species Input Parameters

      real :: tau_s
     !!rho_ref: T_parallel,ref/T_parallel,s.
     !!d_ref: 2*T_parallel,s/(m_ref*v_A,ref^2).

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

      real :: nu_ns
      !! Neutral-species collision frequency
      !! \nu_{n,s}/Omega_{ref} =
      !! n_0 \sigma \sqrt{T_{\parallel,s}/m_s}/Omega_{ref}
      !! where we define Kn^{-1}= n_0 \sigma \rho_{ref}
      !! allowing us to write
      !! \nu_{n,s}/\Omega_{ref} =
      !!
      !! Note: this quantity is derived from
      !! Kn (a global parameter)
      !! and the other species parameters.
      !! The user does not input this value.

   end type specie

   type(specie), dimension(:), allocatable, target :: spec
  !! Dimensionless Species/Component Parameters.

   type(specie), dimension(:), allocatable :: solver_spec

   type(specie), dimension(:, :), allocatable, target :: rad_spec
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

      logical :: eigen_s
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
   type(scanner), dimension(:), allocatable :: scan
  !! Array of scan parameters for all scans to be calculated.

   type(scanner), dimension(:), allocatable :: rad_scan
  !! Array of scan parameters for extended parameter scans (under development).

   real, pointer :: sw, sw2, sw3, sw4
  !!Parameter Sweep parameter values.

   !Susceptibility and elements in tensor form of wave equation

   complex, dimension(:, :, :), allocatable:: susc
  !! Susceptibility tensor.
  !! (1:nspec,1:3,1:3) with the 3x3 subarray arranged as:
  !! (1,1) xx; (1,2) xy; (1,3) xz;
  !! (2,1) yx; (2,2) yy; (2,3) yz;
  !! (3,1) zx; (3,2) zy; (3,3) zz;

   complex, dimension(:, :, :, :), allocatable:: susc_low
  !! low-n components of the susceptibility tensor.
  !! (1:nspec,1:3,1:3,-1:1) with the 3x3 subarray arranged as susc.
  !! The final index contains the n=0 and n=\pm 1 contributions.

   logical :: low_n = .true.
   !!Toggle on low-n susceptibility suplementary calculation.
   !! 25-09-08: Deprecated 'new_low_n' variable, as that is now the only version of the
   !! n=0, \pm 1 calculation we wish to perform.

   complex, dimension(3, 3) :: lam
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

   logical :: positive_roots = .false.
  !! Consider all solutions (false) or only solutions with positive real frequencies (true).

  !integer, parameter :: nr=128
  integer :: nr=128
  !!Number of grid points along real frequency axis
  !! Now read in by [[map_read(subroutine)]]
  
  !integer, parameter :: ni=128
  integer :: ni=128
  !!Number of grid points along imaginary frequency axis
  !! Now read in by [[map_read(subroutine)]]
  
  integer, parameter :: numroots=500

  !!Maximum number of minima to keep for a further refinement from a map search.

   integer, parameter :: nbesmax = 15
   ! !maximum bessel sum counter 1 (from -nbesmax, to nbesmax) used when calcuating fs1 in JET-PLUME. No change for dispersion relation calculation
   !This has no impact on how the dispersion relation is calculated!!!
   !rule of thumb: As j_n(b) is small for b<n/2, nbesmax/2 should be greater than or equal to b_s,max = |(kperp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)| for all species

   !Variables for radial scan of solar wind models: See Option 5.

   integer :: nRad
   !!Number of points to scan in radial models.

   integer :: nRad_sub=1
   !! substepping between nRad points.
   !! Invoke Linear Interpolation between points

   character(100) :: modelName
  !!Input file name for radial model

   real, dimension(:), allocatable :: radius
  !!Radial distance from the Sun, in Rs.

   real, dimension(:), allocatable :: beta_rad, vtp_rad
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

   !variables for fpc
   real    :: vperpmin, vperpmax, vparmin, vparmax
   ! !upper and lowerbounds of normalized velocity (v/vts) space samples (gyro coords)

   real    :: vxmin, vxmax, vymin, vymax, vzmin, vzmax
   ! !upper and lowerbounds of normalized velocity (v/vts) space samples (cart coords)

   real    :: delv
   ! !normalized velocity (delv/vts) space grid spacing

   real    :: elecdircontribution
   ! !Sets components of Electric field (0 (DEFAULT) (or any other value) = Do not modify, 1=Keep only Ex(i.e.Eperp1), 2=Keep only Ey(i.e.Eperp2), 3=Keep only Ez(i.e.Epar))

   logical :: computemoment = .true.
  !! Enables computation of fs1 moments; not recommended as it is prone to inaccuracies and is unnecessary give that analytical moments are computed. Used for verifying self consistency of the code instead.
  !! Note this integrals are *very* challenging numerically using our current integration scheme. It is common to be off by *at least* half an order of magnitude in certain domains.
  !! This is for debug only. Please use the analytic moments.

   real :: EpsilonSokhotski_Plemelj = 0.01 !should be positive!
  !! Epsilon in the Sokhotski–Plemelj theorem, which states int f(x)/(x-a) dx can be approximated using eps->0 int f(x)/(x-a+i eps) dx to 'better' handle the singularity  numerically. This value should be left as zero, only be used by advanced users, and only when the user is computing moments for comparison to the analytic form (because it requires a *very* small delta v to 'work') (Remember to recompile!)

   real :: vxshift = 0.0
  !! shifts grid along vx- used to test stiffness of fs1 in JETPLUME (DEBUG ONLY- keep as zero otherwise)

   real :: vyshift = 0.0
  !! shifts grid along vy- used to test stiffness of fs1 in JETPLUME (DEBUG ONLY- keep as zero otherwise)

   real :: vzshift = 0.0
  !! shifts grid along vz- used to test stiffness of fs1 in JETPLUME (DEBUG ONLY- keep as zero otherwise)

   logical :: useOnlyReferenceWpar = .false.
  !! compute the reference thermal velocity for reference and then 'convert' using temperature ratios, or for each species individually.
  !! Only used for numerical moments of fs1, which are only for debug. For most users, this value is irrelevant.
  !! Choice is approximatley equivalent in domain where error in integrals of fs1 for all species are small, but those domains are rare with physical values
  !! Instead, I recommend setting this to false to allow the error to be 'absorbed' in the calculation of wparS

   integer :: nroots
  !!Number of roots found.

   real, dimension(1:2, 1:numroots) :: wroots
  !!Real and Imaginary components of each solution.

   complex :: omega_val
  !!var used to pass omega to rtsec in fpc.f90

   integer :: nroot_max
  !!Input specified nroots to follow.

   public :: betap, kperp, kpar, vtp, nspec, spec, susc, lam, option, writeOut
   public :: loggridw, loggridg, omi, omf, gami, gamf, nr, ni, numroots, nroot_max, nbesmax
   public :: nroots, wroots, dataName, outputName, use_map, print_Name
   public :: nscan, scan, sw, sw2, sw3, sw4, k_scan, rad_scan, positive_roots
   public :: nRad, nRad_sub
   public :: modelName, rad_spec, radius, beta_rad, vtp_rad
   public :: radial_heating, radial_eigen, pi
   public :: low_n, susc_low
   public :: vperpmin, vperpmax, vparmin, vparmax, delv
   public :: vxmin, vxmax, vymin, vymax, vzmin, vzmax
   public :: elecdircontribution
   public :: computemoment, EpsilonSokhotski_Plemelj, omega_val
   public :: vxshift, vyshift, vzshift
   public :: collision_type, Kn
   public :: normalization, solver_betap, solver_kperp, solver_kpar
   public :: solver_vtp, solver_spec, reference_scan_source
   public :: sync_solver_state, write_normalization_metadata
   public :: has_cold_species, all_species_cold, has_cold_reference, validate_cold_reference, dispersion_status

contains

   logical function has_cold_species()
      has_cold_species = .false.
      if (allocated(spec)) has_cold_species = any(spec%tau_s == 0.)
   end function has_cold_species

   logical function all_species_cold()
      all_species_cold = .false.
      if (allocated(spec)) all_species_cold = size(spec) > 0 .and. all(spec%tau_s == 0.)
   end function all_species_cold

   logical function has_cold_reference()
      has_cold_reference=.false.
      if (allocated(spec)) then
         if (size(spec)>0) has_cold_reference=spec(1)%tau_s==0.
      endif
   end function has_cold_reference

   subroutine validate_cold_reference
      !! Native inertial state: never construct or allocate thermal solver state.
      implicit none
      integer :: i
      if (normalization /= 'd_ref') &
           error stop 'ERROR: cold-fluid species require normalization=d_ref.'
      if (betap /= 0.) &
           error stop 'ERROR: cold-fluid reference d_ref requires betap=0, equal to species_1 tauS.'
      if (vtp <= 0.) error stop 'ERROR: cold-fluid reference d_ref requires positive vtp=v_A,ref/c.'
      if (nscan>0 .and. kpar==0.) &
           error stop 'ERROR: k_parallel=0: exact perpendicular scans are outside Stage 4A.'
      if (.not. all_species_cold() .and. kpar==0.) &
           error stop 'ERROR: warm-species calculations require nonzero kpar; the warm perpendicular limit is deferred.'
      do i=1,nspec
         if (spec(i)%tau_s<0.) error stop 'ERROR: negative tauS is invalid.'
         if (spec(i)%tau_s>0.) then
            if (spec(i)%alph_s<=0.) error stop 'ERROR: warm species require positive alphS.'
         endif
         if (spec(i)%mu_s <= 0. .or. spec(i)%Q_s == 0.) &
              error stop 'ERROR: cold-fluid species require positive muS and nonzero Qs.'
         if (spec(i)%D_s < 0.) &
              error stop 'ERROR: cold-fluid species require nonnegative density Ds.'
      enddo
      if (collision_type /= 0) &
           error stop 'ERROR: cold-fluid calculations require collision_type=0; cold collisions are unsupported.'
      if (option == 6 .or. option == 7) &
           error stop 'ERROR: cold-fluid FPC diagnostics are unavailable.'
      if (option == 5 .or. (nscan > 0 .and. option/=1 .and. option/=2 .and. option/=4)) &
           error stop 'ERROR: cold-fluid parameter scans (including eigenfunctions/heating) are deferred.'
      if (option /= -1 .and. option /= 0 .and. option /= 1 .and. &
           .not.(nscan>0 .and. (option==2 .or. option==4))) &
           error stop 'ERROR: cold-fluid calculations support only option=-1, 0, or 1 without parameter scans.'
   end subroutine validate_cold_reference

   subroutine sync_solver_state
      implicit none
      integer :: i
      real :: length_ratio

      if (normalization == 'd_ref') then
         if (reference_scan_source == 1) spec(1)%tau_s = betap
         if (reference_scan_source == 2) betap = spec(1)%tau_s
      endif
      if (normalization == 'd_ref' .and. spec(1)%tau_s == 0.) &
           error stop 'ERROR: cold-fluid reference must bypass the thermal adapter.'
      if (betap <= 0. .or. vtp <= 0.) then
         write(*,'(a)') 'ERROR: positive betap and vtp are required by the thermal adapter; cold-fluid runs must bypass it.'
         error stop 1
      endif
      do i=1,nspec
         if (spec(i)%tau_s < 0.) then
            write(*,'(a,i0,a)') 'ERROR: species ',i,' has negative tauS; use positive temperature or zero for cold-fluid.'
            error stop 1
         endif
         if (spec(i)%mu_s <= 0. .or. spec(i)%Q_s == 0.) then
            write(*,'(a,i0,a)') 'ERROR: species ',i,' requires positive alphS/muS and nonzero Qs.'
            error stop 1
         endif
         if (spec(i)%tau_s>0.) then
            if (spec(i)%alph_s<=0.) error stop 'ERROR: warm species require positive alphS.'
         endif
         if (spec(i)%tau_s == 0. .and. spec(i)%D_s < 0.) &
              error stop 'ERROR: cold-fluid species require nonnegative density Ds.'
      enddo
      if (has_cold_species()) then
         if (normalization /= 'd_ref') &
              error stop 'ERROR: cold-fluid species require normalization=d_ref in Stage 2.'
         if (spec(1)%tau_s == 0.) &
              error stop 'ERROR: cold-fluid reference species unsupported; species_1 must have positive tauS.'
         if (collision_type /= 0) &
              error stop 'ERROR: cold-fluid mixtures require collision_type=0; cold collisions are unsupported.'
         if (option == 6 .or. option == 7) &
              error stop 'ERROR: cold-fluid FPC diagnostics are unavailable in Stage 2.'
         if (option == 5 .or. (nscan > 0 .and. option/=1 .and. option/=2 .and. option/=4)) &
              error stop 'ERROR: cold-fluid parameter scans (including eigenfunctions/heating) are deferred beyond Stage 2.'
         if (option /= -1 .and. option /= 0 .and. option /= 1 .and. &
           .not.(nscan>0 .and. (option==2 .or. option==4))) &
              error stop 'ERROR: mixed Stage 2 supports only option=-1, 0, or 1 without parameter scans.'
         if (nscan>0 .and. kpar==0.) &
              error stop 'ERROR: k_parallel=0: exact perpendicular scans are outside Stage 4A.'
         if (kpar == 0.) &
              error stop 'ERROR: mixed calculations require nonzero kpar; the warm perpendicular limit is not implemented.'
      endif
      if (normalization == 'd_ref') then
         if (abs(spec(1)%tau_s-betap) > 1.e-8*max(abs(betap),abs(spec(1)%tau_s))) then
            write(*,'(a)') 'ERROR: d_ref requires species_1 tauS = betap.'
            error stop 1
         endif
      endif
      if (.not. allocated(solver_spec)) allocate(solver_spec(nspec))
      solver_betap = betap
      solver_spec = spec
      if (normalization == 'd_ref') then
         length_ratio = sqrt(betap*spec(1)%alph_s)
         solver_kperp = kperp*length_ratio
         solver_kpar = kpar*length_ratio
         solver_vtp = vtp*sqrt(betap)
         do i=1,nspec
            if (spec(i)%tau_s == 0.) cycle ! Never convert a cold temperature.
            solver_spec(i)%tau_s = betap/spec(i)%tau_s
            if (collision_type > 0 .and. Kn <= 1.e25) then
               solver_spec(i)%nu_ns = sqrt(spec(i)%mu_s*spec(i)%tau_s)/(sqrt(2.)*Kn)
            endif
         enddo
      else
         solver_kperp = kperp
         solver_kpar = kpar
         solver_vtp = vtp
      endif
   end subroutine sync_solver_state

   subroutine write_normalization_metadata
      implicit none
      integer :: unit_no, i
      logical :: exists
      character(240) :: metadata_file
      metadata_file = 'data/'//trim(dataName)//'/normalization_'//trim(outputName)//'.txt'
      if (normalization /= 'd_ref') then
         ! Reusing an output name must not leave a d_ref label on rho_ref data.
         inquire(file=trim(metadata_file),exist=exists)
         if (exists) then
            open(newunit=unit_no,file=trim(metadata_file),status='old')
            close(unit_no,status='delete')
         endif
         return
      endif
      open(newunit=unit_no,file=trim(metadata_file), &
           status='replace',action='write')
      write(unit_no,'(a)') 'normalization=d_ref'
      write(unit_no,'(a)') 'kperp,kpar: k*d_ref; vtp: v_A,ref/c; tauS: 2*T_parallel,s/(m_ref*v_A,ref^2)'
      write(unit_no,'(a)') 'betap=tauS(species_1); omega: omega/Omega_ref'
      if (has_cold_species()) then
         write(unit_no,'(a)') 'EM fields: Ex=1; omega*B=(k_d/(v_A/c)) cross E where the gauge is available.'
         write(unit_no,'(a)') 'Species moments/energy/heating/FPC unavailable; eigen moment slots contain NaN.'
      endif
      if (has_cold_reference()) then
         write(unit_no,'(a)') 'reference_species=1; cold reference; no thermal reference scale is constructed.'
         write(unit_no,'(a)') 'Direct evaluation prints raw det(Lambda); thermal diagnostics are unavailable.'
         do i=1,nspec
            if (spec(i)%tau_s<=0.) cycle
            write(unit_no,'(a,i0,a,es24.16)') 'auxiliary_first_warm_species=',i,'; b_star=',spec(i)%tau_s
            exit
         enddo
      else
         write(unit_no,'(a)') 'reference_species=1; rho_ref=v_t_perp,ref/abs(Omega_ref)'
         write(unit_no,'(a)') 'FPC file headers and velocity grids retain the derived rho_ref solver convention.'
      endif
      close(unit_no)
   end subroutine write_normalization_metadata

end module vars
