!=============================================================================!
!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kgklein@arizona.edu                                                        !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!!*MISC. FUNCTIONS                                                          *!!
!=============================================================================!
!=============================================================================!
module functions
  !!Calculates Misc. functions, esp. I/O operations.
   implicit none
   private

   integer :: is
  !!Index for scan and species loops.

   integer :: unit
  !! Custom index for I/O.

   integer, parameter :: stdout_unit = 6
  !! Standard index for I/O.

   integer, save :: input_unit_no
  !! Index for reading in files.

   integer, save ::error_unit_no = stdout_unit
  !! Index for outputing error messages.

   character(50) :: runname
  !! String for input file parameters.

   public :: read_in_params, read_map_input, read_scan_input, read_guess_input
   public :: get_unused_unit, read_radial_input
   public :: solve_norm_facs

contains

   !-=-=-=-=-
   !-=-=-=-=-
   subroutine read_in_params
    !!Read in system parameters.
    !!Input file is argument after executable:
    !!$ ./plume.e example.in
      use vars, only: betap, kperp, kpar, vtp, nspec, spec, susc, option, writeOut
      use vars, only: dataName, nscan, nroot_max, use_map, outputName, pi
      use vars, only: low_n, susc_low, new_low_n
      use vars, only: vperpmin, vperpmax, vparmin, vparmax, delv
      use vars, only: vxmin, vxmax, vymin, vymax, vzmin, vzmax, elecdircontribution
      implicit none

      real :: sum_nq
    !!For Testing Quasineutrality.

      real :: sum_nqv
    !!For zero net current.

      integer :: is
    !!Index for scan and species loops.

      !Read in the basic parameter list.
      nameList /params/ &
         betap, kperp, kpar, vtp, nspec, nscan, option, nroot_max, &
         use_map, low_n, new_low_n, &
         writeOut, dataName, outputName

      nameList /fpc/ &
         vperpmin, vperpmax, vparmin, vparmax, delv, &
         vxmin, vxmax, vymin, vymax, vzmin, vzmax, elecdircontribution

      !default values
      elecdircontribution = 0.

      call get_unused_unit(input_unit_no)
      call get_runname(runname)
      runname = trim(runname)//".in"
      unit = input_unit_no
      open (unit=unit, file=runname, status='old', action='read')
      read (unit=unit, nml=params)

      !Allocate the species variable to have nspec indicies
      allocate (spec(1:nspec))

      !Allocate the susceptibility tensor to have nspec indicies
      allocate (susc(1:nspec, 3, 3))
      if (low_n) &
         allocate (susc_low(1:nspec, 3, 3, 0:1))

      !initialized quasineutrality and current conservation check.
      sum_nq = 0.; sum_nqv = 0.

      !Read in species parameters
      !This is a bit of FORTRAN black magic borrowed from AGK.
      !     which allows us to loop over iterative nml/groupnames.
      do is = 1, nspec
         call get_indexed_namelist_unit(unit, "species", is)
         call spec_read(is)
         !quasineutrality check.
         sum_nq = sum_nq + spec(is)%D_s/spec(is)%Q_s
         !net current check.
         sum_nqv = sum_nqv + spec(is)%vv_s*spec(is)%D_s/spec(is)%Q_s
         close (unit)
      end do

      !Check that ref mass ratio is 1 (in the python wrapper, this is used to help determine the output format to load in sweeps, which is why we check it)
      if (spec(1)%mu_s - 1 > 0.0000001) then
         write (*, '(a,es11.4)') &
            'ERROR: reference mass ratio is not 1:     spec(1)%mu_s =', spec(1)%mu_s
      end if

      !read in fpc params
      if (option == 6 .or. option == 7) then

         rewind (input_unit_no)
         read (unit=input_unit_no, nml=fpc)
      end if

      close (input_unit_no)

      !Notify if not Quasineutral,
      if (abs(sum_nq) .ge. 1.E-6) &
         write (*, '(a,es11.4)') &
         'ERROR: Plasma not quasineutral:     Sum [n_s Q_s] =', sum_nq
      !or not Current-free (not strictly required by linear theory but plasma is likely unstable and thus likely poorly described by linear theory)
      if (abs(sum_nqv) .ge. 1.E-6) &
         write (*, '(a,es11.4)') &
         'ERROR?/WARNING: Non zero Current:        Sum [n_s q_s V_s] =', sum_nqv
      !Take note User: Make sure the above sums are zero.
      !Note if:
      if (spec(1)%vv_s .ne. 0.) &
         write (*, '(a,es11.4)') &
         'ERROR?/WARNING: Not in refernece species rest frame:    v_par drift ref =', spec(1)%vv_s

      !You always need pi. Always. And sometimes pie too.
      pi = 4.*atan(1.)
   end subroutine read_in_params

!-=-=-=-=-
!-=-=-=-=-
   subroutine spec_read(is)
    !!Subroutine for reading in species/component parameters.
      use vars, only: spec
      implicit none
      !Passed
      integer :: is
    !!Species index.

      !Local
      !Dummy values for reading in species parameters
      real :: tauS
    !!Parallel Temperature Ratio.
    !!\(T_{ref}/T_{s}|_{\parallel}\)

      real :: muS
    !!Mass Ratio.
    !!\(m_{ref}/m_{s}\)

      real :: alphS
    !!Temperature Anisotropy.
    !!\(T_{\perp}/T_{\parallel}_s\)

      real :: Qs
    !!Relative charge ratio.
    !!\(q_{ref}/q_{s}\)

      real :: Ds
    !!Density Ratio.
    !!\(n_{s}/n_{ref}\)

      real :: vvS
    !!Relative Drift, normalized to reference Alfven velocity
    !!\(v_{drift}/v_{A,ref}\)
    !! with \(v_{A,ref} = B/\sqrt{4 \pi n_{ref} m_{ref}}\).

      nameList /species/ &
         tauS, muS, alphS, Qs, Ds, vvS
      read (unit=unit, nml=species)

      !Assign dummy values to global component array.
      spec(is)%tau_s = tauS
      spec(is)%mu_s = muS
      spec(is)%alph_s = alphS
      spec(is)%Q_s = Qs
      spec(is)%D_s = Ds
      spec(is)%vv_s = vvS

   end subroutine spec_read

!-=-=-=-=
!-=-=-=-=
   subroutine read_map_input
  !!Read in parameters for bounds on mapping dispersion roots.
  !! Invokes [[map_read(subroutine)]].

  use vars, only: loggridw,loggridg,omi,omf,gami,gamf,nr,ni
  implicit none
  !Append the .in file as first argument in executable
  call get_runname(runname)
  !ie ./plume.e system.in
  runname=trim(runname)//".in"
   
  call get_unused_unit (input_unit_no)
  unit=input_unit_no
  open (unit=unit,file=runname,status='old',action='read')
  call map_read
  close (unit)
  
end subroutine read_map_input

!-=-=-=-=-
!-=-=-=-=-
   subroutine map_read
  !!Subroutine for reading in frequency limits for map search of
  !! complex frequency solution space.

  use vars, only: loggridw,loggridg,omi,omf
  use vars, only: gami,gamf,positive_roots,nr,ni
  implicit none

  nameList /maps/ loggridw,loggridg,omi,omf,gami,gamf,&
       positive_roots,nr,ni

      read (unit=unit, nml=maps)

   end subroutine map_read

!-=-=-=-=
!-=-=-=-=
   subroutine read_scan_input
  !!Read in limits for scans in plasma parameter space.
  !! Invokes [[scan_read(subroutine)]].
      use vars, only: nscan, scan
      implicit none
      integer :: is
  !!Scan index.

      !Allocate scan type object
      allocate (scan(1:nscan))

      call get_unused_unit(input_unit_no)
      call get_runname(runname)
      runname = trim(runname)//".in"
      unit = input_unit_no
      open (unit=unit, file=runname, status='old', action='read')

      do is = 1, nscan
         call get_indexed_namelist_unit(unit, "scan_input", is)
         call scan_read(is)
         close (unit)
      end do
      close (input_unit_no)

   end subroutine read_scan_input

!-=-=-=-=-
!-=-=-=-=-
   subroutine scan_read(is)
  !!Subroutine for reading in scan parameters.
      use vars, only: scan
      implicit none
      !Passed
      integer :: is
  !!Scan index.

      !Local
      integer :: scan_type
  !!Defines kind of parameter scans.

      integer :: scan_style
  !!Defines number of components of scan.

      integer :: ns
  !!Number of output steps.

      integer :: nres
  !!Scan resolution between output steps.

      real    :: swi
  !!Initial value of scanned parameter.

      real    :: swf
  !!Final value of scanned parameter.

      logical :: swlog
  !! Linear or Logarithmic scan.

      logical :: heating
  !! Controls supplementary heating calculation.

      logical :: eigen
  !! Controls supplementary eigenfunction calculation.

      logical :: tensor
  !! Controls supplementary output of susceptibility tensor.

      nameList /scan_input/ &
         scan_type, swi, swf, swlog, scan_style, ns, nres, &
         heating, eigen, tensor
      read (unit=unit, nml=scan_input)

      !Assign dummy values to scan array.
      scan(is)%range_i = swi
      scan(is)%range_f = swf
      scan(is)%log_scan = swlog
      scan(is)%type_s = scan_type
      scan(is)%style_s = scan_style
      scan(is)%n_scan = ns
      scan(is)%n_res = nres
      scan(is)%eigen_s = eigen
      scan(is)%heat_s = heating
      scan(is)%tensor_s = tensor

   end subroutine scan_read

!-=-=-=-=
!-=-=-=-=
   subroutine read_guess_input
  !!Read in initial values for solutions.
  !! Invokes [[om_read(subroutine)]].
      use vars, only: nroot_max, wroots
      implicit none

      integer :: is
  !!Solution index

      call get_unused_unit(input_unit_no)
      call get_runname(runname)
      runname = trim(runname)//".in"
      unit = input_unit_no
      open (unit=unit, file=runname, status='old', action='read')

      do is = 1, nroot_max
         call get_indexed_namelist_unit(unit, "guess", is)
         call om_read(is)

         write (*, '(a,i0,a,2es12.4)') 'Input for Root ', is, ': ', wroots(1:2, is)
         close (unit)
      end do
      close (input_unit_no)

   end subroutine read_guess_input

!-=-=-=-=-
!-=-=-=-=-
   subroutine om_read(is)
  !!Subroutine for reading in initial guesses for complex frequencies of solutions.
      use vars, only: wroots
      implicit none
      !Passed
      integer :: is
  !!Solution Index.

      !Local
      real    :: g_om
  !!Dummy values for real frequency component of solution.

      real    :: g_gam
  !!Dummy values for imaginary frequency component of solution.

      nameList /guess/ &
         g_om, g_gam
      read (unit=unit, nml=guess)

      !Assign dummy values to frequency of solution.
      wroots(1, is) = g_om
      wroots(2, is) = g_gam

   end subroutine om_read

!-=-=-=-=
   subroutine read_radial_input
  !!Read in parameters for radial solar wind model scan.
  !!(in development).
      use vars, only: nRad, modelName, nspec, rad_spec, radius, pi
      use vars, only: beta_rad, vtp_rad, radial_heating, radial_eigen
      use vars, only: betap, vtp, spec, kperp, kpar, k_scan, rad_scan
      implicit none

      !Local
      integer :: is
  !!Species/component index.

      integer :: ir
  !!Radial index.

      character(100) :: readName
  !!String for reading parameters for radial scan.

      !Dummy values for reading in:
      real :: tau_in
  !!Parallel Temperature Ratio.
  !!\(T_{ref}/T_{s}|_{\parallel}\)

      real :: mu_in
  !!Mass Ratio.
  !!\(m_{ref}/m_{s}\)

      real :: alph_in
  !!Temperature Anisotropy.
  !!\(T_{\perp}/T_{\parallel}_s\)

      real :: Q_in
  !!Relative charge ratio.
  !!\(q_{ref}/q_{s}\)

      real :: D_in
  !!Density Ratio.
  !!\(n_{s}/n_{ref}\)

      real :: vv_in
  !!Relative Drift, normalized to reference Alfven velocity
  !!\(v_{drift}/v_{A,ref}\)
  !! with \(v_{A,ref} = B/\sqrt{4 \pi n_{ref} m_{ref}}\).

      !Append the .in file as first argument in executable
      call get_runname(runname)
      !ie ./plume.e system.in
      runname = trim(runname)//".in"

      call get_unused_unit(input_unit_no)
      unit = input_unit_no
      open (unit=unit, file=runname, status='old', action='read')
      call radial_read
      close (unit)

      write (*, *) 'Read in Radial Parameters'
      !Allocate the radial parameter scan for each species
      allocate (rad_spec(1:nspec, 0:nRad))
      !Allocate Radial Distance from Sun
      allocate (radius(0:nRad)); radius = 0.

  !Read in radial parameters for each species
  do is = 1, nspec
     write(readName,'(a,i0,a)')&
          trim(modelName),is,".rad"
     call get_unused_unit (input_unit_no)
     unit=input_unit_no
     open(unit = unit, file=trim(readName) , status='old', action='read')
     do ir = 0,nRad-1
        read(unit,*)&
             radius(ir),&
             tau_in, mu_in, alph_in, Q_in, D_in, vv_in
        !Assign Dummy Variables.
        rad_spec(is,ir)%tau_s = tau_in
        rad_spec(is,ir)%mu_s = mu_in
        rad_spec(is,ir)%alph_s = alph_in
        rad_spec(is,ir)%Q_s = Q_in
        rad_spec(is,ir)%D_s = D_in
        rad_spec(is,ir)%vv_s = vv_in
     enddo
     close(unit)
  enddo

      !Allocate Global Plasma Parameters
      allocate (beta_rad(0:nRad)); beta_rad = 0.
      allocate (vtp_rad(0:nRad)); vtp_rad = 0.

  !Read in Global Parameters, beta_||p and vtp
  write(readName,'(a,i0,a)')&
       trim(modelName),0,".rad"
  call get_unused_unit (input_unit_no)
  unit=input_unit_no
  open(unit = unit, file=trim(readName) , status='old', action='read')
  do ir = 0,nRad-1
     read(unit,*)&
          radius(ir),beta_rad(ir), vtp_rad(ir)          
  enddo
  close(unit)

      !Set parameters to values at intial radial position
      betap = beta_rad(0)
      vtp = vtp_rad(0)

  do is = 1, nspec
     spec(is)%tau_s  = rad_spec(is,0)%tau_s
     spec(is)%mu_s   = rad_spec(is,0)%mu_s
     spec(is)%alph_s = rad_spec(is,0)%alph_s
     spec(is)%Q_s    = rad_spec(is,0)%Q_s
     spec(is)%D_s    = rad_spec(is,0)%D_s
     spec(is)%vv_s    = rad_spec(is,0)%vv_s
  enddo
  
  call get_runname(runname)
  runname=trim(runname)//".in"
  
  !Determine K range which will be explored
  select case(k_scan)
  case(0)
     !simple case of constant kperp, kpar
     write(*,'(a)')'Parametric Scan with fixed (kperp,kpar)'
     call radial_read_0
  case(1)
     !fixed kperp, scan over kpar
     write(*,'(a)')'Parametric Scan with fixed kperp, varying kpar'
     call radial_read_1
  case(2)
     !fixed kpar,  scan over kperp
     write(*,'(a)')'Parametric Scan with varying kperp, fixed kpar'
     call radial_read_2
  case(3)
     !fixed theta, scan over k
     write(*,'(a)')'Parametric Scan with varying |k|, fixed theta'
     call radial_read_3
  case(4)
     !fixed k, scan over theta
     write(*,'(a)')'Parametric Scan with fixed |k|, varying theta'
     call radial_read_4
  case(5)
     !plane scan over (kperp, kpar)
     write(*,'(a)')'Parametric Scan with varying kperp and kpar'
     call radial_read_5
  case(6)
     !plane scan over (k, theta)
     write(*,'(a)')'Parametric Scan with varying |k| and theta'
     call radial_read_6
  end select

      call get_runname(runname)
      runname = trim(runname)//".in"

      !Determine K range which will be explored
      select case (k_scan)
      case (0)
         !simple case of constant kperp, kpar
         call radial_read_0
      case (1)
         !fixed kperp, scan over kpar
         call radial_read_1
      case (2)
         !fixed kpar,  scan over kperp
         call radial_read_2
      case (3)
         !fixed theta, scan over k
         call radial_read_3
      case (4)
         !fixed k, scan over theta
         call radial_read_4
      case (5)
         !plane scan over (kperp, kpar)
         call radial_read_5
      case (6)
         !plane scan over (k, theta)
         call radial_read_6
      end select

   end subroutine read_radial_input

!-=-=-=-=-
!-=-=-=-=-
   subroutine radial_read
  !!Subroutine for reading in radial scan global parameters
      use vars, only: nRad, modelName, radial_heating, radial_eigen, k_scan
      implicit none

  nameList /radial_input/ nRad, modelName, &
       radial_heating, radial_eigen, k_scan
  read (unit=unit,nml=radial_input)

   end subroutine radial_read

!-=-=-=-=-
   subroutine radial_read_0
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(k_{\perp}\) and \(k_{\parallel}\).
      use vars, only: kperp, kpar
      implicit none

      real :: kperp_1
  !! Dummy variable for \(k_{\perp}\) readin.

      real :: kpar_1
  !! Dummy variable for \(k_{\parallel}\) readin.

      nameList /k_range/ kperp_1, kpar_1
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)

      !Set initial wavevector.
      kperp = kperp_1
      kpar = kpar_1

   end subroutine radial_read_0

!-=-=-=-=-
   subroutine radial_read_1
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(k_{\perp}\) and varying \(k_{\parallel}\).
      use vars, only: kperp, kpar, rad_scan
      implicit none
      real :: kperp_1
  !! Dummy variable for \(k_{\perp}\) readin.

      real :: kpar_1
  !! Dummy variable for \(k_{\parallel}\) initial value.

      real :: kpar_2
  !! Dummy variable for \(k_{\parallel}\) final value.

      integer :: nK
  !! Number of output \(k_{\parallel}\) values.

      integer :: kres
  !! Number of steps between output \(k_{\parallel}\) values.

      logical :: rad_log
  !! Logirithmic or linear \(k_{\parallel}\) scan.

      logical :: rad_heat
  !! True turns on supplemental heating calculation.

      logical :: rad_eigen
  !! True turns on supplemental eigenfunction calculation.

      !Read in parameters.
      nameList /k_range/ kperp_1, kpar_1, kpar_2, nK, kres, &
         rad_log, rad_heat, rad_eigen
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)
      !Set initial wavevector
      kperp = kperp_1
      kpar = kpar_1
      !Set \(k_{\parallel}\) range.
      allocate (rad_scan(1))
      rad_scan(1)%range_i = kpar_1
      rad_scan(1)%range_f = kpar_2
      rad_scan(1)%n_scan = nk
      rad_scan(1)%n_res = kres
      rad_scan(1)%heat_s = rad_heat
      rad_scan(1)%eigen_s = rad_eigen
      rad_scan(1)%log_scan = rad_log
      if (rad_log) then
         rad_scan(1)%diff = &
            (log10(kpar_2) - log10(kpar_1))/ &
            real(nk*kres)
      else
         rad_scan(1)%diff = &
            ((kpar_2) - (kpar_1))/ &
            real(nk*kres)
      end if

   end subroutine radial_read_1

!-=-=-=-=-
   subroutine radial_read_2
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(k_{\parallel}\) and varying \(k_{\perp}\).
      use vars, only: kperp, kpar, rad_scan
      implicit none
      real :: kpar_1
  !! Dummy variable for fixed \(k_{\parallel}\).

      real :: kperp_1
  !! Dummy variable for \(k_{\perp}\) initial value.

      real :: kperp_2
  !! Dummy variable for \(k_{\perp}\) final value.

      integer :: nK
  !! Number of output \(k_{\perp}\) values.

      integer :: kres
  !! Number of steps between output \(k_{\perp}\) values.

      logical :: rad_log
  !! Logirithmic or linear \(k_{\perp}\) scan.

      logical :: rad_heat
  !! True turns on supplemental heating calculation.

      logical :: rad_eigen
  !! True turns on supplemental eigenfunction calculation.

      !fixed kpar,  scan over kperp
      nameList /k_range/ kperp_1, kperp_2, kpar_1, nK, kres, &
         rad_log, rad_heat, rad_eigen
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)
      kperp = kperp_1
      kpar = kpar_1
      allocate (rad_scan(1))
      rad_scan(1)%range_i = kperp_1
      rad_scan(1)%range_f = kperp_2
      rad_scan(1)%n_scan = nk
      rad_scan(1)%n_res = kres
      rad_scan(1)%heat_s = rad_heat
      rad_scan(1)%eigen_s = rad_eigen
      rad_scan(1)%log_scan = rad_log
      if (rad_log) then
         rad_scan(1)%diff = &
            (log10(kperp_2) - log10(kperp_1))/ &
            real(nk*kres)
      else
         rad_scan(1)%diff = &
            ((kperp_2) - (kperp_1))/ &
            real(nk*kres)
      end if

   end subroutine radial_read_2

!-=-=-=-=-
   subroutine radial_read_3
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(\theta\) and varying \(|k|\).
      use vars, only: kperp, kpar, rad_scan, pi
      implicit none
      real :: k_1
  !! Dummy variable for \(|k|\) initial value.

      real :: k_2
  !! Dummy variable for \(|k|\) final value.

      real :: theta_1
  !! Dummy variable for \(\theta\) fixed value (in deg.)

      integer :: nK
  !! Number of output \(|k|\) values.

      integer :: kres
  !! Number of steps between output \(|k|\) values.

      logical :: rad_log
  !! Logirithmic or linear \(|k|\) scan.

      logical :: rad_heat
  !! True turns on supplemental heating calculation.

      logical :: rad_eigen
  !! True turns on supplemental eigenfunction calculation.

      !Fixed theta, scan over wavevector amplitude.
      nameList /k_range/ k_1, k_2, theta_1, nK, kres, &
         rad_log, rad_heat, rad_eigen
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)
      !Set initial wavevector.
      kperp = k_1*sin(theta_1*pi/180.)
      kpar = k_1*cos(theta_1*pi/180.)
      allocate (rad_scan(1))
      rad_scan(1)%range_i = k_1
      rad_scan(1)%range_f = k_2
      rad_scan(1)%n_scan = nk
      rad_scan(1)%n_res = kres
      rad_scan(1)%heat_s = rad_heat
      rad_scan(1)%eigen_s = rad_eigen
      rad_scan(1)%log_scan = rad_log
      if (rad_log) then
         rad_scan(1)%diff = &
            (log10(k_2) - log10(k_1))/ &
            real(nk*kres)
      else
         rad_scan(1)%diff = &
            ((k_2) - (k_1))/ &
            real(nk*kres)
      end if

   end subroutine radial_read_3

!-=-=-=-=-
   subroutine radial_read_4
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(|k|\) and varying \(\theta\).
      use vars, only: kperp, kpar, rad_scan, pi
      implicit none

      real :: k_1
  !! Dummy variable for \(|k|\) fixed value.

      real :: theta_1
  !! Dummy variable for \(\theta\) initial value (in deg.).

      real :: theta_2
  !! Dummy variable for \(\theta\) final value (in deg.).

      integer :: nK
  !! Number of output \(k_{\parallel}\) values.

      integer :: kres
  !! Number of steps between output \(k_{\parallel}\) values.

      logical :: rad_log
  !! Logirithmic or linear \(k_{\parallel}\) scan.

      logical :: rad_heat
  !! True turns on supplemental heating calculation.

      logical :: rad_eigen
  !! True turns on supplemental eigenfunction calculation.

      !fixed k, scan over theta
      nameList /k_range/ k_1, theta_1, theta_2, nK, kres, &
         rad_log, rad_heat, rad_eigen
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)
      !set initial k
      kperp = k_1*sin(theta_1*pi/180.)
      kpar = k_1*cos(theta_1*pi/180.)
      allocate (rad_scan(1))
      rad_scan(1)%range_i = theta_1
      rad_scan(1)%range_f = theta_2
      rad_scan(1)%n_scan = nk
      rad_scan(1)%n_res = kres
      rad_scan(1)%heat_s = rad_heat
      rad_scan(1)%eigen_s = rad_eigen
      rad_scan(1)%log_scan = rad_log
      if (rad_log) then
         rad_scan(1)%diff = &
            (log10(theta_2) - log10(theta_1))/ &
            real(nk*kres)
      else
         rad_scan(1)%diff = &
            ((theta_2) - (theta_1))/ &
            real(nk*kres)
      end if

   end subroutine radial_read_4

!-=-=-=-=-
   subroutine radial_read_5
  !!Subroutine for reading in radial scan parameters
  !!for 2D scan over \(k_{\perp}\) and \(k_{\parallel}\).
      use vars, only: kperp, kpar, rad_scan
      implicit none

      real :: kperp_1
  !! Dummy variable for \(k_{\perp}\) initial value.

      real :: kperp_2
  !! Dummy variable for \(k_{\perp}\) final value.

      real :: kpar_1
  !! Dummy variable for \(k_{\parallel}\) initial value.

      real :: kpar_2
  !! Dummy variable for \(k_{\parallel}\) final value.

      integer :: nK
  !! Number of output \(k_{\perp}\) values.

      integer :: kres
  !! Number of steps between output \(k_{\perp}\) values.

      integer :: nK2
  !! Number of output \(k_{\parallel}\) values.

      integer :: kres2
  !! Number of steps between output \(k_{\parallel}\) values.

      logical :: rad_log_perp
  !! Logirithmic or linear \(k_{\perp}\) scan.

      logical :: rad_log_par
  !! Logirithmic or linear \(k_{\parallel}\) scan.

      logical :: rad_heat
  !! True turns on supplemental heating calculation.

      logical :: rad_eigen
  !! True turns on supplemental eigenfunction calculation.

      nameList /k_range/ kperp_1, kperp_2, kpar_1, kpar_2, &
         nk, kres, nk2, kres2, rad_log_perp, rad_log_par, rad_heat, rad_eigen
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)

      !Set initial wavevector.
      kperp = kperp_1
      kpar = kpar_1
      allocate (rad_scan(1:2))
      !Perpendicular scan parameters.
      rad_scan(1)%range_i = kperp_1
      rad_scan(1)%range_f = kperp_2
      rad_scan(1)%n_scan = nk
      rad_scan(1)%n_res = kres
      rad_scan(1)%heat_s = rad_heat
      rad_scan(1)%eigen_s = rad_eigen
      rad_scan(1)%log_scan = rad_log_perp
      if (rad_log_perp) then
         rad_scan(1)%diff = &
            (log10(kperp_2) - log10(kperp_1))/ &
            real(nk*kres)
      else
         rad_scan(1)%diff = &
            ((kperp_2) - (kperp_1))/ &
            real(nk*kres)
      end if

      !Parallel scan parameters.
      rad_scan(2)%range_i = kpar_1
      rad_scan(2)%range_f = kpar_2
      rad_scan(2)%n_scan = nk2
      rad_scan(2)%n_res = kres2
      rad_scan(2)%heat_s = rad_heat
      rad_scan(2)%eigen_s = rad_eigen
      rad_scan(2)%log_scan = rad_log_par
      if (rad_log_par) then
         rad_scan(2)%diff = &
            (log10(kpar_2) - log10(kpar_1))/ &
            real(nk2*kres2)
      else
         rad_scan(2)%diff = &
            ((kpar_2) - (kpar_1))/ &
            real(nk2*kres2)
      end if

   end subroutine radial_read_5

!-=-=-=-=-
   subroutine radial_read_6
  !!Subroutine for reading in radial scan parameters
  !!for 2D scan over \(|k|\) and \(\theta\).
      use vars, only: kperp, kpar, rad_scan, pi
      implicit none

      real :: k_1
  !! Dummy variable for \(|k|\) initial value.

      real :: k_2
  !! Dummy variable for \(|k|\) final value.

      real :: theta_1
  !! Dummy variable for \(\theta\) initial value.

      real :: theta_2
  !! Dummy variable for \(\theta\) final value.

      integer :: nK
  !! Number of output \(|k|\) values.

      integer :: kres
  !! Number of steps between output \(|k|\) values.

      integer :: ntheta
  !! Number of output \(\theta\) values.

      integer :: thetares
  !! Number of steps between output \(\theta\) values.

      logical :: rad_log_k
  !! Logirithmic or linear \(|k|\) scan.

      logical :: rad_log_theta
  !! Logirithmic or linear \(\theta\) scan.

      logical :: rad_heat
  !! True turns on supplemental heating calculation.

      logical :: rad_eigen
  !! True turns on supplemental eigenfunction calculation.

      !Read in parameters.
      nameList /k_range/ k_1, k_2, theta_1, theta_2, &
         nk, kres, ntheta, thetares, rad_log_k, rad_log_theta, rad_heat, rad_eigen
      open (unit=4, file=runname, status='old', action='read')
      read (unit=4, nml=k_range)
      close (4)

      !Set initial wavevectors.
      kperp = k_1*sin(theta_1*pi/180.)
      kpar = k_1*cos(theta_1*pi/180.)
      allocate (rad_scan(1:2))
      !theta scan parameters.
      rad_scan(1)%range_i = theta_1
      rad_scan(1)%range_f = theta_2
      rad_scan(1)%n_scan = ntheta
      rad_scan(1)%n_res = thetares
      rad_scan(1)%heat_s = rad_heat
      rad_scan(1)%eigen_s = rad_eigen
      rad_scan(1)%log_scan = rad_log_theta
      if (rad_log_theta) then
         rad_scan(1)%diff = &
            (log10(theta_2) - log10(theta_1))/ &
            real(ntheta*thetares)
      else
         rad_scan(1)%diff = &
            ((theta_2) - (theta_1))/ &
            real(ntheta*thetares)
      end if

      !k scan
      rad_scan(2)%range_i = k_1
      rad_scan(2)%range_f = k_2
      rad_scan(2)%n_scan = nk
      rad_scan(2)%n_res = kres
      rad_scan(2)%heat_s = rad_heat
      rad_scan(2)%eigen_s = rad_eigen
      rad_scan(2)%log_scan = rad_log_k
      if (rad_log_k) then
         rad_scan(2)%diff = &
            (log10(k_2) - log10(k_1))/ &
            real(nk*kres)
      else
         rad_scan(2)%diff = &
            ((k_2) - (k_1))/ &
            real(nk*kres)
      end if

   end subroutine radial_read_6

!-=-=-=-=
!-=-=-=-=
   subroutine get_runname(runname)
  !! Get runname for output files from input argument
  !! by trimming '.in'.
      implicit none
      integer       :: l
    !! Length of argument.
      character(50) :: arg
    !! Argument after executable.
      character(50), intent(out) :: runname
    !! Argument trimmed of '.in' string.

      !Get the first argument of the program execution command
      call getarg(1, arg)

      !Check if this is the input file and trim .in extension to get runname
      l = len_trim(arg)
      if (l > 3 .and. arg(l - 2:l) == ".in") then
         runname = arg(1:l - 3)
      end if
   end subroutine get_runname
!------------------------------------------------------------------------------

!-=-=-=-=-=-
!The following routines:
!    get_indexed_namelist_unit
!    input_unit_exist
!    get_unused_unit
!    input_unit
!were all adopted from the Astrophysical Gyrokinetic Code (AGK)
!as a means of allowing arbitrary namelist group name input.
!A bit of hassle, but worth the effort.
!-=-=-=-=-=-
   subroutine get_indexed_namelist_unit(unit, nml, index_in)
    !!Extract namelist.
      implicit none
      integer, intent(out) :: unit
      character(*), intent(in) :: nml
      integer, intent(in) :: index_in
      character(500) :: line
      integer :: iunit, iostat, in_file
      integer :: ind_slash
      logical :: exist

      call get_unused_unit(unit)
      ind_slash = index(runname, "/", .True.)
      if (ind_slash .EQ. 0) then !No slash in name
         !Original behaviour
         open (unit=unit, file="."//trim(runname)//".scratch")
      else
         !General behaviour
         open (unit=unit, file=trim(runname(1:ind_slash))//"."//trim(runname(ind_slash + 1:))//".scratch")
      end if

      write (line, *) index_in
      line = nml//"_"//trim(adjustl(line))
      in_file = input_unit_exist(trim(line), exist)

      if (exist) then
         iunit = input_unit(trim(line))
      else
         write (6, *) "get_indexed_namelist: following namelist not found ", trim(line)
         return
      end if

      read (unit=iunit, fmt="(a)") line
      write (unit=unit, fmt="('&',a)") nml

      do
         read (unit=iunit, fmt="(a)", iostat=iostat) line
         if (iostat /= 0 .or. trim(adjustl(line)) == "/") exit
         write (unit=unit, fmt="(a)") trim(line)
      end do
      write (unit=unit, fmt="('/')")
      rewind (unit=unit)
   end subroutine get_indexed_namelist_unit

   !-=-=-=-=-=-
   !-=-=-=-=-=-
   function input_unit_exist(nml, exist)
    !!Is a namelist already open?
      implicit none
      character(*), intent(in) :: nml
      logical, intent(out) :: exist
      integer :: input_unit_exist, iostat
      character(500) :: line
      intrinsic adjustl, trim
      input_unit_exist = input_unit_no
      exist = .true.
      if (input_unit_no > 0) then
         rewind (unit=input_unit_no)
         do
            read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
            if (iostat /= 0) then
               rewind (unit=input_unit_no)
               exit
            end if
            if (trim(adjustl(line)) == "&"//nml) then
               backspace (unit=input_unit_no)
               return
            end if
         end do
      end if
      exist = .false.
   end function input_unit_exist

   !-=-=-=-=-=-
   !-=-=-=-=-=-
   function input_unit(nml)
    !!Determine I/O unit.
      implicit none
      character(*), intent(in) :: nml
      integer :: input_unit, iostat
      character(500) :: line
      intrinsic adjustl, trim
      input_unit = input_unit_no
      if (input_unit_no > 0) then
         rewind (unit=input_unit_no)
         do
            read (unit=input_unit_no, fmt="(a)", iostat=iostat) line
            if (iostat /= 0) then
               rewind (unit=input_unit_no)
               exit
            end if
            if (trim(adjustl(line)) == "&"//nml) then
               backspace (unit=input_unit_no)
               return
            end if
         end do
      end if
      write (unit=error_unit_no, fmt="('Couldn''t find namelist: ',a)") nml
      write (unit=*, fmt="('Couldn''t find namelist: ',a)") nml
   end function input_unit

   subroutine get_unused_unit(unit)
    !!Find a I/O unit that is not currently open.
      implicit none
      integer, intent(out) :: unit
      logical :: od
      unit = 50
      do
         inquire (unit=unit, opened=od)
         if (.not. od) return
         unit = unit + 1
      end do
   end subroutine get_unused_unit
!-=-=-=-=-=-

   !Due to our definition of dimensionless, f1s (which is not the same as performing u sub on v to write in terms of vprime/v, although its close) we
   !   get a different susc tensor and eigenmodes. They are related and can be converted. This function does that.
   ! ======================================================================
   subroutine solve_norm_facs(sigma, a, nval, fac, Efac, ierr)
      ! ======================================================================
      ! Solves the system of equations:
      !   fac(i,j) * Efac(k) = sigma(i,j)  (where k depends on j)
      ! Subject to constraints:
      !   fac(2,1) = -fac(1,2)  (fyx = -fxy)
      !   fac(3,1) =  fac(1,3)  (fzx =  fxz)
      !   fac(3,2) = -fac(2,3)  (fzy = -fyz)
      ! And the normalization condition:
      !   a(1)*Efac(1) + a(2)*Efac(2) + a(3)*Efac(3) = nval
      !
      ! Args:
      !   sigma(3,3) (input, complex(dp)): dielectric tensor ratios (sigmaPLUMEij/sigmaJETPLUMEij)
      !   a(3)       (input, complex(dp)): values related to density moments contribution by each field (set Ex = 1/0/0 Ey = 0/1/0 Ez=0/0/1 and compute JP's density moment to get a1/a2/a3 repsectively)
      !   nval       (input, complex(dp)): The known scalar 'n1s/n0' from plume.
      !   fac(3,3)   (output, real(dp)): The calculated fac matrix. These are the scalar term that the sigmaJETPLUMEij differs by the sigmaPLUMEij
      !   Efac(3)    (output, complex(dp)): The calculated Efac vector. These are the complex scalar tears
      !   ierr       (output, integer): Error code:
      !                                  0 = Success
      !                                  1 = Inconsistent system (no solution)
      !                                  2 = Ambiguous solution (scaling arbitrary, nval=0, S=0)
      !                                  3 = Undetermined fac (Efac=0 or all sigma=0)
      !                                  4 = Calculation failure (e.g., rank(M) < 2)
      !
      ! Uses double precision defined by 'dp'.
      ! ======================================================================
      implicit none

      ! Define double precision kind
      integer, parameter :: dp = kind(1.0d0)

      ! Arguments
      complex(dp), intent(in)  :: sigma(3, 3), a(3), nval
      complex(dp), intent(out)    :: fac(3, 3) !!TODO: this should eventually be real!
      complex(dp), intent(out) :: Efac(3)
      integer, intent(out)     :: ierr

      ! Local variables
      complex(dp) :: M(3, 3)         ! Constraint matrix for Efac
      complex(dp) :: detM           ! Determinant of M
      complex(dp) :: adjM(3, 3)      ! Adjugate of M
      complex(dp) :: facvec(3)      ! Non-normalized vector in null space of M
      complex(dp) :: col_norm(3)    ! Norms of columns of adjM
      complex(dp) :: norm_facvec    ! Norm of chosen facvec
      complex(dp) :: S              ! Dot product a . facvec
      complex(dp) :: c              ! Scaling factor
      real(dp) :: epsilon        ! Tolerance for zero checks
      integer :: i, j, max_col   ! Loop counters, index
      logical :: all_sigma_zero

      ! --- Initialization ---
      epsilon = 1.0e-12_dp
      fac = 0.0_dp
      Efac = 0.0_dp
      ierr = 0 ! Assume success initially

      ! --- Check if all sigma are zero ---
      all_sigma_zero = .true.
      do i = 1, 3
         do j = 1, 3
            if (abs(sigma(i, j)) > epsilon) then
               all_sigma_zero = .false.
               exit ! Exit inner loop
            end if
         end do
         if (.not. all_sigma_zero) exit ! Exit outer loop
      end do

      if (all_sigma_zero) then
         ! If all sigma are zero, then fac(i,j)*Efac(k) = 0.
         ! This means either Efac=0 or some fac=0.
         ! Check the normalization condition: a . Efac = nval
         ! If nval is non-zero, Efac cannot be zero.
         ! If nval is zero, Efac could be zero, or non-zero if a . Efac = 0.
         if (abs(nval) > epsilon) then
            ! Need non-zero Efac such that a . Efac = nval.
            ! M is the zero matrix, so any vector is a solution to M*v=0.
            ! Pick a direction, e.g., facvec = a / norm(a) if a is non-zero
            if (sqrt(abs(dot_product(a, a))) > epsilon) then !TODO: had to convert this logic to make work with complex- check it
               facvec = a/sqrt(dot_product(a, a))
               S = dot_product(a, facvec) ! Should be norm(a)
               c = nval/S
               Efac = c*facvec
               ierr = 3 ! fac is undetermined (0/non_zero)
               return
            else ! a is zero vector
               ierr = 1 ! Inconsistent: 0 = nval (non-zero)
               return
            end if
         else ! nval is zero
            ! a . Efac = 0. Efac could be zero, or any vector orthogonal to a.
            Efac = 0.0_dp ! Simplest solution
            ierr = 3 ! fac is undetermined
            return
         end if
      end if

      ! --- Construct the matrix M ---
      ! M = | sigma(1,2)  sigma(2,1)     0      |
      !     | sigma(1,3)     0      -sigma(3,1) |
      !     |    0       sigma(2,3)  sigma(3,2) |
      M(1, 1) = sigma(1, 2); M(1, 2) = sigma(2, 1); M(1, 3) = 0.0_dp
      M(2, 1) = sigma(1, 3); M(2, 2) = 0.0_dp; M(2, 3) = -sigma(3, 1)
      M(3, 1) = 0.0_dp; M(3, 2) = sigma(2, 3); M(3, 3) = sigma(3, 2)

      ! --- Calculate the determinant of M ---
      detM = M(1, 1)*(M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2)) - &
             M(1, 2)*(M(2, 1)*M(3, 3) - M(2, 3)*M(3, 1)) + &
             M(1, 3)*(M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1))

      ! --- Check if M is singular ---
      if (abs(detM) > epsilon) then
         ! If detM is non-zero, the only solution to M*v = 0 is v = 0.
         ! This implies Efac = 0.
         ! But we established not all sigma are zero, which requires non-zero Efac.
         ! Therefore, the system is inconsistent with the assumptions.
         ierr = 1 ! Inconsistent system
         return
      end if

      ! --- M is singular (or close to it), find null space vector ---
      ! Calculate Adjugate(M) = Transpose(Cofactor(M))
      ! Any non-zero column of Adjugate(M) is in the null space.
      adjM(1, 1) = M(2, 2)*M(3, 3) - M(2, 3)*M(3, 2) ! C11
      adjM(2, 1) = M(1, 3)*M(3, 2) - M(1, 2)*M(3, 3) ! C12 -> adj(2,1)
      adjM(3, 1) = M(1, 2)*M(2, 3) - M(1, 3)*M(2, 2) ! C13 -> adj(3,1)

      adjM(1, 2) = M(2, 3)*M(3, 1) - M(2, 1)*M(3, 3) ! C21 -> adj(1,2)
      adjM(2, 2) = M(1, 1)*M(3, 3) - M(1, 3)*M(3, 1) ! C22
      adjM(3, 2) = M(1, 3)*M(2, 1) - M(1, 1)*M(2, 3) ! C23 -> adj(3,2)

      adjM(1, 3) = M(2, 1)*M(3, 2) - M(2, 2)*M(3, 1) ! C31 -> adj(1,3)
      adjM(2, 3) = M(1, 2)*M(3, 1) - M(1, 1)*M(3, 2) ! C32 -> adj(2,3)
      adjM(3, 3) = M(1, 1)*M(2, 2) - M(1, 2)*M(2, 1) ! C33

      ! Find the column of adjM with the largest norm
      col_norm(1) = sqrt(adjM(1, 1)**2 + adjM(2, 1)**2 + adjM(3, 1)**2)
      col_norm(2) = sqrt(adjM(1, 2)**2 + adjM(2, 2)**2 + adjM(3, 2)**2)
      col_norm(3) = sqrt(adjM(1, 3)**2 + adjM(2, 3)**2 + adjM(3, 3)**2)

      max_col = 1
      if (abs(col_norm(2)) > abs(col_norm(max_col))) max_col = 2 !TODO: had to convert this logic to make work with complex- check it
      if (abs(col_norm(3)) > abs(col_norm(max_col))) max_col = 3 !TODO: had to convert this logic to make work with complex- check it

      facvec(:) = adjM(:, max_col)
      norm_facvec = col_norm(max_col)

      ! Check if a non-zero null vector was found
      if (abs(norm_facvec) < epsilon) then
         ! This happens if rank(M) < 2. The adjugate is the zero matrix.
         ! Since M is not the zero matrix (checked earlier), rank is 1.
         ! The null space is a plane, cannot determine unique direction this way.
         ierr = 4 ! Calculation failure / Rank < 2
         return
      end if

      ! Normalize the null space vector
      facvec = facvec/norm_facvec

      ! --- Use the normalization condition to find scaling factor c ---
      S = dot_product(a, facvec)

      if (abs(S) < epsilon) then
         ! a is orthogonal to the null space vector facvec
         if (abs(nval) < epsilon) then
            ! a . Efac = 0 is required, and nval is 0.
            ! The scaling factor c is arbitrary.
            ierr = 2 ! Ambiguous solution
            c = 1.0_dp ! Return solution for c=1
         else
            ! a . Efac = 0 is required, but nval is non-zero.
            ierr = 1 ! Inconsistent system
            return
         end if
      else
         ! Unique scaling factor c exists
         ierr = 0 ! Unique solution found (so far)
         c = nval/S
      end if

      ! --- Calculate final Efac ---
      Efac = c*facvec

      ! --- Calculate fac components ---
      do i = 1, 3 ! Calculate diagonal elements
         if (abs(Efac(i)) > epsilon) then
            fac(i, i) = sigma(i, i)/Efac(i)
         else
            ! If Efac(i) is zero, sigma(i,i) must also be zero for consistency
            if (abs(sigma(i, i)) > epsilon) then
               ierr = 1 ! Inconsistent
               return
            else
               fac(i, i) = 0.0_dp ! Or could be considered undetermined, but 0 is simplest
            end if
         end if
      end do

      ! Calculate independent off-diagonal elements (fxy, fxz, fyz)
      ! fxy = fac(1,2) = sigma(1,2) / Efac(2)
      if (abs(Efac(2)) > epsilon) then
         fac(1, 2) = sigma(1, 2)/Efac(2)
      else
         if (abs(sigma(1, 2)) > epsilon) then; ierr = 1; return; end if
         fac(1, 2) = 0.0_dp
      end if
      ! fxz = fac(1,3) = sigma(1,3) / Efac(3)
      if (abs(Efac(3)) > epsilon) then
         fac(1, 3) = sigma(1, 3)/Efac(3)
      else
         if (abs(sigma(1, 3)) > epsilon) then; ierr = 1; return; end if
         fac(1, 3) = 0.0_dp
      end if
      ! fyz = fac(2,3) = sigma(2,3) / Efac(3)
      if (abs(Efac(3)) > epsilon) then
         fac(2, 3) = sigma(2, 3)/Efac(3)
      else
         if (abs(sigma(2, 3)) > epsilon) then; ierr = 1; return; end if
         fac(2, 3) = 0.0_dp
      end if

      ! Set dependent off-diagonal elements using symmetry
      fac(2, 1) = -fac(1, 2) ! fyx = -fxy
      fac(3, 1) = fac(1, 3) ! fzx =  fxz
      fac(3, 2) = -fac(2, 3) ! fzy = -fyz

      ! Final consistency checks (already implicitly done by ensuring Efac is in null space
      ! and checking divisions by zero, but could add explicit checks if desired)
      ! e.g., check if abs(fac(2,1)*Efac(1) - sigma(2,1)) < epsilon * abs(sigma(2,1))

   end subroutine solve_norm_facs

end module functions
