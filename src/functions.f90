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
  
  integer, parameter :: stdout_unit=6
  !! Standard index for I/O.
  
  integer, save :: input_unit_no
  !! Index for reading in files.
  
  integer, save ::error_unit_no=stdout_unit
  !! Index for outputing error messages.

  character(50) :: runname
  !! String for input file parameters.

  public :: read_in_params, read_map_input, read_scan_input, read_guess_input
  public :: get_unused_unit, read_radial_input

contains

  !-=-=-=-=-
  !-=-=-=-=-
  subroutine read_in_params
    !!Read in system parameters.
    !!Input file is argument after executable:
    !!$ ./plume.e example.in
    use vars, only : betap,kperp,kpar,vtp,nspec,spec,susc,option,writeOut
    use vars, only : dataName,nscan,nroot_max,use_map,outputName, pi
    use vars, only : low_n, susc_low, new_low_n
    use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv
    use vars, only : vxmin,vxmax,vymin,vymax,vzmin,vzmax,elecdircontribution
    implicit none
    
    real :: sum_nq
    !!For Testing Quasineutrality.

    real :: sum_nqv
    !!For zero net current.

    integer :: is
    !!Index for scan and species loops.

    !Read in the basic parameter list.
    nameList /params/ &
         betap,kperp,kpar,vtp,nspec,nscan,option,nroot_max,&
         use_map,low_n, new_low_n, &
         writeOut,dataName,outputName
    
    nameList /fpc/ &
         vperpmin,vperpmax,vparmin,vparmax,delv,&
         vxmin,vxmax,vymin,vymax,vzmin,vzmax,elecdircontribution

    !default values
    elecdircontribution = 0.

    call get_unused_unit (input_unit_no)
    call get_runname(runname)
    runname=trim(runname)//".in"
    unit=input_unit_no
    open (unit=unit,file=runname,status='old',action='read')
    read (unit=unit,nml=params)

    !Allocate the species variable to have nspec indicies
    allocate (spec(1:nspec))

    !Allocate the susceptibility tensor to have nspec indicies 
    allocate(susc(1:nspec,3,3))
    if (low_n) &
         allocate(susc_low(1:nspec,3,3,0:1))

    !initialized quasineutrality and current conservation check.
    sum_nq = 0.; sum_nqv = 0.

    !Read in species parameters
    !This is a bit of FORTRAN black magic borrowed from AGK.
    !     which allows us to loop over iterative nml/groupnames.
    do is = 1,nspec 
       call get_indexed_namelist_unit (unit, "species", is)
       call spec_read(is)
       !quasineutrality check.
       sum_nq = sum_nq + spec(is)%D_s/spec(is)%Q_s
       !net current check.
       sum_nqv = sum_nqv + spec(is)%vv_s*spec(is)%D_s/spec(is)%Q_s
       close (unit)
    enddo

    !Check that ref mass ratio is 1 (in the python wrapper, this is used to help determine the output format to load in sweeps, which is why we check it)
    if(spec(1)%mu_s-1 > 0.0000001) then 
        write(*,'(a,es11.4)')&
         'ERROR: reference mass ratio is not 1:     spec(1)%mu_s =',spec(1)%mu_s
    end if

    !read in fpc params
    if(option == 6 .or. option == 7) then

      rewind(input_unit_no)
      read (unit=input_unit_no,nml=fpc)
    end if


    close (input_unit_no)

    !Notify if not Quasineutral,
    if (abs(sum_nq).ge.1.E-6) &
         write(*,'(a,es11.4)')&
         'ERROR: Plasma not quasineutral:     Sum [n_s Q_s] =',sum_nq
    !or not Current-free (not strictly required by linear theory but plasma is likely unstable and thus likely poorly described by linear theory)
    if (abs(sum_nqv).ge.1.E-6) &
         write(*,'(a,es11.4)')&
         'ERROR?/WARNING: Non zero Current:        Sum [n_s q_s V_s] =',sum_nqv
    !Take note User: Make sure the above sums are zero.
    !Note if:
    if (spec(1)%vv_s.ne.0.) &
         write(*,'(a,es11.4)')&
         'ERROR?/WARNING: Not in refernece species rest frame:    v_par drift ref =',spec(1)%vv_s

    !You always need pi. Always. And sometimes pie too.
    pi = 4.*atan(1.)
  end subroutine read_in_params

!-=-=-=-=-
!-=-=-=-=-
  subroutine spec_read(is)
    !!Subroutine for reading in species/component parameters.
    use vars, only : spec
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
       tauS,muS,alphS,Qs,Ds,vvS
  read (unit=unit,nml=species)

  !Assign dummy values to global component array.
  spec(is)%tau_s=tauS
  spec(is)%mu_s=muS
  spec(is)%alph_s=alphS
  spec(is)%Q_s=Qs
  spec(is)%D_s=Ds
  spec(is)%vv_s=vvS
  
end subroutine spec_read


!-=-=-=-=
!-=-=-=-=
subroutine read_map_input
  !!Read in parameters for bounds on mapping dispersion roots.
  !! Invokes [[map_read(subroutine)]].
  use vars, only: loggridw,loggridg,omi,omf,gami,gamf
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
  use vars, only: loggridw,loggridg,omi,omf,gami,gamf,positive_roots
  implicit none

  nameList /maps/ loggridw,loggridg,omi,omf,gami,gamf,positive_roots

  read (unit=unit,nml=maps)
  
end subroutine map_read

!-=-=-=-=
!-=-=-=-=
subroutine read_scan_input
  !!Read in limits for scans in plasma parameter space.
  !! Invokes [[scan_read(subroutine)]].
  use vars, only: nscan,scan
  implicit none
  integer :: is
  !!Scan index.

  !Allocate scan type object
  allocate(scan(1:nscan))
  
  call get_unused_unit (input_unit_no)
  call get_runname(runname)
  runname=trim(runname)//".in"
  unit=input_unit_no
  open (unit=unit,file=runname,status='old',action='read')  

  do is = 1,nscan
     call get_indexed_namelist_unit (unit, "scan_input", is)
     call scan_read(is)
     close (unit)
  enddo
  close (input_unit_no)
  
end subroutine read_scan_input

!-=-=-=-=-
!-=-=-=-=-
subroutine scan_read(is)
  !!Subroutine for reading in scan parameters.
  use vars, only : scan
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
       scan_type,swi,swf,swlog,scan_style,ns,nres,&
       heating,eigen,tensor
  read (unit=unit,nml=scan_input)

  !Assign dummy values to scan array.
  scan(is)%range_i =swi
  scan(is)%range_f =swf
  scan(is)%log_scan=swlog
  scan(is)%type_s  =scan_type
  scan(is)%style_s =scan_style
  scan(is)%n_scan  =ns
  scan(is)%n_res   =nres
  scan(is)%eigen_s =eigen
  scan(is)%heat_s  =heating
  scan(is)%tensor_s  =tensor
  
end subroutine scan_read

!-=-=-=-=
!-=-=-=-=
subroutine read_guess_input
  !!Read in initial values for solutions.
  !! Invokes [[om_read(subroutine)]].
  use vars, only: nroot_max,wroots
  implicit none
  
  integer :: is
  !!Solution index
  
  call get_unused_unit (input_unit_no)
  call get_runname(runname)
  runname=trim(runname)//".in"
  unit=input_unit_no
  open (unit=unit,file=runname,status='old',action='read')  

  do is = 1,nroot_max
     call get_indexed_namelist_unit (unit, "guess", is)
     call om_read(is)
     
     write(*,'(a,i0,a,2es12.4)')'Input for Root ',is,': ',wroots(1:2,is)
     close(unit)
  enddo
  close (input_unit_no)
  
end subroutine read_guess_input

!-=-=-=-=-
!-=-=-=-=-
subroutine om_read(is)
  !!Subroutine for reading in initial guesses for complex frequencies of solutions.
  use vars, only : wroots
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
       g_om,g_gam
  read (unit=unit,nml=guess)

  !Assign dummy values to frequency of solution.
  wroots(1,is)=g_om
  wroots(2,is)=g_gam
  
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
  runname=trim(runname)//".in"

  call get_unused_unit (input_unit_no)
  unit=input_unit_no
  open (unit=unit,file=runname,status='old',action='read')
  call radial_read
  close (unit)

  write(*,*)'Read in Radial Parameters'
  !Allocate the radial parameter scan for each species 
  allocate (rad_spec(1:nspec,0:nRad))
  !Allocate Radial Distance from Sun
  allocate (radius(0:nRad)); radius = 0.

  !Read in radial parameters for each species
  do is = 1, nspec
     write(readName,'(a,i0,a)')&
          trim(modelName),is,".rad"
     call get_unused_unit (input_unit_no)
     unit=input_unit_no
     open(unit = unit, file=trim(readName) , status='old', action='read')
     do ir = 0,nRad
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
  allocate (beta_rad(0:nRad)); beta_rad =0.
  allocate (vtp_rad(0:nRad)); vtp_rad =0.

  !Read in Global Parameters, beta_||p and vtp
  write(readName,'(a,i0,a)')&
       trim(modelName),0,".rad"
  call get_unused_unit (input_unit_no)
  unit=input_unit_no
  open(unit = unit, file=trim(readName) , status='old', action='read')
  do ir = 0,nRad
     read(5,*)&
          radius(ir),beta_rad(ir), vtp_rad(ir)          
  enddo
  close(unit)

  !Set parameters to values at intial radial position
  betap = beta_rad(0)
  vtp   = vtp_rad(0)

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
     call radial_read_0
  case(1)
     !fixed kperp, scan over kpar
     call radial_read_1
  case(2)
     !fixed kpar,  scan over kperp
     call radial_read_2
  case(3)
     !fixed theta, scan over k
     call radial_read_3
  case(4)
     !fixed k, scan over theta
     call radial_read_4
  case(5)
     !plane scan over (kperp, kpar)
     call radial_read_5
  case(6)
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
  read (unit=4,nml=radial_input)

  
end subroutine radial_read

!-=-=-=-=-
subroutine radial_read_0
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(k_{\perp}\) and \(k_{\parallel}\).
  use vars, only : kperp, kpar
  implicit none
  
  real :: kperp_1
  !! Dummy variable for \(k_{\perp}\) readin.
  
  real :: kpar_1
  !! Dummy variable for \(k_{\parallel}\) readin.

  nameList /k_range/ kperp_1,kpar_1     
  open (unit=4,file=runname,status='old',action='read')
  read (unit=4,nml=k_range)
  close (4)
  
  !Set initial wavevector.
  kperp = kperp_1
  kpar = kpar_1

end subroutine radial_read_0

!-=-=-=-=-
subroutine radial_read_1
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(k_{\perp}\) and varying \(k_{\parallel}\).
  use vars, only : kperp, kpar, rad_scan
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
  nameList /k_range/ kperp_1,kpar_1,kpar_2,nK,kres,&
       rad_log,rad_heat,rad_eigen
  open (unit=4,file=runname,status='old',action='read')
  read (unit=4,nml=k_range)
  close (4)
  !Set initial wavevector
  kperp = kperp_1
  kpar = kpar_1
  !Set \(k_{\parallel}\) range.
  allocate(rad_scan(1))
  rad_scan(1)%range_i = kpar_1
  rad_scan(1)%range_f = kpar_2
  rad_scan(1)%n_scan  = nk
  rad_scan(1)%n_res   = kres
  rad_scan(1)%heat_s   = rad_heat
  rad_scan(1)%eigen_s  = rad_eigen
  rad_scan(1)%log_scan = rad_log
  if (rad_log) then
     rad_scan(1)%diff =&
          (log10(kpar_2)-log10(kpar_1))/&
          real(nk*kres)
  else
     rad_scan(1)%diff =&
          ((kpar_2)-(kpar_1))/&
          real(nk*kres)
  endif

end subroutine radial_read_1

!-=-=-=-=-
subroutine radial_read_2
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(k_{\parallel}\) and varying \(k_{\perp}\).
  use vars, only : kperp,kpar,rad_scan
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
  nameList /k_range/ kperp_1,kperp_2,kpar_1,nK,kres,&
       rad_log,rad_heat,rad_eigen
  open (unit=4,file=runname,status='old',action='read')
  read (unit=4,nml=k_range)
  close (4)
  kperp = kperp_1
  kpar = kpar_1
  allocate(rad_scan(1))
  rad_scan(1)%range_i = kperp_1
  rad_scan(1)%range_f = kperp_2
  rad_scan(1)%n_scan  = nk
  rad_scan(1)%n_res   = kres
  rad_scan(1)%heat_s   = rad_heat
  rad_scan(1)%eigen_s  = rad_eigen
  rad_scan(1)%log_scan = rad_log
  if (rad_log) then
     rad_scan(1)%diff =&
          (log10(kperp_2)-log10(kperp_1))/&
          real(nk*kres)
  else
     rad_scan(1)%diff =&
          ((kperp_2)-(kperp_1))/&
          real(nk*kres)
  endif
  
end subroutine radial_read_2

!-=-=-=-=-
subroutine radial_read_3
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(\theta\) and varying \(|k|\).
  use vars, only : kperp, kpar, rad_scan, pi
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
  nameList /k_range/ k_1,k_2,theta_1,nK,kres,&
       rad_log,rad_heat,rad_eigen
  open (unit=4,file=runname,status='old',action='read')
  read (unit=4,nml=k_range)
  close (4)
  !Set initial wavevector.
  kperp = k_1*sin(theta_1*pi/180.)
  kpar = k_1*cos(theta_1*pi/180.)
  allocate(rad_scan(1))
  rad_scan(1)%range_i = k_1
  rad_scan(1)%range_f = k_2
  rad_scan(1)%n_scan  = nk
  rad_scan(1)%n_res   = kres
  rad_scan(1)%heat_s   = rad_heat
  rad_scan(1)%eigen_s  = rad_eigen
  rad_scan(1)%log_scan = rad_log
  if (rad_log) then
     rad_scan(1)%diff =&
          (log10(k_2)-log10(k_1))/&
          real(nk*kres)
  else
     rad_scan(1)%diff =&
          ((k_2)-(k_1))/&
          real(nk*kres)
  endif

end subroutine radial_read_3

!-=-=-=-=-
subroutine radial_read_4
  !!Subroutine for reading in radial scan parameters
  !!with fixed \(|k|\) and varying \(\theta\).
  use vars, only : kperp, kpar, rad_scan, pi
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
  nameList /k_range/ k_1,theta_1,theta_2,nK,kres,&
       rad_log,rad_heat,rad_eigen
     open (unit=4,file=runname,status='old',action='read')
     read (unit=4,nml=k_range)
     close (4)
     !set initial k
     kperp = k_1*sin(theta_1*pi/180.)
     kpar = k_1*cos(theta_1*pi/180.)
     allocate(rad_scan(1))
     rad_scan(1)%range_i = theta_1
     rad_scan(1)%range_f = theta_2
     rad_scan(1)%n_scan  = nk
     rad_scan(1)%n_res   = kres
     rad_scan(1)%heat_s   = rad_heat
     rad_scan(1)%eigen_s  = rad_eigen
     rad_scan(1)%log_scan = rad_log
     if (rad_log) then
        rad_scan(1)%diff =&
             (log10(theta_2)-log10(theta_1))/&
             real(nk*kres)
     else
        rad_scan(1)%diff =&
             ((theta_2)-(theta_1))/&
             real(nk*kres)
     endif
     
end subroutine radial_read_4

!-=-=-=-=-
subroutine radial_read_5
  !!Subroutine for reading in radial scan parameters
  !!for 2D scan over \(k_{\perp}\) and \(k_{\parallel}\).
  use vars, only : kperp, kpar, rad_scan
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
  
  nameList /k_range/ kperp_1,kperp_2,kpar_1,kpar_2,&
       nk,kres,nk2,kres2,rad_log_perp,rad_log_par,rad_heat,rad_eigen
  open (unit=4,file=runname,status='old',action='read')
  read (unit=4,nml=k_range)
  close (4)
  
  !Set initial wavevector.
  kperp = kperp_1
  kpar = kpar_1
  allocate(rad_scan(1:2))
  !Perpendicular scan parameters.
  rad_scan(1)%range_i = kperp_1
  rad_scan(1)%range_f = kperp_2
  rad_scan(1)%n_scan  = nk
  rad_scan(1)%n_res   = kres
  rad_scan(1)%heat_s   = rad_heat
  rad_scan(1)%eigen_s  = rad_eigen
  rad_scan(1)%log_scan = rad_log_perp
  if (rad_log_perp) then
     rad_scan(1)%diff =&
          (log10(kperp_2)-log10(kperp_1))/&
          real(nk*kres)
  else
     rad_scan(1)%diff =&
          ((kperp_2)-(kperp_1))/&
          real(nk*kres)
  endif
  
  !Parallel scan parameters.
  rad_scan(2)%range_i = kpar_1
  rad_scan(2)%range_f = kpar_2
  rad_scan(2)%n_scan  = nk2
  rad_scan(2)%n_res   = kres2
  rad_scan(2)%heat_s   = rad_heat
  rad_scan(2)%eigen_s  = rad_eigen
  rad_scan(2)%log_scan = rad_log_par
  if (rad_log_par) then
     rad_scan(2)%diff =&
          (log10(kpar_2)-log10(kpar_1))/&
          real(nk2*kres2)
  else
     rad_scan(2)%diff =&
          ((kpar_2)-(kpar_1))/&
          real(nk2*kres2)
  endif

end subroutine radial_read_5

!-=-=-=-=-
subroutine radial_read_6
  !!Subroutine for reading in radial scan parameters
  !!for 2D scan over \(|k|\) and \(\theta\).
  use vars, only : kperp, kpar, rad_scan, pi
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
  nameList /k_range/ k_1,k_2,theta_1,theta_2,&
       nk,kres,ntheta,thetares,rad_log_k,rad_log_theta,rad_heat,rad_eigen
  open (unit=4,file=runname,status='old',action='read')
  read (unit=4,nml=k_range)
  close (4)

  !Set initial wavevectors.
  kperp = k_1*sin(theta_1*pi/180.)
  kpar = k_1*cos(theta_1*pi/180.)
  allocate(rad_scan(1:2))
  !theta scan parameters.
  rad_scan(1)%range_i = theta_1
  rad_scan(1)%range_f = theta_2
  rad_scan(1)%n_scan  = ntheta
  rad_scan(1)%n_res   = thetares
  rad_scan(1)%heat_s   = rad_heat
  rad_scan(1)%eigen_s  = rad_eigen
  rad_scan(1)%log_scan = rad_log_theta
  if (rad_log_theta) then
     rad_scan(1)%diff =&
          (log10(theta_2)-log10(theta_1))/&
          real(ntheta*thetares)
  else
     rad_scan(1)%diff =&
          ((theta_2)-(theta_1))/&
          real(ntheta*thetares)
  endif

  !k scan
  rad_scan(2)%range_i = k_1
  rad_scan(2)%range_f = k_2
  rad_scan(2)%n_scan  = nk
  rad_scan(2)%n_res   = kres
  rad_scan(2)%heat_s   = rad_heat
  rad_scan(2)%eigen_s  = rad_eigen
  rad_scan(2)%log_scan = rad_log_k
  if (rad_log_k) then
     rad_scan(2)%diff =&
          (log10(k_2)-log10(k_1))/&
          real(nk*kres)
  else
     rad_scan(2)%diff =&
          ((k_2)-(k_1))/&
          real(nk*kres)
  endif

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
    call getarg(1,arg)

    !Check if this is the input file and trim .in extension to get runname
    l = len_trim (arg)
    if (l > 3 .and. arg(l-2:l) == ".in") then
       runname = arg(1:l-3)
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
  subroutine get_indexed_namelist_unit (unit, nml, index_in)
    !!Extract namelist.
    implicit none
    integer, intent (out) :: unit
    character (*), intent (in) :: nml
    integer, intent (in) :: index_in
    character(500) :: line
    integer :: iunit, iostat, in_file
    integer :: ind_slash
    logical :: exist

    call get_unused_unit (unit)
    ind_slash=index(runname,"/",.True.)
    if (ind_slash.EQ.0) then !No slash in name
        !Original behaviour
        open (unit=unit, file="."//trim(runname)//".scratch")
    else
        !General behaviour
        open (unit=unit, file=trim(runname(1:ind_slash))//"."//trim(runname(ind_slash+1:))//".scratch")
    endif

    write (line, *) index_in
    line = nml//"_"//trim(adjustl(line))
    in_file = input_unit_exist(trim(line), exist)

    if (exist) then
       iunit = input_unit(trim(line))
    else
       write(6,*) "get_indexed_namelist: following namelist not found ",trim(line)
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
  function input_unit_exist (nml,exist)
    !!Is a namelist already open?
    implicit none
    character(*), intent (in) :: nml
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
  function input_unit (nml)
    !!Determine I/O unit.
    implicit none
    character(*), intent (in) :: nml
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

  subroutine get_unused_unit (unit)
    !!Find a I/O unit that is not currently open.
    implicit none
    integer, intent (out) :: unit
    logical :: od
    unit = 50
    do
       inquire (unit=unit, opened=od)
       if (.not.od) return
       unit = unit + 1
    end do
  end subroutine get_unused_unit
!-=-=-=-=-=-


end module functions
