!=============================================================================!
!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kris.klein@gmail.com                                                       !!
!!Lunar and Planetary Laboratory, University of Arizona                      !!
!!                                                                           !!
!!*JET-PLUME                                                                 !!
!!Judging Energy Transfer in a Plasma in a Lin. Uni. Mag. Environment        !!
!!collin.crbrown@gmail.com                                                   !!
!!University of Iowa                                                         !!
!!*MAIN PROGRAM                                                             *!!
!=============================================================================!
!=============================================================================!

! NOTE: This code uses an F90 adaptation (by Greg Howes) of the Hot Plasma 
!       Dispersion Relation originally by Eliot Quataert.

! PLUME calculates the hot plasma dispersion relation for a plasma with 
!       an arbitrary number of ion and electron species with relative drifts
!       and bi-Maxwellian velocity distributions.
!       The calculation follows Stix Chapter 10 eqn 66-73.

! The Dispersion relation for omega/Omega_reference 
!     is dependent on four global parameters:
!
!       betap: Plasma Reference Beta:               8 pi n_ref T_ref /B^2
!       kperp: Perpendicular wavelength:         kperp rho_ref
!       kpar : Parallel wavelength:              kparallel rho_ref
!       vtp  : Parallel proton thermal velocity: sqrt(2 T_||ref/m_ref)
!
!     and six species dependent parameters:
!
!       tau_s : T_ref/T_s
!       mu_s  : m_ref/m_s
!       alph_s: T_perp/T_parallel|s
!       Q_s   : q_ref/q_s
!       D_s   : n_s/n_ref
!       vv_s  : v_s,drift/v_Aref
!
!The values for these parameters are extracted from *.in file, appended after
!    the executable program call.

program plume
  use vars, only: option,nscan,use_map,nroot_max,nroots,wroots,writeOut,scan
  use functions, only: read_in_params,read_map_input,read_scan_input
  use functions, only: read_guess_input,read_radial_input
  use disprels, only: map_search,refine_guess,om_scan,om_double_scan,map_scan
  use disprels, only: test_disp, radial_scan
  use fpc, only: compute_fpc_gyro, compute_fpc_cart, write_fs0

  implicit none

  integer :: is,i,j,k,iroot !Integers for code structure

  !Read in global and species parameters from file.in
  !     (specified in command line execution)
  !     Also allocates spec and susc.
  call read_in_params

  !The value of option (extracted from *.in file) determines 
  !    the nature of the dispersion relation calculations
  !    The case options are listed in README.cases
  select case(option)

  case(-1) !Calculate disp(om) at a single (omega, gamma)
        call read_guess_input

        call test_disp

  case(0) !Calculate Roots for input plasma parameters
     !Read in root mapping bounds
     call read_map_input
         
     !Calculate complex roots of the dispersion function (Saved as wroots(1:2,1:nroots))
     call map_search

  case(1) !Calculate Roots for input plasma parameters
          ! OR
          !Read in root values 
          ! THEN
          !Scan over plasma parameters, with range and type specified in *.in file
     
     if (use_map) then
        !Read in root mapping bounds
        call read_map_input
        
        !Calculate complex roots of the dispersion function (Saved as wroots(1:2,1:nroots))
        call map_search

     else !Read in nroot_max (om,gam) inputs from *.in file
        call read_guess_input

        !Take nroot_max inputs and refine guesses
        call refine_guess
     endif
     
     !Read in parameter scan bounds
     call read_scan_input

     do is =1,nscan 
        call om_scan(is)
     enddo

     write(*,*)'Done...'

  case(2)
     !Calculate Roots for input plasma parameters
     ! OR
     !Read in root values 
     ! THEN
     !Scan over two dimensional plasma parameter space
     !     with range and type specified in *.in file
     !Stored in a single file

     if (use_map) then
        !Read in root mapping bounds
        call read_map_input
        
        !Calculate complex roots of the dispersion function
        !    Saved as wroots(1:2,1:nroots)
        call map_search

     else!Read in nroot_max (om,gam) inputs from *.in file
        call read_guess_input

        !Take nroot_max inputs and refine guesses
        call refine_guess

     endif
     
     !Read in parameter scan bounds
     call read_scan_input

     !Scan over two parameters, given in scan_input_1 and scan_input_2
     call om_double_scan 

  case(3)
     !Replicating 
     !SAGA scan
     !from Gullveig (the precursor of this code)
     !A hardwired scan of 
     !  (k, theta) 
     !at a particular value of
     !  (betap, alph_p)

     !start at 
     !(beta,alph_p,kperp, kpar)=
     !(1.0, 1.0,   1.E-3,1.E-3)    

     !scan_1, scan_2 should be over theta, k_fixed, at desired resolution
     !scan_3, scan_4, scan_5  (beta_p, alph_p, (k1,k2))
     !-=-=-=-=-=
     !scan_1: theta (desired resolution) range_i:theta_final
              !Style: -1; Type: 1
     !scan_2: k_fixed angle (desired resolution) range_f:
              !Style: -1; Type: 2
     !scan_3: beta_p (1.0 -> desired beta_p)
              !Style:  0; Type: 2
     !scan_4: alph_p (1.0 -> desired alph_p)
              !Style:  1; Type: 2
     !scan_5: alph_p (1.E-3,1.E-3) -> 
              !(1.E-3 sin(theta_initial), 1.E-3 cos(theta_initial))
              !(range_i,range_f)
              !Style:  0; Type: 0   

     !Initial Roots
     !&guess_1
     !g_om=9.9973E-04   
     !g_gam=-2.2572E-10 
     !/

     !&guess_2
     !g_om=2.0304E-03   
     !g_gam=-5.4273E-05
     !/

     !&guess_3
     !g_om=0.0000E+00   
     !g_gam=-7.2110E-04
     !/

     !&guess_4
     !g_om=1.1830E-03
     !g_gam=-7.3333E-04
     !/  
     
     call read_guess_input
     
     !Take nroot_max inputs and refine guesses
     call refine_guess
     
     !Read in parameter scan bounds
     call read_scan_input

     do is =3,nscan 
        call om_scan(is)
     enddo

     !Scan over two parameters, given in scan_input_1 and scan_input_2
     call om_double_scan    
     
  case(4)!Make multiple maps of complex frequency space

     !Read in root mapping bounds
     call read_map_input
     
     !Read in parameter scan bounds
     !Only scan a single parameter at a time
     call read_scan_input

     !Calculate complex roots of the dispersion function
     
     call map_scan
 
  case(5)!Find roots for parameters along a prescribed path

     !Read in radial scan parameters
     call read_radial_input

     !Determine inital roots
     if (use_map) then
        !Read in root mapping bounds
        call read_map_input
        
        !Calculate complex roots of the dispersion function
        !    Saved as wroots(1:2,1:nroots)
        call map_search

     else !Read in nroot_max (om,gam) inputs from *.in file
        call read_guess_input


        !Take nroot_max inputs and refine guesses
        call refine_guess

     endif     

     !Scan roots over radius
     call radial_scan

   case(6) !calculate field particle correlation as a function of vperp vpar
     write(*,*)'Predicting FPC (gyro coords)...'

     if (use_map) then
        !Read in root mapping bounds
        call read_map_input

        !Calculate complex roots of the dispersion function (Saved as wroots(1:2,1:nroots))
        call map_search

     else!Read in nroot_max (om,gam) inputs from *.in file
        call read_guess_input

        !Take nroot_max inputs and refine guesses
        call refine_guess

     endif

     !Read in parameter scan bounds
     call read_scan_input

     do iroot=1,nroot_max
         call compute_fpc_gyro(iroot)
     enddo

     call write_fs0()

  case(7) !calculate field particle correlation as a function of vx vy vz (i.e. vperp1, vperp2, vpar)
     write(*,*)'Predicting FPC (cart coords)...'

     if (use_map) then
        !Read in root mapping bounds
        call read_map_input

        !Calculate complex roots of the dispersion function (Saved as wroots(1:2,1:nroots))
        call map_search

     else!Read in nroot_max (om,gam) inputs from *.in file
        call read_guess_input

        !Take nroot_max inputs and refine guesses
        call refine_guess

     endif

     !Read in parameter scan bounds
     call read_scan_input

     do iroot=1,nroot_max
         call compute_fpc_cart(iroot)
     enddo
  end select

end program plume
