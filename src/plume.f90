!=============================================================================!
!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kgklein@arizona.edu                                                        !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!!*MAIN PROGRAM                                                             *!!
!=============================================================================!
!=============================================================================!

program plume
  use vars, only: option,nscan,use_map,nroot_max,wroots,writeOut,scan
  use functions, only: read_in_params,read_map_input,read_scan_input
  use functions, only: read_guess_input,read_radial_input
  use disprels, only: map_search,refine_guess,om_scan,om_double_scan,map_scan
  use disprels, only: test_disp, radial_scan
  implicit none
  
  !-=-=-=-=-=-=
  !Variable Declaration:
  !-=-=-=-=-=-=  
  
  integer :: is
  !!Index for looping through parametric variations with [[om_scan(subroutine)]].
  
  !-=-=-=-=-=-=
  !-=-=-=-=-=-=
  !Read in global and species parameters from input file *.in,
  !     specified in command line execution.
  !     Also allocates spec and susc.
  call read_in_params !functions.f90
       

  !The value of option (extracted from *.in file) determines 
  !    the nature of the dispersion relation calculations
  !    The case options are listed in README.cases
  select case(option)

  case(-1)
     !Calculate disp(om) at a single (omega, gamma)

        call read_guess_input !functions.f90
        
        call test_disp !disprels.f90

  case(0)   
     !Calculate Roots for input plasma parameters

     !Read in root mapping bounds
     call read_map_input !functions.f90

     !Calculate complex roots of the dispersion function
     !    At specified wavevector and plasma parameter.
     !    Saved as wroots(1:2,1:nroots)
     call map_search !disprels.f90
          
  case(1)   
     !Calculate Roots for input plasma parameters
     ! OR
     !Read in root values 
     ! THEN
     !Scan over plasma parameters, with range and type specified in *.in file

     if (use_map) then        
        !Read in root mapping bounds        
        call read_map_input !functions.f90
                
        !Calculate complex roots of the dispersion function
        !    At specified wavevector and plasma parameter.
        !    Saved as wroots(1:2,1:nroots)
        call map_search !disprels.f90
        
     else
        !Read in nroot_max complex frequency inputs from guess_N namelists in *.in file
        call read_guess_input !functions.f90
        
        !Take nroot_max inputs and refine guesses
        call refine_guess !disprels.f90
        
     endif
     
     !Read in parameter scan bounds
     call read_scan_input !functions.f90

     !Proceed serially through nscan parameter or wavevector scans.
     do is =1,nscan 
        call om_scan(is) !disprels.f90        
     enddo

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
        call read_map_input !functions.f90
               
        !Calculate complex roots of the dispersion function
        !    At specified wavevector and plasma parameter.
        !    Saved as wroots(1:2,1:nroots)
        call map_search !disprels.f90
        
     else
        !Read in nroot_max complex frequency inputs from guess_N namelists in *.in file
        call read_guess_input !functions.f90        

        !Take nroot_max inputs and refine guesses.
        call refine_guess !disprels.f90
     endif
     
     !Read in parameter scan bounds
     call read_scan_input !functions.f90
          
     !Scan over two parameters, read from scan_input_1 and scan_input_2 namelists.
     call om_double_scan !disprels.f90     
        
  case(3)
     !Replicating 
     !SAGA scan
     !from Gullveig (the precursor of this code)
     !A hardwired scan of 
     !  (k, theta) 
     !at a particular value of
     !  (betap, alph_p)
     !
     !DEPRECATED.
     write(*,'(a)')'Option 3 Deprecated. HALTING.'
     stop
          
  case(4)
     !Make multiple maps of complex frequency space

     !Only scan a single parameter at a time
     if (nscan.gt.1) then
        write(*,'(a)')'Option 4 only accepts one parameter scan. HALTING'
        write(*,'(a)')'Set nscan=1'
        write(*,'(a)')'HALTING'
        stop
     endif
     
     !Read in root mapping bounds
     call read_map_input !functions.f90
         
     !Read in parameter scan bounds
     call read_scan_input !functions.f90
     
     !Calculate complex frequency maps along the prescribed parameter scan.
     call map_scan !disprels.f90
          
  case(5)!Find roots for parameters along a prescribed path
     !Path is set by solar wind models, with values calculated and
     !output by helper function: (under development...)
     
     write(*,'(a)')'Option 5 is under development.'
     write(*,'(a)')'Use with caution...Here be Dragons...'

     !Read in radial scan parameters
     call read_radial_input !functions.f90
     
     !Determine inital roots
     if (use_map) then

        !Read in root mapping bounds
        call read_map_input !functions.f90
                
        !Calculate complex roots of the dispersion function
        !    At specified wavevector and plasma parameter.
        !    Saved as wroots(1:2,1:nroots)
        call map_search !disprels.f90
        
     else!Read in nroot_max (om,gam) inputs from *.in file
        call read_guess_input !functions.f90
        
        !Take nroot_max inputs and refine guesses
        call refine_guess !disprels.f90
        
     endif     

     !Scan roots over radius
     call radial_scan !disprels.f90
     
  end select

end program plume
