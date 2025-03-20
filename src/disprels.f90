!!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kgklein@arizona.edu                                                        !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!!*DISPERSION FUNCTIONS                                                     *!!
!=============================================================================!
!=============================================================================!
module disprels
  !! Dispersion relation function and associated eigenfunction and minimization
  !! calculations.
  implicit none
  private

  integer, parameter :: nbrack=128
  !!Number of Bessel functions to sum

  real, parameter :: tol=1.0E-13
  !!Root Search Tolerance 

  real, parameter :: prec=1.E-7
  !!Root Finding precision  

  
  private :: find_minima,disp,rtsec,zet_in,zetout,bessel,get_out_name
  private :: calc_eigen

  public :: map_search, refine_guess, om_scan, om_double_scan, map_scan, test_disp
  public :: radial_scan

  contains
!-=-=-=-=-

!-=-=-=-=-
!------------------------------------------------------------------------------
!Greg Howes, 2006; Kristopher Klein, 2015
!------------------------------------------------------------------------------
    subroutine map_search
      !! Identifies minima in a prescribed region of complex frequency space.
      use vars, only : loggridw,loggridg,omi,omf,gami,gamf,nr,ni,nroot_max
      use vars, only : nroots,wroots,writeOut,dataName,numroots,scan,print_Name
      use vars, only : betap,kperp,kpar,vtp,nspec,spec,option,outputName,positive_roots
      implicit none
      
      real :: dr
      !!Spacing for real frequency array.

      real :: di
      !!Spacing for imaginary frequency array.

      real :: wr
      !!Real frequency value.

      real :: wi
      !!Imaginary frequency value.

      complex, dimension(:,:), pointer :: om
      !!Complex frequency array.

      complex, dimension(:,:), pointer :: dal
      !!(Complex) value of dispersion relation on frequency array.
      
      real, dimension(:,:), pointer :: val
      !!Amplitude of dispersion relation on frequency array.
      
      complex :: omega
      !!Complex Frequency
      
      logical, parameter :: outmap=.true.
      !!Output ASCII map file.
      
      character(100) :: mapName
      !!Output file name.
      
      integer, dimension(1:2,1:numroots) :: iroots
      !!Indices of roots
      
      logical, parameter :: refine=.true.
      !!T->refine solutions below frequency grid resolution.
      
      complex :: om1
      !!Lower bracket value of solution refinement.
      
      complex :: om2
      !!Upper bracket value of solution refinement.
      
      integer :: iflag
      !!Flag for Root search.
      
      real, dimension(1:6,1:nspec) :: params
      !!Species/component parameters.

      integer :: is
      !! Species/component index.

      integer :: ir
      !! Real Frequency Index.

      integer :: ii
      !! Imaginary Frequency Index.

      integer :: j
      !! Solution index.

      integer :: k
      !! Solution index; removing redundant solutions.

      integer :: ig
      !! Solution index; sorting solution.

      character (50) :: fmt
      !! Output format for solution.
      
      character (50) :: fmt_tnsr
      !! Output format for susceptibility tensor.

      real :: D_gap=1.d-5
      !! Threshold for identifying redundant solutions.

      logical :: unique_root
      !! Testing for redundant solutions.

      logical :: sort_solutions = .false.
      !! Sort solutions from least to most damped.

      logical :: large_solutions = .false.
      !! Eliminate LARGE solutions.
      
      real, dimension(:,:), allocatable :: temp_om
      !! Array for sorting solutions.
      
      real, dimension(:,:), allocatable :: temp_om2
      !! Array for sorting solutions.

      !Allocate array for map values
      allocate(val(0:nr,0:ni)); val(:,:)=0.
      allocate(dal(0:nr,0:ni)); dal(:,:)=cmplx(0.,0.)
      !Allocate array for complex frequencies
      allocate(om(0:nr,0:ni));  om(:,:) =cmplx(0.,0.)

      if (writeOut) then
         write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
         write(*,'(a)')'Root Search:'
         write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
         write(*,'(a)')      'Global Plasma Parameters:'
         write(*,'(a,g14.6)')'k_perp rho_ref   = ',kperp
         write(*,'(a,g14.6)')'k_par  rho_ref   = ',kpar
         write(*,'(a,g14.6)')'Beta_ref,par         = ',betap
         write(*,'(a,g14.6)')'v_t,ref,par/c          = ',vtp
         do is = 1, nspec
            write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
            write(*,'(a,i3)')      'Parameters for Species :',is
            write(*,'(a,g14.6)')'T_||p/T_||s =    ',spec(is)%tau_s
            write(*,'(a,g14.6)')'m_p/m_s =        ',spec(is)%mu_s
            write(*,'(a,g14.6)')'T_perp/T_par|s = ',spec(is)%alph_s
            write(*,'(a,g14.6)')'q_p/q_s =        ',spec(is)%Q_s
            write(*,'(a,g14.6)')'n_s/n_p =        ',spec(is)%D_s
            write(*,'(a,g14.6)')'v_drift s/c =    ',spec(is)%vv_s
         enddo
         write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
         write(*,'(a)')'Searching over:'
         write(*,'(a,es10.3,a,es10.3,a)')'om  \in [',omi,',',omf,']'
         write(*,'(a,es10.3,a,es10.3,a)')'gam \in [',gami,',',gamf,']'
         write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
      endif
    
    write(fmt,'(a,i0,a)')'(6es14.4,',6*nspec,'es14.4)'

    !Determine spacing in complex omega space (Normal or log)
    if (loggridw) then
       dr=(log10(omf)-log10(omi))/real(nr)
    else
       dr=(omf-omi)/real(nr)
    endif

    if (loggridg) then
       di=(log10(abs(gamf))-log10(abs(gami)))/real(ni)
    else
       di=(gamf-gami)/real(ni)
    endif


     !Scan over complex omega space and calculate dispersion relation
     do ir=0,nr
        do ii=0,ni
           if (loggridw) then
              wr=10.**(omi+dr*real(ir))
           else
              wr=omi+dr*real(ir)
           endif
           if (loggridg) then
              wi=-1*(10.**(abs(gami)+di*real(ii)))
           else
              wi=gami+di*real(ii)
           endif

           omega=cmplx(wr,wi)
           om(ir,ii)=omega
           
           dal(ir,ii)=disp(omega)
           val(ir,ii)=abs(dal(ir,ii))           
        enddo
     enddo

     write(*,*)'finding minima'
     
     !Find Local minima (roots) in map
     call find_minima(val,numroots,iroots,nroots)

     write(*,*)'minima found'
     
     if (writeOut) then
        write(*,'(i3,a)')nroots,'  possible local minima found:'
        do j=1,nroots
          write(*,'(a,i4,a,i4)')'ir = ',iroots(1,j),'    ii = ',iroots(2,j)
          write(*,'(2es14.4)')&
               real(om(iroots(1,j),iroots(2,j))),aimag(om(iroots(1,j),iroots(2,j)))
        enddo
     endif
     !Calculate complex frequency value at root and output
     wroots(:,:)=0.
     do j=1,nroots
        wroots(1,j)=real(om(iroots(1,j),iroots(2,j)))
        wroots(2,j)=aimag(om(iroots(1,j),iroots(2,j)))
     enddo

     !Refine roots
     if (refine) then
        if (writeOut) write(*,'(a)')'Refining roots:'
        k=0
        do j=1,nroots
           if ((wroots(1,j) .ne. 0.) .or. (wroots(2,j) .ne. 0.)) then

              omega=cmplx(wroots(1,j),wroots(2,j))
              
              om1=omega*(1.-prec)
              om2=omega*(1.+prec)
              
              omega=rtsec(disp,om1,om2,tol,iflag)
              
              !check to see if the found root is already
              !an identified solution
              if (j.eq.1) then
                 k=k+1
                 wroots(1,k)=real(omega)
                 wroots(2,k)=aimag(omega)
              else
                 unique_root=.true.
                 do ig=1,j-1 !loop through previous solutions
                    if ( (abs(omega-cmplx(wroots(1,ig),wroots(2,ig)))).lt.D_gap) then
                       !repeated mode identified
                       unique_root=.false.
                       exit
                    endif
                 enddo
                 if (unique_root) then
                    k=k+1
                    wroots(1,k)=real(omega)
                    wroots(2,k)=aimag(omega)
                 endif
             endif !
          endif !end omega !=0 loop
       enddo !end nroots loop
    endif !end refine loop

    !Eliminate all roots with a (LARGE) positive damping rate
    if (large_solutions) then
       k=0
       do j=1,nroots
          !LARGE positive Damping rate
          if (wroots(2,j) .le. 1.) then
             k=k+1
             wroots(1:2,k)=wroots(1:2,j)
          endif
       enddo
       if ((k .lt. nroots).and.writeOut)&
            write(*,'(i3,a)')nroots-k,' roots with gamma/Omega_ref>1 eliminated'
       nroots=k
    endif

    !OUTPUT File with values 
    if (outmap) then       
       
       if (option==4) then
          do ii = 1,nspec
             params(1,ii) = spec(ii)%tau_s
             params(2,ii) = spec(ii)%mu_s
             params(3,ii) = spec(ii)%alph_s
             params(4,ii) = spec(ii)%q_s
             params(5,ii) = spec(ii)%D_s
             params(6,ii) = spec(ii)%vv_s
          enddo
          !If constructing mulitple (om,gamma) space maps
          !   for a variety of parameter values
          write(mapName,'(7a)')&
               'data/',trim(dataName),&
               '/dispersion_',trim(outputName),'_',trim(print_Name),&
               '.map'
       else
          !If only constructing a single .map file
          !   for a given Run
          write(mapName,'(6a)')&
               'data/',trim(dataName),'/',&
               'dispersion_',trim(outputName),'.map'
       endif
       
       open(unit=21,file=trim(mapName),status='replace')
       do ir=0,nr
          do ii=0,ni
             !kgk 250108: change in output format: supressing indices,
             !and outputting val and dal directly
             write(21,'(3es14.6,2es14.4)')&
                  om(ir,ii),log10(val(ir,ii)),&                  
                  dal(ir,ii)
          enddo
          write(21,*)
       enddo
       close(21)
        
        if (option==4) then
           !If constructing mulitple (om,gamma) space maps
           !   for a variety of parameter values
           write(mapName,'(7a)')&
                'data/',trim(dataName),&
                '/dispersion_',trim(outputName),'_',trim(print_Name),&
                '.roots'
        else
           !If only constructing a single .root file
           !   for a given Run
           write(mapName,'(6a)')&
                'data/',trim(dataName),'/',&
                'dispersion_',trim(outputName),'.roots'
        endif

        open(unit=21,file=trim(mapName),status='replace')
        do j=1,nroots
           write(21,fmt)&
                kperp,kpar,betap,vtp,wroots(1:2,j),params(1:6,1:nspec)
        enddo
        close(21)        
     endif

    deallocate(val)

    if (sort_solutions) then
       !Sorting from least to most damped mode
       allocate(temp_om(1:2,1:nroot_max));temp_om = -10.
       allocate(temp_om2(1:2,1:nroot_max));temp_om2 = -10.
       
       !Loop through located roots
       do j = 1, nroots
          do is = 1,nroot_max
             if (wroots(2,j).ge.temp_om(2,is)) then
                do ig = 1,is-1
                   temp_om2(1:2,ig)=temp_om(1:2,ig)
                enddo
                temp_om2(1:2,is) = wroots(1:2,j)
                do ig = is+1,nroot_max-1
                   temp_om2(1:2,ig)=temp_om(1:2,ig-1)
                enddo
                temp_om=temp_om2
                exit !exits 'is' loop
             endif
          enddo
       enddo
       
    endif

    !gamma less than 1.
       k=0
       do j=1,nroots
          !if(wroots(1,j).gt.(0.0))then
          if(wroots(2,j).lt.(1.0e0))then!
             k=k+1
             wroots(1:2,k)=wroots(1:2,j)
          endif
       enddo
       nroots=k

       !omega greater than Zero
    if (positive_roots) then
       write(*,'(a)')'Selected Roots with omega > -1.E-10'
       !Only choose roots with om .ge. 0 (or small)!
       k=0
       do j=1,nroots
          !if(wroots(1,j).gt.(0.0))then
          if(wroots(1,j).gt.(-1.0e-10))then!
             k=k+1
             wroots(1:2,k)=wroots(1:2,j)
          endif
       enddo
       nroots=k
    endif
   

    if (writeOut) then
       !WRITE out roots 
       write(*,'(a)')'Dispersion Solutions '
       do is = 1,nroot_max
          write(*,'(i3,2es14.4)')is,wroots(1:2,is)
       enddo
    endif


  end subroutine map_search

  !-=-=-=-=-
  !-=-=-=-=-
  subroutine test_disp
    !! Testing routine for single evaluation of dispersion relation.   
    use vars, only : wroots, nroot_max, writeOut
    use vars, only : kperp, kpar, vtp
    implicit none
    
    complex :: omega
    !! Complex Frequency input.
    
    complex :: D
    !! Output of dispersion relation function.
    
    omega= cmplx(wroots(1,1),wroots(2,1))
    
    D = disp(omega)
    
    write(*,'(4es14.4)') omega,D*vtp**6.
    write(*,'(2es14.4)')kperp,kpar

  end subroutine test_disp

!-=-=-=-=-

!-=-=-=-=-
subroutine refine_guess
  !!Refine solutions for roots of dispersion relation.
  use vars, only : wroots, nroot_max, writeOut
  implicit none
  
  complex :: omega
  !!Complex Frequency.

  complex :: om1
  !!Lower bracket value of solution refinement.

  complex :: om2
  !!Upper bracket value of solution refinement.
  
  integer :: iflag
  !!Flag for Root search
  
  integer :: j
  !!Solution index.

  do j=1,nroot_max
     if (wroots(1,j) .ne. 0. .or. wroots(2,j) .ne. 0.) then
        
        omega=cmplx(wroots(1,j),wroots(2,j))
        
        om1=omega*(1.-prec)
        om2=omega*(1.+prec)
        
        omega=rtsec(disp,om1,om2,tol,iflag)
        
        wroots(1,j)=real(omega)
        wroots(2,j)=aimag(omega)
        
     endif
  enddo

    if (writeOut) then
       !WRITE out roots 
       write(*,'(a)')'Dispersion Solutions: '
       do j = 1,nroot_max
          write(*,'(i3,2es14.4)')j,wroots(1:2,j)
       enddo
    endif

end subroutine refine_guess

!-=-=-=-=-
!-=-=-=-=-=
subroutine map_scan
  !! Calculate complex frequency maps for scanned parameters.
  !! Intervatively calls [[map_search(subroutine)]]
  !! over a fixed range of parameters.
  use vars, only : scan, spec, print_Name, sw, sw2, kperp, kpar
  use vars, only : writeOut, pi
  use vars, only : gami, gamf, omi, omf ,betap
  !Local

  integer :: jj
  !! Parameter step index.
  
  character(150) :: outName
  !! String for I/O identification.
    
  real :: diff
  !! Spacing for scan of primary parameter.

  real :: diff2
  !! Spacing for scan of secondary parameter.
  
  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)

  real :: theta_q
  !! Scanned values of \( \theta \).

  real :: ki
  !! Initial value of \( |k| \).

  real :: kperpi
  !! Initial value of \( k_{\perp} \).

  real :: kpari
  !! Initial value of \( k_{\parallel} \).

  !Assign output name for scan is
  call set_map_pointers(outName,diff,diff2)
  !pointer sw, (and sw2 if needed) is assigned in GET_OUT_NAME 

    if ((scan(1)%style_s)==-1) then
     !Scans with multiple components
     if ((scan(1)%type_s)==0) then
        !k_0 -> k_1; 
        !k_0= current k
        !kperp_1 = range_i, kpar_1 = range_f
        kpari=kpar
        kperpi=kperp
     elseif ((scan(1)%type_s)==1) then
        !theta_0 -> theta_1
        !atan(kperp/kpar) -> theta_0
        !swf ->  theta_1 in degrees
        theta=atan(kperp/kpar)!initial theta value
        ki=kpar/cos(theta)
     elseif ((scan(1)%type_s)==2) then
        !k_0-> k_1 along fixed (current) angle
        !k_0 -> k_1; 
        !k_0= current k
        !|k_1|=swf
        theta = atan(kperp/kpar)
        ki=kpar/cos(theta)
     endif
  endif
  
  !Scan over chosen parameter

  !Output %n_scan steps, with %n_res steps inbetween each output
  !%n_res should be set to 1 for these type of runs, as there
  !is no information passed between each map calculation.
  do jj = 0, (scan(1)%n_scan*scan(1)%n_res) 
     !Advance scanned parameter values
     if ((scan(1)%style_s)==-1)then
        if ((scan(1)%type_s)==0) then
           if (scan(1)%log_scan) then
              !k0->k1
              sw=10.**(log10(kperpi)+diff*real(jj))    
              sw2=10.**(log10(kpari)+diff2*real(jj))    
           else
              sw=(kperpi)+diff*real(jj)    
              sw2=(kpari)+diff2*real(jj)    
           endif
        elseif ((scan(1)%type_s)==1) then
           !theta_0->theta_1
           if (scan(1)%log_scan) then
              theta_q=10.**(log10((theta))+diff*real(jj))    
           else
              theta_q=(((theta))+diff*real(jj))    
           endif
           sw=(ki*sin(theta_q))!kperp
           sw2=(ki*cos(theta_q))!kpar
        elseif ((scan(1)%type_s)==2) then
           !k along contant theta
           if (scan(1)%log_scan) then
              sw=10.**(log10(ki*sin(theta))+diff*real(jj))    
              sw2=10.**(log10(ki*cos(theta))+diff2*real(jj))    
           else
              sw=(ki*sin(theta))+diff*real(jj)    
              sw2=(ki*cos(theta))+diff2*real(jj)    
           endif
        endif
     else
        if (scan(1)%log_scan) then
           sw=10.**(log10(scan(1)%range_i)+diff*real(jj))    
        else
           sw=(scan(1)%range_i)+diff*real(jj)    
        endif
     endif

     if ((scan(1)%style_s)==-1) then
        !Scans with multiple components
        if ((scan(1)%type_s)==0) then
           !k_0 -> k_1; 
           !k_0= current k
           !kperp_1 = range_i, kpar_1 = range_f
           kpari=kpar
           kperpi=kperp
           write(print_Name,'(a,i0,a,i0)')&
                trim(outName),int(10000*sw),&
                '_',int(10000*sw2)
        elseif ((scan(1)%type_s)==1) then
           !theta_0 -> theta_1
           !atan(kperp/kpar) -> theta_0
           !swf ->  theta_1 in degrees
           write(print_Name,'(a,i0)')&
                trim(outName),int(10000*atan(sw/sw2)*180./pi)
        elseif ((scan(1)%type_s)==2) then
           !k_0-> k_1 along fixed (current) angle
           !k_0 -> k_1; 
           !k_0= current k
           !|k_1|=swf
           write(print_Name,'(a,i0)')&
                trim(outName),int(10000*sw/sin(theta))
        endif
     else
        write(print_Name,'(a,i0)')&
             trim(outName),int(10000*sw)
        write(*,*)trim(print_Name)
     endif

     !Edit map span. Hard code in the desired range.
     !Will add an improved interface for the future.
     
     write(*,'(a)')'Modifying Span of Complex Frequency Space:'
     omi=-3.*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
     omf=3.*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
     gami=-3.*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
     gamf=0.5*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
     
     call map_search

  enddo  
  !End Parameter Scan
  
end subroutine map_scan

!-=-=-=-=-=
!-=-=-=-=-=
subroutine om_scan(is)
  !! Calculate solutions as a function of the variation of
  !! a single parameter.
  use vars, only : wroots,scan,nroot_max,nspec,sw,sw2, low_n, pi
  use vars, only : new_low_n   !>>>GGH: 1/18/23
  use vars, only : kperp,kpar,betap,vtp,spec,writeOut,susc
  use functions, only : get_unused_unit
  
  implicit none
  !Passed
  integer :: is
  !!Scan index.

  !Local
  integer :: ii
  !! Solution index.
  
  integer :: jj
  !! Parameter step index.

  character(150) :: outName
  !! String for identifying output files.

  character(150) :: writeName
  !! Full output name.

  character(150) :: tensorName
  !! Output name for susceptibility tensor.
  
  character(100)  :: fmt
  !! Output format for dispersion relation.
  
  character(100)  :: fmt_tnsr
  !! Output format for susceptibility tensor.
  
  integer ,dimension(1:nroot_max) :: out_unit
  !! Main output unit.

  integer ,dimension(1:nroot_max) :: out_unit_2
  !! Suplementary  output unit.
  
  real :: diff
  !! Spacing for parameter scan.

  real :: diff2
  !! Supplemental Spacing for two-component parameter scan.
  
  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)
  
  real :: theta_q              
  !! Scanned values of \( \theta \).

  real :: ki
  !! Initial value of \( |k| \).

  real :: kf
  !! Final value of \( |k| \).

  real :: kperpi
  !! Initial value of \( k_{\perp} \).

  real :: kpari
  !! Initial value of \( k_{\parallel} \).
  
  integer :: out_type
  !! logic for outputing supplementary eigenfunction and heating
  !! calculations.
  !! 0: Frequency, Eigenfuction, Heating.
  !! 1: Frequency, Eigenfuction.
  !! 2: Frequency, Heating.
  !! 3: Frequency.

  complex,dimension(:),allocatable :: omlast
  !!Arrays with complex frequency for each solution.

  complex :: om1
  !!Lower bracket value of solution refinement.
  
  complex :: om2
  !!Upper bracket value of solution refinement.
  
  complex :: omold
  !!Previous frequency value.
  
  complex :: omega
  !!Complex Frequency under evaluation.
  
  real    :: val
  !!Dispersion Relation Value
  
  integer :: iflag
  !!Flag for Root search

  !Eigenfunctions
  complex, dimension(1:3)       :: ef
  !! Electric Field eigenfluctuation.
  
  complex, dimension(1:3)       :: bf
  !! Magnetic Field eigenfluctuation.

  complex, dimension(1:nspec)     :: ns
  !! Density eigenfluctuation.
  
  complex, dimension(1:3,1:nspec) :: Us
  !! Velocity eigenfluctuation.
  
  !Heating
  real, dimension(1:nspec) :: Ps
  !! Power into/out of species/component in one wave period.
  
  real, dimension(1:4,1:nspec) :: Ps_split
  !! Power into/out of species/components, broken into different contributions.
  !! Deprecated.
  
  !>>>GGH: 1/18/23
  real, dimension(1:6,1:nspec) :: Ps_split_new
  !!Power into/out of species/componets.
  !!Corrected LD/TTD calculation (GGH).

  real :: Ew
  !!Wave energy.

  real, dimension(1:6,1:nspec) :: params
  !!Species parameters:

  !Assign output name for scan is
  call get_out_name(outName,tensorName,fmt,fmt_tnsr,out_type,is,diff,diff2)
  !pointer sw, (and sw2 if needed) is assigned in GET_OUT_NAME 

  if (scan(is)%tensor_s) &
       write(*,'(2a)')' =>',trim(tensorName)

  if ((scan(is)%style_s)==-1) then
     !Scans with multiple components
     if ((scan(is)%type_s)==0) then
        !k_0 -> k_1; 
        !k_0= current k
        !kperp_1 = range_i, kpar_1 = range_f
        kpari=kpar
        kperpi=kperp
     elseif ((scan(is)%type_s)==1) then
        !theta_0 -> theta_1
        !atan(kperp/kpar) -> theta_0
        !swf ->  theta_1 in degrees
        theta=atan(kperp/kpar)!initial theta value
        ki=kpar/cos(theta)
     elseif ((scan(is)%type_s)==2) then
        !k_0-> k_1 along fixed (current) angle
        !k_0 -> k_1; 
        !k_0= current k
        !|k_1|=swf
        theta = atan(kperp/kpar)
        ki=kpar/cos(theta)
     endif
  endif


  !Open the (nroot_max) output files
  do ii = 1,nroot_max
     call get_unused_unit(out_unit(ii))
     write(writeName,'(2a,i0)')trim(outName),'.mode',ii
     open(unit=out_unit(ii),file=trim(writeName),status='replace')
     if (scan(is)%tensor_s) then
        call get_unused_unit(out_unit_2(ii))
        write(writeName,'(2a,i0)')trim(tensorName),'.mode',ii
        open(unit=out_unit_2(ii),file=trim(writeName),status='replace')
     endif
  enddo
  
  !Allocate variable for last solution and copy in initial values
  allocate(omlast(nroot_max))
  do ii=1,nroot_max
     omlast(ii)=cmplx(wroots(1,ii),wroots(2,ii))
     if (writeOut) &
          write(*,'(a,i3,a,2es14.4)')&
          'Root ',ii,' In: ',wroots(1,ii),wroots(2,ii)
  enddo

  !Scan over chosen parameter

  !Output %n_scan steps, with %n_res steps inbetween each output
  do jj = 0, (scan(is)%n_scan*scan(is)%n_res) 
     !Advance scanned parameter values
     if ((scan(is)%style_s)==-1)then
        if ((scan(is)%type_s)==0) then
           if (scan(is)%log_scan) then
              !k0->k1
              sw=10.**(log10(kperpi)+diff*real(jj))    
              sw2=10.**(log10(kpari)+diff2*real(jj))    
           else
              sw=(kperpi)+diff*real(jj)    
              sw2=(kpari)+diff2*real(jj)    
           endif
        elseif ((scan(is)%type_s)==1) then
           !theta_0->theta_1
           if (scan(is)%log_scan) then
              theta_q=10.**(log10((theta))+diff*real(jj))    
           else
              theta_q=(((theta))+diff*real(jj))    
           endif
           sw=(ki*sin(theta_q))!kperp
           sw2=(ki*cos(theta_q))!kpar
        elseif ((scan(is)%type_s)==2) then
           !k along contant theta
           if (scan(is)%log_scan) then
              sw=10.**(log10(ki*sin(theta))+diff*real(jj))    
              sw2=10.**(log10(ki*cos(theta))+diff2*real(jj))    
           else
              sw=(ki*sin(theta))+diff*real(jj)    
              sw2=(ki*cos(theta))+diff2*real(jj)    
           endif
        endif
     else
        if (scan(is)%log_scan) then
           sw=10.**(log10(scan(is)%range_i)+diff*real(jj))    
        else
           sw=(scan(is)%range_i)+diff*real(jj)    
        endif
     endif

     !construct params array for formated output
     if (mod(jj,scan(is)%n_res)==0) then
        do ii = 1,nspec
           params(1,ii) = spec(ii)%tau_s
           params(2,ii) = spec(ii)%mu_s
           params(3,ii) = spec(ii)%alph_s
           params(4,ii) = spec(ii)%q_s
           params(5,ii) = spec(ii)%D_s
           params(6,ii) = spec(ii)%vv_s
        enddo
     endif

     !Root Scan....
     do ii=1,nroot_max
        !Bracket values for root search
        omold=omlast(ii)
        om1=omold*(1.-prec)
        om2=omold*(1.+prec)

        iflag=0
        !Find New Omega Value
        omega=rtsec(disp,om1,om2,tol,iflag)
        
        !Save Root
        omlast(ii)=omega        
        
        !Output value every %n_res steps
        if (mod(jj,scan(is)%n_res)==0) then
           !Calculate eigenfns, heating only every %n_res steps
           
           if ((scan(is)%eigen_s).or.((scan(is)%heat_s))) then
              val=abs(disp(omega))
              !>>>GGH: 1/18/23
              call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,&
                   Ps_split_new,scan(is)%eigen_s,scan(is)%heat_s)
              !<<<GGH: 1/18/23
              omega=omlast(ii)
              if (abs(real(omega)).lt.1.E-15) then
                 Ps=0.1
                 Ps_split=-0.1
                 !>>>GGH: 1/18/23
                 Ps_split_new=-0.1
                 !<<<GGH: 1/18/23
              endif
           endif

           if (abs(real(omega)).lt.1.E-15) omega=cmplx(0.,aimag(omega))

           !1, 2, 3, 4

           !5, 6

           !7,8 9,10 11,12 : 13,14 15,16 17,18 : 
           !19,20 21,22 23,24 : 25,26 27,28 29,30 :
           !31,32 33,34
           !35, 36

           select case(out_type)
           case(0)!Om, Eigen, Heating
              if (low_n) then
                               !>>>GGH: 1/18/23
              if (new_low_n) then  !Greg's New low_n calculation to separate LD/TTD
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                      Ps(1:nspec),Ps_split_new(1:6,1:nspec),&
                      params(1:6,1:nspec)
              else  !Old version of low_n using Stix/Quataert version of LD/TTD 
              !<<<GGH: 1/18/23
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                      Ps(1:nspec),Ps_split(1:4,1:nspec),&
                      params(1:6,1:nspec)
              !>>>GGH: 1/18/23
              endif
              !<<<GGH: 1/18/23
              else
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                      Ps(1:nspec),&
                      params(1:6,1:nspec),Ew
              endif
           case(1) !Om, Eigen
              write(out_unit(ii),fmt)&
                   kperp,kpar,betap,vtp,&
                   omega,&            
                   bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                   params(1:6,1:nspec)
           case(2) !Om, Heating
              if (low_n) then
                 !>>>GGH: 1/18/23
              if (new_low_n) then  !Greg's New low_n calculation to separate LD/TTD
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      Ps(1:nspec),Ps_split_new(1:6,1:nspec),&
                      params(1:6,1:nspec)
              else  !Old version of low_n using Stix/Quataert version of LD/TTD 
              !<<<GGH: 1/18/23
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      Ps(1:nspec),Ps_split(1:4,1:nspec),&
                      params(1:6,1:nspec)
              !>>>GGH: 1/18/23
              endif
              !<<<GGH: 1/18/23
              else
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      Ps(1:nspec),&
                      params(1:6,1:nspec),Ew
              endif
           case(3) !Om
              write(out_unit(ii),fmt)&
                   kperp,kpar,betap,vtp,&
                   omega,&            
                   params(1:6,1:nspec)                        
           case default
              write(*,'(a)')&
                   'Out_type outside of allowable parameter range'
              stop
           end select
           if (scan(is)%tensor_s) then
              val=abs(disp(omega))
              write(out_unit_2(ii),fmt_tnsr)&
                   kperp,kpar,betap,vtp,&
                   susc(1:nspec,1,1),susc(1:nspec,1,2),susc(1:nspec,1,3),&
                   susc(1:nspec,2,2),susc(1:nspec,2,3),&
                   susc(1:nspec,3,3),&
                   params(1:6,1:nspec)
              !(1,1) xx; (1,2) xy; (1,3) xz;
              !(2,1) yx; (2,2) yy; (2,3) yz;
              !(3,1) zx; (3,2) zy; (3,3) zz;

           endif
        endif

        !Save Roots for further scans:        
        if (jj==scan(is)%n_scan*scan(is)%n_res) then
           wroots(1,ii)=real(omega)
           wroots(2,ii)=aimag(omega)
           if (writeOut) &
                write(*,'(a,i3,a,2es14.4)')&
                'Root ',ii,' Out: ',wroots(1,ii),wroots(2,ii)
        endif

     enddo
  enddo
  !End Parameter Scan

  !Close output files
  do ii=1,nroot_max
     close(out_unit(ii))
     if (scan(is)%tensor_s) &
          close(out_unit_2(ii))
  enddo

end subroutine om_scan
!-=-=-=-=-

!-=-=-=-=-
subroutine om_double_scan
  !! Calculate solutions as a function of the variation of
  !! two parameters, creating a surface in parameter space.
  use vars, only: scan, spec, betap, kpar, kperp, vtp, sw, sw2, sw3, sw4, pi
  use vars, only: nroot_max, outputName, wroots, nspec, writeOut, susc, low_n
  use vars, only : new_low_n   !>>>GGH: 1/18/23
  use functions, only : get_unused_unit
  implicit none
  
  !Local
  real, dimension (1:2,1:2) :: diff
  !!Parameter steps for both components.

  integer :: is
  !! Scan index.
  
  integer :: ii
  !! Solution index.
  
  integer :: jj
  !! First Parameter step index.

  integer :: kk
  !! Second Parameter step index.

  character(150) :: outName
  !! String for identifying output files.

  character(150) :: writeName
  !! Full output name.

  character(150) :: tensorName
  !! Output name for susceptibility tensor.
  
  character(100)  :: fmt
  !! Output format for dispersion relation.
  
  character(100)  :: fmt_tnsr
  !! Output format for susceptibility tensor.

  integer ,dimension(1:nroot_max) :: out_unit
  !! Main output unit.

  integer ,dimension(1:nroot_max) :: out_unit_2
  !! Suplementary  output unit.

  integer :: out_type
  !! logic for outputing supplementary eigenfunction and heating
  !! calculations.
  !! 0: Frequency, Eigenfuction, Heating.
  !! 1: Frequency, Eigenfuction.
  !! 2: Frequency, Heating.
  !! 3: Frequency.

  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)
  
  real :: theta_q              
  !! Scanned values of \( \theta \).

  real :: ki
  !! Initial value of \( |k| \).

  real :: kf
  !! Final value of \( |k| \).

  real :: kperpi
  !! Initial value of \( k_{\perp} \).

  real :: kpari
  !! Initial value of \( k_{\parallel} \).
  
  !Frequency

  complex,dimension(:),allocatable :: omlast
  !!Arrays with complex frequency for each solution.

  complex :: om1
  !!Lower bracket value of solution refinement.
  
  complex :: om2
  !!Upper bracket value of solution refinement.

  complex :: omold
  !!Previous frequency value.

  complex :: omega
  !!Complex Frequency under evaluation.
  
  real    :: val
  !!Dispersion Relation Value
  
  integer :: iflag
  !!Flag for Root search
  
  complex,dimension(:),allocatable :: omSafe
  !!Saved frequency Arrays for each root for second parameter scan.

  !Eigenfunctions
  complex, dimension(1:3)       :: ef
  !! Electric Field eigenfluctuation.
  
  complex, dimension(1:3)       :: bf
  !! Magnetic Field eigenfluctuation.

  complex, dimension(1:nspec)     :: ns
  !! Density eigenfluctuation.
  
  complex, dimension(1:3,1:nspec) :: Us
  !!Velocity eigenfluctuation.

  !Heating
  real, dimension(1:nspec) :: Ps
  !!Power into/out of species/component in one wave period.
  
  real, dimension(1:4,1:nspec) :: Ps_split
  !!Power into/out of species/components, broken into different contributions.
  !! Deprecated.
  
  !>>>GGH: 1/18/23
  real, dimension(1:6,1:nspec) :: Ps_split_new
  !!Power into/out of species/componets.
  !!Corrected LD/TTD calculation (GGH).

  real :: Ew
  !!Wave energy.

  real, dimension(1:6,1:nspec) :: params
  !!Species parameters:
  
  !Assign output name for scan is
  call get_double_out_name(outName,tensorName,fmt,fmt_tnsr,out_type,diff)
  !pointer sw, sw3 (and sw2, sw4 if needed) is assigned in GET_DOUBLE_OUT_NAME 

  do is = 1,2
     if ((scan(is)%style_s)==-1) then
        !Scans with multiple components
        if ((scan(is)%type_s)==0) then
           !k_0 -> k_1; 
           !k_0= current k
           !kperp_1 = range_i, kpar_1 = range_f
           kpari=kpar
           kperpi=kperp
        elseif ((scan(is)%type_s)==1) then
           !theta_0 -> theta_1
           !atan(kperp/kpar) -> theta_0
           !swf ->  theta_1 in degrees
           theta=atan(kperp/kpar)!initial theta value
           ki=kpar/cos(theta)
        elseif ((scan(is)%type_s)==2) then
           !k_0-> k_1 along fixed (current) angle
           !k_0 -> k_1; 
           !k_0= current k
           !|k_1|=swf
           theta = atan(kperp/kpar)
           ki=kpar/cos(theta)
        endif
     endif
  enddo

  !Open the (nroot_max) output files

  !May want to have a netCDF/h5 file output option.
  !Look into AGK I/O module
  do ii = 1,nroot_max
     call get_unused_unit(out_unit(ii))
     write(writeName,'(2a,i0)')trim(outName),'.mode',ii
     open(unit=out_unit(ii),file=trim(writeName),status='replace')
     if ((scan(1)%tensor_s).and.(scan(2)%tensor_s)) then
        call get_unused_unit(out_unit_2(ii))
        write(writeName,'(2a,i0)')trim(tensorName),'.mode',ii
        open(unit=out_unit_2(ii),file=trim(writeName),status='replace')
     endif
  enddo
  
  !Allocate variable for last solution and copy in initial values
  allocate(omlast(nroot_max))
  allocate(omSafe(nroot_max)); omSafe=cmplx(0.,0.)
  do ii=1,nroot_max
     omlast(ii)=cmplx(wroots(1,ii),wroots(2,ii))
  enddo 

  !Parameter 1 scan
  !Output %n_scan steps, with %n_res steps inbetween each output
  do kk = 0, (scan(1)%n_scan*scan(1)%n_res) 

     !Advance scanned parameter values
     if ((scan(1)%style_s)==-1)then
        if ((scan(1)%type_s)==0) then
           if (scan(1)%log_scan) then
              !k0->k1
              sw=10.**(log10(kperpi)+diff(1,1)*real(kk))    
              sw2=10.**(log10(kpari)+diff(1,2)*real(kk))    
           else
              sw=(kperpi)+diff(1,1)*real(kk)    
              sw2=(kpari)+diff(1,2)*real(kk)    
           endif
        elseif ((scan(1)%type_s)==1) then
           !theta_0->theta_1
           if (scan(1)%log_scan) then
              theta_q=10.**(log10((theta))+diff(1,1)*real(kk))    
           else
              theta_q=(((theta))+diff(1,1)*real(kk))    
           endif
           sw=(ki*sin(theta_q))!kperp
           sw2=(ki*cos(theta_q))!kpar
        elseif ((scan(1)%type_s)==2) then
           !k along contant theta
           if (scan(1)%log_scan) then
              sw=10.**(log10(ki*sin(theta))+diff(1,1)*real(kk))    
              sw2=10.**(log10(ki*cos(theta))+diff(1,1)*real(kk))    
           else
              sw=(ki*sin(theta))+diff(1,1)*real(kk)    
              sw2=(ki*cos(theta))+diff(1,1)*real(kk)    
           endif
        endif
     else
        if (scan(1)%log_scan) then
           sw=10.**(log10(scan(1)%range_i)+diff(1,1)*real(kk))    
        else
           sw=(scan(1)%range_i)+diff(1,1)*real(kk)    
        endif
     endif

     !Root Scan....
     do ii=1,nroot_max
        !Bracket values for root search
        omold=omlast(ii)
        om1=omold*(1.-prec)
        om2=omold*(1.+prec)

        iflag=0
        !Find New Omega Value
        omega=rtsec(disp,om1,om2,tol,iflag)
        
        !Save Root
        omlast(ii)=omega      

     enddo

     !Transition from parameter 1 scan to parameter 2 scan
     if (mod(kk,scan(1)%n_res)==0) then
        if (writeOut) &
             write(*,'(a,es11.4)')&
             'parameter 1: ',sw
        !Save roots
        do ii = 1,nroot_max
           omSafe(ii) = omlast(ii)
        enddo
    
        !Parameter 2 scan
        !Output %n_scan steps, with %n_res steps inbetween each output
        do jj = 0, (scan(2)%n_scan*scan(2)%n_res) 
           !Advance scanned parameter values
           if ((scan(2)%style_s)==-1)then
              if ((scan(2)%type_s)==0) then
                 if (scan(2)%log_scan) then
                    !k0->k1
                    sw3=10.**(log10(kperpi)+diff(2,1)*real(jj))    
                    sw4=10.**(log10(kpari)+diff(2,2)*real(jj))    
                 else
                    sw3=(kperpi)+diff(2,1)*real(jj)    
                    sw4=(kpari)+diff(2,2)*real(jj)    
                 endif
              elseif ((scan(2)%type_s)==1) then
                 !theta_0->theta_1
                 if (scan(2)%log_scan) then
                    theta_q=10.**(log10((theta))+diff(2,1)*real(jj))    
                 else
                    theta_q=(((theta))+diff(2,1)*real(jj))    
                 endif
                 sw3=(ki*sin(theta_q))!kperp
                 sw4=(ki*cos(theta_q))!kpar
              elseif ((scan(2)%type_s)==2) then
                 !k along contant theta
                 if (scan(2)%log_scan) then
                    sw3=10.**(log10(ki*sin(theta_q))+diff(2,1)*real(jj))    
                    sw4=10.**(log10(ki*cos(theta_q))+diff(2,1)*real(jj))    
                 else
                    sw3=(ki*sin(theta_q))+diff(2,1)*real(jj)    
                    sw4=(ki*cos(theta_q))+diff(2,1)*real(jj)    
                 endif
              endif
           else
              if (scan(2)%log_scan) then
                 sw3=10.**(log10(scan(2)%range_i)+diff(2,1)*real(jj))    
              else
                 sw3=(scan(2)%range_i)+diff(2,1)*real(jj)    
              endif
           endif

           !construct params array for formated output
           if (mod(jj,scan(2)%n_res)==0) then
              do ii = 1,nspec
                 params(1,ii) = spec(ii)%tau_s
                 params(2,ii) = spec(ii)%mu_s
                 params(3,ii) = spec(ii)%alph_s
                 params(4,ii) = spec(ii)%q_s
                 params(5,ii) = spec(ii)%D_s
                 params(6,ii) = spec(ii)%vv_s
              enddo
           endif
           
           !Root Scan....
           do ii=1,nroot_max
              !Bracket values for root search
              omold=omlast(ii)
              om1=omold*(1.-prec)
              om2=omold*(1.+prec)
              
              iflag=0
              !Find New Omega Value
              omega=rtsec(disp,om1,om2,tol,iflag)
              
              !Save Root
              omlast(ii)=omega        
              
              !Output value every %n_res steps
              if (mod(jj,scan(2)%n_res)==0) then
                 !Calculate eigenfns, heating only every %n_res steps
                 
                 if ((scan(2)%eigen_s).or.((scan(2)%heat_s))) then
                    val=abs(disp(omega))
                    !>>>GGH: 1/18/23
                    call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,Ps_split_new,scan(2)%eigen_s,scan(2)%heat_s)
                    !<<<GGH: 1/18/23
                    if (abs(real(omega)).lt.1.E-7) then
                       Ps=-0.1
                       Ps_split=1.0
                       !>>>GGH: 1/18/23
                       Ps_split_new=1.0
                       !<<<GGH: 1/18/23
                    endif

                 endif
                 
                 select case(out_type)
                 case(0) !Om, Eigen, Heating                    
                    if (low_n) then
                       !>>>GGH: 1/18/23
                       if (new_low_n) then  !Greg's New low_n calculation to separate LD/TTD
                          write(out_unit(ii),fmt)&
                               kperp,kpar,betap,vtp,&
                               omega,&            
                               bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                               Ps(1:nspec),Ps_split_new(1:6,1:nspec),&
                               params(1:6,1:nspec)
                       else  !Old version of low_n using Stix/Quataert version of LD/TTD 
                          !<<<GGH: 1/18/23
                          write(out_unit(ii),fmt)&
                               kperp,kpar,betap,vtp,&
                               omega,&            
                               bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                               Ps(1:nspec),Ps_split(1:4,1:nspec),&
                               params(1:6,1:nspec)
                          !>>>GGH: 1/18/23
                       endif
                       !<<<GGH: 1/18/23
                    else
                       write(out_unit(ii),fmt)&
                            kperp,kpar,betap,vtp,&
                            omega,&            
                            bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                            Ps(1:nspec),&
                            params(1:6,1:nspec),Ew
                    endif
                 case(1) !Om, Eigen
                    write(out_unit(ii),fmt)&
                         kperp,kpar,betap,vtp,&
                         omega,&            
                         bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                         params(1:6,1:nspec)
                 case(2) !Om, Heating
                    if (low_n) then
                       !>>>GGH: 1/18/23
                       if (new_low_n) then  !Greg's New low_n calculation to separate LD/TTD
                          write(out_unit(ii),fmt)&
                               kperp,kpar,betap,vtp,&
                               omega,&            
                               Ps(1:nspec),Ps_split_new(1:6,1:nspec),&
                               params(1:6,1:nspec)
                       else  !Old version of low_n using Stix/Quataert version of LD/TTD 
                          !<<<GGH: 1/18/23
                          write(out_unit(ii),fmt)&
                               kperp,kpar,betap,vtp,&
                               omega,&            
                               Ps(1:nspec),Ps_split(1:4,1:nspec),&
                               params(1:6,1:nspec)
                          !>>>GGH: 1/18/23
                       endif
                       !<<<GGH: 1/18/23
                    else
                       write(out_unit(ii),fmt)&
                            kperp,kpar,betap,vtp,&
                            omega,&            
                            Ps(1:nspec),&
                            params(1:6,1:nspec),Ew
                    endif
                 case(3) !Om
                    write(out_unit(ii),fmt)&
                         kperp,kpar,betap,vtp,&
                         omega,&            
                         params(1:6,1:nspec)             
                 end select
              endif

              if ((scan(1)%tensor_s).and.(scan(2)%tensor_s)) then
                 val=abs(disp(omega))
                 write(out_unit_2(ii),fmt_tnsr)&
                      kperp,kpar,betap,vtp,&
                      susc(1:nspec,1,1),susc(1:nspec,1,2),susc(1:nspec,1,3),&
                      susc(1:nspec,2,2),susc(1:nspec,2,3),&
                      susc(1:nspec,3,3),&
                      params(1:6,1:nspec)
                 !(1,1) xx; (1,2) xy; (1,3) xz;
                 !(2,1) yx; (2,2) yy; (2,3) yz;
                 !(3,1) zx; (3,2) zy; (3,3) zz;

              endif

           enddo
        enddo     !End Parameter 2 Scan
        jj=0
        !Return scanned parameter 2 value
           if ((scan(2)%style_s)==-1)then
              if ((scan(2)%type_s)==0) then
                 if (scan(2)%log_scan) then
                    !k0->k1
                    sw3=10.**(log10(kperpi)+diff(2,1)*real(jj))    
                    sw4=10.**(log10(kpari)+diff(2,2)*real(jj))    
                 else
                    sw3=(kperpi)+diff(2,1)*real(jj)    
                    sw4=(kpari)+diff(2,2)*real(jj)    
                 endif
              elseif ((scan(2)%type_s)==1) then
                 !theta_0->theta_1
                 if (scan(2)%log_scan) then
                    theta_q=10.**(log10((theta))+diff(2,1)*real(jj))    
                 else
                    theta_q=(((theta))+diff(2,1)*real(jj))    
                 endif
                 sw3=(ki*sin(theta_q))!kperp
                 sw4=(ki*cos(theta_q))!kpar
              elseif ((scan(2)%type_s)==2) then
                 !k along contant theta
                 if (scan(2)%log_scan) then
                    sw3=10.**(log10(ki*sin(theta))+diff(2,1)*real(jj))    
                    sw4=10.**(log10(ki*cos(theta))+diff(2,1)*real(jj))    
                 else
                    sw3=(ki*sin(theta))+diff(2,1)*real(jj)    
                    sw4=(ki*cos(theta))+diff(2,1)*real(jj)    
                 endif
              endif
           else
              if (scan(2)%log_scan) then
                 sw3=10.**(log10(scan(2)%range_i)+diff(2,1)*real(jj))    
              else
                 sw3=(scan(2)%range_i)+diff(2,1)*real(jj)    
              endif
           endif
           !Recall Saved roots
           do ii = 1,nroot_max
              omlast(ii)=omsafe(ii)
              write(out_unit(ii),*)
              if ((scan(1)%tensor_s).and.(scan(2)%tensor_s)) &
                   write(out_unit_2(ii),*)
              write(*,'(a,i3,a,2es14.4)')'Root ',ii,': ',omlast(ii)
           enddo
        endif
     
     enddo  !End Parameter 1 Scan

end subroutine om_double_scan

!-=-=-=-=-
!Kristopher Klein, 2015
!----------------------
!-=-=-=-=-
subroutine radial_scan
  !!Calculates normal modes of the system
  !!for a specified radial solar wind model
  !! (under development)
  use vars, only : nRad, radius, writeOut, nroot_max
  use vars, only : wroots, spec, rad_spec, nspec, k_scan
  use vars, only : vtp_rad, beta_rad, sw, sw2, sw3, sw4
  use vars, only : kperp, kpar, betap, vtp, radial_eigen
  use vars, only : dataName, outputName, radial_heating
  use vars, only : rad_scan, pi
  use functions, only : get_unused_unit
  implicit none
  
  !Local
  integer :: ir
  !!Radial index.
  
  integer :: ii
  !!Mode index.

  integer :: is
  !!Species/component index.
  
  integer :: ij
  !!Primary wavevector index.

  integer :: ik
  !!Secondary wavevector index.
  
  character(150) :: outName
  !! String for identifying output files.

  character(150) :: writeName
  !! Full output name.
  
  character(100)  :: fmt
  !! Output format for dispersion relation.
  
  integer ,dimension(1:nroot_max) :: out_unit
  !! Main output unit.

  integer :: out_type
  !! logic for outputing supplementary eigenfunction and heating
  !! calculations.
  !! 0: Frequency, Eigenfuction, Heating.
  !! 1: Frequency, Eigenfuction.
  !! 2: Frequency, Heating.
  !! 3: Frequency.
  
  complex,dimension(:),allocatable :: omlast
  !!Arrays with complex frequency for each solution.

  complex,dimension(:,:),allocatable :: omSafe
  !!Saved frequency Arrays for each root for second parameter scan.
  
  logical :: mod_write =.FALSE.
  !!For only writing out every n_scan, not n_scan * n_res steps
  !!species parameters: useful for outputing data (INVESTIGATE THIS...)
  
  real, dimension(1:6,1:nspec) :: params
  !!Species/component parameters:

  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)

  real :: ki
  !! Current value of \( |k| \).
  
  !Assign outname
  write(outName,'(4a)')&
       'data/',trim(dataName),'/',&
       trim(outputName)

  !Open the (nroot_max) output files
  do ii = 1,nroot_max
     call get_unused_unit(out_unit(ii))
     write(writeName,'(2a,i0)')trim(outName),'.mode',ii
     open(unit=out_unit(ii),file=trim(writeName),status='replace')
  enddo

  !!ADD split power calculations here... (KGK)
  
  !Specify output formatting
  if (radial_eigen) then
     if (radial_heating) then
        write(fmt,'(a,i0,a)')'(7es14.4,12es14.4,',15*nspec,'es14.4)'
        out_type=0
     else
        write(fmt,'(a,i0,a)')'(7es14.4,12es14.4,',14*nspec,'es14.4)'
        out_type=1
     endif
  else
     if (radial_heating) then
        write(fmt,'(a,i0,a)')'(7es14.4,',7*nspec,'es14.4)'
        out_type=2
     else
        write(fmt,'(a,i0,a)')'(7es14.4,',6*nspec,'es14.4)'
        out_type=3
     endif
  endif

  !Allocate variable for last solution and copy in initial values
  allocate(omlast(nroot_max))
  allocate(omSafe(1:2,nroot_max)); omSafe=cmplx(0.,0.)
  do ii=1,nroot_max
     omlast(ii)=cmplx(wroots(1,ii),wroots(2,ii))
     omSafe(1,ii) = omlast(ii);omSafe(2,ii) = omlast(ii)
     if (writeOut) &
          write(*,'(a,i3,a,2es14.4)')&
          'Root ',ii,' In: ',wroots(1,ii),wroots(2,ii)
  enddo

  !Scan through radial trajectory in parameter space
  do ir = 0, nRad
     !Set global parameters
     betap=beta_rad(ir)
     vtp=vtp_rad(ir)
     !Set Species Parameters
     do is = 1,nspec        
        spec(is)%tau_s  = rad_spec(is,ir)%tau_s
        spec(is)%mu_s   = rad_spec(is,ir)%mu_s
        spec(is)%alph_s = rad_spec(is,ir)%alph_s
        spec(is)%q_s    = rad_spec(is,ir)%q_s
        spec(is)%D_s    = rad_spec(is,ir)%D_s
        spec(is)%vv_s   = rad_spec(is,ir)%vv_s

        params(1,is) = spec(is)%tau_s
        params(2,is) = spec(is)%mu_s
        params(3,is) = spec(is)%alph_s
        params(4,is) = spec(is)%q_s
        params(5,is) = spec(is)%D_s
        params(6,is) = spec(is)%vv_s
     enddo

     if (.true.) &
          write(*,'(a,i0,a,es14.4,a)')&
          'Radius(',ir,') = ',radius(ir),' Rs'

     !Scan (or don't) through kspace
     select case(k_scan)
     case(0)
        !simple case of constant kperp, kpar

        !Root Finder in separate subroutine
        call om_radial(omlast,params,out_unit,fmt,out_type,ir,.true.)
        
     case(1)
        !fixed kperp, scan over kpar

        !reset input roots for new scan of k
        do ii = 1,nroot_max
           omLast(ii)=omSafe(1,ii)
        enddo

        do ij = 0,rad_scan(1)%n_scan*rad_scan(1)%n_res

           !Determine write status
           if (mod(ij,rad_scan(1)%n_res)==0) then
              mod_write=.true.
           else
              mod_write=.false.
           endif

           !set kpar
           if (rad_scan(1)%log_scan) then
              kpar=10.**(log10(rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           else
              kpar=((rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           endif

           !Root Finder in separate subroutine
           call om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
           
           if (ij==0) then
              !Save roots for next scan of k
              do ii=1,nroot_max
                 omSafe(1,ii)=omLast(ii)
              enddo
           endif

        enddo
        
        !Spacing between radial steps for contour plotting
        do ii=1,nroot_max
           write(out_unit(ii),*);!write(out_unit(ii),*)
        enddo

     case(2)
        !fixed kpar, scan over kperp

        !reset input roots for new scan of k
        do ii = 1,nroot_max
           omLast(ii)=omSafe(1,ii)
        enddo

        do ij = 0,rad_scan(1)%n_scan*rad_scan(1)%n_res
           !Determine write status
           if (mod(ij,rad_scan(1)%n_res)==0) then
              mod_write=.true.
           else
              mod_write=.false.
           endif

           !set kperp
           if (rad_scan(1)%log_scan) then
              kperp=10.**(log10(rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           else
              kperp=((rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           endif

           !Root Finder in separate subroutine
           call om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
           
           if (ij==0) then
              !Save roots for next scan of k
              do ii=1,nroot_max
                 omSafe(1,ii)=omLast(ii)
              enddo
           endif

        enddo
        
        !Spacing between radial steps for contour plotting
        do ii=1,nroot_max
           write(out_unit(ii),*);!write(out_unit(ii),*)
        enddo

     case(3)
        !fixed theta, scan over k
        theta = atan(kperp/kpar)
        !ki=kpar/cos(theta)

        !reset input roots for new scan of k
        do ii = 1,nroot_max
           omLast(ii)=omSafe(1,ii)
        enddo

        do ij = 0,rad_scan(1)%n_scan*rad_scan(1)%n_res
           !Determine write status
           if (mod(ij,rad_scan(1)%n_res)==0) then
              mod_write=.true.
           else
              mod_write=.false.
           endif

           !set k
           if (rad_scan(1)%log_scan) then
              ki=10.**(log10(rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           else
              ki=((rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           endif          
           kperp = ki*sin(theta)
           kpar  = ki*cos(theta)

           !Root Finder in separate subroutine
           call om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
           
           if (ij==0) then
              !Save roots for next scan of k
              do ii=1,nroot_max
                 omSafe(1,ii)=omLast(ii)
              enddo
           endif

        enddo
        
        !Spacing between radial steps for contour plotting
        do ii=1,nroot_max
           write(out_unit(ii),*);!write(out_unit(ii),*)
        enddo

     case(4)
        !fixed k, scan over theta
        
        theta = atan(kperp/kpar)
        ki=kpar/cos(theta)

        !reset input roots for new scan of k
        do ii = 1,nroot_max
           omLast(ii)=omSafe(1,ii)
        enddo

        do ij = 0,rad_scan(1)%n_scan*rad_scan(1)%n_res
           !Determine write status
           if (mod(ij,rad_scan(1)%n_res)==0) then
              mod_write=.true.
           else
              mod_write=.false.
           endif

           !set k
           if (rad_scan(1)%log_scan) then
              theta=10.**(log10(rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           else
              theta=((rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
           endif          
           kperp = ki*sin(theta*pi/180.)
           kpar  = ki*cos(theta*pi/180.)

           !Root Finder in separate subroutine
           call om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
           
           if (ij==0) then
              !Save roots for next scan of k
              do ii=1,nroot_max
                 omSafe(1,ii)=omLast(ii)
              enddo
           endif

        enddo
        
        !Spacing between radial steps for contour plotting
        do ii=1,nroot_max
           write(out_unit(ii),*);!write(out_unit(ii),*)
        enddo

     case(5)
        !plane scan over (kperp, kpar)
        
        !reset input roots for new scan of kplane
        do ii = 1,nroot_max
           omLast(ii)=omSafe(1,ii)
        enddo
        
        !kpar loop
        do ik = 0,rad_scan(2)%n_scan*rad_scan(2)%n_res
           !set kpar
           if (rad_scan(2)%log_scan) then
              kpar=10.**(log10(rad_scan(2)%range_i)+rad_scan(2)%diff*real(ik))
           else
              kpar=((rad_scan(2)%range_i)+rad_scan(2)%diff*real(ik))
           endif          

           write(*,'(a,es14.4)')&
                'kpar rho_p :',kpar

           !reset input roots for new scan of kperp
           if (ik.ne.0) then
              do ii = 1,nroot_max
                 omLast(ii)=omSafe(2,ii)
              enddo
           endif

           !kperp loop
           do ij = 0,rad_scan(1)%n_scan*rad_scan(1)%n_res
              !set kperp
              if (rad_scan(1)%log_scan) then
                 kperp=10.**(log10(rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
              else
                 kperp=((rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
              endif
              
              !Determine write status
              if ((mod(ij,rad_scan(1)%n_res)==0).and.&
                   (mod(ik,rad_scan(2)%n_res)==0))then
                 mod_write=.true.
              else
                 mod_write=.false.
              endif

              !Root Finder in separate subroutine
              call om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
              
              if ((ij==0).and.(ik==0)) then
                 !Save roots for next scan of kplane
                 do ii=1,nroot_max
                    omSafe(1,ii)=omLast(ii)
                 enddo
              endif

              if (ij==0) then
                 !Save roots for next scan of kpar
                 do ii=1,nroot_max
                    omSafe(2,ii)=omLast(ii)
                 enddo
              endif
              
           enddo
           !Spacing between radial steps for contour plotting
           if ((mod(ik,rad_scan(2)%n_res)==0))then
              do ii=1,nroot_max
                 write(out_unit(ii),*)
              enddo
           endif
        enddo
        
        !Spacing between radial steps for contour plotting
        do ii=1,nroot_max
           write(out_unit(ii),*)
        enddo

     case(6)
        !plane scan over (k, theta)
        
        !reset input roots for new scan of kplane
        do ii = 1,nroot_max
           omLast(ii)=omSafe(1,ii)
        enddo
        
        !k loop
        do ik = 0,rad_scan(2)%n_scan*rad_scan(2)%n_res
           !set k
           if (rad_scan(2)%log_scan) then
              ki=10.**(log10(rad_scan(2)%range_i)+rad_scan(2)%diff*real(ik))
           else
              ki=((rad_scan(2)%range_i)+rad_scan(2)%diff*real(ik))
           endif

           write(*,'(a,es14.4)')&
                'k rho_p :',ki

           !reset input roots for new scan of kperp
           if (ik.ne.0) then
              do ii = 1,nroot_max
                 omLast(ii)=omSafe(2,ii)
              enddo
           endif

           !theta loop
           do ij = 0,rad_scan(1)%n_scan*rad_scan(1)%n_res
              !set theta
              if (rad_scan(1)%log_scan) then
                 theta=10.**(log10(rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
              else
                 theta=((rad_scan(1)%range_i)+rad_scan(1)%diff*real(ij))
              endif              
              kperp = ki*sin(theta*pi/180.)
              kpar  = ki*cos(theta*pi/180.)

              !Determine write status
              if ((mod(ij,rad_scan(1)%n_res)==0).and.&
                   (mod(ik,rad_scan(2)%n_res)==0))then
                 mod_write=.true.
              else
                 mod_write=.false.
              endif

              !Root Finder in separate subroutine
              call om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
              
              if ((ij==0).and.(ik==0)) then
                 !Save roots for next scan of kplane
                 do ii=1,nroot_max
                    omSafe(1,ii)=omLast(ii)
                 enddo
              endif

              if (ij==0) then
                 !Save roots for next scan of kpar
                 do ii=1,nroot_max
                    omSafe(2,ii)=omLast(ii)
                 enddo
              endif
              
           enddo
           !Spacing between radial steps for contour plotting
           if ((mod(ik,rad_scan(2)%n_res)==0))then
              do ii=1,nroot_max
                 write(out_unit(ii),*)
              enddo
           endif
        enddo
        
        !Spacing between radial steps for contour plotting
        do ii=1,nroot_max
           write(out_unit(ii),*)
        enddo

        !set initial k
        !kperp = k_1*sin(theta_1*pi/180.)
        !kpar = k_1*cos(theta_1*pi/180.)
        !Save ranges
        !sw => k_1; sw2 => k_2
        !sw3 => theta_1; sw4 => theta_2
        !Root Finder in separate subroutine
  end select

enddo

end subroutine radial_scan

!-=-=-=-=-
!-=-=-=-=-
subroutine om_radial(omlast,params,out_unit,fmt,out_type,ir,mod_write)
  !!Scans roots at a given location in parameter space.
  !!Used in conjunction with [[radial_scan(subroutine)]].
  use vars, only: nspec, nroot_max, radial_heating, radial_eigen
  use vars, only: radius, betap, vtp, kperp, kpar
  !Passed

  complex,dimension(1:nspec) :: omlast
  !!Arrays with complex frequency for each solution.

  real, dimension(1:6,1:nspec) :: params
  !!Species parameters:

  integer ,dimension(1:nroot_max) :: out_unit
  !! Main output unit.

  character(100)  :: fmt
  !! Output format for dispersion relation.

  integer :: out_type
  !! logic for outputing supplementary eigenfunction and heating
  !! calculations.
  !! 0: Frequency, Eigenfuction, Heating.
  !! 1: Frequency, Eigenfuction.
  !! 2: Frequency, Heating.
  !! 3: Frequency.

  integer :: ir
  !!Radial index.

  logical :: mod_write
  !!For only writing out every n_scan, not n_scan * n_res steps
  !!species parameters: useful for outputing data (INVESTIGATE THIS...)
  
  !Local
  integer :: ii
  !!Solution index.

  complex :: om1
  !!Lower bracket value of solution refinement.
  
  complex :: om2
  !!Upper bracket value of solution refinement.

  complex :: omold
  !!Previous frequency value.

  complex :: omega
  !!Complex Frequency under evaluation.
  
  real    :: val
  !!Dispersion Relation Value
  
  integer :: iflag
  !!Flag for Root search

  !Eigenfunctions
  complex, dimension(1:3)       :: ef
  !! Electric Field eigenfluctuation.
  
  complex, dimension(1:3)       :: bf
  !! Magnetic Field eigenfluctuation.

  complex, dimension(1:nspec)     :: ns
  !! Density eigenfluctuation.
  
  complex, dimension(1:3,1:nspec) :: Us
  !!Velocity eigenfluctuation.

  !Heating
  real, dimension(1:nspec) :: Ps
  !!Power into/out of species/component in one wave period.
  
  real, dimension(1:4,1:nspec) :: Ps_split
  !!Power into/out of species/components, broken into different contributions.
  !! Deprecated.
  
  !>>>GGH: 1/18/23
  real, dimension(1:6,1:nspec) :: Ps_split_new
  !!Power into/out of species/componets.
  !!Corrected LD/TTD calculation (GGH).

  real :: Ew
  !!Wave energy.
  
  !Root Scan....
  do ii=1,nroot_max
     !Bracket values for root search
        omold=omlast(ii)
        om1=omold*(1.-prec)
        om2=omold*(1.+prec)

        iflag=0
        !Find New Omega Value
        omega=rtsec(disp,om1,om2,tol,iflag)
        
        !Save Root
        omlast(ii)=omega        

        if (mod_write) then
           !Calculate eigenfunctions and heating rates
           if ((radial_heating).or.(radial_eigen)) then
              val=abs(disp(omega))
              call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,Ps_split_new,radial_eigen,radial_heating)
           endif
           
           !Output results
           if (out_type==0) & !Om, Eigen, Heating
                write(out_unit(ii),fmt)&
                radius(ir),kperp,kpar,betap,vtp,&
                omega,&            
                bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                Ps(1:nspec),&
                params(1:6,1:nspec),Ew
           if (out_type==1) & !Om, Eigen
                write(out_unit(ii),fmt)&
                radius(ir),kperp,kpar,betap,vtp,&
                omega,&            
                bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                params(1:6,1:nspec)
           if (out_type==2) & !Om, Heating
                write(out_unit(ii),fmt)&
                radius(ir),kperp,kpar,betap,vtp,&
                omega,&            
                Ps(1:nspec),&
                params(1:6,1:nspec),Ew
           if (out_type==3) & !Om
                write(out_unit(ii),fmt)&
                radius(ir),kperp,kpar,betap,vtp,&
                omega,&            
                params(1:6,1:nspec)                        
        endif
     enddo

end subroutine om_radial

!-=-=-=-=-
!------------------------------------------------------------------------------
!Greg Howes, 2010; Kristopher Klein, 2015
!------------------------------------------------------------------------------
!-=-=-=-=-=-
subroutine calc_eigen(omega,electric,magnetic,vmean,ns,Ps,&
     Ps_split,Ps_split_new,eigen_L,heat_L)
  !!  Calculates the electric and magnetic fields as well as species
  !!     velocities and density fluctuations for identified wave 
  !!     as well as the power emission or absorption.
  use vars, only : new_low_n
!<<<GGH: 1/18/23
  use vars, only : spec,betap,vtp,kperp,kpar,nspec,susc,lam,low_n,susc_low
  implicit none
  !!Passed!!
  complex :: omega
  !!Complex Frequency of identified wave.
  
  !Eigenfunctions.
  complex, dimension(1:3), intent(out)       :: electric
  !! Complex electric field fluctuations
  !! (Ex, Ey, Ez)

  complex, dimension(1:3), intent(out)       :: magnetic
  !! Complex magnetic field fluctuations
  !! (Bx, By, Bz)

  complex, dimension(1:3) :: electric_xy
  !! Components of electric field:
  !! (Ex, Ey, 0)
  
  complex, dimension(1:3) ::  electric_z
  !! Components of electric field:
  !! (0, 0, Ez)
  
  complex, dimension(1:nspec), intent(out)     :: ns
  !!Complex Density fluctuations.
  
  complex, dimension(1:3,1:nspec), intent(out) :: vmean
  !!Complex Velocity fluctuations:
  !! (Uxs, Uys, Uzs)
  
  !Heating
  real, dimension(1:nspec), intent(out) :: Ps
  !!Power into/out of species/components.
  
  real, dimension(1:4,1:nspec), intent(out) :: Ps_split
  !!Power into/out of species

  real, dimension(1:6,1:nspec), intent(out) :: Ps_split_new
  !!<<<GGH: 1/18/23
  !!Power into/out of species/components from LD, TTD, CD.

  logical :: eigen_L
  !!Logical for calculating eigenvalues
  
  logical :: heat_L
  !!Logical for calculating heating

  !Local
  real :: sume
  !!Temporary sum.

  integer :: ii
  !! 1st Tensor index.

  integer :: j
  !! 2nd Tensor index.

  integer :: jj
  !! Species index.
  
  complex :: temp1
  !! Temp. value for complex frequency.
  
  complex, dimension(nspec,3,3) :: susca
  !!The Anti-Hermitian Susceptibility Tensor.
  
  complex, dimension(3,3) :: susch
  !!The Hermitian Susceptibility Tensor.

  complex, dimension(3,3) :: suschold
  !!The Saved Hermitian Susceptibility Tensor.
  
  complex, dimension(3,3) :: dsusch
  !!The Difference in the Hermitian Susceptibility Tensors.
  
  complex, dimension(nspec,3) :: term
  !! Summation terms for heating calculation.
  
  complex, dimension(3) :: term1
  !! Further summation terms for heating calculation.

  real :: ewave
  !! Wave Energy.

  !CALCULATE FIELDS FLUCTUATIONS==========================================
  !Calculate Electric Fields
  !This seems to be correct --GGH 20 JUN 2008
  electric(1) = cmplx(1.,0.)
  electric(3)=-electric(1)*(lam(2,1)*lam(3,2)-lam(3,1)*lam(2,2))
  electric(3)= electric(3)/(lam(2,3)*lam(3,2)-lam(3,3)*lam(2,2))
  electric(2) = -electric(3)*lam(3,3) - electric(1)*lam(3,1)
  electric(2) = electric(2)/lam(3,2)

  !Calculate Magnetic Fields
  !NOTE: There appears to be an overall sign error here --GGH 20 JUN 2008
  !NOTE: Sign error corrected by GGH on 2 SEP 2010 (Verified by EQ 6 JUN 2009)
  !KGK: need to include a factor of sqrt(T_perp p/T_par p) in the case of 
  !proton temperature anisotropy, as the thermal speed 
  !since vtp is for the parallel thermal speed (19 DEC 2013)
  magnetic(1) = -1.* kpar*electric(2)/(omega*vtp*sqrt(spec(1)%alph_s))
  magnetic(2) = -1.* (kperp*electric(3) - kpar*electric(1))/(omega*vtp*sqrt(spec(1)%alph_s))
  magnetic(3) = kperp*electric(2)/(omega*vtp*sqrt(spec(1)%alph_s))
  
  !If (scan(is)%eigen_s) loop
  if (eigen_L) then
     !CALCULATE VELOCITY FLUCTUATIONS========================================
     !vmean is the velocity perturbutation due to the wave for each species
     !    normalized to the electron velocity perturbation
     ! This is the mean velocity normalized to the Alfven velocity
     !KGK: 1-28-2015: Added the drift velocity effects
     vmean(:,:)=0.
     do j=1,3!x,y,z
        do jj = 1,nspec !Species velocity fluctuations
           vmean(j,jj) = -(spec(jj)%Q_s/spec(jj)%D_s)*cmplx(0.,1.)*&
                omega/sqrt(betap)*vtp**2. *sqrt(spec(1)%alph_s)* &
                sum(electric(:)*susc(jj,j,:))
        enddo
     enddo
     
     !CALCULATE DENSITY FLUCTUATIONS========================================
     ! This is ns/ns0
     !KGK: 8-28-2013: switched from real(omega) to the complex valued omega.
     !KGK: 12-19-2013: Added the sqrt(T_perp/T_par) factor
     !KGK: 1-28-2015: Added the drift velocity Dopplar shift
     do jj=1,nspec
        ns(jj) = (vmean(1,jj)*kperp+vmean(3,jj)*kpar)/&
             ((omega-kpar * spec(jj)%vv_s/sqrt(betap*spec(1)%alph_s))&
             *sqrt(betap*spec(1)%alph_s))
     enddo
     
     !EndIf (scan(is)%eigen_s) loop
  endif

  !If (scan(is)%heat_s) loop
  !Greg Howes, 2006; Kristopher Klein, 2015
  if (heat_L) then
     !CALCULATE ELECTRON AND ION HEATING======================================
     temp1 = cmplx(real(omega),0.)
     temp1 = disp(temp1)
          
     do ii = 1, 3 !tensor index
        do j = 1, 3 !tensor index
           do jj = 1, nspec !species index
              susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                   (susc(jj,ii,j) - conjg(susc(jj,j,ii)))
           enddo
           suschold(ii,j) = 0.5*(sum(susc(:,ii,j)) + &
                sum(conjg(susc(:,j,ii))))
        enddo
     enddo
     
     term(:,:)=0.
     term1(:)=0.
     do ii = 1, 3
        do jj = 1, nspec
           term(jj,ii) = sum(conjg(electric(:))*susca(jj,:,ii))     
        enddo
     enddo

     Ps = 0.
     do jj = 1, nspec
        Ps(jj) = sum(term(jj,:)*electric(:))
     enddo
     
     temp1 = disp(cmplx(real(omega*1.000001),0.))
     
     do ii = 1, 3
        do j = 1, 3
           susch(ii,j) = 0.5*(sum(susc(:,ii,j)) + &
                sum(conjg(susc(:,j,ii))))
           dsusch(ii,j)=(1.000001*susch(ii,j)-suschold(ii,j))/0.000001
        enddo
     enddo
          
     ewave = 0.
     do ii = 1, 3
        term1(ii) = sum(conjg(electric(:))*dsusch(:,ii))
     enddo
     
     ewave = sum(term1(:)*electric(:)) + sum(magnetic(:)*conjg(magnetic(:)))

     Ps = Ps/ewave
  
     !LD, TTD, and CD calculation
     if (low_n) then
     !>>>GGH: 1/18/23
     if (new_low_n) then  !Greg's New low_n calculation to separate LD/TTD
       !N=0
        do ii = 1, 3 !tensor index
           do j = 1, 3 !tensor index
              do jj = 1, nspec !species index
                 susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                      (susc_low(jj,ii,j,0) - conjg(susc_low(jj,j,ii,0)))
              enddo
           enddo
        enddo

        !Initialize Ps_split_new
        Ps_split_new(:,:) = 0.
        
        !chi_yy  (TTD term 1)
        Ps_split_new(1,:) =-0.5*cmplx(0.,1.)*&
             conjg(electric(2))*electric(2)* &
             (susc_low(:,2,2,0)-conjg(susc_low(:,2,2,0)))
        !chi_yz  (TTD term 2)
        Ps_split_new(2,:) =-0.5*cmplx(0.,1.)*&
             (electric(3)*conjg(electric(2))*susc_low(:,2,3,0) - &
             conjg(electric(3))*electric(2)*conjg(susc_low(:,2,3,0)))
        !chi_zy  (LD term 1)
        Ps_split_new(3,:) =-0.5*cmplx(0.,1.)*&
             (electric(2)*conjg(electric(3))*susc_low(:,3,2,0) - &
             conjg(electric(2))*electric(3)*conjg(susc_low(:,3,2,0)))
        !chi_zz  (LD term 2)
        Ps_split_new(4,:) =-0.5*cmplx(0.,1.)*&
             conjg(electric(3))*electric(3)* &
             (susc_low(:,3,3,0)-conjg(susc_low(:,3,3,0)))

        !Total n=0 terms
        term(:,:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = sum(conjg(electric(:))*susca(jj,:,ii))     
           enddo
        enddo        
        Ps_split_new(5,:) = 0.
        do jj = 1, nspec
           Ps_split_new(5,jj) = sum(term(jj,:)*electric(:))
        enddo

        !N=1
        do ii = 1, 3 !tensor index
           do j = 1, 3 !tensor index
              do jj = 1, nspec !species index
                 susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                      (susc_low(jj,ii,j,1) - conjg(susc_low(jj,j,ii,1)))
              enddo
           enddo
        enddo

        !Total n=1 terms, Eperp
        electric_xy=electric; electric_xy(3)=cmplx(0.,0.)
        term(:,:)=0.
        term1(:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = sum(conjg(electric_xy(:))*susca(jj,:,ii))     
           enddo
        enddo        
        Ps_split_new(6,:) = 0.
        do jj = 1, nspec
           Ps_split_new(6,jj) = sum(term(jj,:)*electric_xy(:))
        enddo

        !Normalization             
        Ps_split_new = Ps_split_new/ewave
        
     else  !Old version of low_n using Stix/Quataert version of LD/TTD separation
        !<<<GGH: 1/18/23        
        !N=0
        do ii = 1, 3 !tensor index
           do j = 1, 3 !tensor index
              do jj = 1, nspec !species index
                 susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                      (susc_low(jj,ii,j,0) - conjg(susc_low(jj,j,ii,0)))
              enddo
           enddo
        enddo
        
        !LANDAU DAMPING
        term(:,:)=0.
        term1(:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = conjg(electric(3))*susca(jj,3,ii)
           enddo
        enddo        
        Ps_split(1,:) = 0.
        do jj = 1, nspec
           Ps_split(1,jj) = term(jj,3)*electric(3)
        enddo

        !Transit Time Damping
        term(:,:)=0.
        term1(:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = conjg(electric(2))*susca(jj,2,ii)
           enddo
        enddo        
        Ps_split(2,:) = 0.
        do jj = 1, nspec
           Ps_split(2,jj) = term(jj,2)*electric(2)
        enddo
     
        !Total n=0 terms
        term(:,:)=0.
        term1(:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = sum(conjg(electric(:))*susca(jj,:,ii))     
           enddo
        enddo        
        Ps_split(3,:) = 0.
        do jj = 1, nspec
           Ps_split(3,jj) = sum(term(jj,:)*electric(:))
        enddo

        !N=1
        do ii = 1, 3 !tensor index
           do j = 1, 3 !tensor index
              do jj = 1, nspec !species index
                 susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                      (susc_low(jj,ii,j,1) - conjg(susc_low(jj,j,ii,1)))
              enddo
           enddo
        enddo

        !Total n=1 terms, Eperp
        electric_xy=electric; electric_xy(3)=cmplx(0.,0.)
        term(:,:)=0.
        term1(:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = sum(conjg(electric_xy(:))*susca(jj,:,ii))     
           enddo
        enddo        
        Ps_split(4,:) = 0.
        do jj = 1, nspec
           Ps_split(4,jj) = sum(term(jj,:)*electric_xy(:))
        enddo

        !Normalization             
        Ps_split = Ps_split/ewave

        !>>>GGH: 1/18/23
     endif
     !<<<GGH: 1/18/23
     endif

     !EndIf (scan(is)%heat_s) loop
  endif

end subroutine calc_eigen
!-=-=-=-=-=-
!-=-=-=-=-=-
subroutine set_map_pointers(outName,diff,diff2)
  !! Sets step sizes for all classes of parameter scans
  !! called in for [[map_scan(subroutine)]].
  use vars, only : scan,nspec,dataName,writeOut,pi
  use vars, only : spec,sw,sw2,betap,vtp,kperp,kpar
  implicit none
  !Passed
  real :: diff
  !!Spacing for scanned parameter.

  real :: diff2
  !!Suplemental spacing for scanned parameter.

  character(150) :: outName
  !! String for I/O identification.

  !Local
  character(15)  :: param
  !! String to identify scanned parameter.

  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)

  real :: theta_q
  !! Scanned values of \( \theta \).

  real :: ki
  !! Initial value of \( |k| \).

  real :: kf
  !! Final value of \( |k| \).

  real :: kperpi
  !! Initial value of \( k_{\perp} \).

  real :: kpari
  !! Initial value of \( k_{\parallel} \).
  
  if((scan(1)%style_s).ge.0) then
     if (scan(1)%log_scan) then
        !Log spacing
        diff=(log10(scan(1)%range_f)-log10(scan(1)%range_i))/&
             real(scan(1)%n_scan*scan(1)%n_res)
     else
        !Linear spacing
        diff=((scan(1)%range_f)-(scan(1)%range_i))/&
             real(scan(1)%n_scan*scan(1)%n_res)
     endif
  endif

  if (scan(1)%style_s.eq.-1) then
     !Global parameter scan- 2 components
     if (scan(1)%type_s.eq.0) then
        !k_0 -> k_1; 
        !k_0= current k
        !kperp_1 = range_i, kpar_1 = range_f
        write(param,'(a)')'k'
        sw=>kperp
        sw2=>kpar
        kpari=kpar
        kperpi=kperp
        if (writeOut) &          
             write(*,'(4a,es14.4,a,es14.4,a,es14.4,a,es14.4,a)') &
             'Scan over ',trim(param),' from ',&
             '(kperp,kpar) = (',kperp,',',kpar,') to (',&
             scan(1)%range_i,',',scan(1)%range_f,')'

        if (scan(1)%log_scan) then
           !Log spacing
           diff=(log10(scan(1)%range_i)-log10(kperp))/&
                real(scan(1)%n_scan*scan(1)%n_res)
           diff2=(log10(scan(1)%range_f)-log10(kpar))/&
                real(scan(1)%n_scan*scan(1)%n_res)
        else
           !Linear spacing
           diff=((scan(1)%range_i)-(kperp))/&
                real(scan(1)%n_scan*scan(1)%n_res)
           diff2=((scan(1)%range_f)-(kpar))/&
                real(scan(1)%n_scan*scan(1)%n_res)
        endif
        

     elseif (scan(1)%type_s.eq.1) then
        write(param,'(a)')'theta'
        !theta_0 -> theta_1
        !atan(kperp/kpar) -> theta_0
        !swi ->  theta_1 in degrees
        theta=atan(kperp/kpar)!initial theta value
        ki=kpar/cos(theta)
        sw=>kperp
        sw2=>kpar
        if (writeOut) &          
             write(*,'(4a,es14.4,a,es14.4,a,es14.4,a,es14.4,a)') &
             'Scan over ',trim(param),' from ',&
             'theta = ',theta*180./pi,' to ',&
             scan(1)%range_i

        if (scan(1)%log_scan) then
           !Log spacing
           diff=(log10(pi*(scan(1)%range_i)/180.)-log10(theta))/&
                real(scan(1)%n_scan*scan(1)%n_res)
        else
           !Linear spacing
           diff=((pi*(scan(1)%range_i)/180.)-(theta))/&
                real(scan(1)%n_scan*scan(1)%n_res)
        endif

     elseif (scan(1)%type_s.eq.2) then
        write(param,'(a)')'k_fixed_theta'
        !k_0-> k_1 along fixed (current) angle
        !k_0 -> k_1; 
        !k_0= current k
        !|k_1|=swf
        theta = atan(kperp/kpar)
        ki=kpar/cos(theta)
        sw=>kperp
        sw2=>kpar

        if (writeOut) &          
             write(*,'(4a,es14.4,a,es14.4,a,es14.4)') &
             'Scan over ',trim(param),' from ',&
             '|k| = ',ki,' to ',scan(1)%range_f,&
             ' at fixed theta = ',theta * 180./pi

        if (scan(1)%log_scan) then
           !Log spacing
           diff=(log10(sin(theta)*scan(1)%range_f)-log10(kperp))/&
                real(scan(1)%n_scan*scan(1)%n_res)
           diff2=(log10(cos(theta)*scan(1)%range_f)-log10(kpar))/&
                real(scan(1)%n_scan*scan(1)%n_res)
        else
           !Linear spacing
           diff=((sin(theta)*scan(1)%range_f)-(kperp))/&
                real(scan(1)%n_scan*scan(1)%n_res)
           diff2=((cos(theta)*scan(1)%range_f)-(kpar))/&
                real(scan(1)%n_scan*scan(1)%n_res)
        endif

     endif
  
  elseif (scan(1)%style_s.eq.0) then
     !Global parameter scan
     if (scan(1)%type_s.eq.0) then
        write(param,'(a)')'kperp'
        sw=>kperp
     elseif (scan(1)%type_s.eq.1) then
        write(param,'(a)')'kpar'
        sw=>kpar
     elseif (scan(1)%type_s.eq.2) then
        write(param,'(a)')'betap'
        sw=>betap
     elseif (scan(1)%type_s.eq.3) then
        write(param,'(a)')'vtp'
        sw=>vtp
     endif

     if (writeOut) &          
          write(*,'(3a,es14.4,a,es14.4)') &
          'Scan over ',trim(param),' from ',&
          scan(1)%range_i,' to ',scan(1)%range_f

  else
     if (scan(1)%style_s.le.nspec) then
        !Species parameter scan
        if (scan(1)%type_s.eq.0) then
           write(param,'(a)')'tauS'
           sw=>spec(scan(1)%style_s)%tau_s
        elseif (scan(1)%type_s.eq.1) then
           write(param,'(a)')'muS'
           sw=>spec(scan(1)%style_s)%mu_s
        elseif (scan(1)%type_s.eq.2) then
           write(param,'(a)')'alphS'
           sw=>spec(scan(1)%style_s)%alph_s
        elseif (scan(1)%type_s.eq.3) then
           write(param,'(a)')'qs'
           sw=>spec(scan(1)%style_s)%Q_s
        elseif (scan(1)%type_s.eq.4) then
           write(param,'(a)')'ns'
           sw=>spec(scan(1)%style_s)%D_s
        elseif (scan(1)%type_s.eq.5) then
           write(param,'(a)')'Vs'
           sw=>spec(scan(1)%style_s)%vv_s
        endif
     else
        write(*,'(a,i0,a,i0,a)')&
             'Attempting to scan parameters for species ',scan(1)%style_s,&
             ' in a ',nspec,' species plasma.'
     endif

     if (writeOut) &          
          write(*,'(3a,i3,a,es14.4,a,es14.4)') &
          'Scan over ',trim(param),' for species ',&
          scan(1)%style_s,' from ',&
          scan(1)%range_i,' to ',scan(1)%range_f

  endif

  !Scans over n_s, q_s, and vv_s for species s require adjustment of
  !n_s',q_s', and vv_s' for species s' (or, potentially multiple species)
  !Such a variation is not general, and a user wishing to perform such 
  !a scan will have to choose how to vary the s' parameters for their
  !particular case.
  if ((scan(1)%style_s)==1) then
     if ((scan(1)%type_s)==3) then !qp/qs
        write(*,'(a)') &
             'Scanning over q_p/q_s... '
        write(*,'(2a)') &
             'such a scan requires a non general variation of parameters ',&
             'to insure quasineutrality and a current free plasma'
        write(*,'(2a)')&
             'The user should construct an appropriate module inside',&
             'disprels.f90%om_scan module'
     elseif ((scan(1)%type_s)==4) then !ns/np
        write(*,'(a)') &
             'Scanning over n_s/n_s... '
        write(*,'(2a)') &
             'such a scan requires a non general variation of parameters ',&
             'to insure quasineutrality and a current free plasma'
        write(*,'(2a)')&
             'The user should construct an appropriate module inside',&
             'disprels.f90%om_scan module'
     elseif ((scan(1)%type_s)==5) then !vdrift_s/vAp
        write(*,'(a)') &
             'Scanning over v__drift s/v_Ap... '
        write(*,'(2a)') &
             'such a scan requires a non general variation of parameters ',&
             'to insure quasineutrality and a current free plasma'
        write(*,'(2a)')&
             'The user should construct an appropriate module inside',&
             'disprels.f90%om_scan module'
     endif
  endif

  write(outName,'(2a)')&
       trim(param),'_'

end subroutine set_map_pointers

!-=-=-=-=-=-
subroutine get_out_name(outName,tensorName,fmt,fmt_tnsr,out_type,is,diff,diff2)
  !!Sets I/O strings, step sizes, and formats for [[om_scan(subroutine)]].
  use vars, only : scan,nspec,dataName,writeOut,outputName,pi
  use vars, only : spec,sw,sw2,betap,vtp,kperp,kpar,low_n
  !>>>GGH: 1/18/23
  use vars, only : new_low_n
  !<<<GGH: 1/18/23
  implicit none

  !Passed
  character(150) :: outName
  !! String for I/O identification.

  character(150) :: tensorName
  !! String for Susceptability output.

  character (100) :: fmt
  !! Output format for solution.
  
  character (100) :: fmt_tnsr
  !! Output format for susceptibility tensor.
  
  integer :: is
  !!Scan index

  integer :: out_type
  !! logic for outputing supplementary eigenfunction and heating
  !! calculations.
  !! 0: Frequency, Eigenfuction, Heating.
  !! 1: Frequency, Eigenfuction.
  !! 2: Frequency, Heating.
  !! 3: Frequency.
  
  real :: diff
  !!Spacing for scanned parameter.

  real :: diff2
  !!Suplemental spacing for scanned parameter.
  
  !Local
  character(15)  :: param
  !! String to identify scanned parameter.

  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)
  
  real :: theta_q              
  !! Scanned values of \( \theta \).

  real :: ki
  !! Initial value of \( |k| \).

  real :: kf
  !! Final value of \( |k| \).

  real :: kperpi
  !! Initial value of \( k_{\perp} \).

  real :: kpari
  !! Initial value of \( k_{\parallel} \).
  

  if((scan(is)%style_s).ge.0) then
     if (scan(is)%log_scan) then
        !Log spacing
        diff=(log10(scan(is)%range_f)-log10(scan(is)%range_i))/&
             real(scan(is)%n_scan*scan(is)%n_res)
     else
        !Linear spacing
        diff=((scan(is)%range_f)-(scan(is)%range_i))/&
             real(scan(is)%n_scan*scan(is)%n_res)
     endif
  endif


  if (scan(is)%style_s.eq.-1) then
     !Global parameter scan- 2 components
     if (scan(is)%type_s.eq.0) then
        !k_0 -> k_1; 
        !k_0= current k
        !kperp_1 = range_i, kpar_1 = range_f
        write(param,'(a)')'k'
        sw=>kperp
        sw2=>kpar
        kpari=kpar
        kperpi=kperp
        if (writeOut) &          
             write(*,'(4a,es14.4,a,es14.4,a,es14.4,a,es14.4,a)') &
             'Scan over ',trim(param),' from ',&
             '(kperp,kpar) = (',kperp,',',kpar,') to (',&
             scan(is)%range_i,',',scan(is)%range_f,')'

        write(outName,'(7a,i0,a,i0,a,i0,a,i0)')&
             'data/',trim(dataName),'/',&
             trim(outputName),'_',&
             trim(param),'_',int(1000.*kperp),&
             '_',int(1000.*kpar),&                
             '_',int(1000.*scan(is)%range_i),&
             '_',int(1000.*scan(is)%range_f)                

        if (scan(is)%tensor_s) then
           write(tensorName,'(7a,i0,a,i0,a,i0,a,i0)')&
                'data/',trim(dataName),'/',&
                trim(outputName),'_tensor_',&
                trim(param),'_',int(1000.*kperp),&
                '_',int(1000.*kpar),&                
                '_',int(1000.*scan(is)%range_i),&
                '_',int(1000.*scan(is)%range_f)                
        endif

        if (scan(is)%log_scan) then
           !Log spacing
           diff=(log10(scan(is)%range_i)-log10(kperp))/&
                real(scan(is)%n_scan*scan(is)%n_res)
           diff2=(log10(scan(is)%range_f)-log10(kpar))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        else
           !Linear spacing
           diff=((scan(is)%range_i)-(kperp))/&
                real(scan(is)%n_scan*scan(is)%n_res)
           diff2=((scan(is)%range_f)-(kpar))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        endif
        

     elseif (scan(is)%type_s.eq.1) then
        write(param,'(a)')'theta'
        !theta_0 -> theta_1
        !atan(kperp/kpar) -> theta_0
        !swi ->  theta_1 in degrees
        theta=atan(kperp/kpar)!initial theta value
        ki=kpar/cos(theta)
        sw=>kperp
        sw2=>kpar
        if (writeOut) &          
             write(*,'(4a,es14.4,a,es14.4,a,es14.4,a,es14.4,a)') &
             'Scan over ',trim(param),' from ',&
             'theta = ',theta*180./pi,' to ',&
             scan(is)%range_i

        write(outName,'(7a,i0,a,i0)')&
             'data/',trim(dataName),'/',&
             trim(outputName),'_',&
             trim(param),'_',int(10.*theta*180./pi),&
             '_',int(10.*scan(is)%range_i)                

        if (scan(is)%tensor_s) then
           write(tensorName,'(7a,i0,a,i0)')&
                'data/',trim(dataName),'/',&
                trim(outputName),'_tensor_',&
                trim(param),'_',int(10.*theta*180./pi),&
                '_',int(10.*scan(is)%range_i)                
        endif

        if (scan(is)%log_scan) then
           !Log spacing
           diff=(log10(pi*(scan(is)%range_i)/180.)-log10(theta))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        else
           !Linear spacing
           diff=((pi*(scan(is)%range_i)/180.)-(theta))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        endif

     elseif (scan(is)%type_s.eq.2) then
        write(param,'(a)')'k_fixed_theta'
        !k_0-> k_1 along fixed (current) angle
        !k_0 -> k_1; 
        !k_0= current k
        !|k_1|=swf
        theta = atan(kperp/kpar)
        ki=kpar/cos(theta)
        sw=>kperp
        sw2=>kpar

        if (writeOut) &          
             write(*,'(4a,es14.4,a,es14.4,a,es14.4)') &
             'Scan over ',trim(param),' from ',&
             '|k| = ',ki,' to ',scan(is)%range_f,&
             ' at fixed theta = ',theta * 180./pi

        write(outName,'(7a,i0,a,i0,a,i0)')&
             'data/',trim(dataName),'/',&
             trim(outputName),'_',&
             trim(param),'_',int(1000.*ki),&
             '_',int(1000.*scan(is)%range_f),'_theta',&
             int(1000 * theta * 180./pi)

        if (scan(is)%tensor_s) then
           write(tensorName,'(7a,i0,a,i0,a,i0)')&
                'data/',trim(dataName),'/',&
                trim(outputName),'_tensor_',&
                trim(param),'_',int(1000.*ki),&
                '_',int(1000.*scan(is)%range_f),'_theta',&
                int(1000 * theta * 180./pi)
        endif

        if (scan(is)%log_scan) then
           !Log spacing
           diff=(log10(sin(theta)*scan(is)%range_f)-log10(kperp))/&
                real(scan(is)%n_scan*scan(is)%n_res)
           diff2=(log10(cos(theta)*scan(is)%range_f)-log10(kpar))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        else
           !Linear spacing
           diff=((sin(theta)*scan(is)%range_f)-(kperp))/&
                real(scan(is)%n_scan*scan(is)%n_res)
           diff2=((cos(theta)*scan(is)%range_f)-(kpar))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        endif

     endif
  
  elseif (scan(is)%style_s.eq.0) then
     !Global parameter scan
     if (scan(is)%type_s.eq.0) then
        write(param,'(a)')'kperp'
        sw=>kperp
     elseif (scan(is)%type_s.eq.1) then
        write(param,'(a)')'kpar'
        sw=>kpar
     elseif (scan(is)%type_s.eq.2) then
        write(param,'(a)')'betap'
        sw=>betap
     elseif (scan(is)%type_s.eq.3) then
        write(param,'(a)')'vtp'
        sw=>vtp
     endif

     if (writeOut) &          
          write(*,'(3a,es14.4,a,es14.4)') &
          'Scan over ',trim(param),' from ',&
          scan(is)%range_i,' to ',scan(is)%range_f

     write(outName,'(7a,i0,a,i0)')&
          'data/',trim(dataName),'/',&
          trim(outputName),'_',&
          trim(param),'_',int(1000.*scan(is)%range_i),&
          '_',int(1000.*scan(is)%range_f)

        if (scan(is)%tensor_s) then
           write(tensorName,'(7a,i0,a,i0)')&
                'data/',trim(dataName),'/',&
                trim(outputName),'_tensor_',&
                trim(param),'_',int(1000.*scan(is)%range_i),&
                '_',int(1000.*scan(is)%range_f)
        endif

  else
     if (scan(is)%style_s.le.nspec) then
        !Species parameter scan
        if (scan(is)%type_s.eq.0) then
           write(param,'(a)')'tauS'
           sw=>spec(scan(is)%style_s)%tau_s
        elseif (scan(is)%type_s.eq.1) then
           write(param,'(a)')'muS'
           sw=>spec(scan(is)%style_s)%mu_s
        elseif (scan(is)%type_s.eq.2) then
           write(param,'(a)')'alphS'
           sw=>spec(scan(is)%style_s)%alph_s
        elseif (scan(is)%type_s.eq.3) then
           write(param,'(a)')'qs'
           sw=>spec(scan(is)%style_s)%Q_s
        elseif (scan(is)%type_s.eq.4) then
           write(param,'(a)')'ns'
           sw=>spec(scan(is)%style_s)%D_s
        elseif (scan(is)%type_s.eq.5) then
           write(param,'(a)')'Vs'
           sw=>spec(scan(is)%style_s)%vv_s
        endif
        write(outName,'(7a,i0,a,i0,a,i0)')&
             'data/',trim(dataName),'/',&
             trim(outputName),'_',&
             trim(param),'_s',scan(is)%style_s,'_',int(1000.*scan(is)%range_i),&
             '_',int(1000.*scan(is)%range_f)

        if (scan(is)%tensor_s) then
           write(tensorName,'(7a,i0,a,i0,a,i0)')&
                'data/',trim(dataName),'/',&
                trim(outputName),'_tensor_',&
                trim(param),'_s',scan(is)%style_s,'_',int(1000.*scan(is)%range_i),&
                '_',int(1000.*scan(is)%range_f)
        endif

     else
        write(*,'(a,i0,a,i0,a)')&
             'Attempting to scan parameters for species ',scan(is)%style_s,&
             ' in a ',nspec,' species plasma.'
     endif

     if (writeOut) &          
          write(*,'(3a,i3,a,es14.4,a,es14.4)') &
          'Scan over ',trim(param),' for species ',&
          scan(is)%style_s,' from ',&
          scan(is)%range_i,' to ',scan(is)%range_f

  endif

  !Format output structure
  !Global (kperp,kpar,beta,vtp)
  !frequency (om,gam)
  !eigen (Bx,By,Bz,
  !       Ex,Ey,Ez,
  !       Ux,Uy,Uz|spec
  !       n|spec) *optional
  !power (P|spec)
  !Species (tau, mu, alph, Q, n, v_drift)|s
  if (scan(is)%eigen_s) then
     if (scan(is)%heat_s) then
        if (low_n) then
           if (new_low_n) then  !Greg's New low_n calculation
              write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',21*nspec,'es15.6e3)'
           else  !Old version of low_n 
              !<<<GGH: 1/18/23
              write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',19*nspec,'es15.6e3)'
              !>>>GGH: 1/18/23
           endif
           !<<<GGH: 1/18/23
        else
           write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',15*nspec,'es15.6e3,es15.6e3)'
        endif
        out_type=0
     else
        write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',14*nspec,'es15.6e3)'
        out_type=1
     endif
  else
     if (scan(is)%heat_s) then
        if (low_n) then
           !>>>GGH: 1/18/23
           if (new_low_n) then  !Greg's New low_n calculation
              write(fmt,'(a,i0,a)')'(6es15.6e3,',13*nspec,'es15.6e3)'
           else  !Old version of low_n 
              !<<<GGH: 1/18/23
              write(fmt,'(a,i0,a)')'(6es15.6e3,',11*nspec,'es15.6e3)'
              !>>>GGH: 1/18/23
           endif
           !<<<GGH: 1/18/23
        else
           write(fmt,'(a,i0,a)')'(6es15.6e3,',7*nspec,'es15.6e3)'
        endif
        out_type=2
     else
        write(fmt,'(a,i0,a)')'(6es15.6e3,',6*nspec,'es15.6e3)'
        out_type=3
     endif
  endif
  if (scan(is)%tensor_s) &
       write(fmt_tnsr,'(a,i0,a)')'(4es15.6e3,',18*nspec,'es15.6e3)'

  !Scans over n_s, q_s, and vv_s for species s require adjustment of
  !n_s',q_s', and vv_s' for species s' (or, potentially multiple species)
  !Such a variation is not general, and a user wishing to perform such 
  !a scan will have to choose how to vary the s' parameters for their
  !particular case.
  if ((scan(is)%style_s)==1) then
     if ((scan(is)%type_s)==3) then !qp/qs
        write(*,'(a)') &
             'Scanning over q_p/q_s... '
        write(*,'(2a)') &
             'such a scan requires a non general variation of parameters ',&
             'to insure quasineutrality and a current free plasma'
        write(*,'(2a)')&
             'The user should construct an appropriate module inside',&
             'disprels.f90%om_scan module'
     elseif ((scan(is)%type_s)==4) then !ns/np
        write(*,'(a)') &
             'Scanning over n_s/n_s... '
        write(*,'(2a)') &
             'such a scan requires a non general variation of parameters ',&
             'to insure quasineutrality and a current free plasma'
        write(*,'(2a)')&
             'The user should construct an appropriate module inside',&
             'disprels.f90%om_scan module'
     elseif ((scan(is)%type_s)==5) then !vdrift_s/vAp
        write(*,'(a)') &
             'Scanning over v__drift s/v_Ap... '
        write(*,'(2a)') &
             'such a scan requires a non general variation of parameters ',&
             'to insure quasineutrality and a current free plasma'
        write(*,'(2a)')&
             'The user should construct an appropriate module inside',&
             'disprels.f90%om_scan module'
     endif
  endif

end subroutine get_out_name
!-=-=-=-

!-=-=-=-
subroutine get_double_out_name(outName,tensorName,fmt,fmt_tnsr,out_type,diff)
  !!Sets I/O strings, step sizes, and formats for [[om_double_scan(subroutine)]].
  use vars, only : scan,nspec,dataName,writeOut,outputName,pi
  use vars, only : spec,sw,sw2,sw3,sw4,betap,vtp,kperp,kpar, low_n
  !>>>GGH: 1/18/23
  use vars, only : new_low_n
  !<<<GGH: 1/18/23
  implicit none
  !Passed

  character(150) :: outName
  !! String for I/O identification.

  character(150) :: tensorName
  !! String for Susceptability output.

  character (100) :: fmt
  !! Output format for solution.
  
  character (100) :: fmt_tnsr
  !! Output format for susceptibility tensor.

  integer :: out_type
  !! logic for outputing supplementary eigenfunction and heating
  !! calculations.
  !! 0: Frequency, Eigenfuction, Heating.
  !! 1: Frequency, Eigenfuction.
  !! 2: Frequency, Heating.
  !! 3: Frequency.
  
  real,dimension(1:2,1:2) :: diff
  !!Spacing for scanned parameter.
  
  !Local
  character(15), dimension(1:2)  :: param
  !! String to identify scanned parameters.


  real :: theta
  !! \( \atan(k_{\perp}/k_{\parallel}) \)
  
  real :: theta_q              
  !! Scanned values of \( \theta \).

  real :: ki
  !! Initial value of \( |k| \).

  real :: kf
  !! Final value of \( |k| \).

  real :: kperpi
  !! Initial value of \( k_{\perp} \).

  real :: kpari
  !! Initial value of \( k_{\parallel} \).
  
  integer :: is
  !! Scan index.

  if (((scan(1)%style_s)==(scan(2)%style_s)).and.&
       ((scan(1)%type_s)==(scan(2)%type_s))) then
     write(*,'(a)')&
          'Scaning the same parameters in double scan...'
     write(*,'(a)')&
          'Change input file values...'
     stop
  endif

  do is=1,2!iterate over first and second scan information

     if((scan(is)%style_s).ge.0) then
        if (scan(is)%log_scan) then
           diff(is,1)=(log10(scan(is)%range_f)-log10(scan(is)%range_i))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        else
           diff(is,1)=((scan(is)%range_f)-(scan(is)%range_i))/&
                real(scan(is)%n_scan*scan(is)%n_res)
        endif
     endif

     if (scan(is)%style_s.eq.-1) then
        !Global parameter scan- 2 components
        if (scan(is)%type_s.eq.0) then
           !k_0 -> k_1; 
           !k_0= current k
           !kperp_1 = range_i, kpar_1 = range_f
           write(param(is),'(a)')'k'
           if (is==1) then
              sw=>kperp
              sw2=>kpar
           else
              sw3=>kperp
              sw4=>kpar
           endif
           kpari=kpar
           kperpi=kperp
           if (writeOut) &          
                write(*,'(a,i0,4a,es15.6e3,a,es15.6e3,a,es15.6e3,a,es15.6e3,a)') &
                'Scan ',is,' over ',trim(param(is)),' from ',&
                '(kperp,kpar) = (',kperp,',',kpar,') to (',&
                scan(is)%range_i,',',scan(is)%range_f,')'
           
           if (scan(is)%log_scan) then
              !Log spacing
              diff(is,1)=(log10(scan(is)%range_i)-log10(kperp))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              diff(is,2)=(log10(scan(is)%range_f)-log10(kpar))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
           else
              !Linear spacing
              diff(is,1)=((scan(is)%range_i)-(kperp))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              diff(is,2)=((scan(is)%range_f)-(kpar))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
           endif
                      
        elseif (scan(is)%type_s.eq.1) then
           write(param(is),'(a)')'theta'
           !theta_0 -> theta_1
           !atan(kperp/kpar) -> theta_0
           !swi ->  theta_1 in degrees
           theta=atan(kperp/kpar)!initial theta value
           ki=kpar/cos(theta)
           if (is==1) then
              sw=>kperp
              sw2=>kpar
           else
              sw3=>kperp
              sw4=>kpar
           endif
           if (writeOut) &          
                write(*,'(a,i0,4a,es15.6e3,a,es15.6e3,a,es15.6e3,a,es15.6e3,a)') &
                'Scan ',is,' over ',trim(param(is)),' from ',&
                'theta = ',theta*180./pi,' to ',&
                scan(is)%range_i
           
           if (scan(is)%log_scan) then
              !Log spacing
              diff(is,1)=(log10(pi*(scan(is)%range_i)/180.)-log10(theta))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
           else
              !Linear spacing
              diff(is,1)=((pi*(scan(is)%range_i)/180.)-(theta))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
           endif
           
        elseif (scan(is)%type_s.eq.2) then
           write(param(is),'(a)')'k_fixed_theta'
           !k_0-> k_1 along fixed (current) angle
           !k_0 -> k_1; 
           !k_0= current k
           !|k_1|=swf
           theta = atan(kperp/kpar)
           ki=kpar/cos(theta)
           if (is==1) then
              sw=>kperp
              sw2=>kpar
           else
              sw3=>kperp
              sw4=>kpar
           endif
           
           if (writeOut) &          
                write(*,'(a,i0,4a,es15.6e3,a,es15.6e3,a,es15.6e3)') &
                'Scan ',is,' over ',trim(param(is)),' from ',&
                '|k| = ',ki,' to ',scan(is)%range_f,&
                ' at fixed theta = ',theta * 180./pi
                      
           if (scan(is)%log_scan) then
              !Log spacing
              diff(is,1)=(log10(sin(theta)*scan(is)%range_f)-log10(kperp))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              diff(is,2)=(log10(cos(theta)*scan(is)%range_f)-log10(kpar))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
           else
              !Linear spacing
              diff(is,1)=((sin(theta)*scan(is)%range_f)-(kperp))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              diff(is,2)=((cos(theta)*scan(is)%range_f)-(kpar))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
           endif
           
        endif
        
     elseif (scan(is)%style_s.eq.0) then
        !Global parameter scan
        if (scan(is)%type_s.eq.0) then
           write(param(is),'(a)')'kperp'
           if (is==1) sw=>kperp
           if (is==2) sw3=>kperp
        elseif (scan(is)%type_s.eq.1) then
           write(param(is),'(a)')'kpar'
           if (is==1) sw=>kpar
           if (is==2) sw3=>kpar
        elseif (scan(is)%type_s.eq.2) then
           write(param(is),'(a)')'betap'
           if (is==1) sw=>betap
           if (is==2) sw3=>betap
        elseif (scan(is)%type_s.eq.3) then
           write(param(is),'(a)')'vtp'
           if (is==1) sw=>vtp
           if (is==2) sw3=>vtp
        endif
        
        if (writeOut) &          
             write(*,'(3a,es15.6e3,a,es15.6e3)') &
             'Scan over ',trim(param(is)),' from ',&
             scan(is)%range_i,' to ',scan(is)%range_f   
     else
        if (scan(is)%style_s.le.nspec) then
           !Species parameter scan
           if (scan(is)%type_s.eq.0) then
              write(param(is),'(a)')'tauS'
              if (is==1) sw=>spec(scan(is)%style_s)%tau_s
              if (is==2) sw3=>spec(scan(is)%style_s)%tau_s
           elseif (scan(is)%type_s.eq.1) then
              write(param(is),'(a)')'muS'
              if (is==1) sw=>spec(scan(is)%style_s)%mu_s
              if (is==2) sw3=>spec(scan(is)%style_s)%mu_s
           elseif (scan(is)%type_s.eq.2) then
              write(param(is),'(a)')'alphS'
              if (is==1) sw=>spec(scan(is)%style_s)%alph_s
              if (is==2) sw3=>spec(scan(is)%style_s)%alph_s
           elseif (scan(is)%type_s.eq.3) then
              write(param(is),'(a)')'qs'
              if (is==1) sw=>spec(scan(is)%style_s)%Q_s
              if (is==2) sw3=>spec(scan(is)%style_s)%Q_s
           elseif (scan(is)%type_s.eq.4) then
              write(param(is),'(a)')'ns'
              if (is==1) sw=>spec(scan(is)%style_s)%D_s
              if (is==2) sw3=>spec(scan(is)%style_s)%D_s
           elseif (scan(is)%type_s.eq.5) then
              write(param(is),'(a)')'Vs'
              if (is==1) sw=>spec(scan(is)%style_s)%vv_s
              if (is==2) sw3=>spec(scan(is)%style_s)%vv_s
           endif
        else
           write(*,'(a,i0,a,i0,a)')&
                'Attempting to scan parameters for species ',scan(is)%style_s,&
                ' in a ',nspec,' species plasma.'
        endif
        
        if (writeOut) &          
             write(*,'(3a,i3,a,es15.6e3,a,es15.6e3)') &
             'Scan over ',trim(param(is)),' for species ',&
             scan(is)%style_s,' from ',&
             scan(is)%range_i,' to ',scan(is)%range_f
        
     endif
          
     !Scans over n_s, q_s, and vv_s for species s require adjustment of
     !n_s',q_s', and vv_s' for species s' (or, potentially multiple species)
     !Such a variation is not general, and a user wishing to perform such 
     !a scan will have to choose how to vary the s' parameters for their
     !particular case.
     if ((scan(is)%style_s)==1) then
        if ((scan(is)%type_s)==3) then !qp/qs
           write(*,'(a)') &
                'Scanning over q_p/q_s... '
           write(*,'(2a)') &
                'such a scan requires a non general variation of parameters ',&
                'to insure quasineutrality and a current free plasma'
           write(*,'(2a)')&
                'The user should construct an appropriate module inside',&
                'disprels.f90%om_scan module'
        elseif ((scan(is)%type_s)==4) then !ns/np
           write(*,'(a)') &
                'Scanning over n_s/n_s... '
           write(*,'(2a)') &
                'such a scan requires a non general variation of parameters ',&
                'to insure quasineutrality and a current free plasma'
           write(*,'(2a)')&
                'The user should construct an appropriate module inside',&
                'disprels.f90%om_scan module'
        elseif ((scan(is)%type_s)==5) then !vdrift_s/vAp
           write(*,'(a)') &
                'Scanning over v__drift s/v_Ap... '
           write(*,'(2a)') &
                'such a scan requires a non general variation of parameters ',&
                'to insure quasineutrality and a current free plasma'
           write(*,'(2a)')&
                'The user should construct an appropriate module inside',&
                'disprels.f90%om_scan module'
        endif
     endif
     
  enddo!End (is) loop for first and second scan information

     !Format output structure
     !Global (kperp,kpar,beta,vtp)
     !frequency (om,gam)
     !eigen (Bx,By,Bz,
     !       Ex,Ey,Ez,
     !       Ux,Uy,Uz|spec
     !       n|spec) *optional
     !power (P|spec)
     !Species (tau, mu, alph, Q, n, v_drift)|s

     !Scan 2's choice for eigen, heating output will
     !specify the output for the parameter plane
  !ADD PS_SPLIT!
     if (scan(2)%eigen_s) then
        if (scan(2)%heat_s) then
           if (low_n) then
              !>>>GGH: 1/18/23
              if (new_low_n) then  !Greg's New low_n calculation
                 write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',21*nspec,'es15.6e3)'
              else  !Old version of low_n 
                 !<<<GGH: 1/18/23
                 write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',19*nspec,'es15.6e3)'
                 !>>>GGH: 1/18/23
              endif
              !<<<GGH: 1/18/23
           else
              write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',15*nspec,'es15.6e3,es15.6e3)'
           endif
           out_type=0
        else
           write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',14*nspec,'es15.6e3)'
           out_type=1
        endif
     else
        if (scan(2)%heat_s) then
           if (low_n) then
              !>>>GGH: 1/18/23
              if (new_low_n) then  !Greg's New low_n calculation
                 write(fmt,'(a,i0,a)')'(6es15.6e3,',13*nspec,'es15.6e3)'
              else  !Old version of low_n 
                 !<<<GGH: 1/18/23
                 write(fmt,'(a,i0,a)')'(6es15.6e3,',11*nspec,'es15.6e3)'
                 !>>>GGH: 1/18/23
              endif
              !<<<GGH: 1/18/23
           else
              write(fmt,'(a,i0,a)')'(6es15.6e3,',7*nspec,'es15.6e3,es15.6e3)'
           endif
           out_type=2
        else
           write(fmt,'(a,i0,a)')'(6es15.6e3,',6*nspec,'es15.6e3)'
           out_type=3
        endif
     endif

     if ((scan(1)%tensor_s).and.(scan(2)%tensor_s)) &
          write(fmt_tnsr,'(a,i0,a)')'(4es15.6e3,',18*nspec,'es15.6e3)'


     write(outName,'(8a)')&
          'data/',trim(dataName),'/',trim(outputName),'_',&
          trim(param(1)),'_',trim(param(2))

     if ((scan(1)%tensor_s).and.(scan(2)%tensor_s)) then
        write(tensorName,'(8a)')&
             'data/',trim(dataName),'/',trim(outputName),'_tensor_',&
             trim(param(1)),'_',trim(param(2))
     endif

end subroutine get_double_out_name
!-=-=-=-


!THE BEATING HEART OF THE CODE....
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!Collisionless Vlasov-Maxwell Dispersion Relation for an 
!Arbitrarily Large Set of ions and electrons, each with a
!Bi-Maxwellian Velocity Distribution.
!
!Based upon the Proton-Electron Maxwellian Vlasov-Maxwell Dispersion Solver 
!Written by Greg Howes and Elliot Quataert
!
!             Copyright 2015
!             Kristopher G Klein
!-=-=-=-=-=-=
!  The parameters upon which this dispersion relation depends are:
!   1) betap   = 8nT_parallel ref/B^2
!   2) kperp   = k_perp rho_ref
!   3) kpar    = k_parallel rho_ref (allows for +/- k_parallel values)
!   4) vtp     = v_t,ref,parallel/c
!  As well as the species dependent ratios, compared to the selected
!  Reference Species ref:
!   5) tau_s     = T_ref/T_s(both temperatures parallel)
!   6) mu_s      = m_ref/m_s
!   7) alph_s    = T_perp/T_par|_s
!   8) Q_s       = q_ref/q_s
!   9) D_s       = n_s/n_ref
!   10)vv_s      = v_drift s/v_Aref
!   and returns om=omega/Omega_ref
!!!!
!NOTE: REFERENCE SPECIES ASSUMED TO BE FIRST SPECIES
!!!!
   complex function disp(om)
     !!Returns Plasma Dispersion Function for input complex frequency
     !! and global values of the dimensionless plasma parameters.
     use vars,  only: betap,kperp,kpar,vtp,nspec,spec,susc,lam 
     use vars,  only: low_n, susc_low, pi
     !>>>GGH: 1/18/23
     use vars, only : new_low_n
     !<<<GGH: 1/18/23
     implicit none
     complex :: om
     !!Complex Frequency

     complex :: enx2
     !! \(n_x^2= \left(k_{x} c/\omega \right)^2 \)
     complex :: enz2
     !! \(n_z^2= \left(k_{z} c/\omega \right)^2 \)
     complex :: enxnz
     !! \(n_x n_z = k_{x} k_{z}c^2/\omega^2 \)

     complex :: eps_xx
     !!\(\epsilon_{xx} \) element of dielectric tensor

     complex :: eps_yy
     !!\(\epsilon_{yy} \) element of dielectric tensor
     
     complex :: eps_zz
     !!\(\epsilon_{zz} \) element of dielectric tensor

     complex :: eps_xy
     !!\(\epsilon_{xy} \) element of dielectric tensor

     complex :: eps_xz
     !!\(\epsilon_{xz} \) element of dielectric tensor

     complex :: eps_yz
     !!\(\epsilon_{yz} \) element of dielectric tensor

     complex :: eps_xx_t
     !!Temp. \(\epsilon_{xx} \) element of dielectric tensor

     complex :: eps_yy_t
     !!Temp. \(\epsilon_{yy} \) element of dielectric tensor
     
     complex :: eps_zz_t
     !!Temp. \(\epsilon_{zz} \) element of dielectric tensor

     complex :: eps_xy_t
     !!Temp. \(\epsilon_{xy} \) element of dielectric tensor

     complex :: eps_xz_t
     !!Temp. \(\epsilon_{xz} \) element of dielectric tensor

     complex :: eps_yz_t
     !!Temp. \(\epsilon_{yz} \) element of dielectric tensor
     
     real :: lambdap
     !! \( \lambda_{ref}=(k_{\perp} \rho_{ref})^2/2 \)

     real :: lambdas
     !! \( \lambda_{s}=(k_{\perp} \rho_{s})^2/2 \)

     real :: alphp
     !! \( T_{\perp}/T_{\parallel}_{ref} \)

     real :: Vdrifts
     !! \( k_\parallel \rho_{ref} \Delta v_s /
     !! \sqrt{\beta_{parallel,ref} T_{\perp}/T_{\parallel}_{ref}} \)

     !Parameters
     complex, parameter :: c_i =(0.0,1.0)
     !!Imaginary unit i

     integer :: n
     !!Counter for Bessel sums
     
     real, allocatable, dimension(:,:), save :: jn
     !!Modified Bessel functions
     
     real, allocatable, dimension(:,:), save :: jpn
     !!Derivative of Modified Bessel functions
     
     complex, dimension(-nbrack:nbrack), save :: tsi
     !!Argument for Plasma Dispersion Function Z
     
     complex, dimension(-nbrack:nbrack), save :: zz
     !!Value fo Plasma Dispersion Function Z
     
     complex, dimension(1:6) :: norm
     !!normalization for \( \chi_s \)
     !!1: xx, 2: yy, 3: zz, 4: xy, 5: xz, 6: yz
     
     complex, dimension(3,3) :: a
     !!Temp. Matrix for readability
     
     integer :: is
     !! Species/component index.
     
     complex ::temp
     !! Temp. Scalar for readability.

     !Species Parameters called locally
     real :: disp_tau     
     !!Relative Temperature ratio.
     !!\(T_{ref}/T_{s}|_{\parallel}\)
     
     real :: disp_mu
     !!Relative Mass ratio.
     !!\(m_{ref}/m_{s}\)

     real :: disp_alph
     !!Temperature Anisotropy.
     !!\(T_{\perp}/T_{\parallel}_s\)

     real :: disp_Q
     !!Relative charge ratio.
     !!\(q_{ref}/q_{s}\)

     real :: disp_D
     !!Density Ratio.
     !!\(n_{s}/n_{ref}\)

     real :: disp_v
     !!Relative Drift, normalized to reference Alfven velocity
     !!\(v_{drift}/v_{A,ref}\)
     !! with \(v_{A,ref} = B/\sqrt{4 \pi n_{ref} m_{ref}}\).
     
     !For computational efficiency, save Bessel functions between calls==========
     real, save :: kperp_last
     !! Previous value of \(k_{\perp} \rho_{ref} \)
     
     real, allocatable, dimension(:), save :: tau_last
     !!Previous \(T_{ref}/T_{s}|_{\parallel}\) value.

     real, allocatable, dimension(:), save :: alph_last
     !!Previous

     real, allocatable, dimension(:), save :: mu_last
     !!Previous

     real, allocatable, dimension(:), save :: Q_last
     !!Previous
     
     logical, allocatable, dimension(:) :: reuse_bessel
     !! Either use previous Bessel array, or recalculate,
     !! based on if argument for a particular species/component changes.

     if (.not. allocated(jn)) allocate(jn(-nbrack-1:nbrack+1,nspec))          
     if (.not. allocated(jpn)) allocate(jpn(-nbrack:nbrack,nspec))          
     if (.not. allocated(tau_last)) allocate(tau_last(nspec))          
     if (.not. allocated(alph_last)) allocate(alph_last(nspec))          
     if (.not. allocated(mu_last)) allocate(mu_last(nspec))          
     if (.not. allocated(Q_last)) allocate(Q_last(nspec))          
     if (.not. allocated(reuse_bessel)) allocate(reuse_bessel(nspec))

     reuse_bessel=.false.
     !IF
     !lambdas=lambdap*(disp_Q**2. * disp_alph)/(disp_mu * disp_tau * alphp)
     !changes, the bessel functions need to be recalculated...
     if ((kperp .eq. kperp_last ).and.(spec(1)%alph_s .eq. alph_last(1))) then
        do is=1, nspec 
           if ((spec(is)%tau_s .eq. tau_last(is)) .and. &
               (spec(is)%alph_s .eq. alph_last(is)) .and. &
               (spec(is)%mu_s .eq. mu_last(is)) .and. &
               (spec(is)%q_s .eq. q_last(is))) &                
                reuse_bessel(is)=.true.
        enddo
     endif

     !Arguments of Bessel functions for reference species
     lambdap = kperp**2./2.
     
     !temperature anisotropy for reference species
     alphp = spec(1)%alph_s

     !Indices of refraction from the dispersion relation 
     !the dispersion relation  is multiplied by [omega/Omega_R]^2
     !to remove singularities in the function.
     enx2=((kperp/(om*vtp))**2.)/alphp !n_x^2
     enz2=((kpar/(om*vtp))**2.)/alphp  !n_z^2
     enxnz=(kperp*kpar/(om*vtp)**2.)/alphp !n_x n_z

     !Initialize eps factors
     eps_xx = cmplx(0.,0.) 
     eps_yy = cmplx(0.,0.) 
     eps_zz = cmplx(0.,0.) 
     eps_xy = cmplx(0.,0.) 
     eps_xz = cmplx(0.,0.) 
     eps_yz = cmplx(0.,0.) 

     !Calculate the Susceptibility Tensor for Each Species:
     do is=1,nspec
        disp_tau  =spec(is)%tau_s
        disp_mu   =spec(is)%mu_s
        disp_alph =spec(is)%alph_s
        disp_Q    =spec(is)%Q_s
        disp_D    =spec(is)%D_s
        disp_v    =spec(is)%vv_s
        
        !frequently employed v drift normalization
        Vdrifts = kpar * disp_v / sqrt(betap*alphp)

        lambdas=lambdap*(disp_Q**2. * disp_alph)/(disp_mu * disp_tau * alphp)
        
        !Clear the Temporary Susceptability
        eps_xx_t = cmplx(0.,0.) 
        eps_xy_t = cmplx(0.,0.) 
        eps_xz_t = cmplx(0.,0.) 
        eps_yy_t = cmplx(0.,0.) 
        eps_yz_t = cmplx(0.,0.) 
        eps_zz_t = cmplx(0.,0.) 

     !Compute all necessary plasma dispersion functions and bessel functions
        tsi=0. ;  zz=0. ;
        if (.not. reuse_bessel(is)) then 
           do n = -nbrack-1, nbrack+1
              jn(n,is)=bessel(abs(n),lambdas)
           enddo
           do n = -nbrack, nbrack
              jpn(n,is)=0.5*(jn(n+1,is)+jn(n-1,is))
           enddo
        endif

        do n = -nbrack, nbrack
           tsi(n) = sqrt((alphp * disp_tau)/(disp_mu))*&
                ( om - Vdrifts - dble(n)* disp_mu/disp_Q)/kpar
           zz(n)=zet_in(kpar,tsi(n))
        enddo

        !Set normalizations for susceptibilites
        !1- xx, 2- yy, 3- zz, 4- xy, 5- xz, 6- yz
        !-=-=-=-=-=-=-=-=-=-=
        !XX
        norm(1)=( betap * disp_mu * disp_D  )/&
             ( om *om * vtp**2. * disp_Q**2. * lambdas  )

        !-=-=-=-=-=-=-=-=-=-=
        !XY
        norm(4)=((-c_i * betap * disp_mu * disp_D  )/&
             ( om * om* vtp**2. * disp_Q**2. * kpar ))* &
             sqrt(alphp * disp_tau / disp_mu)

        !-=-=-=-=-=-=-=-=-=-=
        !XZ
        norm(5)=( betap * disp_D * kperp )/&
             ( om * om *vtp**2. * disp_Q * lambdas * kpar  )

        !-=-=-=-=-=-=-=-=-=-=
        !YY
        norm(2)=( betap * disp_mu * disp_D  )/&
             ( om * om* vtp**2. * disp_Q**2. * lambdas  )

        !-=-=-=-=-=-=-=-=-=-=
        !YZ
        norm(6)=( c_i * betap * disp_D * kperp )/&
             ( om * om * vtp**2. * disp_Q * kpar  )

        !-=-=-=-=-=-=-=-=-=-=
        !ZZ
        norm(3)=(2. * betap * alphp* disp_tau  * disp_D)/&
             ( om * om* vtp**2. * disp_Q**2. * kpar * kpar * disp_alph)

        !Sum over the Susceptibilities
        do n = 1, nbrack
           !-=-=-=-=-=-=-=-=-=-=
           !XX
           !-=-=-=-=-=-=-=-=-=-=
           eps_xx_t  = eps_xx_t  + &
                n*n*jn(n,is)*( 2.*(disp_alph-1.) &
                + (sqrt(alphp * disp_tau/ disp_mu )/ kpar )*(&
                disp_alph* (om - Vdrifts) * (zz(n)+zz(-n)) &
                + (zz(n)-zz(-n))*(1.-disp_alph)*(n*disp_mu/(disp_Q ))))

           !-=-=-=-=-=-=-=-=-=-=
           !YY
           !-=-=-=-=-=-=-=-=-=-=
           eps_yy_t = eps_yy_t + &
                (n*n*jn(n,is) + 2.*lambdas*lambdas*(jn(n,is)-jpn(n,is)))*&
                (2.*(disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                disp_alph*(om - Vdrifts)*(zz(n)+zz(-n)) + (n*disp_mu/(disp_Q))*&
                (1.-disp_alph)*(zz(n)-zz(-n))))
           
           !-=-=-=-=-=-=-=-=-=-=
           !ZZ
           !-=-=-=-=-=-=-=-=-=-=
           eps_zz_t = eps_zz_t + &
                jn(n,is)*(2.*om* (om*disp_alph - Vdrifts) - &
                2.*n*n*disp_mu**2.* (1.-disp_alph)/(disp_Q**2.) &
                +(sqrt(alphp*disp_tau/disp_mu)/kpar) * (&
                (zz(n)+zz(-n))*(om*om*disp_alph * (om - Vdrifts) &
                + n*n*disp_mu**2.*(om*(3.*disp_alph-2.) -disp_alph*Vdrifts)/disp_Q**2.) &
                +(n * disp_mu/disp_Q )*(zz(n)-zz(-n))*& !changed sign
                (om*(om+2.*disp_alph *Vdrifts -3.*disp_alph*om)+&
                n*n*disp_mu**2.*(1.-disp_alph)/( disp_Q**2.)) ) )!?
           
           
           !-=-=-=-=-=-=-=-=-=-=
           !XY
           !-=-=-=-=-=-=-=-=-=-=
           eps_xy_t = eps_xy_t + &
                n*(jn(n,is)-jpn(n,is))*( &
                (zz(n)-zz(-n))*disp_alph*(om-Vdrifts) &
                +(zz(n)+zz(-n))*(1.-disp_alph)*n*disp_mu/(disp_Q ))

           !-=-=-=-=-=-=-=-=-=-=
           !XZ
           !-=-=-=-=-=-=-=-=-=-=
           eps_xz_t = eps_xz_t + &
                n*jn(n,is)*(2.*n*(1.-disp_alph)*disp_mu/(disp_Q) + &
                (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (zz(n) -zz(-n))*(disp_alph*om*(om-Vdrifts)-&
                (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) + &
                (zz(n) +zz(-n))*(disp_mu*n*(om-2.*om*disp_alph+Vdrifts*disp_alph)/disp_Q) ))

           !-=-=-=-=-=-=-=-=-=-=
           !YZ
           !-=-=-=-=-=-=-=-=-=-=
           eps_yz_t = eps_yz_t + &
                (jn(n,is)-jpn(n,is))*(2.*(om*disp_alph-Vdrifts) &
                +(sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (zz(n)+zz(-n))*(disp_alph*om*(om-Vdrifts) -&
                (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) &
                +(zz(n)-zz(-n))*(n*disp_mu*(om-disp_alph*(2.*om-Vdrifts))/disp_Q) ))
           
           !-=-=-=-=-=-=-=-=-=-=
           !n=\pm 1 susceptibility
           if ((low_n).and.(n==1)) then
              susc_low(is,1,1,1) = norm(1)*&
                   n*n*jn(n,is)*( 2.*(disp_alph-1.) &
                   + (sqrt(alphp * disp_tau/ disp_mu )/ kpar )*(&
                   disp_alph* (om - Vdrifts) * (zz(n)+zz(-n)) &
                   + (zz(n)-zz(-n))*(1.-disp_alph)*(n*disp_mu/(disp_Q ))))

              susc_low(is,1,2,1) = norm(4)*&
                   n*(jn(n,is)-jpn(n,is))*( &
                   (zz(n)-zz(-n))*disp_alph*(om-Vdrifts) &
                   +(zz(n)+zz(-n))*(1.-disp_alph)*n*disp_mu/(disp_Q ))

              susc_low(is,2,1,1) = -susc_low(is,1,2,1)

              susc_low(is,1,3,1) = norm(5)*&
                   n*jn(n,is)*(2.*n*(1.-disp_alph)*disp_mu/(disp_Q) + &
                   (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                   (zz(n) -zz(-n))*(disp_alph*om*(om-Vdrifts)-&
                   (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) + &
                   (zz(n) +zz(-n))*(disp_mu*n*(om-2.*om*disp_alph+Vdrifts*disp_alph)/disp_Q) ))

              susc_low(is,3,1,1) = susc_low(is,1,3,1)

              susc_low(is,2,2,1) = norm(2)*&
                   (n*n*jn(n,is) + 2.*lambdas*lambdas*(jn(n,is)-jpn(n,is)))*&
                   (2.*(disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                   disp_alph*(om - Vdrifts)*(zz(n)+zz(-n)) + (n*disp_mu/(disp_Q))*&
                   (1.-disp_alph)*(zz(n)-zz(-n))))

              susc_low(is,2,3,1) = norm(6)*&
                   (jn(n,is)-jpn(n,is))*(2.*(om*disp_alph-Vdrifts) &
                   +(sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                   (zz(n)+zz(-n))*(disp_alph*om*(om-Vdrifts) -&
                   (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) &
                   +(zz(n)-zz(-n))*(n*disp_mu*(om-disp_alph*(2.*om-Vdrifts))/disp_Q) ))

              susc_low(is,3,2,1) = -susc_low(is,2,3,1)  

              susc_low(is,3,3,1) = norm(3)*&
                   jn(n,is)*(2.*om* (om*disp_alph - Vdrifts) - &
                   2.*n*n*disp_mu**2.* (1.-disp_alph)/(disp_Q**2.) &
                   +(sqrt(alphp*disp_tau/disp_mu)/kpar) * (&
                   (zz(n)+zz(-n))*(om*om*disp_alph * (om - Vdrifts) &
                   + n*n*disp_mu**2.*(om*(3.*disp_alph-2.) -disp_alph*Vdrifts)/disp_Q**2.) &
                   +(n * disp_mu/disp_Q )*(zz(n)-zz(-n))*& !changed sign
                   (om*(om+2.*disp_alph *Vdrifts -3.*disp_alph*om)+&
                   n*n*disp_mu**2.*(1.-disp_alph)/( disp_Q**2.)) ) )!?
           endif

        enddo

        !Add in the n=zero term

        !No n=0 XX term

        !No n=0 XY term

        !No n=0 XZ term

        !YY term
           eps_yy_t = eps_yy_t + &
                (2.*lambdas*lambdas*(jn(0,is)-jpn(0,is)))*&
                ((disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*&
                (disp_alph*(om-Vdrifts))*(zz(0)))

        !YZ term
           eps_yz_t = eps_yz_t + &
                (jn(0,is)-jpn(0,is))*(om*disp_alph-Vdrifts+&
                disp_alph*(om-Vdrifts)*om*(sqrt(alphp*disp_tau/disp_mu)/kpar)*zz(0))

        !ZZ term
           eps_zz_t = eps_zz_t + &
                jn(0,is)*om*(&
                om*disp_alph-Vdrifts +&
                zz(0) * om * (sqrt(alphp*disp_tau/disp_mu)/kpar) &
                *disp_alph*(om-Vdrifts))

           if (low_n) then
              susc_low(is,2,2,0) =  (2.*lambdas*lambdas*(jn(0,is)-jpn(0,is)))*&
                ((disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*&
                (disp_alph*(om-Vdrifts))*(zz(0))) * norm(2)
              susc_low(is,2,3,0) = (jn(0,is)-jpn(0,is))*(om*disp_alph-Vdrifts+&
                disp_alph*(om-Vdrifts)*om*&
                (sqrt(alphp*disp_tau/disp_mu)/kpar)*zz(0)) * norm(6)
              susc_low(is,3,2,0) = -susc_low(is,2,3,0)
              susc_low(is,3,3,0) = jn(0,is)*om*(&
                om*disp_alph-Vdrifts +&
                zz(0) * om * (sqrt(alphp*disp_tau/disp_mu)/kpar) &
                *disp_alph*(om-Vdrifts))* norm(3)
           endif
           
        !Apply the correct species normalization
           !1- xx, 2- yy, 3- zz, 4- xy, 5- xz, 6- yz
        !-=-=-=-=-=-=-=-=-=-=
        !XX
        eps_xx_t=eps_xx_t*norm(1)
        !XY
        eps_xy_t=eps_xy_t*norm(4)
        !XZ
        eps_xz_t=eps_xz_t*norm(5)
        !YY
        eps_yy_t=eps_yy_t*norm(2)
        !YZ
        eps_yz_t=eps_yz_t*norm(6)
        !ZZ
        eps_zz_t=eps_zz_t*norm(3)

        !Add the drift term
        eps_zz_t = eps_zz_t + &
             ((2.*disp_v)/(kpar*vtp**2.)) * &
             (disp_D*disp_tau*sqrt(betap*alphp)/(om*disp_alph*disp_Q**2.))

        !Save the temporary susceptibility to the 'susc' array
        susc(is,1,1) = eps_xx_t
        susc(is,1,2) = eps_xy_t; susc(is,2,1) = -eps_xy_t        
        susc(is,1,3) = eps_xz_t; susc(is,3,1) = eps_xz_t        
        susc(is,2,2) = eps_yy_t
        susc(is,2,3) = eps_yz_t; susc(is,3,2) = -eps_yz_t	       
        susc(is,3,3) = eps_zz_t

        !-=-=-=-=-=-=-=-=-=-=
        !Add to the total susceptibility
        eps_xx = eps_xx + eps_xx_t
        eps_xy = eps_xy + eps_xy_t
        eps_xz = eps_xz + eps_xz_t
        eps_yy = eps_yy + eps_yy_t
        eps_yz = eps_yz + eps_yz_t
        eps_zz = eps_zz + eps_zz_t

        !-=-=-=-=-=-=-=-=-=-=
        !Save values of lambda_s components to save on Bessel function calculations
        tau_last(is) =disp_tau
        alph_last(is)=disp_alph
        Q_last(is) =disp_Q
        mu_last(is)=disp_mu
     enddo
     
     !------------------------------------------------------------------------------
     !------------------------------------------------------------------------------

     !Assign Values of matrix elements in final wave equation
     !NORM!
     a(1,1) = 1. + eps_xx - enz2
     a(1,2) = eps_xy
     a(1,3) = eps_xz + enxnz
     a(2,1) = -a(1,2)
     a(2,2) = 1. + eps_yy - enx2 - enz2 
     a(2,3) = eps_yz
     a(3,1) = a(1,3)
     a(3,2) = -a(2,3)
     a(3,3) = 1. + eps_zz - enx2

     !Set matrix elements into public variable
     lam(:,:)=a(:,:)

     !Calculate the dispersion relation:
     !------------------------------------------------------------------------------
     disp= a(1,1)*( a(2,2)*a(3,3) + a(2,3)**2. ) + &
          2.*a(1,2)*a(2,3)*a(1,3) - a(1,3)**2.*a(2,2) +a(1,2)**2.*a(3,3)

     !------------------------------------------------------------------------------
     !Save values of lambda_s components to save on Bessel function calculations
     kperp_last=kperp

     return

   end function disp

   !------------------------------------------------------------------------------
   !                           Greg Howes, 2005
   !------------------------------------------------------------------------------
   !     NOTE: This routine was adapted from f77 routine by Eliot Quataert
   complex function rtsec(func,x1,x2,xacc,iflag)
     !!Secant routine to find minimum value of dispersion relation function.
     integer, parameter :: maxit=75
     !!Maximum number of iterations.
     
     complex :: func
     !! Function to evaluate.
     
     complex :: x1
     !!Lower bracket value of solution refinement.

     complex :: x2
     !!Upper bracket value of solution refinement.

     complex :: xl
     !!Holding parameter for refined solution.

     complex :: fl
     !!Function evaluated at x1.
     
     complex :: f
     !!Function evaluated at x2.

     complex :: swap
     !!Function swap value for iteration.

     complex :: dx
     !!Step size for further iteration.

     real    :: xacc
     !! Solution tolerance.

     integer :: iflag
     !! flag for tracking iteration.

     integer :: j
     !! Iteration index.

     complex :: maxr
     !!Maximum value for refined solution.

     complex :: minr
     !!Minimum value for refined solution.
     
     logical :: limits=.false.
     !!Enforce a limited range for refined solution.
     
     !Set limits for solution
     if (limits) then
        maxr=x2*10.
        minr=x1*0.1
     endif
     
     fl=func(x1)
     f=func(x2)
          
     if(abs(fl).lt.abs(f))then
        rtsec=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
     else
        xl=x1
        rtsec=x2
     endif
     do  j=1,maxit
        iflag = j
        if (abs(f-fl) .gt. 1.0E-37) then
           dx=(xl-rtsec)*f/(f-fl)
	else
           dx = (x2-x1)/25.0
        end if

        xl=rtsec
        fl=f
        rtsec=rtsec+dx/2.

        !Enforce limits
        if (limits) then
           !Real limit
           if (real(rtsec) .gt. real(maxr) .or. real(rtsec) .lt. real(minr)) &
              rtsec = cmplx(sqrt(real(minr)*real(maxr)), Aimag(rtsec))
           !NOTE: aimag(rtsec) should be negative!
           if (aimag(rtsec) .lt. aimag(maxr) .or. aimag(rtsec) .gt. aimag(minr)) &
              rtsec = cmplx(real(rtsec), -sqrt(aimag(minr)*aimag(maxr)))
        endif

        f=func(rtsec)

        if(abs(dx).lt.xacc.or. Abs(f) .eq. 0.)return
     enddo
     return
   end function rtsec
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!                           Greg Howes, 2006
!------------------------------------------------------------------------------
   subroutine find_minima(val,numroots,iroots,nroots)
     !! Finding minimum values on complex frequency grid.
     use vars, only : ni,nr
     implicit none
     !Passed
     real, dimension(:,:), pointer :: val
     !!Value of Dispersion relation.
     
     integer :: numroots
     !!Maximum number of roots.
     
     integer, dimension(1:2,1:numroots) :: iroots
     !!Indices of roots.
     
     integer, intent(out) :: nroots
     !!Number of roots found.
     
     !Local
     integer :: ir
     !! Real Frequency index.
     
     integer :: ii
     !! Imaginary Frequency index.
     
     !Find local minima in map
     iroots=0
     nroots=0
     do ii=ni,0,-1
        do ir=0,nr
           select case(ir)
           case(0)
              if (val(ir,ii) .lt. val(ir+1,ii)) then
                 select case(ii)
                 case(0)
                    if (val(ir,ii) .lt. val(ir,ii+1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 case(1:ni-1)
                    if (val(ir,ii) .lt. val(ir,ii-1) .and.  &
                         val(ir,ii) .lt. val(ir,ii+1))then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 case(ni)
                    if (val(ir,ii) .lt. val(ir,ii-1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 end select
              endif
           case(1:nr-1)
              if (val(ir,ii) .lt. val(ir-1,ii) .and.  &
                   val(ir,ii) .lt. val(ir+1,ii))then
                 select case(ii)
                 case(0)
                    if (val(ir,ii) .lt. val(ir,ii+1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 case(1:ni-1)
                    if (val(ir,ii) .lt. val(ir,ii-1) .and.  &
                         val(ir,ii) .lt. val(ir,ii+1))then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 case(ni)
                    if (val(ir,ii) .lt. val(ir,ii-1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 end select

              endif
           case(nr)
              if (val(ir,ii) .lt. val(ir-1,ii)) then
                select case(ii)
                 case(0)
                    if (val(ir,ii) .lt. val(ir,ii+1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 case(1:ni-1)
                    if (val(ir,ii) .lt. val(ir,ii-1) .and.  &
                         val(ir,ii) .lt. val(ir,ii+1))then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 case(ni)
                    if (val(ir,ii) .lt. val(ir,ii-1)) then
                       nroots=nroots+1
                       iroots(1,nroots)=ir
                       iroots(2,nroots)=ii
                    endif
                 end select
               endif
           end select
        enddo
     enddo

   end subroutine find_minima

!-------------------------------------------------------------------------------
   function zet_in(kpar,zin)
     !!Evaluate Plasma Dispersion Function (Z) for kparallel >/< 0
     !! Calls [[zetout(function)]] with appropriate argument.
     implicit none
     !passed
     real :: kpar
     !! \( k_{\parallel} \rho_{ref} \)
     
     complex :: zin
     !! Argument of plasma dispersion function.

     complex :: zet_in
     !! Value of plasma dispersion function.
     
     if (kpar.gt.0.) zet_in=zetout(zin)
     if (kpar.lt.0.) zet_in=-zetout(-zin)
     return

   end function zet_in

!------------------------------------------------------------------------------
!Greg Howes, 2004
   complex function zetout(zin)
     !! This is the subroutine used by Linsker to evaluate his Z-function.
     complex :: z,dzetaz,term,fmult,terme,an1,bn1,zsquar
     complex :: hold,temp1,temp2,ddzeta,dddzet,zin,zetaoz
     real :: imagte,imagmu,imagse,imagsu
     real :: error,x,y,fn,realte,realmu,realsu,realse
     logical :: fail
     integer :: n

     error=1.e-7
     z=zin
     zsquar=z*z
     x=real(z)
     y=aimag(z)
     fn=real(zsquar)
     if (y.gt.0.) go to 99
     if (abs(fn) < 174. .and. abs(aimag(zsquar)) < 5.e4) go to 98
     if (fn.gt.0.) go to 97
11    format (' argument wp of subroutine zetvec has too large a negative imaginary part, wp = '/2e14.7)
! GGH: Modification to return with a fail=.true.
     fail=.true.
     return

97   hold=(0.,0.)
     go to 99
98   hold=(0.,1.77245385090551603)*cexp(-zsquar)
99   if (x*x+y*y > 16.) go to 200
     if (abs(y) >= 1.) go to 300
     realte=-2.*x
     imagte=-2.*y
     realmu=.5*(imagte*imagte-realte*realte)
     imagmu=-imagte*realte
     realsu=realte
     imagsu=imagte
     if (x == 0. .and. y == 0.) go to 103
     fn=3.
100  realse=realte
     imagse=imagte
     realte=(realse*realmu-imagse*imagmu)/fn
     imagte=(realse*imagmu+imagse*realmu)/fn
     realse=realsu
     imagse=imagsu
     realsu=realsu+realte
     imagsu=imagsu+imagte
     fn=fn+2.
     if (abs(realse-realsu) > error .or. abs(imagse-imagsu) > error) go to 100
103  x=realsu
     fn=imagsu
     if (y > 0.) hold=(0.,1.77245385090551603)*cexp(-zsquar)
     zetaoz=cmplx(x,fn)+hold
     go to 401
200  fn=5.
     dddzet=6.
     term=dddzet
     fmult=.5/zsquar
201  terme=term
     term=term*fmult*fn*(fn-1.)/(fn-3.)
     zetaoz=term/terme
     if (abs(real(zetaoz))+abs(aimag(zetaoz)) > 1.) go to 250
     zetaoz=dddzet
     dddzet=dddzet+term
     fn=fn+2.
     if (cabs(zetaoz-dddzet) > error) go to 201
250  dddzet=dddzet/(zsquar*zsquar)
     if (y > 0.) go to 260
     fn=1.
     if (y < 0.) fn=2.
     dddzet=dddzet-4.*fn*hold*z*(2.*zsquar-3.)
260  ddzeta=-(4.+(zsquar-.5)*dddzet)/(z*(2.*zsquar-3.))
     dzetaz=(2.-z*ddzeta)/(2.*zsquar-1.)
     zetaoz=-(1.+.5*dzetaz)/z
     go to 401
300  if (y < 0.) z=conjg(z)
     terme=(1.,0.)
     term=(0.,0.)
     dzetaz=term
     fmult=terme
     n=0
     an1=z
     bn1=-z*z+.5
301  temp1=bn1*term+an1*terme
     temp2=bn1*fmult+an1*dzetaz
     zetaoz=temp1/temp2
     dzetaz=(zetaoz-term/fmult)/zetaoz
     if (abs(real(dzetaz)) < error .and. abs(aimag(dzetaz)) < error) go to 302
     bn1=bn1+2.
     n=n+1
     an1=-.5*float(n*(n+n-1))
     terme=term
     dzetaz=fmult
     term=temp1
     fmult=temp2
     if (n < 30) go to 301
302  if (y >= 0.) go to 401
     zetaoz=conjg(zetaoz)+2.*hold
401  zetout=zetaoz
9999 continue

   end function zetout

!------------------------------------------------------------------------------
!                           Greg Howes, 2006
!------------------------------------------------------------------------------
   real function bessel(n,x)
     !!  Call the correct Numerical Recipes Bessel function
     !!     bessels is a function which calls the numerical recipes
     !!     routine of the appropriate order -- this is to avoid
     !!     to many if - then statements in the program.
     !! Note that we are determining I_n(x) e^(-x) to avoid large argument errors.
     use bessels, only: bessim0,bessim1,bessim
     implicit none
     
     real :: x
     !!Bessel function argument.
     integer :: n
     !!Bessel function order.

     select case(n)
     case(0) 
        bessel=bessim0(x)
     case(1)
        bessel = bessim1(x)
     case(2:)
        bessel = bessim(n,x)
     case default
        write(*,'(a,i4,a)') 'ERROR: Bessel function of order ',n,' unexpected'
        stop
     end select

     return 
   end function bessel

end module disprels
