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
  use plasma_z_backend, only: zetout
  !! Dispersion relation function and associated eigenfunction and minimization
  !! calculations.
  implicit none
  private

  include 'kinetic_harmonics.inc'
  !!Number of Bessel functions to sum

  real, parameter :: tol=1.0E-13
  !!Root Search Tolerance 

  real, parameter :: prec=1.E-7
  !!Root Finding precision  

  private :: find_minima,zetout,get_out_name
  public :: calc_eigen, disp, zet_in, rtsec, bessel !! Required for fpc.f90
                                   !! (Consider reducing scope with submodule or
                                   !! other approach- fine for now)

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
      use vars, only : normalization, has_cold_species, dispersion_status,z_backend
      use cold_susceptibility, only : finite_complex, invalid_complex
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

      integer :: invalid_points = 0
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
         if (normalization == 'd_ref') then
            write(*,'(a,g14.6)')'k_perp d_ref     = ',kperp
            write(*,'(a,g14.6)')'k_par  d_ref     = ',kpar
         else
            write(*,'(a,g14.6)')'k_perp rho_ref   = ',kperp
            write(*,'(a,g14.6)')'k_par  rho_ref   = ',kpar
         endif
         write(*,'(a,g14.6)')'Beta_ref,par         = ',betap
         if (normalization == 'd_ref') then
            write(*,'(a,g14.6)')'v_A,ref/c            = ',vtp
         else
            write(*,'(a,g14.6)')'v_t,ref,par/c          = ',vtp
         endif
         do is = 1, nspec
            write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
            write(*,'(a,i3)')      'Parameters for Species :',is
            if (normalization == 'd_ref') then
               write(*,'(a,g14.6)')'2T_||s/(m_ref v_A^2) = ',spec(is)%tau_s
            else
               write(*,'(a,g14.6)')'T_||,ref/T_||s =    ',spec(is)%tau_s
            endif
            write(*,'(a,g14.6)')'m_ref/m_s =        ',spec(is)%mu_s
            write(*,'(a,g14.6)')'T_perp/T_par|s = ',spec(is)%alph_s
            write(*,'(a,g14.6)')'q_ref/q_s =        ',spec(is)%Q_s
            write(*,'(a,g14.6)')'n_s/n_ref =        ',spec(is)%D_s
            write(*,'(a,g14.6)')'v_drift s/v_A,ref =    ',spec(is)%vv_s
         enddo
         write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
         write(*,'(a)')'Searching over:'
         write(*,'(a,es10.3,a,es10.3,a)')'omega/Omega_ref  \in [',omi,',',omf,']'
         write(*,'(a,es10.3,a,es10.3,a)')'gamma/Omega_ref \in [',gami,',',gamf,']'
         write(*,'(a)')'-=-=-=-=-=-=-=-=-=-'
      endif
    
    write(fmt,'(a,i0,a)')'(6es14.4,',6*nspec,'es14.4)'

    invalid_points = 0
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
           if ((has_cold_species() .or. z_backend=='faddeeva') .and. dispersion_status /= 0) then
              val(ir,ii)=huge(1.) ! Search mask, never a substitute tensor/determinant.
              invalid_points=invalid_points+1
           endif
        enddo
     enddo

     if (invalid_points > 0) then
        if (has_cold_species()) then
           write(*,'(a,i0,a)') &
                'WARNING: excluded ',invalid_points,' singular/nonfinite mixed map points from root search.'
        else
           write(*,'(a,i0,a)') &
                'WARNING: excluded ',invalid_points,' singular/nonfinite map points from root search.'
        endif
     endif
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
              if ((has_cold_species() .or. z_backend=='faddeeva') .and. iflag < 0) cycle
              
              !check to see if the found root is already
              !an identified solution
              if ((j.eq.1 .and. .not. (has_cold_species() .or. z_backend=='faddeeva')) .or. &
                   ((has_cold_species() .or. z_backend=='faddeeva') .and. k==0)) then
                 k=k+1
                 wroots(1,k)=real(omega)
                 wroots(2,k)=aimag(omega)
              else
                 unique_root=.true.
                 do ig=1,merge(k,j-1,(has_cold_species() .or. z_backend=='faddeeva')) ! Accepted mixed roots; legacy behavior otherwise.
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
       if ((has_cold_species() .or. z_backend=='faddeeva')) nroots=k
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
             if ((has_cold_species() .or. z_backend=='faddeeva')) then
                if (.not. finite_complex(dal(ir,ii))) val(ir,ii)=real(invalid_complex())
             endif
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

        if (normalization == 'd_ref' .and. option /= 4) then
           do ii=1,nspec
              params(1,ii) = spec(ii)%tau_s
              params(2,ii) = spec(ii)%mu_s
              params(3,ii) = spec(ii)%alph_s
              params(4,ii) = spec(ii)%q_s
              params(5,ii) = spec(ii)%D_s
              params(6,ii) = spec(ii)%vv_s
           enddo
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
       do is = 1,merge(nroots,nroot_max,(has_cold_species() .or. z_backend=='faddeeva'))
          write(*,'(i3,2es14.4)')is,wroots(1:2,is)
       enddo
    endif


  end subroutine map_search

  !-=-=-=-=-
  !-=-=-=-=-
  subroutine test_disp
    !! Testing routine for single evaluation of dispersion relation.   
    use vars, only : wroots, nroot_max, writeOut, nspec, scan
    use vars, only : kperp, kpar, solver_vtp, has_cold_species, has_cold_reference, dispersion_status,z_backend

    use electromagnetic_fields, only: calc_em_eigen
    implicit none

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
     
     real, dimension(1:7,1:nspec) :: Ps_split
     !! Power into/out of species/components, broken into different contributions.
     !! calc_eigen writes all seven channels; retain the first four in stdout.
     

    
    
    complex :: omega
    !! Complex Frequency input.
    
    complex :: D
    !! Output of dispersion relation function.

    integer :: i,j,em_status
    
    

    omega= cmplx(wroots(1,1),wroots(2,1))
    
    D = disp(omega)

    if (has_cold_reference()) then
       write(*,'(a)') 'Cold reference: omega and raw det(Lambda)'
       write(*,'(4es14.4)') omega,D
    else
       write(*,'(4es14.4)') omega,D*solver_vtp**6.
    endif
    write(*,'(2es14.4)')kperp,kpar

    if (has_cold_species()) then
       write(*,'(a,i0)') 'Mixed dispersion status (0=valid, 1=input, 2-4=pole, 5=numerical): ',dispersion_status
       call calc_em_eigen(omega,ef,bf,em_status,disp)
       if (em_status==0) then
          write(*,*) 'ef ='
          write(*,'(3(1x,es14.6))') ef
          write(*,*) 'bf ='
          write(*,'(3(1x,es14.6))') bf
          write(*,'(a)') 'Electromagnetic fields only; species moments, wave energy, heating and FPC unavailable.'
       else
          write(*,'(a,i0)') 'Electromagnetic eigenfunction unavailable; determinant only. Field status: ',em_status
       endif
       return
    endif
    if (z_backend=='faddeeva' .and. dispersion_status/=0) then
       write(*,'(a,i0)') 'Dispersion evaluation failed; diagnostics unavailable. Status: ',dispersion_status
       return
    endif
    call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,.true.,.true.)


    ! --- Print diagnostics ---
    write(*,*) 'ef ='
    write(*,'(3(1x,es14.6))') ef

    write(*,*) 'bf ='
    write(*,'(3(1x,es14.6))') bf

    write(*,*) 'ns ='
    write(*,'(100(1x,es14.6))') (ns(i), i=1,nspec)

    write(*,*) 'Us ='
    do i=1,nspec
       write(*,'(3(1x,es14.6))') (Us(j,i), j=1,3)
    end do

    write(*,*) 'Ps ='
    write(*,'(100(1x,es14.6))') (Ps(i), i=1,nspec)

    write(*,*) 'Ps_split ='
    do i=1,nspec
       write(*,'(4(1x,es14.6))') (Ps_split(j,i), j=1,4)
    end do
    

  end subroutine test_disp

!-=-=-=-=-

!-=-=-=-=-
subroutine refine_guess
  !!Refine solutions for roots of dispersion relation.
  use vars, only : wroots, nroot_max, writeOut, has_cold_species,z_backend
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
     if ((has_cold_species() .or. z_backend=='faddeeva') .and. all(wroots(1:2,j)==0.)) &
          error stop 'ERROR: root guess is at the zero-frequency pole.'
     if (wroots(1,j) .ne. 0. .or. wroots(2,j) .ne. 0.) then
        
        omega=cmplx(wroots(1,j),wroots(2,j))
        
        om1=omega*(1.-prec)
        om2=omega*(1.+prec)
        
        omega=rtsec(disp,om1,om2,tol,iflag)
        if ((has_cold_species() .or. z_backend=='faddeeva') .and. iflag < 0) &
             error stop 'ERROR: root refinement failed: singular/nonfinite trial or iteration limit.'

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
  use scan_support, only: scan_frame,initialize_scan_frame,scan_coordinates,check_scan_kpar, &
       wavevector_scan,refine_scan_point
  use vars, only: has_cold_reference,has_cold_species

  !! Calculate complex frequency maps for scanned parameters.
  !! Intervatively calls [[map_search(subroutine)]]
  !! over a fixed range of parameters.
  use vars, only : scan, spec, print_Name, sw, sw2, kperp, kpar
  use vars, only : writeOut, pi
  use vars, only : gami, gamf, omi, omf ,betap
  use vars, only : collision_type, Kn, nspec, normalization, reference_scan_source
  !Local

  integer :: jj
  !! Parameter step index.

  integer :: iq
  !! species index
  
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

  type(scan_frame) :: geometry
  if (has_cold_species()) call initialize_scan_frame(geometry,kperp,kpar)
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
     if (has_cold_species()) then
        call scan_coordinates(geometry,scan(1),jj,kperp,kpar)
     else
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

     endif
     if (wavevector_scan(scan(1))) call check_scan_kpar()
     do iq=1, nspec
        if (collision_type.gt.0 .and. normalization=='rho_ref') &
             spec(iq)%nu_ns = (sqrt(2.d0)*Kn)**(-1.d0)* &
             sqrt(spec(iq)%mu_s/(spec(iq)%tau_s*spec(1)%alph_s))
     enddo

     if ((scan(1)%style_s)==-1) then
        !Scans with multiple components
        if ((scan(1)%type_s)==0) then
           !k_0 -> k_1; 
           !k_0= current k
           !kperp_1 = range_i, kpar_1 = range_f
           if (.not.has_cold_species()) then
              kpari=kpar
              kperpi=kperp
           endif
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
                trim(outName),int(10000*sign(sqrt(kperp**2+kpar**2),kpar))
        endif
     else
        write(print_Name,'(a,i0)')&
             trim(outName),int(10000*sw)
        write(*,*)trim(print_Name)
     endif

     !Edit map span. Hard code in the desired range.
     !Will add an improved interface for the future.
     
     if (.not.has_cold_species()) then
     write(*,'(a)')'Modifying Span of Complex Frequency Space:'
     if (normalization=='d_ref') then
        omi=-3.*sqrt(kperp**2.+kpar**2.)*sqrt(spec(1)%alph_s)
        omf=-omi
        gami=omi
        gamf=0.5*sqrt(kperp**2.+kpar**2.)*sqrt(spec(1)%alph_s)
     else
        omi=-3.*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
        omf=3.*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
        gami=-3.*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
        gamf=0.5*sqrt(kperp**2.+kpar**2.)/sqrt(betap)
     endif
     
     endif
     call map_search

  enddo  
  !End Parameter Scan
  
end subroutine map_scan

!-=-=-=-=-=
!-=-=-=-=-=
subroutine om_scan(is)
  use electromagnetic_fields, only: calc_mixed_eigen
  use scan_support, only: scan_frame,initialize_scan_frame,scan_coordinates,check_scan_kpar, &
       wavevector_scan,refine_scan_point
  use vars, only: has_cold_reference,has_cold_species

  !! Calculate solutions as a function of the variation of
  !! a single parameter.
  use vars, only : wroots,scan,nroot_max,nspec,sw,sw2, low_n, pi
  use vars, only : kperp,kpar,betap,vtp,spec,writeOut,susc
  use vars, only : collision_type, Kn, normalization, sync_solver_state
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
  
  logical :: root_valid
  integer :: em_status
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
    
  !>>>GGH: 1/18/23
  real, dimension(1:7,1:nspec) :: Ps_split
  !!Power into/out of species/componets.
  !!Corrected LD/TTD calculation (GGH).
  !!Split \pm 1 CD (KGK)

  real :: Ew
  !!Wave energy.

  real, dimension(1:6,1:nspec) :: params
  !!Species parameters:

  type(scan_frame) :: geometry
  if (has_cold_species()) call initialize_scan_frame(geometry,kperp,kpar)
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
     if (has_cold_species()) then
        call scan_coordinates(geometry,scan(is),jj,kperp,kpar)
     else
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

     endif
     if (wavevector_scan(scan(is))) call check_scan_kpar()
     if (.not.has_cold_reference()) call sync_solver_state
     !construct params array for formated output
     do ii = 1,nspec
        !calculate neutral-charged collision frequency
        if (collision_type.gt.0 .and. normalization=='rho_ref') &
             spec(ii)%nu_ns = (sqrt(2.d0)*Kn)**(-1.d0)* &
             sqrt(spec(ii)%mu_s/(spec(ii)%tau_s*spec(1)%alph_s))
        
        if (mod(jj,scan(is)%n_res)==0) then
           params(1,ii) = spec(ii)%tau_s
           params(2,ii) = spec(ii)%mu_s
           params(3,ii) = spec(ii)%alph_s
           params(4,ii) = spec(ii)%q_s
           params(5,ii) = spec(ii)%D_s
           params(6,ii) = spec(ii)%vv_s
        endif
     enddo

     !Root Scan....
     do ii=1,nroot_max
        !Bracket values for root search
        omold=omlast(ii)
        om1=omold*(1.-prec)
        om2=omold*(1.+prec)

        iflag=0
        !Find New Omega Value
        call refine_scan_point(om1,om2,tol,omega,iflag,rtsec,disp)
        root_valid=iflag>0
        if (.not.root_valid) then
           ef=omega;bf=omega;ns=omega;Us=omega
           Ps=real(omega);Ps_split=real(omega);Ew=real(omega)
        endif
        
        !Save Root
        if (root_valid) omlast(ii)=omega
        
        !Output value every %n_res steps
        if (mod(jj,scan(is)%n_res)==0) then
           !Calculate eigenfns, heating only every %n_res steps
           
           if (root_valid .and. (scan(is)%eigen_s.or.scan(is)%heat_s)) then
              if (root_valid) val=abs(disp(omega))
              !>>>GGH: 1/18/23
              if (has_cold_species()) then
                 call calc_mixed_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,scan(is)%heat_s,em_status,disp)
              else
              call calc_eigen(omega,ef,bf,Us,ns,Ps,&
                   Ps_split,scan(is)%eigen_s,scan(is)%heat_s)
              endif
              !<<<GGH: 1/18/23
              omega=omlast(ii)
              if (.not.has_cold_species()) then
              if (abs(real(omega)).lt.1.E-15) then
                 Ps=0.1
                 Ps_split=-0.1
                 !<<<GGH: 1/18/23
              endif
              endif
           endif

           if (root_valid) then
              if (abs(real(omega)).lt.1.E-15) omega=cmplx(0.,aimag(omega))
           endif

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
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                      Ps(1:nspec),Ps_split(1:7,1:nspec),&
                      params(1:6,1:nspec)
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
                 write(out_unit(ii),fmt)&
                      kperp,kpar,betap,vtp,&
                      omega,&            
                      Ps(1:nspec),Ps_split(1:7,1:nspec),&
                      params(1:6,1:nspec)
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
              if (root_valid) val=abs(disp(omega))
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
           wroots(1,ii)=real(omlast(ii))
           wroots(2,ii)=aimag(omlast(ii))
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
  use electromagnetic_fields, only: calc_mixed_eigen
  use scan_support, only: scan_frame,initialize_scan_frame,scan_coordinates,check_scan_kpar, &
       wavevector_scan,refine_scan_point
  use vars, only: has_cold_reference,has_cold_species

  !! Calculate solutions as a function of the variation of
  !! two parameters, creating a surface in parameter space.
  use vars, only: scan, spec, betap, kpar, kperp, vtp, sw, sw2, sw3, sw4, pi
  use vars, only: nroot_max, outputName, wroots, nspec, writeOut, susc, low_n
  use vars, only : collision_type, Kn, normalization, sync_solver_state
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

  integer :: iq
  !! species index

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

  real :: ktmp
  !! Temporary value of \( |k| \), for spliting into \( \perp )\ and \( \parallel \) components.
  
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
  
  logical :: root_valid
  integer :: em_status
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
    
  !>>>GGH: 1/18/23
  real, dimension(1:7,1:nspec) :: Ps_split
  !!Power into/out of species/componets.
  !!Corrected LD/TTD calculation (GGH).
  !!Split \pm 1 CD (KGK)

  real :: Ew
  !!Wave energy.

  real, dimension(1:6,1:nspec) :: params
  !!Species parameters:
  
  type(scan_frame) :: geometry
  if (has_cold_species()) call initialize_scan_frame(geometry,kperp,kpar)
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
           write(*,'(a,2es14.4)')'Initial Value (theta,|k|rho_ref)',theta*180/pi,ki
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

     if (has_cold_species()) then
        call scan_coordinates(geometry,scan(1),kk,kperp,kpar,scan(2),0)
     else
     !Advance scanned parameter values

     !Global Two-Component Scans
     if ((scan(1)%style_s)==-1)then

        !Scan from (kpar_0;kperp_0) to (kpar_1;kperp_1)
        if ((scan(1)%type_s)==0) then           
           if (scan(1)%log_scan) then
              sw=10.**(log10(kperpi)+diff(1,1)*real(kk)) !kperp
              sw2=10.**(log10(kpari)+diff(1,2)*real(kk)) !kpar
           else
              sw=(kperpi)+diff(1,1)*real(kk) !kperp
              sw2=(kpari)+diff(1,2)*real(kk) !kpar
           endif

        !Scan from (theta_0) to (theta_1)
        elseif ((scan(1)%type_s)==1) then
           if (scan(1)%log_scan) then
              theta_q=10.**(log10((theta))+diff(1,1)*real(kk))    
           else
              theta_q=(((theta))+diff(1,1)*real(kk))    
           endif
           sw=(ki*sin(theta_q))!kperp
           sw2=(ki*cos(theta_q))!kpar

        !Scan from (|k_0|) to (|k_1|) at constant theta
        elseif ((scan(1)%type_s)==2) then
           if (scan(1)%log_scan) then
              ktmp=10.**(log10(ki)+diff(1,1)*real(kk)) !k
              sw= ktmp*sin(theta)!kperp
              sw2= ktmp*cos(theta)!kpar
           else                                        
              ktmp=ki+diff(1,1)*real(kk) !k
              sw= ktmp*sin(theta)!kperp
              sw2= ktmp*cos(theta)!kpar
           endif
        endif

     !Single Component Scans
     else
        if (scan(1)%log_scan) then
           sw=10.**(log10(scan(1)%range_i)+diff(1,1)*real(kk))    
        else
           sw=(scan(1)%range_i)+diff(1,1)*real(kk)    
        endif
     endif

     endif
     if (wavevector_scan(scan(1))) call check_scan_kpar()
     !calculate neutral-charged collision frequency
     do iq=1, nspec
        if (collision_type.gt.0 .and. normalization=='rho_ref') &
             spec(iq)%nu_ns = (sqrt(2.d0)*Kn)**(-1.d0)* &
             sqrt(spec(iq)%mu_s/(spec(iq)%tau_s*spec(1)%alph_s))
     enddo
     
     !Root Scan....
     do ii=1,nroot_max
        !Bracket values for root search
        omold=omlast(ii)
        om1=omold*(1.-prec)
        om2=omold*(1.+prec)

        iflag=0
        !Find New Omega Value
        call refine_scan_point(om1,om2,tol,omega,iflag,rtsec,disp)
        root_valid=iflag>0
        if (.not.root_valid) then
           ef=omega;bf=omega;ns=omega;Us=omega
           Ps=real(omega);Ps_split=real(omega);Ew=real(omega)
        endif
        
        !Save Root
        if (root_valid) omlast(ii)=omega

     enddo

     !Transition from parameter 1 scan to parameter 2 scan
     if (mod(kk,scan(1)%n_res)==0) then
        if (writeOut) then
           if ((scan(1)%style_s)==-1)then
              write(*,'(a,2es11.4)')&
                   '(kperp,kpar): ',sw,sw2
           else
              write(*,'(a,es11.4)')&
                   'Parameter 1: ',sw
           endif
        endif
        !Save roots
        do ii = 1,nroot_max
           omSafe(ii) = omlast(ii)
        enddo
    
        !Parameter 2 scan
        !Output %n_scan steps, with %n_res steps inbetween each output
        do jj = 0, (scan(2)%n_scan*scan(2)%n_res) 
           if (has_cold_species()) then
              call scan_coordinates(geometry,scan(1),kk,kperp,kpar,scan(2),jj)
           else
           !Advance scanned parameter values

           !Global Two-Component Scans
           if ((scan(2)%style_s)==-1)then
              
              !Scan from (kpar_0;kperp_0) to (kpar_1;kperp_1)
              if ((scan(2)%type_s)==0) then
                 if (scan(2)%log_scan) then
                    sw3=10.**(log10(kperpi)+diff(2,1)*real(jj)) !kperp
                    sw4=10.**(log10(kpari)+diff(2,2)*real(jj))  !kpar
                 else
                    sw3=(kperpi)+diff(2,1)*real(jj) !kperp
                    sw4=(kpari)+diff(2,2)*real(jj)  !kpar
                 endif

              !Scan from (theta_0) to (theta_1)
              elseif ((scan(2)%type_s)==1) then
                 if (scan(2)%log_scan) then
                    theta_q=10.**(log10((theta))+diff(2,1)*real(jj))    
                 else
                    theta_q=(((theta))+diff(2,1)*real(jj))    
                 endif
                 sw3=(ki*sin(theta_q))!kperp
                 sw4=(ki*cos(theta_q))!kpar

                 !Scan from (|k_0|) to (|k_1|) at constant theta
              elseif ((scan(2)%type_s)==2) then
                 if (scan(2)%log_scan) then                    
                    ktmp=10.**(log10(ki)+diff(2,1)*real(jj)) !k
                    sw3= ktmp*sin(theta_q)!kperp
                    sw4= ktmp*cos(theta_q)!kpar
                 else                                        
                    ktmp=ki+diff(2,1)*real(jj) !k
                    sw3= ktmp*sin(theta_q)!kperp
                    sw4= ktmp*cos(theta_q)!kpar
                 endif
              endif

              !Single Component Scans
           else
              if (scan(2)%log_scan) then
                 sw3=10.**(log10(scan(2)%range_i)+diff(2,1)*real(jj))    
              else
                 sw3=(scan(2)%range_i)+diff(2,1)*real(jj)    
              endif
           endif

           endif
           if (wavevector_scan(scan(2))) call check_scan_kpar()
           if (.not.has_cold_reference()) call sync_solver_state
           !construct params array for formated output
           do ii = 1,nspec
              !calculate neutral-charged collision frequency
              if (collision_type.gt.0 .and. normalization=='rho_ref') &
                   spec(ii)%nu_ns = (sqrt(2.d0)*Kn)**(-1.d0)* &
                   sqrt(spec(ii)%mu_s/(spec(ii)%tau_s*spec(1)%alph_s))
              if (mod(jj,scan(2)%n_res)==0) then
                 params(1,ii) = spec(ii)%tau_s
                 params(2,ii) = spec(ii)%mu_s
                 params(3,ii) = spec(ii)%alph_s
                 params(4,ii) = spec(ii)%q_s
                 params(5,ii) = spec(ii)%D_s
                 params(6,ii) = spec(ii)%vv_s                 
              endif
           enddo
           
           !Root Scan....
           do ii=1,nroot_max
              !Bracket values for root search
              omold=omlast(ii)
              om1=omold*(1.-prec)
              om2=omold*(1.+prec)
              
              iflag=0
              !Find New Omega Value
              call refine_scan_point(om1,om2,tol,omega,iflag,rtsec,disp)
              root_valid=iflag>0
              if (.not.root_valid) then
                 ef=omega;bf=omega;ns=omega;Us=omega
                 Ps=real(omega);Ps_split=real(omega);Ew=real(omega)
              endif
              
              !Save Root
              if (root_valid) omlast(ii)=omega
              
              !Output value every %n_res steps
              if (mod(jj,scan(2)%n_res)==0) then
                 !Calculate eigenfns, heating only every %n_res steps
                 
                 if (root_valid .and. (scan(2)%eigen_s.or.scan(2)%heat_s)) then
                    if (root_valid) val=abs(disp(omega))
                    !>>>GGH: 1/18/23
                    if (has_cold_species()) then
                       call calc_mixed_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,scan(2)%heat_s,em_status,disp)
                    else
                    call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,scan(2)%eigen_s,scan(2)%heat_s)
                    endif
                    !<<<GGH: 1/18/23
                    if (.not.has_cold_species()) then
                    if (abs(real(omega)).lt.1.E-7) then
                       Ps=-0.1
                       Ps_split=1.0
                       !>>>GGH: 1/18/23
                    endif
                    endif

                 endif
                 
                 select case(out_type)
                 case(0) !Om, Eigen, Heating                    
                    if (low_n) then
                       !>>>GGH: 1/18/23
                          write(out_unit(ii),fmt)&
                               kperp,kpar,betap,vtp,&
                               omega,&            
                               bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                               Ps(1:nspec),Ps_split(1:7,1:nspec),&
                               params(1:6,1:nspec)
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
                          write(out_unit(ii),fmt)&
                               kperp,kpar,betap,vtp,&
                               omega,&            
                               Ps(1:nspec),Ps_split(1:7,1:nspec),&
                               params(1:6,1:nspec)
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
                 if (root_valid) val=abs(disp(omega))
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
           if (has_cold_species()) then
              call scan_coordinates(geometry,scan(1),kk,kperp,kpar,scan(2),0)
           else
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
                 if (scan(2)%log_scan) then                    
                    ktmp=10.**(log10(ki)+diff(2,1)*real(jj)) !k
                    sw3= ktmp*sin(theta)!kperp
                    sw4= ktmp*cos(theta)!kpar
                 else                                        
                    ktmp=ki+diff(2,1)*real(jj) !k
                    sw3= ktmp*sin(theta)!kperp
                    sw4= ktmp*cos(theta)!kpar
                 endif
                 !k along contant theta
                 !if (scan(2)%log_scan) then
                 !   sw3=10.**(log10(ki*sin(theta))+diff(2,1)*real(jj))    
                 !   sw4=10.**(log10(ki*cos(theta))+diff(2,2)*real(jj))    
                 !else
                 !   sw3=(ki*sin(theta))+diff(2,1)*real(jj)    
                 !   sw4=(ki*cos(theta))+diff(2,2)*real(jj)    
                 !endif
              endif
           else
              if (scan(2)%log_scan) then
                 sw3=10.**(log10(scan(2)%range_i)+diff(2,1)*real(jj))    
              else
                 sw3=(scan(2)%range_i)+diff(2,1)*real(jj)    
              endif
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
  use vars, only : nRad, nRad_sub, radius, writeOut, nroot_max
  use vars, only : wroots, spec, rad_spec, nspec, k_scan
  use vars, only : vtp_rad, beta_rad, sw, sw2, sw3, sw4
  use vars, only : kperp, kpar, betap, vtp, radial_eigen
  use vars, only : dataName, outputName, radial_heating
  use vars, only : rad_scan, pi
  use vars, only : collision_type, Kn, normalization, sync_solver_state
  use functions, only : get_unused_unit
  implicit none
  
  !Local
  integer :: ir
  !!Radial index.

  integer :: ir_sub
  !!Radial substep index.

  integer :: ir_out
  !!Output radial index, including interpolated points.
  
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

  real :: frac
  !! Radial interpolation fraction.

  real :: radius_current
  !! Current radial position, including interpolated points.

  real :: kperp_start
  !! Starting perpendicular wavevector for radial continuation.

  real :: kpar_start
  !! Starting parallel wavevector for radial continuation.

  nRad_sub=max(1,nRad_sub)
  kperp_start=kperp
  kpar_start=kpar
  
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
        write(fmt,'(a,i0,a)')'(7es14.4,12es14.4,',15*nspec+1,'es14.4)'
        out_type=0
     else
        write(fmt,'(a,i0,a)')'(7es14.4,12es14.4,',14*nspec+1,'es14.4)'
        out_type=1
     endif
  else
     if (radial_heating) then
        write(fmt,'(a,i0,a)')'(7es14.4,',7*nspec+1,'es14.4)'
        out_type=2
     else
        write(fmt,'(a,i0,a)')'(7es14.4,',6*nspec+1,'es14.4)'
        out_type=3
     endif
  endif

  write(*,'(2a)')'Number of Columns:',trim(fmt)
  
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

  !Scan through radial trajectory in parameter space.
  ir_out=0
  do ir = 0, nRad-2
     do ir_sub = 0, nRad_sub-1
        frac=real(ir_sub)/real(nRad_sub)
        call set_radial_state(ir,ir+1,frac,radius_current)
        call report_radius(ir_out,radius_current)
        call run_radial_k_scan(ir_sub==0,radius_current)
        ir_out=ir_out+1
     enddo
  enddo

  call set_radial_state(nRad-1,nRad-1,0.0,radius_current)
  call report_radius(ir_out,radius_current)
  call run_radial_k_scan(.true.,radius_current)

contains

  subroutine set_radial_state(ir0,ir1,frac_in,radius_out)
    integer, intent(in) :: ir0,ir1
    real, intent(in) :: frac_in
    real, intent(out) :: radius_out

    radius_out=(1.-frac_in)*radius(ir0)+frac_in*radius(ir1)
    betap=(1.-frac_in)*beta_rad(ir0)+frac_in*beta_rad(ir1)
    vtp=(1.-frac_in)*vtp_rad(ir0)+frac_in*vtp_rad(ir1)

    do is = 1,nspec
       spec(is)%tau_s  = (1.-frac_in)*rad_spec(is,ir0)%tau_s  + &
            frac_in*rad_spec(is,ir1)%tau_s
       spec(is)%mu_s   = (1.-frac_in)*rad_spec(is,ir0)%mu_s   + &
            frac_in*rad_spec(is,ir1)%mu_s
       spec(is)%alph_s = (1.-frac_in)*rad_spec(is,ir0)%alph_s + &
            frac_in*rad_spec(is,ir1)%alph_s
       spec(is)%q_s    = (1.-frac_in)*rad_spec(is,ir0)%q_s    + &
            frac_in*rad_spec(is,ir1)%q_s
       spec(is)%D_s    = (1.-frac_in)*rad_spec(is,ir0)%D_s    + &
            frac_in*rad_spec(is,ir1)%D_s
       spec(is)%vv_s   = (1.-frac_in)*rad_spec(is,ir0)%vv_s   + &
            frac_in*rad_spec(is,ir1)%vv_s

       if (collision_type.gt.0 .and. normalization=='rho_ref') &
            spec(is)%nu_ns = (sqrt(2.d0)*Kn)**(-1.d0)* &
            sqrt(spec(is)%mu_s/(spec(is)%tau_s*spec(1)%alph_s))

       params(1,is) = spec(is)%tau_s
       params(2,is) = spec(is)%mu_s
       params(3,is) = spec(is)%alph_s
       params(4,is) = spec(is)%q_s
       params(5,is) = spec(is)%D_s
       params(6,is) = spec(is)%vv_s
    enddo
    call sync_solver_state
  end subroutine set_radial_state

  subroutine report_radius(ir_print,radius_in)
    integer, intent(in) :: ir_print
    real, intent(in) :: radius_in

    write(*,'(a,i0,a,es14.4)')&
         'Parameter (',ir_print,') = ',radius_in
  end subroutine report_radius

  subroutine set_start_wavevector
    kperp=kperp_start
    kpar=kpar_start
  end subroutine set_start_wavevector

  subroutine save_start_roots
    do ii=1,nroot_max
       omSafe(1,ii)=omLast(ii)
    enddo
  end subroutine save_start_roots

  subroutine run_start_wavevector(radius_in)
    real, intent(in) :: radius_in

    call set_start_wavevector
    do ii = 1,nroot_max
       omLast(ii)=omSafe(1,ii)
    enddo
    call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,.false.)
    call save_start_roots
  end subroutine run_start_wavevector

  subroutine run_radial_k_scan(full_scan,radius_in)
    logical, intent(in) :: full_scan
    real, intent(in) :: radius_in

    call set_start_wavevector

    if ((.not.full_scan).and.(k_scan.ne.0)) then
       call run_start_wavevector(radius_in)
       return
    endif

     !Scan (or don't) through kspace
     select case(k_scan)
     case(0)
        !simple case of constant kperp, kpar

        !Root Finder in separate subroutine
        call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,full_scan)
        
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
           call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,mod_write)
           
           if (ij==0) then
              !Save roots for next scan of k
              do ii=1,nroot_max
                 omSafe(1,ii)=omLast(ii)
                 write(*,*)'root saved',ii,omlast(ii)
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
           call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,mod_write)
           
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
           call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,mod_write)
           
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
           call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,mod_write)
           
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
              call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,mod_write)
              
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
              call om_radial(omlast,params,out_unit,fmt,out_type,radius_in,mod_write)
              
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
  end subroutine run_radial_k_scan

end subroutine radial_scan

!-=-=-=-=-
!-=-=-=-=-
subroutine om_radial(omlast,params,out_unit,fmt,out_type,radius_current,mod_write)
  !!Scans roots at a given location in parameter space.
  !!Used in conjunction with [[radial_scan(subroutine)]].
  use vars, only: nspec, nroot_max, radial_heating, radial_eigen
  use vars, only: betap, vtp, kperp, kpar
  !Passed

  complex,dimension(1:nroot_max) :: omlast
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

  real :: radius_current
  !!Current radial position, including interpolated points.

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
    
  !>>>GGH: 1/18/23
  real, dimension(1:7,1:nspec) :: Ps_split
  !!Power into/out of species/componets.
  !!Corrected LD/TTD calculation (GGH).
  !!Split \pm 1 CD (KGK)

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
              call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,radial_eigen,radial_heating)
           endif
           
           !Output results
           if (out_type==0) & !Om, Eigen, Heating
                write(out_unit(ii),fmt)&
                kperp,kpar,betap,vtp,&
                omega,&            
                bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                Ps(1:nspec),&
                params(1:6,1:nspec),Ew, radius_current
           if (out_type==1) & !Om, Eigen
                write(out_unit(ii),fmt)&
                kperp,kpar,betap,vtp,&
                omega,&            
                bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
                params(1:6,1:nspec), radius_current
           if (out_type==2) & !Om, Heating
                write(out_unit(ii),fmt)&
                kperp,kpar,betap,vtp,&
                omega,&            
                Ps(1:nspec),&
                params(1:6,1:nspec),Ew, radius_current
           if (out_type==3) & !Om
                write(out_unit(ii),fmt)&
                kperp,kpar,betap,vtp,&
                omega,&            
                params(1:6,1:nspec),radius_current
        endif
     enddo

end subroutine om_radial

!-=-=-=-=-
!------------------------------------------------------------------------------
!Greg Howes, 2010; Kristopher Klein, 2015
!------------------------------------------------------------------------------
!-=-=-=-=-=-
subroutine calc_eigen(omega,electric,magnetic,vmean,ns,Ps,&
     Ps_split,eigen_L,heat_L)
  !!  Calculates the electric and magnetic fields as well as species
  !!     velocities and density fluctuations for identified wave 
  !!     as well as the power emission or absorption.
!<<<GGH: 1/18/23
  use vars, only : spec=>solver_spec,betap=>solver_betap,vtp=>solver_vtp
  use vars, only : kperp=>solver_kperp,kpar=>solver_kpar,nspec,susc,lam,low_n,susc_low
  use vars, only : has_cold_species
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
  
  real, dimension(1:7,1:nspec), intent(out) :: Ps_split
  !!<<<GGH: 1/18/23
  !!Power into/out of species/components from LD, TTD, CD.
  !!Split \pm 1 CD (KGK)

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

  if (has_cold_species()) &
       error stop 'ERROR: cold-fluid eigenfunctions and heating diagnostics are unavailable in Stage 2.'

  !CALCULATE FIELDS FLUCTUATIONS==========================================
  !Calculate Electric Fields
  !This seems to be correct --GGH 20 JUN 2008
  include 'electromagnetic_electric.inc'

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
     !    normalized to (c Eperp1/B0)
     !KGK: 1-28-2015: Added the drift velocity effects
     vmean(:,:)=0.
     do jj = 1,nspec !Species velocity fluctuations
        if (spec(jj)%vv_s.eq.0.) then
           do j=1,3!x,y,z
              !vmean(j,jj) = -(spec(jj)%Q_s/spec(jj)%D_s)*cmplx(0.,1.)*&
              !     omega/sqrt(betap)*vtp**2. *sqrt(spec(1)%alph_s)* &
              !     sum(electric(:)*susc(jj,j,:))
              !  9 AUG 2023: GGH/KGK Fixed U calculation
              !        vmean(j,jj) = -(spec(jj)%Q_s/spec(jj)%D_s)*cmplx(0.,1.)*&
              !             omega/sqrt(betap)*vtp* &
              !             sum(electric(:)*susc(jj,j,:))
              ! 28 AUG 2023: Corrected Version:  Normalize Us1 to (c Eperp1/B0)
              vmean(j,jj) = -cmplx(0.,1.)*omega*sum(electric(:)*susc(jj,j,:)) &
                   *(spec(jj)%Q_s/spec(jj)%D_s)*vtp**2./betap 
           enddo
        else
           !In the non-zero drift velocity case, the density and velocity fluctuations are coupled.
           !This requires solving U_x, which is then used in the solution of U_z.
           !U_x and U_z are then subsequently used in the determination of n1/n0).
           do j=1,2!x,y
              vmean(j,jj) = -cmplx(0.,1.)*omega*sum(electric(:)*susc(jj,j,:)) &
                   *(spec(jj)%Q_s/spec(jj)%D_s)*vtp**2./betap
           enddo
           j=3!z

           vmean(j,jj) = -cmplx(0.,1.)*omega*sum(electric(:)*susc(jj,j,:)) &
                *(spec(jj)%Q_s/spec(jj)%D_s)*vtp**2./betap &
                -vmean(1,jj)*(kperp*spec(jj)%vv_s)/&
                (sqrt(betap*spec(1)%alph_s)*omega-kpar*spec(jj)%vv_s)
           vmean(j,jj) = vmean(j,jj)/&
                (1+(kpar*spec(jj)%vv_s)/&
                (omega*sqrt(betap*spec(1)%alph_s)-kpar*spec(jj)%vv_s))

        endif
        enddo

     !CALCULATE DENSITY FLUCTUATIONS========================================
     ! This is ns/ns0
     !KGK: 8-28-2013: switched from real(omega) to the complex valued omega.
     !KGK: 12-19-2013: Added the sqrt(T_perp/T_par) factor
     !KGK: 1-28-2015: Added the drift velocity Doppler shift
     do jj=1,nspec
        !ns(jj) = (vmean(1,jj)*kperp+vmean(3,jj)*kpar)/&
        !     ((omega-kpar * spec(jj)%vv_s/sqrt(betap*spec(1)%alph_s))&
        !     *sqrt(betap*spec(1)%alph_s))

        !This is (n1s/n0s)/(Ex/B0)
        ns(jj) = (vmean(1,jj)*kperp+vmean(3,jj)*kpar)/&
             ((omega-kpar * spec(jj)%vv_s/sqrt(betap*spec(1)%alph_s))&
             *sqrt(spec(1)%alph_s)*vtp)
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
     !Greg's New low_n calculation to separate LD/TTD
       !N=0
        do ii = 1, 3 !tensor index
           do j = 1, 3 !tensor index
              do jj = 1, nspec !species index
                 susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                      (susc_low(jj,ii,j,0) - conjg(susc_low(jj,j,ii,0)))
              enddo
           enddo
        enddo

        !Initialize Ps_split
        Ps_split(:,:) = 0.
        
        !chi_yy  (TTD term 1)
        Ps_split(1,:) =-0.5*cmplx(0.,1.)*&
             conjg(electric(2))*electric(2)* &
             (susc_low(:,2,2,0)-conjg(susc_low(:,2,2,0)))
        !chi_yz  (TTD term 2)
        Ps_split(2,:) =-0.5*cmplx(0.,1.)*&
             (electric(3)*conjg(electric(2))*susc_low(:,2,3,0) - &
             conjg(electric(3))*electric(2)*conjg(susc_low(:,2,3,0)))
        !chi_zy  (LD term 1)
        Ps_split(3,:) =-0.5*cmplx(0.,1.)*&
             (electric(2)*conjg(electric(3))*susc_low(:,3,2,0) - &
             conjg(electric(2))*electric(3)*conjg(susc_low(:,3,2,0)))
        !chi_zz  (LD term 2)
        Ps_split(4,:) =-0.5*cmplx(0.,1.)*&
             conjg(electric(3))*electric(3)* &
             (susc_low(:,3,3,0)-conjg(susc_low(:,3,3,0)))

        !Total n=0 terms
        term(:,:)=0.
        do ii = 1, 3
           do jj = 1, nspec
              term(jj,ii) = sum(conjg(electric(:))*susca(jj,:,ii))     
           enddo
        enddo        
        Ps_split(5,:) = 0.
        do jj = 1, nspec
           Ps_split(5,jj) = sum(term(jj,:)*electric(:))
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
        Ps_split(6,:) = 0.
        do jj = 1, nspec
           Ps_split(6,jj) = sum(term(jj,:)*electric_xy(:))
        enddo

        !N=-1
        do ii = 1, 3 !tensor index
           do j = 1, 3 !tensor index
              do jj = 1, nspec !species index
                 susca(jj,ii,j) = -0.5*cmplx(0.,1.)* &
                      (susc_low(jj,ii,j,-1) - conjg(susc_low(jj,j,ii,-1)))
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
        Ps_split(7,:) = 0.
        do jj = 1, nspec
           Ps_split(7,jj) = sum(term(jj,:)*electric_xy(:))
        enddo

        !Normalization             
        Ps_split = Ps_split/ewave
        
     endif

     !EndIf (scan(is)%heat_s) loop
  endif

end subroutine calc_eigen
!-=-=-=-=-=-
!-=-=-=-=-=-
subroutine set_map_pointers(outName,diff,diff2)
  use vars, only: has_cold_species
  !! Sets step sizes for all classes of parameter scans
  !! called in for [[map_scan(subroutine)]].
  use vars, only : scan,nspec,dataName,writeOut,pi
  use vars, only : spec,sw,sw2,betap,vtp,kperp,kpar,reference_scan_source
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
  
  reference_scan_source = 0
  if (scan(1)%style_s==0 .and. scan(1)%type_s==2) reference_scan_source = 1
  if (scan(1)%style_s==1 .and. scan(1)%type_s==0) reference_scan_source = 2

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

        if (has_cold_species()) then
           diff=0.;diff2=0. ! Computed by the native coordinate frame instead.
        else
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
  use vars, only: has_cold_species
  !!Sets I/O strings, step sizes, and formats for [[om_scan(subroutine)]].
  use vars, only : scan,nspec,dataName,writeOut,outputName,pi
  use vars, only : spec,sw,sw2,betap,vtp,kperp,kpar,low_n,reference_scan_source
  !>>>GGH: 1/18/23
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
  

  reference_scan_source = 0
  if (scan(is)%style_s==0 .and. scan(is)%type_s==2) reference_scan_source = 1
  if (scan(is)%style_s==1 .and. scan(is)%type_s==0) reference_scan_source = 2

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

        if (has_cold_species()) then
           diff=0.;diff2=0. ! Computed by the native coordinate frame instead.
        else
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
              write(fmt,'(a,i0,a)')'(6es15.6e3,12es15.6e3,',22*nspec,'es15.6e3)'
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
              write(fmt,'(a,i0,a)')'(6es15.6e3,',14*nspec,'es15.6e3)'
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
  use vars, only : normalization, reference_scan_source
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

  reference_scan_source = 0
  do is=1,2
     if (scan(is)%style_s==0 .and. scan(is)%type_s==2) then
        if (reference_scan_source==2 .and. normalization=='d_ref') then
           write(*,'(a)') 'ERROR: d_ref cannot scan betap and species_1 tauS independently.'
           error stop 1
        endif
        reference_scan_source=1
     endif
     if (scan(is)%style_s==1 .and. scan(is)%type_s==0) then
        if (reference_scan_source==1 .and. normalization=='d_ref') then
           write(*,'(a)') 'ERROR: d_ref cannot scan betap and species_1 tauS independently.'
           error stop 1
        endif
        reference_scan_source=2
     endif
  enddo

  if (((scan(1)%style_s)==(scan(2)%style_s)).and.&
       ((scan(1)%type_s)==(scan(2)%type_s))) then
     write(*,'(a)')&
          'Scaning the same parameters in double scan...'
     write(*,'(a)')&
          'Change input file values...'
     stop
  endif

  if ( ((scan(1)%style_s)==-1).and.&
       ((scan(2)%style_s)==-1).and.&
       ((scan(1)%type_s)==2) .and. &
       ((scan(2)%type_s)==1) &
       ) then
     write(*,'(a)')&
          'For (theta,|k| \rho) double scan, set theta as scan_input_1'
     write(*,'(a)')&
          'and |k| as scan_input_2.'
     write(*,'(a)')&
          'HALTING...'
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
              if (writeOut) &
                   write(*,'(a,es15.6e3,a,es15.6e3)')&
                   'Log-spaced scan: d_kperp: ',diff(is,1),' d_kpar: ',diff(is,2)

           else
              !Linear spacing
              diff(is,1)=((scan(is)%range_i)-(kperp))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              diff(is,2)=((scan(is)%range_f)-(kpar))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              if (writeOut) &
                   write(*,'(a,es15.6e3,a,es15.6e3)')&
                   'Linear scan: d_kperp: ',diff(is,1),' d_kpar: ',diff(is,2)
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

              if (writeOut) &
                   write(*,'(a,es15.6e3)')&
                   'Log-spaced scan: d_theta: ',diff(is,1)
              
           else
              !Linear spacing
              diff(is,1)=((pi*(scan(is)%range_i)/180.)-(theta))/&
                   real(scan(is)%n_scan*scan(is)%n_res)

              if (writeOut) &
                   write(*,'(a,es15.6e3)')&
                   'Linear scan: d_theta: ',diff(is,1)
              
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
              diff(is,1)=(log10(scan(is)%range_f)-log10(ki))/&
                   real(scan(is)%n_scan*scan(is)%n_res)
              
              if (writeOut) &
                   write(*,'(a,es15.6e3)')&
                   'Log-spaced scan: d_|k|: ',diff(is,1)
              
           else
              !Linear spacing
              diff(is,1)=((scan(is)%range_f)-(ki))/&
                   real(scan(is)%n_scan*scan(is)%n_res)                      
              
              if (writeOut) &
                   write(*,'(a,es15.6e3)')&
                   'Linear scan: d_|k|: ',diff(is,1)
              
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
              ! eigen + heat + low_n  (MATCH single: 22*nspec)
              write(fmt,'(a,i0,a)') '(6es15.6e3,12es15.6e3,', 22*nspec, 'es15.6e3)'
           else
              ! eigen + heat (non-low_n)
              write(fmt,'(a,i0,a)') '(6es15.6e3,12es15.6e3,', 15*nspec, 'es15.6e3,es15.6e3)'
           endif
           out_type = 0
        else
           ! eigen only
           write(fmt,'(a,i0,a)') '(6es15.6e3,12es15.6e3,', 14*nspec, 'es15.6e3)'
           out_type = 1
        endif
     else
        if (scan(2)%heat_s) then
           if (low_n) then
              ! heat only + low_n  (MATCH single: 14*nspec)
              write(fmt,'(a,i0,a)') '(6es15.6e3,', 14*nspec, 'es15.6e3)'
           else
              ! heat only (non-low_n)  (remove extra trailing es field)
              write(fmt,'(a,i0,a)') '(6es15.6e3,', 7*nspec, 'es15.6e3)'
           endif
           out_type = 2
        else
          ! frequency only
          write(fmt,'(a,i0,a)') '(6es15.6e3,', 6*nspec, 'es15.6e3)'
          out_type = 3
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
!Our method for selecting different dispersion relations
complex function disp(om)
  use vars, only : collision_type, sync_solver_state, has_cold_species, dispersion_status
  use vars, only : has_cold_reference, validate_cold_reference
  
  complex :: om
  !!Complex Frequency

  if (has_cold_reference()) then
     call validate_cold_reference
     disp = disp_mixed(om)
     return
  endif
  call sync_solver_state
  dispersion_status = 0
  if (has_cold_species()) then
     disp = disp_mixed(om)
     return
  endif
  select case (collision_type)

  case (0)
     !! Hot, collisionless dispersion relation
     disp = disp0(om)
  case (1)
     !! Hot plasma dispersion relation
     !! with weak neutral-charged collisions
     disp = disp1(om)
  case default
     write(*,'(a)')'An unknown collision type was selected.'
     write(*,'(a)')'HALTING...'
     stop
  end select

  return
  
end function disp

complex function disp_mixed(om)
  !! Shared cold preflight/status handling; keep outside disp for legacy numerics.
  use plasma_z_backend, only: z_evaluation_guard,begin_z_evaluation,end_z_evaluation
  use vars, only : has_cold_species, all_species_cold, has_cold_reference, dispersion_status,z_backend
  use cold_reference, only : calc_warm_auxiliary
  use kinetic_core, only : guarded_parallel
  use vars, only : nspec, spec, kpar, kperp, vtp, susc, lam, low_n, susc_low
  use cold_susceptibility, only : calc_chi_cold, finite_complex, invalid_complex

  complex :: om
  !!Complex Frequency
  integer :: is, row, column
  complex :: cold_chi(3,3)

  type(z_evaluation_guard) :: z_guard
  logical :: z_reject

  call begin_z_evaluation(om,z_guard,disp_mixed,z_reject)
  if (z_reject) return
  dispersion_status = 0
  if (has_cold_species()) then
     ! Preflight all cold components BEFORE any warm thermal or wave expression.
     ! Inputs here are native d_ref; disp0 uses the derived warm rho_ref state.
     if (.not. finite_complex(om)) then
        dispersion_status = 1
     else if (om == cmplx(0.,0.)) then
        dispersion_status = 2
     else
        do is=1,nspec
           if (spec(is)%tau_s /= 0.) cycle
           call calc_chi_cold(spec(is),om,kpar,kperp,vtp,cold_chi,dispersion_status)
           susc(is,:,:) = cold_chi
           if (dispersion_status /= 0) exit
           if (low_n) susc_low(is,:,:,:) = invalid_complex()
        enddo
     endif
     if (dispersion_status /= 0) then
        susc = invalid_complex()
        lam = invalid_complex()
        if (low_n) susc_low = invalid_complex()
        disp_mixed = invalid_complex()
        call end_z_evaluation(z_guard,disp_mixed)
        return
     endif
  endif
  if (all_species_cold()) then
     disp_mixed = disp_all_cold(om)
  else if (kperp == 0.) then
     if (has_cold_reference()) then
        call guarded_parallel(om,disp_mixed,dispersion_status,disp_mixed_parallel)
     else
        disp_mixed = disp_mixed_parallel(om)
     endif
  else if (has_cold_reference()) then
     call calc_warm_auxiliary(spec,om,kpar,kperp,vtp,susc,dispersion_status,bessel,zet_in)
     if (dispersion_status==0) then
        if (low_n) susc_low=invalid_complex()
        disp_mixed=disp_all_cold(om) ! Shared native wave assembly, now including warm chi.
     else
        disp_mixed=invalid_complex()
     endif
  else
     disp_mixed = disp0(om)
  endif

  if (has_cold_species()) then
     if (dispersion_status==0 .and. .not. finite_complex(disp_mixed)) dispersion_status = 5
     if (dispersion_status==0) then
        do row=1,3
           do column=1,3
              if (.not. finite_complex(lam(row,column))) dispersion_status = 5
           enddo
        enddo
     endif
     if (dispersion_status /= 0) then
        susc = invalid_complex()
        lam = invalid_complex()
        if (low_n) susc_low = invalid_complex()
        disp_mixed = invalid_complex()
     endif
  endif
  call end_z_evaluation(z_guard,disp_mixed)
  return

end function disp_mixed

complex function disp_all_cold(om)
  !! Cold susceptibilities have passed disp_mixed preflight in native d_ref units.
  !! No beta, anisotropy, solver_* state or thermal functions are accessed here.
  use vars, only: kpar, kperp, vtp, susc, lam
  implicit none
  complex, intent(in) :: om
  complex :: nx, nz
  integer :: j
  lam=sum(susc,dim=1)
  do j=1,3
     lam(j,j)=lam(j,j)+1.
  enddo
  nx=kperp/(vtp*om)
  nz=kpar/(vtp*om)
  lam(1,1)=lam(1,1)-nz**2
  lam(2,2)=lam(2,2)-nx**2-nz**2
  lam(3,3)=lam(3,3)-nx**2
  lam(1,3)=lam(1,3)+nx*nz
  lam(3,1)=lam(3,1)+nx*nz
  ! Same Cartesian determinant convention as disp0; transpose, not conjugate.
  disp_all_cold=lam(1,1)*(lam(2,2)*lam(3,3)+lam(2,3)**2) &
       +2.*lam(1,2)*lam(2,3)*lam(1,3)-lam(1,3)**2*lam(2,2)+lam(1,2)**2*lam(3,3)
end function disp_all_cold

complex function disp_mixed_parallel(om)
  !! Exact kperp=0 limit of the same bi-Maxwellian harmonic sum, mixed runs only.
  !! Only n=+/-1 contributes transversely; n=0 gives the longitudinal response.
  !! Cold susceptibilities were evaluated by disp_mixed in native d_ref units.
  use vars, only: nspec, spec, kpar, vtp, susc, lam, low_n, susc_low
  use cold_susceptibility, only: invalid_complex
  implicit none
  complex, intent(in) :: om
  complex :: delta, zeta, response, nz2
  complex, parameter :: imaginary=(0.,1.)
  real :: gyro, plasma, thermal
  integer :: is, n, j

  do is=1,nspec
     if (spec(is)%tau_s == 0.) cycle
     thermal=sqrt(spec(is)%mu_s*spec(is)%tau_s)
     gyro=spec(is)%mu_s/spec(is)%Q_s
     plasma=spec(is)%D_s*spec(is)%mu_s/(spec(is)%Q_s**2*vtp**2)
     delta=om-kpar*spec(is)%vv_s
     susc(is,:,:)=cmplx(0.,0.)
     do n=-1,1,2
        zeta=(delta-n*gyro)/(kpar*thermal)
        response=spec(is)%alph_s-1. + &
             (spec(is)%alph_s*delta+(1.-spec(is)%alph_s)*n*gyro)*zet_in(kpar,zeta)/(kpar*thermal)
        susc(is,1,1)=susc(is,1,1)+plasma*response/(2.*om**2)
        susc(is,1,2)=susc(is,1,2)+imaginary*n*plasma*response/(2.*om**2)
     enddo
     susc(is,2,2)=susc(is,1,1)
     susc(is,2,1)=-susc(is,1,2)
     zeta=delta/(kpar*thermal)
     susc(is,3,3)=2.*plasma*(1.+zeta*zet_in(kpar,zeta))/(kpar*thermal)**2
  enddo
  if (low_n) susc_low=invalid_complex() ! Resonance heating is unavailable for mixtures.
  lam=sum(susc,dim=1)
  do j=1,3
     lam(j,j)=lam(j,j)+1.
  enddo
  nz2=(kpar/(vtp*om))**2
  lam(1,1)=lam(1,1)-nz2
  lam(2,2)=lam(2,2)-nz2
  disp_mixed_parallel=(lam(1,1)*lam(2,2)+lam(1,2)**2)*lam(3,3)
end function disp_mixed_parallel

!-=-=-=-=-=-=
!  The parameters upon which this dispersion relation depends are:
!   1) betap   = 8nT_parallel ref/B^2
!   2) kperp   = k_perp rho_ref
!   3) kpar    = k_parallel rho_ref (allows for +/- k_parallel values)
!   4) vtp     = v_t,ref,parallel/c
!
!  As well as the species dependent ratios, compared to the selected
!  Reference Species ref:
!   5) tau_s     = T_ref/T_s(both temperatures parallel)
!   6) mu_s      = m_ref/m_s
!   7) alph_s    = T_perp/T_par|_s
!   8) Q_s       = q_ref/q_s
!   9) D_s       = n_s/n_ref
!   10)vv_s      = v_drift s/v_Aref
!   11)nu_n,s    = nu_ns/Omega_ref
!
!   and returns om=omega/Omega_ref
!!!!
!NOTE: REFERENCE SPECIES ASSUMED TO BE FIRST SPECIES
!!!!
complex function disp1(om)
  use plasma_z_backend, only: z_evaluation_guard,begin_z_evaluation,end_z_evaluation
  !!Returns Plasma Dispersion Function for input complex frequency
  !! and global values of the dimensionless plasma parameters.
  !! updated to include weak neutral-charged collision
  use adaptive_kinetic, only : dynamic_harmonics, disp_adaptive
  use vars, only: betap=>solver_betap,kperp=>solver_kperp,kpar=>solver_kpar
  use vars, only: vtp=>solver_vtp,nspec,spec=>solver_spec,susc,lam
  use vars,  only: low_n, susc_low, pi
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

  real :: disp_nuns
  !! neutral-charged collision frequency
  !! \(nu_{n,s}/\Omega_{ref}\)
     
  !For computational efficiency, save Bessel functions between calls==========
  real, save :: kperp_last
  !! Previous value of \(k_{\perp} \rho_{ref} \)
     
  real, allocatable, dimension(:), save :: tau_last
  !!Previous \(T_{ref}/T_{s}|_{\parallel}\) value.

  real, allocatable, dimension(:), save :: alph_last
  !!Previous set of temperature anisotropies
  
  real, allocatable, dimension(:), save :: mu_last
  !!Previous set of mass ratios
  
  real, allocatable, dimension(:), save :: Q_last
  !!Previous set of charge ratios
     
  logical, allocatable, dimension(:) :: reuse_bessel
  !! Either use previous Bessel array, or recalculate,
  !! based on if argument for a particular species/component changes.


  type(z_evaluation_guard) :: z_guard
  logical :: z_reject

  call begin_z_evaluation(om,z_guard,disp1,z_reject)
  if (z_reject) return
  if (dynamic_harmonics()) then
     disp1=disp_adaptive(om,bessel,zet_in)
     call end_z_evaluation(z_guard,disp1)
     return
  endif
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
     disp_nuns =spec(is)%nu_ns
     
     include 'kinetic_collisional_species.inc'

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
     enddo !species loop
     

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
     disp1= a(1,1)*( a(2,2)*a(3,3) + a(2,3)**2. ) + &
          2.*a(1,2)*a(2,3)*a(1,3) - a(1,3)**2.*a(2,2) +a(1,2)**2.*a(3,3)

     !------------------------------------------------------------------------------
     !Save values of lambda_s components to save on Bessel function calculations
     kperp_last=kperp

     call end_z_evaluation(z_guard,disp1)

     return
end function disp1
  
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
   complex function disp0(om)
  use plasma_z_backend, only: z_evaluation_guard,begin_z_evaluation,end_z_evaluation
     !!Returns Plasma Dispersion Function for input complex frequency
     !! and global values of the dimensionless plasma parameters.
  use adaptive_kinetic, only : dynamic_harmonics, disp_adaptive
     use vars, only: betap=>solver_betap,kperp=>solver_kperp,kpar=>solver_kpar
     use vars, only: vtp=>solver_vtp,nspec,spec=>solver_spec,susc,lam
     use vars,  only: low_n, susc_low, pi
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


  type(z_evaluation_guard) :: z_guard
  logical :: z_reject

  call begin_z_evaluation(om,z_guard,disp0,z_reject)
  if (z_reject) return
     if (dynamic_harmonics()) then
       disp0=disp_adaptive(om,bessel,zet_in)
       call end_z_evaluation(z_guard,disp0)
       return
     endif
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
           if (spec(is)%tau_s == 0.) cycle
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
        if (spec(is)%tau_s == 0.) then
           ! Ordinary chi_s precomputed in disp(); no thermal conversion/Bessel/Z.
           eps_xx = eps_xx + susc(is,1,1)
           eps_xy = eps_xy + susc(is,1,2)
           eps_xz = eps_xz + susc(is,1,3)
           eps_yy = eps_yy + susc(is,2,2)
           eps_yz = eps_yz + susc(is,2,3)
           eps_zz = eps_zz + susc(is,3,3)
           cycle
        endif
        disp_tau  =spec(is)%tau_s
        disp_mu   =spec(is)%mu_s
        disp_alph =spec(is)%alph_s
        disp_Q    =spec(is)%Q_s
        disp_D    =spec(is)%D_s
        disp_v    =spec(is)%vv_s
        
        include 'kinetic_species.inc'

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
     disp0= a(1,1)*( a(2,2)*a(3,3) + a(2,3)**2. ) + &
          2.*a(1,2)*a(2,3)*a(1,3) - a(1,3)**2.*a(2,2) +a(1,2)**2.*a(3,3)

     !------------------------------------------------------------------------------
     !Save values of lambda_s components to save on Bessel function calculations
     kperp_last=kperp

     call end_z_evaluation(z_guard,disp0)

     return
   end function disp0

   !------------------------------------------------------------------------------
   !                           Greg Howes, 2005
   !------------------------------------------------------------------------------
   !     NOTE: This routine was adapted from f77 routine by Eliot Quataert
   complex function rtsec(func,x1,x2,xacc,iflag)
     !!Secant routine to find minimum value of dispersion relation function.
     use vars, only : has_cold_species, lam,z_backend
     use cold_susceptibility, only : finite_complex, invalid_complex
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
     if (has_cold_species() .or. z_backend=='faddeeva') then
        if (.not. finite_complex(fl)) then
           iflag=-1
           rtsec=invalid_complex()
           return
        endif
     endif
     f=func(x2)
     if (has_cold_species() .or. z_backend=='faddeeva') then
        if (.not. finite_complex(f)) then
           iflag=-1
           rtsec=invalid_complex()
           return
        endif
     endif
          
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

        if (has_cold_species() .or. z_backend=='faddeeva') then
           if (.not. finite_complex(f)) then
              iflag=-1
              rtsec=invalid_complex()
              return
           endif
        endif
        if (has_cold_species()) then
           if (abs(dx)<xacc .or. abs(f)==0.) then
              if (abs(f)/maxval(abs(lam))**3 < 1.e-8) return
           endif
        else
           if(abs(dx).lt.xacc.or. Abs(f) .eq. 0.)return
        endif
     enddo
     if (has_cold_species() .or. z_backend=='faddeeva') then
        iflag=-2
        rtsec=invalid_complex()
     endif
     return
   end function rtsec
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!                           Greg Howes, 2006
!------------------------------------------------------------------------------
   subroutine find_minima(val, numroots, iroots, nroots)
     !! Finding minimum values on complex frequency grid.
     use vars, only : ni, nr
     implicit none
     !Passed

     real,    dimension(:,:), pointer :: val
     !!Value of Dispersion relation.
     !!Covers 0:nr, 0:ni

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
          
     logical :: is_left, is_right, is_bottom, is_top, mid_r, mid_i
     !! Replace Parameter/constant logic that required fixed nr and ni.

     iroots = 0
     nroots = 0

     !Find local minima in map
     do ii = ni, 0, -1
        do ir = 0, nr
           !Real checks
           is_left   = (ir == 0)
           is_right  = (ir == nr)
           mid_r     = (ir > 0  .and. ir <  nr)
           !Imaginary checks
           is_bottom = (ii == 0)
           is_top    = (ii == ni)
           mid_i     = (ii > 0  .and. ii <  ni)
           
           ! --- handle left / middle / right in real index
           if (is_left) then
              if (val(ir,ii) < val(ir+1,ii)) then
                 ! bottom / middle / top in imaginary index
                 if (is_bottom) then
                    if (val(ir,ii) < val(ir,ii+1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 else if (mid_i) then
                    if (val(ir,ii) < val(ir,ii-1) .and. &
                         val(ir,ii) < val(ir,ii+1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 else if (is_top) then
                    if (val(ir,ii) < val(ir,ii-1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 end if
              end if
              
           else if (mid_r) then
              if (val(ir,ii) < val(ir-1,ii) .and. &
                   val(ir,ii) < val(ir+1,ii)) then
                 if (is_bottom) then
                    if (val(ir,ii) < val(ir,ii+1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 else if (mid_i) then
                    if (val(ir,ii) < val(ir,ii-1) .and. &
                         val(ir,ii) < val(ir,ii+1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 else if (is_top) then
                    if (val(ir,ii) < val(ir,ii-1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 end if
              end if
              
           else if (is_right) then
              if (val(ir,ii) < val(ir-1,ii)) then
                 if (is_bottom) then
                    if (val(ir,ii) < val(ir,ii+1)) then
                       nroots = nroots + 1
                       iroots(1,nroots) = ir
                       iroots(2,nroots) = ii
                    end if
                 else if (mid_i) then
                    if (val(ir,ii) < val(ir,ii-1) .and. &
                         val(ir,ii) < val(ir,ii+1)) then
                       nroots = nroots + 1
                    iroots(1,nroots) = ir
                    iroots(2,nroots) = ii
                 end if
              else if (is_top) then
                 if (val(ir,ii) < val(ir,ii-1)) then
                    nroots = nroots + 1
                    iroots(1,nroots) = ir
                    iroots(2,nroots) = ii
                 end if
              end if
           end if
        end if

     end do
  end do
end subroutine find_minima

  
!-------------------------------------------------------------------------------
   function zet_in(kpar,zin)
     !!Evaluate Plasma Dispersion Function (Z) for kparallel >/< 0
     !! Default: binary64 Faddeeva Landau continuation; legacy is selectable.
     !! kpar<0 returns -Z(-zin); failed accurate calls return assigned quiet NaNs.
     !! Evaluation guards mark a nonfinite susceptibility/determinant as failed.
     use vars, only: z_backend
     use plasma_z_faddeeva, only: z_faddeeva_value
     implicit none
     !passed
     real :: kpar
     !! \( k_{\parallel} \rho_{ref} \)
     
     complex :: zin
     !! Argument of plasma dispersion function.

     complex :: zet_in
     !! Value of plasma dispersion function.
     
     if (z_backend=='faddeeva') then
        zet_in=z_faddeeva_value(kpar,zin)
        return
     endif
     if (kpar.gt.0.) zet_in=zetout(zin)
     if (kpar.lt.0.) zet_in=-zetout(-zin)
     return

   end function zet_in

!------------------------------------------------------------------------------
!Greg Howes, 2004


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
