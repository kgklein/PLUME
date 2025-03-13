!!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein (main author of PLUME)                                    !!
!!kris.klein@gmail.com                                                       !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!!*lin fpc routines                                                          !!
!!Collin Brown (author of FPC routine in PLUME)                              !!
!!collin.crbrown@gmail.com or collbrown@uiowa.edum                            !!
!!University of Iowa                                                         !!
!=============================================================================!
!=============================================================================!
module fpc
  !! Perturbed distribution function calculation from eigenfunctions and
  !! associated routines to calculate the FPC from it

  !CB Note: due to an effective sign difference between the textbooks by Stix/Swanson (used to compute eigen functions and fs1 respectively),
  !      we drop the minus sign in the correlation C_{E_i} = -q/2 v_i^2 d fs1/d v_i, so it is "q/2 v_i^2 d fs1/d v_i" when implemented...
  !      This minus sign difference is a result in Swanson 'suppressing the exp[i(k dot r - om t)]' (above 4.180) rather than taking the fourier transform of the product
  !      Note that the fourier transform of this term is F{exp[i(k dot r - om t)]} \proportional_to delta(kx`+kx) delta(ky`+ky) delta(kz`+kz) delta(om`-om)
  !      Effectively swaping the sign in kx ky kz (equivalently changing the frame). Typically, this does not matter as one can still get the same eigen functions, 
  !      but since we are correlating quantities in two diffferent frames, we must account for it somewhere... 
  !!!!!!!(This might be dated now as of mar 12 2025, TODO check this statement and maybe remove this comment....?)!!!!!!

  implicit none
  private :: calc_fs1, calc_exbar
  public :: compute_fpc_gyro, compute_fpc_cart

  real :: bs_last=0.0       
  !! Last Bessel function argument, used to store repeated calculations for efficiency
  real, allocatable :: jbess(:)  
  !! Regular (i.e. not the modified used elsewhere) Bessel functions, stored for efficiency 


  contains

    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    subroutine compute_fpc_cart(wrootindex)
      !! Computes the FPC associated with the linear fs1 and eigenfunction 
      !! response on a square cartesian grid and writes FPC to file
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv,nbesmax
      use vars, only : elecdircontribution
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              
      !! Index of the root to compute fs1 and FPC for

      character(1000) :: filename                      
      !! Output File name
      
      character(1000) :: outputPath                    
      !! Output folder

      character(1000) :: cmd
      !! Varaible to store command line commands
      
      real    :: vxi, vyi, vzi                       
      !! normalized velocity space current value in loop (note: vx, vy, vz corresponds to vperp1,vperp2,vpar, but we use 'x','y','z' as convention)
      
      integer :: vxindex, vyindex, vzindex
      !! loop counters
      
      real    :: vmax3rdval
      !! sampled range when computing projection
      
      complex :: omega
      !!Complex Frequency

      real    :: Cor_par_s, Cor_perp1_s, Cor_perp2_s 
      !! normalized correlation value TODO: remove these and other unused variables...
     
      real    :: wi,gi                               
      !! Freq and Damping of initial guess

      complex :: ominit                              
      !! Complex Frequency initial guess

      complex :: om1,om2                             
      !! Bracket Values
      
      integer :: iflag                               
      !! Flag for Root search

      real, parameter :: tol=1.0E-13                 
      !! Root Search Tolerance

      real, parameter :: prec=1.E-7                  
      !! Root Finding precision

      integer :: numstepvx, numstepvy, numstepvz     
      !! total number of steps in loop

      logical :: ex
      !! used to check if results directory exists
      
      complex, dimension(1:3)       :: ef, bf 
      !! E, B eigenfunction values (all 3 components)

      complex, dimension(1:nspec)     :: ns     
      !! density eigenfunction (all species)
      complex, dimension(1:3,1:nspec) :: Us     
      !! velocity eigenfunction (all species; all 3 componets per specie)

      !Heating (Required parameters of calc eigen)
      real, dimension(1:nspec) :: Ps 
      !! Power into/out of species
      
      real, dimension(1:4,1:nspec) :: Ps_split 
      !! Power into/out of species (Tensor that holds different channels (TTD, LD, CD))

      real, dimension(1:6,1:nspec) :: Ps_split_new 
      !! Power into/out of species updated by Greg G Howes to include off diagnal components
      
      real :: Ew 
      !! wave energy

      !loop counter/ loop parameters
      integer :: is                     
      !! species counter

      integer :: idir                   
      !! debug component contribution counter

      complex    :: tempux1val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempux2val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempux3val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempuy1val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempuy2val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempuy3val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempuz1val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempuz2val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: tempuz3val          
      !! debug temp val that holds susc tensor current contribution to renormalize

      complex    :: A1                     
      !! debug factor missing for fs1 "U" term (see Brown thesis appendx) that we compute via 'brute force'

      complex    :: B1                     
      !! debug factor missing for fs1 "W" term (see Brown thesis appendx) that we compute via 'brute force'

      integer :: unit_s                 
      !! out file unit counter

      real :: start, finish 
      !! debug/test to measure runtime of function

      !arrays that hold values on grid for efficiency
      real, allocatable, dimension(:) :: vvx,vvy,vvz  
      !! Velocity grid values (norm: w_par_s)

      integer :: ivx,ivy,ivz                    
      !! Index for each velocity dimension

      integer :: ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax  
      !! Index limits

      real,  allocatable, dimension(:,:,:,:) :: fs0   
      !! Dimensionless equilibrium fs0

      complex, allocatable, dimension(:,:,:,:) :: fs1 
      !! Perturbed Dist for all species

      real,  allocatable, dimension(:)  :: hatV_s
      !!Flow normalized to wpar_s

      real :: vperp                                   
      !! Perpendicular velocity
      
      real :: phi                                     
      !! Gyrophase angle (azimuthal)
      
      complex, allocatable, dimension(:,:,:,:) :: dfs1dvx,dfs1dvy,dfs1dvz 
      !! Derivatives
      
      real, allocatable, dimension(:,:,:,:) :: corex,corey,corez 
      !! 3V Correlations
      
      real, allocatable, dimension(:,:,:) :: corex_xy,corex_xz,corex_zy 
      !! 2V Corrs
      
      real, allocatable, dimension(:,:,:) :: corey_xy,corey_xz,corey_zy 
      !! 2V Corrs
      
      real, allocatable, dimension(:,:,:) :: corez_xy,corez_xz,corez_zy 
      !! 2V Corrs
      
      complex, allocatable, dimension(:,:,:) :: fs1_xy,fs1_xz,fs1_zy 
      !! 2V fs1 
      
      complex, allocatable, dimension(:) :: ns1 
      !! Density Fluctuation
      
      complex, allocatable, dimension(:,:) :: us1 
      !! Fluid Velocity Fluctuation
      
      real, allocatable, dimension(:) :: jxex,jyey,jzez 
      !! int_v 3V Correlations
      
      integer :: jj                                   
      !! Index
      
      character(100)  :: fmt                          
      !! Eigenfunction Output Format
      
      real :: delv3                                  
      !! delv^3
      
      real :: pi                                     
      !! 3.1415....
      
      character(100)  :: fmt_dbg1,fmt_dbg2           
      !! Eigenfunction Output Format

      complex    :: exbar               
      !! amplitude factor of fs1

      pi = 4.0*ATAN(1.0)
      A1 = 1. !Temporary scalar to fix coeff 
      B1 = 1. !Temporary scalar to fix coeff

      !check if results directory exists
      ! INQUIRE (DIRECTORY='data', EXIST=ex)
      ex = .true. !TODO: make this work for gfortran compiler
      if(ex) then
        write(*,*)"Assuming data folder already exists..."
      else
        write(*,*)"Creating data folder for output..."
        write(*,*)'mkdir data'
        CALL system('mkdir data')
        write(*,*)"Saving output to data folder..."
      endif

      write(outputPath,*) 'data/', trim(dataName) !!TODO: use more general pathing
      ! INQUIRE (DIRECTORY=trim(dataName), EXIST=ex)
      ex = .true.
      if(ex) then
        write(*,*)"assuming subfolder ", trim(dataName), " already exists"
      else
        write(*,*)"Creating data subfolder ", trim(dataName)
        write(cmd,*)'mkdir ',trim(outputPath)
        write(*,*)cmd
        CALL system(cmd)
        write(*,*)"Saving to data subfolder ", trim(dataName)
      endif

      !Grab dispersion relation solution
      wi = wroots(1,wrootindex)
      gi = wroots(2,wrootindex)
      ominit=cmplx(wi,gi)
      om1=ominit*(1.-prec)
      om2=ominit*(1.+prec)

      ! Refine Omega Value
      iflag=0
      omega=rtsec(disp,om1,om2,tol,iflag)
      
      call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,Ps_split_new,.true.,.true.)
      
      
      !==============================================================================
      ! Create Velocity grid and allocate variables
      !==============================================================================

      ! Create velocity grid variables: vvx,vvy,vvz
      ! NOTE: All (x,y,z) velocity values are normalized to species parallel thermal velocity
      ! NOTE: Code below assumes max.min values are integer multiples of delv (if not, throw warning flag!)
      ivxmin=int(vxmin/delv); if (real(ivxmin) .ne. vxmin/delv) write(*,*)'WARNING: vxmin not integer multiple of delv'
      ivxmax=int(vxmax/delv); if (real(ivxmax) .ne. vxmax/delv) write(*,*)'WARNING: vxmax not integer multiple of delv'
      ivymin=int(vymin/delv); if (real(ivymin) .ne. vymin/delv) write(*,*)'WARNING: vymin not integer multiple of delv'
      ivymax=int(vymax/delv); if (real(ivymax) .ne. vymax/delv) write(*,*)'WARNING: vymax not integer multiple of delv'
      ivzmin=int(vzmin/delv); if (real(ivzmin) .ne. vzmin/delv) write(*,*)'WARNING: vzmin not integer multiple of delv'
      ivzmax=int(vzmax/delv); if (real(ivzmax) .ne. vzmax/delv) write(*,*)'WARNING: vzmax not integer multiple of delv'

      !Allocate velocity grid variables
      allocate(vvx(ivxmin:ivxmax)); vvx(:)=0.
      allocate(vvy(ivymin:ivymax)); vvy(:)=0.
      allocate(vvz(ivzmin:ivzmax)); vvz(:)=0.
      !Populate velocity grid values
      do ivx=ivxmin,ivxmax
         vvx(ivx)=real(ivx)*delv
      enddo
      do ivy=ivymin,ivymax
         vvy(ivy)=real(ivy)*delv
      enddo
      do ivz=ivzmin,ivzmax
         vvz(ivz)=real(ivz)*delv
      enddo
      
      !Allocate fs1 and fs0 variables 
      allocate(fs0(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec))
      fs0(:,:,:,:)=0.
      allocate(fs1(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec))
      fs1(:,:,:,:)=0.

      !Allocate Bessel function values
      allocate(jbess(-nbesmax-1:nbesmax+1)); jbess(:)=0.
      
      !=============================================================================
      ! END Create Velocity grid and allocate variables
      !=============================================================================

      !=============================================================================
      ! Calculate fs0 and fs1 on (vx,vy,vz) grid
      !=============================================================================
      allocate(hatV_s(nspec))
      !Loop over (vx,vy,vz) grid and compute fs0 and fs1
      call calc_exbar(omega,ef,bf,exbar)
      do is = 1, nspec
         !Create variable for parallel flow velocity normalized to
         !       species parallel thermal speed
         hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(is)%mu_s*betap))
         do ivx=ivxmin,ivxmax
            do ivy=ivymin,ivymax
               vperp=sqrt(vvx(ivx)*vvx(ivx)+vvy(ivy)*vvy(ivy))
               do ivz=ivzmin,ivzmax
                  !Compute dimensionless equilibrium Distribution value, fs0
                  fs0(ivx,ivy,ivz,is)=fs0hat(vperp,vvz(ivz),hatV_s(is),spec(is)%alph_s)
                  !Compute perturbed  Distribution value, fs1
                  phi = ATAN2(vvy(ivy),vvx(ivx))
                  call calc_fs1(omega,vperp,vvz(ivz),phi,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,elecdircontribution,&
                                    A1,B1,exbar,fs0(ivx,ivy,ivz,is),fs1(ivx,ivy,ivz,is))
               enddo
            enddo
         enddo
      enddo
      !=============================================================================
      ! End Calculate fs0 and fs1 on (vx,vy,vz) grid
      !=============================================================================

      !=============================================================================
      ! Calculate velocity derivatives of fs1
      !=============================================================================
      allocate(dfs1dvx(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec))
      dfs1dvx=0.
      allocate(dfs1dvy(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec))
      dfs1dvy=0.
      allocate(dfs1dvz(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec))
      dfs1dvz=0.

      do is = 1, nspec
         !dfs1/dvx-------------------------------------------------------
         do ivz=ivzmin,ivzmax
            do ivy=ivymin,ivymax
               !First point: 1st order Finite difference
               dfs1dvx(ivxmin,ivy,ivz,is)=(fs1(ivxmin+1,ivy,ivz,is)-fs1(ivxmin,ivy,ivz,is))/(vvx(ivxmin+1)-vvx(ivxmin));
               !Loop: 2nd order Centered Finite difference
               do ivx=ivxmin+1,ivxmax-1
                  dfs1dvx(ivx,ivy,ivz,is)=(fs1(ivx+1,ivy,ivz,is)-fs1(ivx-1,ivy,ivz,is))/(vvx(ivx+1)-vvx(ivx-1));
               end do
               !Last point: 1st order Finite difference
               dfs1dvx(ivxmax,ivy,ivz,is)=(fs1(ivxmax,ivy,ivz,is)-fs1(ivxmax-1,ivy,ivz,is))/(vvx(ivxmax)-vvx(ivxmax-1));
            end do
         end do
         !dfs1/dvy-------------------------------------------------------
         do ivz=ivzmin,ivzmax
            do ivx=ivxmin,ivxmax
               !First point: 1st order Finite difference
               dfs1dvy(ivx,ivymin,ivz,is)=(fs1(ivx,ivymin+1,ivz,is)-fs1(ivx,ivymin,ivz,is))/(vvy(ivymin+1)-vvy(ivymin));
               !Loop: 2nd order Centered Finite difference
               do ivy=ivymin+1,ivymax-1
                  dfs1dvy(ivx,ivy,ivz,is)=(fs1(ivx,ivy+1,ivz,is)-fs1(ivx,ivy-1,ivz,is))/(vvy(ivy+1)-vvy(ivy-1));
               end do
               !Last point: 1st order Finite difference
               dfs1dvy(ivx,ivymax,ivz,is)=(fs1(ivx,ivymax,ivz,is)-fs1(ivx,ivymax-1,ivz,is))/(vvy(ivymax)-vvy(ivymax-1));
            end do
         end do
         !dfs1/dvz-------------------------------------------------------
         do ivy=ivymin,ivymax
            do ivx=ivxmin,ivxmax
               !First point: 1st order Finite difference
               dfs1dvz(ivx,ivy,ivzmin,is)=(fs1(ivx,ivy,ivzmin+1,is)-fs1(ivx,ivy,ivzmin,is))/(vvz(ivzmin+1)-vvz(ivzmin));
               !Loop: 2nd order Centered Finite difference
               do ivz=ivzmin+1,ivzmax-1
                  dfs1dvz(ivx,ivy,ivz,is)=(fs1(ivx,ivy,ivz+1,is)-fs1(ivx,ivy,ivz-1,is))/(vvz(ivz+1)-vvz(ivz-1));
               end do
               !Last point: 1st order Finite difference
               dfs1dvz(ivx,ivy,ivzmax,is)=(fs1(ivx,ivy,ivzmax,is)-fs1(ivx,ivy,ivzmax-1,is))/(vvz(ivzmax)-vvz(ivzmax-1));
            end do
         end do
      enddo
      !=============================================================================
      ! END Calculate velocity derivatives of fs1
      !=============================================================================

      !=============================================================================
      ! Calculate 3V CEx, CEy, CEz Correlations
      !=============================================================================
      allocate(corex(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec)); corex=0.
      allocate(corey(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec)); corey=0.
      allocate(corez(ivxmin:ivxmax,ivymin:ivymax,ivzmin:ivzmax,1:nspec)); corez=0.
      
      !Loop over (vx,vy,vz) grid and compute CEx, CEy, CEz Correlations
      do is = 1, nspec
         do ivx=ivxmin,ivxmax
            do ivy=ivymin,ivymax
               do ivz=ivzmin,ivzmax
                  ! CEx
                  corex(ivx,ivy,ivz,is)=-spec(is)%q_s*0.5*vvx(ivx)*vvx(ivx)* 0.5*(CONJG(dfs1dvx(ivx,ivy,ivz,is))*ef(1) &
                                       +dfs1dvx(ivx,ivy,ivz,is)*CONJG(ef(1)) )
                  
                  ! CEy
                  corey(ivx,ivy,ivz,is)=-spec(is)%q_s*0.5*vvy(ivy)*vvy(ivy)* 0.5*(CONJG(dfs1dvy(ivx,ivy,ivz,is))*ef(2) &
                                        +dfs1dvy(ivx,ivy,ivz,is)*CONJG(ef(2)) )
                  
                  ! CEz
                  corez(ivx,ivy,ivz,is)=-spec(is)%q_s*0.5*vvz(ivz)*vvz(ivz)* 0.5*(CONJG(dfs1dvz(ivx,ivy,ivz,is))*ef(3) &
                                        +dfs1dvz(ivx,ivy,ivz,is)*CONJG(ef(3)) )
                  
               enddo
            enddo
         enddo
      enddo

      !=============================================================================
      ! END Calculate 3V CEx, CEy, CEz Correlations
      !=============================================================================

      !=============================================================================
      ! Reduce to 2V fs1 and 2V CEx, CEy, CEz Correlations
      !=============================================================================
      !Perturbed Distribution------------------------------------------------
      allocate(fs1_xy(ivxmin:ivxmax,ivymin:ivymax,1:nspec)); fs1_xy=0.
      allocate(fs1_xz(ivxmin:ivxmax,ivzmin:ivzmax,1:nspec)); fs1_xz=0.
      allocate(fs1_zy(ivzmin:ivzmax,ivymin:ivymax,1:nspec)); fs1_zy=0.

      !Reduce Distribution
      fs1_xy(:,:,:)=sum(fs1(:,:,:,:),3)*delv
      fs1_xz(:,:,:)=sum(fs1(:,:,:,:),2)*delv
      fs1_zy(:,:,:)=sum(fs1(:,:,:,:),1)*delv !NOTE: Assumes nvy=nvz!
      
      !Transpose zy reductions to correct array ordering: (vz,vy) instead of (vy,vz)
      do is=1,nspec
         fs1_zy(:,:,is)=transpose(fs1_zy(:,:,is))
      enddo
      
      !Correlations---------------------------------------------------------
      allocate(corex_xy(ivxmin:ivxmax,ivymin:ivymax,1:nspec)); corex_xy=0.
      allocate(corex_xz(ivxmin:ivxmax,ivzmin:ivzmax,1:nspec)); corex_xz=0.
      allocate(corex_zy(ivzmin:ivzmax,ivymin:ivymax,1:nspec)); corex_zy=0.
      allocate(corey_xy(ivxmin:ivxmax,ivymin:ivymax,1:nspec)); corey_xy=0.
      allocate(corey_xz(ivxmin:ivxmax,ivzmin:ivzmax,1:nspec)); corey_xz=0.
      allocate(corey_zy(ivzmin:ivzmax,ivymin:ivymax,1:nspec)); corey_zy=0.
      allocate(corez_xy(ivxmin:ivxmax,ivymin:ivymax,1:nspec)); corez_xy=0.
      allocate(corez_xz(ivxmin:ivxmax,ivzmin:ivzmax,1:nspec)); corez_xz=0.
      allocate(corez_zy(ivzmin:ivzmax,ivymin:ivymax,1:nspec)); corez_zy=0.

      !Reduce Correlations
      corex_xy(:,:,:)=sum(corex(:,:,:,:),3)*delv
      corex_xz(:,:,:)=sum(corex(:,:,:,:),2)*delv
      corex_zy(:,:,:)=sum(corex(:,:,:,:),1)*delv !NOTE: Assumes nvy=nvz! TODO: remove this requirement (and fix comments everywhere)
      corey_xy(:,:,:)=sum(corey(:,:,:,:),3)*delv
      corey_xz(:,:,:)=sum(corey(:,:,:,:),2)*delv
      corey_zy(:,:,:)=sum(corey(:,:,:,:),1)*delv !NOTE: Assumes nvy=nvz!
      corez_xy(:,:,:)=sum(corez(:,:,:,:),3)*delv
      corez_xz(:,:,:)=sum(corez(:,:,:,:),2)*delv
      corez_zy(:,:,:)=sum(corez(:,:,:,:),1)*delv !NOTE: Assumes nvy=nvz!
      
      !Transpose zy reductions to correct array ordering: (vz,vy) instead of (vy,vz)
                 !NOTE: Assumes nvy=nvz!
      do is=1,nspec    
         corex_zy(:,:,is)=transpose(corex_zy(:,:,is))
         corey_zy(:,:,is)=transpose(corey_zy(:,:,is))
         corez_zy(:,:,is)=transpose(corez_zy(:,:,is))
      enddo
      !=============================================================================
      ! END Reduce to 2V CEx, CEy, CEz Correlations
      !=============================================================================

      !=============================================================================
      ! Integrate Distribution and Correlations over velocity 
      !=============================================================================

      delv3=delv*delv*delv
      allocate(ns1(nspec)); ns1=0.
      allocate(us1(3,nspec)); us1=0.

      ! if(1==0) then !Force renormalization
      !    !Compute B1
      !    do is = 1, nspec
      !       tempux1val = 0.
      !       tempuy1val = 0.
      !       tempuz1val = 0.
      !       tempux2val = 0.
      !       tempuy2val = 0.
      !       tempuz2val = 0.
      !       tempux3val = 0.
      !       tempuy3val = 0.
      !       tempuz3val = 0.
      !       do idir =1, 3
      !          hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(is)%mu_s*betap))
      !          do ivx=ivxmin,ivxmax
      !             do ivy=ivymin,ivymax
      !                vperp=sqrt(vvx(ivx)*vvx(ivx)+vvy(ivy)*vvy(ivy))
      !                do ivz=ivzmin,ivzmax
      !                   !Compute dimensionless equilibrium Distribution value, fs0
      !                   fs0(ivx,ivy,ivz,is)=fs0hat(vperp,vvz(ivz),hatV_s(is),spec(is)%alph_s)
      !                   !Compute perturbed  Distribution value, fs1
      !                   phi = ATAN2(vvy(ivy),vvx(ivx))
      !                   call calc_fs1(omega,vperp,vvz(ivz),phi,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
      !                                     spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,&
      !                                     real(idir),A1,B1,fs0(ivx,ivy,ivz,is),fs1(ivx,ivy,ivz,is))
      !                enddo
      !             enddo
      !          enddo
      !          if(idir == 1) then
      !             do ivx=ivxmin,ivxmax
      !                tempux1val=tempux1val+vvx(ivx)*sum(sum(fs1(ivx,:,:,is),2),1)*delv3
      !             enddo
      !             do ivy=ivymin,ivymax
      !                tempuy1val=tempuy1val+vvy(ivy)*sum(sum(fs1(:,ivy,:,is),2),1)*delv3
      !             enddo
      !             do ivz=ivzmin,ivzmax
      !                tempuz1val=tempuz1val+vvz(ivz)*sum(sum(fs1(:,:,ivz,is),2),1)*delv3
      !             enddo
      !             tempux1val = tempux1val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !             tempuy1val = tempuy1val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !             tempuz1val = tempuz1val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !          end if

      !          if(idir == 2) then
      !             do ivx=ivxmin,ivxmax
      !                tempux2val=tempux2val+vvx(ivx)*sum(sum(fs1(ivx,:,:,is),2),1)*delv3
      !             enddo
      !             do ivy=ivymin,ivymax
      !                tempuy2val=tempuy2val+vvy(ivy)*sum(sum(fs1(:,ivy,:,is),2),1)*delv3
      !             enddo
      !             do ivz=ivzmin,ivzmax
      !                tempuz2val=tempuz2val+vvz(ivz)*sum(sum(fs1(:,:,ivz,is),2),1)*delv3
      !             enddo
      !             tempux2val = tempux2val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !             tempuy2val = tempuy2val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !             tempuz1val = tempuz1val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !          endif

      !          if(idir == 3) then
      !             do ivx=ivxmin,ivxmax
      !                tempux3val=tempux3val+vvx(ivx)*sum(sum(fs1(ivx,:,:,is),2),1)*delv3
      !             enddo
      !             do ivy=ivymin,ivymax
      !                tempuy3val=tempuy3val+vvy(ivy)*sum(sum(fs1(:,ivy,:,is),2),1)*delv3
      !             enddo
      !             do ivz=ivzmin,ivzmax
      !                tempuz3val=tempuz3val+vvz(ivz)*sum(sum(fs1(:,:,ivz,is),2),1)*delv3
      !             enddo
      !             tempux3val = tempux3val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !             tempuy3val = tempuy3val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !             tempuz1val = tempuz1val*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s !norm to alfven speed
      !          endif

      !       enddo

      !       !try A+B+A
      !       ! A1 = (ABS(Us(2,is))*ABS(tempuz2val)-ABS(Us(3,is))*ABS(tempuy2val)) / ((ABS(tempuy1val)+ABS(tempuy3val))*ABS(tempuz2val)-(ABS(tempuz1val)+ABS(tempuz3val))*ABS(tempuy2val))
      !       ! B1 = (ABS(Us(3,is))*(ABS(tempuy1val)+ABS(tempuy3val))-ABS(Us(2,is))*(ABS(tempuz1val)+ABS(tempuz3val))) / ((ABS(tempuy1val)+ABS(tempuy3val))*ABS(tempuz2val)-(ABS(tempuz1val)+ABS(tempuz3val))*ABS(tempuy2val))

            
      !       ! A1 = (ABS(Us(2,is))*ABS(tempuz3val)-ABS(Us(3,is))*ABS(tempuy3val)) / ((ABS(tempuy1val)+ABS(tempuy2val))*ABS(tempuz3val)-(ABS(tempuz1val)+ABS(tempuz2val))*ABS(tempuy3val))
      !       ! B1 = (ABS(Us(3,is))*(ABS(tempuy1val)+ABS(tempuy2val))-ABS(Us(2,is))*(ABS(tempuz1val)+ABS(tempuz2val))) / ((ABS(tempuy1val)+ABS(tempuy2val))*ABS(tempuz3val)-(ABS(tempuz1val)+ABS(tempuz2val))*ABS(tempuy3val))



      !       ! A1 = (ABS(Us(1,is))*ABS(tempuz3val)-ABS(Us(3,is))*ABS(tempux3val)) / ((ABS(tempux1val)+ABS(tempux2val))*ABS(tempuz3val)-(ABS(tempuz1val)+ABS(tempuz2val))*ABS(tempux3val))
      !       ! B1 = (ABS(Us(3,is))*(ABS(tempux1val)+ABS(tempux2val))-ABS(Us(1,is))*(ABS(tempuz1val)+ABS(tempuz2val))) / ((ABS(tempux1val)+ABS(tempux2val))*ABS(tempuz3val)-(ABS(tempuz1val)+ABS(tempuz2val))*ABS(tempux3val))

      !       ! A1 = ABS(A1)
      !       ! B1 = ABS(B1)

      !       A1 = (Us(1,is)*tempuz3val-Us(3,is)*tempux3val) / (((tempux1val)+(tempux2val))*(tempuz3val)-((tempuz1val)+(tempuz2val))*(tempux3val))
      !       B1 = (Us(3,is)*(tempux1val+(tempux2val))-(Us(1,is))*((tempuz1val)+(tempuz2val))) / (((tempux1val)+(tempux2val))*(tempuz3val)-((tempuz1val)+(tempuz2val))*(tempux3val))


      !       write(*,*)'tempuxvals',tempux1val,tempux2val,tempux3val
      !       write(*,*)'tempuyvals',tempuy1val,tempuy2val,tempuy3val
      !       write(*,*)'tempuzvals',tempuz1val,tempuz2val,tempuz3val
      !       write(*,*)'Debug A1',is,A1, ((tempux1val+tempux2val)*tempuz3val-(tempuz1val+tempuz2val)*tempux3val)
      !       write(*,*)'Debug B1',is,B1, ((tempux1val+tempux2val)*tempuz3val-(tempuz1val+tempuz2val)*tempux3val)

      !       !renormalize values
      !       !renormalize ux_spec
      !       A1 = 1. !hacky way to turn off renormalization before I clean things up.
      !       B1 = 1.
      !       us1(1,is)=A1*tempux1val+A1*tempux2val+B1*tempux3val
      !       us1(2,is)=A1*tempuy1val+A1*tempuy2val+B1*tempuy3val
      !       us1(3,is)=A1*tempuz1val+A1*tempuz2val+B1*tempuz3val


      !       if (elecdircontribution == 1) then
      !          us1(1,is)=A1*tempux1val
      !          us1(2,is)=A1*tempuy1val
      !          us1(3,is)=A1*tempuz1val
      !       else if (elecdircontribution == 2)then
      !          us1(1,is)=A1*tempux2val
      !          us1(2,is)=A1*tempuy2val
      !          us1(3,is)=A1*tempuz2val
      !       else if (elecdircontribution == 3)then
      !          us1(1,is)=B1*tempux3val
      !          us1(2,is)=B1*tempuy3val
      !          us1(3,is)=B1*tempuz3val
      !       else
      !          us1(1,is)=A1*tempux1val+A1*tempux2val+B1*tempux3val
      !          us1(2,is)=A1*tempuy1val+A1*tempuy2val+B1*tempuy3val
      !          us1(3,is)=A1*tempuz1val+A1*tempuz2val+B1*tempuz3val
      !       end if
      !    enddo

      !    !******TODO: renormalize correlation*****
      ! end if

      call calc_exbar(omega,ef,bf,exbar)
      do is = 1, nspec
         !Create variable for parallel flow velocity normalized to
         !       species parallel thermal speed
         hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(is)%mu_s*betap))
         do ivx=ivxmin,ivxmax
            do ivy=ivymin,ivymax
               vperp=sqrt(vvx(ivx)*vvx(ivx)+vvy(ivy)*vvy(ivy))
               do ivz=ivzmin,ivzmax
                  !Compute dimensionless equilibrium Distribution value, fs0
                  fs0(ivx,ivy,ivz,is)=fs0hat(vperp,vvz(ivz),hatV_s(is),spec(is)%alph_s)
                  !Compute perturbed  Distribution value, fs1
                  phi = ATAN2(vvy(ivy),vvx(ivx))
                  call calc_fs1(omega,vperp,vvz(ivz),phi,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,elecdircontribution,A1,B1,exbar,fs0(ivx,ivy,ivz,is),fs1(ivx,ivy,ivz,is))
               enddo
            enddo
         enddo
      enddo



      !Integrate 0th moment
      do is = 1, nspec
      ! Density Fluctuation: Zeroth Moment of delta f
         ns1(is)=sum(sum(sum(fs1(:,:,:,is),3),2),1)*delv3

         !Correct Normalization to n_0R
         ns1(is)=ns1(is)*sqrt(spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))*spec(is)%D_s !TODO: this is not correct and should have an additional factor of Ex/B0 in it (compute with calc exbar rountine)
      

      ! Fluid Velocity: First Moment of total f = delta f (since int v f_0=0)
         !x-component
         do ivx=ivxmin,ivxmax
            us1(1,is)=us1(1,is)+vvx(ivx)*sum(sum(fs1(ivx,:,:,is),2),1)*delv3 !TODO: this is not correct and should have an additional factor of Ex/B0 in it (compute with calc exbar rountine)
         enddo
         !y-component
         do ivy=ivymin,ivymax
            us1(2,is)=us1(2,is)+vvy(ivy)*sum(sum(fs1(:,ivy,:,is),2),1)*delv3 !TODO: this is not correct and should have an additional factor of Ex/B0 in it (compute with calc exbar rountine)
         enddo
         !z-component
         do ivz=ivzmin,ivzmax
            us1(3,is)=us1(3,is)+vvz(ivz)*sum(sum(fs1(:,:,ivz,is),2),1)*delv3 !TODO: this is not correct and should have an additional factor of Ex/B0 in it (compute with calc exbar rountine)
         enddo
      enddo

      !Integrate Correlations
      allocate(jxex(nspec)); jxex=0.
      allocate(jyey(nspec)); jyey=0.
      allocate(jzez(nspec)); jzez=0.
      do is = 1, nspec
         jxex(is)=sum(sum(sum(corex(:,:,:,is),3),2),1)*delv3 !int_v CEx= jxEx
         jyey(is)=sum(sum(sum(corey(:,:,:,is),3),2),1)*delv3 !int_v CEy= jyEy    
         jzez(is)=sum(sum(sum(corez(:,:,:,is),3),2),1)*delv3 !int_v CEz= jzEz
      enddo
 
      !=============================================================================
      ! END Integrate Distribution and Correlations over velocity
      !=============================================================================

      !=============================================================================
      ! Output Distributions and Correlations to Files
      !=============================================================================
      do is = 1, nspec
        call check_nbesmax(MAX(ABS(vxmin),ABS(vxmax),ABS(vymin),ABS(vymax),ABS(vzmin),ABS(vzmax)),&
          spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s)

        !make file to store result
        !TODO: used "get unused unit" to get unit_s to pick correct 'number' to write to
        unit_s = 12+5*is !note: unit = 5,6 are reserved by standard fortran for input form keyboard/ writing to screen
        !TODO: fix formating in the write statements here and in gyro...
        !write(*,'(5A,A,1A,A,16A,I0.2,5A,I0.2)')&
        !'data/',trim(dataName),'/',trim(outputName),'.cparcart.specie',(is),'.mode',wrootindex
        write(filename,'(5A,I0.2,1A,I0.2)')&
        'data/',trim(dataName),'/',trim(outputName),'.cparcart.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating (cart is for cartesian)
        open(unit=unit_s,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cperp1.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+1,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cperp2.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+2,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.dfs.real.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+3,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.dfs.imag.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+4,file=trim(filename),status='replace')

        write(*,*)'Calculating fpc for species ',is
        write(*,*)'Writing omega/kpar V_a normalization to file...'

        !Write header information to output all 4 files
        do jj=0,4
           write(unit_s+jj,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
           write(unit_s+jj,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
           write(unit_s+jj, '(7a22)')'vxmin','vxmax','vymin','vymax','vzmin','vzmax','delv'
           write(unit_s+jj, '(7es22.7)')vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
           write(unit_s+jj, *) '-------------'
        enddo
        
        !CEi(vx,vy)--------------------------------------------------------------
        do ivx=ivxmin,ivxmax
           do ivy=ivymin,ivymax
              if(ABS(corez_xy(ivx,ivy,is)) .lt. 9.999E-99) then
                write(unit_s,'(es17.5)',advance='no') 0.
              else
                write(unit_s,'(es17.5)',advance='no') corez_xy(ivx,ivy,is)
              endif
              if(ABS(corex_xy(ivx,ivy,is)) .lt. 9.999E-99) then
                write(unit_s+1,'(es17.5)',advance='no') 0.
              else
                write(unit_s+1,'(es17.5)',advance='no') corex_xy(ivx,ivy,is)
              endif
              if(ABS(corey_xy(ivx,ivy,is)) .lt. 9.999E-99) then
                write(unit_s+2,'(es17.5)',advance='no') 0. 
              else
                write(unit_s+2,'(es17.5)',advance='no') corey_xy(ivx,ivy,is)
              endif
              if(ABS(real(fs1_xy(ivx,ivy,is))) .lt. 9.999E-99) then
                write(unit_s+3,'(2es17.5)',advance='no') 0. 
              else
                write(unit_s+3,'(2es17.5)',advance='no') real(fs1_xy(ivx,ivy,is))
              endif
              if(ABS(aimag(fs1_xy(ivx,ivy,is))) .lt. 9.999E-99) then
                write(unit_s+4,'(2es17.5)',advance='no') 0. 
              else
                write(unit_s+4,'(2es17.5)',advance='no') aimag(fs1_xy(ivx,ivy,is))
              endif
           enddo
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
          write(unit_s+4,*)
        enddo
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'
        write(unit_s+3,*)'---'
        write(unit_s+4,*)'---'
        
        !CEi(vx,vz)--------------------------------------------------------------
         do ivx=ivxmin,ivxmax
            do ivz=ivzmin,ivzmax
              if(ABS(corez_xz(ivx,ivz,is)) .lt. 9.999E-99) then
                write(unit_s,'(es17.5)',advance='no') 0.
              else
                write(unit_s,'(es17.5)',advance='no') corez_xz(ivx,ivz,is)
              endif
              if(ABS(corex_xz(ivx,ivz,is)) .lt. 9.999E-99) then
                write(unit_s+1,'(es17.5)',advance='no') 0.
              else
                write(unit_s+1,'(es17.5)',advance='no') corex_xz(ivx,ivz,is)
              endif
              if(ABS(corey_xz(ivx,ivz,is)) .lt. 9.999E-99) then
                write(unit_s+2,'(es17.5)',advance='no') 0.
              else
                write(unit_s+2,'(es17.5)',advance='no') corey_xz(ivx,ivz,is)
              endif
              if(ABS(real(fs1_xz(ivx,ivz,is))) .lt. 9.999E-99) then
                write(unit_s+3,'(2es17.5)',advance='no') 0.
              else
                write(unit_s+3,'(2es17.5)',advance='no') real(fs1_xz(ivx,ivz,is))
              endif
              if(ABS(aimag(fs1_xz(ivx,ivz,is))) .lt. 9.999E-99) then
                write(unit_s+4,'(2es17.5)',advance='no') 0.
              else
                write(unit_s+4,'(2es17.5)',advance='no') aimag(fs1_xz(ivx,ivz,is))
              endif
           enddo
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
          write(unit_s+4,*)
        enddo
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'
        write(unit_s+3,*)'---'
        write(unit_s+4,*)'---'

        !CEi(vy,vz)-------------------------------------------------------------
        do ivz=ivzmin,ivzmax
           do ivy=ivymin,ivymax
              if(ABS(corez_zy(ivz,ivy,is)) .lt. 9.999E-99) then
                write(unit_s,'(es17.5)',advance='no') 0.
              else
                write(unit_s,'(es17.5)',advance='no') corez_zy(ivz,ivy,is)
              endif
              if(ABS(corex_zy(ivz,ivy,is)) .lt. 9.999E-99) then
                write(unit_s+1,'(es17.5)',advance='no') 0.
              else
                write(unit_s+1,'(es17.5)',advance='no') corex_zy(ivz,ivy,is)
              endif
              if(ABS(corey_zy(ivz,ivy,is)) .lt. 9.999E-99) then
                write(unit_s+2,'(es17.5)',advance='no') 0.
              else
                write(unit_s+2,'(es17.5)',advance='no') corey_zy(ivz,ivy,is)
              endif
              if(ABS(real(fs1_zy(ivz,ivy,is))) .lt. 9.999E-99) then
                write(unit_s+3,'(2es17.5)',advance='no') 0.
              else
                write(unit_s+3,'(2es17.5)',advance='no') real(fs1_zy(ivz,ivy,is))
              endif
              if(ABS(aimag(fs1_zy(ivz,ivy,is))) .lt. 9.999E-99) then
                write(unit_s+4,'(2es17.5)',advance='no') 0.
              else
                write(unit_s+4,'(2es17.5)',advance='no') aimag(fs1_zy(ivz,ivy,is))
              endif
           enddo
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
          write(unit_s+4,*)
        enddo
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'
        write(unit_s+3,*)'---'
        write(unit_s+4,*)'---'
      end do

      close(unit_s)
      close(unit_s+1)
      close(unit_s+2)
      close(unit_s+3)
      close(unit_s+4)

      !Write Complex Eigenfunction-------------------------------------------
      !The below output is for debugging and is not intended to be used by end users...
      write(filename,'(5A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.eigen.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
      open(unit=unit_s+5,file=trim(filename),status='replace')

      !Write format (consistent with usual PLUME output)
      write(fmt,'(a,i0,a)')'(6es15.6,12es15.6,',15*nspec,'es15.6)'
      write(unit_s+5,fmt)&
           kperp,kpar,betap,vtp,&
           omega,&            
           bf(1:3),ef(1:3),Us(1:3,1:nspec),ns(1:nspec),&
           Ps(1:nspec),Ps_split_new(1:6,1:nspec) !,params(1:6,1:nspec)
      close(unit_s+5)

      !Write Velocity Integrated Moments-------------------------------------------
      !The below output is for debugging and is not intended to be used by end users...
      write(filename,'(5A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.mom.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
      open(unit=unit_s+5,file=trim(filename),status='replace')

      !Write format 
      !If really small (less than 10^-99 in magnitude)- round down to zero for formatting
      do is = 1, nspec
         if(ABS(ns1(is)) .lt. 9.999E-99) then
            ns1(is) = 0.
         endif
         if(ABS(us1(1,is)) .lt. 9.999E-99) then
            us1(1,is) = 0.
         endif
         if(ABS(us1(2,is)) .lt. 9.999E-99) then
            us1(2,is) = 0.
         endif 
         if(ABS(us1(3,is)) .lt. 9.999E-99) then
            us1(3,is) = 0.
         endif 
         if(ABS(jxex(is)) .lt. 9.999E-99) then
            jxex(is) = 0.
         endif
         if(ABS(jyey(is)) .lt. 9.999E-99) then
            jyey(is) = 0.
         endif
         if(ABS(jzez(is)) .lt. 9.999E-99) then
            jzez(is) = 0.
         endif
      enddo

      write(fmt,'(a,i0,a)')'(',11*nspec,'es15.6)'
      write(unit_s+5,fmt)&
           ns1(1:nspec),us1(:,1:nspec),jxex(1:nspec),jyey(1:nspec),jzez(1:nspec)
      close(unit_s+5)
     
      !=============================================================================
      ! END Output Distributions and Correlations to Files
      !=============================================================================

      !=============================================================================
      ! DEBUG: Moment comparisons
      !=============================================================================
      if (1 == 1) then !DEBUG
         !Density
         write(*,'(a,2es12.2,5X,f6.3)')'nir=  ',real(ns(1)),real(ns1(1)), &
              real(ns(1))/real(ns1(1))
         write(*,'(a,2es12.2,5X,f6.3)')'nii=  ',aimag(ns(1)),aimag(ns1(1)), &
              aimag(ns(1))/aimag(ns1(1))
         write(*,'(a,2es12.2,5X,f6.3)')'ner=  ',real(ns(2)),real(ns1(2)), &
              real(ns(2))/real(ns1(2))
         write(*,'(a,2es12.2,5X,f6.3)')'nei=  ',aimag(ns(2)),aimag(ns1(2)), &
              aimag(ns(2))/aimag(ns1(2))
         write(*,*)
         
         !Fluid Moment
         fmt_dbg1='(a,2es12.2,5X,es12.2,5X,2es12.2)'
         fmt_dbg2='(a,2es12.2,5X,es12.2,5X,2f12.3)'
         write(*,fmt_dbg1)'uxir= ',real(Us(1,1)),real(us1(1,1)), &
              real(Us(1,1))/real(us1(1,1)), abs(Us(1,1)), abs(us1(1,1))
         write(*,fmt_dbg2)'uxii= ',aimag(Us(1,1)),aimag(us1(1,1)), &
              aimag(Us(1,1))/aimag(us1(1,1)),atan2(aimag(Us(1,1)),real(Us(1,1)))/pi, &
              atan2(aimag(us1(1,1)),real(us1(1,1)))/pi
         write(*,fmt_dbg1)'uyir= ',real(Us(2,1)),real(us1(2,1)), &
              real(Us(2,1))/real(us1(2,1)), abs(Us(2,1)), abs(us1(2,1))
         write(*,fmt_dbg2)'uyii= ',aimag(Us(2,1)),aimag(us1(2,1)), &
              aimag(Us(2,1))/aimag(us1(2,1)),atan2(aimag(Us(2,1)),real(Us(2,1)))/pi, &
              atan2(aimag(us1(2,1)),real(us1(2,1)))/pi
         write(*,fmt_dbg1)'uzir= ',real(Us(3,1)),real(us1(3,1)), &
              real(Us(3,1))/real(us1(3,1)), abs(Us(3,1)), abs(us1(3,1))
         write(*,fmt_dbg2)'uzii= ',aimag(Us(3,1)),aimag(us1(3,1)), &
              aimag(Us(3,1))/aimag(us1(3,1)),atan2(aimag(Us(3,1)),real(Us(3,1)))/pi, &
              atan2(aimag(us1(3,1)),real(us1(3,1)))/pi
         write(*,*)
         write(*,fmt_dbg1)'uxer= ',real(Us(1,2)),real(us1(1,2)), &
              real(Us(1,2))/real(us1(1,2)), abs(Us(1,2)), abs(us1(1,2))
         write(*,fmt_dbg2)'uxei= ',aimag(Us(1,2)),aimag(us1(1,2)), &
              aimag(Us(1,2))/aimag(us1(1,2)),atan2(aimag(Us(1,2)),real(Us(1,2)))/pi, &
              atan2(aimag(us1(1,2)),real(us1(1,2)))/pi
         write(*,fmt_dbg2)'uyei= ',aimag(Us(2,2)),aimag(us1(2,2)), &
              aimag(Us(2,2))/aimag(us1(2,2)),atan2(aimag(Us(2,2)),real(Us(2,2)))/pi, &
              atan2(aimag(us1(2,2)),real(us1(2,2)))/pi
         write(*,fmt_dbg1)'uzer= ',real(Us(3,2)),real(us1(3,2)), &
              real(Us(3,2))/real(us1(3,2)), abs(Us(3,2)), abs(us1(3,2))
         write(*,fmt_dbg2)'uzei= ',aimag(Us(3,2)),aimag(us1(3,2)), &
              aimag(Us(3,2))/aimag(us1(3,2)),atan2(aimag(Us(3,2)),real(Us(3,2)))/pi, &
              atan2(aimag(us1(3,2)),real(us1(3,2)))/pi
      endif
      !=============================================================================
      ! END DEBUG: Moment comparisons
      !=============================================================================

      !Deallocate variables
      deallocate(ns1,us1,jxex,jyey,jzez)
      deallocate(corex_xy,corex_xz,corex_zy)
      deallocate(corey_xy,corey_xz,corey_zy)
      deallocate(corez_xy,corez_xz,corez_zy)
      deallocate(corex,corey,corez)
      deallocate(dfs1dvx,dfs1dvy,dfs1dvz)
      deallocate(fs0,fs1)
      deallocate(hatV_s)
      deallocate(vvx,vvy,vvz)

    end subroutine compute_fpc_cart


    !------------------------------------------------------------------------------
    !                           Collin Brown, 2020
    !------------------------------------------------------------------------------
    subroutine compute_fpc_gyro(wrootindex)
      !! Computes the FPC associated with the linear fs1 and eigenfunction 
      !! response on a square cartesian grid and writes fs1 and FPC to file

      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv,nbesmax
      use vars, only : elecdircontribution
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              
      !! index of selected root

      character(1000) :: filename                      
      !! Output File name
      
      character(1000) :: outputPath                    
      !! Output folder
      
      character(1000) :: cmd                           
      !! Varaible to store command line commands
      
      real    :: vperpi, vpari                       
      !! normalized velocity space current value in loop
      
      integer :: vperpindex, vparindex               
      !! loop counters
      
      complex :: omega                               
      !! Complex Frequency
      
      real    :: Cor_par_s, Cor_perp_s                
      !! normalized correlation value
      
      real    :: wi,gi                               
      !! Freq and Damping of initial guess
      
      complex :: ominit                              
      !! Complex Frequency initial guess
      
      complex :: om1,om2                             
      !! Bracket Values
      
      integer :: iflag                               
      !! Flag for Root search
      
      real, parameter :: tol=1.0E-13                 
      !! Root Search Tolerance
      
      real, parameter :: prec=1.E-7                  
      !! Root Finding precision
      
      integer :: numstepvperp, numstepvpar           
      !! total number of steps in loop
      
      logical :: ex                                  
      !! used to check if results directory exists
      
      complex, dimension(1:3)       :: ef, bf
      !! E, B
      
      complex, dimension(1:nspec)     :: ns     
      !! density
      
      complex, dimension(1:3,1:nspec) :: Us     
      !! Velocity
      
      !Heating (Required parameters of calc eigen) ----------
      real, dimension(1:nspec) :: Ps 
      !! Power into/out of species
      
      real, dimension(1:4,1:nspec) :: Ps_split 
      !! Power into/out of species

      real, dimension(1:6,1:nspec) :: Ps_split_new 
      !! Power into/out of species (includes all diagonal terms added by GGH in 2023)

      real :: Ew 
      !! wave energy

      !loop counter/ loop parameters -------
      
      integer :: is                     
      !! species counter
      
      integer :: unit_s                 
      !! out file unit counter
      
      character(100)  :: fmt                          
      !! Eigenfunction Output Format

      real :: start, finish 
      !! debug/test to measure runtime of function

      real, allocatable, dimension(:) :: vvperp,vvpar,vvphi  
      !! Velocity grid values (norm: w_par_s)
      
      integer :: ivperp,ivpar,ivphi                    
      !! Index for each velocity dimension
      
      integer :: ivperpmin,ivperpmax,ivparmin,ivparmax,ivphimin,ivphimax  
      !! Index limits
      
      real :: vvperp1temp, vvperp2temp 
      !! Temporary velocity space coordinate variable
      
      real,  allocatable, dimension(:,:,:,:) :: fs0   
      !! Dimensionless equilibrium fs0
      
      complex, allocatable, dimension(:,:,:,:) :: fs1 
      !! Perturbed Dist for all species
      
      real :: fs0_temp   
      !! Temporary Dimensionless equilibrium fs0 at adjacent location of fs0 array in vperp1/vperp2 direction (used for derivatives)
      
      complex, allocatable, dimension(:,:,:,:) :: fs1_plus_delvperp1,fs1_plus_delvperp2,fs1_minus_delvperp1,fs1_minus_delvperp2   
      !! perturbed dist at adjacent location of fs0 array in vperp1/vperp2 direction (used for derivatives) !Perturbed Dist for all species
      
      real :: phi_adjacent 
      !! var used to compute fs0/fs1 at adjacent location in vperp1/vperp2 direction
      
      real :: vperp_adjacent,vperp1_adjacent,vperp2_adjacent 
      !! vars used to compute fs0/fs1 at adjacent location in vperp1/vperp2 direction
      
      real,  allocatable, dimension(:)  :: hatV_s     
      !! Flow normalized to wpar_s
      
      complex, allocatable, dimension(:,:,:,:) :: dfs1dvpar,dfs1dvperp1,dfs1dvperp2
      !! derivative of fs1 on the 3d grid for each species

      real, allocatable, dimension(:,:,:,:) :: corepar,coreperp 
      !! 3V Correlations (in cylindrical coords)

      real, allocatable, dimension(:,:,:) :: corepar_cyln,coreperp_cyln 
      !! 2V Corrs
      
      complex, allocatable, dimension(:,:,:) :: fs1_cyln     
      !! 2V fs1 
      
      complex, allocatable, dimension(:) :: ns1 
      !! Density Fluctuation
      
      complex, allocatable, dimension(:,:) :: us1 
      !! Fluid Velocity Fluctuation
      
      real, allocatable, dimension(:) :: jparepar,jperpeperp 
      !! int_v 3V Correlations (j_i times E_i) - currently not used TODO: implement or remove
      
      integer :: jj                                   
      !! Index
      
      real :: delv3                                  
      !! delv^3
      
      real :: pi
      !! 3.14159...

      complex    :: exbar               
      !! amplitude factor of fs1 which is \propto Ex/B0 (B0 is external B)

      character(100)  :: fmt_dbg1,fmt_dbg2           
      !! Eigenfunction Output Format used for debug

      pi = 4.0*ATAN(1.0)

      !check if results directory exists
      ! INQUIRE (DIRECTORY='data', EXIST=ex)
      ex = .true. !TODO: make this work for gfortran compiler
      if(ex) then
        write(*,*)"Assuming data folder already exists..."
      else
        write(*,*)"Creating data folder for output..."
        write(*,*)'mkdir data'
        CALL system('mkdir data')
        write(*,*)"Saving output to data folder..."
      endif

      write(outputPath,*) 'data/', trim(dataName) !!TODO: use more general pathing
      ! INQUIRE (DIRECTORY=trim(dataName), EXIST=ex)
      ex = .true.
      if(ex) then
        write(*,*)"assuming subfolder ", trim(dataName), "already exists"
      else
        write(*,*)"Creating data subfolder ", trim(dataName)
        write(cmd,*)'mkdir ',trim(outputPath)
        write(*,*)cmd
        CALL system(cmd)
        write(*,*)"Saving to data subfolder ", trim(dataName)
      endif

      !Grab dispersion relation solution
      wi = wroots(1,wrootindex)
      gi = wroots(2,wrootindex)
      ominit=cmplx(wi,gi)
      om1=ominit*(1.-prec)
      om2=ominit*(1.+prec)

      ! Refine Omega Value
      iflag=0
      omega=rtsec(disp,om1,om2,tol,iflag)
      
      call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,Ps_split_new,.true.,.true.)
  
      !==============================================================================
      ! Create Velocity grid and allocate variables
      !==============================================================================

      ! Create velocity grid variables: vvx,vvy,vvz
      ! NOTE: All (x,y,z) velocity values are normalized to species parallel thermal velocity
      ! NOTE: Code below assumes max.min values are integer multiples of delv (if not, throw warning flag!)
      ivparmin=int(vparmin/delv); if (real(ivparmin) .ne. vparmin/delv) &
      write(*,*)'WARNING: vparmin not integer multiple of delv'
      ivparmax=int(vparmax/delv); if (real(ivparmax) .ne. vparmax/delv) &
      write(*,*)'WARNING: vparmax not integer multiple of delv'
      ivperpmin=int(vperpmin/delv); if (real(ivperpmin) .ne. vperpmin/delv) &
      write(*,*)'WARNING: vperpmin not integer multiple of delv'
      ivperpmax=int(vperpmax/delv); if (real(ivperpmax) .ne. vperpmax/delv) &
      write(*,*)'WARNING: vperpmax not integer multiple of delv'
      ivphimin=0 !TODO: load this in input
      ivphimax=15 !WARNING: this is aggressively small

      !Allocate velocity grid variables
      allocate(vvperp(ivperpmin:ivperpmax)); vvperp(:)=0.
      allocate(vvpar(ivparmin:ivparmax)); vvpar(:)=0.
      allocate(vvphi(ivphimin:ivphimax)); vvphi(:)=0.
      !Populate velocity grid values
      do ivperp=ivperpmin,ivperpmax
         vvperp(ivperp)=real(ivperp)*delv
      enddo
      do ivpar=ivparmin,ivparmax
         vvpar(ivpar)=real(ivpar)*delv
      enddo
      do ivphi=ivphimin,ivphimax
         vvphi(ivphi)=real(ivphi)*2*pi/ivphimax
      enddo

      !Allocate fs1 and fs0 variables 
      allocate(fs0(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      fs0(:,:,:,:)=0.
      allocate(fs1(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      fs1(:,:,:,:)=0.
      allocate(fs1_plus_delvperp1(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      fs1_plus_delvperp1(:,:,:,:)=0.
      allocate(fs1_minus_delvperp1(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      fs1_minus_delvperp1(:,:,:,:)=0.
      allocate(fs1_plus_delvperp2(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      fs1_plus_delvperp2(:,:,:,:)=0.
      allocate(fs1_minus_delvperp2(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      fs1_minus_delvperp2(:,:,:,:)=0.

      !Allocate Bessel function values
      allocate(jbess(-nbesmax-1:nbesmax+1)); jbess(:)=0.

      !=============================================================================
      ! END Create Velocity grid and allocate variables
      !=============================================================================

      !=============================================================================
      ! Calculate fs0 and fs1 on (vperp,vpar,vphi) grid
      !=============================================================================
      allocate(hatV_s(nspec))

      !Loop over (vperp,vpar,vphi) grid and compute fs0 and fs1
      call calc_exbar(omega,ef,bf,exbar)
      do is = 1, nspec
        call check_nbesmax(MAX(ABS(vparmin),ABS(vparmax),ABS(vperpmin),ABS(vperpmax)),spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s)

        !Create variable for parallel flow velocity normalized to
        !       species parallel thermal speed
        hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(is)%mu_s*betap))
        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
            do ivphi=ivphimin,ivphimax
                  !Compute dimensionless equilibrium Distribution value, fs0
                  fs0(ivperp,ivpar,ivphi,is)=fs0hat(vvperp(ivperp),vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  !Compute perturbed  Distribution value, fs1
                  call calc_fs1(omega,vvperp(ivperp),vvpar(ivpar),vvphi(ivphi),ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,&
                                    elecdircontribution,(1.,0),(1.,0),exbar,fs0(ivperp,ivpar,ivphi,is),fs1(ivperp,ivpar,ivphi,is))

                  !compute fs1 at adjacent locations in vperp1/vperp2 direction to take derivatives with later
                  !Note: delv may not be the best choice here when it is large. Consider using a separate variable to determine locations that we approximate derivative at
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))+delv
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,elecdircontribution,(1.,0),(1.,0),exbar,fs0_temp,fs1_plus_delvperp1(ivperp,ivpar,ivphi,is))
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))-delv
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,elecdircontribution,(1.,0),(1.,0),exbar,fs0_temp,fs1_minus_delvperp1(ivperp,ivpar,ivphi,is))
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))+delv
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,elecdircontribution,(1.,0),(1.,0),exbar,fs0_temp,fs1_plus_delvperp2(ivperp,ivpar,ivphi,is))
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))-delv
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,elecdircontribution,(1.,0),(1.,0),exbar,fs0_temp,fs1_minus_delvperp2(ivperp,ivpar,ivphi,is))
              enddo
            enddo
          enddo
      enddo

      !=============================================================================
      ! End Calculate fs0 and fs1 on (vx,vy,vz) grid
      !=============================================================================

      !=============================================================================
      ! Calculate velocity derivatives of fs1
      !=============================================================================
      allocate(dfs1dvperp1(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      dfs1dvperp1=0.
      allocate(dfs1dvperp2(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      dfs1dvperp2=0.
      allocate(dfs1dvpar(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec))
      dfs1dvpar=0.

      do is = 1, nspec
        !dfs1/dvpar-------------------------------------------------------
        do ivperp=ivperpmin,ivperpmax
          do ivphi=ivphimin,ivphimax
            !First point: 1st order Finite difference
            dfs1dvpar(ivperp,ivpar,ivphi,is)=&
            (fs1(ivperp,ivpar+1,ivphi,is)-fs1(ivperp,ivpar,ivphi,is))/(vvpar(ivparmin+1)-vvpar(ivparmin));
            !Loop: 2nd order Centered Finite difference
            do ivpar=ivparmin+1,ivparmax-1
              dfs1dvpar(ivperp,ivpar,ivphi,is)=&
              (fs1(ivperp,ivpar+1,ivphi,is)-fs1(ivperp,ivpar-1,ivphi,is))/(vvpar(ivpar+1)-vvpar(ivpar-1));
            enddo
            !Last point: 1st order Finite difference !TODO: ivpar might need to be reset after this-> similar issue may exist in compute_fpc_cart_new. Fix both
            dfs1dvpar(ivperp,ivpar,ivphi,is)=& 
            (fs1(ivperp,ivpar,ivphi,is)-fs1(ivperp,ivpar-1,ivphi,is))/(vvpar(ivparmax)-vvpar(ivparmax-1));
          enddo
        enddo

        !dfs1/dvperp1-------------------------------------------------------
        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
            do ivphi=ivphimin,ivphimax
              dfs1dvperp1(ivperp,ivpar,ivphi,is)=&
              (fs1_plus_delvperp1(ivperp,ivpar,ivphi,is)-fs1_minus_delvperp1(ivperp,ivpar,ivphi,is))/(delv); !Note: this assumes a square grid with even spacings between all points
            enddo
          enddo
        enddo

        !dfs1/dvperp2-------------------------------------------------------
        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
            do ivphi=ivphimin,ivphimax
              dfs1dvperp2(ivperp,ivpar,ivphi,is)=&
              (fs1_plus_delvperp2(ivperp,ivpar,ivphi,is)-fs1_minus_delvperp2(ivperp,ivpar,ivphi,is))/(delv); !Note: this assumes a square grid with even spacings between all points
            enddo
          enddo
        enddo
      enddo

      !=============================================================================
      ! END Calculate velocity derivatives of fs1
      !=============================================================================

      !=============================================================================
      ! Calculate 3V CEpar, CEperp Correlations
      !=============================================================================
      allocate(coreperp(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec)); coreperp=0.
      allocate(corepar(ivperpmin:ivperpmax,ivparmin:ivparmax,ivphimin:ivphimax,1:nspec)); corepar=0.

      !Loop over (vx,vy,vz) grid and compute CEx, CEy, CEz Correlations
      do is = 1, nspec
         do ivperp=ivperpmin,ivperpmax
            do ivpar=ivparmin,ivparmax
               do ivphi=ivphimin,ivphimax
                  ! CEpar
                  corepar(ivperp,ivpar,ivphi,is)=&
                  -spec(is)%q_s*0.5*vvpar(ivpar)*vvpar(ivpar)* 0.5*(CONJG(dfs1dvpar(ivperp,ivpar,ivphi,is))*ef(3) &
                  +dfs1dvpar(ivperp,ivpar,ivphi,is)*CONJG(ef(3)))
                  
                  ! CEperp
                  vvperp1temp = vvperp(ivperp)*COS(vvphi(ivphi))
                  vvperp2temp = vvperp(ivperp)*SIN(vvphi(ivphi))
                  coreperp(ivperp,ivpar,ivphi,is)=&
                  (-spec(is)%q_s*0.5*vvperp1temp*vvperp1temp*(0.5*(CONJG(dfs1dvperp1(ivperp,ivpar,ivphi,is))*ef(1) &
                  +dfs1dvperp1(ivperp,ivpar,ivphi,is)*CONJG(ef(1))))) &
                  -(spec(is)%q_s*0.5*vvperp2temp*vvperp2temp*(0.5*(CONJG(dfs1dvperp2(ivperp,ivpar,ivphi,is))*ef(2) &
                  +dfs1dvperp2(ivperp,ivpar,ivphi,is)*CONJG(ef(2)))))
               enddo
            enddo
         enddo
      enddo   

      !=============================================================================
      ! END Calculate 3V CEpar, CEperp Correlations
      !=============================================================================

      !=============================================================================
      ! Reduce to 2V fs1 and 2V CEpar, CEperp Correlations
      !=============================================================================
      !Perturbed Distribution------------------------------------------------
      allocate(fs1_cyln(ivperpmin:ivperpmax,ivparmin:ivparmax,1:nspec)); fs1_cyln=0.

      !Reduce Distribution
      do is = 1, nspec
        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
            do ivphi=ivphimin,ivphimax
              fs1_cyln(ivperp,ivpar,is) =&
              fs1_cyln(ivperp,ivpar,is)+vvperp(ivperp)*fs1(ivperp,ivpar,ivphi,is)*2*pi/ivphimax
            enddo
          enddo
        enddo
      enddo


      !Reduce Correlations
      allocate(corepar_cyln(ivperpmin:ivperpmax,ivparmin:ivparmax,1:nspec)); corepar_cyln=0.
      allocate(coreperp_cyln(ivperpmin:ivperpmax,ivparmin:ivparmax,1:nspec)); coreperp_cyln=0.

      !Reduce Distribution
      do is = 1, nspec
        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
            do ivphi=ivphimin,ivphimax
              corepar_cyln(ivperp,ivpar,is) =& 
              corepar_cyln(ivperp,ivpar,is)+vvperp(ivperp)*corepar(ivperp,ivpar,ivphi,is)*2*pi/ivphimax
              coreperp_cyln(ivperp,ivpar,is) =&
              coreperp_cyln(ivperp,ivpar,is)+vvperp(ivperp)*coreperp(ivperp,ivpar,ivphi,is)*2*pi/ivphimax
            enddo
          enddo
        enddo
      enddo

      !=============================================================================
      ! END Reduce to 2V fs1 and 2V CEpar, CEperp Correlations
      !=============================================================================
      
      !=============================================================================
      ! Output Distributions and Correlations to Files
      !=============================================================================
      do is = 1, nspec
        !make file to store result
        !TODO: used "get unused unit" to get unit_s to pick correct 'number' to write to
        unit_s = 10+5*is !note: unit = 5,6 are reserved by standard fortran for input form keyboard/ writing to screen
        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cpar.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cperp.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+1,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),&
              '.df1gyro.real.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+2,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),&
              '.df1gyro.imag.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+3,file=trim(filename),status='replace')

        write(*,*)'Calculating fpc for species ',is
        write(*,*)'Writing omega/kpar V_a normalization to file...'

        write(unit_s,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s, '(5a22)')'vperpmin','vperpmax','vparmin','vparmax','delv'
        write(unit_s, '(5es22.7)')vperpmin,vperpmax,vparmin,vparmax,delv
        write(unit_s, *) '-------------'
        write(unit_s+1,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+1,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+1, '(5a22)')'vperpmin','vperpmax','vparmin','vparmax','delv'
        write(unit_s+1, '(5es22.7)')vperpmin,vperpmax,vparmin,vparmax,delv
        write(unit_s+1, *) '-------------'
        write(unit_s+2,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+2,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+2, '(5a22)')'vperpmin','vperpmax','vparmin','vparmax','delv'
        write(unit_s+2, '(5es22.7)')vperpmin,vperpmax,vparmin,vparmax,delv
        write(unit_s+2, *) '-------------'
        write(unit_s+3,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+3,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+3, '(5a22)')'vperpmin','vperpmax','vparmin','vparmax','delv'
        write(unit_s+3, '(5es22.7)')vperpmin,vperpmax,vparmin,vparmax,delv
        write(unit_s+3, *) '-------------'

        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
              if(ABS(corepar_cyln(ivperp,ivpar,is)) .lt. 9.999E-99) then
                write(unit_s,'(es17.5)',advance='no') 0.
              else
                write(unit_s,'(es17.5)',advance='no') corepar_cyln(ivperp,ivpar,is)
              endif
              if(ABS(coreperp_cyln(ivperp,ivpar,is)) .lt. 9.999E-99) then
                write(unit_s+1,'(es17.5)',advance='no') 0.
              else
                write(unit_s+1,'(es17.5)',advance='no') coreperp_cyln(ivperp,ivpar,is)
              endif
              if(ABS(real(fs1_cyln(ivperp,ivpar,is))) .lt. 9.999E-99) then
                write(unit_s+2,'(es17.5)',advance='no') 0.
              else
                write(unit_s+2,'(es17.5)',advance='no') real(fs1_cyln(ivperp,ivpar,is))
              endif
              if(ABS(aimag(fs1_cyln(ivperp,ivpar,is))) .lt. 9.999E-99) then
                write(unit_s+3,'(es17.5)',advance='no') 0.
              else
                write(unit_s+3,'(es17.5)',advance='no') aimag(fs1_cyln(ivperp,ivpar,is))
              endif
           enddo
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
        enddo

        close(unit_s)
        close(unit_s+1)
        close(unit_s+2)
        close(unit_s+3)

      end do

      !Deallocate variables
      ! deallocate(ns1,us1,jxex,jyey,jzez)
      deallocate(corepar_cyln,coreperp_cyln)
      deallocate(corepar,coreperp)
      deallocate(dfs1dvperp1,dfs1dvperp2,dfs1dvpar)
      deallocate(fs1_plus_delvperp1,fs1_minus_delvperp1,fs1_plus_delvperp2,fs1_minus_delvperp2)
      deallocate(fs0,fs1)
      deallocate(hatV_s)
      deallocate(vvperp,vvpar,vvphi)


    end subroutine compute_fpc_gyro

    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    real function fs0hat(vperp,vpar,hatV_s,aleph_s) 
      !! Determines dimensionless species equilibrium VDF \hat{fs0} at (vperp,vpar)
      
      !input
      real         :: vperp, vpar      
      !! normalized velocity space
      
      real         :: hatV_s           
      !! normalized species drift velocity by wpar_s
      
      real         :: q_s              
      !! normalized species charge
      
      real         :: aleph_s          
      !! T_perp/T_parallel_s
      
      real         :: tau_s            
      !! T_ref/T_s_parallel
      
      real         :: mu_s             
      !! m_ref/m_s

      fs0hat = exp(-1.*(vpar-hatV_s)**2.-vperp**2./aleph_s)

    end function fs0hat

    subroutine calc_exbar(omega,ef,bf,exbar)
      !! Computes the amplitdue factor of the peturb distribution
      !! which is propto Ex/B0, a value which is an input but can be related to 
      !! solved and input dimensionless quantities using the linearized lorentz force and amperes law
      !! and a lot of algebra

      use vars, only : betap,kperp,kpar,vtp,nspec,spec

      complex, intent(in)   :: omega            
      !! Complex Frequency
      
      complex, dimension(1:3), intent(in)   :: ef, bf           
      !! E, B
      
      complex, intent(out) :: exbar
      !! Amplitude factor of fs1

      real,  allocatable, dimension(:)  :: hatV_s     
      !! Flow normalized to wpar_s

      complex                           :: omega_temp 
      !! holds the fixed omege with sign of gamma term flipped as PLUME returns it with a minus sign 
      
      complex                           :: numerator  
      !! numerator term of amp term relate to exbar
      
      complex                           :: sumterm    
      !! term that holds running sum of summation
      
      complex                           :: runningterm
      !! term that holds numerator term in summation so we can break it up into many lines for readability
      
      integer                           :: is         
      !! species counter
      
      complex                           :: ii= (0,1.) 
      !! Imaginary unit: 0+1i
      
      real :: kpar_temp 
      !! Debug value that 'fixes' sign due to differences in sign of definitions between the textbooks by stix and swanson
      
      integer :: unit_number
      !! Holds unit number of file used for debug writing to file

      allocate(hatV_s(nspec))

      omega_temp = real(omega)-ii*aimag(omega)
      kpar_temp = kpar

      numerator = -(0,1.)*kpar_temp*bf(2)-(0,1.)*vtp*sqrt(spec(1)%alph_s)*omega_temp !Warning: assumes first species is reference species

      sumterm = (0.,0.)
      runningterm = (0.,0.)
      do is = 1, nspec
        !Create variable for parallel flow velmocity normmalized to
        !       species parallel thermal speed
        hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(1)%mu_s*betap))

        runningterm = -(0,1.)*omega_temp*spec(is)%q_s/spec(is)%mu_s
        runningterm = runningterm + (0.,1.)*omega_temp*spec(is)%q_s/spec(is)%mu_s*hatV_s(is)/spec(is)%tau_s*vtp
        runningterm = runningterm+ef(2)
        runningterm = runningterm+bf(1)*vtp*hatV_s(is)/spec(is)%tau_s
        runningterm = (spec(is)%D_s/spec(is)%q_s)*runningterm/(1-omega_temp**2*(spec(is)%q_s)**2/(spec(is)%mu_s)**2)

        sumterm = sumterm + runningterm
        write(*,*)'debug spec',is,'rt',runningterm,'1-om2',(1-omega_temp**2*(spec(is)%q_s)**2/(spec(is)%mu_s)**2)
        write(*,*)'ef(2)',ef(2),'Ds/Qs',(spec(is)%D_s/spec(is)%q_s),'4thom',omega_temp*spec(is)%q_s/spec(is)%mu_s
        write(*,*)'spec debug',is,spec(is)%mu_s,spec(is)%q_s
        write(*,*)'omega_temp and sqrd',omega_temp,omega_temp**2,'omega and sqrd',omega,omega**2
        write(*,*)''
        runningterm = (0.,0.)
      enddo

      exbar = numerator/((sqrt(betap)*spec(1)%D_s/(vtp**2*sqrt(spec(1)%alph_s)))*sumterm) !compute exbar/B0 (wperp/vAR) !Warning: assumes first species is reference species
      exbar = exbar/(vtp*sqrt(spec(1)%alph_s))

      write(*,*)'--------'
      write(*,*)''

      ! Write to file
      open(newunit=unit_number, file="exbar_output.dat", status="replace", action="write")
      write(unit_number,*) exbar
      close(unit_number)
      
      exbar = (1.,0.)

    end subroutine calc_exbar


    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    
    subroutine calc_fs1(omega,vperp,vpar,phi,ef,bf,hatV_s,q_s,aleph_s,tau_s,mu_s,aleph_r,elecdircontribution,A1,B1,exbar,fs0,fs1)
      !! Determine species perturbed VDF fs1 at given (vperp,vpar,phi)

      use vars, only : betap,kperp,kpar,vtp,pi
      use vars, only : nbesmax
      use bessels, only : bessj_s, bess0_s_prime
      USE nrtype, only: SP, I4B

      complex, intent(in)   :: omega            
      !! Complex Frequency
      
      real, intent(in)      :: vperp, vpar      
      !! normalized velocity space
      
      real, intent(in)      :: phi              
      !! azimuthal angle in velocity space
      
      complex, dimension(1:3), intent(in)   :: ef, bf           
      !! E, B eigenfunctions in all 3 directions

      real, intent(in)      :: hatV_s           
      !! Drift velocity norm to wpar_s
      
      real, intent(in)      :: q_s              
      !! normalized species charge
      
      real, intent(in)      :: aleph_s          
      !! T_perp/T_parallel_s
      
      real, intent(in)      :: tau_s            
      !! T_ref/T_s|_parallel
      
      real, intent(in)      :: mu_s             
      !! m_ref/m_s
      
      real, intent(in)      :: aleph_r          
      !! T_perp/T_parallel_R
      
      real, intent(in)      :: elecdircontribution 
      !! Sets components of Electric field (0 (DEFAULT) (or any other value) = Do not modify, 1=Keep only Ex(i.e.Eperp1), 2=Keep only Ey(i.e.Eperp2), 3=Keep only Ez(i.e.Epar))
      
      complex, intent(in)   :: A1, B1           
      !! Temporary scalars to fix coeff error!
      
      complex, intent(out)  :: exbar            
      !! normalizaiton factor
      
      real, intent(in)      :: fs0              
      !! normalized zero order distribution
      
      complex, intent(out)  :: fs1              
      !! first order distribution

      integer :: n 
      !! bessel sum counter
      
      integer :: m 
      !! bessel sum counter
      
      complex :: ii= (0,1.) 
      !! Imaginary unit: 0+1i
      
      real :: b_s                           
      !! Bessel function argument

      !intermediate values (see calculation notes for definitions)
      real :: kpar_temp
      real :: kperp_temp
      complex :: denom
      complex :: emult
      complex :: Wbar_s
      complex :: Ubar_s

      complex :: ef1,ef2,ef3
      complex :: omega_temp
      real :: phi_temp
      real :: vpar_temp

      phi_temp = phi

      !Used to 'turn off' contributions by other electric fields...
      if (elecdircontribution == 1.) then
         ef1 = ef(1)
         ef2 = 0.
         ef3 = 0
      else if (elecdircontribution == 2.) then
         ef1 = 0.
         ef2 = ef(2)
         ef3 = 0.
      else if (elecdircontribution == 3.) then
         ef1 = 0.
         ef2 = 0.
         ef3 = ef(3)
      else
         ef1 = ef(1)
         ef2 = ef(2)
         ef3 = ef(3)
      end if

      !fix sign definition difference between swanson/ stix
      !Note, this sign difference causes for a strange mixture of signs in the terms (namely in Ubar_s and Wbar_s) but this has been tested!
      if (q_s .gt. 0.) then 
          omega_temp = -real(omega)-ii*aimag(omega) 
          kpar_temp = -kpar !want to keep the same direction
          kperp_temp = kperp 
          vpar_temp = vpar
          ef3 = ef3
          ef2 = ef2
          ef1 = ef1
      else
        omega_temp = real(omega)-ii*aimag(omega)
        kpar_temp = kpar
        vpar_temp = vpar
        ef3 = ef3
        ef2 = ef2
        ef1 = ef1
      end if

      !Compute Bessel functions (if necessary)
      b_s = (kperp_temp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)
      if (b_s .ne. bs_last) then
         !Populate jbess with current values
         do n= -nbesmax-1,nbesmax+1
             jbess(n) = bessj_s(n,b_s)
         enddo
         bs_last=b_s 
      end if

      !Double Bessel Sum to calculate fs1=========================================
      fs1 = (0.,0.)
      !Calculate all parts of solution that don't depend on m or n
      Ubar_s= -2.*vperp/aleph_s*(1.+kpar*sqrt(mu_s/(tau_s*aleph_r))/(real(omega)-ii*aimag(omega))*((aleph_s-1)*vpar_temp-aleph_s*hatV_s))

      do n = -nbesmax,nbesmax
       !Calculate all parts of solution that dosn't depend on m
       denom=(omega_temp-kpar_temp*vpar_temp*sqrt(mu_s/(tau_s*aleph_r))-n*mu_s/q_s)
       Wbar_s=2.*(n*mu_s/(q_s*(real(omega)-ii*aimag(omega)))-1.)*(vpar_temp-hatV_s) - 2.*(n*mu_s/(q_s*(real(omega)-ii*aimag(omega))*aleph_s))*vpar_temp
       if (b_s .ne. 0.) then  !Handle division of first term if b_s=0 (U_bar_s also =0)
          emult=n*jbess(n)*Ubar_s/(b_s)*ef1
          if(n .ne. 0) then
            emult=emult+ii*0.5*(jbess(n-1)-jbess(n+1))*Ubar_s*ef2
          else
            emult=emult+ii*bess0_s_prime(b_s)*Ubar_s*ef2
          end if
          emult=emult+jbess(n)*Wbar_s*ef3
       else
          if(n .ne. 0) then
            emult=+ii*0.5*(jbess(n-1)-jbess(n+1))*Ubar_s*ef2 
          else
            emult=+ii*bess0_s_prime(b_s)*Ubar_s*ef2 
          end if 
          emult=emult + jbess(n)*Wbar_s*ef3
       end if
       
       do m = -nbesmax,nbesmax
          fs1=fs1+jbess(m)*exp(ii*(m-n)*phi_temp)*emult/denom
       enddo
       enddo
      fs1 = -1.*ii*sqrt(mu_s*tau_s/betap)*(exbar/q_s)*fs1*fs0

    end subroutine calc_fs1


   subroutine check_nbesmax(vperpmax,tau_s,mu_s,aleph_r)
      !! We approximate the infitine sums, which include terms like j_n(b), as finite sums from n=-Nlarge to Nlarge.
      !! Since j_n(b) is pretty small when b = n/2, we just make sure that n is large enough for all b (which is related to vperp) values

      !As j_n(b) is small for b<n/2, nbesmax/2 should be greater than or equal to b_s,max = |(kperp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)| for all species
      !TODO: if using wrapper, the user will not see this error- find a way to warn user if using wrapper (or atleast leave note to check output!)
      !TODO: alternatively, automatically set nbesmax to be some multiple of the minimum value calculated using this formula...
      use vars, only : kperp
      use vars, only : nbesmax

      real, intent(in)                      :: vperpmax         
      !! max value to be computed in normalized velocity space
      
      real, intent(in)                      :: tau_s            
      !! T_ref/T_s|_parallel
      
      real, intent(in)                      :: mu_s             
      !! m_ref/m_s
      
      real, intent(in)                      :: aleph_r          
      !! T_perp/T_parallel_R

      real                                  :: b_s              
      !! argument to bessel functions (see calc_fs1)

      b_s = ABS(kperp*vperpmax/sqrt(mu_s*tau_s*aleph_r))

      if(nbesmax/2 < b_s) then
        write(*,*)"***************************************"
        write(*,*)"Warning in check_nbesmax in fpc.f90!!! (see for more info)"
        write(*,*)"nbesmax is too small!!!"
        write(*,*)"nbesmax should be great than b_s = ABS(kperp*vperpmax/sqrt(mu_s*tau_s*aleph_r)) for all species!!!"
        write(*,*)"***************************************"
      endif

   end subroutine check_nbesmax
    
end module fpc
