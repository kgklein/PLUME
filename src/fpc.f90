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

  !CB Note: due to an effective sign difference between the textbooks by Stix/Swanson (used to compute eigen functions and fs1 respectively),
  !      we drop the minus sign in the correlation C_{E_i} = -q/2 v_i^2 d fs1/d v_i, so it is "q/2 v_i^2 d fs1/d v_i" when implemented...
  !      This minus sign difference is a result in Swanson 'suppressing the exp[i(k dot r - om t)]' (above 4.180) rather than taking the fourier transform of the product
  !      Note that the fourier transform of this term is F{exp[i(k dot r - om t)]} \proportional_to delta(kx`+kx) delta(ky`+ky) delta(kz`+kz) delta(om`-om)
  !      Effectively swaping the sign in kx ky kz (equivalently changing the frame). Typically, this does not matter as one can still get the same eigen functions, 
  !      but since we are correlating quantities in two diffferent frames, we must account for it somewhere...

  implicit none
  private :: calc_correlation_par_gyro, calc_correlation_perp_gyro, calc_correlation_i_car, calc_fs0, calc_fs1, &
  calc_fs1_new, calc_fs1_new2

  public :: compute_fpc_gyro, compute_fpc_cart, write_fs0

  !GGH: Optimized Routines
  public :: compute_fpc_cart_new

  !CRB: Optimized Routine
  public :: compute_fpc_gyro_new !TODO: remove old slow routines

  !NEW GLOBAL VARIABLES=========================================================
  real :: bs_last=0.0       !Last Bessel function argument
  real, allocatable :: jbess(:)  !Regular Bessel function values
  real :: eperp1_bar=1.0    !Overall amplitude of the linear mode TODO: move to calc_fs1 and make it local
  !END NEW GLOBAL VARIABLES=====================================================

  contains
    
    !------------------------------------------------------------------------------
    !                           Collin Brown, 2020
    !------------------------------------------------------------------------------
    subroutine compute_fpc_gyro(wrootindex)
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              !index of selected root

      character(1000) :: filename                      !Output File name
      character(1000) :: outputPath                    !Output folder
      character(1000) :: cmd                           !Varaible to store command line commands
      real    :: vperpi, vpari                       !normalized velocity space current value in loop
      integer :: vperpindex, vparindex               !loop counters
      complex :: omega                               !Complex Frequency
      real    :: Cor_par_s, Cor_perp_s                !normalized correlation value
      complex :: fs1, dfs1perp, dfs1z                !normalized distribution function and first partial derivatives
      real    :: wi,gi                               !Freq and Damping of initial guess
      complex :: ominit                              !Complex Frequency initial guess
      complex :: om1,om2                             !Bracket Values
      integer :: iflag                               !Flag for Root search
      real, parameter :: tol=1.0E-13                 !Root Search Tolerance
      real, parameter :: prec=1.E-7                  !Root Finding precision
      integer :: numstepvperp, numstepvpar           !total number of steps in loop
      logical :: ex                                  !used to check if results directory exists
      complex, dimension(1:3)       :: ef, bf !E, B
      complex, dimension(1:nspec)     :: ns     !density
      complex, dimension(1:3,1:nspec) :: Us     !Velocity
      !Heating (Required parameters of calc eigen)
      real, dimension(1:nspec) :: Ps !Power into/out of species
      real, dimension(1:4,1:nspec) :: Ps_split !Power into/out of species
      real, dimension(1:6,1:nspec) :: Ps_split_new !Power into/out of species (GGH)
      real :: Ew !wave energy
      !loop counter/ loop parameters
      integer :: is                     !species counter
      integer :: unit_s                 !out file unit counter
      character(100)  :: fmt                          !Eigenfunction Output Format

      real :: start, finish !debug/test to measure runtime of function

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
   
      do is = 1, nspec
        call check_nbesmax(MAX(ABS(vparmin),ABS(vparmax),ABS(vperpmin),ABS(vperpmax)),spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s)

        !make file to store result
        !TODO: used "get unused unit" to get unit_s to pick correct 'number' to write to
        unit_s = 10+5*is !note: unit = 5,6 are reserved by standard fortran for input form keyboard/ writing to screen
        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cpar.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cperp.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+1,file=trim(filename),status='replace')

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

        !setup loop variables
        numstepvperp = int((vperpmax-vperpmin)/delv)
        numstepvpar = int((vparmax-vparmin)/delv)
        vperpi = vperpmin
        vpari = vparmin

        do vperpindex = 0, numstepvperp
          write(*,*)'vperp = ',vperpi !quick debug 'progress bar'
          do vparindex = 0, numstepvpar
            call calc_correlation_par_gyro(omega,ef,bf,vperpi,vpari,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s,fs1,dfs1z)
            call calc_correlation_perp_gyro(omega,ef,bf,vperpi,vpari,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp_s,fs1,dfs1perp)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp_s) .lt. 9.999E-99) Cor_perp_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp_s
            vpari = vpari+delv
            !call cpu_time(finish) !debug/ performance benchmarking
          end do
          vpari = vparmin
          vperpi = vperpi+delv
          write(unit_s,*)
          write(unit_s+1,*)
        end do
        vpari = vparmin
        vperpi = vperpmin
      end do

    end subroutine compute_fpc_gyro

    subroutine compute_fpc_cart(wrootindex)
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              !index of selected root

      character(1000) :: filename                      !Output File name
      character(1000) :: outputPath                    !Output folder
      character(1000) :: cmd                           !Varaible to store command line commands
      integer :: ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax  !Index limits 
      real    :: vxi, vyi, vzi                       !normalized velocity space current value in loop (note: vx, vy, vz corresponds to vperp1,vperp2,vpar, but we use 'x','y','z' as convention)
      integer :: vxindex, vyindex, vzindex           !loop counters
      real    :: vmax3rdval                          !sampled range when computing projection
      complex :: omega                               !Complex Frequency
      real    :: Cor_par_s, Cor_perp1_s, Cor_perp2_s !normalized correlation value
      complex :: fs1                                 !normalized distribution function
      real    :: wi,gi                               !Freq and Damping of initial guess
      complex :: ominit                              !Complex Frequency initial guess
      complex :: om1,om2                             !Bracket Values
      integer :: iflag                               !Flag for Root search
      real, parameter :: tol=1.0E-13                 !Root Search Tolerance
      real, parameter :: prec=1.E-7                  !Root Finding precision
      integer :: numstepvx, numstepvy, numstepvz     !total number of steps in loop
      logical :: ex                                  !used to check if results directory exists
      complex, dimension(1:3)       :: ef, bf !E, B
      complex, dimension(1:nspec)     :: ns     !density
      complex, dimension(1:3,1:nspec) :: Us     !Velocity
      !Heating (Required parameters of calc eigen)
      real, dimension(1:nspec) :: Ps !Power into/out of species
      real, dimension(1:4,1:nspec) :: Ps_split !Power into/out of species
      real, dimension(1:6,1:nspec) :: Ps_split_new !Power into/out of species (GGH)
      real :: Ew !wave energy
      !loop counter/ loop parameters
      integer :: is                     !species counter
      integer :: unit_s                 !out file unit counter
      character(100)  :: fmt                          !Eigenfunction Output Format

      real :: start, finish !debug/test to measure runtime of function

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

      do is = 1, nspec
        call check_nbesmax(MAX(ABS(vxmin),ABS(vxmax),ABS(vymin),ABS(vymax),ABS(vzmin),ABS(vzmax))&
          ,spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s)

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

        write(unit_s,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s, '(7a22)')'vxmin','vxmax','vymin','vymax','vzmin','vzmax','delv'
        write(unit_s, '(7es22.7)')vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
        write(unit_s, *) '-------------'
        write(unit_s+1,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+1,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+1, '(7a22)')'vxmin','vxmax','vymin','vymax','vzmin','vzmax','delv'
        write(unit_s+1, '(7es22.7)')vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
        write(unit_s+1, *) '-------------'
        write(unit_s+2,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+2,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
                          real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+2, '(7a22)')'vxmin','vxmax','vymin','vymax','vzmin','vzmax','delv'
        write(unit_s+2, '(7es22.7)')vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
        write(unit_s+2, *) '-------------'
        write(unit_s+3,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+3,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
             real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+3, '(7a22)')'vxmin','vxmax','vymin','vymax','vzmin','vzmax','delv'
        write(unit_s+3, '(7es22.7)')vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
        write(unit_s+3, *) '-------------'
        write(unit_s+4,'(8a22)')'tau','bi','kpar','kperp','vti','mu','omega.r','omega.i'
        write(unit_s+4,'(8es22.7)')spec(is)%tau_s,betap,kpar,kperp,vtp,spec(is)%mu_s,&
             real(omega*sqrt(betap)/kpar),aimag(omega*sqrt(betap)/kpar)
        write(unit_s+4, '(7a22)')'vxmin','vxmax','vymin','vymax','vzmin','vzmax','delv'
        write(unit_s+4, '(7es22.7)')vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv
        write(unit_s+4, *) '-------------'

        !setup loop variables
        ivxmin=int(vxmin/delv); if (real(ivxmin) .ne. vxmin/delv) write(*,*)'WARNING: vxmin not integer multiple of delv'
        ivxmax=int(vxmax/delv); if (real(ivxmax) .ne. vxmax/delv) write(*,*)'WARNING: vxmax not integer multiple of delv'
        ivymin=int(vymin/delv); if (real(ivymin) .ne. vymin/delv) write(*,*)'WARNING: vymin not integer multiple of delv'
        ivymax=int(vymax/delv); if (real(ivymax) .ne. vymax/delv) write(*,*)'WARNING: vymax not integer multiple of delv'
        ivzmin=int(vzmin/delv); if (real(ivzmin) .ne. vzmin/delv) write(*,*)'WARNING: vzmin not integer multiple of delv'
        ivzmax=int(vzmax/delv); if (real(ivzmax) .ne. vzmax/delv) write(*,*)'WARNING: vzmax not integer multiple of delv'
      
        numstepvx = (ivxmax-ivxmin)
        numstepvy = (ivymax-ivymin)
        numstepvz = (ivzmax-ivzmin)

        !CEi(vx,vy)----------------------------------------------------------------------------
        vxi = vxmin
        vyi = vymin
        vzi = vzmin

        do vxindex = 0, numstepvx
          do vyindex = 0, numstepvy
            vmax3rdval = vzmax
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,3,3,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s,fs1)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,3,1,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp1_s,fs1)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,3,2,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp2_s,fs1)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp1_s) .lt. 9.999E-99) Cor_perp1_s = 0. !file formating bug fix
            if(ABS(Cor_perp2_s) .lt. 9.999E-99) Cor_perp2_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp1_s
            write(unit_s+2,'(es17.5)',advance='no')Cor_perp2_s
            write(unit_s+3,'(2es17.5)',advance='no')real(fs1)
            write(unit_s+4,'(2es17.5)',advance='no')aimag(fs1)
            vyi = vyi+delv
          end do
          vyi = vymin
          vxi = vxi+delv
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
          write(unit_s+4,*)
        end do
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'
        write(unit_s+3,*)'---'
        write(unit_s+4,*)'---'

        !CEi(vx,vz)----------------------------------------------------------------------------
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        do vxindex = 0, numstepvx
          do vzindex = 0, numstepvz
            vmax3rdval = vymax
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,2,3,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s,fs1)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,2,1,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp1_s,fs1)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,2,2,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp2_s,fs1)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp1_s) .lt. 9.999E-99) Cor_perp1_s = 0. !file formating bug fix
            if(ABS(Cor_perp2_s) .lt. 9.999E-99) Cor_perp2_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp1_s
            write(unit_s+2,'(es17.5)',advance='no')Cor_perp2_s
            write(unit_s+3,'(2es17.5)',advance='no')real(fs1)
            write(unit_s+4,'(2es17.5)',advance='no')aimag(fs1)   
            vzi = vzi+delv
          end do
          vzi = vzmin
          vxi = vxi+delv
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
          write(unit_s+4,*)
        end do    
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'
        write(unit_s+3,*)'---'
        write(unit_s+4,*)'---'

        !CEi(vy,vz)----------------------------------------------------------------------------
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        do vyindex = 0, numstepvy
          do vzindex = 0, numstepvz
            vmax3rdval = vxmax
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,1,3,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s,fs1)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,1,1,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp1_s,fs1)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,1,2,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp2_s,fs1)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp1_s) .lt. 9.999E-99) Cor_perp1_s = 0. !file formating bug fix
            if(ABS(Cor_perp2_s) .lt. 9.999E-99) Cor_perp2_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp1_s
            write(unit_s+2,'(es17.5)',advance='no')Cor_perp2_s
            write(unit_s+3,'(2es17.5)',advance='no')real(fs1)
            write(unit_s+4,'(2es17.5)',advance='no')aimag(fs1)  
            vzi = vzi+delv
          end do
          vzi = vzmin
          vyi = vyi+delv
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
          write(unit_s+3,*)
          write(unit_s+4,*)
        end do
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
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
      
    end subroutine compute_fpc_cart


    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    !TODO: add write/load fs1 to gyro
    subroutine compute_fpc_cart_new(wrootindex)
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vxmin,vxmax,vymin,vymax,vzmin,vzmax,delv,nbesmax
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              !index of selected root

      character(1000) :: filename                      !Output File name
      character(1000) :: outputPath                    !Output folder
      character(1000) :: cmd                           !Varaible to store command line commands
      real    :: vxi, vyi, vzi                       !normalized velocity space current value in loop (note: vx, vy, vz corresponds to vperp1,vperp2,vpar, but we use 'x','y','z' as convention)
      integer :: vxindex, vyindex, vzindex           !loop counters
      real    :: vmax3rdval                          !sampled range when computing projection
      complex :: omega                               !Complex Frequency
      real    :: Cor_par_s, Cor_perp1_s, Cor_perp2_s !normalized correlation value TODO: remove these and other unused variables...
!      complex :: fs1                                 !normalized distribution function
      real    :: wi,gi                               !Freq and Damping of initial guess
      complex :: ominit                              !Complex Frequency initial guess
      complex :: om1,om2                             !Bracket Values
      integer :: iflag                               !Flag for Root search
      real, parameter :: tol=1.0E-13                 !Root Search Tolerance
      real, parameter :: prec=1.E-7                  !Root Finding precision
      integer :: numstepvx, numstepvy, numstepvz     !total number of steps in loop
      logical :: ex                                  !used to check if results directory exists
      complex, dimension(1:3)       :: ef, bf !E, B
      complex, dimension(1:nspec)     :: ns     !density
      complex, dimension(1:3,1:nspec) :: Us     !Velocity
      !Heating (Required parameters of calc eigen)
      real, dimension(1:nspec) :: Ps !Power into/out of species
      real, dimension(1:4,1:nspec) :: Ps_split !Power into/out of species
      real, dimension(1:6,1:nspec) :: Ps_split_new !Power into/out of species (GGH)
      real :: Ew !wave energy
      !loop counter/ loop parameters
      integer :: is                     !species counter
      integer :: unit_s                 !out file unit counter

      real :: start, finish !debug/test to measure runtime of function
      !NEW VARIABLES================================================================
      real, allocatable, dimension(:) :: vvx,vvy,vvz  !Velocity grid values (norm: w_par_s)
      integer :: ivx,ivy,ivz                    !Index for each velocity dimension
      integer :: ivxmin,ivxmax,ivymin,ivymax,ivzmin,ivzmax  !Index limits 
      real,  allocatable, dimension(:,:,:,:) :: fs0   !Dimensionless equilibrium fs0
      complex, allocatable, dimension(:,:,:,:) :: fs1 !Perturbed Dist for all species
      real,  allocatable, dimension(:)  :: hatV_s     !Flow normalized to wpar_s
      real :: vperp                                   !Perpendicular velocity
      real :: phi                                     !Gyrophase angle (azimuthal)
      complex, allocatable, dimension(:,:,:,:) :: dfs1dvx,dfs1dvy,dfs1dvz !Derivatives
      real, allocatable, dimension(:,:,:,:) :: corex,corey,corez !3V Correlations
      real, allocatable, dimension(:,:,:) :: corex_xy,corex_xz,corex_zy !2V Corrs
      real, allocatable, dimension(:,:,:) :: corey_xy,corey_xz,corey_zy !2V Corrs
      real, allocatable, dimension(:,:,:) :: corez_xy,corez_xz,corez_zy !2V Corrs
      complex, allocatable, dimension(:,:,:) :: fs1_xy,fs1_xz,fs1_zy !2V fs1 
      complex, allocatable, dimension(:) :: ns1 !Density Fluctuation
      complex, allocatable, dimension(:,:) :: us1 !Fluid Velocity Fluctuation
      real, allocatable, dimension(:) :: jxex,jyey,jzez !int_v 3V Correlations
      integer :: jj                                   !Index
      character(100)  :: fmt                          !Eigenfunction Output Format
      real :: delv3                                  !delv^3
      real :: pi                                     !3.1415....
      character(100)  :: fmt_dbg1,fmt_dbg2           !Eigenfunction Output Format

      !END NEW VARIABLES============================================================

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
      do is = 1, nspec
         !Create variable for parallel flow velocity normalized to
         !       species parallel thermal speed
         hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(is)%mu_s*betap))
         do ivx=ivxmin,ivxmax
            do ivy=ivymin,ivymax
               vperp=sqrt(vvx(ivx)*vvx(ivx)+vvy(ivy)*vvy(ivy))
               do ivz=ivzmin,ivzmax
                  !Compute dimensionless equilibrium Distribution value, fs0
                  fs0(ivx,ivy,ivz,is)=fs0hat_new(vperp,vvz(ivz),hatV_s(is),spec(is)%alph_s)
                  !Compute perturbed  Distribution value, fs1
                  phi = ATAN2(vvy(ivy),vvx(ivx))
                  call calc_fs1_new(omega,vperp,vvz(ivz),phi,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,fs0(ivx,ivy,ivz,is),fs1(ivx,ivy,ivz,is))
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
      !Integrate 0th and 1st Fluid moments of fs1
      delv3=delv*delv*delv
      allocate(ns1(nspec)); ns1=0.
      allocate(us1(3,nspec)); us1=0.
      do is = 1, nspec
      ! Density Fluctuation: Zeroth Moment of delta f
         ns1(is)=sum(sum(sum(fs1(:,:,:,is),3),2),1)*delv3
         !Correct Normalization to n_0R
         ns1(is)=ns1(is)*sqrt(spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))*spec(is)%D_s
      ! Fluid Velocity: First Moment of total f = delta f (since int v f_0=0)
         !x-component
         do ivx=ivxmin,ivxmax
            us1(1,is)=us1(1,is)+vvx(ivx)*sum(sum(fs1(ivx,:,:,is),2),1)*delv3
         enddo
         !y-component
         do ivy=ivymin,ivymax
            us1(2,is)=us1(2,is)+vvy(ivy)*sum(sum(fs1(:,ivy,:,is),2),1)*delv3
         enddo
         !z-component
         do ivz=ivzmin,ivzmax
            us1(3,is)=us1(3,is)+vvz(ivz)*sum(sum(fs1(:,:,ivz,is),2),1)*delv3
         enddo

         !Correct Normalization to v_AR
         us1(:,is)=us1(:,is)*sqrt(betap*spec(is)%mu_s/(pi*pi*pi*spec(is)%tau_s))/spec(is)%alph_s
      enddo

      !Integrate Correlations
      allocate(jxex(nspec)); jxex=0.
      allocate(jyey(nspec)); jyey=0.
      allocate(jzez(nspec)); jzez=0.
      do is = 1, nspec
      !int_v CEx= jxEx
         jxex(is)=sum(sum(sum(corex(:,:,:,is),3),2),1)*delv3
      !int_v CEy= jyEy
         jyey(is)=sum(sum(sum(corey(:,:,:,is),3),2),1)*delv3
      !int_v CEz= jzEz
         jzez(is)=sum(sum(sum(corez(:,:,:,is),3),2),1)*delv3
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
        do jj=0,5
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
      write(filename,'(5A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.mom.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
      open(unit=unit_s+5,file=trim(filename),status='replace')

      !Write format 
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
      
      !==========================================================================
      !==========================================================================
      !==========================================================================
    end subroutine compute_fpc_cart_new


    !------------------------------------------------------------------------------
    !                           Collin Brown, 2020
    !------------------------------------------------------------------------------
    subroutine compute_fpc_gyro_new(wrootindex)
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv,nbesmax
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              !index of selected root

      character(1000) :: filename                      !Output File name
      character(1000) :: outputPath                    !Output folder
      character(1000) :: cmd                           !Varaible to store command line commands
      real    :: vperpi, vpari                       !normalized velocity space current value in loop
      integer :: vperpindex, vparindex               !loop counters
      complex :: omega                               !Complex Frequency
      real    :: Cor_par_s, Cor_perp_s                !normalized correlation value
      real    :: wi,gi                               !Freq and Damping of initial guess
      complex :: ominit                              !Complex Frequency initial guess
      complex :: om1,om2                             !Bracket Values
      integer :: iflag                               !Flag for Root search
      real, parameter :: tol=1.0E-13                 !Root Search Tolerance
      real, parameter :: prec=1.E-7                  !Root Finding precision
      integer :: numstepvperp, numstepvpar           !total number of steps in loop
      logical :: ex                                  !used to check if results directory exists
      complex, dimension(1:3)       :: ef, bf !E, B
      complex, dimension(1:nspec)     :: ns     !density
      complex, dimension(1:3,1:nspec) :: Us     !Velocity
      !Heating (Required parameters of calc eigen)
      real, dimension(1:nspec) :: Ps !Power into/out of species
      real, dimension(1:4,1:nspec) :: Ps_split !Power into/out of species
      real, dimension(1:6,1:nspec) :: Ps_split_new !Power into/out of species (GGH)
      real :: Ew !wave energy
      !loop counter/ loop parameters
      integer :: is                     !species counter
      integer :: unit_s                 !out file unit counter
      character(100)  :: fmt                          !Eigenfunction Output Format

      real :: start, finish !debug/test to measure runtime of function

      real, allocatable, dimension(:) :: vvperp,vvpar,vvphi  !Velocity grid values (norm: w_par_s)
      integer :: ivperp,ivpar,ivphi                    !Index for each velocity dimension
      integer :: ivperpmin,ivperpmax,ivparmin,ivparmax,ivphimin,ivphimax  !Index limits
      real :: vvperp1temp, vvperp2temp !Temporary velocity space coordinate variable
      real,  allocatable, dimension(:,:,:,:) :: fs0   !Dimensionless equilibrium fs0
      complex, allocatable, dimension(:,:,:,:) :: fs1 !Perturbed Dist for all species
      real :: fs0_temp   !Temporary Dimensionless equilibrium fs0 at adjacent location of fs0 array in vperp1/vperp2 direction (used for derivatives)
      complex, allocatable, dimension(:,:,:,:) :: fs1_plus_delvperp1,fs1_plus_delvperp2,fs1_minus_delvperp1,fs1_minus_delvperp2   !perturbed dist at adjacent location of fs0 array in vperp1/vperp2 direction (used for derivatives) !Perturbed Dist for all species
      real :: phi_adjacent !var used to compute fs0/fs1 at adjacent location in vperp1/vperp2 direction
      real :: vperp_adjacent,vperp1_adjacent,vperp2_adjacent !vars used to compute fs0/fs1 at adjacent location in vperp1/vperp2 direction
      real,  allocatable, dimension(:)  :: hatV_s     !Flow normalized to wpar_s
      complex, allocatable, dimension(:,:,:,:) :: dfs1dvpar,dfs1dvperp1,dfs1dvperp2
      real, allocatable, dimension(:,:,:,:) :: corepar,coreperp !3V Correlations (in cylindrical coords)
      real, allocatable, dimension(:,:,:) :: corepar_cyln,coreperp_cyln !2V Corrs
      complex, allocatable, dimension(:,:,:) :: fs1_cyln     !2V fs1 
      complex, allocatable, dimension(:) :: ns1 !Density Fluctuation
      complex, allocatable, dimension(:,:) :: us1 !Fluid Velocity Fluctuation
      real, allocatable, dimension(:) :: jparepar,jperpeperp !int_v 3V Correlations (j_i times E_i) TODO: implement...
      integer :: jj                                   !Index
      real :: delv3                                  !delv^3
      real :: pi
      character(100)  :: fmt_dbg1,fmt_dbg2           !Eigenfunction Output Format

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
      do is = 1, nspec
        call check_nbesmax(MAX(ABS(vparmin),ABS(vparmax),ABS(vperpmin),ABS(vperpmax)),spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s)

        !Create variable for parallel flow velocity normalized to
        !       species parallel thermal speed
        hatV_s(is)=spec(is)%vv_s*sqrt(spec(is)%tau_s/(spec(is)%mu_s*betap))
        do ivperp=ivperpmin,ivperpmax
          do ivpar=ivparmin,ivparmax
            do ivphi=ivphimin,ivphimax
                  !Compute dimensionless equilibrium Distribution value, fs0
                  fs0(ivperp,ivpar,ivphi,is)=fs0hat_new(vvperp(ivperp),vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  !Compute perturbed  Distribution value, fs1
                  call calc_fs1_new(omega,vvperp(ivperp),vvpar(ivpar),vvphi(ivphi),ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,&
                                    fs0(ivperp,ivpar,ivphi,is),fs1(ivperp,ivpar,ivphi,is))

                  !compute fs1 at adjacent locations in vperp1/vperp2 direction to take derivatives with later
                  !Note: delv may not be the best choice here when it is large. Consider using a separate variable to determine locations that we approximate derivative at
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))+delv
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat_new(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1_new(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,fs0_temp,fs1_plus_delvperp1(ivperp,ivpar,ivphi,is))
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))-delv
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat_new(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1_new(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,fs0_temp,fs1_minus_delvperp1(ivperp,ivpar,ivphi,is))
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))+delv
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat_new(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1_new(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,fs0_temp,fs1_plus_delvperp2(ivperp,ivpar,ivphi,is))
                  vperp1_adjacent = vvperp(ivperp)*COS(vvphi(ivphi))
                  vperp2_adjacent = vvperp(ivperp)*SIN(vvphi(ivphi))-delv
                  phi_adjacent = ATAN2(vperp2_adjacent,vperp1_adjacent) 
                  vperp_adjacent = SQRT(vperp1_adjacent**2+vperp2_adjacent**2)
                  fs0_temp=fs0hat_new(vperp_adjacent,vvpar(ivpar),hatV_s(is),spec(is)%alph_s)
                  call calc_fs1_new(omega,vperp_adjacent,vvpar(ivpar),phi_adjacent,ef,bf,hatV_s(is),spec(is)%q_s,spec(is)%alph_s,&
                                    spec(is)%tau_s,spec(is)%mu_s,spec(1)%alph_s,fs0_temp,fs1_minus_delvperp2(ivperp,ivpar,ivphi,is))
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

        !TODO: write out fs1

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
           enddo
          write(unit_s,*)
          write(unit_s+1,*)
        enddo

        close(unit_s)
        close(unit_s+1)

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


    end subroutine compute_fpc_gyro_new

    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    ! Determine dimensionless species equilibrium VDF \hat{fs0} at (vperp,vpar)
    real function fs0hat_new(vperp,vpar,hatV_s,aleph_s) 
      !input
      real         :: vperp, vpar      !normalized velocity space
      real         :: hatV_s           !normalized species drift velocity by wpar_s
      real         :: q_s              !normalized species charge
      real         :: aleph_s          !T_perp/T_parallel_s
      real         :: tau_s            !T_ref/T_s_parallel
      real         :: mu_s             !m_ref/m_s

      fs0hat_new = exp(-1.*(vpar-hatV_s)**2.-vperp**2./aleph_s)

    end function fs0hat_new

    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    ! Determine species perturbed VDF fs1 at (vperp,vpar)
    subroutine calc_fs1_new(omega,vperp,vpar,phi,ef,bf,hatV_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)
      use vars, only : betap,kperp,kpar,vtp,pi
      use vars, only : nbesmax
      use bessels, only : bessj_s, bess0_s_prime
      USE nrtype, only: SP, I4B

      complex, intent(in)   :: omega            !Complex Frequency
      real, intent(in)      :: vperp, vpar      !normalized velocity space
      real, intent(in)      :: phi              !azimuthal angle in velocity space
      complex, dimension(1:3), intent(in)   :: ef, bf           !E, B
      real, intent(in)      :: hatV_s           !Drift velocity norm to wpar_s
      real, intent(in)      :: q_s              !normalized species charge
      real, intent(in)      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)      :: mu_s             !m_ref/m_s
      real, intent(in)      :: aleph_r          !T_perp/T_parallel_R
      real, intent(in)      :: fs0                 !normalized zero order distribution
      complex, intent(out)                  :: fs1               !first order distribution

      integer :: n !bessel sum counter
      integer :: m !bessel sum counter
      complex :: ii= (0,1.) !Imaginary unit: 0+1i
      real :: b_s                           !Bessel function argument

      !intermediate values (see calculation notes for definitions)
      complex :: denom
      complex :: emult
      complex :: Wbar_s
      complex :: Ubar_s

      complex :: ef1,ef2,ef3
      complex :: omega_temp
      real :: phi_temp

      phi_temp = phi
      ef1 = ef(1) !Used to 'turn off' contributes due to ef1/2/3
      ef2 = ef(2)
      ef3 = ef(3)
      !ef1 = 0
      !ef2 = 0
      !ef3 = 0

      omega_temp = real(omega)+ii*aimag(omega) !`fix` sign as PLUME computes wr and gam such that omega=wr-i*gam but stores wr and gam in the var 'omega'


      !Compute Bessel functions (if necessary)
      ! NOTE: If bs is same as last time, no need to recompute Bessels
      b_s = (kperp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)
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
      Ubar_s= -2.*vperp/aleph_s*(1.+kpar*sqrt(mu_s/(tau_s*aleph_r))/omega_temp*((aleph_s-1)*vpar-aleph_s*hatV_s))

      write(*,*)'denoms (vpar,mu_s)', vpar, mu_s,omega_temp-kpar*vpar*sqrt(mu_s/(tau_s*aleph_r))&
      -mu_s/q_s,omega_temp-kpar*vpar*sqrt(mu_s/(tau_s*aleph_r)),omega_temp-kpar*vpar*sqrt(mu_s/(tau_s*aleph_r))+mu_s/q_s

      write(*,*)jbess(-1)*Ubar_s/b_s*ef1,jbess(0)*Ubar_s/b_s*ef1,jbess(1)*Ubar_s/b_s*ef1
      write(*,*)ii*0.5*(jbess(-1-1)-jbess(-1+1))*Ubar_s*ef2,&
      ii*0.5*bess0_s_prime(b_s)*Ubar_s*ef2,ii*0.5*(jbess(1-1)-jbess(1+1))*Ubar_s*ef2

      do n = -nbesmax,nbesmax
       !Calculate all parts of solution that don't depend on m
       denom=omega_temp-kpar*vpar*sqrt(mu_s/(tau_s*aleph_r))-n*mu_s/q_s
       Wbar_s=2.*(n*mu_s/(q_s*omega_temp)-1.)*(vpar-hatV_s) - 2.*(n*mu_s/(q_s*omega_temp*aleph_s))*vpar
       if (b_s .ne. 0.) then  !Handle division of first term if b_s=0 (U_bar_s also =0)
          emult=n*jbess(n)*Ubar_s/b_s*ef1
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
      fs1 = -1.*ii*sqrt(mu_s*tau_s/betap)/q_s*eperp1_bar*fs1*fs0
      
    end subroutine calc_fs1_new

    !============================================================================
    !============================================================================
    !============================================================================
    !------------------------------------------------------------------------------
    !                           Collin Brown and Greg Howes, 2023
    !------------------------------------------------------------------------------
    ! Verscharen 2016 formula for delta f eq A.2.12???
    ! Determine species perturbed VDF fs1 at (vperp,vpar)
    subroutine calc_fs1_new2(omega,vperp,vpar,phi,ef,bf,hatV_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)
      use vars, only : betap,kperp,kpar,vtp,pi
      use vars, only : nbesmax
      use bessels, only : bessj_s
      USE nrtype, only: SP, I4B

      complex, intent(in)   :: omega            !Complex Frequency
      real, intent(in)      :: vperp, vpar      !normalized velocity space
      real, intent(in)      :: phi              !azimuthal angle in velocity space
      complex, dimension(1:3), intent(in)   :: ef, bf           !E, B
      real, intent(in)      :: hatV_s           !Drift velocity norm to wpar_s
      real, intent(in)      :: q_s              !normalized species charge
      real, intent(in)      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)      :: mu_s             !m_ref/m_s
      real, intent(in)      :: aleph_r          !T_perp/T_parallel_R
      real, intent(in)      :: fs0                 !normalized zero order distribution
      complex, intent(out)                  :: fs1               !first order distribution

      integer :: n !bessel sum counter
      integer :: m !bessel sum counter
      complex :: ii= (0,1.) !Imaginary unit: 0+1i
      real :: b_s                           !Bessel function argument

      !intermediate values (see calculation notes for definitions)
      complex :: denom
      complex :: emult
      complex :: Vbar_s
      complex :: Ubar_s
      complex :: a_ns
                                                                                        
      !Compute Bessel functions (if necessary)
      ! NOTE: If bs is same as last time, no need to recompute Bessels
      b_s = (kperp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)
      if (b_s .ne. bs_last) then
         !Populate jbess with current values
         do n= -nbesmax-1,nbesmax+1
             jbess(n) = bessj_s(n,b_s)
         enddo
         bs_last=b_s 
      end if
      
      !Double Bessel Sum to calculate fs1=========================================
      fs1 = 0.
      !Calculate all parts of solution that don't depend on m or n
      Ubar_s= -2.*vperp/aleph_s*(1.+kpar*sqrt(mu_s/(tau_s*aleph_r))/omega * ( (aleph_s-1)*vpar - aleph_s*hatV_s ))
      Vbar_s= -2.*kperp*vperp*sqrt(mu_s/(tau_s*aleph_r))/(aleph_s*omega)* ( (aleph_s-1)*vpar - aleph_s*hatV_s )

      do n = -nbesmax,nbesmax
         a_ns=omega-kpar*vpar*sqrt(mu_s/(tau_s*aleph_r)) - n*mu_s/q_s
         !Calculate all parts of solution that don't depend on m
         denom=a_ns*a_ns-(mu_s/q_s)**2.
         emult=(ef(1)*Ubar_s-ef(3)*Vbar_s)*(ii*a_ns*cos(phi)+mu_s/q_s*sin(phi)) &
              +ef(2)*Ubar_s*(ii*a_ns*sin(phi)-mu_s/q_s*cos(phi)) &
              -ef(3)*ii*2*(vpar-hatV_s)/a_ns         
         do m = -nbesmax,nbesmax
            fs1=fs1+jbess(m)*jbess(n)*exp(ii*(m-n)*phi)*emult/denom
         enddo
      enddo
      fs1 = -1.*sqrt(mu_s*tau_s/betap)/q_s*eperp1_bar*fs1*fs0
      
    end subroutine calc_fs1_new2

    !============================================================================
    !============================================================================
    !============================================================================

    !------------------------------------------------------------------------------
    !                           Collin Brown, 2020
    !------------------------------------------------------------------------------
    subroutine write_fs0()
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName

      character(1000) :: filename                      !Output File name
      character(1000) :: outputPath                    !Output folder
      character(1000) :: cmd                           !Varaible to store command line commands
      real    :: vperpi, vpari                       !normalized velocity space current value in loop
      integer :: vperpindex, vparindex               !loop counters
      complex :: omega                               !Complex Frequency
      real    :: fs0                                 !normalized distribution function and first partial derivatives
      integer :: numstepvperp, numstepvpar           !total number of steps in loop
      logical :: ex                                  !used to check if results directory exists

      integer :: unit_s                 !out file unit counter

      integer :: is                     !species counter

      !check if results directory exists
      ! INQUIRE (DIRECTORY='data', EXIST=ex) !TODO: make this work for gfortran compiler
      ex = .true.
      if(ex) then
        write(*,*)"Assuming data folder already exists..."
      else
        write(*,*)"Creating data folder for output..."
        write(*,*)'mkdir data'
        CALL system('mkdir data') 
        write(*,*)"Saving output to data folder..."
      endif

      write(outputPath,*) 'data/', trim(dataName)
      ! INQUIRE (DIRECTORY=trim(dataName), EXIST=ex)
      ex = .true.
      if(ex) then
        write(*,*)"Assuming subfolder ", trim(dataName), " already exists"
      else
        write(*,*)"Creating data subfolder ", trim(dataName)
        write(cmd,*)'mkdir ',trim(outputPath)
        write(*,*)cmd
        CALL system(cmd)
        write(*,*)"Saving to data subfolder ", trim(dataName)
      endif

      !loop over each species and write to different files
      do is = 1, nspec
        !make file to store result
        !TODO: used "get unused unit" to get unit_s to pick correct 'number' to write to
        unit_s = 10+is !note: unit = 5,6 are reserved by standard fortran for input form keyboard/ writing to screen
        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.specie',(is),'.fs0'
        open(unit=unit_s,file=trim(filename),status='replace')

        write(*,*)'Calculating fs0 for species ',is

        !setup loop variables
        numstepvperp = int((vperpmax-vperpmin)/delv)
        numstepvpar = int((vparmax-vparmin)/delv)
        vperpi = vperpmin
        vpari = vparmin

        do vperpindex = 0, numstepvperp
          do vparindex = 0, numstepvpar
            call calc_fs0(vperpi,vpari,spec(is)%vv_s,spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,spec(is)%mu_s,fs0)
            if(ABS(real(fs0)) .lt. 9.999E-99) fs0 = (0.,0.) !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')real(fs0)
            vpari = vpari+delv
          end do
          vpari = vparmin
          vperpi = vperpi+delv
          write(unit_s,*)
        end do
        vpari = vparmin !reset counter
        vperpi = vperpmin !reset counter
      end do
    end subroutine write_fs0

    subroutine calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0) !TODO: remove q_s input as it is not used
      use vars, only : betap,vtp,pi
      !input
      real, intent(in)                      :: vperp, vpar      !normalized velocity space
      real, intent(in)                      :: V_s              !normalized species drift velocity by v_alfven
      real, intent(in)                      :: q_s              !normalized species charge
      real, intent(in)                      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)                      :: tau_s            !T_ref/T_s_parallel
      real, intent(in)                      :: mu_s             !m_ref/m_s

      !locals
      real                                  :: hatV_s           !normalized drift velocity by parallel thermal velocity

      !output
      real, intent(out)                  :: fs0              !zero order distribution

      
      hatV_s = V_s*sqrt(tau_s/(mu_s*betap))  !(GGH: Error: Fixed 13 JUL 2023)
      

      fs0 = EXP(-1.*(vpar-hatV_s)**2.-vperp**2./aleph_s)
    end subroutine calc_fs0

    subroutine calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)
      use vars, only : betap,kperp,kpar,vtp,pi
      use vars, only : nbesmax
      use bessels, only : bessj_s
      USE nrtype, only: SP, I4B

      complex, intent(in)                   :: omega            !Complex Frequency
      real, intent(in)                      :: vperp, vpar      !normalized velocity space
      real, intent(in)                      :: phi              !azimuthal angle in velocity space
      complex, dimension(1:3), intent(in)   :: ef, bf           !E, B
      real, intent(in)                      :: V_s              !normalized species drift velocity
      real, intent(in)                      :: q_s              !normalized species charge
      real, intent(in)                      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)                      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)                      :: mu_s             !m_ref/m_s
      real, intent(in)                      :: aleph_r          !T_perp/T_parallel_R
      real, intent(in)                      :: fs0              !normalized zero order distribution

      integer :: n !bessel sum counter
      integer :: m !bessel sum counter
      complex :: i !0+1i
      complex :: omega_temp !omega with correct gamma sign

      !intermediate values (see calculation notes for definitions)
      complex :: UExB
      complex :: A
      complex :: D
      complex :: Wbar_s
      complex :: Ubar_s
      real :: b_s
      complex :: fsi !intermediate variable used to compute summation total
      real, dimension(200) :: jbesselvals !array used to prevent repeated calculations 
                                                                                        

      real                                  :: hatV_s           !normalized drift velocity by parallel thermal velocity


      complex, intent(out)                  :: fs1               !first order distribution

      i = (0,1.)
      omega_temp = real(omega)+i*aimag(omega) !`fix` sign as PLUME computes wr and gam such that omega=wr-i*gam but stores wr and gam in the var 'omega'

      hatV_s = V_s*(tau_s/(mu_s*betap))**(.5)

      !UExB = -ef(1)/sqrt(bf(1)**2+bf(2)**2+bf(3)**2) !-Ex/|\mathbf{B}|

      !GGH ERROR 7 JUL 2023: ExB should use the equilibrium B, not the perturbed B!
!      UExB = -ef(1)/sqrt(abs(bf(1))**2+abs(bf(2))**2+abs(bf(3))**2) !-Ex/|\mathbf{B}|  !GGH: This line removed 13 JUL 2023 to correct amplitude calculation

      !UExB = -(ef(1)*bf(3)-ef(3)*bf(1))/(sqrt(bf(1)**2+bf(2)**2+bf(3)**2))

      ! UExB = (ef(2)*bf(3)-ef(3)*bf(2))**(2.)
      ! UExB = (ef(1)*bf(3)-ef(3)*bf(1))**(2.) + UExB
      ! UExB = (ef(1)*bf(2)-ef(2)*bf(1))**(2.) + UExB
      ! UExB = (UExB)**(0.5)
      ! UExB = -UExB/sqrt(abs(bf(1))**2+abs(bf(2))**2+abs(bf(3))**2)

      A=  mu_s/(q_s*sqrt(betap))* eperp1_bar !NEW version to agree with GH version
! OLD VERSION      
!      A = mu_s**1.5/(q_s*tau_s**0.5)*UExB/vtp
      Ubar_s = -2.*vperp*( (tau_s/mu_s)**(.5)/aleph_s+(kpar/(omega_temp*(aleph_r**0.5)))*(vpar-hatV_s-vpar/aleph_s))
      b_s = (kperp*q_s*vperp)/(mu_s*tau_s*aleph_r)**0.5

      !quick check to see if we take nbesmax to be large enough to approximate our infinite sum
      !note as a rule of thumb, bessj_s(n,x) is 'small' when both n<x and n is large
      !thus, if b_s is sufficiently less than n, we assume that nbesmax is large enough
      !we arbitrarily pick b_s < n/2
      !if(b_s .gt. real(n)/2.) write(*,*) 'Warning, chosen nbesmax may not be large enough to approximate the bessel expansion used for estimating fs1... b_s: ',b_s,' n: ',n,' real(n)/2.: ', real(n)/2.

      !populate jbesselvals to prevent redundant caculations
      n = -nbesmax-1 !need values for [-nbesmax-1,nbesmax+1]
      do while(n <= nbesmax+1)
        jbesselvals(n+1+nbesmax+1) = bessj_s(n,b_s) !here we use the 'unmodified' bessel function rather than the modified bessel function that 'bessel(n,x)' returns
                                                  !note as n starts at a negative value, we must offset our indexing
        n = n+1
      end do

      fs1 = 0.
      n = -nbesmax
      m = -nbesmax

      !Corrected and slightly optimized new version (GGH 13 JUL 2023)
      do while(n <= nbesmax)
           Wbar_s = -2.*((tau_s/mu_s)**(.5)*(vpar-hatV_s)-(n/omega_temp)*(((mu_s*tau_s)**(0.5))/q_s)*(vpar-hatV_s-vpar/aleph_s))
           !GGH ERROR: Should be aleph_r below (FIXED!)
           D = (omega_temp-n*mu_s/q_s-kpar*(mu_s/(aleph_r*tau_s))**(.5)*vpar)
           
          fsi = (0.,0.)
          fsi = (jbesselvals(n+1+nbesmax+1))*Wbar_s*ef(3) !here we use the 'unmodified' bessel function rather than the modified bessel function that 'bessel(n,x)' returns
          fsi = i*((jbesselvals(n+nbesmax+1)-jbesselvals(n+2+nbesmax+1))/2.)*Ubar_s*ef(2)+fsi
          fsi = ( n*(jbesselvals(n+1+nbesmax+1)/b_s)*Ubar_s*ef(1)+fsi)/D
        do while(m <= nbesmax)
          fs1 = fs1 + jbesselvals(m+1+nbesmax+1)*(cexp(i*(m-n)*phi))*fsi
          m = m + 1
        end do
        m = -nbesmax
        n = n + 1
      end do
      fs1 = -1.*i*A*fs0*fs1
    end subroutine calc_fs1


    !TODO: the calc_correlation routines are very similar, and should be combined
    !      That is there should only be 1 calc_correlation_*_gyro func
    subroutine calc_correlation_par_gyro(omega,ef,bf,vperpin,vparin,delv,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,Cor_s,fs1,dfs1z)
      use vars, only : betap,vtp,pi
      !input
      complex, intent(in)                   :: omega            !Complex Frequency
      real, intent(in)                      :: vperpin, vparin  !normalized velocity space current value input
      real, intent(in)                      :: delv             !normalized velocity grid point spacing
      real, intent(in)                      :: V_s              !normalized species drift velocity
      real, intent(in)                      :: q_s              !normalized species charge
      real, intent(in)                      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)                      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)                      :: mu_s             !m_ref/m_s
      real, intent(in)                      :: aleph_r          !T_perp/T_parallel_R
      complex, dimension(1:3),intent(in)    :: ef, bf           !E, B fields from eigen funcs

      !locals
      real                          :: fs0              !zero-order distribution function * c^3 (times c^3 due to normalization)
      complex                       :: fs1_1            !numerical derivation var
      complex                       :: fs1_0            !numerical derivation var
      real                          :: piconst          !3.14159...
      real                          :: vperp, vpar      !normalized velocity space current value in loop

      real                          :: phi              !azimuthal angle (we integrate over this)
      real                          :: delphi
      integer                       :: num_phi          !number over integration sample points
      integer                       :: n_phi            !integration counter

      !output
      real, intent(out)             :: Cor_s            !normalized correlation value
      complex, intent(inout)          :: fs1              !first order distribution * c^3 (times c^3 due to normalization)
      complex, intent(out)          :: dfs1z            !first order distribution partial wrt z * c^4 B_0 (times c^4 B_0 due to normalization)

      num_phi = 20

      piconst = 4.0*ATAN(1.0)
      phi = 0.
      delphi = 2.0*piconst/num_phi
      n_phi = 0
      Cor_s = 0.

      vperp = vperpin
      vpar = vparin !TODO: remove redundant variables
      do while(n_phi < num_phi)
        !simple finite central difference method for derivative
        call calc_fs0(vperp,vpar+delv,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar+delv,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
        call calc_fs0(vperp,vpar-delv,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar-delv,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)
        call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)

        dfs1z = (fs1_1-fs1_0)/(2.*delv)
        Cor_s = vperp*0.5*(-1.*q_s*(vpar**2./2.)*dfs1z*CONJG(ef(3))&
          -1.*q_s*(vpar**2./2.)*CONJG(dfs1z)*ef(3))*2.0*piconst/num_phi+Cor_s

        n_phi = n_phi+1
        phi = phi + delphi
      end do
    end subroutine calc_correlation_par_gyro

    subroutine calc_correlation_perp_gyro(omega,ef,bf,vperpin,vparin,delv,&
                                          V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,Cor_perp_s,fs1,dfs1perp)
      use vars, only : betap,vtp,pi
      !input
      complex, intent(in)                   :: omega            !Complex Frequency
      real, intent(in)                      :: vperpin, vparin  !normalized velocity space input
      real, intent(in)                      :: delv             !normalized velocity grid point spacing
      real, intent(in)                      :: V_s              !normalized species drift velocity
      real, intent(in)                      :: q_s              !normalized species charge
      real, intent(in)                      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)                      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)                      :: mu_s             !m_ref/m_s
      real, intent(in)                      :: aleph_r          !T_perp/T_parallel_R
      complex, dimension(1:3),intent(in)    :: ef, bf           !E, B fields from eigen funcs

      !locals
      real                          :: fs0              !zero-order distribution function * c^3 (times c^3 due to normalization)
      complex                       :: fs1_1            !numerical derivation var
      complex                       :: fs1_0            !numerical derivation var
      real                          :: piconst          !3.14159...
      real                          :: vperp, vpar      !normalized velocity space current value in loop

      real                          :: phi              !azimuthal angle (we integrate over this)
      real                          :: delphi
      integer                       :: num_phi          !number over integration sample points
      integer                       :: n_phi            !integration counter

      real                          :: vxtemp,vytemp,vztemp  !temp vars to convert between gyro and cartesian velocity space when taking derivative of f1 (dist func)
      real                          :: vpartemp,phitemp,vperptemp !temp vars used to convert back from cart to gyro to help compute derivative in the correct direction

      real                          :: Cor_ex_s,Cor_ey_s!intermediate correlation values

      !output
      real, intent(out)             :: Cor_perp_s       !total normalized correlation value
      complex, intent(inout)          :: fs1              !first order distribution * c^3 (times c^3 due to normalization)
      complex, intent(out)          :: dfs1perp         !first order distribution partial wrt z * c^4 B_0 (times c^4 B_0 due to normalization)

      num_phi = 30 !TODO: dont hard code this

      vperp = vperpin
      vpar = vparin
 
      piconst = 4.0*ATAN(1.0)
      phi = 0.
      delphi = 2.0*piconst/num_phi
      n_phi = 0
      Cor_perp_s = 0.
      do while(n_phi < num_phi)
        !TODO: derive analytical form when calculating location to compute f1 for computing df1
        !simple finite central difference method for derivative
        vxtemp = vperp*COS(phi)+delv
        vytemp = vperp*SIN(phi)
        vztemp = vpar
        phitemp = ATAN2(vytemp,vxtemp) 
        vperptemp = SQRT((vxtemp)**2+(vytemp)**2)
        vpartemp = vztemp
        call calc_fs0(vperptemp,vpartemp,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperptemp,vpartemp,phitemp,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
        vxtemp = vperp*COS(phi)-delv
        vytemp = vperp*SIN(phi) 
        vztemp = vpar
        phitemp = ATAN2(vytemp,vxtemp) 
        vperptemp = SQRT((vxtemp)**2+(vytemp)**2)
        vpartemp = vztemp
        call calc_fs0(vperptemp,vpartemp,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperptemp,vpartemp,phitemp,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)
        dfs1perp = (fs1_1-fs1_0)/(2.*delv)
        vxtemp = vperp*COS(phi)
        Cor_ex_s = 0.5*(-1.*q_s*(vxtemp**2./2.)*dfs1perp*CONJG(ef(1))&
          -1.*q_s*(vxtemp**2./2.)*CONJG(dfs1perp)*ef(1))

        !simple finite central difference method for derivative
        vxtemp = vperp*COS(phi) 
        vytemp = (vperp*SIN(phi)+delv)
        vztemp = vpar
        phitemp = ATAN2(vytemp,vxtemp) 
        vperptemp = SQRT((vxtemp)**2+(vytemp)**2)
        vpartemp = vztemp
        call calc_fs0(vperptemp,vpartemp,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperptemp,vpartemp,phitemp,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
        vxtemp = vperp*COS(phi)
        vytemp = (vperp*SIN(phi)-delv)
        vztemp = vpar
        phitemp = ATAN2(vytemp,vxtemp) 
        vperptemp = SQRT((vxtemp)**2+(vytemp)**2)
        vpartemp = vztemp
        call calc_fs0(vperptemp,vpartemp,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperptemp,vpartemp,phitemp,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)
        dfs1perp = (fs1_1-fs1_0)/(2.*delv) !TODO: use different variable name for dfs1perp here for clarity
        vytemp = vperp*SIN(phi)
        Cor_ey_s = 0.5*(-1.*q_s*(vytemp**2./2.)*dfs1perp*CONJG(ef(2))&
          -1.*q_s*(vytemp**2./2.)*CONJG(dfs1perp)*ef(2))

        Cor_perp_s = vperp*(Cor_ex_s+Cor_ey_s)*2.0*piconst/num_phi + Cor_perp_s

        n_phi = n_phi+1
        phi = phi + delphi
      end do
    end subroutine calc_correlation_perp_gyro

    !TODO: account for projection when computing normalization
    subroutine calc_correlation_i_car(omega,ef,bf,vxin,vyin,vzin,vmax3rdval,vmax3rdindex,ceiindex,delv, &
                                      V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,Cor_i_s,fs1)
      !calculates projection of CEi(vx,vy,vz) (int CEi(vx,vy,vz) dxi) where vmax3rdindex species which direction the rountine integrates over

      use vars, only : betap,vtp,pi
      !input
      complex, intent(in)                   :: omega            !Complex Frequency
      complex, dimension(1:3),intent(in)    :: ef, bf           !E, B fields from eigen funcs
      real, intent(in)                      :: vxin, vyin, vzin !normalized velocity space current value in loop
      real, intent(in)                      :: vmax3rdval       !integration or range (-vmax3rdindex,+vmax3rdindex)
      integer, intent(in)                   :: vmax3rdindex     !integration dir index (vx<->1, vx<->2, vz<->3)
      integer, intent(in)                   :: ceiindex         !select which val to compute (cex=cperp1<->1,cey=cperp2<->2,cez=cpar<->3)
      real, intent(in)                      :: delv             !normalized velocity grid point spacing
      real, intent(in)                      :: V_s              !normalized species drift velocity
      real, intent(in)                      :: q_s              !normalized species charge
      real, intent(in)                      :: aleph_s          !T_perp/T_parallel_s
      real, intent(in)                      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)                      :: mu_s             !m_ref/m_s
      real, intent(in)                      :: aleph_r          !T_perp/T_parallel_R

      !locals
      real                          :: vx, vy, vz       !velocity coordinates

      real                          :: fs0              !zero-order distribution function * c^3 (times c^3 due to normalization)
      complex                       :: sum_fs1              !perturbed distri
      complex, intent(out)                       :: fs1              !perturbed distribution function 
      complex                       :: fs1_1            !numerical derivation var
      complex                       :: fs1_0            !numerical derivation var
      complex                       :: dfs1i            !derivative in direction controlled by vmax3rdindex
      real                          :: piconst          !3.14159...

      integer                       :: intmax           !max number of integration sample points
      integer                       :: inti             !integration counter

      real                          :: vperp,vpar,phi   !gyro coordinates vars
      real                          :: vival            !integration coordinate value tracker (specified by vmax3rdindex)
      real                          :: vicor            !temp var that holds velocity term withinn FPC. Is defined by ceiindex


      !output
      real, intent(out)             :: Cor_i_s        !total normalized correlation value

      vx = vxin
      vy = vyin
      vz = vzin

      intmax = int((vmax3rdval+vmax3rdval)/delv)
      inti = 0
      vival = -1*vmax3rdval

      piconst = 4.0*ATAN(1.0)
      Cor_i_s = 0.

      sum_fs1=0.
      do while(inti < intmax)
        if(vmax3rdindex == 1)then 
          vx = vival
        endif
        if(vmax3rdindex == 2)then 
          vy = vival
        endif
        if(vmax3rdindex == 3)then 
          vz = vival
        endif

        !simple finite central difference method for derivative
        if(ceiindex == 1)then 
          phi = ATAN2(vy,vx+delv) 
          vperp = SQRT((vx+delv)**2+vy**2)
          vpar = vz
          call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
          call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
          phi = ATAN2(vy,vx-delv) 
          vperp = SQRT((vx-delv)**2+vy**2)
          vpar = vz
          call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
          call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)
          dfs1i = (fs1_1-fs1_0)/(2.*delv)
        endif
        if(ceiindex == 2)then 
          phi = ATAN2(vy+delv,vx) 
          vperp = SQRT(vx**2+(vy+delv)**2)
          vpar = vz
          call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
          call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
          phi = ATAN2(vy-delv,vx) 
          vperp = SQRT(vx**2+(vy-delv)**2)
          vpar = vz
          call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
          call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)
          dfs1i = (fs1_1-fs1_0)/(2.*delv) 
        endif
        if(ceiindex == 3)then 
          phi = ATAN2(vy,vx) 
          vperp = SQRT(vx**2+vy**2)
          vpar = vz+delv
          call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
          call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
          phi = ATAN2(vy,vx) 
          vperp = SQRT(vx**2+vy**2)
          vpar = vz-delv
          call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
          call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)
          dfs1i = (fs1_1-fs1_0)/(2.*delv)
        endif

        phi = ATAN2(vy,vx) 
        vperp = SQRT(vx**2+vy**2)
        vpar = vz
        call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)

        if(ceiindex == 1)then 
          vicor = vx
        endif
        if(ceiindex == 2)then 
          vicor = vy
        endif
        if(ceiindex == 3)then 
          vicor = vz
        endif

        !Cor_i_s = (-1.*q_s*(vicor**2./2.)*dfs1i*ef(ceiindex))*delv+Cor_i_s
        Cor_i_s = (q_s*(vicor**2./2.)*dfs1i*ef(ceiindex))*delv+Cor_i_s

        sum_fs1=sum_fs1+fs1*delv        
        inti = inti + 1
        vival = vival + delv
     end do
     fs1=sum_fs1
   end subroutine calc_correlation_i_car

   subroutine check_nbesmax(vperpmax,tau_s,mu_s,aleph_r)
      !As j_n(b) is small for b<n/2, nbesmax/2 should be greater than or equal to b_s,max = |(kperp*q_s*vperp)/sqrt(mu_s*tau_s*aleph_r)| for all species
      !TODO: if using wrapper, the user will not see this error- find a way to warn user if using wrapper (or atleast leave note to check output!)
      !TODO: alternatively, automatically set nbesmax to be some multiple of the minimum value calculated using this formula...
      use vars, only : kperp
      use vars, only : nbesmax

      real, intent(in)                      :: vperpmax         !max value to be computed in normalized velocity space
      real, intent(in)                      :: tau_s            !T_ref/T_s|_parallel
      real, intent(in)                      :: mu_s             !m_ref/m_s
      real, intent(in)                      :: aleph_r          !T_perp/T_parallel_R

      real                                  :: b_s              !argument to bessel functions (see calc_fs1)


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
