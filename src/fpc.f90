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
!!collin.crbrown@gmail.com or collbrown@uiowa.edu                            !!
!!University of Iowa                                                         !!
!=============================================================================!
!=============================================================================!
module fpc
  implicit none
  private :: calc_correlation_par_gyro, calc_correlation_perp_gyro, calc_correlation_i_car, calc_fs0, calc_fs1

  public :: compute_fpc_gyro, compute_fpc_cart, write_fs0

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

      ! Refine Omega Value (again)
      iflag=0
      omega=rtsec(disp,om1,om2,tol,iflag)
      
      call calc_eigen(omega,ef,bf,Us,ns,Ps,Ps_split,Ps_split_new,.true.,.true.)
      ef(1) = ef(1)
      ef(2) = ef(2) !note: 'coord trans': The routine used to calculate the eigen functions (calc_eigen) and the routine used to calculate fs0/fs1/CEi(i.e. FPC related functions) use subtly different coordinates, so we fix them here
      ef(3) = ef(3) !   In field aligned coordinates, Epar (aka Ez) is fixed in the direction of the guiding magnetic field, but Eperp1/Eperp2 (aka Ex/Ey) have freedom in what direction they can be in
                    !   In this code, we use two different sources that take different coordinates, and we must account for them here.
                    !TODO: move this note after fixing code^^^^^
                    
      do is = 1, nspec
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
      ef(1) = ef(1)
      ef(2) = ef(2)
      ef(3) = ef(3)

      ! !DEBUG TEST: hard code needed quants
      ! ef(1) = (1,0)
      ! ef(2) = (2,0)
      ! ef(3) = (0,0)
      ! omega = (.2,-.2)


      do is = 1, nspec
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

        !DEBUG
        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.debugcperp1.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+3,file=trim(filename),status='replace')
        write(unit_s+3,'(16A)')'CEperp1 vx vy vz'
        !END DEBUG

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

        !setup loop variables
        numstepvx = int((vxmax-vxmin)/delv)
        numstepvy = int((vymax-vymin)/delv)
        numstepvz = int((vzmax-vzmin)/delv)

        !CEi(vx,vy)----------------------------------------------------------------------------
        vxi = vxmin
        vyi = vymin
        vzi = vzmin

        do vxindex = 0, numstepvx
          do vyindex = 0, numstepvy
            vmax3rdval = vzmax
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,3,3,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,3,1,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp1_s)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,3,2,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp2_s)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp1_s) .lt. 9.999E-99) Cor_perp1_s = 0. !file formating bug fix
            if(ABS(Cor_perp2_s) .lt. 9.999E-99) Cor_perp2_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp1_s
            write(unit_s+2,'(es17.5)',advance='no')Cor_perp2_s

            !DEBUG
            write(unit_s+3,'(es17.5,es17.5,es17.5,es17.5)')Cor_perp1_s,vxi,vyi,vzi
            !END DEBUG

            vyi = vyi+delv
          end do
          vyi = vymin
          vxi = vxi+delv
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
        end do
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'

        !CEi(vx,vz)----------------------------------------------------------------------------
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        do vxindex = 0, numstepvx
          do vzindex = 0, numstepvz
            vmax3rdval = vymax
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,2,3,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,2,1,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp1_s)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,2,2,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp2_s)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp1_s) .lt. 9.999E-99) Cor_perp1_s = 0. !file formating bug fix
            if(ABS(Cor_perp2_s) .lt. 9.999E-99) Cor_perp2_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp1_s
            write(unit_s+2,'(es17.5)',advance='no')Cor_perp2_s
            vzi = vzi+delv
          end do
          vzi = vzmin
          vxi = vxi+delv
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
        end do    
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'

        !CEi(vy,vz)----------------------------------------------------------------------------
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        do vyindex = 0, numstepvy
          do vzindex = 0, numstepvz
            vmax3rdval = vxmax
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,1,3,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,1,1,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp1_s)
            call calc_correlation_i_car(omega,ef,bf,vxi,vyi,vzi,vmax3rdval,1,2,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_perp2_s)
            if(ABS(Cor_par_s) .lt. 9.999E-99) Cor_par_s = 0. !file formating bug fix
            if(ABS(Cor_perp1_s) .lt. 9.999E-99) Cor_perp1_s = 0. !file formating bug fix
            if(ABS(Cor_perp2_s) .lt. 9.999E-99) Cor_perp2_s = 0. !file formating bug fix
            write(unit_s,'(es17.5)',advance='no')Cor_par_s
            write(unit_s+1,'(es17.5)',advance='no')Cor_perp1_s
            write(unit_s+2,'(es17.5)',advance='no')Cor_perp2_s
            vzi = vzi+delv
          end do
          vzi = vzmin
          vyi = vyi+delv
          write(unit_s,*)
          write(unit_s+1,*)
          write(unit_s+2,*)
        end do
        vxi = vxmin
        vyi = vymin
        vzi = vzmin
        write(unit_s,*)'---'
        write(unit_s+1,*)'---'
        write(unit_s+2,*)'---'

      end do

    end subroutine compute_fpc_cart

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
      complex :: fs0                                 !normalized distribution function and first partial derivatives
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
          write(*,*)'vperp = ',vperpi !quick debug 'progress bar'
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
      complex, intent(out)                  :: fs0              !zero order distribution

      hatV_s = V_s*(tau_s/(mu_s*betap))

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
      complex, intent(in)                   :: fs0              !normalized zero order distribution

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

      UExB = -ef(1)/sqrt(abs(bf(1))**2+abs(bf(2))**2+abs(bf(3))**2) !-Ex/|\mathbf{B}|

      !UExB = -(ef(1)*bf(3)-ef(3)*bf(1))/(sqrt(bf(1)**2+bf(2)**2+bf(3)**2))

      ! UExB = (ef(2)*bf(3)-ef(3)*bf(2))**(2.)
      ! UExB = (ef(1)*bf(3)-ef(3)*bf(1))**(2.) + UExB
      ! UExB = (ef(1)*bf(2)-ef(2)*bf(1))**(2.) + UExB
      ! UExB = (UExB)**(0.5)
      ! UExB = -UExB/sqrt(abs(bf(1))**2+abs(bf(2))**2+abs(bf(3))**2)

      A = mu_s**1.5/(q_s*tau_s**0.5)*UExB/vtp
      Ubar_s = -2.*vperp*((tau_s/mu_s)**(.5)/aleph_s+(kpar/(omega_temp*(aleph_r**0.5)))*(vpar-hatV_s-vpar/aleph_s))
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

      do while(n <= nbesmax)
        do while(m <= nbesmax)
          Wbar_s = -2.*((tau_s/mu_s)**(.5)*(vpar-hatV_s)-(n/omega_temp)*(((mu_s*tau_s)**(0.5))/q_s)*(vpar-hatV_s-vpar/aleph_s))
          D = (omega_temp-n*mu_s/q_s-kpar*(mu_s/(aleph_s*tau_s))**(.5)*vpar)
          fsi = (0.,0.)
          fsi = (jbesselvals(n+1+nbesmax+1))*Wbar_s*ef(3) !here we use the 'unmodified' bessel function rather than the modified bessel function that 'bessel(n,x)' returns
          fsi = i*((jbesselvals(n+nbesmax+1)-jbesselvals(n+2+nbesmax+1))/2.)*Ubar_s*ef(2)+fsi
          fsi = n*(jbesselvals(n+1+nbesmax+1)/b_s)*Ubar_s*ef(1)+fsi
          fsi = jbesselvals(m+1+nbesmax+1)*(cexp(i*(m-n)*phi)/D)*fsi
          fs1 = fs1 + fsi
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
      complex                       :: fs0              !zero-order distribution function * c^3 (times c^3 due to normalization)
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
      complex, intent(out)          :: fs1              !first order distribution * c^3 (times c^3 due to normalization)
      complex, intent(out)          :: dfs1z            !first order distribution partial wrt z * c^4 B_0 (times c^4 B_0 due to normalization)

      num_phi = 20

      piconst = 4.0*ATAN(1.0)
      phi = 0.
      delphi = 2.0*piconst/num_phi
      n_phi = 0
      Cor_s = 0.

      vperp = vperpin
      vpar = -vparin !see note 'coord trans'
      do while(n_phi < num_phi)
        !simple finite central difference method for derivative
        call calc_fs0(vperp,vpar+delv,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar+delv,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
        call calc_fs0(vperp,vpar-delv,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar-delv,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)

        call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0) !TODO: update how fs0/ fs1 out is computed!
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
      complex                       :: fs0              !zero-order distribution function * c^3 (times c^3 due to normalization)
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
      complex, intent(out)          :: fs1              !first order distribution * c^3 (times c^3 due to normalization)
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
          -1.*q_s*(vperp**2./2.)*CONJG(dfs1perp)*ef(1))

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
          -1.*q_s*(vperp**2./2.)*CONJG(dfs1perp)*ef(2))

        Cor_perp_s = vperp*(Cor_ex_s+Cor_ey_s)*2.0*piconst/num_phi + Cor_perp_s

        n_phi = n_phi+1
        phi = phi + delphi
      end do
    end subroutine calc_correlation_perp_gyro

    !TODO: account for projection when computing normalization
    subroutine calc_correlation_i_car(omega,ef,bf,vxin,vyin,vzin,vmax3rdval,vmax3rdindex,ceiindex,delv, &
                                      V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,Cor_i_s)
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

      complex                       :: fs0              !zero-order distribution function * c^3 (times c^3 due to normalization)
      complex                       :: fs1              !perturbed distribution function 
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
        Cor_i_s = 0.5*(-1.*q_s*(vicor**2./2.)*dfs1i*CONJG(ef(ceiindex))&
          -1.*q_s*(vicor**2./2.)*CONJG(dfs1i)*ef(ceiindex))*delv+Cor_i_s

        inti = inti + 1
        vival = vival + delv
      end do
    end subroutine calc_correlation_i_car
    
end module fpc
