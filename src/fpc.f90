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
  private :: calc_correlation_par, calc_fs0, calc_fs1 

  public :: compute_fpc, write_fs0

  contains
    
    !------------------------------------------------------------------------------
    !                           Collin Brown, 2020
    !------------------------------------------------------------------------------
    subroutine compute_fpc(wrootindex) !TODO: move this and all related routines to separate f90 file
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName
      
      use disprels, only : calc_eigen, rtsec, disp

      integer, intent(in) :: wrootindex              !index of selected root

      character(60) :: filename                      !Output File name
      character(60) :: outputPath                    !Output folder
      character(60) :: cmd                           !Varaible to store command line commands
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
      !>>>KGK: 1/31/23; allow GGH's new LD/TTD calculation
      real, dimension(1:6,1:nspec) :: Ps_split_new !Power into/out of species (GGH)
      !<<<KGG: 1/31/23
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
        write(*,*)"assuming subfolder ", trim(dataName), "alreay exists"
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
        !make file to store result
        !TODO: used "get unused unit" to get unit_s to pick correct 'number' to write to
        !TODO: update sample input file to include all input varaibles (see vars.f90. Ex: we are missing values for tensor_s, thus a 'random' logical is assigned here)
        unit_s = 10+is !note: unit = 5,6 are reserved by standard fortran for input form keyboard/ writing to screen
        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cpar.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s,file=trim(filename),status='replace')

        write(filename,'(5A,I0.2,1A,I0.2)')'data/',trim(dataName),'/',trim(outputName),'.cperp.specie',(is),'.mode',wrootindex !Assumes nspec,nroots < 100 for filename formating
        open(unit=unit_s+1,file=trim(filename),status='replace')

        write(*,*)'Calculating fpc for species ',is
        write(*,*)'Writing omega/kpar V_a normalization to file. WIP'

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
            call calc_correlation_par(omega,ef,bf,vperpi,vpari,delv,spec(is)%vv_s,&
                                      spec(is)%Q_s,spec(is)%alph_s,spec(is)%tau_s,&
                                      spec(is)%mu_s,spec(1)%alph_s,Cor_par_s,fs1,dfs1z)
            call calc_correlation_perp(omega,ef,bf,vperpi,vpari,delv,spec(is)%vv_s,&
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

    end subroutine compute_fpc

    !------------------------------------------------------------------------------
    !                           Collin Brown, 2020
    !------------------------------------------------------------------------------
    subroutine write_fs0()
      use vars, only : betap,kperp,kpar,vtp,nspec,spec
      use vars, only : vperpmin,vperpmax,vparmin,vparmax,delv
      use vars, only : wroots, nroots
      use vars, only : outputName, dataName

      character(60) :: filename                      !Output File name
      character(60) :: outputPath                    !Output folder
      character(60) :: cmd                           !Varaible to store command line commands
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
        write(*,*)"Assuming subfolder ", trim(dataName), "alreay exists"
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
        !TODO: update sample input file to include all input varaibles (see vars.f90. Ex: we are missing values for tensor_s, thus a 'random' logical is assigned here)
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
      real, dimension(1000) :: jbesselvals !array used to prevent repeated calculations 
                                                                                        

      real                                  :: hatV_s           !normalized drift velocity by parallel thermal velocity


      complex, intent(out)                  :: fs1               !first order distribution

      !real :: start, finish !debug/test to measure runtime of function

      i = (0,1)
      omega_temp = real(omega)+i*aimag(omega) !`fix` sign as people PLUME returns omega as omega=wr-i*gam

      i = (0,1.)
      hatV_s = V_s*(tau_s/(mu_s*betap))**(.5)

      UExB = -ef(1)/sqrt(bf(1)**2+bf(2)**2+bf(3)**2) !-Ex/|\mathbf{B}|
      A = mu_s**1.5/(q_s*tau_s**0.5)*UExB/vtp
      Ubar_s = -2.*vperp*((tau_s/mu_s)**(.5)/aleph_s+(kpar/(omega_temp*(aleph_r**0.5)))*(vpar-hatV_s-vpar/aleph_s))
      b_s = (kperp*q_s*vperp)/(mu_s*tau_s*aleph_s)**0.5

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
                                                  !not sure how to index this jbesselvals array 'cleanly'. Just note the potential for off by one errors as we need 2*nbesmax+2 values here
        n = n+1
      end do

      fs1 = 0.
      n = -nbesmax
      m = -nbesmax

      do while(n <= nbesmax)
        do while(m <= nbesmax)
          Wbar_s = -2.*((tau_s/mu_s)**(.5)*(vpar-hatV_s)-(n/omega_temp)*(((mu_s*tau_s)**(0.5))/q_s)*(vpar-hatV_s-vpar/aleph_s))
          D = (omega_temp-n*mu_s/q_s-kpar*(mu_s/(aleph_s*tau_s))**(.5)*vpar)
          fsi = 0.
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

    subroutine calc_correlation_par(omega,ef,bf,vperp,vpar,delv,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,Cor_s,fs1,dfs1z)
      use vars, only : betap,vtp,pi
      !input
      complex, intent(in)                   :: omega            !Complex Frequency
      real, intent(in)                      :: vperp, vpar      !normalized velocity space current value in loop
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

      real                          :: phi              !azimuthal angle (we integrate over this)
      real                          :: delphi
      integer                       :: num_phi          !number over integration sample points
      integer                       :: n_phi            !integration counter

      !output
      real, intent(out)             :: Cor_s            !normalized correlation value
      complex, intent(out)          :: fs1              !first order distribution * c^3 (times c^3 due to normalization)
      complex, intent(out)          :: dfs1z            !first order distribution partial wrt z * c^4 B_0 (times c^4 B_0 due to normalization)

      num_phi = 20 !TODO: dont hard code this

      piconst = 4.0*ATAN(1.0)
      phi = 0.
      delphi = 2.0*piconst/num_phi
      n_phi = 0
      Cor_s = 0.
      do while(n_phi < num_phi)
        !simple finite central difference method for derivative
        call calc_fs0(vperp,vpar+delv,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar+delv,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
        call calc_fs0(vperp,vpar-delv,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar-delv,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)

        call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)

        dfs1z = (fs1_1-fs1_0)/(2.*delv)
        Cor_s = vperp*-1.*q_s*(vpar**2./2.)*dfs1z*ef(3)*2.0*piconst/num_phi+Cor_s !Note: this implicitly is a real cast due to the type of Cor_s
        n_phi = n_phi+1
        phi = phi + delphi
      end do
    end subroutine calc_correlation_par

    subroutine calc_correlation_perp(omega,ef,bf,vperp,vpar,delv,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,Cor_perp_s,fs1,dfs1perp)
      use vars, only : betap,vtp,pi
      !input
      complex, intent(in)                   :: omega            !Complex Frequency
      real, intent(in)                      :: vperp, vpar      !normalized velocity space current value in loop
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

      real                          :: phi              !azimuthal angle (we integrate over this)
      real                          :: delphi
      integer                       :: num_phi          !number over integration sample points
      integer                       :: n_phi            !integration counter

      real                          :: Cor_ex_s,Cor_ey_s!intermediate correlation values

      !output
      real, intent(out)             :: Cor_perp_s       !total normalized correlation value
      complex, intent(out)          :: fs1              !first order distribution * c^3 (times c^3 due to normalization)
      complex, intent(out)          :: dfs1perp         !first order distribution partial wrt z * c^4 B_0 (times c^4 B_0 due to normalization)

      num_phi = 20 !TODO: dont hard code this

      piconst = 4.0*ATAN(1.0)
      phi = 0.
      delphi = 2.0*piconst/num_phi
      n_phi = 0
      Cor_perp_s = 0.
      do while(n_phi < num_phi)
        !simple finite central difference method for derivative
        call calc_fs0(vperp+delv,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp+delv,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_1)
        call calc_fs0(vperp-delv,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp-delv,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1_0)

        call calc_fs0(vperp,vpar,V_s,q_s,aleph_s,tau_s,mu_s,fs0)
        call calc_fs1(omega,vperp,vpar,phi,ef,bf,V_s,q_s,aleph_s,tau_s,mu_s,aleph_r,fs0,fs1)

        !
        dfs1perp = (fs1_1-fs1_0)/(2.*delv)
        Cor_ex_s = -1.*q_s*(vperp**2./2.)*dfs1perp*ef(1)
        Cor_ey_s = -1.*q_s*(vperp**2./2.)*dfs1perp*ef(2)
        Cor_perp_s = vperp*(Cor_ex_s+Cor_ey_s)*2.0*piconst/num_phi + Cor_perp_s

        n_phi = n_phi+1
        phi = phi + delphi
      end do
    end subroutine calc_correlation_perp
    
end module fpc
