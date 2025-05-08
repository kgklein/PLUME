!=============================================================================!
!=============================================================================!
!!*PLUME                                                                    *!!
!!Plasma in a Linear Uniform Magnetized Environment                          !!
!!                                                                           !!
!!Kristopher Klein                                                           !!
!!kgklein@arizona.edu                                                        !!
!!Lunar and Planetary Laboratory, University of Arizona
!!                                                                           !!
!                             Bessel Functions
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2005
!                                Greg Howes
!------------------------------------------------------------------------------
!  From Numerical Recipes
!------------------------------------------------------------------------------
module bessels
  !!Calculates Bessel Functions necessary for the dispersion calculation.
  implicit none
  private

  public :: bessim0,bessim1,bessim
  public :: bessj0, bessj1
  public :: bessj_s,bess0_s_prime

 contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
   FUNCTION bessim0(x)
     !! Calculates the Modified Bessel Function of order zero.
     !! Determines I_0(x) e^(-x) instead of I_0(x)
     !! to avoid large argument problems.
     USE nrtype; USE nrutil_trim, ONLY : poly
     IMPLICIT NONE
     
     REAL(SP), INTENT(IN) :: x
     !! Argument of Modified Bessel Function.
     
     REAL(SP) :: bessim0
     !! Returned value of the Bessel Function calculation.
     
     REAL(SP) :: ax
     !! Absolute value of the argument.
     
     REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
          3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
          0.45813e-2_dp/)
     !!Parameters for calculation.
     
     REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
          0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
          -0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
          0.392377e-2_dp/)
     !!Parameters for calculation.
     
     ax=abs(x)
     if (ax < 3.75) then
        bessim0=poly(real((x/3.75_sp)**2,dp),p)*exp(-ax)
     else
        bessim0=(1./sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
     end if
     
   END FUNCTION bessim0
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
   FUNCTION bessim1(x)
     !! Calculates the Modified Bessel Function of order one.
     !! Determines I_1(x) e^(-x) instead of I_1(x) to
     !! avoid large argument problems.
     USE nrtype; USE nrutil_trim, ONLY : poly
     IMPLICIT NONE

     REAL(SP), INTENT(IN) :: x
     !! Argument of Modified Bessel Function.

     REAL(SP) :: bessim1
     !! Returned value of the Bessel Function calculation.

     REAL(SP) :: ax
     !! Absolute value of the argument.

     REAL(DP), DIMENSION(7) :: p = (/0.5_dp,0.87890594_dp,&
          0.51498869_dp,0.15084934_dp,0.2658733e-1_dp,&
          0.301532e-2_dp,0.32411e-3_dp/)
     !!Parameters for calculation.

     REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,-0.3988024e-1_dp,&
          -0.362018e-2_dp,0.163801e-2_dp,-0.1031555e-1_dp,&
          0.2282967e-1_dp,-0.2895312e-1_dp,0.1787654e-1_dp,&
          -0.420059e-2_dp/)
     !!Parameters for calculation.

     ax=abs(x)
     if (ax < 3.75) then
        bessim1=ax*poly(real((x/3.75_sp)**2,dp),p)*exp(-ax)
     else
        bessim1=(1./sqrt(ax))*poly(real(3.75_sp/ax,dp),q)
     end if
     if (x < 0.0) bessim1=-bessim1
     
   END FUNCTION bessim1
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
   FUNCTION bessim(n,x)
     !! Calculates the Modified Bessel Function of arbitrary order.
     !! Determines I_n(x) e^(-x) instead of I_n(x) to
     !! avoid large argument problems.
     USE nrtype; USE nrutil_trim, ONLY : assert
     IMPLICIT NONE
     
     INTEGER(I4B), INTENT(IN) :: n
     !!Order of Modified Bessel Function.
     
     REAL(SP), INTENT(IN) :: x
     !!Order of Modified Bessel Function.
     
     REAL(SP) :: bessim
     !! Returned value of the Bessel Function calculation.
     
     INTEGER(I4B), PARAMETER :: IACC=40,IEXP=maxexponent(x)/2
     !! Parameters for Calculation.
     
     INTEGER(I4B) :: j,m
     !! Index for Calculation.
     
     REAL(SP) :: bi,bim,bip,tox
     !! Dummy Variables for Calculation.
     
     call assert(n >= 2, 'bessim args')
     bessim=0.0
     if (x*x <= 8.0_sp*tiny(x)) RETURN
     tox=2.0_sp/abs(x)
     bip=0.0
     bi=1.0
     m=2*((n+int(sqrt(real(IACC*n,sp)))))
     do j=m,1,-1
        bim=bip+j*tox*bi
        bip=bi
        bi=bim
        if (exponent(bi) > IEXP) then
           bessim=scale(bessim,-IEXP)
           bi=scale(bi,-IEXP)
           bip=scale(bip,-IEXP)
        end if
        if (j == n) bessim=bip
     end do
     bessim=bessim*bessim0(x)/bi
     if (x < 0.0 .and. mod(n,2) == 1) bessim=-bessim
   END FUNCTION bessim
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
   FUNCTION bessj0(x)
     USE nrtype; USE nrutil_trim, ONLY : poly
     IMPLICIT NONE
     REAL(SP), INTENT(IN) :: x
     REAL(SP) :: bessj0
     REAL(SP) :: ax,xx,z
     REAL(DP) :: y
     REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
          0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
     REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
          0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
          -0.934945152e-7_dp/)
     REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
          651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
          -184.9052456_dp/)
     REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
          9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
     if (abs(x) < 8.0) then
        y=x**2
        bessj0=poly(y,r)/poly(y,s)
     else
        ax=abs(x)
        z=8.0_sp/ax
        y=z**2
        xx=ax-0.785398164_sp
        bessj0=sqrt(0.636619772_sp/ax)*(cos(xx)*&
             poly(y,p)-z*sin(xx)*poly(y,q))
     end if
   END FUNCTION bessj0
!------------------------------------------------------------------------------
!                           Greg Howes, 2005 & Collin Brown 2022
!------------------------------------------------------------------------------
   FUNCTION bessj1(x)
     USE nrtype; USE nrutil_trim, ONLY : poly
     IMPLICIT NONE
     REAL(SP), INTENT(IN) :: x
     REAL(SP) :: bessj1
     REAL(SP) :: ax,xx,z
     REAL(DP) :: y
     REAL(DP), DIMENSION(6) :: r = (/72362614232.0_dp,&
          -7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
          15704.48260_dp,-30.16036606_dp/)
     REAL(DP), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
          18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
     REAL(DP), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
          -0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
     REAL(DP), DIMENSION(5) :: q = (/0.04687499995_dp,&
          -0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
          0.105787412e-6_dp/)
     if (abs(x) < 8.0) then
        y=x**2
        bessj1=x*(poly(y,r)/poly(y,s))
     else
        ax=abs(x)
        z=8.0_sp/ax
        y=z**2
        xx=ax-2.356194491_sp
        bessj1=sqrt(0.636619772_sp/ax)*(cos(xx)*&
             poly(y,p)-z*sin(xx)*poly(y,q))*sign(1.0_sp,x)
     end if
   END FUNCTION bessj1

   FUNCTION bessj_s(n, x)
      USE nrtype; USE nrutil_trim, ONLY : assert
      INTEGER(I4B) :: n
      REAL(SP), INTENT(IN) :: x
      REAL(SP) :: bessj_s
      INTEGER(I4B), PARAMETER :: IACC=40, IEXP=maxexponent(x)/2
      INTEGER(I4B) :: j, jsum, m
      REAL(DP) :: ax, bj, bjm, bjp, summ, tox
      LOGICAL :: n_is_negative !fix to handle negative n's

      if(n < 0) then
       n_is_negative = .true.
       n = n*(-1)
      else
       n_is_negative = .false.
      end if

      ax = abs(x)
      if(n == 0) then
        bessj_s = bessj0(REAL(ax,SP))
      else if (n == 1) then
        bessj_s = bessj1(REAL(ax,SP))
      else if (ax*ax <= 8.0D0*tiny(x)) then
         bessj_s = 0D0
      else if (ax > real(n, DP)) then
         tox = 2.0D0/ax
         bjm = bessj0(REAL(ax,SP))
         bj = bessj1(REAL(ax,SP))
         do j = 1, n-1
            bjp = j*tox*bj-bjm
            bjm = bj
            bj = bjp
         end do
         bessj_s = bj
      else
         tox = 2D0/ax
         m = 2*((n+int(sqrt(real(IACC*n, DP))))/2)
         bessj_s = 0D0
         jsum = 0
         summ = 0D0
         bjp = 0D0
         bj = 1D0
         do j = m, 1, -1
            bjm = j*tox*bj-bjp
            bjp = bj
            bj = bjm
            if (exponent(bj) > IEXP) then
               bj = scale(bj, -IEXP)
               bjp = scale(bjp, -IEXP)
               bessj_s = scale(bessj_s, -IEXP)
               summ = scale(summ,-IEXP)
            end if
            if (jsum /= 0)  summ = summ+bj
            jsum = 1 - jsum
            if (j == n)     bessj_s = bjp
         end do
         summ = 2D0*summ-bj
         bessj_s = bessj_s/summ
      end if
      if (x < 0D0 .and. mod(n, 2) == 1 .and. .not. n_is_negative) bessj_s = -bessj_s
      if (x > 0D0 .and. mod(n, 2) == 1 .and. n_is_negative) bessj_s = -bessj_s
      if (n_is_negative) n = n*(-1) !return n to being negative. 
   END FUNCTION bessj_s

   FUNCTION bess0_s_prime(x)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL :: bess0_s_prime
      integer :: n 

      n = 1
      bess0_s_prime = -bessj_s(n,x)
   END FUNCTION bess0_s_prime
 end module bessels
