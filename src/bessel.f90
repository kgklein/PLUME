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

 contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
   FUNCTION bessim0(x)
     !! Calculates the Modified Bessel Function of order zero.
     !! Determines I_0(x) e^(-x) instead of I_0(x) to large argument problems.
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
     !! Determines I_1(x) e^(-x) instead of I_1(x) to large argument problems.
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
     !! Determines I_n(x) e^(-x) instead of I_n(x) to large argument problems.
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
 end module bessels
