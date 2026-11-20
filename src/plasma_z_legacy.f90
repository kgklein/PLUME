module plasma_z_legacy
  implicit none
  private
  public :: zetout_legacy
contains
   pure complex function zetout_legacy(zin)
     !! This is the subroutine used by Linsker to evaluate his Z-function.
     complex :: z,dzetaz,term,fmult,terme,an1,bn1,zsquar
     complex :: hold,temp1,temp2,ddzeta,dddzet,zetaoz
     complex,intent(in) :: zin
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
99   if (x*x+y*y > 9.) go to 200
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
401  zetout_legacy=zetaoz
9999 continue

   end function zetout_legacy
end module
