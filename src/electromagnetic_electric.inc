! Historical Ex=1 cofactor elimination; shared without edited duplicates.
  electric(1) = cmplx(1.,0.)
  electric(3)=-electric(1)*(lam(2,1)*lam(3,2)-lam(3,1)*lam(2,2))
  electric(3)= electric(3)/(lam(2,3)*lam(3,2)-lam(3,3)*lam(2,2))
  electric(2) = -electric(3)*lam(3,3) - electric(1)*lam(3,1)
  electric(2) = electric(2)/lam(3,2)

