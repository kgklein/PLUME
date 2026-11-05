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
        !here, we add i nu_ns/Omega_ref to the numerator
        tsi(n) = sqrt((alphp * disp_tau)/(disp_mu))*&
             ( om - Vdrifts - dble(n)* disp_mu/disp_Q + c_i * disp_nuns)/kpar
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
        !Minimal change; addition of (alph_s-1)*i*\nu_ns
        !-=-=-=-=-=-=-=-=-=-=
        eps_xx_t  = eps_xx_t  + &
             n*n*jn(n,is)*( 2.*(disp_alph-1.) &
             + (sqrt(alphp * disp_tau/ disp_mu )/ kpar )*(&
             (disp_alph* (om - Vdrifts) + (disp_alph-1.d0)*c_i*disp_nuns)* &
             (zz(n)+zz(-n)) &
             + (zz(n)-zz(-n))*(1.-disp_alph)*(n*disp_mu/(disp_Q ))))

        !-=-=-=-=-=-=-=-=-=-=
        !YY
        !Minimal change; addition of (alph_s-1)*i*\nu_ns
        !-=-=-=-=-=-=-=-=-=-=
        eps_yy_t = eps_yy_t + &
             (n*n*jn(n,is) + 2.*lambdas*lambdas*(jn(n,is)-jpn(n,is)))*&
             (2.*(disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
             (disp_alph*(om - Vdrifts)+(disp_alph-1.d0)*c_i*disp_nuns)*&
             (zz(n)+zz(-n)) + (n*disp_mu/(disp_Q))*&
             (1.-disp_alph)*(zz(n)-zz(-n))))

        !-=-=-=-=-=-=-=-=-=-=
        !ZZ
        !Paired longitudinal Krook response; see tests/COLLISIONAL_ZZ_CORRECTION.md.
        !-=-=-=-=-=-=-=-=-=-=
        eps_zz_t = eps_zz_t + &
             jn(n,is)*(2.*om* &
             (om*disp_alph - Vdrifts+c_i *disp_nuns*disp_alph) - &
             2.*n*n*disp_mu**2.* (1.-disp_alph)/(disp_Q**2.) &
             +(sqrt(alphp*disp_tau/disp_mu)/kpar) * (&
             (zz(n)+zz(-n))*&
             (om*(om+c_i*disp_nuns)*&
             (disp_alph * (om - Vdrifts) +c_i*disp_nuns*disp_alph)&
             + n*n*disp_mu**2.*&
             (om*(3.*disp_alph-2.) -disp_alph*Vdrifts+&
             2.d0*c_i*disp_nuns*(disp_alph-1.d0))/disp_Q**2.) &
             +(n * disp_mu/disp_Q )*(zz(n)-zz(-n))*&
             (om*om*(1.d0-3.d0*disp_alph) +&
             (2.d0*om + c_i *disp_nuns)*disp_alph*Vdrifts+&
             (disp_alph-1.d0)*(disp_nuns*disp_nuns-2.d0*c_i*disp_nuns*om)&
             -2.d0*c_i*disp_nuns*om*disp_alph+&
             n*n*disp_mu**2.*(1.-disp_alph)/( disp_Q**2.)) ) )

        !-=-=-=-=-=-=-=-=-=-=
        !XY
        !Minimal change; addition of (alph_s-1)*i*\nu_ns
        !-=-=-=-=-=-=-=-=-=-=
        eps_xy_t = eps_xy_t + &
             n*(jn(n,is)-jpn(n,is))*( &
             (zz(n)-zz(-n))*(disp_alph*(om-Vdrifts) +(disp_alph-1)*c_i*disp_nuns)&
             +(zz(n)+zz(-n))*(1.-disp_alph)*n*disp_mu/(disp_Q ))

        !-=-=-=-=-=-=-=-=-=-=
        !XZ
        !Moderate changes; addition of several (alph_s-1)*i*\nu_ns terms
        !-=-=-=-=-=-=-=-=-=-=
        eps_xz_t = eps_xz_t + &
             n*jn(n,is)*(2.*n*(1.-disp_alph)*disp_mu/(disp_Q) + &
             (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
             (zz(n) -zz(-n))*&
             ((om+c_i*disp_nuns)*((om-Vdrifts)*disp_alph+c_i *disp_nuns*(disp_alph-1.d0) )-&
             (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) + &
             (zz(n) +zz(-n))*&
             (disp_mu*n*(om-2.*om*disp_alph+Vdrifts*disp_alph&
             +2.d0*c_i*disp_nuns*(1.d0-disp_alph))/disp_Q) ))

        !-=-=-=-=-=-=-=-=-=-=
        !YZ
        !Moderate changes; addition of several (alph_s-1)*i*\nu_ns terms
        !-=-=-=-=-=-=-=-=-=-=
        eps_yz_t = eps_yz_t + &
             (jn(n,is)-jpn(n,is))*(2.*(om*disp_alph-Vdrifts+c_i*disp_nuns*(disp_alph-1.d0)) &
             +(sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
             (zz(n)+zz(-n))*&
             ((om+c_i*disp_nuns)*(disp_alph*(om-Vdrifts) +c_i*disp_nuns*(disp_alph-1.d0) )  -&
             (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph) ) &
             +(zz(n)-zz(-n))*&
             (n*disp_mu*(om-disp_alph*(2.*om-Vdrifts) +&
             2.d0*c_i*disp_nuns*(1.d0-disp_alph))/disp_Q) ))

        !-=-=-=-=-=-=-=-=-=-=
        !n=\pm 1 susceptibility
        if ((low_n).and.(n==1)) then

           !\chi_xx
           susc_low(is,1,1,1) = norm(1)*&
                n*n*jn(n,is)*( (disp_alph-1.) &
                + (sqrt(alphp * disp_tau/ disp_mu )/ kpar )*(&
                (disp_alph* (om - Vdrifts) + (disp_alph-1.d0)*c_i*disp_nuns)* &
                (zz(n)) &
                + (zz(n))*(1.-disp_alph)*(n*disp_mu/(disp_Q ))))

           susc_low(is,1,1,-1) = norm(1)*&
                n*n*jn(n,is)*( (disp_alph-1.) &
                + (sqrt(alphp * disp_tau/ disp_mu )/ kpar )*(&
                (disp_alph* (om - Vdrifts) + (disp_alph-1.d0)*c_i*disp_nuns)* &
                (zz(-n)) &
                + (-zz(-n))*(1.-disp_alph)*(n*disp_mu/(disp_Q ))))

           !\chi_xy
           susc_low(is,1,2,1) = norm(4)*&
                n*(jn(n,is)-jpn(n,is))*( &
                (zz(n))*(disp_alph*(om-Vdrifts) +(disp_alph-1.d0)*c_i*disp_nuns)&
                +(zz(n))*(1.-disp_alph)*n*disp_mu/(disp_Q ))

           susc_low(is,1,2,-1) = norm(4)*&
                n*(jn(n,is)-jpn(n,is))*( &
                (-zz(-n))*(disp_alph*(om-Vdrifts) +(disp_alph-1.d0)*c_i*disp_nuns) &
                +(zz(-n))*(1.-disp_alph)*n*disp_mu/(disp_Q ))

           !\chi_yx
           susc_low(is,2,1,1) = -susc_low(is,1,2,1)

           susc_low(is,2,1,-1) = -susc_low(is,1,2,-1)

           !\chi_xz
           susc_low(is,1,3,1) = norm(5)*&
                n*jn(n,is)*(n*(1.-disp_alph)*disp_mu/(disp_Q) + &
                (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (zz(n))*&
                ((om+c_i*disp_nuns)*((om-Vdrifts)*disp_alph+c_i *disp_nuns*(disp_alph-1.d0) ) -&
                (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) + &
                (zz(n))*&
                (disp_mu*n*(om-2.*om*disp_alph+Vdrifts*disp_alph&
                +2.d0*c_i*disp_nuns*(1.d0-disp_alph))/disp_Q) ))

           susc_low(is,1,3,-1) = norm(5)*&
                n*jn(n,is)*(n*(1.-disp_alph)*disp_mu/(disp_Q) + &
                (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (-zz(-n))*&
                ((om+c_i*disp_nuns)*((om-Vdrifts)*disp_alph+c_i *disp_nuns*(disp_alph-1.d0) ) -&
                (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph)) + &
                (zz(-n))*&
                (disp_mu*n*(om-2.*om*disp_alph+Vdrifts*disp_alph&
                +2.d0*c_i*disp_nuns*(1.d0-disp_alph))/disp_Q) ))

           !\chi_zx
           susc_low(is,3,1,1) = susc_low(is,1,3,1)

           susc_low(is,3,1,-1) = susc_low(is,1,3,-1)

           !\chi_yy
           susc_low(is,2,2,1) = norm(2)*&
                (n*n*jn(n,is) + 2.*lambdas*lambdas*(jn(n,is)-jpn(n,is)))*&
                ((disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (disp_alph*(om - Vdrifts) + (disp_alph-1.d0)*c_i*disp_nuns)*(zz(n)) + &
                (n*disp_mu/(disp_Q))*&
                (1.-disp_alph)*(zz(n))))

           susc_low(is,2,2,-1) = norm(2)*&
                (n*n*jn(n,is) + 2.*lambdas*lambdas*(jn(n,is)-jpn(n,is)))*&
                ((disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (disp_alph*(om - Vdrifts) + (disp_alph-1.d0)*c_i*disp_nuns)*(zz(-n)) + &
                (n*disp_mu/(disp_Q))*&
                (1.-disp_alph)*(-zz(-n))))

           !\chi_yz
           susc_low(is,2,3,1) = norm(6)*&
                (jn(n,is)-jpn(n,is))*((om*disp_alph-Vdrifts) &
                +(sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (zz(n))*&
                ((om+c_i*disp_nuns)*(disp_alph*(om-Vdrifts) +c_i*disp_nuns*(disp_alph-1.d0) )  -&
                (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph) ) &
                +(zz(n))*&
                (n*disp_mu*(om-disp_alph*(2.*om-Vdrifts) +&
                2.d0*c_i*disp_nuns*(1.d0-disp_alph))/disp_Q) ))

           susc_low(is,2,3,-1) = norm(6)*&
                (jn(n,is)-jpn(n,is))*((om*disp_alph-Vdrifts) &
                +(sqrt(alphp*disp_tau/disp_mu)/kpar)*(&
                (zz(-n))*&
                ((om+c_i*disp_nuns)*(disp_alph*(om-Vdrifts) +c_i*disp_nuns*(disp_alph-1.d0) )  -&
                (n*n*disp_mu**2./(disp_Q**2.))*(1.-disp_alph) ) &
                +(-zz(-n))*&
                (n*disp_mu*(om-disp_alph*(2.*om-Vdrifts) +&
                2.d0*c_i*disp_nuns*(1.d0-disp_alph))/disp_Q) ))

           !\chi_zy
           susc_low(is,3,2,1) = -susc_low(is,2,3,1)

           susc_low(is,3,2,-1) = -susc_low(is,2,3,-1)

           ! Longitudinal signed harmonics include their full nonresonant
           ! terms. The bulk sum retains its separate drift correction.
           block
              integer :: zz_harmonic
              real :: zz_gyro
              complex :: zz_shift,zz_resonance,zz_coefficient
              zz_shift=om+c_i*disp_nuns
              do zz_harmonic=-1,1,2
                 zz_gyro=zz_harmonic*n*disp_mu/disp_Q
                 zz_resonance=zz_shift-zz_gyro
                 zz_coefficient=disp_alph*om+(1.-disp_alph)*zz_gyro
                 susc_low(is,3,3,zz_harmonic)=norm(3)*jn(n,is)*(&
                      zz_coefficient*zz_resonance+zz_gyro*Vdrifts+&
                      (sqrt(alphp*disp_tau/disp_mu)/kpar)*zz(zz_harmonic*n)*(&
                      zz_coefficient*zz_resonance**2+&
                      disp_alph*(zz_gyro-om)*Vdrifts*zz_resonance))
              enddo
           end block

        endif

     enddo !end n loop

     !Add in the n=zero term

     !No n=0 XX term

     !No n=0 XY,YX terms

     !No n=0 XZ,ZX terms

     !YY term
     eps_yy_t = eps_yy_t + &
          (2.*lambdas*lambdas*(jn(0,is)-jpn(0,is)))*&
          ((disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*&
          (disp_alph*(om-Vdrifts)+(disp_alph-1.d0)*c_i*disp_nuns)*(zz(0)))

     !YZ term
     eps_yz_t = eps_yz_t + &
          (jn(0,is)-jpn(0,is))*(om*disp_alph-Vdrifts+c_i*disp_nuns*(disp_alph-1.d0)+&
          zz(0)*(sqrt(alphp*disp_tau/disp_mu)/kpar)*&
          (disp_alph*(om-Vdrifts)+c_i*disp_nuns*(disp_alph-1.d0))*(om+c_i*disp_nuns))

     !ZZ term
     eps_zz_t = eps_zz_t + &
          jn(0,is)*om*(&
          om*disp_alph-Vdrifts +c_i*disp_nuns*disp_alph+&
          zz(0) * (om +c_i*disp_nuns)* (sqrt(alphp*disp_tau/disp_mu)/kpar) &
          *(disp_alph*(om-Vdrifts)+c_i*disp_nuns*disp_alph) )


     if (low_n) then

        !xx
        susc_low(is,1,1,0)=cmplx(0.d0,0.d0)

        !yy
        susc_low(is,2,2,0) =  (2.*lambdas*lambdas*(jn(0,is)-jpn(0,is)))*&
             ((disp_alph-1.) + (sqrt(alphp*disp_tau/disp_mu)/kpar)*&
             (disp_alph*(om-Vdrifts)+(disp_alph-1.d0)*c_i*disp_nuns)*(zz(0))) * norm(2)

        !xy
        susc_low(is,1,2,0)=cmplx(0.d0,0.d0)

        !yx
        susc_low(is,2,1,0)=cmplx(0.d0,0.d0)

        !xz
        susc_low(is,1,3,0)=cmplx(0.d0,0.d0)

        !zx
        susc_low(is,3,1,0)=cmplx(0.d0,0.d0)

        !yz
        susc_low(is,2,3,0) = &
             (jn(0,is)-jpn(0,is))*(om*disp_alph-Vdrifts+c_i*disp_nuns*(disp_alph-1.d0)+&
             zz(0)*(sqrt(alphp*disp_tau/disp_mu)/kpar)*&
             (disp_alph*(om-Vdrifts)+c_i*disp_nuns*(disp_alph-1.d0))*(om+c_i*disp_nuns))*&
             norm(6)

        !zy
        susc_low(is,3,2,0) = -susc_low(is,2,3,0)

        ! Full n=0 longitudinal harmonic; no separate drift term is added
        ! to this diagnostic. Unlike the bulk regrouping, it contains no -w*V.
        susc_low(is,3,3,0)=norm(3)*jn(0,is)*disp_alph*om*(om+c_i*disp_nuns)*(&
             1.+zz(0)*(sqrt(alphp*disp_tau/disp_mu)/kpar)*(om+c_i*disp_nuns-Vdrifts))


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

