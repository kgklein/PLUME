!=============================================================================!
!=============================================================================!
!                Gyrokinetic Dispersion Relation Functions
!=============================================================================!
!=============================================================================!
!
!------------------------------------------------------------------------------
!                             Copyright,  2005
!                                Greg Howes
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
module gk_disp_gf
   implicit none
   private

   real, target :: bi                              !Ion Plasma Beta
   real, target :: bipar                              !Ion Parallel Plasma Beta
   real, target :: tr                              !Temp Ratio = Te/Ti
   real, target :: kp                              !k_perp rho_i
   real, target :: alphi                           !(Tperp / Tpar)_ion
   real, target :: alphe                           !(Tperp / Tpar)_electron
   real, target :: kapi, kape                             !kappa
   real    :: krat                                 !k_par/k_perp
   real    :: vti                                  !Ion Thermal Velocity
   real    :: nuii, nuee                            !Like-particle collisionality
   real    :: nui, nue                              !Ion/electron collisionality
   real    :: nuh                                  !Hyperviscosity
   real    :: nuh_exp                              !Hyperviscous exponent (del^(nuh_exp))
   real    :: mr                                   !Mass Ratio mi/me
   logical :: krookcoll                            !T=Krook operator
   logical :: collisions                           !T=collisions on
   logical :: hypervisc                            !T=Hyperviscosity on
   logical :: fail
   integer :: choice                               !Choice for dispersion relation
   integer, parameter :: standard = 0                !Collisionless Inviscid
   integer, parameter :: collisional = 1             !Collisional
   integer, parameter :: hyperviscous = 2            !Hyperviscous
   integer, parameter :: krook = 3                   !Krook Collisions
   integer, parameter :: bimkap = 4                   !Bi-Maxwellian/Bi-kappa
   real, parameter :: c = 1.0                        !Speed of light (c=1)
   complex, parameter :: ii = (0., 1.)                !i
   real, parameter :: pi = 3.1415927                 !Argument k_perp rho_i
   complex :: ta, tc, td                             !Intermediate calculations
   complex :: tap, tcp, tcpp, tdp, tf                  !Intermediate calculations
   real    :: tb, te, tz                             !Intermediate calculations
   complex :: tai, tae, tci, tce, tdi, tde           !Intermediate calculations
   real    :: tbi, tbe, tei, tee                   !Intermediate calculations
   real :: alphi_tw                                !Auxilliary parameter A14f in Kunz 2015
   real :: alphi_twi, alphi_twe                              !Auxilliary parameter A14f in Kunz 2015
   real    :: taa                                  !New intermediates
   complex :: tbb                                  !New intermediates
   integer, parameter :: ne = 50                     !Number of v-space energies
   integer, parameter :: nl = 50                     !Number of v-space pitch angles
   real, parameter :: ecut = 64.                   !Max energy in v-space

   public :: bi, bipar, tr, alphi, alphe, kp, vti, mr, pi, krat, kapi, kape
   public :: nui, nue, collisions, col_disp
   public :: nuii, nuee, krook_disp, krookcoll
   public :: nuh, nuh_exp, hypervisc, hv_disp
   public :: choice, standard, collisional, hyperviscous, krook, bimkap
   public :: disp, dispbimkap, dispanabikap
   public :: disp_prime
   public :: epbp
   public :: heating
   public :: get_eigenfuncs, get_eigenfuncs_bimkap
   public :: fail
   public :: dist, ne, nl
   public :: kpara

contains
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the Gyrokinetic Dispersion relation as a function of frequency om
   complex function disp(om)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex :: om                                !Complex Frequency
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments

      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kp*kp*tr/mr

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)

      !Find Bessel Functions
      g0i = bessim0(xi)
      g1i = bessim0(xi) - bessim1(xi)
      g2i = 2.*g1i
      g0e = bessim0(xe)
      g1e = bessim0(xe) - bessim1(xe)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions
      call zetvec(xii, zi, fail)
      if (fail) return
      call zetvec(xie, ze, fail)
      if (fail) return

!     zi=Z_o(xii)
!     ze=Z_o(xie)

      !Calculate intermediate values
      ta = 1 + g0i*xii*zi + (1.+g0e*xie*ze)/tr
      tb = 1.-g0i + (1.-g0e)/tr
      tc = g1i*xii*zi - g1e*xie*ze
      td = g2i*xii*zi + g2e*xie*ze*tr
      te = g1i - g1e

      !Calculate Value of Dispersion relation D(om)
!     disp=(1./(om*om) - tb/xi*(1.-tb/ta) )*(1.-0.5*bi*td+0.5*bi*tc*tc/ta) -&
!          0.5*bi/xi*(te+tb*tc/ta)**2.
      !Rewritten for better numerical behavior
!     disp=(xi*ta/(om*om) - ta*tb + tb*tb)*(2.*ta/bi - ta*td + tc*tc) - &
!          (ta*te + tb*tc)**2.
      disp = (xi*ta - ta*tb*(om*om) + tb*tb*(om*om))*(2.*ta/bi - ta*td + tc*tc) - &
             (om*om)*(ta*te + tb*tc)**2.

!     disp=(xi*ta/(om*om) - ta*tb + tb*tb)*(2./bi-td)

      !DEBUG
      if (0 .eq. 1) then
         write (*, '(5es14.6)') om, disp, abs(disp)
      end if

      return
   end function disp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Jason Tenbarge 2016
   !Analytical ki-kappa solution for k_perp rho_i, me / mi << 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   complex function dispanabikap(om)
      use gk_funcs, only: modzetvec
      complex :: om
      complex :: xii, z
      real    :: Tpi, Tpe                           !(Tperp/T)_s
      real    :: Tbi, Tbe                           !(Tpar/T)_s
      complex :: C0bi, C1bi, C2bi
      real :: C0be, kap
      complex :: A, B, C

      Tpi = 3.*alphi/(2.*alphi + 1.)
      Tpe = 3.*alphe/(2.*alphe + 1.)
      Tbi = 3./(2.*alphi + 1.)
      Tbe = 3./(2.*alphe + 1.)

      xii = om/sqrt(bi*Tbi)
      z = modzetvec(xii, kap)
      z = z*2.*kap/(2.*kap - 1.)
      C0be = (2.*kap - 1.)/(2.*kap - 3.)
      C0bi = C0be*(1.+xii*z)
      C1bi = 1./(2.*kap) + (1.-1./(2.*kap))*(1.+xii*xii/kap)*(1.+xii*z)
      C2bi = (1.-(3./(4.*kap)) + (1.-3./(2.*kap))*xii*xii/(2.*kap))*(1./(kap - 1.)) + &
             (1.+xii*xii/kap)*(1.+xii*xii/kap)*(1.-3./(2.*kap))*(1.-1./(2.*kap))*(kap/(kap - 1.))* &
             (1.+xii*z)

      A = C0bi*2./(bi*Tbi) + C0be*2./(bi*Tbe*tr)
      B = bi*Tpi*(C2bi*alphi - 1.) + bi*Tpe*tr*(alphe - 1.) - 1.
      C = C1bi*alphi - alphe

      dispanabikap = A*B - C*C

   end function dispanabikap
!------------------------------------------------------------------------------
!                           Jason TenBarge, 2016
!------------------------------------------------------------------------------
!  Calculates the Gyrokinetic Dispersion relation as a function of frequency om
   complex function dispbimkap(om)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs
      use DEIntegrator, only: recursiveintegrate
      implicit none
      complex :: om                                !Complex Frequency
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: kpe                               !k rho_e
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (electron)
      real    :: g00i, g00pi, g10pi, g20pi            !Perp velocity integrals (ion)
      complex :: g00bi, g10bi, g20bi               !Parallel velocity integrals (ion)
      real    :: g00e, g00pe, g10pe, g20pe            !Perp velocity integrals (elec)
      complex :: g00be, g10be, g20be               !Parallel velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments
      real    :: C0pi, C0pe                               !C0perp
      real    :: Tperpi, Tperpe                           !(Tperp/T)_s
      real    :: Tpari, Tpare                           !(Tpar/T)_s

      real :: x0 = 0., xcut = 100.                 !Double exponential v_perp / v_th_perp integral bounds
      integer :: evals
      real :: error, repart, impart

      Tperpi = 3.*alphi/(2.*alphi + 1.)
      Tperpe = 3.*alphe/(2.*alphe + 1.)
      Tpari = 3./(2.*alphi + 1.)
      Tpare = 3./(2.*alphe + 1.)

      kpe = kp*sqrt((Tperpe/Tperpi)*tr/mr)
      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kpe*kpe

      !Calculate arguments of plasma dispersion function
      !Omega / sqrt(beta_parallel_s)
      xii = om/sqrt(bi*Tpari)
      xie = om/sqrt(bi*Tpare*mr*tr)

      if (kapi .ne. 1.) then
         xii = xii*sqrt(2.*(kapi + 1.)/(2.*kapi - 3.))

         C0pi = (2.*kapi - 1.)/(2.*kapi - 3.)

         !Branch necessary because the g00 and g00p integrals converge to 1+eps for k rho_i << 1, leading to highly
         !unphysical results
         if (kp < 0.1) then
            g00i = 1.-xi
            g00pi = C0pi - xi
            g10pi = 1.-1.5*xi
            g20pi = 2.*(1.-1.5*xi/C0pi)
         else
            call recursiveintegrate(func00, x0, xcut, kp, xii, kapi, g00i, evals, error, .true.)
            call recursiveintegrate(func00perp, x0, xcut, kp, xii, kapi, g00pi, evals, error, .true.)
            call recursiveintegrate(func10perp, x0, xcut, kp, xii, kapi, g10pi, evals, error, .true.)
            call recursiveintegrate(func20perp, x0, xcut, kp, xii, kapi, g20pi, evals, error, .true.)
         end if

         call recursiveintegrate(refunc00par, 1., xcut + 1., kp, xii, kapi, repart, evals, error, .true.)
         call recursiveintegrate(imfunc00par, 1., xcut + 1., kp, xii, kapi, impart, evals, error, .true.)
         g00bi = cmplx(repart, impart)
         call recursiveintegrate(refunc10par, 1., xcut + 1., kp, xii, kapi, repart, evals, error, .true.)
         call recursiveintegrate(imfunc10par, 1., xcut + 1., kp, xii, kapi, impart, evals, error, .true.)
         g10bi = cmplx(repart, impart)
         call recursiveintegrate(refunc20par, 1., xcut + 1., kp, xii, kapi, repart, evals, error, .true.)
         call recursiveintegrate(imfunc20par, 1., xcut + 1., kp, xii, kapi, impart, evals, error, .true.)
         g20bi = cmplx(repart, impart)

         tai = C0pi - g00pi + g00bi*alphi
         tbi = (C0pi - g00pi)
         tci = alphi*g10bi - g10pi
         tdi = alphi*g20bi - g20pi
         tei = g10pi
         alphi_twi = xi + 0.5*bi*Tpari*(alphi - 1.)*(1.-g00i)

      else
         !Find Bessel Functions for bimax
         g0i = bessim0(xi)
         g1i = bessim0(xi) - bessim1(xi)
         g2i = 2.*g1i

         !Calculate Plasma Dispersion functions
         call zetvec(xii, zi, fail)
         if (fail) then
            !write(*,*) 'Warning! Ion disp function failed!'
            fail = .false.
         end if

         tai = 1 + g0i*(alphi - 1.+xii*zi*alphi)
         tbi = 1.-g0i
         tci = g1i*(alphi - 1.+xii*zi*alphi)
         tdi = g2i*(alphi - 1.+xii*zi*alphi)
         tei = g1i
         alphi_twi = xi + 0.5*bi*Tpari*(alphi - 1.)*(1.-g0i)
      end if

      if (kape .ne. 1.) then
         xie = xie*sqrt(2.*(kape + 1.)/(2.*kape - 3.))
         C0pe = (2.*kape - 1.)/(2.*kape - 3.)
         !Branch necessary because the g00 and g00p integrals converge to 1+eps for k rho_i << 1, leading to highly
         !unphysical results
         if (kp < 0.1) then
            g00e = 1.-xe
            g00pe = C0pe - xe
            g10pe = 1.-1.5*xe
            g20pe = 2.*(1.-1.5*xe/C0pe)
         else
            call recursiveintegrate(func00, x0, xcut, kpe, xie, kape, g00e, evals, error, .true.)
            call recursiveintegrate(func00perp, x0, xcut, kpe, xie, kape, g00pe, evals, error, .true.)
            call recursiveintegrate(func10perp, x0, xcut, kpe, xie, kape, g10pe, evals, error, .true.)
            call recursiveintegrate(func20perp, x0, xcut, kpe, xie, kape, g20pe, evals, error, .true.)
         end if

         call recursiveintegrate(refunc00par, 1., xcut + 1., kpe, xie, kape, repart, evals, error, .true.)
         call recursiveintegrate(imfunc00par, 1., xcut + 1., kpe, xie, kape, impart, evals, error, .true.)
         g00be = cmplx(repart, impart)
         call recursiveintegrate(refunc10par, 1., xcut + 1., kpe, xie, kape, repart, evals, error, .true.)
         call recursiveintegrate(imfunc10par, 1., xcut + 1., kpe, xie, kape, impart, evals, error, .true.)
         g10be = cmplx(repart, impart)
         call recursiveintegrate(refunc20par, 1., xcut + 1., kpe, xie, kape, repart, evals, error, .true.)
         call recursiveintegrate(imfunc20par, 1., xcut + 1., kpe, xie, kape, impart, evals, error, .true.)
         g20be = cmplx(repart, impart)

         tae = (C0pe - g00pe + g00be*alphe)*Tperpi/(Tperpe*tr)
         tbe = (C0pe - g00pe)*Tperpi/(Tperpe*tr)
         tce = -(alphe*g10be - g10pe)
         tde = (alphe*g20be - g20pe)*tr*Tperpe/Tperpi
         tee = -g10pe
         alphi_twe = 0.5*bi*Tpare*tr*(alphe - 1.)*(1.-g00e)*mr*Tperpi/(Tperpe*tr)
      else
         !Find Bessel Functions for bimax
         g0e = bessim0(xe)
         g1e = bessim0(xe) - bessim1(xe)
         g2e = 2.*g1e
         ze = 0.
         !Calculate Plasma Dispersion functions
         call zetvec(xie, ze, fail)
         if (fail) then
            !write(*,*) 'Warning! Electron disp function failed!'
            fail = .false.
         end if
         !Calculate intermediate values
         tae = (1.+g0e*(alphe - 1.+xie*ze*alphe))*Tperpi/(Tperpe*tr)
         tbe = (1.-g0e)*Tperpi/(Tperpe*tr)
         tce = -g1e*(alphe - 1.+xie*ze*alphe)
         tde = g2e*(alphe - 1.+xie*ze*alphe)*tr*Tperpe/Tperpi
         tee = -g1e
         alphi_twe = 0.5*bi*Tpare*tr*(alphe - 1.)*(1.-g0e)*mr*Tperpi/(Tperpe*tr)
      end if

      !Calculate intermediate values
      ta = tai + tae
      tb = tbi + tbe
      tc = tci + tce
      td = tdi + tde
      te = tei + tee
      alphi_tw = alphi_twi + alphi_twe

      !Calculate Value of Dispersion relation D(om)
      dispbimkap = (alphi_tw*ta - ta*tb*(om*om) + tb*tb*(om*om))*(2.*ta/(bi*Tperpi) - ta*td + tc*tc) - &
                   (om*om)*(ta*te + tb*tc)**2.

      !DEBUG
      if (.false.) then
         ! write(*,*) 'xi= ',xi,' alphi_twi= ',0.5*bi*Tbi*(alphi-1.)*(1.-g0i),' alphi_twe= ',0.5*bi*Tbe*tr*(alphe-1.)*(1.-g00e)*mr*Tpi/(Tpe*tr)
         ! write(*,*) 'g0i= ',g0i, '; g1i= ',g1i,'; g2i=',g2i
         ! write(*,*) 'g00e= ',g00e, '; g00pe= ',g00pe,'; g10pe=',g10pe,'; g20pe=',g20pe
         ! write(*,*) 'g00be= ',g00be, '; g10be= ',g10be,'; g20be=',g20be
         write (*, *) 'A= ', ta, '; B= ', tb, '; C= ', tc, '; D= ', td, '; E= ', te
         write (*, *) 'alphi_tw= ', alphi_tw
         !write(*,'(5es14.6)')om,dispbikap,abs(dispbikap)
         STOP
      end if

      return
   end function dispbimkap
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the Gyrokinetic Dispersion relation as a function of frequency om
   complex function disp_prime(om)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex :: om                                !Complex Frequency
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments

      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kp*kp*tr/mr

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)

      !Find Bessel Functions
      g0i = bessim0(xi)
      g1i = bessim0(xi) - bessim1(xi)
      g2i = 2.*g1i
      g0e = bessim0(xe)
      g1e = bessim0(xe) - bessim1(xe)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions
      call zetvec(xii, zi, fail)
      if (fail) return
      call zetvec(xie, ze, fail)
      if (fail) return

!     zi=Z_o(xii)
!     ze=Z_o(xie)

      !Calculate intermediate values
      ta = 1 + g0i*xii*zi + (1.+g0e*xie*ze)/tr
      tb = 1.-g0i + (1.-g0e)/tr
      tc = g1i*xii*zi - g1e*xie*ze
      td = g2i*xii*zi + g2e*xie*ze*tr
      te = g1i - g1e

      !Calculate Value of Dispersion relation D'(om)
      disp_prime = om*om/(xi*ta)*(tb*(ta - tb) + &
                                  (ta*te + tb*tc)**2./(2.*ta/bi - ta*td + tc*tc)) - 1.

      return
   end function disp_prime
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the Gyrokinetic Dispersion relation as a function of frequency om
   complex function krook_disp(om)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex :: om                                !Frequency om/(k_par v_A)
      real    :: ai, ae                             !ai=kp^2/2,ae=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !om/(k_par vt_s)
      complex :: zetai, zetae                       !i nu_ss /(k_par vt_s)
      complex :: psii, psie                         !psi_s = xi_s+zeta_s
      complex :: qa, qb, qc, qd, qe, qf                 !Intermediate calculations
      complex :: qg1i, qg2i, qg3i, qg4i               !Intermediate calculations
      complex :: qh1i, qh2i, qh3i                    !Intermediate calculations
      complex :: qg1e, qg2e, qg3e, qg4e               !Intermediate calculations
      complex :: qh1e, qh2e, qh3e                    !Intermediate calculations
      complex :: cta, ctb, ctc, ctd, cte, ctf           !Intermediate calculations
      complex :: gai, gbi, gci                       !Intermediate calculations
      complex :: gae, gbe, gce                       !Intermediate calculations
      complex :: hai, hbi, hci                       !Intermediate calculations
      complex :: hae, hbe, hce                       !Intermediate calculations

      !Set up complex arguments
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)
      zetai = ii*nuii/sqrt(bi)
      zetae = ii*nuee/sqrt(bi*mr*tr)
      psii = xii + zetai
      psie = xie + zetae

      !Calculate arguments of Bessel function
      ai = 0.5*kp*kp
      ae = 0.5*kp*kp*tr/mr

      !Find Bessel Functions
      g0i = bessim0(ai)
      g1i = bessim0(ai) - bessim1(ai)
      g2i = 2.*g1i
      g0e = bessim0(ae)
      g1e = bessim0(ae) - bessim1(ae)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions (functions of psii & psie)
      call zetvec(psii, zi, fail)
      if (fail) return
      call zetvec(psie, ze, fail)
      if (fail) return

      !Calculate Intermediate values:
      qg1i = (g0i*g0i - g1i*g1i*ai)*zetai*zi
      qg1e = (g0e*g0e - g1e*g1e*ae)*zetae*ze/tr
      qg2i = (g0i*g1i - g1i*g2i*ai)*zetai*zi
      qg2e = -(g0e*g1e - g1e*g2e*ae)*zetae*ze
      qg3i = (g1i*g1i - g2i*g2i*ai)*zetai*zi
      qg3e = (g1e*g1e - g2e*g2e*ae)*zetae*ze*tr
      qg4i = (g0i*g0i - g1i*g1i*ai)
      qg4e = (g0e*g0e - g1e*g1e*ae)/tr

      qh1i = 2.*g0i*g0i*(1.+psii*zi)**2.
      qh1e = 2.*g0e*g0e*(1.+psie*ze)**2./tr
      qh2i = 2.*g0i*g1i*(1.+psii*zi)**2.
      qh2e = -2.*g0e*g1e*(1.+psie*ze)**2.
      qh3i = 2.*g1i*g1i*(1.+psii*zi)**2.
      qh3e = 2.*g1e*g1e*(1.+psie*ze)**2.*tr

      gai = xii*zi
      gae = xie*ze
      gbi = 1.+zetai*zi
      gbe = 1.+zetae*ze
      gci = ii*nuii/om*(1.+psii*zi)
      gce = ii*nuee/om*(1.+psie*ze)

      hai = zetai*xii
      hae = zetae*xie
      hbi = zetai*zetai
      hbe = zetae*zetae
      hci = ii*nuii/om*zetai*psii
      hce = ii*nuee/om*zetae*psie

      cta = 1 + g0i*xii*zi + (1.+g0e*xie*ze)/tr
      ctb = 1.-g0i*(1.+zetai*zi) + (1.-g0e*(1.+zetae*ze))/tr
      ctc = g1i*xii*zi - g1e*xie*ze
      ctd = g2i*xii*zi + g2e*xie*ze*tr
      cte = g1i*(1.+zetai*zi) - g1e*(1.+zetae*ze)
      ctf = g0i*ii*nuii/om*(1.+psii*zi) + g0e*ii*nuee/om*(1.+psie*ze)/tr

      qa = cta - (qg1i*gai + qg1e*gae + qh1i*hai + qh1e*hae)
      qb = ctb + (qg1i*gbi + qg1e*gbe + qh1i*hbi + qh1e*hbe)
      qc = ctc - (qg2i*gai + qg2e*gae + qh2i*hai + qh2e*hae)
      qd = ctd - 2./bi - (qg3i*gai + qg3e*gae + qh3i*hai + qh3e*hae)
      qe = cte - (qg2i*gbi + qg2e*gbe + qh2i*hbi + qh2e*hbe)
      qf = ctf - (qg1i*gci + qg1e*gce + qg4i*gci + qg4e*gce + qh1i*hci + qh1e*hce)

      !Calculate Value of Dispersion relation D(om)
!     disp=(1./(om*om) - tb/xi*(1.-tb/ta) )*(1.-0.5*bi*td+0.5*bi*tc*tc/ta) -&
!          0.5*bi/xi*(te+tb*tc/ta)**2.
      !Rewritten for better numerical behavior
!     disp=(xi*ta/(om*om) - ta*tb + tb*tb)*(2.*ta/bi - ta*td + tc*tc - &
!          (ta*te + tb*tc)**2.
!     disp=(ai*ta - ta*tb*(om*om) + tb*tb*(om*om))*(2.*ta/bi - ta*td + tc*tc) - &
!          (om*om)*(ta*te + tb*tc)**2.
      krook_disp = (ai*qa - qa*qf*om*om - qa*qb*om*om + qb*qb*om*om)*(-qa*qd + qc*qc) - &
                   om*om*(qa*qe + qb*qc)**2.

      return
   end function krook_disp
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the Collisional Gyrokinetic Dispersion relation as a function of
!      frequency om using a Krook Collision operator
   complex function col_disp(om)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex :: om                                !Complex Frequency
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments
      complex :: psi, pse                           !Plasma Dispersion arguments

      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kp*kp*tr/mr

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)
      psi = (om + ii*nui)/sqrt(bi)
      pse = (om + ii*nue)/sqrt(bi*mr*tr)

      !Find Bessel Functions
      g0i = bessim0(xi)
      g1i = bessim0(xi) - bessim1(xi)
      g2i = 2.*g1i
      g0e = bessim0(xe)
      g1e = bessim0(xe) - bessim1(xe)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions
      call zetvec(psi, zi, fail)
      if (fail) return
      call zetvec(pse, ze, fail)
      if (fail) return
!     zi=Z_o(psi)
!     ze=Z_o(pse)

      !Calculate intermediate values
      tap = 1 + g0i*psi*zi + (1.+g0e*pse*ze)/tr
      tb = 1.-g0i + (1.-g0e)/tr
      tcp = g1i*xii*zi - g1e*xie*ze
      tcpp = g1i*psi*zi - g1e*pse*ze
      tdp = g2i*xii*zi + g2e*xie*ze*tr
      te = g1i - g1e
      tf = g0i*(1.+psi*zi)*psi + g0e*(1.+pse*ze)*pse*sqrt(mr/tr)

      !Calculate Value of Dispersion relation D(om)
!     disp=(1./(om*om) - tb/xi*(1.-tb/ta) )*(1.-0.5*bi*td+0.5*bi*tc*tc/ta) -&
!          0.5*bi/xi*(te+tb*tc/ta)**2.
      !Rewritten for better numerical behavior
!     disp=(xi*ta/(om*om) - ta*tb + tb*tb)*(2.*ta/bi - ta*td + tc*tc) - &
!          (ta*te + tb*tc)**2.
      col_disp = (xi*tap - om*sqrt(bi)*tb*tf)*(2.*tap/bi - tap*tdp + tcp*tcpp) - &
                 om*(tap*te + tb*tcpp)*(om*tap*(tcpp + te) - sqrt(bi)*tcp*tf)

!     disp=(xi*ta/(om*om) - ta*tb + tb*tb)*(2./bi-td)

      !DEBUG
      if (0 .eq. 1) then
         write (*, '(5es14.6)') om, col_disp, abs(col_disp)
      end if

      return
   end function col_disp
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the Hyperviscous Gyrokinetic Dispersion relation as a function of frequency om
   complex function hv_disp(om)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex :: om                                !Complex Frequency
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments
      complex :: psi, pse                           !Plasma Dispersion arguments
      complex :: psibar                            !Plasma Dispersion arguments
      real    :: nu                                !Hyperviscous coefficient (func of k_perp)

      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kp*kp*tr/mr

      !Hyperviscosity
      nu = nuh*kp**nuh_exp

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)
      psi = (om + ii*nu)/sqrt(bi)    !NOTE: This might be dodgy!
      pse = (om + ii*nu)/sqrt(bi*mr*tr)
      psibar = om + ii*nu

      !Find Bessel Functions
      g0i = bessim0(xi)
      g1i = bessim0(xi) - bessim1(xi)
      g2i = 2.*g1i
      g0e = bessim0(xe)
      g1e = bessim0(xe) - bessim1(xe)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions
      call zetvec(psi, zi, fail)
      if (fail) return
      call zetvec(pse, ze, fail)
      if (fail) return
!     zi=Z_o(psi)
!     ze=Z_o(pse)

      !Calculate intermediate values
      ta = 1 + g0i*xii*zi + (1.+g0e*xie*ze)/tr
      tb = 1.-g0i + (1.-g0e)/tr
      tc = g1i*xii*zi - g1e*xie*ze
      td = g2i*xii*zi + g2e*xie*ze*tr
      te = g1i - g1e

      taa = 1 + 1./tr
      tbb = g0i*xii*zi + (g0e*xie*ze)/tr

      !Calculate Value of Hyperviscous Dispersion relation D(om)
      hv_disp = (xi*ta - ta*tb*(om*om) + tb*tb*(om*om) &
                 + ii*nu*om*((tb - ta)*taa - tb*tbb) + nu*nu*taa*tbb) &
                *(2.*ta/bi - ta*td + tc*tc) &
                - (ta*te*om + tb*tc*om - ii*nu*tc*taa)**2.

!     hv_disp=(xi*(ta+tb) +(tz*om + tb*(om+ii*nu))*(tz*om - ta*(om+ii*nu))) * &
!          ((ta+tb)*(td-2./bi)-tc*tc) + &
!          (c*(tz*om-ta*(om+ii*nu/om)) - te*om*(ta+tb))**2.

!     ta=1+ 1./tr
!     tap=1+g0i*xii*zi+(1.+g0e*xie*ze)/tr
!     tb=g0i*xii*zi + (g0e*xie*ze)/tr
!     tz=g0i + g0e/tr
!     tc=g1i*xii*zi - g1e*xie*ze
!     td=g2i*xii*zi + g2e*xie*ze*tr
!     te=g1i-g1e

!     hv_disp=(xi*tap + (tz*om + tb*psibar)*(tz*om - ta*psibar)) * &
!          (2.*tap/bi - tap*td + tc*tc) - &
!          (tc*(tz*om-ta*psibar) -te*om*tap)**2.

      !DEBUG
      if (0 .eq. 1) then
         write (*, '(5es14.6)') om, hv_disp, abs(hv_disp)
      end if

      return
   end function hv_disp
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the E_perp and B_perp as a function of frequency om
   subroutine epbp(om, eprho, bprho)
      implicit none
      complex, intent(in) :: om                    !Complex Frequency
      real, intent(out) :: eprho, bprho             !|E_perp|rho_i, |B_perp|rho_i
      complex :: x, y, z                             !Solutions

      !Set y and solve for x and z
      y = (1., 0)

      x = (tc*te - tb*(2./bi - td))/(tc*tc + ta*(2./bi - td))
      z = -1.*(ta*te + tb*tc)/(tc*tc + ta*(2./bi - td))

      !NOTE: Changed abs() to real() below
!     eprho=kp*abs(x+y)
!     bprho=kp/abs(om)*sqrt(bi)*c/vti*abs(y)
!     bprho=abs(kp/om*sqrt(bi)*c/vti*y)

      eprho = kp*real(x + y)
      bprho = real(kp/om*sqrt(bi)*c/vti*y)

   end subroutine epbp

!------------------------------------------------------------------------------
!                           Jason TenBarge, 2016
!------------------------------------------------------------------------------
!  Calculates eigenfunctions
   subroutine get_eigenfuncs_bimkap(om, bx, by, bz, ex, ey, ez, dni, dne, vxi, vyi, vzi, vxe, vye, vze)
      use gk_funcs
      use DEIntegrator, only: recursiveintegrate

      implicit none
      complex, intent(in) :: om                    !Complex Frequency
      complex, intent(out) :: bx, by, bz             !B
      complex, intent(out) :: ex, ey, ez             !E
      complex, intent(out) :: dni, dne                !density
      complex, intent(out) :: vxi, vyi, vzi                !ion velocity
      complex, intent(out) :: vxe, vye, vze                !electron velocity
      complex :: x, y, z                             !Solutions
      real :: Tpi, Tpe                                  !(Tperp/T)_s
      real :: Tib, Teb                                  !(Tpar/T)_s

      Tpi = 3.*alphi/(2.*alphi + 1.)
      Tpe = 3.*alphe/(2.*alphe + 1.)
      Tib = 3./(2.*alphi + 1.)
      Teb = 3./(2.*alphe + 1.)

      !Set y and solve for x and z
      !omega Apar / kpar c
      y = (1., 0)
      !phi - y
      x = (tc*te - tb*(2./(bi*Tpi) - td))/(tc*tc + ta*(2./(bi*Tpi) - td))
      !Tperp^i Bpar / qi B0
      z = -1.*(ta*te + tb*tc)/(tc*tc + ta*(2./(bi*Tpi) - td))
      !hat(delta B) = Ti delta B / qi B0
      bx = -krat*z/Tpi
      by = -ii*kp*sqrt(bi/Tpi)*y/(2.*om)
      bz = z/Tpi
      !hat(delta E) = E / kpar
      ex = -ii*(x + y)/krat
      ey = 2.*om*z/(sqrt(bi*Tpi)*kp)
      ez = -ii*x

      !Calculate fluctuations
      !hat(delta ns) = Ti delta ns / qi ns
      dni = -(tai*x + tbi*y + tci*z)/Tpi
      dne = (tae*x + tbe*y + tce*z)/Tpi
      !hat(delta v) = Ti delta vs / qi vA
      vxi = ii*kp*sqrt(bi/Tpi)/2.*(tei*y - tci*x - tdi*z)
      vxe = -ii*kp*sqrt(bi/Tpi)/2.*(tee*y - tce*x - tde*z)
      vyi = 0.
      vye = 0.
      vzi = -om*((tai - tbi)*x + alphi_twi*y/om**2.+(tci + tei)*z)/Tpi
      vze = om*((tae - tbe)*x + alphi_twe*y/om**2.+(tce + tee)*z)/Tpi

   end subroutine get_eigenfuncs_bimkap
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the ion to electron heating ratio P_i/P_e
   subroutine get_eigenfuncs(om, epar, apar, bpar, epar2, apar2, bpar2)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex, intent(in) :: om                    !Complex Frequency
      real, intent(out) ::  epar, apar, bpar         !E_par, A_par, B_par
      real, intent(out) ::  epar2, apar2, bpar2      !E_par, A_par, B_par
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments
      complex :: x, y, z                             !Solutions
      complex :: x2, y2, z2                             !Solutions

      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kp*kp*tr/mr

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)

      !Find Bessel Functions
      g0i = bessim0(xi)
      g1i = bessim0(xi) - bessim1(xi)
      g2i = 2.*g1i
      g0e = bessim0(xe)
      g1e = bessim0(xe) - bessim1(xe)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions
      call zetvec(xii, zi, fail)
      if (fail) return
      call zetvec(xie, ze, fail)
      if (fail) return
!     zi=Z_o(xii)
!     ze=Z_o(xie)

      !Calculate intermediate values
      ta = 1 + g0i*xii*zi + (1.+g0e*xie*ze)/tr
      tb = 1.-g0i + (1.-g0e)/tr
      tc = g1i*xii*zi - g1e*xie*ze
      td = g2i*xii*zi + g2e*xie*ze*tr
      te = g1i - g1e

      !Set y and solve for x and z
      y = (1., 0)
      x = (tc*(xi*ta/om**2.-ta*tb + tb*tb) - tb*(ta*te + tc*tb))/(ta*(ta*te + tb*tc))*y
      z = (xi*ta/om**2.-ta*tb + tb*tb)/(-1.*(ta*te + tb*tc))*y

      y2 = (1., 0)
      x2 = (tc*te - tb*(2./bi - td))/(tc*tc + ta*(2./bi - td))*y2
      z2 = -1.*(ta*te + tb*tc)/(tc*tc + ta*(2./bi - td))*y2

      epar = sqrt(x*conjg(x))
      apar = sqrt(y*conjg(y))
      bpar = sqrt(z*conjg(z))

      epar2 = sqrt(x2*conjg(x2))
      apar2 = sqrt(y2*conjg(y2))
      bpar2 = sqrt(z2*conjg(z2))

   end subroutine get_eigenfuncs
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates the ion to electron heating ratio P_i/P_e
   subroutine heating(om, piti, pete, pipe)
      use gk_bessel, only: bessim0, bessim1
      use gk_funcs, only: Z_o, zetvec
      implicit none
      complex, intent(in) :: om                    !Complex Frequency
      real, intent(out) :: piti, pete, pipe          !P_i*T_i, P_e*T_e,P_i/P_e
      real    :: xi, xe                             !xi=kp^2/2,xe=kp^2/2*tr
      real    :: g0i, g1i, g2i                       !Velocity integrals (ion)
      real    :: g0e, g1e, g2e                       !Velocity integrals (elec)
      complex :: zi, ze                             !Ion/electron Plasma Dispersion
      complex :: xii, xie                           !Plasma Dispersion arguments
      complex :: x, y, z                             !Solutions
      complex :: xc, yc, zc                          !Solutions

      !Calculate arguments of Bessel function
      xi = 0.5*kp*kp
      xe = 0.5*kp*kp*tr/mr

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)

      !Find Bessel Functions
      g0i = bessim0(xi)
      g1i = bessim0(xi) - bessim1(xi)
      g2i = 2.*g1i
      g0e = bessim0(xe)
      g1e = bessim0(xe) - bessim1(xe)
      g2e = 2.*g1e

      !Calculate Plasma Dispersion functions
      call zetvec(xii, zi, fail)
      if (fail) return
      call zetvec(xie, ze, fail)
      if (fail) return
!     zi=Z_o(xii)
!     ze=Z_o(xie)

      !Calculate intermediate values
      ta = 1 + g0i*xii*zi + (1.+g0e*xie*ze)/tr
      tb = 1.-g0i + (1.-g0e)/tr
      tc = g1i*xii*zi - g1e*xie*ze
      td = g2i*xii*zi + g2e*xie*ze*tr
      te = g1i - g1e

      !Set y and solve for x and z
      y = (1., 0)
      x = (tc*te - tb*(2./bi - td))/(tc*tc + ta*(2./bi - td))
      z = -1.*(ta*te + tb*tc)/(tc*tc + ta*(2./bi - td))
!     x=(tc*te - tb*(td-2./bi))/(tc*tc+ta*(td-2./bi))
!     z=-1.*(ta*te+tb*tc)/(tc*tc+ta*(td-2./bi))

      xc = conjg(x) + conjg(y)*(1.-om/conjg(om))
      yc = conjg(y)*om/conjg(om)
!     xc=conjg(x)
!     yc= conjg(y)
      zc = conjg(z)

      piti = sqrt(pi)*om*conjg(om)*exp(-om*conjg(om)/bi)* &
             (g0i*x*conjg(x) + g1i*(x*conjg(z) + z*conjg(x)) + g2i*z*conjg(z))
      pete = sqrt(pi)*om*conjg(om)*exp(-om*conjg(om)/(bi*tr*mr))* &
             (g0e*x*conjg(x)/tr - g1e*(x*conjg(z) + z*conjg(x)) + tr*g2e*z*conjg(z))

!     piti=real(ii*om*((xc*y+x*yc+y*yc)*g0i + (y*zc + yc*z)*g1i - &
!          xii*zi*(x*xc*g0i + (xc*z+x*zc)*g1i + z*zc*g2i)))
!     pete=real(ii*om*((xc*y+x*yc+y*yc)*g0e - tr*(y*zc + yc*z)*g1e - &
!          xie*ze*(x*xc*g0e - tr*(xc*z+x*zc)*g1e + tr*tr*z*zc*g2e)))

!     piti=real(-ii*om*((x*y*2.*g0i + y*y*g0i + y*z*g1i) - &
!          xii*zi*(x*x*g0i + x*z*2.*g1i + z*z*g2i)))
!     pete=real(-ii*om*((x*y*2.*g0e + y*y*g0e -tr*y*z*g1e) - &
!          xie*ze*(x*x*g0e - tr*x*z*2.*g1e + tr*tr*z*z*g2e)))

!     piti=abs(((x*y*2.*g0i + y*y*g0i + y*z*g1i) - &
!          xii*zi*(x*x*g0i + x*z*2.*g1i + z*z*g2i)))
!     pete=abs(((x*y*2.*g0e + y*y*g0e -tr*y*z*g1e) - &
!          xie*ze*(x*x*g0e - tr*x*z*2.*g1e + tr*tr*z*z*g2e)))

      pipe = piti/pete*sqrt(tr*mr)

!     pipe=tr*real(ii*((xc*y+x*yc+y*yc)*g0i + (y*zc + yc*z)*g1i - &
!          xii*zi*(x*xc*g0i + (xc*z+x*zc)*g1i + z*zc*g2i)))/ &
!          real(ii*((xc*y+x*yc+y*yc)*g0e - tr*(y*zc + yc*z)*g1e - &
!          xie*ze*(x*xc*g0e - tr*(xc*z+x*zc)*g1e + tr*tr*z*zc*g2e)))

!     pipe=(-ii*om*((x*y*2.*g0i + y*y*g0i + y*z*g1i) - &
!          xii*zi*(x*x*g0i + x*z*2.*g1i + z*z*g2i))) / &
!          (-ii*om*((x*y*2.*g0e + y*y*g0e -tr*y*z*g1e) - &
!          xie*ze*(x*x*g0e - tr*x*z*2.*g1e + tr*tr*z*z*g2e)))

   end subroutine heating
!------------------------------------------------------------------------------
!                           Greg Howes, 2005
!------------------------------------------------------------------------------
!  Calculates distribution functions g(s,vperp,vpar)
   subroutine dist(om, vperp, vpar, gdist, hdist)
      use gk_bessel, only: bessj0, bessj1
      implicit none
      complex, intent(in) :: om                    !Complex Frequency
      complex, dimension(1:, 1:, -nl:), intent(out) :: gdist, hdist  !Distribution function g_s
      real, dimension(1:, -nl:), intent(out) :: vperp, vpar !V-space grid values
      complex :: x, y, z                             !Solutions
      complex :: xii, xie                           !Plasma Dispersion arguments
      real :: de, dl                              !Energy, pitch angle intervals
      real :: e, l                                !Energy, pitch angle
      real :: gri, gre                            !Bessel Function arguments
      real :: j0i, j1i, j0e, j1e                    !Bessel Functions
      integer :: i, j

      !Set up vspace grid values
      de = ecut/real(ne)
      dl = 2/real(2*nl)
      !Determine vperp and vpar for polar (e,l) grid (l=vpar/v)
      do i = 1, ne
         e = de*real(i)
         do j = -nl, nl
            l = dl*real(j)
            vperp(i, j) = sqrt(e*(1 - l*l))
            vpar(i, j) = l*sqrt(e)
         end do
      end do

      !Set y and solve for x and z
      y = (1., 0)
      x = (tc*te - tb*(2./bi - td))/(tc*tc + ta*(2./bi - td))
      z = -1.*(ta*te + tb*tc)/(tc*tc + ta*(2./bi - td))

      !NOTE: REMOVE THIS!
!     z=(0.0,0.0)

      !Calculate arguments of Bessel function J_n
      gri = kp
      gre = kp*sqrt(tr/mr)

      !Calculate arguments of plasma dispersion function
      xii = om/sqrt(bi)
      xie = om/sqrt(bi*mr*tr)

      !Calculate distribution functions for each v-space point
      !NOTE: ions=1, electrons=2
      gdist = 0.
      hdist = 0.
      do i = 1, ne
         do j = -nl, nl

            !Find Bessel Functions
            j0i = bessj0(gri*vperp(i, j))
            j1i = bessj1(gri*vperp(i, j))
            j0e = bessj0(gre*vperp(i, j))
            j1e = bessj1(gre*vperp(i, j))

            !Ions
            gdist(1, i, j) = vpar(i, j)/(xii - cmplx(vpar(i, j), 0.))*(j0i*x + &
                                                                       j1i/kp*2.*vperp(i, j)*z)
            !Electrons
            gdist(2, i, j) = -1./tr*vpar(i, j)/(xie - cmplx(vpar(i, j), 0.))*(j0e*x + &
                                                                              j1e/kp*sqrt(tr*mr)*2.*vperp(i, j)*z)

            !Ions
            hdist(1, i, j) = j0i*y + (xii/(xii - cmplx(vpar(i, j), 0.)))*(j0i*x + &
                                                                          j1i/kp*2.*vperp(i, j)*z)
            !Electrons
            hdist(2, i, j) = -1./tr*(j0e*y + xie/(xie - cmplx(vpar(i, j), 0.))*(j0e*x + &
                                                                                j1e/kp*sqrt(tr*mr)*2.*vperp(i, j)*z))

         end do
      end do

   end subroutine dist
!------------------------------------------------------------------------------
!                           Greg Howes, 2007
!------------------------------------------------------------------------------
!  Calculate the parallel wavenumber as a function of kperp for the critically
!  balanced MHD Alfven and kinetic Alfven wave cascades
   real function kpara(tk0, tkp, tbi, ttite)
      implicit none
      !Passed
      real :: tk0                                   !Driving scale k_0 rho_i
      real :: tkp                                   !k_perp rho_i
      real :: tbi                                   !beta_i
      real :: ttite                                 !Ti/Te
      !Local
      real, parameter :: kcm1 = 2.5                   !Kolmogorov Constant ~6
      real, parameter :: kcm2 = 2.2                   !Kolmogorov Constant ~6
      real, parameter :: kck1 = 2.5                   !Kolmogorov Constant ~6
      real, parameter :: kck2 = 2.2                   !Kolmogorov Constant ~6
      real :: d                                    !sqrt(bi+2/(1+1/tite))

      d = sqrt(tbi + 2./(1.+1./ttite))

      kpara = tk0**(1./3.)*(kcm2/kcm1**0.5*tkp**(2./3.) + &
                            kck2*d/(kck1**0.5*tbi**0.5)*tkp**(7./3.))/(1.+tkp**2.)

      return
   end function kpara
!------------------------------------------------------------------------------
end module gk_disp_gf

