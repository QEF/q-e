!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
! ========================================================
!
!
! Routine to evaluate the derivative of the exchange potential for exchange current
!
subroutine pbex_current(rho, grho, iflag, sx)
   USE kinds, ONLY: DP
   USE constants, ONLY: pi
   use exch_gga, only: pbex
   implicit none
   real(DP) :: rho, grho(1:3) !nb, qui compare tutto il gradiente
   real(DP) :: sx(1:3) !
   integer  :: iflag
   real(DP) :: kf, agrho, s1, s2, exunif, dsg, fx
   ! (3*pi2*|rho|)^(1/3)
   ! |grho|
   ! |grho|/(2*kf*|rho|)
   ! s^2
   ! exchange energy LDA part
   real(DP) :: f1, f2, f3
   ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
   real(DP) :: third, c1, c2
   parameter(third=1.d0/3.d0, c1=0.75d0/pi, &
             c2=3.093667726280136d0)
   ! parameters of the functional
   real(DP) :: k(3), mu(3)
   data k/0.804d0, 1.2450D0, 0.804d0/, &
      mu/0.21951d0, 0.21951d0, 0.12345679012345679012d0/
   !
   agrho = sqrt(grho(1)**2 + grho(2)**2 + grho(3)**2)

   kf = c2*rho**third

   dsg = 0.5d0/kf

   s1 = agrho*dsg/rho
   s2 = s1*s1

   !
   !   current
   !
   f1 = s2*mu(iflag)/k(iflag)
   f2 = 1.d0 + f1
   f3 = 2.d0*mu(iflag)/(f2*f2)
   fx = f3/((2.d0*kf*rho)*(2.d0*kf*rho))
!
   exunif = -c1*kf
!
   sx(1:3) = exunif*fx*grho(1:3)

end subroutine pbex_current
!
!
!=============================================================================
!
! Only for testing, the previous routine using a finite difference approach.
!
subroutine pbex_current_numerical(rho, grho, iflag, sx)
   USE kinds, ONLY: DP
   USE constants, ONLY: pi
   use exch_gga, only: pbex
   implicit none
!
   real(DP) :: rho, grho(1:3) !nb, qui compare tutto il gradiente
   real(DP) :: sx(1:3) !
   integer  :: iflag
   real(DP) :: dgrho, agrho, sx_one, sx_two, der, delta
   real(DP) :: v1x, v2x !servono solo per chiamare la routine
!
   delta = 1.0E-10_DP
!
   agrho = grho(1)**2 + grho(2)**2 + grho(3)**2
   dgrho = (grho(1) + delta)**2 + grho(2)**2 + grho(3)**2

   call pbex(rho, agrho, iflag, sx_one, v1x, v2x)

   call pbex(rho, dgrho, iflag, sx_two, v1x, v2x)

   der = (sx_two - sx_one)/delta

   sx(1) = der/rho

   agrho = grho(1)**2 + grho(2)**2 + grho(3)**2
   dgrho = grho(1)**2 + (grho(2) + delta)**2 + grho(3)**2

   call pbex(rho, agrho, iflag, sx_one, v1x, v2x)

   call pbex(rho, dgrho, iflag, sx_two, v1x, v2x)

   der = (sx_two - sx_one)/delta

   sx(2) = der/rho

   agrho = grho(1)**2 + grho(2)**2 + grho(3)**2
   dgrho = grho(1)**2 + grho(2)**2 + (grho(3) + delta)**2

   call pbex(rho, agrho, iflag, sx_one, v1x, v2x)

   call pbex(rho, dgrho, iflag, sx_two, v1x, v2x)

   der = (sx_two - sx_one)/delta

   sx(3) = der/rho

end subroutine pbex_current_numerical
!
!
!=============== Routine for correlation part of the current  ================================================
!
subroutine pbec_current(rho, grho, iflag, sc)
   !---------------------------------------------------------------
   !
   ! PBE correlation (without LDA part)
   ! iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
   ! iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
   !
   USE kinds, ONLY: DP
   use corr_gga, only: pbec, pw
   implicit none
   integer, intent(in) :: iflag
   real(DP) :: rho, grho(1:3), sc(1:3)
!
   real(DP) :: ga, be(2), ec, vc
   parameter(ga=0.031091d0)
   data be/0.066725d0, 0.046d0/
   real(DP) :: third, pi34, xkf, xks
   parameter(third=1.d0/3.d0, pi34=0.6203504908994d0)
   parameter(xkf=1.919158292677513d0, xks=1.128379167095513d0)
   ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
   real(DP) :: kf, ks, rs, t, expe, af, bf, y, ay, by, cy
   !
   rs = pi34/rho**third
   call pw(rs, 1, ec, vc)
   kf = xkf/rs
   ks = xks*sqrt(kf)
   t = sqrt(grho(1)**2 + grho(2)**2 + grho(3)**2)/(2.d0*ks*rho)

   expe = exp(-ec/ga)

   af = be(iflag)/ga*(1.d0/(expe - 1.d0))

   y = af*t*t
   ay = (1.d0 + 2.d0*y)/(1.d0 + y + y*y)
   by = (y + y*y)*(be(iflag)/af + ga) + ga
   cy = 2.d0*be(iflag)*ga*ay/by

!
   sc(1:3) = grho(1:3)/((2.d0*ks*rho)*(2.d0*ks*rho))*cy

   !
   return
end subroutine pbec_current

!====================================================================================
! Only for testing, previous routine with a finite different approach

subroutine pbec_current_numerical(rho, grho, iflag, sc)
   !---------------------------------------------------------------
   !
   ! PBE correlation (without LDA part)
   ! iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
   ! iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
   !
   USE kinds, ONLY: DP
   USE constants, ONLY: pi
   use corr_gga, only: pbec
   implicit none
!
   real(DP) :: rho, grho(1:3) !nb, qui compare tutto il gradiente
   real(DP) :: sc(1:3) !
   integer  :: iflag
   real(DP) :: dgrho, agrho, sc_one, sc_two, der, delta
   real(DP) :: v1c, v2c !servono solo per chiamare la routine
!
   delta = 1.0E-10_DP
!
   agrho = grho(1)**2 + grho(2)**2 + grho(3)**2
   dgrho = (grho(1) + delta)**2 + grho(2)**2 + grho(3)**2

   call pbec(rho, agrho, iflag, sc_one, v1c, v2c)

   call pbec(rho, dgrho, iflag, sc_two, v1c, v2c)

   der = (sc_two - sc_one)/delta

   sc(1) = der/rho

   agrho = grho(1)**2 + grho(2)**2 + grho(3)**2
   dgrho = grho(1)**2 + (grho(2) + delta)**2 + grho(3)**2

   call pbec(rho, agrho, iflag, sc_one, v1c, v2c)

   call pbec(rho, dgrho, iflag, sc_two, v1c, v2c)

   der = (sc_two - sc_one)/delta

   sc(2) = der/rho

   agrho = grho(1)**2 + grho(2)**2 + grho(3)**2
   dgrho = grho(1)**2 + grho(2)**2 + (grho(3) + delta)**2

   call pbec(rho, agrho, iflag, sc_one, v1c, v2c)

   call pbec(rho, dgrho, iflag, sc_two, v1c, v2c)

   der = (sc_two - sc_one)/delta

   sc(3) = der/rho

end subroutine pbec_current_numerical

!================================================================================
