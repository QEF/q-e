!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dgcxc_spin (rup, rdw, gup, gdw, vrrxup, vrrxdw, vrsxup, &
 vrsxdw, vssxup, vssxdw, vrrcup, vrrcdw, vrscup, vrscdw, vssc, &
 vrzcup, vrzcdw)
!-----------------------------------------------------------------------
!
!    This routine computes the derivative of the exchange and correlatio
!    potentials with respect to the density, the gradient and zeta
!
use parameters, only : DP
implicit none  
real(kind=DP) :: rup, rdw, gup (3), gdw (3), vrrxup, vrrxdw, vrsxup, &
 vrsxdw, vssxup, vssxdw, vrrcup, vrrcdw, vrscup, vrscdw, vssc, &
 vrzcup, vrzcdw
                                   ! input: the charges and the gradient
                                   ! output: derivatives of the exchange
                                   ! output: derivatives of the exchange
                                   ! output: derivatives of the correlat
                                   ! output: derivatives of the correlat
                                   ! output: derivatives of the correlat
!
!    local variables
!
real(kind=DP) :: r, zeta, sup2, sdw2, s2, s, sup, sdw, dr, dzeta, ds, &
 drup, drdw, dsup, dsdw, sx, sc, v1xupp, v1xdwp, v2xupp, v2xdwp, &
 v1xupm, v1xdwm, v2xupm, v2xdwm, v1cupp, v1cdwp, v2cp, v1cupm, &
 v1cdwm, v2cm
                                        ! charge densities and square gr
                                        ! gradients
                                        ! delta charge densities and gra
                                        ! delta gradients
                                        ! energies
                                          ! exchange potentials
                                          ! exchange potentials
                                          ! coorelation potentials
                                          ! coorelation potentials
real(kind=DP) :: eps  

parameter (eps = 1.d-6)  
r = rup + rdw  
if (r.gt.eps) then  
   zeta = (rup - rdw) / r  
else  
   zeta = 2.d0  

endif  
sup2 = gup (1) **2 + gup (2) **2 + gup (3) **2  
sdw2 = gdw (1) **2 + gdw (2) **2 + gdw (3) **2  

s2 = (gup (1) + gdw (1) ) **2 + (gup (2) + gdw (2) ) **2 + &
 (gup (3) + gdw (3) ) **2
sup = sqrt (sup2)  
sdw = sqrt (sdw2)  
s = sqrt (s2)  
!
!     up part of exchange
!

if (rup.gt.eps.and.sup.gt.eps) then  
   drup = min (1.d-4, 1.d-2 * rup)  
   dsup = min (1.d-4, 1.d-2 * sdw)  
!
!    derivatives of exchange: up part
!
   call gcx_spin (rup + drup, rdw, sup2, sdw2, sx, v1xupp, v1xdwp, &
    v2xupp, v2xdwp)

   call gcx_spin (rup - drup, rdw, sup2, sdw2, sx, v1xupm, v1xdwm, &
    v2xupm, v2xdwm)
   vrrxup = 0.5d0 * (v1xupp - v1xupm) / drup  


   vrsxup = 0.25d0 * (v2xupp - v2xupm) / drup  
   call gcx_spin (rup, rdw, (sup + dsup) **2, sdw2, sx, v1xupp, &
    v1xdwp, v2xupp, v2xdwp)

   call gcx_spin (rup, rdw, (sup - dsup) **2, sdw2, sx, v1xupm, &
    v1xdwm, v2xupm, v2xdwm)
   vrsxup = vrsxup + 0.25d0 * (v1xupp - v1xupm) / dsup / sup  
   vssxup = 0.5d0 * (v2xupp - v2xupm) / dsup / sup  
else  
   vrrxup = 0.d0  
   vrsxup = 0.d0  
   vssxup = 0.d0  

endif  

if (rdw.gt.eps.and.sdw.gt.eps) then  
   drdw = min (1.d-4, 1.d-2 * rdw)  
   dsdw = min (1.d-4, 1.d-2 * sdw)  
!
!    derivatives of exchange: down part
!
   call gcx_spin (rup, rdw + drdw, sup2, sdw2, sx, v1xupp, v1xdwp, &
    v2xupp, v2xdwp)

   call gcx_spin (rup, rdw - drdw, sup2, sdw2, sx, v1xupm, v1xdwm, &
    v2xupm, v2xdwm)
   vrrxdw = 0.5d0 * (v1xdwp - v1xdwm) / drdw  

   vrsxdw = 0.25d0 * (v2xdwp - v2xdwm) / drdw  
   call gcx_spin (rup, rdw, sup2, (sdw + dsdw) **2, sx, v1xupp, &
    v1xdwp, v2xupp, v2xdwp)

   call gcx_spin (rup, rdw, sup2, (sdw - dsdw) **2, sx, v1xupm, &
    v1xdwm, v2xupm, v2xdwm)
   vrsxdw = vrsxdw + 0.25d0 * (v1xdwp - v1xdwm) / dsdw / sdw  
   vssxdw = 0.5d0 * (v2xdwp - v2xdwm) / dsdw / sdw  
else  
   vrrxdw = 0.d0  
   vrsxdw = 0.d0  
   vssxdw = 0.d0  
endif  
!
!     derivatives of correlation
!

if (r.gt.eps.and.abs (zeta) .le.1.d0.and.s.gt.eps) then  

   dr = min (1.d-4, 1.d-2 * r)  
   call gcc_spin (r + dr, zeta, s2, sc, v1cupp, v1cdwp, v2cp)  

   call gcc_spin (r - dr, zeta, s2, sc, v1cupm, v1cdwm, v2cm)  
   vrrcup = 0.5d0 * (v1cupp - v1cupm) / dr  

   vrrcdw = 0.5d0 * (v1cdwp - v1cdwm) / dr  

   ds = min (1.d-4, 1.d-2 * s)  
   call gcc_spin (r, zeta, (s + ds) **2, sc, v1cupp, v1cdwp, v2cp)  

   call gcc_spin (r, zeta, (s - ds) **2, sc, v1cupm, v1cdwm, v2cm)  
   vrscup = 0.5d0 * (v1cupp - v1cupm) / ds / s  
   vrscdw = 0.5d0 * (v1cdwp - v1cdwm) / ds / s  

   vssc = 0.5d0 * (v2cp - v2cm) / ds / s  
   dzeta = min (1.d-4, 1.d-2 * abs (zeta) )  

   if (dzeta.lt.1.d-7) dzeta = 1.d-7  
   call gcc_spin (r, zeta + dzeta, s2, sc, v1cupp, v1cdwp, v2cp)  

   call gcc_spin (r, zeta - dzeta, s2, sc, v1cupm, v1cdwm, v2cm)  
   vrzcup = 0.5d0 * (v1cupp - v1cupm) / dzeta  
   vrzcdw = 0.5d0 * (v1cdwp - v1cdwm) / dzeta  
else  
   vrrcup = 0.d0  
   vrrcdw = 0.d0  
   vrscup = 0.d0  
   vrscdw = 0.d0  
   vssc = 0.d0  
   vrzcup = 0.d0  
   vrzcdw = 0.d0  

endif  
return  
end subroutine dgcxc_spin
