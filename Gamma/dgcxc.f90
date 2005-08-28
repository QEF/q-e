!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dgcxc (r, s2, vrrx, vsrx, vssx, vrrc, vsrc, vssc)
!-----------------------------------------------------------------------
USE kinds, only : DP
implicit none
real(DP) :: r, s2, vrrx, vsrx, vssx, vrrc, vsrc, vssc
real(DP) :: dr, s, ds

real(DP) :: sx, sc, v1xp, v2xp, v1cp, v2cp, v1xm, v2xm, v1cm, &
 v2cm
s = sqrt (s2)
dr = min (1.d-4, 1.d-2 * r)

ds = min (1.d-4, 1.d-2 * s)
call gcxc (r + dr, s2, sx, sc, v1xp, v2xp, v1cp, v2cp)

call gcxc (r - dr, s2, sx, sc, v1xm, v2xm, v1cm, v2cm)
vrrx = 0.5d0 * (v1xp - v1xm) / dr

vrrc = 0.5d0 * (v1cp - v1cm) / dr
vsrx = 0.25d0 * (v2xp - v2xm) / dr

vsrc = 0.25d0 * (v2cp - v2cm) / dr
call gcxc (r, (s + ds) **2, sx, sc, v1xp, v2xp, v1cp, v2cp)

call gcxc (r, (s - ds) **2, sx, sc, v1xm, v2xm, v1cm, v2cm)
vsrx = vsrx + 0.25d0 * (v1xp - v1xm) / ds / s

vsrc = vsrc + 0.25d0 * (v1cp - v1cm) / ds / s
vssx = 0.5d0 * (v2xp - v2xm) / ds / s

vssc = 0.5d0 * (v2cp - v2cm) / ds / s
return
end subroutine dgcxc
