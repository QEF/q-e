!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine constrain (dist0, g, dg, dg2, dist, nat, tau, alat)
  !
  !     this is a user supplied routine which define the constrain equation
  !
  !     g(tau,dist) = 0
  !
  !     where tau are the atomic positions ( in alat units ) and dist is
  !     the distance of two atoms which is in this case the one
  !     dimensional constrain, dg is in output the value of the gradient
  !     of g and dg2 is its square modulus.
  !     DA 1997
  !
#include "machine.h"
  USE kinds, ONLY: DP
  implicit none
  !
  integer, INTENT(IN) :: nat
  real(kind=DP), INTENT(IN) :: tau (3, nat), alat, dist
  real(kind=DP), INTENT(OUT):: dg (3, nat), dg2, g, dist0
  !
  real(kind=DP) :: x1, x2, y1, y2, z1, z2
  real(kind=DP), external :: DDOT
  !
  x1 = tau (1, 1)
  y1 = tau (2, 1)
  z1 = tau (3, 1)
  x2 = tau (1, 2)
  y2 = tau (2, 2)
  z2 = tau (3, 2)
  !
  dist0 = (x1 - x2) **2 + (y1 - y2) **2 + (z1 - z2) **2
  !
  g = dist0 - dist
  !
  dg (:,:) = 0.d0
  dg (1, 1) = 2 * (x1 - x2)
  dg (1, 2) = 2 * (x2 - x1)
  dg (2, 1) = 2 * (y1 - y2)
  dg (2, 2) = 2 * (y2 - y1)
  dg (3, 1) = 2 * (z1 - z2)
  dg (3, 2) = 2 * (z2 - z1)
  !
  dg2 = DDOT (3 * nat, dg, 1, dg, 1)
  !
  return
end subroutine constrain

