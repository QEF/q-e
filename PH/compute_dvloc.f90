!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine compute_dvloc (mode, dvlocin)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "machine.h"

  use pwcom
  USE kinds, only : DP
  use phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: mode
  ! input: the actual perturbation

  complex(kind=DP) :: dvlocin (nrxxs)
  ! output: the change of the local potential
  !
  !   And the local variables
  !
  integer :: na, nt, mu, ig, ibnd, ir
  ! counters
  complex(kind=DP) :: gtau, gu, fact, u1, u2, u3, gu0
  ! auxiliary variables

  call start_clock ('com_dvloc')
  dvlocin (:) = (0.d0, 0.d0)
  do na = 1, nat
     fact = tpiba * (0.d0, - 1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if ( abs (u (mu + 1, mode) ) + abs (u (mu + 2, mode) ) + &
          abs (u (mu + 3, mode) ) > 1.0d-12) then
        nt = ityp (na)
        u1 = u (mu + 1, mode)
        u2 = u (mu + 2, mode)
        u3 = u (mu + 3, mode)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) * eigts3 ( &
                ig3 (ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           dvlocin (nls (ig) ) = dvlocin (nls (ig) ) + vlocq (ig, nt) &
                * gu * fact * gtau
        enddo
     endif
  enddo
  !
  ! Now we compute dV_loc/dtau in real space
  !

  call cft3s (dvlocin, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)

  call stop_clock ('com_dvloc')
  return
end subroutine compute_dvloc
