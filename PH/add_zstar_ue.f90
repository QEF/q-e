!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine add_zstar_ue (imode0, npe)
  !-----------------------------------------------------------------------
  ! add the contribution of the modes imode0+1 -> imode+npe
  ! to the effective charges Z(Us,E) (Us=scf,E=bare)
  !
  ! trans =.true. is needed for this calculation to be meaningful
  !
#include "f_defs.h"
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  USE io_files, ONLY: iunigk
  use phcom
  implicit none

  integer :: imode0, npe

  integer :: ibnd, jpol, ipert, nrec, mode, ik
  ! counter on bands
  ! counter on polarization
  ! counter on pertubations
  ! counter on records
  ! counter on modes
  ! counter on k points

  real(DP) :: weight

  complex(DP) :: ZDOTC

  call start_clock('add_zstar_ue')
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) read (iunigk) npw, igk
     npwq = npw
     weight = wk (ik)
     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     do jpol = 1, 3
        !
        ! read/compute DeltaV*psi(bare) for electric field
        !
        call dvpsi_e (ik, jpol)
        !
        do ipert = 1, npe
           mode = imode0 + ipert
           nrec = (ipert - 1) * nksq + ik
           !
           ! read DeltaV*psi(scf) for phonon mode # mode
           !

           call davcio (dpsi, lrdwf, iudwf, nrec, -1)
           do ibnd = 1, nbnd_occ(ik)
              zstarue0 (mode, jpol) = zstarue0 (mode, jpol) - 2.d0 * weight * &
                   ZDOTC (npw, dpsi (1, ibnd), 1, dvpsi (1, ibnd), 1)
           enddo
        enddo
     enddo

  enddo
  call stop_clock('add_zstar_ue')
  return
end subroutine add_zstar_ue
