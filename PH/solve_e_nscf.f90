!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e_nscf( avg_iter, thresh, ik, ipol, dvscfs, &
                         auxg, spsi, auxr )
  !-----------------------------------------------------------------------
  !
  !   Solve the linear system which defines the change of the wavefunctions
  !   due to the electric field for a given k_point in a non self-consistent
  !   way. The self-consistent variation of the potential has been computed
  !   previously and is in dvscfs.
  !
#include "f_defs.h"
  use kinds, only : DP
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  use becmod
  use phcom
  implicit none

  !
  !  Input variables
  !
  integer :: ik, ipol
  ! input: k-point under consideration
  ! input: polarization of the electric field

  real(kind=DP) :: thresh, avg_iter
  ! input: convergence threshold
  ! in/out: # of diagonalization iterations

  complex(kind=DP) :: dvscfs (nrxxs, 3), auxg (npwx), spsi (npwx), &
                      auxr(nrxxs)
  ! input: potential on the smooth grid
  ! auxiliary space
  ! auxiliary space

  !
  !  Local variables
  !
  integer :: ibnd, ir, ig, nrec
  ! counter on bands
  ! counter on mesh points
  ! counter on G-points
  ! the record number

  !
  ! Calculates [H,x]*psi_kpoint
  !
  dpsi (:,:) = (0.d0, 0.d0)
  this_pcxpsi_is_on_file(:,:)=.false.
  call dvpsi_e (ik, ipol)

  do ig = 1, npw
     g2kin (ig) = ( (xk (1, ik) + g(1, igk (ig)) ) **2 + &
                    (xk (2, ik) + g(2, igk (ig)) ) **2 + &
                    (xk (3, ik) + g(3, igk (ig)) ) **2 ) *tpiba2
  enddo
  !
  ! Calculates dvscf*psi_k in G_space,
  !
  do ibnd = 1, nbnd_occ (ik)
     auxr (:) = (0.d0, 0.d0)
     do ig = 1, npw
        auxr (nls (igk (ig))) = evc (ig, ibnd)
     end do
     call cft3s (auxr, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
     do ir = 1, nrxxs
        auxr (ir) = auxr(ir) * dvscfs(ir, ipol)
     end do
     call cft3s (auxr, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
     do ig = 1, npwq
        dvpsi (ig, ibnd) = dvpsi(ig, ibnd) + auxr(nls (igkq (ig)))
     enddo
  enddo
  !
  ! starting value for  delta_psi is read from iudwf
  !
  nrec = (ipol - 1) * nksq + ik
  call davcio (dpsi, lrdwf, iudwf, nrec, -1)
  call pcgreen (avg_iter, thresh, ik, et (1, ik), auxg, spsi)

  return
end subroutine solve_e_nscf
