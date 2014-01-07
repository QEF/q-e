!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine drhodvloc (nu_i0, npe, drhoscf, wdyn)
  !-----------------------------------------------------------------------
  !
  !    This subroutine computes the contribution of the local
  !    potential to the electronic term <psi|dv-e ds|dpsi> of the dynamical
  !    matrix. It can be used both for KB and for US pseudopotentials.
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : nat
  USE fft_base, ONLY : dfftp, dffts
  USE cell_base, ONLY : omega
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag

  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  implicit none

  integer :: npe, nu_i0
  ! input: the number of perturbation of this representations
  ! input: the initial position of the mode
  complex(DP) :: drhoscf (dfftp%nnr, nspin_mag, npe), wdyn (3 * nat, 3 * nat)
  ! the change of density due to perturbations
  ! auxiliary matrix where drhodv is stored

  integer :: ipert, is, nu_i, nu_j
  ! counter on perturbations
  ! counter on spin polarizations
  ! counter on the i modes
  ! counter on the j modes

  complex(DP) :: zdotc, dynwrk (3 * nat, 3 * nat)
  complex(DP), allocatable :: dvloc (:)
  ! d Vloc / dtau

  allocate (dvloc( dffts%nnr))
  dynwrk (:,:) = (0.d0, 0.d0)
  !
  ! We need a sum over all perturbations
  !
  do nu_j = 1, 3 * nat
     call compute_dvloc (nu_j, dvloc)
     do ipert = 1, npe
        nu_i = nu_i0 + ipert
        do is = 1, nspin_lsda
           dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + &
                zdotc (dffts%nnr, drhoscf (1, is, ipert), 1, dvloc, 1) * &
                  omega / (dffts%nr1 * dffts%nr2 * dffts%nr3)
        enddo
     enddo

  enddo
  !
  ! collect contributions from nodes of a pool (sum over G & R space)
  !
  call mp_sum ( dynwrk, intra_bgrp_comm )

  wdyn(:,:) = wdyn(:,:) + dynwrk(:,:)
  deallocate(dvloc)
  return
end subroutine drhodvloc
