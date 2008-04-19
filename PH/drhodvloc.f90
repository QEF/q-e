!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine drhodvloc (nu_i0, nper, drhoscf, wdyn)
  !-----------------------------------------------------------------------
  ! following comment is obsolete
  !    This subroutine computes the electronic term
  !    <psi|dv|dpsi> of the dynamical matrix. It can be used both for KB
  !    and for US pseudopotentials.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat
  use pwcom
  USE kinds, only : DP
  use phcom
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  implicit none

  integer :: nper, nu_i0
  ! input: the number of perturbation of this representations
  ! input: the initial position of the mode
  complex(DP) :: drhoscf (nrxx, nspin, npertx), wdyn (3 * nat, 3 * nat)
  ! the change of density due to perturbations
  ! auxiliary matrix where drhodv is stored

  integer :: ipert, is, nu_i, nu_j, nspin0
  ! counter on perturbations
  ! counter on spin polarizations
  ! counter on the i modes
  ! counter on the j modes

  complex(DP) :: ZDOTC, dynwrk (3 * nat, 3 * nat)
  complex(DP), allocatable :: dvloc (:)
  ! d Vloc / dtau

  nspin0=nspin
  if (nspin==4) nspin0=1
  allocate (dvloc( nrxxs))    
  dynwrk (:,:) = (0.d0, 0.d0)
  !
  ! We need a sum over all perturbations
  !
  do nu_j = 1, 3 * nat
     call compute_dvloc (nu_j, dvloc)
     do ipert = 1, nper
        nu_i = nu_i0 + ipert
        do is = 1, nspin0
           dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + &
                ZDOTC (nrxxs, drhoscf (1, is, ipert), 1, dvloc, 1) * &
                  omega / (nr1s * nr2s * nr3s)
        enddo
     enddo

  enddo
#ifdef __PARA
  !
  ! collect contributions from nodes of a pool (sum over G & R space)
  !
  call mp_sum ( dynwrk, intra_pool_comm )
#endif

  wdyn(:,:) = wdyn(:,:) + dynwrk(:,:)
  deallocate(dvloc)
  return
end subroutine drhodvloc
