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
#include "machine.h"

  use pwcom
  use parameters, only : DP
  use phcom
  implicit none

  integer :: nper, nu_i0
  ! input: the number of perturbation of this repres
  ! input: the initial position of the mode
  complex(kind=DP) :: drhoscf (nrxx, nspin, npertx), wdyn (3 * nat, 3 * &
       nat)
  ! the change of density due to
  ! perturbations
  ! auxiliary matrix where drhodv is stored

  integer :: ipert, is, nu_i, nu_j
  ! counter on perturbations
  ! counter on spin polarizations
  ! counter on the i modes
  ! counter on the j modes

  complex(kind=DP) :: ZDOTC, dynwrk (3 * nat, 3 * nat)
  complex(kind=DP), allocatable :: dvloc (:)
  ! the scalar product functions
  ! auxiliary dynamical matrix
  ! d Vloc / dtau

  allocate (dvloc( nrxxs))    

  call setv (2 * 3 * nat * 3 * nat, 0.d0, dynwrk, 1)
  !
  ! We need a sum over all perturbations
  !

  do nu_j = 1, 3 * nat
     call compute_dvloc (nu_j, dvloc)
     do ipert = 1, nper
        nu_i = nu_i0 + ipert
        do is = 1, nspin
           dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + ZDOTC (nrxxs, drhoscf &
                (1, is, ipert), 1, dvloc, 1) * omega / (nr1s * nr2s * nr3s)
        enddo
     enddo

  enddo
#ifdef PARA
  !
  ! collect contributions from nodes of a pool (sum over G & R space)
  !

  call reduce (18 * nat * nat, dynwrk)
#endif

  call DAXPY (18 * nat * nat, 1.d0, dynwrk, 1, wdyn, 1)
  deallocate(dvloc)
  return
end subroutine drhodvloc
