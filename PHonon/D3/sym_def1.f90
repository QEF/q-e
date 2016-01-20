!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine sym_def1 (def, irr)
  !---------------------------------------------------------------------
  ! Symmetrizes the first order changes of the Fermi energies of an
  ! irreducible representation. These objects are defined complex because
  ! perturbations may be complex
  !
  ! Used in the q=0 metallic case only.
  !
  USE kinds, only : DP
  use pwcom
  use phcom
  use d3com

  USE lr_symm_base, ONLY : nsymq, irgq

  implicit none

  integer :: irr
  ! input: the representation under consideration

  complex (DP) :: def (npertx)
  ! inp/out: the fermi energy changes

  integer :: ipert, jpert, isym, irot
  ! counter on perturbations
  ! counter on perturbations
  ! counter on symmetries
  ! the rotation

  complex (DP) :: w_def (npertx)
  ! the fermi energy changes (work array)

  do ipert = 1, npertg0 (irr)
     def (ipert) =  DBLE (def (ipert) )
  enddo
  if (nsymq == 1) return
  !
  ! Here we symmetrize with respect to the small group of q
  !
  w_def (:) = (0.d0, 0.d0)
  do ipert = 1, npertg0 (irr)
     do isym = 1, nsymg0
        irot = irgq (isym)
        do jpert = 1, npertg0 (irr)
           w_def (ipert) = w_def (ipert) + tg0 (jpert, ipert, irot, irr) &
                * def (jpert)
        enddo
     enddo
  enddo
  !
  ! normalize and exit
  !
  def (:) = w_def(:) / DBLE(nsymq)

  return
end subroutine sym_def1
