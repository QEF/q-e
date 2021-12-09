!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
MODULE sym_def_module
CONTAINS
subroutine sym_def (def, irr)
  !---------------------------------------------------------------------
  !! Symmetrizes the first order changes of the Fermi energies of an
  !! irreducible representation. These objects are defined complex because
  !! perturbations may be complex.
  !
  !! Used in the q=0 metallic case only.
  !
  USE kinds, only : DP
  USE modes,   ONLY : npert, t, tmq
  USE control_ph,           ONLY : lgamma_gamma

  USE lr_symm_base, ONLY : minus_q, nsymq

  implicit none

  integer :: irr
  !! input: the representation under consideration
  complex(DP) :: def(3)
  !! inp/out: the fermi energy changes.  
  !! NB: def(3) should be def(npertx), but it is used only at Gamma
  !!     where the dimension of irreps never exceeds 3.
  !
  ! ... local variables
  !
  integer :: ipert, jpert, isym, irot
  ! counter on perturbations
  ! counter on perturbations
  ! counter on symmetries
  ! the rotation

  complex(DP) :: w_def(3)
  ! the fermi energy changes (work array)

  IF (lgamma_gamma) RETURN
  if (nsymq == 1 .and. (.not.minus_q) ) return
  if (npert(irr) > 3) CALL errore("sym_def", "npert(irr) exceeds 3", 1)
  !
  ! first the symmetrization   S(irotmq)*q = -q + Gi if necessary
  !
  if (minus_q) then
     w_def = (0.d0, 0.d0)
     do ipert = 1, npert (irr)
        do jpert = 1, npert (irr)
           w_def (ipert) = w_def (ipert) + tmq (jpert, ipert, irr) &
                * def (jpert)
        enddo
     enddo
     do ipert = 1, npert (irr)
        def (ipert) = 0.5d0 * (def (ipert) + CONJG(w_def (ipert) ) )
     enddo
  endif
  !
  ! Here we symmetrize with respect to the small group of q
  !
  w_def = (0.d0, 0.d0)
  do ipert = 1, npert (irr)
     do isym = 1, nsymq
        irot = isym
        do jpert = 1, npert (irr)
           w_def (ipert) = w_def (ipert) + t (jpert, ipert, irot, irr) &
                * def (jpert)
        enddo
     enddo
  enddo
  !
  ! normalize and exit
  !
  def = w_def / DBLE(nsymq)

  return
end subroutine sym_def
END MODULE sym_def_module
