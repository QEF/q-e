!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE sym_def(def)
  !---------------------------------------------------------------------
  !! Symmetrizes the first order changes of the Fermi energies of an
  !! irreducible representation. These objects are defined complex because
  !! perturbations may be complex.
  !
  !! Used in the q=0 metallic case only.
  !
  USE kinds,        ONLY : DP
  USE control_lr,   ONLY : lgamma_gamma
  USE lr_symm_base, ONLY : minus_q, nsymq, lr_npert, upert, upert_mq
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout) :: def(3)
  !! inp/out: the fermi energy changes.
  !! NB: def(3) should be def(npertx), but it is used only at Gamma
  !!     where the dimension of irreps never exceeds 3.
  !
  ! ... local variables
  !
  INTEGER :: ipert, jpert, isym
  ! counter on perturbations
  ! counter on perturbations
  ! counter on symmetries
  !
  COMPLEX(DP) :: w_def(3)
  ! the fermi energy changes (work array)
  !
  IF (lgamma_gamma) RETURN
  if (nsymq == 1 .and. (.not.minus_q) ) return
  if (lr_npert > 3) CALL errore("sym_def", "lr_npert cannot exceed 3 at q=0", 1)
  !
  ! first the symmetrization   S(irotmq)*q = -q + Gi if necessary
  !
  if (minus_q) then
     w_def = (0.d0, 0.d0)
     do ipert = 1, lr_npert
        do jpert = 1, lr_npert
           w_def(ipert) = w_def(ipert) + upert_mq(jpert, ipert) * def(jpert)
        enddo
     enddo
     do ipert = 1, lr_npert
        def(ipert) = 0.5d0 * (def(ipert) + CONJG(w_def(ipert)) )
     enddo
  endif
  !
  ! Here we symmetrize with respect to the small group of q
  !
  w_def = (0.d0, 0.d0)
  do ipert = 1, lr_npert
     do isym = 1, nsymq
        do jpert = 1, lr_npert
           w_def(ipert) = w_def(ipert) + upert(jpert, ipert, isym) * def(jpert)
        enddo
     enddo
  enddo
  !
  ! normalize and exit
  !
  def = w_def / DBLE(nsymq)
  !
END SUBROUTINE sym_def
