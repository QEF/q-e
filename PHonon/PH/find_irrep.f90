!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE find_irrep()
  !---------------------------------------------------------------------
  !
  !  Computes the variables needed to pass to the pattern representation
  !     u      the patterns
  !     nirr   the number of irreducible representation
  !     npert  the dimension of each irreducible representation
  !
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat
  USE symm_base,     ONLY : nsym
  USE control_ph,    ONLY : lgamma_gamma
  USE modes,         ONLY : u, npert, nirr
  USE qpoint,        ONLY : xq
  USE control_flags, ONLY : modenum

  IMPLICIT NONE

  REAL(DP) :: w2(3*nat)

  IF (nsym > 1.AND..NOT.lgamma_gamma.AND.modenum==0) THEN
     CALL set_irr_new (xq, u, npert, nirr, w2)
  ELSE
     CALL set_irr_nosym_new (u, npert, nirr)
  ENDIF

  RETURN
  END SUBROUTINE find_irrep

!-----------------------------------------------------------------------
SUBROUTINE find_irrep_sym()
  !-----------------------------------------------------------------------
  !
  !   Computes the variables needed to symmetrize in the pattern representation
  !     t      the matrices of the small group of q on the pattern basis
  !     tmq    the matrix of the symmetry which sends q -> -q + G
  !
  !
  USE kinds,         ONLY : DP
  USE control_ph,    ONLY : lgamma_gamma
  USE symm_base,     ONLY : nsym
  USE modes,         ONLY : npertx, npert, nirr, t, tmq

  IMPLICIT NONE

  INTEGER :: irr
  ! counters

  IF (lgamma_gamma) RETURN

  npertx = 0
  DO irr = 1, nirr
     npertx = max (npertx, npert (irr) )
  ENDDO
  CALL allocate_pert()
  CALL set_irr_sym_new (t, tmq, npertx )

  RETURN
END SUBROUTINE find_irrep_sym
