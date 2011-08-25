!
! Copyright (C) 2001-2009 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_pert()
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities depending on the
  ! maximum number of perturbations npertx
  !
  USE ions_base, ONLY : nat

  USE modes, ONLY : npertx, t, tmq

  implicit none
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  ALLOCATE ( t ( npertx, npertx, 48, 3 * nat ) )
  ALLOCATE ( tmq ( npertx, npertx, 3 * nat ) )

  RETURN
END SUBROUTINE allocate_pert

!-----------------------------------------------------------------------
subroutine deallocate_pert()
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities depending on the
  ! maximum number of perturbations npertx
  !
  USE modes, ONLY : t, tmq

  IMPLICIT NONE
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  IF (ASSOCIATED(t)) DEALLOCATE ( t )
  IF (ASSOCIATED(tmq)) DEALLOCATE ( tmq )

  RETURN
END SUBROUTINE deallocate_pert
