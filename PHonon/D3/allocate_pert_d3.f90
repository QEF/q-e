!
! Copyright (C) 2001-2009 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_pert_d3()
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities depending on the
  ! maximum number of perturbations
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat

  USE modes, ONLY : npertx, t, tmq
  USE modesg0, ONLY : tg0

  USE control_lr, ONLY : lgamma

  implicit none
  !
  !  allocate space for the quantities with dimensions that depend
  !  on the maximum number of perturbations
  !
  ALLOCATE (t (npertx, npertx, 48, 3*nat))
  ALLOCATE (tmq (npertx, npertx, 3*nat))
  IF (lgamma) THEN
     tg0 => t
  ELSE
     allocate (tg0( npertx, npertx, 48, 3*nat))
  ENDIF
  RETURN
END SUBROUTINE allocate_pert_d3
