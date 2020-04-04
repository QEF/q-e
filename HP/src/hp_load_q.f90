!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hp_load_q()
  !-----------------------------------------------------------------------
  !
  ! This is a driver to the HP initialization routines.
  !
  USE klist,        ONLY : nks
  USE io_global,    ONLY : stdout
  USE qpoint,       ONLY : nksq, ikks, ikqs
  USE control_lr,   ONLY : lgamma
  USE ldaU_hp,      ONLY : code
  !
  IMPLICIT NONE
  INTEGER :: ik
  !
  ! ... nksq is the number of k-points, NOT including k+q points
  !
  IF ( lgamma ) THEN
     !
     nksq = nks
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik=1,nksq
        ikks(ik) = ik
        ikqs(ik) = ik
     ENDDO
     !
  ELSE
     !
     nksq = nks / 2
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik=1,nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  END IF
  !
  ! Allocate various arrays
  !
  CALL hp_allocate_q()
  !
  ! Setup various control variables
  !
  CALL hp_setup_q()
  !
  ! Output summary of the main variables
  !
  CALL hp_summary_q()
  !
  ! Open all necessary files
  !
  CALL hp_openfil_q()
  !
  ! Initialize all quantities which do not depend on the
  ! linear response to the perturbation
  !
  CALL hp_init_q()
  !
  WRITE( stdout, '(/5x,"Total time spent up to now is:")')
  !
  CALL print_clock (code)
  !
  RETURN
  !
END SUBROUTINE hp_load_q
