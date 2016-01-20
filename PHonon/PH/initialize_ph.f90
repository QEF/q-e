!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE initialize_ph()
  !-----------------------------------------------------------------------
  !
  ! This is a driver to the phonon initialization routines.
  !
  USE klist,  ONLY : nks, nkstot
  !
  USE qpoint, ONLY : nksq, nksqtot, ikks, ikqs
  USE control_lr, ONLY : lgamma
  !
  IMPLICIT NONE
  INTEGER :: ik
  !
  ! ... nksq is the number of k-points, NOT including k+q points
  !
  IF ( lgamma ) THEN
     !
     nksq = nks
     nksqtot = nkstot
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik=1,nksq
        ikks(ik) = ik
        ikqs(ik) = ik
     ENDDO
     !
  ELSE
     !
     nksq = nks / 2
     nksqtot = nkstot / 2
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik=1,nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  END IF
  !
  !  Allocate the phonon variables
  !
  CALL allocate_phq()
  !
  !  Set the main control variable of the phonon code
  !
  CALL phq_setup()
  !
  !  Recover the status if available
  !
  CALL phq_recover()
  !
  !  Output summary of the main variables of the phonon code
  !
  CALL phq_summary()
  !
  !  Open the files of the phonon code
  !
  CALL openfilq()
  !
  !  Initialize all quantities which do not depend on the
  !  linear response to the perturbation
  !
  CALL phq_init()
  !
  CALL print_clock( 'PHONON' )
  !
  RETURN

END SUBROUTINE initialize_ph
