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
  !! This is a driver to the phonon initialization routines.
  !
  USE klist,  ONLY : nks, nkstot
  !
  USE qpoint, ONLY : nksq, nksqtot, ikks, ikqs
  USE qpoint_aux, ONLY : ikmks, ikmkmqs
  USE control_lr, ONLY : lgamma
  USE noncollin_module, ONLY : noncolin, domag
  !
  IMPLICIT NONE
  INTEGER :: ik
  !
  ! ... nksq is the number of k-points, NOT including k+q points
  !
  IF ( lgamma ) THEN
     !
     IF (noncolin.AND.domag) THEN
        nksq = nks/2
        nksqtot = nkstot/2
        ALLOCATE(ikks(nksq), ikqs(nksq))
        ALLOCATE(ikmks(nksq), ikmkmqs(nksq))
        DO ik=1,nksq
           ikks(ik) = 2*ik-1
           ikqs(ik) = 2*ik-1
           ikmks(ik) = 2*ik
           ikmkmqs(ik) = 2*ik
        ENDDO
     ELSE
        nksq = nks
        nksqtot = nkstot
        ALLOCATE(ikks(nksq), ikqs(nksq))
        DO ik=1,nksq
           ikks(ik) = ik
           ikqs(ik) = ik
        ENDDO
     END IF
     !
  ELSE
     !
     IF (noncolin.AND.domag) THEN
        nksq = nks / 4
        nksqtot = nkstot / 4
        ALLOCATE(ikks(nksq), ikqs(nksq))
        ALLOCATE(ikmks(nksq), ikmkmqs(nksq))
        DO ik=1,nksq
           ikks(ik) = 4 * ik - 3
           ikqs(ik) = 4 * ik - 2
           ikmks(ik) = 4 * ik - 1
           ikmkmqs(ik) = 4 * ik 
        ENDDO
     ELSE
        nksq = nks / 2
        nksqtot = nkstot / 2
        ALLOCATE(ikks(nksq), ikqs(nksq))
        DO ik=1,nksq
           ikks(ik) = 2 * ik - 1
           ikqs(ik) = 2 * ik
        ENDDO
     END IF
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
