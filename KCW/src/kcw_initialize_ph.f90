!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kcw_initialize_ph()
  !-----------------------------------------------------------------------
  !
  !! This routine initialize all the relevant quantity for the LR calculation
  !! of the screening parameter. 
  !
  USE klist,             ONLY : nks, nkstot
  USE qpoint,            ONLY : nksq, ikks, ikqs, nksqtot
  USE control_lr,        ONLY : lgamma
  USE noncollin_module,  ONLY : domag, noncolin
  USE qpoint_aux,        ONLY : ikmks, ikmkmqs
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
  !  Allocate the variables for the LR calculation at q
  CALL kcw_allocate_q()
  !
  !  Set the main control variable of the LR calculation at q
  CALL kcw_q_setup()
  !
  !  Open the files for the LR calculation
  CALL kcw_openfilq()
  !
  !  Initialize all quantities which do not depend on the 
  !  linear response to the perturbation
  CALL kcw_init_q()
  !
  !
  CALL print_clock( 'NSCF' )
  !
  RETURN
  
END SUBROUTINE kcw_initialize_ph
