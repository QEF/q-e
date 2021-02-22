!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kc_initialize_ph()
  !-----------------------------------------------------------------------
  !
  !! This routine initialize all the relevant quantity for the LR calculation
  !! of the screening parameter. 
  !
  USE klist,           ONLY : nks
  USE qpoint,          ONLY : nksq, ikks, ikqs
  USE control_ph,      ONLY : all_done
  USE control_lr,      ONLY : lgamma
  USE acfdtest,        ONLY : skip_ph 
  USE ph_restart,      ONLY : ph_writefile
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
  !  Allocate the phonon variables
  CALL allocate_phq()
  !
  !  Set the main control variable of the phonon code
  skip_ph =.TRUE. 
  CALL kc_phq_setup()
  skip_ph =.FALSE. 
  !
  all_done =.FALSE. 
  !
  !  Open the files of the phonon code
  CALL kc_openfilq()
  !
  !  Initialize all quantities which do not depend on the 
  !  linear response to the perturbation
  CALL phq_init()
  !
  !
  CALL print_clock( 'NSCF' )
  !
  RETURN
  
END SUBROUTINE kc_initialize_ph
