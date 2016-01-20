!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE clean_pw_ph(iq)
  !-----------------------------------------------------------------------
  !
  ! This routine deallocate all the variables of pwscf and of the
  ! phonon code, and reset the same variables as after reading input in
  ! phq_readin, so that it is possible to start a calculation at
  ! a new q.
  !
  USE kinds,           ONLY : DP
  USE modes,           ONLY : nirr
  USE partial,         ONLY : done_irr
  USE disp,            ONLY : done_iq
  USE control_ph,      ONLY : rec_code_read
  USE save_ph,         ONLY : restore_ph_input_variables

  USE lr_symm_base,    ONLY : nsymq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !
  INTEGER :: irr
  !
  done_iq(iq)=.TRUE.
  DO irr=1,nirr
     IF (.NOT.done_irr(irr)) done_iq(iq)=.FALSE.
  ENDDO
  CALL clean_pw( .FALSE. )
  CALL deallocate_phq()
  rec_code_read=-1000
  nsymq=0
  !
  ! ... Close the files
  !
  CALL close_phq( .TRUE. )
  !
  CALL restore_ph_input_variables()
  !
RETURN
END SUBROUTINE clean_pw_ph
