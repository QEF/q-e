!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE ions()
  !----------------------------------------------------------------------------
  !
  USE control_flags, ONLY : conv_ions, restart, lscf, lmd, lbfgs
  USE force_mod,     ONLY : lforce, lstres
  !
  IMPLICIT NONE
  !
  CALL start_clock( 'ions' )
  !
  conv_ions = .TRUE.
  !
  ! ... recover from a previous run, if appropriate
  !
  IF ( restart .AND. lscf ) CALL restart_in_ions()
  !
  IF ( lforce ) CALL forces()
  !
  IF ( lstres ) CALL stress()
  !
  IF ( lmd .OR. lbfgs ) THEN
     !
     ! ... first we move the ions
     !
     CALL move_ions()
     !
     ! ... then we save restart information for the new configuration
     !
     CALL punch( 'config' )
     !
     CALL save_in_ions()
     !
  END IF
  !
  CALL stop_clock( 'ions' )
  !
  RETURN
  !
END SUBROUTINE ions
