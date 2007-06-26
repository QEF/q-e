!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE davcio( vect, nword, unit, nrec, io )
  !----------------------------------------------------------------------------
  !
  ! ... direct-access vector input/output
  ! ... read/write nword words starting from the address specified by vect
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword, unit, nrec, io
    ! input: the dimension of vect
    ! input: the unit where to read/write
    ! input: the record where to read/write
    ! input: flag if < 0 reading if > 0 writing
  REAL(DP), INTENT(INOUT) :: vect(nword)
   ! input/output: the vector to read/write
  !
  INTEGER :: ios
    ! integer variable for I/O control
  LOGICAL :: opnd
  !
  !
  CALL start_clock( 'davcio' )
  !
  INQUIRE( UNIT = unit )
  !
  IF ( unit  <= 0 ) CALL errore(  'davcio', 'wrong unit', 1 )
  IF ( nrec  <= 0 ) CALL errore(  'davcio', 'wrong record number', 2 )
  IF ( nword <= 0 ) CALL errore(  'davcio', 'wrong record length', 3 )
  IF ( io    == 0 ) CALL infomsg( 'davcio', 'nothing to do?' )
  !
  INQUIRE( UNIT = unit, OPENED = opnd )
  !
  IF ( .NOT. opnd ) &
     CALL errore(  'davcio', 'unit is not opened', unit )
  !
  ios = 0
  !
  IF ( io < 0 ) THEN
     !
     READ( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     !
  ELSE IF ( io > 0 ) THEN
     !
     WRITE( UNIT = unit, REC = nrec, IOSTAT = ios ) vect
     !
  END IF
  !
  IF ( ios /= 0 ) THEN
     !
     WRITE( stdout, '(1X,"IOS = ",I2)' ) ios
     !
     CALL errore( 'davcio', 'i/o error in davcio', unit )
     !
  END IF
  !
  CALL stop_clock( 'davcio' )
  !
  RETURN
  !
END SUBROUTINE davcio
