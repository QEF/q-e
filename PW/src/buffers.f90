!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE buffers
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  PRIVATE
  PUBLIC :: open_buffer, get_buffer, save_buffer, close_buffer
  !
  SAVE
  !
  ! ... global variables: buffer1 is the memory buffer
  !                       nword_ the record length
  !                       buffered is true when buffer has been written
  !                       at least once - read from file otherwise
  COMPLEX(DP), ALLOCATABLE :: buffer1(:,:)
  INTEGER :: nword_
  LOGICAL :: buffered = .FALSE.
  !
  CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE open_buffer (unit, extension, nword, maxrec, exst)
  !-----------------------------------------------------------------------
  !
  !     unit > 6 : connect unit "unit" to file $wfc_dir/$prefix."extension" 
  !     for direct I/O access, with record length = nword complex numbers;
  !     maxrec is ignored, exst=T(F) if the file (does not) exists
  !
  !     unit =-10: in addition to opening unit 10 as above, allocate a buffer
  !     for storing up to maxrec records of length nword complex numbers
  !
  USE io_files,  ONLY : diropn, wfc_dir
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: extension
  INTEGER, INTENT(IN) :: unit, nword, maxrec
  LOGICAL, INTENT(OUT) :: exst
  !
  INTEGER :: ierr
  !
  exst = .FALSE.
  IF ( unit == -10) THEN
     exst =  ALLOCATED ( buffer1 )
     IF ( exst ) THEN
        CALL infomsg ('open_buffer', 'buffer already allocated')
     ELSE
        nword_ = nword
        ALLOCATE ( buffer1 ( nword, maxrec ) )
        buffered = .FALSE.
     END IF
  END IF
  !
  IF ( unit == -10 .OR. unit > 6 ) THEN
     CALL diropn ( ABS(unit), extension, 2*nword, exst, wfc_dir )
  ELSE
     CALL errore ('open_buffer', 'incorrect unit specified', ABS(unit))
  END IF
  !
  RETURN
  !
END SUBROUTINE open_buffer
!
!----------------------------------------------------------------------------
SUBROUTINE save_buffer( vect, nword, unit, nrec )
  !----------------------------------------------------------------------------
  !
  ! ... copy vect(1:nword) into the "nrec"-th record of
  ! ... - a previously allocated buffer, if unit = -10
  ! ... - a previously opened direct-access file with unit > 6
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword, unit, nrec
  COMPLEX(DP), INTENT(IN) :: vect(nword)
  !
  IF ( unit == -10 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        !
        IF ( nrec > SIZE ( buffer1, 2) )  &
           CALL errore ('save_buffer', 'too many records', ABS(nrec))
        IF ( nword /= SIZE ( buffer1, 1) )  &
           CALL errore ('save_buffer', 'record length mismatch', ABS(nword))
        !
        buffer1(:,nrec) = vect(:)
        buffered = .TRUE.
        !
     ELSE
        !
        CALL errore ('save_buffer', 'buffer not allocated', ABS(unit))
        !
     END IF
     !
  ELSE IF ( unit > 6 ) THEN
     !
     CALL davcio ( vect, 2*nword, unit, nrec, +1 )
     !
  ELSE
     !
     CALL errore ('save_buffer', 'incorrect unit specified', ABS(unit))
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE save_buffer
!
!----------------------------------------------------------------------------
SUBROUTINE get_buffer( vect, nword, unit, nrec )
  !----------------------------------------------------------------------------
  !
  ! ... copy vect(1:nword) from the "nrec"-th record of
  ! ... - a previously allocated buffer, if unit = -10 ;
  ! ...   if buffer never written: read it from unit 10
  ! ... - a previously opened direct-access file with unit > 6
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nword, unit, nrec
  COMPLEX(DP), INTENT(OUT) :: vect(nword)
  !
  IF ( unit == -10 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        IF ( nrec > SIZE ( buffer1, 2) )  &
           CALL errore ('get_buffer', 'no such record', ABS(nrec))
        IF ( nword /= SIZE ( buffer1, 1) )  &
           CALL errore ('get_buffer', 'record length mismatch', ABS(nword))
        IF ( buffered ) THEN
           vect(:) = buffer1(:,nrec)
        ELSE
           CALL davcio ( vect, 2*nword, ABS(unit), nrec, -1 )
        END IF
     ELSE
        CALL errore ('get_buffer', 'buffer not allocated', ABS(unit))
     END IF
     !
  ELSE IF ( unit > 6 ) THEN
     !
     CALL davcio ( vect, 2*nword, unit, nrec, -1 )
     !
  ELSE
     !
     CALL errore ('get_buffer', 'incorrect unit specified', ABS(unit))
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE get_buffer
!
SUBROUTINE close_buffer ( unit, status )
  !
  !     unit > 6 : close unit with status "status" ('keep' or 'delete')
  !     unit =-10: deallocate buffer
  !                if status='keep', save to previosly opened file
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: unit
  CHARACTER(LEN=*), INTENT(IN) :: status
  !
  INTEGER :: i
  LOGICAL :: exst, opnd
  !
  IF ( unit == -10 ) THEN
     !
     IF ( ALLOCATED ( buffer1 ) ) THEN
        !
        IF ( TRIM(status) == 'KEEP' .OR. TRIM(status) == 'keep') THEN
           DO i = 1, SIZE (buffer1, 2)
              CALL davcio ( buffer1(1,i), 2*nword_, ABS(unit), i, +1 )
           END DO
           CLOSE( UNIT = ABS(unit), STATUS = status )
        END IF
        !
        DEALLOCATE (buffer1)
        buffered = .FALSE.
        !
     ELSE
        !
        CALL infomsg ('close_buffer', 'buffer not allocated')
        !
     END IF
     !
  END IF
  !
  IF ( unit == -10 .OR. unit > 6 ) THEN
     !
     INQUIRE( UNIT = ABS(unit), OPENED = opnd )
     IF ( opnd ) CLOSE( UNIT = ABS(unit), STATUS = status )
     !
  ELSE
     !
     CALL infomsg ('get_buffer', 'incorrect unit specified')
     !
  END IF
  !
END SUBROUTINE close_buffer
!
END MODULE buffers
