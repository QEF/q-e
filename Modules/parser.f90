!
! Copyright (C) 2001-2004 Carlo Cavazzoni and PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! ... SUBROUTINE field_count:   accepts two string (one of them is optional) 
!                               and one integer and count the number of fields
!                               in the string separated by a blank or a tab 
!                               character. If the optional string is specified
!                               (it has anyway len=1) it is assumed as the 
!                               separator character.
!                               Ignores any character following the exclamation
!                               mark (fortran comment)
!
! ... SUBROUTINE con_cam:       counts the number of fields in a string 
!                               separated by the optional character
!
! ... SUBROUTINE field_compare: accepts two strings and one integer. Counts the
!                               fields contained in the first string and 
!                               compares it with the integer. 
!                               If they are less than the integer calls the 
!                               routine error and show by the second string the
!                               name of the field where read-error occurred.
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE parser
  !----------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE kinds
  !
  PRIVATE
  !
  PUBLIC :: parse_unit, field_count, read_line
  !
  INTEGER :: parse_unit = 5 ! normally 5, but can be set otherwise
  !
  CONTAINS
  !
  !
  !--------------------------------------------------------------------------
  PURE SUBROUTINE field_count( num, line, car )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER,                    INTENT(OUT) :: num
    CHARACTER(LEN=*),           INTENT(IN)  :: line
    CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: car
#if defined (__XLF)
    ! ... with the IBM xlf compiler some combination of flags lead to
    ! ... variables being defined as static, hence giving a conflict
    ! ... with PURE function. We then force the variable to be AUTOMATIC
    CHARACTER(LEN=1), AUTOMATIC             :: sep1, sep2    
    INTEGER, AUTOMATIC                      :: j
#else
    CHARACTER(LEN=1)                        :: sep1, sep2    
    INTEGER                                 :: j
#endif
    !
    !
    num = 0
    !
    IF ( .NOT. present(car) ) THEN
       !
       sep1 = char(32)  ! ... blank character
       sep2 = char(9)   ! ... tab character
       !
       DO j = 2, MAX( LEN( line ), 256 )
          !
          IF ( line(j:j) == '!' .OR. line(j:j) == char(0) ) THEN
             !
             IF ( line(j-1:j-1) /= sep1 .AND. line(j-1:j-1) /= sep2 ) THEN
                !
                num = num + 1
                !
             END IF   
             !
             EXIT
             !
          END IF
          !
          IF ( ( line(j:j) == sep1 .OR. line(j:j) == sep2 ) .AND. &
               ( line(j-1:j-1) /= sep1 .AND. line(j-1:j-1) /= sep2 ) ) THEN
             !
             num = num + 1
             !
          END IF
          !
       END DO
       !
    ELSE
       !
       sep1 = car
       !
       DO j = 2, MAX( LEN( line ), 256 )
          ! 
          IF ( line(j:j) == '!' .OR. &
               line(j:j) == char(0) .OR. line(j:j) == char(32) ) THEN
             !
             IF ( line(j-1:j-1) /= sep1 ) num = num + 1
             !
             EXIT
             !
          END IF
          !
          IF ( line(j:j) == sep1 .AND. line(j-1:j-1) /= sep1 ) num = num + 1
          !
       END DO
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE field_count
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_line( line, nfield, field, end_of_file )
    !--------------------------------------------------------------------------
    !
    USE mp,        ONLY : mp_bcast
    USE mp_global, ONLY : world_comm
    USE io_global, ONLY : ionode, ionode_id
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),           INTENT(OUT) :: line
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: field
    INTEGER,          OPTIONAL, INTENT(IN)  :: nfield
    LOGICAL,          OPTIONAL, INTENT(OUT) :: end_of_file
    LOGICAL                                 :: tend
    !
    !
    IF( LEN( line ) < 256 ) THEN
       CALL errore(' read_line ', ' input line too short ', MAX(LEN(line),1) )
    END IF
    !
    IF ( ionode ) THEN
30     READ (parse_unit, fmt='(A256)', ERR=10, END=10) line
       IF( line == ' ' .OR. line(1:1) == '#' ) GO TO 30
       tend = .FALSE.
       GO TO 20
10     tend = .TRUE.
20     CONTINUE
    END IF
    !
    CALL mp_bcast( tend, ionode_id, world_comm )
    CALL mp_bcast( line, ionode_id, world_comm )
    !
    IF( PRESENT(end_of_file) ) THEN
       end_of_file = tend
    ELSE IF( tend ) THEN
       CALL infomsg(' read_line ', ' end of file ', -1 )
    ELSE
       IF( PRESENT(field) ) CALL field_compare( line, nfield, field )
    END IF
    !
    RETURN
    !
  END SUBROUTINE read_line
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE field_compare( str, nf, var )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER,          INTENT(IN) :: nf
    CHARACTER(LEN=*), INTENT(IN) :: str
    INTEGER                      :: nc
    !
    CALL field_count( nc, str )
    !
    IF( nc < nf ) &
      CALL errore( ' field_compare ', &
                 & ' wrong number of fields: ' // TRIM( var ), 1 )
    !
    RETURN
    !
  END SUBROUTINE field_compare
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE con_cam(num, line, car)
    !--------------------------------------------------------------------------
    CHARACTER(LEN=*) :: line
    CHARACTER(LEN=1) :: sep
    CHARACTER(LEN=1), OPTIONAL :: car
    INTEGER :: num, j

    num = 0
    IF (len(line) .GT. 256 ) THEN
       WRITE( stdout,*) 'riga ', line
       WRITE( stdout,*) 'lunga ', len(line)
       num = -1
       RETURN
    END IF

    WRITE( stdout,*) '1riga ', line
    WRITE( stdout,*) '1lunga ', len(line)
    IF ( .NOT. present(car) ) THEN
       sep=char(32)             !char(32) is the blank character
    ELSE
       sep=car
    END IF

    DO j=2, MAX(len(line),256)
       IF ( line(j:j) == '!' .OR. line(j:j) == char(0)) THEN
          RETURN
       END IF
       IF ( (line(j:j) .EQ. sep) .AND. &
            (line(j-1:j-1) .NE. sep) )  THEN
          num = num + 1
       END IF
    END DO
    RETURN
  END SUBROUTINE con_cam
  !
END MODULE parser
