!
! Copyright (C) 2001-2004 Carlo Cavazzoni and PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE parser
  !----------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  PURE FUNCTION int_to_char( int )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: int
    CHARACTER (LEN=6)   :: int_to_char
    !
    !   
    IF ( int < 10 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I1)" ) int
       !
    ELSE IF ( int < 100 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I2)" ) int
       !
    ELSE IF ( int < 1000 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I3)" ) int
       !
    ELSE IF ( int < 10000 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I4)" ) int
       !
    ELSE      
       ! 
       WRITE( UNIT = int_to_char , FMT = "(I5)" ) int     
       !
    END IF    
    !
    RETURN
    !
  END FUNCTION int_to_char
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE delete_if_present( filename )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*) :: filename
    LOGICAL          :: exst, opnd
    INTEGER          :: iunit
    !
    INQUIRE( FILE = filename, EXIST = exst )
    !
    IF ( exst ) THEN
       !
       unit_loop: DO iunit = 99, 1, - 1
          !
          INQUIRE( UNIT = iunit, OPENED = opnd )
          !
          IF ( .NOT. opnd ) THEN
             !
             OPEN(  UNIT = iunit, FILE = filename , STATUS = 'OLD' )
             CLOSE( UNIT = iunit, STATUS = 'DELETE' )
             WRITE( UNIT = stdout,         &
                    FMT = '(/,5X,"WARNING: ",A," file was present; old file deleted")' ) filename
             !
             RETURN
             !
          END IF
          !
       END DO unit_loop
       !
       CALL errore( 'delete_if_present', 'free unit not found ?!?', 1 )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE  
  !
  !
  !-----------------------------------------------------------------------
  logical function matches (string1, string2)  
    !-----------------------------------------------------------------------
    !
    ! .true. if string 1 is contained in string2, .false. otherwise
    !
    implicit none  
    character (len=*) :: string1, string2  
    integer :: len1, len2, l  


    len1 = len_trim(string1)  
    len2 = len_trim(string2)  
    do l = 1, len2 - len1 + 1  
       if (string1 (1:len1) .eq.string2 (l:l + len1 - 1) ) then  
          matches = .true.  
          return  
       endif
    enddo

    matches = .false.  
    return  
  end function matches
  !
  !-----------------------------------------------------------------------
  function capital (character)  
    !-----------------------------------------------------------------------
    !
    !   converts character to capital if lowercase
    !   copy character to output in all other cases
    !
    implicit none  
    character (len=1) :: capital, character
    !
    character(len=26) :: minuscole='abcdefghijklmnopqrstuvwxyz', &
                         maiuscole='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: i
    !
    do i=1,26
       if (character.eq.minuscole(i:i)) then
          capital=maiuscole(i:i)
          return
       end if
    end do
    capital = character  
    !
    return  
  end function capital

  SUBROUTINE field_count(num, line, car)
    CHARACTER(LEN=*) :: line
    CHARACTER(LEN=1) :: sep1, sep2
    CHARACTER(LEN=1), OPTIONAL :: car
    INTEGER :: num, j

    num = 0
    IF ( .NOT. present(car) ) THEN

      sep1=char(32)  !blank character
      sep2=char(9)   !tab character
      DO j=2, MAX(len(line),256)
       IF ( line(j:j) == '!' .OR. line(j:j) == char(0)) THEN
         IF ( (line(j-1:j-1) .NE. sep1) .AND. &
              (line(j-1:j-1) .NE. sep2) ) THEN
            num = num + 1
         END IF
         EXIT
       END IF
       IF ( ( (line(j:j) .EQ. sep1) .OR. &
              (line(j:j) .EQ. sep2) ) .AND. &
            ( (line(j-1:j-1) .NE. sep1) .AND. &
              (line(j-1:j-1) .NE. sep2) ) )  THEN
          num = num + 1
        END IF
      END DO
    ELSE
      sep1=car
      DO j=2, MAX(len(line),256)
        IF ( line(j:j) == '!' .OR. line(j:j) == char(0) .OR. &
          line(j:j) == char(32)) THEN
          IF( line(j-1:j-1) .NE. sep1 ) THEN
            num = num + 1
          END IF
          EXIT
        END IF
        IF ( (line(j:j) .EQ. sep1) .AND. &
             (line(j-1:j-1) .NE. sep1) )  THEN
          num = num + 1
        END IF
      END DO
    END IF

    RETURN

  END SUBROUTINE field_count


        SUBROUTINE read_line( line, nfield, field, end_of_file )
          USE mp, ONLY: mp_bcast
          USE mp_global, ONLY: group
          USE io_global, ONLY: ionode, ionode_id
          CHARACTER(LEN=*), INTENT(OUT) :: line
          CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field
          INTEGER, OPTIONAL, INTENT(IN) :: nfield
          LOGICAL, OPTIONAL, INTENT(OUT) :: end_of_file
          LOGICAL :: tend

          IF( LEN( line ) < 256 ) THEN
            CALL errore(' read_line ', ' input line too short ', LEN( line ) )
          END IF

          IF ( ionode ) THEN
            READ (5, fmt='(A256)', END=10) line
            tend = .FALSE.
            GO TO 20
 10         tend = .TRUE.
 20         CONTINUE
          END IF
          CALL mp_bcast(tend, ionode_id, group)
          CALL mp_bcast(line, ionode_id, group)

          IF( PRESENT(end_of_file) ) THEN
            end_of_file = tend
          ELSE IF( tend ) THEN
            CALL errore(' read_line ', ' end of file ', 0 )
          ELSE
            IF( PRESENT(field) ) CALL field_compare(line, nfield, field)
          END IF

          RETURN
        END SUBROUTINE

!SUBROUTINE con_cam:       count the number of fields in a string separated by
!                          the optional character
!SUBROUTINE field_count:   accept two string (one of them is optional) and one integer
!                          and count the number of fields in the string separated by a
!                          blank or a tab character. If the optional string is specified
!                          (it has anyway len=1) it is assumed as the separator character
!                          Ignore any charcter following the exclamation mark
!                          (fortran comment)
!SUBROUTINE field_compare: accept two strings and one integer. Count the
!                          fields contained in the first string and compare
!                          it with the integer. If they are less than the
!                          integer call the routine error and show by the
!                          second string the name of the field where read-error
!                          occurred.
!SUBROUTINE p_err(I,R,S,L): call the appropriate error subroutine (the same as
!                          flib/error.f90) depending on which kind of parameter
!                          has to be printed


  SUBROUTINE field_compare(str, nf, var)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER, INTENT(IN) :: nf

    CHARACTER(LEN=*), INTENT(OUT) :: str
    INTEGER :: nc

    CALL field_count(nc, str)
    IF(nc .LT. nf) THEN
      CALL errore(' field_compare ', ' wrong number of fields: ' // TRIM(var), 1 )
    END IF
  END SUBROUTINE field_compare

  SUBROUTINE con_cam(num, line, car)
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

END MODULE parser
