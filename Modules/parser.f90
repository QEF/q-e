!
! Copyright (C) 2001 Carlo Cavazzoni and PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module parser

  USE kinds

  INTERFACE parser_error
     MODULE PROCEDURE p_err_I, p_err_R, p_err_S, p_err_L
  END INTERFACE

contains

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
            CALL error(' read_line ', ' input line too short ', LEN( line ) )
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
            CALL error(' read_line ', ' end of file ', 0 )
          ELSE
            IF( PRESENT(field) ) CALL field_compare(line, nfield, field)
          END IF

          RETURN
        END SUBROUTINE

!SUBROUTINE cpitoa: turn an integer into a string
!SUBROUTINE unitname: ***
!SUBROUTINE myunitname: given the processor num. return the output unit num.
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


  SUBROUTINE cpitoa(n,str)
!   trasforma un numero intero in una stringa
    integer, intent(in) :: n
    character(LEN=*) str
    integer i, npow, j, nq, ntmp
    logical :: lzero

    j     = 1
    IF( n .eq. 0) THEN
      str(j:j) = '0'
      j = j + 1
    ELSE
      IF( n .lt. 0 ) THEN
        str(j:j) = '-'
        j = j + 1
        ntmp  = -n
      ELSE
        ntmp  =  n
      END IF
      lzero = .FALSE.
      do i = 9, 0, -1
        npow = 10**i
        nq   = ntmp / npow
        ntmp = MOD(ntmp,npow)
        if( lzero .or. (nq.ne.0)) then
          lzero = .TRUE.
          select case(nq)
            case (9)
              str(j:j) = '9'
            case (8)
              str(j:j) = '8'
            case (7)
              str(j:j) = '7'
            case (6)
              str(j:j) = '6'
            case (5)
              str(j:j) = '5'
            case (4)
              str(j:j) = '4'
            case (3)
              str(j:j) = '3'
            case (2)
              str(j:j) = '2'
            case (1)
              str(j:j) = '1'
            case (0)
              str(j:j) = '0'
            case default
              str(j:j) = '0'
          end select
          j = j + 1
        end if
      end do
    END IF
    str(j:j) = ' '
    return
  END SUBROUTINE


  SUBROUTINE unitname(n,name)
    character(len=*) :: name
    character(len=5) :: num
    integer :: l
    call cpitoa(n,num)
    l    = index(num,' ') - 1
    name = 'fort.'//num(1:l)
    return
  END SUBROUTINE

  SUBROUTINE myunitname(me,name)
    character(len=*) :: name
    character(len=5) :: num
    integer :: l
    call cpitoa(me,num)
    l    = index(num,' ') - 1
    name = 'fort_6.'//num(1:l)
    return
  END SUBROUTINE

  SUBROUTINE field_compare(str, nf, var)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: var
    INTEGER, INTENT(IN) :: nf

    CHARACTER(LEN=*), INTENT(OUT) :: str
    INTEGER :: nc

    CALL field_count(nc, str)
    IF(nc .LT. nf) THEN
      CALL parser_error(' field_compare ', ' wrong number of fields ', var)
    END IF
  END SUBROUTINE field_compare

  SUBROUTINE con_cam(num, line, car)
    CHARACTER(LEN=*) :: line
    CHARACTER(LEN=1) :: sep
    CHARACTER(LEN=1), OPTIONAL :: car
    INTEGER :: num, j

    num = 0
    IF (len(line) .GT. 256 ) THEN
       WRITE(*,*) 'riga ', line
       WRITE(*,*) 'lunga ', len(line)
       num = -1
       RETURN
    END IF

    WRITE(*,*) '1riga ', line
    WRITE(*,*) '1lunga ', len(line)
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

  SUBROUTINE p_err_I(a,b,n)
      USE mp
      USE mp_global, ONLY: mpime
      USE parameters

      IMPLICIT NONE

! ... declare subroutine arguments
      CHARACTER(LEN=*) a, b
      INTEGER n

! ... declare function

! ... print the error message
      WRITE (6,100) mpime, a, b, n

! ... terminate the program
      OPEN(UNIT=15, FILE='CRASH', POSITION='append', STATUS='unknown')
      WRITE (15,100) mpime, a, b, n
      CALL cpflush  ! flush output streams
      CALL mp_end   ! terminate MPI

100   FORMAT (/,' *** from PE : ',I3,'  *** in routine ',A, &
              /,' *** error msg. : ',A,' *** code ',I5, &
              /,' *** aborting ***', /)
      STOP
  END SUBROUTINE p_err_I

  SUBROUTINE p_err_R(a,b,r)

      USE mp
      USE mp_global, ONLY: mpime
      USE parameters

      IMPLICIT NONE

! ... declare subroutine arguments
      CHARACTER(LEN=*) a, b
      REAL(DBL) r

! ... declare function

! ... print the error message
      WRITE (6,100) mpime, a, b, r

! ... terminate the program
      OPEN(UNIT=15, FILE='CRASH', POSITION='append', STATUS='unknown')
      WRITE (15,100) mpime, a, b, r
      CALL cpflush  ! flush output streams
      CALL mp_end   ! terminate MPI


100   FORMAT (/,' *** from PE : ',I3,'  *** in routine ',A, &
              /,' *** error msg. : ',A,' *** code ',F16.8, &
              /,' *** aborting ***', /)
      STOP
  END SUBROUTINE p_err_R

  SUBROUTINE p_err_S(a,b,c)

      USE mp
      USE mp_global, ONLY: mpime
      USE parameters

      IMPLICIT NONE

! ... declare subroutine arguments
      CHARACTER(LEN=*) a, b, c

! ... print the error message
      WRITE (6,100) mpime, a, b, c

! ... terminate the program
      OPEN(UNIT=15, FILE='CRASH', POSITION='append', STATUS='unknown')
      WRITE (15,100) mpime, a, b, c
      CALL cpflush  ! flush output streams
      CALL mp_end   ! terminate MPI


100   FORMAT (/,' *** from PE : ',I3,'  *** in routine ',A, &
              /,' *** error msg. : ',A,' *** code ',A, &
              /,' *** aborting ***', /)
      STOP
  END SUBROUTINE p_err_S

  SUBROUTINE p_err_L(a,b,t)

      USE mp
      USE mp_global, ONLY: mpime
      USE parameters

      IMPLICIT NONE

! ... declare subroutine arguments
      CHARACTER(LEN=*) a, b
      LOGICAL t

! ... declare function

! ... print the error message
      WRITE (6,100) mpime, a, b, t

! ... terminate the program
      OPEN(UNIT=15, FILE='CRASH', POSITION='append', STATUS='unknown')
      WRITE (15,100) mpime, a, b, t
      CALL cpflush  ! flush output streams
      CALL mp_end   ! terminate MPI

100   FORMAT (/,' *** from PE : ',I3,'  *** in routine ',A, &
              /,' *** error msg. : ',A,' *** code ',L2, &
              /,' *** aborting ***', /)
      STOP
  END SUBROUTINE p_err_L


END MODULE parser
