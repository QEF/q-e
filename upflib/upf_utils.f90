!
! Copyright (C) 2020 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE upf_utils
  !
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: capital, lowercase, isnumeric, matches, version_compare

! ... FUNCTION capital : converts a lowercase letter to uppercase
!                        returns input character if not a lowercase letter
! ... FUNCTION lowercase : as above, in reverse
! ... FUNCTION isnumeric : returns .true. if input character is a digit
!
! ... FUNCTION matches   : returns .true. if string1 matches string2  
!
! ... FUNCTION version_compare: Compare two version strings; the result can be
!                               "newer", "equal", "older", ""
! 
  CHARACTER(LEN=26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz', &
                                  upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

CONTAINS
  
  !-----------------------------------------------------------------------
FUNCTION capital( in_char )  
  !-----------------------------------------------------------------------
  !
  ! ... converts character to capital if lowercase
  ! ... copy character to output in all other cases
  !
  IMPLICIT NONE  
  !
  CHARACTER(LEN=1), INTENT(IN) :: in_char
  CHARACTER(LEN=1)             :: capital
  INTEGER                      :: i
  !
  DO i=1, 26
     IF ( in_char == lower(i:i) ) THEN
        capital = upper(i:i)
        RETURN
     END IF
  END DO
  capital = in_char
  !
END FUNCTION capital
!
!-----------------------------------------------------------------------
FUNCTION lowercase( in_char )  
  !-----------------------------------------------------------------------
  !
  ! ... converts character to lowercase if capital
  ! ... copy character to output in all other cases
  !
  IMPLICIT NONE  
  !
  CHARACTER(LEN=1), INTENT(IN) :: in_char
  CHARACTER(LEN=1)             :: lowercase
  INTEGER                      :: i
  !
  DO i=1, 26
     IF ( in_char == upper(i:i) ) THEN
        lowercase = lower(i:i)
        RETURN
     END IF
  END DO
  lowercase = in_char
  !
END FUNCTION lowercase
!
!-----------------------------------------------------------------------
LOGICAL FUNCTION isnumeric ( in_char )  
  !-----------------------------------------------------------------------
  !
  ! ... check if a character is a number
  !
  IMPLICIT NONE  
  !
  CHARACTER(LEN=1), INTENT(IN) :: in_char
  CHARACTER(LEN=10), PARAMETER :: numbers = '0123456789'
  INTEGER                      :: i
  !
  DO i=1, 10
     isnumeric = ( in_char == numbers(i:i) )
     IF ( isnumeric ) RETURN
  END DO
  !
END FUNCTION isnumeric
!
!-----------------------------------------------------------------------
FUNCTION matches( string1, string2 )  
  !-----------------------------------------------------------------------
  !
  ! ... .TRUE. if string1 is contained in string2, .FALSE. otherwise
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*), INTENT(IN) :: string1, string2
  LOGICAL                       :: matches
  INTEGER                       :: len1, len2, l  
  !
  !
  len1 = LEN_TRIM( string1 )  
  len2 = LEN_TRIM( string2 )  
  !
  DO l = 1, ( len2 - len1 + 1 )  
     IF ( string1(1:len1) == string2(l:(l+len1-1)) ) THEN  
        matches = .TRUE.  
        RETURN  
     END IF
  END DO
  matches = .FALSE.
  !
END FUNCTION matches
  !
  !--------------------------------------------------------------------------
  SUBROUTINE version_parse(str, major, minor, patch, ierr)
    !--------------------------------------------------------------------------
    !   
    ! Determine the major, minor and patch numbers from 
    ! a version string with the fmt "i.j.k"
    ! 
    ! The ierr variable assumes the following values
    !
    ! ierr < 0     emtpy string
    ! ierr = 0     no problem
    ! ierr > 0     fatal error
    !   
    IMPLICIT NONE
    CHARACTER(*),     INTENT(in)    :: str
    INTEGER,          INTENT(out)   :: major, minor, patch, ierr
    !
    INTEGER       :: i1, i2, length
    INTEGER       :: ierrtot
    CHARACTER(10) :: num(3)

    !
    major = 0
    minor = 0
    patch = 0

    length = LEN_TRIM( str )
    !
    IF ( length == 0 ) THEN
       !
       ierr = -1
       RETURN
       !
    ENDIF

    i1 = SCAN( str, ".")
    i2 = SCAN( str, ".", BACK=.TRUE.)
    !
    IF ( i1 == 0 .OR. i2 == 0 .OR. i1 == i2 ) THEN
       !
       ierr = 1
       RETURN
       !
    ENDIF
    !
    num(1) = str(    1 : i1-1 )
    num(2) = str( i1+1 : i2-1 )
    num(3) = str( i2+1 : )
    !
    ierrtot = 0
    !
    READ( num(1), *, IOSTAT=ierr ) major
    IF (ierr/=0) RETURN
    !
    READ( num(2), *, IOSTAT=ierr ) minor
    IF (ierr/=0) RETURN
    !
    READ( num(3), *, IOSTAT=ierr ) patch
    IF (ierr/=0) RETURN
    !
  END SUBROUTINE version_parse
  !
  !--------------------------------------------------------------------------
  FUNCTION version_compare(str1, str2)
    !--------------------------------------------------------------------------
    !   
    ! Compare two version strings; the result is
    !   
    ! "newer":   str1 is newer that str2    
    ! "equal":   str1 is equal   to str2    
    ! "older":   str1 is older than str2    
    ! " ":       str1 or str2 has a wrong format
    !   
    IMPLICIT NONE
    CHARACTER(*)  :: str1, str2
    CHARACTER(10) :: version_compare
    !   
    INTEGER   :: version1(3), version2(3)
    INTEGER   :: basis, icheck1, icheck2
    INTEGER   :: ierr
    !   

    version_compare = " " 
    !   
    CALL version_parse( str1, version1(1), version1(2), version1(3), ierr) 
    IF ( ierr/=0 ) RETURN
    !   
    CALL version_parse( str2, version2(1), version2(2), version2(3), ierr) 
    IF ( ierr/=0 ) RETURN
    !   
    !   
    basis = 1000
    !   
    icheck1 = version1(1) * basis**2 + version1(2)* basis + version1(3) 
    icheck2 = version2(1) * basis**2 + version2(2)* basis + version2(3) 
    !   
    IF ( icheck1 > icheck2 ) THEN
       !   
       version_compare = 'newer'
       !   
    ELSEIF( icheck1 == icheck2 ) THEN
       !   
       version_compare = 'equal'
       !   
    ELSE
       !   
       version_compare = 'older'
       !   
    ENDIF
    !   
  END FUNCTION version_compare

END MODULE upf_utils
