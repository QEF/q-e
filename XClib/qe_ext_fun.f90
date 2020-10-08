!
! Copyright (C) 2001-2004 Carlo Cavazzoni and PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
     !   
     IF ( string1(1:len1) == string2(l:(l+len1-1)) ) THEN  
        !
        matches = .TRUE.  
        !
        RETURN  
        !
     END IF
     !
  END DO
  !
  matches = .FALSE.
  ! 
  RETURN
  !
END FUNCTION matches
!

! !-----------------------------------------------------------------------
! FUNCTION capital( in_char )  
!   !-----------------------------------------------------------------------
!   !
!   ! ... converts character to capital if lowercase
!   ! ... copy character to output in all other cases
!   !
!   IMPLICIT NONE  
!   !
!   CHARACTER(LEN=1), INTENT(IN) :: in_char
!   CHARACTER(LEN=1)             :: capital
!   CHARACTER(LEN=26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz', &
!                                   upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!   INTEGER                      :: i
!   !
!   !
!   DO i=1, 26
!      !
!      IF ( in_char == lower(i:i) ) THEN
!         !
!         capital = upper(i:i)
!         !
!         RETURN
!         !
!      END IF
!      !
!   END DO
!   !
!   capital = in_char
!   !
!   RETURN 
!   !
! END FUNCTION capital
!
