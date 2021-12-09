!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION capital( in_char )  
  !-----------------------------------------------------------------------
  !! Converts character to capital if lowercase.
  !! Copy character to output in all other cases.
  !
  IMPLICIT NONE  
  !
  CHARACTER(LEN=1), INTENT(IN) :: in_char
  CHARACTER(LEN=1)             :: capital
  CHARACTER(LEN=26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz', &
                                  upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  INTEGER                      :: i
  !
  !
  DO i=1, 26
     !
     IF ( in_char == lower(i:i) ) THEN
        !
        capital = upper(i:i)
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  capital = in_char
  !
  RETURN 
  !
END FUNCTION capital
!
!-----------------------------------------------------------------------
FUNCTION lowercase( in_char )  
  !-----------------------------------------------------------------------
  !! Converts character to lowercase if capital.
  !! Copy character to output in all other cases.
  !
  IMPLICIT NONE  
  !
  CHARACTER(LEN=1), INTENT(IN) :: in_char
  CHARACTER(LEN=1)             :: lowercase
  CHARACTER(LEN=26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz', &
                                  upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  INTEGER                      :: i
  !
  !
  DO i=1, 26
     !
     IF ( in_char == upper(i:i) ) THEN
        !
        lowercase = lower(i:i)
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  lowercase = in_char
  !
  RETURN 
  !
END FUNCTION lowercase
