!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
 !-----------------------------------------------------------------------
  FUNCTION trimcheck ( directory )
    !-----------------------------------------------------------------------
    !
    ! ... verify if directory ends with / or \ for Windows, add one if needed
    ! ... trim white spaces and put the result in trimcheck
    !
    IMPLICIT NONE
    !
    CHARACTER (LEN=*), INTENT(IN) :: directory
    CHARACTER (LEN=256) :: trimcheck
#if defined (_WIN32)
#if defined (__PGI)
    CHARACTER (LEN=1) :: separator = '/'
#else
    CHARACTER (LEN=1) :: separator = '\'
#endif
#else
    CHARACTER (LEN=1) :: separator = '/'
#endif
    INTEGER  :: l
    !
    l = LEN_TRIM( directory )
    IF ( l == 0 ) CALL errore( 'trimcheck', ' input name empty', 1)
    !
    IF ( directory(l:l) == separator ) THEN
       trimcheck = TRIM ( ADJUSTL(directory) )
    ELSE
       trimcheck = TRIM ( ADJUSTL(directory) )
       IF ( l < LEN( trimcheck ) ) THEN
          trimcheck = TRIM ( ADJUSTL(directory) ) // separator
       ELSE
          CALL errore(  'trimcheck', ' input name too long', l )
       END IF
    END IF
    !
    RETURN
    !
  END FUNCTION trimcheck
  !
