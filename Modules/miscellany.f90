!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
MODULE miscellany
  !---------------------------------------------------------------------------
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  CONTAINS
     ! 
     PURE FUNCTION norm( vect )
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), DIMENSION(:), INTENT(IN)  :: vect
       REAL (KIND=DP)                            :: norm
       !
       !
       norm = SQRT( DOT_PRODUCT( vect , vect ) )
       !
       RETURN
       !
     END FUNCTION norm
     !
     !
     PURE FUNCTION int_to_char( int )
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
END MODULE miscellany
