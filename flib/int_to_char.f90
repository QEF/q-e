  !-----------------------------------------------------------------------
  FUNCTION int_to_char( int )
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
