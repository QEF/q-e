  !-----------------------------------------------------------------------
  FUNCTION int_to_char( i )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: i
    CHARACTER (LEN=6)   :: int_to_char
    !
    !   
    IF ( i < 10 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I1)" ) i
       !
    ELSE IF ( i < 100 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I2)" ) i
       !
    ELSE IF ( i < 1000 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I3)" ) i
       !
    ELSE IF ( i < 10000 ) THEN
       !
       WRITE( UNIT = int_to_char , FMT = "(I4)" ) i
       !
    ELSE      
       ! 
       WRITE( UNIT = int_to_char , FMT = "(I5)" ) i     
       !
    END IF    
    !
    RETURN
    !
  END FUNCTION int_to_char
