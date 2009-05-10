
  !-----------------------------------------------------------------------
  FUNCTION int_to_char( i )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: i
    CHARACTER (LEN=6)   :: int_to_char
    CHARACTER :: c
    INTEGER   :: n, j
    !   
    IF( i < 0 ) &
       CALL errore( ' int_to_char ', ' negative input not allowed ', 1 )
    !   
    n = i
    j = 1
    DO WHILE( j < 7 ) 
       int_to_char(j:j) = CHAR( MOD( n, 10 ) + ICHAR( '0' ) )
       n = n / 10
       IF( n == 0 ) EXIT
       j = j + 1
    END DO
    !
    IF( j < 7 ) THEN
       DO n = 1, j/2
          c = int_to_char( n : n )
          int_to_char( n : n ) = int_to_char( j-n+1 : j-n+1 )
          int_to_char( j-n+1 : j-n+1 ) = c
       END DO
       IF( j < 6 ) int_to_char(j+1:6) = ' '
    ELSE
       int_to_char(1:6) = '*'
    END IF
    !
    RETURN
    !
  END FUNCTION int_to_char
