!
SUBROUTINE mem_counter(valore, sumsign, label)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: valore
    INTEGER, INTENT(IN) :: sumsign
    CHARACTER (LEN=*), INTENT(IN) :: label
    !
    INTEGER, SAVE :: current_memory = 0
    INTEGER, SAVE :: last_print_val = 0
    INTEGER, PARAMETER :: MB = 1024*1024
    !
    current_memory = current_memory + valore*sumsign
    IF ( ( current_memory - last_print_val ) > 2*MB) THEN
      WRITE(*, '("Max allocated memory: ", I6, " MB, variable name: ",a)') &
              current_memory/MB, trim(label)
      last_print_val = current_memory
    END IF
    !
END SUBROUTINE
