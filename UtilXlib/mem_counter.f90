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
#if defined (__DEBUG)
    if ( sumsign == 1) then
       WRITE(*, '("  allocating ", I6, " kB, variable: ",a)') &
       valore/1024, trim(label)
    else
       WRITE(*, '("deallocating ", I6, " kB, variable: ",a)') &
       valore/1024, trim(label)
    end if
#endif
    current_memory = current_memory + valore*sumsign
    IF ( ( current_memory - last_print_val ) > 2*MB) THEN
      WRITE(*, '("Max allocated memory: ", I6, " MB, variable: ",a)') &
              current_memory/MB, trim(label)
      last_print_val = current_memory
    END IF
    !
END SUBROUTINE
