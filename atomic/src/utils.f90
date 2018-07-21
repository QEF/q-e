module utils
implicit none

private 
public assert

contains

    subroutine assert(condition)
    ! If condition == .false., it aborts the program.
    !
    ! Arguments
    ! ---------
    !
    logical, intent(in) :: condition
    !
    ! Example
    ! -------
    !
    ! call assert(a == 5)

    if (.not. condition) call stop_error("Assert failed.")
    end subroutine
end module utils
