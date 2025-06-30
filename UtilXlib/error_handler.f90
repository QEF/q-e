!
! Copyright (C) 2002-2025 Quantum ESPRESSO Foundation
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE errore( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output: 
  ! ... if ierr <= 0 it does nothing,
  ! ... if ierr  > 0 it stops.
  !
  ! ...          **** Important note for parallel execution ***
  !
  ! ... in parallel execution unit 6 is written only by the first node;
  ! ... all other nodes have unit 6 redirected to nothing (/dev/null).
  ! ... We write to the "*" unit instead, that appears on all nodes.
  ! ... Effective but annoying!
  !
  USE util_param
#if defined(__PTRACE) && defined(__INTEL_COMPILER)
  USE ifcore,    ONLY : tracebackqq
#endif
  USE mp,        ONLY : mp_abort, mp_rank
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routine
    ! the output message
  INTEGER,          INTENT(IN) :: ierr
    ! the error flag
  INTEGER :: crashunit, mpime
  CHARACTER(LEN=6) :: cerr
  !
  IF( ierr <= 0 ) RETURN
  !
  ! ... the error message is written on the "*" unit
  !
  WRITE( cerr, FMT = '(I6)' ) ierr
  WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
  WRITE( UNIT = *, FMT = '(5X,"Error in routine ",A," (",A,"):")' ) &
        TRIM(calling_routine), TRIM(ADJUSTL(cerr))
  WRITE( UNIT = *, FMT = '(5X,A)' ) TRIM(message)
  WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
  !
  WRITE( *, '("     stopping ...")' )
  !
  FLUSH( stdout )
  !
#if defined(__PTRACE)
#if defined(__INTEL_COMPILER)
    call tracebackqq(user_exit_code=-1)
#elif __GFORTRAN__
    call backtrace
#else
    WRITE( UNIT = 0, FMT = '(5X,A)' ) "Printing strace..."
    CALL ptrace()
#endif
#endif
  !
  !  .. write the message to a file and close it before exiting
  !  .. this will prevent loss of information on systems that
  !  .. do not flush the open streams
  !  .. added by C.C.
  !
  OPEN( NEWUNIT = crashunit, FILE = crash_file, &
        POSITION = 'APPEND', STATUS = 'UNKNOWN' )
  !
  WRITE( UNIT = crashunit, FMT = '(/,1X,78("%"))' )
#if defined(__MPI)
  mpime = mp_rank(MPI_COMM_WORLD)
  WRITE( UNIT = crashunit, FMT = '(5X,"task #",I10)' ) mpime
#endif
  WRITE( UNIT = crashunit, &
         FMT = '(5X,"from ",A," : error #",I10)' ) TRIM(calling_routine), ierr
  WRITE( UNIT = crashunit, FMT = '(5X,A)' ) TRIM(message)
  WRITE( UNIT = crashunit, FMT = '(1X,78("%"),/)' )
  !
  CLOSE( UNIT = crashunit )
  !
#if defined(__MPI)
  ! ... try to exit in a smooth way
  !
  CALL mp_abort(1,MPI_COMM_WORLD)
  !
#endif
  !
  STOP 1
  !
END SUBROUTINE errore
!
!----------------------------------------------------------------------
SUBROUTINE infomsg( routine, message )
  !----------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an info message
  ! ... from a given routine to output.
  !
  USE util_param
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*) :: routine, message
  ! the name of the calling routine
  ! the output message
  !
!  IF ( ionode ) THEN   !if not ionode it is redirected to /dev/null anyway
     !
     WRITE( stdout , '(5X,"Message from routine ",A,":")' ) routine
     WRITE( stdout , '(5X,A)' ) message
     !
!  END IF
  !
  RETURN
  !
END SUBROUTINE infomsg
!
module error_handler
  implicit none
  private

  public :: init_error, add_name, chop_name, error_mem, warning

  type chain
   character (len=35)   :: routine_name
   type(chain), pointer :: previous_link
  end type chain

  type(chain), pointer :: routine_chain

contains

  subroutine init_error(routine_name)
    implicit none
    character (len=*), intent(in) :: routine_name

    allocate(routine_chain)

    routine_chain%routine_name  =  routine_name
    nullify(routine_chain%previous_link)

    return
  end subroutine init_error

  subroutine add_name(routine_name)
    implicit none
    character (len=*), intent(in) :: routine_name
    type(chain), pointer          :: new_link

    allocate(new_link)
    new_link%routine_name  =  routine_name
    new_link%previous_link => routine_chain
    routine_chain          => new_link

    return
  end subroutine add_name

  subroutine chop_name
    implicit none
    type(chain), pointer :: chopped_chain

    chopped_chain => routine_chain%previous_link
    deallocate(routine_chain)
    routine_chain => chopped_chain

    return
  end subroutine chop_name

  recursive subroutine trace_back(error_code)

    implicit none
    integer :: error_code

    write(unit=*,fmt=*) "   Called by ", routine_chain%routine_name
    if (.not.associated(routine_chain%previous_link)) then
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
       write(unit=*,fmt=*) " "
       if( error_code > 0 ) then
          stop
       else
          return
       end if
    end if

    routine_chain => routine_chain%previous_link
    call trace_back(error_code)

  end subroutine trace_back

  subroutine error_mem(message,error_code)
    character (len=*), intent(in) :: message
    integer, intent(in), optional :: error_code
    integer                       :: action_code
    type(chain), pointer          :: save_chain

    if (present(error_code)) then
       action_code = error_code
    else
       action_code = 1
    end if

    if( action_code /= 0 ) then
       write(unit=*,fmt=*) " "
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"

       if( action_code > 0 ) then
          write(unit=*,fmt=*) "   Fatal error in routine `", &
               trim(routine_chain%routine_name),"': ",message
       else
          write(unit=*,fmt=*) "   Warning from routine `", &
               trim(routine_chain%routine_name),"': ",message
          save_chain => routine_chain
       end if
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
       routine_chain => routine_chain%previous_link
       call trace_back(action_code)
       routine_chain => save_chain
    end if

    return
  end subroutine error_mem

  subroutine warning(message)
    character (len=*), intent(in) :: message
    call error_mem(message,-1)
    return
  end subroutine warning

end module error_handler
