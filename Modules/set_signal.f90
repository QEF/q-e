MODULE set_signal
! This module is a Fortran 2003 interface to the customize_signals.c C file
! Compatible with Intel/PGI/Gcc(>=4.3) compilers

! This module is compiled only if the following preprocessing option
! is enabled
#if defined(__TRAP_SIGUSR1) || defined(__TERMINATE_GRACEFULLY)

USE iso_c_binding
USE io_global, ONLY : stdout
USE mp_world,  ONLY : root, world_comm, mpime
USE mp, ONLY : mp_bcast

IMPLICIT NONE

LOGICAL,VOLATILE::signal_trapped
INTEGER(kind=c_int),PARAMETER :: SIGINT = 2_c_int

INTERFACE 
   FUNCTION init_signal_USR1(new_handler) BIND(c, name = "init_signal_USR1")
     USE iso_c_binding
     TYPE(C_FUNPTR),VALUE,INTENT(IN):: new_handler
     INTEGER(C_INT)::init_signal_USR1
   END FUNCTION init_signal_USR1

   FUNCTION init_TERMINATE_GRACEFULLY(new_handler) BIND(c, name = "init_TERMINATE_GRACEFULLY")
     USE iso_c_binding
     TYPE(C_FUNPTR),VALUE,INTENT(IN):: new_handler
     INTEGER(C_INT)::init_TERMINATE_GRACEFULLY
   END FUNCTION init_TERMINATE_GRACEFULLY

   FUNCTION init_signal(signum, new_handler) BIND(c, name = "init_signal")
     USE iso_c_binding
     INTEGER(C_INT),VALUE :: signum
     TYPE(C_FUNPTR), VALUE,INTENT(IN) :: new_handler
     INTEGER(C_INT)::init_signal
   END FUNCTION init_signal

END INTERFACE

CONTAINS

#if defined(__TRAP_SIGUSR1)
SUBROUTINE set_signal_USR1(routine)
  USE iso_c_binding
  TYPE(C_FUNPTR),TARGET::ptr
  INTERFACE 
     SUBROUTINE routine(signal) bind(C)
       USE iso_c_binding
       INTEGER(C_INT),VALUE, INTENT(IN)::signal
     END SUBROUTINE routine

  END INTERFACE
       
  ptr = C_FUNLOC(routine)
  
  IF (init_signal_USR1(ptr) .NE. 0) THEN
     CALL errore("set_signal_USR1", "The association of signal USR1 failed!", 1)
  ENDIF
 
END SUBROUTINE set_signal_USR1
#endif

#if defined(__TERMINATE_GRACEFULLY)
SUBROUTINE set_TERMINATE_GRACEFULLY(routine)
  USE iso_c_binding
  TYPE(C_FUNPTR),TARGET::ptr
  INTERFACE
     SUBROUTINE routine(signal) bind(C)
       USE iso_c_binding
       INTEGER(C_INT),VALUE, INTENT(IN)::signal
     END SUBROUTINE routine

  END INTERFACE

  ptr = C_FUNLOC(routine)

  IF (init_TERMINATE_GRACEFULLY(ptr) .NE. 0) THEN
     CALL errore("set_TERMINATE_GRACEFULLY", "The association of signals INT or TERM failed!", 1)
  ENDIF

END SUBROUTINE set_TERMINATE_GRACEFULLY
#endif

! Unused. Here for possible future developments
SUBROUTINE set_signal_action(signal, routine)
  USE iso_c_binding
  INTEGER::signal
  TYPE(C_FUNPTR),TARGET::ptr
  INTERFACE 
     SUBROUTINE routine(signal) bind(C)
       USE iso_c_binding
       INTEGER(C_INT),VALUE::signal
     END SUBROUTINE routine
  END INTERFACE

  ptr = C_FUNLOC(routine)
       
  IF (init_signal(signal, ptr) .NE. 0) THEN
     CALL errore("set_signal", "The association of the signal failed!", 1)
  ENDIF
END SUBROUTINE set_signal_action


! Sets the signal_trapped flag on all nodes/processors
! Only the master will use the signal, though
SUBROUTINE custom_handler(signum) BIND(c)
  USE iso_c_binding
#if defined(__MPI)
  USE mp_world,      ONLY : world_comm
  USE mp,            ONLY : mp_abort
#endif
  INTEGER(C_INT),VALUE,INTENT(IN):: signum
  ! Double CTRL-C will stop immediately;
  ! This cannot be done with any signal because some implementation of MPI 
  ! send SIGTERM to every process when SIGINT (aka CTRL-C) is received
  IF(signal_trapped.and.signum==SIGINT) THEN
    WRITE(stdout, '(/,5x,a)') "**** SIGNAL ALREADY TRAPPED: terminating immediately!!", signum
#if defined(__MPI)
    CALL mp_abort(signum, world_comm)
#else
    STOP 1
#endif
  ELSE
    WRITE(stdout, '(/,5x,a)') "**** Trapped signal: trying to terminate gracefully", signum
    IF(signum==SIGINT) &
      WRITE(stdout, '(5x,a)') "**** press CTRL-C again to terminate immediately (no restart possible!)", signum
    !
    signal_trapped = .TRUE.
  ENDIF
  !
END SUBROUTINE custom_handler


! Set the signal handler for SIGUSR1 to 'custom_handler' 
! Every processor will trap the signal, howver only 0 will actually
! use the result (required since the default action for SIGUSR1 is exit)
SUBROUTINE signal_trap_init
  USE iso_c_binding
#if defined(__TRAP_SIGUSR1)
  WRITE(stdout, FMT='(5x,a)') "signal trapping enabled: kill the code with -SIGUSR1 to stop cleanly the simulation "
  CALL set_signal_USR1(custom_handler)
#endif
#if defined(__TERMINATE_GRACEFULLY)
  WRITE(stdout, FMT='(/,5x,a)') "Signal trapping enabled: code will terminate cleanly with SIGINT, SIGTERM, SIGUSR1, SIGUSR2, SIGXCPU"
  WRITE(stdout, FMT='(5x,a)') "Type CTRL-C twice to terminate immediately (no restart possible!)"
  CALL set_TERMINATE_GRACEFULLY(custom_handler)
#endif
END SUBROUTINE signal_trap_init

FUNCTION signal_detected()
  LOGICAL::signal_detected
  ! If the signal is trapped, set the exit status and broadcast it
  ! DO NOT broadcast the signal_trapped variable or you will be Very
  ! Sorry
  signal_detected = signal_trapped
     
  CALL mp_bcast(signal_detected, root, world_comm)

END FUNCTION signal_detected

#else

USE io_global, ONLY : stdout

CONTAINS

! Place holders to employ when the signal trapping feature is disabled
SUBROUTINE signal_trap_init
  WRITE(stdout, FMT=*) "signal trapping disabled: compile with "
  WRITE(stdout, FMT=*) "-D__TRAP_SIGUSR1 to enable this feature"
END SUBROUTINE signal_trap_init

FUNCTION signal_detected()
  LOGICAL::signal_detected
  signal_detected = .FALSE.
END FUNCTION signal_detected
  
#endif

END MODULE set_signal
