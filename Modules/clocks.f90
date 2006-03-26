!
! Copyright (C) 2001-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mytime
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER,        PARAMETER :: maxclock = 100
  REAL (DP), PARAMETER :: notrunning = - 1.D0
  ! 
  REAL (DP)     :: myclock(maxclock), t0(maxclock)
  CHARACTER (LEN=12) :: clock_label(maxclock)
  INTEGER            :: called(maxclock)
  !
  INTEGER :: nclock = 0
  LOGICAL :: no
  !
END MODULE mytime
!
!----------------------------------------------------------------------------
SUBROUTINE init_clocks( go )
  !----------------------------------------------------------------------------
  !
  ! flag = .TRUE.  : clocks will run
  ! flag = .FALSE. : only clock #1 will run
  !
  USE kinds,  ONLY : DP
  USE mytime, ONLY : called, t0, myclock, no, notrunning, maxclock
  !
  IMPLICIT NONE
  !
  LOGICAL :: go
  INTEGER :: n
  !
  !
  no = .NOT. go
  !
  DO n = 1, maxclock
     !
     called(n)  = 0
     myclock(n) = 0.D0
     t0(n)      = notrunning
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE init_clocks
!
!----------------------------------------------------------------------------
SUBROUTINE start_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : nclock, clock_label, notrunning, no, maxclock, t0
  !
  IMPLICIT NONE
  !
  REAL (DP)    :: scnds
  CHARACTER (LEN=*) :: label
  INTEGER           :: n
  !
  !
  IF ( no .AND. ( nclock == 1 ) ) RETURN
  !
  DO n = 1, nclock
     !
     IF ( label == clock_label(n) ) THEN
        !
        ! ... found previously defined clock: check if not already started,
        ! ... store in t0 the starting time
        !
        IF ( t0(n) /= notrunning ) THEN
           WRITE( stdout, '("start_clock: clock # ",I2," for ",A12, &
                          & " already started")' ) n, label
        ELSE
           t0(n) = scnds()
        END IF
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  ! ... clock not found : add new clock for given label
  !
  IF ( nclock == maxclock ) THEN
     !
     WRITE( stdout, '("start_clock: Too many clocks! call ignored")' )
     !
  ELSE
     !
     nclock              = nclock + 1
     clock_label(nclock) = label
     t0(nclock)          = scnds()
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE start_clock
!
!----------------------------------------------------------------------------
SUBROUTINE stop_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : no, nclock, clock_label, myclock, &
                        notrunning, called, t0
  !
  IMPLICIT NONE
  !
  REAL (DP)    :: scnds
  CHARACTER (LEN=*) :: label
  INTEGER           :: n
  !
  !
  IF ( no ) RETURN
  !
  DO n = 1, nclock
     !
     IF ( label == clock_label(n) ) THEN
        !
        ! ... found previously defined clock : check if properly initialised,
        ! ... add elapsed time, increase the counter of calls
        !
        IF ( t0(n) == notrunning ) THEN
           !
           WRITE( stdout, '("stop_clock: clock # ",I2," for ",A12, &
                          & " not running")' ) n, label
           !
        ELSE
           !
           myclock(n) = myclock(n) + scnds() - t0(n)
           t0(n)      = notrunning
           called(n)  = called(n) + 1
           !
        END IF
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  ! ... clock not found
  !
  WRITE( stdout, '("stop_clock: no clock for ",A12," found !")' ) label
  !
  RETURN
  !
END SUBROUTINE stop_clock
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : nclock, clock_label
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*) :: label
  INTEGER           :: n
  !
  !
  IF ( label == ' ' ) THEN
     !
     WRITE( stdout, * )
     !
     DO n = 1, nclock
        !
        CALL print_this_clock( n )
        !
     END DO
     !
  ELSE
     !
     DO n = 1, nclock
        !
        IF ( label == clock_label(n) ) THEN
           !
           CALL print_this_clock( n )
           !
           RETURN
           !
        END IF
        !
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE print_clock
!
!----------------------------------------------------------------------------
SUBROUTINE print_this_clock( n )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : no, nclock, clock_label, myclock, &
                        notrunning, called, t0
  USE mp,        ONLY : mp_max, mp_min
  USE mp_global, ONLY : intra_image_comm, my_image_id
  !
  IMPLICIT NONE
  !
  REAL(DP) :: scnds
  INTEGER       :: n
  REAL(DP) :: elapsed_cpu_time, nsec
  INTEGER       :: nday, nhour, nmin
  !
  !
  IF ( t0(n) == notrunning ) THEN
     !
     ! ... clock stopped, print the stored value for the cpu time
     !
     elapsed_cpu_time = myclock(n)
     !
  ELSE
     !
     ! ... clock not stopped, print the current value of the cpu time
     !
     elapsed_cpu_time = myclock(n) + scnds() - t0(n)
     !
  END If
  !
#if defined (__PARA)
  !
  ! ... In the parallel case it is far from clear which value to print
  ! ... The following is the maximum over all nodes and pools. NOTA BENE:
  ! ... some trouble could arise if a clock is not started on all nodes
  !
  ! ... by uncommenting the following line the extreme operation is removed
  ! ... may be useful for testing purposes :
  !
!#define DEBUG
  !
#  ifndef DEBUG
  !
  CALL mp_max( elapsed_cpu_time, intra_image_comm )
  !
#  endif
#endif
  !
  IF ( n == 1 ) THEN
     !
     ! ... The first clock is written as days/hour/min/sec
     !
     nday  = elapsed_cpu_time / 86400
     nsec  = elapsed_cpu_time - 86400 * nday
     nhour = nsec / 3600 
     nsec  = nsec - 3600 * nhour
     nmin  = nsec / 60
     nsec  = nsec - 60 * nmin
     !
     IF ( nday > 0 ) THEN
        !    
        WRITE( stdout, &
               '(5X,A12," : ",3X,I2,"d",3X,I2,"h",I2, "m CPU time"/)' ) &
             clock_label(n), nday, nhour, nmin
        !
     ELSE IF ( nhour > 0 ) THEN
        !
        WRITE( stdout, &
               '(5X,A12," : ",3X,I2,"h",I2,"m CPU time"/)' ) &
             clock_label(n), nhour, nmin
        !
     ELSE IF ( nmin > 0 ) THEN
        !
        WRITE( stdout, &
               '(5X,A12," : ",I2,"m",F5.2,"s CPU time"/)' ) &
             clock_label(n), nmin, nsec
        !
     ELSE
        !
        WRITE( stdout, &
               '(5X,A12," : ",3X,F5.2,"s CPU time"/)' ) &
             clock_label(n), nsec
        !
     END IF
     !
  ELSE IF ( called(n) == 1 .OR. t0(n) /= notrunning ) THEN
     !
     ! ... for clocks that have been called only once
     !
     WRITE( stdout, &
            '(5X,A12," :",F9.2,"s CPU")') clock_label(n), elapsed_cpu_time
     !
  ELSE IF ( called(n) == 0 ) THEN
     !
     ! ... for clocks that have never been called
     !
     WRITE( stdout, &
            '("print_this: clock # ",I2," for ",A12," never called !")' ) &
         n, clock_label(n)
     !
  ELSE
     !
     ! ... for all other clocks
     !
     WRITE( stdout, &
            '(5X,A12," :",F9.2,"s CPU (",I8," calls,",F8.3," s avg)")' ) &
         clock_label(n), elapsed_cpu_time, &
         called(n), ( elapsed_cpu_time / called(n) )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE print_this_clock
!
!----------------------------------------------------------------------------
FUNCTION get_clock( label )
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mytime,    ONLY : no, nclock, clock_label, myclock, &
                        notrunning, called, t0
  USE mp,        ONLY : mp_max, mp_min
  USE mp_global, ONLY : intra_image_comm 
  !
  IMPLICIT NONE
  !
  REAL(DP)     :: get_clock
  REAL(DP)     :: scnds
  CHARACTER (LEN=*) :: label
  INTEGER           :: n
  !
  !
  IF ( no ) THEN
     !
     IF ( label == clock_label(1) ) THEN
        !
        get_clock = scnds()
        !
     ELSE
        !
        get_clock = notrunning
        !
     END IF
     !
     RETURN
     !
  END IF
  !
  DO n = 1, nclock
     !
     IF ( label == clock_label(n) ) THEN
        !
        IF ( t0(n) == notrunning ) THEN
           !
           get_clock = myclock(n)
           !
        ELSE
           !
           get_clock = myclock(n) + scnds() - t0(n)
           !
        END IF
        !
        ! ... In the parallel case, use the maximum over all nodes and pools
        !
        CALL mp_max( get_clock, intra_image_comm )
        !
        RETURN
        !
     END IF
     !
  END DO
  !
  ! ... clock not found
  !
  get_clock = notrunning
  !
  WRITE( stdout, '("get_clock: no clock for ",A12," found !")') label
  !
  RETURN
  !
END FUNCTION get_clock
