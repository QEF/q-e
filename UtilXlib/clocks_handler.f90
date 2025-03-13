!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... Time-printing utilities - Contains the following subroutines:
!     init_clocks( go )    initialization - must be called first
!                          go = .TRUE. : up to "maxclock" clocks can be started
!                          go = .FALSE.: only clock #1 can be started
!     start_clock( label )   starts clock "label" (max 12 characters)
!                            if "label" has never been started, initializes it
!                            issues warning if "label" already started
!     stop_clock( label )    stops  clock "label"
!                            issues warning if "label" is either not running
!                            or has never been started
!     print_clock( label )   print cpu and wall time measured by clock "label"
!                            clock "label" may be running or stopped 
!                            and remains in the same state
!                            issues warning if "label" has never been started
! ... and the following function (real(kind=dp):
!     get_clock( label )     return wall time measured by clock "label"
!                            returns -1 if "label" has never been started
! ... All output and warnings are written to stdout
! ... Clocks should be started, read, stopped either on all processors, or 
! ... only on one, but not half and half! For parallel debugging, uncomment:
!#define __TRACE
! ... See also comments in subroutine print_this_clock about parallel case
!
!----------------------------------------------------------------------------
MODULE mytime
  !----------------------------------------------------------------------------
  !
  USE util_param, ONLY : DP
  USE parallel_include
#if defined(__CUDA)
  USE cudafor
#endif
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER,  PARAMETER :: maxclock = 128
  REAL(DP), PARAMETER :: notrunning = - 1.0_DP
  !
  REAL(DP)          :: cputime(maxclock), t0cpu(maxclock)
  REAL(DP)          :: walltime(maxclock), t0wall(maxclock)
  REAL(DP)          :: gputime(maxclock)
  CHARACTER(len=12) :: clock_label(maxclock)
  INTEGER           :: called(maxclock)
  INTEGER           :: gpu_called(maxclock)
  INTEGER           :: max_print_depth = maxclock  ! used to gauge the amount of output. default: a very deep depth
  !
  REAL(DP)          :: mpi_per_thread = 1.0_DP

  INTEGER :: nclock = 0
  LOGICAL :: no
#if defined (__TRACE)
  INTEGER :: trace_depth = 0
  INTEGER :: mpime
#endif
#if defined(__CUDA)
  type(cudaEvent) :: gpu_starts(maxclock), gpu_stops(maxclock)
#else
  ! dummy values never used!
  INTEGER :: gpu_starts(maxclock), gpu_stops(maxclock)
#endif
  INTERFACE
     FUNCTION f_wall ( ) BIND(C,name="cclock") RESULT(t)
       USE ISO_C_BINDING
       REAL(kind=c_double) :: t
     END FUNCTION f_wall
     FUNCTION f_tcpu ( ) BIND(C,name="scnds") RESULT(t)
       USE ISO_C_BINDING
       REAL(kind=c_double) :: t
     END FUNCTION f_tcpu
  END INTERFACE
  !
END MODULE mytime
!
!----------------------------------------------------------------------------
SUBROUTINE init_clocks( go )
  !----------------------------------------------------------------------------
  !
  ! ... go = .TRUE.  : clocks will run
  ! ... go = .FALSE. : only clock #1 will run
  !
  USE util_param,  ONLY : DP, stdout
  USE mytime, ONLY : called, t0cpu, cputime, no, notrunning, maxclock, &
       clock_label, walltime, t0wall, nclock, mpi_per_thread
  ! ... GPU related timers
  USE mytime, ONLY : gpu_starts, gpu_stops, gpu_called, gputime
#if defined (__TRACE)
  USE mytime, ONLY : mpime, max_print_depth, MPI_COMM_WORLD
#endif
#if defined(__CUDA)
  USE cudafor
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: go
  INTEGER :: n, ierr
  !
#if defined(_OPENMP)
  INTEGER, EXTERNAL :: omp_get_max_threads
  mpi_per_thread = 1.0_DP/omp_get_max_threads()
#endif
  no = .not. go
  nclock = 0
  !
  DO n = 1, maxclock
     !
     called(n)      = 0
     gpu_called(n)  = 0
     cputime(n)     = 0.0_DP
     t0cpu(n)       = notrunning
     walltime(n)    = 0.0_DP
     t0wall(n)      = notrunning
     gputime(n)     = 0.0_DP
     clock_label(n) = ' '
#if defined(__CUDA)
     ierr = cudaEventCreate(gpu_starts(n))
     ierr = cudaEventCreate(gpu_stops(n))
#endif
     !
  ENDDO
#if defined (__TRACE)
  write(stdout,*) '*** Code flow traced exploiting clocks calls ***'
  write(stdout,*) '--- Code flow traced down to depth ',max_print_depth
  mpime = 0
#if defined(__MPI)
  ierr = 0
  CALL mpi_comm_rank(MPI_COMM_WORLD,mpime,ierr)
  IF (ierr/=0) then
     WRITE( stdout, fmt='( "*** error in init_clocks call to mpi_comm_rank ***")' )
     WRITE( stdout, fmt='( "*** error code: ",I5)' ) ierr
     ! abort with extreme prejudice across the entire MPI set of tasks
     CALL mpi_abort(MPI_COMM_WORLD,ierr, ierr)
  END IF
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE init_clocks
!
!----------------------------------------------------------------------------
SUBROUTINE set_trace_max_depth(max_depth_)
  USE mytime, ONLY: max_print_depth
  IMPLICIT NONE
  INTEGER  :: max_depth_
  ! 
  max_print_depth = max_depth_   
END SUBROUTINE set_trace_max_depth 
! 
!----------------------------------------------------------------------------
SUBROUTINE start_clock( label )
  !----------------------------------------------------------------------------
  !
  USE util_param,     ONLY : DP, stdout
#if defined (__TRACE)
  USE mytime,    ONLY : trace_depth, mpime, max_print_depth
#endif
  USE mytime,    ONLY : nclock, clock_label, notrunning, no, maxclock, &
                        t0cpu, t0wall, f_wall, f_tcpu
  USE nvtx
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: label
  !
  CHARACTER(len=12):: label_
  INTEGER          :: n
  !
#if defined (__TRACE)
  if (trace_depth <= max_print_depth ) &  ! used to gauge the ammount of output
  WRITE( stdout,'(I3," depth=",I2," start_clock ",A )') mpime,trace_depth, label ; FLUSH(stdout)
  !WRITE( stdout, '("mpime = ",I2,", TRACE (depth=",I2,") Start: ",A12)') mpime, trace_depth, label
  trace_depth = trace_depth + 1
#endif
  !
  IF ( no .and. ( nclock == 1 ) ) RETURN
  !
  ! ... prevent trouble if label is longer than 12 characters
  !
  label_ = trim ( label )
  !
  DO n = 1, nclock
     !
     IF ( clock_label(n) == label_ ) THEN
        !
        ! ... found previously defined clock: check if not already started,
        ! ... store in t0cpu the starting time
        !
        IF ( t0cpu(n) /= notrunning ) THEN
!            WRITE( stdout, '("start_clock: clock # ",I2," for ",A12, &
!                           & " already started")' ) n, label_
        ELSE
           t0cpu(n) = f_tcpu()
           t0wall(n)= f_wall()

           call nvtxStartRange(label_, n)
        ENDIF
        !
        RETURN
        !
     ENDIF
     !
  ENDDO
  !
  ! ... clock not found : add new clock for given label
  !
  IF ( nclock == maxclock ) THEN
     !
     WRITE( stdout, '("start_clock(",A,"): Too many clocks! call ignored")' ) label
     !
  ELSE
     !
     nclock              = nclock + 1
     clock_label(nclock) = label_
     t0cpu(nclock)       = f_tcpu()
     t0wall(nclock)      = f_wall()
     call nvtxStartRange(label_, n)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE start_clock
!
!----------------------------------------------------------------------------
SUBROUTINE start_clock_gpu( label )
  !----------------------------------------------------------------------------
  !
  USE util_param,     ONLY : DP, stdout
#if defined (__TRACE)
  USE mytime,    ONLY : trace_depth, mpime, max_print_depth
#endif
  USE mytime,    ONLY : nclock, clock_label, notrunning, no, maxclock, &
                        t0cpu, t0wall, f_wall, f_tcpu, &
                        gputime, gpu_starts, gpu_stops
  USE nvtx,      ONLY : nvtxStartRange
#if defined(__CUDA)
  USE cudafor
#endif
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: label
  !
  CHARACTER(len=12):: label_
  INTEGER          :: n, ierr
  !
#if defined (__TRACE)
  if (trace_depth <= max_print_depth ) &  ! used to gauge the ammount of output
  WRITE( stdout,'(I3," depth=",I2," start_clock ",A )') mpime,trace_depth, label ; FLUSH(stdout)
  !WRITE( stdout, '("mpime = ",I2,", TRACE (depth=",I2,") Start: ",A12)') mpime, trace_depth, label
  trace_depth = trace_depth + 1
#endif
  !
  IF ( no .and. ( nclock == 1 ) ) RETURN
  !
  ! ... prevent trouble if label is longer than 12 characters
  !
  label_ = trim ( label )
  !
  DO n = 1, nclock
     !
     IF ( clock_label(n) == label_ ) THEN
        !
        ! ... found previously defined clock: check if not already started,
        ! ... store in t0cpu the starting time
        !
        IF ( t0cpu(n) /= notrunning ) THEN
!            WRITE( stdout, '("start_clock: clock # ",I2," for ",A12, &
!                           & " already started")' ) n, label_
        ELSE
#if defined(__CUDA)
           ierr = cudaEventRecord(gpu_starts(n),0)
#endif
           t0cpu(n) = f_tcpu()
           t0wall(n)= f_wall()
           call nvtxStartRange(label_, n)
        ENDIF
        !
        RETURN
        !
     ENDIF
     !
  ENDDO
  !
  ! ... clock not found : add new clock for given label
  !
  IF ( nclock == maxclock ) THEN
     !
     WRITE( stdout, '("start_clock(",A,"): Too many clocks! call ignored")' ) label
     !
  ELSE
     !
     nclock              = nclock + 1
     clock_label(nclock) = label_
#if defined(__CUDA)
     ierr = cudaEventRecord(gpu_starts(nclock),0)
#endif
     t0cpu(nclock)       = f_tcpu()
     t0wall(nclock)      = f_wall()
     call nvtxStartRange(label_, n)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE start_clock_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE stop_clock( label )
  !----------------------------------------------------------------------------
  !
  USE util_param,     ONLY : DP, stdout
#if defined (__TRACE)
  USE mytime,    ONLY : trace_depth, mpime, max_print_depth
#endif
  USE mytime,    ONLY : no, nclock, clock_label, cputime, walltime, &
                        notrunning, called, t0cpu, t0wall, f_wall, f_tcpu
  USE nvtx
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: label
  !
  CHARACTER(len=12):: label_
  INTEGER          :: n
  !
#if defined (__TRACE)
  trace_depth = trace_depth - 1
  if (trace_depth <= max_print_depth ) &  ! used to gauge the ammount of output
  WRITE( stdout,'(I3," depth=",I2," stop_clock ",A )') mpime, trace_depth, label ; FLUSH(stdout)
  !WRITE( *, '("mpime = ",I2,", TRACE (depth=",I2,") End: ",A12)') mpime, trace_depth, label
#endif
  !
  IF ( no ) RETURN
  !
  ! ... prevent trouble if label is longer than 12 characters
  !
  label_ = trim ( label )
  !
  DO n = 1, nclock
     !
     IF ( clock_label(n) == label_ ) THEN
        !
        ! ... found previously defined clock : check if properly initialised,
        ! ... add elapsed time, increase the counter of calls
        !
        IF ( t0cpu(n) == notrunning ) THEN
           !
           WRITE( stdout, '("stop_clock: clock # ",I2," for ",A12, " not running")' ) n, label
           !
        ELSE
           !
           cputime(n)   = cputime(n) + f_tcpu() - t0cpu(n)
           walltime(n)  = walltime(n)+ f_wall() - t0wall(n)
           t0cpu(n)     = notrunning
           t0wall(n)    = notrunning
           called(n)    = called(n) + 1

           call nvtxEndRange
           !
        ENDIF
        !
        RETURN
        !
     ENDIF
     !
  ENDDO
  !
  ! ... clock not found
  !
  WRITE( stdout, '("stop_clock: no clock for ",A12," found !")' ) label
  !
  RETURN
  !
END SUBROUTINE stop_clock
!
SUBROUTINE stop_clock_gpu( label )
  !----------------------------------------------------------------------------
  !
  USE util_param,     ONLY : DP, stdout
#if defined (__TRACE)
  USE mytime,    ONLY : trace_depth, mpime, max_print_depth
#endif
  USE mytime,    ONLY : no, nclock, clock_label, cputime, walltime, &
                        notrunning, called, t0cpu, t0wall, f_wall, f_tcpu, &
                        gpu_called, gputime, gpu_starts, gpu_stops
  USE nvtx,      ONLY:  nvtxEndRange
#if defined(__CUDA)
  USE cudafor
#endif
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: label
  !
  CHARACTER(len=12):: label_
  INTEGER          :: n, ierr
  REAL             :: time
  !
#if defined (__TRACE)
  trace_depth = trace_depth - 1
  if (trace_depth <= max_print_depth ) &  ! used to gauge the ammount of output
  WRITE( stdout,'(I3," depth=",I2," stop_clock ",A )') mpime, trace_depth, label ; FLUSH(stdout)
  !WRITE( *, '("mpime = ",I2,", TRACE (depth=",I2,") End: ",A12)') mpime, trace_depth, label
#endif
  !
  IF ( no ) RETURN
  !
  ! ... initialize time used in CUDA APIs if __CUDA is present.
  !
  time = 0.0
  !
  ! ... prevent trouble if label is longer than 12 characters
  !
  label_ = trim ( label )
  !
  DO n = 1, nclock
     !
     IF ( clock_label(n) == label_ ) THEN
        !
        ! ... found previously defined clock : check if properly initialised,
        ! ... add elapsed time, increase the counter of calls
        !
        IF ( t0cpu(n) == notrunning ) THEN
           !
           WRITE( stdout, '("stop_clock: clock # ",I2," for ",A12, " not running")' ) n, label
           !
        ELSE
           !
           cputime(n)   = cputime(n) + f_tcpu() - t0cpu(n)
#if defined(__CUDA)
           ierr         = cudaEventRecord(gpu_stops(n),0)
           ierr         = cudaEventSynchronize(gpu_stops(n))
           ierr         = cudaEventElapsedTime(time, gpu_starts(n), gpu_stops(n))
#endif           
           gputime(n)   = gputime(n) + time
           gpu_called(n)= gpu_called(n) + 1
           !
           walltime(n)  = walltime(n)+ f_wall() - t0wall(n)
           t0cpu(n)     = notrunning
           t0wall(n)    = notrunning
           called(n)    = called(n) + 1
           call nvtxEndRange()
           !
        ENDIF
        !
        RETURN
        !
     ENDIF
     !
  ENDDO
  !
  ! ... clock not found
  !
  WRITE( stdout, '("stop_clock_gpu: no clock for ",A12," found !")' ) label
  !
  RETURN
  !
END SUBROUTINE stop_clock_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock( label )
  !----------------------------------------------------------------------------
  !
  USE util_param, ONLY : stdout
  USE mytime,     ONLY : nclock, clock_label, gpu_called
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*) :: label
  !
  CHARACTER(len=12) :: label_
  INTEGER          :: n
  LOGICAL          :: print_gpu
  !
  print_gpu = ANY(gpu_called > 0)
  !
  IF ( label == ' ' ) THEN
     !
     WRITE( stdout, * )
     !
     DO n = 1, nclock
        !
        CALL print_this_clock( n )
        IF(print_gpu) CALL print_this_clock_gpu( n )
        !
     ENDDO
     !
  ELSE
     !
     ! ... prevent trouble if label is longer than 12 characters
     !
     label_ = trim ( label )
     !
     DO n = 1, nclock
        !
        IF ( clock_label(n) == label_ ) THEN
           !
           CALL print_this_clock( n )
           IF(print_gpu) CALL print_this_clock_gpu( n )
           !
           exit
           !
        ENDIF
        !
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE print_clock
!
!----------------------------------------------------------------------------
FUNCTION get_cpu_and_wall( n) result (t) 
  !--------------------------------------------------------------------------
  !
  USE util_param, ONLY : DP 
  USE mytime,     ONLY : clock_label, cputime, walltime, mpi_per_thread, &
                         notrunning, called, t0cpu, t0wall, f_wall, f_tcpu
  IMPLICIT NONE 
  ! 
  INTEGER  :: n 
  REAL(DP) :: t(2)  
  !
  IF (t0cpu(n) == notrunning ) THEN 
     t(1) = cputime(n)
     t(2)  = walltime(n)
   ELSE 
     t(1)   = cputime(n) + f_tcpu() - t0cpu(n)
     t(2)   = walltime(n)+ f_wall() - t0wall(n)
   END IF 
#if defined(PRINT_AVG_CPU_TIME_PER_THREAD)
  ! rescale the elapsed cpu time on a per-thread basis
  t(1)    = t(1)  * mpi_per_thread
#endif
END FUNCTION get_cpu_and_wall      
!----------------------------------------------------------------------------
SUBROUTINE print_this_clock( n )
  !----------------------------------------------------------------------------
  !
  USE util_param, ONLY : DP, stdout
  USE mytime,     ONLY : clock_label, cputime, walltime, mpi_per_thread, &
                         notrunning, called, t0cpu, t0wall, f_wall, f_tcpu
  !
  IMPLICIT NONE
  !
  INTEGER  :: n
  REAL(DP) :: elapsed_cpu_time, elapsed_wall_time, nsec, msec
  INTEGER  :: nday, nhour, nmin, nmax, mday, mhour, mmin
  !
  !
  IF ( t0cpu(n) == notrunning ) THEN
     !
     ! ... clock stopped, print the stored value for the cpu time
     !
     elapsed_cpu_time = cputime(n)
     elapsed_wall_time= walltime(n)
     !
  ELSE
     !
     ! ... clock not stopped, print the current value of the cpu time
     !
     elapsed_cpu_time   = cputime(n) + f_tcpu() - t0cpu(n)
     elapsed_wall_time  = walltime(n)+ f_wall() - t0wall(n)
     called(n)  = called(n) + 1
     !
  ENDIF
  !
!#define PRINT_AVG_CPU_TIME_PER_THREAD
#if defined(PRINT_AVG_CPU_TIME_PER_THREAD)
  ! rescale the elapsed cpu time on a per-thread basis
  elapsed_cpu_time   = elapsed_cpu_time * mpi_per_thread
#endif
  !
  nmax = called(n)
  !
  ! ... In the parallel case there are several possible approaches
  ! ... The safest one is to leave each clock independent from the others
  ! ... Another possibility is to print the maximum across all processors
  ! ... This is done by uncommenting the following lines
  !
  ! CALL mp_max( elapsed_cpu_time, intra_image_comm )
  ! CALL mp_max( elapsed_wall_time, intra_image_comm )
  ! CALL mp_max( nmax, intra_image_comm )
  !
  ! ... In the last line we assume that the maximum cpu time
  ! ... is associated to the maximum number of calls
  ! ... NOTA BENE: by uncommenting the above lines you may run into
  ! ... serious trouble if clocks are not started on all nodes
  !
  IF ( n == 1 ) THEN
     !
     ! ... The first clock is written as days/hour/min/sec
     !
#if defined(__CLOCK_SECONDS)
     !
     WRITE( stdout, &
        '(5X,A12," : ",F9.2,"s CPU ",F9.2,"s WALL"/)' ) &
        clock_label(n), elapsed_cpu_time, elapsed_wall_time
     !
#else
     !
     nday  = elapsed_cpu_time / 86400
     nsec  = elapsed_cpu_time - 86400 * nday
     nhour = nsec / 3600
     nsec  = nsec - 3600 * nhour
     nmin  = nsec / 60
     nsec  = nsec - 60 * nmin
     !
     ! ... The first clock writes elapsed (wall) time as well
     !
     mday  = elapsed_wall_time / 86400
     msec  = elapsed_wall_time - 86400 * mday
     mhour = msec / 3600
     msec  = msec - 3600 * mhour
     mmin  = msec / 60
     msec  = msec - 60 * mmin
     !
     IF ( nday > 0 ) THEN
        !
        WRITE( stdout, ADVANCE='no', &
               FMT='(5X,A12," : ",1X,I2,"d",I2,"h",I2,"m CPU ")' ) &
             clock_label(n), nday, nhour, nmin
        !
     ELSEIF ( nhour > 0 ) THEN
        !
        WRITE( stdout, ADVANCE='no', &
                FMT='(5X,A12," : ",4X,I2,"h",I2,"m CPU ")' ) &
             clock_label(n), nhour, nmin
        !
     ELSEIF ( nmin > 0 ) THEN
        !
        WRITE( stdout, ADVANCE='no', &
               FMT='(5X,A12," : ",1X,I2,"m",F5.2,"s CPU ")' ) &
             clock_label(n), nmin, nsec
        !
     ELSE
        !
        WRITE( stdout, ADVANCE='no', &
               FMT='(5X,A12," : ",4X,F5.2,"s CPU ")' )&
             clock_label(n), nsec
        !
     ENDIF
     IF ( mday > 0 ) THEN
        !
        WRITE( stdout, '(1X,I2,"d",I2,"h",I2,"m WALL"/)' ) &
             mday, mhour, mmin
        !
     ELSEIF ( mhour > 0 ) THEN
        !
        WRITE( stdout, '(4X,I2,"h",I2,"m WALL"/)' ) &
             mhour, mmin
        !
     ELSEIF ( mmin > 0 ) THEN
        !
        WRITE( stdout, '(1X,I2,"m",F5.2,"s WALL"/)' ) &
             mmin, msec
        !
     ELSE
        !
        WRITE( stdout, '(4X,F5.2,"s WALL"/)' ) &
             msec
        !
     ENDIF
#endif
     !
  ELSEIF ( nmax == 1 .or. t0cpu(n) /= notrunning ) THEN
     !
     ! ... for clocks that have been called only once
     !
     WRITE( stdout, &
            '(5X,A12," : ",F9.2,"s CPU ",F9.2,"s WALL (",I8," calls)")' ) &
                clock_label(n), elapsed_cpu_time, elapsed_wall_time, nmax
     !
  ELSEIF ( nmax == 0 ) THEN
     !
     ! ... for clocks that have never been called
     !
     WRITE( stdout, &
            '("print_this: clock # ",I2," for ",A12," never called !"/)' ) &
                n, clock_label(n)
     !
  ELSE
     !
     ! ... for all other clocks
     !
     WRITE( stdout, &
        '(5X,A12," : ",F9.2,"s CPU ",F9.2,"s WALL (",I8," calls)")' ) &
        clock_label(n), elapsed_cpu_time, elapsed_wall_time, nmax
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE print_this_clock
!----------------------------------------------------------------------------
SUBROUTINE print_this_clock_gpu( n )
  !----------------------------------------------------------------------------
  !
  USE util_param, ONLY : DP, stdout
  USE mytime,     ONLY : clock_label, gputime, called, gpu_called, mpi_per_thread                
  !
  IMPLICIT NONE
  !
  INTEGER  :: n
  REAL(DP) :: elapsed_gpu_time
  INTEGER  :: nday, nhour, nmin, nmax, mday, mhour, mmin
  !
  !
  !
  nmax = gpu_called(n)
  elapsed_gpu_time = gputime(n) / 1000.d0 ! GPU times are stored in ms
  if (nmax == 0) return
  !
  ! ... In the parallel case there are several possible approaches
  ! ... The safest one is to leave each clock independent from the others
  ! ... Another possibility is to print the maximum across all processors
  ! ... This is done by uncommenting the following lines
  !
  ! CALL mp_max( elapsed_cpu_time, intra_image_comm )
  ! CALL mp_max( elapsed_wall_time, intra_image_comm )
  ! CALL mp_max( nmax, intra_image_comm )
  !
  ! ... In the last line we assume that the maximum cpu time
  ! ... is associated to the maximum number of calls
  ! ... NOTA BENE: by uncommenting the above lines you may run into
  ! ... serious trouble if clocks are not started on all nodes
  !
  IF ( n == 1 ) THEN
     !
     ! ... The first clock is written as days/hour/min/sec
     !
     WRITE( stdout, &
        '(5X,A12," : ",F9.2,"s GPU "/)' ) &
        clock_label(n), elapsed_gpu_time
     !
     !
  ELSE
     !
     ! ... for all other clocks
     !
     WRITE( stdout, &
        '(35X,F9.2,"s GPU  (",I8," calls)")' ) &
        elapsed_gpu_time, nmax
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE print_this_clock_gpu
!
!----------------------------------------------------------------------------
FUNCTION get_clock( label )
  !----------------------------------------------------------------------------
  !
  USE util_param, ONLY : DP
  USE mytime,     ONLY : no, nclock, clock_label, walltime, &
                         notrunning, t0wall, t0cpu, f_wall
  !
  IMPLICIT NONE
  !
  REAL(DP)         :: get_clock
  CHARACTER(len=*) :: label
  INTEGER          :: n
  !
  !
  IF ( no ) THEN
     !
     IF ( label == clock_label(1) ) THEN
        !
        get_clock = f_wall()
        !
     ELSE
        !
        get_clock = notrunning
        !
     ENDIF
     !
     RETURN
     !
  ENDIF
  !
  DO n = 1, nclock
     !
     IF ( label == clock_label(n) ) THEN
        !
        IF ( t0cpu(n) == notrunning ) THEN
           !
           get_clock = walltime(n)
           !
        ELSE
           !
           get_clock = walltime(n) + f_wall() - t0wall(n)
           !
        ENDIF
        !
        ! ... See comments in subroutine print_this_clock about parallel case
        !
        ! CALL mp_max( get_clock, intra_image_comm )
        !
        RETURN
        !
     ENDIF
     !
  ENDDO
  !
  ! ... clock not found
  !
  get_clock = notrunning
  !
  RETURN
  !
END FUNCTION get_clock
