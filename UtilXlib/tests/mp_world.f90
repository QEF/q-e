!----------------------------------------------------------------------------
MODULE mp_world
  !----------------------------------------------------------------------------
  !
#if defined(__MPI)
  USE mpi
#endif
  USE mp, ONLY : mp_barrier, mp_start, mp_end, mp_stop 
  USE mp, ONLY : mp_count_nodes
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... World group - all QE routines using mp_world_start to start MPI
  ! ... will work in the communicator passed as input to mp_world_start
  !
  INTEGER :: nnode = 1 ! number of nodes
  INTEGER :: nproc = 1  ! number of processors
  INTEGER :: mpime = 0  ! processor index (starts from 0 to nproc-1)
  INTEGER :: root  = 0  ! index of the root processor
  INTEGER :: world_comm = 0  ! communicator
  !
  ! ... library_mode =.true. if QE is called as a library by an external code
  ! ... if true, MPI_Init()     is not called when starting MPI,
  ! ...          MPI_Finalize() is not called when stopping MPI
  !
  !
#if defined(__MPI)
  LOGICAL :: library_mode = .FALSE.
#endif
  PRIVATE
  PUBLIC :: nnode, nproc, mpime, root, world_comm, mp_world_start, mp_world_end
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_world_start ( my_world_comm )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: my_world_comm
#if defined(__MPI)
    INTEGER :: ierr
#endif
#if defined(_OPENMP)
    INTEGER :: PROVIDED
#endif
    !
    world_comm = my_world_comm
    !
    ! ... check if mpi is already initialized (library mode) or not
    ! 
#if defined(__MPI)
    CALL MPI_Initialized ( library_mode, ierr)
    IF (ierr/=0) CALL mp_stop( 8000 )
    IF (.NOT. library_mode ) THEN
#if defined(_OPENMP)
       CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
       CALL MPI_Init(ierr)
#endif
       IF (ierr/=0) CALL mp_stop( 8001 )
    END IF
#endif
    !
    CALL mp_start( nproc, mpime, world_comm )
    !CALL mp_count_nodes ( nnode, world_comm )
    !
    RETURN
    !
  END SUBROUTINE mp_world_start
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_world_end ( )
    !-----------------------------------------------------------------------
#if defined(__MPI)
    INTEGER :: ierr
#endif
    !
    CALL mp_barrier( world_comm )
    CALL mp_end ( world_comm, .true. )
#if defined(__MPI)
    CALL mpi_finalize(ierr)
    IF (ierr/=0) CALL mp_stop( 8002 )
#endif
    !
  END SUBROUTINE mp_world_end
  !
END MODULE mp_world
