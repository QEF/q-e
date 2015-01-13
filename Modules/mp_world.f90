!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_world
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_start, mp_end, mp_stop
  USE io_global, ONLY : meta_ionode_id, meta_ionode
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... World group - all QE routines using mp_world_start to start MPI
  ! ... will work in the communicator passed as input to mp_world_start
  !
  INTEGER :: nproc = 1  ! number of processors
  INTEGER :: mpime = 0  ! processor index (starts from 0 to nproc-1)
  INTEGER :: root  = 0  ! index of the root processor
  INTEGER :: world_comm = 0  ! communicator
  !
  ! ... library_mode =.true. if QE is called as a library by an external code
  ! ... if true, MPI_Init()     is not called when starting MPI,
  ! ...          MPI_Finalize() is not called when stopping MPI
  !
  LOGICAL :: library_mode = .FALSE.
  !
  PRIVATE
  PUBLIC ::nproc, mpime, root, world_comm, mp_world_start, mp_world_end
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_world_start ( my_world_comm )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: my_world_comm
    INTEGER :: ierr

    INTEGER :: PROVIDED    
    INTEGER :: MPI_THREAD_FUNNELED
    !
    world_comm = my_world_comm
    !
    ! ... check if mpi is already initialized (library mode) or not
    ! 
#if defined(__MPI)
    CALL mpi_initialized ( library_mode, ierr)
    IF (ierr/=0) CALL mp_stop( 8000 )
    IF (.NOT. library_mode ) THEN
#if defined(__OPENMP)
       CALL mpi_init_thread(MPI_THREAD_FUNNELED,PROVIDED,ierr)
#else
       CALL mpi_init(ierr)
#endif
       IF (ierr/=0) CALL mp_stop( 8001 )
    END IF
#endif
    !
    CALL mp_start( nproc, mpime, world_comm )
    !
    ! ... meta_ionode is true if this processor is the root processor
    ! ... of the world group - "ionode_world" would be a better name
    ! ... meta_ionode_id is the index of such processor
    !
    meta_ionode = ( mpime == root )
    meta_ionode_id = root
    !
    RETURN
    !
  END SUBROUTINE mp_world_start
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_world_end ( )
    !-----------------------------------------------------------------------
    INTEGER :: ierr
    !
    CALL mp_barrier( world_comm )
    CALL mp_end ( world_comm )
#if defined(__MPI)
    IF (.NOT. library_mode ) THEN
       CALL mpi_finalize(ierr)
       IF (ierr/=0) CALL mp_stop( 8002 )
    END IF
#endif
    !
  END SUBROUTINE mp_world_end
  !
END MODULE mp_world
