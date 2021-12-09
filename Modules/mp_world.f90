!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_world
  !----------------------------------------------------------------------------
  !! World group - all QE routines using \(\texttt{mp_world_start}\) to 
  !! start MPI will work in the communicator passed as input to 
  !! \(\texttt{mp_world_start}\)
  !
  USE mp, ONLY : mp_barrier, mp_start, mp_end, mp_stop, mp_count_nodes 
  USE io_global, ONLY : meta_ionode_id, meta_ionode
#if defined(__CUDA)
  use cudafor, ONLY : cudaSetDevice, cudaGetDeviceCount, cudaDeviceSynchronize
#endif
  !
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... World group - all QE routines using mp_world_start to start MPI
  ! ... will work in the communicator passed as input to mp_world_start
  !
  INTEGER :: nnode = 1
  !! number of nodes
  INTEGER :: nproc = 1
  !! number of processors
  INTEGER :: mpime = 0
  !! processor index (starts from 0 to nproc-1)
  INTEGER :: root  = 0
  !! index of the root processor
  INTEGER :: world_comm = 0
  !! communicator
  !
#if defined(__MPI)
  LOGICAL :: library_mode = .FALSE.
  !! library\_mode = TRUE if QE is called as a library by an external code.  
  !! If true:  
  !! - \(\texttt{MPI_Init()}\)    is not called when starting MPI;  
  !! - \(\texttt{MPI_Finalize()}\) is not called when stopping MPI.
#endif
  !
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
    INTEGER :: color, key, ndev
#if defined(__MPI) || defined(__CUDA)
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
       CALL MPI_Init_thread(MPI_THREAD_MULTIPLE, PROVIDED, ierr)
#else
       CALL MPI_Init(ierr)
#endif
       IF (ierr/=0) CALL mp_stop( 8001 )
    END IF
#endif
    !
    CALL mp_count_nodes ( nnode, color, key, world_comm )
    !
#if defined(__CUDA)
    ierr = cudaGetDeviceCount( ndev )
    IF (ierr/=0) CALL mp_stop( 9000 + ierr )
    !
    ! WARNING: the OpenMP standard does not guarantee that the thread
    !          pool remains the same in different parallel regions.
    !          However, this is apparently what all major implementations do.
!$omp parallel firstprivate(key, ndev) reduction(max:ierr)
    ierr = cudaSetDevice(mod(key, ndev))
!$omp end parallel
    IF (ierr/=0) CALL mp_stop( 9100 + ierr )
    ierr = cudaDeviceSynchronize()
    IF (ierr/=0) CALL mp_stop( 9200 + ierr )
#if defined(__DEBUG)
    write(*,*) "MPI ", key, " on node ", color, " is using GPU: ", mod(key, ndev)
#endif
#endif
    !
    CALL mp_start( nproc, mpime, world_comm )
    !
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
#if defined(__MPI)
    INTEGER :: ierr
#endif
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
