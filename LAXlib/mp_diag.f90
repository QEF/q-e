!

! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE laxlib_processors_grid
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... linear-algebra group (also known as "ortho" or "diag" group). 
  ! ... Used for parallelization of dense-matrix diagonalization used in
  ! ... iterative diagonalization/orthonormalization, matrix-matrix products
  !
  INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
  INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
  INTEGER :: me_ortho1   = 0  ! task id for the ortho group
  INTEGER :: nproc_ortho = 1  ! size of the ortho group:
  INTEGER :: leg_ortho   = 1  ! the distance in the father communicator
                              ! of two neighbour processors in ortho_comm
  INTEGER :: ortho_comm  = 0  ! communicator for the ortho group
  INTEGER :: ortho_row_comm  = 0  ! communicator for the ortho row group
  INTEGER :: ortho_col_comm  = 0  ! communicator for the ortho col group
  INTEGER :: ortho_comm_id= 0 ! id of the ortho_comm
  INTEGER :: ortho_parent_comm  = 0  ! parent communicator from which ortho group has been created
  !
  INTEGER :: ortho_cntx = -1  ! BLACS context for ortho_comm
  !
  LOGICAL :: do_distr_diag_inside_bgrp = .true. ! whether the distributed diagoalization should be performed
                                                ! at the band group level (bgrp) or at its parent level
  !
  LOGICAL :: lax_is_initialized = .false.
  !
CONTAINS
  !
  SUBROUTINE laxlib_end_drv ( )
    !  
    !  free resources associated to the communicator
    !
    IF( .not.  lax_is_initialized ) THEN
    !   CALL lax_error__( ' laxlib_end ', ' laxlib was not initialized ', 1 )
       WRITE(*,"(' laxlib_end: laxlib was not initialized ')")
       RETURN
    END IF
    !
    CALL laxlib_comm_free( ortho_comm )
    IF(  ortho_comm_id > 0  ) THEN
       CALL laxlib_comm_free( ortho_col_comm )
       CALL laxlib_comm_free( ortho_row_comm )
    ENDIF
#if defined __SCALAPACK
    IF(  ortho_cntx /= -1 ) CALL BLACS_GRIDEXIT( ortho_cntx )
    ortho_cntx = -1
#endif
    !
    lax_is_initialized = .false.
    !
    np_ortho(2) = 1
    me_ortho(2) = 0
    me_ortho1   = 0
    nproc_ortho = 1
    leg_ortho   = 1
    ortho_comm  = 0
    ortho_row_comm  = 0
    ortho_col_comm  = 0
    ortho_comm_id= 0
    ortho_parent_comm  = 0
    ortho_cntx = -1  ! BLACS context for ortho_comm
    do_distr_diag_inside_bgrp = .true.
    !
  END SUBROUTINE laxlib_end_drv
  !
!------------------------------------------------------------------------------!
      FUNCTION laxlib_rank( comm )
        USE laxlib_parallel_include
        IMPLICIT NONE
        INTEGER :: laxlib_rank
        INTEGER, INTENT(IN) :: comm
        INTEGER :: ierr, taskid

        ierr = 0
        taskid = 0
#if defined(__MPI)
        CALL mpi_comm_rank(comm,taskid,ierr)
        IF (ierr/=0) CALL lax_error__( ' laxlib_rank ', ' problem getting MPI rank ', 1 )
#endif
        laxlib_rank = taskid
      END FUNCTION laxlib_rank

!------------------------------------------------------------------------------!
      FUNCTION laxlib_size( comm )
        USE laxlib_parallel_include
        IMPLICIT NONE
        INTEGER :: laxlib_size
        INTEGER, INTENT(IN) :: comm
        INTEGER :: ierr, numtask

        ierr = 0
        numtask = 1
#if defined(__MPI)
        CALL mpi_comm_size(comm,numtask,ierr)
        IF (ierr/=0) CALL lax_error__( ' laxlib_size ', ' problem getting MPI size ', 1 )
#endif
        laxlib_size = numtask
      END FUNCTION laxlib_size

      SUBROUTINE laxlib_comm_split( old_comm, color, key, new_comm )
         USE laxlib_parallel_include
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: old_comm
         INTEGER, INTENT (IN) :: color, key
         INTEGER, INTENT (OUT) :: new_comm
         INTEGER :: ierr
         ierr = 0
#if defined(__MPI)
         CALL MPI_COMM_SPLIT( old_comm, color, key, new_comm, ierr )
         IF (ierr/=0) CALL lax_error__( ' laxlib_comm_split ', ' problem splitting MPI communicator ', 1 )
#else
         new_comm = old_comm
#endif
      END SUBROUTINE  laxlib_comm_split

      SUBROUTINE laxlib_comm_free( comm )
         USE laxlib_parallel_include
         IMPLICIT NONE
         INTEGER, INTENT (INOUT) :: comm
         INTEGER :: ierr
         ierr = 0
#if defined(__MPI)
         IF( comm /= MPI_COMM_NULL ) THEN
            CALL mpi_comm_free( comm, ierr )
            IF (ierr/=0) CALL lax_error__( ' laxlib_comm_free ', ' problem freeing MPI communicator ', 1 )
         END IF
#endif
         RETURN
      END SUBROUTINE laxlib_comm_free

  !
END MODULE laxlib_processors_grid
