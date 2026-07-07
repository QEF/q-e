!
! Copyright (C) 2003-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_rdiagh( n, m, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  !!----------------------------------------------------------------------------
  !!
  !! Called by diagh interface.
  !! Calculates eigenvalues and eigenvectors of the standard problem.
  !! Solve Hv = ev, with H symmetric matrix.
  !! real matrices version.
  !! On output H matrix is unchanged.
  !!
  !! LAPACK version - uses both DSYEVD and DSYEVX
  !
  USE laxlib_parallel_include
  !
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated (m <= n)
  REAL(DP), INTENT(INOUT) :: h(n,n)
  !! matrix to be diagonalized
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues (only the first m entries are written when m < n)
  REAL(DP), INTENT(OUT) :: v(n,m)
  !! eigenvectors (column-wise)
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  !
  INTEGER               :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)              :: abstol
  INTEGER,  ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP), ALLOCATABLE :: work(:), v_temp(:,:)
  LOGICAL               :: all_eigenvalues
  INTEGER,  EXTERNAL    :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  CALL start_clock( 'rdiagh' )
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     all_eigenvalues = ( m == n )
     !
     ! ... allocate workspace for all eigenvectors
     !
     ALLOCATE( v_temp(n,n) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ! ... use DSYEVD (divide-and-conquer) for all eigenvalues
        !
        ! ... check for optimal block size
        !
        nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
        !
        IF ( nb < 5 .OR. nb >= n ) THEN
           lwork = 8*n
        ELSE
           lwork = ( nb + 3 )*n
        END IF
        !
        ! ... estimate workspace size for DSYEVD
        !
        lwork = MAX( lwork, 1 + 6*n + 2*n*n )
        !
        ALLOCATE( work( lwork ) )
        ALLOCATE( iwork( 3 + 5*n ) )
        !
        ! ... copy H to v_temp (will be overwritten with eigenvectors)
        !
        !$omp parallel do
        DO i = 1, n
           v_temp(1:n,i) = h(1:n,i)
        END DO
        !$omp end parallel do
        !
        CALL DSYEVD( 'V', 'U', n, v_temp, n, e, work, lwork, &
                     iwork, SIZE(iwork), info )
        !
        DEALLOCATE( iwork )
        !
     ELSE
        !
        ! ... use DSYEVX for subset of eigenvalues
        !
        nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
        !
        IF ( nb < 5 .OR. nb >= n ) THEN
           lwork = 8*n
        ELSE
           lwork = ( nb + 3 )*n
        END IF
        !
        ALLOCATE( work( lwork ) )
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... calculate only m lowest eigenvalues
        !
        abstol = 0.D0
        !
        ! ... copy H to v_temp (will be overwritten)
        !
        !$omp parallel do
        DO i = 1, n
           v_temp(1:n,i) = h(1:n,i)
        END DO
        !$omp end parallel do
        !
        CALL DSYEVX( 'V', 'I', 'U', n, v_temp, n, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v_temp, n, &
                     work, lwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
     END IF
     !
     ! ... copy first m eigenvectors to output
     !
     !$omp parallel do
     DO i = 1, m
        v(1:n,i) = v_temp(1:n,i)
     END DO
     !$omp end parallel do
     !
     DEALLOCATE( v_temp )
     DEALLOCATE( work )
     !
     IF ( info > 0 ) THEN
        CALL lax_error__( 'rdiagh', 'eigenvectors failed to converge', ABS( info ) )
     ELSE IF ( info < 0 ) THEN
        CALL lax_error__( 'rdiagh', 'incorrect call to DSYEV*', ABS( info ) )
     END IF
     !
     ! Note: H matrix is preserved (we diagonalized v_temp, not h)
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiagh', 'error broadcasting array e', ABS( info ))
  CALL MPI_BCAST( v, SIZE(v), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'rdiagh', 'error broadcasting array v', ABS( info ))
#endif
  !
  CALL stop_clock( 'rdiagh' )
  !
  RETURN
  !
END SUBROUTINE laxlib_rdiagh

!----------------------------------------------------------------------------
SUBROUTINE laxlib_prdiagh( n, h, e, v, idesc )
  !----------------------------------------------------------------------------
  !
  !! Called by pdiagh interface.
  !! Calculates eigenvalues and eigenvectors of the standard problem.
  !! Solve Hv = ev, with H symmetric matrix.
  !! real matrices version.
  !! On output H matrix is unchanged.
  !!
  !! Parallel version with full data distribution.
  !!
  !! Communicators: ortho_comm is the subset of ortho_parent_comm whose ranks
  !! form the ScaLAPACK/BLACS grid used for the parallel linear algebra.
  !! Output is replicated on all processors of ortho_parent_comm.
  !!
  !
  USE laxlib_parallel_include
  USE laxlib_descriptor,      ONLY : la_descriptor, laxlib_intarray_to_desc
  USE laxlib_processors_grid, ONLY : ortho_parent_comm
#if defined __SCALAPACK
  USE laxlib_processors_grid, ONLY : ortho_cntx, np_ortho, me_ortho, ortho_comm
  USE dspev_module,     ONLY : pdsyevd_drv
#endif
  !
  IMPLICIT NONE
  !
  include 'laxlib_kinds.fh'
  include 'laxlib_param.fh'
  include 'laxlib_mid.fh'
  include 'laxlib_low.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
  REAL(DP), INTENT(INOUT) :: h(n,n)
  !! matrix to be diagonalized; replicated input on the processes of ortho_comm
  !! (only read on ranks with desc%active_node > 0). On output H is unchanged.
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues; replicated output on the processes of ortho_parent_comm
  REAL(DP), INTENT(OUT) :: v(n,n)
  !! eigenvectors (column-wise); replicated output on the processes of
  !! ortho_parent_comm.
  INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
  !! laxlib descriptor
  !
  TYPE(la_descriptor) :: desc
  !
  INTEGER, PARAMETER  :: root = 0
  INTEGER             :: nx, info
#if defined __SCALAPACK
  INTEGER             :: ortho_cntx_loc, ortho_comm_loc
#endif
    ! local block size
  REAL(DP), ALLOCATABLE :: h_dist(:,:)
    ! work space used only in parallel diagonalization
  !
  ! ... input h is copied and distributed for ScaLAPACK diagonalization
  !
  CALL start_clock( 'rdiagh' )
  !
  CALL laxlib_intarray_to_desc(desc,idesc)
  !
  IF( desc%active_node > 0 ) THEN
     !
     nx   = desc%nrcx
     !
     ! Note: Unlike pdiaghg, h here is the full n x n replicated matrix; the
     ! routine distributes it to block-cyclic format internally.
     !
     ! ... allocate distributed matrix
     !
     ALLOCATE( h_dist( nx, nx ) )
     !
     ! ... distribute replicated H matrix to block-cyclic format
     !
     CALL laxlib_dsqmdst_x( n, h, n, h_dist, nx, desc )
     !
     ! ... diagonalize using ScaLAPACK/ELPA
     !
#if defined(__SCALAPACK)
     !
     ortho_cntx_loc = ortho_cntx
     ortho_comm_loc = ortho_comm
     !
     CALL pdsyevd_drv( .true., n, nx, h_dist, SIZE(h_dist,1), e, ortho_cntx_loc, ortho_comm_loc )
     !
#else
     !
     CALL laxlib_pdsyevd( .true., n, idesc, h_dist, SIZE( h_dist, 1 ), e )
     !
#endif
     !
     ! ... collect distributed eigenvectors back to replicated format
     !
     CALL laxlib_dsqmcll_x( n, h_dist, nx, v, n, desc, desc%comm )
     !
     DEALLOCATE( h_dist )
     !
  END IF
  !
  ! ... broadcast e and v: ranks outside the ortho pool
  ! ... (desc%active_node == 0) skip the collect in laxlib_dsqmcll_x.
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root, ortho_parent_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'prdiagh', 'error broadcasting array e', ABS( info ) )
  CALL MPI_BCAST( v, SIZE(v), MPI_DOUBLE_PRECISION, root, ortho_parent_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'prdiagh', 'error broadcasting array v', ABS( info ) )
#endif
  !
  CALL stop_clock( 'rdiagh' )
  !
  RETURN
  !
END SUBROUTINE laxlib_prdiagh
!
