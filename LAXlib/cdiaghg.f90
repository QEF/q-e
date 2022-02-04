!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_cdiaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  !----------------------------------------------------------------------------
  !
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem.
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! LAPACK version - uses both ZHEGV and ZHEGVX
  !!
  !
  USE laxlib_parallel_include
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n)
  !! matrix to be diagonalized
  COMPLEX(DP), INTENT(INOUT) :: s(ldh,n)
  !! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
  !! eigenvectors (column-wise)
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  !
  INTEGER                  :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 1 .OR. nb >= n) THEN
        !
        lwork = 2*n
        !
     ELSE
        !
        lwork = ( nb + 1 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ALLOCATE( rwork( 3*n - 2 ) )
        !
        ! ... calculate all eigenvalues (overwritten to v)
        !
        v(:,:) = h(:,:)
        !
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
        !
     ELSE
        !
        ALLOCATE( rwork( 7*n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) )
        END DO
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... calculate only m lowest eigenvalues
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        ! ... the following commented lines calculate optimal lwork
        !
        !lwork = -1
        !
        !CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
        !             0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
        !             work, lwork, rwork, iwork, ifail, info )
        !
        !lwork = INT( work(1) ) + 1
        !
        !IF( lwork > SIZE( work ) ) THEN
        !   DEALLOCATE( work )
        !   ALLOCATE( work( lwork ) )
        !END IF
        !
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     !
     DEALLOCATE( rwork )
     DEALLOCATE( work )
     !
     IF ( info > n ) THEN
        CALL lax_error__( 'cdiaghg', 'S matrix not positive definite', ABS( info ) )
     ELSE IF ( info > 0 ) THEN
        CALL lax_error__( 'cdiaghg', 'eigenvectors failed to converge', ABS( info ) )
     ELSE IF ( info < 0 ) THEN
        CALL lax_error__( 'cdiaghg', 'incorrect call to ZHEGV*', ABS( info ) )
     END IF
     !
     ! ... restore input S matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO
     !
     DEALLOCATE( sdiag )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e', ABS( info ) )
  CALL MPI_BCAST( v, SIZE(v), MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v', ABS( info ) )
#endif
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_cdiaghg
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_cdiaghg_gpu( n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm)
  !----------------------------------------------------------------------------
  !!
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! GPU VERSION.
  !
#if defined(_OPENMP)
  USE omp_lib
#endif
  !
#if defined(__CUDA)
  USE cudafor
  !
  USE cusolverdn
#endif
  !
  USE laxlib_parallel_include
  !
  ! NB: the flag below can be used to decouple LAXlib from devXlib.
  !     This will make devXlib an optional dependency of LAXlib when
  !     the library will be decoupled from QuantumESPRESSO.
#define __USE_GLOBAL_BUFFER
#if defined(__USE_GLOBAL_BUFFER) && defined(__CUDA)
  USE device_fbuff_m,        ONLY : dev=>dev_buf, pin=>pin_buf
#define VARTYPE POINTER
#else
#define VARTYPE ALLOCATABLE
#endif
  !
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h_d(ldh,n)
  !! matrix to be diagonalized, allocated on the GPU
  COMPLEX(DP), INTENT(INOUT) :: s_d(ldh,n)
  !! overlap matrix, allocated on the GPU
  REAL(DP), INTENT(OUT) :: e_d(n)
  !! eigenvalues, , allocated on the GPU
  COMPLEX(DP),  INTENT(OUT) :: v_d(ldh,n)
  !! eigenvectors (column-wise), , allocated on the GPU
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  !
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) :: h_d, s_d, e_d, v_d
#endif
  !
  INTEGER              :: lwork, info
  !
  REAL(DP)             :: abstol
  INTEGER, ALLOCATABLE :: ifail(:)
  INTEGER, VARTYPE     :: iwork(:)
  REAL(DP), VARTYPE    :: rwork(:)
  COMPLEX(DP), VARTYPE :: work(:)
  !
  COMPLEX(DP), VARTYPE :: v_h(:,:)
  REAL(DP), VARTYPE    :: e_h(:)
#if (! defined(__USE_GLOBAL_BUFFER)) && defined(__CUDA)
  ATTRIBUTES( PINNED ) :: work, iwork, rwork, v_h, e_h
#endif
  !
  INTEGER              :: lwork_d, lrwork_d, liwork, lrwork
  REAL(DP), VARTYPE    :: rwork_d(:)
  COMPLEX(DP), VARTYPE :: work_d(:)
  ! various work space
  !
  ! Temp arrays to save H and S.
  REAL(DP), VARTYPE    :: h_diag_d(:), s_diag_d(:)
#if defined(__CUDA)
  ATTRIBUTES( DEVICE ) :: work_d, rwork_d, h_diag_d, s_diag_d
  INTEGER                      :: devInfo_d, h_meig
  ATTRIBUTES( DEVICE )         :: devInfo_d
  TYPE(cusolverDnHandle), SAVE :: cuSolverHandle
  LOGICAL, SAVE                :: cuSolverInitialized = .FALSE.
  !
  COMPLEX(DP), VARTYPE   :: h_bkp_d(:,:), s_bkp_d(:,:)
  ATTRIBUTES( DEVICE )   :: h_bkp_d, s_bkp_d
#endif
  INTEGER :: i, j
#undef VARTYPE
  !
  !
  !
  !
  CALL start_clock_gpu( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
      !
      ! Keeping compatibility for both CUSolver and CustomEigensolver, CUSolver below
      !
#if defined(__CUDA)

#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(h_bkp_d(n,n), s_bkp_d(n,n), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate h_bkp_d or s_bkp_d ', ABS( info ) )
#else
      CALL dev%lock_buffer( h_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate h_bkp_d ', ABS( info ) )
      CALL dev%lock_buffer( s_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate s_bkp_d ', ABS( info ) )
#endif
      !
!$cuf kernel do(2)
      DO j=1,n
         DO i=1,n
            h_bkp_d(i,j) = h_d(i,j)
            s_bkp_d(i,j) = s_d(i,j)
         ENDDO
      ENDDO
      !
#if defined(_OPENMP)
      IF (omp_get_num_threads() > 1) CALL lax_error__( ' cdiaghg_gpu ', 'cdiaghg_gpu is not thread-safe',  ABS( info ) )
#endif
      IF ( .NOT. cuSolverInitialized ) THEN
         info = cusolverDnCreate(cuSolverHandle)
         IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnCreate',  ABS( info ) )
         cuSolverInitialized = .TRUE.
      ENDIF
      !
      info = cusolverDnZhegvdx_bufferSize(cuSolverHandle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                               n, h_d, ldh, s_d, ldh, 0.D0, 0.D0, 1, m, h_meig, e_d, lwork_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnZhegvdx_bufferSize failed ', ABS( info ) )
      !
#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(work_d(1*lwork_d), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
#else
      CALL dev%lock_buffer( work_d,  lwork_d, info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
#endif
      !
      info = cusolverDnZhegvdx(cuSolverHandle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                  n, h_d, ldh, s_d, ldh, 0.D0, 0.D0, 1, m, h_meig, e_d, work_d, lwork, devInfo_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnZhegvdx failed ', ABS( info ) )
!$cuf kernel do(2)
      DO j=1,n
         DO i=1,n
            IF(j <= m) v_d(i,j) = h_d(i,j)
            h_d(i,j) = h_bkp_d(i,j)
            s_d(i,j) = s_bkp_d(i,j)
         ENDDO
      ENDDO
      !
      !
      ! Do not destroy the handle to save the (re)creation time on each call.
      !
      !info = cusolverDnDestroy(cuSolverHandle)
      !IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnDestroy failed ', ABS( info ) )
      !
#if ! defined(__USE_GLOBAL_BUFFER)
      DEALLOCATE(work_d)
      DEALLOCATE(h_bkp_d, s_bkp_d)
#else
      CALL dev%release_buffer( work_d,  info )
      CALL dev%release_buffer( h_bkp_d, info )
      CALL dev%release_buffer( s_bkp_d, info )
#endif
      !
      ! Keeping compatibility for both CUSolver and CustomEigensolver, CustomEigensolver below
      !
#else
     CALL lax_error__( 'cdiaghg', 'Called GPU eigensolver without GPU support', 1 )
#endif
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
#if defined __GPU_MPI
  info = cudaDeviceSynchronize()
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (first)', ABS( info ) )
  CALL MPI_BCAST( e_d, n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
  CALL MPI_BCAST( v_d, ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
  info = cudaDeviceSynchronize() ! this is probably redundant...
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (second)', ABS( info ) )
#else
  ALLOCATE(e_h(n), v_h(ldh,m))
  e_h(1:n) = e_d(1:n)
  v_h(1:ldh, 1:m) = v_d(1:ldh, 1:m)
  CALL MPI_BCAST( e_h, n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
  CALL MPI_BCAST( v_h, ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
  e_d(1:n) = e_h(1:n)
  v_d(1:ldh, 1:m) = v_h(1:ldh, 1:m)
  DEALLOCATE(e_h, v_h)
#endif
#endif
  !
  CALL stop_clock_gpu( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_cdiaghg_gpu
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE laxlib_pcdiaghg( n, h, s, ldh, e, v, idesc )
  !----------------------------------------------------------------------------
  !
  !! Called by pdiaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem.
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! Parallel version with full data distribution
  !!
  !
  USE laxlib_parallel_include
  USE laxlib_descriptor,      ONLY : la_descriptor, laxlib_intarray_to_desc
  USE laxlib_processors_grid, ONLY : ortho_parent_comm
#if defined __SCALAPACK
  USE laxlib_processors_grid, ONLY : ortho_cntx, np_ortho, me_ortho, ortho_comm
  USE zhpev_module,     ONLY : pzheevd_drv
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
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh)
  !! matrix to be diagonalized
  COMPLEX(DP), INTENT(INOUT) :: s(ldh,ldh)
  !! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
  !! eigenvectors (column-wise)
  INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
  !! laxlib descriptor
  !  
  TYPE(la_descriptor) :: desc
  !
  INTEGER, PARAMETER  :: root = 0
  INTEGER             :: nx, info
#if defined __SCALAPACK
  INTEGER             :: descsca( 16 )
#endif
    ! local block size
  COMPLEX(DP), ALLOCATABLE :: ss(:,:), hh(:,:), tt(:,:)
    ! work space used only in parallel diagonalization
  !
  ! ... input s and h are copied so that they are not destroyed
  !
  CALL start_clock( 'cdiaghg' )
  !
  CALL laxlib_intarray_to_desc(desc,idesc)
  !
  IF( desc%active_node > 0 ) THEN
     !
     nx   = desc%nrcx
     !
     IF( nx /= ldh ) &
        CALL lax_error__(" pcdiaghg ", " inconsistent leading dimension ", ldh )
     !
     ALLOCATE( hh( nx, nx ) )
     ALLOCATE( ss( nx, nx ) )
     !
     hh(1:nx,1:nx) = h(1:nx,1:nx)
     ss(1:nx,1:nx) = s(1:nx,1:nx)
     !
  END IF

  CALL start_clock( 'cdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of sl ( L is stored in sl )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined __SCALAPACK
     CALL descinit( descsca, n, n, desc%nrcx, desc%nrcx, 0, 0, ortho_cntx, SIZE( ss, 1 ) , info )
     !
     IF( info /= 0 ) CALL lax_error__( ' cdiaghg ', ' desckinit ', ABS( info ) )
#endif
     !
#if defined __SCALAPACK

     CALL pzpotrf( 'L', n, ss, 1, 1, descsca, info )

     IF( info /= 0 ) CALL lax_error__( ' cdiaghg ', ' problems computing cholesky ', ABS( info ) )
#else
     CALL laxlib_pzpotrf( ss, nx, n, idesc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:choldc' )
  !
  ! ... L is inverted ( sl = L^-1 )
  !
  CALL start_clock( 'cdiaghg:inversion' )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined __SCALAPACK
     !CALL clear_upper_tr( ss )
     ! set to zero the upper triangle of ss
     !
     CALL sqr_setmat( 'U', n, ZERO, ss, size(ss,1), idesc )
     !
     CALL pztrtri( 'L', 'N', n, ss, 1, 1, descsca, info )
     !
     IF( info /= 0 ) CALL lax_error__( ' cdiaghg ', ' problems computing inverse ', ABS( info ) )
#else
     CALL laxlib_pztrtri( ss, nx, n, idesc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:inversion' )
  !
  ! ... vl = L^-1*H
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, idesc )
     !
  END IF
  !
  ! ... hl = ( L^-1*H )*(L^-1)^T
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'C', n, ONE, v, nx, ss, nx, ZERO, hh, nx, idesc )
     !
     ! ensure that "hh" is really Hermitian, it is sufficient to set the diagonal
     ! properly, because only the lower triangle of hh will be used
     ! 
     CALL sqr_setmat( 'H', n, ZERO, hh, size(hh,1), idesc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  !
  IF ( desc%active_node > 0 ) THEN
     ! 
#ifdef TEST_DIAG
     CALL test_drv_begin()
#endif

#if defined(__SCALAPACK)
     !
     CALL pzheevd_drv( .true., n, desc%nrcx, hh, e, ortho_cntx, ortho_comm )
     !
#else
     !
     CALL laxlib_pzheevd( .true., n, idesc, hh, SIZE( hh, 1 ), e )
     !
#endif
     !
#ifdef TEST_DIAG
     CALL test_drv_end()
#endif
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'C', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, idesc )
     !
  END IF
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root, ortho_parent_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'pcdiaghg', 'error broadcasting array e', ABS( info ) )
#endif
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     DEALLOCATE( ss, hh )
  END IF
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE test_drv_begin()
     ALLOCATE( tt( n, n ) )
     CALL laxlib_zsqmcll_x( n, hh, nx, tt, n, desc, desc%comm )
     RETURN
  END SUBROUTINE test_drv_begin
  !
  SUBROUTINE test_drv_end()
     !
     INTEGER :: i, j, k
     COMPLEX(DP), ALLOCATABLE :: diag(:,:)
     !
     IF( desc%myc == 0 .AND. desc%myr == 0 ) THEN

        write( 100, fmt="(A20,2D18.10)" ) ' e code = ', e( 1 ), e( n )
        ALLOCATE( diag( n*(n+1)/2, 1 ) )
        k = 1
        !write( 100, fmt="(I5)" ) n
        DO j = 1, n
           DO i = j, n
              diag( k, 1 ) = tt( i, j )
              !write( 100, fmt="(2I5,2D18.10)" ) i, j, tt( i, j )
              k = k + 1
           END DO
        END DO
        call zhpev_drv( 'V', 'L', N, diag(:,1), e, tt, n )
        write( 100, fmt="(A20,2D18.10)" ) ' e test = ', e( 1 ), e( n )
        !write( 100, * ) 'eigenvalues and eigenvectors'
        DO j = 1, n
           !write( 100, fmt="(1I5,1D18.10,A)" ) j, e( j )
           DO i = 1, n
              !write( 100, fmt="(2I5,2D18.10)" ) i, j, tt( i, j )
           END DO
        END DO
        close(100)
        DEALLOCATE( diag )
     END IF
#if defined __MPI
     CALL MPI_BCAST( tt, SIZE(tt), MPI_DOUBLE_COMPLEX, 0, desc%comm, info )
     IF ( info /= 0 ) &
        CALL lax_error__( 'test_drv_end', 'error broadcasting array e', ABS( info ) )
#endif
     CALL laxlib_zsqmdst_x( n, tt, n, hh, nx, desc )
     DEALLOCATE( tt )
     CALL lax_error__('cdiaghg','stop serial',1)
     RETURN
  END SUBROUTINE test_drv_end
  !
END SUBROUTINE laxlib_pcdiaghg
!
