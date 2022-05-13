!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

INTERFACE laxlib_start
  SUBROUTINE laxlib_start_drv( ndiag_, parent_comm, do_distr_diag_inside_bgrp_  )
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: ndiag_  ! (IN) input number of procs in the diag group, (OUT) actual number
    INTEGER, INTENT(IN) :: parent_comm ! parallel communicator inside which the distributed linear algebra group
                                       ! communicators are created
    LOGICAL, INTENT(IN) :: do_distr_diag_inside_bgrp_  ! comme son nom l'indique
  END SUBROUTINE 
END INTERFACE laxlib_start

INTERFACE laxlib_getval
SUBROUTINE laxlib_getval_ ( nproc_ortho, leg_ortho, np_ortho, me_ortho, ortho_comm, ortho_row_comm, ortho_col_comm, &
  ortho_comm_id, ortho_parent_comm, ortho_cntx, do_distr_diag_inside_bgrp  )
  IMPLICIT NONE
  INTEGER, OPTIONAL, INTENT(OUT) :: nproc_ortho
  INTEGER, OPTIONAL, INTENT(OUT) :: leg_ortho
  INTEGER, OPTIONAL, INTENT(OUT) :: np_ortho(2)
  INTEGER, OPTIONAL, INTENT(OUT) :: me_ortho(2)
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_row_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_col_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_comm_id
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_parent_comm
  INTEGER, OPTIONAL, INTENT(OUT) :: ortho_cntx
  LOGICAL, OPTIONAL, INTENT(OUT) :: do_distr_diag_inside_bgrp
END SUBROUTINE
END INTERFACE

INTERFACE dspev_drv
   SUBROUTINE dspev_drv_x( JOBZ, UPLO, N, AP, W, Z, LDZ )
     IMPLICIT NONE
     include 'laxlib_kinds.fh'
     CHARACTER ::       JOBZ, UPLO
     INTEGER   ::       LDZ, N
     REAL(DP) ::  AP( * ), W( * ), Z( LDZ, * )
   END SUBROUTINE
   SUBROUTINE pdspev_drv_x ( jobz, ap, lda, w, z, ldz, nrl, n, nproc, mpime, comm )
     IMPLICIT NONE
     include 'laxlib_kinds.fh'
     CHARACTER, INTENT(IN) :: JOBZ
     INTEGER, INTENT(IN) :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, INTENT(IN) :: comm
     REAL(DP) :: ap( lda, * ), w( * ), z( ldz, * )
   END SUBROUTINE
END INTERFACE

INTERFACE zhpev_drv
   SUBROUTINE zhpev_drv_x( JOBZ, UPLO, N, AP, W, Z, LDZ )
     IMPLICIT NONE
     include 'laxlib_kinds.fh'
     CHARACTER ::       JOBZ, UPLO
     INTEGER   ::       LDZ, N
     COMPLEX(DP) ::  AP( * ), Z( LDZ, * )
     REAL(DP) ::  W( * )
   END SUBROUTINE 
   SUBROUTINE pzhpev_drv_x( jobz, ap, lda, w, z, ldz, nrl, n, nproc, mpime, comm )
     IMPLICIT NONE
     include 'laxlib_kinds.fh'
     CHARACTER :: JOBZ
     INTEGER, INTENT(IN) :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, INTENT(IN) :: comm
     COMPLEX(DP) :: ap( lda, * ), z( ldz, * )
     REAL(DP) :: w( * )
   END SUBROUTINE 
END INTERFACE

   INTERFACE distribute_lambda
      SUBROUTINE distribute_lambda_x( lambda_repl, lambda_dist, idesc )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         include 'laxlib_param.fh'
         REAL(DP), INTENT(IN)  :: lambda_repl(:,:)
         REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
         INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE distribute_lambda_x
   END INTERFACE

   INTERFACE collect_lambda
      SUBROUTINE collect_lambda_x( lambda_repl, lambda_dist, idesc )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         include 'laxlib_param.fh'
         REAL(DP), INTENT(OUT) :: lambda_repl(:,:)
         REAL(DP), INTENT(IN)  :: lambda_dist(:,:)
         INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE collect_lambda_x
   END INTERFACE

   INTERFACE setval_lambda
      SUBROUTINE setval_lambda_x( lambda_dist, i, j, val, idesc )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         include 'laxlib_param.fh'
         REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
         INTEGER,  INTENT(IN)  :: i, j
         REAL(DP), INTENT(IN)  :: val
         INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE setval_lambda_x
   END INTERFACE

   INTERFACE distribute_zmat
      SUBROUTINE distribute_zmat_x( zmat_repl, zmat_dist, idesc )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         include 'laxlib_param.fh'
         REAL(DP), INTENT(IN)  :: zmat_repl(:,:)
         REAL(DP), INTENT(OUT) :: zmat_dist(:,:)
         INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE distribute_zmat_x
   END INTERFACE

   INTERFACE collect_zmat
      SUBROUTINE collect_zmat_x( zmat_repl, zmat_dist, idesc )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         include 'laxlib_param.fh'
         REAL(DP), INTENT(OUT) :: zmat_repl(:,:)
         REAL(DP), INTENT(IN)  :: zmat_dist(:,:)
         INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE collect_zmat_x
   END INTERFACE

   INTERFACE laxlib_init_desc
      SUBROUTINE laxlib_init_desc_x( idesc, n, nx, np, me, comm, cntx, includeme )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
         INTEGER, INTENT(IN)  :: n   !  the size of this matrix
         INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing this descriptor or the same data distribution
         INTEGER, INTENT(IN)  :: np(2), me(2), comm, cntx
         INTEGER, INTENT(IN)  :: includeme
      END SUBROUTINE
      SUBROUTINE laxlib_multi_init_desc_x( idesc, idesc_ip, rank_ip, n, nx  )
		IMPLICIT NONE
		include 'laxlib_param.fh'
		INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
		INTEGER, INTENT(OUT) :: idesc_ip(:,:,:)
		INTEGER, INTENT(OUT) :: rank_ip(:,:)
		INTEGER, INTENT(IN)  :: n   !  the size of this matrix
		INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing this descriptor or the same data distribution
      END SUBROUTINE
   END INTERFACE

   INTERFACE laxlib_local_dims
      SUBROUTINE descla_local_dims( i2g, nl, n, nx, np, me )
         IMPLICIT NONE
         INTEGER, INTENT(OUT) :: i2g  !  global index of the first local element
         INTEGER, INTENT(OUT) :: nl   !  local number of elements
         INTEGER, INTENT(IN)  :: n    !  number of actual element in the global array
         INTEGER, INTENT(IN)  :: nx   !  dimension of the global array (nx>=n) to be distributed
         INTEGER, INTENT(IN)  :: np   !  number of processors
         INTEGER, INTENT(IN)  :: me   !  taskid for which i2g and nl are computed
      END SUBROUTINE
   END INTERFACE

   INTERFACE blk2cyc_redist
      SUBROUTINE blk2cyc_redist_x( n, a, lda, nca, b, ldb, ncb, idesc )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         INTEGER, INTENT(IN) :: n
         INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
         REAL(DP) :: a( lda, nca ), b( ldb, ncb )
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE
      SUBROUTINE blk2cyc_zredist_x( n, a, lda, nca, b, ldb, ncb, idesc )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         INTEGER, INTENT(IN) :: n
         INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
         COMPLEX(DP) :: a( lda, nca ), b( ldb, ncb )
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE
   END INTERFACE

   INTERFACE cyc2blk_redist
      SUBROUTINE cyc2blk_redist_x( n, a, lda, nca, b, ldb, ncb, idesc )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         INTEGER, INTENT(IN) :: n
         INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
         REAL(DP) :: a( lda, nca ), b( ldb, ncb )
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE
      SUBROUTINE cyc2blk_zredist_x( n, a, lda, nca, b, ldb, ncb, idesc )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         INTEGER, INTENT(IN) :: n
         INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
         COMPLEX(DP) :: a( lda, nca ), b( ldb, ncb )
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE
   END INTERFACE


INTERFACE laxlib_dsqmred
SUBROUTINE laxlib_dsqmred_x( na, a, lda, idesca, nb, b, ldb, idescb )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT(IN) :: na
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,lda)  !  matrix to be redistributed into b
   INTEGER, INTENT(IN) :: idesca(LAX_DESC_SIZE)
   INTEGER, INTENT(IN) :: nb
   INTEGER, INTENT(IN) :: ldb
   REAL(DP)            :: b(ldb,ldb)
   INTEGER, INTENT(IN) :: idescb(LAX_DESC_SIZE)
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_zsqmred
SUBROUTINE laxlib_zsqmred_x( na, a, lda, idesca, nb, b, ldb, idescb )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT(IN) :: na
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP)         :: a(lda,lda)  !  matrix to be redistributed into b
   INTEGER, INTENT(IN) :: idesca(LAX_DESC_SIZE)
   INTEGER, INTENT(IN) :: nb
   INTEGER, INTENT(IN) :: ldb
   COMPLEX(DP)         :: b(ldb,ldb)
   INTEGER, INTENT(IN) :: idescb(LAX_DESC_SIZE)
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_dsqmsym
SUBROUTINE laxlib_dsqmsym_x( n, a, lda, idesc )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#if defined (__CUDA)
SUBROUTINE laxlib_dsqmsym_gpu_x( n, a, lda, idesc )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   REAL(DP), INTENT(INOUT), DEVICE    :: a(:,:)
   ATTRIBUTES(DEVICE)  :: a
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#endif
END INTERFACE

INTERFACE laxlib_zsqmher
SUBROUTINE laxlib_zsqmher_x( n, a, lda, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP)         :: a(lda,lda)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
END INTERFACE

INTERFACE sqr_setmat
SUBROUTINE sqr_dsetmat_x( what, n, alpha, a, lda, idesc )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   CHARACTER(LEN=1), INTENT(IN) :: what
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: alpha
   INTEGER, INTENT(IN) :: lda
   REAL(DP) :: a(lda,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
SUBROUTINE sqr_zsetmat_x( what, n, alpha, a, lda, idesc )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   CHARACTER(LEN=1), INTENT(IN) :: what
   INTEGER, INTENT(IN) :: n
   COMPLEX(DP), INTENT(IN) :: alpha
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP) :: a(lda,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
END INTERFACE


INTERFACE sqr_mm_cannon
SUBROUTINE sqr_dmm_cannon_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#if defined (__CUDA)
SUBROUTINE sqr_dmm_cannon_gpu_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   REAL(DP), DEVICE :: a(:,:), b(:,:), c(:,:)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#endif
SUBROUTINE sqr_smm_cannon_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   REAL(SP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   REAL(SP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
SUBROUTINE sqr_zmm_cannon_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   COMPLEX(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   COMPLEX(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
END INTERFACE

INTERFACE sqr_tr_cannon
SUBROUTINE sqr_tr_cannon_x( n, a, lda, b, ldb, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP)            :: a(lda,*), b(ldb,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
SUBROUTINE sqr_tr_cannon_sp_x( n, a, lda, b, ldb, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(SP)            :: a(lda,*), b(ldb,*)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#if defined (__CUDA)
SUBROUTINE sqr_tr_cannon_gpu_x( n, a, lda, b, ldb, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP), INTENT(IN),  DEVICE :: a(:,:)
   REAL(DP), INTENT(OUT), DEVICE :: b(:,:)
   INTEGER, INTENT(IN) :: idesc(:)
END SUBROUTINE
#endif
END INTERFACE

INTERFACE redist_row2col
SUBROUTINE redist_row2col_x( n, a, b, ldx, nx, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: ldx, nx
   REAL(DP)            :: a(ldx,nx), b(ldx,nx)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#if defined (__CUDA)
SUBROUTINE redist_row2col_gpu_x( n, a, b, ldx, nx, idesc )
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: ldx, nx
   REAL(DP), DEVICE    :: a(:,:)
   REAL(DP), DEVICE    :: b(:,:)
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE
#endif
END INTERFACE
