!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
INTERFACE diaghg
SUBROUTINE laxlib_rdiaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  INTEGER, INTENT(IN) :: n, m, ldh
  REAL(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
  REAL(DP), INTENT(OUT) :: e(n)
  REAL(DP), INTENT(OUT) :: v(ldh,m)
  INTEGER,  INTENT(IN)  :: me_bgrp, root_bgrp, intra_bgrp_comm
END SUBROUTINE
#ifdef __CUDA
SUBROUTINE laxlib_rdiaghg_gpu( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  INTEGER, INTENT(IN) :: n, m, ldh
  REAL(DP), DEVICE, INTENT(INOUT) :: h(ldh,n), s(ldh,n)
  REAL(DP), DEVICE, INTENT(OUT) :: e(n)
  REAL(DP), DEVICE, INTENT(OUT) :: v(ldh,m)
  INTEGER,  INTENT(IN)  :: me_bgrp, root_bgrp, intra_bgrp_comm
END SUBROUTINE
#endif
SUBROUTINE laxlib_cdiaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  INTEGER, INTENT(IN) :: n, m, ldh
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
  REAL(DP), INTENT(OUT) :: e(n)
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
  INTEGER, INTENT(IN) :: me_bgrp, root_bgrp, intra_bgrp_comm
END SUBROUTINE
#ifdef __CUDA
SUBROUTINE laxlib_cdiaghg_gpu( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  INTEGER, INTENT(IN) :: n, m, ldh
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: h(ldh,n), s(ldh,n)
  REAL(DP), DEVICE, INTENT(OUT) :: e(n)
  COMPLEX(DP), DEVICE, INTENT(OUT) :: v(ldh,m)
  INTEGER, INTENT(IN) :: me_bgrp, root_bgrp, intra_bgrp_comm
END SUBROUTINE
#endif
END INTERFACE
!----------------------------------------------------------------------------
INTERFACE pdiaghg
SUBROUTINE laxlib_pcdiaghg( n, h, s, ldh, e, v, idesc )
  IMPLICIT NONE
  include 'laxlib_param.fh'
  include 'laxlib_kinds.fh'
  INTEGER, INTENT(IN) :: n, ldh
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
  REAL(DP), INTENT(OUT) :: e(n)
  COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
  INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE 
SUBROUTINE laxlib_prdiaghg( n, h, s, ldh, e, v, idesc )
  IMPLICIT NONE
  include 'laxlib_param.fh'
  include 'laxlib_kinds.fh'
  INTEGER, INTENT(IN) :: n, ldh
  REAL(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
  REAL(DP), INTENT(OUT) :: e(n)
  REAL(DP), INTENT(OUT) :: v(ldh,ldh)
  INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
END SUBROUTINE 
END INTERFACE
!----------------------------------------------------------------------------
   INTERFACE laxlib_print_matrix
      SUBROUTINE print_lambda_x( lambda, idesc, n, nshow, nudx, ccc, ionode, iunit )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         REAL(DP), INTENT(IN) :: lambda(:,:,:), ccc
         INTEGER, INTENT(IN) :: idesc(:,:)
         INTEGER, INTENT(IN) :: n, nshow, nudx
         LOGICAL, INTENT(IN) :: ionode
         INTEGER, INTENT(IN) :: iunit
      END SUBROUTINE
   END INTERFACE
!----------------------------------------------------------------------------
   INTERFACE laxlib_diagonalize
      SUBROUTINE diagonalize_parallel_x( n, rhos, rhod, s, idesc )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         REAL(DP), INTENT(IN)  :: rhos(:,:) !  input symmetric matrix
         REAL(DP)              :: rhod(:)   !  output eigenvalues
         REAL(DP)              :: s(:,:)    !  output eigenvectors
         INTEGER,  INTENT(IN) :: n         !  size of the global matrix
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
      END SUBROUTINE
      SUBROUTINE diagonalize_serial_x( n, rhos, rhod )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         INTEGER,  INTENT(IN)  :: n
         REAL(DP)              :: rhos(:,:)
         REAL(DP)              :: rhod(:)
      END SUBROUTINE
#ifdef __CUDA
      SUBROUTINE diagonalize_serial_gpu( m, rhos, rhod, s, info )
         IMPLICIT NONE
         include 'laxlib_kinds.fh'
         INTEGER, INTENT(IN) :: m
         REAL(DP), DEVICE, INTENT(IN) :: rhos(:,:)
         REAL(DP), DEVICE, INTENT(OUT) :: rhod(:)
         REAL(DP), DEVICE, INTENT(OUT) :: s(:,:)
         INTEGER, INTENT(OUT) :: info
      END SUBROUTINE
#endif
   END INTERFACE
!----------------------------------------------------------------------------
   INTERFACE desc_init
      !
      SUBROUTINE laxlib_desc_init1( nsiz, nx, la_proc, idesc, rank_ip, idesc_ip )
         IMPLICIT NONE
	 include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
	 INTEGER, INTENT(IN)  :: nsiz
	 INTEGER, INTENT(OUT) :: nx
         LOGICAL, INTENT(OUT) :: la_proc
         INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
         INTEGER, INTENT(OUT), ALLOCATABLE :: rank_ip( :, : )
         INTEGER, INTENT(OUT), ALLOCATABLE :: idesc_ip(:,:,:)
      END SUBROUTINE
      !
      SUBROUTINE laxlib_desc_init2( nsiz, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
         IMPLICIT NONE
	 include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
	 INTEGER, INTENT(IN)  :: nsiz
	 INTEGER, INTENT(OUT) :: nx
         LOGICAL, INTENT(OUT) :: la_proc
         INTEGER, INTENT(OUT) :: idesc(LAX_DESC_SIZE)
         INTEGER, INTENT(OUT), ALLOCATABLE :: rank_ip( :, : )
	 INTEGER, INTENT(OUT), ALLOCATABLE :: irc_ip(:)
         INTEGER, INTENT(OUT), ALLOCATABLE :: nrc_ip(:)
      END SUBROUTINE
   END INTERFACE

