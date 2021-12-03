MODULE LAXlib
#ifdef __CUDA
  USE cudafor
#endif
  IMPLICIT NONE
  !
  INTERFACE diaghg
     MODULE PROCEDURE cdiaghg_cpu_, rdiaghg_cpu_
#ifdef __CUDA
     MODULE PROCEDURE cdiaghg_gpu_, rdiaghg_gpu_
#endif
  END INTERFACE
  !
  INTERFACE pdiaghg
     MODULE PROCEDURE pcdiaghg_, prdiaghg_
#ifdef __CUDA
     MODULE PROCEDURE pcdiaghg__gpu, prdiaghg__gpu
#endif
  END INTERFACE
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE cdiaghg_cpu_( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by diaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! complex matrices version.
    !! On output both matrix are unchanged.
    !!
    !! LAPACK version - uses both ZHEGV and ZHEGVX
    !
#if defined (__CUDA)
    USE cudafor
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
    LOGICAL, OPTIONAL ::  offload
    !! optionally solve the eigenvalue problem on the GPU
    LOGICAL :: loffload
      !
#if defined(__CUDA)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: v_d(:,:), h_d(:,:), s_d(:,:)
    REAL(DP),    ALLOCATABLE, DEVICE :: e_d(:)
    INTEGER :: info
#endif
    !
    loffload = .false.
    !
    ! the following ifdef ensures no offload if not compiling from GPU 
#if defined(__CUDA)
    IF (PRESENT(offload)) loffload = offload
#endif
    !
    ! ... always false when compiling without CUDA support
    !
    IF ( loffload ) THEN
#if defined(__CUDA)
      !
      ALLOCATE(s_d, source=s); ALLOCATE(h_d, source=h)
      ALLOCATE(e_d(n), v_d(ldh,n))
      !
      CALL laxlib_cdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm)
      !
      e = e_d
      v(1:ldh,1:m) = v_d(1:ldh,1:m)
      !
      DEALLOCATE(h_d, s_d, e_d, v_d)
#endif
    ELSE
      CALL laxlib_cdiaghg(n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm)
    END IF
    !
    RETURN
    !
  END SUBROUTINE cdiaghg_cpu_
  !
#if defined(__CUDA)
  !----------------------------------------------------------------------------
  SUBROUTINE cdiaghg_gpu_( n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, onhost )
    !----------------------------------------------------------------------------
    !
    !! Called by diaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! complex matrices version.
    !! On output both matrix are unchanged.
    !!
    !! GPU version
    !!
    !
    USE cudafor
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculate
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), DEVICE, INTENT(INOUT) :: h_d(ldh,n)
    !! matrix to be diagonalized
    COMPLEX(DP), DEVICE, INTENT(INOUT) :: s_d(ldh,n)
    !! overlap matrix
    REAL(DP), DEVICE, INTENT(OUT) :: e_d(n)
    !! eigenvalues
    COMPLEX(DP), DEVICE, INTENT(OUT) :: v_d(ldh,n)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    LOGICAL, OPTIONAL ::  onhost
    !! optionally solve the eigenvalue problem on the CPU
    LOGICAL :: lonhost
      !
    COMPLEX(DP), ALLOCATABLE :: v(:,:), h(:,:), s(:,:)
    REAL(DP),    ALLOCATABLE :: e(:)
    !
    INTEGER :: info
    !
    lonhost = .false.
    !
    IF (PRESENT(onhost)) lonhost = onhost
    !
    !
    IF ( lonhost ) THEN
      !
      ALLOCATE(s, source=s_d); ALLOCATE(h, source=h_d)
      ALLOCATE(e(n), v(ldh,m))
      !
      CALL laxlib_cdiaghg(n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm)
      !
      e_d = e
      v_d(1:ldh,1:m) = v(1:ldh,1:m)
      !
      DEALLOCATE(h, s, e, v)
    ELSE
      CALL laxlib_cdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm)
    END IF
    !
    RETURN
    !
  END SUBROUTINE cdiaghg_gpu_
#endif
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rdiaghg_cpu_( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by diaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! real matrices version.
    !! On output both matrix are unchanged.
    !!
    !! LAPACK version - uses both DSYGV and DSYGVX
    !!
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculate
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,n)
    !! matrix to be diagonalized
    REAL(DP), INTENT(INOUT) :: s(ldh,n)
    !! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,m)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    LOGICAL, OPTIONAL ::  offload
    !! optionally solve the eigenvalue problem on the GPU   
    LOGICAL :: loffload
      !
#if defined(__CUDA)
    REAL(DP), ALLOCATABLE, DEVICE :: v_d(:,:), h_d(:,:), s_d(:,:)
    REAL(DP), ALLOCATABLE, DEVICE :: e_d(:)
    INTEGER :: info
#endif
    !
    loffload = .false.
    !
    ! the following ifdef ensures no offload if not compiling from GPU 
#if defined(__CUDA)
    IF (PRESENT(offload)) loffload = offload
#endif
    !
    ! ... always false when compiling without CUDA support
    !
    IF ( loffload ) THEN
#if defined(__CUDA)
      !
      ALLOCATE(s_d, source=s); ALLOCATE(h_d, source=h)
      ALLOCATE(e_d(n), v_d(ldh,n))
      !
      CALL laxlib_rdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm)
      !
      e = e_d
      v(1:ldh,1:m) = v_d(1:ldh,1:m)
      !
      DEALLOCATE(h_d, s_d, e_d, v_d)
#endif
    ELSE
      CALL laxlib_rdiaghg(n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm)
    END IF
    !
    RETURN
    !
  END SUBROUTINE rdiaghg_cpu_
  !
#if defined(__CUDA)
  !----------------------------------------------------------------------------
  SUBROUTINE rdiaghg_gpu_( n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, onhost )
    !----------------------------------------------------------------------------
    !
    !! Called by diaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! real matrices version.
    !! On output both matrix are unchanged.
    !!
    !! GPU version 
    !!
    !
    USE cudafor
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculate
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), DEVICE, INTENT(INOUT) :: h_d(ldh,n)
    !! matrix to be diagonalized
    REAL(DP), DEVICE, INTENT(INOUT) :: s_d(ldh,n)
    !! overlap matrix
    REAL(DP), DEVICE, INTENT(OUT) :: e_d(n)
    !! eigenvalues
    REAL(DP), DEVICE, INTENT(OUT) :: v_d(ldh,n)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    LOGICAL, OPTIONAL ::  onhost
    !! optionally solve the eigenvalue problem on the CPU
    LOGICAL :: lonhost
    !
    REAL(DP), ALLOCATABLE :: v(:,:), h(:,:), s(:,:)
    REAL(DP),    ALLOCATABLE :: e(:)
    !
    INTEGER :: info
    !
    lonhost = .false.
    !
    IF (PRESENT(onhost)) lonhost = onhost
    !
    !
    IF ( lonhost ) THEN
      !
      ALLOCATE(s, source=s_d); ALLOCATE(h, source=h_d)
      ALLOCATE(e(n), v(ldh,m))
      !
      CALL laxlib_rdiaghg(n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm)
      !
      e_d = e
      v_d(1:ldh,1:m) = v(1:ldh,1:m)
      !
      DEALLOCATE(h, s, e, v)
    ELSE
      CALL laxlib_rdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm)
    END IF
    !
    RETURN
    !
  END SUBROUTINE rdiaghg_gpu_
#endif
  !
  !  === Parallel diagonalization interface subroutines
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE prdiaghg_( n, h, s, ldh, e, v, idesc, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! real matrices version.
    !! On output both matrix are unchanged.
    !!
    !! Parallel version with full data distribution
    !!
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,ldh)
    !! matrix to be diagonalized
    REAL(DP), INTENT(INOUT) :: s(ldh,ldh)
    !! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    LOGICAL, OPTIONAL ::  offload
    !! place-holder, offloading on GPU not implemented yet
    LOGICAL :: loffload

    CALL laxlib_prdiaghg( n, h, s, ldh, e, v, idesc)
      
  END SUBROUTINE
  !----------------------------------------------------------------------------
  SUBROUTINE pcdiaghg_( n, h, s, ldh, e, v, idesc, offload )
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
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
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
    LOGICAL, OPTIONAL ::  offload
    !! place-holder, offloading on GPU not implemented yet
    LOGICAL :: loffload

    CALL laxlib_pcdiaghg( n, h, s, ldh, e, v, idesc)
      
  END SUBROUTINE
  !
#if defined(__CUDA)
  !----------------------------------------------------------------------------
  SUBROUTINE prdiaghg__gpu( n, h_d, s_d, ldh, e_d, v_d, idesc, onhost )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! real matrices version.
    !! On output both matrix are unchanged.
    !!
    !! Parallel GPU version with full data distribution
    !!
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT), DEVICE :: h_d(ldh,ldh)
    !! matrix to be diagonalized
    REAL(DP), INTENT(INOUT), DEVICE :: s_d(ldh,ldh)
    !! overlap matrix
    REAL(DP), INTENT(OUT), DEVICE :: e_d(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT), DEVICE :: v_d(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    LOGICAL, OPTIONAL ::  onhost
    !! place-holder, prdiaghg on GPU not implemented yet
    LOGICAL :: lonhost
    !
    REAL(DP), ALLOCATABLE :: v(:,:), h(:,:), s(:,:)
    REAL(DP), ALLOCATABLE :: e(:)
    
    ALLOCATE(h(ldh,ldh), s(ldh,ldh), e(n), v(ldh,ldh))
    h = h_d; s = s_d;
    CALL laxlib_prdiaghg( n, h, s, ldh, e, v, idesc)
    e_d = e; v_d = v
    DEALLOCATE(h,s,v,e)
    ! 
  END SUBROUTINE
  !----------------------------------------------------------------------------
  SUBROUTINE pcdiaghg__gpu( n, h_d, s_d, ldh, e_d, v_d, idesc, onhost )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! complex matrices version.
    !! On output both matrix are unchanged.
    !!
    !! Parallel GPU version with full data distribution
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT), DEVICE :: h_d(ldh,ldh)
    !! matrix to be diagonalized
    COMPLEX(DP), INTENT(INOUT), DEVICE :: s_d(ldh,ldh)
    !! overlap matrix
    REAL(DP), INTENT(OUT), DEVICE :: e_d(n)
    !! eigenvalues
    COMPLEX(DP), INTENT(OUT), DEVICE :: v_d(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    LOGICAL, OPTIONAL ::  onhost
    !! place-holder, pcdiaghg on GPU not implemented yet
    LOGICAL :: lonhost
      !
    COMPLEX(DP), ALLOCATABLE :: v(:,:), h(:,:), s(:,:)
    REAL(DP), ALLOCATABLE :: e(:)
    
    ALLOCATE(h(ldh,ldh), s(ldh,ldh), e(n), v(ldh,ldh))
    h = h_d; s = s_d;
    CALL laxlib_pcdiaghg( n, h, s, ldh, e, v, idesc)
    e_d = e; v_d = v
    DEALLOCATE(h,s,v,e)
    !
  END SUBROUTINE
#endif
END MODULE LAXlib
