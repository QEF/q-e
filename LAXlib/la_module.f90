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
  SUBROUTINE cdiaghg_cpu_( n, m, h, s, ldh, e, v, offload )
    !----------------------------------------------------------------------------
    !
    ! ... calculates eigenvalues and eigenvectors of the generalized problem
    ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
    ! ... On output both matrix are unchanged
    !
    ! ... LAPACK version - uses both ZHEGV and ZHEGVX
    !
    USE la_param,          ONLY : DP
#if defined (__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, m, ldh
      ! dimension of the matrix to be diagonalized
      ! number of eigenstates to be calculate
      ! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
      ! actually intent(in) but compilers don't know and complain
      ! matrix to be diagonalized
      ! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
      ! eigenvalues
    COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
      ! eigenvectors (column-wise)
    LOGICAL, OPTIONAL ::  offload
      ! optionally solve the eigenvalue problem on the GPU
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
      CALL cdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d)
      !
      e = e_d
      v(1:ldh,1:m) = v_d(1:ldh,1:m)
      !
      DEALLOCATE(h_d, s_d, e_d, v_d)
#endif
    ELSE
      CALL cdiaghg(n, m, h, s, ldh, e, v)
    END IF
    !
    RETURN
    !
  END SUBROUTINE cdiaghg_cpu_
  !
#if defined(__CUDA)
  !----------------------------------------------------------------------------
  SUBROUTINE cdiaghg_gpu_( n, m, h_d, s_d, ldh, e_d, v_d, onhost )
    !----------------------------------------------------------------------------
    !
    ! ... calculates eigenvalues and eigenvectors of the generalized problem
    ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
    ! ... On output both matrix are unchanged
    !
    ! ... LAPACK version - uses both ZHEGV and ZHEGVX
    !
    USE la_param,          ONLY : DP
    USE cudafor
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, m, ldh
      ! dimension of the matrix to be diagonalized
      ! number of eigenstates to be calculate
      ! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), DEVICE, INTENT(INOUT) :: h_d(ldh,n), s_d(ldh,n)
      ! actually intent(in) but compilers don't know and complain
      ! matrix to be diagonalized
      ! overlap matrix
    REAL(DP), DEVICE, INTENT(OUT) :: e_d(n)
      ! eigenvalues
    COMPLEX(DP), DEVICE, INTENT(OUT) :: v_d(ldh,n)
      ! eigenvectors (column-wise)
    LOGICAL, OPTIONAL ::  onhost
      ! optionally solve the eigenvalue problem on the CPU
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
      CALL cdiaghg(n, m, h, s, ldh, e, v)
      !
      e_d = e
      v_d(1:ldh,1:m) = v(1:ldh,1:m)
      !
      DEALLOCATE(h, s, e, v)
    ELSE
      CALL cdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d)
    END IF
    !
    RETURN
    !
  END SUBROUTINE cdiaghg_gpu_
#endif
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rdiaghg_cpu_( n, m, h, s, ldh, e, v, offload )
    !----------------------------------------------------------------------------
    !
    ! ... general interface for rdiaghg
    !
    USE la_param,          ONLY : DP
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, m, ldh
      ! dimension of the matrix to be diagonalized
      ! number of eigenstates to be calculate
      ! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
      ! actually intent(in) but compilers don't know and complain
      ! matrix to be diagonalized
      ! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
      ! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,m)
      ! eigenvectors (column-wise)
    LOGICAL, OPTIONAL ::  offload
      ! optionally solve the eigenvalue problem on the GPU
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
      CALL rdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d)
      !
      e = e_d
      v(1:ldh,1:m) = v_d(1:ldh,1:m)
      !
      DEALLOCATE(h_d, s_d, e_d, v_d)
#endif
    ELSE
      CALL rdiaghg(n, m, h, s, ldh, e, v)
    END IF
    !
    RETURN
    !
  END SUBROUTINE rdiaghg_cpu_
  !
#if defined(__CUDA)
  !----------------------------------------------------------------------------
  SUBROUTINE rdiaghg_gpu_( n, m, h_d, s_d, ldh, e_d, v_d, onhost )
    !----------------------------------------------------------------------------
    !
    ! ... General interface to rdiaghg_gpu
    !
    USE la_param,          ONLY : DP
    USE cudafor
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, m, ldh
      ! dimension of the matrix to be diagonalized
      ! number of eigenstates to be calculate
      ! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), DEVICE, INTENT(INOUT) :: h_d(ldh,n), s_d(ldh,n)
      ! actually intent(in) but compilers don't know and complain
      ! matrix to be diagonalized
      ! overlap matrix
    REAL(DP), DEVICE, INTENT(OUT) :: e_d(n)
      ! eigenvalues
    REAL(DP), DEVICE, INTENT(OUT) :: v_d(ldh,n)
      ! eigenvectors (column-wise)
    LOGICAL, OPTIONAL ::  onhost
      ! optionally solve the eigenvalue problem on the CPU
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
      CALL rdiaghg(n, m, h, s, ldh, e, v)
      !
      e_d = e
      v_d(1:ldh,1:m) = v(1:ldh,1:m)
      !
      DEALLOCATE(h, s, e, v)
    ELSE
      CALL rdiaghg_gpu(n, m, h_d, s_d, ldh, e_d, v_d)
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
  SUBROUTINE prdiaghg_( n, h, s, ldh, e, v, desc, offload )
    !----------------------------------------------------------------------------
    !
    ! ... calculates eigenvalues and eigenvectors of the generalized problem
    ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
    ! ... On output both matrix are unchanged
    !
    ! ... Parallel version with full data distribution
    !
    USE la_param,    ONLY : DP
    USE descriptors, ONLY : la_descriptor
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, ldh
      ! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
      ! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
      ! matrix to be diagonalized
      ! overlap matrix
    !
    REAL(DP), INTENT(OUT) :: e(n)
      ! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,ldh)
      ! eigenvectors (column-wise)
    TYPE(la_descriptor), INTENT(IN) :: desc
      !
    LOGICAL, OPTIONAL ::  offload
      ! place-holder, offloading on GPU not implemented yet
    LOGICAL :: loffload

    CALL prdiaghg( n, h, s, ldh, e, v, desc)
      
  END SUBROUTINE
  !----------------------------------------------------------------------------
  SUBROUTINE pcdiaghg_( n, h, s, ldh, e, v, desc, offload )
    !----------------------------------------------------------------------------
    !
    ! ... calculates eigenvalues and eigenvectors of the generalized problem
    ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
    ! ... On output both matrix are unchanged
    !
    ! ... Parallel version with full data distribution
    !
    USE la_param,    ONLY : DP
    USE descriptors, ONLY : la_descriptor
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, ldh
      ! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
      ! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
      ! matrix to be diagonalized
      ! overlap matrix
    !
    REAL(DP), INTENT(OUT) :: e(n)
      ! eigenvalues
    COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
      ! eigenvectors (column-wise)
    TYPE(la_descriptor), INTENT(IN) :: desc
      !
    LOGICAL, OPTIONAL ::  offload
      ! place-holder, offloading on GPU not implemented yet
    LOGICAL :: loffload

    CALL pcdiaghg( n, h, s, ldh, e, v, desc)
      
  END SUBROUTINE
  !
#if defined(__CUDA)
  !----------------------------------------------------------------------------
  SUBROUTINE prdiaghg__gpu( n, h_d, s_d, ldh, e_d, v_d, desc, onhost )
    !----------------------------------------------------------------------------
    !
    ! ... calculates eigenvalues and eigenvectors of the generalized problem
    ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
    ! ... On output both matrix are unchanged
    !
    ! ... Parallel version with full data distribution
    !
    USE la_param,    ONLY : DP
    USE descriptors, ONLY : la_descriptor    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, ldh
      ! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
      ! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT), DEVICE :: h_d(ldh,ldh), s_d(ldh,ldh)
      ! matrix to be diagonalized
      ! overlap matrix
    !
    REAL(DP), INTENT(OUT), DEVICE :: e_d(n)
      ! eigenvalues
    REAL(DP), INTENT(OUT), DEVICE :: v_d(ldh,ldh)
      ! eigenvectors (column-wise)
    TYPE(la_descriptor), INTENT(IN) :: desc
      !
    LOGICAL, OPTIONAL ::  onhost
      ! place-holder, prdiaghg on GPU not implemented yet
    LOGICAL :: lonhost
      !
    REAL(DP), ALLOCATABLE :: v(:,:), h(:,:), s(:,:)
    REAL(DP), ALLOCATABLE :: e(:)
    
    ALLOCATE(h(ldh,ldh), s(ldh,ldh), e(n), v(ldh,ldh))
    h = h_d; s = s_d;
    CALL prdiaghg( n, h, s, ldh, e, v, desc)
    e_d = e; v_d = v
    DEALLOCATE(h,s,v,e)
    ! 
  END SUBROUTINE
  !----------------------------------------------------------------------------
  SUBROUTINE pcdiaghg__gpu( n, h_d, s_d, ldh, e_d, v_d, desc, onhost )
    !----------------------------------------------------------------------------
    !
    ! ... calculates eigenvalues and eigenvectors of the generalized problem
    ! ... Hv=eSv, with H symmetric matrix, S overlap matrix.
    ! ... On output both matrix are unchanged
    !
    ! ... Parallel version with full data distribution
    !
    USE la_param,    ONLY : DP
    USE descriptors, ONLY : la_descriptor    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n, ldh
      ! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
      ! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT), DEVICE :: h_d(ldh,ldh), s_d(ldh,ldh)
      ! matrix to be diagonalized
      ! overlap matrix
    !
    REAL(DP), INTENT(OUT), DEVICE :: e_d(n)
      ! eigenvalues
    COMPLEX(DP), INTENT(OUT), DEVICE :: v_d(ldh,ldh)
      ! eigenvectors (column-wise)
    TYPE(la_descriptor), INTENT(IN) :: desc
      !
    LOGICAL, OPTIONAL ::  onhost
      ! place-holder, pcdiaghg on GPU not implemented yet
    LOGICAL :: lonhost
      !
    COMPLEX(DP), ALLOCATABLE :: v(:,:), h(:,:), s(:,:)
    REAL(DP), ALLOCATABLE :: e(:)
    
    ALLOCATE(h(ldh,ldh), s(ldh,ldh), e(n), v(ldh,ldh))
    h = h_d; s = s_d;
    CALL pcdiaghg( n, h, s, ldh, e, v, desc)
    e_d = e; v_d = v
    DEALLOCATE(h,s,v,e)
    !
  END SUBROUTINE
#endif
END MODULE LAXlib
