!
MODULE test_helpers
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: solve_with_zhegvd, solve_with_dsygvd, &
            verify_generalized_eigenpairs, &
            hermitian, symmetric

  include 'laxlib_kinds.fh'

  INTERFACE verify_generalized_eigenpairs
    MODULE PROCEDURE verify_generalized_eigenpairs_real
    MODULE PROCEDURE verify_generalized_eigenpairs_complex
  END INTERFACE

CONTAINS

SUBROUTINE solve_with_zhegvd(n, v, s, ldh, e)
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  complex(DP) :: v(ldh,n)
  complex(DP) :: s(ldh,n)
  real(DP) :: e(n)
  INTEGER                  :: n
  !
  INTEGER                  :: lwork, lrwork, liwork, info, ldh
  !
  INTEGER,     ALLOCATABLE :: iwork(:)
  REAL(DP),    ALLOCATABLE :: rwork(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
  ! various work space

  !
  ALLOCATE(work(1), rwork(1), iwork(1))
  CALL ZHEGVD( 1, 'V', 'U', n, v, ldh, &
              s, ldh, e, work, -1, rwork, -1, iwork, -1, info )

  IF (info /= 0) print *, "Workspace not computed!"

  lwork = work(1)
  lrwork = rwork(1)
  liwork = iwork(1)

  DEALLOCATE(work, rwork, iwork)
  ALLOCATE(work(lwork), rwork(lrwork), iwork(liwork))

  CALL ZHEGVD( 1, 'V', 'U', n, v, ldh, &
              s, ldh, e, work, lwork, rwork, lrwork, iwork, liwork, info )

  DEALLOCATE(work, rwork, iwork)
  !
END SUBROUTINE solve_with_zhegvd
!
SUBROUTINE solve_with_dsygvd(n, v, s, ldh, e)
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  REAL(DP) :: v(ldh,n)
  REAL(DP) :: s(ldh,n)
  real(DP) :: e(n)
  INTEGER                  :: n
  !
  INTEGER                  :: lwork, liwork, info, ldh
  !
  INTEGER,     ALLOCATABLE :: iwork(:)
  REAL(DP),    ALLOCATABLE :: work(:)
  ! various work space

  !
  ALLOCATE(work(1), iwork(1))
  CALL dsygvd( 1, 'V', 'U', n, v, ldh, &
              s, ldh, e, work, -1, iwork, -1, info )

  IF (info /= 0) print *, "Workspace not computed!"

  lwork = work(1)
  liwork = iwork(1)

  DEALLOCATE(work, iwork)
  ALLOCATE(work(lwork), iwork(liwork))
  !
  CALL dsygvd( 1, 'V', 'U', n, v, ldh, &
              s, ldh, e, work, lwork, iwork, liwork, info )
  !
  DEALLOCATE(work, iwork)
  !
END SUBROUTINE solve_with_dsygvd
!
!----------------------------------------------------------------------------
SUBROUTINE verify_generalized_eigenpairs_real(n, m, h, s, ldh, e, v, max_residual)
!----------------------------------------------------------------------------
! Verify H*v = e*S*v for real symmetric generalized eigenvalue problem
! Returns maximum normalized residual: max_i ||H*v_i - e_i*S*v_i|| / (|e_i| * ||S*v_i||)
  IMPLICIT NONE
  include 'laxlib_kinds.fh'

  INTEGER, INTENT(IN) :: n, m, ldh
  REAL(DP), INTENT(IN) :: h(ldh,n), s(ldh,n), e(m), v(ldh,m)
  REAL(DP), INTENT(OUT) :: max_residual

  REAL(DP), ALLOCATABLE :: hv(:), sv(:)
  REAL(DP) :: svnorm, residual
  REAL(DP), EXTERNAL :: DNRM2
  INTEGER :: i

  ALLOCATE(hv(n), sv(n))
  max_residual = 0.0_DP

  DO i = 1, m
    ! Compute H*v_i using BLAS-2 symmetric matrix-vector product
    CALL DSYMV('U', n, 1.0_DP, h, ldh, v(1,i), 1, 0.0_DP, hv, 1)

    ! Compute S*v_i
    CALL DSYMV('U', n, 1.0_DP, s, ldh, v(1,i), 1, 0.0_DP, sv, 1)

    ! Compute residual: H*v_i - e_i*S*v_i
    hv(1:n) = hv(1:n) - e(i) * sv(1:n)

    ! Compute normalized residual
    svnorm = DNRM2(n, sv, 1)
    residual = DNRM2(n, hv, 1) / (ABS(e(i)) * svnorm)

    max_residual = MAX(max_residual, residual)
  END DO

  DEALLOCATE(hv, sv)
END SUBROUTINE verify_generalized_eigenpairs_real

!----------------------------------------------------------------------------
SUBROUTINE verify_generalized_eigenpairs_complex(n, m, h, s, ldh, e, v, max_residual)
!----------------------------------------------------------------------------
! Verify H*v = e*S*v for complex Hermitian generalized eigenvalue problem
! Returns maximum normalized residual: max_i ||H*v_i - e_i*S*v_i|| / (|e_i| * ||S*v_i||)
  IMPLICIT NONE
  include 'laxlib_kinds.fh'

  INTEGER, INTENT(IN) :: n, m, ldh
  COMPLEX(DP), INTENT(IN) :: h(ldh,n), s(ldh,n), v(ldh,m)
  REAL(DP), INTENT(IN) :: e(m)
  REAL(DP), INTENT(OUT) :: max_residual

  COMPLEX(DP), ALLOCATABLE :: hv(:), sv(:)
  COMPLEX(DP), PARAMETER :: one = (1.0_DP, 0.0_DP), zero = (0.0_DP, 0.0_DP)
  REAL(DP) :: svnorm, residual
  REAL(DP), EXTERNAL :: DZNRM2
  INTEGER :: i

  ALLOCATE(hv(n), sv(n))
  max_residual = 0.0_DP

  DO i = 1, m
    ! Compute H*v_i using BLAS-2 Hermitian matrix-vector product
    CALL ZHEMV('U', n, one, h, ldh, v(1,i), 1, zero, hv, 1)

    ! Compute S*v_i
    CALL ZHEMV('U', n, one, s, ldh, v(1,i), 1, zero, sv, 1)

    ! Compute residual: H*v_i - e_i*S*v_i
    hv(1:n) = hv(1:n) - CMPLX(e(i), 0.0_DP, KIND=DP) * sv(1:n)

    ! Compute normalized residual
    svnorm = DZNRM2(n, sv, 1)
    residual = DZNRM2(n, hv, 1) / (ABS(e(i)) * svnorm)

    max_residual = MAX(max_residual, residual)
  END DO

  DEALLOCATE(hv, sv)
END SUBROUTINE verify_generalized_eigenpairs_complex

!----------------------------------------------------------------------------
SUBROUTINE hermitian(mSize, M)
!----------------------------------------------------------------------------
! Generate a random Hermitian matrix for testing
  IMPLICIT NONE
  include 'laxlib_kinds.fh'

  INTEGER, INTENT(IN) :: mSize
  COMPLEX(DP), INTENT(OUT) :: M(mSize,mSize)
  !
  REAL(DP), ALLOCATABLE :: rnd(:)
  INTEGER :: h, k, j
  !
  ALLOCATE(rnd(mSize*(mSize+1)))
  CALL RANDOM_NUMBER(rnd)
  rnd = 1.d0*rnd - 5.d-1
  !
  M = (0.d0, 0.d0)
  j = 1
  DO k=1,mSize
    DO h=1,mSize
      IF(h>k) THEN
        M(h,k) = CMPLX(rnd(j), rnd(j+1))
        M(k,h) = CONJG(M(h,k))
        j=j+2;
      ELSE IF(k == h) THEN
        M(k,h) = CMPLX(mSize, 0.d0, kind=DP)
      END IF
    END DO
  END DO
  !
  DEALLOCATE(rnd)
  !
END SUBROUTINE hermitian

!----------------------------------------------------------------------------
SUBROUTINE symmetric(mSize, M)
!----------------------------------------------------------------------------
! Generate a random symmetric matrix for testing
  IMPLICIT NONE
  include 'laxlib_kinds.fh'

  INTEGER, INTENT(IN) :: mSize
  REAL(DP), INTENT(OUT) :: M(mSize,mSize)
  !
  REAL(DP), ALLOCATABLE :: rnd(:)
  INTEGER :: h, k, j
  !
  ALLOCATE(rnd(mSize*(mSize+1)/2))
  CALL RANDOM_NUMBER(rnd)
  rnd = 1.d0*rnd - 5.d-1
  !
  M = 0.d0
  j = 1
  DO k=1,mSize
    DO h=1,mSize
      IF(h>k) THEN
        M(h,k) = rnd(j)
        M(k,h) = M(h,k)
        j=j+1;
      ELSE IF(k == h) THEN
        M(k,h) = REAL(mSize, kind=DP)
      END IF
    END DO
  END DO
  !
  DEALLOCATE(rnd)
  !
END SUBROUTINE symmetric

END MODULE test_helpers
!