!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE mdiis
  !--------------------------------------------------------------------------
  !
  ! ... this module keeps data for MDIIS algorism.
  ! ... (A.Kovalenko et al., J. Comput. Chem. 1999, 20, 928-936)
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_sum, mp_bcast, mp_rank
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define data of MDIIS
  TYPE mdiis_type
    INTEGER           :: mbox       ! maximum size of box
    INTEGER           :: nbox       ! number of saved vectors in box
    INTEGER,  POINTER :: ibox(:)    ! index of vectors in box
    INTEGER           :: nvec       ! dimension of vector
    REAL(DP), POINTER :: vbox(:,:)  ! vectors in box
    REAL(DP), POINTER :: vres(:,:)  ! residual vectors
    REAL(DP), POINTER :: rmat(:,:)  ! matrix of dot-products
    REAL(DP), POINTER :: coef(:)    ! coefficients of DIIS
    REAL(DP)          :: eta        ! step radius
    INTEGER           :: next       ! number of extrapolation
  END TYPE mdiis_type
  !
  ! ... public components
  PUBLIC :: mdiis_type
  PUBLIC :: allocate_mdiis
  PUBLIC :: deallocate_mdiis
  PUBLIC :: reset_mdiis
  PUBLIC :: update_by_mdiis
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_mdiis(mdiist, mbox, nvec, eta, next)
    !--------------------------------------------------------------------------
    !
    ! ... initialize mdiis_type
    !
    IMPLICIT NONE
    !
    TYPE(mdiis_type), INTENT(INOUT) :: mdiist
    INTEGER,          INTENT(IN)    :: mbox
    INTEGER,          INTENT(IN)    :: nvec
    REAL(DP),         INTENT(IN)    :: eta
    INTEGER,          INTENT(IN)    :: next
    !
    mdiist%mbox = mbox
    mdiist%nbox = 0
    mdiist%nvec = nvec
    mdiist%eta  = eta
    mdiist%next = next
    !
    ALLOCATE(mdiist%ibox(mbox))
    ALLOCATE(mdiist%rmat(mbox, mbox))
    ALLOCATE(mdiist%coef(mbox))
    IF (nvec > 0) THEN
      ALLOCATE(mdiist%vbox(nvec, mbox))
      ALLOCATE(mdiist%vres(nvec, mbox))
    END IF
    !
  END SUBROUTINE allocate_mdiis
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_mdiis(mdiist)
    !--------------------------------------------------------------------------
    !
    ! ... finalize mdiis_type
    !
    IMPLICIT NONE
    !
    TYPE(mdiis_type), INTENT(INOUT) :: mdiist
    !
    mdiist%mbox = 0
    mdiist%nbox = 0
    mdiist%nvec = 0
    mdiist%eta  = 0.0_DP
    !
    IF(ASSOCIATED(mdiist%ibox)) DEALLOCATE(mdiist%ibox)
    IF(ASSOCIATED(mdiist%rmat)) DEALLOCATE(mdiist%rmat)
    IF(ASSOCIATED(mdiist%coef)) DEALLOCATE(mdiist%coef)
    IF(ASSOCIATED(mdiist%vbox)) DEALLOCATE(mdiist%vbox)
    IF(ASSOCIATED(mdiist%vres)) DEALLOCATE(mdiist%vres)
    !
  END SUBROUTINE deallocate_mdiis
  !
  !--------------------------------------------------------------------------
  SUBROUTINE reset_mdiis(mdiist, keep1)
    !--------------------------------------------------------------------------
    !
    ! ... reset mdiis_type
    ! ...
    ! ... Variables:
    ! ...   keep1: if true, keep the latest vector, and delete the others.
    ! ...          if false, delete all vectors (default).
    !
    IMPLICIT NONE
    TYPE(mdiis_type),  INTENT(INOUT) :: mdiist
    LOGICAL, OPTIONAL, INTENT(IN)    :: keep1
    !
    INTEGER :: ibox1
    LOGICAL :: keep1_
    !
    EXTERNAL :: dcopy
    !
    IF (PRESENT(keep1)) THEN
      keep1_ = keep1
    ELSE
      keep1_ = .FALSE.
    END IF
    !
    IF (keep1_) THEN
      ibox1 = mdiist%ibox(mdiist%nbox)
      mdiist%nbox       = 1
      mdiist%ibox(1)    = 1
      mdiist%rmat(1, 1) = mdiist%rmat(ibox1, ibox1)
      mdiist%coef(1)    = 1.0_DP
      IF (ibox1 /= 1 .AND. mdiist%nvec > 0) THEN
        CALL dcopy(mdiist%nvec, mdiist%vbox(1, ibox1), 1, mdiist%vbox(1, 1), 1)
        CALL dcopy(mdiist%nvec, mdiist%vres(1, ibox1), 1, mdiist%vres(1, 1), 1)
      END IF
      !
    ELSE
      mdiist%nbox = 0
    END IF
    !
  END SUBROUTINE reset_mdiis
  !
  !--------------------------------------------------------------------------
  SUBROUTINE update_by_mdiis(mdiist, vbox1, vres1, comm)
    !--------------------------------------------------------------------------
    !
    ! ... save vector and solve MDIIS-equation to update vector
    ! ...
    ! ... Variables:
    ! ...   vbox1: vector of the latest iteration  (for in)
    ! ...          vector modified by MDIIS method (for out)
    ! ...   vres1: residual of the latest iteration
    ! ...   comm:  MPI's communicator, to sum up dot-products
    !
    IMPLICIT NONE
    !
    TYPE(mdiis_type),  INTENT(INOUT) :: mdiist
    REAL(DP),          INTENT(INOUT) :: vbox1(1:*)
    REAL(DP),          INTENT(IN)    :: vres1(1:*)
    INTEGER, OPTIONAL, INTENT(IN)    :: comm
    !
    INTEGER :: irank
    INTEGER :: ierr
    LOGICAL :: lmpi
    !
    IF (PRESENT(comm)) THEN
      lmpi = .TRUE.
    ELSE
      lmpi = .FALSE.
    END IF
    !
    CALL save_vbox1_vres1()
    !
    CALL make_rmat()
    !
    IF(mdiist%mbox <= 2) THEN
      CALL update_vbox1_extpol()
      RETURN
    END IF
    !
    IF(mdiist%nbox < mdiist%mbox .AND. mdiist%nbox < MIN(mdiist%mbox, mdiist%next)) THEN
      CALL update_vbox1_extpol()
      RETURN
    END IF
    !
    ierr  = 0
    irank = 0
#if defined (__MPI)
    IF (lmpi) THEN
      irank = mp_rank(comm)
    END IF
#endif
    IF (irank == 0) THEN
      CALL solve_mdiis(ierr)
    END IF
    IF (lmpi) THEN
      CALL mp_bcast(ierr, 0, comm)
    END IF
    !
    IF (ierr == 0) THEN
      IF (lmpi) THEN
        CALL mp_bcast(mdiist%coef, 0, comm)
      END IF
      CALL update_vbox1_mdiis()
      !
    ELSE
      CALL update_vbox1_extpol()
      CALL reset_mdiis(mdiist, .TRUE.)
    END IF
    !
  CONTAINS
    !
    SUBROUTINE save_vbox1_vres1()
      IMPLICIT NONE
      INTEGER :: i
      INTEGER :: ibox1
      !
      EXTERNAL :: dcopy
      !
      IF (mdiist%nbox < mdiist%mbox) THEN
        IF (mdiist%nvec > 0) THEN
          CALL dcopy(mdiist%nvec, vbox1(1), 1, mdiist%vbox(1, mdiist%nbox + 1), 1)
          CALL dcopy(mdiist%nvec, vres1(1), 1, mdiist%vres(1, mdiist%nbox + 1), 1)
        END IF
        mdiist%ibox(mdiist%nbox + 1) = mdiist%nbox + 1
        mdiist%nbox = mdiist%nbox + 1
        !
      ELSE
        ibox1 = mdiist%ibox(1)
        IF (mdiist%nvec > 0) THEN
          CALL dcopy(mdiist%nvec, vbox1(1), 1, mdiist%vbox(1, ibox1), 1)
          CALL dcopy(mdiist%nvec, vres1(1), 1, mdiist%vres(1, ibox1), 1)
        END IF
        DO i = 2, mdiist%nbox
          mdiist%ibox(i - 1) = mdiist%ibox(i)
        END DO
        mdiist%ibox(mdiist%nbox) = ibox1
      END IF
    END SUBROUTINE save_vbox1_vres1
    !
    SUBROUTINE make_rmat()
      IMPLICIT NONE
      INTEGER :: i
      INTEGER :: ibox1
      INTEGER :: ibox2
      !
      REAL(DP), EXTERNAL :: ddot
      !
      ibox1 = mdiist%ibox(mdiist%nbox)
      DO i = 1, mdiist%nbox
        ibox2 = mdiist%ibox(i)
        IF (mdiist%nvec > 0) THEN
          mdiist%rmat(ibox2, ibox1) = &
          & ddot(mdiist%nvec, mdiist%vres(1, ibox2), 1, mdiist%vres(1, ibox1), 1)
        ELSE
          mdiist%rmat(ibox2, ibox1) = 0.0_DP
        END IF
      END DO
      !
      IF (lmpi) THEN
        CALL mp_sum(mdiist%rmat(:, ibox1), comm)
      END IF
      !
      DO i = 1, mdiist%nbox
        ibox2 = mdiist%ibox(i)
        mdiist%rmat(ibox1, ibox2) = mdiist%rmat(ibox2, ibox1)
      END DO
    END SUBROUTINE make_rmat
    !
    SUBROUTINE inverse(n, amat, ierr)
      IMPLICIT NONE
      INTEGER,  INTENT(IN)    :: n
      REAL(DP), INTENT(INOUT) :: amat(n, n)
      INTEGER,  INTENT(OUT)   :: ierr
      !
      INTEGER               :: i
      INTEGER               :: lwork
      REAL(DP), ALLOCATABLE :: work(:)
      REAL(DP), ALLOCATABLE :: evec1(:,:)
      REAL(DP), ALLOCATABLE :: evec2(:,:)
      REAL(DP), ALLOCATABLE :: eval(:)
      !
      REAL(DP), PARAMETER :: SMALL_EVAL = 1.0E-30_DP
      !
      EXTERNAL :: dsyev
      EXTERNAL :: dgemm
      !
      ierr = 0
      lwork = 3 * n
      ALLOCATE(work(lwork))
      ALLOCATE(evec1(n, n))
      ALLOCATE(evec2(n, n))
      ALLOCATE(eval(n))
      !
      evec1 = amat
      CALL dsyev('V', 'U', n, evec1, n, eval, work, lwork, ierr)
      IF (ierr /= 0) THEN
        ierr = ABS(ierr)
        GOTO 100
      END IF
      !
      ierr = 0
      DO i = 1, n
        IF (ABS(eval(i)) > SMALL_EVAL) THEN
          evec2(:, i) = evec1(:, i) / eval(i)
        ELSE
          ierr = ierr + 1
          evec2(:, i) = 0.0_DP
        END IF
      END DO
      !
      IF (ierr <= (n / 2)) THEN
        ierr = 0
      END IF
      !
      CALL dgemm('N', 'T', n, n, n, 1.0_DP, evec2, n, evec1, n, 0.0_DP, amat, n)
      !
100   CONTINUE
      DEALLOCATE(work)
      DEALLOCATE(evec1)
      DEALLOCATE(evec2)
      DEALLOCATE(eval)
    END SUBROUTINE inverse
    !
    SUBROUTINE solve_mdiis(ierr)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: ierr
      !
      INTEGER               :: i1
      INTEGER               :: i2
      INTEGER               :: ibox1
      INTEGER               :: ibox2
      REAL(DP), ALLOCATABLE :: dmat(:,:)
      REAL(DP), ALLOCATABLE :: dvec(:)
      !
      ALLOCATE(dmat(mdiist%nbox + 1, mdiist%nbox + 1))
      ALLOCATE(dvec(mdiist%nbox + 1))
      !
      ! ... make dmat
      DO i2 = 1, mdiist%nbox
        ibox2 = mdiist%ibox(i2)
        DO i1 = 1, mdiist%nbox
          ibox1 = mdiist%ibox(i1)
          dmat(i1, i2) = mdiist%rmat(ibox1, ibox2)
        END DO
      END DO
      !
      dmat(              :, mdiist%nbox + 1) = 1.0_DP
      dmat(mdiist%nbox + 1,               :) = 1.0_DP
      dmat(mdiist%nbox + 1, mdiist%nbox + 1) = 0.0_DP
      !
      ! ... solve linear equation
      CALL inverse(mdiist%nbox + 1, dmat, ierr)
      IF (ierr /= 0) THEN
        GOTO 100
      END IF
      !
      dvec(:) = dmat(:, mdiist%nbox + 1)
      dvec(1:mdiist%nbox) = dvec(1:mdiist%nbox) / SUM(dvec(1:mdiist%nbox))
      DO i1 = 1, mdiist%nbox
        ibox1 = mdiist%ibox(i1)
        mdiist%coef(ibox1) = dvec(i1)
      END DO
      !
100   CONTINUE
      DEALLOCATE(dmat)
      DEALLOCATE(dvec)
    END SUBROUTINE solve_mdiis
    !
    SUBROUTINE update_vbox1_mdiis()
      IMPLICIT NONE
      INTEGER               :: i
      INTEGER               :: ibox1
      REAL(DP), ALLOCATABLE :: vadd(:)
      !
      EXTERNAL :: dcopy
      EXTERNAL :: daxpy
      !
      IF (mdiist%nvec <= 0) THEN
        RETURN
      END IF
      !
      ALLOCATE(vadd(mdiist%nvec))
      !
      vbox1(1:mdiist%nvec) = 0.0_DP
      DO i = 1, mdiist%nbox
        ibox1 = mdiist%ibox(i)
        CALL dcopy(mdiist%nvec, mdiist%vbox(1, ibox1), 1, vadd, 1)
        CALL daxpy(mdiist%nvec, mdiist%eta, mdiist%vres(1, ibox1), 1, vadd, 1)
        CALL daxpy(mdiist%nvec, mdiist%coef(ibox1), vadd, 1, vbox1, 1)
      END DO
      !
      DEALLOCATE(vadd)
    END SUBROUTINE update_vbox1_mdiis
    !
    SUBROUTINE update_vbox1_extpol()
      IMPLICIT NONE
      INTEGER               :: ibox1
      INTEGER               :: ibox2
      REAL(DP), ALLOCATABLE :: vadd(:)
      !
      EXTERNAL :: dcopy
      EXTERNAL :: daxpy
      !
      IF (mdiist%nvec <= 0) THEN
        RETURN
      END IF
      !
      IF (mdiist%nbox > 1) THEN
        ALLOCATE(vadd(mdiist%nvec))
        !
        ibox1 = mdiist%ibox(mdiist%nbox - 1)
        ibox2 = mdiist%ibox(mdiist%nbox)
        CALL dcopy(mdiist%nvec, mdiist%vres(1, ibox2), 1, vadd, 1)
        CALL daxpy(mdiist%nvec, +1.0_DP, mdiist%vbox(1, ibox2), 1, vadd, 1)
        CALL daxpy(mdiist%nvec, -1.0_DP, mdiist%vbox(1, ibox1), 1, vadd, 1)
        CALL daxpy(mdiist%nvec, mdiist%eta, vadd, 1, vbox1, 1)
        !
        DEALLOCATE(vadd)
        !
      ELSE
        ibox2 = mdiist%ibox(mdiist%nbox)
        CALL daxpy(mdiist%nvec, mdiist%eta, mdiist%vres(1, ibox2), 1, vbox1, 1)
      END IF
    END SUBROUTINE update_vbox1_extpol
    !
  END SUBROUTINE update_by_mdiis
  !
END MODULE mdiis
