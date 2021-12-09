  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE low_lvl
  !----------------------------------------------------------------------
  !!
  !! This module contains low level routines that are used throughout EPW.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    LOGICAL FUNCTION hslt(a, b, eps)
    !-----------------------------------------------------------------------
    !!
    !! Compare two real number and return the result
    !!
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: a
    !! Input number a
    REAL(KIND = DP), INTENT(in) :: b
    !! Input number b
    REAL(KIND = DP), INTENT(in) :: eps
    !! Tolerence
    !
    IF (ABS(a - b) < eps) THEN
      hslt = .FALSE.
    ELSE
      hslt = (a < b )
    ENDIF
    !-----------------------------------------------------------------------
    END FUNCTION hslt
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION s()
    !-----------------------------------------------------------------------
    !!
    !! s-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP) :: s

    s = 1.d0 / DSQRT(fpi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION s
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION px(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! p-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! Phi
    REAL(KIND = DP) :: px
    !! Output
    !
    ! Local variable
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    px =  DSQRT(3.d0 / fpi) * sint * cos(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION px
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION py(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! p-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! Cos
    REAL(KIND = DP), INTENT(in) :: phi
    !! Phi
    REAL(KIND = DP) :: py
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! Sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    py =  DSQRT(3.d0 / fpi) * sint * sin(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION py
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION p_z(cost)
    !-----------------------------------------------------------------------
    !!
    !! p-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP) :: p_z
    !! Output
    !
    p_z =  DSQRT(3.d0 / fpi) * cost
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION p_z
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dz2(cost)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds, ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP) :: dz2
    !! Output
    !
    dz2 =  DSQRT(1.25d0 / fpi) * (3.d0 * cost * cost - 1.d0)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dz2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dxz(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dxz
    !! Output
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dxz =  DSQRT(15.d0 / fpi) * sint * cost * COS(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dxz
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dyz(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dyz
    !! output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dyz =  DSQRT(15.d0 / fpi) * sint * cost * SIN(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dyz
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dx2my2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dx2my2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dx2my2 =  DSQRT(3.75d0 / fpi) * sint * sint * COS(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dx2my2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dxy(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dxy
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dxy =  DSQRT(3.75d0 / fpi) * sint * sint * SIN(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dxy
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fz3(cost)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP) :: fz3
    !! Output
    !
    fz3 = 0.25d0 * DSQRT(7.d0 / pi) * (5.d0 * cost * cost - 3.d0) * cost
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fz3
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fxz2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fxz2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fxz2 = 0.25d0 * DSQRT(10.5d0 / pi) * (5.d0 * cost * cost - 1.d0) * sint * COS(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fxz2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fyz2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fyz2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fyz2 = 0.25d0 * DSQRT(10.5d0 / pi) * (5.d0 * cost * cost - 1.d0) * sint * SIN(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fyz2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fzx2my2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fzx2my2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fzx2my2 = 0.25d0 * DSQRT(105d0/pi) * sint * sint * cost * COS(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fzx2my2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fxyz(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fxyz
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fxyz = 0.25d0 * DSQRT(105d0 / pi) * sint * sint * cost * SIN(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fxyz
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fxx2m3y2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fxx2m3y2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fxx2m3y2 = 0.25d0 * DSQRT(17.5d0 / pi) * sint * sint * sint * COS(3.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fxx2m3y2
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION fy3x2my2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fy3x2my2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fy3x2my2 = 0.25d0 * DSQRT(17.5d0 / pi) * sint * sint * sint * SIN(3.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fy3x2my2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE init_random_seed()
    !-----------------------------------------------------------------------
    !!
    !! Create seeds for random number generation
    !!
    !
    IMPLICIT NONE
    !
    INTEGER :: i
    !! Division by number running from 1 to n
    INTEGER :: n
    !! Random number
    INTEGER :: clock
    !! Clock count
    INTEGER :: ierr
    !! Error status
    INTEGER, ALLOCATABLE :: seed(:)
    !! Seeds
    !
    CALL RANDOM_SEED(SIZE = n)
    ALLOCATE(seed(n), STAT = ierr)
    IF (ierr /= 0) CALL errore('init_random_seed', 'Error allocating seed', 1)
    !
    CALL SYSTEM_CLOCK(COUNT = clock)
    !
    seed = clock + 37 * (/(i - 1, i = 1, n)/)
    CALL RANDOM_SEED(PUT = seed)
    !
    DEALLOCATE(seed, STAT = ierr)
    IF (ierr /= 0) CALL errore('init_random_seed', 'Error deallocating seed', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE init_random_seed
    !-----------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    SUBROUTINE hpsort_eps_epw(n, ra, ind, eps)
    !---------------------------------------------------------------------
    !! This routine is adapted from flib/hpsort_eps
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! and considering two elements being equal if their values differ
    !! for less than "eps".
    !! n is input, ra is replaced on output by its sorted rearrangement.
    !! create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (ra).
    !! in case of equal values in the data array (ra) the values in the
    !! index array (ind) are used to order the entries.
    !! if on input ind(1)  = 0 then indices are initialized in the routine,
    !! if on input ind(1) != 0 then indices are assumed to have been
    !!                initialized before entering the routine and these
    !!                indices are carried around during the sorting process
    !!
    !! no work space needed !
    !! free us from machine-dependent sorting-routines !
    !!
    !! adapted from Numerical Recipes pg. 329 (new edition)
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: n
    !! Size of the array
    INTEGER, INTENT(inout) :: ind(n)
    !! Array
    REAL(KIND = DP), INTENT(inout) :: ra(n)
    !! Sorted array
    REAL(KIND = DP), INTENT(in) :: eps
    !! Tolerence
    !
    ! Local variables
    INTEGER :: i
    !!
    INTEGER :: ir
    !!
    INTEGER :: j
    !!
    INTEGER :: l
    !!
    INTEGER :: iind
    !!
    REAL(KIND = DP) :: rra
    !! Input array
    !
    ! initialize index array
    IF (ind(1) == 0) THEN
      DO i = 1, n
        ind(i) = i
      ENDDO
    ENDIF
    ! nothing to order
    IF (n < 2) RETURN
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1
    !
    ir = n
    !
    SORTING: DO
      ! still in hiring phase
      IF (l > 1) THEN
        l = l - 1
        rra  = ra(l)
        iind = ind(l)
        ! in retirement-promotion phase.
      ELSE
        ! clear a space at the end of the array
        rra  = ra(ir)
        iind = ind(ir)
        ! retire the top of the heap into it
        ra(ir) = ra(1)
        !
        ind(ir) = ind(1)
        ! decrease the size of the corporation
        ir = ir - 1
        ! done with the last promotion
        IF (ir == 1) THEN
          ! the least competent worker at all !
          ra(1) = rra
          ind(1) = iind
          EXIT sorting
        ENDIF
      ENDIF
      ! Wheter in hiring or promotion phase, we
      i = l
      ! Set up to place rra in its proper level
      j = l + l
      !
      DO WHILE(j <= ir)
        IF (j < ir) THEN
          ! compare to better underling
          IF (hslt(ra(j), ra(j + 1), eps)) THEN
            j = j + 1
          ENDIF
        ENDIF
        ! demote rra
        IF (hslt(rra, ra(j), eps)) THEN
          ra(i) = ra(j)
          ind(i) = ind(j)
          i = j
          j = j + j
        ELSE
          ! set j to terminate do-while loop
          j = ir + 1
        ENDIF
      ENDDO
      ra(i) = rra
      ind(i) = iind
      !
    ENDDO sorting
    !
    !----------------------------------------------------------------------
    END SUBROUTINE hpsort_eps_epw
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------
    SUBROUTINE set_ndnmbr(pool, proc, procp, npool, ndlab)
    !----------------------------------------------------------------
    !!
    !!  create ndlab label from pool and proc numbers
    !!
    !!  The rule for deciding the node number is based on
    !!  the restriction set in startup.f90 that every pool
    !!  has the same number of procs.
    !!
    IMPLICIT NONE
    !
    CHARACTER(LEN = 4), INTENT(out) :: ndlab
    !! Label
    INTEGER, INTENT(in) :: pool
    !! Pool = 1,..., npool
    INTEGER, INTENT(in) :: proc
    !! Processor = 0,..., nproc_pool-1
    INTEGER, INTENT(in) :: procp
    !!
    INTEGER, INTENT(in) :: npool
    !!
    ! Local variables
    INTEGER :: node
    !! Number of nodes
    INTEGER :: nprocs
    !!
    !
    nprocs = npool * procp
    !
    node = (pool - 1) * procp + proc + 1
    !
    ndlab = '    '
    IF (nprocs < 10) THEN
      WRITE(ndlab(1:1), '(i1)') node
    ELSEIF (nprocs < 100 ) then
      IF (node < 10) THEN
        WRITE(ndlab(1:1), '(i1)') node
      ELSE
        WRITE(ndlab(1:2), '(i2)') node
      ENDIF
    ELSEIF (nprocs < 100) THEN
      IF (node < 10) THEN
        WRITE(ndlab(1:1), '(i1)') node
      ELSEIF (node < 100) THEN
        WRITE(ndlab(1:2), '(i2)') node
      ELSE
        WRITE(ndlab(1:3), '(i3)') node
      ENDIF
    ELSE
      IF (node < 10) THEN
        WRITE(ndlab(1:1), '(i1)') node
      ELSEIF (node < 100) THEN
        WRITE(ndlab(1:2), '(i2)') node
      ELSEIF (node < 1000) THEN
        WRITE(ndlab(1:3), '(i3)') node
      ELSE
        WRITE(ndlab(1:4), '(i4)') node
      ENDIF
    ENDIF
    !
    !----------------------------------------------------------------
    END SUBROUTINE set_ndnmbr
    !----------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE system_mem_usage(valueRSS)
    !----------------------------------------------------------------------
    !
    !! Report the memory usage ( VIRT and REAL ) from the current PID
    !! process ( so it will be master only in case of MPI ).
    !! Memory is reported from the /proc/PID_NUMBER/status file
    !
#ifdef __INTEL_COMPILER
    USE ifport !if on intel compiler
#endif
    USE io_global,   ONLY : stdout
    USE io_var,      ONLY : iunimem
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: valueRSS(2)
    !! Contains the value of the memory in kB
    !
    ! Local variables
    CHARACTER(LEN = 200) :: filename = ' '
    !! Name of the file
    CHARACTER(LEN = 80) :: line
    !! Line in the file
    CHARACTER(LEN = 8) :: pid_char = ' '
    !!
    LOGICAL :: ifxst
    !! Does the file exists
#if defined(__PGI) || defined(__CRAY) || defined(__XLF) || defined(__FLANG)
    INTEGER, EXTERNAL :: getpid
    !! PID of the process
#endif
    INTEGER :: pid
    !! PID of the process
    !
    valueRSS = -1    ! return negative number if not found
    !
    ! Get process ID
    !
    pid = getpid()
    WRITE(pid_char, '(I8)') pid
    filename = '/proc/' // TRIM(ADJUSTL(pid_char)) // '/status'
    !
    ! Read system file
    !
    INQUIRE(FILE = filename, EXIST = ifxst)
    IF (.NOT. ifxst) THEN
      WRITE(stdout, '(a)') 'System file does not exist'
      RETURN
    ENDIF
    !
    OPEN(UNIT = iunimem, FILE = filename, ACTION = 'read')
    !
    DO
      READ(iunimem, '(a)', END = 120) line
      ! Peak virtual memory usage
      IF (line(1:7) == 'VmPeak:') THEN
        READ(line(8:), *) valueRSS(1)
      ENDIF
      ! Peak resident set size
      IF (line(1:6) == 'VmHWM:') THEN
        READ(line(7:), *) valueRSS(2)
        CLOSE(UNIT = iunimem, STATUS = 'keep')
        EXIT
      ENDIF
    ENDDO
    120 CONTINUE
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE system_mem_usage
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    FUNCTION utility_zdotu(a, b)
    !--------------------------------------------------------------------------
    !!
    !! Dot product function
    !!
    USE kinds, ONLY: DP
    !
    COMPLEX(KIND = DP), INTENT(in)  :: a(:)
    !! Input vector
    COMPLEX(KIND = DP), INTENT(in)  :: b(:)
    !!
    COMPLEX(KIND = DP) :: utility_zdotu
    !! Output
    !
    utility_zdotu = SUM(a * b)
    !
    RETURN
    !--------------------------------------------------------------------------
    END FUNCTION utility_zdotu
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE degen_sort(input_array, sizes, output, repeat_list)
    !--------------------------------------------------------------------------
    !!
    !! Find degenerate values using a bubble sorting algorithms
    !!
    !! From: https://stackoverflow.com/questions/7502489/bubble-sort-algorithm-javascript/37901176
    !!
    !! On exititing, repeat_list contains 0 for bands that are non-degenerate and
    !! a group index for the one that are.
    !! Example: the following set of eigenenergies from array = [0,0.1,0.1,0.1,0.2,0.3,0.3]
    !!          gives repeat_list = [0,1,1,1,0,2,2]
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : eps8, eps20, eps6
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(out) :: output
    !! Return true on return if degenercies found
    INTEGER, INTENT(in) :: sizes
    !! Sizes of the array
    INTEGER, INTENT(out) :: repeat_list(sizes)
    !! Array containing the degeneracices
    REAL(KIND = DP), INTENT(in) :: input_array(sizes)
    !! Input array
    !
    ! Local variables
    INTEGER :: j
    !! Index of size
    INTEGER :: degen_label
    !! Degen index
    !
    output         = .FALSE.
    degen_label    = 0
    repeat_list(:) = 0
    !
    DO j = 1, sizes - 1
      IF (0.5d0 * (ABS(input_array(j) - input_array(j + 1)) /&
          (ABS(input_array(j)) + ABS(input_array(j + 1)) + eps20)) < eps6) THEN
        IF (j == 1) THEN
          degen_label = 1
        ELSE
          IF (0.5d0 * (ABS(input_array(j) - input_array(j - 1)) /&
              (ABS(input_array(j)) + ABS(input_array(j - 1)) + eps20)) > eps6) THEN
            degen_label = degen_label + 1
          ENDIF
        ENDIF
        repeat_list(j)     = degen_label
        repeat_list(j + 1) = degen_label
        output             = .TRUE.
      ELSE
        repeat_list(j + 1) = 0
      ENDIF
    ENDDO
    !
    RETURN
    !--------------------------------------------------------------------------
    END SUBROUTINE degen_sort
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    FUNCTION matinv3(A) RESULT(B)
    !-----------------------------------------------------------------------
    !!
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : eps160
    !
    REAL(KIND = DP), INTENT(in) :: A(3, 3)
    !! Matrix
    !
    ! Local variable
    REAL(KIND = DP) :: detinv
    !! Inverse of the determinant
    REAL(KIND = DP) :: B(3, 3)
    !! Inverse matrix
    !
    ! Calculate the inverse determinant of the matrix
    detinv = 1 / (A(1, 1) * A(2, 2) * A(3, 3) - A(1, 1) * A(2, 3) * A(3, 2) &
                - A(1, 2) * A(2, 1) * A(3, 3) + A(1, 2) * A(2, 3) * A(3, 1) &
                + A(1, 3) * A(2, 1) * A(3, 2) - A(1, 3) * A(2, 2) * A(3, 1))
    !
    IF (detinv < eps160) THEN
      CALL errore('matinv3', 'Inverse does not exist ', 1)
    ENDIF
    !
    ! Calculate the inverse of the matrix
    B(1, 1) = +detinv * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2))
    B(2, 1) = -detinv * (A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1))
    B(3, 1) = +detinv * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))
    B(1, 2) = -detinv * (A(1, 2) * A(3, 3) - A(1, 3) * A(3, 2))
    B(2, 2) = +detinv * (A(1, 1) * A(3, 3) - A(1, 3) * A(3, 1))
    B(3, 2) = -detinv * (A(1, 1) * A(3, 2) - A(1, 2) * A(3, 1))
    B(1, 3) = +detinv * (A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2))
    B(2, 3) = -detinv * (A(1, 1) * A(2, 3) - A(1, 3) * A(2, 1))
    B(3, 3) = +detinv * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1))
    !-----------------------------------------------------------------------
    END FUNCTION matinv3
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    PURE FUNCTION find_minimum(grid, grid_dim) RESULT(minpos)
    !-----------------------------------------------------------------------
    !!
    !! Return the position of a minimum and its value in a grid
    !!
    USE kinds,         ONLY : DP
    !
    INTEGER, INTENT(in) :: grid_dim
    !! Grid dimension
    REAL(KIND = DP), INTENT(in) :: grid(grid_dim)
    !! Actual grid
    !
    ! Local variable
    INTEGER :: i
    !! Index of the grid
    REAL(KIND = DP) :: minvalore
    !! Minimum value
    INTEGER :: minpos
    !! Return the minimum position
    !
    minvalore = grid(1)
    minpos = 1
    DO i = 2, grid_dim
      IF (grid(i) < minvalore) THEN !supposing only one minimum
        minpos = i
        minvalore = grid(i)
      ENDIF
    ENDDO
    !-----------------------------------------------------------------------
    END FUNCTION find_minimum
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    LOGICAL FUNCTION eqvect_strict(x, y, accep)
    !-----------------------------------------------------------------------
    !!
    !! This function test if two tridimensional vectors are equal
    !!
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: x(3)
    !! input: input vector
    REAL(KIND = DP), INTENT(in) :: y(3)
    !! input: second input vector
    REAL(KIND = DP), INTENT(in) :: accep
    !! acceptance parameter
    !
    eqvect_strict = ABS(x(1) - y(1)) < accep .AND. &
                    ABS(x(2) - y(2)) < accep .AND. &
                    ABS(x(3) - y(3)) < accep
    !
    !----------------------------------------------------------------------
    END FUNCTION eqvect_strict
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE fractrasl(npw, igk, evc, eigv1, eig0v)
    !----------------------------------------------------------------------
    !!
    !! Routine to compute fractional translations
    !!
    USE kinds, ONLY : DP
    USE wvfct, ONLY : nbnd, npwx
    USE noncollin_module, ONLY : noncolin, npol
    USE elph2, ONLY : ngxxf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: npw
    !! Number of plane-waves
    INTEGER, INTENT(in) :: igk(npw)
    !! G mapping
    COMPLEX(KIND = DP), INTENT(inout) :: evc(npwx * npol, nbnd)
    !!
    COMPLEX(KIND = DP), INTENT(in) :: eigv1(ngxxf)
    !! Eigenvalues
    COMPLEX(KIND = DP), INTENT(in) :: eig0v
    !! Eigenvalues
    !
    INTEGER :: ig
    !! Counter on G-vectors
    INTEGER :: ibnd
    !! Counter on bands
    !
    DO ibnd = 1, nbnd
      DO ig = 1, npw
        evc(ig, ibnd) = evc(ig, ibnd) * eigv1(igk(ig)) * eig0v
        IF (noncolin) THEN
          evc(ig + npwx, ibnd) = evc(ig + npwx, ibnd) * eigv1(igk(ig)) * eig0v
        ENDIF
      ENDDO
    ENDDO
    !
    !------------------------------------------------------------
    END SUBROUTINE fractrasl
    !------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rotate_cart(x, s, sx)
    !-----------------------------------------------------------------------
    !!
    !! A simple symmetry operation in cartesian coordinates
    !! ( s is INTEGER and in crystal coord!)
    !!
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: s(3, 3)
    !! Symmetry matrix
    REAL(KIND = DP), INTENT(in) :: x(3)
    !! Input x
    REAL(KIND = DP), INTENT(out) :: sx(3)
    !! Output rotated x
    !
    ! Local variables
    INTEGER :: i
    !! Cartesian direction
    REAL(KIND = DP) :: xcrys(3)
    !! x in cartesian coords
    !
    xcrys = x
    CALL cryst_to_cart(1, xcrys, at, -1)
    DO i = 1, 3
       sx(i) = DBLE(s(i,1)) * xcrys(1) &
             + DBLE(s(i,2)) * xcrys(2) &
             + DBLE(s(i,3)) * xcrys(3)
    ENDDO
    CALL cryst_to_cart(1, sx, bg, +1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE rotate_cart
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE mem_size(nmodes, nkf)
    !--------------------------------------------------------------------------
    !!
    !! This routine estimates the amount of memory taken up by
    !! the $$<k+q| dV_q,nu |k>$$ on the fine meshes and prints
    !! out a useful(?) message
    !!
    USE io_global, ONLY : stdout
    USE kinds,     ONLY : DP
    USE elph2,     ONLY : nbndfst
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Number of modes
    INTEGER, INTENT(in) :: nkf
    !! Number of k-points in pool
    !
    ! Local variables
    CHARACTER(LEN = 256) :: chunit
    !! Unit name
    INTEGER :: imelt
    !! Size in number of elements
    REAL(KIND = DP) :: rmelt
    !! Size in byte
    !
    imelt = (nbndfst**2) * nmodes * nkf
    rmelt = imelt * 8 / 1048576.d0 ! 8 bytes per number, value in Mb
    IF (rmelt < 1000.0) THEN
      chunit =  ' Mb '
      IF (rmelt < 1.0) THEN
        chunit = ' Kb '
        rmelt  = rmelt * 1024.d0
      ENDIF
    ELSE
      rmelt = rmelt / 1024.d0
      chunit = ' Gb '
    ENDIF
    WRITE(stdout, '(/, 5x, a, i13, a, f7.2, a, a)') "Number of ep-matrix elements per pool :", &
         imelt, " ~= ", rmelt, TRIM(chunit), " (@ 8 bytes/ DP)"
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE mem_size
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mem_size_eliashberg(vmelt, imelt)
    !-----------------------------------------------------------------------
    !
    !  This routine estimates the amount of memory taken up or
    !  released by different arrays
    !
    USE io_global,     ONLY : stdout
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : max_memlt
    USE eliashbergcom, ONLY : memlt_pool
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: vmelt
    !! 1 for integer variables and 2 for real variables
    INTEGER, INTENT(in) :: imelt
    !! > 0 memory added or < 0 memory subtracted
    !
    REAL(KIND = DP) :: rmelt
    !! change in memory
    !
    ! This is only a quick fix since the routine was written for parallel
    ! execution - FG June 2014
#if !defined(__MPI)
    my_pool_id = 0
#endif
    !
    rmelt = zero
    rmelt = DBLE(imelt) * 4.d0 / 1073741824.d0 ! 4 bytes per number, value in Gb
    IF (vmelt == 2) &
      rmelt = 2.d0 * rmelt ! 8 bytes per number, value in Gb
    rmelt = rmelt + memlt_pool(my_pool_id + 1)
    !
    memlt_pool(:) = zero
    memlt_pool(my_pool_id + 1) = rmelt
    !
    ! collect contributions from all pools
    CALL mp_sum(memlt_pool, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (MAXVAL(memlt_pool(:)) > max_memlt) THEN
      WRITE(stdout, '(/, 5x, a, a, f9.4, a)') "Size of required memory per pool:", &
            " ~= ", MAXVAL(memlt_pool(:)), " Gb"
      CALL errore('mem_size_eliashberg', 'Size of required memory exceeds max_memlt', 1)
    ELSEIF(MAXVAL(memlt_pool(:)) > 0.5d0 * max_memlt) THEN
      WRITE(stdout, '(/, 5x, a, a, f9.4, a)') "Size of allocated memory per pool:", &
            " ~= ", MAXVAL(memlt_pool(:)), " Gb"
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mem_size_eliashberg
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE memlt_eliashberg(itemp, cname)
    !-----------------------------------------------------------------------
    !!
    !! Estimate the memory requirements for anisotropic Eliashberg equations
    !! on imaginary axis or real axis (analytic continuation)
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : max_memlt, nqstep
    USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, nqfs, limag_fly, &
                              lacon_fly, memlt_pool
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    CHARACTER(LEN = 4), INTENT(in) :: cname
    !! calculation type
    !
    !Local variables
    INTEGER :: imelt
    !! size array
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k parallelization
    !
    REAL(KIND = DP) :: rmelt
    !! change in memory
    !
    ! This is only a quick fix since the routine was written for parallel
    ! execution - FG June 2014
#if !defined(__MPI)
    my_pool_id = 0
#endif
    !
    limag_fly = .FALSE.
    lacon_fly = .FALSE.
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    imelt = (upper_bnd - lower_bnd + 1) * MAXVAL(nqfs(:)) * nbndfs**2
    IF (cname == 'imag') THEN
      ! get the size of the akeri that needa to be stored in each pool
      imelt = imelt * (2 * nsiw(itemp))
    ELSEIF (cname == 'acon') THEN
      ! get the size of a2fij that needs to be stored in each pool
      imelt = imelt * nqstep
    ENDIF
    rmelt = DBLE(imelt) * 8.d0 / 1073741824.d0 ! 8 bytes per number, value in Gb
    rmelt = rmelt + memlt_pool(my_pool_id + 1)
    !
    memlt_pool(:) = zero
    memlt_pool(my_pool_id + 1) = rmelt
    !
    ! collect contributions from all pools
    CALL mp_sum(memlt_pool, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (MAXVAL(memlt_pool(:)) > max_memlt) THEN
      WRITE(stdout, '(/, 5x, a, a, f9.4, a)') "Size of required memory per pool:", &
            " ~= ", MAXVAL(memlt_pool(:)), " Gb"
      IF (cname == 'imag') limag_fly = .TRUE.
      IF (cname == 'acon') lacon_fly = .TRUE.
      !
      ! remove memory required for akeri or a2fij
      CALL mem_size_eliashberg(2, -imelt)
      !
    ENDIF
    !
    IF (limag_fly) THEN
      WRITE(stdout, '(/, 5x, a/)') "akeri is calculated on the fly since its size exceedes max_memlt"
    ELSEIF (lacon_fly) THEN
      WRITE(stdout, '(/, 5x, a/)') "a2fij is calculated on the fly since its size exceedes max_memlt"
    ELSE
      WRITE(stdout, '(/, 5x, a, a, f9.4, a)') "Size of allocated memory per pool:", &
            " ~= ", MAXVAL(memlt_pool(:)), " Gb"
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE memlt_eliashberg
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE create_interval(size_v, v, n_intval, val_intval, pos_intval)
    !--------------------------------------------------------------------------
    !!
    !! This routine divide an integer vector "v" of size "size_v" into "n_intval" intervals
    !! The routine then store the first element of each intervals into "val_intval"
    !! and store their position in "pos_intval".
    !!
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: size_v
    !! Size of vector to bisect
    INTEGER, INTENT(in) :: v(size_v)
    !! Vector to bisect
    INTEGER, INTENT(in) :: n_intval
    !! Number of intervals indexes (so there are (n_intval - 1) number of intervals)
    INTEGER, INTENT(out) :: val_intval(n_intval)
    !! Value of the first element of each intervals
    INTEGER, INTENT(out) :: pos_intval(n_intval)
    !! Position of the first element of each intervals
    !
    !  Local variables
    INTEGER :: nkl
    !! Size of each intervals
    INTEGER :: nkr
    !! Remaining of the division of the full vector by nkl
    INTEGER :: i
    !! Index looping on the invervals
    !
    ! Initialize
    val_intval(:) = 0
    pos_intval(:) = 0
    !
    ! We want to have n_intval indexes, i.e. n_intval-1 intervals
    nkl = size_v / (n_intval - 1)
    nkr = size_v - nkl * (n_intval - 1)
    !
    ! The reminder goes to the first nkr intervals (0...nkr-1)
    !
    DO i = 1, n_intval
      pos_intval(i) = nkl * (i - 1)
      IF (i < nkr) pos_intval(i) = pos_intval(i) + i
      IF (i >= nkr) pos_intval(i) = pos_intval(i) + nkr
    ENDDO
    !
    ! In case the reminder is 0
    IF (nkr == 0) THEN
      DO i = 1, (n_intval - 1)
        pos_intval(i) = pos_intval(i) + 1
      ENDDO
    ENDIF
    !
    DO i = 1, n_intval
      val_intval(i) = v(pos_intval(i))
    ENDDO
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE create_interval
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE bisection(size_v, v, ind, n_intval, val_intval, pos_intval)
    !--------------------------------------------------------------------------
    !!
    !! The subroutine first perform a pre-serach to determine in which intervals
    !! the index "ind" lies.
    !! Then the subroutines find the position of the index "ind" within that interaval.
    !!
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: size_v
    !! Size of vector to bisect
    INTEGER, INTENT(in) :: v(size_v)
    !! Vector to bisect
    INTEGER, INTENT(in) :: n_intval
    !! Number of intervals indexes (so there are (n_intval - 1) number of intervals
    INTEGER, INTENT(in) :: val_intval(n_intval)
    !! Value of the first element of each intervals
    INTEGER, INTENT(in) :: pos_intval(n_intval)
    !! Position of the first element of each intervals
    INTEGER, INTENT(inout) :: ind
    !! In entry, index to find in the vector; on exit, position of the element;
    !! If 0 on RETURN, the element is not present in the list
    !
    ! Local variable
    INTEGER :: it
    !! Iteration index
    INTEGER :: lp
    !! Left pointer
    INTEGER :: rp
    !! Right pointer
    INTEGER :: cp
    !! Center pointer
    INTEGER :: nit
    !! Number of iterations
    INTEGER :: prod
    !! Product to determine on which side of the interval we are
    INTEGER :: prod_2
    !! Product
    INTEGER :: ierr
    !! Error number
    INTEGER :: subsize
    !! Size of the chunk of v we are considering
    INTEGER :: pre_minloc
    !! Minimum location of presearch
    INTEGER :: pre_minval
    !! Minimum value of presearch
    INTEGER, ALLOCATABLE :: v_red(:)
    !! Reduced vector for bisection
    INTEGER, ALLOCATABLE :: v_sgn(:)
    !! v - ind
    !
    ! 1) Deal with special cases
    !
    ! The point is not within the fsthick window ==> return 0
    IF (ind > v(size_v) .OR. (ind < v(1))) THEN
      ind = 0
      RETURN
    ENDIF
    !
    ! Look at the first element of each intervals.
    pre_minloc = MINLOC(ABS(val_intval(:) - ind), 1)
    pre_minval = MINVAL(ABS(val_intval(:) - ind), 1)
    !
    ! The ind is one of those first element (special case) - we are done
    IF (pre_minval == 0) THEN
      ind = pos_intval(pre_minloc)
      RETURN
    ENDIF
    !
    ! The ind is the second or one before last element (special case treated here because the
    ! algo below cannot deal with the two extreme cases.
    ! We compute the lower left bound and upper right bound. The ind is between those two bounds.
    IF (pre_minloc == 1) THEN
      lp = pos_intval(1)
      rp = pos_intval(2)
    ELSEIF (pre_minloc == n_intval) THEN
      lp = pos_intval(n_intval - 1)
      rp = pos_intval(n_intval)
    ELSE
      !
      ! 2) Do a pre-search to deterine in which intervals the index lies.
      !
      ! The sign of the product tels us on which side of the interval does the index ind lies.
      IF ( val_intval(pre_minloc) - ind /= 0 .AND. val_intval(pre_minloc - 1) - ind /= 0) THEN
        prod = (val_intval(pre_minloc) - ind) / ABS(val_intval(pre_minloc) - ind) &
             * (val_intval(pre_minloc - 1) - ind) / ABS(val_intval(pre_minloc - 1) - ind)
      ELSE
        prod = 0
      ENDIF
      !
      IF (prod < 0) THEN
        ! Then the correspondence is in between pos_intval(pre_minloc-1) and pos_intval(pre_minloc)
        lp = pos_intval(pre_minloc - 1)
        rp = pos_intval(pre_minloc)
      ELSEIF (prod > 0) THEN
        ! Then the correspondence is in between pos_intval(pre_minloc) and pos_intval(pre_minloc+1)
        lp = pos_intval(pre_minloc)
        rp = pos_intval(pre_minloc + 1)
      ELSEIF (prod == 0) THEN
        ! It cannot be that val_intval(pre_minloc)-ind) == 0 because we already checked
        ind = pos_intval(pre_minloc - 1)
        RETURN
      ENDIF
    ENDIF
    !
    subsize = rp - lp + 1
    ALLOCATE(v_red(subsize), STAT = ierr)
    IF (ierr /= 0) CALL errore('bisection', 'Error allocating v_red', 1)
    ALLOCATE(v_sgn(subsize), STAT = ierr)
    IF (ierr /= 0) CALL errore('bisection', 'Error allocating v_sgn', 1)
    !
    !Shift lp to 1 and rp to subsize
    IF (pre_minloc == 1) THEN
      v_red(:) = v(pos_intval(1):pos_intval(2))
      lp = lp - pos_intval(1) + 1
      rp = rp - pos_intval(1) + 1
    ELSEIF (pre_minloc == n_intval) THEN
      v_red(:) = v(pos_intval(n_intval - 1):pos_intval(n_intval))
      lp = lp - pos_intval(n_intval - 1) + 1
      rp = rp - pos_intval(n_intval - 1) + 1
    ELSE
      IF (prod < 0) THEN
        v_red(:) = v(pos_intval(pre_minloc - 1):pos_intval(pre_minloc))
        lp = lp - pos_intval(pre_minloc - 1) + 1
        rp = rp - pos_intval(pre_minloc - 1) + 1
      ELSEIF (prod > 0) THEN
        v_red(:) = v(pos_intval(pre_minloc):pos_intval(pre_minloc + 1))
        lp = lp - pos_intval(pre_minloc) + 1
        rp = rp - pos_intval(pre_minloc) + 1
      ENDIF
    ENDIF
    !
    ! 3) Do a bisection on the interval
    !
    v_sgn(:) = v_red(:) - ind
    !
    ! Each iteration of the bisection method divide the range by 2 so that the
    ! maximum number of iteration is log2(subsize). We multiply by 2 for safety.
    nit = (NINT(LOG(REAL(subsize, KIND = DP)) / LOG(2.0d0)) + 1) * 2 !Safe range
    ! Add call errore if it == nit
    !
    bisec : DO it = 1, nit
      !
      IF (it == nit) CALL errore('bisection', 'Maximum number of iteration reached in bisection', 1)
      !
      ! Check if we found a correspondence
      IF (v_sgn(lp) == 0) THEN
        ind = lp
        IF (pre_minloc == 1) THEN
          ind = ind + pos_intval(1) - 1
          EXIT bisec
        ELSEIF (pre_minloc == n_intval) THEN
          ind = ind + pos_intval(n_intval - 1) - 1
          EXIT bisec
        ELSE
          IF (prod < 0) THEN
            ind = ind + pos_intval(pre_minloc - 1) - 1
          ELSEIF (prod > 0) THEN
            ind = ind + pos_intval(pre_minloc) - 1
          ENDIF
          EXIT bisec
        ENDIF
      ELSEIF (v_sgn(rp) == 0) THEN
        ind = rp
        IF (pre_minloc == 1) THEN
          ind = ind + pos_intval(1) - 1
          EXIT bisec
        ELSEIF (pre_minloc == n_intval) THEN
          ind = ind + pos_intval(n_intval - 1) - 1
          EXIT bisec
        ELSE
          IF (prod < 0) THEN
            ind = ind + pos_intval(pre_minloc - 1) - 1
          ELSEIF (prod > 0) THEN
            ind = ind + pos_intval(pre_minloc) - 1
          ENDIF
          EXIT bisec
        ENDIF
      ENDIF
      ! v_sgn(lp) and v_sgn(rp) cannot be 0 at this point
      prod_2 = (v_sgn(lp) / ABS(v_sgn(lp))) * (v_sgn(rp) / ABS(v_sgn(rp)))
      !
      ! If prod not 0 but they are contiguous, ind = 0
      IF (lp + 1 == rp) THEN
        ind = 0
        EXIT bisec
      ENDIF
      ! Check for the case prod > 0
      IF (prod_2 > 0) THEN
        CALL errore('bisection', 'Error in bisection algorithm', 1)
      ENDIF
      ! After the checks, move lp or rp
      IF (MOD(subsize, 2) == 0) THEN
        cp = lp + subsize / 2 - 1
      ELSE
        cp = lp + (subsize + 1) / 2 - 1
      ENDIF
      IF (v_sgn(cp) == 0) THEN
        ind = cp
        IF (pre_minloc == 1) THEN
          ind = ind + pos_intval(1) - 1
          EXIT bisec
        ELSEIF (pre_minloc == n_intval) THEN
          ind = ind + pos_intval(n_intval - 1) - 1
          EXIT bisec
        ELSE
          IF (prod < 0) THEN
            ind = ind + pos_intval(pre_minloc - 1) - 1
          ELSEIF (prod > 0) THEN
            ind = ind + pos_intval(pre_minloc) - 1
          ENDIF
          EXIT bisec
        ENDIF
      ELSE
        ! v_sgn(lp) and v_sgn(cp) cannot be 0 at this point
        prod_2 = (v_sgn(lp) / ABS(v_sgn(lp))) * (v_sgn(cp) / ABS(v_sgn(cp)))
      ENDIF
      IF (prod_2 < 0) THEN
        rp = cp
      ELSE
        lp = cp
      ENDIF
      subsize = rp - lp + 1
      IF (rp < lp) CALL errore('bisection', 'Error in bisection algorithm', 1)
    ENDDO bisec
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE bisection
    !--------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE s_crystocart(s, sr, at, bg)
    !----------------------------------------------------------------------
    !!
    !! This routine transform a symmetry matrix expressed in the
    !! basis of the crystal axis in the cartesian basis.
    !!
    !! SP - Feb 2020
    !! Routine taken from PP/src/sym_band.f90 and adapted for EPW.
    !!
    USE kinds,    ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: s(3, 3)
    !! Matrix in crystal axis
    REAL(KIND = DP), INTENT(in) :: at(3, 3)
    !! Direct lattice vectors
    REAL(KIND = DP), INTENT(in) :: bg(3, 3)
    !! Reciprocal lattice vectors
    REAL(KIND = DP), INTENT(out) :: sr(3, 3)
    ! Output matrix in cartesian axis
    !
    ! Local variables
    REAL(KIND = DP) :: sa(3, 3)
    !! Temporary matrix
    REAL(KIND = DP) :: sb(3, 3)
    !! Temporary matrix
    !
    sa(:, :) = DBLE(s(:, :))
    sb = MATMUL(bg, sa)
    sr(:, :) = MATMUL(at, TRANSPOSE(sb))
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE s_crystocart
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    INTEGER FUNCTION copy_sym_epw(nrot_, sym, indsym)
    !-----------------------------------------------------------------------
    !!
    !! Imported and adapted from copy_sym in symm_base.f90 in QE
    !!
    USE kinds,     ONLY : DP
    USE symm_base, ONLY : irt, s, ft, sname, t_rev
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nrot_
    !! number of rotations
    INTEGER, INTENT(inout) :: indsym(48)
    !! mapping between the original sym. indices and the new indices
    LOGICAL, INTENT(inout) :: sym(48)
    !! .TRUE. if rotation isym is a sym.op. of the crystal
    !! (i.e. not of the bravais lattice only)
    !
    ! Local variables
    CHARACTER(LEN = 45) :: nametemp
    !! Temporary name of the rotation
    INTEGER :: stemp(3, 3)
    !! Temporary s matrix for that rotation
    INTEGER :: ttemp
    !! Temporary time-reversal for that rotation
    INTEGER :: ierr
    !! Error index
    INTEGER :: irot
    !! Rotation index
    INTEGER :: jrot
    !! Rotation index
    INTEGER, ALLOCATABLE :: irtemp(:)
    !! Temporary irt.
    REAL(KIND = DP) :: ft_(3)
    !! Fractional translation
    !
    ALLOCATE(irtemp(SIZE(irt, 2)), STAT = ierr)
    IF (ierr /= 0) CALL errore('copy_sym_epw', 'Error allocating irtemp', 1)
    !
    jrot = 0
    !
    DO irot = 1, nrot_
      IF (sym(irot)) THEN
        jrot = jrot + 1
        IF (irot > jrot) THEN
          stemp = s(:, :, jrot)
          s(:, :, jrot) = s(:, :, irot)
          s(:, :, irot) = stemp
          ft_(:) = ft(:, jrot)
          ft(:, jrot) = ft(:, irot)
          ft(:, irot) = ft_(:)
          irtemp(:) = irt(jrot, :)
          irt(jrot, :) = irt(irot, :)
          irt(irot, :) = irtemp(:)
          nametemp = sname(jrot)
          sname(jrot) = sname(irot)
          sname(irot) = nametemp
          ttemp = t_rev(jrot)
          t_rev(jrot) = t_rev(irot)
          t_rev(irot) = ttemp
          ttemp = indsym(jrot)
          indsym(jrot) = indsym(irot)
          indsym(irot) = ttemp
        ENDIF
      ENDIF
    ENDDO
    !
    sym(1:jrot) = .TRUE.
    sym(jrot + 1:nrot_) = .FALSE.
    !
    DEALLOCATE(irtemp, STAT = ierr)
    IF (ierr /= 0) CALL errore('copy_sym_epw', 'Error deallocating irtemp', 1)
    !
    copy_sym_epw = jrot
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END FUNCTION copy_sym_epw
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE fix_sym(lcase)
    !----------------------------------------------------------------------
    !!
    !! This routine tries to fix some symmetry-related issues in EPW.
    !! This is just a temporary fix before restructuring the code.
    !!
    !! HL - July 2020
    !!
    USE kinds,         ONLY : DP
    USE symm_base,     ONLY : t_rev, time_reversal, nrot, nsym, invsym, &
                              nofrac, nosym_evc, find_sym
    USE ions_base,     ONLY : nat, tau, ityp
    USE epwcom,        ONLY : epw_no_t_rev, epw_tr, epw_nosym, epw_crysym, &
                              mp_mesh_k, mp_mesh_q
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lcase
    !! .FALSE.: before interpolation
    !! .TRUE.: after interpolation
    !
    REAL(DP) :: mdum(3,nat)
    !! Dummy values for local magnitization
    !
    IF (epw_no_t_rev) t_rev = 0
    !
    time_reversal = epw_tr
    !
    IF (lcase) THEN
      IF (epw_crysym) THEN
        IF (mp_mesh_k .OR. mp_mesh_q) THEN
          !
          nosym_evc = .FALSE.
          nofrac = .FALSE.
          !
          CALL find_sym(nat, tau, ityp, .FALSE., mdum)
          WRITE(stdout, '(a,i3)') 'fix_sym: nrot=', nrot
          WRITE(stdout, '(a,i3)') 'fix_sym: nsym=', nsym
          nrot = nsym
          !
        ENDIF
      ENDIF
    ENDIF
    !
    IF (epw_nosym) THEN
      nrot = 1
      nsym = 1
      invsym = .FALSE.
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE fix_sym
    !-----------------------------------------------------------------------
  !-------------------------------------------------------------------------
  END MODULE low_lvl
  !-------------------------------------------------------------------------
