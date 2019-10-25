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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi
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
    !-----------------------------------------------------------------------
    SUBROUTINE fermiwindow() 
    !-----------------------------------------------------------------------
    !!
    !! Find the band indices of the first
    !! and last state falling within the window e_fermi+-efermithickness
    !! 
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf
    USE epwcom,        ONLY : fsthick, nbndsub
    USE pwcom,         ONLY : ef
    USE mp,            ONLY : mp_max, mp_min
    USE mp_global,     ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER :: ik
    !! Counter on k-points in the pool
    INTEGER :: ibnd
    !! Counter on bands
    REAL(KIND = DP) :: ebnd
    !! Eigenvalue at etf(ibnd, ik)
    REAL(KIND = DP) :: ebndmin
    !! Minimum eigenvalue
    REAL(KIND = DP) :: ebndmax
    !! Maximum eigenvalue
    REAL(KIND = DP) :: tmp
    !
    !
    ibndmin = 100000
    ibndmax = 0
    ebndmin =  1.d8
    ebndmax = -1.d8
    !
    DO ik = 1, nkqf
      DO ibnd = 1, nbndsub
        ebnd = etf(ibnd, ik)
        !
        IF (ABS(ebnd - ef) < fsthick) THEN
          ibndmin = MIN(ibnd, ibndmin)
          ibndmax = MAX(ibnd, ibndmax)
          ebndmin = MIN(ebnd, ebndmin)
          ebndmax = MAX(ebnd, ebndmax)
        ENDIF
        !
      ENDDO
    ENDDO
    !
    tmp = DBLE(ibndmin)
    CALL mp_min(tmp, inter_pool_comm)
    ibndmin = NINT(tmp)
    CALL mp_min(ebndmin, inter_pool_comm)
    !
    tmp = DBLE(ibndmax)
    CALL mp_max(tmp, inter_pool_comm)
    ibndmax = NINT(tmp)
    CALL mp_max(ebndmax, inter_pool_comm)
    !
    WRITE(stdout,'(/14x,a,i5,2x,a,f9.3)') 'ibndmin = ', ibndmin, 'ebndmin = ', ebndmin
    WRITE(stdout,'(14x,a,i5,2x,a,f9.3/)') 'ibndmax = ', ibndmax, 'ebndmax = ', ebndmax
    !
    !----------------------------------------------------------------------
    END SUBROUTINE fermiwindow
    !---------------------------------------------------------------------
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
    CHARACTER(LEN = 3), INTENT(out) :: ndlab 
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
    ndlab = '   '
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
        WRITE(ndlab, '(i4)') node
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
#if defined(__PGI) || defined(__CRAY) || defined(__XLF)
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
    PURE FUNCTION matinv3(A) RESULT(B)
    !-----------------------------------------------------------------------
    !!
    !! Performs a direct calculation of the inverse of a 3×3 matrix. 
    !! 
    USE kinds, ONLY : DP
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
    !---------------------------------------------------------------------------
    SUBROUTINE read_modes(iunpun, current_iq, ierr)
    !---------------------------------------------------------------------------
    !!
    !! This routine reads the displacement patterns.
    !!
    USE modes,        ONLY : nirr, npert, u
    USE lr_symm_base, ONLY : minus_q, nsymq
    USE iotk_module,  ONLY : iotk_index, iotk_scan_dat, iotk_scan_begin, &
                             iotk_scan_end
    USE io_global,    ONLY : meta_ionode, meta_ionode_id
    USE mp,           ONLY : mp_bcast
    USE mp_global,    ONLY : world_comm
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: current_iq
    !! Current q-point 
    INTEGER, INTENT(in) :: iunpun
    !! Current q-point 
    INTEGER, INTENT(out) :: ierr
    !! Error
    !
    ! Local variables
    INTEGER :: imode0, imode
    !! Counter on modes
    INTEGER :: irr
    !! Counter on irreducible representations
    INTEGER :: ipert
    !! Counter on perturbations at each irr
    INTEGER :: iq
    !! Current q-point 
    !
    ierr = 0
    IF (meta_ionode) THEN
      CALL iotk_scan_begin(iunpun, "IRREPS_INFO")
      !
      CALL iotk_scan_dat(iunpun, "QPOINT_NUMBER", iq)
    ENDIF
    CALL mp_bcast(iq,  meta_ionode_id, world_comm)
    IF (iq /= current_iq) CALL errore('read_modes', ' Problems with current_iq', 1)
    ! 
    IF (meta_ionode) THEN
      !
      CALL iotk_scan_dat(iunpun, "QPOINT_GROUP_RANK", nsymq)
      CALL iotk_scan_dat(iunpun, "MINUS_Q_SYM", minus_q)
      CALL iotk_scan_dat(iunpun, "NUMBER_IRR_REP", nirr)
      imode0 = 0
      DO irr = 1, nirr
        CALL iotk_scan_begin(iunpun, "REPRESENTION" // TRIM(iotk_index(irr)))
        CALL iotk_scan_dat(iunpun, "NUMBER_OF_PERTURBATIONS", npert(irr))
        DO ipert = 1, npert(irr)
          imode = imode0 + ipert
          CALL iotk_scan_begin(iunpun, "PERTURBATION" // TRIM(iotk_index(ipert)))
          CALL iotk_scan_dat(iunpun, "DISPLACEMENT_PATTERN", u(:, imode))
          CALL iotk_scan_end(iunpun, "PERTURBATION" // TRIM(iotk_index(ipert)))
        ENDDO
        imode0 = imode0 + npert(irr)
        CALL iotk_scan_end(iunpun, "REPRESENTION" // TRIM(iotk_index(irr)))
      ENDDO
      !
      CALL iotk_scan_end(iunpun, "IRREPS_INFO")
      !
    ENDIF
    !
    CALL mp_bcast(nirr   , meta_ionode_id, world_comm)
    CALL mp_bcast(npert  , meta_ionode_id, world_comm)
    CALL mp_bcast(nsymq  , meta_ionode_id, world_comm)
    CALL mp_bcast(minus_q, meta_ionode_id, world_comm)
    CALL mp_bcast(u      , meta_ionode_id, world_comm)
    !
    RETURN
    ! 
    !---------------------------------------------------------------------------
    END SUBROUTINE read_modes
    !---------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------
    SUBROUTINE fractrasl(npw, igk, evc, eigv1, eig0v)
    !------------------------------------------------------------
    !!
    !! Routine to compute fractional translations
    !! 
    USE kinds, ONLY : DP
    USE wvfct, ONLY : nbnd, npwx
    USE gvect, ONLY : ngm
    USE noncollin_module, ONLY : noncolin, npol
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: npw
    !! Number of plane-waves
    INTEGER, INTENT(in) :: igk(npw)
    !! G mapping
    COMPLEX(KIND = DP), INTENT(inout) :: evc(npwx * npol, nbnd)
    !! 
    COMPLEX(KIND = DP), INTENT(in) :: eigv1(ngm)
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
    USE kinds, ONLY : DP
    USE cell_base, ONLY : at, bg
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: x(3)
    !! Input x
    INTEGER, INTENT(in) :: s(3,3)
    !! Symmetry matrix
    REAL(KIND = DP), INTENT(out) :: sx(3)
    !! Output rotated x
    !
    REAL(KIND = DP) :: xcrys(3)
    !! x in cartesian coords
    INTEGER :: i
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
    !-----------------------------------------------------------------------
    SUBROUTINE compute_dos(itemp, ef0, dos)
    !-----------------------------------------------------------------------
    !!
    !! This routine computes the density of states at a given fermi level.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : two, eps16, ryd2mev
    USE epwcom,        ONLY : ngaussw, nstemp, nbndsub, degaussw
    USE elph2,         ONLY : etf, nkqf, wkf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(inout) :: dos(nstemp)
    !! DOS to compute for the temperature itemp.
    !
    ! Local variables
    REAL(KIND = DP), EXTERNAL :: dos_ef
    ! 
    ! divide by two to have DOS/spin
    IF (ABS(degaussw) < eps16) THEN
      ! use 1 meV instead
      dos(itemp) = dos_ef(ngaussw, 1.0d0 / ryd2mev, ef0(itemp), etf, wkf, nkqf, nbndsub) / two
    ELSE
      dos(itemp) = dos_ef(ngaussw, degaussw, ef0(itemp), etf, wkf, nkqf, nbndsub) / two
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_dos
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE fermicarrier(itemp, etemp, ef0, efcb, ctype)
    !-----------------------------------------------------------------------
    !!
    !!  This routine computes the Fermi energy associated with a given 
    !!  carrier concentration using bissection for insulators or
    !!  semi-conductors.
    !!
    !-----------------------------------------------------------------------
    USE cell_base, ONLY : omega, alat, at
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : etf, nkf, wkf, efnew, nkqf
    USE constants_epw, ONLY : ryd2ev, bohr2ang, ang2cm, eps5, kelvin2eV, zero, eps80
    USE noncollin_module, ONLY : noncolin
    USE pwcom,     ONLY : nelec
    USE epwcom,    ONLY : int_mob, nbndsub, ncarrier, nstemp, fermi_energy, &
                          system_2d, carrier, efermi_read, assume_metal, ngaussw
    USE klist_epw, ONLY : isk_dummy
    USE mp,        ONLY : mp_barrier, mp_sum, mp_max, mp_min
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    INTEGER, INTENT(out) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in kBT [Ry] unit.
    REAL(KIND = DP), INTENT(inout) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(inout) :: efcb(nstemp)
    !! Second fermi level for the temperature itemp
    REAL(KIND = DP), EXTERNAL :: efermig
    !! External function to calculate the fermi energy
    ! 
    ! Local variables 
    INTEGER :: i
    !! Index for the bisection iteration
    INTEGER :: ik
    !! k-point index per pool
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ivbm
    !! Index of the VBM
    INTEGER :: icbm
    !! Index of the CBM
    INTEGER, PARAMETER :: maxiter = 500 ! 300
    !! Maximum interation
    REAL(KIND = DP) :: fermi
    !! Fermi level returned
    REAL(KIND = DP) :: fermicb
    !! Fermi level returned for second Fermi level
    REAL(KIND = DP) :: fnk
    !! Fermi-Diract occupation
    REAL(KIND = DP) :: ks_exp(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT
    REAL(KIND = DP) :: ks_expcb(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT for CB
    REAL(KIND = DP) :: fermi_exp
    !! Fermi level in exponential format
    REAL(KIND = DP) :: rel_err
    !! Relative error
    REAL(KIND = DP) :: factor
    !! Factor that goes from number of carrier per unit cell to number of
    !! carrier per cm^-3
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: evbm
    !! Energy of the VBM
    REAL(KIND = DP) :: ecbm
    !! Energy of the CBM
    REAL(KIND = DP) :: Ef
    !! Energy of the current Fermi level for the bisection method
    REAL(KIND = DP) :: Elw
    !! Energy lower bound for the bisection method
    REAL(KIND = DP) :: Eup
    !! Energy upper bound for the bisection method
    REAL(KIND = DP) :: hole_density
    !! Hole carrier density
    REAL(KIND = DP) :: electron_density 
    !! Electron carrier density
    REAL(KIND = DP), PARAMETER :: maxarg = 200.d0
    !! Maximum value for the argument of the exponential
    !
    IF (assume_metal) THEN
      !! set conduction band chemical potential to 0 since it is irrelevent
      ctype = -1  ! act like it's for holes
      efcb(itemp) = 0.0
      ef0(itemp) = efermig(etf, nbndsub, nkqf, nelec, wkf, etemp, ngaussw, 0, isk_dummy)
      RETURN
    ENDIF
    Ef      = zero
    fermi   = zero
    fermicb = zero
    inv_cell = 1.0d0 / omega
    ! 
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF (system_2d) inv_cell = inv_cell * at(3, 3) * alat
    ! vbm index
    IF (noncolin) THEN
      ivbm = FLOOR(nelec / 1.0d0)
    ELSE
      ivbm = FLOOR(nelec / 2.0d0)
    ENDIF  
    icbm = ivbm + 1 ! Nb of bands
    !
    ! Initialization value. Should be large enough ...
    evbm = -10000
    ecbm = 10000 ! In Ry
    ! 
    DO ik = 1, nkf
      ikk = 2 * ik - 1
      DO ibnd = 1, nbndsub
        IF (ibnd < ivbm + 1) THEN
          IF (etf(ibnd, ikk) > evbm) THEN
            evbm = etf (ibnd, ikk)
          ENDIF
        ENDIF
        ! Find cbm index 
        IF (ibnd > ivbm) THEN
          IF (etf(ibnd, ikk) < ecbm) THEN
            ecbm = etf (ibnd, ikk)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !
    ! Find max and min across pools
    !         
    CALL mp_max(evbm, inter_pool_comm)
    CALL mp_min(ecbm, inter_pool_comm)
    !    
    IF (itemp == 1) THEN
      WRITE(stdout, '(5x,"Valence band maximum    = ",f10.6," eV")') evbm * ryd2ev
      WRITE(stdout, '(5x,"Conduction band minimum = ",f10.6," eV")') ecbm * ryd2ev
    ENDIF
    ! 
    ! Store e^(e_nk/kbT) on each core
    DO ik = 1, nkf
      DO ibnd = 1, nbndsub
        ikk = 2 * ik - 1
        ! Because the number are so large. It does lead to instabilities
        ! Therefore we rescale everything to the VBM
        IF (ABS(etemp) < eps80) THEN
          CALL errore('fermicarrier', 'etemp cannot be 0', 1)
        ELSE
          arg = (etf(ibnd, ikk) - evbm) / etemp 
        ENDIF
        !
        IF (arg < - maxarg) THEN
          ks_exp(ibnd, ik) = 0.0d0
        ELSE
          ks_exp(ibnd, ik) = EXP(arg)
        ENDIF
      ENDDO
    ENDDO
    !
    ! Store e^(e_nk/kbT) on each core for the electrons (CBM only)
    DO ik = 1, nkf
      DO ibnd = 1, nbndsub
        ikk = 2 * ik - 1
        ! Because the number are so large. It does lead to instabilities
        ! Therefore we rescale everything to the CBM
        arg = (etf(ibnd, ikk) - ecbm) / etemp
        !
        IF (arg > maxarg) THEN
          ks_expcb(ibnd, ik) = 1.0d200
        ELSE
          ks_expcb(ibnd, ik) = EXP(arg)
        ENDIF
      ENDDO
    ENDDO
    !
    ! Case 1 : Intrinsic mobilities (electron and hole concentration are the same)   
    ! Starting bounds energy for the biscection method. The energies are rescaled to the VBM
    Elw = 1.0d0  ! This is e^0 = 1.0 
    Eup = 1d-160 ! This is e^(-large) = 0.0 (small)
    IF (int_mob .AND. .NOT. carrier) THEN
      ! Use bisection method
      DO i = 1, maxiter
        !
        !WRITE(stdout,*),'Iteration ',i
        ! We want Ef = (Eup + Elw) / 2.d0 but the variables are exp therefore:
        Ef = DSQRT(Eup) * DSQRT(Elw)     
        ! 
        !WRITE(stdout,*),'Ef ', - log (Ef) * etemp * ryd2ev
        hole_density = 0.0
        electron_density = 0.0
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ! Compute hole carrier concentration
          DO ibnd = 1, ivbm 
            ! Discard very large numbers
            IF (ks_exp(ibnd, ik) * Ef > 1d60) THEN
              fnk = 0.0d0
            ELSE
              fnk = 1.0d0 / (ks_exp(ibnd, ik) * Ef  + 1.0d0)  
            ENDIF
            ! The wkf(ikk) already include a factor 2
            hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk)
          ENDDO
          ! Compute electron carrier concentration
          DO ibnd = icbm, nbndsub            
            ! Discard very large numbers
            IF (ks_exp(ibnd, ik) * Ef > 1d60) THEN
              fnk = 0.0d0
            ELSE
              fnk = 1.0d0 / (ks_exp(ibnd, ik) * Ef  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            electron_density = electron_density + wkf(ikk) * fnk
          ENDDO    
          ! 
        ENDDO 
        !
        CALL mp_sum(hole_density, inter_pool_comm)
        CALL mp_sum(electron_density, inter_pool_comm)
        ! 
        ! WRITE(stdout,*),'hole_density ',hole_density * (1.0d0/omega) * ( bohr2ang * ang2cm  )**(-3)
        ! WRITE(stdout,*),'electron_density ',electron_density * (1.0d0/omega) * (bohr2ang * ang2cm  )**(-3)
        ! CALL FLUSH(stdout)
        IF (ABS(hole_density) < eps80) THEN 
          rel_err = -1000d0
        ELSE
          rel_err = (hole_density - electron_density) / hole_density
        ENDIF
        !
        IF (ABS(rel_err) < eps5) THEN
          fermi_exp = Ef
          fermi = evbm - (LOG(fermi_exp) * etemp)
          EXIT
        ELSEIF ((rel_err) > eps5) THEN
          Elw = Ef                           
        ELSE                                   
          Eup = Ef 
        ENDIF
      ENDDO ! iteration
    ENDIF 
    ! 
    ! Case 2 :
    ! Hole doped mobilities (Carrier concentration should be larger than 1E5 cm^-3)   
    factor = inv_cell * (bohr2ang * ang2cm)**(-3)
    Eup = 1d-160 ! e^(-large) = 0.0 (small)
    Elw = 1.0d0 ! e^0 = 1
    IF (ncarrier < -1E5 .OR. (int_mob .AND. carrier)) THEN
      IF (int_mob .AND. carrier) ncarrier = - ABS(ncarrier)
      ! Use bisection method
      DO i = 1, maxiter
        ! We want Ef = (Eup + Elw) / 2.d0 but the variables are exp therefore:
        Ef = DSQRT(Eup) * DSQRT(Elw)
        ! 
        hole_density = 0.0
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ! Compute hole carrier concentration
          DO ibnd = 1, ivbm
            ! Discard very large numbers
            IF (ks_exp(ibnd, ik) * Ef > 1d60) THEN
              fnk = 0.0d0
            ELSE
              fnk = 1.0d0 / (ks_exp(ibnd, ik) * Ef  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
          ENDDO
          ! 
        ENDDO
        !
        CALL mp_sum(hole_density, inter_pool_comm)
        !
        ! WRITE(stdout,*),'hole_density ',hole_density * (1.0d0/omega) * ( bohr2ang * ang2cm  )**(-3)
        ! CALL FLUSH(stdout)
        IF (ABS(hole_density) < eps80) THEN
          rel_err = -1000.0d0
        ELSE
          ! In this case ncarrier is a negative number
          rel_err = (hole_density - ABS(ncarrier)) / hole_density
        ENDIF
        !
        IF (ABS(rel_err) < eps5) THEN
          fermi_exp = Ef
          fermi = evbm - (LOG(fermi_exp) * etemp)
          EXIT
        ELSEIF ((rel_err) > eps5) THEN
          Elw = Ef
        ELSE
          Eup = Ef
        ENDIF
      ENDDO ! iteration
    ENDIF
    ! 
    ! Case 3 : Electron doped mobilities (Carrier concentration should be larger than 1E5 cm^-3)   
    Eup = 1.0d0 ! e^(0) =1
    Elw = 1.0d80 ! e^large yields fnk = 1
    IF (ncarrier > 1E5 .OR. (int_mob .AND. carrier)) THEN
      IF (int_mob .AND. carrier) ncarrier = ABS(ncarrier)
      ! Use bisection method
      DO i = 1, maxiter
        ! We want Ef = (Eup + Elw) / 2.d0 but the variables are exp therefore:
        Ef = DSQRT(Eup) * DSQRT(Elw)
        ! 
        electron_density = 0.0
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ! Compute electron carrier concentration
          DO ibnd = icbm, nbndsub
            ! Discard very large numbers
            IF (ks_expcb(ibnd, ik) * Ef > 1d60) THEN
              fnk = 0.0d0
            ELSE
              fnk = 1.0d0 / (ks_expcb(ibnd, ik) * Ef  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            electron_density = electron_density + wkf(ikk) * fnk * factor
          ENDDO
          ! 
        ENDDO
        !
        CALL mp_sum(electron_density, inter_pool_comm)
        ! 
        IF (ABS(electron_density) < eps80) THEN
          rel_err = 1000.0d0
        ELSE
          ! In this case ncarrier is a negative number
          rel_err = (electron_density - ncarrier) / electron_density
        ENDIF
        !
        IF (ABS(rel_err) < eps5) THEN
          fermi_exp = Ef
          fermicb = ecbm - (LOG(fermi_exp) * etemp)
          EXIT
        ELSEIF ((rel_err) > eps5) THEN
          Eup = Ef
        ELSE
          Elw = Ef
        ENDIF
      ENDDO ! iteration
    ENDIF
    ! 
    IF (i == maxiter) THEN
      WRITE(stdout, '(5x,"Warning: too many iterations in bisection"/ &
           &      5x,"Ef = ",f10.6)' ) fermi * ryd2ev
    ENDIF
    ! 
    ! Print results
    !  
    WRITE(stdout, '(/5x,"Temperature ",f8.3," K")' ) etemp * ryd2ev / kelvin2eV
    ! 
    ! Small gap semiconductor. Computes intrinsic mobility by placing 
    ! the Fermi level such that carrier density is equal for electron and holes
    IF (int_mob .AND. .NOT. carrier) THEN
      !
      ef0(itemp) = fermi
      WRITE(stdout, '(5x,"Mobility Fermi level = ",f10.6," eV")' )  ef0(itemp) * ryd2ev
      ! We only compute 1 Fermi level so we do not need the other
      efcb(itemp) = 0
      ctype = -1
      !   
    ENDIF
    ! 
    ! Large bandgap semiconductor. Place the gap at the value ncarrier.
    ! The user want both VB and CB mobilities. 
    IF (int_mob .AND. carrier) THEN
      ! 
      ef0(itemp) = fermi
      WRITE(stdout, '(5x,"Mobility VB Fermi level = ",f10.6," eV")' )  ef0(itemp) * ryd2ev
      ! 
      efcb(itemp) = fermicb
      WRITE(stdout, '(5x,"Mobility CB Fermi level = ",f10.6," eV")' )  efcb(itemp) * ryd2ev
      ctype = 0
      !  
    ENDIF
    ! 
    ! User decide the carrier concentration and choose to only look at VB or CB  
    IF (.NOT. int_mob .AND. carrier) THEN
      ! 
      ! VB only
      IF (ncarrier < 0.0) THEN
        ef0(itemp) = fermi
        WRITE(stdout, '(5x,"Mobility VB Fermi level = ",f10.6," eV")' )  ef0(itemp) * ryd2ev
        ! We only compute 1 Fermi level so we do not need the other
        efcb(itemp) = 0
        ctype = -1
      ELSE ! CB 
        efcb(itemp) = fermicb
        WRITE(stdout, '(5x,"Mobility CB Fermi level = ",f10.6," eV")' )  efcb(itemp) * ryd2ev
        ! We only compute 1 Fermi level so we do not need the other
        ef0(itemp) = 0
        ctype = 1
      ENDIF
    ENDIF
    ! 
    ! In the case were we do not want mobility (just scattering rates)
    IF (.NOT. int_mob .AND. .NOT. carrier) THEN
      IF (efermi_read) THEN
        !
        ef0(itemp) = fermi_energy
        !
      ELSE !SP: This is added for efficiency reason because the efermig routine is slow
        ef0(itemp) = efnew
      ENDIF
      ! We only compute 1 Fermi level so we do not need the other
      efcb(itemp) = 0
      ctype = -1
      !  
    ENDIF
    !
    RETURN   
    !-----------------------------------------------------------------------
    END SUBROUTINE fermicarrier
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    FUNCTION sumkg_seq(et, nbnd, nks, wk, degauss, ngauss, e, is, isk)
    !-----------------------------------------------------------------------
    !!
    !!  This function computes the number of states under a given energy e
    !!
    USE kinds, ONLY : DP
    USE mp,    ONLY : mp_sum
    ! 
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: nks
    !! the total number of K points
    INTEGER, INTENT(in) :: nbnd
    !! the number of bands
    INTEGER, INTENT(in) :: ngauss
    !! the type of smearing
    INTEGER, INTENT(in) :: is
    !!
    INTEGER, INTENT(in) :: isk(nks)
    !!
    REAL(KIND = DP), INTENT(in) :: wk (nks)
    !! the weight of the k points
    REAL(KIND = DP), INTENT(in) :: et (nbnd, nks)
    !! the energy eigenvalues
    REAL(KIND = DP), INTENT(in) :: degauss
    !! gaussian broadening
    REAL(KIND = DP), INTENT(in) :: e
    !! the energy to check
    REAL(KIND = DP)  :: sumkg_seq
    !! Output of the function
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k points
    INTEGER :: ibnd
    !! Counter on the band energy
    REAL(KIND = DP) ::sum1
    !! Sum
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Function which compute the smearing 
    !
    sumkg_seq = 0.d0
    DO ik = 1, nks
      sum1 = 0.d0
      IF (is /= 0) THEN
        IF (isk(ik) /= is) CYCLE
      ENDIF
      DO ibnd = 1, nbnd
        sum1 = sum1 + wgauss((e - et(ibnd, ik)) / degauss, ngauss)
      ENDDO
      sumkg_seq = sumkg_seq + wk (ik) * sum1
    ENDDO
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END FUNCTION sumkg_seq
    !--------------------------------------------------------------------------
    ! 
    !--------------------------------------------------------------------
    FUNCTION efermig_seq(et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk)
    !--------------------------------------------------------------------
    !!
    !! Finds the Fermi energy - Gaussian Broadening
    !! (see Methfessel and Paxton, PRB 40, 3616 (1989 )
    !!
    USE io_global,     ONLY : stdout
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : ryd2ev, eps10
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nks
    !! Number of k-points per pool
    INTEGER, INTENT(in) :: nbnd
    !! Number of band
    INTEGER, INTENT(in) :: Ngauss
    !! 
    INTEGER, INTENT(in) :: is
    !! 
    INTEGER, INTENT(in) :: isk(nks)
    !! 
    REAL(KIND = DP), INTENT(in) :: wk(nks)
    !!
    REAL(KIND = DP), INTENT(in) :: et(nbnd, nks)
    !!
    REAL(KIND = DP), INTENT(in) :: Degauss
    !!
    REAL(KIND = DP), INTENT(in) :: nelec
    !! Number of electron (charge)
    ! 
    ! Local variables
    INTEGER :: i
    !! Iteration number. 
    INTEGER :: kpoint
    !! k-point index
    INTEGER, PARAMETER :: maxiter = 300
    !! Maximum iteration number 
    REAL(KIND = DP) :: efermig_seq 
    !! Fermi energy (seq)
    REAL(KIND = DP) :: Ef
    !! Fermi energy
    REAL (KIND = DP) :: Eup
    !! Upper bound
    REAL (KIND = DP) :: Elw
    !! Lower bound 
    REAL (KIND = DP) :: sumkup 
    !! Sum upper one
    REAL (KIND = DP) :: sumklw
    !! Sum lower one
    REAL (KIND = DP) :: sumkmid
    !! Sum final
    !
    ! Find bounds for the Fermi energy. Very safe choice!
    !
    Elw = et(1, 1)
    Eup = et(nbnd, 1)
    DO kpoint = 2, nks
      Elw = MIN(Elw, et(1, kpoint))
      Eup = MAX(Eup, et(nbnd, kpoint))
    ENDDO
    Eup = Eup + 2 * Degauss
    Elw = Elw - 2 * Degauss
    !
    ! Bisection method
    !
    sumkup = sumkg_seq(et, nbnd, nks, wk, Degauss, Ngauss, Eup, is, isk)
    sumklw = sumkg_seq(et, nbnd, nks, wk, Degauss, Ngauss, Elw, is, isk)
    IF ((sumkup - nelec) < -eps10 .OR. (sumklw - nelec) > eps10) THEN
      CALL errore ('efermig_seq', 'internal error, cannot bracket Ef', 1)
    ENDIF
    DO i = 1, maxiter
      Ef = (Eup + Elw) / 2.d0
      sumkmid = sumkg_seq(et, nbnd, nks, wk, Degauss, Ngauss, Ef, is, isk)
      IF (ABS(sumkmid-nelec) < eps10) THEN
        efermig_seq = Ef
        RETURN
      ELSEIF ((sumkmid - nelec) < -eps10) THEN
        Elw = Ef
      ELSE
        Eup = Ef
      ENDIF
    ENDDO
    IF (is /= 0) WRITE(stdout, '(5x, "Spin Component #", i3)') is
    WRITE(stdout, '(5x,"Warning: too many iterations in bisection"/ &
         &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) Ef * ryd2ev, sumkmid
    !
    efermig_seq = Ef
    RETURN
    !
    !--------------------------------------------------------------------------
    END FUNCTION efermig_seq
    !--------------------------------------------------------------------------
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
    WRITE(stdout, '(/,5x,a, i13, a,f7.2,a,a)') "Number of ep-matrix elements per pool :", &
         imelt, " ~= ", rmelt, TRIM(chunit), " (@ 8 bytes/ DP)"
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE mem_size
    !--------------------------------------------------------------------------
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE broadening(ik, ikk, ikq, w, vmefp, eta)
    !--------------------------------------------------------------------------
    !!
    !! This routine computes the adaptative broadening
    !! It requires electronic and phononic velocities
    !! The implemented equation is Eq. 18 of Computer Physics Communications 185, 1747 (2014)
    !! Samuel Ponce & Francesco Macheda
    !!
    USE cell_base,     ONLY : alat, bg
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : nbndfst, nkf, dmef, vmef, ibndmin, etf
    USE epwcom,        ONLY : vme, nqf1, nqf2, nqf3
    USE phcom,         ONLY : nmodes
    USE constants_epw, ONLY : eps40, ryd2mev, twopi, zero, eps6, eps8, eps4
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Current k-point on that core
    INTEGER, INTENT(in) :: ikk
    !! Current k point on that core (ikk = 2 * ik + 1)
    INTEGER, INTENT(in) :: ikq
    !! k+q point on that core
    REAL(KIND = DP), INTENT(in) :: w(nmodes)
    !! Phonon frequencies
    REAL(KIND = DP), INTENT(out) :: eta(nmodes, nbndfst, nkf)
    !! Adaptative smearing value
    COMPLEX(KIND = DP), INTENT(in) :: vmefp(3, nmodes, nmodes)
    !! Phonon velocity
    !
    ! Local variables
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: imode, jmode
    !! Mode index
    INTEGER :: n_av
    !! To average eta_av
    REAL(KIND = DP) :: vel_diff(3)
    !! Velocity difference when computed adaptative broadening
    REAL(KIND = DP) :: eta_tmp(3)
    !! Temporary adaptative broadening
    REAL(KIND = DP) :: eta_deg(nmodes, nbndfst)
    !! Average eta over degenerate states
    REAL(KIND = DP) :: e_1
    !! Eigenvalue 1 for deg. testing
    REAL(KIND = DP) :: e_2
    !! Eigenvalue 2 for deg. testing
    REAL(KIND = DP) :: w_1
    !! Phonon frequency for degeneracy checking
    REAL(KIND = DP) :: w_2
    !! Phonon frequency for degeneracy checking
    REAL(KIND = DP) :: vmeq(3, nmodes)
    !! Local phonon velocity
    REAL(KIND = DP) :: vmek(3, nbndfst)
    !! Local electron velocity
    REAL(KIND = DP) :: vmeq_av(3)
    !! Average phonon velocity
    REAL(KIND = DP) :: vmek_av(3)
    !! Average phonon velocity
    !
    eta_deg(:, :) = zero
    vmeq(:, :) = zero 
    vmek(:, :) = zero 
    ! 
    ! First average the phonon velocities
    DO imode = 1, nmodes
      w_1 = w(imode)
      vmeq_av(:) = zero 
      n_av = 0
      DO jmode = 1, nmodes
        w_2 = w(jmode)
        IF (ABS(w_2 - w_1) < eps6) THEN
          n_av   = n_av + 1
          vmeq_av(:) = vmeq_av(:) + REAL(vmefp(:, jmode, jmode), KIND = DP) 
        ENDIF
      ENDDO
      vmeq(:, imode) = vmeq_av(:) / FLOAT(n_av) 
    ENDDO
    ! 
    ! Average electron velocity
    DO ibnd = 1, nbndfst
      e_1 = etf(ibndmin - 1 + ibnd, ikk)
      vmek_av(:) = zero
      n_av   = 0
      DO jbnd = 1, nbndfst
        e_2 = etf(ibndmin - 1 + jbnd, ikk)
        IF (ABS(e_2 - e_1) < eps4) THEN
          n_av = n_av + 1
          IF (vme) THEN
            vmek_av(:) = vmek_av(:) + REAL(vmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq), KIND = DP)
          ELSE 
            vmek_av(:) = vmek_av(:) + REAL(dmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq), KIND = DP) 
          ENDIF 
        ENDIF
      ENDDO
      vmek(:, ibnd) = vmek_av(:) / FLOAT(n_av)      
    ENDDO  
    ! 
    ! vmefp and vmef are obtained using irvec, which are without alat; therefore I multiply them to bg without alat
    DO ibnd = 1, nbndfst
      DO imode = 1, nmodes
        IF (w(imode) > 0) THEN
          vel_diff(:) = vmeq(:, imode) / (2d0 * w(imode)) - vmek(:, ibnd)
          !IF (vme) THEN
          !  vel_diff(:) = REAL(vmefp(:, imode, imode) / &
          !                    (2d0 * w(imode)) - vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikq))
          !ELSE
          !  vel_diff(:) = REAL(vmefp(:, imode ,imode) / &
          !                    (2d0 * w(imode)) - dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikq))
          !ENDIF
          IF (SQRT(DOT_PRODUCT(vel_diff, vel_diff)) < eps40) THEN
            eta(imode, ibnd, ik) = 1.0d0 / ryd2mev
          ELSE
            eta_tmp(1) = (twopi / alat) * ABS(DOT_PRODUCT(vel_diff(:), bg(:, 1)) / DBLE(nqf1))
            eta_tmp(2) = (twopi / alat) * ABS(DOT_PRODUCT(vel_diff(:), bg(:, 2)) / DBLE(nqf2))
            eta_tmp(3) = (twopi / alat) * ABS(DOT_PRODUCT(vel_diff(:), bg(:, 3)) / DBLE(nqf3))
            !eta(imode, ibnd, ik) = MAXVAL(eta_tmp) !Eq. (24) of PRB 97 075405 (2015)
            !eta(imode, ibnd, ik) = DSQRT(eta_tmp(1)**2+eta_tmp(2)**2+eta_tmp(3)**2)/DSQRT(12d0) !Eq. (18) of Computer Physics Communications 185 (2014) 1747–1758
            ! The prefactor 0.5 is arbitrary and is to speedup convergence
            eta(imode, ibnd, ik) = 0.5d0 * DSQRT(eta_tmp(1)**2+eta_tmp(2)**2+eta_tmp(3)**2) / SQRT(12d0)
          ENDIF
          ! 
          ! If the smearing is too small, set 1 meV. Too small value are numerically unstable. 
          IF (eta(imode, ibnd, ik) * ryd2mev < 1.0d0) THEN
            eta(imode, ibnd, ik) = 1.0d0 / ryd2mev
          ENDIF
        ELSE
          ! Fixed value 1 meV
          eta(imode, ibnd, ik) = 1.0d0 / ryd2mev
        ENDIF
      ENDDO
    ENDDO
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE broadening
    !--------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  END MODULE low_lvl
  !-------------------------------------------------------------------------
