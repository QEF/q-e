!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform memcpy and memset on the device with CUDA Fortran
! cuf_memXXX contains a CUF KERNEL to perform the selected operation
!
SUBROUTINE cuf_memcpy_r1d(array, array_in, l)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(l)
  REAL(DP), INTENT(IN) :: array_in(l)
  INTEGER, INTENT(IN) :: l
#if defined(__CUDA)
  attributes(DEVICE) :: array, array_in
#endif
  !
  INTEGER :: i
  !
  !$cuf kernel do
  DO i=1, l
     array(i) = array_in(i)
  ENDDO
  !
END SUBROUTINE cuf_memcpy_r1d
!
SUBROUTINE cuf_memcpy_r2d(array, array_in, l, m)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array( l, m )
  REAL(DP), INTENT(IN) :: array_in( l, m )
  INTEGER, INTENT(IN) :: l, m
#if defined(__CUDA)
  attributes(DEVICE) :: array, array_in
#endif

  !
  INTEGER :: i, j
  !
  !$cuf kernel do (2)
  DO j=1, m
     DO i=1, l
        array(i, j) = array_in(i, j)
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memcpy_r2d
!
SUBROUTINE cuf_memcpy_r3d(array, array_in, l, m, n)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array( l, m , n )
  REAL(DP), INTENT(IN) :: array_in( l, m , n )
  INTEGER, INTENT(IN) :: l, m, n
#if defined(__CUDA)
  attributes(DEVICE) :: array, array_in
#endif
  !
  INTEGER :: i, j, k
  !
  !$cuf kernel do(3)
  DO k=1, n
     DO j=1, m
        DO i=1, l
           array(i, j, k) = array_in(i, j, k)
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memcpy_r3d
!
SUBROUTINE cuf_memset_r1d(array, val, l)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(l)
  REAL(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: l
#if defined(__CUDA)
  attributes(DEVICE) :: array
#endif
  !
  INTEGER :: i
  !
  !$cuf kernel do
  DO i=1, l
     array(i) = val
  ENDDO
  !
END SUBROUTINE cuf_memset_r1d
!
SUBROUTINE cuf_memset_r2d(array, val, l, m)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array( l, m )
  REAL(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: l, m
#if defined(__CUDA)
  attributes(DEVICE) :: array
#endif

  !
  INTEGER :: i, j
  !
  !$cuf kernel do (2)
  DO j=1, m
     DO i=1, l
        array(i, j) = val
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memset_r2d
!
SUBROUTINE cuf_memset_r3d(array, val, l, m , n)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array( l, m , n )
  REAL(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: l, m , n
#if defined(__CUDA)
  attributes(DEVICE) :: array
#endif
  !
  INTEGER :: i, j, k
  !
  !$cuf kernel do(3)
  DO k=1, n
     DO j=1, m
        DO i=1, l
           array(i, j, k) = val
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memset_r3d
!
SUBROUTINE cuf_memcpy_c1d(array, array_in, l)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array(l)
  COMPLEX(DP), INTENT(IN) :: array_in(l)
  INTEGER, INTENT(IN) :: l
#if defined(__CUDA)
  attributes(DEVICE) :: array, array_in
#endif
  !
  INTEGER :: i
  !
  !$cuf kernel do
  DO i=1, l
     array(i) = array_in(i)
  ENDDO
  !
END SUBROUTINE cuf_memcpy_c1d
!
SUBROUTINE cuf_memcpy_c2d(array, array_in, l, m)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array( l, m )
  COMPLEX(DP), INTENT(IN) :: array_in( l, m )
  INTEGER, INTENT(IN) :: l, m
#if defined(__CUDA)
  attributes(DEVICE) :: array, array_in
#endif

  !
  INTEGER :: i, j
  !
  !$cuf kernel do (2)
  DO j=1, m
     DO i=1, l
        array(i, j) = array_in(i, j)
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memcpy_c2d
!
SUBROUTINE cuf_memcpy_c3d(array, array_in, l, m , n)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array( l, m , n )
  COMPLEX(DP), INTENT(IN) :: array_in( l, m , n )
  INTEGER, INTENT(IN) :: l, m , n
#if defined(__CUDA)
  attributes(DEVICE) :: array, array_in
#endif
  !
  INTEGER :: i, j, k
  !
  !$cuf kernel do(3)
  DO k=1, n
     DO j=1, m
        DO i=1, l
           array(i, j, k) = array_in(i, j, k)
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memcpy_c3d
!
SUBROUTINE cuf_memset_c1d(array, val, l)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array(l)
  COMPLEX(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: l
#if defined(__CUDA)
  attributes(DEVICE) :: array
#endif
  !
  INTEGER :: i
  !
  !$cuf kernel do
  DO i=1, l
     array(i) = val
  ENDDO
  !
END SUBROUTINE cuf_memset_c1d
!
SUBROUTINE cuf_memset_c2d(array, val, l, m)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array( l, m )
  COMPLEX(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: l, m
#if defined(__CUDA)
  attributes(DEVICE) :: array
#endif

  !
  INTEGER :: i, j
  !
  !$cuf kernel do (2)
  DO j=1, m
     DO i=1, l
        array(i, j) = val
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memset_c2d
!
SUBROUTINE cuf_memset_c3d(array, val, l, m , n)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array( l, m , n )
  COMPLEX(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: l, m , n
#if defined(__CUDA)
  attributes(DEVICE) :: array
#endif
  !
  INTEGER :: i, j, k
  !
  !$cuf kernel do(3)
  DO k=1, n
     DO j=1, m
        DO i=1, l
           array(i, j, k) = val
        ENDDO
     ENDDO
  ENDDO
  !
END SUBROUTINE cuf_memset_c3d
!
