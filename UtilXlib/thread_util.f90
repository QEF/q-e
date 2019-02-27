!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform threaded memcpy and memset
! threaded_memXXX contains a parallel do region
! threaded_barrier_memXXX contains a do region without parallel
! threaded_nowait_memXXX contains a do region without parallel and a nowait at the end do
!
SUBROUTINE threaded_memcpy(array, array_in, length)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: array_in(length)
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
  !$omp parallel do
  DO i=1, length
     array(i) = array_in(i)
  ENDDO
  !$omp end parallel do
  !
END SUBROUTINE threaded_memcpy

SUBROUTINE threaded_barrier_memcpy(array, array_in, length)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: array_in(length)
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
  !$omp do
  DO i=1, length
     array(i) = array_in(i)
  ENDDO
  !$omp end do
  !
END SUBROUTINE threaded_barrier_memcpy

SUBROUTINE threaded_nowait_memcpy(array, array_in, length)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: array_in(length)
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
  !$omp do
  DO i=1, length
     array(i) = array_in(i)
  ENDDO
  !$omp end do nowait
  !
END SUBROUTINE threaded_nowait_memcpy

SUBROUTINE threaded_memset(array, val, length)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
  !$omp parallel do
  DO i=1, length
     array(i) = val
  ENDDO
  !$omp end parallel do
  !
END SUBROUTINE threaded_memset

SUBROUTINE threaded_barrier_memset(array, val, length)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
  !$omp do
  DO i=1, length
     array(i) = val
  ENDDO
  !$omp end do
  !
END SUBROUTINE threaded_barrier_memset

SUBROUTINE threaded_nowait_memset(array, val, length)
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: val
  INTEGER, INTENT(IN) :: length
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
  !$omp do
  DO i=1, length
     array(i) = val
  ENDDO
  !$omp end do nowait
  !
END SUBROUTINE threaded_nowait_memset
