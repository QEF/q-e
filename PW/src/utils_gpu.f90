!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! General-purpose routines for scalar products, printing,
! linear-algebra operators for exact-exchange and localization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE matcalc_k_gpu (label, DoE, PrtMat, ik, ninner, n, m, U, V, mat, ee)
  !
  USE kinds,                ONLY : dp
  USE io_global,ONLY : stdout
  USE wvfct,                ONLY : wg, npwx
  USE wvfct_gpum,           ONLY : using_wg_d,wg_d
  USE becmod_subs_gpum,     ONLY : calbec_gpu
  USE noncollin_module,     ONLY : noncolin, npol
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: PrtMat, ik, ninner, n, m
  COMPLEX(dp), INTENT(IN) :: U(ninner,n), V(ninner,m)
  COMPLEX(dp), INTENT(OUT):: mat(n,m)
  REAL(DP), INTENT(OUT) :: ee
  CHARACTER(len=*), INTENT(IN) :: label
#if defined(__CUDA)
  attributes(DEVICE) :: U, V, mat
#endif
  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock_gpu('matcalc')

  string = 'M-'
  mat = (0.0_dp, 0.0_dp)
  IF(noncolin) THEN
    noncolin = .false.
    CALL calbec_gpu(ninner, U, V, mat, m)
    noncolin = .true.
  ELSE
    CALL calbec_gpu(ninner, U, V, mat, m)
  ENDIF

  IF(DoE) THEN
    CALL using_wg_d(0)
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    IF( PrtMat > 1 ) CALL errore("matcalc_k_gpu", "matcalc_k_gpu cannot print matrix", 1)
    string = 'E-'
    ee = 0.0_dp
    !$cuf kernel do (1) 
    DO i = 1,n
      ee = ee + wg_d(i,ik)*DBLE(mat(i,i))
    ENDDO
    IF ( PrtMat > 0 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock_gpu('matcalc')

END SUBROUTINE matcalc_k_gpu

