!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fermisurfer_common
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & b_low, &
  & b_high
  !
  CONTAINS
  !
!----------------------------------------------------------------------------
SUBROUTINE rotate_k_fs(equiv)
  !--------------------------------------------------------------------------
  !
  ! This routine find the equivalent k-point in irr-BZ for the whole BZ
  ! Also compute the max. and min. band index containing Fermi surfaces.
  !
  USE kinds,     ONLY : DP
  USE klist,     ONLY : xk, nks, two_fermi_energies
  USE symm_base, ONLY : nsym, s, time_reversal, t_rev
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : at
  USE ener,      ONLY : ef, ef_up, ef_dw
  USE start_k,   ONLY : k1, k2, k3, nk1, nk2, nk3
  USE wvfct,     ONLY : nbnd, et
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(OUT) :: equiv(nk1,nk2,nk3)
  !
  INTEGER :: isym, ik, ikv(3), nk, ibnd
  REAL(8) :: xk_frac(3), kv(3), ef1, ef2
  logical :: ldone(nk1,nk2,nk3)
  !
  WRITE(stdout,*)
  WRITE(stdout,'(5x,a,i7)') "Number of bands : ", nbnd
  WRITE(stdout,'(5x,a,i7)') "Number of k times spin : ", nks
  WRITE(stdout,'(5x,a,i7)') "Number of symmetries : ", nsym
  !
  ! Which band contains Fermi level ?
  !
  IF(two_fermi_energies) THEN
     ef1 = ef_up
     ef2 = ef_dw
  ELSE
     ef1 = ef
     ef2 = ef
  END IF
  !
  DO ibnd = 1, nbnd
     !
     IF(MINVAL(et(       ibnd    ,1:nks)) < MAX(ef1, ef2)) b_high = ibnd
     IF(MAXVAL(et(nbnd - ibnd + 1,1:nks)) > MIN(ef1, ef2)) b_low = nbnd - ibnd + 1
     !
  END DO
  !
  WRITE(stdout,'(5x,a,i7)') "Lowest band which contains FS : ", b_low
  WRITE(stdout,'(5x,a,i7)') "Highest band which contains FS : ", b_high
  !
  IF(nspin == 2) THEN
     nk = nks / 2
  ELSE
     nk = nks
  END IF
  ldone(1:nk1, 1:nk2, 1:nk3) = .FALSE.
  !
  DO ik = 1, nk
     !
     xk_frac(1:3) = matmul(xk(1:3,ik), at(1:3,1:3)) * REAL((/nk1, nk2, nk3/), DP)
     !
     DO isym = 1, nsym
        !
        kv(1:3) = MATMUL(REAL(s(1:3,1:3,isym), DP), xk_frac(1:3))
        IF(t_rev(isym)==1) kv(1:3) = - kv(1:3)
        !
        kv(1:3) = kv(1:3) - 0.5_dp * REAL((/k1, k2, k3/), DP)
        ikv(1:3) = NINT(kv(1:3))
        !
        IF(ANY(ABS(kv(1:3) - REAL(ikv(1:3), DP)) > 1d-8)) CYCLE
        !
        ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/)) + 1
        !
        equiv(ikv(1), ikv(2), ikv(3)) = ik
        ldone(ikv(1), ikv(2), ikv(3)) = .TRUE.
        !
        ! Time-Reversal
        !
        IF (time_reversal) THEN
           !
           ikv(1:3) = - (ikv(1:3) - 1) - (/k1, k2, k3/)
           ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/)) + 1
           !
           equiv(ikv(1), ikv(2), ikv(3)) = ik
           ldone(ikv(1), ikv(2), ikv(3)) = .TRUE.
           !
        END IF
        !
     END DO ! END isym
     !
  END DO ! End ik
  !
  ! Check
  !
  IF(COUNT(.NOT. ldone) /= 0) &
  &     WRITE(stdout,*)  "  # of elements that are not done : ", COUNT(.NOT. ldone)
  !
END SUBROUTINE rotate_k_fs
!
!----------------------------------------------------------------------------
SUBROUTINE write_fermisurfer(eig, mat, filename)
  !----------------------------------------------------------------------------
  !
  ! This routine output a matrix element on the Fermi surface
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE start_k,   ONLY : nk1, nk2, nk3, k1, k2, k3
  USE io_global, ONLY : stdout, ionode
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: eig(b_low:b_high,nk1,nk2,nk3), &
  &                      mat(b_low:b_high,nk1,nk2,nk3)
  CHARACTER(*),INTENT(IN) :: filename
  !
  INTEGER :: ibnd, i1, i2, i3, fo
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  WRITE(stdout,'(5x,a,f18.8,5x,a,f18.8)') &
  &  "Max : ", MAXVAL(mat(b_low:b_high,1:nk1,1:nk2,1:nk3)), &
  &  "Min : ", MINVAL(mat(b_low:b_high,1:nk1,1:nk2,1:nk3))
  !
  IF(ionode) THEN
     !
     fo = find_free_unit()
     OPEN(fo, file = TRIM(filename))
     !
     WRITE(fo,'(3i6)') nk1, nk2, nk3
     !
     WRITE(fo,'(i6)') 1 + k1
     !
     WRITE(fo,'(i6)') b_high - b_low + 1
     !
     ! Write with single-precision
     !
     WRITE(fo,*) REAL(bg(1:3,1))
     WRITE(fo,*) REAL(bg(1:3,2))
     WRITE(fo,*) REAL(bg(1:3,3))
     !
     DO ibnd = b_low, b_high
        DO i1 = 1, nk1
           DO i2 = 1, nk2
              DO i3 = 1, nk3
                 WRITE(fo,*) REAL(eig(ibnd,i1,i2,i3))
              END DO
           END DO
        END DO
     END DO
     !
     DO ibnd = b_low, b_high
        DO i1 = 1, nk1
           DO i2 = 1, nk2
              DO i3 = 1, nk3
                 WRITE(fo,*) REAL(mat(ibnd,i1,i2,i3)) 
              END DO
           END DO
        END DO
     END DO
     !
     CLOSE(fo)
     !
  END IF
  !
END SUBROUTINE write_fermisurfer
!
END MODULE fermisurfer_common
