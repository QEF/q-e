! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)
!
!----------------------------------------------------------------------
SUBROUTINE wannier_occupancies( occ )
  !----------------------------------------------------------------------
  !! This routine computes occupation of each wannier. It is assumed that
  !! WF generated already and stored if the buffer.
  !
  USE kinds,        ONLY: DP
  USE wannier_new,  ONLY: nwan, pp
  USE io_global,    ONLY: stdout 
  USE wvfct,        ONLY: nbnd, et, wg
  USE klist,        ONLY: nks
  USE lsda_mod,     ONLY: current_spin, lsda, nspin, isk
  USE io_files
  USE buffers
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: occ(nwan, nwan, nspin)
  !! the occupation of each wannier
  !
  ! ... local variables
  !
  INTEGER :: i,j,k,ik
  !
  occ = ZERO
  current_spin = 1
  !
  DO ik = 1, nks
     IF (lsda) current_spin  = isk(ik)
     CALL get_buffer( pp, nwordwpp, iunwpp, ik)
     DO i = 1, nwan
        DO j = 1,nwan
           DO k = 1, nbnd
              occ(i,j,current_spin) = occ(i,j,current_spin) + pp(i,k) * &
                                                      CONJG(pp(j,k))*wg(k,ik)
           END DO
        END DO
     END DO
  END DO
  !
  IF (nspin==1) occ = occ * 0.5_DP
  !
  RETURN
  !
END SUBROUTINE wannier_occupancies
!
