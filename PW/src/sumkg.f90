!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION sumkg( et, nbnd, nks, wk, degauss, ngauss, e, is, isk )
  !-----------------------------------------------------------------------
  !! This function computes the number of states under a given 
  !! energy \( e \).
  !
  USE kinds
  USE mp_pools,  ONLY: inter_pool_comm
  USE mp,        ONLY: mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sumkg
  !! output: 
  INTEGER, INTENT(IN) :: nks
  !! the total number of K points
  INTEGER, INTENT(IN) :: nbnd
  !! the number of bands
  INTEGER, INTENT(IN) :: ngauss
  !! the type of smearing
  REAL(DP), INTENT(IN) :: wk(nks)
  !! the weight of the k points
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! the energy eigenvalues
  REAL(DP), INTENT(IN) :: degauss
  !! gaussian broadening
  REAL(DP), INTENT(IN) :: e
  !! the energy to check
  INTEGER, INTENT(IN) :: is
  !! the spin label
  INTEGER, INTENT(IN) :: isk(nks)
  !! the spin index for each k-point
  !
  ! ... local variables
  !
  REAL(DP), EXTERNAL :: wgauss
  ! function which compute the smearing
  REAL(DP) :: sum1
  INTEGER :: ik, ibnd
  ! counter on k points
  ! counter on the band energy
  !
  !
  sumkg = 0.d0
  !
  DO ik = 1, nks
     !
     sum1 = 0.d0
     IF (is /= 0) THEN
        IF (isk(ik) /= is) CYCLE
     ENDIF
     DO ibnd = 1, nbnd
        sum1 = sum1 + wgauss( (e-et(ibnd,ik))/degauss, ngauss )
     ENDDO
     sumkg = sumkg + wk (ik) * sum1
     !
  ENDDO
  !
#if defined(__MPI)
  CALL mp_sum( sumkg, inter_pool_comm )
#endif
  !
  RETURN
  !
END FUNCTION sumkg

