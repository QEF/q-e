!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE iweights( nks, wk, nbnd, nelec, et, Ef, wg, is, isk )
  !--------------------------------------------------------------------
  !! Calculates weights for semiconductors and insulators (bands
  !! are either empty or filled).
  !! On output, Ef is the highest occupied Kohn-Sham level.
  !
  !! NOTE: wg must be (INOUT) and not (OUT) because if is/=0 only terms for
  !! spin=is are initialized; the remaining terms should be kept, not lost.
  !
  USE kinds
  USE noncollin_module,    ONLY: noncolin
  USE mp,                  ONLY: mp_max
  USE mp_pools,            ONLY: inter_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  INTEGER, INTENT(IN) :: is
  !! spin label (0 or 1,2)
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  REAL(DP), INTENT(IN) :: wk(nks)
  !! weight of k points
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! eigenvalues of the hamiltonian
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  REAL(DP), INTENT(OUT) :: Ef
  !! the Fermi energy
  !
  ! ... local variables
  !
  INTEGER :: kpoint, ibnd
  !
  CALL iweights_only( nks, wk, is, isk, nbnd, nelec, wg )
  !
  Ef = -1.0d+20
  !
  DO kpoint = 1, nks
     !
     IF (is /= 0) THEN
        IF (isk(kpoint) /=  is) CYCLE
     ENDIF
     !
     DO ibnd = 1, nbnd
        IF ( wg(ibnd, kpoint) > 0.0_DP ) Ef = MAX( Ef, et(ibnd,kpoint) )
     ENDDO
     !
  ENDDO
  !
  ! find max across pools
  !
  CALL mp_max( ef, inter_pool_comm )
  !
  RETURN
  !
END SUBROUTINE iweights
!
!
!--------------------------------------------------------------------
SUBROUTINE iweights_only( nks, wk, is, isk, nbnd, nelec, wg )
  !--------------------------------------------------------------------
  !! Calculates weights for semiconductors and insulators (bands
  !! are either empty or filled).
  !
  USE kinds
  USE noncollin_module, ONLY: noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  INTEGER, INTENT(IN) :: is
  !! spin label (0 or 1,2)
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  REAL(DP), INTENT(IN) :: wk(nks)
  !! weight of k points
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP), INTENT(OUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  !
  ! ... local variables
  !
  REAL(DP) :: degspin 
  INTEGER :: kpoint, ibnd
  !
  degspin = 2.0_DP
  IF (noncolin) degspin = 1.0_DP
  IF (is /= 0)  degspin = 1.0_DP
  !
  DO kpoint = 1, nks
     !
     IF (is /= 0) THEN
        IF (isk(kpoint) /=  is) CYCLE
     ENDIF
     !
     DO ibnd = 1, nbnd
        !
        IF ( ibnd <= NINT(nelec)/degspin ) THEN
           wg(ibnd,kpoint) = wk(kpoint)
        ELSE
           wg(ibnd,kpoint) = 0.0_DP
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE iweights_only
