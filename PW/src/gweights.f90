!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE gweights( nks, wk, nbnd, nelec, degauss, ngauss, &
                     et, ef, demet, wg, is, isk )
  !--------------------------------------------------------------------
  !! Calculates Ef and weights with the gaussian spreading technique.
  !! Wrapper routine: computes first Ef, then the weights.
  !
  !! NOTE: wg must be (INOUT) and not (OUT) because if is/=0 only terms for
  !! spin=is are initialized; the remaining terms should be kept, not lost.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  INTEGER, INTENT(IN) :: ngauss
  !! type of smearing technique
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
  REAL(DP), INTENT(IN) :: degauss
  !! smearing parameter
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  REAL(DP), INTENT(OUT) :: ef
  !! the Fermi energy
  REAL(DP), INTENT(OUT) :: demet
  !! variational correction ("-TS") for metals
  !
  REAL(DP), EXTERNAL :: efermig
  !
  ! Calculate the Fermi energy ef
  !
  ef = efermig( et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk )
  !
  ! Calculate weights
  !
  CALL gweights_only( nks, wk, is, isk, nbnd, nelec, degauss, &
     ngauss, et, ef, demet, wg )
  !
  RETURN
  !
END SUBROUTINE gweights
!
!--------------------------------------------------------------------
subroutine gweights_mix (nks, wk, nbnd, nelec, degauss, ngauss, &
     et, ef, demet, wg, is, isk, beta)
  !--------------------------------------------------------------------
  !     calculates Ef and weights with the gaussian spreading technique
  ! ... Wrapper routine: computes first Ef, then the weights
  ! ... Also IN and OUT Ef are mixed.
  ! ... This routine is called in weight with GC-SCF calculation. 
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks, nbnd, ngauss, is, isk(nks)
  REAL(DP), INTENT(IN) :: wk (nks), et (nbnd, nks), nelec, degauss
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  REAL(DP), INTENT(INOUT) :: wg (nbnd, nks)
  REAL(DP), INTENT(INOUT) :: ef
  REAL(DP), INTENT(OUT) :: demet
  REAL(DP), INTENT(IN) :: beta
  !
  REAL(DP) :: ef_by_n
  REAL(DP), EXTERNAL :: efermig
  !
  ! Calculate the Fermi energy ef
  !
  ef_by_n = efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk)
  !
  ! Mixing the Fermi energy ef
  !
  ef = beta * ef + (1.0_DP - beta) * ef_by_n
  !
  ! Calculate weights
  !
  CALL gweights_only (nks, wk, is, isk, nbnd, nelec, degauss, &
     ngauss, et, ef, demet, wg)
  !
  RETURN
  !
end subroutine gweights_mix
!
!--------------------------------------------------------------------
SUBROUTINE gweights_only( nks, wk, is, isk, nbnd, nelec, degauss, &
                          ngauss, et, ef, demet, wg )
  !--------------------------------------------------------------------
  !! Calculates weights with the gaussian spreading technique.
  !! Fermi energy is provided in input.
  !
  !! NOTE: wg must be (inout) and not (out) because if is/=0 only terms for
  !! spin=is are initialized; the remaining terms should be kept, not lost.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  INTEGER, INTENT(IN) :: ngauss
  !! type of smearing technique
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
  REAL(DP), INTENT(IN) :: degauss
  !! smearing parameter
  REAL(DP), INTENT(IN) :: ef
  !! the Fermi energy
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  REAL(DP), INTENT(OUT) :: demet
  !! variational correction ("-TS") for metals
  !
  ! ... local variables
  !
  INTEGER :: kpoint, ibnd
  REAL(DP) , EXTERNAL :: wgauss, w1gauss
  !
  demet = 0._DP
  !
  DO kpoint = 1, nks
     !
     IF (is /= 0) THEN
        IF (isk(kpoint) /= is) CYCLE
     ENDIF
     !
     DO ibnd = 1, nbnd
        ! Calculate the gaussian weights
        wg(ibnd,kpoint) = wk(kpoint) * &
                            wgauss( (ef-et(ibnd,kpoint))/degauss, ngauss )
        !
        ! The correct (i.e. variational) form of the band energy is 
        !    Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        ! This differs by the term "demet" from the sum of KS eigenvalues:
        !    Eks = \sum wg(n,k) et(n,k)
        ! which is non variational. When a Fermi-Dirac function is used
        ! for a given T, the variational energy is REALly the free energy F,
        ! and F = E - TS , with E = non variational energy, -TS = demet
        !
        demet = demet + wk(kpoint) * &
                 degauss * w1gauss( (ef-et(ibnd,kpoint))/degauss, ngauss )
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE gweights_only
