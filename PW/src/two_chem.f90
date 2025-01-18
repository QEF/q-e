!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE two_chem
  !
  !! This module contains subroutines for the calculation using two chemical potentials
  !
  IMPLICIT NONE
  !
  LOGICAL :: twochem
  !determines whether the calculation is performing with two chemical potentials,
  ! one for the holes and one for the electrons calculation
  !
  PUBLIC :: gweights_twochem, gweights_only_twochem
  SAVE
  !
CONTAINS 
!----------------------------------------------------------------------
SUBROUTINE init_twochem()
  !----------------------------------------------------------------------
  !! This routine checks that the two chemical potential calculation 
  !! is inizialized correctly. 
  !
  USE input_parameters,     ONLY : occupations
  USE kinds
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : degauss_cond, nelec, nelec_cond, two_fermi_energies
  USE wvfct,                ONLY : nbnd, nbnd_cond
  USE noncollin_module,     ONLY : noncolin
  USE io_files,             ONLY : restart_dir
  USE io_global,            ONLY : ionode
  USE control_flags,        ONLY : use_gpu
 !
 !
  IMPLICIT NONE

  !
  IF (nbnd_cond==0) THEN
     !
     IF (noncolin) THEN
        !
        nbnd_cond = nbnd - NINT(nelec)
        !
     ELSE
        ! 
        nbnd_cond = nbnd - NINT(nelec)/2
        !
     END IF
     !
  END IF
  !
  !initializes the number of conduction band if the user didn't
  !
  WRITE(stdout, *) "---------------------------------2CHEM----------------------------------"
  WRITE(stdout, *) " You are performing a constrained density-functional perturbation theory"
  WRITE(stdout, *) " employing two chemical potentials, one for electrons and one for holes."
  WRITE(stdout, *) " Please refer to: "
  WRITE(stdout, *) " Giovanni Marini, Matteo Calandra "
  WRITE(stdout, *) " Lattice dynamics of photoexcited insulators" 
  WRITE(stdout, *) " constrained density-functional perturbation theory" 
  WRITE(stdout, *) " Phys. Rev. B 104, 144103 (2021)"
  WRITE(stdout, *) " doi:10.1103/PhysRevB.104.144103"
  WRITE(stdout, *) 
  WRITE(stdout, 9060) nbnd_cond
  WRITE(stdout, 9061) nelec_cond
  WRITE(stdout, *) "---------------------------------2CHEM----------------------------------"
  !
  ! message for two-chemical potential calculation 
  !
  IF (use_gpu) CALL errore('init_twochem', 'twochem with GPU not present in this version',1)  
  IF (trim(occupations) /= 'smearing') CALL errore( 'init_twochem', &
             & 'two chemical potential calculation requires smearing', 1 )
  !
  ! checks that smearing is being employed 
  ! 
  IF (noncolin) THEN
     !
     IF (nbnd_cond.GT.(nbnd - NINT(nelec))) CALL errore( 'init_twochem', &
            & 'non collinear calculation and nbnd_cond > nbnd - NINT(nelec)', 1 )
     !
  ELSE
     ! 
     IF (nbnd_cond.GT.(nbnd - NINT(nelec)/2)) CALL errore( 'init_twochem', &
            & 'collinear calculation and nbnd_cond > nbnd - NINT(nelec)/2', 1 )
     !
  END IF
  ! checks that nbnd_cond <= nbnd - NINT(nelec) (non collinear)
  !     or      nbnd_cond <= nbnd - NINT(nelec)/2 (collinear)

  IF (nelec_cond.GE.nelec) CALL errore( 'init_twochem', &
            & 'nelec_cond greater than nelec', 1 )
  ! checks that the number of electrons in the conduction band 
  ! is less than total number of electrons
  !
  IF (two_fermi_energies) CALL errore( 'init_twochem', &
            & 'fixed total magnetization with twochem not implemented', 1 )
  ! checks that total magnetization is not fixed 
  !
9060 FORMAT( '     The conduction manifold is constituted by',I3, ' bands' )
9061 FORMAT( '    ', F8.4, ' electrons are placed in the conduction manifold' )
  !
  RETURN
  !
END SUBROUTINE init_twochem
!

!--------------------------------------------------------------------
SUBROUTINE gweights_twochem( nks, wk, nbnd,nbnd_cond, nelec,nelec_cond, &
                    degauss,degauss_cond,ngauss, et, ef,ef_cond, demet, wg, is, isk )
  !--------------------------------------------------------------------
  !! Calculates Ef and weights with the gaussian spreading technique in the case of twochem = .TRUE.
  !! with two fermi levels, one for the valence bands and one for the conduction bands 
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
  INTEGER, INTENT(IN) :: nbnd_cond
  !! number of conduction bands
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
  REAL(DP), INTENT(IN) :: nelec_cond
  !! number of conduction electrons
  REAL(DP), INTENT(IN) :: degauss
  !! smearing parameter
  REAL(DP), INTENT(IN) :: degauss_cond
  !! degauss in the conduction band
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  REAL(DP), INTENT(OUT) :: ef,ef_cond
  !! the Fermi energy
  REAL(DP), INTENT(OUT) :: demet
  !! variational correction ("-TS") for metals
  !
  ! Calculate the Fermi energy ef
  !
  ef = efermig_twochem( et,nbnd,1,nbnd-nbnd_cond, nks, nelec-nelec_cond, wk, degauss, ngauss, is, isk )
  ef_cond = efermig_twochem( et,nbnd,nbnd-nbnd_cond+1,nbnd, nks, nelec_cond, wk, degauss_cond, ngauss, is, isk )
  !
  ! Calculate weights
  !
  CALL gweights_only_twochem( nks, wk, is, isk, nbnd,nbnd_cond, nelec,nelec_cond, degauss, &
     degauss_cond,ngauss, et, ef,ef_cond, demet, wg )
  !
  RETURN
  !
END SUBROUTINE gweights_twochem
!
!--------------------------------------------------------------------

SUBROUTINE gweights_only_twochem( nks, wk, is, isk, nbnd,nbnd_cond, nelec,nelec_cond, degauss, &
                          degauss_cond,ngauss, et, ef,ef_cond, demet, wg )
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
    INTEGER, INTENT(IN) :: nbnd_cond
  !! number of conduction bands
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
  REAL(DP), INTENT(IN) :: nelec_cond
  !! number of conduction electrons
  REAL(DP), INTENT(IN) :: degauss
  !! smearing parameter
    REAL(DP), INTENT(IN) :: degauss_cond
  !! degauss in the conduction band
  REAL(DP), INTENT(IN) :: ef,ef_cond
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
     DO ibnd = 1, nbnd-nbnd_cond
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
     DO ibnd = nbnd-nbnd_cond+1,nbnd
        ! Calculate the gaussian weights
        wg(ibnd,kpoint) = wk(kpoint) * &
                            wgauss( (ef_cond-et(ibnd,kpoint))/degauss_cond, ngauss )
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
                 degauss_cond * w1gauss( (ef_cond-et(ibnd,kpoint))/degauss_cond, ngauss )
     ENDDO

     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE gweights_only_twochem
!
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
FUNCTION efermig_twochem( et,nbnd,nbnd1,nbnd2, nks, nelec, wk, Degauss, Ngauss, is, isk )
  !--------------------------------------------------------------------
  !! Finds the Fermi energy in a subset of bands - Gaussian Broadening. 
  !! (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  USE mp,        ONLY : mp_max, mp_min
  USE mp_pools,  ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd,nbnd1,nbnd2
  !! number of bands
  INTEGER, INTENT(IN) :: Ngauss
  !! type of smearing technique
  INTEGER, INTENT(IN) :: is
  !! spin label (0 or 1,2)
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  REAL(DP), INTENT(IN) :: wk(nks)
  !! weight of k points
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! eigenvalues of the hamiltonian
  REAL(DP), INTENT(IN) :: Degauss
  !! smearing parameter
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP) :: efermig_twochem
  !! the Fermi energy
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: eps = 1.0d-10
  INTEGER, PARAMETER :: maxiter = 300
  !
  REAL(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  INTEGER :: i, kpoint, Ngauss_
  !
  !  ... find (very safe) bounds for the Fermi energy:
  !  Elw = lowest, Eup = highest energy among all k-points.
  !  Works with distributed k-points, also if nks=0 on some processor
  !
  Elw = 1.0E+8
  Eup =-1.0E+8
  DO kpoint = 1, nks
     Elw = MIN( Elw, et(nbnd1,kpoint) )
     Eup = MAX( Eup, et(nbnd2,kpoint) )
  ENDDO
  Eup = Eup + 5 * Degauss 
  Elw = Elw - 5 * Degauss
  !
  ! ... find min and max across pools
  !
  CALL mp_max( eup, inter_pool_comm )
  CALL mp_min( elw, inter_pool_comm )
  !
  ! ... Bisection method
  !
  ! ... perform a preliminary determination with the Gaussian broadening
  ! to safely locate Ef mid-gap in the insulating case
  !
  !   Ngauss_ = 0 ! currently disabled
  Ngauss_ = Ngauss
  !
1 CONTINUE
  !
  !write(*,*) nelec
  sumkup = sumkg_twochem( et,nbnd,nbnd1,nbnd2, nks, wk, Degauss, Ngauss_, Eup, is, isk )
  sumklw = sumkg_twochem( et,nbnd,nbnd1,nbnd2, nks, wk, Degauss, Ngauss_, Elw, is, isk )
  !write(*,*) sumkup,sumklw
  IF ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
       CALL errore( 'efermig twochem', 'internal error, cannot bracket Ef', 1 )
  DO i = 1, maxiter
     Ef = (Eup + Elw) / 2.d0
     sumkmid = sumkg_twochem( et,nbnd,nbnd1,nbnd2, nks, wk, Degauss, Ngauss_, Ef, is, isk )
     IF ( ABS( sumkmid-nelec ) < eps) THEN
        efermig_twochem = Ef
        ! refine the search with the input Ngauss value if not already done
        IF (Ngauss /= Ngauss_) THEN
           Elw = Ef - Degauss ; Eup = Ef + Degauss ; Ngauss_ = Ngauss
           GOTO 1
        ENDIF
        RETURN
     ELSEIF ( (sumkmid-nelec) < -eps) THEN
        Elw = Ef
     ELSE
        Eup = Ef
     ENDIF
  ENDDO
  IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
       Ef * rytoev, sumkmid
  !
  efermig_twochem = Ef
  !
  !
  RETURN
  !
END FUNCTION efermig_twochem

!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION sumkg_twochem( et,nbnd,nbnd1,nbnd2, nks, wk, degauss, ngauss, e, is, isk )
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
  REAL(DP) :: sumkg_twochem
  !! output: 
  INTEGER, INTENT(IN) :: nks
  !! the total number of K points
  INTEGER, INTENT(IN) :: nbnd,nbnd1,nbnd2
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
  sumkg_twochem = 0.d0
  !
  DO ik = 1, nks
     !
     sum1 = 0.d0
     IF (is /= 0) THEN
        IF (isk(ik) /= is) CYCLE
     ENDIF
     DO ibnd = nbnd1, nbnd2
        sum1 = sum1 + wgauss( (e-et(ibnd,ik))/degauss, ngauss )
     ENDDO
     sumkg_twochem = sumkg_twochem + wk (ik) * sum1
     !
  ENDDO
  !
#if defined(__MPI)
  CALL mp_sum( sumkg_twochem, inter_pool_comm )
#endif
  !
  RETURN
  !
END FUNCTION sumkg_twochem
END MODULE two_chem
