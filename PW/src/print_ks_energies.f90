!
! Copyright (C) 2007-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!----------------------------------------------------------------------------
SUBROUTINE print_ks_energies()
  !----------------------------------------------------------------------------
  !
  ! ... printout of Kohn-Sham eigenvalues
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : xk, ngk, nks, nkstot, wk, lgauss, ltetra, &
                                   two_fermi_energies
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ener,                 ONLY : ef, ef_up, ef_dw, eband 
  USE lsda_mod,             ONLY : lsda, nspin
  USE spin_orb,             ONLY : lforcet
  USE wvfct,                ONLY : nbnd, et, wg
  USE control_flags,        ONLY : conv_elec, lbands, iverbosity
  USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE mp_pools,             ONLY : inter_pool_comm 

  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
  INTEGER, ALLOCATABLE :: &
      ngk_g(:)       ! number of plane waves summed on all nodes
  REAL(DP) :: &
      ehomo, elumo   ! highest occupied and lowest unoccupied levels
  INTEGER :: &
      i,            &! counter on polarization
      ik,           &! counter on k points
      ibnd           ! counter on bands
  !
  IF (nkstot >= 100 .and. iverbosity <= 0 ) THEN
     WRITE( stdout, '(/,5x,a)') &
     "Number of k-points >= 100: set verbosity='high' to print the bands."
  ELSE
     !
     ALLOCATE ( ngk_g (nkstot) ) 
     !
     ngk_g(1:nks) = ngk(:)
     CALL mp_sum( ngk_g(1:nks), intra_bgrp_comm )
     CALL ipoolrecover( ngk_g, 1, nkstot, nks )
     CALL mp_bcast( ngk_g, root_bgrp, intra_bgrp_comm )
     CALL mp_bcast( ngk_g, root_bgrp, inter_bgrp_comm )
!
! band energy is not available in non-scf calculations (AlexS)
!
     IF (lforcet) THEN
        eband = 0.0_dp
        DO ik = 1, nks
           DO i = 1, nbnd
              eband = eband + et(i,ik) * wg(i,ik)
           END DO
        END DO
        CALL mp_sum( eband, inter_pool_comm )
        WRITE (stdout,'(/,"------")')
        WRITE (stdout,*)  'eband, Ef (eV) = ',eband*rytoev,ef*rytoev
        WRITE (stdout,'("------",/)')
     ENDIF
     !
     DO ik = 1, nkstot
        !
        IF ( lsda ) THEN
           !
           IF ( ik == 1 ) WRITE( stdout, 9015)
           IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016)
           !
        END IF
        !
        IF ( conv_elec ) THEN
           WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngk_g(ik)
        ELSE
           WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
        END IF
        !
        WRITE( stdout, 9030 ) ( et(ibnd,ik) * rytoev, ibnd = 1, nbnd )
        !
        IF( iverbosity > 0 .AND. .NOT. lbands ) THEN
           !
           WRITE( stdout, 9032 )
           IF (ABS(wk(ik))>1.d-10) THEN
              WRITE( stdout, 9030 ) ( wg(ibnd,ik)/wk(ik), ibnd = 1, nbnd )
           ELSE
              WRITE( stdout, 9030 ) ( wg(ibnd,ik), ibnd = 1, nbnd )
           ENDIF
           !
        END IF
        !
     END DO
     !
     DEALLOCATE ( ngk_g )
     !
  ENDIF
  !
  ! ... print HOMO/Top of the VB and LUMO/Bottom of the CB, or E_Fermi
  !
  IF ( .NOT. lbands ) THEN
     !
     IF ( lgauss .OR. ltetra ) THEN
        !
        ! ... presumably a metal: print Fermi energy
        !
        IF ( two_fermi_energies ) THEN
           WRITE( stdout, 9041 ) ef_up*rytoev, ef_dw*rytoev
        ELSE
           WRITE( stdout, 9040 ) ef*rytoev
        END IF
        !
     ELSE IF ( .NOT. one_atom_occupations ) THEN
        !
        ! ... presumably not a metal: print HOMO (and LUMO if available)
        !
        CALL get_homo_lumo (ehomo, elumo)
        !
        IF ( elumo < 1d+6) THEN
           WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
        ELSE
           WRITE( stdout, 9043 ) ehomo*rytoev
        END IF
        !
     END IF
     !
  END IF
  !
  FLUSH( stdout )
  RETURN
  !
  ! ... formats
  !
9015 FORMAT(/' ------ SPIN UP ------------'/ )
9016 FORMAT(/' ------ SPIN DOWN ----------'/ )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
9021 FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   bands (ev):'/ )
9030 FORMAT( '  ',8F9.4 )
9032 FORMAT(/'     occupation numbers ' )
9043 FORMAT(/'     highest occupied level (ev): ',F10.4 )
9042 FORMAT(/'     highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9041 FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev' )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
!
  !
END SUBROUTINE print_ks_energies
!
!----------------------------------------------------------------------------
SUBROUTINE get_homo_lumo ( ehomo, elumo )
  !----------------------------------------------------------------------------
  !
  ! ... Compute estimated HOMO and LUMO from occupations
  ! ... HOMO = largest  E_k with occupation > eps (set to 0.001, see below)
  ! ... LUMO = smallest E_k with occupation < eps
  ! ... Can be done also for metals, in order to check if there is a gap
  ! ... If LUMO is not available, a large 1.0D+6 value is returned
  ! ... In parallel execution, only "ionode" returns the correct values
  !
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nkstot, wk
  USE wvfct,                ONLY : nbnd, et, wg
  USE io_global,            ONLY : ionode
  !
  IMPLICIT NONE
  !
  REAL(dp), PARAMETER :: eps = 0.001_dp ! threshold for zero occupancy
  REAL(DP), INTENT(OUT) :: &
      ehomo, elumo   ! highest occupied and lowest unoccupied levels

  INTEGER :: &
      kbnd,         &! possible position of HOMO
      ibnd, ik       ! counters on bands and k-points
  !
  ehomo=-1D+6
  elumo=+1D+6
  !
  IF ( .NOT. ionode ) RETURN
  !
  k_loop: DO ik = 1, nkstot
     ! exclude states with zero weight (present in phonon calculation)
     IF ( ABS(wk(ik)) > 1.d-10) THEN
        kbnd = nbnd
        band_loop: DO ibnd = 1, nbnd
           ! Allow for negative occupancies (arising in Methfessel-Paxton)
           IF ( ABS(wg(ibnd,ik)) / wk(ik) < eps ) THEN
              kbnd = ibnd - 1
              EXIT band_loop
           END IF
        END DO band_loop
        ehomo = MAX ( ehomo, et(kbnd,ik) )
        IF ( kbnd < nbnd ) THEN
           elumo = MIN ( elumo, et(kbnd+1,ik) )
        END IF
     END IF
  END DO k_loop
  !
END SUBROUTINE get_homo_lumo
