!
! Copyright (C) 2007-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!----------------------------------------------------------------------------
SUBROUTINE print_ks_energies_nonscf ( ef_scf, ef_scf_up, ef_scf_dw )
  !----------------------------------------------------------------------------
  !
  ! ... printout of Kohn-Sham eigenvalues, Fermi energies, HOMO, LUMO if available
  ! ... The input Fermi energy from scf is passed as argument and printed as well
  !
  USE kinds,                ONLY : dp
  USE control_flags,        ONLY : lbands
  IMPLICIT NONE
  !
  REAL(dp), INTENT(in) :: ef_scf, ef_scf_up, ef_scf_dw
  !
  CALL print_ks_only ( )
  !
  IF ( .NOT. lbands) CALL print_ks_ef_homolumo ( .true., ef_scf, ef_scf_up, ef_scf_dw )
  !
END SUBROUTINE print_ks_energies_nonscf
!
!----------------------------------------------------------------------------
SUBROUTINE print_ks_energies( )
  !----------------------------------------------------------------------------
  !
  ! ... printout of Kohn-Sham eigenvalues, Fermi energies, HOMO, LUMO if available
  !
  USE kinds,                ONLY : dp
  USE control_flags,        ONLY : lbands
  !
  CALL print_ks_only ( )
  !
  IF ( .NOT. lbands) CALL print_ks_ef_homolumo ( .false., 0.0_dp, 0.0_dp, 0.0_dp )
  !
END SUBROUTINE print_ks_energies
!
!----------------------------------------------------------------------------
SUBROUTINE print_ks_only( )
  !----------------------------------------------------------------------------
  !
  ! ... printout of Kohn-Sham eigenvalues
  !
  USE kinds,                ONLY : dp
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : xk, ngk, nks, nkstot, wk
  USE ener,                 ONLY : ef, eband 
  USE lsda_mod,             ONLY : lsda, nspin
  USE noncollin_module,     ONLY : lforcet
  USE wvfct,                ONLY : nbnd, et, wg
  USE control_flags,        ONLY : conv_elec, lbands, iverbosity
  USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE mp_pools,             ONLY : inter_pool_comm 
  USE add_dmft_occ,         ONLY : dmft_updated
  !
  USE wvfct_gpum,           ONLY : using_et
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
  INTEGER, ALLOCATABLE :: &
      ngk_g(:)       ! number of plane waves summed on all nodes
  INTEGER :: &
      i,            &! counter on polarization
      ik,           &! counter on k points
      ibnd           ! counter on bands
  !
  CALL using_et(0)
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
        ELSEIF ( dmft_updated ) THEN
           WRITE( stdout, 9019 ) ( xk(i,ik), i = 1, 3 )
        ELSE
           WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
        END IF
        !
        IF ( .NOT. dmft_updated ) WRITE( stdout, 9030 ) ( et(ibnd,ik) * rytoev, ibnd = 1, nbnd )
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
  ! ... formats
  !
9015 FORMAT(/' ------ SPIN UP ------------'/ )
9016 FORMAT(/' ------ SPIN DOWN ----------'/ )
9019 FORMAT(/'          k =',3F7.4,':' )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
9021 FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   bands (ev):'/ )
9030 FORMAT( '  ',8F9.4 )
9032 FORMAT(/'     occupation numbers ' )
  !
END SUBROUTINE print_ks_only
!
!----------------------------------------------------------------------------
SUBROUTINE print_ks_ef_homolumo ( print_ef_scf, ef_scf, ef_scf_up, ef_scf_dw )
  !----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : dp
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout
  USE fixed_occ,            ONLY : one_atom_occupations
  USE klist,                ONLY : two_fermi_energies, lgauss, ltetra
  USE ener,                 ONLY : ef, ef_up, ef_dw
  !
  IMPLICIT NONE
  LOGICAL, INTENT(in) :: print_ef_scf
  REAL(dp), INTENT(in) :: ef_scf, ef_scf_up, ef_scf_dw
  REAL(dp) :: ehomo, elumo ! highest occupied and lowest unoccupied levels
  !
  ! ... print HOMO/Top of the VB and LUMO/Bottom of the CB, or E_Fermi
  !
  IF ( lgauss .OR. ltetra ) THEN
     !
     ! ... presumably a metal: print Fermi energy
     !
     IF ( two_fermi_energies ) THEN
        WRITE( stdout, 9041 ) ef_up*rytoev, ef_dw*rytoev
        IF ( print_ef_scf ) &
             WRITE( stdout, 9051 ) ef_scf_up*rytoev, ef_scf_dw*rytoev
     ELSE
        WRITE( stdout, 9040 ) ef*rytoev
        IF ( print_ef_scf ) WRITE( stdout, 9050 ) ef_scf*rytoev
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
  FLUSH( stdout )
  !
  ! ... formats
  !
9043 FORMAT(/'     highest occupied level (ev): ',F10.4 )
9042 FORMAT(/'     highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9041 FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev' )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
9051 FORMAT( '     (compare with: ',2F10.4,' eV, computed in scf)' )
9050 FORMAT( '     (compare with: ', F10.4,' eV, computed in scf)' )
  !
END SUBROUTINE print_ks_ef_homolumo
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
  USE wvfct_gpum,           ONLY : using_et
  !
  IMPLICIT NONE
  !
  REAL(dp), PARAMETER :: eps = 0.001_dp ! threshold for zero occupancy
  REAL(dp), INTENT(OUT) :: &
      ehomo, elumo   ! highest occupied and lowest unoccupied levels

  INTEGER :: &
      kbnd,         &! possible position of HOMO
      ibnd, ik       ! counters on bands and k-points
  !
  CALL using_et(0)
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
        IF ( kbnd > 0 ) ehomo = MAX ( ehomo, et(kbnd,ik) )
        !
        IF ( kbnd < nbnd ) THEN
           elumo = MIN ( elumo, et(kbnd+1,ik) )
        END IF
     END IF
  END DO k_loop
  !
END SUBROUTINE get_homo_lumo
