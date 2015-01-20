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
  USE io_global,            ONLY : stdout, ionode
  USE klist,                ONLY : xk, ngk, nks, nkstot, wk
  USE lsda_mod,             ONLY : lsda, nspin
  USE wvfct,                ONLY : nbnd, et, wg
  USE control_flags,        ONLY : conv_elec, lbands, iverbosity
  USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
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
  IF ( .NOT. lbands ) CALL get_homo_lumo (ehomo, elumo)
  !
  CALL flush_unit( stdout )
  !
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
  !
END SUBROUTINE print_ks_energies
!
!----------------------------------------------------------------------------
SUBROUTINE get_homo_lumo ( ehomo, elumo )
  !----------------------------------------------------------------------------
  !
  ! ... printout of Kohn-Sham eigenvalues: HOMO and LUMO, or Fermi Energy
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout, ionode
  USE ener,                 ONLY : ef, ef_up, ef_dw 
  USE klist,                ONLY : nelec, nkstot, lgauss, two_fermi_energies,&
                                   nelup, neldw
  USE lsda_mod,             ONLY : nspin
  USE ktetra,               ONLY : ltetra
  USE wvfct,                ONLY : nbnd, et
  USE fixed_occ,            ONLY : f_inp, tfixed_occ, one_atom_occupations
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: &
      ehomo, elumo   ! highest occupied and lowest unoccupied levels

  ! ... a few local variables

  INTEGER :: &
      i,            &! counter on polarization
      ik,           &! counter on k points
      kbnd,         &! counter on bands
      ibnd_up,      &! counter on bands
      ibnd_dw,      &! counter on bands
      ibnd         
  !
  ehomo=-1E+6
  elumo=+1E+6
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
  ELSE
     !
     ! ... presumably not a metal: store in ibnd the position of HOMO
     ! ... (or in ibnd_up, ibnd_dw for LSDA calculations)
     !
     IF ( tfixed_occ ) THEN
        ibnd    = 0
        ibnd_up = 0
        ibnd_dw = 0
        DO kbnd = 1, nbnd
           IF ( nspin == 1 .OR. nspin == 4 ) THEN
              IF ( f_inp(kbnd,1) > 0.D0 ) ibnd = kbnd
           ELSE
              IF ( f_inp(kbnd,1) > 0.D0 ) ibnd_up = kbnd
              IF ( f_inp(kbnd,2) > 0.D0 ) ibnd_dw = kbnd
              ibnd = MAX(ibnd_up, ibnd_dw)
           END IF
        END DO
     ELSE
        IF ( nspin == 1 ) THEN
           ibnd = NINT( nelec ) / 2
        ELSE
           ibnd    = NINT( nelec )
           ibnd_up = NINT( nelup )
           ibnd_dw = NINT( neldw )
        END IF
     END IF
     !
     ! ... print HOMO and LUMO (or just the HOMO if LUMO is not there)
     !
     IF ( ionode .AND. .NOT. one_atom_occupations ) THEN
        !
        IF ( nspin == 1 .OR. nspin == 4 ) THEN
           ehomo = MAXVAL( et(ibnd,  1:nkstot) )
           IF ( nbnd > ibnd ) THEN
              elumo = MINVAL( et(ibnd+1,1:nkstot) )
              WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
           ELSE
              WRITE( stdout, 9043 ) ehomo*rytoev
           ENDIF
        ELSE
           IF ( ibnd_up == 0 ) THEN
              !
              ehomo = MAXVAL( et(ibnd_dw,1:nkstot/2) )
              !
           ELSE IF ( ibnd_dw == 0 ) THEN
              !
              ehomo = MAXVAL( et(ibnd_up,1:nkstot/2) )
              !
           ELSE
              !
              ehomo = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                   MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
              !
           END IF
           IF ( nbnd > ibnd_up .AND. nbnd > ibnd_dw ) THEN
              elumo = MIN( MINVAL( et(ibnd_up+1,1:nkstot/2) ), &
                   MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) ) )
              WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
           ELSE
              WRITE( stdout, 9043 ) ehomo*rytoev
           ENDIF
        END IF
        !
     END IF
  END IF
  !
  RETURN
  !
  ! ... formats
  !
9043 FORMAT(/'     highest occupied level (ev): ',F10.4 )
9042 FORMAT(/'     highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9041 FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev' )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
  !
END SUBROUTINE get_homo_lumo
