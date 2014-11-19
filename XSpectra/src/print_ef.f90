!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE get_homo_lumo(ehomo,elumo)
  !----------------------------------------------------------------------------
  !
  ! ... printout of ef for a metal
  ! ... printout HOMO and LUMO for an insulator
  ! .... HOMO --> Energy of the Highest Occupied state  
  ! .... LUMO --> Energy of the Lowest Unoccupied state  
  !
  ! adapted from subroutine print_ks_energies in PW/src 
  ! DC June 2014
  !
  !----------------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout, ionode
  USE ener,                 ONLY : ef, ef_up, ef_dw
  USE klist,                ONLY : xk, nelec, ngk, nks, nkstot, &
                                   lgauss, two_fermi_energies, nelup, neldw, &
                                   wk
  USE lsda_mod,             ONLY : lsda, nspin, isk
  USE ktetra,               ONLY : ltetra
  USE wvfct,                ONLY : nbnd, et, wg
  USE fixed_occ,            ONLY : f_inp, tfixed_occ, one_atom_occupations
  USE control_flags,        ONLY : conv_elec, lbands, iverbosity
  USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast


  IMPLICIT NONE
  !
  ! Arguments
  REAL(DP), INTENT(OUT) :: &
      ehomo, elumo   ! highest occupied and lowest unoccupied levels
  !
  ! Local variables
  INTEGER :: &
      i,            &! counter on polarization
      ik,           &! counter on k points
      kbnd,         &! counter on bands
      ibnd_up,      &! counter on bands
      ibnd_dw,      &! counter on bands
      ibnd

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
              IF ( nbnd > ibnd ) elumo = MINVAL( et(ibnd+1,1:nkstot) )
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
              IF ( nbnd > ibnd ) &
                 !
                 elumo = MIN( MINVAL( et(ibnd_up+1,1:nkstot/2) ), &
                              MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) ) )
           END IF
           !
           IF ( nbnd > ibnd ) THEN
              WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
           ELSE
              WRITE( stdout, 9043 ) ehomo*rytoev
           END IF
           !
        END IF
     END IF
  !
  ! ... formats
  !
9043 FORMAT(/'     highest occupied level (ev): ',F10.4 )
9042 FORMAT(/'     highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9041 FORMAT(/'     the spin up/dw Fermi energies are ',2F10.4,' ev' )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
     !
END SUBROUTINE get_homo_lumo

