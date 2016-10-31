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
  USE klist,                ONLY : xk, ngk, nks, nkstot, wk, lgauss, &
                                   two_fermi_energies
  USE ktetra,               ONLY : ltetra
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ener,                 ONLY : ef, ef_up, ef_dw 
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
      band_energy, & ! band energy
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


!----------
! To calculate the Band Energy (AlexS)
!
     if(lforcet) then
      call weights()
      band_energy = 0.d0
      do ik = 1, nks
       do i = 1, nbnd
        band_energy = band_energy + et(i,ik) * wg(i,ik)
       enddo
      enddo
      CALL mp_sum( band_energy, inter_pool_comm )
      write(6,*) 
      write(6,*) '------'
      write(6,*)  'eband, Ef (eV) = ',band_energy*rytoev,ef*rytoev
      write(6,*) '------'
      write(6,*)
     endif
!----------

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
     CALL get_homo_lumo (ehomo, elumo)
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
  ! ... Compute estimated HOMO and LUMO from electron counting
  ! ... This is done also for metals, in order to check if there is a gap
  ! ... If LUMO is not available, a large 1.0D+6 value is returned
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nelec, nkstot, nelup, neldw
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : nbnd, et
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: &
      ehomo, elumo   ! highest occupied and lowest unoccupied levels

  INTEGER :: &
      kbnd,         &! counter on bands
      ibnd_up,      &! position of HOMO (spin-up bands)
      ibnd_dw,      &! position of HOMO (spin-down bands)
      ibnd           ! position of HOMO (non-LSDA case)
  !
  ehomo=-1D+6
  elumo=+1D+6
  !
  ibnd    = 0
  ibnd_up = 0
  ibnd_dw = 0
  !
  ! ... store in ibnd the position of the presumed HOMO
  ! ... (or in ibnd_up, ibnd_dw for LSDA calculations)
  !
  IF ( tfixed_occ ) THEN
     DO kbnd = 1, nbnd
        IF ( nspin == 1 .OR. nspin == 4 ) THEN
           IF ( f_inp(kbnd,1) > 0.D0 ) ibnd = kbnd
        ELSE
           IF ( f_inp(kbnd,1) > 0.D0 ) ibnd_up = kbnd
           IF ( f_inp(kbnd,2) > 0.D0 ) ibnd_dw = kbnd
        END IF
     END DO
  ELSE
     IF ( nspin == 1 ) THEN
        ibnd = NINT( nelec ) / 2
     ELSE IF ( nspin == 4 ) THEN
        ibnd = NINT( nelec )
     ELSE
        ibnd_up = NINT( nelup )
        ibnd_dw = NINT( neldw )
     END IF
  END IF
  !
  ! ... estimate HOMO and LUMO (or just the HOMO if LUMO is not there)
  !
  IF ( nspin == 1 .OR. nspin == 4 ) THEN
     IF ( ibnd > 0 .AND. ibnd <= nbnd ) THEN
        ehomo = MAXVAL( et(ibnd,  1:nkstot) )
     END IF
     IF ( ibnd > 0 .AND. ibnd < nbnd ) THEN
        elumo = MINVAL( et(ibnd+1,1:nkstot) )
     END IF
  ELSE
     IF ( ibnd_up == 0 .AND. ibnd_dw > 0 .AND. ibnd_dw <= nbnd) THEN
        ehomo = MAXVAL( et(ibnd_dw,1:nkstot/2) )
     ELSE IF ( ibnd_dw == 0 .AND. ibnd_up > 0 .AND. ibnd_up <= nbnd) THEN
        ehomo = MAXVAL( et(ibnd_up,1:nkstot/2) )
     ELSE IF ( ibnd_dw > 0 .AND. ibnd_up > 0 .AND. &
               ibnd_dw <= nbnd .AND. ibnd_up <= nbnd) THEN
        ehomo = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                     MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
     END IF
     IF ( ibnd_up < nbnd .AND. ibnd_dw < nbnd ) THEN
        elumo = MIN( MINVAL( et(ibnd_up+1,1:nkstot/2) ), &
                     MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) ) )
     ENDIF
  END IF
  !
END SUBROUTINE get_homo_lumo
