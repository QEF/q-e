!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
#include "f_defs.h"
!
!-----------------------------------------------------------------------
  SUBROUTINE non_scf (ik_)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield, lberry, nberrycyc, fact_hepsi
  USE constants,            ONLY : rytoev
  USE control_flags,        ONLY : lbands
  USE ener,                 ONLY : ef
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunwfc, nwordwfc, iunefield
  USE klist,                ONLY : xk, wk, degauss, nelec, nks, nkstot, &
                                   lgauss, ngauss
  USE ktetra,               ONLY : ltetra, ntetra, tetra  
  USE lsda_mod,             ONLY : lsda, nspin, isk
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE wvfct,                ONLY : nbnd, et, npwx
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: ik_
  !
  ! ... local variables
  !
  INTEGER :: iter = 1, i, ik, inberry, ibnd, ibnd_up, ibnd_dw, kbnd
  !
  REAL(DP) :: dr2 = 0.d0,  ehomo, elumo ! highest occupied and lowest unoccupied levels
  !
  REAL(DP), EXTERNAL :: efermit, efermig, get_clock
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  !
  !
  CALL start_clock( 'electrons' )
  !
  WRITE( stdout, 9002 )
  !
  CALL flush_unit( stdout )
  !
  ! ... diagonalization of the KS hamiltonian
  !
  IF ( lelfield ) THEN
     !
     ALLOCATE(fact_hepsi(nks))
     !
     DO inberry = 1, nberrycyc
        !
        ALLOCATE( psi( npwx, nbnd ) )
        !
        DO ik=1,nks
           !
           CALL davcio( psi, nwordwfc, iunwfc,    ik, -1 )
           CALL davcio( psi, nwordwfc, iunefield, ik,  1 )
           !
        END DO
        !
        DEALLOCATE( psi )
        !
        !...set up electric field hermitean operator
        !
        CALL h_epsi_her_set ( )
        !
        CALL c_bands( iter, ik_, dr2 )
        !
     END DO
     !
     DEALLOCATE(fact_hepsi)
     !
  ELSE
     !
     CALL c_bands( iter, ik_, dr2 )
     !
  END IF
  !
  CALL poolrecover( et, nbnd, nkstot, nks )
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  WRITE( stdout, 9102 )
  !
  ! ... write band eigenvalues
  !
  DO ik = 1, nkstot
     !
     IF ( lsda ) THEN
        !   
        IF ( ik == 1 ) WRITE( stdout, 9015 )
        IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016 )
        !
     END IF
     !
     WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
     WRITE( stdout, 9030 ) ( et(ibnd,ik)*rytoev, ibnd = 1, nbnd )
     !
  END DO
  !
  IF ( lgauss .AND. .NOT. lbands ) THEN
     !
     ef = efermig( et, nbnd, nks, nelec, wk, degauss, ngauss, 0, isk )
     !
     WRITE( stdout, 9040 ) ef*rytoev
     !
  ELSE IF ( ltetra .AND. .NOT. lbands ) THEN
     !
     ef = efermit( et, nbnd, nks, nelec, nspin, ntetra, tetra, 0, isk )
     !
     WRITE( stdout, 9040 ) ef*rytoev
     !
  ELSE IF (.NOT. lbands) THEN
     !
     IF ( tfixed_occ ) THEN
        ibnd = 0
        DO kbnd = 1, nbnd
           IF ( nspin == 1 .OR. nspin == 4 ) THEN
              IF ( f_inp(kbnd,1) > 0.D0 ) ibnd = kbnd
           ELSE
              IF ( f_inp(kbnd,1) > 0.D0 ) ibnd_up   = kbnd
              IF ( f_inp(kbnd,2) > 0.D0 ) ibnd_dw = kbnd
           END IF
        END DO
     ELSE
        IF ( nspin == 1 ) THEN
           ibnd = NINT( nelec ) / 2
        ELSE
           ibnd = NINT( nelec )
        END IF
     END IF
     !
     IF ( ionode .AND. nbnd > ibnd ) THEN
        !
        IF ( nspin == 1 .OR. nspin == 4 ) THEN
           ehomo = MAXVAL( et(ibnd  ,1:nkstot) )
           elumo = MINVAL( et(ibnd+1,1:nkstot) )
           !
        ELSE
           ehomo = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
           elumo = MIN( MINVAL( et(ibnd_up+1,1:nkstot/2) ), &
                MINVAL( et(ibnd_dw+1,nkstot/2+1:nkstot) ) )
        END IF
        IF ( ehomo < elumo ) &
             WRITE( stdout, 9042 ) ehomo*rytoev, elumo*rytoev
        !
     END IF
     !
  END IF
  !
  CALL flush_unit( stdout )
  !
  ! ... do a Berry phase polarization calculation if required
  !
  IF ( lberry ) CALL c_phase()
  !
  CALL stop_clock( 'electrons' )
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F9.2,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9015 FORMAT(/' ------ SPIN UP ------------'/ )
9016 FORMAT(/' ------ SPIN DOWN ----------'/ )
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
9030 FORMAT( '  ',8F9.4 )
9040 FORMAT(/'     the Fermi energy is ',F10.4,' ev' )
9042 FORMAT(/'     Highest occupied, lowest unoccupied level (ev): ',2F10.4 )
9102 FORMAT(/'     End of band structure calculation' )
  !
END SUBROUTINE non_scf
