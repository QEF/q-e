!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE wfcinit()
  !----------------------------------------------------------------------------
  !
  !     This routine computes an estimate of the starting wavefunctions
  !     from superposition of atomic wavefunctions.
  !
  USE parameters,            ONLY :  DP
  USE io_global,             ONLY :  stdout
  USE constants,             ONLY :  tpi, rytoev
  USE brilz,                 ONLY :  tpiba2
  USE basis,                 ONLY :  natomwfc, startingwfc
  USE gvect,                 ONLY :  g
  USE klist,                 ONLY :  xk, nks, nkstot
  USE lsda_mod,              ONLY :  lsda, current_spin, isk
  USE varie,                 ONLY :  isolve, iprint, reduce_io
  USE wvfct,                 ONLY :  nbnd, npw, npwx, igk, g2kin, et
  USE us,                    ONLY :  vkb, okvan
  USE ldaU,                  ONLY :  swfcatom, lda_plus_u
  USE io_files,              ONLY :  iunat, nwordwfc, iunwfc, iunigk, &
                                     nwordatwfc
  USE wavefunctions_module,  ONLY :  evc
  !
  IMPLICIT NONE
  !
  INTEGER :: &
      ik,               &!  counter on k points
      ibnd,             &!     "     "   bands
      ig,               &!     "     "   plane waves
      ipol,             &!     "     "   polarization
      n_starting_wfc     !  number of starting wavefunctions
  COMPLEX(KIND=DP), ALLOCATABLE :: & 
      wfcatom(:,:)       !  atomic wfcs for initialization
  REAL (KIND=DP) :: &
      rr,               &!
      arg                !
  REAL (KIND=DP), EXTERNAL :: &
      rndm               !  random function generation
  !
  !
  ! ... state what is going to happen
  !
  IF ( startingwfc == 'file' ) THEN
     !
     WRITE( stdout, '(5X,A)' ) 'Starting wfc from file'
     !
     ! ... read the wavefunction into memory (if it is not done in c_bands)
     !
     IF ( nks == 1 .AND. reduce_io ) &
        CALL davcio( evc, nwordwfc, iunwfc, 1, -1 )
     !
     RETURN
     !
  END IF
  !
  CALL start_clock( 'wfcinit' )
  !
  IF ( startingwfc == 'atomic' ) THEN
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,A)' ) 'Starting wfc are atomic'
     ELSE
        WRITE( stdout, '(5X,A,I3,A)') &
            'Starting wfc are atomic + ',nbnd-natomwfc, ' random wfc'
     END IF
     n_starting_wfc = MAX( natomwfc, nbnd )
  ELSE
     WRITE( stdout, '(5X,A)' ) 'Starting wfc are random'
     n_starting_wfc = nbnd
  END IF
  !
  ! ... Needed for LDA+U
  !
  IF ( lda_plus_u ) CALL orthoatwfc
  IF ( nks > 1 ) REWIND( iunigk )
  !
  ! ... we start a loop on k points
  !
  ALLOCATE( wfcatom(npwx,n_starting_wfc) )
  !
  DO ik = 1, nks
     IF ( lsda ) current_spin = isk(ik)
     IF ( nks > 1 ) READ( UNIT = iunigk ) npw, igk
     !
     ! ... here we compute the kinetic energy
     !
     DO ig = 1, npw
        g2kin(ig) = ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                    ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                    ( xk(3,ik) + g(3,igk(ig)) )**2
     END DO
     !
     ! ... Put the correct units on the kinetic energy
     !
     g2kin(:) = g2kin(:) * tpiba2
     !
     IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunat, ik, - 1 )
     !
     IF ( startingwfc == 'atomic' ) THEN
        !
        CALL atomic_wfc( ik, wfcatom )
        !
        ! ... if not enough atomic wfc are available, complete with random wfcs
        !
        DO ibnd = ( natomwfc + 1 ), nbnd
           DO ig = 1, npw
              rr = rndm()
              arg = tpi * rndm()
              wfcatom(ig,ibnd) = CMPLX( rr * COS( arg ), rr * SIN( arg ) ) / &
                                 ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                   ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                   ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.0D0 )
           END DO
        END DO
        !
     ELSE
        DO ibnd = 1, nbnd
           DO ig = 1, npw
              rr = rndm()
              arg = tpi * rndm()
              wfcatom(ig,ibnd) = CMPLX( rr * COS( arg ), rr * SIN( arg ) ) / &
                                 ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                   ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                   ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.0D0 )
           END DO
        END DO
        !
     END IF
     !
     CALL init_us_2( npw, igk, xk(1,ik), vkb )
     !
     ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
     !
     IF ( isolve == 1 ) THEN
        CALL cinitcgg( npwx, npw, n_starting_wfc, nbnd, &
                       wfcatom, evc, et(1,ik) )
     ELSE
        CALL rotate_wfc( npwx, npw, n_starting_wfc, nbnd, wfcatom, &
                         okvan, evc, et(1,ik) )
     END IF
     !
     DO ibnd = 1, nbnd
        DO ig = npw + 1, npwx
           evc(ig,ibnd) = ( 0.D0, 0.D0 )
        END DO
     END DO
     !
     IF ( nks > 1 .OR. .NOT. reduce_io ) &
        CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
     !
  END DO
  !
  DEALLOCATE( wfcatom )
  !
  IF ( iprint == 1 ) THEN
#ifdef __PARA
     CALL poolrecover( et, nbnd, nkstot, nks )
#endif
     DO ik = 1, nkstot
        WRITE( stdout, 9010 ) ( xk(ipol,ik), ipol = 1, 3 )
        WRITE( stdout, '(2X,8F9.4)' ) ( et(ibnd,ik) * rytoev, ibnd = 1, nbnd )
     END DO
  END IF
#ifdef FLUSH
  CALL flush( 6 )
#endif
  CALL stop_clock( 'wfcinit' )
  !
  RETURN
  !
9010 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/)
  !
END SUBROUTINE wfcinit

