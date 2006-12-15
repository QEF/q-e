!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
! 
!----------------------------------------------------------------------------
SUBROUTINE wfcinit()
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes an estimate of the starting wavefunctions
  ! ... from superposition of atomic wavefunctions.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : gamma_only
  USE constants,            ONLY : tpi
  USE cell_base,            ONLY : tpiba2
  USE basis,                ONLY : natomwfc, startingwfc
  USE gvect,                ONLY : g
  USE klist,                ONLY : xk, nks
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE control_flags,        ONLY : isolve, reduce_io
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et,&
                                   wg, current_k
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_files,             ONLY : iunsat, nwordwfc, iunwfc, iunigk, &
                                   nwordatwfc
  USE wavefunctions_module, ONLY : evc
  USE random_numbers,       ONLY : rndm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, is, ibnd
  !
  !
  CALL start_clock( 'wfcinit' )
  !
  ! ... Needed for LDA+U (not yet in gamma)
  !
  IF ( lda_plus_u ) CALL orthoatwfc()
  !
  IF ( gamma_only ) THEN
     !
     CALL wfcinit_gamma()
     !
  ELSE
     !
     CALL wfcinit_k()
     !
  END IF  
  !
  CALL stop_clock( 'wfcinit' )  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE wfcinit_gamma()
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version       
       !
       USE gvect,  ONLY : gstart
       USE becmod, ONLY : rbecp
       !
       IMPLICIT NONE
       !       
       INTEGER :: ibnd, ig, n_starting_wfc
         ! counter on  bands
         !    "     "  plane waves
         ! number of starting wavefunctions
       COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)
         ! atomic wfcs for initialization
       REAL(DP), ALLOCATABLE :: etatom(:)
         ! atomic eigenvalues
       REAL(DP) :: rr, arg
       !
       !
       ! ... state what is going to happen
       !
       IF ( startingwfc == 'file' ) THEN
          !
          WRITE( stdout, '(5X,"Starting wfc from file")' )
          !
          ! ... read the wavefunction into memory (if it is not done in c_bands)
          !
          IF ( nks == 1 .AND. reduce_io ) &
             CALL davcio( evc, nwordwfc, iunwfc, 1, -1 )
          !
          RETURN
          !
       ELSE IF ( startingwfc == 'atomic' ) THEN
          !
          IF ( natomwfc >= nbnd ) THEN
             !
             WRITE( stdout, '(5X,"Starting wfc are atomic")' ) 
             !
          ELSE
             !
             WRITE( stdout, &
                  & '(5X,"Starting wfc are atomic + ",I3," random wfc")' ) &
                 nbnd-natomwfc
             !
          END IF
          !
          n_starting_wfc = MAX( natomwfc, nbnd )
          !
       ELSE
          !
          WRITE( stdout, '(5X,"Starting wfc are random")' )
          !
          n_starting_wfc = nbnd
          !
       END IF
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       ! ... we start a loop on k points
       !
       ! ... rbecp contains <beta|psi> - used in h_psi and s_psi
       !
       ALLOCATE( rbecp( nkb, n_starting_wfc ) )
       ALLOCATE( wfcatom( npwx, n_starting_wfc ) )
       ALLOCATE( etatom( n_starting_wfc ) )
       !
       wfcatom (:,:) = (0.D0, 0.d0)
       !       
       DO ik = 1, nks
          !
          current_k = ik
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          IF ( nks > 1 ) READ( iunigk ) npw, igk
          !
          ! ... here we compute the kinetic energy
          !
          DO ig = 1, npw
             !
             g2kin(ig) = ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                         ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                         ( xk(3,ik) + g(3,igk(ig)) )**2
             !
          END DO
          !
          ! ... Put the correct units on the kinetic energy
          !
          g2kin(:) = g2kin(:) * tpiba2
          !
          IF ( lda_plus_u ) &
             CALL davcio( swfcatom, nwordatwfc, iunsat, ik, - 1 )
          !
          IF ( startingwfc == 'atomic' ) THEN
             !
             CALL atomic_wfc( ik, wfcatom )
             !
             ! ... if not enough atomic wfc are available, complete 
             ! ... with random wfcs
             !
             DO ibnd = ( natomwfc + 1 ), nbnd
                !
                DO ig = 1, npw
                   !
                   rr  = rndm()
                   arg = tpi * rndm()
                   !
                   wfcatom(ig,ibnd) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                      ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                        ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                        ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.D0 )
                   !
                END DO
                !
             END DO
             !
          ELSE
             !
             DO ibnd = 1, nbnd
                !
                DO ig = 1, npw
                   !
                   rr  = rndm()
                   arg = tpi * rndm()
                   !
                   wfcatom(ig,ibnd) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                      ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                        ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                        ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.D0 )
                   !
                END DO
                !
             END DO
             !
          END IF
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
          !
          IF ( isolve == 1 ) THEN
             !
             CALL rinitcgg( npwx, npw, n_starting_wfc, &
                            nbnd, wfcatom, wfcatom, etatom )
             !
          ELSE
             !
             CALL rotate_wfc_gamma( npwx, npw, n_starting_wfc, gstart, &
                                    nbnd, wfcatom, okvan, wfcatom, etatom )
          END IF
          !
          ! ... the first nbnd wavefunctions and eigenvalues are copied
          !
          et(1:nbnd,ik) = etatom(1:nbnd)
          !
          evc(:,1:nbnd) = wfcatom(:,1:nbnd)
          !
          IF ( nks > 1 .OR. .NOT. reduce_io ) &
             CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
          !
       END DO
       !
       DEALLOCATE( rbecp )
       DEALLOCATE( wfcatom )
       DEALLOCATE( etatom )
       !
       RETURN
       !
     END SUBROUTINE wfcinit_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE wfcinit_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version       
       !
       USE becmod, ONLY : becp, becp_nc
       !
       IMPLICIT NONE
       !
       INTEGER :: ibnd, ig, ipol, ig2, n_starting_wfc
         ! counter on  bands
         !    "     "  plane waves
         !    "     "  polarization
         ! number of starting wavefunctions
       COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)
         !  atomic wfcs for initialization
       REAL(DP), ALLOCATABLE :: etatom(:)
         ! atomic eigenvalues
       REAL (DP) :: rr, arg
       !
       !
       ! ... state what is going to happen
       !
       IF ( startingwfc == 'file' ) THEN
          !
          WRITE( stdout, '(5X,"Starting wfc from file")' )
          !
          ! ... read the wavefunction into memory (if it is not done in c_bands)
          !
          IF ( nks == 1 .AND. reduce_io ) THEN
             !
             CALL davcio( evc, nwordwfc, iunwfc, 1, -1 )
             !
          ENDIF
          !
          RETURN
          !
       ELSE IF ( startingwfc == 'atomic' ) THEN
          !
          IF ( natomwfc >= nbnd ) THEN
             !
             WRITE( stdout, '(5X,"Starting wfc are atomic")' ) 
             !
          ELSE
             !
             WRITE( stdout, &
                  & '(5X,"Starting wfc are atomic + ",I3," random wfc")' ) &
                 nbnd-natomwfc
             !
          END IF
          !
          n_starting_wfc = MAX( natomwfc, nbnd )
          !
       ELSE
          !
          WRITE( stdout, '(5X,"Starting wfc are random")' )
          !
          n_starting_wfc = nbnd
          !
       END IF
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       ! ... we start a loop on k points
       !
       ! ... becp contains <beta|psi> - used in h_psi and s_psi
       !
       IF ( noncolin ) THEN
          !
          ALLOCATE( becp_nc( nkb, npol, n_starting_wfc ) )
          !
       ELSE
          !
          ALLOCATE( becp( nkb, n_starting_wfc ) )
          !
       END IF
       !
       ALLOCATE( wfcatom( npwx*npol, n_starting_wfc ) )
       !
       ALLOCATE( etatom( n_starting_wfc ) )
       !
       DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          IF ( nks > 1 ) READ( UNIT = iunigk ) npw, igk
          !
          ! ... here we compute the kinetic energy
          !
          DO ig = 1, npw
             !
             g2kin(ig) = ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                         ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                         ( xk(3,ik) + g(3,igk(ig)) )**2
             !
          END DO
          !
          ! ... Put the correct units on the kinetic energy
          !
          g2kin(:) = g2kin(:) * tpiba2
          !
          IF ( lda_plus_u ) &
             CALL davcio( swfcatom, nwordatwfc, iunsat, ik, -1 )
          !
          IF ( startingwfc == 'atomic' ) THEN
             !
             IF ( noncolin ) THEN
                !
                CALL atomic_wfc_nc( ik, wfcatom )
                !
             ELSE
                !
                CALL atomic_wfc( ik, wfcatom )
                !
             END IF
             !
             ! ... if not enough atomic wfc are available, complete 
             ! ... with random wfcs
             !
             DO ibnd = ( natomwfc + 1 ), nbnd
                !
                DO ipol = 1, npol
                   !
                   DO ig = 1, npw
                      !
                      ig2 = ig + (ipol-1)*npwx
                      rr  = rndm()
                      arg = tpi * rndm()
                      !
                      wfcatom(ig2,ibnd) = &
                              CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                    ( ( xk(1,ik) + g(1,igk(ig) ) )**2 + &
                                      ( xk(2,ik) + g(2,igk(ig) ) )**2 + &
                                      ( xk(3,ik) + g(3,igk(ig) ) )**2 + 1.0D0 )
                   END DO
                   !
                END DO
                !
             END DO
             !
          ELSE
             !
             DO ibnd = 1, nbnd
                !
                DO ipol = 1, npol
                   !
                   DO ig = 1, npw
                      !
                      ig2 = ig + (ipol-1)*npwx
                      rr  = rndm()
                      arg = tpi * rndm()
                      !
                      wfcatom(ig2,ibnd) = &
                              CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                     ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                       ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                       ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.D0 )
                   END DO
                   !
                END DO
                !
             END DO
             !
          END IF
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
          !
          IF ( isolve == 1 ) THEN
             !
             CALL cinitcgg( npwx, npw, n_starting_wfc, &
                            nbnd, wfcatom, evc, etatom, .true. )
             !
          ELSE
             !
             IF ( noncolin ) THEN
                !
                CALL rotate_wfc_nc( npwx, npw, n_starting_wfc, nbnd, &
                                   wfcatom, npol, okvan, evc, etatom )
                !
             ELSE
                !
                CALL rotate_wfc( npwx, npw, n_starting_wfc, &
                                 nbnd, wfcatom, okvan, evc, etatom )
                !
             END IF
             !
             evc(npw+1:npwx,1:nbnd) = (0.d0, 0.d0)
             !
          END IF
          !
          ! ... the first nbnd wavefunctions and eigenvalues are copied
          !
          et(1:nbnd,ik) = etatom(1:nbnd)
          !
          evc(npw+1:npwx,1:nbnd) = (0.d0, 0.d0)
          IF ( noncolin ) THEN
             evc(npwx+npw+1:2*npwx,1:nbnd) = (0.d0,0.d0)
          END IF
          !
          IF ( nks > 1 .OR. .NOT. reduce_io ) &
             CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
          !
       END DO
       !
       IF ( noncolin ) THEN
          !
          DEALLOCATE( becp_nc )
          !
       ELSE
          !
          DEALLOCATE( becp )
          !
       END IF
       !
       DEALLOCATE( wfcatom )
       DEALLOCATE( etatom )
       !
       RETURN
       !
     END SUBROUTINE wfcinit_k
     !
END SUBROUTINE wfcinit
