!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
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
  USE ener,                 ONLY : ef , demet  
  USE constants,            ONLY : tpi, rytoev
  USE cell_base,            ONLY : tpiba2
  USE basis,                ONLY : natomwfc, startingwfc
  USE gvect,                ONLY : g
  USE klist,                ONLY : lgauss, degauss, ngauss, xk, nks, &
                                   wk, nkstot, nelec
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE control_flags,        ONLY : isolve, iprint, reduce_io
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, wg
  USE us,                   ONLY : okvan
  USE fixed_occ,            ONLY : f_inp, tfixed_occ
  USE ktetra,               ONLY : ltetra, ntetra, tetra
  USE uspp,                 ONLY : nkb, vkb
  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE io_files,             ONLY : iunat, nwordwfc, iunwfc, iunigk, &
                                   nwordatwfc
  USE wavefunctions_module, ONLY : evc
  USE mp_global,            ONLY : intra_image_comm, me_image, root_image
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER  :: ik, is, ibnd
  !
  !
  CALL start_clock( 'wfcinit' )
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
  demet = 0.D0
  !
  CALL poolrecover( et, nbnd, nkstot, nks )
  !
  ! ... occupations are computed here
  !
  IF ( .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) THEN
     !
     ! ... calculate weights for the insulator case
     !
     CALL iweights( nks, wk, nbnd, nelec, et, ef, wg )
     !
  ELSE IF ( ltetra ) THEN
     !
     ! ... calculate weights for the metallic case
     !
     CALL poolrecover( et, nbnd, nkstot, nks )
     !
     IF ( me_image == root_image ) THEN
        !
        CALL tweights( nkstot, nspin, nbnd, nelec, ntetra, tetra, et, ef, wg )
        !
     END IF
     !
     CALL poolscatter( nbnd, nkstot, wg, nks, wg )
     !
     CALL mp_bcast( ef, root_image, intra_image_comm )
     !
  ELSE IF ( lgauss ) THEN
     !
     CALL gweights( nks, wk, nbnd, nelec, degauss, ngauss, et, ef, demet, wg )
     !
  ELSE IF ( tfixed_occ ) THEN
     !
     ef = - 1.0D+20
     !
     wg = f_inp
     !
     DO is = 1, nspin
        !
        DO ibnd = 1, nbnd
           !
           IF ( wg(ibnd,is) > 0.D0 ) ef = MAX( ef, et(ibnd,is) )
           !
        END DO
        !
     END DO
     !
  END IF  
  !
  IF ( iprint == 1 ) THEN
     !
     DO ik = 1, nkstot
        !
        WRITE( stdout, &
               '(/,10X,"k =",3F7.4,5X,"band energies (ev):"/)' ) xk(:,ik)
        WRITE( stdout, '(2X,8F9.4)') et(:,ik) * rytoev
        !
     END DO
     !
  END IF
  !
#if defined (FLUSH)
       CALL flush( stdout )
#endif
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
       INTEGER                        :: ibnd, ig, ipol, n_starting_wfc
         ! counter on k points
         !    "     "   bands
         !    "     "  plane waves
         !    "     "  polarization
         ! number of starting wavefunctions
       COMPLEX(KIND=DP), ALLOCATABLE  :: wfcatom(:,:)
         ! atomic wfcs for initialization
       REAL(KIND=DP)                  :: rr, arg
       REAL(KIND=DP)                  :: rndm
       EXTERNAL                          rndm
         ! random function generation
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
       IF ( startingwfc == 'atomic' ) THEN
          !
          IF ( natomwfc >= nbnd ) THEN
             WRITE( stdout, '(5X,A)' ) 'Starting wfc are atomic'
          ELSE
             WRITE( stdout, '(5X,A,I3,A)' ) &
                 'Starting wfc are atomic + ', nbnd-natomwfc, ' random wfc'
          END IF
          !
          n_starting_wfc = MAX( natomwfc, nbnd )
          !
       ELSE
          !
          WRITE( stdout, '(5X,A)' ) 'Starting wfc are random'
          !
          n_starting_wfc = nbnd
          !
       END IF
       !
       ! ... Needed for LDA+U (not yet in gamma)
       !
       IF ( lda_plus_u ) CALL orthoatwfc()
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       ! ... we start a loop on k points
       !
       ! ... rbecp contains <beta|psi> - used in h_psi and s_psi
       !
       ALLOCATE( rbecp( nkb, n_starting_wfc ) )
       ALLOCATE( wfcatom( npwx, n_starting_wfc ) )
       !
       DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          IF ( nks > 1 ) READ( iunigk ) npw, igk
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
          IF ( lda_plus_u ) &
               CALL davcio( swfcatom, nwordatwfc, iunat, ik, - 1 )
          !
          IF ( startingwfc == 'atomic' ) THEN
             !
             CALL atomic_wfc( ik, wfcatom )
             !
             ! ... if not enough atomic wfc are available, complete 
             ! ... with random wfcs
             !
             DO ibnd = ( natomwfc + 1 ), nbnd
                DO ig = 1, npw
                   rr = rndm()
                   arg = tpi * rndm()
                   wfcatom(ig,ibnd) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                      ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                        ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                        ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.D0 )
                END DO
             END DO
          ELSE
             DO ibnd = 1, nbnd
                DO ig = 1, npw
                   rr = rndm ()
                   arg = tpi * rndm ()
                   wfcatom(ig,ibnd) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                      ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                        ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                        ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.D0 )
                END DO
             END DO
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
          !
          IF ( isolve == 1 ) THEN
             !
             CALL cinitcgg( npwx, npw, n_starting_wfc, &
                            nbnd, wfcatom, evc, et(1,ik) )
             !
          ELSE
             !
             CALL rotate_wfc_gamma( npwx, npw, n_starting_wfc, gstart, &
                                    nbnd, wfcatom, okvan, evc, et(1,ik) )
          END IF
          !
          evc(npw+1:npwx,:) = ( 0.D0, 0.D0 )
          !
          IF ( nks > 1 .OR. .NOT. reduce_io ) &
             CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
          !
       END DO
       !
       DEALLOCATE( rbecp )
       DEALLOCATE( wfcatom )
       !
       RETURN
       !
     END SUBROUTINE wfcinit_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE wfcinit_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version       
       !
       USE becmod, ONLY : becp
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: &
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
       ! ... becp contains <beta|psi> - used in h_psi and s_psi
       !
       ALLOCATE( becp( nkb, n_starting_wfc ) )
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
          IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunat, ik, -1 )
          !
          IF ( startingwfc == 'atomic' ) THEN
             !
             CALL atomic_wfc( ik, wfcatom )
             !
             ! ... if not enough atomic wfc are available, complete 
             ! ... with random wfcs
             !
             DO ibnd = ( natomwfc + 1 ), nbnd
                DO ig = 1, npw
                   rr = rndm()
                   arg = tpi * rndm()
                   wfcatom(ig,ibnd) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                      ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                        ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                        ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.D0 )
                END DO
             END DO
             !
          ELSE
             DO ibnd = 1, nbnd
                DO ig = 1, npw
                   rr = rndm()
                   arg = tpi * rndm()
                   wfcatom(ig,ibnd) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
                                      ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                                        ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                                        ( xk(3,ik) + g(3,igk(ig)) )**2 + 1.0D0 )
                END DO
             END DO
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... Diagonalize the Hamiltonian on the basis of atomic wfcs
          !
          IF ( isolve == 1 ) THEN
             CALL cinitcgg( npwx, npw, n_starting_wfc, &
                            nbnd, wfcatom, evc, et(1,ik) )
          ELSE
             CALL rotate_wfc( npwx, npw, n_starting_wfc, &
                              nbnd, wfcatom, okvan, evc, et(1,ik) )
          END IF
          !
          evc(npw+1:npwx,:) = ( 0.D0, 0.D0 )
          !
          IF ( nks > 1 .OR. .NOT. reduce_io ) &
               CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
          !
       END DO
       !
       DEALLOCATE( becp )
       DEALLOCATE( wfcatom )
       !
       RETURN
       !
     END SUBROUTINE wfcinit_k
     !
END SUBROUTINE wfcinit
