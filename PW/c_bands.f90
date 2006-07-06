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
SUBROUTINE c_bands( iter, ik_, dr2 )
  !----------------------------------------------------------------------------
  !
  ! ... this is a wrapper to specific calls
  !
  ! ... internal procedures :
  !
  ! ... c_bands_gamma()   : for gamma sampling of the BZ (optimized algorithms)
  ! ... c_bands_k()       : for arbitrary BZ sampling (general algorithm)
  ! ... test_exit_cond()  : the test on the iterative diagonalization
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : gamma_only
  USE io_files,             ONLY : iunigk, nwordatwfc, iunat, iunwfc, nwordwfc, iunefield
  USE cell_base,            ONLY : tpiba2 
  USE klist,                ONLY : nkstot, nks, wk, xk, nelec
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE gvect,                ONLY : g, gstart, ecfixed, qcutz, q2sigma, nrxx, &
                                   nr1, nr2, nr3  
  USE wvfct,                ONLY : g2kin, wg, nbndx, et, nbnd, npwx, igk, &
                                   npw, current_k, btype
  USE control_flags,        ONLY : diis_ndim, ethr, lscf, max_cg_iter, &
                                   isolve, istep, reduce_io, conv_ions
  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE scf,                  ONLY : vltot
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc, evc_nc  
  USE g_psi_mod,            ONLY : h_diag, s_diag, h_diag_nc, s_diag_nc
  USE bp,                   ONLY : lelfield, evcel
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variables
  !
  INTEGER :: ik_, iter
    ! k-point already done
    ! current iterations
  REAL(DP) :: dr2
    ! current accuracy of self-consistency
  !
  ! ... local variables
  !
  REAL(DP) :: avg_iter, cg_iter, v_of_0
    ! average number of iterations
    ! number of iterations in Conjugate-Gradient
    ! the average of the potential
  INTEGER :: ik, ig, ibnd, dav_iter, diis_iter, ntry, notconv
    ! counter on k points
    ! counter on G vectors
    ! counter on bands
    ! number of iterations in Davidson
    ! number of iterations in DIIS
    ! number or repeated call to diagonalization in case of non convergence
    ! number of notconverged elements
  LOGICAL :: lrot
    ! .TRUE. if the wfc have already be rotated
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: erf
    ! error function  
  !
  IF ( ik_ == nks ) THEN
     !
     ik_ = 0
     !
     RETURN
     !
  END IF
  !
  CALL start_clock( 'c_bands' )
  !
  IF ( noncolin ) THEN
     !
     ALLOCATE( h_diag_nc( npwx, npol ) )
     ALLOCATE( s_diag_nc( npwx, npol ) )
     !
  ELSE
     !
     ALLOCATE( h_diag( npwx ) )
     ALLOCATE( s_diag( npwx ) )
     !
  END IF
  !
  IF ( lelfield ) ALLOCATE( evcel( npwx, nbnd ) )
  !
  IF ( gamma_only ) THEN
     !
     CALL c_bands_gamma()
     !
  ELSE
     !
     CALL c_bands_k()
     !
  END IF
  !
  IF ( noncolin ) THEN
     !
     DEALLOCATE( s_diag_nc )
     DEALLOCATE( h_diag_nc )
     !
  ELSE
     !
     DEALLOCATE( s_diag )
     DEALLOCATE( h_diag )
     !
  END IF
  !
  IF ( lelfield ) DEALLOCATE( evcel )
  !
  CALL stop_clock( 'c_bands' )  
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE c_bands_gamma()
       !-----------------------------------------------------------------------
       !  
       ! ... This routine is a driver for the diagonalization routines of the
       ! ... total Hamiltonian at Gamma point only
       ! ... It reads the Hamiltonian and an initial guess of the wavefunctions
       ! ... from a file and computes initialization quantities for the
       ! ... diagonalization routines.
       ! ... There are two types of iterative diagonalization:
       ! ... a) Davidson algorithm (all-band)
       ! ... b) Conjugate Gradient (band-by-band)
       ! ... c) DIIS algorithm (all-band) 
       !
       USE becmod,           ONLY : rbecp
       USE real_diis_module, ONLY : rdiisg
       !
       IMPLICIT NONE
       !     
       !
       ! ... becp contains <beta|psi> - used in h_psi and s_psi
       !
       ALLOCATE( rbecp( nkb, nbnd ) )
       !
       IF ( isolve == 0 ) THEN
          !
          WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
          !
       ELSE IF ( isolve == 1 ) THEN
          !
          WRITE( stdout, '(5X,"CG style diagonalization")')
          !
       ELSE IF ( isolve == 2 ) THEN
          !
          WRITE( stdout, '(5X,"DIIS style diagonalization")')
          !
       END IF
       !
       avg_iter = 0.D0
       !
       ! ... v_of_0 is (Vloc)(G=0)
       !
       v_of_0 = SUM( vltot(1:nrxx) ) / DBLE( nr1 * nr2 * nr3 )
       !
       CALL reduce( 1, v_of_0 )
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       ! ... For each k point diagonalizes the hamiltonian
       !
       k_loop: DO ik = 1, nks
          !
          current_k = ik
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          ! ... Reads the Hamiltonian and the list k+G <-> G of this k point
          !
          IF ( nks > 1 ) READ( iunigk ) npw, igk
          !
          ! ... do not recalculate k-points if restored from a previous run
          !
          IF ( ik <= ik_ ) THEN
             !
             CALL save_in_cbands( iter, ik, dr2 )
             !
             CYCLE k_loop
             !
          END IF          
          !
          ! ... various initializations
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... read in wavefunctions from the previous iteration
          !
          IF ( nks > 1 .OR. .NOT. reduce_io ) &
             CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
          !
          ! ... Needed for LDA+U
          !
          IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunat, ik, -1 )
          !
          ! ... sets the kinetic energy
          !
          g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk(1:npw)) )**2 + &
                           ( xk(2,ik) + g(2,igk(1:npw)) )**2 + &
                           ( xk(3,ik) + g(3,igk(1:npw)) )**2 ) * tpiba2
          !
          IF ( qcutz > 0.D0 ) THEN
             !
             DO ig = 1, npw
                !
                g2kin(ig) = g2kin(ig) + qcutz * &
                            ( 1.D0 + erf( (g2kin(ig) - ecfixed ) / q2sigma ) )
                !
             END DO
             !
          END IF
          !
          IF ( isolve == 1 ) THEN
             !
             ! ... Conjugate-Gradient diagonalization
             !
             ! ... h_diag is the precondition matrix
             !
             FORALL( ig = 1 : npw )
                !
                h_diag(ig) = 1.D0 + g2kin(ig) + &
                             SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
                !
             END FORALL
             !
             ntry = 0
             !
             CG_loop : DO
                !
                IF ( iter > 1 .OR. istep > 0 .OR. ntry > 0 ) THEN
                   !
                   CALL rinitcgg( npwx, npw, nbnd, nbnd, evc, evc, et(1,ik) )
                   !
                   avg_iter = avg_iter + 1.D0
                   !
                END IF
                !
                CALL rcgdiagg( npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
                               h_diag, ethr, max_cg_iter, .NOT. lscf, &
                               notconv, cg_iter )
                !
                avg_iter = avg_iter + cg_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
                IF ( nks > 1 .OR. .NOT. reduce_io ) &
                   CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
                !
                ntry = ntry + 1                
                !
                ! ... exit condition
                !
                IF ( test_exit_cond() ) EXIT  CG_loop
                !
             END DO CG_loop
             !
          ELSE IF ( isolve == 2 ) THEN
             !
             ! ... RMM-DIIS method
             !
             h_diag(1:npw) = g2kin(1:npw) + v_of_0
             !
             CALL usnldiag( h_diag, s_diag )
             !
             ntry = 0
             !
             RMMDIIS_loop: DO
                !
                CALL rdiisg( npw, npwx, nbnd, evc, &
                             et(1,ik), btype(1,ik), notconv, diis_iter, iter )
                !  
                avg_iter = avg_iter + diis_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
                IF ( nks > 1 .OR. .NOT. reduce_io ) &
                   CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
                !
                ntry = ntry + 1                
                !
                ! ... exit condition
                !
                IF ( test_exit_cond() ) EXIT  RMMDIIS_loop
                !
             END DO RMMDIIS_loop
             !
          ELSE
             !
             ! ... Davidson diagonalization
             !
             ! ... h_diag are the diagonal matrix elements of the 
             ! ... hamiltonian used in g_psi to evaluate the correction 
             ! ... to the trial eigenvectors
             !
             h_diag(1:npw) = g2kin(1:npw) + v_of_0
             !
             CALL usnldiag( h_diag, s_diag )
             !
             ntry = 0
             !
             david_loop: DO
                !
                lrot = ( iter == 1 )
                !
                CALL regterg( npw, npwx, nbnd, nbndx, evc, ethr, &
                              okvan, gstart, et(1,ik), btype(1,ik), &
                              notconv, lrot, dav_iter )
                !
                avg_iter = avg_iter + dav_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
                IF ( nks > 1 .OR. .NOT. reduce_io ) &
                   CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
                !   
                ntry = ntry + 1
                !
                ! ... exit condition
                !
                IF ( test_exit_cond() ) EXIT  david_loop
                !
             END DO david_loop
             !
          END IF
          !
          IF ( notconv > MAX( 5, nbnd / 4 ) ) THEN
             !
             CALL errore( 'c_bands', &
                        & 'too many bands are not converged', 1 )
             !
          END IF
          !
          ! ... save restart information
          !
          CALL save_in_cbands( iter, ik, dr2 )
          !
       END DO k_loop
       !
       ik_ = 0
       !
       CALL poolreduce( 1, avg_iter )
       !
       avg_iter = avg_iter / nkstot
       !
       WRITE( stdout, &
              '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
           ethr, avg_iter
       !
       ! ... deallocate work space
       !
       DEALLOCATE( rbecp )
       !
       RETURN
       !
     END SUBROUTINE c_bands_gamma
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE c_bands_k()
       !-----------------------------------------------------------------------
       !
       ! ... This routine is a driver for the diagonalization routines of the
       ! ... total Hamiltonian at each k-point.
       ! ... It reads the Hamiltonian and an initial guess of the wavefunctions
       ! ... from a file and computes initialization quantities for the
       ! ... diagonalization routines.
       ! ... There are three types of iterative diagonalization:
       ! ... a) Davidson algorithm (all-band)
       ! ... b) Conjugate Gradient (band-by-band)
       ! ... c) DIIS algorithm
       !
       USE becmod,              ONLY : becp, becp_nc
       USE complex_diis_module, ONLY : cdiisg
       USE control_flags,       ONLY : lbands
       USE check_stop,          ONLY : check_stop_now
       !
       IMPLICIT NONE
       !
       ! ... here the local variables
       !
       INTEGER :: ipol
       !
       ! ... becp contains <beta|psi> - used in h_psi and s_psi
       !
       IF ( noncolin ) THEN
          !
          ALLOCATE( becp_nc( nkb, npol, nbnd ) )
          !
       ELSE
          !
          ALLOCATE( becp( nkb, nbnd ) )
          !
       END IF
       !
       IF ( isolve == 0 ) THEN
          !
          WRITE( stdout, '(5X,"Davidson diagonalization with overlap")')
          !
       ELSE IF ( isolve == 1 ) THEN
          !
          WRITE( stdout, '(5X,"CG style diagonalization")')
          !
       ELSE IF ( isolve == 2 ) THEN
          !
          WRITE( stdout, '(5X,"DIIS style diagonalization")')
          !
       ELSE
          !
          CALL errore( 'c_bands', 'isolve not implemented', 1 )
          !
       END IF
       !
       avg_iter = 0.D0
       !
       ! ... v_of_0 is (Vloc)(G=0)
       !
       v_of_0 = SUM( vltot(1:nrxx) ) / DBLE( nr1 * nr2 * nr3 )
       !
       CALL reduce( 1, v_of_0 )
       !
       if ( nks > 1 ) REWIND( iunigk )
       !
       ! ... For each k point diagonalizes the hamiltonian
       !
       k_loop: DO ik = 1, nks
          !
          ! ... read wave function for electric field
          !
          IF ( lelfield ) CALL davcio( evcel, nwordwfc, iunefield, ik, -1 )
          !
          current_k = ik
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          ! ... Reads the Hamiltonian and the list k+G <-> G of this k point
          !
          IF ( nks > 1 ) READ( iunigk ) npw, igk
          !
          ! ... do not recalculate k-points if restored from a previous run
          !
          IF ( ik <= ik_ ) THEN
             !
             CALL save_in_cbands( iter, ik, dr2 )
             !
             CYCLE k_loop
             !
          END IF
          !
          ! ... various initializations
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... read in wavefunctions from the previous iteration
          !
          IF ( noncolin ) THEN
             !
             IF ( nks > 1 .OR. .NOT. reduce_io.OR. lelfield ) &
                CALL davcio( evc_nc, nwordwfc, iunwfc, ik, -1 )
             !
          ELSE
             !
             IF ( nks > 1 .OR. .NOT. reduce_io .OR. lelfield) &
                CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
             !
          END IF
          !   
          ! ... Needed for LDA+U
          !
          IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunat, ik, -1 )
          !
          ! ... sets the kinetic energy
          !
          g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk(1:npw)) )**2 + &
                           ( xk(2,ik) + g(2,igk(1:npw)) )**2 + &
                           ( xk(3,ik) + g(3,igk(1:npw)) )**2 ) * tpiba2
          !
          IF ( qcutz > 0.D0 ) THEN
             !
             DO ig = 1, npw
                !
                g2kin(ig) = g2kin(ig) + qcutz * &
                            ( 1.D0 + erf( ( g2kin(ig) - ecfixed ) / q2sigma ) )
                !
             END DO
             !
          END IF
          !
          IF ( isolve == 1 ) THEN
             !
             ! ... Conjugate-Gradient diagonalization
             !
             ! ... h_diag is the precondition matrix
             !
             IF ( noncolin ) THEN
                !
                h_diag_nc = 1.D0
                !
                FORALL( ig = 1 : npwx )
                   !
                   h_diag_nc(ig,:) = 1.D0 + g2kin(ig) + &
                                     SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
                   !
                END FORALL
                !
             ELSE
                !
                h_diag = 1.D0
                !
                FORALL( ig = 1 : npw )
                   !
                   h_diag(ig) = 1.D0 + g2kin(ig) + &
                                SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
                   !
                END FORALL
                !
             END IF
             !
             ntry = 0
             !
             CG_loop : DO
                !
                IF ( iter > 1 .OR. istep > 0 .OR. ntry > 0 ) THEN
                   !
                   IF ( noncolin ) THEN
                      !
                      CALL cinitcgg( npwx, npw, nbnd, &
                                     nbnd, evc_nc, evc_nc, et(1,ik) )
                      !
                   ELSE
                      !
                      CALL cinitcgg( npwx, npw, nbnd, nbnd, evc, evc, et(1,ik) )
                      !
                   END IF
                   !
                   avg_iter = avg_iter + 1.D0
                   !
                END IF
                !
                IF ( noncolin ) THEN
                   !
                   CALL ccgdiagg( npwx, npw, nbnd, evc_nc, et(1,ik), &
                                  btype(1,ik), h_diag_nc, ethr, max_cg_iter, &
                                  .NOT. lscf, notconv, cg_iter )
                   !
                ELSE
                   !
                   CALL ccgdiagg( npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
                                  h_diag, ethr, max_cg_iter, .NOT. lscf, &
                                  notconv, cg_iter )
                   !
                END IF
                !
                avg_iter = avg_iter + cg_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
                IF ( noncolin ) THEN
                   !
                   IF ( nks > 1 .OR. .NOT. reduce_io .OR. lelfield ) &
                      CALL davcio( evc_nc, nwordwfc, iunwfc, ik, 1 )
                   !
                ELSE
                   !
                   IF ( nks > 1 .OR. .NOT. reduce_io .OR. lelfield ) &
                      CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
                   !
                END IF
                !
                ntry = ntry + 1                
                !
                ! ... exit condition
                !
                IF ( test_exit_cond() ) EXIT  CG_loop
                !
             END DO CG_loop
             !
          ELSE IF ( isolve == 2 ) THEN
             !
             ! ... RMM-DIIS method
             !
             IF ( noncolin ) THEN
                !
                DO ipol = 1, npol
                   !
                   h_diag_nc(1:npw,ipol) = g2kin(1:npw) + v_of_0
                   !
                END DO
                !
             ELSE
                !
                h_diag(1:npw) = g2kin(1:npw) + v_of_0
                !
             END IF
             !
             CALL usnldiag( h_diag_nc, s_diag_nc )
             !
             ntry = 0
             !
             RMMDIIS_loop: DO
                !
                IF ( noncolin ) THEN
                   !
                   CALL cdiisg_nc( npw, npwx, nbnd, diis_ndim, &
                                   evc_nc, et(1,ik), ethr, btype(1,ik), &
                                   notconv, diis_iter, iter, npol )
                   !
                ELSE
                   !
                   CALL cdiisg( npw, npwx, nbnd, evc, et(1,ik), &
                                btype(1,ik), notconv, diis_iter, iter )
                   !
                END IF
                !  
                avg_iter = avg_iter + diis_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
                IF ( noncolin ) THEN
                   !
                   IF ( nks > 1 .OR. .NOT. reduce_io ) &
                      CALL davcio( evc_nc, nwordwfc, iunwfc, ik, 1 )
                   !
                ELSE
                   !
                   IF ( nks > 1 .OR. .NOT. reduce_io ) &
                      CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
                   !
                END IF
                !
                ntry = ntry + 1                
                !
                ! ... exit condition
                !
                IF ( test_exit_cond() ) EXIT  RMMDIIS_loop
                !
             END DO RMMDIIS_loop
             !
          ELSE
             !
             ! ... Davidson diagonalization
             !
             ! ... h_diag are the diagonal matrix elements of the
             ! ... hamiltonian used in g_psi to evaluate the correction 
             ! ... to the trial eigenvectors
             !
             IF ( noncolin ) THEN
                !
                DO ipol = 1, npol
                   !
                   h_diag_nc(1:npw,ipol) = g2kin(1:npw) + v_of_0
                   !
                END DO
                !
                CALL usnldiag_nc( h_diag_nc, s_diag_nc )
                !
             ELSE
                !
                h_diag(1:npw) = g2kin(1:npw) + v_of_0
                !
                CALL usnldiag( h_diag, s_diag )
                !
             END IF
             !
             ntry = 0
             !
             david_loop: DO
                !
                lrot = ( iter == 1 )
                !
                IF ( noncolin ) THEN
                   !
                   CALL cegterg( npw, npwx, nbnd, nbndx, evc_nc, &
                                 ethr, okvan, et(1,ik), btype(1,ik), &
                                 notconv, lrot, dav_iter )
                   !
                ELSE
                   !
                   CALL cegterg( npw, npwx, nbnd, nbndx, evc, ethr, &
                                 okvan, et(1,ik), btype(1,ik), notconv, &
                                 lrot, dav_iter )
                   !
                END IF
                !
                avg_iter = avg_iter + dav_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
                IF ( noncolin ) THEN
                   !
                   IF ( nks > 1 .OR. .NOT. reduce_io .OR. lelfield ) &
                      CALL davcio( evc_nc, nwordwfc, iunwfc, ik, 1 )
                   !
                ELSE
                   !
                   IF ( nks > 1 .OR. .NOT. reduce_io .OR. lelfield) &
                      CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
                   !
                END IF
                !
                ntry = ntry + 1                
                !
                ! ... exit condition
                !
                IF ( test_exit_cond() ) EXIT david_loop                
                !
             END DO david_loop 
             !
          END IF
          !
          IF ( notconv > MAX( 5, nbnd / 4 ) ) THEN
             !
             CALL errore( 'c_bands', &
                        & 'too many bands are not converged', 1 )
             !
          END IF
          !
          ! ... save restart information
          !
          CALL save_in_cbands( iter, ik, dr2 )
          !
          IF (lbands) THEN
             IF (check_stop_now() )  call stop_run(.FALSE.)
          ENDIF
       END DO k_loop
       !
       ik_ = 0
       !
       CALL poolreduce( 1, avg_iter )
       !
       avg_iter = avg_iter / nkstot
       !
       WRITE( stdout, &
              '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
           ethr, avg_iter
       !
       ! ... deallocate work space
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
       RETURN
       !
     END SUBROUTINE c_bands_k
     !
     !-----------------------------------------------------------------------
     FUNCTION test_exit_cond()
       !-----------------------------------------------------------------------
       !
       ! ... this logical function is .TRUE. when iterative diagonalization
       ! ... is converged
       !
       IMPLICIT NONE
       !
       LOGICAL :: test_exit_cond
       !
       !
       test_exit_cond = .NOT. ( ( ntry <= 5 ) .AND. &
                                ( ( .NOT. lscf .AND. ( notconv > 0 ) ) .OR. &
                                  (       lscf .AND. ( notconv > 5 ) ) ) )
       !                          
     END FUNCTION test_exit_cond
     !     
END SUBROUTINE c_bands
