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
SUBROUTINE check_v_eff ( veff, charge )
  !----------------------------------------------------------------------------
  !
  ! ... this is a wrapper to specific calls
  !
  ! ... internal procedures :
  !
  ! ... diag_v_eff()      : for diagonalizing effective potential
  ! ... test_exit_cond()  : the test on the iterative diagonalization
  !
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps4
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : gamma_only
  USE io_files,             ONLY : iunigk, nwordatwfc, iunat, iunwfc, nwordwfc
  USE cell_base,            ONLY : tpiba2 
  USE klist,                ONLY : nkstot, nks, xk, nelec
  USE uspp,                 ONLY : okvan
  USE cell_base,            ONLY : omega
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g, gg, gstart, ecfixed, qcutz, q2sigma, nrxx, &
                                   nr1, nr2, nr3, nrx1, nrx2, nrx3,ngm, ecutwfc, nl  
  USE wvfct,                ONLY : g2kin, wg, nbndx, et, nbnd, npwx, igk, &
                                   npw
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid  
  USE control_flags,        ONLY : diis_ndim, ethr, lscf, max_cg_iter, &
                                   isolve, reduce_io
  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE scf,                  ONLY : rho, vltot, vrs
  USE lsda_mod,             ONLY : nspin, current_spin, lsda, isk
  USE wavefunctions_module, ONLY : psic , evc 
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE eff_v,                ONLY : rho_fft, rho_veff, evc_veff, nelecr
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variables
  !
  REAL(KIND=DP)              :: veff (nrxx, nspin)        ! in: effective potential
  REAL(KIND=DP), ALLOCATABLE :: vrs_ (:, :)        ! to keep the local potential
  REAL(KIND=DP) ::    charge    ! out: the charge difference  between  rho_check & rho-fft
  !
  ! ... local variables
  !
  REAL(KIND=DP) :: avg_iter, v_of_0
    ! average number of iterations
    ! the average of the potential
  REAL(KIND=DP), ALLOCATABLE :: k_gamma(:)
    ! gamma point
  INTEGER :: ik, ig, ibnd, dav_iter, ntry, notconv
    ! counter on k points
    ! counter on G vectors
    ! counter on bands
    ! number of iterations in Davidson
    ! number or repeated call to diagonalization in case of non convergence
    ! number of notconverged elements
  INTEGER, ALLOCATABLE :: btype(:)
    ! type of band: conduction (1) or valence (0)  
  COMPLEX (KIND=DP), ALLOCATABLE :: evc_(:,:)
    !  evc_   contains the  refined estimates of the eigenvectors
    !
    ! ... external functions
    !
  REAL(KIND=DP), EXTERNAL :: erf
    ! error function  
    !
    !
  CALL start_clock( 'c_bands' )
  !
  ! ... allocate arrays
  !
  ALLOCATE( vrs_ ( nrxx, nspin ) )
  ALLOCATE( h_diag( npwx ) )    
  ALLOCATE( s_diag( npwx ) )   
  ALLOCATE( btype(  nbnd ) )       
  ALLOCATE( evc_(npwx,nbnd ) )
  !
  CALL diag_v_eff()
  !
  ! ... deallocate arrays
  !
  DEALLOCATE( s_diag )
  DEALLOCATE( h_diag )
  DEALLOCATE( btype )
  DEALLOCATE( evc_ )
  DEALLOCATE( vrs_  )
  !       
  CALL stop_clock( 'c_bands' )  
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedure
     !-----------------------------------------------------------------------
     SUBROUTINE diag_v_eff()
       !-----------------------------------------------------------------------
       !
       ! ... This routine is a driver for the diagonalization routines of the
       ! ... total Hamiltonian at g-point using Davidson algorithm.
       !
       IMPLICIT NONE
       !
       ! ... here the local variables
       !
       INTEGER :: ir
       !
       REAL(KIND=DP) :: w1           ! weights
       !       
!       WRITE( stdout, '(5X,"Davidson diagonalization (with overlap)")')
       !
       avg_iter = 0.D0
       !
       ! ... v_of_0 is (Vloc)(G=0)
       !
       v_of_0 = SUM( vltot(1:nrxx) ) / REAL( nr1 * nr2 * nr3 )
       !
       CALL reduce( 1, v_of_0 )
       !
       nks = 1  ! for TF+vW
       !
       ! ... For each k point diagonalizes the hamiltonian
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          ! ... generates the Hamiltonian and the 
          !     list k+G <-> G of this k point
          !
          allocate( k_gamma(3) )
          k_gamma = 0.d0
          call gk_sort (k_gamma , ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
          !
          ! ... various initializations
          !
          nkb = 0   ! for TF+vW
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... read in wavefunctions from the previous iteration
          !
!          IF ( nks > 1 .OR. .NOT. reduce_io ) &
!             CALL davcio( evc_, nwordwfc, iunwfc, ik, -1 )
          !
          ! trial wave function for V_eff
          !
          DO ibnd = 1, nbnd
             !
             psic(1:nrxx) = sqrt(abs(rho_fft(1:nrxx,1)))
             !
             CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
             !
             evc_(1:npw,ibnd) = psic(nls(igk(1:npw)))
             !
          ENDDO
  !
          ! ... sets the kinetic energy
          !
          xk(1:3,ik) = k_gamma(1:3)
          g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk(1:npw)) )**2 + &
                           ( xk(2,ik) + g(2,igk(1:npw)) )**2 + &
                           ( xk(3,ik) + g(3,igk(1:npw)) )**2 ) * tpiba2
          !
          !
          IF ( qcutz > 0.D0 ) THEN
             DO ig = 1, npw
                g2kin (ig) = g2kin(ig) + qcutz * &
                             ( 1.D0 + erf( ( g2kin(ig) - ecfixed ) / q2sigma ) )
             END DO
          END IF
          !
          btype(:) = 0
          !
          ! ... a band is considered empty when its occupation is less 
          ! ... than 1.0 %
          !   
          WHERE( wg(:,ik) < 0.01D0 ) btype(:) = 0
          !
          IF ( isolve == 0 ) THEN
             !
             ! ... Davidson diagonalization
             !
             ! ... h_diag are the diagonal matrix elements of the
             ! ... hamiltonian used in g_psi to evaluate the correction 
             ! ... to the trial eigenvectors
             !
             h_diag(1:npw) = g2kin(1:npw) + v_of_0
             !
!             CALL usnldiag( h_diag, s_diag )
             s_diag(:) = 1.d0
             !
             ntry = 0
             !
             !
             ! set input value for TF+vW
             !
             ethr = 1.D-12
             okvan = .false.
             btype(:) = 0
             !
             david_loop: DO
                !
                ! pass the effective potential for TF+vW
                !
                vrs_ = vrs
                vrs = veff
                !
                CALL cegterg( npw, npwx, nbnd, nbndx, evc_, ethr, &
                              okvan, et(1,ik), btype, notconv, dav_iter )
                !
                avg_iter = avg_iter + dav_iter
                !
                ! ... save wave-functions to be used as input for the
                ! ... iterative diagonalization of the next scf iteration 
                ! ... and for rho calculation
                !
!                IF ( nks > 1 .OR. .NOT. reduce_io ) &
!                   CALL davcio( evc_, nwordwfc, iunwfc, ik, 1 )
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
       END DO k_loop
       !
       CALL poolreduce( 1, avg_iter )
       !
       avg_iter = avg_iter / nkstot
       !
!       WRITE( stdout, &
!              '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
!           ethr, avg_iter
       !
       ! compute the charge density from v_eff 
       !
       rho_veff =   0.D0
       !
       ik = 1  ! for TF+vW
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       DO ibnd = 1, nbnd
          !
          psic(:) = ( 0.D0, 0.D0 )
          !
          psic(nls(igk(1:npw))) = evc_(1:npw,ibnd)
          !
          CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
          !
          ! compute the weight
          ! 
          nelecr = sum(rho_fft) * omega / (nr1*nr2*nr3)
#ifdef __PARA
          call reduce(1,nelecr)
#endif
          w1 =  nelecr  / omega 
          !
          ! ... increment the charge density ...
          !
          DO ir = 1, nrxxs
             !
             rho_veff(ir,current_spin) = rho_veff(ir,current_spin) + &
                                          w1 * ( REAL( psic(ir) )**2 + &
                                                 AIMAG( psic(ir) )**2 )
             !
          END DO
          nelecr = sum(rho_veff) * omega / (nr1*nr2*nr3)
#ifdef __PARA
          call reduce(1,nelecr)
#endif
          !
       END DO
       !
       !  compute the charge difference
       !
       charge = 0.d0
       DO ir = 1, nrxx
          charge = charge + abs( rho_fft(ir,nspin) - rho_veff(ir,nspin) )
       END DO
       charge = charge * omega / (nr1*nr2*nr3) / nelecr
#ifdef __PARA
          call reduce(1,charge)
#endif
       !
       ! return the value for vrs and keep evc_ in evc_eff
       !
       vrs = vrs_
       evc_veff = evc_
       !
       deallocate( k_gamma )
       !
       RETURN
       !
     END SUBROUTINE diag_v_eff
     !
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
END SUBROUTINE check_v_eff
