!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
  USE io_files,             ONLY : iunigk, nwordatwfc, iunsat, iunwfc, nwordwfc
  USE cell_base,            ONLY : tpiba2
  USE klist,                ONLY : nkstot, nks, xk, nelec
  USE uspp,                 ONLY : okvan
  USE cell_base,            ONLY : omega
  USE uspp,                 ONLY : vkb, nkb
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : g, gg, gstart, ngm, nl
  USE wvfct,                ONLY : g2kin, wg, nbndx, et, nbnd, npwx, igk, &
                                   ecutwfc, npw
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE control_flags,        ONLY : ethr, lscf, isolve
  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE scf,                  ONLY : vltot, vrs, v_of_0
  USE lsda_mod,             ONLY : nspin, current_spin, lsda, isk
  USE wavefunctions_module, ONLY : psic , evc
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE eff_v,                ONLY : rho_fft, rho_veff, evc_veff, nelecr
  USE mp_global,            ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variables
  !
  REAL(kind=DP)              :: veff (dfftp%nnr, nspin)        ! in: effective potential
  REAL(kind=DP), ALLOCATABLE :: vrs_ (:, :)        ! to keep the local potential
  REAL(kind=DP) ::    charge    ! out: the charge difference  between  rho_check & rho-fft
  !
  ! ... local variables
  !
  REAL(kind=DP) :: avg_iter
    ! average number of iterations
    ! the average of the potential
  REAL(kind=DP), ALLOCATABLE :: k_gamma(:)
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
  COMPLEX (kind=DP), ALLOCATABLE :: evc_(:,:)
    !  evc_   contains the  refined estimates of the eigenvectors
    !
    ! ... external functions
    !
  REAL(kind=DP), EXTERNAL :: qe_erf
    ! error function
    !
    !
  CALL start_clock( 'c_bands' )
  !
  ! ... allocate arrays
  !
  ALLOCATE( vrs_ ( dfftp%nnr, nspin ) )
  ALLOCATE( h_diag( npwx,1 ) )
  ALLOCATE( s_diag( npwx,1 ) )
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
       REAL(kind=DP) :: w1           ! weights
       !
!       WRITE( stdout, '(5X,"Davidson diagonalization (with overlap)")')
       !
       avg_iter = 0.D0
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
          ALLOCATE( k_gamma(3) )
          k_gamma = 0.d0
          CALL gk_sort (k_gamma , ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
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
             psic(1:dfftp%nnr) = sqrt(abs(rho_fft(1:dfftp%nnr,1)))
             !
             CALL fwfft ('Wave', psic, dffts)
             !
             evc_(1:npw,ibnd) = psic(nls(igk(1:npw)))
             !
          ENDDO
  !
          ! ... sets the kinetic energy
          !
          xk(1:3,ik) = k_gamma(1:3)
          !
          CALL g2_kin ( ik )
          !
          ! ... a band is considered empty when its occupation is less
          ! ... than 1.0 %
          !
          btype(:) = 0
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
             h_diag(1:npw,1) = g2kin(1:npw) + v_of_0
             !
!             CALL usnldiag( h_diag, s_diag )
             s_diag(:,1) = 1.d0
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
                CALL cegterg_vdw( npw, npwx, nbnd, nbndx, evc_, ethr, &
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
                IF ( test_exit_cond() ) exit david_loop
                !
             ENDDO david_loop
             !
          ENDIF
          !
          IF ( notconv > max( 5, nbnd / 4 ) ) THEN
             !
             CALL errore( 'c_bands', &
                        & 'too many bands are not converged', 1 )
             !
          ENDIF
          !
       ENDDO k_loop
       !
       CALL mp_sum( avg_iter, inter_pool_comm )
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
          CALL invfft ('Wave', psic, dffts)
          !
          ! compute the weight
          !
          nelecr = sum(rho_fft) * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
#ifdef __MPI
          CALL mp_sum( nelecr, intra_pool_comm )
#endif
          w1 =  nelecr  / omega
          !
          ! ... increment the charge density ...
          !
          DO ir = 1, dffts%nnr
             !
             rho_veff(ir,current_spin) = rho_veff(ir,current_spin) + &
                                          w1 * ( REAL( psic(ir) )**2 + &
                                                 aimag( psic(ir) )**2 )
             !
          ENDDO
          nelecr = sum(rho_veff) * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
#ifdef __MPI
          CALL mp_sum( nelecr, intra_pool_comm )
#endif
          !
       ENDDO
       !
       !  compute the charge difference
       !
       charge = 0.d0
       DO ir = 1, dfftp%nnr
          charge = charge + abs( rho_fft(ir,nspin) - rho_veff(ir,nspin) )
       ENDDO
       charge = charge * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3) / nelecr
#ifdef __MPI
          CALL mp_sum( charge, intra_pool_comm )
#endif
       !
       ! return the value for vrs and keep evc_ in evc_eff
       !
       vrs = vrs_
       evc_veff = evc_
       !
       DEALLOCATE( k_gamma )
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
       test_exit_cond = .not. ( ( ntry <= 5 ) .and. &
                                ( ( .not. lscf .and. ( notconv > 0 ) ) .or. &
                                  (       lscf .and. ( notconv > 5 ) ) ) )
       !
     END FUNCTION test_exit_cond
     !
END SUBROUTINE check_v_eff
