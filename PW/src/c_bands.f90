!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE c_bands( iter )
  !----------------------------------------------------------------------------
  !! Driver routine for the Hamiltonian diagonalization ones.
  !! It reads the Hamiltonian and an initial guess of the wavefunctions
  !! from a file and computes initialization quantities for the
  !! diagonalization routines.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE klist,                ONLY : nkstot, nks, ngk, igk_k, xk
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, isolve, restart, use_gpu, iverbosity
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_projectors, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions,        ONLY : evc
  USE bp,                   ONLY : lelfield
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE check_stop,           ONLY : check_stop_now
  USE gcscf_module,         ONLY : lgcscf
  USE add_dmft_occ,         ONLY : dmft, dmft_updated
  USE uspp_init,            ONLY : init_us_2
  USE device_fbuff_m,       ONLY : dev_buf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  !! iteration index
  !
  ! ... local variablems
  !
  REAL(DP) :: avg_iter
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios
  ! ik : counter on k points
  ! ik_: k-point already done in a previous run
  LOGICAL :: exst
  LOGICAL,EXTERNAL :: rmm_use_davidson, rmm_use_paro
  !
  INTEGER :: ierr
  !
  !
  CALL start_clock( 'c_bands' ); !write (*,*) 'start c_bands' ; FLUSH(6)
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands( ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  ! ... (not needed for a single k-point: this is done in wfcinit, 
  ! ...  directly from file, in order to avoid wasting memory)
  !
  DO ik = 1, ik_
     IF ( nks > 1 .OR. lelfield ) THEN
        CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
        !$acc update device(evc)
     END IF
  ENDDO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSEIF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")')
  ELSEIF ( isolve == 3 ) THEN
     WRITE( stdout, '(5X,"ParO style diagonalization")')
  ELSEIF ( isolve == 4 ) THEN
     IF (rmm_use_davidson(iter)) THEN 
       WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
     ELSE IF (rmm_use_paro(iter)) THEN 
      WRITE( stdout, '(5X,"ParO style diagonalization")')
     ELSE 
       WRITE( stdout, '(5X,"RMM-DIIS diagonalization")')
     END IF 
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
  ENDIF
  !
  if (iverbosity > 0) CALL print_mem_usage(stdout, 'c_bands before calling an iterative solver')
  !
  ! ... For each k point diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) CALL phase_factor(ik)
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     CALL g2_kin( ik )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb, .true. )
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF ( nks > 1 .OR. lelfield ) THEN
          CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !$acc update device(evc)
     END IF
     !
     ! ... Needed for DFT+Hubbard
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (Hubbard_projectors.NE.'pseudo') ) THEN
        CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
        !$acc update device(wfcU)
     END IF
     !
     ! ... diagonalization of bands for k-point ik
     ! ... (skip only in charge self-consistent DFT+DMFT calculations)
     !
     IF (.NOT. ( dmft .AND. .NOT. dmft_updated ) ) THEN
        call diag_bands ( iter, ik, avg_iter )
        !sync evc here to allow later use of converged wavefunction on host
        !$acc update self(evc)
     END IF
     !
     ! ... save wave-functions to be used as input for the
     ! ... iterative diagonalization of the next scf iteration
     ! ... and for rho calculation
     !
     IF ( nks > 1 .OR. lelfield ) CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     !
     IF (ik <= nkdum) THEN
        IF (check_stop_now()) THEN
           CALL save_in_cbands( ik, ethr, avg_iter, et )
           RETURN
        ENDIF
     ENDIF
     !
     CALL dev_buf%reinit( ierr )
     IF ( ierr .ne. 0 ) CALL infomsg( 'c_bands', 'Cannot reset GPU buffers! Some buffers still locked.' )
     !
  ENDDO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, &
       '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
       ethr, avg_iter
  !
  CALL stop_clock( 'c_bands' ); !write (*,*) 'stop c_bands' ; FLUSH(6)
  !
  RETURN
  !
END SUBROUTINE c_bands
!
!----------------------------------------------------------------------------
SUBROUTINE diag_bands( iter, ik, avg_iter )
  !----------------------------------------------------------------------------
  !! Driver routine for diagonalization at each k-point. Types of iterative
  !! diagonalizations currently in use:
  !
  !! * Davidson algorithm (all-band);
  !! * Conjugate Gradient (band-by-band);
  !! * Projected Preconditioned Conjugate Gradient (block);
  !! * Parallel Orbital update (all-band);
  !! * RMM-DIIS algorithm (all-band).
  !
  !! Internal procedures:
  !
  !! * \(\textrm{diag_bands_gamma}\)(): optimized algorithms for gamma sampling
  !!                                    of the BZ (real Hamiltonian);
  !! * \(\textrm{diag_bands_k}\)(): general algorithm for arbitrary BZ sampling
  !!                                (complex Hamiltonian);
  !! * \(\textrm{test_exit_cond}\)(): the test on the iterative diagonalization.
  !
  USE kinds,                ONLY : DP
  USE buffers,              ONLY : get_buffer
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfc, iunefieldp, iunefieldm
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE gvect,                ONLY : gstart
  USE wvfct,                ONLY : g2kin, nbndx, et, nbnd, npwx, btype
  USE control_flags,        ONLY : ethr, lscf, max_cg_iter, isolve, &
                                   rmm_ndim, rmm_conv, gs_nblock, &
                                   gamma_only, use_para_diag, use_gpu
  USE ldaU,                 ONLY : hub_pot_fix
  USE noncollin_module,     ONLY : npol
  USE wavefunctions,        ONLY : evc
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE bp,                   ONLY : lelfield, evcel, evcelp, evcelm, bec_evcel, &
                                   gdir, l3dstring, efield, efield_cry
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type, &
                                   allocate_bec_type_acc, deallocate_bec_type_acc
  USE klist,                ONLY : nks, ngk
  USE gcscf_module,         ONLY : lgcscf
  USE mp_bands,             ONLY : nproc_bgrp, intra_bgrp_comm, inter_bgrp_comm, &
                                   my_bgrp_id, nbgrp
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE xc_lib,               ONLY : exx_is_active
  USE gcscf_module,         ONLY : lgcscf
  !
  USE control_flags,        ONLY : scissor
  USE sci_mod,              ONLY : evcc
#if defined (__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_h_diag
#endif
  !
  IMPLICIT NONE
  !
  ! please do not capitalize (FORD rules)
  !  
  INTEGER, INTENT(IN) :: iter
  !! iteration index
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  REAL(KIND=DP), INTENT(INOUT) :: avg_iter
  !! average number of H*psi products
  !
  ! ... local variables
  !
  REAL(KIND=DP) :: cg_iter, rmm_iter
  ! (weighted) number of iterations in RMM-DIIS
  INTEGER :: npw, dav_iter, ntry, notconv, nhpsi
  ! number of iterations in Davidson
  ! number or repeated call to diagonalization in case of non convergence
  ! number of notconverged elements
  INTEGER :: ierr, ipw, ibnd, ibnd_start, ibnd_end
  !
  LOGICAL :: lrot
  ! .TRUE. if the wfc have already be rotated
  !
  COMPLEX (DP), POINTER :: hevc(:,:), sevc(:,:)
  !
  ! Davidson and RMM-DIIS diagonalization uses these external routines on groups of nvec bands
  EXTERNAL h_psi, s_psi, g_psi
  EXTERNAL h_psi_gpu, s_psi_acc
  ! subroutine h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
  ! subroutine s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
  ! subroutine g_psi(npwx,npw,nvec,psi,eig)   computes G*psi -> psi
  !------------------------------------------------------------------------
  ! CG diagonalization uses these external routines on a single band
  EXTERNAL hs_1psi, s_1psi, hs_psi
  EXTERNAL hs_psi_gpu
  EXTERNAL hs_1psi_gpu, s_1psi_gpu
  LOGICAL, EXTERNAL   :: rmm_use_davidson, rmm_use_paro
  ! subroutine hs_1psi(npwx,npw,psi,hpsi,spsi)  computes H*psi and S*psi
  ! subroutine s_1psi(npwx,npw,psi,spsi)        computes S*psi (if needed)
  ! In addition to the above the initial wfc rotation uses h_psi, and s_psi
  !------------------------------------------------------------------------
  ! subroutine h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
  ! subroutine s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
  !------------------------------------------------------------------------
  ! ParO diagonalization uses these external routines on a single band
  ! subroutine hs_1psi(npwx,npw,psi,hpsi,spsi)  computes H*psi and S*psi
  ! In addition to the above the initial wfc rotation uses h_psi, and s_psi
  !
  ALLOCATE( h_diag( npwx, npol ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' diag_bands ', ' cannot allocate h_diag ', ABS(ierr) )
  !
  ALLOCATE( s_diag( npwx, npol ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' diag_bands ', ' cannot allocate s_diag ', ABS(ierr) )
  !$acc enter data create (h_diag, s_diag)
  !
  ipw=npwx
  CALL mp_sum(ipw, intra_bgrp_comm)
  IF ( nbndx > ipw ) &
     CALL errore ( 'diag_bands', 'too many bands, or too few plane waves',1)
  !
  ! ... allocate space for <beta_i|psi_j> - used in h_psi and s_psi
  !
  CALL allocate_bec_type_acc( nkb, nbnd, becp, intra_bgrp_comm )
  !
  npw = ngk(ik)
  IF ( gamma_only ) THEN
     !
     CALL diag_bands_gamma()
     !
  ELSE
     !
     CALL diag_bands_k()
     !
  ENDIF
  !
  ! ... deallocate work space
  !
  CALL deallocate_bec_type_acc( becp )
  !$acc exit data delete (h_diag, s_diag)
  DEALLOCATE( s_diag )
  DEALLOCATE( h_diag )
  !
  IF ( notconv > MAX( 5, nbnd / 4 ) ) THEN
     !
     IF ( hub_pot_fix ) THEN
        ! If perturbing a Hubbard manifold using Hubbard_alpha, 
        ! need to diagonalize with very tight convergence threshholds
        ! even during the first iteration.
        ! Thus, c_bands should not throw an error even though many
        ! bands might not achieve this convergence criterion.
        !
        WRITE( stdout, '(5X,"c_bands: ",I2, " eigenvalues not converged")' ) notconv
        WRITE( stdout, '(5X,"WARNING: c_bands: not aborting due to active Hubbard alpha")' ) notconv
     ELSE
        CALL errore( 'c_bands', 'too many bands are not converged', 1 )
     ENDIF
     !
  ELSEIF ( notconv > 0 ) THEN
     !
     WRITE( stdout, '(5X,"c_bands: ",I2, " eigenvalues not converged")' ) notconv
     !
  ENDIF
  !
  RETURN
  !
 CONTAINS
  !
  ! ... internal procedures
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_bands_gamma()
    !-----------------------------------------------------------------------
    !
    ! ... Diagonalization of a real Hamiltonian
    !
    IMPLICIT NONE
    !
    INTEGER :: ig
    !
    IF ( isolve == 1 .OR. isolve == 2 .OR. isolve == 3 .OR. rmm_use_paro(iter))   THEN
       !
       ! ... (Projected Preconditioned) Conjugate-Gradient diagonalization
       !
       ! ... h_diag is the precondition matrix
       !
       IF ( isolve == 1 .OR. isolve == 2 ) THEN
          !$acc parallel loop present(g2kin)
          DO ig = 1, npw 
             h_diag(ig,1) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
          END DO
       ELSE
          CALL usnldiag( npw, npol, h_diag, s_diag )
       END IF
       !
       ntry = 0
       !
       CG_loop : DO
          !
          IF ( isolve == 1 .OR. isolve == 2 ) THEN
             lrot = ( iter == 1 .AND. ntry == 0 )
             !
             IF ( .NOT. lrot ) THEN
                !
                IF (.not. use_gpu) THEN
                   CALL rotate_wfc( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
                ELSE
                   CALL rotate_wfc_gpu( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
                END IF
                !
                avg_iter = avg_iter + 1.D0
                !
             ENDIF
          ENDIF
          !
          IF ( isolve == 1 ) THEN
             IF (.not. use_gpu) THEN
                CALL rcgdiagg( hs_1psi, s_1psi, h_diag, &
                         npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
                         ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
             ELSE
                CALL rcgdiagg( hs_1psi_gpu, s_1psi_gpu, h_diag, &
                         npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
                         ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
                !
             END IF
             !
             avg_iter = avg_iter + cg_iter
             !
          ELSE
             !
             IF (.not. use_gpu ) THEN
               CALL paro_gamma_new( h_psi, s_psi, hs_psi, g_psi, okvan, &
                          npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
               !
               avg_iter = avg_iter + nhpsi/float(nbnd) 
               ! write (6,*) ntry, avg_iter, nhpsi
               !
             ELSE
               !$acc host_data use_device(et)
               CALL paro_gamma_new( h_psi_gpu, s_psi_acc, hs_psi_gpu, g_psi, okvan, &
                          npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
               !$acc end host_data
               !
               avg_iter = avg_iter + nhpsi/float(nbnd) 
               ! write (6,*) ntry, avg_iter, nhpsi
               !
             ENDIF  
          ENDIF
          !
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  CG_loop
          !
       ENDDO CG_loop
       !$acc update self(et)
       !
    ELSE IF ( isolve == 4 .AND. .NOT. rmm_use_davidson(iter)) THEN
       !
       ! ... RMM-DIIS diagonalization
       !
       ALLOCATE( hevc  ( npwx*npol, nbnd ) )
       !$acc enter data create(hevc)
       IF ( okvan ) THEN
          ALLOCATE( sevc( npwx*npol, nbnd ) )
          !$acc enter data create(sevc)
       ELSE
          sevc => evc
       END IF
       !
       ntry = 0
       !
       RMM_loop : DO
          !
          lrot = ( iter == 1 .AND. ntry == 0 )
          !
!edp
!          IF ( .NOT. lrot ) THEN
          IF (lrot .AND. .NOT. lscf ) THEN
              !!
              !$acc parallel loop present(g2kin)
              DO ig = 1, npw
                 h_diag(ig,1) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
              END DO
              !
              IF (.not. use_gpu ) THEN
                CALL paro_gamma_new( h_psi, s_psi, hs_psi, g_psi, okvan, &
                           npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
                !
                avg_iter = avg_iter + nhpsi/float(nbnd) 
                ! write (6,*) ntry, avg_iter, nhpsi
                !
              ELSE
                !$acc host_data use_device(et)
                CALL paro_gamma_new( h_psi_gpu, s_psi_acc, hs_psi_gpu, g_psi, okvan, &
                           npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
                !$acc end host_data
                !$acc update self(et)
                !
                avg_iter = avg_iter + nhpsi/float(nbnd) 
                ! write (6,*) ntry, avg_iter, nhpsi
                !
              ENDIF  
               !
          ELSE IF ( .NOT. lrot ) THEN
             !
             IF (.not. use_gpu) THEN
                CALL rotate_xpsi_driver( h_psi, s_psi, h_psi, s_psi, npwx, npw, nbnd, nbnd, evc, npol, okvan, &
                               evc, hevc, sevc, et(:,ik), use_para_diag, .TRUE. )
#if defined(__CUDA)
             ELSE
                CALL rotate_xpsi_driver( h_psi, s_psi, h_psi_gpu, s_psi_acc, npwx, npw, nbnd, nbnd, evc, npol, okvan, &
                               evc, hevc, sevc, et(:,ik), use_para_diag, .TRUE.)
#endif
             END IF
             !
             avg_iter = avg_iter + 1.D0
             !
          END IF
          !
          !
          IF (.not. use_gpu) THEN
            CALL rrmmdiagg( h_psi, s_psi, npwx, npw, nbnd, evc, hevc, sevc, &
                         et(1,ik), g2kin(1), btype(1,ik), ethr, rmm_ndim, &
                         okvan, lrot, exx_is_active(), notconv, rmm_iter )
          ELSE
             CALL rrmmdiagg( h_psi_gpu, s_psi_acc, npwx, npw, nbnd, evc, hevc, sevc, &
                          et(1,ik), g2kin, btype(1,ik), ethr, rmm_ndim, &
                          okvan, lrot, exx_is_active(), notconv, rmm_iter )
          END IF
          !
          !
          IF ( lscf .AND. ( .NOT. rmm_conv ) ) notconv = 0
          !
          avg_iter = avg_iter + rmm_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  RMM_loop
          !
       END DO RMM_loop
       !
       ! ... Gram-Schmidt orthogonalization
       !
       CALL gram_schmidt_gamma( npwx, npw, nbnd, evc, hevc, sevc, et(1,ik), &
                       okvan, .TRUE., .TRUE., gs_nblock )
       !
       avg_iter = avg_iter + 0.5D0
       !
       !$acc exit data delete(hevc)
       DEALLOCATE( hevc )
       IF ( okvan ) THEN
          !$acc exit data delete(sevc)
          DEALLOCATE( sevc )
       ELSE
          NULLIFY( sevc )
       END IF
       !
    ELSE
       !
       ! ... Davidson diagonalization
       !
       ! ... h_diag are the diagonal matrix elements of the
       ! ... hamiltonian used in g_psi to evaluate the correction
       ! ... to the trial eigenvectors
       !
       CALL usnldiag( npw, npol, h_diag, s_diag )
#if defined (__OSCDFT)
       IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_h_diag(oscdft_ctx, h_diag)
#endif
       !
       ntry = 0
       !
       david_loop: DO
          !
          lrot = ( iter == 1 )
          !
          IF (.not. use_gpu) THEN
             IF ( use_para_diag ) THEN
!                ! make sure that all processors have the same wfc
                CALL pregterg( h_psi, s_psi, okvan, g_psi, &
                            npw, npwx, nbnd, nbndx, evc, ethr, &
                            et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
             ELSE
                CALL regterg (  h_psi, s_psi, okvan, g_psi, &
                         npw, npwx, nbnd, nbndx, evc, ethr, &
                         et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
             END IF
             ! 
          ELSE
             IF ( use_para_diag ) THEN
                !$acc host_data use_device(et)
                CALL pregterg_gpu( h_psi_gpu, s_psi_acc, okvan, g_psi, &
                            npw, npwx, nbnd, nbndx, evc, ethr, &
                            et(1, ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                !$acc end host_data
                !
             ELSE
                !
                !$acc host_data use_device(et)
                CALL regterg (  h_psi_gpu, s_psi_acc, okvan, g_psi, &
                         npw, npwx, nbnd, nbndx, evc, ethr, &
                         et(1, ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                !$acc end host_data
             END IF
             !$acc update self(et)
          END IF
          !
          avg_iter = avg_iter + dav_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  david_loop
          !
       ENDDO david_loop
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE diag_bands_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_bands_k()
    !-----------------------------------------------------------------------
    !! Complex Hamiltonian diagonalization.
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ipol
    REAL(DP) :: eps=0.000001d0
    ! --- Define a small number ---
    INTEGER :: ig
    !
    !write (*,*) ' enter diag_bands_k'; FLUSH(6)
    IF ( lelfield ) THEN
       !
       ! ... save wave functions from previous iteration for electric field
       !
       evcel = evc
       !
       !... read projectors from disk
       !
       IF (.NOT.l3dstring .AND. ABS(efield)>eps ) THEN
          CALL get_buffer (evcelm(:,:,gdir), nwordwfc, iunefieldm, ik+(gdir-1)*nks)
          CALL get_buffer (evcelp(:,:,gdir), nwordwfc, iunefieldp, ik+(gdir-1)*nks)
       ELSE
          DO ipol = 1, 3
             IF ( ABS(efield_cry(ipol))>eps ) THEN
                CALL get_buffer( evcelm(:,:,ipol), nwordwfc, iunefieldm, ik+(ipol-1)*nks )
                CALL get_buffer( evcelp(:,:,ipol), nwordwfc, iunefieldp, ik+(ipol-1)*nks )
             ENDIF
          ENDDO
       ENDIF
       !
       IF ( okvan ) THEN
          !
          CALL allocate_bec_type( nkb, nbnd, bec_evcel )
          !
          !$acc update self(vkb)
          CALL calbec( npw, vkb, evcel, bec_evcel )
          !
       ENDIF
       !
    ENDIF
    !
    IF (scissor) evcc = evc
    !
    !write (*,*) ' current isolve value ( 0 Davidson, 1 CG, 3 PARO, 4 RMM)', isolve; FLUSH(6)
    IF ( isolve == 1 .OR. isolve == 2 .OR. isolve == 3 .or. rmm_use_paro(iter)) THEN
       !
       ! ... (Projected Preconditioned) Conjugate-Gradient diagonalization
       !
       ! ... h_diag is the precondition matrix
       !
       !write (*,*) ' inside CG solver branch '
       !
       IF ( isolve == 1 .OR. isolve == 2) THEN
          !$acc parallel loop present(g2kin)
          DO ig = 1, npwx
             h_diag(ig,:) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
          END DO
       ELSE
          CALL usnldiag( npw, npol, h_diag, s_diag )
       ENDIF
       !
       ntry = 0
       !
       CG_loop : DO
          !
          IF ( isolve == 1 .OR. isolve == 2 ) THEN
             lrot = ( iter == 1 .AND. ntry == 0 )
             !
             IF ( .NOT. lrot ) THEN
                !
                IF ( .not. use_gpu ) THEN
                   CALL rotate_wfc( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
                ELSE
                   CALL rotate_wfc_gpu( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
                END IF
                !
                avg_iter = avg_iter + 1.D0
             ENDIF
          ENDIF
          !
          IF ( isolve == 1) then
             IF ( .not. use_gpu ) THEN
                CALL ccgdiagg( hs_1psi, s_1psi, h_diag, &
                         npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), &
                         ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
             ELSE
                CALL ccgdiagg( hs_1psi_gpu, s_1psi_gpu, h_diag, &
                         npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), &
                         ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
             END IF
             !
             avg_iter = avg_iter + cg_iter
             !
          ELSE 
             !
             IF ( .not. use_gpu ) THEN
               CALL paro_k_new( h_psi, s_psi, hs_psi, g_psi, okvan, &
                        npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
               !
               avg_iter = avg_iter + nhpsi/float(nbnd) 
               ! write (6,*) ntry, avg_iter, nhpsi
             ELSE
               !$acc host_data use_device(et)
               CALL paro_k_new( h_psi_gpu, s_psi_acc, hs_psi_gpu, g_psi, okvan, &
                        npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
               !$acc end host_data
               !
               avg_iter = avg_iter + nhpsi/float(nbnd) 
               ! write (6,*) ntry, avg_iter, nhpsi
               !
             END IF
          ENDIF
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  CG_loop
          !
       ENDDO CG_loop
       !$acc update self(et)
       !
    ELSE IF ( isolve == 4 .AND. .NOT. rmm_use_davidson(iter) )  THEN
       !
       ! ... RMM-DIIS diagonalization
       !
       ALLOCATE( hevc( npwx*npol, nbnd ) )
       !$acc enter data create(hevc)
       IF ( okvan ) THEN
          ALLOCATE( sevc( npwx*npol, nbnd ) )
          !$acc enter data create(sevc)
       ELSE
          sevc => evc
       END IF
       !
       ntry = 0
       !
       !
       RMM_loop : DO
          !
          lrot = ( iter == 1 .AND. ntry == 0 )
          !
!edp
!          IF ( .NOT. lrot ) THEN
          IF (lrot .AND. .NOT. lscf ) THEN
              CALL usnldiag(npw, npol, h_diag, s_diag )
              !
              IF ( .not. use_gpu ) THEN
                CALL paro_k_new( h_psi, s_psi, hs_psi, g_psi, okvan, &
                         npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
                !
                avg_iter = avg_iter + nhpsi/float(nbnd) 
                ! write (6,*) ntry, avg_iter, nhpsi
              ELSE
                !$acc host_data use_device(et)
                CALL paro_k_new( h_psi_gpu, s_psi_acc, hs_psi_gpu, g_psi, okvan, &
                         npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
                !$acc end host_data
                !$acc update self(et)
                !
                avg_iter = avg_iter + nhpsi/float(nbnd) 
                ! write (6,*) ntry, avg_iter, nhpsi
                !
              END IF
              !
          ELSE IF ( .NOT. lrot ) THEN
             !
             IF ( .not. use_gpu ) THEN
                CALL rotate_xpsi_driver( h_psi, s_psi, h_psi, s_psi, npwx, npw, nbnd, nbnd, evc, npol, okvan, &
                                  evc, hevc, sevc, et(:,ik), & 
                                  use_para_diag, gamma_only )
#if defined(__CUDA)
             ELSE
                CALL rotate_xpsi_driver( h_psi, s_psi, h_psi_gpu, s_psi_acc, npwx, npw, nbnd, nbnd, evc, npol, okvan, &
                                  evc, hevc, sevc, et(:,ik), &
                                  use_para_diag, gamma_only )
#endif
             END IF
             !
             avg_iter = avg_iter + 1.D0
             !
          END IF
          !
          IF ( .not. use_gpu ) THEN
             CALL crmmdiagg( h_psi, s_psi, npwx, npw, nbnd, npol, evc, hevc, sevc, &
                             et(1,ik), g2kin(1), btype(1,ik), ethr, rmm_ndim, &
                             okvan, lrot, exx_is_active(), notconv, rmm_iter )
          ELSE
             CALL crmmdiagg( h_psi_gpu, s_psi_acc, npwx, npw, nbnd, npol, evc, hevc, sevc, &
                             et(1,ik), g2kin(1), btype(1,ik), ethr, rmm_ndim, &
                             okvan, lrot, exx_is_active(), notconv, rmm_iter )
          END IF
          !
          IF ( lscf .AND. ( .NOT. rmm_conv ) ) notconv = 0
          !
          avg_iter = avg_iter + rmm_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  RMM_loop
          !
       END DO RMM_loop
       !
       ! ... Gram-Schmidt orthogonalization
       !
       CALL gram_schmidt_k( npwx, npw, nbnd, npol, evc, hevc, sevc, et(1,ik), &
                          okvan, .TRUE., .TRUE., gs_nblock )
       !
       avg_iter = avg_iter + 0.5D0
       !
       !$acc exit data delete(hevc)
       DEALLOCATE( hevc )
       IF ( okvan ) THEN
          !$acc exit data delete(sevc)
          DEALLOCATE( sevc )
       ELSE
          NULLIFY( sevc )
       END IF
       !
    ELSE
       !
       ! ... Davidson diagonalization
       !
       ! ... h_diag are the diagonal matrix elements of the
       ! ... hamiltonian used in g_psi to evaluate the correction
       ! ... to the trial eigenvectors
       !
       CALL usnldiag( npw, npol, h_diag, s_diag )
#if defined (__OSCDFT)
       IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_h_diag(oscdft_ctx, h_diag)
#endif
       !
       ntry = 0
       !
       david_loop: DO
          !
          lrot = ( iter == 1 )
          !
          IF (.not. use_gpu ) THEN
             IF ( use_para_diag ) then
                !
                CALL pcegterg( h_psi, s_psi, okvan, g_psi, &
                               npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                               et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                !
             ELSE
                !
                CALL cegterg ( h_psi, s_psi, okvan, g_psi, &
                               npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                               et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
             END IF
          ELSE
             IF ( use_para_diag ) then
                !
                !$acc host_data use_device(et)
                CALL pcegterg_gpu( h_psi_gpu, s_psi_acc, okvan, g_psi, &
                               npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                               et(1, ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                !$acc end host_data
                !
             ELSE
                !
                !$acc host_data use_device(et)
                CALL cegterg ( h_psi_gpu, s_psi_acc, okvan, g_psi, &
                               npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                               et(1, ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                !$acc end host_data 
             END IF
             !$acc update self(et)
          END IF
          !
          avg_iter = avg_iter + dav_iter
          !
          ! ... save wave-functions to be used as input for the
          ! ... iterative diagonalization of the next scf iteration
          ! ... and for rho calculation
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT david_loop
          !
       ENDDO david_loop
       !
    ENDIF
    !
    IF ( lelfield .AND. okvan ) CALL deallocate_bec_type( bec_evcel )
    !
    RETURN
    !
  END SUBROUTINE diag_bands_k
  !
  !-----------------------------------------------------------------------
  FUNCTION test_exit_cond()
    !-----------------------------------------------------------------------
    !! This logical function is .TRUE. when iterative diagonalization
    !! is converged.
    !
    IMPLICIT NONE
    !
    LOGICAL :: test_exit_cond
    !
    
    IF ( lscf .AND. lgcscf ) THEN
       !
       ! ... tight condition for GC-SCF
       !
       test_exit_cond = .NOT. ( ( ntry <= 8 ) .AND. ( notconv > 0 ) )
       !
    ELSE
       !
       test_exit_cond = .NOT. ( ( ntry <= 5 ) .AND. &
            ( ( .NOT. lscf .AND. ( notconv > 0 ) ) .OR. &
            (       lscf .AND. ( notconv > 5 ) ) ) )
       !
    END IF
    !
  END FUNCTION test_exit_cond
  !
END SUBROUTINE diag_bands
!
!----------------------------------------------------------------------------
SUBROUTINE c_bands_efield( iter )
  !----------------------------------------------------------------------------
  !! Driver routine for Hamiltonian diagonalization under an electric field.
  !
  USE noncollin_module,     ONLY : npol
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : nberrycyc, fact_hepsi, &
                                   evcel, evcelp, evcelm, gdir, l3dstring,&
                                   efield, efield_cry
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : nbnd, npwx
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  !! iteration index
  !
  ! ... local variables
  !
  INTEGER :: inberry, ipol, ierr
  !
  !
  ALLOCATE( evcel ( npol*npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcel ', ABS( ierr ) )
  ALLOCATE( evcelm( npol*npwx, nbnd, 3  ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcelm ', ABS( ierr ) )
  ALLOCATE( evcelp( npol*npwx, nbnd, 3 ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcelp ', ABS( ierr ) )
  ALLOCATE( fact_hepsi(nks, 3), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate fact_hepsi ', ABS( ierr ) )
  !
  DO inberry = 1, nberrycyc
     !
     !...set up electric field hermitean operator
     !
     FLUSH(stdout)
     IF (.NOT.l3dstring) THEN
        CALL h_epsi_her_set (gdir, efield)
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_set(ipol, efield_cry(ipol))
        ENDDO
     ENDIF
     FLUSH(stdout)
     !
     CALL c_bands( iter )
     !
  ENDDO
  !
  DEALLOCATE( fact_hepsi )
  DEALLOCATE( evcelp )
  DEALLOCATE( evcelm )
  DEALLOCATE( evcel  )
  !
  RETURN
  !
END SUBROUTINE c_bands_efield
!
!------------------------------------------------------------------------------
SUBROUTINE c_bands_nscf( )
  !----------------------------------------------------------------------------
  !! Driver routine for Hamiltonian diagonalization routines
  !! specialized to non-self-consistent calculations (no electric field).
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE starting_scf,         ONLY : starting_wfc
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, restart, isolve, io_level, iverbosity, use_gpu
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_projectors, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions,        ONLY : evc
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE check_stop,           ONLY : check_stop_now
  USE uspp_init,            ONLY : init_us_2
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(DP) :: avg_iter, ethr_
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios
  ! ik_: k-point already done in a previous run
  ! ik : counter on k points
  LOGICAL :: exst
  !
  REAL(DP), EXTERNAL :: get_clock
  !
  !
  CALL start_clock( 'c_bands' )
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands( ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  !
  DO ik = 1, ik_
     CALL get_buffer( evc, nwordwfc, iunwfc, ik )
  ENDDO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSEIF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")' )
  ELSEIF ( isolve == 3 ) THEN
     WRITE( stdout, '(5X,"ParO style diagonalization")')
  ELSEIF ( isolve == 4 ) THEN
     WRITE( stdout, '(5X,"RMM-DIIS diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve )
  ENDIF
  !
  ! ... For each k point (except those already calculated if restarting)
  ! ... diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) CALL phase_factor(ik)
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     CALL g2_kin( ik )
     ! 
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb , .true.)
     !
     ! ... Needed for DFT+Hubbard
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (Hubbard_projectors.NE.'pseudo') ) THEN
        CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
        !$acc update device(wfcU)
     END IF
     !
     ! ... calculate starting  wavefunctions
     !
     IF ( iverbosity > 0 .AND. npool == 1 ) THEN
        WRITE( stdout, 9001 ) ik, nks
     ELSEIF ( iverbosity > 0 .AND. npool > 1 ) THEN
        WRITE( stdout, 9002 ) ik, nks
     ENDIF
     !
     IF ( TRIM(starting_wfc) == 'file' ) THEN
        !
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        !
     ELSE
        !
        CALL init_wfc( ik )
        !
     ENDIF
     !
     ! ... diagonalization of bands for k-point ik
     !
     CALL diag_bands( 1, ik, avg_iter )
     !$acc update self(evc)
     !
     ! ... save wave-functions (unless disabled in input)
     !
     IF ( io_level > -1 ) CALL save_buffer( evc, nwordwfc, iunwfc, ik )
     !
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     IF (ik <= nkdum) THEN
        !
        ! ... stop requested by user: save restart information,
        ! ... save wavefunctions to file
        !
        IF ( check_stop_now() ) THEN
           CALL save_in_cbands( ik, ethr, avg_iter, et )
           RETURN
        ENDIF
        !
     ENDIF
     !
     ! report about timing
     !
     IF ( iverbosity > 0 ) THEN
        WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
        FLUSH( stdout )
     ENDIF
     !
  ENDDO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, '(/,5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)' ) &
       ethr, avg_iter
  !
  CALL stop_clock( 'c_bands' )
  !
  RETURN
  !
  ! formats
  !
9002 FORMAT(/'     Computing kpt #: ',I5, '  of ',I5,' on this pool' )
9001 FORMAT(/'     Computing kpt #: ',I5, '  of ',I5 )
9000 FORMAT( '     total cpu time spent up to now is ',F10.1,' secs' )
  !
END SUBROUTINE c_bands_nscf

FUNCTION rmm_use_davidson(iter_) RESULT (res)
  USE control_flags, ONLY: rmm_with_davidson
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iter_ 
  LOGICAL :: res 
  res = (rmm_with_davidson) .AND. ( iter_ < 3 .OR. MOD(iter_,5) == 0) 
END FUNCTION rmm_use_davidson

FUNCTION rmm_use_paro(iter_) RESULT (res)
  USE control_flags, ONLY: rmm_with_davidson
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iter_ 
  LOGICAL  :: res 
  res = (.NOT. rmm_with_davidson) .AND.  (MOD(iter_,5) == 1) 
END FUNCTION rmm_use_paro
