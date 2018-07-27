!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band_gpu()
  !----------------------------------------------------------------------------
  !
  ! ... Calculates the symmetrized charge density and related quantities
  ! ... Also computes the occupations and the sum of occupied eigenvalues.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : eband
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft, fft_interpolate
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k, igk_k_d
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions,        ONLY : evc, psic
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE spin_orb,             ONLY : lspinorb, domag
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum
  USE funct,                ONLY : dft_is_meta
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   becp
  USE wavefunctions_gpum, ONLY : evc_d, using_evc, using_evc_d
  USE wvfct_gpum,                ONLY : using_et
  USE uspp_gpum,                 ONLY : becsum_d, ebecsum_d, vkb_d, using_vkb, &
                                        using_vkb_d, using_becsum_d, using_ebecsum_d, &
                                        using_becsum, using_ebecsum
  USE becmod_subs_gpum,          ONLY : using_becp_auto
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ir,   &! counter on 3D r points
             is,   &! counter on spin polarizations
             ig,   &! counter on g vectors
             ibnd, &! counter on bands
             ik,   &! counter on k points
             npol_,&! auxiliary dimension for noncolin case
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  REAL (DP), ALLOCATABLE :: kplusg (:)
  !
  !
  CALL start_clock( 'sum_band' )
  !
  CALL using_becsum_d(2)
  !
  becsum_d(:,:,:) = 0.D0
  if (tqr) CALL using_ebecsum_d(2)
  if (tqr) ebecsum_d(:,:,:) = 0.D0
  rho%of_r(:,:)      = 0.D0
  rho%of_g(:,:)      = (0.D0, 0.D0)
  if ( dft_is_meta() .OR. lxdm ) then
     rho%kin_r(:,:)      = 0.D0
     rho%kin_g(:,:)      = (0.D0, 0.D0)
  end if
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL weights ( )
  !
  IF (one_atom_occupations) CALL new_evc()
  !
  ! ... btype, used in diagonalization, is set here: a band is considered empty
  ! ... and computed with low accuracy only when its occupation is < 0.01, and
  ! ... only if option diago_full_acc is false; otherwise, use full accuracy
  !
  btype(:,:) = 1
  IF ( .NOT. diago_full_acc ) THEN
     !
     FORALL( ik = 1:nks, wk(ik) > 0.D0 )
        WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
     END FORALL
     !
  END IF
  !
  ! ... Needed for LDA+U: compute occupations of Hubbard states
  !
  IF (lda_plus_u) THEN
     IF(noncolin) THEN
        CALL new_ns_nc(rho%ns_nc)
     ELSE
        CALL new_ns(rho%ns)
     ENDIF
  ENDIF
  !
  ! ... for band parallelization: set band computed by this processor
  !
  call divide ( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type (nkb,nbnd, becp,intra_bgrp_comm)
  IF ( okvan ) CALL using_becp_auto(2)
  IF (dft_is_meta() .OR. lxdm) ALLOCATE (kplusg(npwx))
  !
  ! ... specialized routines are called to sum at Gamma or for each k point 
  ! ... the contribution of the wavefunctions to the charge
  ! ... The band energy contribution eband is computed together with the charge
  !
  eband         = 0.D0
  !
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma_gpu()
     !
  ELSE
     !
     CALL sum_band_k_gpu()
     !
  END IF
  CALL mp_sum( eband, inter_pool_comm )
  CALL mp_sum( eband, inter_bgrp_comm )
  !
  IF (dft_is_meta() .OR. lxdm) DEALLOCATE (kplusg)
  IF ( okvan ) CALL deallocate_bec_type ( becp )
  IF ( okvan ) CALL using_becp_auto(2)
  !
  ! ... sum charge density over pools (distributed k-points) and bands
  !
  CALL mp_sum( rho%of_r, inter_pool_comm )
  CALL mp_sum( rho%of_r, inter_bgrp_comm )
  IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
  !
  ! ... bring the unsymmetrized rho(r) to G-space (use psic as work array)
  !
  DO is = 1, nspin
     psic(1:dffts%nnr) = rho%of_r(1:dffts%nnr,is)
     psic(dffts%nnr+1:) = 0.0_dp
     CALL fwfft ('Rho', psic, dffts)
     rho%of_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
     rho%of_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
  END DO

  IF( okvan )  THEN
     !
     ! ... becsum is summed over bands (if bgrp_parallelization is done)
     ! ... and over k-points (but it is not symmetrized)
     !
     CALL using_becsum(1)   ! use host copy to do the comunication. This avoids going back an forth GPU data
     CALL mp_sum(becsum, inter_bgrp_comm )
     CALL mp_sum(becsum, inter_pool_comm )
     CALL using_becsum_d(0)
     !
     ! ... same for ebecsum, a correction to becsum (?) in real space
     !
     IF (tqr) CALL using_ebecsum(1)
     IF (tqr) CALL mp_sum(ebecsum, inter_pool_comm )
     IF (tqr) CALL mp_sum(ebecsum, inter_bgrp_comm )
     IF (tqr) CALL using_ebecsum_d(0)
     !
     ! ... PAW: symmetrize becsum and store it
     ! ... FIXME: the same should be done for USPP as well
     !
     IF ( okpaw ) THEN
        rho%bec(:,:,:) = becsum(:,:,:)
        CALL PAW_symmetrize(rho%bec)
     END IF
     !
     ! ... Here we add the (unsymmetrized) Ultrasoft contribution to the charge
     !
     CALL addusdens_gpu(rho%of_g(:,:))
     !
  ENDIF
  !
  ! ... symmetrize rho(G) 
  !
  CALL sym_rho ( nspin_mag, rho%of_g )
  !
  ! ... synchronize rho%of_r to the calculated rho%of_g (use psic as work array)
  !
  DO is = 1, nspin_mag
     psic(:) = ( 0.D0, 0.D0 )
     psic(dfftp%nl(:)) = rho%of_g(:,is)
     IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho%of_g(:,is) )
     CALL invfft ('Rho', psic, dfftp)
     rho%of_r(:,is) = psic(:)
  END DO
  !
  ! ... rho_kin(r): sum over bands, k-points, bring to G-space, symmetrize,
  ! ... synchronize with rho_kin(G)
  !
  IF ( dft_is_meta() .OR. lxdm) THEN
     !
     CALL mp_sum( rho%kin_r, inter_pool_comm )
     CALL mp_sum( rho%kin_r, inter_bgrp_comm )
     DO is = 1, nspin
        psic(1:dffts%nnr) = rho%kin_r(1:dffts%nnr,is)
        psic(dffts%nnr+1:) = 0.0_dp
        CALL fwfft ('Rho', psic, dffts)
        rho%kin_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
        rho%of_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
     END DO
     !
     IF (.NOT. gamma_only) CALL sym_rho( nspin, rho%kin_g )
     !
     DO is = 1, nspin
        psic(:) = ( 0.D0, 0.D0 )
        psic(dfftp%nl(:)) = rho%kin_g(:,is)
        IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho%kin_g(:,is) )
        CALL invfft ('Rho', psic, dfftp)
        rho%kin_r(:,is) = psic(:)
     END DO
     !
  END IF
  !
  CALL stop_clock( 'sum_band' )
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_gamma_gpu()
       !-----------------------------------------------------------------------
       !
       ! ... gamma version
       !
       USE becmod,        ONLY : becp
       USE mp_bands,      ONLY : me_bgrp
       USE mp,            ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, ioff, ioff_tg, nxyp, incr, v_siz, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi_d(:)
       COMPLEX(DP), ALLOCATABLE :: psic_d(:)
       REAL(DP),    ALLOCATABLE :: tg_rho_d(:), tg_rho_h(:)
       REAL(DP),    ALLOCATABLE :: rho_d(:,:)
       INTEGER,     POINTER     :: dffts_nl_d(:), dffts_nlm_d(:)
       LOGICAL :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp, ierr
#if defined(__CUDA)
       attributes(device) :: psic_d, tg_psi_d, tg_rho_d, rho_d
       attributes(device) :: dffts_nl_d, dffts_nlm_d
       attributes(pinned) :: tg_rho_h
#endif
       !
       CALL using_evc_d(0); CALL using_et(0)
       dffts_nl_d  => dffts%nl_d
       dffts_nlm_d => dffts%nlm_d
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (dft_is_meta() .OR. lxdm) )
       !
       incr = 2

       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg 
          !
          ALLOCATE( tg_psi_d( v_siz ) )
          ALLOCATE( tg_rho_d( v_siz ) )
          ALLOCATE( tg_rho_h( v_siz ) )
          !
          incr  = 2 *  fftx_ntgrp(dffts)
          !
       ELSE
          ALLOCATE( rho_d, MOLD=rho%of_r ) ! OPTIMIZE HERE, use buffers (and batched FFT)
          ALLOCATE(psic_d(dfftp%nnr))
          rho_d = 0.0_DP
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( use_tg ) tg_rho_d = 0.0_DP
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          IF ( nks > 1 ) CALL using_evc(2) ! get_buffer(evc, ...) evc is updated (intent out)
          IF ( nks > 1 ) CALL using_evc_d(0) ! sync on the GPU
          !
          IF ( nkb > 0 ) CALL using_vkb_d(2)
          IF ( nkb > 0 ) CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             ! ... the sum of eband and demet is the integral for
             ! ... e < ef of e n(e) which reduces for degauss=0 to the sum of
             ! ... the eigenvalues.
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
          END DO
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             IF( use_tg ) THEN
                !
                tg_psi_d(:) = ( 0.D0, 0.D0 )
                ioff   = 0
                !
                CALL tg_get_nnr( dffts, right_nnr )
                ntgrp = fftx_ntgrp(dffts)
                !
                DO idx = 1, 2*ntgrp, 2
                   !
                   ! ... 2*ntgrp ffts at the same time
                   !
                   IF( idx + ibnd - 1 < ibnd_end ) THEN
!$cuf kernel do(1) <<<*,*>>>
                      DO j = 1, npw
                         tg_psi_d(dffts_nl_d (j)+ioff)=     evc_d(j,idx+ibnd-1)+&
                              (0.0d0,1.d0) * evc_d(j,idx+ibnd)
                         tg_psi_d(dffts_nlm_d(j)+ioff)=CONJG(evc_d(j,idx+ibnd-1) -&
                              (0.0d0,1.d0) * evc_d(j,idx+ibnd) )
                      END DO
                   ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
!$cuf kernel do(1) <<<*,*>>>
                      DO j = 1, npw
                         tg_psi_d(dffts_nl_d (j)+ioff)=       evc_d(j,idx+ibnd-1)
                         tg_psi_d(dffts_nlm_d(j)+ioff)=CONJG( evc_d(j,idx+ibnd-1) )
                      END DO
                   END IF

                   ioff = ioff + right_nnr

                END DO
                !
                CALL invfft ('tgWave', tg_psi_d, dffts )
                !
                ! Now the first proc of the group holds the first two bands
                ! of the 2*ntgrp bands that we are processing at the same time,
                ! the second proc. holds the third and fourth band
                ! and so on
                !
                ! Compute the proper factor for each band
                !
                idx = fftx_tgpe(dffts) + 1
                !
                ! Remember two bands are packed in a single array :
                ! proc 0 has bands ibnd   and ibnd+1
                ! proc 1 has bands ibnd+2 and ibnd+3
                ! ....
                !
                idx = 2 * idx - 1
                !
                IF( idx + ibnd - 1 < ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = wg( idx + ibnd    , ik) / omega
                ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = w1
                ELSE
                   w1 = 0.0d0
                   w2 = w1
                END IF
                !
                CALL tg_get_group_nr3( dffts, right_nr3 )
                !
                CALL get_rho_gamma_gpu(tg_rho_d, dffts%nr1x*dffts%nr2x*right_nr3, w1, w2, tg_psi_d)
                !
             ELSE
                !
                psic_d(:) = ( 0.D0, 0.D0 )
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   !$cuf kernel do(1)
                   DO j=1,npw
                      psic_d(dffts_nl_d(j))  = evc_d(j,ibnd) + &
                                              ( 0.D0, 1.D0 ) * evc_d(j,ibnd+1)
                      psic_d(dffts_nlm_d(j)) = CONJG( evc_d(j,ibnd) - &
                                              ( 0.D0, 1.D0 ) * evc_d(j,ibnd+1) )
                   END DO
                   !
                ELSE
                   !
                   !$cuf kernel do(1)
                   DO j=1,npw
                      psic_d(dffts_nl_d (j))  = evc_d(j,ibnd)
                      psic_d(dffts_nlm_d(j)) = CONJG( evc_d(j,ibnd) )
                   END DO
                   !
                END IF
                !
                CALL invfft ('Wave', psic_d, dffts)
                !
                w1 = wg(ibnd,ik) / omega
                !
                ! ... increment the charge density ...
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   w2 = wg(ibnd+1,ik) / omega
                   !
                ELSE
                   !
                   w2 = w1
                   !
                END IF
                !
                CALL get_rho_gamma_gpu(rho_d(:,current_spin), dffts%nnr, w1, w2, psic_d(:))
                !
             END IF
             !
             IF (dft_is_meta() .OR. lxdm) THEN
                DO j=1,3
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   kplusg (1:npw) = (xk(j,ik)+g(j,1:npw)) * tpiba

                   IF ( ibnd < ibnd_end ) THEN
                      ! ... two ffts at the same time
                      psic(dffts%nl (1:npw))=CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                            ( evc(1:npw,ibnd) + &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                      psic(dffts%nlm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) - &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   ELSE
                      psic(dffts%nl(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      psic(dffts%nlm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) )
                   END IF
                   !
                   CALL invfft ('Wave', psic, dffts)
                   !
                   ! ... increment the kinetic energy density ...
                   !
                   DO ir = 1, dffts%nnr
                      rho%kin_r(ir,current_spin) = &
                                           rho%kin_r(ir,current_spin) + &
                                           w1 *  DBLE( psic(ir) )**2 + &
                                           w2 * AIMAG( psic(ir) )**2
                   END DO
                   !
                END DO
             END IF
             !
             !
          END DO
          !
          IF( use_tg ) THEN
             tg_rho_h = tg_rho_d
             CALL tg_reduce_rho( rho%of_r, tg_rho_h, current_spin, dffts )
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec_gpu ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       IF( .not. use_tg ) THEN
          rho%of_r = rho_d
       END IF
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL using_becsum(1)
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL mp_sum( becsum, becp%comm )
       IF( okvan .AND. becp%comm /= mp_get_comm_null() .and. tqr ) CALL using_ebecsum(1)
       IF( okvan .AND. becp%comm /= mp_get_comm_null() .and. tqr ) CALL mp_sum( ebecsum, becp%comm )
       !
       IF( use_tg ) THEN
          DEALLOCATE( tg_psi_d )
          DEALLOCATE( tg_rho_d )
          DEALLOCATE( tg_rho_h )
       ELSE
          DEALLOCATE(rho_d)
          DEALLOCATE(psic_d)
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma_gpu
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k_gpu()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE wavefunctions_gpum, ONLY : psic_nc_d
       USE mp_bands,     ONLY : me_bgrp
       USE mp,           ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, ioff, ioff_tg, nxyp, incr, v_siz, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi_d(:), tg_psi_nc_d(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho_d(:), tg_rho_nc_d(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho_h(:), tg_rho_nc_h(:,:)
       REAL(DP),    ALLOCATABLE :: rho_d(:,:)
       COMPLEX(DP), ALLOCATABLE :: psic_d(:)
       INTEGER,     POINTER     :: dffts_nl_d(:)
       LOGICAL  :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp, ierr
       !
#if defined(__CUDA)
       attributes(device) :: psic_d, tg_psi_d, tg_rho_d, tg_psi_nc_d, tg_rho_nc_d
       attributes(device) :: rho_d, dffts_nl_d
       attributes(pinned) :: tg_rho_h, tg_rho_nc_h
#endif
       !
       CALL using_evc(0); CALL using_et(0)
       dffts_nl_d => dffts%nl_d
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (dft_is_meta() .OR. lxdm) )
       !
       incr = 1
       !
       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg
          !
          IF (noncolin) THEN
             ALLOCATE( tg_psi_nc_d( v_siz, npol ) )
             ALLOCATE( tg_rho_nc_d( v_siz, nspin_mag ) )
             ALLOCATE( tg_rho_nc_h( v_siz, nspin_mag ) )
          ELSE
             ALLOCATE( tg_psi_d( v_siz ) )
             ALLOCATE( tg_rho_d( v_siz ) )
             ALLOCATE( tg_rho_h( v_siz ) )
          ENDIF
          !
          incr  = fftx_ntgrp(dffts)
          !
       ELSE
          ALLOCATE(rho_d, MOLD=rho%of_r) ! OPTIMIZE HERE, use buffers!
          ALLOCATE(psic_d(dffts%nnr))
          ! This is used as reduction variable on the device
          rho_d = 0.0_DP
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF( use_tg ) THEN
            IF (noncolin) THEN
               tg_rho_nc_d = 0.0_DP
            ELSE
               tg_rho_d = 0.0_DP
            ENDIF
          ENDIF

          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          IF ( nks > 1 ) CALL using_evc(2)
          IF ( nks > 1 ) CALL using_evc_d(0)  ! sync evc on GPU, OPTIMIZE (use async here)
          !
          IF ( nkb > 0 ) CALL using_vkb_d(2)
          IF ( nkb > 0 ) &
             CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             IF( use_tg ) THEN
                DO idx = 1, fftx_ntgrp(dffts)
                   IF( idx + ibnd - 1 <= ibnd_end ) eband = eband + et( idx + ibnd - 1, ik ) * wg( idx + ibnd - 1, ik )
                END DO
             ELSE
                eband = eband + et( ibnd, ik ) * wg( ibnd, ik )
             END IF
             !
             ! ... the sum of eband and demet is the integral for e < ef of
             ! ... e n(e) which reduces for degauss=0 to the sum of the
             ! ... eigenvalues
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
                IF( use_tg ) THEN
                   !
                   tg_psi_nc_d = ( 0.D0, 0.D0 )
                   !
                   CALL tg_get_nnr( dffts, right_nnr )
                   ntgrp = fftx_ntgrp( dffts )
                   !
                   ioff   = 0
                   !
                   DO idx = 1, ntgrp
                      !
                      ! ... ntgrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= ibnd_end ) THEN
!$cuf kernel do(1)
                         DO j = 1, npw
                            tg_psi_nc_d( dffts_nl_d(igk_k_d(j,ik) ) + ioff, 1 ) = &
                                                       evc_d( j, idx+ibnd-1 )
                            tg_psi_nc_d( dffts_nl_d(igk_k_d(j,ik) ) + ioff, 2 ) = &
                                                       evc_d( j+npwx, idx+ibnd-1 )
                         END DO
                      END IF

                      ioff = ioff + right_nnr

                   END DO
                   !
                   CALL invfft ('tgWave', tg_psi_nc_d(:,1), dffts)
                   CALL invfft ('tgWave', tg_psi_nc_d(:,2), dffts)
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the ntgrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   idx = fftx_tgpe(dffts) + 1 
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   END IF
                   !
                   CALL tg_get_group_nr3( dffts, right_nr3 )
                   !
                   ! OPTIMIZE HERE : this is a sum of all densities in first spin channel
                   DO ipol=1,npol
                      CALL get_rho_gpu(tg_rho_nc_d(:,1), dffts%nr1x * dffts%nr2x* right_nr3, w1, tg_psi_nc_d(:,ipol))
                   ENDDO
                   !
                   IF (domag) CALL get_rho_domag_gpu(tg_rho_nc_d(:,:), dffts%nr1x*dffts%nr2x*dffts%my_nr3p, w1, tg_psi_nc_d(:,:))
                   !
                ELSE
!
!     Noncollinear case without task groups
!
                   psic_nc_d = (0.D0,0.D0)

!$cuf kernel do(1)
                   DO j = 1, npw
                      psic_nc_d( dffts_nl_d(igk_k_d(j,ik) ), 1 ) = &
                                                 evc_d( j, ibnd )
                      psic_nc_d( dffts_nl_d(igk_k_d(j,ik) ), 2 ) = &
                                                 evc_d( j+npwx, ibnd )
                   END DO
                   !
                   CALL invfft ('Wave', psic_nc_d(:,1), dffts)
                   CALL invfft ('Wave', psic_nc_d(:,2), dffts)
                   !
                   ! increment the charge density ...
                   !
                   DO ipol=1,npol
                      CALL get_rho_gpu(rho_d(:,1), dffts%nnr, w1, psic_nc_d(:,ipol))
                   END DO
                   !
                   ! In this case, calculate also the three
                   ! components of the magnetization (stored in rho%of_r(ir,2-4))
                   !
                   IF (domag) THEN
                      CALL get_rho_domag_gpu(rho_d(1:,1:), dffts%nnr, w1, psic_nc_d(1:,1:))
                   ELSE
                      rho_d(:,2:4)=0.0_DP  ! OPTIMIZE HERE: this memset can be avoided
                   END IF
                   !
                END IF
                !
             ELSE
                !
                IF( use_tg ) THEN
                   !
!$cuf kernel do(1)
                   DO j = 1, SIZE( tg_psi_d )
                      tg_psi_d(j) = ( 0.D0, 0.D0 )
                   END DO
                   !
                   ioff   = 0
                   !
                   CALL tg_get_nnr( dffts, right_nnr )
                   ntgrp = fftx_ntgrp( dffts )
                   !
                   DO idx = 1, ntgrp
                      !
                      ! ... ntgrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= ibnd_end ) THEN
!$cuf kernel do(1)
                         DO j = 1, npw
                            tg_psi_d( dffts_nl_d(igk_k_d(j,ik))+ioff ) = evc_d(j,idx+ibnd-1)
                         END DO
                      END IF

                      ioff = ioff + right_nnr

                   END DO
                   !
                   CALL invfft ('tgWave', tg_psi_d, dffts)
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the ntgrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   idx = fftx_tgpe(dffts) + 1
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   END IF
                   !
                   CALL tg_get_group_nr3( dffts, right_nr3 )
                   !
                   CALL get_rho_gpu(tg_rho_d, dffts%nr1x * dffts%nr2x * right_nr3, w1, tg_psi_d)
                   !
                ELSE
                   !
                   psic_d(:) = ( 0.D0, 0.D0 )
                   !
                   !$cuf kernel do(1)
                   DO j = 1, npw
                      psic_d(dffts_nl_d(igk_k_d(j,ik))) = evc_d(j,ibnd)
                   END DO
                   !
                   CALL invfft ('Wave', psic_d, dffts)
                   !
                   ! ... increment the charge density ...
                   !
                   CALL get_rho_gpu(rho_d(:,current_spin), dffts%nnr, w1, psic_d(:))

                END IF
                !
                IF (dft_is_meta() .OR. lxdm) THEN
                   DO j=1,3
                      psic(:) = ( 0.D0, 0.D0 )
                      !
                      kplusg (1:npw) = (xk(j,ik)+g(j,igk_k(1:npw,ik))) * tpiba
                      psic(dffts_nl_d(igk_k(1:npw,ik)))=CMPLX(0d0,kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      !
                      CALL invfft ('Wave', psic, dffts)
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      CALL get_rho(rho%kin_r(:,current_spin), dffts%nnr, w1, psic)
                   END DO
                END IF
                !
             END IF
             !
          END DO
          !
          IF( use_tg ) THEN
             !
             ! reduce the charge across task group
             !
             IF (noncolin)       tg_rho_nc_h = tg_rho_nc_d
             IF (.not. noncolin) tg_rho_h    = tg_rho_d
             CALL tg_reduce_rho( rho%of_r, tg_rho_nc_h, tg_rho_h, current_spin, noncolin, domag, dffts )
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec_gpu ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       IF (.not. use_tg ) THEN
          rho%of_r = rho_d
       END IF
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL using_becsum(1)
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL mp_sum( becsum, becp%comm )
       IF( okvan .AND. becp%comm /= mp_get_comm_null() .and. tqr ) CALL using_ebecsum(1)
       IF( okvan .AND. becp%comm /= mp_get_comm_null() .and. tqr ) CALL mp_sum( ebecsum, becp%comm )
       !
       IF( use_tg ) THEN
          IF (noncolin) THEN
             DEALLOCATE( tg_psi_nc_d )
             DEALLOCATE( tg_rho_nc_d )
             DEALLOCATE( tg_rho_nc_h )
          ELSE
             DEALLOCATE( tg_psi_d )
             DEALLOCATE( tg_rho_d )
             DEALLOCATE( tg_rho_h )
          END IF
       ELSE
          DEALLOCATE(rho_d) ! OPTIMIZE HERE, use buffers!
          DEALLOCATE(psic_d) ! OPTIMIZE HERE, use buffers!
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k_gpu
     !
     !
     SUBROUTINE get_rho_gpu(rho_loc_d, nrxxs_loc, w1_loc, psic_loc_d)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(:)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_d(:)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir

!$cuf kernel do(1)
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir) = rho_loc_d(ir) + &
                         w1_loc * ( DBLE( psic_loc_d(ir) )**2 + &
                                   AIMAG( psic_loc_d(ir) )**2 )
           !
        END DO

     END SUBROUTINE get_rho_gpu
     !
     SUBROUTINE get_rho(rho_loc_h, nrxxs_loc, w1_loc, psic_loc_h)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_h(nrxxs_loc)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_h(nrxxs_loc)
        INTEGER :: ir

        DO ir = 1, nrxxs_loc
           !
           rho_loc_h(ir) = rho_loc_h(ir) + &
                         w1_loc * ( DBLE( psic_loc_h(ir) )**2 + &
                                   AIMAG( psic_loc_h(ir) )**2 )
           !
        END DO

     END SUBROUTINE get_rho

     SUBROUTINE get_rho_gamma_gpu(rho_loc_d, nrxxs_loc, w1_loc, w2_loc, psic_loc_d)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP) :: psic_loc_d(nrxxs_loc)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir

!$cuf kernel do(1)
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir) = rho_loc_d(ir) + &
                         w1_loc * DBLE( psic_loc_d(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc_d(ir) )**2
           !
        END DO

     END SUBROUTINE get_rho_gamma_gpu


     SUBROUTINE get_rho_domag_gpu(rho_loc_d, nrxxs_loc, w1_loc, psic_loc_d)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(:, :)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_d(:, :)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir

!$cuf kernel do(1)
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir,2) = rho_loc_d(ir,2) + w1_loc*2.D0* &
                          (DBLE(psic_loc_d(ir,1))* DBLE(psic_loc_d(ir,2)) + &
                          AIMAG(psic_loc_d(ir,1))*AIMAG(psic_loc_d(ir,2)))
 
           rho_loc_d(ir,3) = rho_loc_d(ir,3) + w1_loc*2.D0* &
                          (DBLE(psic_loc_d(ir,1))*AIMAG(psic_loc_d(ir,2)) - &
                           DBLE(psic_loc_d(ir,2))*AIMAG(psic_loc_d(ir,1)))

           rho_loc_d(ir,4) = rho_loc_d(ir,4) + w1_loc* &
                          (DBLE(psic_loc_d(ir,1))**2+AIMAG(psic_loc_d(ir,1))**2 &
                          -DBLE(psic_loc_d(ir,2))**2-AIMAG(psic_loc_d(ir,2))**2)
           !
        END DO

     END SUBROUTINE get_rho_domag_gpu

END SUBROUTINE sum_band_gpu

!----------------------------------------------------------------------------
SUBROUTINE sum_bec_gpu ( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
  !----------------------------------------------------------------------------
  !
  ! This routine computes the sum over bands
  !     \sum_i <\psi_i|\beta_l>w_i<\beta_m|\psi_i>
  ! for point "ik" and, for LSDA, spin "current_spin" 
  ! Calls calbec to compute "becp"=<beta_m|psi_i> 
  ! Output is accumulated (unsymmetrized) into "becsum", module "uspp"
  !
  ! Routine used in sum_band (if okvan) and in compute_becsum, called by hinit1 (if okpaw)
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE kinds,         ONLY : DP
  USE becmod,        ONLY : becp, calbec, allocate_bec_type
  USE control_flags, ONLY : gamma_only, tqr
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp,          ONLY : nkb, becsum, ebecsum, indv_ijkb0
  USE uspp_param,    ONLY : upf, nh, nhm
  USE wvfct,         ONLY : nbnd, wg, et, current_k
  USE klist,         ONLY : ngk, nkstot
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions, ONLY : evc
  USE realus,        ONLY : real_space, &
                            invfft_orbital_gamma, calbec_rs_gamma, &
                            invfft_orbital_k, calbec_rs_k
  USE us_exx,        ONLY : store_becxx0
  USE mp_bands,      ONLY : nbgrp,inter_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  ! GPU modules here
  USE wavefunctions_gpum, ONLY : evc_d, using_evc, using_evc_d
  USE wvfct_gpum,                ONLY : et_d, using_et, using_et_d
  USE uspp_gpum,                 ONLY : using_indv_ijkb0, &
                                        using_becsum, using_ebecsum
  USE becmod_subs_gpum,          ONLY : calbec_gpu, using_becp_auto, using_becp_d_auto
  USE becmod_gpum,               ONLY : becp_d
  USE uspp_gpum,                 ONLY : vkb_d, becsum_d, ebecsum_d, indv_ijkb0_d
  USE uspp_gpum,                 ONLY : using_vkb_d, using_becsum_d, using_ebecsum_d, using_indv_ijkb0_d
  !
  ! Used to avoid unnecessary memcopy
  USE funct,                     ONLY : dft_is_hybrid
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1_d(:,:), auxk2_d(:,:), aux_nc_d(:,:)
  REAL(DP), ALLOCATABLE    :: auxg_d(:,:), aux_gk_d(:,:), aux_egk_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: becsum_nc_d, auxk1_d, auxk2_d, aux_nc_d
  attributes(DEVICE) :: auxg_d, aux_gk_d, aux_egk_d
#endif
  INTEGER :: ibnd, ibnd_loc, nbnd_loc  ! counters on bands
  INTEGER :: npw, ikb, jkb, ih, jh, ijh, na, np, is, js, nhnt
  ! counters on beta functions, atoms, atom types, spin, and auxiliary vars
  !
  REAL(DP), POINTER :: becp_d_r_d(:,:)
  COMPLEX(DP), POINTER :: becp_d_k_d(:,:), becp_d_nc_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: becp_d_r_d, becp_d_k_d, becp_d_nc_d
#endif
  REAL(DP), ALLOCATABLE :: wg_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: wg_d
#endif
  !
  ALLOCATE(wg_d, SOURCE=wg)   !!! OPTIMIZE HERE: move this to duplicated modules?
  CALL using_indv_ijkb0_d(0)
  CALL using_becsum_d(1)
  if (tqr) CALL using_ebecsum_d(1)
  !
  npw = ngk(ik)
  IF ( .NOT. real_space ) THEN
     CALL using_evc_d(0); CALL using_vkb_d(0); CALL using_becp_d_auto(2)
     ! calbec computes becp = <vkb_i|psi_j>
     CALL calbec_gpu( npw, vkb_d, evc_d, becp_d )
  ELSE
     CALL using_evc(0); CALL using_becp_auto(2)
     if (gamma_only) then
        do ibnd = ibnd_start, ibnd_end, 2
           call invfft_orbital_gamma(evc,ibnd,ibnd_end) 
           call calbec_rs_gamma(ibnd,ibnd_end,becp%r)
        enddo
        call mp_sum(becp%r,inter_bgrp_comm)
     else
        current_k = ik
        becp%k = (0.d0,0.d0)
        do ibnd = ibnd_start, ibnd_end
           call invfft_orbital_k(evc,ibnd,ibnd_end) 
           call calbec_rs_k(ibnd,ibnd_end)
        enddo
       call mp_sum(becp%k,inter_bgrp_comm)
     endif
  ENDIF
  !
  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  IF(dft_is_hybrid()) THEN       ! This if condition is not present in the CPU code!! Add it?
     CALL using_becp_auto(0)
     CALL store_becxx0(ik, becp)
  ENDIF
  !
  CALL start_clock( 'sum_band:becsum' )
  CALL using_becp_d_auto(0)
  !
  DO np = 1, ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ! allocate work space used to perform GEMM operations
        !
        IF ( gamma_only ) THEN
           nbnd_loc = becp_d%nbnd_loc
           ALLOCATE( auxg_d( nbnd_loc, nh(np) ) )
        ELSE
           ALLOCATE( auxk1_d( ibnd_start:ibnd_end, nh(np)*npol ), &
                     auxk2_d( ibnd_start:ibnd_end, nh(np)*npol ) )
        END IF
        IF ( noncolin ) THEN
           ALLOCATE ( aux_nc_d( nh(np)*npol,nh(np)*npol ) ) 
        ELSE
           ALLOCATE ( aux_gk_d( nh(np),nh(np) ) ) 
           if (tqr) ALLOCATE ( aux_egk_d( nh(np),nh(np) ) ) 
        END IF
        !
        !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
        !   run from index i=indv_ijkb0(na)+1 to i=indv_ijkb0(na)+nh(nt)
        !
        nhnt = nh(np)
        DO na = 1, nat
           !
           IF (ityp(na)==np) THEN
              !
              ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
              ! copy into aux1, aux2 the needed data to perform a GEMM
              !
              IF ( noncolin ) THEN
                 !
                 becp_d_nc_d => becp_d%nc_d
!$cuf kernel do(2)
                 DO is = 1, npol
                    DO ih = 1, nhnt
                       ikb = indv_ijkb0_d(na) + ih
                       DO ibnd = ibnd_start, ibnd_end
                          auxk1_d(ibnd,ih+(is-1)*nhnt)= becp_d_nc_d(ikb,is,ibnd)
                          auxk2_d(ibnd,ih+(is-1)*nhnt)= wg_d(ibnd,ik) * &
                                                        becp_d_nc_d(ikb,is,ibnd)
                       END DO
                    END DO
                 END DO
                 !
                 CALL cublasZgemm ( 'C', 'N', npol*nhnt, npol*nhnt, this_bgrp_nbnd, &
                      (1.0_dp,0.0_dp), auxk1_d, this_bgrp_nbnd, auxk2_d, this_bgrp_nbnd, &
                      (0.0_dp,0.0_dp), aux_nc_d, npol*nhnt )
                 !
              ELSE IF ( gamma_only ) THEN
                 !
                 becp_d_r_d => becp_d%r_d
!$cuf kernel do(2)
                DO ih = 1, nhnt
                    DO ibnd_loc = 1, nbnd_loc
                       ikb = indv_ijkb0_d(na) + ih
                       ibnd = ibnd_loc + becp_d%ibnd_begin - 1
                       auxg_d(ibnd_loc,ih) = becp_d_r_d(ikb,ibnd_loc) * wg_d(ibnd,ik)
                    END DO
                 END DO
                 !
                 ! NB: band parallelizazion has not been performed in this case because 
                 !     bands were already distributed across R&G processors.
                 !     Contribution to aux_gk is scaled by 1.d0/nbgrp so that the becsum
                 !     summation across bgrps performed outside will gives the right value.
                 !
                 CALL cublasDgemm ( 'N', 'N', nhnt, nhnt, nbnd_loc, &
                      1.0_dp/nbgrp, becp_d%r_d(indv_ijkb0(na)+1,1), nkb,    &
                      auxg_d, nbnd_loc, 0.0_dp, aux_gk_d, nhnt )
               if (tqr) then
                 CALL using_et_d(0)
!$cuf kernel do(1)
                 DO ih = 1, nhnt
                    ikb = indv_ijkb0_d(na) + ih
                    DO ibnd_loc = 1, nbnd_loc
                    auxg_d(ibnd_loc,ih) = et_d(ibnd_loc,ik) * auxg_d(ibnd_loc,ih)
                    END DO
                 END DO

                 CALL cublasDgemm ( 'N', 'N', nhnt, nhnt, nbnd_loc, &
                      1.0_dp/nbgrp, becp_d%r_d(indv_ijkb0(na)+1,1), nkb,    &
                      auxg_d, nbnd_loc, 0.0_dp, aux_egk_d, nhnt )
               end if
                 !
              ELSE
                 !
                 becp_d_k_d => becp_d%k_d
!$cuf kernel do(2) <<<*,*>>>
                 DO ih = 1, nhnt
                    DO ibnd = ibnd_start, ibnd_end
                       ikb = indv_ijkb0_d(na) + ih
                       auxk1_d(ibnd,ih) = becp_d_k_d(ikb,ibnd) 
                       auxk2_d(ibnd,ih) = wg_d(ibnd,ik)*becp_d_k_d(ikb,ibnd)
                    END DO
                 END DO
                 !
                 ! only the real part is computed
                 !
                 CALL cublasDgemm ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1_d, 2*this_bgrp_nbnd, auxk2_d, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_gk_d, nhnt )
                 !
               if (tqr) then
                 CALL using_et_d(0)
!$cuf kernel do(2)
                 DO ih = 1, nhnt
                    DO ibnd = ibnd_start, ibnd_end
                       ikb = indv_ijkb0_d(na) + ih
                       auxk2_d(ibnd,ih) = et_d(ibnd,ik)*auxk2_d(ibnd,ih)
                    END DO
                 END DO

                 CALL cublasDgemm ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1_d, 2*this_bgrp_nbnd, auxk2_d, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_egk_d, nhnt )
               end if

              END IF
              !
              ! copy output from GEMM into desired format
              !
              IF (noncolin .AND. .NOT. upf(np)%has_so) THEN
                 CALL add_becsum_nc_gpu (na, np, aux_nc_d, becsum_d )
              ELSE IF (noncolin .AND. upf(np)%has_so) THEN
                 CALL add_becsum_so_gpu (na, np, aux_nc_d, becsum_d )
              ELSE
                 !
!$cuf kernel do(2) <<<*,*>>>
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ijh = jh + ((ih-1)*(2*nhnt-ih))/2  ! or use  ijtoh_d(ih,jh,np) ?  OPTIMIZE !!
                       !
                       ! nondiagonal terms summed and collapsed into a
                       ! single index (matrix is symmetric wrt (ih,jh))
                       !
                       IF ( jh == ih ) THEN
                          becsum_d(ijh,na,current_spin) = &
                               becsum_d(ijh,na,current_spin) + aux_gk_d (ih,jh)
                          if (tqr) ebecsum_d(ijh,na,current_spin) = &
                             ebecsum_d(ijh,na,current_spin) + aux_egk_d (ih,jh)
                       ELSE IF ( jh > ih ) THEN
                          becsum_d(ijh,na,current_spin) = &
                               becsum_d(ijh,na,current_spin) + aux_gk_d(ih,jh)*2.0_dp
                          if (tqr) ebecsum_d(ijh,na,current_spin) = &
                             ebecsum_d(ijh,na,current_spin) + aux_egk_d(ih,jh)*2.0_dp
                       END IF
                    END DO
                 END DO
                 !
              END IF
           END IF
           !
        END DO
        !
        IF ( noncolin ) THEN
           DEALLOCATE ( aux_nc_d )
        ELSE
           DEALLOCATE ( aux_gk_d  ) 
           if (tqr) DEALLOCATE ( aux_egk_d  ) 
        END IF
        IF ( gamma_only ) THEN
           DEALLOCATE( auxg_d )
        ELSE
           DEALLOCATE( auxk2_d, auxk1_d )
        END IF
        !
     END IF
     !
  END DO
  !
  DEALLOCATE(wg_d)
  !
  CALL stop_clock( 'sum_band:becsum' )
  !
END SUBROUTINE sum_bec_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_nc_gpu ( na, np, becsum_nc_d, becsum_d )
!----------------------------------------------------------------------------
  !
  ! This routine multiplies becsum_nc by the identity and the Pauli matrices,
  ! saves it in becsum for the calculation of augmentation charge and
  ! magnetization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE spin_orb,             ONLY : domag
  !
  USE uspp_gpum,            ONLY : ijtoh_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc_d(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum_d(nhm*(nhm+1)/2,nat,nspin_mag)
#if defined(__CUDA)
  attributes(DEVICE) :: becsum_nc_d, becsum_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, ijh, ipol, jpol, nhnp
  REAL(DP) :: fac
  !
  nhnp = nh(np)

!$cuf kernel do(2) <<<*,*>>>
  DO ih = 1, nhnp
     DO jh = 1, nhnp
        IF ( jh >= ih ) THEN
           !ijh = jh + ((ih-1)*(2*nhnp-ih))/2  is this faster? Does it matter?
           ijh=ijtoh_d(ih,jh,np)
           IF ( ih == jh ) THEN
              fac = 1.0_dp
           ELSE
              fac = 2.0_dp
           END IF
           becsum_d(ijh,na,1)= becsum_d(ijh,na,1) + fac * &
                   DBLE( becsum_nc_d(ih,1,jh,1) + becsum_nc_d(ih,2,jh,2) )
           IF (domag) THEN
              becsum_d(ijh,na,2)= becsum_d(ijh,na,2) + fac *  &
                   DBLE( becsum_nc_d(ih,1,jh,2) + becsum_nc_d(ih,2,jh,1) )
              becsum_d(ijh,na,3)= becsum_d(ijh,na,3) + fac * DBLE( (0.d0,-1.d0)* &
                  (becsum_nc_d(ih,1,jh,2) - becsum_nc_d(ih,2,jh,1)) )
              becsum_d(ijh,na,4)= becsum_d(ijh,na,4) + fac * &
                   DBLE( becsum_nc_d(ih,1,jh,1) - becsum_nc_d(ih,2,jh,2) )
           END IF
        END IF
     END DO
  END DO
  
END SUBROUTINE add_becsum_nc_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_so_gpu( na, np, becsum_nc_d, becsum_d )
  !----------------------------------------------------------------------------
  !
  ! This routine multiplies becsum_nc by the identity and the Pauli matrices,
  ! rotates it as appropriate for the spin-orbit case, saves it in becsum
  ! for the calculation of augmentation charge and magnetization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE spin_orb,             ONLY : domag
  !
  USE uspp_gpum,            ONLY : ijtoh_d, nhtol_d, nhtoj_d, indv_d
  USE spin_orb_gpum,        ONLY : fcoef_d
  !
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc_d(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum_d(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, lh, kh, ijh, is1, is2, nhnt
  COMPLEX(DP) :: fac

#if defined(__CUDA)
  attributes (DEVICE) :: becsum_nc_d, becsum_d
#endif
  !
  nhnt = nh(np)
  !
  !$cuf kernel do(1)
  DO ih = 1, nhnt
     DO jh = 1, nhnt
        ijh=ijtoh_d(ih,jh,np)
        DO kh = 1, nhnt
           IF ( (nhtol_d(kh,np)==nhtol_d(ih,np)).AND. &
                (ABS(nhtoj_d(kh,np)-nhtoj_d(ih,np))<1.d8).AND. &
                (indv_d(kh,np)==indv_d(ih,np)) ) THEN ! same_lj(kh,ih,np)
              DO lh=1,nhnt
                 IF ( (nhtol_d(lh,np)==nhtol_d(jh,np)).AND. &
                      (ABS(nhtoj_d(lh,np)-nhtoj_d(jh,np))<1.d8).AND. &
                      (indv_d(lh,np)==indv_d(jh,np)) ) THEN   !same_lj(lh,jh,np)) THEN
                    DO is1=1,npol
                       DO is2=1,npol
                          fac=becsum_nc_d(kh,is1,lh,is2)
                          becsum_d(ijh,na,1)=becsum_d(ijh,na,1) + DBLE( fac * &
                               (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,1,is2,np) + &
                                fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,2,is2,np)  ) )
                          IF (domag) THEN
                            becsum_d(ijh,na,2)=becsum_d(ijh,na,2) + DBLE( fac * &
                                (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,2,is2,np) +&
                                 fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,1,is2,np)  ) )
                            becsum_d(ijh,na,3)=becsum_d(ijh,na,3) + DBLE( fac*(0.d0,-1.d0)*&
                               (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,2,is2,np) - &
                                fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,1,is2,np)  ))
                            becsum_d(ijh,na,4)=becsum_d(ijh,na,4) + DBLE(fac * &
                               (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,1,is2,np) - &
                                fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,2,is2,np)  ) )
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END DO

END SUBROUTINE add_becsum_so_gpu
