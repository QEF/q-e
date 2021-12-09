!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band()
  !----------------------------------------------------------------------------
  !! Calculates the symmetrized charge density and related quantities.  
  !! Also computes the occupations and the sum of occupied eigenvalues.
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : eband
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho, rhoz_or_updw
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, nhtol, nhtoj, indv, okvan, &
                                   becsum_d, ebecsum_d
  USE uspp_param,           ONLY : nh, nhm
  USE wavefunctions,        ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum
  USE xc_lib,               ONLY : xclib_dft_is
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   becp
  USE gcscf_module,         ONLY : lgcscf, gcscf_calc_nelec
  USE io_global,            ONLY : stdout
  USE add_dmft_occ,         ONLY : dmft, dmft_updated, v_dmft
  USE wavefunctions_gpum,   ONLY : using_evc
  USE wvfct_gpum,           ONLY : using_et
  USE becmod_subs_gpum,     ONLY : using_becp_auto
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
             nt,   &! counter on atomic types
             npol_,&! auxiliary dimension for noncolin case
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  REAL (DP), ALLOCATABLE :: kplusg (:)
  !
  !
  CALL start_clock( 'sum_band' )
  !
  if ( nhm > 0 ) then
     becsum(:,:,:) = 0.D0
     if (tqr) ebecsum(:,:,:) = 0.D0
  end if
  rho%of_r(:,:)      = 0.D0
  rho%of_g(:,:)      = (0.D0, 0.D0)
  if ( xclib_dft_is('meta') .OR. lxdm ) then
     rho%kin_r(:,:)      = 0.D0
     rho%kin_g(:,:)      = (0.D0, 0.D0)
  end if
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL start_clock( 'sum_band:weights' )
  !
  ! ... for DMFT skip weights in the first iteration since they were loaded from file
  ! ... and are manipulated elsewhere
  !
  IF (.NOT. ( dmft .AND. .NOT. dmft_updated) ) THEN
     CALL weights ( )
  ENDIF
  !
  CALL stop_clock( 'sum_band:weights' )
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
  ! ... Needed for DFT+U(+V): compute occupations of Hubbard states
  !
  IF (lda_plus_u) THEN
    IF (lda_plus_u_kind.EQ.0) THEN
       !
       CALL new_ns(rho%ns)
       !
       DO nt = 1, ntyp
          IF (is_hubbard_back(nt)) CALL new_nsb(rho%nsb)
       ENDDO
       !
    ELSEIF (lda_plus_u_kind.EQ.1) THEN
       !
       IF (noncolin) THEN
          CALL new_ns_nc(rho%ns_nc)
       ELSE
          CALL new_ns(rho%ns)
       ENDIF
       !
    ELSEIF (lda_plus_u_kind.EQ.2) THEN 
       !
       CALL new_nsg()
       !
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
  IF ( okvan ) CALL allocate_bec_type (nkb, this_bgrp_nbnd, becp, intra_bgrp_comm)
  IF ( okvan ) CALL using_becp_auto(2)
  IF (xclib_dft_is('meta') .OR. lxdm) ALLOCATE (kplusg(npwx))
  !
  ! ... specialized routines are called to sum at Gamma or for each k point 
  ! ... the contribution of the wavefunctions to the charge
  ! ... The band energy contribution eband is computed together with the charge
  !
  eband         = 0.D0
  !
  CALL start_clock( 'sum_band:loop' )
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma()
     !
  ELSE
     !
     CALL sum_band_k()
     !
  END IF
  CALL stop_clock( 'sum_band:loop' )
  CALL mp_sum( eband, inter_pool_comm )
  CALL mp_sum( eband, inter_bgrp_comm )
  !
  IF (xclib_dft_is('meta') .OR. lxdm) DEALLOCATE (kplusg)
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
     CALL mp_sum(becsum, inter_bgrp_comm )
     CALL mp_sum(becsum, inter_pool_comm )
     !
     ! ... same for ebecsum, a correction to becsum (?) in real space
     !
     IF (tqr) CALL mp_sum(ebecsum, inter_pool_comm )
     IF (tqr) CALL mp_sum(ebecsum, inter_bgrp_comm )
     !
#if defined __CUDA
     if (nhm>0) then
        becsum_d=becsum
        if (tqr) ebecsum_d=ebecsum
     endif
#endif
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
     CALL addusdens(rho%of_g(:,:))
     !
  ENDIF
  !
  ! ... symmetrize rho(G) 
  !
  CALL start_clock( 'sum_band:sym_rho' )
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
  IF ( xclib_dft_is('meta') .OR. lxdm) THEN
     !
     CALL mp_sum( rho%kin_r, inter_pool_comm )
     CALL mp_sum( rho%kin_r, inter_bgrp_comm )
     DO is = 1, nspin
        psic(1:dffts%nnr) = rho%kin_r(1:dffts%nnr,is)
        psic(dffts%nnr+1:) = 0.0_dp
        CALL fwfft ('Rho', psic, dffts)
        rho%kin_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
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
  CALL stop_clock( 'sum_band:sym_rho' )
  !
  ! ... if LSDA rho%of_r and rho%of_g are converted from (up,dw) to
  ! ... (up+dw,up-dw) format.
  !
  IF ( nspin == 2 ) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
  !
  ! ... sum number of electrons, for GC-SCF
  !
  IF ( lgcscf ) CALL gcscf_calc_nelec()
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
     SUBROUTINE sum_band_gamma()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for gamma version.
       !
       USE becmod,        ONLY : becp
       USE mp_bands,      ONLY : me_bgrp
       USE mp,            ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines, ONLY : fftx_ntgrp, fftx_tgpe, &
                          tg_reduce_rho, tg_get_nnr, tg_get_group_nr3
       USE uspp_init,     ONLY : init_us_2
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, ioff, ioff_tg, nxyp, incr, v_siz, j, ir3
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:)
       LOGICAL :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp
       !
       CALL using_evc(0); CALL using_et(0)
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (xclib_dft_is('meta') .OR. lxdm) )
       !
       incr = 2

       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg 
          !
          ALLOCATE( tg_psi( v_siz ) )
          ALLOCATE( tg_rho( v_siz ) )
          !
          incr  = 2 *  fftx_ntgrp(dffts)
          !
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( use_tg ) tg_rho = 0.0_DP
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          CALL start_clock( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )

          IF ( nks > 1 ) CALL using_evc(1) ! get_buffer(evc, ...) evc is updated (intent out)

          CALL stop_clock( 'sum_band:buffer' )
          !
          CALL start_clock( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          CALL stop_clock( 'sum_band:init_us_2' )
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
                tg_psi(:) = ( 0.D0, 0.D0 )
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
                      DO j = 1, npw
                         tg_psi(dffts%nl (j)+ioff)=     evc(j,idx+ibnd-1)+&
                              (0.0d0,1.d0) * evc(j,idx+ibnd)
                         tg_psi(dffts%nlm(j)+ioff)=CONJG(evc(j,idx+ibnd-1) -&
                              (0.0d0,1.d0) * evc(j,idx+ibnd) )
                      END DO
                   ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
                      DO j = 1, npw
                         tg_psi(dffts%nl (j)+ioff)=       evc(j,idx+ibnd-1)
                         tg_psi(dffts%nlm(j)+ioff)=CONJG( evc(j,idx+ibnd-1) )
                      END DO
                   END IF

                   ioff = ioff + right_nnr

                END DO
                !
                CALL invfft ('tgWave', tg_psi, dffts )
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
                CALL get_rho_gamma(tg_rho, dffts%nr1x*dffts%nr2x*right_nr3, w1, w2, tg_psi)
                !
             ELSE
                !
                psic(:) = ( 0.D0, 0.D0 )
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   psic(dffts%nl(1:npw))  = evc(1:npw,ibnd) + &
                                           ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1)
                   psic(dffts%nlm(1:npw)) = CONJG( evc(1:npw,ibnd) - &
                                           ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   !
                ELSE
                   !
                   psic(dffts%nl (1:npw))  = evc(1:npw,ibnd)
                   psic(dffts%nlm(1:npw)) = CONJG( evc(1:npw,ibnd) )
                   !
                END IF
                !
                CALL invfft ('Wave', psic, dffts)
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
                CALL get_rho_gamma(rho%of_r(:,current_spin), dffts%nnr, w1, w2, psic)
                !
             END IF
             !
             IF (xclib_dft_is('meta') .OR. lxdm) THEN
                CALL using_evc(0)
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

             CALL tg_reduce_rho( rho%of_r, tg_rho, current_spin, dffts )
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL mp_sum( becsum, becp%comm )
       IF( okvan .AND. becp%comm /= mp_get_comm_null() .and. tqr ) CALL mp_sum( ebecsum, becp%comm )
       !
       IF( use_tg ) THEN
          DEALLOCATE( tg_psi )
          DEALLOCATE( tg_rho )
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for k-points version
       !
       USE mp_bands,     ONLY : me_bgrp
       USE mp,           ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines, ONLY : fftx_ntgrp, fftx_tgpe, &
                          tg_reduce_rho, tg_get_nnr, tg_get_group_nr3
       USE uspp_init,    ONLY : init_us_2
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, ioff, ioff_tg, nxyp, incr, v_siz, j, ir3
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_psi_nc(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:), tg_rho_nc(:,:)
       LOGICAL  :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp
       !
       ! chunking parameters
       INTEGER, PARAMETER :: blocksize = 256
       INTEGER :: numblock
       !
       CALL using_evc(0); CALL using_et(0)
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (xclib_dft_is('meta') .OR. lxdm) )
       !
       incr = 1
       !
       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg
          !
          IF (noncolin) THEN
             ALLOCATE( tg_psi_nc( v_siz, npol ) )
             ALLOCATE( tg_rho_nc( v_siz, nspin_mag ) )
          ELSE
             ALLOCATE( tg_psi( v_siz ) )
             ALLOCATE( tg_rho( v_siz ) )
          ENDIF
          !
          incr  = fftx_ntgrp(dffts)
          !
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF( use_tg ) THEN
            IF (noncolin) THEN
               tg_rho_nc = 0.0_DP
            ELSE
               tg_rho = 0.0_DP
            ENDIF
          ENDIF

          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          CALL start_clock( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          IF ( nks > 1 ) CALL using_evc(1)

          CALL stop_clock( 'sum_band:buffer' )
          !
          CALL start_clock( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
          CALL stop_clock( 'sum_band:init_us_2' )
          !
          ! ... for DMFT the eigenvectors are updated using v_dmft from add_dmft_occ.f90
          !
          IF ( dmft .AND. .NOT. dmft_updated) THEN
             ! 
             DO j = 1, npw
                CALL ZGEMM('C', 'N', nbnd, 1, nbnd, (1.d0,0.d0), v_dmft(:,:,ik), nbnd, evc(j,:), nbnd, (0.d0,0.d0), evc(j,:), nbnd)
             ENDDO
             !
             IF ( nks > 1 ) &
                CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
             !
          ENDIF
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
                   tg_psi_nc = ( 0.D0, 0.D0 )
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
                         DO j = 1, npw
                            tg_psi_nc( dffts%nl(igk_k(j,ik) ) + ioff, 1 ) = &
                                                       evc( j, idx+ibnd-1 )
                            tg_psi_nc( dffts%nl(igk_k(j,ik) ) + ioff, 2 ) = &
                                                       evc( j+npwx, idx+ibnd-1 )
                         END DO
                      END IF

                      ioff = ioff + right_nnr

                   END DO
                   !
                   CALL invfft ('tgWave', tg_psi_nc(:,1), dffts )
                   CALL invfft ('tgWave', tg_psi_nc(:,2), dffts)
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
                   DO ipol=1,npol
                      CALL get_rho(tg_rho_nc(:,1), dffts%nr1x * dffts%nr2x* right_nr3, w1, tg_psi_nc(:,ipol))
                   ENDDO
                   !
                   IF (domag) CALL get_rho_domag(tg_rho_nc(:,:), dffts%nr1x*dffts%nr2x*dffts%my_nr3p, w1, tg_psi_nc(:,:))
                   !
                ELSE
!
!     Noncollinear case without task groups
!
                   psic_nc = (0.D0,0.D0)
                   DO ig = 1, npw
                      psic_nc(dffts%nl(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
                      psic_nc(dffts%nl(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
                   END DO
                   CALL invfft ('Wave', psic_nc(:,1), dffts)
                   CALL invfft ('Wave', psic_nc(:,2), dffts)
                   !
                   ! increment the charge density ...
                   !
                   DO ipol=1,npol
                      CALL get_rho(rho%of_r(:,1), dffts%nnr, w1, psic_nc(:,ipol))
                   END DO
                   !
                   ! In this case, calculate also the three
                   ! components of the magnetization (stored in rho%of_r(ir,2-4))
                   !
                   IF (domag) THEN
                      CALL get_rho_domag(rho%of_r(:,:), dffts%nnr, w1, psic_nc(:,:))
                   ELSE
                      rho%of_r(:,2:4)=0.0_DP
                   END IF
                   !
                END IF
                !
             ELSE
                !
                IF( use_tg ) THEN
                   !
                   CALL tg_get_nnr( dffts, right_nnr )
                   !
                   ntgrp = fftx_ntgrp( dffts )
                   !
                   ! compute the number of chuncks
                   numblock  = (npw+blocksize-1)/blocksize
                   !
                !$omp parallel
                   CALL threaded_barrier_memset(tg_psi, 0.D0, ntgrp*right_nnr*2)
                   !
                   ! ... ntgrp ffts at the same time
                   !
                   !$omp do collapse(2)
                   DO idx = 0, MIN(ntgrp-1, ibnd_end-ibnd)
                      DO j = 1, numblock
                         tg_psi( dffts%nl(igk_k((j-1)*blocksize+1:MIN(j*blocksize, npw),ik))+right_nnr*idx ) = &
                            evc((j-1)*blocksize+1:MIN(j*blocksize, npw),idx+ibnd)
                      END DO
                   END DO
                   !$omp end do nowait
                !$omp end parallel
                   !
                   CALL invfft ('tgWave', tg_psi, dffts)
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
                   CALL get_rho(tg_rho, dffts%nr1x * dffts%nr2x * right_nr3, w1, tg_psi)
                   !
                ELSE
                   !
                !$omp parallel
                   CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
                   !
                   !$omp do
                   DO j = 1, npw
                      psic(dffts%nl(igk_k(j,ik))) = evc(j,ibnd)
                   ENDDO
                   !$omp end do nowait
                !$omp end parallel
                   !
                   CALL invfft ('Wave', psic, dffts)
                   !
                   ! ... increment the charge density ...
                   !
                   CALL get_rho(rho%of_r(:,current_spin), dffts%nnr, w1, psic)

                END IF
                !
                IF (xclib_dft_is('meta') .OR. lxdm) THEN
                   DO j=1,3
                      psic(:) = ( 0.D0, 0.D0 )
                      !
                      kplusg (1:npw) = (xk(j,ik)+g(j,igk_k(1:npw,ik))) * tpiba
                      psic(dffts%nl(igk_k(1:npw,ik)))=CMPLX(0d0,kplusg(1:npw),kind=DP) * &
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
             CALL tg_reduce_rho( rho%of_r, tg_rho_nc, tg_rho, current_spin, noncolin, domag, dffts )
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF( okvan .AND. becp%comm /= mp_get_comm_null() ) CALL mp_sum( becsum, becp%comm )
       IF( okvan .AND. becp%comm /= mp_get_comm_null() .and. tqr ) CALL mp_sum( ebecsum, becp%comm )
       !
       IF( use_tg ) THEN
          IF (noncolin) THEN
             DEALLOCATE( tg_psi_nc )
             DEALLOCATE( tg_rho_nc )
          ELSE
             DEALLOCATE( tg_psi )
             DEALLOCATE( tg_rho )
          END IF
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
     !
     SUBROUTINE get_rho(rho_loc, nrxxs_loc, w1_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc(nrxxs_loc)

        INTEGER :: ir

        !$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * ( DBLE( psic_loc(ir) )**2 + &
                                   AIMAG( psic_loc(ir) )**2 )
           !
        END DO
        !$omp end parallel do

     END SUBROUTINE get_rho

     SUBROUTINE get_rho_gamma(rho_loc, nrxxs_loc, w1_loc, w2_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP) :: psic_loc(nrxxs_loc)

        INTEGER :: ir

        !$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * DBLE( psic_loc(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc(ir) )**2
           !
        END DO
        !$omp end parallel do

     END SUBROUTINE get_rho_gamma


     SUBROUTINE get_rho_domag(rho_loc, nrxxs_loc, w1_loc, psic_loc)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc(:, :)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc(:, :)

        INTEGER :: ir

        !$omp parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir,2) = rho_loc(ir,2) + w1_loc*2.D0* &
                          (DBLE(psic_loc(ir,1))* DBLE(psic_loc(ir,2)) + &
                          AIMAG(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)))
 
           rho_loc(ir,3) = rho_loc(ir,3) + w1_loc*2.D0* &
                          (DBLE(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)) - &
                           DBLE(psic_loc(ir,2))*AIMAG(psic_loc(ir,1)))

           rho_loc(ir,4) = rho_loc(ir,4) + w1_loc* &
                          (DBLE(psic_loc(ir,1))**2+AIMAG(psic_loc(ir,1))**2 &
                          -DBLE(psic_loc(ir,2))**2-AIMAG(psic_loc(ir,2))**2)
           !
        END DO
        !$omp end parallel do

     END SUBROUTINE get_rho_domag

END SUBROUTINE sum_band

!----------------------------------------------------------------------------
SUBROUTINE sum_bec ( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
  !----------------------------------------------------------------------------
  !! This routine computes the sum over bands:
  !
  !! \[ \sum_i \langle\psi_i|\beta_l\rangle w_i \langle\beta_m|\psi_i\rangle \]
  !
  !! for point "ik" and, for LSDA, spin "current_spin".  
  !! Calls calbec to compute \(\text{"becp"}=\langle \beta_m|\psi_i \rangle\).  
  !! Output is accumulated (unsymmetrized) into "becsum", module "uspp".
  !
  !! Routine used in sum_band (if okvan) and in compute_becsum, called by hinit1 (if okpaw).
  !
  USE kinds,              ONLY : DP
  USE becmod,             ONLY : becp, calbec, allocate_bec_type
  USE control_flags,      ONLY : gamma_only, tqr
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp,               ONLY : nkb, vkb, becsum, ebecsum, ofsbeta
  USE uspp_param,         ONLY : upf, nh, nhm
  USE wvfct,              ONLY : nbnd, wg, et, current_k
  USE klist,              ONLY : ngk
  USE noncollin_module,   ONLY : noncolin, npol
  USE wavefunctions,      ONLY : evc
  USE realus,             ONLY : real_space, &
                                 invfft_orbital_gamma, calbec_rs_gamma, &
                                 invfft_orbital_k, calbec_rs_k
  USE us_exx,             ONLY : store_becxx0
  USE mp_bands,           ONLY : nbgrp,inter_bgrp_comm
  USE mp,                 ONLY : mp_sum
  USE wavefunctions_gpum, ONLY : using_evc
  USE wvfct_gpum,         ONLY : using_et
  USE becmod_subs_gpum,   ONLY : using_becp_auto
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the point the sum is computed for
  INTEGER, INTENT(IN) :: current_spin
  !! the input spin
  INTEGER, INTENT(IN) :: ibnd_start
  !! first band of the parallel block
  INTEGER, INTENT(IN) :: ibnd_end
  !! last band of the parallel block
  INTEGER, INTENT(IN) :: this_bgrp_nbnd
  !! number of bands of the parallel block
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1(:,:), auxk2(:,:), aux_nc(:,:)
  REAL(DP), ALLOCATABLE :: auxg(:,:), aux_gk(:,:), aux_egk(:,:)
  INTEGER :: ibnd, ibnd_loc, nbnd_loc, kbnd  ! counters on bands
  INTEGER :: npw, ikb, jkb, ih, jh, ijh, na, np, is, js
  ! counters on beta functions, atoms, atom types, spin
  !
  CALL using_evc(0) ! calbec->in ; invfft_orbital_gamma|k -> in
  CALL using_et(0)
  CALL using_becp_auto(2)
  !
  CALL start_clock( 'sum_band:calbec' )
  npw = ngk(ik)
  IF ( .NOT. real_space ) THEN
     ! calbec computes becp = <vkb_i|psi_j>
     CALL calbec( npw, vkb, evc(:,ibnd_start:ibnd_end), becp )
  ELSE
     if (gamma_only) then
        do kbnd = 1, this_bgrp_nbnd, 2 !  ibnd_start, ibnd_end, 2
           ibnd = ibnd_start + kbnd - 1
           call invfft_orbital_gamma(evc,ibnd,ibnd_end) 
           call calbec_rs_gamma(kbnd,this_bgrp_nbnd,becp%r)
        enddo
     else
        current_k = ik
        becp%k = (0.d0,0.d0)
        do kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
           ibnd = ibnd_start + kbnd - 1
           call invfft_orbital_k(evc,ibnd,ibnd_end) 
           call calbec_rs_k(kbnd,this_bgrp_nbnd)
        enddo
     endif
  ENDIF
  CALL stop_clock( 'sum_band:calbec' )
  !
  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  CALL store_becxx0(ik, becp)
  !
  CALL start_clock( 'sum_band:becsum' )
  !
  DO np = 1, ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ! allocate work space used to perform GEMM operations
        !
        IF ( gamma_only ) THEN
           nbnd_loc = becp%nbnd_loc
           ALLOCATE( auxg( nbnd_loc, nh(np) ) )
        ELSE
           ALLOCATE( auxk1( ibnd_start:ibnd_end, nh(np)*npol ), &
                     auxk2( ibnd_start:ibnd_end, nh(np)*npol ) )
        END IF
        IF ( noncolin ) THEN
           ALLOCATE ( aux_nc( nh(np)*npol,nh(np)*npol ) ) 
        ELSE
           ALLOCATE ( aux_gk( nh(np),nh(np) ) ) 
           if (tqr) ALLOCATE ( aux_egk( nh(np),nh(np) ) ) 
        END IF
        !
        !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
        !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
        !
        DO na = 1, nat
           !
           IF (ityp(na)==np) THEN
              !
              ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
              ! copy into aux1, aux2 the needed data to perform a GEMM
              !
              IF ( noncolin ) THEN
                 !
                 !$omp parallel do default(shared), private(is,ih,ikb,ibnd)
                 DO is = 1, npol
                    DO ih = 1, nh(np)
                       ikb = ofsbeta(na) + ih
                       DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                          ibnd = ibnd_start + kbnd - 1
                          auxk1(ibnd,ih+(is-1)*nh(np))= becp%nc(ikb,is,kbnd)
                          auxk2(ibnd,ih+(is-1)*nh(np))= wg(ibnd,ik) * &
                                                        becp%nc(ikb,is,kbnd)
                       END DO
                    END DO
                 END DO
                 !$omp end parallel do
                 !
                 CALL ZGEMM ( 'C', 'N', npol*nh(np), npol*nh(np), this_bgrp_nbnd, &
                      (1.0_dp,0.0_dp), auxk1, this_bgrp_nbnd, auxk2, this_bgrp_nbnd, &
                      (0.0_dp,0.0_dp), aux_nc, npol*nh(np) )
                 !
              ELSE IF ( gamma_only ) THEN
                 !
                 !$omp parallel do default(shared), private(ih,ikb,ibnd,ibnd_loc)
                 DO ih = 1, nh(np)
                    ikb = ofsbeta(na) + ih
                    DO ibnd_loc = 1, nbnd_loc
                       ibnd = (ibnd_start - 1) + ibnd_loc + becp%ibnd_begin - 1
                       auxg(ibnd_loc,ih)= wg(ibnd,ik)*becp%r(ikb,ibnd_loc) 
                    END DO
                 END DO
                 !$omp end parallel do
                 CALL DGEMM ( 'N', 'N', nh(np), nh(np), nbnd_loc, &
                      1.0_dp, becp%r(ofsbeta(na)+1,1), nkb,    &
                      auxg, nbnd_loc, 0.0_dp, aux_gk, nh(np) )
               if (tqr) then
                 !$omp parallel do default(shared), private(ih,ikb,ibnd,ibnd_loc)
                 DO ih = 1, nh(np)
                    ikb = ofsbeta(na) + ih
                    DO ibnd_loc = 1, nbnd_loc
                       ibnd = (ibnd_start - 1) + ibnd_loc + becp%ibnd_begin - 1
                       auxg(ibnd_loc,ih) = et(ibnd,ik) * auxg(ibnd_loc,ih)
                    END DO
                 END DO
                 !$omp end parallel do
                 CALL DGEMM ( 'N', 'N', nh(np), nh(np), nbnd_loc, &
                      1.0_dp, becp%r(ofsbeta(na)+1,1), nkb,    &
                      auxg, nbnd_loc, 0.0_dp, aux_egk, nh(np) )
               end if
                 !
              ELSE
                 !
                 !$omp parallel do default(shared), private(ih,ikb,ibnd)
                 DO ih = 1, nh(np)
                    ikb = ofsbeta(na) + ih
                    DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                       ibnd = ibnd_start + kbnd - 1
                       auxk1(ibnd,ih) = becp%k(ikb,kbnd) 
                       auxk2(ibnd,ih) = wg(ibnd,ik)*becp%k(ikb,kbnd)
                    END DO
                 END DO
                 !$omp end parallel do
                 !
                 ! only the real part is computed
                 !
                 CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_gk, nh(np) )
                 !
                 if (tqr) then
                   !$omp parallel do default(shared), private(ih,ikb,ibnd)
                   DO ih = 1, nh(np)
                      ikb = ofsbeta(na) + ih
                      DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                         ibnd = ibnd_start + kbnd - 1
                         auxk2(ibnd,ih) = et(ibnd,ik)*auxk2(ibnd,ih)
                      END DO
                   END DO
                   !$omp end parallel do
                   !
                   CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
                        1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
                        0.0_dp, aux_egk, nh(np) )
                 end if
                 !
              END IF
              !
              ! copy output from GEMM into desired format
              !
              IF (noncolin .AND. .NOT. upf(np)%has_so) THEN
                 CALL add_becsum_nc (na, np, aux_nc, becsum )
              ELSE IF (noncolin .AND. upf(np)%has_so) THEN
                 CALL add_becsum_so (na, np, aux_nc,becsum )
              ELSE
                 ijh = 0
                 DO ih = 1, nh(np)
                    DO jh = ih, nh(np)
                       ijh = ijh + 1
                       !
                       ! nondiagonal terms summed and collapsed into a
                       ! single index (matrix is symmetric wrt (ih,jh))
                       !
                       IF ( jh == ih ) THEN
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk (ih,jh)
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                               ebecsum(ijh,na,current_spin) + aux_egk (ih,jh)
                       ELSE
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk(ih,jh)*2.0_dp
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                               ebecsum(ijh,na,current_spin) + aux_egk(ih,jh)*2.0_dp
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
           DEALLOCATE ( aux_nc )
        ELSE
           DEALLOCATE ( aux_gk  ) 
           if (tqr) DEALLOCATE ( aux_egk  ) 
        END IF
        IF ( gamma_only ) THEN
           DEALLOCATE( auxg )
        ELSE
           DEALLOCATE( auxk2, auxk1 )
        END IF
        !
     END IF
     !
  END DO
  !
  CALL stop_clock( 'sum_band:becsum' )
  !
END SUBROUTINE sum_bec
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_nc ( na, np, becsum_nc, becsum )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the
  !! Pauli matrices, saves it in \(\text{becsum}\) for the calculation of 
  !! augmentation charge and magnetization.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, ijh
  REAL(dp) :: fac
  !
  ijh=0
  DO ih = 1, nh(np)
     DO jh = ih, nh(np)
        ijh=ijh+1
        IF ( ih == jh ) THEN
           fac = 1.0_dp
        ELSE
           fac = 2.0_dp
        END IF
        becsum(ijh,na,1)= becsum(ijh,na,1) + fac * &
                DBLE( becsum_nc(ih,1,jh,1) + becsum_nc(ih,2,jh,2) )
        IF (domag) THEN
           becsum(ijh,na,2)= becsum(ijh,na,2) + fac *  &
                DBLE( becsum_nc(ih,1,jh,2) + becsum_nc(ih,2,jh,1) )
           becsum(ijh,na,3)= becsum(ijh,na,3) + fac * DBLE( (0.d0,-1.d0)* &
               (becsum_nc(ih,1,jh,2) - becsum_nc(ih,2,jh,1)) )
           becsum(ijh,na,4)= becsum(ijh,na,4) + fac * &
                DBLE( becsum_nc(ih,1,jh,1) - becsum_nc(ih,2,jh,2) )
        END IF
     END DO
  END DO
  
END SUBROUTINE add_becsum_nc
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_so( na, np, becsum_nc, becsum )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the Pauli
  !! matrices, rotates it as appropriate for the spin-orbit case, saves it in 
  !! \(\text{becsum}\) for the calculation of augmentation charge and magnetization.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE uspp,                 ONLY : ijtoh, nhtol, nhtoj, indv
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
  USE upf_spinorb,          ONLY : fcoef
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na
  !! atom index
  INTEGER, INTENT(IN) :: np
  !! atomic type
  COMPLEX(DP), INTENT(IN) :: becsum_nc(nh(np),npol,nh(np),npol)
  !! becsum - noncolin
  REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
  !! sum over bands
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, lh, kh, ijh, is1, is2
  COMPLEX(DP) :: fac
  !
  DO ih = 1, nh(np)
     DO jh = 1, nh(np)
        ijh=ijtoh(ih,jh,np)
        DO kh = 1, nh(np)
           IF (same_lj(kh,ih,np)) THEN
              DO lh=1,nh(np)
                 IF (same_lj(lh,jh,np)) THEN
                    DO is1=1,npol
                       DO is2=1,npol
                          fac=becsum_nc(kh,is1,lh,is2)
                          becsum(ijh,na,1)=becsum(ijh,na,1) + fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) + &
                               fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                          IF (domag) THEN
                            becsum(ijh,na,2)=becsum(ijh,na,2)+fac * &
                                (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) +&
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                            becsum(ijh,na,3)=becsum(ijh,na,3)+fac*(0.d0,-1.d0)*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                           becsum(ijh,na,4)=becsum(ijh,na,4) + fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END DO
!
CONTAINS
   LOGICAL FUNCTION same_lj(ih,jh,np)
   INTEGER :: ih, jh, np
   !
   same_lj = ((nhtol(ih,np)==nhtol(jh,np)).AND. &
             (ABS(nhtoj(ih,np)-nhtoj(jh,np))<1.d8).AND. &
             (indv(ih,np)==indv(jh,np)) )
   !
   END FUNCTION same_lj

END SUBROUTINE add_becsum_so
