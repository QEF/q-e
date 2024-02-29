!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band_gpu()
  !----------------------------------------------------------------------------
  !! Calculates the symmetrized charge density and related quantities.  
  !! Also computes the occupations and the sum of occupied eigenvalues.
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
  USE fft_interfaces,       ONLY : invfft
  USE fft_rho,              ONLY : rho_g2r, rho_r2g
  USE fft_wave,             ONLY : wave_g2r, tgwave_g2r
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho, rhoz_or_updw
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, nhtol, nhtoj, indv, okvan, &
                                   becsum_d, ebecsum_d
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions,        ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum
  USE xc_lib,               ONLY : xclib_dft_is
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                                   becp
  USE gcscf_module,         ONLY : lgcscf, gcscf_calc_nelec
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
  REAL(DP), ALLOCATABLE :: kplusg(:)
  COMPLEX(DP), ALLOCATABLE :: psicd(:), kplusg_evc(:,:)
  !
  !
  CALL start_clock_gpu( 'sum_band' )
  !
  IF ( nhm > 0 ) THEN
     becsum(:,:,:) = 0.D0
     IF (tqr) ebecsum(:,:,:) = 0.D0
     becsum_d(:,:,:) = 0.D0
     IF (tqr) ebecsum_d(:,:,:) = 0.D0
  ENDIF
  rho%of_r(:,:) = 0.D0
  rho%of_g(:,:) = (0.D0, 0.D0)
  IF ( xclib_dft_is('meta') .OR. lxdm ) THEN
     rho%kin_r(:,:) = 0.D0
     rho%kin_g(:,:) = (0.D0, 0.D0)
  ENDIF
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL start_clock_gpu( 'sum_band:weights' )
  CALL weights()
  CALL stop_clock_gpu( 'sum_band:weights' )
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
  ENDIF
  !
  ! ... Needed for DFT+Hubbard: compute occupations of Hubbard states
  !
  IF (lda_plus_u) THEN
    IF (lda_plus_u_kind==0) THEN
       !
       IF (noncolin) THEN
          CALL new_ns_nc(rho%ns_nc)
       ELSE
          CALL new_ns(rho%ns)
       ENDIF
       !
       DO nt = 1, ntyp
          IF (is_hubbard_back(nt)) CALL new_nsb( rho%nsb )
       ENDDO
       !
    ELSEIF (lda_plus_u_kind==1) THEN
       !
       IF (noncolin) THEN
          CALL new_ns_nc( rho%ns_nc )
       ELSE
          CALL new_ns( rho%ns )
       ENDIF
       !
    ELSEIF (lda_plus_u_kind==2) THEN 
       !
       IF (noncolin) THEN
          CALL new_nsg_nc()
       ELSE
          CALL new_nsg()
       ENDIF
       !
    ENDIF
  ENDIF
  !
  ! ... for band parallelization: set band computed by this processor
  !
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type_acc( nkb, this_bgrp_nbnd, becp, intra_bgrp_comm )
  IF (xclib_dft_is('meta') .OR. lxdm) ALLOCATE( kplusg(npwx), kplusg_evc(npwx,2) )
  !
  ! ... specialized routines are called to sum at Gamma or for each k point 
  ! ... the contribution of the wavefunctions to the charge
  ! ... The band energy contribution eband is computed together with the charge
  !
  eband = 0.D0
  !
  CALL start_clock_gpu( 'sum_band:loop' )
  !
  !$acc enter data create(evc)
  IF (noncolin) THEN
    !$acc enter data create(psic_nc)
  ELSE
    !$acc enter data create(psic)
  ENDIF
  !
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma_gpu()
     !
  ELSE
     !
     CALL sum_band_k_gpu()
     !
  ENDIF
  !
  !$acc exit data delete(evc)
  IF (noncolin) THEN
    !$acc exit data delete(psic_nc)
  ELSE
    !$acc exit data delete(psic)
  ENDIF
  !
  CALL stop_clock_gpu( 'sum_band:loop' )
  CALL mp_sum( eband, inter_pool_comm )
  CALL mp_sum( eband, inter_bgrp_comm )
  !
  IF (xclib_dft_is('meta') .OR. lxdm) THEN
    DEALLOCATE( kplusg, kplusg_evc )
  ENDIF
  IF ( okvan ) CALL deallocate_bec_type_acc ( becp )
  !
  ! ... sum charge density over pools (distributed k-points) and bands
  !
  CALL mp_sum( rho%of_r, inter_pool_comm )
  CALL mp_sum( rho%of_r, inter_bgrp_comm )
  IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
  !
  ! ... bring the unsymmetrized rho(r) to G-space (use psic as work array)
  !
  CALL rho_r2g( dffts, rho%of_r, rho%of_g )
  !
  IF( okvan )  THEN
     !
     ! ... becsum is summed over bands (if bgrp_parallelization is done)
     ! ... and over k-points (but it is not symmetrized)
     !
     ! ... use host copy to do the comunication. This avoids going back an forth GPU data
     ! ... becsum=becsum_d     not needed 
     ! ... since becsum is already uptodate, see sum_band*gpu
     !
     CALL mp_sum(becsum, inter_bgrp_comm )
     CALL mp_sum(becsum, inter_pool_comm )
     becsum_d=becsum
     !
     ! ... same for ebecsum, a correction to becsum (?) in real space
     !
     IF (tqr) THEN
        !ebecsum=ebecsum_d   not needed as above
        CALL mp_sum(ebecsum, inter_pool_comm )
        CALL mp_sum(ebecsum, inter_bgrp_comm )
        ebecsum_d=ebecsum
     ENDIF
     !
     ! ... PAW: symmetrize becsum and store it
     ! ... FIXME: the same should be done for USPP as well
     !
     IF ( okpaw ) THEN
        rho%bec(:,:,:) = becsum(:,:,:)
        CALL PAW_symmetrize( rho%bec )
     ENDIF
     !
     ! ... Here we add the (unsymmetrized) Ultrasoft contribution to the charge
     !
     CALL addusdens( rho%of_g )
     !
  ENDIF
  !
  ! ... symmetrize rho(G) 
  !
  CALL start_clock_gpu( 'sum_band:sym_rho' )
  CALL sym_rho( nspin_mag, rho%of_g )
  !
  ! ... synchronize rho%of_r to the calculated rho%of_g (use psic as work array)
  !
  CALL rho_g2r( dfftp, rho%of_g, rho%of_r )
  !
  ! ... rho_kin(r): sum over bands, k-points, bring to G-space, symmetrize,
  ! ... synchronize with rho_kin(G)
  !
  IF ( xclib_dft_is('meta') .OR. lxdm) THEN
     !
     CALL mp_sum( rho%kin_r, inter_pool_comm )
     CALL mp_sum( rho%kin_r, inter_bgrp_comm )
     !
     CALL rho_r2g( dffts, rho%kin_r, rho%kin_g )
     !
     IF (.NOT. gamma_only) CALL sym_rho( nspin, rho%kin_g )
     !
     CALL rho_g2r( dfftp, rho%kin_g, rho%kin_r )
     !
  ENDIF
  CALL stop_clock_gpu( 'sum_band:sym_rho' )
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
     SUBROUTINE sum_band_gamma_gpu()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for gamma version.
       !
       USE becmod,                 ONLY : becp
       USE mp_bands,               ONLY : me_bgrp
       USE mp,                     ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       USE uspp_init,              ONLY : init_us_2
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, ioff, ioff_tg, nxyp, incr, v_siz, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:)
       REAL(DP),    ALLOCATABLE :: tg_rho_d(:), tg_rho_h(:)
       REAL(DP),    ALLOCATABLE :: rho_d(:,:)
       LOGICAL :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp, ierr, ebnd, i, brange
       REAL(DP) :: kplusgi
#if defined(__CUDA)
       attributes(device) :: tg_rho_d, rho_d
       attributes(pinned) :: tg_rho_h
#endif
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (xclib_dft_is('meta') .OR. lxdm) )
       !
       incr = 2
       !
       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg 
          !
          ALLOCATE( tg_psi( v_siz ) )
          !$acc enter data create(tg_psi)
          ALLOCATE( tg_rho_d( v_siz ) )
          ALLOCATE( tg_rho_h( v_siz ) )
          !
          incr = 2 * fftx_ntgrp(dffts)
          !
       ELSE
          ALLOCATE( rho_d, MOLD=rho%of_r ) ! OPTIMIZE HERE, use buffers (and batched FFT)
          rho_d = 0.0_DP
       ENDIF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( use_tg ) tg_rho_d = 0.0_DP
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          CALL start_clock_gpu( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !$acc update device(evc)
          !
          CALL stop_clock_gpu( 'sum_band:buffer' )
          !
          CALL start_clock_gpu( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
          !
          CALL stop_clock_gpu( 'sum_band:init_us_2' )
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
          ENDDO
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             IF( use_tg ) THEN
                !
                CALL tgwave_g2r( evc(1:npw,ibnd:ibnd_end), tg_psi, dffts, npw )
                !
                ! ... Now the first proc of the group holds the first two bands
                ! ... of the 2*ntgrp bands that we are processing at the same time,
                ! ... the second proc. holds the third and fourth band and so on.
                !
                ! ... Compute the proper factor for each band
                !
                idx = fftx_tgpe(dffts) + 1
                !
                ! ... Remember two bands are packed in a single array :
                !    - proc 0 has bands ibnd   and ibnd+1
                !    - proc 1 has bands ibnd+2 and ibnd+3
                !    - ....
                !
                idx = 2*idx - 1
                !
                IF( idx + ibnd - 1 < ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = wg( idx + ibnd    , ik) / omega
                ELSEIF( idx + ibnd - 1 == ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = w1
                ELSE
                   w1 = 0.0d0
                   w2 = w1
                ENDIF
                !
                CALL tg_get_group_nr3( dffts, right_nr3 )
                !
                !$acc host_data use_device(tg_psi)
                CALL get_rho_gamma_gpu( tg_rho_d, dffts%nr1x*dffts%nr2x*right_nr3, &
                                        w1, w2, tg_psi )
                !$acc end host_data
                !
             ELSE
                !
                ebnd = ibnd
                IF ( ibnd < ibnd_end ) ebnd = ebnd + 1
                !
                CALL wave_g2r( evc(1:npw,ibnd:ebnd), psic, dffts )
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
                ENDIF
                !
                !$acc host_data use_device(psic)
                CALL get_rho_gamma_gpu( rho_d(:,current_spin), dffts%nnr, w1, w2, psic )
                !$acc end host_data
                !
             ENDIF
             !
             IF (xclib_dft_is('meta') .OR. lxdm) THEN
                !
                DO j = 1, 3
                   DO i = 1, npw
                     kplusgi = (xk(j,ik)+g(j,i)) * tpiba
                     kplusg_evc(i,1) = CMPLX(0._DP,kplusgi,KIND=DP) * evc(i,ibnd)
                     IF ( ibnd < ibnd_end ) kplusg_evc(i,2) = CMPLX(0._DP,kplusgi,KIND=DP) * evc(i,ibnd+1)
                   ENDDO
                   !
                   ebnd = ibnd
                   IF ( ibnd < ibnd_end ) ebnd = ebnd + 1
                   brange = ebnd-ibnd+1
                   !
                   CALL wave_g2r( kplusg_evc(1:npw,1:brange), psic, dffts )
                   !$acc update self(psic)
                   !
                   ! ... increment the kinetic energy density ...
                   !
                   DO ir = 1, dffts%nnr
                      rho%kin_r(ir,current_spin) = rho%kin_r(ir,current_spin) + &
                                                    w1 *  DBLE( psic(ir) )**2 + &
                                                    w2 * AIMAG( psic(ir) )**2
                   ENDDO
                ENDDO
                !
             ENDIF
             !
          ENDDO
          !
          IF( use_tg ) THEN
             tg_rho_h = tg_rho_d
             CALL tg_reduce_rho( rho%of_r, tg_rho_h, current_spin, dffts )
          ENDIF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec_gpu( ik, current_spin, ibnd_start, ibnd_end, &
                                         this_bgrp_nbnd )
          !
       ENDDO k_loop
       !
       IF( .NOT. use_tg ) THEN
          rho%of_r = rho_d
       ENDIF
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. nhm>0) THEN
          !becsum=becsum_d      not needed, since already updated in sum_bec_gpu
          CALL mp_sum( becsum, becp%comm )
          becsum_d=becsum
       ENDIF
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. tqr .AND. nhm>0) THEN
          !ebecsum=ebecsum_d    as above
          CALL mp_sum( ebecsum, becp%comm )
          ebecsum_d=ebecsum
       ENDIF
       !
       IF( use_tg ) THEN
          !$acc exit data delete(tg_psi)
          DEALLOCATE( tg_psi )
          DEALLOCATE( tg_rho_d )
          DEALLOCATE( tg_rho_h )
       ELSE
          DEALLOCATE(rho_d)
       ENDIF
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma_gpu
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k_gpu()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for k-points version
       !
       USE mp_bands,               ONLY : me_bgrp
       USE mp,                     ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       USE uspp_init,              ONLY : init_us_2
       USE control_flags,          ONLY : many_fft
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, ioff, ioff_tg, nxyp, incr, v_siz
       COMPLEX(DP), ALLOCATABLE :: psicd(:)
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_psi_nc(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho_d(:), tg_rho_nc_d(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho_h(:), tg_rho_nc_h(:,:)
       REAL(DP),    ALLOCATABLE :: rho_d(:,:)
       LOGICAL  :: use_tg
       INTEGER :: nnr, right_nnr, right_nr3, right_inc, ntgrp, ierr
       INTEGER :: i, j, group_size, hm_vec(3)
       REAL(DP) :: kplusgi
       !
#if defined(__CUDA)
       attributes(device) :: tg_rho_d, tg_rho_nc_d
       attributes(device) :: rho_d
       attributes(pinned) :: tg_rho_h, tg_rho_nc_h
#endif
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (xclib_dft_is('meta') .OR. lxdm) )
       !
       incr = 1
       nnr  = dffts%nnr
       !
       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg
          !
          IF (noncolin) THEN
             ALLOCATE( tg_psi_nc( v_siz, npol ) )
             !$acc enter data create(tg_psi_nc)
             ALLOCATE( tg_rho_nc_d( v_siz, nspin_mag ) )
             ALLOCATE( tg_rho_nc_h( v_siz, nspin_mag ) )
          ELSE
             ALLOCATE( tg_psi( v_siz ) )
             !$acc enter data create(tg_psi)
             ALLOCATE( tg_rho_d( v_siz ) )
             ALLOCATE( tg_rho_h( v_siz ) )
          ENDIF
          !
          incr = fftx_ntgrp(dffts)
          !
       ELSE
          ALLOCATE( rho_d, MOLD=rho%of_r ) ! OPTIMIZE HERE, use buffers!
          IF (noncolin .OR. (xclib_dft_is('meta') .OR. lxdm)) THEN
            ALLOCATE( psicd(dffts%nnr) )
            incr = 1
          ELSE
            ALLOCATE( psicd(dffts%nnr*many_fft) )
            incr = many_fft
          ENDIF
          !$acc enter data create(psicd)
          ! ... This is used as reduction variable on the device
          rho_d = 0.0_DP
       ENDIF
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
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          CALL start_clock_gpu( 'sum_band:buffer' )
          IF ( nks > 1 ) THEN
             CALL get_buffer( evc, nwordwfc, iunwfc, ik )
          ENDIF
          !$acc update device(evc)
          CALL stop_clock_gpu( 'sum_band:buffer' )
          !
          CALL start_clock_gpu( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
          !
          CALL stop_clock_gpu( 'sum_band:init_us_2' )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             !IF( use_tg ) THEN
             DO idx = 1, incr
                IF( idx+ibnd-1 <= ibnd_end ) eband = eband + et(idx+ibnd-1,ik) * &
                                                             wg(idx+ibnd-1,ik)
             ENDDO
             !ELSE
             !   eband = eband + et( ibnd, ik ) * wg( ibnd, ik )
             !END IF
             !
             ! ... the sum of eband and demet is the integral for e < ef of
             ! ... e n(e) which reduces for degauss=0 to the sum of the
             ! ... eigenvalues
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
                IF( use_tg ) THEN
                   !
                   CALL tgwave_g2r( evc(1:npw,ibnd:ibnd_end), tg_psi_nc(:,1), dffts, &
                                    npw, igk_k(:,ik) )
                   CALL tgwave_g2r( evc(npwx+1:npwx+npw,ibnd:ibnd_end), tg_psi_nc(:,2), &
                                    dffts, npw, igk_k(:,ik) )
                   !
                   ! ... Now the first proc of the group holds the first band
                   ! ... of the ntgrp bands that we are processing at the same time,
                   ! ... the second proc. holds the second and so on
                   !
                   ! ... Compute the proper factor for each band
                   !
                   idx = fftx_tgpe(dffts) + 1 
                   !
                   ! ... Remember
                   ! ... proc 0 has bands ibnd
                   ! ... proc 1 has bands ibnd+1
                   ! ...    ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   ENDIF
                   !
                   CALL tg_get_group_nr3( dffts, right_nr3 )
                   !
                   ! OPTIMIZE HERE : this is a sum of all densities in first spin channel
                   !
                   !$acc host_data use_device(tg_psi_nc)
                   DO ipol = 1, npol
                      CALL get_rho_gpu( tg_rho_nc_d(:,1), dffts%nr1x*dffts%nr2x* &
                                        right_nr3, w1, tg_psi_nc(:,ipol) )
                   ENDDO
                   !
                   IF (domag) CALL get_rho_domag_gpu( tg_rho_nc_d(:,:), dffts%nr1x* &
                                      dffts%nr2x*dffts%my_nr3p, w1, tg_psi_nc(:,:) )
                   !$acc end host_data
                   !
                ELSE
                   !
                   ! ... Noncollinear case without task groups
                   !
                   CALL wave_g2r( evc(1:npw,ibnd:ibnd), psic_nc(:,1), &
                                  dffts, igk=igk_k(:,ik) )
                   CALL wave_g2r( evc(npwx+1:npwx+npw,ibnd:ibnd), psic_nc(:,2), &
                                  dffts, igk=igk_k(:,ik) )
                   !
                   ! ... Increment the charge density
                   !
                   DO ipol = 1, npol
                      !$acc host_data use_device(psic_nc)
                      CALL get_rho_gpu( rho_d(:,1), dffts%nnr, w1, psic_nc(:,ipol) )
                      !$acc end host_data
                   ENDDO
                   !
                   ! ... In this case, calculate also the three
                   ! ... components of the magnetization (stored in rho%of_r(ir,2-4))
                   !
                   IF (domag) THEN
                      !$acc host_data use_device(psic_nc)
                      CALL get_rho_domag_gpu( rho_d(1:,1:), dffts%nnr, w1, psic_nc(1:,1:) )
                      !$acc end host_data
                   ELSE
                      rho_d(:,2:4) = 0.0_DP  ! OPTIMIZE HERE: this memset can be avoided
                   ENDIF
                   !
                ENDIF
                !
             ELSE
                !
                IF( use_tg ) THEN
                   !
                   CALL tgwave_g2r( evc(1:npw,ibnd:ibnd_end), tg_psi, dffts, npw, &
                                    igk_k(:,ik) )
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
                   IF( idx+ibnd-1 <= ibnd_end ) THEN
                      w1 = wg(idx+ibnd-1,ik) / omega
                   ELSE
                      w1 = 0.0d0
                   ENDIF
                   !
                   CALL tg_get_group_nr3( dffts, right_nr3 )
                   !
                   !$acc host_data use_device(tg_psi)
                   CALL get_rho_gpu( tg_rho_d, dffts%nr1x*dffts%nr2x*right_nr3, w1, tg_psi )
                   !$acc end host_data
                   !
                ELSEIF (many_fft > 1 .AND. (.NOT. (xclib_dft_is('meta') .OR. lxdm))) THEN
                   !
                   group_size = MIN(many_fft,ibnd_end-(ibnd-1))
                   hm_vec(1)=group_size ; hm_vec(2)=npw ; hm_vec(3)=group_size
                   !
                   CALL wave_g2r( evc(:,ibnd:ibnd+group_size-1), psicd, &
                                  dffts, igk=igk_k(:,ik), howmany_set=hm_vec )
                   !
                   ! ... increment the charge density ...
                   !
                   DO i = 0, group_size-1
                     w1 = wg(ibnd+i,ik) / omega
                     !$acc host_data use_device(psicd)
                     CALL get_rho_gpu( rho_d(:,current_spin), nnr, w1, psicd(i*nnr+1:) )
                     !$acc end host_data
                   ENDDO
                   !
                ELSE
                   !
                   CALL wave_g2r( evc(1:npw,ibnd:ibnd), psicd, dffts, igk=igk_k(:,ik) )
                   !
                   ! ... increment the charge density ...
                   !
                   !$acc host_data use_device(psicd)
                   CALL get_rho_gpu( rho_d(:,current_spin), dffts%nnr, w1, psicd )
                   !$acc end host_data
                   !
                ENDIF
                !
                IF (xclib_dft_is('meta') .OR. lxdm) THEN
                   DO j=1,3
                      DO i = 1, npw
                        kplusgi = (xk(j,ik)+g(j,igk_k(i,ik))) * tpiba
                        kplusg_evc(i,1) = CMPLX(0._DP,kplusgi,KIND=DP) * evc(i,ibnd)
                      ENDDO
                      !
                      CALL wave_g2r( kplusg_evc(1:npw,1:1), psic, dffts, igk=igk_k(:,ik) )
                      !$acc update self(psic)
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      CALL get_rho( rho%kin_r(:,current_spin), dffts%nnr, w1, psic )
                   ENDDO
                ENDIF
                !
             ENDIF
             !
          ENDDO
          !
          IF ( use_tg ) THEN
             !
             ! ... reduce the charge across task group
             !
             IF (noncolin)       tg_rho_nc_h = tg_rho_nc_d
             IF (.NOT. noncolin) tg_rho_h    = tg_rho_d
             CALL tg_reduce_rho( rho%of_r, tg_rho_nc_h, tg_rho_h, current_spin, &
                                 noncolin, domag, dffts )
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
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. nhm>0 ) THEN 
          !becsum=becsum_d     not needed, since already updated in sum_bec_gpu
          CALL mp_sum( becsum, becp%comm )
          becsum_d=becsum
       ENDIF
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. tqr .AND. nhm> 0) THEN
          !ebecsum=ebecsum_d    not needed as above
          CALL mp_sum( ebecsum, becp%comm )
          ebecsum_d=ebecsum
       ENDIF
       !
       IF( use_tg ) THEN
          IF (noncolin) THEN
             !$acc exit data delete(tg_psi_nc)
             DEALLOCATE( tg_psi_nc )
             DEALLOCATE( tg_rho_nc_d )
             DEALLOCATE( tg_rho_nc_h )
          ELSE
             !$acc exit data delete(tg_psi)
             DEALLOCATE( tg_psi )
             DEALLOCATE( tg_rho_d )
             DEALLOCATE( tg_rho_h )
          END IF
       ELSE
          DEALLOCATE( rho_d ) ! OPTIMIZE HERE, use buffers!
          !$acc exit data delete(psicd)
          DEALLOCATE( psicd )
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k_gpu
     !
     !---------------
     SUBROUTINE get_rho_gpu(rho_loc_d, nrxxs_loc, w1_loc, psic_loc_d)
        !------------
        !
        IMPLICIT NONE
        !
        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(:)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_d(:)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir
        !
        !$cuf kernel do(1)
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir) = rho_loc_d(ir) + &
                         w1_loc * ( DBLE( psic_loc_d(ir) )**2 + &
                                   AIMAG( psic_loc_d(ir) )**2 )
        END DO
        !
     END SUBROUTINE get_rho_gpu
     !
     !-------------
     SUBROUTINE get_rho( rho_loc_h, nrxxs_loc, w1_loc, psic_loc_h )
        !----------
        !
        IMPLICIT NONE
        !
        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_h(nrxxs_loc)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_h(nrxxs_loc)
        INTEGER :: ir
        !
        DO ir = 1, nrxxs_loc
           !
           rho_loc_h(ir) = rho_loc_h(ir) + &
                         w1_loc * ( DBLE( psic_loc_h(ir) )**2 + &
                                   AIMAG( psic_loc_h(ir) )**2 )
           !
        ENDDO
        !
     END SUBROUTINE get_rho
     !
     !-----------------
     SUBROUTINE get_rho_gamma_gpu( rho_loc_d, nrxxs_loc, w1_loc, w2_loc, psic_loc_d )
        !--------------
        !
        IMPLICIT NONE
        !
        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP) :: psic_loc_d(nrxxs_loc)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir
        !
        !$cuf kernel do(1)
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir) = rho_loc_d(ir) + &
                         w1_loc * DBLE( psic_loc_d(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc_d(ir) )**2
           !
        END DO
        !
     END SUBROUTINE get_rho_gamma_gpu
     !
     !--------------
     SUBROUTINE get_rho_domag_gpu( rho_loc_d, nrxxs_loc, w1_loc, psic_loc_d )
        !-----------
        !
        IMPLICIT NONE
        !
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
     !
END SUBROUTINE sum_band_gpu

!----------------------------------------------------------------------------
SUBROUTINE sum_bec_gpu ( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd ) 
  !----------------------------------------------------------------------------
  !
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
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#define cublasZgemm zgemm
#define cublasDgemm dgemm
#endif
  USE kinds,              ONLY : DP
  USE becmod,             ONLY : becp, calbec
  USE control_flags,      ONLY : gamma_only, tqr, offload_type 
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp,               ONLY : nkb, becsum, ebecsum, ofsbeta, &
                                 becsum_d, ebecsum_d, vkb
  USE uspp_param,         ONLY : upf, nh, nhm
  USE wvfct,              ONLY : nbnd, wg, et, current_k
  USE klist,              ONLY : ngk, nkstot
  USE noncollin_module,   ONLY : noncolin, npol
  USE wavefunctions,      ONLY : evc
  USE realus,             ONLY : real_space, &
                                 invfft_orbital_gamma, calbec_rs_gamma, &
                                 invfft_orbital_k, calbec_rs_k
  USE us_exx,             ONLY : store_becxx0
  USE mp_bands,           ONLY : nbgrp,inter_bgrp_comm
  USE mp,                 ONLY : mp_sum
  USE upf_spinorb,        ONLY : fcoef
  !
  ! Used to avoid unnecessary memcopy
  USE xc_lib,             ONLY : xclib_dft_is
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1_d(:,:), auxk2_d(:,:), aux_nc_d(:,:)
  REAL(DP), ALLOCATABLE    :: auxg1_d(:,:), auxg2_d(:,:), aux_gk_d(:,:), aux_egk_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: auxk1_d, auxk2_d, aux_nc_d
  attributes(DEVICE) :: auxg1_d, auxg2_d, aux_gk_d, aux_egk_d
#endif
  INTEGER :: ibnd, kbnd, ibnd_loc, nbnd_loc, ibnd_begin  ! counters on bands
  INTEGER :: npw, ikb, jkb, ih, jh, ijh, na, np, is, js, nhnt, offset
  ! counters on beta functions, atoms, atom types, spin, and auxiliary vars
  !
  CALL start_clock_gpu( 'sum_band:calbec' )
  npw = ngk(ik)
  IF ( .NOT. real_space ) THEN
     !$acc data present_or_copyin(evc) 
     CAll calbec(offload_type, npw, vkb, evc(:,ibnd_start:ibnd_end), becp )
     !$acc end data
  ELSE
     if (gamma_only) then
        do ibnd = ibnd_start, ibnd_end, 2
           call invfft_orbital_gamma(evc,ibnd,ibnd_end) 
           call calbec_rs_gamma(ibnd,ibnd_end,becp%r)
        enddo
        call mp_sum(becp%r,inter_bgrp_comm)
        !$acc update device(becp%r)
     else
        current_k = ik
        becp%k = (0.d0,0.d0)
        do ibnd = ibnd_start, ibnd_end
           call invfft_orbital_k(evc,ibnd,ibnd_end)
           call calbec_rs_k(ibnd,ibnd_end)
        enddo
        call mp_sum(becp%k,inter_bgrp_comm)
        !$acc update device(becp%k)
     endif
  ENDIF
  CALL stop_clock_gpu( 'sum_band:calbec' )
  !
  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  IF(xclib_dft_is('hybrid')) THEN       ! This if condition is not present in the CPU code!! Add it?
     if(allocated(becp%r)) then
       !$acc update self(becp%r)
     elseif(allocated(becp%k)) then
       !$acc update self(becp%k)
     elseif(allocated(becp%nc)) then
       !$acc update self(becp%nc)
     endif
     CALL store_becxx0(ik, becp)
  ENDIF
  !
  CALL start_clock_gpu( 'sum_band:becsum' )
  !
  !$acc data copyin(wg)
  DO np = 1, ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ! allocate work space used to perform GEMM operations
        !
        IF ( gamma_only ) THEN
           nbnd_loc = becp%nbnd_loc
           ALLOCATE( auxg1_d( nbnd_loc, nh(np) ) )
           ALLOCATE( auxg2_d( nbnd_loc, nh(np) ) )
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
        !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
        !
        nhnt = nh(np)
        DO na = 1, nat
           !
           IF (ityp(na)==np) THEN
              !
              ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
              ! copy into aux1, aux2 the needed data to perform a GEMM
              !
              offset = ofsbeta(na)
              IF ( noncolin ) THEN
                 !
                 !$acc parallel loop collapse(2)
                 DO is = 1, npol
                    DO ih = 1, nhnt
                       ikb = offset + ih
                       DO kbnd = 1, this_bgrp_nbnd 
                          ibnd = ibnd_start + kbnd -1 
                          auxk1_d(ibnd,ih+(is-1)*nhnt)= becp%nc(ikb,is,kbnd)
                          auxk2_d(ibnd,ih+(is-1)*nhnt)= wg(ibnd,ik) * &
                                                        becp%nc(ikb,is,kbnd)
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
                 ibnd_begin = becp%ibnd_begin
                 !$acc parallel loop collapse(2)
                 DO ih = 1, nhnt
                    DO ibnd_loc = 1, nbnd_loc
                       ikb = offset + ih
                       ibnd = (ibnd_start -1) + ibnd_loc + ibnd_begin - 1
                       auxg1_d(ibnd_loc,ih) = becp%r(ikb,ibnd_loc)
                       auxg2_d(ibnd_loc,ih) = becp%r(ikb,ibnd_loc) * wg(ibnd,ik)
                    END DO
                 END DO
                 CALL cublasDgemm ( 'T', 'N', nhnt, nhnt, nbnd_loc, &
                      1.0_dp, auxg1_d, nbnd_loc,    &
                      auxg2_d, nbnd_loc, 0.0_dp, aux_gk_d, nhnt )
                 !
                 if (tqr) then
                   !$acc parallel loop collapse(2)
                   DO ih = 1, nhnt
                      DO ibnd_loc = 1, nbnd_loc
                      ibnd = (ibnd_start -1) + ibnd_loc + ibnd_begin - 1
                      auxg2_d(ibnd_loc,ih) = et(ibnd,ik) * auxg2_d(ibnd_loc,ih)
                      END DO
                   END DO

                   CALL cublasDgemm ( 'T', 'N', nhnt, nhnt, nbnd_loc, &
                        1.0_dp, auxg1_d, nbnd_loc,    &
                        auxg2_d, nbnd_loc, 0.0_dp, aux_egk_d, nhnt )
                 end if
                 !
              ELSE
                 !
                 !$acc parallel loop collapse(2)
                 DO ih = 1, nhnt
                    DO kbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
                       ibnd = ibnd_start + kbnd -1 
                       ikb = offset + ih
                       auxk1_d(ibnd,ih) = becp%k(ikb,kbnd) 
                       auxk2_d(ibnd,ih) = wg(ibnd,ik)*becp%k(ikb,kbnd)
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
                   !$acc parallel loop collapse(2)
                   DO ih = 1, nhnt
                      DO ibnd = ibnd_start, ibnd_end
                         auxk2_d(ibnd,ih) = et(ibnd,ik)*auxk2_d(ibnd,ih)
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
!$acc host_data use_device(fcoef)
                 CALL add_becsum_so_gpu (na, np, fcoef, aux_nc_d, becsum_d )
!$acc end host_data
              ELSE
                 !
                 !$acc parallel loop collapse(2)
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
           DEALLOCATE( auxg2_d, auxg1_d )
        ELSE
           DEALLOCATE( auxk2_d, auxk1_d )
        END IF
        !
     END IF
     !
  END DO
  !$acc end data
  !
  ! sync 
  if (nhm > 0) then
     becsum=becsum_d
     if (tqr) ebecsum=ebecsum_d
  endif
  !
  CALL stop_clock_gpu( 'sum_band:becsum' )
  !
END SUBROUTINE sum_bec_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_nc_gpu ( na, np, becsum_nc_d, becsum_d )
!----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the
  !! Pauli matrices, saves it in \(\text{becsum}\) for the calculation of 
  !! augmentation charge and magnetization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE uspp,                 ONLY : ijtoh_d
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
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
SUBROUTINE add_becsum_so_gpu( na, np, fcoef_d, becsum_nc_d, becsum_d )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the Pauli
  !! matrices, rotates it as appropriate for the spin-orbit case, saves it in 
  !! \(\text{becsum}\) for the calculation of augmentation charge and magnetization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
  USE uspp,                 ONLY : ijtoh_d, nhtol_d, nhtoj_d, indv_d
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: fcoef_d(nhm,nhm,2,2,ntyp)
  !! function needed to account for spinors.
!$acc declare deviceptr(fcoef_d)
  COMPLEX(DP), INTENT(IN) :: becsum_nc_d(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum_d(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, lh, kh, ijh, is1, is2, nhnt
  COMPLEX(DP) :: fac
#if defined(__CUDA)
  attributes (DEVICE) :: fcoef_d, becsum_nc_d, becsum_d
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

