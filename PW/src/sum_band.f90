!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
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
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr, sic
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_rho,              ONLY : rho_g2r, rho_r2g
  USE fft_wave,             ONLY : wave_g2r
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho, rhoz_or_updw
  USE sic_mod,              ONLY : isp, pol_type
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, okvan
  USE uspp_param,           ONLY : nh, nhm
  USE wavefunctions,        ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  USE xc_lib,               ONLY : xclib_dft_is
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                                   becp
  USE gcscf_module,         ONLY : lgcscf, gcscf_calc_nelec
  USE add_dmft_occ,         ONLY : dmft, dmft_updated, v_dmft
#if defined (__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_sum_band
#endif
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
  REAL(dp), EXTERNAL :: e_band
  COMPLEX(DP), ALLOCATABLE :: psic(:,:)
  !! Work space used for FFTs in this routine
  !
  CALL start_clock( 'sum_band' )
  !
  ALLOCATE( psic(dfftp%nnr,npol) )
  !$acc enter data create(psic)
  !
  IF ( nhm > 0 ) THEN
     ! Note: becsum and ebecsum are computed on GPU and copied to CPU
     !$acc kernels
     becsum(:,:,:) = 0.D0
     !$acc end kernels
     IF (tqr) THEN
        !$acc kernels
        ebecsum(:,:,:) = 0.D0
        !$acc end kernels
     END IF
  ENDIF
  ! FIXME: next line  needed for old (21.7 or so) NVIDIA compilers
  !$acc enter data create(rho)
  ! FIXME: rho should be on GPU always, not just in this routine
  !$acc enter data create(rho%of_r)
  !$acc kernels
  rho%of_r(:,:) = 0.0_DP
  !$acc end kernels
  rho%of_g(:,:) = (0.0_dp, 0.0_dp)
  IF ( xclib_dft_is('meta') .OR. lxdm ) THEN
     !$acc enter data create(rho%kin_r)
     !$acc kernels
     rho%kin_r(:,:) = 0.0_dp
     !$acc end kernels
     rho%kin_g(:,:) = (0.0_dp, 0.0_dp)
  ENDIF
  IF (sic) THEN
     rho%pol_r(:,:) = 0.0_dp
     rho%pol_g(:,:) = (0.0_dp,0.0_dp)
  END IF
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL start_clock( 'sum_band:weights' )
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
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_sum_band(oscdft_ctx)
#endif
  !
  ! ... for band parallelization: set band computed by this processor
  !
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type_acc( nkb, this_bgrp_nbnd, becp, intra_bgrp_comm )
  !
  ! ... specialized routines are called to sum at Gamma or for each k point 
  ! ... the contribution of the wavefunctions to the charge
  !
  CALL start_clock( 'sum_band:loop' )
  !
  IF ( ( dffts%has_task_groups ) .AND. .NOT. &
       ( xclib_dft_is('meta') .OR. lxdm .OR. dmft .OR. sic ) ) THEN
     ! Task groups: no GPU, no special cases, etc.
     IF ( gamma_only ) THEN
        !
        CALL sum_band_gamma_tg ()
        !
     ELSE
        !
        CALL sum_band_k_tg()
        !
     ENDIF
  ELSE
     ! General case, but no task groups
     IF ( gamma_only ) THEN
        !
        CALL sum_band_gamma()
        !
     ELSE
        !
        CALL sum_band_k()
        !
     ENDIF
  END IF
  !
  CALL stop_clock( 'sum_band:loop' )
  !
  !$acc exit data delete(psic)
  DEALLOCATE(psic)
  !
  ! ... Compute here the sum of band eigenvalues "eband"
  !
  eband =  e_band( )
  !
  ! ... sum charge density over pools (distributed k-points) and bands
  !
  !$acc update host(rho%of_r)
  !$acc exit data delete(rho%of_r)
  !$acc exit data delete(rho)
  CALL mp_sum( rho%of_r, inter_pool_comm )
  CALL mp_sum( rho%of_r, inter_bgrp_comm )
  IF (sic) then
     CALL mp_sum( rho%pol_r, inter_pool_comm )
     CALL mp_sum( rho%pol_r, inter_bgrp_comm )
  END IF
  !
  ! ... bring the unsymmetrized rho(r) to G-space
  !
  CALL rho_r2g( dffts, rho%of_r, rho%of_g )
  IF(sic) CALL rho_r2g( dffts, rho%pol_r, rho%pol_g )
  !
  IF( okvan )  THEN
     !
     ! ... becsum = \sum_k\sum_i w_{k,i} <\psi_{k,i}|\beta_l><\beta_m|psi_{k,i}> is computed
     ! ... in sum_band_*. Use host copy to do the communications
     !$acc update host(becsum)
     if (tqr) then
        !$acc update host(ebecsum)
     endif
     !
     CALL deallocate_bec_type_acc ( becp )
     !
     ! ... becsums must be summed over bands (with bgrp parallelization)
     ! ... and over k-points (unsymmetrized). Then the CPU and GPU copies are aligned.
     !
     CALL mp_sum(becsum, inter_bgrp_comm )
     CALL mp_sum(becsum, inter_pool_comm )
     !$acc update device(becsum)
     !
     ! ... same for ebecsum, a correction needed for real-space algorithm
     !
     IF (tqr) THEN
        CALL mp_sum(ebecsum, inter_pool_comm )
        CALL mp_sum(ebecsum, inter_bgrp_comm )
        !$acc update device(ebecsum)
     ENDIF
     !
     ! ... Needed for PAW: becsums are stored into rho%bec and symmetrized so that they reflect
     ! ... a real integral in k-space, not only on the irreducible zone. 
     ! ... For USPP there is no need to do this as becsums are only used
     ! ... to compute the density, which is symmetrized later.
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
  CALL start_clock( 'sum_band:sym_rho' )
  CALL sym_rho( nspin_mag, rho%of_g )
  IF (sic) CALL sym_rho ( nspin_mag, rho%pol_g )
  !
  ! ... synchronize rho%of_r to the calculated rho%of_g
  !
  CALL rho_g2r( dfftp, rho%of_g, rho%of_r )
  IF(sic) CALL rho_g2r( dfftp, rho%pol_g, rho%pol_r )
  !
  ! ... rho_kin(r): sum over bands, k-points, bring to G-space, symmetrize,
  ! ... synchronize with rho_kin(G)
  !
  IF ( xclib_dft_is('meta') .OR. lxdm) THEN
     !
     !$acc exit data delete(rho%kin_r)
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
  CALL stop_clock( 'sum_band:sym_rho' )
  !
  ! ... if LSDA rho%of_r and rho%of_g are converted from (up,dw) to
  ! ... (up+dw,up-dw) format.
  !
  IF ( nspin == 2 ) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
  IF (sic) THEN
    rho%pol_r(:,2) = rho%pol_r(:,1) 
    rho%pol_g(:,2) = rho%pol_g(:,1) 
  END IF
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
       USE uspp_init,      ONLY : init_us_2
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, incr, j
       INTEGER ::  ierr, ebnd, i, brange
       REAL(DP) :: kplusgi
       COMPLEX(DP), ALLOCATABLE :: grad_psic(:,:)
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       incr = 2
       IF (xclib_dft_is('meta') .OR. lxdm) THEN
          ALLOCATE( grad_psic(npwx,2) )
          !$acc enter data create(grad_psic)
       ENDIF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          CALL start_clock( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !$acc update device(evc)
          !
          CALL stop_clock( 'sum_band:buffer' )
          !
          CALL start_clock( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
          !
          CALL stop_clock( 'sum_band:init_us_2' )
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             ebnd = ibnd
             IF ( ibnd < ibnd_end ) ebnd = ebnd + 1
             !
             CALL wave_g2r( evc(1:npw,ibnd:ebnd), psic(:,1), dffts )
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
             CALL get_rho_gamma( rho%of_r(:,current_spin), dffts%nnr, w1, w2, psic )
             !
             IF (xclib_dft_is('meta') .OR. lxdm) THEN
                !
                !$acc data present(g,igk_k,evc,psic,rho%kin_r) copyin(xk)
                DO j = 1, 3
                   !$acc parallel loop
                   DO i = 1, npw
                      kplusgi = ( xk(j,ik) + g(j,igk_k(i,ik)) ) * tpiba
                      grad_psic(i,1) = CMPLX(0._DP,kplusgi,KIND=DP) * evc(i,ibnd)
                      IF ( ibnd < ibnd_end ) grad_psic(i,2) = &
                           CMPLX(0._DP,kplusgi,KIND=DP) * evc(i,ibnd+1)
                   ENDDO
                   !
                   ebnd = ibnd
                   IF ( ibnd < ibnd_end ) ebnd = ebnd + 1
                   brange = ebnd-ibnd+1
                   !
                   CALL wave_g2r( grad_psic(1:npw,1:brange), psic(:,1), dffts )
                   !
                   ! ... increment the kinetic energy density ...
                   !  
                   CALL get_rho_gamma( rho%kin_r(:,current_spin), dffts%nnr, w1, w2, psic )
                ENDDO
                !$acc end data
             ENDIF
             !
          ENDDO
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec( ik, current_spin, ibnd_start, ibnd_end, &
               this_bgrp_nbnd )
          !
       ENDDO k_loop
       !
       IF (xclib_dft_is('meta') .OR. lxdm) THEN
          !$acc exit data delete(grad_psic)
          DEALLOCATE( grad_psic )
          !$acc update host(rho%kin_r)
       END IF
       RETURN
       !
     END SUBROUTINE sum_band_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for k-points version
       !
       USE uspp_init,              ONLY : init_us_2
       USE klist,                  ONLY : nelec
       USE control_flags,          ONLY : many_fft
       USE wavefunctions,          ONLY : psic_nc
       !! psic_nc is allocated work space for noncolinear case
       !! FIXME: too many work spaces!
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, incr
       COMPLEX(DP), ALLOCATABLE :: psicd(:)
       COMPLEX(DP), ALLOCATABLE :: grad_psic(:,:)
       ! polaron calculation
       REAL(DP), ALLOCATABLE :: rho_p(:)
       COMPLEX(DP), ALLOCATABLE :: psic_p(:)
       !$acc declare device_resident(rho_p, psic_p)
       INTEGER :: ierr
       INTEGER :: i, j, group_size, hm_vec(3)
       REAL(DP) :: kplusgi
       REAL(DP) :: wg_p
       INTEGER  :: ibnd_p
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF(sic) THEN
          ALLOCATE(rho_p(dffts%nnr))
          ALLOCATE(psic_p(dffts%nnr*2)) ! why *2?
          !$acc kernels
          rho_p=0.0_dp
          !$acc end kernels
          wg_p = 0.0
          ibnd_p = nelec/2+1 
       END IF
       !
       IF (noncolin) THEN
          incr = 1
          !$acc enter data create(psic_nc)
       ELSE IF (xclib_dft_is('meta') .OR. lxdm) THEN
          ! many_fft cannot be used with meta-GGA and XDM
          incr = 1
       ELSE
          incr = many_fft
       ENDIF
       IF (xclib_dft_is('meta') .OR. lxdm) THEN
          ALLOCATE( grad_psic(npwx,incr) )
          !$acc enter data create(grad_psic)
       END IF
       !
       ALLOCATE( psicd(dffts%nnr*incr) )
       !$acc data create(psicd)
       !
       ! ... This is used as reduction variable on the device
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          CALL start_clock( 'sum_band:buffer' )
          IF ( nks > 1 ) THEN
             CALL get_buffer( evc, nwordwfc, iunwfc, ik )
          ENDIF
          !$acc update device(evc)
          CALL stop_clock( 'sum_band:buffer' )
          !
          CALL start_clock( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
          !
          CALL stop_clock( 'sum_band:init_us_2' )
          !
          ! ... for DMFT the eigenvectors are updated using v_dmft from add_dmft_occ.f90
          !
          IF ( dmft .AND. .NOT. dmft_updated) THEN
             ! BEWARE: untested on GPUs - should work if v_dmft is present on host  
             !$acc update host(evc)
             DO j = 1, npw
                CALL ZGEMM( 'T', 'N', nbnd, 1, nbnd, (1.d0,0.d0), v_dmft(:,:,ik), &
                            nbnd, evc(j,:), nbnd, (0.d0,0.d0), evc(j,:), nbnd )
             ENDDO
             !
             IF ( nks > 1 ) THEN
                  CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
             ELSE
                !$acc update device(evc)
             END IF
             !
          END IF
          !
          ! ... calculate polaron density
          !
          IF ( sic .AND. current_spin==isp ) THEN
             CALL wave_g2r( evc(1:npw,ibnd_p:ibnd_p), psic_p, dffts, igk=igk_k(:,ik) )
             !
             CALL get_rho_k(rho_p, dffts%nnr, wg(1,ik)/omega, psic_p)
             !$acc update host(rho_p)
             rho%pol_r(:,1) = rho_p(:)
             wg_p = wg_p + wg(ibnd_p,ik)
          ENDIF
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
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
                   CALL get_rho_k( rho%of_r(:,1), dffts%nnr, w1, psic_nc(:,ipol) )
                ENDDO
                !
                ! ... In this case, calculate also the three
                ! ... components of the magnetization (stored in rho%of_r(ir,2-4))
                !
                IF (domag) &
                   CALL get_rho_domag( rho%of_r(:,:), dffts%nnr, w1, psic_nc(1:,1:) )
                !
             ELSE IF ( incr > 1 ) THEN
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
                   CALL get_rho_k( rho%of_r(:,current_spin), dffts%nnr, w1, psicd(i*dffts%nnr+1:) )
                ENDDO
                !
             ELSE
                !
                CALL wave_g2r( evc(1:npw,ibnd:ibnd), psicd, dffts, igk=igk_k(:,ik) )
                !
                ! ... increment the charge density ...
                !
                CALL get_rho_k( rho%of_r(:,current_spin), dffts%nnr, w1, psicd )
                !
                IF (xclib_dft_is('meta') .OR. lxdm) THEN
                   !$acc data present(g,igk_k,evc,rho%kin_r) copyin(xk)
                   DO j=1,3
                      !$acc parallel loop
                      DO i = 1, npw
                         kplusgi = (xk(j,ik)+g(j,igk_k(i,ik))) * tpiba
                         grad_psic(i,1) = CMPLX(0._DP,kplusgi,KIND=DP) * evc(i,ibnd)
                      ENDDO
                      !
                      CALL wave_g2r( grad_psic(1:npw,1:1), psicd, dffts, igk=igk_k(:,ik) )
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      CALL get_rho_k( rho%kin_r(:,current_spin), dffts%nnr, w1, psicd )
                   ENDDO
                   !$acc end data
                ENDIF
                !
             ENDIF
             !
          ENDDO
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       !$acc end data
       !
       IF(sic .and. pol_type == 'h') THEN
          wg_p = 1.0 - wg_p
          DEALLOCATE(psic_p)
          DEALLOCATE(rho_p)
       END IF
       !
       DEALLOCATE( psicd )
       !
       IF (xclib_dft_is('meta') .OR. lxdm) THEN
          !$acc exit data delete(grad_psic)
          DEALLOCATE( grad_psic )
          !$acc update host(rho%kin_r)
       END IF
       IF ( noncolin) THEN
          !$acc exit data delete(psic_nc)
       ENDIF
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
     !---------------
     SUBROUTINE get_rho_k(rho_loc, nrxxs_loc, w1_loc, psic_loc)
        !------------
        !
        IMPLICIT NONE
        !
        INTEGER :: nrxxs_loc
        REAL(DP), INTENT(INOUT) :: rho_loc(:)
        REAL(DP) :: w1_loc
        COMPLEX(DP), INTENT(IN) :: psic_loc(:)
        INTEGER :: ir
        !
        !$acc data present(rho_loc, psic_loc)
        !$acc parallel loop
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * ( DBLE( psic_loc(ir) )**2 + &
                                   AIMAG( psic_loc(ir) )**2 )
        END DO
        !$acc end data
        !
     END SUBROUTINE get_rho_k
     !
     !-----------------
     SUBROUTINE get_rho_gamma( rho_loc, nrxxs_loc, w1_loc, w2_loc, psic_loc )
        !--------------
        !
        IMPLICIT NONE
        !
        INTEGER :: nrxxs_loc
        REAL(DP), INTENT(INOUT) :: rho_loc(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP), INTENT(IN) :: psic_loc(nrxxs_loc)
        INTEGER :: ir
        !
        !$acc data present(rho_loc, psic_loc)
        !$acc parallel loop
        DO ir = 1, nrxxs_loc
           !
           rho_loc(ir) = rho_loc(ir) + &
                         w1_loc * DBLE( psic_loc(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc(ir) )**2
           !
        END DO
        !$acc end data
        !
     END SUBROUTINE get_rho_gamma
     !
     !--------------
     SUBROUTINE get_rho_domag( rho_loc, nrxxs_loc, w1_loc, psic_loc )
        !-----------
        !
        IMPLICIT NONE
        !
        INTEGER :: nrxxs_loc
        REAL(DP), INTENT(INOUT) :: rho_loc(:, :)
        REAL(DP) :: w1_loc
        COMPLEX(DP), INTENT(IN) :: psic_loc(:, :)
        INTEGER :: ir

        !$acc data present( rho_loc, psic_loc)
        !$acc parallel loop
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
        !$acc end data

     END SUBROUTINE get_rho_domag
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_gamma_tg()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for gamma version.
       !
       USE fft_helper_subroutines, ONLY : fftx_ntgrp, fftx_tgpe, &
                                          tg_reduce_rho, tg_get_group_nr3
       USE fft_wave,               ONLY : tgwave_g2r
       USE uspp_init,              ONLY : init_us_2
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, incr, v_siz, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:)
       INTEGER :: right_nr3, ntgrp, ierr, i
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       incr = 2
       !
       v_siz = dffts%nnr_tg 
       !
       ALLOCATE( tg_psi( v_siz ) )
       ALLOCATE( tg_rho( v_siz ) )
       !
       incr = 2 * fftx_ntgrp(dffts)
       !
       k_loop: DO ik = 1, nks
          !
          tg_rho = 0.0_DP
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          CALL start_clock( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          !
          CALL stop_clock( 'sum_band:buffer' )
          !
          CALL start_clock( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
          !
          CALL stop_clock( 'sum_band:init_us_2' )
          !
          DO ibnd = ibnd_start, ibnd_end, incr
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
             CALL get_rho_gamma( tg_rho, dffts%nr1x*dffts%nr2x*right_nr3, &
                  w1, w2, tg_psi )
             !
          ENDDO
          !
          CALL tg_reduce_rho( rho%of_r, tg_rho, current_spin, dffts )
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec( ik, current_spin, ibnd_start, ibnd_end, &
                                         this_bgrp_nbnd )
          !
       ENDDO k_loop
       !
       !
       DEALLOCATE( tg_psi )
       DEALLOCATE( tg_rho )
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma_tg
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k_tg()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for k-points version
       !
       USE fft_helper_subroutines, ONLY : fftx_ntgrp, fftx_tgpe, &
                                          tg_reduce_rho, tg_get_group_nr3
       USE fft_wave,               ONLY : tgwave_g2r
       USE uspp_init,              ONLY : init_us_2
       USE klist,                  ONLY : nelec
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, incr, v_siz
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_psi_nc(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:), tg_rho_nc(:,:)
       !
       INTEGER :: right_nr3, ntgrp, ierr
       INTEGER :: i, j, group_size, hm_vec(3)
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       incr = 1
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
       incr = fftx_ntgrp(dffts)
       !
       k_loop: DO ik = 1, nks
          !
          IF (noncolin) THEN
             tg_rho_nc = 0.0_DP
          ELSE
             tg_rho = 0.0_DP
          ENDIF
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          CALL start_clock( 'sum_band:buffer' )
          IF ( nks > 1 ) THEN
             CALL get_buffer( evc, nwordwfc, iunwfc, ik )
          ENDIF
          CALL stop_clock( 'sum_band:buffer' )
          !
          CALL start_clock( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
          !
          CALL stop_clock( 'sum_band:init_us_2' )
          !
          ! ... for DMFT the eigenvectors are updated using v_dmft from add_dmft_occ.f90
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
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
                DO ipol = 1, npol
                   CALL get_rho_k( tg_rho_nc(:,1), dffts%nr1x*dffts%nr2x* &
                        right_nr3, w1, tg_psi_nc(:,ipol) )
                ENDDO
                !
                IF (domag) CALL get_rho_domag( tg_rho_nc(:,:), dffts%nr1x* &
                     dffts%nr2x*dffts%my_nr3p, w1, tg_psi_nc(:,:) )
                !
                !
             ELSE
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
                CALL get_rho_k( tg_rho, dffts%nr1x*dffts%nr2x*right_nr3, w1, tg_psi )
             END IF
             !
          ENDDO
          !
          ! ... reduce the charge across task group
          !
          IF (.not.noncolin) THEN
             CALL tg_reduce_rho( rho%of_r, tg_rho_nc, tg_rho, current_spin, &
                  noncolin, domag, dffts )
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd ) 
          !
       END DO k_loop
       !
       IF (noncolin) THEN
          DEALLOCATE( tg_psi_nc )
          DEALLOCATE( tg_rho_nc )
       ELSE
          DEALLOCATE( tg_psi )
          DEALLOCATE( tg_rho )
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k_tg
     !
END SUBROUTINE sum_band
!
!----------------------------------------------------------------------------
FUNCTION e_band ( )
  !----------------------------------------------------------------------------
  !
  !! Calculation of band energy sum - Paolo Giannozzi Oct. 2024
  !! To be called after "weights" 
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : et, wg, nbnd
  !
  IMPLICIT NONE
  !
  REAL(dp) :: e_band
  INTEGER  :: ik, ibnd
  !
  e_band = 0.0_dp
  !
  k_loop: DO ik = 1, nks
     !
     band_loop: DO ibnd = 1, nbnd
        !
        e_band = e_band + wg(ibnd,ik) * et(ibnd,ik)
        !
     END DO band_loop
     !
  END DO k_loop
  !
  CALL mp_sum (e_band, inter_pool_comm)
  !
END FUNCTION e_band
