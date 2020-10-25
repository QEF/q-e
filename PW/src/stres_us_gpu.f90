!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_us_gpu( ik, gk_d, sigmanlc )
  !----------------------------------------------------------------------------
  !! nonlocal (separable pseudopotential) contribution to the stress
  !! NOTICE: sum of partial results over procs is performed in calling routine
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k_d
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf, lmaxkb, nh, nhm
  USE uspp,                 ONLY : nkb
  USE spin_orb,             ONLY : lspinorb
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum, mp_get_comm_null, &
                                   mp_circular_shift_left 
  !
  USE wavefunctions_gpum,   ONLY : using_evc, using_evc_d, evc_d
  USE wvfct_gpum,           ONLY : using_et
  USE uspp_gpum,            ONLY : vkb_d, using_vkb, using_vkb_d, &
                                   deeq_d, using_deeq_d
  USE becmod_gpum,          ONLY : becp_d, bec_type_d
  USE becmod_subs_gpum,     ONLY : using_becp_auto, using_becp_d_auto, &
                                   calbec_gpu
  USE device_fbuff_m,             ONLY : dev_buf
  USE device_memcpy_m,        ONLY : dev_memcpy
  !
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ik
  !! k-point index
  REAL(DP), INTENT(IN)    :: gk_d(npwx,3)
  !! wave function components for fixed k-point
  REAL(DP), INTENT(INOUT) :: sigmanlc(3,3)
  !! stress tensor, non-local contribution
  !
  ! ... local variables
  !
  REAL(DP), POINTER :: qm1_d(:)
  REAL(DP) :: q
  INTEGER  :: npw , iu, np, ierr
  !
  INTEGER :: na1, np1, nh_np1, ijkb01, itot
  LOGICAL :: ismulti_np
  INTEGER, ALLOCATABLE :: shift(:)
  INTEGER, ALLOCATABLE :: ix_d(:,:), ityp_d(:), nh_d(:), shift_d(:)
  LOGICAL, ALLOCATABLE :: is_multinp_d(:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: gk_d, qm1_d, is_multinp_d, ix_d, shift_d, &
                        ityp_d, nh_d
#endif 
  !
  CALL using_evc_d(0)
  CALL using_evc(0)
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  IF ( nks > 1 ) THEN
    CALL using_vkb_d(1)
    CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
  ENDIF
  !
  CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm ) 
  CALL using_becp_auto(2)
  !
  CALL using_evc_d(0)
  CALL using_vkb(0)
  CALL using_vkb_d(0)
  CALL using_becp_d_auto(2)
  !
  CALL calbec_gpu( npw, vkb_d, evc_d, becp_d )
  !
  CALL dev_buf%lock_buffer( qm1_d, npwx, ierr )
  !$cuf kernel do (1) <<<*,*>>>
  DO iu = 1, npw
     q = SQRT( gk_d(iu,1)**2 + gk_d(iu,2)**2 + gk_d(iu,3)**2 )
     IF ( q > eps8 ) THEN
        qm1_d(iu) = 1._DP / q
     ELSE
        qm1_d(iu) = 0._DP
     ENDIF
  ENDDO
  !
  !----------define index arrays (type, atom, etc.) for cuf kernel loops------
  ALLOCATE( is_multinp_d(nat*nhm) )
  ALLOCATE( ix_d(nat*nhm,4), ityp_d(nat), nh_d(ntyp) )
  ityp_d = ityp  ;  nh_d = nh
  !
  ALLOCATE( shift(nat) )
  ijkb01 = 0 
  DO iu = 1, ntyp
    DO na1 = 1, nat
      IF (ityp(na1) == iu ) THEN
         shift(na1) = ijkb01
         ijkb01 = ijkb01 + nh(iu)
      ENDIF
    ENDDO
  ENDDO
  !
  ALLOCATE( shift_d(nat) )
  shift_d = shift
  !
  ijkb01 = 0
  itot=0
  DO np = 1, ntyp
    DO na1 = 1, nat
      np1 = ityp(na1)
      IF (np /= np1) CYCLE
      ijkb01 = shift(na1)
      nh_np1 = nh(np1)
      ismulti_np = upf(np1)%tvanp .OR. upf(np1)%is_multiproj
      IF ( .NOT. ismulti_np .AND. noncolin .AND. lspinorb ) &
                                   CALL errore('stres_us','wrong case',1)
      !$cuf kernel do (1) <<<*,*>>>
      DO iu = itot+1, itot+nh_np1
        ix_d(iu,1) = na1                !na 
        ix_d(iu,2) = iu-itot            !ih
        ix_d(iu,3) = nh_np1             !nh(np)
        ix_d(iu,4) = ijkb01             !ishift
        is_multinp_d(iu) = ismulti_np
      ENDDO 
      itot = itot + nh_np1
    ENDDO
  ENDDO
  !--------------
  !
  !
  IF ( gamma_only ) THEN
     !
     CALL stres_us_gamma_gpu()
     !
  ELSE
     !
     CALL stres_us_k_gpu()
     !
  ENDIF
  !
  CALL dev_buf%release_buffer( qm1_d, ierr )
  !
  DEALLOCATE( is_multinp_d )
  DEALLOCATE( ix_d, ityp_d, nh_d, shift_d )
  DEALLOCATE( shift )
  !
  CALL deallocate_bec_type( becp ) 
  CALL using_becp_auto(2)
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_gamma_gpu()
       !-----------------------------------------------------------------------
       !! nonlocal contribution to the stress - gamma version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: na, np, nt, ibnd, ipol, jpol, l, i, ikb,  &
                   jkb, ih, jh, ibnd_loc,ijkb0,nh_np, nproc, &
                   nbnd_loc, nbnd_begin, icyc, ishift, nhmx
       REAL(DP) :: dot11, dot21, dot31, dot22, dot32, dot33, &
                   qm1i, gk1, gk2, gk3, wg_nk, fac, evps, aux
       COMPLEX(DP) :: worksum, cv, wsum1, wsum2, wsum3, ps, evci
       !
       COMPLEX(DP), POINTER :: ps_d(:), deff_d(:,:,:), dvkb_d(:,:,:)
       REAL(DP),    POINTER :: becpr_d(:,:)
       !
       REAL(DP) :: xyz(3,3)
       INTEGER  :: ierrs(3)
       !
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1._DP, 0._DP, 0._DP, 0._DP, 1._DP, 0._DP, 0._DP, 0._DP, &
                  1._DP /
#if defined(__CUDA)
       attributes(DEVICE) :: deff_d, becpr_d, ps_d, dvkb_d
#endif
       !
       CALL using_becp_auto(0)
       !
       IF( becp%comm /= mp_get_comm_null() ) THEN
          nproc      = becp%nproc
          nbnd_loc   = becp%nbnd_loc
          nbnd_begin = becp%ibnd_begin
          IF( ( nbnd_begin + nbnd_loc - 1 ) > nbnd ) nbnd_loc = nbnd - nbnd_begin + 1
       ELSE
          nproc      = 1
          nbnd_loc   = nbnd
          nbnd_begin = 1
       ENDIF
       !
       CALL dev_buf%lock_buffer( deff_d, (/ nhm,nhm,nat /), ierrs(1) )
       !
       ! ... diagonal contribution - if the result from "calbec" are not 
       ! ... distributed, must be calculated on a single processor
       !
       ! ... for the moment when using_gpu is true becp is always fully present in all processors
       !
       CALL using_et(0) ! compute_deff : intent(in)
       !
       CALL dev_buf%lock_buffer( ps_d, nkb, ierrs(2) )
       IF (ANY(ierrs(1:2) /= 0)) CALL errore( 'stres_us_gpu', 'cannot allocate buffers', -1 )
       !
       becpr_d => becp_d%r_d 
       !
       evps = 0._DP
       !
       compute_evps: IF ( .NOT. (nproc == 1 .AND. me_pool /= root_pool) ) THEN
         !
         DO ibnd_loc = 1, nbnd_loc
            ibnd = ibnd_loc + becp%ibnd_begin - 1 
            CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
            wg_nk = wg(ibnd,ik)
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itot
              !
              ih = ix_d(i,2)       ; na = ix_d(i,1)
              ishift = ix_d(i,4)   ; ikb = ishift + ih
              !
              IF (.NOT. is_multinp_d(i)) THEN
                 aux = wg_nk * DBLE(deff_d(ih,ih,na)) * &
                                       ABS(becpr_d(ikb,ibnd_loc))**2
              ELSE
                 nh_np = ix_d(i,3)
                 !
                 aux = wg_nk * DBLE(deff_d(ih,ih,na))         &
                                     * ABS(becpr_d(ikb,ibnd_loc))**2  &
                             +  becpr_d(ikb,ibnd_loc)* wg_nk * 2._DP  &
                                * SUM( DBLE(deff_d(ih,ih+1:nh_np,na)) &
                                * becpr_d(ishift+ih+1:ishift+nh_np,ibnd_loc))
              ENDIF
              evps = evps + aux
              !
            ENDDO
            !
         ENDDO
         !
       ENDIF compute_evps
       !
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !------------------------------------
       CALL dev_buf%lock_buffer( dvkb_d, (/ npwx,nkb,4 /), ierrs(3) )
       IF (ierrs(3) /= 0) CALL errore( 'stres_us_gpu', 'cannot allocate buffers', -1 )
       !
       CALL gen_us_dj_gpu( ik, dvkb_d(:,:,4) )
       IF ( lmaxkb > 0 ) THEN 
         DO ipol = 1, 3
           CALL gen_us_dy_gpu( ik, xyz(1,ipol), dvkb_d(:,:,ipol))
         ENDDO
       ENDIF
       !
       !
       DO icyc = 0, nproc -1
          !
          DO ibnd_loc = 1, nbnd_loc
             !  
             ibnd = ibnd_loc + becp%ibnd_begin - 1 
             CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
             !
             !$cuf kernel do (1) <<<*,*>>>
             DO i = 1, itot
               !
               ih = ix_d(i,2)     ; na = ix_d(i,1)
               ishift = ix_d(i,4) ; ikb = ishift + ih
               !
               IF (.NOT. is_multinp_d(i)) THEN
                  ps_d(ikb) =  deff_d(ih,ih,na) * CMPLX(becpr_d(ikb,ibnd_loc))
               ELSE
                  nh_np = ix_d(i,3)
                  !
                  ps_d(ikb) = CMPLX( SUM( becpr_d(ishift+1:ishift+nh_np,ibnd_loc) &
                                    * DBLE(deff_d(ih,1:nh_np,na))))
               ENDIF
               !
             ENDDO
             !
             dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
             dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
             !
             !$cuf kernel do(2) <<<*,*>>> 
             DO na =1, nat 
                DO i = 1, npw
                   worksum = (0._DP,0._DP) 
                   np = ityp_d(na) 
                   ijkb0 = shift_d(na)
                   nh_np = nh_d(np)
                   DO ih = 1, nh_np
                      ikb = ijkb0 + ih  
                      worksum = worksum + ps_d(ikb) * dvkb_d(i,ikb,4) 
                   ENDDO
                   evci = evc_d(i,ibnd) 
                   gk1  = gk_d(i,1) 
                   gk2  = gk_d(i,2) 
                   gk3  = gk_d(i,3) 
                   qm1i = qm1_d(i) 
                   !  
                   cv = evci * CMPLX(gk1 * gk1  * qm1i)
                   dot11 = dot11 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv) 
                   !  
                   cv = evci * CMPLX(gk2 * gk1 * qm1i )
                   dot21 = dot21 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   ! 
                   cv = evci * CMPLX(gk3 * gk1 * qm1i)
                   dot31 = dot31 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   ! 
                   cv = evci * CMPLX(gk2 * gk2 * qm1i)
                   dot22 = dot22 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   ! 
                   cv = evci * CMPLX(gk3 * gk2 * qm1i)
                   dot32 = dot32 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk3 * gk3 * qm1i)
                   dot33 = dot33 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                ENDDO
             ENDDO
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             sigmanlc(:,1) = sigmanlc(:,1) - 4._DP * wg(ibnd,ik) * [dot11, dot21, dot31] 
             sigmanlc(:,2) = sigmanlc(:,2) - 4._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
             sigmanlc(:,3) = sigmanlc(:,3) - 4._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]                  

             !
             ! ... non diagonal contribution - derivative of the spherical harmonics
             ! ... (no contribution from l=0)
             !
             IF ( lmaxkb == 0 ) CYCLE 
             !
             !------------------------------------
             !
             dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
             dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
             !
             !$cuf kernel do(2) <<<*,*>>> 
             DO ikb = 1, nkb 
                DO i = 1, npw  
                   !
                   wsum1 = ps_d(ikb)*dvkb_d(i,ikb,1) 
                   wsum2 = ps_d(ikb)*dvkb_d(i,ikb,2) 
                   wsum3 = ps_d(ikb)*dvkb_d(i,ikb,3)      
                   !
                   evci = evc_d(i,ibnd) 
                   gk1 = gk_d(i,1)
                   gk2 = gk_d(i,2)
                   gk3 = gk_d(i,3) 
                   !
                   cv = evci * CMPLX(gk1 )
                   dot11 = dot11 + DBLE( wsum1)* DBLE(cv) + DIMAG(wsum1)*DIMAG(cv)
                   dot21 = dot21 + DBLE( wsum2)* DBLE(cv) + DIMAG(wsum2)*DIMAG(cv)
                   dot31 = dot31 + DBLE( wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                   !
                   cv = evci * CMPLX( gk2)
                   dot22 = dot22 + DBLE( wsum2)* DBLE(cv) + DIMAG(wsum2)*DIMAG(cv) 
                   dot32 = dot32 + DBLE( wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                   ! 
                   cv =  evci * CMPLX( gk3 )
                   dot33 = dot33 + DBLE( wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                ENDDO
             ENDDO 
             !
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             !
             sigmanlc(:,1) = sigmanlc(:,1) -4._DP * wg(ibnd, ik) * [dot11, dot21, dot31]
             sigmanlc(:,2) = sigmanlc(:,2) -4._DP * wg(ibnd, ik) * [0._DP, dot22, dot32]
             sigmanlc(:,3) = sigmanlc(:,3) -4._DP * wg(ibnd, ik) * [0._DP, 0._DP, dot33]
             IF ( nproc > 1 ) THEN
                 CALL errore ('stres_us_gamma_gpu line 303', &
                       'unexpected error nproc be 1 with GPU acceleration', 100) 
                !CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
                !CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
                !CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
             ENDIF
          ENDDO
       ENDDO
       !
10     CONTINUE
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       ENDDO
       !
       !
       CALL dev_buf%release_buffer( deff_d, ierrs(1) )
       CALL dev_buf%release_buffer( ps_d,   ierrs(2) )
       CALL dev_buf%release_buffer( dvkb_d, ierrs(3) )
       !
       !
       RETURN
       !
     END SUBROUTINE stres_us_gamma_gpu
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE stres_us_k_gpu()
       !----------------------------------------------------------------------  
       !! nonlocal contribution to the stress - k-points version       
       !
#if defined(_OPENMP) && defined(__PGI)
       USE omp_lib
#endif
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: na, np, ibnd, ipol, jpol, l, i, nt, &
                   ikb, jkb, ih, jh, is, js, ijs, ishift, nh_np
       REAL(DP) :: fac, evps, dot11, dot21, dot31, dot22, dot32, dot33, aux
       COMPLEX(DP) :: qm1i, gk1, gk2, gk3, ps, ps_nc(2)
       COMPLEX(DP) :: cv, cv1, cv2, worksum, worksum1, worksum2, evci, evc1i, &
                      evc2i, ps1, ps2, ps1d1, ps1d2, ps1d3, ps2d1, ps2d2,     &
                      ps2d3, psd1, psd2, psd3
       !
       COMPLEX(DP), POINTER :: dvkb_d(:,:,:)
       COMPLEX(DP), POINTER :: deff_d(:,:,:), deff_nc_d(:,:,:,:)
       COMPLEX(DP), POINTER :: ps_d(:), ps_nc_d(:,:)
       COMPLEX(DP), POINTER :: becpk_d(:,:), becpnc_d(:,:,:)
       INTEGER :: nhmx, ierrs(3)
       !
       REAL(DP) :: xyz(3,3)
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1._DP, 0._DP, 0._DP, 0._DP, 1._DP, 0._DP, 0._DP, &
                  0._DP, 1._DP /
       !
#if defined(__CUDA) 
       ATTRIBUTES(DEVICE) :: ps_d, ps_nc_d, becpnc_d, becpk_d, dvkb_d, &
                             deff_d, deff_nc_d
#endif
       !
       ! WORKAROUND STARTS ==================================================
       !
       ! There seems to be a bug with the OpenMP code generated by PGI 18.5
       ! for this subroutine.
       !
#if defined(_OPENMP) && defined(__PGI)
       INTEGER :: num_threads
       num_threads=omp_get_max_threads()
       CALL omp_set_num_threads(1)
#endif
       ! WORKAROUND ENDS ====================================================
       !
       !
       !
       evps = 0._DP
       ! ... diagonal contribution
       !
       CALL dev_buf%lock_buffer( dvkb_d, (/ npwx,nkb,4 /), ierrs(1) )
       !
       CALL gen_us_dj_gpu( ik, dvkb_d(:,:,4) )
       IF ( lmaxkb > 0 ) THEN 
         DO ipol = 1, 3
           CALL gen_us_dy_gpu( ik, xyz(1,ipol), dvkb_d(:,:,ipol) )
         ENDDO
       ENDIF
       !
       IF (noncolin) THEN
          CALL dev_buf%lock_buffer( ps_nc_d, (/ nkb,npol /), ierrs(2) )
          CALL dev_buf%lock_buffer( deff_nc_d, (/ nhm,nhm,nat,nspin /), ierrs(3) )
          becpnc_d => becp_d%nc_d
       ELSE 
          CALL dev_buf%lock_buffer ( ps_d, nkb, ierrs(2) )
          CALL dev_buf%lock_buffer( deff_d, (/ nhm,nhm,nat /), ierrs(3) )
          becpk_d => becp_d%k_d
       ENDIF
       IF (ANY(ierrs /= 0)) CALL errore( 'stres_us_gpu', 'cannot allocate buffers', -1 )
       !
       CALL using_deeq_d(0)
       !
       CALL using_et(0)
       !
       ! ... the contribution is calculated only on one processor because
       ! ... partial results are later summed over all processors
       !
       IF ( me_bgrp /= root_bgrp ) GO TO 100
       !
       !
       DO ibnd = 1, nbnd
          fac = wg(ibnd,ik)
          IF (ABS(fac) < 1.d-9) CYCLE
          IF (noncolin) THEN
             CALL compute_deff_nc_gpu( deff_nc_d , et(ibnd,ik) )
          ELSE
             CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
          ENDIF
          !
          !
          IF (noncolin) THEN
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itot
              !
              ih = ix_d(i,2)     ; na = ix_d(i,1)
              ishift = ix_d(i,4) ; ikb = ishift + ih
              aux = 0.d0
              !
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 ijs = 0
                 DO is = 1, npol
                   DO js = 1, npol
                      ijs = ijs + 1
                      aux = aux + fac * DBLE(deff_nc_d(ih,ih,na,ijs) &
                                         * CONJG(becpnc_d(ikb,is,ibnd)) *   &
                                                 becpnc_d(ikb,js,ibnd))
                   ENDDO
                 ENDDO
                 !
              ELSE
                 !
                 nh_np = ix_d(i,3)
                 !
                 ijs = 0
                 DO is = 1, npol
                   DO js = 1, npol
                      ijs = ijs + 1
                      aux = aux + fac * DBLE(deff_nc_d(ih,ih,na,ijs) &
                                         * CONJG(becpnc_d(ikb,is,ibnd))*    &
                                                 becpnc_d(ikb,js,ibnd))
                   ENDDO
                 ENDDO
                 !
                 DO jh = ih+1, nh_np
                    jkb = ishift + jh
                    ijs = 0
                    DO is = 1, npol
                      DO js = 1, npol
                         ijs = ijs + 1
                         aux = aux + 2._DP*fac * &
                                       DBLE( deff_nc_d(ih,jh,na,ijs) * &
                                       (CONJG(becpnc_d(ikb,is,ibnd)) *   &
                                              becpnc_d(jkb,js,ibnd)) )
                      ENDDO
                    ENDDO
                 ENDDO
                 !
              ENDIF
              evps = evps + aux
            ENDDO
            !
          ELSE
            !
            aux = 0.d0
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itot
              !
              ih = ix_d(i,2)     ; na = ix_d(i,1)
              ishift = ix_d(i,4) ; ikb = ishift + ih
              !
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 aux = fac * deff_d(ih,ih,na) * &
                               ABS(becpk_d(ikb,ibnd) )**2
                 !
              ELSE
                 !
                 nh_np = ix_d(i,3)
                 !
                 aux = fac * DBLE(deff_d(ih,ih,na) * &
                                ABS(becpk_d(ikb,ibnd) )**2) + &
                                DBLE(SUM( deff_d(ih,ih+1:nh_np,na) * &
                                fac * 2._DP*DBLE( CONJG(becpk_d(ikb,ibnd)) &
                                * becpk_d(ishift+ih+1:ishift+nh_np,ibnd) ) ))
                 !
              ENDIF
              evps = evps + aux
            ENDDO
            !
          ENDIF
          !
       ENDDO
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       ENDDO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !
       DO ibnd = 1, nbnd
         !
         IF ( noncolin ) THEN
            !
            CALL compute_deff_nc_gpu( deff_nc_d, et(ibnd,ik) )
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itot
              !
              ih = ix_d(i,2)     ; na = ix_d(i,1)
              ishift = ix_d(i,4) ; ikb = ishift + ih
              !
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 DO is = 1, npol
                   ijs = (is-1)*npol
                   ps_nc_d(ikb,is) = SUM( becpnc_d(ikb,1:npol,ibnd) * &
                                          deff_nc_d(ih,ih,na,ijs+1:ijs+npol) )
                 ENDDO
                 !
              ELSE
                 !
                 nh_np = ix_d(i,3)
                 !
                 DO is = 1, npol
                   ijs = (is-1)*npol
                   ps_nc_d(ikb,is) = SUM( becpnc_d(ishift+1:ishift+nh_np,1:npol,ibnd) * &
                                          deff_nc_d(ih,1:nh_np,na,ijs+1:ijs+npol) )
                 ENDDO  
                 !
              ENDIF
            ENDDO
            !
         ELSE
            !
            CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itot
               !
               ih = ix_d(i,2)     ; na = ix_d(i,1)
               ishift = ix_d(i,4) ; ikb = ishift + ih
               !
               IF (.NOT. is_multinp_d(i)) THEN
                  ps_d(ikb) = CMPLX(deeq_d(ih,ih,na,current_spin)) * &
                                    becpk_d(ikb,ibnd)
               ELSE 
                  nh_np = ix_d(i,3)
                  !
                  ps_d(ikb) = SUM( becpk_d(ishift+1:ishift+nh_np,ibnd) * &
                                                 deff_d(ih,1:nh_np,na) )
               ENDIF
               !
            ENDDO
            !
         ENDIF
         !
         dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
         dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
         !
         !
         IF (noncolin) THEN
            !$cuf kernel do(2) <<<*,*>>>    
            DO ikb =1, nkb
               DO i = 1, npw    
                  evc1i = evc_d(i, ibnd)
                  evc2i = evc_d(i+npwx,ibnd)
                  qm1i = CMPLX(qm1_d(i))
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  worksum1 = ps_nc_d(ikb,1) * dvkb_d(i,ikb,4)     
                  worksum2 = ps_nc_d(ikb,2) * dvkb_d(i,ikb,4)   
                  !   
                  cv1 = evc1i * gk1 * gk1 * qm1i
                  cv2 = evc2i * gk1 * gk1 * qm1i
                  dot11 = dot11 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)
                  !   
                  cv1 = evc1i * gk2 * gk1 * qm1i
                  cv2 = evc2i * gk2 * gk1 * qm1i
                  dot21 = dot21 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk3 * gk1 * qm1i
                  cv2 = evc2i * gk3 * gk1 * qm1i
                  dot31 = dot31 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk2 * gk2 * qm1i
                  cv2 = evc2i * gk2 * gk2 * qm1i
                  dot22 = dot22 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk3 * gk2 * qm1i
                  cv2 = evc2i * gk3 * gk2 * qm1i
                  dot32 = dot32 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)
                  !   
                  cv1  = evc1i * gk3 * gk3 * qm1i
                  cv2  = evc2i * gk3 * gk3 * qm1i
                  dot33 = dot33 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)
               ENDDO
            ENDDO
            !
         ELSE   
            !
            !
            !$cuf kernel do(2) <<<*,*>>>   
            DO ikb = 1, nkb
               DO i = 1, npw
                  !
                  worksum = ps_d(ikb) *dvkb_d(i,ikb,4)   
                  !    
                  evci = evc_d(i,ibnd)
                  qm1i = CMPLX(qm1_d(i))
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  !
                  cv = evci * gk1 * gk1 * qm1i
                  dot11 = dot11 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk2 * gk1 * qm1i
                  dot21 = dot21 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk1 * qm1i
                  dot31 = dot31 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk2 * gk2 * qm1i
                  dot22 = dot22 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk2 * qm1i
                  dot32 = dot32 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk3 * qm1i
                  dot33 = dot33 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !   
               ENDDO   
            ENDDO   
            !   
            !IF ( me_bgrp /= root_bgrp ) stop     
            !   
         ENDIF
         !
         sigmanlc(:,1) = sigmanlc(:,1) - 2._DP * wg(ibnd,ik) * [dot11, dot21, dot31]
         sigmanlc(:,2) = sigmanlc(:,2) - 2._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
         sigmanlc(:,3) = sigmanlc(:,3) - 2._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]
         !      
         ! ... non diagonal contribution - derivative of the spherical harmonics
         ! ... (no contribution from l=0)
         !      
         IF ( lmaxkb == 0 ) CYCLE       
         !      
         dot11 = 0._DP ;  dot21 = 0._DP      
         dot31 = 0._DP ;  dot22 = 0._DP      
         dot32 = 0._DP ;  dot33 = 0._DP      
         !
         IF (noncolin) THEN      
            !      
            !$cuf kernel do(2) <<<*,*>>>      
            DO ikb =1, nkb
               DO i = 1, npw
                  !       
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  !
                  ps1 = ps_nc_d(ikb,1)
                  ps2 = ps_nc_d(ikb,2)
                  !
                  ps1d1 = ps1 * dvkb_d(i,ikb,1)
                  ps1d2 = ps1 * dvkb_d(i,ikb,2)       
                  ps1d3 = ps1 * dvkb_d(i,ikb,3)       
                  !
                  ps2d1 = ps2 * dvkb_d(i,ikb,1)
                  ps2d2 = ps2 * dvkb_d(i,ikb,2)
                  ps2d3 = ps2 * dvkb_d(i,ikb,3)
                  !
                  evc1i = evc_d(i,ibnd)       
                  evc2i = evc_d(i+npwx,ibnd)       
                  !      
                  cv1 = evc1i * gk1
                  cv2 = evc2i * gk1
                  dot11 = dot11 + DBLE(ps1d1)*DBLE(cv1)   + &
                                  DIMAG(ps1d1)*DIMAG(cv1) + &
                                  DBLE(ps2d1)*DBLE(cv2)   + &
                                  DIMAG(ps2d1)*DIMAG(cv2)
                  !
                  dot21 = dot21 + DBLE(ps1d2)*DBLE(cv1)   + &
                                  DIMAG(ps1d2)*DIMAG(cv1) + &
                                  DBLE(ps2d2)*DBLE(cv2)   + &
                                  DIMAG(ps2d2)*DIMAG(cv2)
                  !
                  dot31 = dot31 + DBLE(ps1d3)*DBLE(cv1)   + &
                                  DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2)   + &
                                  DIMAG(ps2d3)*DIMAG(cv2)       
                  !
                  cv1 = evc1i * gk2
                  cv2 = evc2i * gk2
                  dot22 = dot22 + DBLE(ps1d2)*DBLE(cv1)   + &
                                  DIMAG(ps1d2)*DIMAG(cv1) + &
                                  DBLE(ps2d2)*DBLE(cv2)   + &
                                  DIMAG(ps2d2)*DIMAG(cv2)
                  !
                  dot32 = dot32 + DBLE(ps1d3)*DBLE(cv1)   + &
                                  DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2)   + &
                                  DIMAG(ps2d3)*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk3
                  cv2 = evc2i * gk3
                  dot33 = dot33 + DBLE(ps1d3)*DBLE(cv1)   + &
                                  DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2)   + &
                                  DIMAG(ps2d3)*DIMAG(cv2)      
                  !
               ENDDO
            ENDDO
            !
         ELSE
            !
            !$cuf kernel do(2) <<<*,*>>>
            DO ikb = 1, nkb
               DO i = 1, npw
                 ps   = ps_d(ikb)
                 psd1 = ps*dvkb_d(i,ikb,1)
                 psd2 = ps*dvkb_d(i,ikb,2)       
                 psd3 = ps*dvkb_d(i,ikb,3)
                 evci = evc_d(i,ibnd)
                 gk1  = CMPLX(gk_d(i,1))
                 gk2  = CMPLX(gk_d(i,2))
                 gk3  = CMPLX(gk_d(i,3))
                 !
                 cv = evci * gk1
                 dot11 = dot11 + DBLE(psd1)*DBLE(cv) + DIMAG(psd1)*DIMAG(cv)
                 dot21 = dot21 + DBLE(psd2)*DBLE(cv) + DIMAG(psd2)*DIMAG(cv)
                 dot31 = dot31 + DBLE(psd3)*DBLE(cv) + DIMAG(psd3)*DIMAG(cv)
                 !
                 cv = evci * gk2
                 dot22 = dot22 + DBLE(psd2)*DBLE(cv) + DIMAG(psd2)*DIMAG(cv)
                 dot32 = dot32 + DBLE(psd3)*DBLE(cv) + DIMAG(psd3)*DIMAG(cv)
                 !
                 cv = evci * gk3
                 dot33 = dot33 + DBLE(psd3)*DBLE(cv) + DIMAG(psd3)*DIMAG(cv)
               ENDDO
            ENDDO
            !   
         ENDIF      
         !       
         sigmanlc(:,1) = sigmanlc(:,1) -2._DP * wg(ibnd,ik) * [dot11, dot21, dot31]
         sigmanlc(:,2) = sigmanlc(:,2) -2._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
         sigmanlc(:,3) = sigmanlc(:,3) -2._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]
         !
       ENDDO 
       !
10     CONTINUE
       !
       CALL dev_buf%release_buffer( dvkb_d, ierrs(1) )
       !
       IF (noncolin) THEN 
          CALL dev_buf%release_buffer( ps_nc_d, ierrs(2) )
          CALL dev_buf%release_buffer( deff_nc_d, ierrs(3) )
       ELSE 
          CALL dev_buf%release_buffer( ps_d, ierrs(2) )
          CALL dev_buf%release_buffer( deff_d, ierrs(3) )
       ENDIF
       !
       ! WORKAROUND STARTS ==================================================
       !
       ! ... and now restore the previous value.
       !
#if defined(_OPENMP) && defined(__PGI)
       CALL omp_set_num_threads(num_threads)
#endif
       ! WORKAROUND ENDS ====================================================
       !
       RETURN
       !
     END SUBROUTINE stres_us_k_gpu
     !
END SUBROUTINE stres_us_gpu
