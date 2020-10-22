!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE stres_knl_gpu( sigmanlc, sigmakin )
  !-----------------------------------------------------------------------
  !! Computes the kinetic + nonlocal contribuition to the stress
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi, e2
  USE cell_base,            ONLY: omega, alat, at, bg, tpiba
  USE gvect,                ONLY: g
  USE gvecw,                ONLY: qcutz, ecfixed, q2sigma
  USE klist,                ONLY: nks, xk, ngk, igk_k_d
  USE io_files,             ONLY: iunwfc, nwordwfc
  USE buffers,              ONLY: get_buffer
  USE symme,                ONLY: symmatrix
  USE wvfct,                ONLY: npwx, nbnd, wg
  USE control_flags,        ONLY: gamma_only, use_gpu
  USE noncollin_module,     ONLY: noncolin, npol
  USE wavefunctions,        ONLY: evc
  USE mp_pools,             ONLY: inter_pool_comm
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp,                   ONLY: mp_sum
  USE wavefunctions_gpum,   ONLY: using_evc, using_evc_d, evc_d  
  USE device_fbuff_m,             ONLY : dev_buf
  USE gvect_gpum,           ONLY : g_d
  USE device_memcpy_m,        ONLY : dev_memcpy
 
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sigmanlc(3,3)
  !! non-local contribution to stress
  REAL(DP) :: sigmakin(3,3)
  !! kinetic contribution to stress
  !
  ! ... local variables
  !
  REAL(DP), POINTER :: gk_d(:,:), kfac_d (:)
  REAL(DP) :: twobysqrtpi, gk2, arg, s11, s21, s31, s22, s32, s33, &
              xk1, xk2, xk3, tmpf, wg_nk 
  INTEGER  :: npw, ik, l, m, i, ibnd, is
  INTEGER  :: ierr(2)
  !
#if defined(__CUDA) 
  ATTRIBUTES(DEVICE)  :: gk_d, kfac_d 
#endif   
  !
  CALL using_evc(0)
  CALL using_evc_d(0)
  !
  CALL dev_buf%lock_buffer( gk_d, [npwx,3], ierr(1) )
  CALL dev_buf%lock_buffer( kfac_d, npwx, ierr(2)  )
  !
  sigmanlc(:,:) = 0._DP
  sigmakin(:,:) = 0._DP
  twobysqrtpi   = 2._DP/SQRT(pi)
  !
  kfac_d(:) = 1._DP
  !
  s11 = 0._DP ;  s22 = 0._DP
  s21 = 0._DP ;  s32 = 0._DP
  s31 = 0._DP ;  s33 = 0._DP
  !
  DO ik = 1, nks
     IF ( nks > 1 ) THEN
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        CALL using_evc(2)
        CALL using_evc_d(0)
     ENDIF 
     !
     npw = ngk(ik)
     !
     xk1 = xk(1,ik)
     xk2 = xk(2,ik)
     xk3 = xk(3,ik)
     !
     !$cuf kernel do(1) <<<*,*>>>
     DO i = 1, npw
        gk_d(i,1) = ( xk1 + g_d(1,igk_k_d(i,ik)) ) * tpiba
        gk_d(i,2) = ( xk2 + g_d(2,igk_k_d(i,ik)) ) * tpiba
        gk_d(i,3) = ( xk3 + g_d(3,igk_k_d(i,ik)) ) * tpiba
        IF (qcutz > 0._DP) THEN
           gk2 = gk_d(i,1)**2 + gk_d(i,2)**2 + gk_d(i,3)**2
           arg = ( (gk2-ecfixed)/q2sigma )**2
           kfac_d(i) = 1._DP + qcutz / q2sigma * twobysqrtpi * EXP(-arg)
        ENDIF
     ENDDO
     !
     ! ... kinetic contribution
     !
     DO ibnd = 1, nbnd
        wg_nk = wg(ibnd,ik)
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, npw
           IF (noncolin) THEN
              tmpf = wg_nk * kfac_d(i) * &
               ( DBLE(CONJG(evc_d(   i  ,ibnd))*evc_d(   i  ,ibnd)) + &
                 DBLE(CONJG(evc_d(i+npwx,ibnd))*evc_d(i+npwx,ibnd)) )
           ELSE
              tmpf = wg_nk * kfac_d(i) * &
                    DBLE( CONJG(evc_d(i, ibnd) ) * evc_d(i, ibnd) )
           ENDIF
           s11 = s11 + tmpf * gk_d(i,1)*gk_d(i,1)
           s21 = s21 + tmpf * gk_d(i,2)*gk_d(i,1)
           s31 = s31 + tmpf * gk_d(i,3)*gk_d(i,1)
           s22 = s22 + tmpf * gk_d(i,2)*gk_d(i,2)
           s32 = s32 + tmpf * gk_d(i,3)*gk_d(i,2)
           s33 = s33 + tmpf * gk_d(i,3)*gk_d(i,3)
        ENDDO
     ENDDO
     !
     !  ... contribution from the  nonlocal part
     !
     CALL stres_us_gpu( ik, gk_d, sigmanlc )
     !
  ENDDO
  !
  sigmakin(:,1) =  [s11,  s21,  s31]
  sigmakin(:,2) =  [0._DP,s22,  s32]
  sigmakin(:,3) =  [0._DP,0._DP,s33]
  !
  CALL dev_buf%release_buffer( kfac_d, ierr(2)  )
  CALL dev_buf%release_buffer( gk_d,   ierr(1)  )
  !
  ! ... the kinetic term must be summed over PW's and over k-points
  !
  CALL mp_sum( sigmakin, intra_bgrp_comm )
  CALL mp_sum( sigmakin, inter_pool_comm )
  !
  ! ... the nonlocal term is summed here only over k-points, because we add
  ! ... to it the US term from augmentation charge derivatives
  !
  CALL mp_sum( sigmanlc, inter_pool_comm )
  !
  ! ... add US term from augmentation charge derivatives, sum result over PW's
  !
  CALL addusstress( sigmanlc )
  CALL mp_sum( sigmanlc, intra_bgrp_comm )
  !
  DO l = 1, 3
     DO m = 1, l-1
        sigmanlc(m,l) = sigmanlc(l,m)
        sigmakin(m,l) = sigmakin(l,m)
     ENDDO
  ENDDO
  !
  IF ( gamma_only ) THEN
     sigmakin(:,:) = 2._DP * e2 / omega * sigmakin(:,:)
  ELSE
     sigmakin(:,:) = e2 / omega * sigmakin(:,:)
  ENDIF
  sigmanlc(:,:) = -1._DP / omega * sigmanlc(:,:)
  !
  ! ... symmetrize stress
  !
  CALL symmatrix( sigmakin )
  CALL symmatrix( sigmanlc )
  !
  RETURN
  !
END SUBROUTINE stres_knl_gpu 

