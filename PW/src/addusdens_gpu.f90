!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusdens_gpu(rho)
  !----------------------------------------------------------------------
  !! Add US contribution to the charge density to \(\text{rho}(G)\).
  !
  USE realus,               ONLY : addusdens_r
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho(dfftp%ngm,nspin_mag)
  !! Charge density in G space
  !
  IF ( tqr ) THEN
     CALL addusdens_r( rho )
  ELSE
     CALL addusdens_g_gpu( rho )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE addusdens_gpu
!
!
!----------------------------------------------------------------------
SUBROUTINE addusdens_g_gpu(rho)
  !----------------------------------------------------------------------
  !! This routine adds to the charge density \(\text{rho}(G)\) in reciprocal space
  !! the part which is due to the US augmentation.
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#define cublasDgemm Dgemm
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, eigts1, eigts2, eigts3, mill
  USE gvect_gpum,           ONLY : gg_d, g_d, &
                                   eigts1_d, eigts2_d, eigts3_d, mill_d
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  USE uspp_gpum,            ONLY : becsum_d, using_becsum_d
  USE device_fbuff_m,             ONLY : dev_buf, pin_buf
  USE device_memcpy_m,        ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: rho(dfftp%ngm,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ngm_s, ngm_e, ngm_l, ngm_s_tmp, ngm_e_tmp, ngm_l_tmp
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ig, na, nt, ih, jh, ijh, is, nab, nb, nij
  ! counters

  REAL(DP), POINTER :: tbecsum_d(:,:,:)
  ! \sum_kv <\psi_kv|\beta_l><beta_m|\psi_kv> for each species of atoms
  REAL(DP), POINTER :: qmod_d (:), ylmk0_d (:,:)
  ! modulus of G, spherical harmonics
  COMPLEX(DP), POINTER :: skk_d(:,:)
  COMPLEX(DP), POINTER :: aux2_d(:,:)
  ! structure factors, US contribution to rho
  COMPLEX(DP), POINTER ::  aux_d (:,:), qgm_d(:)
  COMPLEX(DP), POINTER ::  aux_h (:,:)
  ! work space for rho(G,nspin), Fourier transform of q
  INTEGER :: ij, im, ierr
#if defined(__CUDA)
  attributes(device) :: tbecsum_d, qmod_d, ylmk0_d, skk_d, &
                        aux2_d, aux_d, qgm_d
  attributes(pinned) :: aux_h
#endif
  IF (.not.okvan) RETURN

  CALL start_clock_gpu ('addusdens')
  !
  CALL pin_buf%lock_buffer(aux_h, (/ ngm, nspin_mag /), ierr ) !ALLOCATE (aux_h (ngm, nspin_mag) )
  CALL dev_buf%lock_buffer(aux_d, (/ ngm, nspin_mag /), ierr ) !ALLOCATE (aux_d (ngm, nspin_mag) )
  IF( ierr /= 0 ) &
     CALL errore( ' addusdens_gpu ',' cannot allocate aux_d ', ABS(ierr) )

  !
  CALL dev_memset(aux_d, (0.d0, 0.d0), [ 1, ngm ], 1, [ 1, nspin_mag ], 1)
  !
  ! With k-point/bgrp parallelization, distribute G-vectors across all processors
  ! ngm_s = index of first G-vector for this processor (in the k-point x bgrp pool)
  ! ngm_e = index of last  G-vector for this processor (in the k-point x bgrp pool)
  ! ngm_l = local number of G-vectors 
  !
  CALL divide( inter_pool_comm, ngm, ngm_s_tmp, ngm_e_tmp ) ; ngm_l_tmp = ngm_e_tmp - ngm_s_tmp + 1
  CALL divide( inter_bgrp_comm, ngm_l_tmp, ngm_s, ngm_e ) ; ngm_l = ngm_e - ngm_s + 1 
  ngm_s = ngm_s + ngm_s_tmp - 1 ; ngm_e = ngm_e + ngm_s_tmp -1
  ! for the extraordinary unlikely case of more processors than G-vectors
  IF ( ngm_l <= 0 ) GO TO 10
  !
  ! Sync becsum if needed
  CALL using_becsum_d(0)
  !
  !ALLOCATE (qmod_d(ngm_l), qgm_d(ngm_l) )
  !ALLOCATE (ylmk0_d(ngm_l, lmaxq * lmaxq) )
  CALL dev_buf%lock_buffer(qmod_d, ngm_l, ierr )
  CALL dev_buf%lock_buffer(qgm_d, ngm_l, ierr )
  CALL dev_buf%lock_buffer(ylmk0_d, (/ ngm_l, lmaxq * lmaxq /), ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' addusdens_gpu ',' cannot allocate ylmk0_d ', ABS(ierr) )

  CALL ylmr2_gpu (lmaxq * lmaxq, ngm_l, g_d(1,ngm_s), gg_d(ngm_s), ylmk0_d)
  
!$cuf kernel do(1) <<<*,*>>>
  DO ig = 1, ngm_l
     qmod_d (ig) = SQRT(gg_d(ngm_s+ig-1))*tpiba
  ENDDO
  !
  ! Use largest size for buffer
  nij = nhm*(nhm+1)/2
  CALL dev_buf%prepare_buffer(aux2_d, (/ ngm_l,nij /), ierr )
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ! count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        !
        !ALLOCATE ( skk_d(ngm_l,nab), tbecsum_d(nij,nab,nspin_mag), aux2_d(ngm_l,nij) )
        CALL dev_buf%lock_buffer(skk_d, (/ ngm_l,nab /), ierr)
        CALL dev_buf%lock_buffer(tbecsum_d, (/ nij,nab,nspin_mag /), ierr )
        CALL dev_buf%lock_buffer(aux2_d, (/ ngm_l,nij /), ierr )
        IF( ierr /= 0 ) &
            CALL errore( ' addusdens_gpu ',' cannot allocate aux2_d ', ABS(ierr) )
        !
        call start_clock_gpu( 'addusd:skk')
        nb = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) THEN
              nb = nb + 1
              !tbecsum(:,nb,:) = becsum(1:nij,na,1:nspin_mag)
!$cuf kernel do(2) <<<*,*>>>
              DO im = 1, nspin_mag
                 DO ij = 1, nij
                   tbecsum_d(ij,nb,im) = becsum_d(ij,na,im)
                 ENDDO
              ENDDO

!$cuf kernel do(1) <<<*,*>>>
              DO ig = 1, ngm_l
                 skk_d(ig,nb) = eigts1_d (mill_d (1,ngm_s+ig-1), na) * &
                              eigts2_d (mill_d (2,ngm_s+ig-1), na) * &
                              eigts3_d (mill_d (3,ngm_s+ig-1), na)
              ENDDO
           ENDIF
        ENDDO
        call stop_clock_gpu( 'addusd:skk')

        DO is = 1, nspin_mag
           ! sum over atoms
           CALL cublasDgemm( 'N', 'T', 2*ngm_l, nij, nab, 1.0_dp, skk_d, 2*ngm_l,&
                tbecsum_d(1,1,is), nij, 0.0_dp, aux2_d, 2*ngm_l )
           ! sum over lm indices of Q_{lm}
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL qvan2_gpu (ngm_l, ih, jh, nt, qmod_d, qgm_d, ylmk0_d)
!$cuf kernel do(1) <<<*,*>>>
                 DO ig = 1, ngm_l
                    aux_d(ngm_s+ig-1,is) = aux_d(ngm_s+ig-1,is)+aux2_d(ig,ijh)*qgm_d(ig)
                 ENDDO

             ENDDO
           ENDDO
        ENDDO
        !DEALLOCATE (aux2_d, tbecsum_d, skk_d )
        CALL dev_buf%release_buffer(skk_d, ierr)
        CALL dev_buf%release_buffer(tbecsum_d, ierr)
        CALL dev_buf%release_buffer(aux2_d, ierr)
     ENDIF
  ENDDO
  !
  !DEALLOCATE (ylmk0_d)
  !DEALLOCATE (qgm_d, qmod_d)
  CALL dev_buf%release_buffer(ylmk0_d, ierr)
  CALL dev_buf%release_buffer(qgm_d, ierr)
  CALL dev_buf%release_buffer(qmod_d, ierr)
  !
  10 CONTINUE
  !
  CALL dev_memcpy(aux_h, aux_d)
  CALL mp_sum( aux_h, inter_bgrp_comm )
  CALL mp_sum( aux_h, inter_pool_comm )
  !
  !     add aux to the charge density in reciprocal space
  !
  rho(:,:) = rho(:,:) + aux_h (:,:)
  !
  !DEALLOCATE (aux_h, aux_d)
  CALL pin_buf%release_buffer(aux_h, ierr)
  CALL dev_buf%release_buffer(aux_d, ierr)
  !
  CALL stop_clock_gpu( 'addusdens' )
  !
  RETURN
END SUBROUTINE addusdens_g_gpu

