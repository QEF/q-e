!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! Workaround for missing interface
#if ! defined (__CUDA)
#define tg_gather_gpu tg_gather
#endif
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma_gpu( lda, n, m, psi_d, v_d, hpsi_d )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - Gamma point.
  !
  USE parallel_include
  USE kinds,          ONLY : DP
  USE control_flags,  ONLY : many_fft
  USE mp_bands,       ONLY : me_bgrp
  USE fft_base,       ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines
#if defined(__CUDA)
  USE device_fbuff_m, ONLY : dev_buf
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi_d(lda,m)
  COMPLEX(DP), INTENT(inout):: hpsi_d(lda,m)
  REAL(DP),    INTENT(in)   :: v_d(dffts%nnr)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, v_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr, right_nnr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:), vpsi(:,:)
  ! ... Variables for task groups
  LOGICAL :: use_tg
  REAL(DP), POINTER :: tg_v_d(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: tg_v_d
  INTEGER :: dffts_nnr, idx, ebnd, brange
  INTEGER :: ierr, ioff
  ! ... Variables to handle batched FFT
  INTEGER :: group_size, pack_size, remainder, howmany, hm_vec(3)
  REAL(DP):: fac
  !
  CALL start_clock_gpu( 'vloc_psi' )
  !
  ALLOCATE( psi(lda,m) )
  !$acc data create( psi )
  !$acc kernels
  psi = psi_d
  !$acc end kernels
  !
  incr = 2*many_fft
  !
  use_tg = dffts%has_task_groups
  !
  IF ( use_tg ) THEN
     CALL start_clock_gpu( 'vloc_psi:tg_gather' )
     dffts_nnr = dffts%nnr_tg
     incr = 2*fftx_ntgrp(dffts)
     CALL dev_buf%lock_buffer( tg_v_d, dffts_nnr, ierr )
     ALLOCATE( tg_psic(dffts_nnr), tg_vpsi(dffts_nnr,incr) )
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock_gpu( 'vloc_psi:tg_gather' )
  ELSE
     dffts_nnr = dffts%nnr
     ALLOCATE( psic(dffts_nnr*incr), vpsi(dffts_nnr,incr) )
  ENDIF
  !
  ! ... The local potential V_Loc psi:
  !    - fft to real space;
  !    - product with the potential v on the smooth grid;
  !    - back to reciprocal space.
  !
  IF( use_tg ) THEN
     !
     !$acc data create(tg_psic,tg_vpsi)
     DO ibnd = 1, m, incr
        !
        CALL tgwave_g2r( psi(1:n,ibnd:m), tg_psic, dffts, n )
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        !$acc parallel loop
        DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
           tg_psic(j) = tg_psic(j) * tg_v_d(j)
        ENDDO
        !
        brange = m-ibnd+1
        !
        CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n )
        !
        DO idx = 1, 2*fftx_ntgrp(dffts), 2
           IF ( idx+ibnd-1<m ) THEN
              !$acc parallel loop
              DO j = 1, n
                 hpsi_d(j,ibnd+idx-1) = hpsi_d(j,ibnd+idx-1) + 0.5d0 * tg_vpsi(j,idx)
                 hpsi_d(j,ibnd+idx) = hpsi_d(j,ibnd+idx) + 0.5d0 * tg_vpsi(j,idx+1)
              ENDDO
           ELSEIF ( idx+ibnd-1==m ) THEN
              !$acc parallel loop
              DO j = 1, n
                 hpsi_d(j,ibnd+idx-1) = hpsi_d(j,ibnd+idx-1) + tg_vpsi(j,idx)
              ENDDO
           ENDIF
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ELSEIF (many_fft > 1) THEN
     !
     !$acc data create(psic,vpsi)
     DO ibnd = 1, m, incr
        !
        group_size = MIN(2*many_fft, m-(ibnd-1))
        pack_size = (group_size/2) ! This is FLOOR(group_size/2)
        remainder = group_size - 2*pack_size
        howmany = pack_size + remainder
        hm_vec(1)=group_size ; hm_vec(2)=n ; hm_vec(3)=howmany
        !
        CALL wave_g2r( psi(:,ibnd:ibnd+group_size-1), psic, dffts, howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, howmany-1
          DO j = 1, dffts_nnr
            psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v_d(j)
          ENDDO
        ENDDO
        !
        CALL wave_r2g( psic, vpsi, dffts, howmany_set=hm_vec )
        !
        IF ( pack_size > 0 ) THEN
           !$acc parallel loop collapse(2)
           DO idx = 0, pack_size-1
              DO j = 1, n
                 hpsi_d(j,ibnd+idx*2)   = hpsi_d(j,ibnd+idx*2)   + vpsi(j,idx*2+1)
                 hpsi_d(j,ibnd+idx*2+1) = hpsi_d(j,ibnd+idx*2+1) + vpsi(j,idx*2+2)
              ENDDO
           ENDDO
        ENDIF
        !
        IF (remainder > 0) THEN
           !$acc parallel loop
           DO j = 1, n
              hpsi_d(j,ibnd+group_size-1) = hpsi_d(j,ibnd+group_size-1) + &
                                            vpsi(j,group_size)
           ENDDO
        ENDIF
        !
     ENDDO
     !$acc end data
     !
  ELSE
     !
     !$acc data create(psic,vpsi)
     DO ibnd = 1, m, incr
        !
        ebnd = ibnd
        IF ( ibnd<m ) ebnd = ibnd + 1
        !
        CALL wave_g2r( psi(1:n,ibnd:ebnd), psic, dffts )
        !        
        !$acc parallel loop
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v_d(j)
        ENDDO
        !
        brange=1 ;  fac=1.d0
        IF ( ibnd<m ) THEN
          brange=2 ;  fac=0.5d0
        ENDIF
        !
        CALL wave_r2g( psic, vpsi(:,1:brange), dffts )
        !
        !$acc parallel loop
        DO j = 1, n
          hpsi_d(j,ibnd) = hpsi_d(j,ibnd) + fac*vpsi(j,1)
          IF ( ibnd<m ) hpsi_d(j,ibnd+1) = hpsi_d(j,ibnd+1) + fac*vpsi(j,2)
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psi )
  !
  IF( use_tg ) THEN
     DEALLOCATE( tg_psic, tg_vpsi )
     CALL dev_buf%release_buffer( tg_v_d, ierr )
  ELSE
     DEALLOCATE( psic, vpsi )
  ENDIF
  !
  CALL stop_clock_gpu ('vloc_psi')
#endif
  !
  RETURN
END SUBROUTINE vloc_psi_gamma_gpu
!
!@njs: vloc_psi_k
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k_gpu( lda, n, m, psi_d, v_d, hpsi_d )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points. GPU double.
  !
  !   fft to real space
  !   product with the potential v on the smooth grid
  !   back to reciprocal space
  !   addition to the hpsi
  !
  USE parallel_include
  USE kinds,         ONLY : DP
  USE wvfct,         ONLY : current_k
  USE klist,         ONLY : igk_k
  USE mp_bands,      ONLY : me_bgrp
  USE control_flags, ONLY : many_fft
  USE fft_base,      ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines
#if defined(__CUDA)
  USE device_fbuff_m,    ONLY : dev_buf
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi_d(lda,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi_d(lda,m)
  REAL(DP), INTENT(IN) :: v_d(dffts%nnr)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, v_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ebnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:), vpsi(:,:)
  ! ... Task Groups
  LOGICAL :: use_tg
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:)
  REAL(DP), POINTER :: tg_v_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: tg_v_d
  !
  INTEGER :: dffts_nnr, idx, group_size, hm_vec(3)
  INTEGER :: ierr, brange
  !
  CALL start_clock_gpu ('vloc_psi')
  use_tg = dffts%has_task_groups
  !
  ALLOCATE( psi(lda,m) )
  !$acc data create( psi )
  !$acc kernels
  psi = psi_d
  !$acc end kernels
  !
  incr = many_fft
  !
  IF( use_tg ) THEN
     CALL start_clock_gpu ('vloc_psi:tg_gather')
     dffts_nnr =  dffts%nnr_tg
     incr = fftx_ntgrp(dffts)
     CALL dev_buf%lock_buffer( tg_v_d, dffts_nnr, ierr )
     ALLOCATE( tg_psic(dffts_nnr), tg_vpsi(dffts_nnr,incr) )
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock_gpu ('vloc_psi:tg_gather')
  ELSE
     dffts_nnr = dffts%nnr
     ALLOCATE( psic(dffts_nnr*incr), vpsi(dffts_nnr,incr) )
  ENDIF
  !
  IF( use_tg ) THEN
     !
     CALL tg_get_nnr( dffts, right_nnr )
     !
     !$acc data create(tg_psic,tg_vpsi)
     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !
        CALL tgwave_g2r( psi(1:n,ibnd:m), tg_psic, dffts, n, igk_k(:,current_k) )
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        !$acc parallel loop
        DO j = 1, dffts%nr1x*dffts%nr2x* right_nr3
           tg_psic(j) = tg_psic(j) * tg_v_d(j)
        ENDDO
        !
        brange = m-ibnd+1
        !
        CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n, igk_k(:,current_k) )
        !
        DO idx = 1, fftx_ntgrp(dffts)
           IF ( idx+ibnd-1 <= m ) THEN
              !$acc parallel loop
              DO j = 1, n
                 hpsi_d(j,ibnd+idx-1) = hpsi_d(j,ibnd+idx-1) + tg_vpsi(j,idx-1)
              ENDDO
           ENDIF
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ELSEIF (many_fft > 1) THEN
     !
     !$acc data create(psic,vpsi)
     DO ibnd = 1, m, incr
        !
        group_size = MIN(many_fft,m-(ibnd-1))
        hm_vec(1)=group_size ; hm_vec(2)=n ; hm_vec(3)=group_size
        ebnd = ibnd+group_size-1
        !
        CALL wave_g2r( psi(:,ibnd:ebnd), psic, dffts, igk=igk_k(:,current_k), &
                       howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO j = 1, dffts_nnr
              psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v_d(j)
           ENDDO
        ENDDO
        !
        CALL wave_r2g( psic, vpsi, dffts, igk=igk_k(:,current_k), &
                       howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO j = 1, n
              hpsi_d(j,ibnd+idx) = hpsi_d(j,ibnd+idx) + vpsi(j,idx+1)
           ENDDO
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ELSE
     !
     !$acc data create(psic,vpsi)
     DO ibnd = 1, m
        !
        CALL wave_g2r( psi(1:n,ibnd:ibnd), psic, dffts, igk=igk_k(:,current_k) )
        !
        !$acc parallel loop
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v_d(j)
        ENDDO
        !
        CALL wave_r2g( psic, vpsi(1:n,:), dffts, igk=igk_k(:,current_k) )
        !
        !$acc parallel loop
        DO i = 1, n
           hpsi_d(i,ibnd) = hpsi_d(i,ibnd) + vpsi(i,1)
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psi )
  !
  IF( use_tg ) THEN
     DEALLOCATE( tg_psic, tg_vpsi )
     CALL dev_buf%release_buffer( tg_v_d, ierr )
  ELSE
     DEALLOCATE( psic, vpsi )
  ENDIF
  !
  CALL stop_clock_gpu( 'vloc_psi' )
#endif
  !
99 format ( 20 ('(',2f12.9,')') )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_k_gpu
!
!@njs: vloc_psi_nc
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc_gpu( lda, n, m, psi_d, v_d, hpsi_d )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - non-collinear - 
  !! GPU version.
  !
  USE parallel_include
  USE kinds,               ONLY : DP
  USE wvfct,               ONLY : current_k
  USE klist,               ONLY : igk_k
  USE mp_bands,            ONLY : me_bgrp
  USE fft_base,            ONLY : dffts, dfftp
  USE fft_wave
  USE lsda_mod,            ONLY : nspin
  USE noncollin_module,    ONLY : npol, domag
  USE fft_helper_subroutines
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  REAL(DP), INTENT(IN) :: v_d(dfftp%nnr,4) ! beware dimensions!
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi_d(lda,npol,m)
#if defined(__CUDA)
  attributes(DEVICE) :: v_d, psi_d, hpsi_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), psic_nc(:,:), vpsic_nc(:,:)
  ! ... Variables for task groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:), tg_vpsi(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: tg_v_d
  INTEGER :: dffts_nnr, idx, ioff, ii, ie, brange
  INTEGER :: right_nnr, right_nr3, right_inc
  !
  CALL start_clock_gpu ('vloc_psi')
  !
  incr = 1
  use_tg = dffts%has_task_groups 
  !
  ALLOCATE( psi(lda*npol,m) )
  !$acc data create( psi )
  !$acc kernels
  psi = psi_d
  !$acc end kernels
  !
  IF( use_tg ) THEN
     CALL start_clock_gpu( 'vloc_psi:tg_gather' )
     dffts_nnr = dffts%nnr_tg
     incr = fftx_ntgrp(dffts)
     IF (domag) THEN
        ALLOCATE( tg_v_d(dffts_nnr,4) )
        DO is = 1, nspin
           CALL tg_gather_gpu( dffts, v_d(:,is), tg_v_d(:,is) )
        ENDDO
     ELSE
        ALLOCATE( tg_v_d(dffts_nnr,1) )
        CALL tg_gather_gpu( dffts, v_d(:,1), tg_v_d(:,1) )
     ENDIF
     ALLOCATE( tg_psic(dffts_nnr,npol), tg_vpsi(lda,incr) )
     CALL stop_clock_gpu( 'vloc_psi:tg_gather' )
  ELSE
     dffts_nnr = dffts%nnr
     ALLOCATE( psic_nc(dffts_nnr,npol), vpsic_nc(lda,1) )
  ENDIF
  !
  ! ... the local potential V_Loc psi. First the psi in real space
  !
  IF( use_tg ) THEN
     !
     !$acc data create( tg_psic, tg_vpsi )
     DO ibnd = 1, m, incr
        !
        DO ipol = 1, npol
           ii = lda*(ipol-1)+1
           ie = lda*ipol
           CALL tgwave_g2r( psi(ii:ie,ibnd:m), tg_psic(:,ipol), dffts, n, &
                            igk_k(:,current_k) )
        ENDDO
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        IF (domag) THEN
           !$acc parallel loop
           DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
              sup  = tg_psic(j,1) * (tg_v_d(j,1)+tg_v_d(j,4)) + &
                     tg_psic(j,2) * (tg_v_d(j,2)-(0.d0,1.d0)*tg_v_d(j,3))
              sdwn = tg_psic(j,2) * (tg_v_d(j,1)-tg_v_d(j,4)) + &
                     tg_psic(j,1) * (tg_v_d(j,2)+(0.d0,1.d0)*tg_v_d(j,3))
              tg_psic(j,1) = sup
              tg_psic(j,2) = sdwn
           ENDDO
        ELSE
           !$acc parallel loop
           DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
              tg_psic(j,:) = tg_psic(j,:) * tg_v_d(j,1)
           ENDDO
        ENDIF
        !
        brange = m-ibnd+1
        !
        DO ipol = 1, npol
           CALL tgwave_r2g( tg_psic(:,ipol), tg_vpsi(:,1:brange), dffts, n, igk_k(:,current_k) )
           !
           CALL tg_get_recip_inc( dffts, right_inc )
           !
           ioff = 0
           !
           DO idx = 1, fftx_ntgrp(dffts)
             IF ( idx+ibnd-1<=m ) THEN
               !$acc parallel loop
               DO j = 1, n
                 hpsi_d(j,ipol,ibnd+idx-1) = hpsi_d(j,ipol,ibnd+idx-1) + tg_vpsi(j,idx)
               ENDDO
             ENDIF
             ioff = ioff + right_inc
           ENDDO
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ELSE
     !
     !$acc data create( psic_nc, vpsic_nc )
     DO ibnd = 1, m, incr
        !
        !$acc kernels
        psic_nc = (0.d0,0.d0)
        !$acc end kernels
        DO ipol = 1, npol
           ii = lda*(ipol-1)+1
           ie = lda*(ipol-1)+n
           CALL wave_g2r( psi(ii:ie,ibnd:ibnd), psic_nc(:,ipol), dffts, &
                          igk=igk_k(:,current_k) )
        ENDDO
        !
        IF (domag) THEN
           !$acc parallel loop
           DO j = 1, dffts_nnr
              sup  = psic_nc(j,1) * (v_d(j,1)+v_d(j,4)) + &
                     psic_nc(j,2) * (v_d(j,2)-(0.d0,1.d0)*v_d(j,3))
              sdwn = psic_nc(j,2) * (v_d(j,1)-v_d(j,4)) + &
                     psic_nc(j,1) * (v_d(j,2)+(0.d0,1.d0)*v_d(j,3))
              psic_nc(j,1) = sup
              psic_nc(j,2) = sdwn
           ENDDO
        ELSE
           !$acc parallel loop
           DO j = 1, dffts_nnr
              psic_nc(j,:) = psic_nc(j,:) * v_d(j,1)
           ENDDO
        ENDIF
        !
        DO ipol = 1, npol
           CALL wave_r2g( psic_nc(:,ipol), vpsic_nc(1:n,:), dffts, &
                          igk=igk_k(:,current_k) )
           !
           !$acc parallel loop
           DO j = 1, n
              hpsi_d(j,ipol,ibnd) = hpsi_d(j,ipol,ibnd) + vpsic_nc(j,1)
           ENDDO
        ENDDO
        !
     ENDDO
     !$acc end data
     !
  ENDIF
  !
  IF( use_tg ) THEN
     DEALLOCATE( tg_v_d, tg_psic, tg_vpsi )
  ELSE
     DEALLOCATE( psic_nc, vpsic_nc )
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psi )
  !
  CALL stop_clock_gpu ('vloc_psi')
#endif
  !
  RETURN
  !
END SUBROUTINE vloc_psi_nc_gpu

