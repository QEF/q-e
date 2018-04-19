!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
#if defined(__CUDA)
!@njs: vloc_psi_gamma, psi, v, hpsi, psic
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma_gpu(lda, n, m, psi_d, v_d, hpsi_d)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - Gamma point
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  USE wavefunctions_module,       ONLY: psic_h => psic
  USE wavefunctions_module_gpum,       ONLY: psic_d
  USE fft_helper_subroutines
  USE buffers_module,       ONLY : buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(in)   :: psi_d (lda, m)
  COMPLEX(DP), DEVICE, INTENT(inout):: hpsi_d (lda, m)
  REAL(DP), DEVICE, INTENT(in) :: v_d(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr, right_nnr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  !
  LOGICAL :: use_tg
  ! Variables for task groups
!@njs: tg_v, tg_psic
  REAL(DP),    DEVICE, POINTER :: tg_v_d(:)
  COMPLEX(DP), DEVICE, POINTER :: tg_psic_d(:)
  INTEGER,     DEVICE, POINTER :: dffts_nl_d(:), dffts_nlm_d(:)
  INTEGER :: v_siz, idx, ioff
  INTEGER :: ierr
  !
  CALL start_clock ('vloc_psi')
  incr = 2
  !
  IF (.not. allocated(psic_d)) ALLOCATE(psic_d, SOURCE=psic_h)
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dffts%nnr_tg
     !
     CALL buffer%lock_buffer( tg_v_d, v_siz, ierr )
     CALL buffer%lock_buffer( tg_psic_d, v_siz, ierr )
     !
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
     incr = 2 * fftx_ntgrp(dffts)
     !
  ENDIF
  ! Sync fft data
  CALL buffer%lock_buffer( dffts_nl_d, size(dffts%nl), ierr )
  CALL buffer%lock_buffer( dffts_nlm_d, size(dffts%nlm), ierr )
  dffts_nl_d(1:n) = dffts%nl(1:n)
  dffts_nlm_d(1:n) = dffts%nlm(1:n)
  ! End Sync
  
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        !
        CALL tg_get_nnr( dffts, right_nnr )
        !
        tg_psic_d = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*fftx_ntgrp(dffts), 2
           IF( idx + ibnd - 1 < m ) THEN
              !$cuf kernel do(1) <<<,>>>
              DO j = 1, n
                 tg_psic_d(dffts_nl_d (j)+ioff) =        psi_d(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi_d(j,idx+ibnd)
                 tg_psic_d(dffts_nlm_d(j)+ioff) = conjg( psi_d(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi_d(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              !$cuf kernel do(1) <<<,>>>
              DO j = 1, n
                 tg_psic_d(dffts_nl_d (j)+ioff) =        psi_d(j,idx+ibnd-1)
                 tg_psic_d(dffts_nlm_d(j)+ioff) = conjg( psi_d(j,idx+ibnd-1) )
              ENDDO
           ENDIF

           ioff = ioff + right_nnr

        ENDDO
        !
     ELSE
        !
        psic_d(:) = (0.d0, 0.d0)
        IF (ibnd < m) THEN
           ! two ffts at the same time
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              psic_d(dffts_nl_d (j))=      psi_d(j,ibnd) + (0.0d0,1.d0)*psi_d(j,ibnd+1)
              psic_d(dffts_nlm_d(j))=conjg(psi_d(j,ibnd) - (0.0d0,1.d0)*psi_d(j,ibnd+1))
           ENDDO
        ELSE
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              psic_d (dffts_nl_d (j)) =       psi_d(j, ibnd)
              psic_d (dffts_nlm_d(j)) = conjg(psi_d(j, ibnd))
           ENDDO
        ENDIF
        !
     ENDIF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        CALL invfft ('tgWave', tg_psic_d, dffts )
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        !$cuf kernel do(1) <<<,>>>
        DO j = 1, dffts%nr1x * dffts%nr2x * right_nr3
           tg_psic_d (j) = tg_psic_d (j) * tg_v_d(j)
        ENDDO
        !
        CALL fwfft ('tgWave', tg_psic_d, dffts )
        !
     ELSE
        !
        CALL invfft ('Wave', psic_d, dffts)
        !
        !$cuf kernel do(1) <<<,>>>
        DO j = 1, dffts%nnr
           psic_d (j) = psic_d (j) * v_d(j)
        ENDDO
        !
        CALL fwfft ('Wave', psic_d, dffts)
        !
     ENDIF
     !
     !   addition to the total product
     !
     IF( use_tg ) THEN
        !
        ioff   = 0
        !
        CALL tg_get_recip_inc( dffts, right_inc )
        !
        DO idx = 1, 2*fftx_ntgrp(dffts), 2
           !
           IF( idx + ibnd - 1 < m ) THEN
              !$cuf kernel do(1) <<<,>>>
              DO j = 1, n
                 fp= ( tg_psic_d( dffts_nl_d(j) + ioff ) +  &
                       tg_psic_d( dffts_nlm_d(j) + ioff ) ) * 0.5d0
                 fm= ( tg_psic_d( dffts_nl_d(j) + ioff ) -  &
                       tg_psic_d( dffts_nlm_d(j) + ioff ) ) * 0.5d0
                 hpsi_d (j, ibnd+idx-1) = hpsi_d (j, ibnd+idx-1) + &
                                        cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi_d (j, ibnd+idx  ) = hpsi_d (j, ibnd+idx  ) + &
                                        cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              !$cuf kernel do(1) <<<,>>>
              DO j = 1, n
                 hpsi_d (j, ibnd+idx-1) = hpsi_d (j, ibnd+idx-1) + &
                                         tg_psic_d( dffts_nl_d(j) + ioff )
              ENDDO
           ENDIF
           !
           ioff = ioff + right_inc
           !
        ENDDO
        !
     ELSE
        IF (ibnd < m) THEN
           ! two ffts at the same time
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              fp = (psic_d (dffts_nl_d(j)) + psic_d (dffts_nlm_d(j)))*0.5d0
              fm = (psic_d (dffts_nl_d(j)) - psic_d (dffts_nlm_d(j)))*0.5d0
              hpsi_d (j, ibnd)   = hpsi_d (j, ibnd)   + &
                                 cmplx( dble(fp), aimag(fm),kind=DP)
              hpsi_d (j, ibnd+1) = hpsi_d (j, ibnd+1) + &
                                 cmplx(aimag(fp),- dble(fm),kind=DP)
           ENDDO
        ELSE
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              hpsi_d (j, ibnd)   = hpsi_d (j, ibnd)   + psic_d (dffts_nl_d(j))
           ENDDO
        ENDIF
     ENDIF
     !
  ENDDO
  !
  IF( use_tg ) THEN
     !
     CALL buffer%release_buffer( tg_psic_d, ierr )
     CALL buffer%release_buffer( tg_v_d, ierr )
     !
  ENDIF
  CALL buffer%release_buffer( dffts_nl_d, ierr )
  CALL buffer%release_buffer( dffts_nlm_d, ierr )
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_gamma_gpu
!
!@njs: vloc_psi_k
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k_gpu(lda, n, m, psi_d, v_d, hpsi_d)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - k-points
  !
  !   fft to real space
  !   product with the potential v on the smooth grid
  !   back to reciprocal space
  !   addition to the hpsi
  !
  USE parallel_include
  USE kinds, ONLY : DP
  USE wvfct, ONLY : current_k
  USE klist, ONLY : igk_k_host => igk_k
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  USE fft_helper_subroutines
  USE wavefunctions_module, ONLY: psic_h => psic
  USE wavefunctions_module_gpum, ONLY: psic_d
  USE buffers_module,       ONLY : buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(in)   :: psi_d (lda, m)
  COMPLEX(DP), DEVICE, INTENT(inout):: hpsi_d (lda, m)
  REAL(DP), DEVICE, INTENT(in) :: v_d(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  LOGICAL :: use_tg
  ! Task Groups
  REAL(DP),    DEVICE, POINTER :: tg_v_d(:)
  COMPLEX(DP), DEVICE, POINTER :: tg_psic_d(:)
  INTEGER,     DEVICE, POINTER :: dffts_nl_d(:)
  INTEGER,     DEVICE, POINTER :: igk_k_d(:)
  INTEGER :: v_siz, idx, ioff
  INTEGER :: ierr
  !
  IF (.not. allocated(psic_d)) ALLOCATE(psic_d, SOURCE=psic_h)
  CALL start_clock ('vloc_psi')
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dffts%nnr_tg
     !
     CALL buffer%lock_buffer( tg_v_d,    v_siz, ierr )
     CALL buffer%lock_buffer( tg_psic_d, v_siz, ierr )
     !
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
  ENDIF
  CALL buffer%lock_buffer( dffts_nl_d, size(dffts%nl) , ierr )
  CALL buffer%lock_buffer( igk_k_d, n, ierr )
  !
  dffts_nl_d = dffts%nl
  igk_k_d(1:n) = igk_k_host(1:n, current_k)
  !
  IF( use_tg ) THEN

     CALL tg_get_nnr( dffts, right_nnr )

     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !
        tg_psic_d = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, fftx_ntgrp(dffts)

           IF( idx + ibnd - 1 <= m ) THEN
!$cuf kernel do(1) <<<,>>>
              DO j = 1, n
                 tg_psic_d(dffts_nl_d (igk_k_d(j))+ioff) =  psi_d(j,idx+ibnd-1)
              ENDDO

           ENDIF

        !write (6,*) 'wfc G ', idx+ibnd-1
        !write (6,99) (tg_psic(i+ioff), i=1,400)

           ioff = ioff + right_nnr
        ENDDO
        !
        CALL  invfft ('tgWave', tg_psic_d, dffts )
        !write (6,*) 'wfc R ' 
        !write (6,99) (tg_psic(i), i=1,400)
        !
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
!$cuf kernel do(1) <<<,>>>
        DO j = 1, dffts%nr1x*dffts%nr2x* right_nr3
           tg_psic_d (j) = tg_psic_d (j) * tg_v_d(j)
        ENDDO
!
        !write (6,*) 'v psi R ' 
        !write (6,99) (tg_psic(i), i=1,400)
        !
        CALL fwfft ('tgWave',  tg_psic_d, dffts )
        !
        !   addition to the total product
        !
        ioff   = 0
        !
        CALL tg_get_recip_inc( dffts, right_inc )
        !
        DO idx = 1, fftx_ntgrp(dffts)
           !
           IF( idx + ibnd - 1 <= m ) THEN
!$cuf kernel do(1) <<<,>>>
              DO j = 1, n
                 hpsi_d (j, ibnd+idx-1) = hpsi_d (j, ibnd+idx-1) + &
                    tg_psic_d( dffts_nl_d(igk_k_d(j)) + ioff )
              ENDDO
!
           ENDIF
           !
        !write (6,*) 'v psi G ', idx+ibnd-1
        !write (6,99) (tg_psic(i+ioff), i=1,400)

           ioff = ioff + right_inc
           !
        ENDDO
        !
     ENDDO
  ELSE
     DO ibnd = 1, m
        !
        !!! == OPTIMIZE HERE == (setting to 0 and setting elements!)
        psic_d(:) = (0.d0, 0.d0)
!$cuf kernel do(1) <<<,>>>
        DO j = 1, n
          psic_d (dffts_nl_d (igk_k_d(j))) = psi_d(j, ibnd)
        END DO
!
        !write (6,*) 'wfc G ', ibnd
        !write (6,99) (psic(i), i=1,400)
        !
        CALL invfft ('Wave', psic_d, dffts)
        !write (6,*) 'wfc R ' 
        !write (6,99) (psic(i), i=1,400)
        !
!$cuf kernel do(1) <<<,>>>
        DO j = 1, dffts%nnr
           psic_d (j) = psic_d (j) * v_d(j)
        ENDDO
!
        !write (6,*) 'v psi R ' 
        !write (6,99) (psic(i), i=1,400)
        !
        CALL fwfft ('Wave', psic_d, dffts)
        !
        !   addition to the total product
        !
!$cuf kernel do(1) <<<,>>>
        DO j = 1, n
           hpsi_d (j, ibnd)   = hpsi_d (j, ibnd)   + psic_d (dffts_nl_d(igk_k_d(j)))
        ENDDO
!
        !write (6,*) 'v psi G ', ibnd
        !write (6,99) (psic(i), i=1,400)
        !
     ENDDO
  ENDIF
  !
  IF( use_tg ) THEN
     !
     CALL buffer%release_buffer( tg_psic_d, ierr )
     CALL buffer%release_buffer( tg_v_d, ierr )
     !
  ENDIF
  CALL buffer%release_buffer( dffts_nl_d, ierr )
  CALL buffer%release_buffer( igk_k_d, ierr )
  CALL stop_clock ('vloc_psi')
  !
99 format ( 20 ('(',2f12.9,')') )

  RETURN
END SUBROUTINE vloc_psi_k_gpu
!
!@njs: vloc_psi_nc
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc_gpu (lda, n, m, psi_d, v_d, hpsi_d)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - noncolinear
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE wvfct, ONLY : current_k
  USE klist, ONLY : igk_k_host => igk_k
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE lsda_mod,      ONLY : nspin
  USE spin_orb,      ONLY : domag
  USE noncollin_module,     ONLY: npol
  USE wavefunctions_module,       ONLY: psic_nc_h => psic_nc
  USE wavefunctions_module_gpum,  ONLY: psic_nc_d
  USE fft_helper_subroutines
  USE buffers_module,       ONLY : buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  REAL(DP), DEVICE, INTENT(in) :: v_d(dfftp%nnr,4) ! beware dimensions!
  COMPLEX(DP), DEVICE, INTENT(in)   :: psi_d (lda*npol, m)
  COMPLEX(DP), DEVICE, INTENT(inout):: hpsi_d (lda,npol,m)
  INTEGER,     DEVICE, POINTER :: dffts_nl_d(:)
  INTEGER,     DEVICE, POINTER :: igk_k_d(:)
  !
  INTEGER :: ibnd, j,ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    DEVICE, ALLOCATABLE :: tg_v_d(:,:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: tg_psic_d(:,:)
  INTEGER :: v_siz, idx, ioff
  INTEGER :: right_nnr, right_nr3, right_inc
  INTEGER :: ierr
  !
  CALL start_clock ('vloc_psi')
  !
  IF (.not. allocated(psic_nc_d)) ALLOCATE(psic_nc_d, SOURCE=psic_nc_h)
  incr = 1
  !
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz = dffts%nnr_tg
     IF (domag) THEN
        ALLOCATE( tg_v_d( v_siz, 4 ) )
        DO is=1,nspin
           CALL tg_gather_gpu( dffts, v_d(:,is), tg_v_d(:,is) )
        ENDDO
     ELSE
        ALLOCATE( tg_v_d( v_siz, 1 ) )
        CALL tg_gather_gpu( dffts, v_d(:,1), tg_v_d(:,1) )
     ENDIF
     ALLOCATE( tg_psic_d( v_siz, npol ) )
     CALL stop_clock ('vloc_psi:tg_gather')

     incr = fftx_ntgrp(dffts)
  ENDIF
  CALL buffer%lock_buffer( dffts_nl_d, size(dffts%nl), ierr )
  CALL buffer%lock_buffer( igk_k_d, n, ierr )
  dffts_nl_d = dffts%nl
  igk_k_d(1:n) = igk_k_host(1:n, current_k)
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr

     IF( use_tg ) THEN
        !
        CALL tg_get_nnr( dffts, right_nnr )
        !
        DO ipol = 1, npol
           !
           !  == OPTIMIZE HERE == setting data twice
           tg_psic_d(:,ipol) = ( 0.D0, 0.D0 )
           ioff   = 0
           !
           DO idx = 1, fftx_ntgrp(dffts)
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 !$cuf kernel do(1) <<<,>>>
                 DO j = 1, n
                    tg_psic_d( dffts_nl_d( igk_k_d(j) ) + ioff, ipol ) = &
                       psi_d( j +(ipol-1)*lda, idx+ibnd-1 )
                 ENDDO
              ENDIF

              ioff = ioff + right_nnr

           ENDDO
           !
           CALL invfft ('tgWave', tg_psic_d(:,ipol), dffts )
           !
        ENDDO
        !
     ELSE
        psic_nc_d = (0.d0,0.d0)
        DO ipol=1,npol
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              psic_nc_d(dffts_nl_d(igk_k_d(j)),ipol) = psi_d(j+(ipol-1)*lda,ibnd)
           ENDDO
           CALL invfft ('Wave', psic_nc_d(:,ipol), dffts)
        ENDDO
     ENDIF

     !
     !   product with the potential v = (vltot+vr) on the smooth grid
     !
     IF( use_tg ) THEN
        CALL tg_get_group_nr3( dffts, right_nr3 )
        IF (domag) THEN
           !$cuf kernel do(1) <<<,>>>
           DO j=1, dffts%nr1x*dffts%nr2x*right_nr3
              sup = tg_psic_d(j,1) * (tg_v_d(j,1)+tg_v_d(j,4)) + &
                    tg_psic_d(j,2) * (tg_v_d(j,2)-(0.d0,1.d0)*tg_v_d(j,3))
              sdwn = tg_psic_d(j,2) * (tg_v_d(j,1)-tg_v_d(j,4)) + &
                     tg_psic_d(j,1) * (tg_v_d(j,2)+(0.d0,1.d0)*tg_v_d(j,3))
              tg_psic_d(j,1)=sup
              tg_psic_d(j,2)=sdwn
           ENDDO
        ELSE
           !$cuf kernel do(1) <<<,>>>
           DO j=1, dffts%nr1x*dffts%nr2x*right_nr3
              tg_psic_d(j,:) = tg_psic_d(j,:) * tg_v_d(j,1)
           ENDDO
        ENDIF
     ELSE
        IF (domag) THEN
           !$cuf kernel do(1) <<<,>>>
           DO j=1, dffts%nnr
              sup = psic_nc_d(j,1) * (v_d(j,1)+v_d(j,4)) + &
                    psic_nc_d(j,2) * (v_d(j,2)-(0.d0,1.d0)*v_d(j,3))
              sdwn = psic_nc_d(j,2) * (v_d(j,1)-v_d(j,4)) + &
                     psic_nc_d(j,1) * (v_d(j,2)+(0.d0,1.d0)*v_d(j,3))
              psic_nc_d(j,1)=sup
              psic_nc_d(j,2)=sdwn
           ENDDO
        ELSE
           !$cuf kernel do(1) <<<,>>>
           DO j=1, dffts%nnr
              psic_nc_d(j,:) = psic_nc_d(j,:) * v_d(j,1)
           ENDDO
        ENDIF
     ENDIF
     !
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        DO ipol = 1, npol

           CALL fwfft ('tgWave', tg_psic_d(:,ipol), dffts )
           !
           ioff   = 0
           !
           CALL tg_get_recip_inc( dffts, right_inc )
           !
           DO idx = 1, fftx_ntgrp(dffts)
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 !$cuf kernel do(1) <<<,>>>
                 DO j = 1, n
                    hpsi_d (j, ipol, ibnd+idx-1) = hpsi_d (j, ipol, ibnd+idx-1) + &
                                 tg_psic_d( dffts_nl_d(igk_k_d(j)) + ioff, ipol )
                 ENDDO
              ENDIF
              !
              ioff = ioff + right_inc
              !
           ENDDO

        ENDDO
        !
     ELSE

        DO ipol=1,npol
           CALL fwfft ('Wave', psic_nc_d(:,ipol), dffts)
        ENDDO
        !
        !   addition to the total product
        !
        DO ipol=1,npol
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              hpsi_d(j,ipol,ibnd) = hpsi_d(j,ipol,ibnd) + &
                                  psic_nc_d(dffts_nl_d(igk_k_d(j)),ipol)
           ENDDO
        ENDDO

     ENDIF

  ENDDO

  IF( use_tg ) THEN
     !
     DEALLOCATE(tg_v_d, tg_psic_d)
     !
  ENDIF
  CALL buffer%release_buffer( dffts_nl_d, ierr )
  CALL buffer%release_buffer( igk_k_d, ierr )
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_nc_gpu
!@nje
#endif
