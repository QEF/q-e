!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma_gpu(lda, n, m, psi_d, v_d, hpsi_d)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - Gamma point
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE control_flags, ONLY : many_fft
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  USE fft_helper_subroutines
  USE qe_buffers,    ONLY : qe_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi_d (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi_d (lda, m)
  REAL(DP),    INTENT(in)   :: v_d(dffts%nnr)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, v_d
#endif
  !
  INTEGER :: ibnd, j, incr, right_nnr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  COMPLEX(DP), POINTER :: psic_d(:)
  REAL(DP),    POINTER :: tg_v_d(:)
  COMPLEX(DP), POINTER :: tg_psic_d(:)
  INTEGER,     POINTER :: dffts_nl_d(:), dffts_nlm_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: psic_d, tg_v_d, tg_psic_d, dffts_nl_d, dffts_nlm_d
#endif
  INTEGER :: v_siz, idx, ioff
  INTEGER :: ierr
  ! Variables to handle batched FFT
  INTEGER :: group_size, pack_size, remainder, howmany
  REAL(DP):: v_tmp
  !
  CALL start_clock ('vloc_psi')
  incr = 2 * many_fft
  !
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dffts%nnr_tg
     !
     CALL qe_buffer%lock_buffer( tg_v_d, v_siz, ierr )
     CALL qe_buffer%lock_buffer( tg_psic_d, v_siz, ierr )
     !
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
     incr = 2 * fftx_ntgrp(dffts)
     !
  ELSE
     CALL qe_buffer%lock_buffer( psic_d, dffts%nnr * many_fft, ierr )
  ENDIF
  ! Sync fft data
  dffts_nl_d => dffts%nl_d
  dffts_nlm_d => dffts%nlm_d
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
     ELSE IF (many_fft > 1) THEN
        !
        psic_d(:) = (0.d0, 0.d0)
        !
        ! FFT batching strategy is defined here:
        !  the buffer in psi_c can contain 2*many_fft bands.
        !  * group_size: is the number of bands that will be transformed.
        !  * pack_size:  is the number of slots in psi_c used to store
        !                the (couples) of bands
        !  * remainder:  can be 1 or 0, if 1, a spare band should be added
        !                in the the first slot of psi_c not occupied by
        !                the couples of bands.
        !
        group_size = MIN(2*many_fft, m - (ibnd -1))
        pack_size  = (group_size/2) ! This is FLOOR(group_size/2)
        remainder  = group_size - 2*pack_size
        howmany    = pack_size+remainder
        !
        ! two ffts at the same time
        IF ( pack_size > 0 ) THEN
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              DO idx = 0, pack_size-1
                 psic_d(dffts_nl_d (j) + idx*dffts%nnr)=      psi_d(j,ibnd+2*idx) + (0.0d0,1.d0)*psi_d(j,ibnd+2*idx+1)
                 psic_d(dffts_nlm_d(j) + idx*dffts%nnr)=conjg(psi_d(j,ibnd+2*idx) - (0.0d0,1.d0)*psi_d(j,ibnd+2*idx+1))
              end do
           ENDDO
        END IF
        IF (remainder > 0) THEN
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              psic_d (dffts_nl_d (j) + pack_size*dffts%nnr) =       psi_d(j, ibnd+group_size-1)
              psic_d (dffts_nlm_d(j) + pack_size*dffts%nnr) = conjg(psi_d(j, ibnd+group_size-1))
           ENDDO
        ENDIF
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
     ELSE IF ( many_fft > 1 ) THEN
        !
        CALL invfft ('Wave', psic_d, dffts, howmany=howmany)
        !
        !$cuf kernel do(1) <<<,>>>
        DO j = 1, dffts%nnr
           v_tmp = v_d(j)
           DO idx = 0, howmany-1
              psic_d (j + idx*dffts%nnr) = psic_d (j+ idx*dffts%nnr) * v_tmp
           END DO
        ENDDO
        !
        CALL fwfft ('Wave', psic_d, dffts, howmany=howmany)
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
     ELSE IF ( many_fft > 1 ) THEN
        IF ( pack_size > 0 ) THEN
           ! two ffts at the same time
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              DO idx = 0, pack_size-1
                 ioff = idx*dffts%nnr
                 fp = (psic_d (ioff + dffts_nl_d(j)) + psic_d (ioff + dffts_nlm_d(j)))*0.5d0
                 fm = (psic_d (ioff + dffts_nl_d(j)) - psic_d (ioff + dffts_nlm_d(j)))*0.5d0
                 hpsi_d (j, ibnd + idx*2)   = hpsi_d (j, ibnd + idx*2)   + &
                                    cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi_d (j, ibnd + idx*2 +1) = hpsi_d (j, ibnd + idx*2 +1) + &
                                    cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ENDDO
        END IF
        IF (remainder > 0) THEN
           !$cuf kernel do(1) <<<,>>>
           DO j = 1, n
              hpsi_d (j, ibnd + group_size-1)   = hpsi_d (j, ibnd + group_size-1)   + psic_d (pack_size*dffts%nnr + dffts_nl_d(j))
           ENDDO
        ENDIF
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
     CALL qe_buffer%release_buffer( tg_psic_d, ierr )
     CALL qe_buffer%release_buffer( tg_v_d, ierr )
     !
  ELSE
     CALL qe_buffer%release_buffer( psic_d, ierr )
  END IF
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
  USE klist, ONLY : igk_k_d
  USE mp_bands,      ONLY : me_bgrp
  USE control_flags, ONLY : many_fft
  USE fft_base,      ONLY : dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  USE fft_helper_subroutines
  USE wavefunctions, ONLY: psic_h => psic
  !USE wavefunctions_gpum, ONLY: psic_d
  USE qe_buffers,    ONLY : qe_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi_d (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi_d (lda, m)
  REAL(DP),    INTENT(in)   :: v_d(dffts%nnr)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, v_d
#endif
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  LOGICAL :: use_tg
  ! Task Groups
  COMPLEX(DP), POINTER :: psic_d(:)
  REAL(DP),    POINTER :: tg_v_d(:)
  COMPLEX(DP), POINTER :: tg_psic_d(:)
  INTEGER,     POINTER :: dffts_nl_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: psic_d, tg_v_d, tg_psic_d, dffts_nl_d
#endif
  !
  REAL(DP) :: v_tmp
  INTEGER :: v_siz, idx, ioff
  INTEGER :: ierr
  INTEGER :: group_size
  
  CALL start_clock ('vloc_psi')
  use_tg = dffts%has_task_groups
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dffts%nnr_tg
     !
     CALL qe_buffer%lock_buffer( tg_v_d,    v_siz, ierr )
     CALL qe_buffer%lock_buffer( tg_psic_d, v_siz, ierr )
     !
     CALL tg_gather_gpu( dffts, v_d, tg_v_d )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
  ELSE
     CALL qe_buffer%lock_buffer( psic_d, dffts%nnr*many_fft, ierr )
  ENDIF
  !
  dffts_nl_d => dffts%nl_d
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
                 tg_psic_d(dffts_nl_d (igk_k_d(j, current_k))+ioff) =  psi_d(j,idx+ibnd-1)
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
                    tg_psic_d( dffts_nl_d(igk_k_d(j, current_k)) + ioff )
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
  ELSE IF (many_fft > 1) THEN
     DO ibnd = 1, m, many_fft
        !
        !!! == OPTIMIZE HERE == (setting to 0 and setting elements!)
        psic_d(:) = (0.d0, 0.d0)
        group_size = MIN(many_fft, m - (ibnd -1))

!$cuf kernel do(1) <<<,>>>
        DO j = 1, n
           DO idx = 0, group_size-1
              psic_d (dffts_nl_d (igk_k_d(j, current_k)) + idx*dffts%nnr) = psi_d(j, ibnd+idx)
           END DO
        END DO
        !
        CALL invfft ('Wave', psic_d, dffts, howmany=group_size)
        !
!$cuf kernel do(1) <<<,>>>
        DO j = 1, dffts%nnr
           v_tmp = v_d(j)
           DO idx = 0, group_size-1
              psic_d (j + idx*dffts%nnr) = psic_d (j + idx*dffts%nnr) * v_tmp
           END DO
        ENDDO
        !
        CALL fwfft ('Wave', psic_d, dffts, howmany=group_size)
        !
        !   addition to the total product
        !
!$cuf kernel do(1) <<<*,*>>>
        DO j = 1, n
           DO idx = 0, group_size-1
              hpsi_d (j, ibnd + idx)   = hpsi_d (j, ibnd + idx)   + psic_d (dffts_nl_d(igk_k_d(j, current_k)) + idx * dffts%nnr)
           ENDDO
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
          psic_d (dffts_nl_d (igk_k_d(j, current_k))) = psi_d(j, ibnd)
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
           hpsi_d (j, ibnd)   = hpsi_d (j, ibnd)   + psic_d (dffts_nl_d(igk_k_d(j, current_k)))
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
     CALL qe_buffer%release_buffer( tg_psic_d, ierr )
     CALL qe_buffer%release_buffer( tg_v_d, ierr )
     !
  ELSE
     CALL qe_buffer%release_buffer( psic_d, ierr )
  ENDIF
  !
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
  USE klist, ONLY : igk_k_d
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE lsda_mod,      ONLY : nspin
  USE spin_orb,      ONLY : domag
  USE noncollin_module,     ONLY: npol
  USE wavefunctions_gpum,  ONLY: psic_nc_d
  USE fft_helper_subroutines
  USE qe_buffers,    ONLY : qe_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  REAL(DP),    INTENT(in)   :: v_d(dfftp%nnr,4) ! beware dimensions!
  COMPLEX(DP), INTENT(in)   :: psi_d (lda*npol, m)
  COMPLEX(DP), INTENT(inout):: hpsi_d (lda,npol,m)
#if defined(__CUDA)
  attributes(DEVICE) :: v_d, psi_d, hpsi_d
#endif
  !
  INTEGER :: ibnd, j,ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic_d(:,:)
  INTEGER,     POINTER      :: dffts_nl_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: tg_v_d, tg_psic_d, dffts_nl_d
#endif
  INTEGER :: v_siz, idx, ioff
  INTEGER :: right_nnr, right_nr3, right_inc
  INTEGER :: ierr
  !
  CALL start_clock ('vloc_psi')
  !
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
  !CALL qe_buffer%lock_buffer( dffts_nl_d, size(dffts%nl), ierr )
  dffts_nl_d => dffts%nl_d
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
                    tg_psic_d( dffts_nl_d( igk_k_d(j, current_k) ) + ioff, ipol ) = &
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
              psic_nc_d(dffts_nl_d(igk_k_d(j, current_k)),ipol) = psi_d(j+(ipol-1)*lda,ibnd)
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
                                 tg_psic_d( dffts_nl_d(igk_k_d(j, current_k)) + ioff, ipol )
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
                                  psic_nc_d(dffts_nl_d(igk_k_d(j, current_k)),ipol)
           ENDDO
        ENDDO

     ENDIF

  ENDDO

  IF( use_tg ) THEN
     !
     DEALLOCATE(tg_v_d, tg_psic_d)
     !
  ENDIF
  !CALL qe_buffer%release_buffer( dffts_nl_d, ierr )
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_nc_gpu

