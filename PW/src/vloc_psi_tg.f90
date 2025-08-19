!
! Copyright (C) 2003-2024 QUantum ESPRESS Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_tg_gamma( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique and "task groups"
  !
  USE parallel_include
  USE kinds,                   ONLY : DP
  USE mp_bands,                ONLY : me_bgrp
  USE fft_base,                ONLY : dffts
  USE fft_wave
  USE wavefunctions,           ONLY : psic
  USE fft_helper_subroutines,  ONLY : fftx_ntgrp, tg_get_group_nr3, &
                                      tg_get_recip_inc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  COMPLEX(DP), ALLOCATABLE :: vpsi(:,:) 
  ! ... Variables for task groups
  INTEGER :: v_siz, idx, ebnd, brange
  REAL(DP) :: fac
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:)
  !
  CALL start_clock( 'vloc_psi' )
  !
  IF (.not. dffts%has_task_groups ) CALL errore('vloc_psi','no task groups?',1)
  !
  CALL start_clock( 'vloc_psi:tg_gather' )
  incr = 2*fftx_ntgrp(dffts)     
  v_siz =  dffts%nnr_tg
  ALLOCATE( tg_v(v_siz) )
  ALLOCATE( tg_psic(v_siz) )
  CALL tg_gather( dffts, v, tg_v )
  ALLOCATE( tg_vpsi(n,incr) )
  CALL stop_clock( 'vloc_psi:tg_gather' )
  !
  DO ibnd = 1, m, incr
     !
     ! ... the local potential V_Loc psi. First bring psi to real space
     !
     CALL tgwave_g2r( psi(:,ibnd:m), tg_psic, dffts, n )
     !
     CALL tg_get_group_nr3( dffts, right_nr3 )
     !
     ! ... product with the potential v on the smooth grid
     !
     DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
        tg_psic(j) = tg_psic(j) * tg_v(j)
     ENDDO
     !
     ! ... back to reciprocal space
     ! ... addition to the total product
     !
     brange = m-ibnd+1
     !
     CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n )
     !
     DO idx = 1, 2*fftx_ntgrp(dffts), 2
        IF ( idx+ibnd-1<m ) THEN
           DO j = 1, n
              hpsi(j,ibnd+idx-1) = hpsi(j,ibnd+idx-1) + 0.5d0 * tg_vpsi(j,idx)
              hpsi(j,ibnd+idx)   = hpsi(j,ibnd+idx) + 0.5d0 * tg_vpsi(j,idx+1)
           ENDDO
        ELSEIF ( idx+ibnd-1==m ) THEN
           DO j = 1, n
              hpsi(j,ibnd+idx-1) = hpsi(j,ibnd+idx-1) + tg_vpsi(j,idx)
           ENDDO
        ENDIF
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE( tg_psic )
  DEALLOCATE( tg_v )
  DEALLOCATE( tg_vpsi )
  !
  CALL stop_clock( 'vloc_psi' )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_tg_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_tg_k( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points:
  !
  !! * fft to real space;
  !! * product with the potential v on the smooth grid;
  !! * back to reciprocal space;
  !! * addition to the hpsi.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, tg_get_group_nr3
  USE wavefunctions,          ONLY : psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
  !! Hamiltonian dot psi
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !! the total pot. in real space (smooth grid) for current spin
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr
  INTEGER :: i, iin, right_nnr, right_nr3, right_inc
  COMPLEX(DP), ALLOCATABLE :: vpsi(:,:)
  ! ... chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
  ! ... Task Groups
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:), tg_vpsi(:,:)
  INTEGER :: v_siz, idx, brange
  !
  IF (.not. dffts%has_task_groups ) CALL errore('vloc_psi','no task groups?',2)
  !
  CALL start_clock( 'vloc_psi' )
  !
  CALL start_clock( 'vloc_psi:tg_gather' )
  v_siz = dffts%nnr_tg
  incr = fftx_ntgrp(dffts)
  ALLOCATE( tg_v(v_siz) )
  ALLOCATE( tg_psic(v_siz), tg_vpsi(lda,incr) )
  CALL tg_gather( dffts, v, tg_v )
  CALL stop_clock( 'vloc_psi:tg_gather' )
  !
  CALL tg_get_nnr( dffts, right_nnr )
  !
  ! ... compute the number of chuncks
  numblock = (n+blocksize-1)/blocksize
  !
  DO ibnd = 1, m, fftx_ntgrp(dffts)
     !
     CALL tgwave_g2r( psi(:,ibnd:m), tg_psic, dffts, n, igk_k(:,current_k) )
     !
     !        write (6,*) 'wfc R '
     !        write (6,99) (tg_psic(i), i=1,400)
     !
     CALL tg_get_group_nr3( dffts, right_nr3 )
     !
     !$omp parallel do
     DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
        tg_psic(j) = tg_psic(j) * tg_v(j)
     ENDDO
     !$omp end parallel do
     !
     !        write (6,*) 'v psi R ' 
     !        write (6,99) (tg_psic(i), i=1,400)
     !
     brange = m-ibnd+1
     !
     CALL tgwave_r2g( tg_psic, tg_vpsi(:,1:brange), dffts, n, igk_k(:,current_k) )
     !
     !$omp parallel do collapse(2)
     DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
        DO j = 1, numblock
           DO iin = (j-1)*blocksize+1, MIN(j*blocksize,n)
              hpsi(iin,ibnd+idx) = hpsi(iin,ibnd+idx) + tg_vpsi(iin,idx+1)
           ENDDO
        ENDDO
     ENDDO
     !$omp end parallel do
     !
  ENDDO
  !
  DEALLOCATE( tg_psic, tg_vpsi )
  DEALLOCATE( tg_v )
  !
  CALL stop_clock( 'vloc_psi' )
  !
99 format ( 20 ('(',2f12.9,')') )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_tg_k
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_tg_nc( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - noncollinear case.
  !
  USE parallel_include
  USE kinds,                  ONLY : DP
  USE wvfct,                  ONLY : current_k
  USE klist,                  ONLY : igk_k
  USE mp_bands,               ONLY : me_bgrp
  USE fft_base,               ONLY : dffts, dfftp
  USE fft_wave
  USE lsda_mod,               ONLY : nspin
  USE noncollin_module,       ONLY : npol, domag
  USE wavefunctions,          ONLY : psic_nc
  USE fft_helper_subroutines, ONLY : fftx_ntgrp, tg_get_nnr, &
                                     tg_get_group_nr3, tg_get_recip_inc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  REAL(DP), INTENT(IN) :: v(dfftp%nnr,4) ! beware dimensions!
  !! the total pot. in real space (smooth grid)
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,npol,m)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j,ipol, incr, is, ii, ie
  COMPLEX(DP) :: sup, sdwn
  COMPLEX(DP), ALLOCATABLE :: vpsi(:,:)
  ! ... Variables for task groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:), tg_vpsi(:,:)
  INTEGER :: v_siz, idx, ioff, brange
  INTEGER :: right_nr3, right_inc
  !
  IF (.not. dffts%has_task_groups ) CALL errore('vloc_psi','no task groups?',3)
  !
  CALL start_clock( 'vloc_psi' )
  !
  CALL start_clock( 'vloc_psi:tg_gather' )
  incr = fftx_ntgrp(dffts)
  v_siz = dffts%nnr_tg
  IF (domag) THEN
     ALLOCATE( tg_v(v_siz,4) )
     DO is = 1, nspin
        CALL tg_gather( dffts, v(:,is), tg_v(:,is) )
     ENDDO
  ELSE
     ALLOCATE( tg_v(v_siz,1) )
     CALL tg_gather( dffts, v(:,1), tg_v(:,1) )
  ENDIF
  ALLOCATE( tg_psic(v_siz,npol), tg_vpsi(lda,incr) )
  CALL stop_clock( 'vloc_psi:tg_gather' )
  !
  DO ibnd = 1, m, incr
     !
     DO ipol = 1, npol
        ii = lda*(ipol-1)+1
        ie = lda*(ipol-1)+n
        CALL tgwave_g2r( psi(ii:ie,ibnd:m), tg_psic(:,ipol), dffts, n, &
             igk_k(:,current_k) )
     ENDDO
     !
     CALL tg_get_group_nr3( dffts, right_nr3 )
     !
     IF (domag) THEN
        DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
           sup  = tg_psic(j,1) * (tg_v(j,1)+tg_v(j,4)) + &
                tg_psic(j,2) * (tg_v(j,2)-(0.d0,1.d0)*tg_v(j,3))
           sdwn = tg_psic(j,2) * (tg_v(j,1)-tg_v(j,4)) + &
                tg_psic(j,1) * (tg_v(j,2)+(0.d0,1.d0)*tg_v(j,3))
           tg_psic(j,1) = sup
           tg_psic(j,2) = sdwn
        ENDDO
     ELSE
        DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
           tg_psic(j,:) = tg_psic(j,:) * tg_v(j,1)
        ENDDO
     ENDIF
     !
     brange = m-ibnd+1
     !
     DO ipol = 1, npol
        !
        CALL tgwave_r2g( tg_psic(:,ipol), tg_vpsi(:,1:brange), dffts, n, &
             igk_k(:,current_k) )
        !
        CALL tg_get_recip_inc( dffts, right_inc )
        !
        ioff = 0
        !$omp parallel do
        DO idx = 1, fftx_ntgrp(dffts)
           IF ( idx+ibnd-1<=m ) THEN
              DO j = 1, n
                 hpsi(j,ipol,ibnd+idx-1) = hpsi(j,ipol,ibnd+idx-1) + tg_vpsi(j,idx)
              ENDDO
           ENDIF
           ioff = ioff + right_inc
        ENDDO
        !$omp end parallel do
        !
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE( tg_v )
  DEALLOCATE( tg_psic, tg_vpsi )
  !
  CALL stop_clock ('vloc_psi')
  !
  RETURN
  !
END SUBROUTINE vloc_psi_tg_nc
