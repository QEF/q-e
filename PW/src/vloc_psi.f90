!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - Gamma point.
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
  COMPLEX(DP), ALLOCATABLE :: psi2(:,:)
  !
  !Variables for task groups
  LOGICAL :: use_tg
  INTEGER :: v_siz, idx, ebnd, brange
  REAL(DP) :: fac
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  !
  CALL start_clock( 'vloc_psi' )
  incr = 2
  !
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock( 'vloc_psi:tg_gather' )
     v_siz =  dffts%nnr_tg
     !
     ALLOCATE( tg_v(v_siz) )
     ALLOCATE( tg_psic(v_siz) )
     !
     CALL tg_gather( dffts, v, tg_v )
     CALL stop_clock( 'vloc_psi:tg_gather' )
     !
     incr = 2 * fftx_ntgrp(dffts)
     !
     ALLOCATE( psi2(n,incr) )
  ELSE
     ALLOCATE( psi2(n,2) )
  ENDIF
  !
  ! ... the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     ! ... FFT to real space
     !
     IF ( use_tg ) THEN
        CALL tgwave_g2r( psi, tg_psic, dffts, n, ibnd, m )
     ELSE
        ebnd = ibnd
        IF ( ibnd < m ) ebnd = ebnd + 1
        !
        CALL wave_g2r( psi(1:n,ibnd:ebnd), psic, dffts )
     ENDIF
     !
     ! ... product with the potential v on the smooth grid
     !
     IF ( use_tg ) THEN
        CALL tg_get_group_nr3( dffts, right_nr3 )
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*right_nr3
           tg_psic(j) = tg_psic(j) * tg_v(j)
        ENDDO
     ELSE
        DO j = 1, dffts%nnr
           psic(j) = psic(j) * v(j)
        ENDDO
     ENDIF
     !
     ! ... back to reciprocal space
     ! ... addition to the total product
     !
     IF( use_tg ) THEN
        !
        CALL tgwave_r2g( tg_psic, psi2, dffts, n, 1, m-ibnd+1 )
        !
        DO idx = 1, 2*fftx_ntgrp(dffts), 2
           IF ( idx+ibnd-1<m ) THEN
              DO j = 1, n
                 hpsi(j,ibnd+idx-1) = hpsi(j,ibnd+idx-1) + 0.5d0 * psi2(j,idx)
                 hpsi(j,ibnd+idx)   = hpsi(j,ibnd+idx) + 0.5d0 * psi2(j,idx+1)
              ENDDO
           ELSEIF ( idx+ibnd-1==m ) THEN
              DO j = 1, n
                 hpsi(j,ibnd+idx-1) = hpsi(j,ibnd+idx-1) + psi2(j,idx)
              ENDDO
           ENDIF
        ENDDO
        !
     ELSE
        !
        brange=1 ;  fac=1.d0
        IF ( ibnd<m ) THEN
          brange=2 ;  fac=0.5d0
        ENDIF
        !
        CALL wave_r2g( psic, psi2(:,1:brange), dffts )
        !
        DO j = 1, n
          hpsi(j,ibnd) = hpsi(j,ibnd) + fac*psi2(j,1)
          IF ( ibnd<m ) hpsi(j,ibnd+1) = hpsi(j,ibnd+1) + fac*psi2(j,2)
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF( use_tg ) THEN
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
  ENDIF
  DEALLOCATE( psi2 )
  !
  CALL stop_clock ('vloc_psi')
  !
  RETURN
  !
END SUBROUTINE vloc_psi_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k( lda, n, m, psi, v, hpsi )
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
  COMPLEX(DP), ALLOCATABLE :: psicv(:), psi2(:,:)
  !
  ! chunking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
  !
  ! Task Groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx
  !
  CALL start_clock( 'vloc_psi' )
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock( 'vloc_psi:tg_gather' )
     v_siz =  dffts%nnr_tg
     !
     incr = fftx_ntgrp(dffts)
     ALLOCATE( tg_v(v_siz) )
     ALLOCATE( tg_psic(v_siz) )
     ALLOCATE( psi2(lda,incr) )
     !
     CALL tg_gather( dffts, v, tg_v )
     CALL stop_clock( 'vloc_psi:tg_gather' )
  ELSE
     ALLOCATE( psicv(dffts%nnr) )
     ALLOCATE( psi2(lda,1) )
  ENDIF
  !
  IF( use_tg ) THEN
     !
     CALL tg_get_nnr( dffts, right_nnr )
     !
     ! ... compute the number of chuncks
     numblock = (n+blocksize-1)/blocksize
     !
     DO ibnd = 1, m, fftx_ntgrp(dffts)
        !
        CALL tgwave_g2r( psi, tg_psic, dffts, n, ibnd, m, igk_k(:,current_k) )
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
        CALL tgwave_r2g( tg_psic, psi2, dffts, n, 1, m-ibnd+1, igk_k(:,current_k) )
        !
        !$omp parallel do collapse(2)
        DO idx = 0, MIN(fftx_ntgrp(dffts)-1, m-ibnd)
           DO j = 1, numblock
              DO iin = (j-1)*blocksize+1, MIN(j*blocksize,n)
                 hpsi(iin,ibnd+idx) = hpsi(iin,ibnd+idx) + psi2(iin,idx+1)
              ENDDO
           ENDDO
        ENDDO
        !$omp end parallel do
        !
     ENDDO
  ELSE
     DO ibnd = 1, m
        !
        CALL wave_g2r( psi(1:n,ibnd:ibnd), psic, dffts, igk=igk_k(:,current_k) )
        !
!        write (6,*) 'wfc R '
!        write (6,99) (psic(i), i=1,400)
        !
        !$omp parallel do
        DO j = 1, dffts%nnr
           psicv(j) = psic(j) * v(j)
        ENDDO
        !$omp end parallel do
        !
!        write (6,*) 'v psi R '
!        write (6,99) (psic(i), i=1,400)
        !
        CALL wave_r2g( psicv, psi2(1:n,:), dffts, igk=igk_k(:,current_k) )
        !
        !$omp parallel do
        DO i = 1, n
           hpsi(i,ibnd) = hpsi(i,ibnd) + psi2(i,1)
        ENDDO
        !$omp end parallel do
        !
!        write (6,*) 'v psi G ', ibnd
!        write (6,99) (psic(i), i=1,400)
        !
     ENDDO
  ENDIF
  !
  IF ( use_tg ) THEN
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
  ELSE
     DEALLOCATE( psicv )
  ENDIF
  DEALLOCATE( psi2 )
  !
  CALL stop_clock ('vloc_psi')
  !
99 format ( 20 ('(',2f12.9,')') )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_k
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc( lda, n, m, psi, v, hpsi )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - noncollinear.
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
  COMPLEX(DP), ALLOCATABLE :: psi2(:,:)
  !
  ! Variables for task groups
  LOGICAL :: use_tg
  REAL(DP), ALLOCATABLE :: tg_v(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:)
  INTEGER :: v_siz, idx, ioff
  INTEGER :: right_nr3, right_inc
  !
  CALL start_clock ('vloc_psi')
  !
  incr = 1
  !
  use_tg = dffts%has_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock( 'vloc_psi:tg_gather' )
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
     ALLOCATE( tg_psic(v_siz,npol) )
     CALL stop_clock( 'vloc_psi:tg_gather' )
     !
     incr = fftx_ntgrp(dffts)
     ALLOCATE( psi2(lda,incr) )
  ELSE
     ALLOCATE( psi2(lda,1) )
  ENDIF
  !
  ! ... the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        DO ipol = 1, npol
           ii = lda*(ipol-1)+1
           ie = lda*(ipol-1)+n
           CALL tgwave_g2r( psi(ii:ie,:), tg_psic(:,ipol), dffts, n, ibnd, m, &
                            igk_k(:,current_k) )
        ENDDO
     ELSE
        psic_nc = (0.d0,0.d0)
        !
        DO ipol = 1, npol
           ii = lda*(ipol-1)+1
           ie = lda*(ipol-1)+n
           CALL wave_g2r( psi(ii:ie,ibnd:ibnd), psic_nc(:,ipol), dffts, &
                          igk=igk_k(:,current_k) )
        ENDDO
     ENDIF
     !
     ! ... product with the potential v = (vltot+vr) on the smooth grid
     !
     IF( use_tg ) THEN
        CALL tg_get_group_nr3( dffts, right_nr3 )
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
     ELSE
        IF (domag) THEN
           DO j = 1, dffts%nnr
              sup  = psic_nc(j,1) * (v(j,1)+v(j,4)) + &
                     psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
              sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + &
                     psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
              psic_nc(j,1) = sup
              psic_nc(j,2) = sdwn
           ENDDO
        ELSE
           DO j = 1, dffts%nnr
              psic_nc(j,:) = psic_nc(j,:) * v(j,1)
           ENDDO
        ENDIF
     ENDIF
     !
     ! ... back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        DO ipol = 1, npol
           !
           CALL tgwave_r2g( tg_psic(:,ipol), psi2, dffts, n, 1, m-ibnd+1, &
                            igk_k(:,current_k) )
           !
           CALL tg_get_recip_inc( dffts, right_inc )
           !
           ioff = 0
!$omp parallel do
           DO idx = 1, fftx_ntgrp(dffts)
             IF ( idx+ibnd-1<=m ) THEN
               DO j = 1, n
                 hpsi(j,ipol,ibnd+idx-1) = hpsi(j,ipol,ibnd+idx-1) + psi2(j,idx)
               ENDDO
             ENDIF
             ioff = ioff + right_inc
           ENDDO
!$omp end parallel do
           !
        ENDDO
        !
     ELSE
        !
        DO ipol = 1, npol
           !
           CALL wave_r2g( psic_nc(:,ipol), psi2(1:n,:), dffts, igk=igk_k(:,current_k) )
           !
!$omp parallel do
           DO j = 1, n
              hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + psi2(j,1)
           ENDDO
!$omp end parallel do
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF( use_tg ) THEN
     DEALLOCATE( tg_v )
     DEALLOCATE( tg_psic )
  ENDIF
  DEALLOCATE( psi2 )
  !
  CALL stop_clock ('vloc_psi')
  !
  RETURN
  !
END SUBROUTINE vloc_psi_nc
