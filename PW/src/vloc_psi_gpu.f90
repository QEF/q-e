!
! Copyright (C) 2003-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma_gpu( lda, n, m, psi_d, v, hpsi_d )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - Gamma point.
  !! - fft \psi to real space;
  !! - product with the potential v on the smooth grid;
  !! - back to reciprocal space;
  !! - add to hpsi.
  !
  USE parallel_include
  USE kinds,          ONLY : DP
  USE control_flags,  ONLY : many_fft
  USE mp_bands,       ONLY : me_bgrp
  USE fft_base,       ONLY : dffts
  USE fft_wave
  USE fft_helper_subroutines
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi_d(lda,m)
  COMPLEX(DP), INTENT(inout):: hpsi_d(lda,m)
  REAL(DP),    INTENT(in)   :: v(dffts%nnr)
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr, right_nnr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:)
  INTEGER :: dffts_nnr, idx, ebnd, brange
  INTEGER :: ierr, ioff
  ! ... Variables to handle batched FFT
  INTEGER :: group_size, pack_size, remainder, howmany, hm_vec(3)
  REAL(DP):: fac
  !
  IF ( dffts%has_task_groups ) CALL errore('Vloc_psi_gpu','no task groups!',1)
  !
  CALL start_clock_gpu( 'vloc_psi' )
  !
  incr = 2*many_fft
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( psi(n,incr) )
  ALLOCATE( psic(dffts_nnr*incr) )
  !$acc data deviceptr(psi_d,hpsi_d) present(v) create(psi,psic)
  !
  IF (many_fft > 1) THEN
     !
     DO ibnd = 1, m, incr
        !
        group_size = MIN(2*many_fft, m-(ibnd-1))
        pack_size = (group_size/2) ! This is FLOOR(group_size/2)
        remainder = group_size - 2*pack_size
        howmany = pack_size + remainder
        hm_vec(1)=group_size ; hm_vec(2)=n ; hm_vec(3)=howmany
        !$acc parallel loop collapse(2)
        DO idx = 1, group_size
           DO j = 1, n
              psi(j,idx) = psi_d(j,ibnd+idx-1)
           END DO
        END DO
        !
        CALL wave_g2r( psi(:,1:group_size), psic, dffts, howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, howmany-1
          DO j = 1, dffts_nnr
            psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v(j)
          ENDDO
        ENDDO
        !
        CALL wave_r2g( psic, psi, dffts, howmany_set=hm_vec )
        !
        IF ( pack_size > 0 ) THEN
           !$acc parallel loop collapse(2)
           DO idx = 0, pack_size-1
              DO j = 1, n
                 hpsi_d(j,ibnd+idx*2)   = hpsi_d(j,ibnd+idx*2)   + psi(j,idx*2+1)
                 hpsi_d(j,ibnd+idx*2+1) = hpsi_d(j,ibnd+idx*2+1) + psi(j,idx*2+2)
              ENDDO
           ENDDO
        ENDIF
        !
        IF (remainder > 0) THEN
           !$acc parallel loop
           DO j = 1, n
              hpsi_d(j,ibnd+group_size-1) = hpsi_d(j,ibnd+group_size-1) + &
                                            psi(j,group_size)
           ENDDO
        ENDIF
        !
     ENDDO
     !
  ELSE
     !
     DO ibnd = 1, m, incr
        !
        brange = 1
        IF ( ibnd<m ) brange = 2
        !$acc parallel loop collapse(2)
        DO idx = 1, brange
           DO j = 1, n
              psi(j,idx) = psi_d(j,ibnd+idx-1)
           END DO
        END DO
        !
        CALL wave_g2r( psi(:,1:brange), psic, dffts )
        !        
        !$acc parallel loop
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v(j)
        ENDDO
        !
        CALL wave_r2g( psic, psi(:,1:brange), dffts )
        !
        fac=1.d0
        IF ( ibnd<m ) fac=0.5d0
        !
        !$acc parallel loop
        DO j = 1, n
          hpsi_d(j,ibnd) = hpsi_d(j,ibnd) + fac*psi(j,1)
          IF ( ibnd<m ) hpsi_d(j,ibnd+1) = hpsi_d(j,ibnd+1) + fac*psi(j,2)
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psi )
  DEALLOCATE( psic )
  !
  CALL stop_clock_gpu ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_gamma_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k_gpu( lda, n, m, psi_d, v, hpsi_d )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - k-points. GPU double.
  !! - fft \psi to real space;
  !! - product with the potential v on the smooth grid;
  !! - back to reciprocal space;
  !! - add to hpsi.
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
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi_d(lda,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi_d(lda,m)
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ebnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:)
  !
  INTEGER :: dffts_nnr, idx, group_size, hm_vec(3)
  INTEGER :: ierr, brange
  !
  IF ( dffts%has_task_groups ) CALL errore('Vloc_psi_gpu','no task groups!',2)
  !
  CALL start_clock_gpu ('vloc_psi')
  !
  incr = many_fft
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( psi(n,incr) )
  ALLOCATE( psic(dffts_nnr*incr) )
  !$acc data deviceptr(psi_d,hpsi_d) present(v, igk_k) create(psi,psic)
  !
  IF (many_fft > 1) THEN
     !
     DO ibnd = 1, m, incr
        !
        group_size = MIN(many_fft,m-(ibnd-1))
        hm_vec(1)=group_size ; hm_vec(2)=n ; hm_vec(3)=group_size
        ebnd = ibnd+group_size-1
        !$acc parallel loop collapse(2)
        DO idx = 1, group_size
           DO j = 1, n
              psi(j,idx) = psi_d(j,ibnd+idx-1)
           END DO
        END DO
        !
        CALL wave_g2r( psi(:,1:group_size), psic, dffts, igk=igk_k(:,current_k),&
                       howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2) 
        DO idx = 0, group_size-1
           DO j = 1, dffts_nnr
              psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v(j)
           ENDDO
        ENDDO
        !
        CALL wave_r2g( psic, psi, dffts, igk=igk_k(:,current_k), &
                       howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO j = 1, n
              hpsi_d(j,ibnd+idx) = hpsi_d(j,ibnd+idx) + psi(j,idx+1)
           ENDDO
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
     DO ibnd = 1, m
        !
        idx = 1
        !$acc parallel loop
        DO j = 1, n
           psi(j,idx) = psi_d(j,ibnd+idx-1)
        END DO
        CALL wave_g2r( psi(:,idx:idx), psic, dffts, igk=igk_k(:,current_k) )
        !
        !$acc parallel loop 
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v(j)
        ENDDO
        !
        CALL wave_r2g( psic, psi(:,:), dffts, igk=igk_k(:,current_k) )
        !
        !$acc parallel loop
        DO i = 1, n
           hpsi_d(i,ibnd) = hpsi_d(i,ibnd) + psi(i,idx)
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psic )
  DEALLOCATE( psi )
  !
  CALL stop_clock_gpu( 'vloc_psi' )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_k_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc_gpu( lda, n, m, psi_d, v, hpsi_d )
  !-----------------------------------------------------------------------
  !! Calculation of Vloc*psi using dual-space technique - non-collinear - 
  !! - fft \psi to real space;
  !! - product with the potential v on the smooth grid;
  !! - back to reciprocal space;
  !! - add to hpsi.
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
  REAL(DP), INTENT(IN) :: v(dfftp%nnr,4) ! beware dimensions!
  COMPLEX(DP), INTENT(IN) :: psi_d(lda*npol,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi_d(lda,npol,m)
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), psic_nc(:,:)
  INTEGER :: dffts_nnr, idx, ioff, ii, ie, brange
  INTEGER :: right_nnr, right_nr3, right_inc
  !
  IF ( dffts%has_task_groups ) CALL errore('Vloc_psi_gpu','no task groups!',3)
  !
  CALL start_clock_gpu ('vloc_psi')
  !
  incr = 1
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( psi(n,npol) )
  ALLOCATE( psic_nc(dffts_nnr,npol) )
  !
  !$acc data deviceptr(psi_d,hpsi_d) present(v, igk_k) create(psi,psic_nc)
  !
  ! ... the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr
     !
     !$acc parallel loop collapse(2)
     DO ipol = 1, npol
        DO j = 1, n
           psi(j,ipol) = psi_d(j+lda*(ipol-1),ibnd)
        END DO
     END DO
     DO ipol = 1, npol
        CALL wave_g2r( psi(:,ipol:ipol), psic_nc(:,ipol), dffts, &
             igk=igk_k(:,current_k) )
     ENDDO
     !
     IF (domag) THEN
        !$acc parallel loop
        DO j = 1, dffts_nnr
           sup  = psic_nc(j,1) * (v(j,1)+v(j,4)) + &
                psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
           sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + &
                psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
           psic_nc(j,1) = sup
           psic_nc(j,2) = sdwn
        ENDDO
     ELSE
        !$acc parallel loop collapse(2)
        DO ipol = 1, npol
           DO j = 1, dffts_nnr
              psic_nc(j,ipol) = psic_nc(j,ipol) * v(j,1)
           ENDDO
        ENDDO
     ENDIF
     !
     DO ipol = 1, npol
        CALL wave_r2g( psic_nc(:,ipol), psi(:,1:1), dffts, &
             igk=igk_k(:,current_k) )
        !
        !$acc parallel loop
        DO j = 1, n
           hpsi_d(j,ipol,ibnd) = hpsi_d(j,ipol,ibnd) + psi(j,1)
        ENDDO
     ENDDO
     !
  ENDDO
  !$acc end data
  DEALLOCATE( psic_nc )
  DEALLOCATE( psi )
  !
  CALL stop_clock_gpu ('vloc_psi')
  !
  RETURN
  !
END SUBROUTINE vloc_psi_nc_gpu

