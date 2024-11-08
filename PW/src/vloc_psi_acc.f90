!
! Copyright (C) 2003-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma_acc( lda, n, m, psi, v, hpsi )
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
  COMPLEX(DP), INTENT(in)   :: psi(lda,m)
  COMPLEX(DP), INTENT(inout):: hpsi(lda,m)
  REAL(DP),    INTENT(in)   :: v(dffts%nnr)
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, incr, right_nnr, right_nr3, right_inc
  COMPLEX(DP) :: fp, fm
  !
  COMPLEX(DP), ALLOCATABLE :: psi1(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:)
  INTEGER :: dffts_nnr, idx, ebnd, brange
  INTEGER :: ierr, ioff
  ! ... Variables to handle batched FFT
  INTEGER :: group_size, pack_size, remainder, howmany, hm_vec(3)
  REAL(DP):: fac
  !
  IF ( dffts%has_task_groups ) CALL errore('Vloc_psi_acc','no task groups!',1)
  !
  CALL start_clock_gpu( 'vloc_psi' )
  !
  incr = 2*many_fft
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( psi1(n,incr) )
  ALLOCATE( psic(dffts_nnr*incr) )
  !$acc data present_or_copyout(hpsi) present_or_copyin(psi,v) create(psi1,psic)
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
              psi1(j,idx) = psi(j,ibnd+idx-1)
           END DO
        END DO
        !
        CALL wave_g2r( psi1(:,1:group_size), psic, dffts, howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, howmany-1
          DO j = 1, dffts_nnr
            psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v(j)
          ENDDO
        ENDDO
        !
        CALL wave_r2g( psic, psi1, dffts, howmany_set=hm_vec )
        !
        IF ( pack_size > 0 ) THEN
           !$acc parallel loop collapse(2)
           DO idx = 0, pack_size-1
              DO j = 1, n
                 hpsi(j,ibnd+idx*2)   = hpsi(j,ibnd+idx*2)   + psi1(j,idx*2+1)
                 hpsi(j,ibnd+idx*2+1) = hpsi(j,ibnd+idx*2+1) + psi1(j,idx*2+2)
              ENDDO
           ENDDO
        ENDIF
        !
        IF (remainder > 0) THEN
           !$acc parallel loop
           DO j = 1, n
              hpsi(j,ibnd+group_size-1) = hpsi(j,ibnd+group_size-1) + &
                                            psi1(j,group_size)
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
              psi1(j,idx) = psi(j,ibnd+idx-1)
           END DO
        END DO
        !
        CALL wave_g2r( psi1(:,1:brange), psic, dffts )
        !        
        !$acc parallel loop
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v(j)
        ENDDO
        !
        CALL wave_r2g( psic, psi1(:,1:brange), dffts )
        !
        fac=1.d0
        IF ( ibnd<m ) fac=0.5d0
        !
        !$acc parallel loop
        DO j = 1, n
          hpsi(j,ibnd) = hpsi(j,ibnd) + fac*psi1(j,1)
          IF ( ibnd<m ) hpsi(j,ibnd+1) = hpsi(j,ibnd+1) + fac*psi1(j,2)
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psi1 )
  DEALLOCATE( psic )
  !
  CALL stop_clock_gpu ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_gamma_acc
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k_acc( lda, n, m, psi, v, hpsi )
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
  COMPLEX(DP), INTENT(IN) :: psi(lda,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi(lda,m)
  REAL(DP), INTENT(IN) :: v(dffts%nnr)
  !
  ! ... local variables
  !
  INTEGER :: ibnd, ebnd, j, incr
  INTEGER :: i, right_nnr, right_nr3, right_inc
  !
  COMPLEX(DP), ALLOCATABLE :: psi1(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic(:)
  !
  INTEGER :: dffts_nnr, idx, group_size, hm_vec(3)
  INTEGER :: ierr, brange
  !
  IF ( dffts%has_task_groups ) CALL errore('Vloc_psi_acc','no task groups!',2)
  !
  CALL start_clock_gpu ('vloc_psi')
  !
  incr = many_fft
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( psi1(n,incr) )
  ALLOCATE( psic(dffts_nnr*incr) )
  !$acc data present_or_copyout(hpsi) present_or_copyin(psi,v) create(psi1,psic)
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
              psi1(j,idx) = psi(j,ibnd+idx-1)
           END DO
        END DO
        !
        CALL wave_g2r( psi1(:,1:group_size), psic, dffts, igk=igk_k(:,current_k),&
                       howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2) 
        DO idx = 0, group_size-1
           DO j = 1, dffts_nnr
              psic(idx*dffts_nnr+j) = psic(idx*dffts_nnr+j) * v(j)
           ENDDO
        ENDDO
        !
        CALL wave_r2g( psic, psi1, dffts, igk=igk_k(:,current_k), &
                       howmany_set=hm_vec )
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO j = 1, n
              hpsi(j,ibnd+idx) = hpsi(j,ibnd+idx) + psi1(j,idx+1)
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
           psi1(j,idx) = psi(j,ibnd+idx-1)
        END DO
        CALL wave_g2r( psi1(:,idx:idx), psic, dffts, igk=igk_k(:,current_k) )
        !
        !$acc parallel loop 
        DO j = 1, dffts_nnr
           psic(j) = psic(j) * v(j)
        ENDDO
        !
        CALL wave_r2g( psic, psi1(:,:), dffts, igk=igk_k(:,current_k) )
        !
        !$acc parallel loop
        DO i = 1, n
           hpsi(i,ibnd) = hpsi(i,ibnd) + psi1(i,idx)
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( psic )
  DEALLOCATE( psi1 )
  !
  CALL stop_clock_gpu( 'vloc_psi' )
  !
  RETURN
  !
END SUBROUTINE vloc_psi_k_acc
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc_acc( lda, n, m, psi, v, hpsi )
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
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(INOUT):: hpsi(lda,npol,m)
  !
  ! ... local variables
  !
  INTEGER :: ibnd, j, ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  COMPLEX(DP), ALLOCATABLE :: psi1(:,:), psic(:,:)
  INTEGER :: dffts_nnr, idx, ioff, ii, ie, brange
  INTEGER :: right_nnr, right_nr3, right_inc
  !
  IF ( dffts%has_task_groups ) CALL errore('Vloc_psi_acc','no task groups!',3)
  !
  CALL start_clock_gpu ('vloc_psi')
  !
  incr = 1
  dffts_nnr = dffts%nnr
  !
  ALLOCATE( psi1(n,npol) )
  ALLOCATE( psic(dffts_nnr,npol) )
  !
  !$acc data present_or_copyout(hpsi) present_or_copyin(psi,v) create(psi1,psic)
  !
  ! ... the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr
     !
     !$acc parallel loop collapse(2)
     DO ipol = 1, npol
        DO j = 1, n
           psi1(j,ipol) = psi(j+lda*(ipol-1),ibnd)
        END DO
     END DO
     DO ipol = 1, npol
        CALL wave_g2r( psi1(:,ipol:ipol), psic(:,ipol), dffts, &
             igk=igk_k(:,current_k) )
     ENDDO
     !
     IF (domag) THEN
        !$acc parallel loop
        DO j = 1, dffts_nnr
           sup  = psic(j,1) * (v(j,1)+v(j,4)) + &
                psic(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
           sdwn = psic(j,2) * (v(j,1)-v(j,4)) + &
                psic(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
           psic(j,1) = sup
           psic(j,2) = sdwn
        ENDDO
     ELSE
        !$acc parallel loop collapse(2)
        DO ipol = 1, npol
           DO j = 1, dffts_nnr
              psic(j,ipol) = psic(j,ipol) * v(j,1)
           ENDDO
        ENDDO
     ENDIF
     !
     DO ipol = 1, npol
        CALL wave_r2g( psic(:,ipol), psi1(:,1:1), dffts, &
             igk=igk_k(:,current_k) )
        !
        !$acc parallel loop
        DO j = 1, n
           hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + psi1(j,1)
        ENDDO
     ENDDO
     !
  ENDDO
  !$acc end data
  DEALLOCATE( psic )
  DEALLOCATE( psi1 )
  !
  CALL stop_clock_gpu ('vloc_psi')
  !
  RETURN
  !
END SUBROUTINE vloc_psi_nc_acc

