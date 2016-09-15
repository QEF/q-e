!
! Copyright (C) 2003-2013 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_gamma(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - Gamma point
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, dtgs
  USE fft_parallel,  ONLY : tg_gather
  USE fft_interfaces,ONLY : fwfft, invfft
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr
  COMPLEX(DP) :: fp, fm
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff
  !
  CALL start_clock ('vloc_psi')
  incr = 2
  !
  use_tg = dtgs%have_task_groups 
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dtgs%tg_nnr * dtgs%nogrp
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, dtgs, v, tg_v )
     CALL stop_clock ('vloc_psi:tg_gather')
     !
     incr = 2 * dtgs%nogrp
     !
  ENDIF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*dtgs%nogrp, 2
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 tg_psic(nls (j)+ioff) =        psi(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd)
                 tg_psic(nlsm(j)+ioff) = conjg( psi(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) )
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 tg_psic(nls (j)+ioff) =        psi(j,idx+ibnd-1)
                 tg_psic(nlsm(j)+ioff) = conjg( psi(j,idx+ibnd-1) )
              ENDDO
           ENDIF

           ioff = ioff + dtgs%tg_nnr

        ENDDO
        !
     ELSE
        !
        psic(:) = (0.d0, 0.d0)
        IF (ibnd < m) THEN
           ! two ffts at the same time
           DO j = 1, n
              psic(nls (j))=      psi(j,ibnd) + (0.0d0,1.d0)*psi(j,ibnd+1)
              psic(nlsm(j))=conjg(psi(j,ibnd) - (0.0d0,1.d0)*psi(j,ibnd+1))
           ENDDO
        ELSE
           DO j = 1, n
              psic (nls (j)) =       psi(j, ibnd)
              psic (nlsm(j)) = conjg(psi(j, ibnd))
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
        CALL invfft ('Wave', tg_psic, dffts, dtgs)
        !
        DO j = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO
        !
        CALL fwfft ('Wave', tg_psic, dffts, dtgs)
        !
     ELSE
        !
        CALL invfft ('Wave', psic, dffts)
        !
        DO j = 1, dffts%nnr
           psic (j) = psic (j) * v(j)
        ENDDO
        !
        CALL fwfft ('Wave', psic, dffts)
        !
     ENDIF
     !
     !   addition to the total product
     !
     IF( use_tg ) THEN
        !
        ioff   = 0
        !
        DO idx = 1, 2*dtgs%nogrp, 2
           !
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 fp= ( tg_psic( nls(j) + ioff ) +  &
                       tg_psic( nlsm(j) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(j) + ioff ) -  &
                       tg_psic( nlsm(j) + ioff ) ) * 0.5d0
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                        cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + &
                                        cmplx(aimag(fp),- dble(fm),kind=DP)
              ENDDO
           ELSEIF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                                         tg_psic( nls(j) + ioff )
              ENDDO
           ENDIF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
           !
        ENDDO
        !
     ELSE
        IF (ibnd < m) THEN
           ! two ffts at the same time
           DO j = 1, n
              fp = (psic (nls(j)) + psic (nlsm(j)))*0.5d0
              fm = (psic (nls(j)) - psic (nlsm(j)))*0.5d0
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + &
                                 cmplx( dble(fp), aimag(fm),kind=DP)
              hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + &
                                 cmplx(aimag(fp),- dble(fm),kind=DP)
           ENDDO
        ELSE
           DO j = 1, n
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(j))
           ENDDO
        ENDIF
     ENDIF
     !
  ENDDO
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  ENDIF
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - k-points
  !
  USE parallel_include
  USE kinds, ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE wvfct, ONLY : current_k
  USE klist, ONLY : igk_k
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, dtgs
  USE fft_parallel,  ONLY : tg_gather
  USE fft_interfaces,ONLY : fwfft, invfft
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(dffts%nnr)
  !
  INTEGER :: ibnd, j, incr
  !
  LOGICAL :: use_tg
  ! Task Groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff
  !
  CALL start_clock ('vloc_psi')
  use_tg = dtgs%have_task_groups 
  !
  incr = 1
  !
  IF( use_tg ) THEN
     !
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz =  dtgs%tg_nnr * dtgs%nogrp
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, dtgs, v, tg_v )
     CALL stop_clock ('vloc_psi:tg_gather')

     incr = dtgs%nogrp
     !
  ENDIF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  DO ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, dtgs%nogrp

           IF( idx + ibnd - 1 <= m ) THEN
!$omp parallel do
              DO j = 1, n
                 tg_psic(nls (igk_k(j,current_k))+ioff) =  psi(j,idx+ibnd-1)
              ENDDO
!$omp end parallel do
           ENDIF

           ioff = ioff + dtgs%tg_nnr

        ENDDO
        !
        CALL  invfft ('Wave', tg_psic, dffts, dtgs)
        !
     ELSE
        !
        psic(:) = (0.d0, 0.d0)
        psic (nls (igk_k(1:n,current_k))) = psi(1:n, ibnd)
        !
        CALL invfft ('Wave', psic, dffts)
        !
     ENDIF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
!$omp parallel do
        DO j = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        ENDDO
!$omp end parallel do
        !
        CALL fwfft ('Wave',  tg_psic, dffts, dtgs)
        !
     ELSE
        !
!$omp parallel do
        DO j = 1, dffts%nnr
           psic (j) = psic (j) * v(j)
        ENDDO
!$omp end parallel do
        !
        CALL fwfft ('Wave', psic, dffts)
        !
     ENDIF
     !
     !   addition to the total product
     !
     IF( use_tg ) THEN
        !
        ioff   = 0
        !
        DO idx = 1, dtgs%nogrp
           !
           IF( idx + ibnd - 1 <= m ) THEN
!$omp parallel do
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + &
                    tg_psic( nls(igk_k(j,current_k)) + ioff )
              ENDDO
!$omp end parallel do
           ENDIF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
           !
        ENDDO
        !
     ELSE
!$omp parallel do
        DO j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk_k(j,current_k)))
        ENDDO
!$omp end parallel do
     ENDIF
     !
  ENDDO
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  ENDIF
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_k
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_nc (lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - noncolinear
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gvecs, ONLY : nls, nlsm
  USE wvfct, ONLY : current_k
  USE klist, ONLY : igk_k
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dffts, dfftp, dtgs
  USE fft_parallel,  ONLY : tg_gather
  USE fft_interfaces,ONLY : fwfft, invfft
  USE lsda_mod,      ONLY : nspin
  USE spin_orb,      ONLY : domag
  USE noncollin_module,     ONLY: npol
  USE wavefunctions_module, ONLY: psic_nc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: lda, n, m
  REAL(DP), INTENT(in) :: v(dfftp%nnr,4) ! beware dimensions!
  COMPLEX(DP), INTENT(in)   :: psi (lda*npol, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda,npol,m)
  !
  INTEGER :: ibnd, j,ipol, incr, is
  COMPLEX(DP) :: sup, sdwn
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:)
  INTEGER :: v_siz, idx, ioff
  !
  CALL start_clock ('vloc_psi')
  !
  incr = 1
  !
  use_tg = dtgs%have_task_groups 
  !
  IF( use_tg ) THEN
     CALL start_clock ('vloc_psi:tg_gather')
     v_siz = dtgs%tg_nnr * dtgs%nogrp
     IF (domag) THEN
        ALLOCATE( tg_v( v_siz, 4 ) )
        DO is=1,nspin
           CALL tg_gather( dffts, dtgs, v(:,is), tg_v(:,is) )
        ENDDO
     ELSE
        ALLOCATE( tg_v( v_siz, 1 ) )
        CALL tg_gather( dffts, dtgs, v(:,1), tg_v(:,1) )
     ENDIF
     ALLOCATE( tg_psic( v_siz, npol ) )
     CALL stop_clock ('vloc_psi:tg_gather')

     incr = dtgs%nogrp
  ENDIF
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  DO ibnd = 1, m, incr

     IF( use_tg ) THEN
        !
        DO ipol = 1, npol
           !
           tg_psic(:,ipol) = ( 0.D0, 0.D0 )
           ioff   = 0
           !
           DO idx = 1, dtgs%nogrp
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    tg_psic( nls( igk_k(j,current_k) ) + ioff, ipol ) = &
                       psi( j +(ipol-1)*lda, idx+ibnd-1 )
                 ENDDO
              ENDIF

              ioff = ioff + dtgs%tg_nnr

           ENDDO
           !
           CALL invfft ('Wave', tg_psic(:,ipol), dffts, dtgs)
           !
        ENDDO
        !
     ELSE
        psic_nc = (0.d0,0.d0)
        DO ipol=1,npol
           DO j = 1, n
              psic_nc(nls(igk_k(j,current_k)),ipol) = psi(j+(ipol-1)*lda,ibnd)
           ENDDO
           CALL invfft ('Wave', psic_nc(:,ipol), dffts)
        ENDDO
     ENDIF

     !
     !   product with the potential v = (vltot+vr) on the smooth grid
     !
     IF( use_tg ) THEN
        IF (domag) THEN
           DO j=1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
              sup = tg_psic(j,1) * (tg_v(j,1)+tg_v(j,4)) + &
                    tg_psic(j,2) * (tg_v(j,2)-(0.d0,1.d0)*tg_v(j,3))
              sdwn = tg_psic(j,2) * (tg_v(j,1)-tg_v(j,4)) + &
                     tg_psic(j,1) * (tg_v(j,2)+(0.d0,1.d0)*tg_v(j,3))
              tg_psic(j,1)=sup
              tg_psic(j,2)=sdwn
           ENDDO
        ELSE
           DO j=1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
              tg_psic(j,:) = tg_psic(j,:) * tg_v(j,1)
           ENDDO
        ENDIF
     ELSE
        IF (domag) THEN
           DO j=1, dffts%nnr
              sup = psic_nc(j,1) * (v(j,1)+v(j,4)) + &
                    psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
              sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + &
                     psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
              psic_nc(j,1)=sup
              psic_nc(j,2)=sdwn
           ENDDO
        ELSE
           DO j=1, dffts%nnr
              psic_nc(j,:) = psic_nc(j,:) * v(j,1)
           ENDDO
        ENDIF
     ENDIF
     !
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        DO ipol = 1, npol

           CALL fwfft ('Wave', tg_psic(:,ipol), dffts, dtgs)
           !
           ioff   = 0
           !
           DO idx = 1, dtgs%nogrp
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    hpsi (j, ipol, ibnd+idx-1) = hpsi (j, ipol, ibnd+idx-1) + &
                                 tg_psic( nls(igk_k(j,current_k)) + ioff, ipol )
                 ENDDO
              ENDIF
              !
              ioff = ioff + dffts%nr3x * dffts%nsw( me_bgrp + 1 )
              !
           ENDDO

        ENDDO
        !
     ELSE

        DO ipol=1,npol
           CALL fwfft ('Wave', psic_nc(:,ipol), dffts)
        ENDDO
        !
        !   addition to the total product
        !
        DO ipol=1,npol
           DO j = 1, n
              hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + &
                                  psic_nc(nls(igk_k(j,current_k)),ipol)
           ENDDO
        ENDDO

     ENDIF

  ENDDO

  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_v )
     DEALLOCATE( tg_psic )
     !
  ENDIF
  CALL stop_clock ('vloc_psi')
  !
  RETURN
END SUBROUTINE vloc_psi_nc
