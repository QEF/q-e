!
! Copyright (C) 2003-2009 PWSCF group
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
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE wvfct,   ONLY : igk
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist, use_task_groups
  USE fft_parallel,  ONLY : tg_cft3s
  USE fft_base,      ONLY : dffts
  USE task_groups,   ONLY : tg_gather
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(nrxxs)
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
  !
  incr = 2
  !
  use_tg = ( use_task_groups ) .and. ( m >= nogrp )
  !
  if( use_tg ) then
     !
     v_siz =  dffts%nnrx * nogrp
     !
     allocate( tg_v   ( v_siz ) )
     allocate( tg_psic( v_siz ) )
     !
     call tg_gather( dffts, v, tg_v )
     !
     incr = 2 * nogrp
     !
  endif
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  do ibnd = 1, m, incr
     !
     if( use_tg ) then
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        do idx = 1, 2*nogrp, 2
           if( idx + ibnd - 1 < m ) then
              do j = 1, n
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd)
                 tg_psic(nlsm(igk(j))+ioff) = conjg( psi(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) )
              enddo
           elseif( idx + ibnd - 1 == m ) then
              do j = 1, n
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1)
                 tg_psic(nlsm(igk(j))+ioff) = conjg( psi(j,idx+ibnd-1) )
              enddo
           endif

           ioff = ioff + dffts%nnrx

        enddo
        !
     else
        !
        psic(:) = (0.d0, 0.d0)
        if (ibnd < m) then
           ! two ffts at the same time
           do j = 1, n
              psic (nls (igk(j))) =       psi(j, ibnd) + (0.0d0,1.d0)*psi(j, ibnd+1)
              psic (nlsm(igk(j))) = conjg(psi(j, ibnd) - (0.0d0,1.d0)*psi(j, ibnd+1))
           enddo
        else
           do j = 1, n
              psic (nls (igk(j))) =       psi(j, ibnd)
              psic (nlsm(igk(j))) = conjg(psi(j, ibnd))
           enddo
        endif
        !
     endif
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     if( use_tg ) then
        !
        call tg_cft3s ( tg_psic, dffts, 2, use_tg )
        !
        do j = 1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        enddo
        !
        call tg_cft3s ( tg_psic, dffts, -2, use_tg )
        !
     else
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        !
        do j = 1, nrxxs
           psic (j) = psic (j) * v(j)
        enddo
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
        !
     endif
     !
     !   addition to the total product
     !
     if( use_tg ) then
        !
        ioff   = 0
        !
        do idx = 1, 2*nogrp, 2
           !
           if( idx + ibnd - 1 < m ) then
              do j = 1, n
                 fp= ( tg_psic( nls(igk(j)) + ioff ) +  tg_psic( nlsm(igk(j)) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(igk(j)) + ioff ) -  tg_psic( nlsm(igk(j)) + ioff ) ) * 0.5d0
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + cmplx( dble(fp), aimag(fm),kind=DP)
                 hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + cmplx(aimag(fp),- dble(fm),kind=DP)
              enddo
           elseif( idx + ibnd - 1 == m ) then
              do j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + tg_psic( nls(igk(j)) + ioff )
              enddo
           endif
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
           !
        enddo
        !
     else
        if (ibnd < m) then
           ! two ffts at the same time
           do j = 1, n
              fp = (psic (nls(igk(j))) + psic (nlsm(igk(j))))*0.5d0
              fm = (psic (nls(igk(j))) - psic (nlsm(igk(j))))*0.5d0
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + cmplx( dble(fp), aimag(fm),kind=DP)
              hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + cmplx(aimag(fp),- dble(fm),kind=DP)
           enddo
        else
           do j = 1, n
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
           enddo
        endif
     endif
     !
  enddo
  !
  if( use_tg ) then
     !
     deallocate( tg_psic )
     deallocate( tg_v )
     !
  endif
  !
  return
END SUBROUTINE vloc_psi_gamma
!
!-----------------------------------------------------------------------
SUBROUTINE vloc_psi_k(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - k-points
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE wvfct,   ONLY : igk
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist, use_task_groups
  USE fft_parallel,  ONLY : tg_cft3s
  USE fft_base,      ONLY : dffts
  USE task_groups,   ONLY : tg_gather
  USE wavefunctions_module,  ONLY: psic
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(in)   :: psi (lda, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda, m)
  REAL(DP), INTENT(in) :: v(nrxxs)
  !
  INTEGER :: ibnd, j, incr
  !
  LOGICAL :: use_tg
  ! Task Groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff
  !
  !
  use_tg = ( use_task_groups ) .and. ( m >= nogrp )
  !
  incr = 1
  !
  if( use_tg ) then
     !
     v_siz =  dffts%nnrx * nogrp
     !
     allocate( tg_v   ( v_siz ) )
     allocate( tg_psic( v_siz ) )
     !
     call tg_gather( dffts, v, tg_v )
     incr = nogrp
     !
  endif
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  do ibnd = 1, m, incr
     !
     if( use_tg ) then
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        do idx = 1, nogrp

           if( idx + ibnd - 1 <= m ) then
!$omp parallel do
              do j = 1, n
                 tg_psic(nls (igk(j))+ioff) =  psi(j,idx+ibnd-1)
              enddo
!$omp end parallel do
           endif

           ioff = ioff + dffts%nnrx

        enddo
        !
        call tg_cft3s ( tg_psic, dffts, 2, use_tg )
        !
     else
        !
        psic(:) = (0.d0, 0.d0)
        psic (nls (igk(1:n))) = psi(1:n, ibnd)
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        !
     endif
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     if( use_tg ) then
        !
!$omp parallel do
        do j = 1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        enddo
!$omp end parallel do
        !
        call tg_cft3s ( tg_psic, dffts, -2, use_tg )
        !
     else
        !
!$omp parallel do
        do j = 1, nrxxs
           psic (j) = psic (j) * v(j)
        enddo
!$omp end parallel do
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
        !
     endif
     !
     !   addition to the total product
     !
     if( use_tg ) then
        !
        ioff   = 0
        !
        do idx = 1, nogrp
           !
           if( idx + ibnd - 1 <= m ) then
!$omp parallel do
              do j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + tg_psic( nls(igk(j)) + ioff )
              enddo
!$omp end parallel do
           endif
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
           !
        enddo
        !
     else
!$omp parallel do
        do j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
        enddo
!$omp end parallel do
     endif
     !
  enddo
  !
  if( use_tg ) then
     !
     deallocate( tg_psic )
     deallocate( tg_v )
     !
  endif
  !
  return
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
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE gvect,   ONLY : nrxx
  USE wvfct,   ONLY : igk
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist, use_task_groups
  USE fft_parallel,  ONLY : tg_cft3s
  USE fft_base,      ONLY : dffts
  USE task_groups,   ONLY : tg_gather
  USE noncollin_module,     ONLY: npol
  USE wavefunctions_module, ONLY: psic_nc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  REAL(DP), INTENT(in) :: v(nrxx,4)
  COMPLEX(DP), INTENT(in)   :: psi (lda*npol, m)
  COMPLEX(DP), INTENT(inout):: hpsi (lda,npol,m)
  !
  INTEGER :: ibnd, j,ipol, incr
  COMPLEX(DP) :: sup, sdwn
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:,:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:)
  INTEGER :: v_siz, idx, ioff
  !
  !
  incr = 1
  !
  use_tg = ( use_task_groups ) .and. ( m >= nogrp )
  !
  if( use_tg ) then
     v_siz = dffts%nnrx * nogrp
     allocate( tg_v( v_siz, 4 ) )
     call tg_gather( dffts, v(:,1), tg_v(:,1) )
     call tg_gather( dffts, v(:,2), tg_v(:,2) )
     call tg_gather( dffts, v(:,3), tg_v(:,3) )
     call tg_gather( dffts, v(:,4), tg_v(:,4) )
     allocate( tg_psic( v_siz, npol ) )
     incr = nogrp
  endif
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  do ibnd = 1, m, incr

     if( use_tg ) then
        !
        do ipol = 1, npol
           !
           tg_psic(:,ipol) = ( 0.D0, 0.D0 )
           ioff   = 0
           !
           do idx = 1, nogrp
              !
              if( idx + ibnd - 1 <= m ) then
                 do j = 1, n
                    tg_psic( nls( igk(j) ) + ioff, ipol ) = psi( j +(ipol-1)*lda, idx+ibnd-1 )
                 enddo
              endif

              ioff = ioff + dffts%nnrx

           enddo
           !
           call tg_cft3s ( tg_psic(:,ipol), dffts, 2, use_tg )
           !
        enddo
        !
     else
        psic_nc = (0.d0,0.d0)
        do ipol=1,npol
           do j = 1, n
              psic_nc(nls(igk(j)),ipol) = psi(j+(ipol-1)*lda,ibnd)
           enddo
           call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        enddo
     endif

     !
     !   product with the potential v = (vltot+vr) on the smooth grid
     !
     if( use_tg ) then
        do j=1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           sup = tg_psic(j,1) * (tg_v(j,1)+tg_v(j,4)) + &
                 tg_psic(j,2) * (tg_v(j,2)-(0.d0,1.d0)*tg_v(j,3))
           sdwn = tg_psic(j,2) * (tg_v(j,1)-tg_v(j,4)) + &
                  tg_psic(j,1) * (tg_v(j,2)+(0.d0,1.d0)*tg_v(j,3))
           tg_psic(j,1)=sup
           tg_psic(j,2)=sdwn
        enddo
     else
        do j=1, nrxxs
           sup = psic_nc(j,1) * (v(j,1)+v(j,4)) + &
                 psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
           sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + &
                  psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
           psic_nc(j,1)=sup
           psic_nc(j,2)=sdwn
        enddo
     endif
     !
     !   back to reciprocal space
     !
     if( use_tg ) then
        !
        do ipol = 1, npol

           call tg_cft3s ( tg_psic(:,ipol), dffts, -2, use_tg )
           !
           ioff   = 0
           !
           do idx = 1, nogrp
              !
              if( idx + ibnd - 1 <= m ) then
                 do j = 1, n
                    hpsi (j, ipol, ibnd+idx-1) = hpsi (j, ipol, ibnd+idx-1) + tg_psic( nls(igk(j)) + ioff, ipol )
                 enddo
              endif
              !
              ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
              !
           enddo

        enddo
        !
     else

        do ipol=1,npol
           call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
        enddo
        !
        !   addition to the total product
        !
        do ipol=1,npol
           do j = 1, n
              hpsi(j,ipol,ibnd) = hpsi(j,ipol,ibnd) + psic_nc(nls(igk(j)),ipol)
           enddo
        enddo

     endif

  enddo

  if( use_tg ) then
     !
     deallocate( tg_v )
     deallocate( tg_psic )
     !
  endif
  !
  return
END SUBROUTINE vloc_psi_nc
