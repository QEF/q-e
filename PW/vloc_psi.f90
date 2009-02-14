!
! Copyright (C) 2003-2009 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine vloc_psi(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - Gamma point
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE wvfct,   ONLY : igk
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
  USE fft_parallel,  ONLY : tg_cft3s
  USE fft_base,      ONLY : dffts
  USE control_flags, ONLY : use_task_groups
  USE task_groups,   ONLY : tg_gather
  USE wavefunctions_module,  ONLY: psic
  !
  implicit none
  !
  integer, intent(IN) :: lda, n, m
  complex(DP), intent(in)   :: psi (lda, m)
  complex(DP), intent(inout):: hpsi (lda, m)
  real(DP), intent(in) :: v(nrxxs)
  !
  integer :: i, ibnd, j, incr
  complex(DP) :: fp, fm
  !
  logical :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff, nsiz
  !
  !
  call start_clock ('h_psi:vloc')
  !
  incr = 2
  !
  use_tg = ( use_task_groups ) .AND. ( m >= nogrp )
  !
  IF( use_tg ) THEN
     !
     v_siz =  dffts%nnrx * nogrp
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
     ! 
     incr = 2 * nogrp
     !
  END IF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  do ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, 2*nogrp, 2
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1) + &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) 
                 tg_psic(nlsm(igk(j))+ioff) = CONJG( psi(j,idx+ibnd-1) - &
                                      (0.0d0,1.d0) * psi(j,idx+ibnd) )
              END DO
           ELSE IF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1)
                 tg_psic(nlsm(igk(j))+ioff) = CONJG( psi(j,idx+ibnd-1) )
              END DO
           END IF

           ioff = ioff + dffts%nnrx

        END DO
        !
     ELSE
        !
        psic(:) = (0.d0, 0.d0)
        if (ibnd < m) then
           ! two ffts at the same time
           do j = 1, n
              psic (nls (igk(j))) =       psi(j, ibnd) + (0.0d0,1.d0)*psi(j, ibnd+1)
              psic (nlsm(igk(j))) = CONJG(psi(j, ibnd) - (0.0d0,1.d0)*psi(j, ibnd+1))
           enddo
        else
           do j = 1, n
              psic (nls (igk(j))) =       psi(j, ibnd)
              psic (nlsm(igk(j))) = CONJG(psi(j, ibnd))
           enddo
        end if
        !
     END IF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        call tg_cft3s ( tg_psic, dffts, 2, use_tg )
        ! 
        do j = 1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        enddo
        !
        call tg_cft3s ( tg_psic, dffts, -2, use_tg )
        !
     ELSE
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        !
        do j = 1, nrxxs
           psic (j) = psic (j) * v(j)
        enddo
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
        !
     END IF
     !
     !   addition to the total product
     !
     IF( use_tg ) THEN
        !
        ioff   = 0
        !
        DO idx = 1, 2*nogrp, 2
           !
           IF( idx + ibnd - 1 < m ) THEN
              DO j = 1, n
                 fp= ( tg_psic( nls(igk(j)) + ioff ) +  tg_psic( nlsm(igk(j)) + ioff ) ) * 0.5d0
                 fm= ( tg_psic( nls(igk(j)) + ioff ) -  tg_psic( nlsm(igk(j)) + ioff ) ) * 0.5d0
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + CMPLX( DBLE(fp), AIMAG(fm))
                 hpsi (j, ibnd+idx  ) = hpsi (j, ibnd+idx  ) + CMPLX(AIMAG(fp),- DBLE(fm))
              END DO
           ELSE IF( idx + ibnd - 1 == m ) THEN
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1) + tg_psic( nls(igk(j)) + ioff )
              END DO
           END IF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
           !
        END DO
        !
     ELSE
        if (ibnd < m) then
           ! two ffts at the same time
           do j = 1, n
              fp = (psic (nls(igk(j))) + psic (nlsm(igk(j))))*0.5d0
              fm = (psic (nls(igk(j))) - psic (nlsm(igk(j))))*0.5d0
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + CMPLX( DBLE(fp), AIMAG(fm))
              hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + CMPLX(AIMAG(fp),- DBLE(fm))
           enddo
        else
           do j = 1, n
              hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
           enddo
        end if
     END IF
     !
  enddo
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  END IF
  !
  call stop_clock ('h_psi:vloc')
  !
  return
end subroutine vloc_psi
!
!-----------------------------------------------------------------------
subroutine vloc_psi_k(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - k-points
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE wvfct,   ONLY : igk
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
  USE fft_parallel,  ONLY : tg_cft3s
  USE fft_base,      ONLY : dffts
  USE control_flags, ONLY : use_task_groups
  USE task_groups,   ONLY : tg_gather
  USE wavefunctions_module,  ONLY: psic
  !
  implicit none
  !
  integer, intent(IN) :: lda, n, m
  complex(DP), intent(in)   :: psi (lda, m)
  complex(DP), intent(inout):: hpsi (lda, m)
  real(DP), intent(in) :: v(nrxxs)
  !
  integer :: i, ibnd, j, incr
  !
  logical :: use_tg
  ! Task Groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: v_siz, idx, ioff, nsiz
  !
  !
  call start_clock ('h_psi:vloc')
  !
  use_tg = ( use_task_groups ) .AND. ( m >= nogrp )
  !
  incr = 1
  !
  IF( use_tg ) THEN
     !
     v_siz =  dffts%nnrx * nogrp
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( v_siz ) )
     !
     CALL tg_gather( dffts, v, tg_v )
     incr = nogrp
     !
  END IF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  do ibnd = 1, m, incr
     !
     IF( use_tg ) THEN
        !
        tg_psic = (0.d0, 0.d0)
        ioff   = 0
        !
        DO idx = 1, nogrp

           IF( idx + ibnd - 1 <= m ) THEN
              DO j = 1, n
                 tg_psic(nls (igk(j))+ioff) =  psi(j,idx+ibnd-1)
              END DO
           END IF

           ioff = ioff + dffts%nnrx

        END DO
        !
        call tg_cft3s ( tg_psic, dffts, 2, use_tg )
        !
     ELSE
        !
        psic(:) = (0.d0, 0.d0)
        psic (nls (igk(1:n))) = psi(1:n, ibnd)
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        !
     END IF
     !
     !   fft to real space
     !   product with the potential v on the smooth grid
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        do j = 1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        enddo
        !
        call tg_cft3s ( tg_psic, dffts, -2, use_tg )
        !
     ELSE
        !
        do j = 1, nrxxs
           psic (j) = psic (j) * v(j)
        enddo
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
        !
     END IF
     !
     !   addition to the total product
     !
     IF( use_tg ) THEN
        !
        ioff   = 0
        !
        DO idx = 1, nogrp
           !
           IF( idx + ibnd - 1 <= m ) THEN
              DO j = 1, n
                 hpsi (j, ibnd+idx-1) = hpsi (j, ibnd+idx-1)
              END DO
           END IF
           !
           ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
           !
        END DO
        !
     ELSE
        do j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
        enddo
     END IF
     !
  enddo
  !
  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_psic )
     DEALLOCATE( tg_v )
     !
  END IF
  !
  call stop_clock ('h_psi:vloc')
  !
  return
end subroutine vloc_psi_k
!
!-----------------------------------------------------------------------
subroutine vloc_psi_nc (lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Calculation of Vloc*psi using dual-space technique - noncolinear
  !
  USE parallel_include
  USE kinds,   ONLY : DP
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE gvect,   ONLY : nrxx
  USE wvfct,   ONLY : igk
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
  USE fft_parallel,  ONLY : tg_cft3s
  USE fft_base,      ONLY : dffts
  USE control_flags, ONLY : use_task_groups
  USE task_groups,   ONLY : tg_gather
  USE noncollin_module,     ONLY: npol
  USE wavefunctions_module, ONLY: psic_nc
  !
  implicit none
  !
  integer, intent(IN) :: lda, n, m
  real(DP), intent(in) :: v(nrxx,4)
  complex(DP), intent(in)   :: psi (lda*npol, m)
  complex(DP), intent(inout):: hpsi (lda,npol,m)
  !
  integer :: i, ibnd, j,ipol, incr
  complex(DP) :: sup, sdwn
  !
  LOGICAL :: use_tg
  ! Variables for task groups
  REAL(DP),    ALLOCATABLE :: tg_v(:,:) 
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:,:) 
  INTEGER :: v_siz, idx, ioff, nsiz
  !
  !
  call start_clock ('h_psi:vloc')
  !
  incr = 1    
  !              
  use_tg = ( use_task_groups ) .AND. ( m >= nogrp )
  !        
  IF( use_tg ) THEN
     v_siz = dffts%nnrx * nogrp
     ALLOCATE( tg_v( v_siz, 4 ) )
     CALL tg_gather( dffts, v(:,1), tg_v(:,1) )
     CALL tg_gather( dffts, v(:,2), tg_v(:,2) )
     CALL tg_gather( dffts, v(:,3), tg_v(:,3) )
     CALL tg_gather( dffts, v(:,4), tg_v(:,4) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = nogrp
  END IF
  !
  ! the local potential V_Loc psi. First the psi in real space
  !
  do ibnd = 1, m, incr

     IF( use_tg ) THEN
        !
        DO ipol = 1, npol
           !
           tg_psic(:,ipol) = ( 0.D0, 0.D0 )
           ioff   = 0
           ! 
           DO idx = 1, nogrp
              !
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    tg_psic( nls( igk(j) ) + ioff, ipol ) = psi( j +(ipol-1)*lda, idx+ibnd-1 )
                 END DO
              END IF
  
              ioff = ioff + dffts%nnrx
     
           END DO 
           ! 
           call tg_cft3s ( tg_psic(:,ipol), dffts, 2, use_tg )
           !
        END DO
        !
     ELSE
        psic_nc = (0.d0,0.d0)
        do ipol=1,npol
           do j = 1, n
              psic_nc(nls(igk(j)),ipol) = psi(j+(ipol-1)*lda,ibnd)
           enddo
           call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        enddo
     END IF

     !
     !   product with the potential v = (vltot+vr) on the smooth grid
     !
     IF( use_tg ) THEN
        do j=1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           sup = tg_psic(j,1) * (tg_v(j,1)+tg_v(j,4)) + &
                 tg_psic(j,2) * (tg_v(j,2)-(0.d0,1.d0)*tg_v(j,3))
           sdwn = tg_psic(j,2) * (tg_v(j,1)-tg_v(j,4)) + &
                  tg_psic(j,1) * (tg_v(j,2)+(0.d0,1.d0)*tg_v(j,3))
           tg_psic(j,1)=sup
           tg_psic(j,2)=sdwn
        end do
     ELSE
        do j=1, nrxxs
           sup = psic_nc(j,1) * (v(j,1)+v(j,4)) + &
                 psic_nc(j,2) * (v(j,2)-(0.d0,1.d0)*v(j,3))
           sdwn = psic_nc(j,2) * (v(j,1)-v(j,4)) + &
                  psic_nc(j,1) * (v(j,2)+(0.d0,1.d0)*v(j,3))
           psic_nc(j,1)=sup
           psic_nc(j,2)=sdwn
        end do
     END IF
     !
     !   back to reciprocal space
     !
     IF( use_tg ) THEN
        !
        DO ipol = 1, npol

           call tg_cft3s ( tg_psic(:,ipol), dffts, -2, use_tg )
           !
           ioff   = 0
           !  
           DO idx = 1, nogrp
              !  
              IF( idx + ibnd - 1 <= m ) THEN
                 DO j = 1, n
                    hpsi (j, ipol, ibnd+idx-1) = hpsi (j, ipol, ibnd+idx-1) + tg_psic( nls(igk(j)) + ioff, ipol )
                 END DO
              END IF
              !
              ioff = ioff + dffts%nr3x * dffts%nsw( me_pool + 1 )
              !
           END DO

        END DO
        !
     ELSE

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

     END IF
 
  enddo

  IF( use_tg ) THEN
     !
     DEALLOCATE( tg_v )
     DEALLOCATE( tg_psic )
     !
  END IF
  !
  call stop_clock ('h_psi:vloc')
  !
  return
end subroutine vloc_psi_nc
