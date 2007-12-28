!
! Copyright (C) 2003 PWSCF group
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
  USE parallel_include
  USE kinds, only : DP
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE wvfct,   ONLY : igk
  USE wavefunctions_module,  ONLY: psic
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool, nolist
  USE fft_parallel,  ONLY : tg_cft3s
  USE sticks,        ONLY : dffts
  !
  implicit none
  !
  integer :: lda, n, m
  complex(DP) :: psi (lda, m), hpsi (lda, m)
  real(DP) :: v(nrxxs)
  !
  complex(DP) :: fp, fm
  integer :: i, ibnd, j, incr, ierr, idx, ioff, nsiz
  logical :: use_tg
  ! counters
  !
  ! Task Groups
  REAL(DP),    ALLOCATABLE :: tg_v(:)
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:)
  INTEGER :: recv_cnt( nogrp ), recv_displ( nogrp )
  INTEGER :: v_siz


  call start_clock ('vloc_psi')
  !
  incr = 2
  !
  use_tg = ( dffts%use_task_groups ) .AND. ( m >= nogrp )

  IF( use_tg ) THEN
     !
     ! call errore( ' vloc_psi ', ' task_groups not yet implemented ', 1 )
     !
     v_siz =  dffts%nnrx * ( nogrp + 1 )
     !
     ALLOCATE( tg_v   ( v_siz ) )
     ALLOCATE( tg_psic( dffts%nnrx * ( nogrp + 1 ) ) )
     !
     tg_v = 0.0d0
     !
     !  The potential in v is distributed accros all processors
     !  We need to redistribute it so that it is completely contained in the
     !  processors of an orbital TASK-GROUP
     !
     recv_cnt(1)   = dffts%npp( nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
     recv_displ(1) = 0
     DO i = 2, NOGRP
        recv_cnt(i) = dffts%npp( nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
        recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
     ENDDO

#if defined (__PARA) && defined (__MPI)
     nsiz = dffts%npp( me_pool+1 ) * dffts%nr1x * dffts%nr2x
     
     CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
          tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, ogrp_comm, IERR)
     IF( ierr /= 0 ) &
        call errore( ' vloc_psi ', ' MPI_Allgatherv ', ABS( ierr ) )
#endif

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
                 tg_psic(nls (igk(j))+ioff) =        psi(j,idx+ibnd-1) + (0.0d0,1.d0) * psi(j,idx+ibnd) 
                 tg_psic(nlsm(igk(j))+ioff) = CONJG( psi(j,idx+ibnd-1) - (0.0d0,1.d0) * psi(j,idx+ibnd) )
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
     IF( use_tg ) THEN
        !
        call tg_cft3s ( tg_psic, dffts, 2)

        do j = 1, nrx1s * nrx2s * dffts%tg_npp( me_pool + 1 )
           tg_psic (j) = tg_psic (j) * tg_v(j)
        enddo

        call tg_cft3s ( tg_psic, dffts, -2)
        !
     ELSE
        !
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        !
        !   product with the potential v on the smooth grid
        !
        do j = 1, nrxxs
           psic (j) = psic (j) * v(j)
        enddo
        !
        !   back to reciprocal space
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
  call stop_clock ('vloc_psi')
  !
  return
end subroutine vloc_psi
