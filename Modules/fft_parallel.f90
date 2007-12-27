!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
! This module contains drivers for parallel fft
!
MODULE fft_parallel
!
   IMPLICIT NONE
   SAVE
!
CONTAINS
!
!  Task groups driver
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dffts, sign )
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------
  !TASK GROUPS FFT ROUTINE.
  !Added: C. Bekas, Oct. 2005. Adopted from the CPMD code (A. Curioni)
  !Revised by Carlo Cavazzoni 2007.
  !
  ! ... sign = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT YET IMPLEMENTED WITH TASK GROUPS
  ! ... sign = +-2 : parallel 3d fft for wavefunctions
  !
  ! ... sign = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  ! ...              fft along z using pencils        (cft_1z)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... sign = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (sign>0) or before (sign<0) the fft on z direction
  !
  ! ...  Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE fft_base,   ONLY : fft_scatter
  USE kinds,      ONLY : DP
  USE mp_global,  only : me_pool, nproc_pool, ogrp_comm, npgrp, nogrp, &
                         intra_pool_comm
  USE fft_types,  ONLY : fft_dlay_descriptor
  USE task_groups
  USE parallel_include

  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: f( : ) 
  type (fft_dlay_descriptor), intent(in) :: dffts
  INTEGER,     INTENT(IN)    :: sign
  !
  INTEGER                    :: mc, i, j, ii, iproc, k
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  INTEGER                    :: idx, ierr
  COMPLEX(DP), ALLOCATABLE   :: yf(:), aux (:)
  COMPLEX(DP) :: stmp
  INTEGER, DIMENSION(NOGRP)  :: send_cnt, send_displ, recv_cnt, recv_displ
  INTEGER, SAVE :: ib = 1
  !
  !
  CALL start_clock( 'cft3s' )
  !
  n1  = dffts%nr1
  n2  = dffts%nr2
  n3  = dffts%nr3
  nx1 = dffts%nr1x
  nx2 = dffts%nr2x
  nx3 = dffts%nr3x
  !
  ALLOCATE( aux( (NOGRP+1)*strd ) )
  ALLOCATE( YF ( (NOGRP+1)*strd ) )
  !
  me_p = me_pool + 1
  !
  IF ( sign > 0 ) THEN
     !
     IF ( sign /= 2 ) THEN
        !
        CALL errore( ' tg_cfft ', ' task groups are implemented only for waves ', 1 )
        !
     ELSE
        !
        !
        send_cnt(1)   = nx3 * dffts%nsw( me_p )
        IF( nx3 * dffts%nsw( me_p ) > strd ) THEN
           CALL errore( ' tg_cfft ', ' inconsistent strd ', 1 )
        END IF
        send_displ(1) = 0
        recv_cnt(1)   = nx3 * dffts%nsw( nolist(1) + 1 )
        recv_displ(1) = 0
        DO idx = 2, nogrp
           send_cnt(idx)   = nx3 * dffts%nsw( me_p )
           send_displ(idx) = send_displ(idx-1) + strd
           recv_cnt(idx)   = nx3 * dffts%nsw( nolist(idx) + 1 )
           recv_displ(idx) = recv_displ(idx-1) + recv_cnt(idx-1)
        ENDDO

        IF( recv_displ(nogrp) + recv_cnt(nogrp) > SIZE( yf ) ) THEN
           CALL errore( ' tg_cfft ', ' inconsistent size ', 1 )
        END IF
        IF( send_displ(nogrp) + send_cnt(nogrp) > SIZE( f ) ) THEN
           CALL errore( ' tg_cfft ', ' inconsistent size ', 2 )
        END IF

        CALL start_clock( 'ALLTOALL' )
        !
        !  Collect all the sticks of the different states,
        !  in "yf" processors will have all the sticks of the OGRP

#if defined __MPI

        CALL MPI_ALLTOALLV( f(1), send_cnt, send_displ, MPI_DOUBLE_COMPLEX, yf(1), recv_cnt, &
         &                     recv_displ, MPI_DOUBLE_COMPLEX, ogrp_comm, IERR)
        IF( ierr /= 0 ) THEN
           CALL errore( ' tg_cfft ', ' alltoall error 1 ', ABS(ierr) )
        END IF

#endif

        !
        CALL cft_1z( yf, tmp_nsw( me_p ), n3, nx3, sign, aux )

        !
        CALL fft_scatter( aux, nx3, (nogrp+1)*strd, f, tmp_nsw, tmp_npp, sign, use_tg = .TRUE. )
        !
        f(:)      = ( 0.D0 , 0.D0 )
        ii        = 0
        !
        DO iproc = 1, nproc_pool
           !
           DO i = 1, dffts%nsw( iproc )
              !
              mc = dffts%ismap( i + dffts%iss(iproc) )
              !
              ii = ii + 1
              !
              DO j = 1, tmp_npp (me_p)
                 !
                 f(mc+(j-1)*nx1*nx2) = aux(j+(ii-1)*tmp_npp( me_p ) )
                 !
              END DO
              !
           END DO
           !
        END DO
        !
     END IF
     !
     CALL cft_2xy( f, tmp_npp( me_p ), n1, n2, nx1, nx2, sign, dffts%iplw )

     !
  ELSE
     !
     IF ( sign /= -2 ) THEN
        !
        CALL errore( ' tg_cfft ', ' task groups are implemented only for waves ', 1 )
        !
     ELSE
        !
        CALL cft_2xy( f, tmp_npp(me_p), n1, n2, nx1, nx2, sign, dffts%iplw )
        !
        ii = 0
        !
        DO iproc = 1, nproc_pool
           !
           DO i = 1, dffts%nsw(iproc)
              !
              mc = dffts%ismap( i + dffts%iss(iproc) )
              !
              ii = ii + 1
              !
              DO j = 1, tmp_npp(me_p)
                 !
                 aux(j+(ii-1)*tmp_npp(me_p)) = f( mc + (j-1)*nx1*nx2)
                 !
              END DO
              !
           END DO
           !
        END DO
        !
        CALL fft_scatter( aux, nx3, (NOGRP+1)*strd, f, tmp_nsw, tmp_npp, sign, use_tg = .TRUE. )
        !
        CALL cft_1z( aux, tmp_nsw(me_p), n3, nx3, sign, yf )
        !
        !  Bring pencils back to their original distribution
        !
        send_cnt  (1) = nx3 * dffts%nsw( nolist(1) + 1 )
        send_displ(1) = 0
        recv_cnt  (1) = nx3 * dffts%nsw( me_p )
        recv_displ(1) = 0
        DO idx = 2, NOGRP
           send_cnt  (idx) = nx3 * dffts%nsw( nolist(idx) + 1 )
           send_displ(idx) = send_displ(idx-1) + send_cnt(idx-1)
           recv_cnt  (idx) = nx3 * dffts%nsw( me_p )
           recv_displ(idx) = recv_displ(idx-1) + recv_cnt(idx-1)
        ENDDO

        IF( recv_displ(nogrp) + recv_cnt(nogrp) > SIZE( f ) ) THEN
           CALL errore( ' tg_cfft ', ' inconsistent size ', 3 )
        END IF
        IF( send_displ(nogrp) + send_cnt(nogrp) > SIZE( yf ) ) THEN
           CALL errore( ' tg_cfft ', ' inconsistent size ', 4 )
        END IF

        CALL start_clock( 'ALLTOALL' )

#if defined __MPI
        CALL MPI_Alltoallv( yf(1), &
             send_cnt, send_displ, MPI_DOUBLE_COMPLEX, f(1), &
             recv_cnt, recv_displ, MPI_DOUBLE_COMPLEX, ogrp_comm, IERR)
        IF( ierr /= 0 ) THEN
           CALL errore( ' tg_cfft ', ' alltoall error 2 ', ABS(ierr) )
        END IF
#endif   

        CALL stop_clock( 'ALLTOALL' )

     END IF
     !
  END IF
  !
  DEALLOCATE( aux )
  DEALLOCATE( yf )
  !
  CALL stop_clock( 'cft3s' )
  !
  RETURN
  !
END SUBROUTINE tg_cft3s
!
END MODULE fft_parallel
