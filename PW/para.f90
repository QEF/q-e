!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!
! ... here are all parallel subroutines (wrappers to MPI calls) used 
! ... by the PWscf code
!
! ... "reduce"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE reduce( dim, ps )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM if available, MPI otherwhise
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), INTENT(INOUT) :: ps(dim)
  !
  IF ( dim <= 0 .OR. nproc_pool <= 1 ) RETURN
  !
  CALL start_clock( 'reduce' )

  CALL reduce_base_real( dim, ps, intra_pool_comm, -1 )

  CALL stop_clock( 'reduce' )
  !
  RETURN
  !
END SUBROUTINE reduce
!
!------------------------------------------------------------------------
SUBROUTINE poolreduce( dim, ps )
  !-----------------------------------------------------------------------
  !
  ! ... Sums a distributed variable ps(dim) over the pools.
  ! ... This MPI-only version uses a fixed-length buffer
  !
  USE mp_global, ONLY : inter_pool_comm, intra_image_comm, &
                        my_pool_id, nproc_pool, npool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP), INTENT(INOUT) :: ps(dim)
  !
  IF ( dim <= 0 .OR. npool <= 1 ) RETURN
  !
  CALL start_clock( 'poolreduce' )
  !
  CALL mp_barrier( intra_image_comm )  !  WHY on image? carlo c.
  !
  CALL reduce_base_real( dim, ps, inter_pool_comm, -1 )
  !
  CALL stop_clock( 'poolreduce' )
  !
  RETURN
  !
END SUBROUTINE poolreduce
!
! ... "gather"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE gather( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers nproc_pool distributed data on the first processor of every pool
  !
  ! ... REAL*8  f_in  = distributed variable (nxx)
  ! ... REAL*8  f_out = gathered variable (nrx1*nrx2*nrx3)
  !
  USE fft_base,  ONLY : dfftp
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in( dfftp%nnr ), f_out(*)
  !
#if defined (__PARA)  
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(me_pool), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, root_pool,    &
                    intra_pool_comm, info )
  !
  CALL errore( 'gather', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE gather
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym( f_in, f_out )
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (in phonon code)
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nrx1*nrx2*nrx3)
  !
  USE fft_base,  ONLY : dfftp
  USE mp_global, ONLY : intra_pool_comm, intra_image_comm, &
                        nproc_pool, me_pool
  USE mp,        ONLY : mp_barrier
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in( dfftp%nnr ), f_out(*)
  !
#if defined (__PARA)  
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'cgather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_ALLGATHERV( f_in, recvcount(me_pool), MPI_DOUBLE_PRECISION, &
                       f_out, recvcount, displs, MPI_DOUBLE_PRECISION, &
                       intra_pool_comm, info )
  !
  CALL errore( 'cgather_sym', 'info<>0', info )
  !
  CALL mp_barrier( intra_image_comm )
  !
  CALL stop_clock( 'cgather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym
!
!----------------------------------------------------------------------------
SUBROUTINE cgather_smooth ( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers data on the smooth AND complex fft grid
  !
  ! ... gathers nproc_pool distributed data on the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = distributed variable ( dffts%nnr )
  ! ... COMPLEX*16  f_out = gathered variable (nrx1s*nrx2s*nrx3s)
  !
  USE fft_base,  ONLY : dffts
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in( dffts%nnr ), f_out(*)
  !
#if defined (__PARA)  
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = 2 * dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( intra_pool_comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(me_pool), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, root_pool,    &
                    intra_pool_comm, info )
  !
  CALL errore( 'gather', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_smooth
!
! ... "scatter"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE scatter( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... REAL*8  f_in  = gathered variable (nrx1*nrx2*nrx3)
  ! ... REAL*8  f_out = distributed variable (nxx)
  !
  USE fft_base,  ONLY : dfftp
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in(*), f_out( dfftp%nnr )
  !
#if defined (__PARA)  
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( intra_pool_comm )
  !  
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(me_pool), MPI_DOUBLE_PRECISION, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /= dfftp%nnr ) f_out(sendcount(me_pool)+1:dfftp%nnr) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE scatter

!----------------------------------------------------------------------------
SUBROUTINE cscatter_sym( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nrx1*nrx2*nrx3)
  ! ... COMPLEX*16  f_out = distributed variable (nxx)
  !
  USE fft_base,  ONLY : dfftp
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(*), f_out( dfftp%nnr )
  !
#if defined (__PARA)  
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'cscatter_sym' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( intra_pool_comm )
  !  
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(me_pool), MPI_DOUBLE_PRECISION, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'cscatter_sym', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /=  dfftp%nnr  ) f_out(sendcount(me_pool)+1: dfftp%nnr ) = 0.D0
  !
  CALL stop_clock( 'cscatter_sym' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_sym
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_smooth( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data on the smooth AND complex fft grid
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nrx1s*nrx2s*nrx3s)
  ! ... COMPLEX*16  f_out = distributed variable ( dffts%nnr)
  !
  USE fft_base,  ONLY : dffts
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(*), f_out( dffts%nnr )
  !
#if defined (__PARA)  
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = 2 * dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     END IF
     !
  END DO
  !
  CALL mp_barrier( intra_pool_comm )
  !  
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(me_pool), MPI_DOUBLE_PRECISION, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /=  dffts%nnr  ) f_out(sendcount(me_pool)+1: dffts%nnr ) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_smooth
!
!----------------------------------------------------------------------------
SUBROUTINE poolscatter( nsize, nkstot, f_in, nks, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... This routine scatters a quantity ( typically the eigenvalues )
  ! ... among the pools. 
  ! ... On input, f_in is required only on the first node of the first pool. 
  ! ... f_in and f_out may coincide.
  ! ... Not a smart implementation!
  !
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, &
                        my_pool_id, npool, me_pool, root_pool, kunit
  USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  INTEGER :: nsize, nkstot, nks
    ! first dimension of vectors f_in and f_out
    ! number of k-points per pool
    ! total number of k-points
  REAL(DP) :: f_in(nsize,nkstot), f_out(nsize,nks)
    ! input  ( contains values for all k-point )
    ! output ( only for k-points of mypool )
  !
#if defined (__PARA)  
  !
  INTEGER :: rest, nbase
    ! the rest of the integer division nkstot / npo
    ! the position in the original list
  !
  !
  ! ... copy from the first node of the first pool
  ! ... to the first node of all the other pools
  !
  IF ( me_pool == root_pool ) &
     CALL mp_bcast( f_in, root_pool, inter_pool_comm )
  !
  ! ... distribute the vector on the first node of each pool
  !
  rest = nkstot / kunit - ( nkstot / kunit / npool ) * npool 
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  f_out(:,1:nks) = f_in(:,(nbase+1):(nbase+nks))
  !
  ! ... copy from the first node of every pool
  ! ... to the other nodes of every pool
  !
  CALL mp_bcast( f_out, root_pool, intra_pool_comm )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE poolscatter
!
! ... other parallel subroutines
!
!------------------------------------------------------------------------
SUBROUTINE poolextreme( ps, iflag )
  !-----------------------------------------------------------------------
  !
  ! ... Finds the maximum (iflag.gt.0) or the minimum (iflag.le.0) value o
  ! ... a real variable among the values distributed across the pools and
  ! ... returns this value to all pools.
  !
  USE mp_global, ONLY : inter_pool_comm, intra_image_comm, npool
  USE mp,        ONLY : mp_barrier  
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER  :: iflag
  REAL(DP) :: ps
  !
#if defined (__PARA)  
  !
  IF ( npool <= 1 ) RETURN
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( iflag > 0 ) THEN
     !
     CALL parallel_max_real( 1, ps, inter_pool_comm, -1 )
     !
  ELSE
     !
     CALL parallel_min_real( 1, ps, inter_pool_comm, -1 )
     !
  END IF
  !
#endif
  !
  RETURN
  !
END SUBROUTINE poolextreme
!
!-----------------------------------------------------------------------
SUBROUTINE poolrecover( vec, length, nkstot, nks )
  !----------------------------------------------------------------------- 
  !
  ! ... recovers on the first processor of the first pool a 
  ! ... distributed vector
  !
  USE mp_global, ONLY : inter_pool_comm, intra_image_comm, &
                        npool, me_pool, root_pool, my_pool_id, kunit
  USE mp,        ONLY : mp_barrier  
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER  :: length, nks, nkstot
  REAL(DP) :: vec(length,nkstot)
  !
#if defined (__PARA)  
  !
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: i, nks1, rest, fine, nbase, info
  !
  !
  IF ( npool <= 1 ) RETURN
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) &
       CALL errore( 'poolrecover', 'nkstot/kunit is not an integer', nkstot )
  !
  nks1 = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( me_pool == root_pool .AND. my_pool_id > 0 ) THEN
     !
     CALL MPI_SEND( vec, (length*nks), MPI_DOUBLE_PRECISION, 0, 17, &
                    inter_pool_comm, info )
     !     
     CALL errore( 'poolrecover', 'info<>0 in send', info )
     !
  END IF
  !
  DO i = 2, npool
     !
     IF ( i <= rest ) THEN
        !
        fine = nks1 + kunit
        !
        nbase = ( nks1 + kunit ) * ( i - 1 )
        !
     ELSE
        !
        fine = nks1
        !
        nbase = rest * (nks1 + kunit) + (i - 1 - rest) * nks1
        !
     END IF
     !
     IF ( me_pool == root_pool .AND. my_pool_id == 0 ) THEN
        !
        CALL MPI_RECV( vec(1,nbase+1), (length*fine), MPI_DOUBLE_PRECISION, &
                       (i-1), 17, inter_pool_comm, status, info )
        !
        CALL errore( 'poolrecover', 'info<>0 in recv', info )
        !
     END IF
     !
  END DO
  !
#endif
  !
  RETURN
  !
END SUBROUTINE poolrecover
!
!------------------------------------------------------------------------
SUBROUTINE ipoolrecover( ivec, length, nkstot, nks )
  !------------------------------------------------------------------------
  !
  ! ... as above, for an integer vector
  !
  USE mp_global, ONLY : inter_pool_comm, intra_image_comm, &
                        npool, me_pool, root_pool, my_pool_id, kunit
  USE mp,        ONLY : mp_barrier  
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER :: length, nks, nkstot
  INTEGER :: ivec(length,nkstot)
  !
#if defined (__PARA)  
  !
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: i, nks1, rest, fine, nbase, info
  !
  !
  IF ( npool <= 1 ) RETURN
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) & 
       CALL errore( 'poolrecover', 'nkstot/kunit is not an integer', nkstot )
  !
  nks1 = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( me_pool == root_pool .AND. my_pool_id > 0 ) THEN
     !
     CALL MPI_SEND( ivec, (length*nks), MPI_INTEGER, 0, 17, &
                    inter_pool_comm, info )
     !
     CALL errore( 'ipoolrecover', 'info<>0 in send', info )
     !
  END IF
  !
  DO i = 2, npool
     !
     IF ( i <= rest ) THEN
        !
        fine = nks1 + kunit
        !
        nbase = ( nks1 + kunit ) * ( i - 1 )
        !
     ELSE
        !
        fine = nks1
        !
        nbase = rest * ( nks1 + kunit ) + ( i - 1 - rest ) * nks1
        !
     END IF
     !
     IF ( me_pool == root_pool .AND. my_pool_id == 0 ) THEN
        !
        CALL MPI_RECV( ivec(1,nbase+1), (length*fine), MPI_INTEGER, &
                       (i-1), 17, inter_pool_comm, status, info )
        !
        CALL errore( 'ipoolrecover', 'info<>0 in recv', info )
        !
     END IF
     !
  END DO
  !
#endif
  !
  RETURN
  !
END SUBROUTINE ipoolrecover
!
!-----------------------------------------------------------------------
SUBROUTINE extreme( ps, iflag )
  !-----------------------------------------------------------------------
  !
  ! ... Finds the maximum (iflag.gt.0) or the minimum (iflag.le.0) value
  ! ... of a real variable among the values distributed on a given pool
  !
  USE mp_global, ONLY : intra_image_comm
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: ps
  INTEGER  :: iflag
  !
#if defined (__PARA)  
  !
  IF ( iflag > 0 ) THEN
     !
     CALL parallel_max_real( 1, ps, intra_image_comm, -1 )
     !
  ELSE
     !
     CALL parallel_min_real( 1, ps, intra_image_comm, -1 )
     !
  END IF
  !
#endif
  !
  RETURN
  !
END SUBROUTINE extreme
