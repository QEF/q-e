!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE para_const
  !----------------------------------------------------------------------------
  !
  SAVE
  !
  INTEGER, PARAMETER :: &
#if defined(__QK_USER__)
      maxproc  = 2048     !  maximum number of processors. make it big for 'big iron'.
#else
      maxproc  = 128      !  maximum number of processors
#endif
END MODULE para_const
!
!
!----------------------------------------------------------------------------
MODULE pfft
  !----------------------------------------------------------------------------
  !
  ! ... parallel fft information for the dense grid
  !
  USE para_const
  !
  SAVE
  !
  INTEGER :: &
      npp(maxproc),     &!  number of plane per processor
      ncp(maxproc),     &!  number of (density) columns per proc
      ncp0(maxproc),    &!  starting column for each processor
      ncplane,          &!  number of columns in a plane
      nct,              &!  total number of non-zero columns
      nxx                !  local fft data dim
  !
END MODULE pfft
!
!
!----------------------------------------------------------------------------
MODULE pffts
  !----------------------------------------------------------------------------
  !
  ! ... parallel fft information for the smooth grid
  !
  USE para_const
  !
  SAVE
  !
  INTEGER :: &
       nkcp(maxproc)     !  number of (wfs) columns per processor
  INTEGER :: &
       npps(maxproc),   &!  number of plane per processor
       ncps(maxproc),   &!  number of (density) columns per proc
       ncp0s(maxproc),  &!  starting column for each processor
       ncplanes,        &!  number of columns in a plane
       ncts,            &!  total number of non-zero columns
       nxxs              !  local fft data dim
  !
END MODULE pffts
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
  ! ...              uses SHMEM for the T3D/T3E, MPI otherwhise
  !
  USE mp_global, ONLY : intra_pool_comm, my_pool_id, nproc_pool, npool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER            :: dim
  REAL (DP)     :: ps(dim)
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, nbuf
  INTEGER, PARAMETER :: maxb = 10000
  !
# if (defined __SHMEM && defined __ALTIX) || (defined __SHMEM && defined __ORIGIN)
  INTEGER            :: sym_len
  LOGICAL            :: first
  REAL (DP)     :: buff(*), snd_buff(*)
  POINTER               (buff_p, buff), (snd_buff_p, snd_buff)
  COMMON /sym_heap1/    buff_p, snd_buff_p, sym_len, first
# else
  REAL (DP)     :: buff(maxb)  
# endif
  !
# if defined (__SHMEM)
  !
  ! ... SHMEM specific 
  !
  INCLUDE 'mpp/shmem.fh'
# if defined (__ALTIX) || defined (__ORIGIN)
      INTEGER    :: pWrkSync(SHMEM_REDUCE_SYNC_SIZE),                    &
     &              pWrkData(1024*1024), start
      DATA pWrkSync /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/
      DATA pWrkData / 1048576 * 0 /
# else
  INTEGER :: pWrkSync, pWrkData, start
  COMMON / SH_SYNC / pWrkSync(SHMEM_BARRIER_SYNC_dim)
  COMMON / SH_DATA / pWrkData(1024*1024)
  DATA pWrkData / 1048576 * 0 /
  DATA pWrkSync / SHMEM_BARRIER_SYNC_dim * SHMEM_SYNC_VALUE /
!DIR$ CACHE_ALIGN /SH_SYNC/
!DIR$ CACHE_ALIGN /SH_DATA/
#  endif
  !
#  endif
  !
  !
  IF ( dim <= 0 .OR. nproc_pool <= 1 ) RETURN
  !
  CALL start_clock( 'reduce' )
  !
  ! ... syncronize processes - maybe unneeded on T3D but necessary on T3E !!!
  !
  CALL mp_barrier( intra_pool_comm )
  !
  nbuf = dim / maxb
  !
#  if defined (__SHMEM)
# if defined (__ALTIX) || defined (__ORIGIN)
  IF (dim .GT. sym_len) THEN
     IF (sym_len .NE. 0) THEN
        CALL shpdeallc( snd_buff_p, info, -1 )
     END IF
     sym_len = dim
     CALL shpalloc( snd_buff_p, 2*sym_len, info, -1 )
  END IF
  IF (first .NE. .TRUE.) THEN
     CALL shpalloc( buff_p, 2*maxb, info, -1 )
     first = .TRUE.
  END IF
  snd_buff(1:dim) = ps(1:dim)
# endif
  !
  start = my_pool_id * nproc_pool
  !
#  endif
  !
  DO n = 1, nbuf
     !
#  if defined (__SHMEM)
     !
# if defined (__ALTIX) || defined (__ORIGIN)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, snd_buff(1+(n-1)*maxb), maxb, &
                                  start, 0, nproc_pool, pWrkData, pWrkSync )
# else
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+(n-1)*maxb), maxb, &
                                  start, 0, nproc_pool, pWrkData, pWrkSync )
#endif
     !                             
#  else
     !
     CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_REAL8, &
                         MPI_SUM, intra_pool_comm, info )
     !                    
     CALL errore( 'reduce', 'error in allreduce1', info )
     !
#  endif
     !
     ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
#  if defined (__SHMEM)
     !
# if defined (__ALTIX) || defined (__ORIGIN)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, snd_buff(1+nbuf*maxb),          &
     &                            (dim-nbuf*maxb), start, 0, nproc_pool,&
     &                            pWrkData, pWrkSync )
# else
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+nbuf*maxb), (dim-nbuf*maxb), &
                                  start, 0, nproc_pool, pWrkData, pWrkSync )
# endif
     !                             
#  else
     !
     CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_REAL8, &
                         MPI_SUM, intra_pool_comm, info )
     !
     CALL errore( 'reduce', 'error in allreduce2', info )
     !
#  endif
     !
     ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     !
  END IF
  !
  CALL stop_clock( 'reduce' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE reduce
!
!----------------------------------------------------------------------------
SUBROUTINE ireduce( dim, is )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable is(dim) over the processors.
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, npool
  USE mp,        ONLY : mp_barrier
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER            :: dim, is(dim)
  !
#if defined (__PARA)  
  !
  INTEGER            :: info, n, m, nbuf
  INTEGER, PARAMETER :: maxi = 500
  INTEGER            :: buff(maxi)
  !
  !
  IF ( dim <= 0 .OR. nproc_pool <= 1 ) RETURN
  !
  ! ... syncronize processes
  !
  CALL mp_barrier( intra_pool_comm )
  !
  nbuf = dim / maxi
  !
  DO n = 1, nbuf
     !
     CALL MPI_ALLREDUCE( is(1+(n-1)*maxi), buff, maxi, MPI_INTEGER, &
                         MPI_SUM, intra_pool_comm, info )
     !   
     CALL errore( 'ireduce', 'error in allreduce 1', info )
     !
     is((1+(n-1)*maxi):(n*maxi)) = buff(1:maxi)
     !
  END DO
  !
  ! ... possible remaining elements < maxi
  !
  IF ( ( dim - nbuf * maxi ) > 0 ) THEN
     !
     CALL MPI_ALLREDUCE( is(1+nbuf*maxi), buff, (dim-nbuf*maxi), MPI_INTEGER, &
                         MPI_SUM, intra_pool_comm, info )
     !
     CALL errore( 'reduce', 'error in allreduce 2', info )
     !
     is((1+nbuf*maxi):dim) = buff(1:(dim-nbuf*maxi))
     !
  END IF
  !
#endif
  !
  RETURN
  !
END SUBROUTINE ireduce
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
  INTEGER            :: dim
  REAL (DP)     :: ps(dim)
  !
#if defined (__PARA)  
  !
  INTEGER, PARAMETER :: maxb = 10000
  REAL (DP)     :: buff(maxb)
  INTEGER            :: info, nbuf, n
  !
  !
  IF ( dim <= 0 .OR. npool <= 1 ) RETURN
  !
  CALL start_clock( 'poolreduce' )
  !
  ! ... MPI syncronize processes
  !
  CALL mp_barrier( intra_image_comm )
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     CALL MPI_ALLREDUCE( ps(1+(n-1)*maxb), buff, maxb, MPI_REAL8, &
                         MPI_SUM, inter_pool_comm, info )
     !
     CALL errore( 'poolreduce', 'info<>0 at allreduce1', info )
     !
     ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     !
  END DO
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     CALL MPI_ALLREDUCE( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), MPI_REAL8, &
                         MPI_SUM, inter_pool_comm, info )
     !
     CALL errore( 'poolreduce', 'info<>0 at allreduce2', info )
     !
     ps((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxb))
     !
  END IF
  !
  CALL stop_clock( 'poolreduce' )
  !
#endif
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
  USE pfft,      ONLY : ncplane, npp, nxx   
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  REAL (DP) :: f_in(nxx), f_out(*)
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = ncplane * npp(proc+1)
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
  CALL MPI_GATHERV( f_in, recvcount(me_pool), MPI_REAL8, f_out, &
                    recvcount, displs, MPI_REAL8, root_pool,    &
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
  USE pfft,      ONLY : ncplane, npp, nxx     
  USE mp_global, ONLY : intra_pool_comm, intra_image_comm, &
                        nproc_pool, me_pool
  USE mp,        ONLY : mp_barrier
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(nxx), f_out(*)
  !
#if defined (__PARA)  
  !
  INTEGER          :: proc, info
  INTEGER          :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'cgather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = 2 * ncplane * npp(proc+1)
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
  CALL MPI_ALLGATHERV( f_in, recvcount(me_pool), MPI_REAL8, &
                       f_out, recvcount, displs, MPI_REAL8, &
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
!----------------------------------------------------------------------------
SUBROUTINE cgather_smooth ( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers data on the smooth AND complex fft grid
  !
  ! ... gathers nproc_pool distributed data on the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = distributed variable (nxxs)
  ! ... COMPLEX*16  f_out = gathered variable (nrx1s*nrx2s*nrx3s)
  !
  USE pffts,     ONLY : ncplanes, npps, nxxs   
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX (DP) :: f_in(nxxs), f_out(*)
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:nproc_pool-1), recvcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     recvcount(proc) = 2 * ncplanes * npps(proc+1)
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
  CALL MPI_GATHERV( f_in, recvcount(me_pool), MPI_REAL8, f_out, &
                    recvcount, displs, MPI_REAL8, root_pool,    &
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
  USE pfft,      ONLY : ncplane, npp, nxx 
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  REAL (DP) :: f_in(*), f_out(nxx)
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = ncplane * npp(proc+1)
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
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_REAL8,   &
                     f_out, sendcount(me_pool), MPI_REAL8, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /= nxx ) f_out(sendcount(me_pool)+1:nxx) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE scatter
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_smooth( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data on the smooth AND complex fft grid
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nrx1s*nrx2s*nrx3s)
  ! ... COMPLEX*16  f_out = distributed variable (nxxs)
  !
  USE pffts,     ONLY : ncplanes, npps, nxxs
  USE mp_global, ONLY : intra_pool_comm, nproc_pool, &
                        me_pool, root_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  COMPLEX (DP) :: f_in(*), f_out(nxxs)
  !
#if defined (__PARA)  
  !
  INTEGER        :: proc, info
  INTEGER        :: displs(0:nproc_pool-1), sendcount(0:nproc_pool-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = 2 * ncplanes * npps(proc+1)
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
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_REAL8,   &
                     f_out, sendcount(me_pool), MPI_REAL8, &
                     root_pool, intra_pool_comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(me_pool) /= nxxs ) f_out(sendcount(me_pool)+1:nxxs) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_smooth
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
  REAL (DP) :: f_in(nsize,nkstot), f_out(nsize,nks)
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
!----------------------------------------------------------------------------
SUBROUTINE fft_scatter1( f_in, nrx3, nxx_, f_aux, ncp_, npp_, sign )
  !----------------------------------------------------------------------------
  !
  ! ... transpose the fft grid across nodes
  !
  ! ... a)  From columns to planes (sign > 0)
  !
  ! ...     "columns" (or "pencil") representation:
  !
  ! ...        processor "me" has ncp_(me) contiguous columns along z
  ! ...        Each column has nrx3 elements for a fft of order nr3
  ! ...        nrx3 can be = nr3 + 1 in order to reduce memory conflicts.
  !
  ! ...     The transpose take places in two steps:
  !
  ! ...     1) on each processor the columns are divided into slices along z
  ! ...        that are stored contiguously. On processor "me", slices for
  ! ...        processor "proc" are npp_(proc)*ncp_(me) big
  ! ...     2) all processors communicate to exchange slices
  ! ...        (all columns with z in the slice belonging to "me"
  ! ...        must be received, all the others must be sent to "proc")
  !     
  ! ...     Finally one gets the "planes" representation:
  !
  ! ...        processor "me" has npp_(me) complete xy planes
  !
  ! ... b)  From planes to columns (sign < 0)
  !
  ! ...     Quite the same in the opposite direction
  !
  ! ...  The output is overwritten on f_in ; f_aux is used as work space
  !
  USE para_const, ONLY : maxproc
  USE mp_global,  ONLY : intra_pool_comm, nproc_pool, me_pool
  USE mp,         ONLY : mp_barrier
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER        :: nrx3, nxx_, sign, ncp_(maxproc), npp_(maxproc)
  REAL (DP) :: f_in(2*nxx_), f_aux(2*nxx_)
  !
#if defined (__PARA)  
  !
  INTEGER        :: dest, from, k, proc, ierr
  INTEGER        :: offset1(0:maxproc-1), sendcount(0:maxproc-1), &
                    sdispls(0:maxproc-1), recvcount(0:maxproc-1), &
                    rdispls(0:maxproc-1)
  !
  !
  IF ( nproc_pool == 1 ) RETURN
  !
  CALL start_clock( 'fft_scatter' )
  !
  ! ... sendcount(proc): amount of data processor "me" must send to processor
  ! ... recvcount(proc): amount of data processor "me" must receive from
  !
  DO proc = 0, ( nproc_pool - 1 )
     !
     sendcount(proc) = 2 * npp_(proc+1) * ncp_(me_pool+1)
     recvcount(proc) = 2 * npp_(me_pool+1) * ncp_(proc+1)
     !
  END DO
  !
  ! ... offset1(proc)   is used to locate the slices to be sent to proc
  ! ... sdispls(proc)+1 is the beginning of data that must be sent to proc
  ! ... rdispls(proc)+1 is the beginning of data that must be received from proc
  !
  offset1(0) = 1
  sdispls(0) = 0
  rdispls(0) = 0
  !
  DO proc = 1, ( nproc_pool - 1 )
     !
     offset1(proc) = offset1(proc-1) + 2 * npp_(proc)
     sdispls(proc) = sdispls(proc-1) + sendcount(proc-1)
     rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
     !
  END DO
  !
  ierr = 0
  !
  IF ( sign > 0 ) THEN
     !
     ! ... "forward" scatter from columns to planes
     !
     ! ... step one: store contiguously the slices
     !
     DO proc = 0, ( nproc_pool - 1 )
        !
        from = offset1(proc)
        dest = 1 + sdispls(proc)
        !
        DO k = 1, ncp_(me_pool+1)
           !
           CALL DCOPY( 2*npp_(proc+1), f_in(from+2*(k-1)*nrx3), &
                       1, f_aux(dest+2*(k-1)*npp_(proc+1)), 1 )
           !
        END DO
        !
     END DO
     !
     ! ... maybe useless; ensures that no garbage is present in the output
     !
     f_in(:) = ( 0.D0, 0.D0 )
     !
     ! ... step two: communication
     !
     CALL mp_barrier( intra_pool_comm )
     !
     CALL MPI_ALLTOALLV( f_aux, sendcount, sdispls, MPI_REAL8, &
                         f_in, recvcount, rdispls, MPI_REAL8,  &
                         intra_pool_comm, ierr )
     !     
     CALL errore( 'fft_scatter', 'info<>0', ierr )
     !
  ELSE
     !
     ! ... "backward" scatter from planes to columns
     !
     ! ... step two: communication
     !
     CALL mp_barrier( intra_pool_comm )
     !
     CALL MPI_ALLTOALLV( f_in, recvcount, rdispls, MPI_REAL8,  &
                         f_aux, sendcount, sdispls, MPI_REAL8, &
                         intra_pool_comm, ierr )
     !     
     CALL errore( 'fft_scatter', 'info<>0', ierr )
     !
     ! ... step one: store contiguously the columns
     !
     f_in(:) = ( 0.D0, 0.D0 )
     !
     DO proc = 0,( nproc_pool - 1 ) 
        !
        from = 1 + sdispls(proc)
        !
        dest = offset1(proc)
        !
        DO k = 1, ncp_(me_pool+1)
           !
           CALL DCOPY( 2*npp_(proc+1), f_aux(from+2*(k-1)*npp_(proc+1)), &
                       1, f_in (dest+2*(k-1)*nrx3), 1 )
           !
        END DO
        !
     END DO
     !
  END IF
  !
  CALL stop_clock( 'fft_scatter' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE fft_scatter1
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
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER        :: iflag
  REAL (DP) :: ps
  !
#if defined (__PARA)  
  !
  INTEGER        :: info
  REAL (DP) :: psr 
  !
  !
  IF ( npool <= 1 ) RETURN
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( iflag > 0 ) THEN
     !
     CALL MPI_ALLREDUCE( ps, psr, 1, MPI_REAL8, MPI_MAX, &
                         inter_pool_comm, info )
     !
     CALL errore( 'poolextreme', 'info<>0 in allreduce1', info )
     !
  ELSE
     !
     CALL MPI_ALLREDUCE( ps, psr, 1, MPI_REAL8, MPI_MIN, &
                         inter_pool_comm, info )
     !
     CALL errore( 'poolextreme', 'info<>0 in allreduce2', info )
     !
  END IF
  !
  ps = psr
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
  INTEGER        :: length, nks, nkstot
  REAL (DP) :: vec(length,nkstot)
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
     CALL MPI_SEND( vec, (length*nks), MPI_REAL8, 0, 17, &
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
        CALL MPI_RECV( vec(1,nbase+1), (length*fine), MPI_REAL8, &
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
  USE mp,        ONLY : mp_barrier  
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  REAL (DP) :: ps
  INTEGER        :: iflag
  !
#if defined (__PARA)  
  !
  REAL (DP) :: psr
  INTEGER        :: info
  !
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( iflag > 0 ) THEN
     !
     CALL MPI_ALLREDUCE( ps, psr, 1, MPI_REAL8, MPI_MAX, &
                         intra_image_comm, info )
     !
  ELSE
     !
     CALL MPI_ALLREDUCE( ps, psr, 1, MPI_REAL8, MPI_MIN, &
                         intra_image_comm, info )
     !
  END IF
  !
  ps = psr
  !
#endif
  !
  RETURN
  !
END SUBROUTINE extreme
