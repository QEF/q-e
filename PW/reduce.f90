!
! Copyright (C) 2001-2204 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef SHMEM
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE reduce( dim, ps )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  ! ...              uses SHMEM for the T3D/T3E, MPI otherwhise
  !
#if defined (__PARA)
  !
  USE mp_global, ONLY : intra_pool_comm, my_pool_id, nproc_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'    
  INTEGER            :: dim
  REAL (KIND=DP)     :: ps(dim)
  INTEGER            :: info, n, nbuf
  INTEGER, PARAMETER :: maxb = 10000
  REAL (KIND=DP)     :: buff(maxb)  
  !
#  if defined (SHMEM)
  !
  ! ... SHMEM specific 
  !
  INCLUDE 'mpp/shmem.fh'
  INTEGER :: pWrkSync, pWrkData, start
  COMMON / SH_SYNC / pWrkSync(SHMEM_BARRIER_SYNC_dim)
  COMMON / SH_DATA / pWrkData(1024*1024)
  DATA pWrkData / 1048576 * 0 /
  DATA pWrkSync / SHMEM_BARRIER_SYNC_dim * SHMEM_SYNC_VALUE /
!DIR$ CACHE_ALIGN /SH_SYNC/
!DIR$ CACHE_ALIGN /SH_DATA/
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
#  if defined (SHMEM)
  !
  start = my_pool_id * nproc_pool
  !
#  endif
  !
  DO n = 1, nbuf
     !
#  if defined (SHMEM)
     !
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+(n-1)*maxb), &
                                  maxb, start, 0, nprocp, pWrkData, pWrkSync )
     !                             
#  else
     !
     CALL MPI_allreduce( ps(1+(n-1)*maxb), buff, maxb, &
                         MPI_REAL8, MPI_SUM, intra_pool_comm, info )
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
#  if defined (SHMEM)
     !
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+nbuf*maxb), (dim-nbuf*maxb), &
                                  start, 0, nprocp, pWrkData, pWrkSync )
     !                             
#  else
     !
     CALL MPI_allreduce( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), &
                         MPI_REAL8, MPI_SUM, intra_pool_comm, info )
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
!
!----------------------------------------------------------------------------
SUBROUTINE ireduce( dim, is )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable is(dim) over the processors.
  !
#if defined (__PARA)
  !
  USE mp_global, ONLY : intra_pool_comm, nproc_pool
  USE mp,        ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'
  !  
  INTEGER            :: dim, is(dim)
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
     CALL MPI_allreduce( is(1+(n-1)*maxi), buff, maxi, &
                         MPI_INTEGER, MPI_SUM, intra_pool_comm, info )
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
     CALL MPI_allreduce( is(1+nbuf*maxi), buff, (dim-nbuf*maxi), &
                         MPI_INTEGER, MPI_SUM, intra_pool_comm, info )
     !
     CALL errore( 'reduce', 'error in allreduce 2', info )
     !
     is((1+nbuf*maxb):dim) = buff(1:(dim-nbuf*maxi))
     !
  END IF
  !
#endif
  !
  RETURN
  !
END SUBROUTINE ireduce
