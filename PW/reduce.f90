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
SUBROUTINE reduce( size, ps )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(size) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) size
  ! ...              uses SHMEM for the T3D/T3E, MPI otherwhise
  !
#if defined (__PARA)
  USE para,  ONLY : MPI_COMM_POOL, mypool, nprocp
  USE mp,    ONLY : mp_barrier
#endif
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER        :: size
  REAL (KIND=DP) :: ps(size)
  !
#if defined (__PARA)
  !
  ! ... MPI specific
  !
  INCLUDE 'mpif.h'  
  INTEGER        :: info, n, nbuf
#define MAXB 10000
  REAL (KIND=DP) :: buff(MAXB)  
#if defined (SHMEM)
  !
  ! ... SHMEM specific 
  !
  INCLUDE 'mpp/shmem.fh'
  INTEGER :: pWrkSync, pWrkData, start
  COMMON / SH_SYNC / pWrkSync (SHMEM_BARRIER_SYNC_SIZE)
  COMMON / SH_DATA / pWrkData (1024 * 1024)
  DATA pWrkData / 1048576 * 0 /
  DATA pWrkSync / SHMEM_BARRIER_SYNC_SIZE * SHMEM_SYNC_VALUE /
!DIR$ CACHE_ALIGN /SH_SYNC/
!DIR$ CACHE_ALIGN /SH_DATA/
#endif
  !
  !
  IF ( nprocp <= 1 ) RETURN
  IF ( size <= 0 ) RETURN
  !
  CALL start_clock( 'reduce' )
  !
  ! ... syncronize processes - maybe unneeded on T3D but necessary on T3E !!!
  !
  CALL mp_barrier( MPI_COMM_POOL )
  !
  nbuf = size / MAXB
  !
#if defined (SHMEM)
  !
  start = ( mypool - 1 ) * nprocp
  !
#endif
  !
  DO n = 1, nbuf
     !
#if defined (SHMEM)
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+(n-1)*MAXB), &
                                  MAXB, start, 0, nprocp, pWrkData, pWrkSync )
#else
     CALL MPI_allreduce( ps(1+(n-1)*MAXB), buff, MAXB, &
                         MPI_REAL8, MPI_SUM, MPI_COMM_POOL, info )
     CALL errore( 'reduce', 'error in allreduce1', info )
#endif
     !
     ps((1+(n-1)*MAXB):(n*MAXB)) = buff(1:MAXB)
     !CALL DCOPY( MAXB, buff, 1, ps(1+(n-1)*MAXB), 1 )
     !
  END DO
  !
  ! ... possible remaining elements < maxb
  !
  IF ( ( size - nbuf * MAXB ) > 0 ) THEN
     !
#if defined (SHMEM)
     !
     CALL SHMEM_REAL8_SUM_TO_ALL( buff, ps(1+nbuf*MAXB), (size-nbuf*MAXB), &
                                  start, 0, nprocp, pWrkData, pWrkSync )
#else
     CALL MPI_allreduce( ps(1+nbuf*MAXB), buff, (size-nbuf*MAXB), &
                         MPI_REAL8, MPI_SUM, MPI_COMM_POOL, info )
     !
     CALL errore( 'reduce', 'error in allreduce2', info )
#endif
     !
     ps((1+nbuf*MAXB):(1+size)) = buff(1:(size-nbuf*MAXB))
     !CALL DCOPY( (size-nbuf*MAXB), buff, 1, ps(1+nbuf*MAXB), 1 )
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
SUBROUTINE ireduce( size, is )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable is(size) over the processors.
  !
#if defined (__PARA)
  !
  USE para,   ONLY : MPI_COMM_POOL, nprocp
  USE mp,     ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'
  !
#define MAXI 500  
  !  
  INTEGER :: size, is(size)
  INTEGER :: info, n, m, nbuf
  INTEGER :: buff(MAXI)
  !
  !
  IF ( nprocp <= 1 ) RETURN
  IF ( size <= 0 ) RETURN
  !
  ! ... syncronize processes
  !
  CALL mp_barrier( MPI_COMM_POOL )
  !
  nbuf = size / MAXI
  !
  do n = 1, nbuf
     !
     CALL MPI_allreduce( is(1+(n-1)*MAXI), buff, MAXI, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_POOL, info )
     !   
     CALL errore( 'ireduce', 'error in allreduce 1', info )
     !
     is((1+(n-1)*MAXI):(n*MAXI)) = buff(1:MAXI)
     !do m = 1, MAXI
     !   is (m + (n - 1) * MAXI) = buff (m)
     !enddo
     !
  END DO
  !
  ! ... possible remaining elements < MAXI
  !
  IF ( ( size - nbuf * MAXI ) > 0 ) THEN
     !
     CALL MPI_allreduce( is(1+nbuf*MAXI), buff, (size-nbuf*MAXI), &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_POOL, info )
     !
     CALL errore( 'reduce', 'error in allreduce 2', info )
     !
     is((1+nbuf*MAXB):(1+size)) = buff(1:(size-nbuf*MAXB))
     !do m = 1, size-nbuf * MAXI
     !   is (m + nbuf * MAXI) = buff (m)
     !enddo
     !
  END IF
  !
#endif
  !
  RETURN
  !
END SUBROUTINE ireduce
