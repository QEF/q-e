!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!------------------------------------------------------------------------
SUBROUTINE poolreduce( dim, ps )
  !-----------------------------------------------------------------------
  !
  ! ... Sums a distributed variable ps(dim) over the pools.
  ! ... This MPI-only version uses a fixed-length buffer
  !       
#if defined (__PARA)
  !  
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, nproc_pool
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'
  INTEGER            :: dim
  REAL (KIND=DP)     :: ps(*)
  INTEGER, PARAMETER :: maxb = 10000
  REAL (KIND=DP)     :: buff(maxb)
  INTEGER            :: info, nbuf, n
  !
  !
  IF ( dim <= 0 .OR. nproc_pool <= 1 ) RETURN
  !
  CALL start_clock( 'poolreduce' )
  !
  ! ... MPI syncronize processes
  !
  CALL mp_barrier()
  !
  nbuf = dim / maxb
  !
  DO n = 1, nbuf
     !
     CALL MPI_allreduce( ps(1+(n-1)*maxb), buff, maxb, &
                         MPI_REAL8, MPI_SUM, inter_pool_comm, info )
     !
     CALL errore( 'poolreduce', 'info<>0 at allreduce1', info )
     !
     ps((1+(n-1)*maxb):(n*maxb)) = buff(1:maxb)
     !CALL DCOPY( maxb, buff, 1, ps(1+(n-1)*maxb), 1 )
     !
  END DO
  !
  IF ( ( dim - nbuf * maxb ) > 0 ) THEN
     !
     CALL MPI_allreduce( ps(1+nbuf*maxb), buff, (dim-nbuf*maxb), &
                         MPI_REAL8, MPI_SUM, inter_pool_comm, info )
     !
     CALL errore( 'poolreduce', 'info<>0 at allreduce2', info )
     !
     ps((1+nbuf*maxb):dim) = buff(1:dim-nbuf*maxb)
     !CALL DCOPY( dim-nbuf*maxb, buff, 1, ps(1+nbuf*maxb), 1 )
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
