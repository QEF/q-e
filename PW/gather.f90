!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE gather( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers nprocp distributed data on the first processor of every pool
  !
  ! ... REAL*8  f_in  = distributed variable (nxx)
  ! ... REAL*8  f_out = gathered variable (nrx1*nrx2*nrx3)
  !
#if defined (__PARA)
  !
  USE pfft,      ONLY : ncplane, npp, nxx   
  USE mp_global, ONLY : intra_pool_comm
  USE io_global, ONLY : ionode_id
  USE para,      ONLY : me, nprocp
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'  
  REAL (KIND=DP) :: f_in(nxx), f_out(*)
  INTEGER        :: proc, info, displs(nprocp), recvcount(nprocp)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 1, nprocp
     !
     recvcount(proc) = ncplane * npp(proc)
     !
     IF ( proc == 1 ) THEN
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
  CALL MPI_gatherv( f_in, recvcount(me), MPI_REAL8, f_out, recvcount, &
                    displs, MPI_REAL8, ionode_id, intra_pool_comm, info )
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
