!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!-----------------------------------------------------------------------
subroutine cgather_sym (f_in, f_out)
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (in phonon code)
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nrx1*nrx2*nrx3)
  !
#if defined (__PARA)
  !
  USE mp_global, ONLY : intra_pool_comm
  USE para,      ONLY : me, mypool
  USE mp,        ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  INCLUDE 'mpif.h'
  COMPLEX(KIND=DP) :: f_in(nxx), f_out(*)
  INTEGER          :: proc, info, displs(nprocp), recvcount(nprocp)
  !
  CALL start_clock( 'cgather' )
  !
  DO proc = 1, nprocp
     !
     recvcount(proc) = 2 * ncplane * npp(proc)
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
  CALL MPI_allgatherv( f_in, recvcount(me), MPI_REAL8, f_out, &
                       recvcount, displs, MPI_REAL8, intra_pool_comm, info )
  !
  CALL errore( 'cgather_sym', 'info<>0', info )
  !
  CALL mp_barrier()
  !
  CALL stop_clock( 'cgather' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym

