!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE poolscatter( nsize, nkstot, f_in, nks, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... This routine scatters a quantity (typically the eigenvalues)
  ! ... among the pools. On input, f_in is required only on the
  ! ... first node of the first pool. f_in and f_out may coincide.
  ! ... Not a smart implementation!
  !  
#if defined (__PARA)
  !
  USE kinds,      ONLY : DP
  USE para,       ONLY : me, mypool, npool, kunit, MPI_COMM_POOL, MPI_COMM_ROW
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  INTEGER :: nsize, nkstot, nks
    ! first dimension of vectors f_in and f_out
    ! number of k-points per pool
    ! total number of k-points
  REAL (KIND=DP) :: f_in(nsize,nkstot), f_out(nsize,nks)
    ! input (contains values for all k-point
    ! output(only for k-points of mypool)
  INTEGER :: rest, nbase
    ! the rest of the integer division nkstot/npo
    ! the position in the original list
  !
  !
  ! ... copy from the first node of the first pool
  ! ... to the first node of all the other pools
  !
  IF ( me == 1 ) CALL mp_bcast( f_in, ionode_id, MPI_COMM_ROW )
  !
  ! ... distribute the vector on the first node of each pool
  !
  rest = nkstot / kunit - ( nkstot / kunit / npool ) * npool 
  !
  nbase = nks * ( mypool - 1 )
  !
  IF ( mypool > rest ) nbase = nbase + rest * kunit
  !
  CALL DCOPY( nsize * nks, f_in(1,nbase+1), 1, f_out, 1 )
  !
  ! ... copy from the first node  of every pool
  ! ... to the other nodes of every pool
  !
  CALL mp_bcast( f_out, ionode_id, MPI_COMM_POOL )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE poolscatter
