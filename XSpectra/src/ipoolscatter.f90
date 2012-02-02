!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE ipoolscatter( nsize, nkstot, f_in, nks, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... This routine scatters a quantity ( typically the eigenvalues )
  ! ... among the pools.
  ! ... On input, f_in is required only on the first node of the first
  ! pool.
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
  INTEGER :: f_in(nsize,nkstot), f_out(nsize,nks)
    ! input  ( contains values for all k-point )
    ! output ( only for k-points of mypool )
  !
#if defined (__MPI)
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
END SUBROUTINE ipoolscatter
