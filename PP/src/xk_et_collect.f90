!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE xk_et_collect( xk_collect, et_collect, xk, et, nkstot, nks, nbnd )
  !----------------------------------------------------------------------------
  !
  ! ... This routine collects the k points (with granularity kunit) among 
  ! ... nodes and sets the variable xk_collect and wk_collect with the total 
  ! ... number of k-points
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : my_pool_id, npool, kunit, inter_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: nkstot, nks, nbnd
    ! total number of k-points
    ! number of k-points per pool
  REAL (DP) :: xk(3,nks), et(nbnd,nks)
  REAL (DP) :: xk_collect(3,nkstot), et_collect(nbnd,nkstot)
    ! k-points
    ! k-point weights
  !
#if defined (__MPI)
  !
  INTEGER :: nbase, rest, nks1
  !
  xk_collect=0.d0
  !
  et_collect=0.d0
  !
  nks1    = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1.ne.nks) &
     call errore('xk_et_collect','problems with nks1',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  ! copy the original points in the correct position of the list
  !
  xk_collect(:,nbase+1:nbase+nks) = xk(:,1:nks)
  !
  et_collect(:,nbase+1:nbase+nks)=et(:,1:nks)
  !
  CALL mp_sum( xk_collect, inter_pool_comm )
  !
  CALL mp_sum( et_collect, inter_pool_comm )
  !
#else
  xk_collect=xk
  et_collect=et
#endif
  !
  RETURN
  !
END SUBROUTINE xk_et_collect
!
