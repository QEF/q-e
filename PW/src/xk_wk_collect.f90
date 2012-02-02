!
! Copyright (C) 2007-2012 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE xk_wk_collect( xk_collect, wk_collect, xk, wk, nkstot, nks )
  !----------------------------------------------------------------------------
  !
  ! ... This routine collects the k points (with granularity kunit) among 
  ! ... nodes and sets the variable xk_collect and wk_collect with the total 
  ! ... number of k-points
  !
  USE io_global, only : stdout
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : my_pool_id, npool, kunit
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: nkstot, nks
    ! total number of k-points
    ! number of k-points per pool
  REAL (DP) :: xk(3,nks), wk(nks)
  REAL (DP) :: xk_collect(3,nkstot), wk_collect(nkstot)
    ! k-points
    ! k-point weights
  !
#if defined (__MPI)
  !
  INTEGER :: nbase, rest, nks1
  !
  xk_collect=0.d0
  !
  wk_collect=0.d0
  !
  nks1    = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1.ne.nks) &
     call errore('xk_wk_collect','problems with nks1',1)
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
  wk_collect(nbase+1:nbase+nks)=wk(1:nks)
  !
  CALL mp_sum( xk_collect, inter_pool_comm )
  !
  CALL mp_sum( wk_collect, inter_pool_comm )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE xk_wk_collect
!
!----------------------------------------------------------------------------
SUBROUTINE wg_all(wg_collect, wg, nkstot, nks )
!----------------------------------------------------------------------------
  !
  ! ... This routine collects all the weights and copy them in all pools.
  !
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : my_pool_id, npool, kunit
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_sum
  USE wvfct,     ONLY : nbnd
  !
  IMPLICIT NONE
  !
  INTEGER :: nkstot, nks
    ! total number of k-points
    ! number of k-points per pool
  REAL (DP) :: wg(nbnd, nks)
  REAL (DP) :: wg_collect(nbnd, nkstot)
    ! distributed weights of the k points of this pool
    ! collected weights of all k points
    !
#if defined (__MPI)
  !
  INTEGER :: nbase, rest, nks1
  !
  wg_collect=(0.0_DP, 0.0_DP)
  !
  nks1    = ( nkstot / npool )
  !
  rest = ( nkstot - nks1 * npool ) 
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + 1
  !
  IF (nks1.ne.nks) &
     call errore('wg_all','problems with nks1',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest 
  !
  ! copy the original wavefunctions in the correct position of the list
  !
  wg_collect(:,nbase+1:nbase+nks) = wg(:,1:nks)
  !
  CALL mp_sum( wg_collect, inter_pool_comm )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE wg_all
!
!
INTEGER FUNCTION find_current_k(ik, nkstot, nks)
  !----------------------------------------------------------------------------
  !
  ! ... This function receives the index of a k point in the list 
  ! ... of nks k-points within a pool and gives the index in the 
  ! ... full list of nkstot k-points
  !
  !
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : my_pool_id, npool, kunit
  !
  IMPLICIT NONE
  !
  INTEGER :: nkstot, nks
    ! total number of k-points
    ! number of k-points per pool
  INTEGER :: ik
    ! k-points
  !
#if defined (__MPI)
  !
  INTEGER :: nbase, rest, nks1
  !
  nks1    = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1.ne.nks) &
     call errore('isk_ngk_collect','problems with nks1',1)
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
  !
  find_current_k = nbase+ik
#else
  find_current_k = ik
#endif
  !
  RETURN

END FUNCTION find_current_k
!
