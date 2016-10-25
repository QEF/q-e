!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE divide_et_impera( nkstot, xk, wk, isk, nks )
  !----------------------------------------------------------------------------
  !
  ! ... This routine divides the k points across nodes, sets the variable
  ! ... nks equal to the local (on this processors) number of k-points
  ! ... (nkstot on input is the total number of k-points)
  ! ... The distributed has "granularity kunit", that is, kunit consecutive 
  ! ... points stay on the same processor. Usually kunit=1; kunit=2 is used 
  ! ... in phonon calculations, when one has interspersed k_i and k_i+q and
  ! ... it is needed that they stay on the same processor
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : my_pool_id, npool, kunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nkstot
  ! total number of k-points
  INTEGER, INTENT(INOUT) :: isk(nkstot)
  ! spin index of each kpoint (used in LSDA calculations only)
  REAL (DP), INTENT(INOUT) :: xk(3,nkstot), wk(nkstot)
  ! k-points (on all processors)
  ! k-point weights
  INTEGER, INTENT(OUT)  :: nks
  ! number of k-points per pool
  !
  INTEGER :: ik, nbase, rest
  !
  ! simple case: no pools
  !
  IF ( npool == 1 ) THEN
     nks = nkstot
     RETURN
  END IF
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) &
     CALL errore( 'divide_et_impera', 'nkstot/kunit is not an integer', nkstot )
  !
  nks    = kunit * ( nkstot / kunit / npool )
  !
  IF (nks == 0) CALL errore('divide_et_impera','some nodes have no k-points', 1)
  !
  rest = ( nkstot - nks * npool ) / kunit
  !
  IF ( my_pool_id < rest ) nks = nks + kunit
  !
  ! ... calculates nbase = the position in the list of the first point
  ! ...                    that belongs to this pool, minus one
  !
  nbase = nks * my_pool_id
  IF ( my_pool_id >= rest ) nbase = nbase + rest * kunit
  !
  ! ... displace the nks points in the pool to the first positions of the list
  !
  IF ( nbase > 0 ) THEN
     !
     xk(:,1:nks) = xk(:,nbase+1:nbase+nks)
     wk (1:nks)  = wk(nbase+1:nbase+nks)
     isk(1:nks)  =isk(nbase+1:nbase+nks)
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE divide_et_impera
!----------------------------------------------------------------------------
FUNCTION global_kpoint_index ( nkstot, ik ) RESULT (ik_g)
  !----------------------------------------------------------------------------
  
  ! ... Returns the index in the global list of k-points
  ! ... of k-point "ik" in this pool

  USE mp_pools, ONLY : npool, my_pool_id, kunit

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nkstot, ik
  ! total number of k-points
  ! index of k-point
  INTEGER  :: ik_g
  ! index in global list corresponding to ik in pool
  INTEGER  :: nks
  ! this is actually the number of k-points in this pool
  !
  INTEGER  :: nkbl, rest
  !
  ! ... nkbl = number of blocks of "kunit" k-points
  !
  nkbl = nkstot / kunit
  !
  ! ... nks = k-points per pool
  !
  nks = kunit * ( nkbl / npool )
  !
  ! ... if npool not a divisor of nkstot/kunit, find out the rest
  !
  rest = ( nkstot - nks * npool ) / kunit
  !
  ! ... Assign the remaining k-points to the first "rest" pools
  !
  IF ( my_pool_id < rest ) nks = nks + kunit
  !
  ! ... find out the global index of ik-th k-point in this pool
  !
  ik_g = nks*my_pool_id + ik
  IF ( my_pool_id >= rest ) ik_g = ik_g + rest*kunit
  !
END FUNCTION global_kpoint_index
!----------------------------------------------------------------------------
FUNCTION local_kpoint_index ( nkstot, ik_g ) RESULT (ik)
  !----------------------------------------------------------------------------
  
  ! ... Returns the local index index of a k-point, if it belongs to
  ! ... current pool, or -1 if it does not

  USE mp_pools, ONLY : npool, my_pool_id, kunit

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nkstot, ik_g
  ! total number of k-points
  ! global index of k-point
  INTEGER  :: ik ! return the index if we have it, -1 if we don't
  INTEGER  :: nks
  ! this is actually the number of k-points in this pool
  !
  INTEGER  :: nkbl, rest, nks_before, ik_g0
  !
  ! ... nkbl = number of blocks of "kunit" k-points
  nkbl = nkstot / kunit
  ! ... nks = k-points per pool
  nks = kunit * ( nkbl / npool )
  ! ... if npool not a divisor of nkstot/kunit, find out the rest
  rest = ( nkstot - nks * npool ) / kunit
  ! ... Assign the remaining k-points to the first "rest" pools
  IF ( my_pool_id < rest ) nks = nks + kunit  
  ! find the global index of the first k-point in my pool
  ik_g0 = nks*my_pool_id
  IF ( my_pool_id >= rest ) ik_g0 = ik_g0 + rest*kunit
  ik = ik_g - ik_g0
  IF(ik<=0 .or. ik>nks)  ik = -1
  !
END FUNCTION
!----------------------------------------------------------------------------
 
