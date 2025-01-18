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
  !! Author: Paolo Giannozzi (for recent additions and cleanup), others

  !! This routine divides the k points across nodes, sets the variable
  !! nks equal to the local (on this processors) number of k-points
  !! (nkstot on input is the total number of k-points)
  !! The distributed has "granularity kunit", that is, kunit consecutive 
  !! points stay on the same processor. Usually kunit=1; kunit=2 is used 
  !! in phonon calculations, when one has interspersed k_i and k_i+q and
  !! it is needed that they stay on the same processor.
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : my_pool_id, npool, kunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nkstot
  !! total number of k-points
  INTEGER, INTENT(INOUT) :: isk(nkstot)
  !! spin index of each kpoint (used in LSDA calculations only)
  REAL (DP), INTENT(INOUT) :: xk(3,nkstot)
  !! k-points (on all processors)
  REAL (DP), INTENT(INOUT) :: wk(nkstot)
  !! k-point weights
  INTEGER, INTENT(OUT)  :: nks
  !! number of k-points per pool
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
  IF (nks == 0) CALL infomsg('divide_et_impera', &
          'suboptimal parallelization: some nodes have no k-points')
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
  !! Returns the index in the global list of k-points 
  !! of k-point "ik" in this pool
  !
  USE mp_pools, ONLY : npool, my_pool_id, kunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nkstot
  !! total number of k-points
  INTEGER, INTENT(IN) :: ik
  !! index of k-point
  INTEGER  :: ik_g
  !! index in global list corresponding to ik in pool
  INTEGER  :: nks
  !! this is actually the number of k-points in this pool
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
  !! Returns the local index of a k-point, if it belongs to
  !! current pool, or -1 if it does not
  !
  USE mp_pools,   ONLY: npool, my_pool_id, kunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nkstot
  !! total number of k-points
  INTEGER, INTENT(IN) :: ik_g
  !! global index of k-point
  INTEGER  :: ik 
  !! return the index if we have it, -1 if we don't
  !
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
 
!----------------------------------------------------------------------------
SUBROUTINE pool_and_local_kpoint_index ( nkstot, ik_g, ipool, ik_l)
  !----------------------------------------------------------------------------
  !! For given global k point index ik_g, calculate the pool index where that
  !! k point belongs to, and calculate the local k point index in that pool.
  !! Adapted from kfold.f90 of EPW
  !
  USE mp_pools,   ONLY: npool, kunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nkstot
  !! total number of k-points
  INTEGER, INTENT(IN) :: ik_g
  !! global index of k-point
  INTEGER, INTENT(OUT) :: ipool
  !! index of the pool where ik_g belongs to
  INTEGER, INTENT(OUT) :: ik_l
  !! local index of k-point
  !
  INTEGER  :: jpool, nks, rest, nbase
  !
  IF (ik_g > nkstot) CALL errore("pool_and_local_kpoint_index", &
      "ik_g cannot be greater than nkstot", 1)
  !
  IF (npool == 1) THEN
    ! simple case: no pools
    ipool = 0
    ik_l = ik_g
    RETURN
  ENDIF
  !
  ! Loop over pools and find whether ik_g belongs to this pool
  !
  DO jpool = 0, npool - 1
    !
    nks = kunit * ( nkstot / kunit / npool )
    !
    rest = ( nkstot - nks * npool ) / kunit
    !
    IF ( jpool < rest ) nks = nks + kunit
    !
    ! ... calculates nbase = the position in the list of the first point
    ! ...                    that belongs to this pool, minus one
    !
    nbase = nks * jpool
    IF ( jpool >= rest ) nbase = nbase + rest * kunit
    !
    ! This pool has global k point index nbase+1:nbase+nks
    !
    IF (ik_g >= nbase+1 .AND. ik_g <= nbase+nks) THEN
      ipool = jpool
      ik_l = ik_g - nbase
      RETURN
    ENDIF
    !
  ENDDO
  !
  ! if the subroutine did not return inside the loop, something is wrong.
  !
  CALL errore("pool_and_local_kpoint_index", "ipool not found", 1)
  !
END SUBROUTINE pool_and_local_kpoint_index
!----------------------------------------------------------------------------
