!
! Copyright (C) 2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE el_ph_collect( nmodes, el_ph_mat, el_ph_mat_collect, nksqtot, nksq )
  !----------------------------------------------------------------------------
  !
  ! ... This routine collects the electron-phonon matrix elements.
  !
  USE io_global, only : stdout
  USE kinds,     ONLY : DP
  USE mp_pools, ONLY : my_pool_id, npool, kunit, inter_pool_comm
  USE mp,        ONLY : mp_sum
  USE ions_base, ONLY : nat
  USE wvfct,     ONLY : nbnd
  !
  IMPLICIT NONE
  !
  INTEGER :: nksqtot, nksq, nmodes
    ! total number of k-points
    ! number of k-points per pool
  COMPLEX (DP) :: el_ph_mat(nbnd,nbnd,nksq,nmodes)
  COMPLEX (DP) :: el_ph_mat_collect(nbnd,nbnd,nksqtot,nmodes)
    ! electron-phonon matrix elements
    ! collected electron-phonon matrix elements
  !
#if defined (__MPI)
  !
  INTEGER :: nbase, rest, nks1
  !
  el_ph_mat_collect=(0.0_DP, 0.0_DP)
  !
  nks1    = ( nksqtot / npool )
  !
  rest = ( nksqtot - nks1 * npool ) 
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + 1
  !
  IF (nks1.ne.nksq) &
     call errore('el_ph_collect','problems with nks1',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nksq * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest 
  !
  ! copy the original points in the correct position of the list
  !
  el_ph_mat_collect(:,:,nbase+1:nbase+nksq,:) = el_ph_mat(:,:,1:nksq,:)
  !
  CALL mp_sum( el_ph_mat_collect, inter_pool_comm )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE el_ph_collect

!----------------------------------------------------------------------------
SUBROUTINE el_ph_distribute( nmodes, el_ph_mat, el_ph_mat_collect, nksqtot, nksq )

  !----------------------------------------------------------------------------
  !
  ! ... This routine distributes the electron-phonon matrix elements.
  !
  USE io_global, only : stdout
  USE kinds,     ONLY : DP
  USE mp_pools, ONLY : my_pool_id, npool, kunit, inter_pool_comm
  USE mp,        ONLY : mp_sum
  USE ions_base, ONLY : nat
  USE wvfct,     ONLY : nbnd
  !
  IMPLICIT NONE
  !
  INTEGER :: nksqtot, nksq, nmodes
    ! total number of k-points
    ! number of k-points per pool
    ! number of perturbation
  COMPLEX (DP) :: el_ph_mat(nbnd,nbnd,nksq,nmodes)
  COMPLEX (DP) :: el_ph_mat_collect(nbnd,nbnd,nksqtot,nmodes)
    ! electron-phonon matrix elements
    ! collected electron-phonon matrix elements
  !
#if defined (__MPI)
  !
  INTEGER :: nbase, rest, nks1
  !
  el_ph_mat=(0.0_DP, 0.0_DP)
  !
  nks1    = ( nksqtot / npool )
  !
  rest = ( nksqtot - nks1 * npool ) 
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + 1
  !
  IF (nks1.ne.nksq) &
     call errore('el_ph_distribute','problems with nks1',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nksq * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest 
  !
  ! copy the original points in the correct position of the list
  !
  el_ph_mat(:,:,1:nksq,:) = el_ph_mat_collect(:,:,nbase+1:nbase+nksq,:) 
  !
#endif
  !
  RETURN
  !
END SUBROUTINE el_ph_distribute
