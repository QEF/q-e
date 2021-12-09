!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... parallel subroutines (wrappers to MPI calls) used by PWscf
! ... for k-point parallelization ("pools")
!
!----------------------------------------------------------------------------
SUBROUTINE poolscatter_matrix( length, nkstot, f_in, nks, f_out )
  !----------------------------------------------------------------------------
  !! This routine distributes a real array (e.g. eigenvalues) from the
  !! first processor of the first pool to all other pools.
  !
  !! * On input: f_in(length,nkstot) contains data for all "nkstot" k-points
  !! * On output: f_out(length,nks) contains the data for the "nks" k-point
  !!   belonging to the current pool. f_in and f_out may coincide.
  !
  ! FIXME: The copy from f_in to f_out should be made safer if f_in=f_out
  ! FIXME: Quick-and-dirty implementation: shouldn't broadcast the contents of
  ! FIXME: the first processor, just distribute the content of each processor
  ! FIXME: of the first pool to each corresponding processors of other pools
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : intra_pool_comm, inter_pool_comm, &
                        my_pool_id, npool, me_pool, root_pool, kunit
  USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! first dimension of vectors f_in and f_out
  INTEGER, INTENT(IN) :: nkstot
  !! number of k-points per pool
  INTEGER, INTENT(IN) :: nks
  !! total number of k-points
  COMPLEX(DP), INTENT(IN) :: f_in(length,length,nkstot)
  !! input: contains values for all k-point
  COMPLEX(DP), INTENT(OUT) :: f_out(length,length,nks)
  !! output: only for k-points of mypool
  !
  ! ... local variables
  !
  INTEGER :: rest, nbase
  ! the rest of the integer division nkstot / npool
  ! the position in the original list
  !
  ! ... copy from the first node of the first pool
  ! ... to the first node of all the other pools
  !
#if defined (__MPI)
  IF ( me_pool == root_pool ) &
     CALL mp_bcast( f_in, root_pool, inter_pool_comm )
#endif
  !
  ! ... distribute the vector on the first node of each pool
  !
  rest = nkstot / kunit - ( nkstot / kunit / npool ) * npool 
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  f_out(:,:,1:nks) = f_in(:,:,(nbase+1):(nbase+nks))
  !
  ! ... copy from the first proc of a pool to the other procs of the same pool
  !
#if defined (__MPI)
  CALL mp_bcast( f_out, root_pool, intra_pool_comm )
#endif
  !
  RETURN
  !
END SUBROUTINE poolscatter_matrix
!
!----------------------------------------------------------------------------
SUBROUTINE poolscatter( length, nkstot, f_in, nks, f_out )
  !----------------------------------------------------------------------------
  !! This routine distributes a real array (e.g. eigenvalues) from the
  !! first processor of the first pool to all other pools.
  !
  !! * On input: f_in(length,nkstot) contains data for all "nkstot" k-points
  !! * On output: f_out(length,nks) contains the data for the "nks" k-point
  !!   belonging to the current pool. f_in and f_out may coincide.
  !
  ! FIXME: The copy from f_in to f_out should be made safer if f_in=f_out
  ! FIXME: Quick-and-dirty implementation: shouldn't broadcast the contents of
  ! FIXME: the first processor, just distribute the content of each processor
  ! FIXME: of the first pool to each corresponding processors of other pools
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : intra_pool_comm, inter_pool_comm, &
                        my_pool_id, npool, me_pool, root_pool, kunit
  USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! first dimension of vectors f_in and f_out
  INTEGER, INTENT(IN) :: nkstot
  !! number of k-points per pool
  INTEGER, INTENT(IN) :: nks
  !! total number of k-points
  REAL(DP), INTENT(IN) :: f_in(length,nkstot)
  !! input: contains values for all k-point
  REAL(DP), INTENT(OUT) :: f_out(length,nks)
  !! output: only for k-points of mypool
  !
  ! ... local variables
  !
  INTEGER :: rest, nbase
  ! the rest of the integer division nkstot / npool
  ! the position in the original list
  !
  ! ... copy from the first node of the first pool
  ! ... to the first node of all the other pools
  !
#if defined (__MPI)
  IF ( me_pool == root_pool ) &
     CALL mp_bcast( f_in, root_pool, inter_pool_comm )
#endif
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
  ! ... copy from the first proc of a pool to the other procs of the same pool
  !
#if defined (__MPI)
  CALL mp_bcast( f_out, root_pool, intra_pool_comm )
#endif
  !
  RETURN
  !
END SUBROUTINE poolscatter
!
!----------------------------------------------------------------------------
SUBROUTINE poolcollect( length, nks, f_in, nkstot, f_out )
  !----------------------------------------------------------------------------
  !! Collects a real array f_in, distributed across pools, from all pools,
  !! into a real array f_out.
  !
  !! * On input: f_in(length,nks) contains data for the "nks" k-points
  !!   of the current pool, on all pools;
  !! * On output: f_out(length,nkstot) contains data for all "nkstot" k-points
  !!   on all pools.
  !
  !! f_in and f_out must differ! Honors "kunit"
  !
  USE kinds,     ONLY : DP
  USE mp_pools,  ONLY : my_pool_id, npool, kunit, &
                        inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! first dimension of arrays
  INTEGER, INTENT(IN) :: nks
  !! number of k-points per pool
  INTEGER, INTENT(IN) :: nkstot
  !! total number of k-points
  REAL(DP), INTENT(IN) :: f_in(length,nks)
  !! pool-distributed function
  REAL(DP), INTENT(OUT) :: f_out(length,nkstot)
  !! pool-collected function
  !
  ! ... local variables
  !
  INTEGER :: nbase, rest, nks1
  !
  nks1 = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1 /= nks) CALL errore( 'xk_collect', 'inconsistent number of k-points', 1 )
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
  f_out = 0.0_DP
  f_out(:,nbase+1:nbase+nks) = f_in(:,1:nks)
  !
  CALL mp_sum( f_out, inter_pool_comm )
  !
  !
  RETURN
  !
END SUBROUTINE poolcollect
!
!-----------------------------------------------------------------------
SUBROUTINE poolrecover( vec, length, nkstot, nks )
  !----------------------------------------------------------------------- 
  !! Gathers a real array (e.g. eigenvalues) distributed across pools
  !! from all pools into the first pool.
  !
  !! * differences from "poolcollect": 
  !!    * in-place
  !!    * result available only on first proc of first pool;
  !! * On input: vec(length,nks) contains data for the "nks" k-points
  !!   of the current pool;
  !! * On output: vec(length,nkstot) contains data for all "nkstot" k-points
  !!   on the first processor of the first pool;
  !! * vec(1:length,1:nks) is unchanged on output;
  !! * Opposite of "poolscatter". Honors "kunit".
  !
  USE kinds,            ONLY : DP
  USE mp_images,        ONLY : intra_image_comm
  USE mp_pools,         ONLY : inter_pool_comm, npool, me_pool, root_pool, &
                               my_pool_id, kunit
  USE mp,               ONLY : mp_barrier  
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER  :: length
  !! first dimension of arrays
  INTEGER  :: nks
  !! number of k-points per pool
  INTEGER  :: nkstot
  !! total number of k-points
  REAL(DP) :: vec(length,nkstot)
  !! I/O: see routine comments
  !
#if defined (__MPI)  
  !
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: i, nks1, rest, fine, nbase, info
  !
  !
  IF ( npool <= 1 ) RETURN
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) &
       CALL errore( 'poolrecover', 'nkstot/kunit is not an integer', nkstot )
  !
  nks1 = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( me_pool == root_pool .AND. my_pool_id > 0 .AND. nks > 0 ) THEN
     !
     CALL MPI_SEND( vec, (length*nks), MPI_DOUBLE_PRECISION, 0, 17, &
                    inter_pool_comm, info )
     !     
     CALL errore( 'poolrecover', 'info<>0 in send', info )
     !
  ENDIF
  !
  DO i = 2, npool
     !
     IF ( i <= rest ) THEN
        !
        fine = nks1 + kunit
        !
        nbase = ( nks1 + kunit ) * ( i - 1 )
        !
     ELSE
        !
        fine = nks1
        !
        nbase = rest * (nks1 + kunit) + (i - 1 - rest) * nks1
        !
     ENDIF
     !
     IF ( me_pool == root_pool .AND. my_pool_id == 0 .AND. nbase < nkstot ) THEN
        !
        CALL MPI_RECV( vec(1,nbase+1), (length*fine), MPI_DOUBLE_PRECISION, &
                       (i-1), 17, inter_pool_comm, status, info )
        !
        CALL errore( 'poolrecover', 'info<>0 in recv', info )
        !
     ENDIF
     !
  ENDDO
  !
#endif
  !
  RETURN
  !
END SUBROUTINE poolrecover
!
!------------------------------------------------------------------------
SUBROUTINE ipoolrecover( ivec, length, nkstot, nks )
  !----------------------------------------------------------------------
  !! As poolrecover, for an integer vector.
  !
  USE mp_images,        ONLY : intra_image_comm
  USE mp_pools,         ONLY : inter_pool_comm, npool, me_pool, &
                               root_pool, my_pool_id, kunit
  USE mp,               ONLY : mp_barrier  
  USE parallel_include    
  !
  IMPLICIT NONE
  !
  INTEGER  :: length
  !! first dimension of arrays
  INTEGER  :: nks
  !! number of k-points per pool
  INTEGER  :: nkstot
  !! total number of k-points
  INTEGER :: ivec(length,nkstot)
  !! I/O: see comments in \(\texttt{poolrecover}\) routine.
  !
#if defined (__MPI)  
  !
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: i, nks1, rest, fine, nbase, info
  !
  !
  IF ( npool <= 1 ) RETURN
  !
  IF ( MOD( nkstot, kunit ) /= 0 ) & 
       CALL errore( 'poolrecover', 'nkstot/kunit is not an integer', nkstot )
  !
  nks1 = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  CALL mp_barrier( intra_image_comm )
  !
  IF ( me_pool == root_pool .AND. my_pool_id > 0 .AND. nks > 0 ) THEN
     !
     CALL MPI_SEND( ivec, (length*nks), MPI_INTEGER, 0, 17, &
                    inter_pool_comm, info )
     !
     CALL errore( 'ipoolrecover', 'info<>0 in send', info )
     !
  ENDIF
  !
  DO i = 2, npool
     !
     IF ( i <= rest ) THEN
        !
        fine = nks1 + kunit
        !
        nbase = ( nks1 + kunit ) * ( i - 1 )
        !
     ELSE
        !
        fine = nks1
        !
        nbase = rest * ( nks1 + kunit ) + ( i - 1 - rest ) * nks1
        !
     ENDIF
     !
     IF ( me_pool == root_pool .AND. my_pool_id == 0 .AND. nbase < nkstot ) THEN
        !
        CALL MPI_RECV( ivec(1,nbase+1), (length*fine), MPI_INTEGER, &
                       (i-1), 17, inter_pool_comm, status, info )
        !
        CALL errore( 'ipoolrecover', 'info<>0 in recv', info )
        !
     ENDIF
     !
  ENDDO
  !
#endif
  !
  RETURN
  !
END SUBROUTINE ipoolrecover
