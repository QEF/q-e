!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_pool()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initialize the pool
  !
#if defined (__PARA)
  !
  USE para,      ONLY : me, mypool, npool, nprocp, MPI_COMM_POOL, MPI_COMM_ROW
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : mp_global_group_start
  USE parallel_include   
  !
  IMPLICIT NONE
  !
  ! ... MPI pool division in pools
  !
  INTEGER :: ierr, rank
  !
  !
  ! ... set "mypool" and reset "me"
  !
  rank   = me
  mypool = ( me - 1 ) / nprocp + 1
  me     = me - ( mypool - 1 ) * nprocp
  !
  CALL mp_barrier()
  !
  CALL MPI_comm_split( MPI_COMM_WORLD, mypool, rank, MPI_COMM_POOL, ierr )
  !
  CALL errore( 'init_pool', 'MPI_COMM_POOL is wrong', ierr )
  !
  IF ( npool > 1 ) THEN
     !
     CALL MP_barrier()
     !
     CALL MPI_comm_split( MPI_COMM_WORLD, me, rank, MPI_COMM_ROW, ierr )
     !
     call errore( 'init_pool', 'MPI_COMM_ROW is wrong', ierr )
     !
  END IF                                 
  !
  ! ... Initialize globally accessible pool variables
  !
  ! ... me             =>  me_pool + 1
  ! ... mypool         =>  my_pool_id + 1
  ! ... MPI_COMM_POOL  =>  intra_pool_comm
  ! ... MPI_COMM_ROW   =>  inter_pool_comm
  ! ... nprocp         =>  nproc_pool
  ! 
  CALL mp_global_group_start( ( me - 1 ), ( mypool - 1 ), &
                              MPI_COMM_POOL, MPI_COMM_ROW, nprocp )
  !                            
#endif
  !
  RETURN
  !
END SUBROUTINE init_pool
