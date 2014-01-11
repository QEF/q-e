!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_pools
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Pool groups (processors within a pool of k-points)
  ! ... Subdivision of image group, used for k-point parallelization 
  !
  INTEGER :: npool       = 1  ! number of "k-points"-pools
  INTEGER :: nproc_pool  = 1  ! number of processors within a pool
  INTEGER :: me_pool     = 0  ! index of the processor within a pool 
  INTEGER :: root_pool   = 0  ! index of the root processor within a pool
  INTEGER :: my_pool_id  = 0  ! index of my pool
  INTEGER :: inter_pool_comm  = 0  ! inter pool communicator
  INTEGER :: intra_pool_comm  = 0  ! intra pool communicator
  !
  INTEGER :: kunit = 1  ! granularity of k-point distribution 
                        ! kunit=1 standard case. In phonon k and k+q must
                        ! be on the same pool, so kunit=2.
  ! 
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_pools( npool_, parent_comm )
    !---------------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into "pools"
    ! ... Requires: npool_, read from command line
    ! ...           parent_comm, typically world_comm = group of all processors
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npool_, parent_comm
    !
    INTEGER :: parent_nproc = 1, parent_mype  = 0
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... npool_ must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    npool = npool_
    IF ( npool < 1 .OR. npool > parent_nproc ) CALL errore( 'mp_start_pools',&
                          'invalid number of pools, out of range', 1 )

    IF ( MOD( parent_nproc, npool ) /= 0 ) CALL errore( 'mp_start_pools', &
           'invalid number of pools, parent_nproc /= nproc_pool * npool', 1 )  
    !
    ! ... number of cpus per pool of k-points (created inside each parent group)
    !
    nproc_pool = parent_nproc / npool
    !
    !
    ! ... my_pool_id  =  pool index for this processor    ( 0 : npool - 1 )
    ! ... me_pool     =  processor index within the pool  ( 0 : nproc_pool - 1 )
    !
    my_pool_id = parent_mype / nproc_pool    
    me_pool    = MOD( parent_mype, nproc_pool )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the intra_pool_comm communicator is created
    !
    CALL mp_comm_split ( parent_comm, my_pool_id, parent_mype, intra_pool_comm )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the inter_pool_comm communicator is created
    !
    CALL mp_comm_split ( parent_comm, me_pool, parent_mype, inter_pool_comm )
    !
#endif
    !
    RETURN
  END SUBROUTINE mp_start_pools
  !
END MODULE mp_pools
