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
  ! ... This routine initialize the pool :  MPI division in pools and images
  !
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : mpime, me_image, my_image_id, nproc, nproc_image, &
                        nproc_pool, npool, nimage, me_pool, my_pool_id, &
                        intra_image_comm, inter_image_comm, &
                        intra_pool_comm, inter_pool_comm
  USE parallel_include
  
  USE mp_global, ONLY : mp_global_group_start
  USE para,      ONLY : me, mypool, nprocp, npool_ => npool
  !
  IMPLICIT NONE
  !
  INTEGER :: ierr = 0
  ! 
  !
#if defined (__PARA)
  !  
  ! ... here we set all parallel indeces (defined in mp_global): 
  !
  !
  ! ... number of cpus per image
  !
  nproc_image = nproc / nimage
  !
  IF ( MOD( nproc, nimage ) /= 0 ) &
     CALL errore( 'startup', 'nproc /= nproc_image * nimage', 1 )  
  !
  ! ... my_image_id  =  image index for this processor   ( 0 : nimage - 1 )
  ! ... me_image     =  processor index within the image ( 0 : nproc_image - 1 )
  !
  my_image_id = mpime / nproc_image
  me_image    = MOD( mpime, nproc_image )
  !
  CALL mp_barrier()
  !
  ! ... the intra_image_comm communicator is created
  !
  CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, &
                       my_image_id, mpime, intra_image_comm, ierr )
  !
  CALL errore( 'init_pool', 'intra_image_comm is wrong', ierr )
  !
  CALL mp_barrier()
  !
  ! ... the inter_image_comm communicator is created                     
  !     
  CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, &
                       me_image, mpime, inter_image_comm, ierr )  
  !
  CALL errore( 'init_pool', 'inter_image_comm is wrong', ierr )
  !
  ! ... number of cpus per pool of k-points (they are created inside each image)
  !
  nproc_pool = nproc_image / npool
  !
  IF ( MOD( nproc, npool ) /= 0 ) &
     CALL errore( 'startup', 'nproc /= nproc_pool * npool', 1 )  
  !
  ! ... my_pool_id  =  pool index for this processor    ( 0 : npool - 1 )
  ! ... me_pool     =  processor index within the pool  ( 0 : nproc_pool - 1 )
  !
  my_pool_id = me_image / nproc_pool    
  me_pool    = MOD( me_image, nproc_pool )
  !
  CALL mp_barrier( intra_image_comm )
  !
  ! ... the intra_pool_comm communicator is created
  !
  CALL MPI_COMM_SPLIT( intra_image_comm, &
                       my_pool_id, me_image, intra_pool_comm, ierr )
  !
  CALL errore( 'init_pool', 'intra_pool_comm is wrong', ierr )
  !
  CALL mp_barrier( intra_image_comm )
  !
  ! ... the inter_pool_comm communicator is created
  !
  CALL MPI_COMM_SPLIT( intra_image_comm, &
                       me_pool, me_image, inter_pool_comm, ierr )
  !
  call errore( 'init_pool', 'inter_pool_comm is wrong', ierr )
  !
  ! ... compatibility with old PWscf routines
  !
  ! ... me      =>  me_pool + 1
  ! ... mypool  =>  my_pool_id + 1
  ! ... nprocp  =>  nproc_pool
  !
  me     = me_pool + 1
  mypool = my_pool_id + 1
  nprocp = nproc_pool
  npool_ = npool
  !
#endif
  !
  RETURN
  !
END SUBROUTINE init_pool
