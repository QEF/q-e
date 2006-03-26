!
! Copyright (C) 2002-2004 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_global
  !----------------------------------------------------------------------------
  !
#if defined (__SHMEM)
  USE shmem_include
#endif
  !
  USE parallel_include
  !
  IMPLICIT NONE 
  !
  SAVE
  !
  INTEGER :: mpime = 0  ! absolute processor index starting from 0
  INTEGER :: root  = 0  ! index of the absolute root processor
  INTEGER :: nproc = 1  ! absolute number of processor
  INTEGER :: world_comm = 0  ! communicator of all processor
  INTEGER :: kunit = 1  ! granularity of k-point distribution
  !
  ! ... indeces ( all starting from 0 !!! )
  !
  INTEGER :: me_pool     = 0  ! index of the processor within a pool 
  INTEGER :: me_image    = 0  ! index of the processor within an image
  INTEGER :: root_pool   = 0  ! index of the root processor within a pool
  INTEGER :: root_image  = 0  ! index of the root processor within an image
  INTEGER :: my_pool_id  = 0  ! index of my pool
  INTEGER :: my_image_id = 0  ! index of my image
  !
  INTEGER :: npool       = 1  ! number of "k-points"-pools
  INTEGER :: nimage      = 1  ! number of "path-images"-pools
  INTEGER :: nogrp       = 1  ! number of "task groups" 
  INTEGER :: npgrp       = 1  ! number of processor withing a "task group" 
  INTEGER :: nproc_pool  = 1  ! number of processor within a pool
  INTEGER :: nproc_image = 1  ! number of processor within an image
  !
  ! ... communicators
  !
  INTEGER :: inter_pool_comm  = 0  ! inter pool communicator
  INTEGER :: intra_pool_comm  = 0  ! intra pool communicator
  INTEGER :: inter_image_comm = 0  ! inter image communicator
  INTEGER :: intra_image_comm = 0  ! intra image communicator  
  INTEGER :: me_pgrp          = 0  ! index of the processor in plane-wave group (task grouping)
  INTEGER :: me_ogrp          = 0  ! index of the processor in orbital group (task grouping)
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE mp_global_start( root_i, mpime_i, group_i, nproc_i )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: root_i, mpime_i, group_i, nproc_i
       !
       root             = root_i
       mpime            = mpime_i
       world_comm       = group_i
       nproc            = nproc_i
       nproc_pool       = nproc_i
       nproc_image      = nproc_i
       my_pool_id       = 0
       my_image_id      = 0
       me_pool          = mpime
       me_image         = mpime
       root_pool        = root
       root_image       = root
       inter_pool_comm  = group_i
       intra_pool_comm  = group_i
       inter_image_comm = group_i
       intra_image_comm = group_i
       !
       RETURN
       !
     END SUBROUTINE mp_global_start
     !
     !-----------------------------------------------------------------------     
     SUBROUTINE mp_global_group_start( mep, myp, nprocp, num_of_pools )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !     
       INTEGER, INTENT(IN) :: mep, myp, nprocp, num_of_pools
       !
       me_pool    = mep
       my_pool_id = myp
       nproc_pool = nprocp
       npool      = num_of_pools
       !
       RETURN
       !
     END SUBROUTINE mp_global_group_start
     !
!
!----------------------------------------------------------------------------
SUBROUTINE init_pool()
  !----------------------------------------------------------------------------
  !
  ! ... This routine initialize the pool :  MPI division in pools and images
  !
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE parallel_include
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
  IF ( nproc < nimage ) &
     CALL errore( 'startup', 'nproc < nimage', 1 )
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
#endif
  !
  RETURN
  !
END SUBROUTINE init_pool
     !
END MODULE mp_global
