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
  INTEGER :: group = 0  ! group communicator
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
       group            = group_i
       nproc            = nproc_i
       nproc_pool       = nproc_i
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
END MODULE mp_global
