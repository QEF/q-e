!
! Copyright (C) 2002-2003 PWSCF-FPMD-CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

     MODULE mp_global

       USE shmem_include
       USE parallel_include
       USE kinds
       IMPLICIT NONE
       SAVE

         INTEGER :: mpime = 0  ! absolute processor index starting from 0
         INTEGER :: root  = 0  ! index of the absolute root processor
         INTEGER :: nproc = 1  ! absolute number of processor
         INTEGER :: group = 0  ! group communicator
         INTEGER :: me_pool    = 0  ! index of the processor within a pool (starting from 0 !!! )
         INTEGER :: root_pool  = 0  ! index of the root processor within a pool (starting from 0 !!! )
         INTEGER :: my_pool_id = 0  ! index of my pool (starting from 0 !!! )
         INTEGER :: nproc_pool = 1       ! number of processor within a pool
         INTEGER :: inter_pool_comm = 0  ! inter pool communicator
         INTEGER :: intra_pool_comm = 0  ! intra pool communicator

       CONTAINS

        SUBROUTINE mp_global_start(root_i, mpime_i, group_i, nproc_i )
          INTEGER, INTENT(IN) :: root_i, mpime_i, group_i, nproc_i
          root   = root_i
          mpime  = mpime_i
          group  = group_i
          nproc  = nproc_i
          nproc_pool = nproc_i
          my_pool_id = 0
          me_pool    = mpime
          root_pool  = root
          inter_pool_comm = group_i
          intra_pool_comm = group_i
          RETURN
        END SUBROUTINE

        SUBROUTINE mp_global_group_start( mep, myp, intrap, interp, nprocp )
          INTEGER, INTENT(IN) :: mep, myp, intrap, interp, nprocp
          me_pool    = mep
          my_pool_id = myp
          intra_pool_comm = intrap
          inter_pool_comm = interp
          nproc_pool = nprocp
          RETURN
        END SUBROUTINE

     END MODULE mp_global
