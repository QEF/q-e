!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_global
  !----------------------------------------------------------------------------
  !
  ! ... Wrapper module, for compatibility only, DO NOT USE ANY LONGER
  !
  USE mp_world
  USE mp_images
  USE mp_pools
  USE mp_pots
  USE mp_bands
  USE mp_diag
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... number of processors for the various groups: values read from file
  !
  INTEGER :: nproc_file = 1
  INTEGER :: nproc_image_file = 1
  INTEGER :: nproc_pool_file  = 1
  INTEGER :: nproc_pot_file = 1
  INTEGER :: nproc_ortho_file = 1
  INTEGER :: nproc_bgrp_file  = 1
  INTEGER :: ntask_groups_file= 1
  !
  ! ... "task" groups (for band parallelization of FFT)
  !
  INTEGER :: ntask_groups = 1  ! number of proc. in an orbital "task group" 
  PRIVATE :: ntask_groups
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_startup ( start_images )
    !-----------------------------------------------------------------------
    ! ... This wrapper subroutine initializes all parallelization levels.
    ! ... If option with_images=.true., processes are organized into images,
    ! ... each performing a quasi-indipendent calculation, such as a point
    ! ..  in configuration space (NEB) or a phonon irrep (PHonon)
    ! ... Within each image processes are further subdivided into various
    ! ... groups and parallelization levels
    !
    USE command_line_options, ONLY : get_command_line, &
        nimage_, npool_, npot_, ndiag_, nband_, ntg_
    USE io_global, ONLY : io_global_start
    IMPLICIT NONE
    LOGICAL, INTENT(IN), OPTIONAL :: start_images
    LOGICAL :: do_images
    !
    CALL mp_global_start( )
    CALL get_command_line ( )
    !
    do_images = PRESENT(start_images) 
    IF ( do_images ) do_images = start_images
    IF ( do_images ) THEN
       CALL mp_start_images ( nimage_, world_comm )
    ELSE
       CALL mp_init_image ( world_comm  )
    END IF
    !
    CALL mp_start_pots  ( npot_, intra_image_comm )
    CALL mp_start_pools ( npool_, intra_image_comm )
    CALL mp_start_bands ( nband_, intra_pool_comm )
    ntask_groups = ntg_
    CALL mp_start_diag  ( ndiag_, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE mp_startup
  !
  FUNCTION get_ntask_groups()
     IMPLICIT NONE
     INTEGER :: get_ntask_groups
     get_ntask_groups = ntask_groups
     RETURN
  END FUNCTION get_ntask_groups
  !
END MODULE mp_global
