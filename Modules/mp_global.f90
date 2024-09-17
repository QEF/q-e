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
  !! Wrapper module, for compatibility. Contains a few "leftover" variables
  !! used for checks (all the *_file variables, read from data file),
  !! plus the routine mp_startup initializing MPI and the command line, 
  !! plus the routine mp_global_end stopping MPI.  
  !! Do not use this module to reference variables (e.g. communicators)
  !! belonging to each of the various parallelization levels:
  !! use the specific modules instead. 
  ! ... PLEASE DO NOT ADD NEW STUFF TO THIS MODULE. Removing stuff is ok.
  !
  USE mp_world, ONLY: world_comm, mp_world_start, mp_world_end
  USE mp_images
  USE mp_pools
  USE mp_bands
  USE mp_bands_TDDFPT
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... number of processors for the various groups: values read from file
  !
  INTEGER :: nproc_file = 1
  INTEGER :: nproc_image_file = 1
  INTEGER :: nproc_pool_file  = 1
  INTEGER :: nproc_ortho_file = 1
  INTEGER :: nproc_bgrp_file  = 1
  INTEGER :: ntask_groups_file= 1
  INTEGER :: nyfft_file= 1
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_startup ( my_world_comm, start_images, images_only )
    !-----------------------------------------------------------------------
    !! This wrapper subroutine initializes all parallelization levels.  
    !! If option with_images=.true., processes are organized into images,
    !! each performing a quasi-indipendent calculation, such as a point
    !! in configuration space (NEB) or a phonon irrep (PHonon).  
    !! Within each image processes are further subdivided into various
    !! groups and parallelization levels.
    !
    !  IMPORTANT NOTICE 1: since the command line is read here, it may be
    !                      convenient to call it in serial execution as well
    !  IMPORTANT NOTICE 2: most parallelization levels are initialized here 
    !                      but at least some will be moved to a later stage
    !  SPECIAL CASE: command-line options "-ntg" and "-nyfft", introduced to
    !                improve scaling, coexist and are in part interchangeable.
    !                If task groups are available, -ntg is faster than -nyfft
    !                because it communicates larger data chuncks less frequently
    !                Sometimes task groups are not available as for instance
    !                when metagga is used or realus or for conjugate gradient.
    !                For those cases, -nyfft can be used instead.
    !                You may specify one or another: the same value will be set
    !                for both ntg and nyfft. These variables are kept separated
    !                to help understanding which operation belong to task groups
    !                or to nyfft, allowing to differenciate them if need arises.
    !
    USE command_line_options, ONLY : get_command_line, &
        nimage_, npool_, nband_, ntg_, nyfft_
    USE parallel_include
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN), OPTIONAL :: my_world_comm
    LOGICAL, INTENT(IN), OPTIONAL :: start_images
    LOGICAL, INTENT(IN), OPTIONAL :: images_only
    LOGICAL :: do_images, only_images
    INTEGER :: my_comm
    !
    my_comm = MPI_COMM_WORLD
    IF ( PRESENT(my_world_comm) ) my_comm = my_world_comm
    !
    CALL mp_world_start( my_comm )
    CALL get_command_line ( )
#if defined (__CUDA)
    IF ( ntg_ > 1 ) CALL errore('mp_startup','No task groups for GPUs',ntg_)
#endif
    !
    do_images = .FALSE.
    IF ( PRESENT(start_images) ) do_images = start_images
    IF ( do_images ) THEN
       CALL mp_start_images ( nimage_, world_comm )
    ELSE
       CALL mp_init_image ( world_comm  )
    END IF
    !
    ! if images_only is present and .t., do not initialize all other 
    ! parallelization levels - Temporary solution for backward-compatibility
    !
    only_images = .FALSE.
    IF ( PRESENT(images_only) ) only_images = images_only
    IF ( only_images ) RETURN
    !
    ! npool_ is 0 if not specified in command line
    IF ( npool_== 0 ) npool_ = 1
    CALL mp_start_pools ( npool_, intra_image_comm )
    ! ntg_ is 0 if not specified in command line
    IF ( ntg_== 0 ) ntg_ = 1
#if defined (__CUDA_OPTIMIZED)
    CALL mp_start_bands ( 1 , ntg_, nyfft_, intra_pool_comm )
#else
    CALL mp_start_bands ( nband_, ntg_, nyfft_, intra_pool_comm )
#endif
    !
    RETURN
    !
  END SUBROUTINE mp_startup
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_end ( )
    !-----------------------------------------------------------------------
    !
    USE mp, ONLY : mp_comm_free
    USE mp_orthopools, ONLY: mp_stop_orthopools
    !
    CALL mp_comm_free ( intra_bgrp_comm )
    CALL mp_comm_free ( inter_bgrp_comm )
    CALL mp_comm_free ( intra_pool_comm )
    CALL mp_comm_free ( inter_pool_comm )
    CALL mp_stop_orthopools( ) ! cleans orthopools if used in exx
    CALL mp_world_end( )
    !
    RETURN
    !
  END SUBROUTINE mp_global_end
  !
END MODULE mp_global
