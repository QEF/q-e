!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_world
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_start, mp_end
  USE io_global, ONLY : meta_ionode_id, meta_ionode
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... World group (all processors)
  !
  INTEGER :: nproc = 1  ! number of processors
  INTEGER :: mpime = 0  ! processor index (starts from 0 to nproc-1)
  INTEGER :: root  = 0  ! index of the root processor
  INTEGER :: world_comm = 0  ! communicator
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_start ( my_world_comm )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: my_world_comm
    !
    world_comm = my_world_comm
    ! ... get the basic parameters from communications sub-system
    ! ... to handle processors
    ! ... nproc = number of processors
    ! ... mpime = processor number, starting from 0
    !
    CALL mp_start( nproc, mpime, world_comm )
    !
    ! ... meta_ionode is true if this processor is the root processor
    ! ... of the world group - "ionode_world" would be a better name
    ! ... meta_ionode_id is the index of such processor
    !
    meta_ionode = ( mpime == root )
    meta_ionode_id = root
    !
    RETURN
    !
  END SUBROUTINE mp_global_start
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_global_end ( )
    !-----------------------------------------------------------------------
    !
    CALL mp_barrier( world_comm )
    CALL mp_end ()
    !
  END SUBROUTINE mp_global_end
  !
END MODULE mp_world
