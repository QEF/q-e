!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run_path( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  ! ... Called at the end of the run with flag = .TRUE. (removes 'restart')
  ! ... or during execution with flag = .FALSE. (does not remove 'restart')
  !
  USE io_global,          ONLY : ionode
  USE mp_global,          ONLY : nimage, mp_global_end
  USE environment,        ONLY : environment_end
  USE control_flags,      ONLY : lpath, twfcollect, lconstrain, &
                                 io_level, llondon
  USE io_files,           ONLY : iunwfc, iunigk, iunefield, iunefieldm,&
                                 iunefieldp, iuntmp
  USE buffers,            ONLY : close_buffer
  USE path_variables,     ONLY : path_deallocation
  USE image_io_routines,   ONLY : io_image_stop
  USE london_module,      ONLY : dealloca_london
  USE constraints_module, ONLY : deallocate_constraint
  USE input_parameters,   ONLY : deallocate_input_parameters
  USE bp,                 ONLY : lelfield
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL             :: exst, opnd, flag2
  !
  !
  !
  !
  CALL io_image_stop()
  !
! call pwscf stop run routine, close files and deallocate arrays
  !
  CALL stop_run( flag )
  !
  CALL path_deallocation()
  !
  IF ( .not. flag ) THEN
     !
     STOP 1
     !
  END IF
  !
END SUBROUTINE stop_run_path
!
!----------------------------------------------------------------------------
