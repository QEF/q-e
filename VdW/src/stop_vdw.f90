!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE stop_vdw
  !--------------------------------------------------------------------
  !
  ! Synchronize processes before stopping.
  !
  USE control_flags, ONLY: twfcollect
  USE io_files, ONLY: iunwfc
  USE mp_global, ONLY: mp_global_end
  USE parallel_include
#ifdef __MPI

  INTEGER :: info
  LOGICAL :: op

  CALL print_clock_vdw ()

  INQUIRE ( iunwfc, opened = op )

  IF ( op ) THEN
     IF (twfcollect) THEN
        CLOSE (unit = iunwfc, status = 'delete')
     ELSE
        CLOSE (unit = iunwfc, status = 'keep')
     ENDIF
  ENDIF

#endif

  CALL mp_global_end ()

  STOP
END SUBROUTINE stop_vdw
