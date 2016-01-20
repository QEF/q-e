!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE stop_d3 (flag)
!-----------------------------------------------------------------------
!
!    This routine closes all files before stopping
!    flag is no longer used
!
  USE pwcom
  USE phcom
  USE d3com
  USE control_flags, ONLY : twfcollect
  USE io_files,      ONLY : iunigk
  USE mp_global,     ONLY : me_pool, root_pool, mp_global_end

  use control_lr, ONLY : lgamma

  IMPLICIT NONE

  LOGICAL :: flag

  IF (twfcollect ) THEN
     CLOSE (unit = iuwfc, status = 'delete')
  ELSE
     CLOSE (unit = iuwfc, status = 'keep')
  END IF
  CLOSE (unit = iubar, status = 'keep')
  CLOSE (unit = iudwf, status = 'keep')

  IF ( me_pool == root_pool ) THEN
     !
     CLOSE (unit = iudrho, status = 'keep')
     IF (.NOT.lgamma) CLOSE (unit = iud0rho, status = 'keep')
     IF(nlcc_any) THEN
       CLOSE (unit = iudrho+1000, status = 'keep')
       IF (.NOT.lgamma) CLOSE (unit = iud0rho+1000, status = 'keep')
     ENDIF
     !
  END IF

  CLOSE (unit = iunigk, status = 'delete')
  IF (.NOT.lgamma) THEN
     CLOSE (unit = iud0qwf, status = 'keep')
     CLOSE (unit = iudqwf, status = 'keep')
  ENDIF
  CLOSE (unit = iupdqvp, status = 'keep')
  IF (.NOT.lgamma) CLOSE (unit = iupd0vp, status = 'keep')
  IF (degauss.NE.0.d0) THEN
     CLOSE (unit = iudpdvp_1, status = 'keep')
     IF (.NOT.lgamma) THEN
        CLOSE (unit = iudpdvp_2, status = 'keep')
        CLOSE (unit = iudpdvp_3, status = 'keep')
     ENDIF
  ENDIF
  CALL print_clock_d3

  CALL mp_global_end ()

  STOP
  RETURN
END SUBROUTINE stop_d3
