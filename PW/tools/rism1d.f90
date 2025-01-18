!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM rism1d
  !----------------------------------------------------------------------------
  !
  ! ... Stand-alone program of 1D-RISM calculation.
  ! ... Same input file as in pw.x or cp.x is used.
  !
  USE check_stop,           ONLY : check_stop_init, check_stop_now
  USE command_line_options, ONLY : input_file_
  USE environment,          ONLY : environment_start, environment_end
  USE input_parameters,     ONLY : trism
  USE io_global,            ONLY : stdout
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE solvmol,              ONLY : deallocate_solVs
  USE read_input,           ONLY : read_input_file
  USE rism1d_facade,        ONLY : rism1d_summary, rism1d_finalize, &
                                 & rism1d_prepare, rism1d_run, rism1d_write_to_restart, &
                                 & rism1d_write_to_show, rism1d_print_clock
  !
  IMPLICIT NONE
  !
  LOGICAL :: conv_rism1d
  !
  ! ... initialize environment
  CALL mp_startup()
  CALL environment_start('RISM1D')
  !
  ! ... read input file, and set internal variables
  CALL read_input_file('PW', input_file_)
  IF (.NOT. trism) THEN
    CALL errore('rism1d', 'One must set trism = .TRUE.', 1)
  END IF
  !
  CALL iosys()
  !
  ! ... setup initial conditions
  conv_rism1d = .FALSE.
  CALL check_stop_init()
  !
  ! ... perform calculation
  IF (check_stop_now()) THEN
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"Cannot calculate 1D-RISM")')
    WRITE(stdout, '()')
  ELSE
    CALL rism1d_summary()
    CALL rism1d_prepare()
    CALL rism1d_run(conv_rism1d)
    CALL rism1d_write_to_restart()
    IF (conv_rism1d) THEN
      CALL rism1d_write_to_show()
    END IF
  END IF
  !
  ! ... finalize environment
  CALL rism1d_finalize()
  CALL deallocate_solVs()
  WRITE(stdout, '(/,5X,"1D-RISM routines")')
  CALL rism1d_print_clock()
  CALL environment_end('RISM1D')
  CALL mp_global_end()
  !
  IF (.NOT. conv_rism1d) THEN
    STOP 1
  END IF
  !
END PROGRAM rism1d
