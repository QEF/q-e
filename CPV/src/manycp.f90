!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM manycp
  !----------------------------------------------------------------------------
  !
  ! ... Poor-man cp.x parallel launcher. Usage (for mpirun):
  ! ...    mpirun -np Np manycp.x -ni Ni [other options]
  ! ... or whatever is appropriate for your parallel environment
  ! ... Starts Ni cp.x instances each running on Np/Ni processors
  ! ... Each cp.x instances
  ! ... * reads input data from from cp_N.in, N=0,..,,Ni-1 if no input
  ! ...   file is specified via the -i option; from "input_file"_N
  ! ...   if command-line options -i "input_file" is specified
  ! ... * saves temporary and final data to "outdir"_N/ directory
  ! ...   (or to tmp_N/ if outdir='./')
  ! ... * writes output to cp_N.out in the current directory if no input
  ! ...   file is specified via the -i option; to "input_file"_N.out
  ! ...   if command-line options -i "input_file" is specified
  !
  USE input,             ONLY : iosys_pseudo, iosys
  USE input_parameters,  ONLY : outdir
  USE environment,       ONLY : environment_start, environment_end
  USE io_global,         ONLY : ionode, ionode_id, stdout
  USE mp_global,         ONLY : mp_startup
  USE mp_images,         ONLY : my_image_id
  USE read_input,        ONLY : read_input_file
  USE check_stop,        ONLY : check_stop_init
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  !
  INTEGER :: i
  LOGICAL :: opnd
  CHARACTER(LEN=256) :: filin, filout
  CHARACTER(LEN=7) :: image_label
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( 'MANYCP' )
  !
  ! ... Image-specific input files
  !
  image_label = '_' // int_to_char(my_image_id)
  IF ( TRIM (input_file_) == ' ') THEN
     filin = 'cp' // TRIM(image_label)  // '.in'
  ELSE
     filin = TRIM(input_file_) // TRIM(image_label) // '.in'
  END IF
  !
  ! ... Here open image-specific output files
  !
  IF ( ionode ) THEN
     !
     INQUIRE ( UNIT = stdout, OPENED = opnd )
     IF (opnd) CLOSE ( UNIT = stdout )
     IF ( TRIM (input_file_) == ' ') THEN
        filout = 'cp' // TRIM(image_label)  // '.out'
     ELSE
        filout = TRIM(input_file_) // TRIM(image_label) // '.out'
     END IF
     OPEN( UNIT = stdout, FILE = TRIM(filout), STATUS = 'UNKNOWN' )
     !
  END IF
  !
  CALL start_clock('CP')
  CALL read_input_file ( prog='CP', input_file_=filin )
  !
  ! ... Set image-specific value for "outdir", starting from input value
  ! ... (read in read_input_file)
  !
  DO i=LEN_TRIM(outdir),1,-1
     IF ( outdir(i:i) /= '/' .AND. outdir(i:i) /= '.' ) EXIT
  END DO
  ! ... i = position of last character different from '/' and '.'
  IF ( i == 0 ) THEN
     outdir = 'tmp' // trim(image_label) // '/'
  ELSE
     outdir = outdir(1:i) // trim(image_label) // '/'
  END IF
  !
  ! ... Perform actual calculation
  !
  ! ... read in pseudopotentials files and then
  ! ... copy pseudopotential parameters into internal variables
  !
  CALL iosys_pseudo()
  !
  ! ... copy input parameters from input_parameter module
  !
  CALL iosys()
  !
  CALL check_stop_init()
  !
  CALL cpr_loop( 1 )
  !
  CALL stop_run(  )
  CALL do_stop( .TRUE. )
  !
  STOP
  !
END PROGRAM manycp
