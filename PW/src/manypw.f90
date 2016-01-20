!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM manypw
  !----------------------------------------------------------------------------
  !
  ! ... Poor-man pw.x parallel launcher. Usage (for mpirun):
  ! ...    mpirun -np Np many.x -ni Ni [other options]
  ! ... or whatever is appropriate for your parallel environment
  ! ... Starts Ni pw.x instances each running on Np/Ni processors
  ! ... Each pw.x instances
  ! ... * reads input data from from pw_N.in, N=0,..,,Ni-1 if no input
  ! ...   file is specified via the -i option; from "input_file"_N
  ! ...   if command-line options -i "input_file" is specified
  ! ... * saves temporary and final data to "outdir"_N/ directory
  ! ...   (or to tmp_N/ if outdir='./')
  ! ... * writes output to pw_N.out in the current directory if no input
  ! ...   file is specified via the -i option; to "input_file"_N.out
  ! ...   if command-line options -i "input_file" is specified
  !
  USE input_parameters,  ONLY : outdir
  USE environment,       ONLY : environment_start, environment_end
  USE io_global,         ONLY : ionode, ionode_id, stdout
  USE mp_global,         ONLY : mp_startup
  USE mp_images,         ONLY : my_image_id
  USE mp,                ONLY : mp_bcast
  USE mp_world,          ONLY : root, world_comm
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  !
  INTEGER :: i, exit_status, first_image, ios
  LOGICAL :: opnd
  CHARACTER(LEN=256) :: filin, filout
  CHARACTER(LEN=7) :: image_label, var_first_index
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( 'MANYPW' )
  !
  CALL get_environment_variable( 'FIRST_IMAGE_INDEX', var_first_index )
  READ(var_first_index, *, iostat=ios) first_image
  IF(ios/= 0) first_image = 0
  CALL mp_bcast( first_image   , root, world_comm )
  !
  ! ... Image-specific input files
  !
  image_label = '_' // int_to_char(my_image_id+first_image)
  IF ( TRIM (input_file_) == ' ') THEN
     filin = 'pw' // TRIM(image_label)  // '.in'
  ELSE
     filin = TRIM(input_file_) // TRIM(image_label) // '.in'
  END IF
  !
  ! ... Here open image-specific output files
  !     Note: this operation is also done in environment_start when nimage>1
  !     which causes the creation of annoying out.#_0 empty files
  !     here we delete these files, which is not optimal but a bit more tidy
  IF ( ionode ) THEN
     !
     INQUIRE ( UNIT = stdout, OPENED = opnd )
     IF (opnd) CLOSE ( UNIT = stdout, status="DELETE" )
     IF ( TRIM (input_file_) == ' ') THEN
        filout = 'pw' // TRIM(image_label)  // '.out'
     ELSE
        filout = TRIM(input_file_) // TRIM(image_label) // '.out'
     END IF
     OPEN( UNIT = stdout, FILE = TRIM(filout), STATUS = 'UNKNOWN' )
     !
  END IF
  !
  CALL start_clock('PWSCF')
  CALL read_input_file ( prog='PW', input_file_=filin )
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
  CALL run_pwscf ( exit_status )
  !
  CALL stop_run( exit_status )
  CALL do_stop( exit_status )
  !
  STOP
  !
END PROGRAM manypw
