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
  USE check_stop,        ONLY : check_stop_init
  USE environment,       ONLY : environment_start, environment_end
  USE input_parameters,  ONLY : outdir
  USE io_global,         ONLY : ionode, ionode_id, stdout
  USE mp_global,         ONLY : mp_startup, nimage, my_image_id
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  !
  INTEGER :: i
  LOGICAL :: opnd
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=7) :: image_label
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
#ifdef __MPI
  CALL mp_startup ( start_images=.true. )
#endif
  CALL environment_start ( 'MANYPW' )
  !
  ! ... Image-specific input files
  !
  image_label = '_' // int_to_char(my_image_id)
  IF ( TRIM (input_file_) == ' ') THEN
     filename = 'pw' // TRIM(image_label)  // '.in'
  ELSE
     filename = TRIM(input_file_) // TRIM(image_label) 
  END IF
  CALL read_input_file ( prog='PW', input_file_=filename )
  !
  ! ... Here open image-specific output files
  !
  IF ( ionode ) THEN
     !
     INQUIRE ( UNIT = stdout, OPENED = opnd )
     IF (opnd) CLOSE ( UNIT = stdout )
     IF ( TRIM (input_file_) == ' ') THEN
        filename = 'pw' // TRIM(image_label)  // '.out'
     ELSE
        filename = TRIM(input_file_) // TRIM(image_label) // '.out'
     END IF
     OPEN( UNIT = stdout, FILE = TRIM(filename), STATUS = 'UNKNOWN' )
     !
  END IF
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
  ! ... Here copy data read from input to internal modules
  !
  CALL iosys()
  !
  ! ... Perform actual calculation
  !
  CALL compute_pwscf  ( )
  !
  STOP
  !
END PROGRAM manypw
!
!----------------------------------------------------------------------------
SUBROUTINE compute_pwscf ( ) 
  !----------------------------------------------------------------------------
  !
  ! ... This is just the main of pw.x with the initialization and read section
  ! ... taken out 
  !
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma
  USE environment,      ONLY : environment_start
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_global,        ONLY : mp_global_end, intra_image_comm
  USE xml_io_base,      ONLY : create_directory, change_directory
  USE read_input,       ONLY : read_input_file
  !
  IMPLICIT NONE
  !
  !
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... open, read, close input file ALREADY DONE
  !
  !!! CALL read_input_file ('PW')
  !
  ! ... convert to internal variables ALREADY DONE
  !
  !!! CALL iosys()
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  ! call to void routine for user defined / plugin patches initializations
  !
  CALL plugin_initialization()
  !
  CALL check_stop_init()
  !
  CALL setup ()
  !
  CALL init_run()
  !
  IF ( check_stop_now() ) THEN
     CALL punch( 'all' )
     CALL stop_run( .TRUE. )
  ENDIF
  !
  main_loop: DO
     !
     ! ... electronic self-consistency
     !
     CALL electrons()
     !
     IF ( .NOT. conv_elec ) THEN
       CALL punch( 'all' )
       CALL stop_run( conv_elec )
     ENDIF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' )
     conv_ions = .TRUE.
     !
     ! ... recover from a previous run, if appropriate
     !
     IF ( restart .AND. lscf ) CALL restart_in_ions()
     !
     ! ... file in CASINO format written here if required
     !
     CALL pw2casino()
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( lstres ) CALL stress ( sigma )
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        if (fix_volume) CALL impose_deviatoric_stress(sigma)
        !
        if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions()
        !
        ! ... then we save restart information for the new configuration
        !
        IF ( istep < nstep .AND. .NOT. conv_ions ) THEN
           !
           CALL punch( 'config' )
           CALL save_in_ions()
           !
        END IF
        !
     END IF
     !
     CALL stop_clock( 'ions' )
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
     IF ( lmd .OR. lbfgs ) CALL hinit1()
     !
  END DO main_loop
  !
  ! ... save final data file
  !
  CALL punch('all')
  CALL stop_run( conv_ions )
  !
  !  END IF      
  !
  STOP
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE compute_pwscf
