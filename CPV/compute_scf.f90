!
! Copyright (C) 2002-2004 PWSCF-FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE compute_scf( N_in, N_fin, stat  )
  !-----------------------------------------------------------------------
  !
  ! ... this subroutine is the main scf-driver for all "path" calculations
  ! ... ( called by Modules/path_base.f90/born_oppenheimer() subroutine )
  !
  USE kinds,             ONLY : DP
  USE input_parameters,  ONLY : if_pos, sp_pos, rd_pos, ion_positions
  USE input_parameters,  ONLY : restart_mode
  USE ions_base,         ONLY : nat
  USE constants,         ONLY : e2
  USE control_flags,     ONLY : conv_elec, ethr, ndr, program_name
  USE io_files,          ONLY : iunpath, iunexit, outdir, prefix, scradir
  USE io_global,         ONLY : stdout
  USE path_formats,      ONLY : scf_fmt
  USE path_variables,    ONLY : pos, pes, grad_pes, num_of_images, &
                                dim, suspended_image, frozen
  USE parser,            ONLY : int_to_char
  USE mp_global,         ONLY : mpime, my_pool_id
  USE mp,                ONLY : mp_barrier
  USE check_stop,        ONLY : check_stop_now
  USE restart_file,      ONLY : check_restartfile
  USE main_module,       ONLY : cpmain
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN)   :: N_in, N_fin
  LOGICAL, INTENT(OUT)  :: stat
  !
  ! ... local variables definition
  ! 
  INTEGER               :: image
  REAL (KIND=DP)        :: tcpu 
  CHARACTER (LEN=256)   :: outdir_saved
  CHARACTER (LEN=80)    :: restart_mode_saved
  LOGICAL               :: file_exists, opnd, tstop 
  REAL(DP), ALLOCATABLE :: tau(:,:)
  REAL(DP), ALLOCATABLE :: fion(:,:)
  REAL(DP)              :: etot
  INTEGER               :: ia, is, isa, ipos
  !
  ! ... end of local variables definition
  !
  !
  stat = .TRUE.
  tcpu = 0.0d0
  !
  IF( nat > 0 ) THEN
     ALLOCATE( tau( 3, nat ), fion( 3, nat ) )
  ELSE
     CALL errore( ' cploop ', ' nat less or equal 0 ', 1 )
  END IF
  !
  outdir_saved = outdir
  restart_mode_saved = restart_mode
  ! 
  DO image = N_in, N_fin
     !
     IF ( frozen(image) ) CYCLE
     !
     suspended_image = image
     !
     tstop = check_stop_now()
     stat  = .NOT. tstop
     !
     IF( tstop ) RETURN
     !
     outdir  = TRIM( outdir_saved ) // "/" // TRIM( prefix ) // "_" // &
               TRIM( int_to_char( image ) ) // "/" 
     !          
     scradir = outdir
     !
     WRITE( UNIT = iunpath, FMT = scf_fmt ) tcpu, image
     !
     ! ... unit stdout is connected to the appropriate file
     !
     IF ( mpime == 0 .AND. my_pool_id == 0 ) THEN
        INQUIRE( UNIT = stdout, OPENED = opnd )
        IF ( opnd ) CLOSE( UNIT = stdout )
        OPEN( UNIT = stdout, FILE = TRIM( outdir )//'CP.out', &
             STATUS = 'UNKNOWN', POSITION = 'APPEND' )
     END IF
     !
     DO ia = 1, nat
        !
        ! ... rd_pos already in bohr units
        !
        rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),image) 
        !
     END DO
     !
     ion_positions = 'from_input'
     !
     IF( check_restartfile( outdir, ndr ) ) THEN
        WRITE( iunpath, * ) ' restarting calling readfile '
        restart_mode = 'restart'
     ELSE
        WRITE( iunpath, * ) ' restarting from scratch '
        restart_mode = 'from_scratch'
     END IF
     !
     ! ... perform an electronic minimization using CPMAIN
     !
     IF( program_name == 'CP90' ) THEN
       CALL cprmain( tau, fion, etot )
     ELSE IF( program_name == 'FPMD' ) THEN
       CALL cpmain( tau, fion, etot )
     ELSE
       CALL errore( ' compute_scf ', ' unknown program ', 1 )
     END IF
     !
     IF ( mpime == 0 .AND. my_pool_id == 0 ) THEN
        INQUIRE( UNIT = stdout, OPENED = opnd )
        IF ( opnd ) CLOSE( UNIT = stdout )
     END IF
     !
     IF ( .NOT. conv_elec ) THEN
        !
        WRITE( iunpath, '(/,5X,"WARNING :  scf convergence NOT achieved",/)' )
        !   
        stat = .FALSE.
        !
        RETURN
        !
     END IF
     !
     ! ... gradients already in ( hartree / bohr )
     !
     DO ia = 1, nat
        grad_pes( 1 + ( ia - 1 ) * 3, image ) = - fion( 1, ia )
        grad_pes( 2 + ( ia - 1 ) * 3, image ) = - fion( 2, ia )
        grad_pes( 3 + ( ia - 1 ) * 3, image ) = - fion( 3, ia )
     END DO
     !
     !
     pes(image) = etot  ! energy already in hartree
     !
     ! ... input values are restored at the end of each iteration
     !
     ethr = 0.D0
     !
  END DO
  !
  outdir  = outdir_saved
  restart_mode  = restart_mode_saved
  scradir = './'
  !
  suspended_image = 0
  !
  DEALLOCATE( tau, fion )
  !
  RETURN
  !
END SUBROUTINE compute_scf
