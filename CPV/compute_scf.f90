    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_scf( N_in, N_fin, stat  )
      !-----------------------------------------------------------------------
      !
      USE kinds
      USE input_parameters,  ONLY : if_pos, sp_pos, rd_pos, tscal, ion_positions
      USE input_parameters,  ONLY : outdir, prefix, nat, restart_mode
      USE input_parameters,  ONLY : scradir => tscradir_inp, ndr
      USE constants,         ONLY : e2
      USE control_flags,     ONLY : conv_elec, ethr
      USE io_files,          ONLY : iunneb, iunexit
      USE io_global,         ONLY : stdout
      USE formats,           ONLY : scf_fmt
      USE neb_variables,     ONLY : pos, PES, PES_gradient, num_of_images, &
                                    dim, suspended_image
      USE parser,        ONLY : int_to_char
      USE mp_global,         ONLY : mpime, my_pool_id
      USE mp,                ONLY : mp_barrier
      USE check_stop,       ONLY : check_stop_now
      !USE main_module,       ONLY : cpmain
      USE restart,           ONLY : check_restartfile
      !
      IMPLICIT NONE
      !
      ! ... local variables definition
      !
      INTEGER, INTENT(IN)    :: N_in, N_fin
      LOGICAL, INTENT(OUT)   :: stat
      INTEGER                :: image
      REAL (KIND=DP)         :: tcpu 
      CHARACTER (LEN=80)     :: outdir_saved, restart_mode_saved
      LOGICAL                :: file_exists, opnd, tstop 

      REAL(dbl), ALLOCATABLE :: tau( :, : )
      REAL(dbl), ALLOCATABLE :: fion( :, : )
      REAL(dbl) :: etot

      INTEGER :: ia, is, isa, ipos
      !
      ! ... end of local variables definition
      !
      ! ... external functions definition
      !
      ! ... end of external functions definition
      !

      stat = .TRUE.
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
         suspended_image = image
         !
         tstop = check_stop_now()
         stat  = .NOT. tstop
         !
         IF( tstop ) THEN
           RETURN
         END IF
         !
         !
         outdir  = TRIM( outdir_saved ) // "/" // TRIM( prefix ) // "_" // &
                    TRIM( int_to_char( image ) ) // "/" 
         scradir = outdir
         !
         WRITE( UNIT = iunneb, FMT = scf_fmt ) tcpu, image
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
         !
         DO ia = 1, nat
            !
            ! ... rd_pos already in bohr units
            !
            rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),image) 
            !
         END DO
  
         tscal         = .FALSE.
         ion_positions = 'from_input'
         IF( check_restartfile( outdir, ndr ) ) THEN
           WRITE( 6, * ) ' restarting calling readfile '
           restart_mode = 'restart'
         ELSE
           WRITE( 6, * ) ' restarting from scratch '
           restart_mode = 'from_scratch'
         END IF

         !
         ! perform an electronic minimization using CPMAIN
         !
         ! CALL cpmain( tau, fion, etot )

         call cprmain( tau(1,1), fion(1,1), etot)
 
         !WRITE(*,*) '++++++++++++++++++++++++++++++++++++'
         !WRITE(*,*) tau
         !WRITE(*,*) '++++++++++++++++++++++++++++++++++++'
         !WRITE(*,*) fion
         !WRITE(*,*) '++++++++++++++++++++++++++++++++++++'
         !WRITE(*,*) etot
         !WRITE(*,*) '++++++++++++++++++++++++++++++++++++'
     
         !
         IF ( mpime == 0 .AND. my_pool_id == 0 ) THEN
            INQUIRE( UNIT = stdout, OPENED = opnd )
            IF ( opnd ) CLOSE( UNIT = stdout )
         END IF 
         !
         IF ( .NOT. conv_elec ) THEN
            !
            WRITE( iunneb, '(/,5X,"WARNING :  scf convergence NOT achieved",/, &
                             & 5X,"stopping in compute_scf()...",/)' )
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
           PES_gradient( 1 + ( ia - 1 ) * 3, image ) = - fion( 1, ia )
           PES_gradient( 2 + ( ia - 1 ) * 3, image ) = - fion( 2, ia )
           PES_gradient( 3 + ( ia - 1 ) * 3, image ) = - fion( 3, ia )
           !write(6,*) fion(1,ia), fion(2,ia), fion(3,ia)
         END DO
         !
         !
         PES(image) = etot  ! energy already in hartree
         !
         ! ... input values are restored at the end of each iteration
         !
         ethr = 0.D0
         !
         !
      END DO         
      !
      outdir  = outdir_saved
      restart_mode  = restart_mode_saved
      scradir = './'
      !
      suspended_image = 0
      DEALLOCATE( tau, fion )
      !
      RETURN
      !
    END SUBROUTINE compute_scf  
