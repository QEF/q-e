    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_scf( N_in, N_fin, stat  )
      !-----------------------------------------------------------------------
      !
      USE kinds
      USE input_parameters, ONLY : if_pos, sp_pos, startingwfc, startingpot
      USE constants,        ONLY : e2
      USE control_flags,    ONLY : conv_elec, ethr
      USE check_stop,       ONLY : check_stop_now
      USE brilz,            ONLY : alat
      USE basis,            ONLY : tau, ityp, nat, &
                                   startingwfc_ => startingwfc, &
                                   startingpot_ => startingpot    
      USE ener,             ONLY : etot
      USE force_mod,        ONLY : force
      USE relax,            ONLY : if_pos_ => if_pos
      USE extfield,         ONLY : tefield, forcefield
      USE io_files,         ONLY : prefix, tmp_dir, &
                                   iunneb, iunexit, exit_file
      USE io_global,        ONLY : stdout
      USE formats,          ONLY : scf_fmt
      USE neb_variables,    ONLY : pos, PES, PES_gradient, num_of_images, &
                                   dim, suspended_image
      USE parser,           ONLY : int_to_char
#if defined (__PARA)
      USE para,             ONLY : me, mypool
      USE mp,               ONLY : mp_barrier
#endif        
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN)  :: N_in, N_fin
      LOGICAL, INTENT(OUT) :: stat      
      !
      ! ... local variables definition
      !
      INTEGER              :: image, ia
      REAL (KIND=DP)       :: tcpu 
      CHARACTER (LEN=80)   :: tmp_dir_saved
      LOGICAL              :: file_exists, opnd 
      !
      ! ... end of local variables definition
      !
      ! ... external functions definition
      !
      REAL (KIND=DP), EXTERNAL :: get_clock
      !
      ! ... end of external functions definition
      !
      !
      stat = .TRUE.
      !
      tmp_dir_saved = tmp_dir
      ! 
      DO image = N_in, N_fin
         !
         suspended_image = image
         !
         IF( check_stop_now() ) THEN
            !   
            stat = .FALSE.
            !
            RETURN
            !    
         END IF
         !
         tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                   TRIM( int_to_char( image ) ) // "/" 

         tcpu = get_clock( 'PWSCF' )
         !
         WRITE( UNIT = iunneb, FMT = scf_fmt ) tcpu, image
         !
         CALL clean_pw()
         !
         CALL close_files()
         !
         ! ... unit stdout is connected to the appropriate file
         !
#if defined (__PARA)
         IF ( me == 1 .AND. mypool == 1 ) THEN
#endif  
            INQUIRE( UNIT = stdout, OPENED = opnd )
            IF ( opnd ) CLOSE( UNIT = stdout )
            OPEN( UNIT = stdout, FILE = TRIM( tmp_dir )//'PW.out', &
                  STATUS = 'UNKNOWN', POSITION = 'APPEND' )
#if defined (__PARA)
         END IF 
#endif 
         !
         IF ( .NOT. ALLOCATED( tau ) )      ALLOCATE( tau( 3, nat ) )
         IF ( .NOT. ALLOCATED( ityp ) )     ALLOCATE( ityp( nat ) )
         IF ( .NOT. ALLOCATED( force ) )    ALLOCATE( force( 3, nat ) )  
         IF ( .NOT. ALLOCATED( if_pos_ ) )  ALLOCATE( if_pos_( 3, nat ) )
         IF ( tefield  .AND. .NOT. ALLOCATED( forcefield ) ) &
                                           ALLOCATE( forcefield( 3, nat ) )
         !
         DO ia = 1, nat
            !
            ! ... tau is in alat units ( pos is in bohr )
            !
            tau(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),image) / alat
            !
         END DO
         !
         if_pos_(:,:) = if_pos(:,1:nat) 
         ityp(:)      = sp_pos(1:nat)
         !
         CALL init_run()
         !
         ! ... potential and wavefunctions are extrapolated
         !
         CALL update_pot()
         !
         CALL electrons()
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
         CALL forces()
         !
         ! ... energy is converted from rydberg to hartree
         !
         PES(image) = etot / e2
         !
         ! ... gradients are converted from ( rydberg / bohr ) 
         ! ... to ( hartree / bohr )
         !
         PES_gradient(:,image) = - RESHAPE( SOURCE = force, &
                                            SHAPE = (/ dim /) ) / e2 
         !
         ! ... input values are restored at the end of each iteration
         !
         ethr = 0.D0
         startingpot_ = startingpot
         startingwfc_ = startingwfc
         !
         CALL reset_k_points()
         !
      END DO         
      !
      tmp_dir = tmp_dir_saved
      !
      suspended_image = 0
      !
      RETURN
      !
    END SUBROUTINE compute_scf  
