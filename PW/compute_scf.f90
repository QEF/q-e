!
! Copyright (C) 2002-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE compute_scf( N_in, N_fin, stat  )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : if_pos, sp_pos, startingwfc, startingpot
  USE constants,        ONLY : e2
  USE control_flags,    ONLY : conv_elec, istep, alpha0, beta0
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
                               dim, suspended_image, istep_neb
  USE parser,           ONLY : int_to_char
  USE para,             ONLY : me, mypool
  USE io_global,        ONLY : ionode_id
  USE mp_global,        ONLY : inter_pool_comm, intra_pool_comm    
  USE mp,               ONLY : mp_bcast, mp_barrier 
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
  REAL(KIND=DP), ALLOCATABLE :: tauold(:,:,:)
    ! previous positions of atoms  
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
  istep = istep_neb
  !
  stat = .TRUE.
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( me == 1 .AND. mypool == 1 ) THEN
     !
     ALLOCATE( tauold( 3, nat, 3 ) )
     !
  END IF    
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
     IF ( me == 1 .AND. mypool == 1 ) THEN
        !  
        INQUIRE( UNIT = stdout, OPENED = opnd )
        IF ( opnd ) CLOSE( UNIT = stdout )
        OPEN( UNIT = stdout, FILE = TRIM( tmp_dir ) // 'PW.out', &
              STATUS = 'UNKNOWN', POSITION = 'APPEND' )
        !
     END IF
     !
     IF ( .NOT. ALLOCATED( tau ) )      ALLOCATE( tau( 3, nat ) )
     IF ( .NOT. ALLOCATED( ityp ) )     ALLOCATE( ityp( nat ) )
     IF ( .NOT. ALLOCATED( force ) )    ALLOCATE( force( 3, nat ) )  
     IF ( .NOT. ALLOCATED( if_pos_ ) )  ALLOCATE( if_pos_( 3, nat ) )
     IF ( tefield  .AND. .NOT. ALLOCATED( forcefield ) ) &
                                        ALLOCATE( forcefield( 3, nat ) )
     !
     ! ... tau is in alat units ( pos is in bohr )
     !     
     tau = RESHAPE( SOURCE = pos(:,image), &
                    SHAPE = SHAPE( tau ) ) / alat
     !
     if_pos_(:,:) = if_pos(:,1:nat) 
     ityp(:)      = sp_pos(1:nat)
     !
     CALL init_run()
     !
     IF ( me == 1 .AND. mypool == 1 ) THEN 
        !     
        ! ... the file containing old positions is opened 
        ! ... ( needed for extrapolation )
        !
        CALL seqopn( 4, TRIM( prefix ) // '.update', 'FORMATTED', file_exists ) 
        !
        IF ( file_exists ) THEN
           !
           READ( UNIT = 4, FMT = * ) tauold
           !
        ELSE
           !
           tauold = 0.D0
           !
        END IF
        !
        CLOSE( UNIT = 4, STATUS = 'KEEP' )  
        !
        ! ... find the best coefficients for the extrapolation of the potential
        !
        CALL find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )        
        !
     END IF
     !
     IF ( me == 1 ) CALL mp_bcast( alpha0, ionode_id, inter_pool_comm )
     IF ( me == 1 ) CALL mp_bcast( beta0,  ionode_id, inter_pool_comm )
     !
     CALL mp_bcast( alpha0, ionode_id, intra_pool_comm )
     CALL mp_bcast( beta0,  ionode_id, intra_pool_comm )     
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
     IF ( me == 1 .AND. mypool == 1 ) THEN 
        !
        ! ... save the previous two steps ( a total of three steps is saved )
        !
        tauold(:,:,3) = tauold(:,:,2)
        tauold(:,:,2) = tauold(:,:,1)
        tauold(:,:,1) = tau(:,:)              
        !
        CALL seqopn( 4, TRIM( prefix ) // '.update', 'FORMATTED', file_exists ) 
        !
        WRITE( UNIT = 4, FMT = * ) tauold
        !
        CLOSE( UNIT = 4, STATUS = 'KEEP' )
        !
     END IF   
     !
     ! ... input values are restored at the end of each iteration
     !
     startingpot_ = startingpot
     startingwfc_ = startingwfc
     !
     CALL reset_k_points()
     !
  END DO
  !
  IF ( me == 1 .AND. mypool == 1 ) THEN
     !
     DEALLOCATE( tauold )
     !
  END IF      
  !
  tmp_dir = tmp_dir_saved
  !
  suspended_image = 0
  !
  RETURN
  !
END SUBROUTINE compute_scf
