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
  ! ... this subroutine is the main scf-driver for NEB 
  ! ... ( called by Modules/neb_base.f90/born_oppenheimer() subroutine ) 
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : if_pos, sp_pos, startingwfc, startingpot, &
                               diago_thr_init
  USE constants,        ONLY : e2
  USE control_flags,    ONLY : conv_elec, istep, history, alpha0, beta0, ethr
  USE check_stop,       ONLY : check_stop_now
  USE cell_base,        ONLY : alat
  USE ions_base,        ONLY : tau, ityp, nat
  USE basis,            ONLY : startingwfc_ => startingwfc, &
                               startingpot_ => startingpot    
  USE ener,             ONLY : etot
  USE force_mod,        ONLY : force
  USE ions_base,        ONLY : if_pos_ => if_pos
  USE extfield,         ONLY : tefield, forcefield
  USE io_files,         ONLY : prefix, tmp_dir, &
                               iunneb, iunupdate, exit_file, iunexit
  USE io_global,        ONLY : stdout
  USE formats,          ONLY : scf_fmt, scf_fmt_para
  USE neb_variables,    ONLY : pos, PES, PES_gradient, num_of_images, &
                               dim, suspended_image, istep_neb, frozen
  USE parser,           ONLY : int_to_char
  USE io_global,        ONLY : ionode
  USE mp_global,        ONLY : inter_image_comm, intra_image_comm, &
                               my_image_id, me_image, root_image, nimage
  USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum, mp_min
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
  INTEGER                    :: image, ia, istat
  REAL (KIND=DP)             :: tcpu 
  CHARACTER (LEN=80)         :: tmp_dir_saved
  LOGICAL                    :: file_exists, opnd   
  REAL(KIND=DP), ALLOCATABLE :: tauold(:,:,:)
    ! previous positions of atoms (needed for extrapolation)
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
  ! ... all processes are syncronized (needed to have an ordered output)
  !
  CALL mp_barrier( intra_image_comm )
  CALL mp_barrier( inter_image_comm )  
  !
  istep = istep_neb + 1
  istat = 0 
  !
  ! ... only the first cpu on each image needs the tauold vector
  !
  IF ( me_image == root_image ) ALLOCATE( tauold( 3, nat, 3 ) )  
  !
  tmp_dir_saved = tmp_dir
  !
  ! ... vectors PES and PES_gradient are initalized to zero for all images on
  ! ... all nodes: this is needed for the final mp_sum()
  !
  IF ( my_image_id == root_image ) THEN
     !
     FORALL( image = N_in:N_fin, ( .NOT. frozen(image) ) )
        !
        PES(image)            = 0.D0
        PES_gradient(:,image) = 0.D0
        !
     END FORALL     
     !
  ELSE
     !
     PES(N_in:N_fin)            = 0.D0
     PES_gradient(:,N_in:N_fin) = 0.D0
     !   
  END IF
  !
  ! ... only the first cpu initializes the file needed by parallelization 
  ! ... among images
  !
  IF ( ( nimage > 1 ) .AND. ionode ) CALL new_image_init()
  !
  image = N_in + my_image_id
  ! 
  scf_loop: DO
     !     
     suspended_image = image
     !
     IF ( check_stop_now( iunneb ) ) THEN
        !   
        istat = 1
        !
        ! ... in case of parallelization on images a stop signal
        ! ... is sent via the "EXIT" file
        !
        IF ( nimage > 1 ) CALL stop_other_images()        
        !
        EXIT scf_loop
        !    
     END IF
     !
     ! ... self-consistency only for non-frozen images
     !
     IF ( .NOT. frozen(image) ) THEN
        !
        tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                  TRIM( int_to_char( image ) ) // "/" 
               
        !
        tcpu = get_clock( 'PWSCF' )
        !
        IF ( nimage > 1 ) THEN
           !
           WRITE( UNIT = iunneb, FMT = scf_fmt_para ) my_image_id, tcpu, image
           !
        ELSE
           !
           WRITE( UNIT = iunneb, FMT = scf_fmt ) tcpu, image
           !
        END IF   
        !
        CALL clean_pw(.true.)
        !
        CALL close_files()
        !
        ! ... unit stdout is connected to the appropriate file
        !
        IF ( me_image == root_image ) THEN
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
        IF ( tefield .AND. .NOT. ALLOCATED( forcefield ) ) &
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
        ! ... initialization of the scf calculation
        !
        CALL init_run()
        !
        IF ( me_image == root_image ) THEN 
           !     
           ! ... the file containing old positions is opened 
           ! ... ( needed for extrapolation )
           !
           CALL seqopn( iunupdate, TRIM( prefix ) // '.update', &
                        'FORMATTED', file_exists ) 
           !
           IF ( file_exists ) THEN
              !
              READ( UNIT = iunupdate, FMT = * ) history
              READ( UNIT = iunupdate, FMT = * ) tauold
              !
           ELSE
              !
              history = 0
              tauold  = 0.D0
              !
              WRITE( UNIT = iunupdate, FMT = * ) history
              WRITE( UNIT = iunupdate, FMT = * ) tauold
              !
           END IF
           !
           CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
           !
        END IF
        !
        CALL mp_bcast( history, root_image, intra_image_comm )
        !
        IF ( conv_elec ) THEN
           !
           ! ... potential and wavefunctions are extrapolated only if
           ! ... scf on the previous image was achieved :
           !
           IF ( me_image == root_image ) THEN 
              !
              ! ... find the best coefficients for the extrapolation of 
              ! ... the potential
              !
              CALL find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )        
              !
           END IF
           !
           CALL mp_bcast( alpha0, root_image, intra_image_comm )
           CALL mp_bcast( beta0,  root_image, intra_image_comm )                
           !
           CALL update_pot()
           !
        END IF
        !
        ! ... self-consistency loop
        !
        CALL electrons()
        !
        IF ( .NOT. conv_elec ) THEN
           !   
           istat = 1
           !
           WRITE( UNIT = iunneb, &
                  FMT = '(/,5X,"WARNING :  scf convergence NOT achieved",/)' )
           !
           ! ... in case of parallelization on images a stop signal
           ! ... is sent via the "EXIT" file
           !
           IF ( nimage > 1 ) CALL stop_other_images()
           !
           EXIT scf_loop
           !
        END IF
        !
        ! ... self-consistent forces
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
        IF ( me_image == root_image ) THEN 
           !
           ! ... save the previous two steps 
           ! ... ( a total of three steps is saved )
           !
           tauold(:,:,3) = tauold(:,:,2)
           tauold(:,:,2) = tauold(:,:,1)
           tauold(:,:,1) = tau(:,:)              
           !
           history = MIN( 3, ( history + 1 ) )
           !
           CALL seqopn( iunupdate, &
                      & TRIM( prefix ) // '.update', 'FORMATTED', file_exists ) 
           !
           WRITE( UNIT = iunupdate, FMT = * ) history
           WRITE( UNIT = iunupdate, FMT = * ) tauold
           !
           CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
           !
        END IF
        !
     END IF   
     !
     ! ... the new image is obtained
     !
     IF ( me_image == root_image ) CALL get_new_image( image )
     !
     CALL mp_bcast( image, root_image, intra_image_comm )
     !        
     ! ... input values are restored at the end of each iteration ( they are
     ! ... modified in init_run )
     !
     startingpot_ = startingpot
     startingwfc_ = startingwfc
     !
     ethr = diago_thr_init     
     !
     CALL reset_k_points()
     !
     ! ... exit if finished
     !
     IF ( image > N_fin ) EXIT scf_loop     
     !
  END DO scf_loop
  !
  IF ( me_image == root_image ) DEALLOCATE( tauold )
  !
  tmp_dir = tmp_dir_saved
  !
  IF ( nimage > 1 ) THEN
     !
     WRITE( UNIT = iunneb, &
            FMT = '(/"image ",I2," waiting at the barrier")' ) my_image_id
     !
     CALL mp_barrier( intra_image_comm )
     CALL mp_barrier( inter_image_comm )
     !
     ! ... PES and PES_gradient are communicated among "image" pools
     !
     CALL mp_sum( PES(N_in:N_fin),            inter_image_comm )
     CALL mp_sum( PES_gradient(:,N_in:N_fin), inter_image_comm )
     CALL mp_sum( istat,                      inter_image_comm )
     !
  END IF
  !
  ! ... global status is computed here
  !
  IF ( istat == 0 ) THEN
     !
     stat = .TRUE.
     !
     suspended_image = 0
     !
  ELSE
     !
     stat = .FALSE.
     !
     IF ( nimage > 1 ) THEN
        !
        CALL mp_min( suspended_image, inter_image_comm )
        !
        OPEN(  UNIT = iunexit, FILE = TRIM( exit_file ) )
        CLOSE( UNIT = iunexit, STATUS = 'DELETE' )
        !
     END IF   
     !     
  END IF      
  !
  ! ... afetr the first call to compute_scf the input values of startingpot
  ! ... and startingwfc are both set to 'file'
  !
  startingpot = 'file'
  startingwfc = 'file'
  startingpot_ = startingpot
  startingwfc_ = startingwfc
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE new_image_init()
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine initializes the file needed for the 
       ! ... parallelization among images
       !
       USE io_files,  ONLY : iunnewimage
       USE mp_global, ONLY : nimage
       !
       IMPLICIT NONE       
       !
       !
       OPEN( UNIT = iunnewimage, FILE = TRIM( tmp_dir_saved ) // &
           & TRIM( prefix ) // '.newimage' , STATUS = 'UNKNOWN' )
       !
       WRITE( iunnewimage, * ) N_in + nimage
       ! 
       CLOSE( UNIT = iunnewimage, STATUS = 'KEEP' )       
       !
       RETURN
       !
     END SUBROUTINE new_image_init
     !
     !-----------------------------------------------------------------------
     SUBROUTINE get_new_image( image )
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine is used to get the new image to work on
       ! ... the "prefix.BLOCK" file is needed to avoid (when present) that 
       ! ... other jobs try to read/write on file "prefix.newimage" 
       !
       USE io_files,  ONLY : iunnewimage, iunblock
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(OUT) :: image
       INTEGER              :: ioerr
       LOGICAL              :: opened, exists
       !
       !
       IF ( nimage > 1 ) THEN
          !
          open_loop: DO
             !          
             OPEN( UNIT = iunblock, FILE = TRIM( tmp_dir_saved ) // &
                 & TRIM( prefix ) // '.BLOCK' , IOSTAT = ioerr, STATUS = 'NEW' )
             !
             IF ( ioerr > 0 ) CYCLE open_loop
             !
             INQUIRE( UNIT = iunnewimage, OPENED = opened )
             !
             IF ( .NOT. opened ) THEN
                !
                INQUIRE( FILE = TRIM( tmp_dir_saved ) // &
                       & TRIM( prefix ) // '.newimage', EXIST = exists )
                !
                IF ( exists ) THEN
                   !
                   OPEN( UNIT = iunnewimage, FILE = TRIM( tmp_dir_saved ) // &
                       & TRIM( prefix ) // '.newimage' , STATUS = 'OLD' )
                   !
                   READ( iunnewimage, * ) image
                   !
                   CLOSE( UNIT = iunnewimage, STATUS = 'DELETE' )
                   !
                   OPEN( UNIT = iunnewimage, FILE = TRIM( tmp_dir_saved ) // &
                       & TRIM( prefix ) // '.newimage' , STATUS = 'NEW' )
                   !
                   WRITE( iunnewimage, * ) image + 1
                   ! 
                   CLOSE( UNIT = iunnewimage, STATUS = 'KEEP' )
                   !
                   EXIT open_loop
                   !
                END IF
                !
             END IF
             !
          END DO open_loop
          !
          CLOSE( UNIT = iunblock, STATUS = 'DELETE' )
          !
       ELSE
          !
          image = image + 1
          !
       END IF      
       !
       RETURN
       !
     END SUBROUTINE get_new_image
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stop_other_images()
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine is used to send a stop signal to other images
       ! ... this is done by creating the exit_file on the working directory
       !
       USE io_files,  ONLY : iunexit, exit_file
       !
       IMPLICIT NONE
       !
       !
       OPEN( UNIT = iunexit, FILE = TRIM( exit_file ) )
       CLOSE( UNIT = iunexit, STATUS = 'KEEP' )               
       !
       RETURN       
       !
     END SUBROUTINE stop_other_images
     !
END SUBROUTINE compute_scf
