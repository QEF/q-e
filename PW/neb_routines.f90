!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE neb_routines
  !-----------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the NEB implementation into the PWSCF code
  ! ... Written by Carlo Sbraccia ( 04-11-2003 )
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE constants,  ONLY : AU, BOHR_RADIUS_ANGS, eV_to_kelvin
  USE neb_base,   ONLY : initialize_neb, compute_action, compute_tangent, &
                         elastic_constants, gradient, search_stationary_points, &
                         compute_error, path_tangent_

  !
  PRIVATE
  !
  PUBLIC :: initialize_neb, search_mep
  !   
  CONTAINS
    !
    !    
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE control_flags, ONLY : time_max, istep, nstep
      USE io_files,      ONLY : iunneb, iunexit, exit_file     
      USE formats,       ONLY : run_output, run_output_T_const
      USE neb_variables, ONLY : num_of_images, dim, pos, PES, error,       &
                                climbing, optimization,  CI_scheme,        &
                                Emax_index, temp, Emax, neb_thr, conv_neb, &
                                suspended_image, lsteep_des, lquick_min ,  &
                                ldamped_dyn, lmol_dyn 
      USE io_routines,   ONLY : write_restart, write_dat_files, write_output 
#if defined (__LANGEVIN)   
      USE parser,        ONLY : int_to_char
#endif
#if defined (__PARA)
      USE para,          ONLY : me, mypool
      USE mp,            ONLY : mp_barrier      
#endif
      USE minimization_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL (KIND=DP)      :: err, tcpu
      INTEGER             :: image
      LOGICAL             :: stat
      INTEGER             :: N_in, N_fin
      LOGICAL             :: file_exists
#if defined (__LANGEVIN)
      REAL (KIND=DP)      :: langevin_action( num_of_images )
      INTEGER, PARAMETER  :: iunlangevin = 666
      CHARACTER(LEN=20)   :: langevin_fmt
#endif      
      !
      ! ... external functions
      !
      REAL (kind=DP), EXTERNAL :: get_clock
      !
      !
#if defined (__LANGEVIN)  
      !
      langevin_fmt = "(I3," // TRIM( int_to_char( num_of_images - 2 ) ) // &
                   & "(F8.5))"   
      !
      OPEN( UNIT = iunlangevin, FILE = "langevin", STATUS = "UNKNOWN" )
#endif
      !
      conv_neb = .FALSE.
      !
      IF ( istep == nstep ) THEN
         !
         CALL compute_tangent()
         ! 
         CALL gradient()
         !
         CALL compute_error( err )
         !
         CALL write_dat_files()
         !
         suspended_image = 0
         !
         RETURN 
         !  
      END IF 
      ! 
      IF ( optimization ) THEN
         !
         N_in  = 1
         N_fin = num_of_images
         !     
      ELSE
         !
         N_in  = 2
         N_fin = num_of_images - 1
         ! 
      END IF
      !
      minimization: DO
         !
         ! ... the restart file is written 
         !
         CALL write_restart()
         !
         ! ... minimization step is done only in case of no suspended images
         ! ... when the simulation is started from scratch all gradients are 
         ! ... zero and the elastic band remains in the old position.
         !
         IF ( suspended_image == 0 ) THEN
            !
            first_minimization_loop: DO image = N_in, N_fin
               !
               IF ( lsteep_des ) THEN
                  !
                  CALL steepest_descent( image )
                  !
               ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
                  !
                  CALL velocity_Verlet_first_step( image )
                  !
               END IF
               !
            END DO first_minimization_loop 
            !
         END IF
         !
         ! ... the programs checks if the user has required a soft
         ! ... exit or if if maximum CPU time has been exceeded
         !
         INQUIRE( FILE = TRIM( exit_file ), EXIST = file_exists )
         !
         tcpu = get_clock( 'PWSCF' )
         !
         IF ( file_exists ) THEN
            !
#if defined (__PARA)
            !
            ! ... all jobs are syncronized
            !
            CALL mp_barrier()
            !
            IF ( me == 1 .AND. mypool == 1 ) THEN
#endif
               OPEN( UNIT = iunexit, FILE = TRIM( exit_file ), STATUS = "OLD" )
               CLOSE( UNIT = iunexit, STATUS = "DELETE" )
#if defined (__PARA)
            END IF
#endif       
            !  
            WRITE( UNIT = iunneb, &
                   FMT = '(/,5X,"WARNING :  soft exit required",/, &
                   & 5X,"stopping in searc_mep() ...",/)' )  
            !  
            CALL stop_pw( .FALSE. )
            !
         ELSE IF ( tcpu > time_max ) THEN
            ! 
            WRITE( UNIT = iunneb, &
                   FMT = '(5X,"Maximum CPU time exceeded",2F15.2)') &
                tcpu, time_max
            !  
            CALL stop_pw( .FALSE. )
            !
         END IF
         !
         ! ... energies and gradients acting on each image of the elastic band 
         ! ... are computed calling a driver that perfoms the scf calculations
         !
         IF ( istep == 0 ) THEN
            !
            CALL born_oppenheimer_PES( .TRUE., stat )
            !
         ELSE
            !
            CALL born_oppenheimer_PES( optimization, stat )
            !
         END IF
         !
         IF ( .NOT. stat ) THEN
            !
            EXIT minimization
            !
         END IF
         !
         IF ( CI_scheme == "highest-TS" ) THEN
            ! 
            climbing = .FALSE.
            !  
            climbing(Emax_index) = .TRUE.
            !
         ELSE IF ( CI_scheme == "all-SP" ) THEN
            !
            CALL search_stationary_points()
            !
         END IF
         !
         CALL compute_tangent()
         CALL gradient()
         ! 
         ! ... a second minimization step is needed for those algorithms 
         ! ... based on a velocity Verlet scheme
         !
         IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            second_minimization_loop: DO image = N_in, N_fin 
               !
               IF ( lquick_min ) THEN
                  !
                  CALL quick_min_second_step( image )
                  !
               ELSE IF ( ldamped_dyn .OR. lmol_dyn ) THEN
                  !    
                  CALL velocity_Verlet_second_step( image )
                  !
               END IF
               !
            END DO second_minimization_loop
            !
            IF ( lmol_dyn ) CALL thermalization( N_in , N_fin   )
            !
         END IF
         !
         CALL compute_error( err )
         !
         CALL write_dat_files()
         !
#if defined (__LANGEVIN)         
         CALL compute_action( langevin_action )
#endif                  
         !
         istep = istep + 1
         !
         ! ... informations are written on the standard output
         !
#if defined (__LANGEVIN)         
         WRITE( UNIT = iunlangevin, FMT = langevin_fmt ) &
             istep, &
             langevin_action(2:num_of_images-1) * ( AU * BOHR_RADIUS_ANGS )
         !
         WRITE( UNIT = iunneb, FMT = '(/)' )
         !
         DO image = 2, ( num_of_images - 1 )
            !
            WRITE( UNIT = iunneb, &
                   FMT = '(5X,"image = ", I2, "   action = ",F14.8)' ) &
                image, langevin_action(image)  * ( AU * BOHR_RADIUS_ANGS )
            !
         END DO   
         !
         WRITE( UNIT = iunneb, FMT = '(/,5X,"langevin action = ",F14.8)' ) &
             SUM( langevin_action )  * ( AU * BOHR_RADIUS_ANGS ) 
#endif         
         !
         IF ( lmol_dyn ) THEN
            !
            WRITE( UNIT = iunneb, FMT = run_output_T_const ) &
                istep,                    &
                temp * AU * eV_to_kelvin, &
                err * ( AU / BOHR_RADIUS_ANGS ) 
            !
         ELSE
            !            
            WRITE( UNIT = iunneb, FMT = run_output ) &
                istep,                  &
                ( Emax - PES(1) ) * AU, &
                err * ( AU / BOHR_RADIUS_ANGS )
            !   
         END IF
         !
         CALL write_output()
         !
         ! ... the program checks if the convergence has been achieved
         ! 
         IF ( ( err * AU / BOHR_RADIUS_ANGS ) <= neb_thr )  THEN
            !
            conv_neb = .TRUE.
            !
            EXIT minimization
            !  
         END IF
         !
         ! ... the programs checks if the maximum number of iterations has 
         ! ... been reached or if the user has required a soft exit
         !
         IF ( istep >= nstep ) THEN
            !
            suspended_image = 0
            !
            EXIT minimization
            !  
         END IF 
         !
         suspended_image = 0
         ! 
      END DO minimization
      !
#if defined (__LANGEVIN)  
      CLOSE( UNIT = iunlangevin )
#endif      
      !
      RETURN
      !
    END SUBROUTINE search_mep
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_PES( flag, stat )
      !-----------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, Emax_index, Emin, Emax, &
                                PES, PES_gradient, suspended_image
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      LOGICAL, INTENT(IN)   :: flag
      LOGICAL, INTENT(OUT)  :: stat            
      !
      ! ... local variables
      !
      INTEGER               :: i, image
      INTEGER               :: N_in, N_fin
      ! 
      ! ... end of local variables
      !
      !
      IF ( flag ) THEN
         !
         N_in = 1
         N_fin = num_of_images
         !
      ELSE
         !
         N_in = 2
         N_fin = ( num_of_images - 1 )
         !
      END IF
      !
      IF ( suspended_image /= 0 ) N_in = suspended_image
      !
      Emin = + 1.0D16
      Emax = - 1.0D16
      !
      CALL compute_scf( N_in, N_fin, stat )  
      !
      IF ( .NOT. stat ) RETURN
      !
      DO image = 1, num_of_images
         !  
         IF ( PES(image) <= Emin ) Emin = PES(image)
         !
         IF ( PES(image) >= Emax ) THEN
            !
            Emax = PES(image)         
            !
            Emax_index = image
            !
         END IF
         !
      END DO 
      !
      RETURN
      !
    END SUBROUTINE born_oppenheimer_PES
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_scf( N_in, N_fin, stat  )
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : if_pos, sp_pos, startingwfc, startingpot
      USE constants,        ONLY : e2
      USE control_flags,    ONLY : time_max, conv_elec, ethr
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
         INQUIRE( FILE = TRIM( exit_file ), EXIST = file_exists )
         !
         tcpu = get_clock( 'PWSCF' )
         !
         IF ( file_exists ) THEN
            !
#if defined (__PARA)
            !
            ! ... all jobs are syncronized
            !
            CALL mp_barrier()
            !
            IF ( me == 1 .AND. mypool == 1 ) THEN
#endif    
               OPEN( UNIT = iunexit, FILE = TRIM( exit_file ), STATUS = "OLD" )
               CLOSE( UNIT = iunexit, STATUS = "DELETE" )
#if defined (__PARA)
            END IF
#endif       
            !  
            WRITE( iunneb, '(/,5X,"WARNING :  soft exit required",/, &
                             & 5X,"stopping in compute_scf()...",/)' )
            !   
            stat = .FALSE.
            !
            RETURN
            !    
         ELSE IF ( tcpu > time_max ) THEN
            !
            WRITE( iunneb, '(5X,"Maximum CPU time exceeded",2F15.2)' ) &
                tcpu, time_max
            !
            stat = .FALSE.
            !
            RETURN
            !       
         END IF
         !
         tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                   TRIM( int_to_char( image ) ) // "/" 
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
    !
END MODULE neb_routines
