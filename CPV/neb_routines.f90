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

  INTEGER :: nstep_neb
  !
  PUBLIC :: iosys_neb
  PUBLIC :: initialize_neb
  PUBLIC :: search_mep
  !   
  CONTAINS
    !
    !    
    !-----------------------------------------------------------------------
    SUBROUTINE iosys_neb()
      !-----------------------------------------------------------------------
      USE neb_variables, ONLY : lsteep_des, lquick_min , ldamped_dyn, lmol_dyn, &
                            num_of_images_  => num_of_images, &
                            CI_scheme_      => CI_scheme, &
                            VEC_scheme_     => VEC_scheme, &
                            optimization_   => optimization, &
                            damp_           => damp, &
                            temp_req_       => temp_req, &
                            ds_             => ds, &
                            k_max_          => k_max, &
                            k_min_          => k_min, &
                            neb_thr_        => neb_thr
      !
      USE input_parameters, ONLY : CI_scheme, VEC_scheme, minimization_scheme, &
                            num_of_images, optimization, damp, temp_req, &
                            ds, k_max, k_min, neb_thr
      !
      USE input_parameters, ONLY : outdir, prefix, restart_mode, calculation
      USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc
      USE input_parameters, ONLY : rd_ht, trd_ht, cell_symmetry
      USE input_parameters, ONLY : nstep
      !
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: mpime
      USE mp,        ONLY: mp_bcast, mp_barrier, mp_sum
      USE parser, ONLY: int_to_char
      USE cell_base, ONLY: cell_base_init

      IMPLICIT NONE

      INTEGER :: image
      INTEGER :: ios
      CHARACTER(LEN=256) :: outdir_saved
      CHARACTER(LEN=256) :: filename

      INTEGER :: c_mkdir

      IF ( num_of_images < 2 ) THEN
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': num_of_images must be at least 2', 1 )
      END IF
      !
      IF ( ( CI_scheme /= "no-CI"      ) .AND. &
           ( CI_scheme /= "highest-TS" ) .AND. &
           ( CI_scheme /= "all-SP"     ) .AND. &
           ( CI_scheme /= "manual"     ) ) THEN
         !
         CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown CI_scheme', 1 )
         !
      END IF
      IF ( ( VEC_scheme /= "energy-weighted" )   .AND. &
           ( VEC_scheme /= "gradient-weighted" ) ) THEN
         !
         CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown VEC_scheme', 1 )
        !
      END IF

      SELECT CASE ( minimization_scheme )
      !
      CASE ( "sd" )
         lsteep_des  = .TRUE.
         lquick_min  = .FALSE.
         ldamped_dyn = .FALSE.
         lmol_dyn    = .FALSE.
      CASE ( "quick-min" )
         lsteep_des  = .FALSE.
         lquick_min  = .TRUE.
         ldamped_dyn = .FALSE.
         lmol_dyn    = .FALSE.
      CASE ( "damped-dyn" )
         lsteep_des  = .FALSE.
         lquick_min  = .FALSE.
         ldamped_dyn = .TRUE.
         lmol_dyn    = .FALSE.
      CASE ( "mol-dyn" )
         lsteep_des  = .FALSE.
         lquick_min  = .FALSE.
         ldamped_dyn = .FALSE.
         lmol_dyn    = .TRUE.
         IF ( temp_req == 0 ) &
            WRITE( stdout,'(/,T2,"WARNING: tepm_req has not been set" )')
         !
         temp_req = temp_req / ( eV_to_kelvin * AU )
         !
      CASE default
         !
         CALL errore( ' iosys ','calculation=' // TRIM( calculation ) // &
                   & ': unknown minimization_scheme', 1 )
         !
      END SELECT
 
      CALL cell_base_init( ibrav, celldm, trd_ht, cell_symmetry, rd_ht, &
             a, b, c, cosab, cosac, cosbc )

      num_of_images_ = num_of_images
      CI_scheme_     = CI_scheme
      VEC_scheme_    = VEC_scheme
      optimization_  = optimization
      damp_          = damp
      temp_req_      = temp_req
      ds_            = ds
      k_max_         = k_max
      k_min_         = k_min
      neb_thr_       = neb_thr


      nstep_neb = nstep
      nstep     = 100

      outdir_saved = outdir
      !
      DO image = 1, num_of_images
        !
        ios = 0
        outdir = TRIM( outdir_saved ) // TRIM( prefix ) //"_" // &
                 TRIM( int_to_char( image ) ) // '/'
        WRITE(*,*) 'CREATING:', outdir
        !
        IF ( ionode ) THEN
           !
           ! ... a scratch directory for this image of the elastic band is
           ! ... created ( only by the master node )
           !
           ios = c_mkdir( TRIM( outdir ), LEN_TRIM( outdir ) )
           !
        END IF
        !
        ! ... all jobs are syncronized
        !
        CALL mp_barrier()
        !
        ! ... each job checks whether the scratch directory is accessible
        ! ... or not
        !
        filename = TRIM( outdir ) // 'cp' // TRIM( int_to_char( mpime ) )
        !
        OPEN( UNIT = 4, FILE = TRIM( filename ) , &
              STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
        CALL mp_sum( ios )
        !
        IF ( ios /= 0 ) &
           CALL errore( 'outdir: ', TRIM( outdir ) // &
                      & ' non existent or non writable', 1 )
        !
        ! ... if starting from scratch all temporary files are removed
        ! ... from outdir ( only by the master node )
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           IF ( ionode ) THEN
              !
              ! ... standard output of the self consistency is removed
              !
              OPEN( UNIT = 4, FILE = TRIM( outdir ) // 'FPMD.out', &
                    STATUS = 'UNKNOWN' )
              CLOSE( UNIT = 4, STATUS = 'DELETE' )
              !
           END IF
           !
        END IF
        !
      END DO
      !
      outdir = outdir_saved

      RETURN
    END SUBROUTINE
    !    
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE control_flags,  ONLY : istep
      USE io_files,       ONLY : iunneb, iunexit, exit_file     
      USE formats,        ONLY : run_output, run_output_T_const
      USE neb_variables,  ONLY : num_of_images, dim, pos, PES, error, &
                                 climbing, optimization,  CI_scheme, &
                                 Emax_index, temp, Emax, neb_thr, conv_neb, &
                                 suspended_image, lsteep_des, lquick_min , &
                                 ldamped_dyn, lmol_dyn 
      USE io_routines,    ONLY : write_restart, write_dat_files 
      USE mp_global,      ONLY : mpime, my_pool_id
      USE mp,             ONLY : mp_barrier      
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
      LOGICAL             :: tstop
      !
      ! ... end of local variables
      ! ... external functions
      !
      !
      conv_neb = .FALSE.
      !
      IF ( istep == nstep_neb ) THEN
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
         !!tstop = check_stop_now()
         !!!
         !!IF ( tstop ) THEN
         !!   !
         !!   ! ... all jobs are syncronized
         !!   !
         !!   CALL mp_barrier()
         !!   !
         !!   WRITE( UNIT = iunneb, &
         !!          FMT = '(/,5X,"WARNING :  soft exit required",/, &
         !!          & 5X,"stopping in searc_mep() ...",/)' )  
         !!   RETURN
         !!   !
         !!END IF
         !
         !
         ! ... energies and gradients acting on each image of the elastic band 
         ! ... are computed calling a driver that perfoms the scf calculations
         !
         IF ( istep == 0 ) THEN
            !
            CALL born_oppenhimer_PES( .TRUE., stat )
            !
         ELSE
            !
            CALL born_oppenhimer_PES( optimization, stat )
            !
         END IF
         !
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
         !
         CALL compute_tangent()
         !
         !
         CALL gradient()
         !
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
         !
         CALL compute_error( err )
         !
         !
         CALL write_dat_files()
         !
         istep = istep + 1
         !  
         IF ( lmol_dyn ) THEN
            !
            WRITE( UNIT = iunneb, FMT = run_output_T_const ) &
                istep, &
                temp * AU * eV_to_kelvin, &
                err * ( AU / BOHR_RADIUS_ANGS ) 
            !
         ELSE
            !
            WRITE( UNIT = iunneb, FMT = run_output ) &
                istep, &
                ( Emax - PES(1) ) * AU, &
                err * ( AU / BOHR_RADIUS_ANGS )
            !   
         END IF     
         !
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
         !
         ! ... the programs checks if the maximum number of iterations has 
         ! ... been reached or if the user has required a soft exit
         !
         IF ( istep >= nstep_neb ) THEN
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
      RETURN
      !
    END SUBROUTINE search_mep
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE born_oppenhimer_PES( flag, stat )
      !-----------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, Emax_index, Emin, Emax, &
                                PES, PES_gradient, suspended_image
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      LOGICAL, INTENT(IN)        :: flag
      LOGICAL, INTENT(OUT)       :: stat
      INTEGER                    :: i, image
      INTEGER                    :: N_in, N_fin
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
    END SUBROUTINE born_oppenhimer_PES
    !
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
      USE formats,           ONLY : scf_fmt
      USE neb_variables,     ONLY : pos, PES, PES_gradient, num_of_images, &
                                    dim, suspended_image
      USE parser,        ONLY : int_to_char
      USE mp_global,         ONLY : mpime, my_pool_id
      USE mp,                ONLY : mp_barrier
      !USE environment,       ONLY : check_stop_now
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
         !tstop = check_stop_now()
         tstop = .FALSE.
         stat  = .NOT. tstop
         !
         IF( tstop ) THEN
           RETURN
         END IF
         !
         !
         outdir  = TRIM( outdir_saved ) // TRIM( prefix ) // "_" // &
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
            OPEN( UNIT = stdout, FILE = TRIM( outdir )//'FPMD.out', &
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

    !
END MODULE neb_routines 
