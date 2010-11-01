!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_fes_grads( fii, lii, stat )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : e2
  USE input_parameters,   ONLY : startingwfc, startingpot, diago_thr_init
  USE basis,              ONLY : starting_wfc, starting_pot
  USE metadyn_vars,       ONLY : ncolvar, dfe_acc, new_target, to_target, &
                                 to_new_target, sw_nstep, fe_nstep, eq_nstep
  USE path_variables,     ONLY : grad_fes => grad_pes, &
                                 pos, num_of_images, istep_path, pending_image
  USE constraints_module, ONLY : lagrange, constr_target, init_constraint, &
                                 deallocate_constraint
  USE control_flags,      ONLY : istep, nstep, ethr, conv_ions, ldamped
  USE cell_base,          ONLY : alat, at
  USE ions_base,          ONLY : nat, tau, ityp
  USE path_formats,       ONLY : scf_fmt, scf_fmt_para
  USE io_files,           ONLY : prefix, tmp_dir, iunpath, iunaxsf, &
                                 delete_if_present
  USE constants,          ONLY : bohr_radius_angs
  USE io_global,          ONLY : stdout, ionode, ionode_id, meta_ionode
  USE mp_global,          ONLY : inter_image_comm, intra_image_comm, &
                                 my_image_id, nimage, root_image
  USE mp,                 ONLY : mp_bcast, mp_barrier, mp_sum, mp_min
  USE check_stop,         ONLY : check_stop_now
  USE path_io_routines,   ONLY : new_image_init, get_new_image, &
                                 stop_other_images
  USE metadyn_base,       ONLY : add_domain_potential
  USE metadyn_io,         ONLY : write_axsf_file
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: fii, lii
  LOGICAL, INTENT(OUT) :: stat
  !
  INTEGER            :: i, image
  REAL(DP)           :: tcpu
  CHARACTER(LEN=10)  :: stage
  INTEGER            :: fe_step0, sw_step0
  CHARACTER(LEN=256) :: tmp_dir_saved, filename, basename
  LOGICAL            :: lfirst_scf = .TRUE.
  LOGICAL            :: opnd, file_exists
  LOGICAL            :: ldamped_saved
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  REAL(DP),         EXTERNAL :: get_clock
  !
  !
  CALL flush_unit( iunpath )
  !
  tmp_dir_saved = tmp_dir
  ldamped_saved = ldamped
  !
  ! ... vectors fes and grad_fes are initalized to zero for all images on
  ! ... all nodes: this is needed for the final mp_sum()
  !
  IF ( my_image_id == root_image ) THEN
     !
     grad_fes(:,:) = 0.D0
     !
  ELSE
     !
     grad_fes(:,fii:lii) = 0.D0
     !
  END IF
  !
  ! ... only the first cpu initializes the file needed by parallelization
  ! ... among images
  !
  IF ( meta_ionode ) CALL new_image_init( fii, tmp_dir_saved )
  !
  image = fii + my_image_id
  !
  ! ... all processes are syncronized (needed to have an ordered output)
  !
  CALL mp_barrier()
  !
  fes_loop: DO
     !
     ! ... exit if available images are finished
     !
     IF ( image > lii ) EXIT fes_loop
     !
     pending_image = image
     !
     IF ( check_stop_now( iunpath ) ) THEN
        !
        stat = .FALSE.
        !
        ! ... in case of parallelization on images a stop signal
        ! ... is sent via the "EXIT" file
        !
        IF ( nimage > 1 ) CALL stop_other_images()
        !
        IF ( interrupt_run( stat ) ) RETURN
        !
     END IF
     !
     ! ... calculation of the mean-force
     !
     tcpu = get_clock( 'sm' )
     !
     IF ( nimage > 1 ) THEN
        !
        WRITE( UNIT = iunpath, FMT = scf_fmt_para ) my_image_id, tcpu, image
        !
     ELSE
        !
        WRITE( UNIT = iunpath, FMT = scf_fmt ) tcpu, image
        !
     END IF
     !
     tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // &
               "_" // TRIM( int_to_char( image ) ) // "/"
     !
     basename = TRIM( tmp_dir ) // TRIM( prefix )
     !
     ! ... unit stdout is connected to the appropriate file
     !
     IF ( ionode ) THEN
        !
        INQUIRE( UNIT = stdout, OPENED = opnd )
        IF ( opnd ) CLOSE( UNIT = stdout )
        OPEN( UNIT = stdout, FILE = TRIM( tmp_dir ) // 'PW.out', &
              STATUS = 'UNKNOWN', POSITION = 'APPEND' )
        !
     END IF
     !
     filename = TRIM( tmp_dir ) // "therm_average.restart"
     !
     INQUIRE( FILE = filename, EXIST = file_exists )
     !
     IF ( file_exists ) THEN
        !
        ! ... we read the previous positions, the value of the accumulators,
        ! ... and the number of steps already performed for this image from
        ! ... a restart file
        !
        IF ( ionode ) THEN
           !
           OPEN( UNIT = 1000, FILE = filename )
           !
           READ( 1000, * ) stage
           READ( 1000, * ) tau(:,:)
           READ( 1000, * ) nstep
           READ( 1000, * ) to_target
           READ( 1000, * ) dfe_acc
           !
           CLOSE( UNIT = 1000 )
           !
        END IF
        !
        CALL mp_bcast( stage,     ionode_id, intra_image_comm )
        CALL mp_bcast( tau,       ionode_id, intra_image_comm )
        CALL mp_bcast( nstep,     ionode_id, intra_image_comm )
        CALL mp_bcast( to_target, ionode_id, intra_image_comm )
        CALL mp_bcast( dfe_acc,   ionode_id, intra_image_comm )
        !
     ELSE
        !
        stage = 'tobedone'
        !
     END IF
     !
     CALL clean_pw( .FALSE. )
     CALL deallocate_constraint()
     !
     CALL init_constraint( nat, tau, ityp, alat )
     !
     CALL setup ()
     CALL init_run()
     !
     fe_step0 = 1
     sw_step0 = 1
     !
     SELECT CASE( stage )
     CASE( 'done' )
        !
        ! ... do nothing and recompute the average quantities
        !
     CASE( 'tobedone' )
        !
        new_target(:) = pos(:,image)
        !
        to_target(:) = ( new_target(:) - &
                         constr_target(1:ncolvar) ) / DBLE( sw_nstep )
        !
        dfe_acc = 0.D0
        !
        stage = 'switch'
        !
     CASE( 'switch' )
        !
        dfe_acc = 0.D0
        !
        sw_step0 = nstep
        !
     CASE( 'mean-force' )
        !
        fe_step0 = nstep
        !
     CASE DEFAULT
        !
        CALL errore( 'compute_fes_grads', &
                     'stage ' // TRIM( stage ) // ' unknown', 1 )
        !
     END SELECT
     !
     IF ( stage == 'switch' ) THEN
        !
        ! ... first the collective variables are "adiabatically" changed to
        ! ... the new vales by using MD without damping
        !
        WRITE( stdout, '(/,5X,"adiabatic switch of the system ", &
                            & "to the new coarse-grained positions",/)' )
        !
        CALL delete_if_present( TRIM( basename ) // '.md' )
        CALL delete_if_present( TRIM( basename ) // '.update' )
        !
        ldamped       = .FALSE.
        lfirst_scf    = .TRUE.
        to_new_target = .TRUE.
        !
        nstep = sw_nstep
        !
        DO i = sw_step0, sw_nstep
           !
           CALL electronic_scf( lfirst_scf, stat )
           !
           IF ( interrupt_run( stat ) ) RETURN
           !
           lfirst_scf = .FALSE.
           !
           CALL move_ions()
           !
        END DO
        !
        ldamped = ldamped_saved
        !
        stage = 'mean-force'
        !
        CALL write_restart( 'mean-force', 0 )
        !
     END IF
     !
     IF ( stage == 'mean-force' ) THEN
        !
        ! ... then the free energy gradients are computed
        !
        WRITE( stdout, '(/,5X,"calculation of the mean force",/)' )
        !
        CALL delete_if_present( TRIM( basename ) // '.md' )
        CALL delete_if_present( TRIM( basename ) // '.bfgs' )
        CALL delete_if_present( TRIM( basename ) // '.update' )
        !
        to_new_target = .FALSE.
        !
        nstep = fe_nstep
        !
        DO i = fe_step0, fe_nstep
           !
           CALL electronic_scf( .FALSE., stat )
           !
           IF ( interrupt_run( stat ) ) RETURN
           !
           CALL move_ions()
           !
           IF ( ldamped .AND. conv_ions ) EXIT
           !
        END DO
        !
     END IF
     !
     ! ... the averages are computed here (converted to Hartree a.u.)
     !
     IF ( ldamped ) THEN
        !
        ! ... zero temperature case
        !
        grad_fes(:,image) = - lagrange(1:ncolvar) / e2
        !
     ELSE
        !
        ! ... finite temperature case
        !
        grad_fes(:,image) = dfe_acc(:) / DBLE( fe_nstep - eq_nstep ) / e2
        !
     END IF
     !
     ! ... notice that grad_fes(:,image) have been computed, so far, by
     ! ... ionode only: here we broadcast to all the other cpus
     !
     CALL mp_bcast( grad_fes(:,image), ionode_id, intra_image_comm )
     !
     IF ( ionode ) THEN
        !
        ! ... the restart file is written here
        !
        CALL write_restart( 'done', 0 )
        !
     END IF
     !
     ! ... the new image is obtained (by ionode only)
     !
     CALL get_new_image( image, tmp_dir_saved )
     !
     CALL mp_bcast( image, ionode_id, intra_image_comm )
     !
     ! ... input values are restored at the end of each iteration ( they are
     ! ... modified by init_run )
     !
     starting_pot = startingpot
     starting_wfc = startingwfc
     !
     ethr = diago_thr_init
     !
     CALL close_files()
     !
     CALL reset_k_points ( )
     !
  END DO fes_loop
  !
  CALL mp_barrier()
  !
  IF ( meta_ionode ) THEN
     !
     ! ... when all the images are done the stage is changed from
     ! ... 'done' to 'tobedone'
     !
     DO image = fii, lii
        !
        tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // &
                & "_" // TRIM( int_to_char( image ) ) // "/"
        !
        filename = TRIM( tmp_dir ) // "therm_average.restart"
        !
        OPEN( UNIT = 1000, FILE = filename )
        !
        READ( 1000, * ) stage
        READ( 1000, * ) tau(:,:)
        READ( 1000, * ) nstep
        READ( 1000, * ) to_target
        READ( 1000, * ) dfe_acc
        !
        CLOSE( UNIT = 1000 )
        !
        CALL write_restart( 'tobedone', 0 )
        !
     END DO
     !
     ! ... here the meta_ionode writes the axsf file for this iteration
     ! ... by reading the positions from the restart-file
     !
     filename = TRIM( prefix ) // "_" // &
              & TRIM( int_to_char( istep_path + 1 ) ) // ".axsf"
     !
     OPEN( UNIT = iunaxsf, FILE = filename, ACTION = "WRITE" )
     !
     WRITE( UNIT = iunaxsf, FMT = '(" ANIMSTEPS ",I5)' ) num_of_images
     WRITE( UNIT = iunaxsf, FMT = '(" CRYSTAL ")' )
     WRITE( UNIT = iunaxsf, FMT = '(" PRIMVEC ")' )
     WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
          at(1,1)*alat*bohr_radius_angs, &
          at(2,1)*alat*bohr_radius_angs, &
          at(3,1)*alat*bohr_radius_angs
     WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
          at(1,2)*alat*bohr_radius_angs, &
          at(2,2)*alat*bohr_radius_angs, &
          at(3,2)*alat*bohr_radius_angs
     WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
          at(1,3)*alat*bohr_radius_angs, &
          at(2,3)*alat*bohr_radius_angs, &
          at(3,3)*alat*bohr_radius_angs
     !
     DO image = 1, num_of_images
        !
        tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // &
                & "_" // TRIM( int_to_char( image ) ) // "/"
        !
        filename = TRIM( tmp_dir ) // "therm_average.restart"
        !
        OPEN( UNIT = 1000, FILE = filename )
        !
        READ( 1000, * ) stage
        READ( 1000, * ) tau(:,:)
        !
        CLOSE( UNIT = 1000 )
        !
        CALL write_axsf_file( image, tau, alat )
        !
     END DO
     !
     CLOSE( UNIT = iunaxsf )
     !
  END IF
  !
  CALL add_domain_potential()
  !
  tmp_dir = tmp_dir_saved
  !
  ! ... after the first call to compute_fes_grads the input values of
  ! ... startingpot and startingwfc are both set to 'file'
  !
  startingpot = 'file'
  startingwfc = 'file'
  starting_pot= startingpot
  starting_wfc= startingwfc
  !
  pending_image = 0
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... grad_fes is communicated among "image" pools
     !
     CALL mp_sum( grad_fes(:,fii:lii), inter_image_comm )
     !
  END IF
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_restart( stage, nstep )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN) :: stage
      INTEGER,          INTENT(IN) :: nstep
      !
      OPEN( UNIT = 1000, FILE = filename )
      !
      WRITE( 1000, * ) TRIM( stage )
      WRITE( 1000, * ) tau(:,:)
      WRITE( 1000, * ) nstep
      WRITE( 1000, * ) to_target
      WRITE( 1000, * ) dfe_acc
      !
      CLOSE( UNIT = 1000 )
      !
    END SUBROUTINE write_restart
    !
    !------------------------------------------------------------------------
    FUNCTION interrupt_run( stat )
      !-----------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: stat
      LOGICAL             :: interrupt_run
      !
      interrupt_run = .NOT. stat
      !
      IF ( stat ) RETURN
      !
      pending_image = 1
      !
      filename = TRIM( prefix ) // "_" // &
               & TRIM( int_to_char( istep_path + 1 ) ) // ".axsf"
      !
      CALL delete_if_present( TRIM( filename ) )
      !
      filename = TRIM( tmp_dir ) // "therm_average.restart"
      !
      CALL write_restart( stage, istep - 1 )
      !
      IF ( nimage > 1 ) CALL stop_other_images()
      !
    END FUNCTION interrupt_run
    !
END SUBROUTINE compute_fes_grads
!
!----------------------------------------------------------------------------
SUBROUTINE reset_init_mag()
  !----------------------------------------------------------------------------
  !
  USE dfunct,                 only : newd
  IMPLICIT NONE
  !
  CALL hinit0()
  CALL potinit()
  CALL newd()
  CALL wfcinit()
  !
  RETURN
  !
END SUBROUTINE reset_init_mag
!
!----------------------------------------------------------------------------
SUBROUTINE electronic_scf( lfirst_scf, stat )
  !----------------------------------------------------------------------------
  !
  USE control_flags, ONLY : conv_elec, ethr
  USE io_files,      ONLY : iunpath
  USE io_global,     ONLY : ionode
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN)  :: lfirst_scf
  LOGICAL, INTENT(OUT) :: stat
  !
  !
  IF ( .NOT. lfirst_scf ) THEN
     !
     ethr = 1.D-5
     !
     CALL hinit1()
     !
  END IF
  !
  CALL electrons()
  !
  stat = conv_elec
  !
  IF ( .NOT.conv_elec ) THEN
     !
     IF ( ionode ) &
        WRITE( UNIT = iunpath, &
               FMT = '(/,5X,"WARNING :  scf convergence NOT achieved",/)' )
     !
     RETURN
     !
  END IF
  !
  CALL forces()
  !
  RETURN
  !
END SUBROUTINE electronic_scf
