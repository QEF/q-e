!
! Copyright (C) 2002-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE compute_fes_grads( fii, lii, stat )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : e2
  USE input_parameters,   ONLY : startingwfc, startingpot, diago_thr_init
  USE basis,              ONLY : startingwfc_ => startingwfc, &
                                 startingpot_ => startingpot
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
     tcpu = get_clock( 'PWSCF' )
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
     ! ... the averages are computed here (coverted ot Hartree a.u.)
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
     startingpot_ = startingpot
     startingwfc_ = startingwfc
     !
     ethr = diago_thr_init
     !
     CALL close_files()
     !
     CALL reset_k_points()
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
     ! ... by reading the postions from the restart-file
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
  startingpot_ = startingpot
  startingwfc_ = startingwfc
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
SUBROUTINE metadyn()
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE constraints_module, ONLY : constr_target
  USE ions_base,          ONLY : tau
  USE cell_base,          ONLY : alat
  USE io_files,           ONLY : iunaxsf, iunmeta, prefix, tmp_dir
  USE metadyn_vars,       ONLY : ncolvar, etot_av, fe_grad, metadyn_fmt, &
                                 to_new_target, metadyn_history,         &
                                 max_metadyn_iter, first_metadyn_iter,   &
                                 gaussian_pos
  USE metadyn_base,       ONLY : add_gaussians, add_domain_potential, &
                                 evolve_collective_vars
  USE metadyn_io,         ONLY : write_axsf_file, write_metadyn_restart
  USE io_global,          ONLY : ionode, stdout
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: dirname
  INTEGER            :: iter
  REAL(DP)           :: norm_fe_grad
  LOGICAL            :: lfirst_scf = .TRUE.
  !
  !
  dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
  !
  iter = first_metadyn_iter
  !
  metadyn_loop: DO
     !
     IF ( iter > 0 ) THEN
        !
        CALL add_gaussians( iter )
        !
        CALL add_domain_potential()
        !
        norm_fe_grad = norm( fe_grad )
        !
        CALL evolve_collective_vars( norm_fe_grad )
        !
        WRITE( stdout, '(/,5X,"adiabatic switch of the system ", &
                            & "to the new coarse-grained positions",/)' )
        !
        ! ... the system is "adiabatically" moved to the new constr_target
        !
        CALL move_to_target( lfirst_scf )
        !
     END IF
     !
     iter = iter + 1
     !
     metadyn_history(:,iter) = gaussian_pos(:)
     !
     IF ( ionode ) CALL write_axsf_file( iter, tau, alat )
     !
     WRITE( stdout, '(/,5X,"calculation of the mean force",/)' )
     !
     CALL free_energy_grad( lfirst_scf )
     !
     IF ( ionode ) THEN
        !
        WRITE( UNIT = iunmeta, FMT = metadyn_fmt ) &
            iter, constr_target(1:ncolvar), etot_av, gaussian_pos(:), fe_grad(:)
        !
        CALL flush_unit( iunmeta )
        CALL flush_unit( iunaxsf )
        !
     END IF
     !
     CALL write_metadyn_restart( dirname, iter, tau, etot_av, alat )
     !
     IF ( iter >= max_metadyn_iter ) EXIT metadyn_loop
     !
  END DO metadyn_loop
  !
  IF ( ionode ) THEN
     !
     CLOSE( UNIT = iunaxsf )
     CLOSE( UNIT = iunmeta )
     !
  END IF
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE free_energy_grad( lfirst_scf )
      !------------------------------------------------------------------------
      !
      USE constants,          ONLY : e2
      USE ener,               ONLY : etot
      USE lsda_mod,           ONLY : lsda
      USE control_flags,      ONLY : ldamped, conv_ions, nstep
      USE metadyn_vars,       ONLY : fe_nstep, eq_nstep, dfe_acc, etot_av
      USE constraints_module, ONLY : lagrange
      USE io_files,           ONLY : tmp_dir, prefix, delete_if_present
      USE io_global,          ONLY : ionode_id
      USE mp_global,          ONLY : intra_image_comm
      USE mp,                 ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(INOUT) :: lfirst_scf
      !
      INTEGER :: i
      LOGICAL :: stat
      !
      !
      etot_av = 0.D0
      dfe_acc = 0.D0
      !
      IF ( lsda ) CALL reset_init_mag()
      !
      CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
      CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs' )
      CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.update' )
      !
      to_new_target = .FALSE.
      !
      nstep = fe_nstep
      !
      DO i = 1, fe_nstep
         !
         CALL electronic_scf( lfirst_scf, stat )
         !
         lfirst_scf = .FALSE.
         !
         IF ( .NOT. stat ) CALL stop_run( stat )
         !
         CALL move_ions()
         !
         IF ( ldamped .AND. conv_ions ) EXIT
         !
      END DO
      !
      ! ... the averages are computed here and converted to Hartree
      !
      IF ( ldamped ) THEN
         !
         ! ... zero temperature
         !
         etot_av = etot / e2
         !
         fe_grad(:) = - lagrange(1:ncolvar) / e2
         !
      ELSE
         !
         ! ... finite temperature
         !
         etot_av = etot_av / DBLE( fe_nstep ) / e2
         !
         fe_grad(:) = dfe_acc(:) / DBLE( fe_nstep - eq_nstep ) / e2
         !
      END IF
      !
      ! ... notice that etot_av and fe_grad have been computed, so far, by
      ! ... ionode only: here we broadcast to all the other cpus
      !
      CALL mp_bcast( etot_av, ionode_id, intra_image_comm )
      CALL mp_bcast( fe_grad, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE free_energy_grad
    !
    !------------------------------------------------------------------------
    SUBROUTINE move_to_target( lfirst_scf )
      !------------------------------------------------------------------------
      !
      USE metadyn_vars,  ONLY : sw_nstep
      USE lsda_mod,      ONLY : lsda
      USE control_flags, ONLY : ldamped, nstep
      USE io_files,      ONLY : tmp_dir, prefix, delete_if_present
      !
      LOGICAL, INTENT(INOUT) :: lfirst_scf
      !
      INTEGER :: i
      LOGICAL :: stat, ldamped_saved
      !
      !
      ldamped_saved = ldamped
      !
      IF ( lsda ) CALL reset_init_mag()
      !
      CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
      CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.update' )
      !
      ldamped = .FALSE.
      !
      to_new_target = .TRUE.
      !
      nstep = sw_nstep
      !
      DO i = 1, sw_nstep
         !
         CALL electronic_scf( lfirst_scf, stat )
         !
         lfirst_scf = .FALSE.
         !
         IF ( .NOT. stat ) CALL stop_run( stat )
         !
         CALL move_ions()
         !
      END DO
      !
      ldamped = ldamped_saved
      !
      RETURN
      !
    END SUBROUTINE move_to_target
    !
END SUBROUTINE metadyn
!
!----------------------------------------------------------------------------
SUBROUTINE reset_init_mag()
  !----------------------------------------------------------------------------
  !
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
