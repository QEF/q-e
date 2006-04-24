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
SUBROUTINE compute_fes_grads( N_in, N_fin, stat )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : e2
  USE input_parameters,   ONLY : startingwfc, startingpot, diago_thr_init
  USE basis,              ONLY : startingwfc_ => startingwfc, &
                                 startingpot_ => startingpot
  USE metadyn_vars,       ONLY : new_target, to_target, to_new_target, &
                                 dfe_acc, shake_nstep, fe_nstep
  USE path_variables,     ONLY : pos, pes, grad_pes, frozen, &
                                 num_of_images, istep_path, suspended_image
  USE constraints_module, ONLY : lagrange, target, init_constraint, &
                                 deallocate_constraint
  USE control_flags,      ONLY : istep, nstep, history, ethr, conv_ions, ldamped
  USE cell_base,          ONLY : alat, at
  USE ener,               ONLY : etot
  USE ions_base,          ONLY : nat, tau, ityp
  USE path_formats,       ONLY : scf_fmt, scf_fmt_para
  USE io_files,           ONLY : prefix, tmp_dir, iunpath, iunaxsf, &
                                 delete_if_present
  USE constants,          ONLY : bohr_radius_angs
  USE io_global,          ONLY : stdout, ionode, ionode_id, meta_ionode
  USE mp_global,          ONLY : inter_image_comm, intra_image_comm, &
                                 my_image_id, nimage, root
  USE mp,                 ONLY : mp_bcast, mp_barrier, mp_sum, mp_min
  USE check_stop,         ONLY : check_stop_now
  USE path_io_routines,   ONLY : new_image_init, get_new_image, &
                                 stop_other_images  
  USE metadyn_base,       ONLY : add_domain_potential
  USE metadyn_io,         ONLY : write_axsf_file
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: N_in, N_fin
  LOGICAL, INTENT(OUT) :: stat
  INTEGER              :: image
  REAL(DP)             :: tcpu
  CHARACTER(LEN=256)   :: tmp_dir_saved, filename
  LOGICAL              :: lfirst_scf = .TRUE.
  LOGICAL              :: opnd, file_exists
  LOGICAL              :: ldamped_saved
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  REAL(DP),         EXTERNAL :: get_clock
  !
  !
  CALL flush_unit( iunpath )
  !
  IF ( ionode ) THEN
     !
     OPEN( UNIT = iunaxsf, FILE = TRIM( prefix ) // "_" // &
         & TRIM( int_to_char( istep_path + 1 ) ) // ".axsf", &
           STATUS = "UNKNOWN", ACTION = "WRITE" )
     !
     WRITE( UNIT = iunaxsf, FMT = '(" ANIMSTEPS ",I5)' ) num_of_images
     WRITE( UNIT = iunaxsf, FMT = '(" CRYSTAL ")' )
     WRITE( UNIT = iunaxsf, FMT = '(" PRIMVEC ")' )
     WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
          at(1,1) * alat * bohr_radius_angs, &
          at(2,1) * alat * bohr_radius_angs, &
          at(3,1) * alat * bohr_radius_angs
     WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
          at(1,2) * alat * bohr_radius_angs, &
          at(2,2) * alat * bohr_radius_angs, &
          at(3,2) * alat * bohr_radius_angs
     WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
          at(1,3) * alat * bohr_radius_angs, &
          at(2,3) * alat * bohr_radius_angs, &
          at(3,3) * alat * bohr_radius_angs
     !
  END IF
  !
  tmp_dir_saved = tmp_dir
  ldamped_saved = ldamped
  !
  ! ... vectors pes and grad_pes are initalized to zero for all images on
  ! ... all nodes: this is needed for the final mp_sum()
  !
  IF ( my_image_id == root ) THEN
     !
     FORALL( image = N_in:N_fin, .NOT. frozen(image)   )
        !
        grad_pes(:,image) = 0.D0
        !
     END FORALL     
     !
  ELSE
     !
     grad_pes(:,N_in:N_fin) = 0.D0
     !   
  END IF
  !
  ! ... only the first cpu initializes the file needed by parallelization 
  ! ... among images
  !
  IF ( meta_ionode ) CALL new_image_init( N_in, tmp_dir_saved )
  !
  image = N_in + my_image_id
  !
  ! ... all processes are syncronized (needed to have an ordered output)
  !
  CALL mp_barrier()
  !
  fes_loop: DO
     !
     ! ... exit if available images are finished
     !
     IF ( image > N_fin ) EXIT fes_loop
     !     
     suspended_image = image
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
        EXIT fes_loop
        !    
     END IF
     !
     ! ... free-energy gradient ( for non-frozen images only )
     !
     IF ( .NOT. frozen(image) ) THEN
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
        ! ... we read the previous positions for this image from a restart file
        !
        filename = TRIM( tmp_dir ) // "thermodinamic_average.restart"
        !
        INQUIRE( FILE = filename, EXIST = file_exists )
        !
        IF ( file_exists ) THEN
           !
           IF ( ionode ) THEN
              !
              OPEN( UNIT = 1000, FILE = filename )
              !
              READ( 1000, * ) tau
              !
              CLOSE( UNIT = 1000 )
              !
           END IF
           !
           CALL mp_bcast( tau, ionode_id )
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
        ! ... the new values of the order-parameters are set here
        !
        new_target(:) = pos(:,image)
        !
        to_target(:) = ( new_target(:) - target(:) ) / DBLE( shake_nstep )
        !
        ! ... first the system is "adiabatically" moved to the new target
        ! ... by using MD without damping
        !
        ldamped = .FALSE.
        !
        history = 0
        !
        lfirst_scf = .TRUE.
        !
        CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
        CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.update' )
        !
        to_new_target = .TRUE.
        !
        nstep = shake_nstep
        !
        DO istep = 1, shake_nstep
           !
           CALL electronic_scf( lfirst_scf, stat )
           !
           IF ( .NOT. stat ) RETURN
           !
           lfirst_scf = .FALSE.
           !
           CALL move_ions()
           !
        END DO
        !
        ldamped = ldamped_saved
        !
        ! ... then the free energy gradients are computed
        !
        CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
        CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs' )
        !
        to_new_target = .FALSE.
        !
        nstep = fe_nstep
        !
        DO istep = 1, fe_nstep
           !
           CALL electronic_scf( .FALSE., stat )
           !
           IF ( .NOT. stat )  RETURN
           !
           CALL move_ions()
           !
           IF ( ldamped .AND. conv_ions ) EXIT
           !
        END DO
        !
        ! ... the averages are computed here (coverted ot Hartree a.u.)
        !
        IF ( ldamped ) THEN
           !
           ! ... zero temperature
           !
           grad_pes(:,image) = - lagrange(:) / e2
           !
           pes(image) = etot / e2
           !
        ELSE
           !
           ! ... finite temperature
           !
           grad_pes(:,image) = dfe_acc(:) / DBLE( istep ) / e2
           !
        END IF
        !
     END IF
     !
     IF ( ionode ) CALL write_axsf_file( image, tau, alat )
     !
     ! ... the new image is obtained (by ionode only)
     !
     CALL get_new_image( image, tmp_dir_saved )
     !
     CALL mp_bcast( image, ionode_id, intra_image_comm )
     !        
     ! ... input values are restored at the end of each iteration ( they are
     ! ... modified in init_run )
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
  CLOSE( UNIT = iunaxsf )
  !
  CALL add_domain_potential()
  !
  tmp_dir = tmp_dir_saved
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... grad_pes is communicated among "image" pools
     !
     CALL mp_sum( grad_pes(:,N_in:N_fin), inter_image_comm )
     !
  END IF
  !
  ! ... after the lfirst call to compute_scf the input values of startingpot
  ! ... and startingwfc are both set to 'file'
  !
  startingpot = 'file'
  startingwfc = 'file'
  startingpot_ = startingpot
  startingwfc_ = startingwfc
  !
  stat = .TRUE.
  !
  suspended_image = 0
  !
  RETURN  
  !
END SUBROUTINE compute_fes_grads
!
!------------------------------------------------------------------------
SUBROUTINE metadyn()
  !------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE constraints_module, ONLY : target
  USE ions_base,          ONLY : tau
  USE cell_base,          ONLY : alat
  USE io_files,           ONLY : iunaxsf, iunmeta, prefix, tmp_dir
  USE metadyn_vars,       ONLY : fe_grad, metadyn_fmt, to_new_target, &
                                 metadyn_history, max_metadyn_iter,   &
                                 first_metadyn_iter, gaussian_pos, etot_av
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
        ! ... the system is "adiabatically" moved to the new target
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
            iter, target(:), etot_av, gaussian_pos(:), fe_grad(:)
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
      USE control_flags,      ONLY : istep, ldamped, conv_ions, nstep
      USE metadyn_vars,       ONLY : fe_nstep, dfe_acc, etot_av
      USE constraints_module, ONLY : lagrange
      USE io_files,           ONLY : tmp_dir, prefix, delete_if_present
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(INOUT) :: lfirst_scf
      !
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
      nstep = fe_nstep
      !
      to_new_target = .FALSE.
      !
      DO istep = 1, fe_nstep
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
         fe_grad(:) = - lagrange(:) / e2
         !
      ELSE
         !
         ! ... finite temperature
         !
         etot_av = etot_av / DBLE( istep ) / e2
         !
         fe_grad(:) = dfe_acc(:) / DBLE( istep ) / e2
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE free_energy_grad
    !
    !------------------------------------------------------------------------
    SUBROUTINE move_to_target( lfirst_scf )
      !------------------------------------------------------------------------
      !
      USE metadyn_vars,  ONLY : shake_nstep
      USE lsda_mod,      ONLY : lsda
      USE control_flags, ONLY : istep, ldamped, nstep
      USE io_files,      ONLY : tmp_dir, prefix, delete_if_present
      !
      LOGICAL, INTENT(INOUT) :: lfirst_scf
      !
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
      nstep = shake_nstep
      !
      DO istep = 1, shake_nstep
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
  IF ( .NOT. conv_elec ) THEN
     !
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
