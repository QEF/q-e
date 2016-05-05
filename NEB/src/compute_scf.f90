!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE compute_scf( fii, lii, stat  )
  !----------------------------------------------------------------------------
  !
  ! ... this subroutine is the main scf-driver for all "path" calculations
  ! ... ( called by Modules/path_base.f90/born_oppenheimer() subroutine )
  !
  ! ... for each image in the path, it performs the self-consistent loop
  ! ... computing the energy and the forces
  !
  ! ... Written by Carlo Sbraccia (2003-2006)
  !
  USE basis,            ONLY : starting_wfc, starting_pot
  USE kinds,            ONLY : DP
  USE constants,        ONLY : e2
  USE control_flags,    ONLY : conv_elec
  USE cell_base,        ONLY : alat
  USE ions_base,        ONLY : tau, nat, ityp, zv
  USE ener,             ONLY : etot, ef
  USE force_mod,        ONLY : force
  USE io_files,         ONLY : prefix, tmp_dir, wfc_dir,  iunupdate, seqopn, &
                               exit_file, iunexit, delete_if_present
  USE path_io_units_module, ONLY : iunpath
  USE path_formats,     ONLY : scf_fmt, scf_fmt_para
  USE path_variables,   ONLY : pos, pes, grad_pes, dim1, pending_image, &
                               istep_path, frozen, num_of_images, &
                               first_last_opt
  USE io_global,        ONLY : stdout, ionode, ionode_id, meta_ionode
  USE mp_global,        ONLY : inter_image_comm, intra_image_comm, &
                               my_image_id, nimage, root_image
  USE mp_world,         ONLY : world_comm
  USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum, mp_min
  USE path_io_routines, ONLY : new_image_init, get_new_image, &
                               stop_other_images
  USE fcp_opt_routines, ONLY : fcp_neb_nelec, fcp_neb_ef
  USE fcp_variables,    ONLY : lfcpopt
  USE klist,            ONLY : nelec, tot_charge
  USE extrapolation,    ONLY : update_neb
  USE funct,            ONLY : stop_exx, dft_is_hybrid
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: fii, lii         ! indexes to first and last images
  LOGICAL, INTENT(OUT) :: stat
  !
  INTEGER               :: fii_, lii_      ! local copies of fii and lii
  INTEGER               :: image, istat
  REAL(DP)              :: tcpu
  CHARACTER (LEN=256)   :: tmp_dir_saved
  LOGICAL               :: file_exists, opnd
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
  fii_ = fii
  lii_ = lii
  !
  istat = 0
  !
  FLUSH( iunpath )
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... vectors pes and grad_pes are initialized to zero for all images on
     ! ... all nodes: this is needed for the final mp_sum()
     !
     IF ( my_image_id == root_image ) THEN
        !
        FORALL( image = fii:lii, .NOT.frozen(image) )
           !
           pes(image)        = 0.D0
           grad_pes(:,image) = 0.D0
           !
        END FORALL
        !
     ELSE
        !
        pes(fii:lii)        = 0.D0
        grad_pes(:,fii:lii) = 0.D0
        !
     END IF
     !
     IF ( lfcpopt ) THEN
        !
        IF ( my_image_id == root_image ) THEN
           !
           FORALL( image = fii:lii, .NOT.frozen(image) )
              !
              fcp_neb_ef(image) = 0.D0
              !
           END FORALL
           !
        ELSE
           !
           fcp_neb_ef(fii:lii) = 0.D0
           !
        END IF
        !
     END IF
     !
  END IF
  !
  ! ... all processes are syncronized (needed to have a readable output)
  !
  CALL mp_barrier( world_comm )
  !
  IF ( nimage > 1 .AND. .NOT.first_last_opt ) THEN
     !
     ! ... self-consistency on the first and last images is done separately
     !
     IF ( fii == 1 ) THEN
        !
        IF ( my_image_id == root_image ) THEN
           !
           IF (dft_is_hybrid()) call stop_exx() 
           CALL do_scf( 1, istat )
           !
           IF ( istat /= 0 ) GOTO 1
           !
        END IF
        !
        fii_ = 2
        !
     END IF
     IF ( lii == num_of_images ) THEN
        !
        IF ( my_image_id == root_image + 1 ) THEN
           !
           IF (dft_is_hybrid()) call stop_exx() 
           CALL do_scf( num_of_images, istat )
           !
           IF ( istat /= 0 ) GOTO 1
           !
        END IF
        !
        lii_ = lii - 1
        !
     END IF
     !
  END IF
  !
  ! ... only the first cpu initializes the file needed by parallelization
  ! ... among images
  !
  IF ( meta_ionode ) CALL new_image_init( nimage, fii_, tmp_dir_saved )
  !
  image = fii_ + my_image_id
  !
  scf_loop: DO
     !
     ! ... exit if available images are finished
     !
     IF ( image > lii_ ) EXIT scf_loop
     !
     pending_image = image
     !
     IF (dft_is_hybrid()) call stop_exx() 
     CALL do_scf( image, istat )
     !
     IF ( istat /= 0 ) GOTO 1
     !
     ! ... the new image is obtained (by ionode only)
     !
     CALL get_new_image( nimage, image, tmp_dir_saved )
     !
     CALL mp_bcast( image, ionode_id, intra_image_comm )
     !
  END DO scf_loop
  !
  ! ... after the first call to compute_scf the input values of startingpot
  ! ... and startingwfc are both set to 'file'
  !
  starting_pot = 'file'
  starting_wfc = 'file'
  !
  ! ... finalization of the job (this point is also reached in case of error
  ! ... condition)
  !
1 CALL mp_barrier( world_comm )
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... pes and grad_pes are communicated among "image" pools
     !
     CALL mp_sum( pes(fii:lii),        inter_image_comm )
     CALL mp_sum( grad_pes(:,fii:lii), inter_image_comm )
     IF ( lfcpopt ) CALL mp_sum( fcp_neb_ef(fii:lii), inter_image_comm )
     CALL mp_sum( istat,               inter_image_comm )
     !
  END IF
  !
  ! ... global status is computed here
  !
  IF ( istat == 0 ) THEN
     !
     stat = .TRUE.
     !
     pending_image = 0
     !
  ELSE
     !
     stat = .FALSE.
     !
     IF ( nimage > 1 ) THEN
        !
        CALL mp_min( pending_image, inter_image_comm )
        !
        IF ( meta_ionode ) CALL delete_if_present( exit_file )
        !
     END IF
     !
     IF ( meta_ionode ) THEN
        !
        ! ... some image didn't converge:  extrapolation is no longer
        ! ...                              possible, files are removed
        !
        WRITE( UNIT = iunpath, &
               FMT = '(/,5X,"cleaning-up extrapolation files"/)' )
        !
        DO image = pending_image, lii
           !
           tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                     TRIM( int_to_char( image ) ) // "/"
           !
           CALL delete_if_present( TRIM( tmp_dir ) // &
                                   TRIM( prefix ) // '.update' )
           !
        END DO
        !
     END IF
     !
  END IF
  !
  tmp_dir = tmp_dir_saved
  !
  RETURN
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE do_scf( image, istat )
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : diago_thr_init
      USE control_flags,    ONLY : ethr
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)    :: image
      INTEGER, INTENT(INOUT) :: istat
      !
      REAL(DP), EXTERNAL :: get_clock
      !
      ! ... self-consistency ( for non-frozen images only )
      !
      IF ( frozen(image) ) RETURN
      !
      CALL clean_pw( .FALSE. )
      !
      tcpu = get_clock( 'NEB' )
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
      tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                TRIM( int_to_char( image ) ) // "/"
      wfc_dir = tmp_dir
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
      ! ... tau is in alat units ( pos is in bohr )
      !
      tau = RESHAPE( pos(:,image), SHAPE( tau ) ) / alat
      !
      WRITE( stdout, '(/,5X,"coordinates at iteration ",I3,/)' ) istep_path
      !
      CALL output_tau( .FALSE., .FALSE. )
      !
      ! ... initialization of the scf calculation
      !
      CALL start_clock('PWSCF')
      CALL setup ()
      !
      IF ( lfcpopt ) THEN
         nelec = fcp_neb_nelec(image)
         tot_charge = SUM( zv(ityp(1:nat)) ) - nelec
      ENDIF
      !
      CALL init_run()
      !
      ! ... charge density and wavefunction extrapolation
      !
      CALL update_neb ()
      !
      ! ... self-consistency loop
      !
      CALL electrons()
      !
      CALL punch( 'all' )
      !
      ! ... scf convergence is checked here
      !
      IF ( .NOT.conv_elec ) THEN
         !
         istat = 1
         !
         WRITE( UNIT = iunpath, &
                FMT = '(/,5X,"WARNING :  scf convergence ", &
                       &     "NOT achieved on image ",I3)' ) image
         !
         ! ... in case of parallelization on images a stop signal
         ! ... is sent via the "EXIT" file
         !
         IF ( nimage > 1 ) CALL stop_other_images()
         !
         RETURN
         !
      END IF
      !
      ! ... self-consistent forces
      !
      CALL forces()
      !
      ! ... energy is converted from rydberg to hartree
      !
      pes(image) = etot / e2
      !
      ! ... gradients are converted from rydberg/bohr to hartree/bohr
      !
      grad_pes(:,image) = - RESHAPE( force, (/ dim1 /) ) / e2
      !
      ethr = diago_thr_init
      !
      CALL close_files(.FALSE.)
      !
      IF ( lfcpopt ) fcp_neb_ef(image) = ef
      !
      RETURN
      !
    END SUBROUTINE do_scf
    !
END SUBROUTINE compute_scf
