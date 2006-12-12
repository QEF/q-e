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
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : startingwfc, startingpot, diago_thr_init
  USE constants,        ONLY : e2
  USE control_flags,    ONLY : conv_elec, istep, history, &
                               alpha0, beta0, ethr, pot_order
  USE check_stop,       ONLY : check_stop_now
  USE vlocal,           ONLY : strf
  USE cell_base,        ONLY : bg, alat
  USE gvect,            ONLY : ngm, g, nr1, nr2, nr3, eigts1, eigts2, eigts3
  USE ions_base,        ONLY : tau, nat, nsp, ityp
  USE basis,            ONLY : startingwfc_ => startingwfc, &
                               startingpot_ => startingpot
  USE ener,             ONLY : etot
  USE force_mod,        ONLY : force
  USE io_files,         ONLY : prefix, tmp_dir, iunpath, iunupdate, &
                               exit_file, iunexit, delete_if_present
  USE path_formats,     ONLY : scf_fmt, scf_fmt_para
  USE path_variables,   ONLY : pos, pes, grad_pes, dim, pending_image, &
                               istep_path, frozen, write_save, num_of_images
  USE io_global,        ONLY : stdout, ionode, ionode_id, meta_ionode
  USE mp_global,        ONLY : inter_image_comm, intra_image_comm, &
                               my_image_id, nimage, root_image
  USE mp,               ONLY : mp_bcast, mp_barrier, mp_sum, mp_min
  USE path_io_routines, ONLY : new_image_init, get_new_image, &
                               stop_other_images
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: fii, lii         ! indeces of first and last images
  LOGICAL, INTENT(OUT) :: stat
  !
  INTEGER               :: fii_, lii_      ! local copies of ini and inl
  INTEGER               :: image, istat
  REAL(DP)              :: tcpu 
  CHARACTER (LEN=256)   :: tmp_dir_saved
  LOGICAL               :: file_exists, opnd
  REAL(DP), ALLOCATABLE :: tauold(:,:,:)
   ! previous positions of atoms (needed for extrapolation)
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  REAL(DP),         EXTERNAL :: get_clock
  !
  !
  fii_ = fii
  lii_ = lii
  !
  istep = istep_path
  istat = 0
  !
  CALL flush_unit( iunpath )
  !
  ALLOCATE( tauold( 3, nat, 3 ) )
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... vectors pes and grad_pes are initalized to zero for all images on
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
     ! ... self-consistency on the first and last images is done separately
     !
     IF ( fii == 1 ) THEN
        !
        IF ( my_image_id == root_image ) THEN
           !
           CALL do_scf( 1, istat )
           !
           IF ( interrupt_run( istat ) ) RETURN
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
           CALL do_scf( num_of_images, istat )
           !
           IF ( interrupt_run( istat ) ) RETURN
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
  IF ( meta_ionode ) CALL new_image_init( fii_, tmp_dir_saved )
  !
  image = fii_ + my_image_id
  !
  ! ... all processes are syncronized (needed to have an ordered output)
  !
  CALL mp_barrier()
  !
  scf_loop: DO
     !
     ! ... exit if available images are finished
     !
     IF ( image > lii_ ) EXIT scf_loop
     !
     pending_image = image
     !
     IF ( check_stop_now( iunpath ) ) THEN
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
     CALL do_scf( image, istat )
     !
     IF ( interrupt_run( istat ) ) RETURN
     !
     ! ... the new image is obtained (by ionode only)
     !
     CALL get_new_image( image, tmp_dir_saved )
     !
     CALL mp_bcast( image, ionode_id, intra_image_comm )
     !
  END DO scf_loop
  !
  tmp_dir = tmp_dir_saved
  !
  ! ... after the first call to compute_scf the input values of startingpot
  ! ... and startingwfc are both set to 'file'
  !
  startingpot = 'file'
  startingwfc = 'file'
  startingpot_ = startingpot
  startingwfc_ = startingwfc
  !
  DEALLOCATE( tauold )
  !
  CALL mp_barrier()
  !
  IF ( nimage > 1 ) THEN
     !
     ! ... pes and grad_pes are communicated among "image" pools
     !
     CALL mp_sum( pes(fii:lii),        inter_image_comm )
     CALL mp_sum( grad_pes(:,fii:lii), inter_image_comm )
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
  END IF
  !
  RETURN
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE do_scf( image, istat )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)    :: image
      INTEGER, INTENT(INOUT) :: istat
      !
      ! ... self-consistency ( for non-frozen images only )
      !
      IF ( frozen(image) ) RETURN
      !
      CALL clean_pw( .FALSE. )
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
      tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                TRIM( int_to_char( image ) ) // "/"
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
      WRITE( stdout, '(/,5X,"coordinates at iteration ",I3,/)' ) istep
      !
      CALL output_tau( .FALSE. )
      !
      ! ... initialization of the scf calculation
      !
      CALL init_run()
      !
      IF ( ionode ) THEN
         !
         ! ... the file containing old positions is opened 
         ! ... ( needed for extrapolation )
         !
         CALL seqopn( iunupdate, 'update', 'FORMATTED', file_exists )
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
      CALL mp_bcast( history, ionode_id, intra_image_comm )
      CALL mp_bcast( tauold,  ionode_id, intra_image_comm )
      !
      IF ( conv_elec .AND. history > 0 ) THEN
         !
         ! ... potential and wavefunctions are extrapolated only if
         ! ... we are starting a new self-consistency ( scf on the
         ! ... previous image was achieved )
         !
         IF ( ionode ) THEN 
            !
            ! ... find the best coefficients for the extrapolation of 
            ! ... the potential
            !
            CALL find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )
            !
         END IF
         !
         CALL mp_bcast( alpha0, ionode_id, intra_image_comm )
         CALL mp_bcast( beta0,  ionode_id, intra_image_comm )
         !
         IF ( pot_order > 0 ) THEN
            !
            ! ... structure factors of the old positions are computed 
            ! ... (needed for the old atomic charge)
            !
            CALL struc_fact( nat, tauold(:,:,1), nsp, ityp, ngm, g, bg, &
                             nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
            !
         END IF
         !
         CALL update_pot()
         !
      END IF
      !
      ! ... self-consistency loop
      !
      CALL electrons()
      !
      ! ... scf convergence is checked
      !
      IF ( .NOT.conv_elec ) THEN
         !
         istat = 1
         !
         WRITE( UNIT = iunpath, &
                FMT = '(/,5X,"WARNING :  scf convergence NOT achieved",/)' )
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
      grad_pes(:,image) = - RESHAPE( force, (/ dim /) ) / e2
      !
      IF ( ionode ) THEN
         !
         ! ... save the previous two steps
         ! ... ( a total of three ionic steps is saved )
         !
         tauold(:,:,3) = tauold(:,:,2)
         tauold(:,:,2) = tauold(:,:,1)
         tauold(:,:,1) = tau(:,:)
         !
         history = MIN( 3, ( history + 1 ) )
         !
         CALL seqopn( iunupdate, 'update', 'FORMATTED', file_exists ) 
         !
         WRITE( UNIT = iunupdate, FMT = * ) history
         WRITE( UNIT = iunupdate, FMT = * ) tauold
         !
         CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
         !
      END IF
      !
      ! ... the save file is written ( if required )
      !
      IF ( write_save ) CALL punch( 'all' )
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
      RETURN
      !
    END SUBROUTINE do_scf
    !
    !-----------------------------------------------------------------------
    FUNCTION interrupt_run( istat )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(INOUT) :: istat
      LOGICAL                :: interrupt_run
      !
      interrupt_run = ( istat == 1 )
      !
      IF ( istat == 0 ) RETURN
      !
      DEALLOCATE( tauold )
      !
      CALL mp_barrier()
      !
      IF ( nimage > 1 ) THEN
         !
         ! ... pes and grad_pes are communicated among "image" pools
         !
         CALL mp_sum( pes(fii:lii),        inter_image_comm )
         CALL mp_sum( grad_pes(:,fii:lii), inter_image_comm )
         CALL mp_sum( istat,               inter_image_comm )
         !
      END IF
      !
      ! ... global status is computed here
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
      RETURN
      !
    END FUNCTION interrupt_run
    !
END SUBROUTINE compute_scf
