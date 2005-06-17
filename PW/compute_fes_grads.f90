!
! Copyright (C) 2005 PWSCF-FPMD-CPV groups
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
  USE coarsegrained_vars, ONLY : new_target, to_target, dfe_acc, &
                                 max_shake_iter, max_fe_iter, num_acc, &
                                 fe_grad_thr
  USE path_variables,     ONLY : pos, pes, grad_pes, frozen, &
                                 num_of_images, istep_path, suspended_image
  USE constraints_module, ONLY : lagrange, target, init_constraint, &
                                 check_constrain, deallocate_constraint
  USE dynam,              ONLY : dt
  USE control_flags,      ONLY : conv_elec, istep, history, wg_set, &
                                 alpha0, beta0, ethr, pot_order,    &
                                 conv_ions, ldamped
  USE cell_base,          ONLY : alat, at, bg
  USE gvect,              ONLY : ngm, g, nr1, nr2, nr3, eigts1, eigts2, eigts3
  USE vlocal,             ONLY : strf
  USE ener,               ONLY : etot
  USE ions_base,          ONLY : nat, nsp, tau, ityp, if_pos
  USE path_formats,       ONLY : scf_fmt, scf_fmt_para
  USE io_files,           ONLY : prefix, tmp_dir, iunpath, iunaxsf, &
                                 iunupdate, exit_file, iunexit
  USE parser,             ONLY : int_to_char, delete_if_present
  USE constants,          ONLY : bohr_radius_angs
  USE io_global,          ONLY : stdout, ionode, ionode_id, meta_ionode
  USE mp_global,          ONLY : inter_image_comm, intra_image_comm, &
                                 my_image_id, nimage, root
  USE mp,                 ONLY : mp_bcast, mp_barrier, mp_sum, mp_min
  USE check_stop,         ONLY : check_stop_now
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: N_in, N_fin
  LOGICAL, INTENT(OUT) :: stat
  INTEGER              :: image, iter, counter
  REAL (KIND=DP)       :: tcpu, error
  CHARACTER (LEN=256)  :: tmp_dir_saved, filename
  LOGICAL              :: opnd, file_exists, ldamped_saved
  LOGICAL              :: tfirst
  REAL(KIND=DP), ALLOCATABLE :: tauold(:,:,:)
    ! previous positions of atoms (needed for extrapolation)
  !
  REAL (KIND=DP), EXTERNAL :: get_clock
  !
  !
  CALL flush_unit( iunpath )
  !
  ALLOCATE( tauold( 3, nat, 3 ) )
  !
  OPEN( UNIT = iunaxsf, FILE = TRIM( prefix ) // "_" // &
      & TRIM( int_to_char( istep_path ) ) // ".axsf", &
        STATUS = "UNKNOWN", ACTION = "WRITE" )
  !
  WRITE( UNIT = iunaxsf, FMT = '(" ANIMSTEPS ",I3)' ) num_of_images
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
  tmp_dir_saved = tmp_dir
  ldamped_saved = ldamped
  !
  ! ... vectors pes and grad_pes are initalized to zero for all images on
  ! ... all nodes: this is needed for the final mp_sum()
  !
  IF ( my_image_id == root ) THEN
     !
     FORALL( image = N_in : N_fin, .NOT. frozen(image)   )
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
  ! ... only the tfirst cpu initializes the file needed by parallelization 
  ! ... among images
  !
  IF ( meta_ionode ) CALL new_image_init()
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
     ! ... self-consistency ( for non-frozen images only )
     !
     IF ( .NOT. frozen(image) ) THEN
        !
        CALL clean_pw( .FALSE. )
        !
        CALL deallocate_constraint()
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
           dfe_acc(:,:) = 0.D0
           !
           OPEN( UNIT = 1000, FILE = filename )
           !
           READ( 1000, * ) tau
           READ( 1000, * ) dfe_acc(:,2:num_acc)
           READ( 1000, * ) counter
           !
           CLOSE( UNIT = 1000 )
           !
           counter = MIN( counter + 1, num_acc )
           !
        ELSE
           !
           dfe_acc(:,:) = 0.D0
           !
           counter = 1
           !
        END IF
        !
        ! ... the new value of the order-parameter is set here
        !
        CALL init_constraint( nat, tau, alat, ityp, if_pos )
        !
        new_target(:) = pos(:,image)
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
           ! ... we are starting a new self-consistency (scf on the 
           ! ... previous image was achieved)
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
        wg_set = .FALSE.
        tfirst = .TRUE.
        !
        ! ... tfirst the system is "adiabatically" moved to the new target
        ! ... This is always done with standard MD
        !
        ldamped = .FALSE.
        !
        to_target(:) = new_target(:) - target(:)
        !
        CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
        !
        DO iter = 1, max_shake_iter
           !
           istep = iter
           !
           CALL electronic_scf( tfirst, stat )
           !
           IF ( .NOT. stat ) RETURN
           !
           CALL move_ions()
           !
           target(:) = target(:) + to_target(:) / DBLE( max_shake_iter )
           !
           tfirst = .FALSE.
           !
        END DO
        !
        ldamped = ldamped_saved
        !
        ! ... then the free energy gradients are computed
        !
        CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
        !
        DO iter = 1, max_fe_iter
           !
           istep = iter
           !
           CALL electronic_scf( tfirst, stat )
           !
           IF ( .NOT. stat )  RETURN
           !
           CALL move_ions()
           !
           IF ( ldamped ) THEN
              !
              ! ... zero temperature
              !
              IF ( conv_ions ) EXIT
              !
           ELSE
              !
              ! ... finite temperature
              !
              dfe_acc(:,1) = dfe_acc(:,1) - lagrange(:)
              !
           END IF
           !
           tfirst = .FALSE.
           !
        END DO
        !
        ! ... the averages are computed here
        !
        IF ( ldamped ) THEN
           !
           ! ... zero temperature
           !
           grad_pes(:,image) = - lagrange(:) / e2
           !
           pes(image) = etot
           !
        ELSE
           !
           ! ... finite temperature
           !
           dfe_acc(:,1) = dfe_acc(:,1) / DBLE( istep )
           !
           grad_pes(:,image) = 0.D0
           !
           DO iter = 1, counter
              !
              grad_pes(:,image) = grad_pes(:,image) + dfe_acc(:,iter)
              !
           END DO
           !
           grad_pes(:,image) = grad_pes(:,image) / DBLE( counter ) / e2
           !
        END IF
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
        CALL write_config( image )
        !
        ! ... the restart file is written here
        !
        OPEN( UNIT = 1000, FILE = filename )
        !
        WRITE( 1000, * ) tau
        WRITE( 1000, * ) dfe_acc(:,1:num_acc-1)
        WRITE( 1000, * ) counter
        !
        CLOSE( UNIT = 1000 )
        !
     END IF
     !
     ! ... the new image is obtained (by ionode only)
     !
     CALL get_new_image( image )
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
  DEALLOCATE( tauold )
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
  ! ... after the tfirst call to compute_scf the input values of startingpot
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
       USE io_files,       ONLY : iunnewimage
       USE path_variables, ONLY : tune_load_balance
       !
       IMPLICIT NONE       
       !
       IF ( nimage == 1 .OR. .NOT. tune_load_balance ) RETURN
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
       USE io_files,       ONLY : iunnewimage, iunblock
       USE io_global,      ONLY : ionode
       USE path_variables, ONLY : tune_load_balance
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(INOUT) :: image
       INTEGER                :: ioerr
       CHARACTER (LEN=256)    :: filename
       LOGICAL                :: opened, exists
       !
       !
       IF ( .NOT. ionode ) RETURN
       !
       IF ( nimage > 1 ) THEN
          !
          IF ( tune_load_balance ) THEN
             !
             filename = TRIM( tmp_dir_saved ) // TRIM( prefix ) // '.BLOCK'
             !
             open_loop: DO
                !          
                OPEN( UNIT = iunblock, FILE = TRIM( filename ), &
                    & IOSTAT = ioerr, STATUS = 'NEW' )
                !
                IF ( ioerr > 0 ) CYCLE open_loop
                !
                INQUIRE( UNIT = iunnewimage, OPENED = opened )
                !
                IF ( .NOT. opened ) THEN
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
             END DO open_loop
             !
             CLOSE( UNIT = iunblock, STATUS = 'DELETE' )
             !
          ELSE
             !
             image = image + nimage
             !
          END IF
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
       USE io_global, ONLY : ionode
       !
       IMPLICIT NONE
       !
       !
       IF ( .NOT. ionode ) RETURN
       !
       OPEN( UNIT = iunexit, FILE = TRIM( exit_file ) )
       CLOSE( UNIT = iunexit, STATUS = 'KEEP' )               
       !
       RETURN       
       !
     END SUBROUTINE stop_other_images
     !
END SUBROUTINE compute_fes_grads
!
!----------------------------------------------------------------------------
SUBROUTINE write_config( iter )
  !----------------------------------------------------------------------------
  !
  USE input_parameters, ONLY : atom_label
  USE io_files,         ONLY : iunaxsf
  USE constants,        ONLY : bohr_radius_angs
  USE ions_base,        ONLY : nat, tau, ityp
  USE cell_base,        ONLY : alat
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  INTEGER             :: atom
  !
  !
  WRITE( UNIT = iunaxsf, FMT = '(" PRIMCOORD ",I3)' ) iter
  WRITE( UNIT = iunaxsf, FMT = '(I5,"  1")' ) nat
  !
  DO atom = 1, nat
     !
     WRITE( UNIT = iunaxsf, FMT = '(A2,3(2X,F18.10))' ) &
            TRIM( atom_label(ityp(atom)) ), &
             tau(1,atom) * alat * bohr_radius_angs, &
             tau(2,atom) * alat * bohr_radius_angs, &
             tau(3,atom) * alat * bohr_radius_angs
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE write_config
!
!----------------------------------------------------------------------------
SUBROUTINE electronic_scf( tfirst, stat )
  !----------------------------------------------------------------------------
  !
  USE control_flags, ONLY : conv_elec
  USE io_files,      ONLY : iunpath
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN)  :: tfirst
  LOGICAL, INTENT(OUT) :: stat
  !
  !
  IF ( .NOT. tfirst ) CALL hinit1()
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
