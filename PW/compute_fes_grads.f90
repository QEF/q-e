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
  USE coarsegrained_vars, ONLY : new_target, to_target, dfe_acc, &
                                 max_shake_iter, max_fe_iter, num_acc, &
                                 fe_grad_thr
  USE path_variables,     ONLY : pos, pes, grad_pes, frozen, &
                                 num_of_images, istep_path
  USE constraints_module, ONLY : lagrange, target, init_constraint, &
                                 check_constrain, deallocate_constraint
  USE dynam,              ONLY : dt
  USE control_flags,      ONLY : istep, conv_ions, ldamped
  USE cell_base,          ONLY : alat, at
  USE ener,               ONLY : etot
  USE ions_base,          ONLY : nat, tau, ityp, if_pos
  USE path_formats,       ONLY : scf_fmt
  USE io_files,           ONLY : prefix, tmp_dir, iunpath, iunaxsf
  USE parser,             ONLY : int_to_char, delete_if_present
  USE constants,          ONLY : bohr_radius_angs
  USE io_global,          ONLY : stdout, ionode, ionode_id, meta_ionode
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: N_in, N_fin
  LOGICAL, INTENT(OUT) :: stat
  INTEGER              :: image, iter, counter
  REAL (KIND=DP)       :: tcpu, error
  CHARACTER (LEN=256)  :: tmp_dir_saved, filename
  LOGICAL              :: opnd, exists, ldamped_saved
  LOGICAL              :: first
  !
  REAL (KIND=DP), EXTERNAL :: get_clock
  !
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
  !
  ldamped_saved = ldamped
  !
  image = N_in
  !
  fes_loop: DO
     !
     ! ... exit if available images are finished
     !
     IF ( image > N_fin ) EXIT fes_loop
     !
     IF ( frozen(image) ) THEN
        !
        ! ... the forces are needed only for non-frozen images 
        !
        image = image + 1
        !
        CYCLE fes_loop
        !
     END IF
     !
     tcpu = get_clock( 'PWSCF' )
     !
     CALL clean_pw( .FALSE. )
     !
     CALL deallocate_constraint()
     !
     WRITE( UNIT = iunpath, FMT = scf_fmt ) tcpu, image     
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
     INQUIRE( FILE = filename, EXIST = exists )
     !
     IF ( exists ) THEN
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
     CALL write_config( image )
     !
     ! ... first the system is "adiabatically" moved to the new target
     !
     to_target(:) = new_target(:) - target(:)
     !
     CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) // '.md' )
     !
     ldamped = .FALSE.
     !
     DO iter = 1, max_shake_iter
        !
        istep = iter
        !
        CALL electronic_scf( first )
        !
        CALL move_ions()
        !
        target(:) = target(:) + to_target(:) / DBLE( max_shake_iter )
        !
        first = .FALSE.
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
        CALL electronic_scf( first )
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
        first = .FALSE.
        !
     END DO
     !
     ! ... the averages are computed here
     !
     IF ( ldamped ) THEN
        !
        ! ... zero temperature
        !
        grad_pes(:,image) = - lagrange(:)
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
        grad_pes(:,image) = grad_pes(:,image) / DBLE( counter )
        !
     END IF
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
     ! ... the new image
     !
     image = image + 1
     !
  END DO fes_loop
  !
  CLOSE( UNIT = iunaxsf )
  !
  tmp_dir = tmp_dir_saved
  !
  stat = .TRUE.
  !
END SUBROUTINE compute_fes_grads
!
!------------------------------------------------------------------------
SUBROUTINE write_config( iter )
  !------------------------------------------------------------------------
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
!------------------------------------------------------------------------
SUBROUTINE electronic_scf( first )
  !------------------------------------------------------------------------
  !
  USE control_flags, ONLY : conv_elec
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: first 
  !
  !
  IF ( .NOT. first ) CALL hinit1()
  !
  CALL electrons()
  !
  IF ( .NOT. conv_elec ) CALL stop_pw( conv_elec )
  !
  CALL forces()
  !
  RETURN
  !
END SUBROUTINE electronic_scf
