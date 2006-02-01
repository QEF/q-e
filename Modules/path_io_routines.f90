!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
MODULE path_io_routines
  !----------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines used for I/O in path
  ! ... optimisations
  !
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : au, bohr_radius_angs
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: io_path_start, io_path_stop
  PUBLIC :: read_restart
  PUBLIC :: write_restart, write_dat_files, write_output
  PUBLIC :: new_image_init, get_new_image, stop_other_images
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE io_path_start()
       !-----------------------------------------------------------------------
       !
       USE io_global, ONLY : stdout
       USE io_global, ONLY : ionode, ionode_id
       USE mp_global, ONLY : me_image, root_image
       !
       IMPLICIT NONE
       !
       !
       ! ... the I/O node is set again according to the number of parallel
       ! ... images that have been required: for each parallel image there
       ! ... is only one node that does I/O
       !
       IF ( me_image == root_image ) THEN
          !
          ionode = .TRUE.
          !
       ELSE
          !
          ionode = .FALSE.
          !
       END IF
       !
       ionode_id = root_image
       !
       ! ... stdout is connected to a file ( different for each image ) 
       ! ... via unit 17 ( only root_image performes I/O )
       !
       IF ( me_image == root_image ) stdout = 17
       !
       RETURN
       !
     END SUBROUTINE io_path_start
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE io_path_stop()
       !-----------------------------------------------------------------------
       !
       USE io_global,  ONLY : stdout, io_global_start
       USE mp_global,  ONLY : mpime, root
       !
       IMPLICIT NONE
       !
       !
       ! ... the original I/O node is set again 
       !
       CALL io_global_start( mpime, root )
       !
       ! ... stdout is reconnected to standard output
       !
       stdout = 6
       !
       RETURN
       !
     END SUBROUTINE io_path_stop
     !
     !-----------------------------------------------------------------------
     SUBROUTINE read_restart()
       !-----------------------------------------------------------------------
       !
       USE control_flags,    ONLY : istep, nstep, conv_elec, &
                                    lneb, lsmd, lcoarsegrained
       USE io_files,         ONLY : iunpath, iunrestart, path_file   
       USE input_parameters, ONLY : if_pos
       USE path_variables,   ONLY : istep_path, nstep_path, suspended_image, &
                                    dim, num_of_images, pos, pes, grad_pes,  &
                                    frozen , lquick_min, lmol_dyn,           &
                                    ldamped_dyn, first_last_opt
       USE path_variables,   ONLY : vel, Emax, Emin, Emax_index
       USE path_variables,   ONLY : Nft, ft_coeff, num_of_modes
       USE io_global,        ONLY : meta_ionode, meta_ionode_id
       USE mp,               ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !    
       INTEGER              :: i, j, ia, ierr
       CHARACTER (LEN=256)  :: input_line
       REAL (DP)       :: val
       LOGICAL              :: exists
       LOGICAL, EXTERNAL    :: matches
       !
       !
       IF ( meta_ionode ) THEN
          !
          WRITE( UNIT = iunpath, &
                 FMT = '(/,5X,"reading file ", A,/)') TRIM( path_file )
          !
          INQUIRE( FILE = TRIM( path_file ), EXIST = exists )
          !
          IF ( .NOT. exists ) &
             CALL errore( 'read_restart()', 'restart file not found', 1 )
          !
          OPEN( UNIT = iunrestart, FILE = path_file, STATUS = "OLD", &
                ACTION = "READ" )
          !
          READ( UNIT = iunrestart, FMT = '(256A)' ) input_line
          !
          IF ( matches( "RESTART INFORMATION", input_line ) ) THEN
             !
             READ( UNIT = iunrestart, FMT = * ) istep_path
             READ( UNIT = iunrestart, FMT = * ) nstep_path
             READ( UNIT = iunrestart, FMT = * ) suspended_image
             READ( UNIT = iunrestart, FMT = * ) conv_elec
             !
          ELSE   
             !
             ! ... mandatory fields
             !
             CALL errore( 'read_restart()', 'RESTART INFORMATION missing', 1 )
             !
          END IF
          !
          READ( UNIT = iunrestart, FMT = '(256A)' ) input_line
          !
          IF( matches( "NUMBER OF IMAGES", input_line ) ) THEN
             !
             ! ... optional field
             !
             READ( UNIT = iunrestart, FMT = * ) num_of_images
             !
             IF ( lsmd ) THEN
                !
                ! ... fourier dimensions updated
                !
                Nft = ( num_of_images - 1 )
                !
                num_of_modes = ( Nft - 1 )
                !
                ft_coeff = 2.D0 / DBLE( Nft )
                !
             END IF
             !
             READ( UNIT = iunrestart, FMT = '(256A)' ) input_line
             !
          END IF
          !
          IF ( .NOT. ( matches( "ENERGIES, POSITIONS AND GRADIENTS", &
                                 input_line ) ) ) THEN
             !
             ! ... mandatory fields
             !
             CALL errore( 'read_restart()', &
                          'ENERGIES, POSITIONS AND GRADIENTS missing', 1 )
             !
          END IF
          !
          IF ( lcoarsegrained ) THEN
             !
             DO i = 1, num_of_images
                !
                READ( UNIT = iunrestart, FMT = * )
                !
                DO j = 1, dim
                   !
                   READ( UNIT = iunrestart, FMT = * ) pos(j,i), grad_pes(j,i)
                   !
                END DO
                !
             END DO
             !
          ELSE
             !
             READ( UNIT = iunrestart, FMT = * )
             READ( UNIT = iunrestart, FMT = * ) pes(1)
             !
             ia = 0  
             !
             if_pos = 0
             !
             DO j = 1, dim, 3 
                !
                ia = ia + 1
                !
                READ( UNIT = iunrestart, FMT = * ) &
                    pos(j+0,1),                    &
                    pos(j+1,1),                    &
                    pos(j+2,1),                    &
                    grad_pes(j+0,1),               &
                    grad_pes(j+1,1),               &
                    grad_pes(j+2,1),               &
                    if_pos(1,ia),                  &
                    if_pos(2,ia),                  &
                    if_pos(3,ia) 
                !
                grad_pes(:,1) = grad_pes(:,1) * &
                               DBLE( RESHAPE( if_pos, (/ dim /) ) )
                !
             END DO
             !
             DO i = 2, num_of_images
                !
                READ( UNIT = iunrestart, FMT = * )
                READ( UNIT = iunrestart, FMT = * ) pes(i)
                !
                DO j = 1, dim, 3 
                   !
                   READ( UNIT = iunrestart, FMT = * ) &
                       pos(j+0,i),                    &
                       pos(j+1,i),                    &
                       pos(j+2,i),                    &
                       grad_pes(j+0,i),               &
                       grad_pes(j+1,i),               &
                       grad_pes(j+2,i)
                   !
                END DO
                !
                grad_pes(:,i) = grad_pes(:,i) * &
                                DBLE( RESHAPE( if_pos, (/ dim /) ) )
                !
             END DO
             !
             IF ( suspended_image == 0 ) THEN
                !
                Emin       = MINVAL( pes(:) )
                Emax       = MAXVAL( pes(:) )
                Emax_index = MAXLOC( pes(:), 1 )
                !
             END IF
             !
          END IF
          !
          READ( UNIT = iunrestart, FMT = '(256A)', IOSTAT = ierr ) input_line
          !
          IF ( ( ierr == 0 ) .AND. &
               ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) ) THEN
             !
             IF ( matches( "QUICK-MIN FIELDS", input_line ) ) THEN
                !
                ! ... optional fields
                !
                IF ( lcoarsegrained ) THEN
                   !
                   DO i = 1, num_of_images
                      !
                      READ( UNIT = iunrestart, FMT = * )
                      READ( UNIT = iunrestart, FMT = * ) frozen(i)                      
                      !
                      DO j = 1, dim
                         !
                         READ( UNIT = iunrestart, FMT = * ) vel(j,i)
                         !
                      END DO
                      !
                   END DO                   
                   !
                ELSE
                   !
                   DO i = 1, num_of_images
                      !
                      READ( UNIT = iunrestart, FMT = * )
                      READ( UNIT = iunrestart, FMT = * ) frozen(i)
                      !
                      DO j = 1, dim, 3
                         !
                         READ( UNIT = iunrestart, FMT = * ) &
                             vel(j+0,i),                    & 
                             vel(j+1,i),                    &
                             vel(j+2,i)
                         !
                      END DO
                      !
                      vel(:,i) = vel(:,i) * &
                                 DBLE( RESHAPE( if_pos, (/ dim /) ) )
                      !
                   END DO
                   !
                END IF
                !
             END IF
             !
          END IF
          !
          CLOSE( iunrestart )
          !
       END IF
       !
       ! ... broadcast to all nodes
       !
       CALL mp_bcast( istep_path,      meta_ionode_id )
       CALL mp_bcast( nstep_path,      meta_ionode_id )
       CALL mp_bcast( suspended_image, meta_ionode_id )
       CALL mp_bcast( conv_elec,       meta_ionode_id )
       !
       CALL mp_bcast( pos,      meta_ionode_id )
       CALL mp_bcast( if_pos,   meta_ionode_id )
       CALL mp_bcast( pes,      meta_ionode_id )
       CALL mp_bcast( grad_pes, meta_ionode_id )
       !
       CALL mp_bcast( num_of_images, meta_ionode_id )
       !
       CALL mp_bcast( Emax,       meta_ionode_id )
       CALL mp_bcast( Emin,       meta_ionode_id )
       CALL mp_bcast( Emax_index, meta_ionode_id )
       !
       IF ( lsmd ) THEN
          !          
          CALL mp_bcast( num_of_modes, meta_ionode_id )
          CALL mp_bcast( Nft,          meta_ionode_id )
          CALL mp_bcast( ft_coeff,     meta_ionode_id )
          !
       END IF
       !
       IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
          !
          CALL mp_bcast( frozen, meta_ionode_id )
          CALL mp_bcast( vel,    meta_ionode_id )
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE read_restart
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_restart()
       !-----------------------------------------------------------------------
       !
       USE input_parameters, ONLY : if_pos       
       USE io_files,         ONLY : iunrestart, path_file, tmp_dir 
       USE control_flags,    ONLY : conv_elec, lneb, lsmd, lcoarsegrained
       USE path_variables,   ONLY : istep_path, nstep_path, suspended_image,  &
                                    dim, num_of_images, pos, pes, grad_pes,   &
                                    vel, frozen, Nft, num_of_modes, lmol_dyn, &
                                    lquick_min, ldamped_dyn, first_last_opt
       USE path_formats,     ONLY : energy, restart_first, restart_others, &
                                    quick_min
       USE io_global,        ONLY : meta_ionode
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER             :: i, j, ia
       CHARACTER (LEN=256) :: file
       CHARACTER(LEN=6), EXTERNAL :: int_to_char
       !
       IF ( meta_ionode ) THEN
          !
          ! ... first the restart file is written in the working directory
          !
          OPEN( UNIT = iunrestart, FILE = path_file, STATUS = "UNKNOWN", &
                ACTION = "WRITE" )
          !
          CALL write_common_fields( iunrestart )
          !
          IF (  lquick_min .OR. ldamped_dyn .OR. lmol_dyn  ) THEN
             !
             CALL write_quick_min_fields( iunrestart )
             ! 
          END IF
          !
          CLOSE( iunrestart )
          !
          ! ... then, if suspended_image == 0, it is also written on the 
          ! ... scratch direcoty (a backup copy at each iteration)
          !
          IF ( suspended_image == 0 ) THEN
             !
             file = TRIM( tmp_dir ) // &
                    TRIM( path_file ) // TRIM( int_to_char( istep_path ) )
             !
             OPEN( UNIT = iunrestart, FILE = TRIM( file ), &
                   STATUS = "UNKNOWN",  ACTION = "WRITE" )
             !
             CALL write_common_fields( iunrestart )
             !
             IF (  lquick_min .OR. ldamped_dyn .OR. lmol_dyn  ) THEN
                !
                CALL write_quick_min_fields( iunrestart )
                ! 
             END IF
             !
             CLOSE( iunrestart )          
             !
          END IF
          !
       END IF
       !
       CONTAINS
         !
         !-------------------------------------------------------------------
         SUBROUTINE write_common_fields( in_unit )
           !-------------------------------------------------------------------
           !
           IMPLICIT NONE
           !
           INTEGER, INTENT(IN) :: in_unit
           !
           !
           WRITE( UNIT = in_unit, FMT = '("RESTART INFORMATION")' )
           !
           WRITE( UNIT = in_unit, FMT = '(I4)' ) istep_path
           WRITE( UNIT = in_unit, FMT = '(I4)' ) nstep_path
           WRITE( UNIT = in_unit, FMT = '(I4)' ) suspended_image
           WRITE( UNIT = in_unit, FMT = '(L1)' ) conv_elec
           !
           WRITE( UNIT = in_unit, FMT = '("NUMBER OF IMAGES")' )
           !
           WRITE( UNIT = in_unit, FMT = '(I4)' ) num_of_images
           !
           WRITE( UNIT = in_unit, &
                  FMT = '("ENERGIES, POSITIONS AND GRADIENTS")' )
           !
           DO i = 1, num_of_images
              !
              WRITE( UNIT = in_unit, FMT = '("Image: ",I4)' ) i
              !
              IF ( lcoarsegrained ) THEN
                 !
                 DO j = 1, dim
                    !
                    WRITE( UNIT = in_unit, &
                           FMT = '(2(2X,F18.12))' ) pos(j,i), grad_pes(j,i)
                    !
                 END DO
                 !
              ELSE
                 !
                 WRITE( UNIT = in_unit, FMT = energy ) pes(i)
                 !
                 ia = 0
                 !
                 DO j = 1, dim, 3
                    !
                    ia = ia + 1
                    !
                    IF ( i == 1 ) THEN
                       !
                       WRITE( UNIT = in_unit, FMT = restart_first ) &
                           pos(j+0,i),                              &
                           pos(j+1,i),                              &
                           pos(j+2,i),                              &
                           grad_pes(j+0,i),                         &
                           grad_pes(j+1,i),                         &
                           grad_pes(j+2,i),                         &
                           if_pos(1,ia),                            &
                           if_pos(2,ia),                            &
                           if_pos(3,ia) 
                       !
                    ELSE
                       !
                       WRITE( UNIT = in_unit, FMT = restart_others ) &
                           pos(j+0,i),                               &
                           pos(j+1,i),                               &
                           pos(j+2,i),                               &
                           grad_pes(j+0,i),                          &
                           grad_pes(j+1,i),                          & 
                           grad_pes(j+2,i)
                       !
                    END IF
                    !
                 END DO
                 !
              END IF
              !
           END DO
           !
           RETURN
           !
         END SUBROUTINE write_common_fields
         !
         !-------------------------------------------------------------------
         SUBROUTINE write_quick_min_fields( in_unit )
           !-------------------------------------------------------------------
           !
           IMPLICIT NONE
           !
           INTEGER, INTENT(IN) :: in_unit
           !
           !
           WRITE( UNIT = in_unit, FMT = '("QUICK-MIN FIELDS")' )
           !
           DO i = 1, num_of_images
              !
              WRITE( UNIT = in_unit, FMT = '("Image: ",I4)' ) i
              WRITE( UNIT = in_unit, &
                     FMT = '(2(L1,1X))' ) frozen(i)
              !
              IF ( lcoarsegrained ) THEN
                 !
                 DO j = 1, dim
                    !
                    WRITE( UNIT = in_unit, &
                           FMT = '(2X,F18.12)' ) vel(j,i)
                    !
                 END DO
                 !
              ELSE
                 !
                 DO j = 1, dim, 3
                    !
                    WRITE( UNIT = in_unit, FMT = quick_min ) &
                        vel(j+0,i),                          & 
                        vel(j+1,i),                          &
                        vel(j+2,i)
                    !
                 END DO
                 !
              END IF
              !
           END DO
           !
           RETURN
           !
         END SUBROUTINE write_quick_min_fields
         !
     END SUBROUTINE write_restart
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_dat_files()
       !-----------------------------------------------------------------------
       !
       USE constants,        ONLY : pi
       USE input_parameters, ONLY : atom_label
       USE control_flags,    ONLY : lcoarsegrained
       USE cell_base,        ONLY : alat, at
       USE ions_base,        ONLY : ityp, nat
       USE path_formats,     ONLY : dat_fmt, int_fmt, xyz_fmt, axsf_fmt
       USE path_variables,   ONLY : pos, grad_pes, pes, num_of_images, &
                                    path_length, react_coord, tangent, &
                                    dim, Emax_index, error
       USE io_files,         ONLY : iundat, iunint, iunxyz, iunaxsf, &
                                    dat_file, int_file, xyz_file, axsf_file
       USE io_global,        ONLY : meta_ionode
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL (DP)              :: R, delta_R, x, delta_x
       REAL (DP), ALLOCATABLE :: d_R(:)
       REAL (DP), ALLOCATABLE :: a(:), b(:), c(:), d(:), F(:)
       REAL (DP)              :: ener, ener_0
       INTEGER                     :: i, j, n, atom, image
       INTEGER, PARAMETER          :: max_i = 100
       !
       !
       IF ( .NOT. meta_ionode ) RETURN
       !
       ! ... the *.dat and *.int files are written here
       !
       OPEN( UNIT = iundat, FILE = dat_file, STATUS = "UNKNOWN", &
             ACTION = "WRITE" )
       OPEN( UNIT = iunint, FILE = int_file, STATUS = "UNKNOWN", &
             ACTION = "WRITE" )
       !
       ALLOCATE( d_R( dim ) )
       !
       ALLOCATE( a( num_of_images - 1 ) )
       ALLOCATE( b( num_of_images - 1 ) )
       ALLOCATE( c( num_of_images - 1 ) )
       ALLOCATE( d( num_of_images - 1 ) )
       ALLOCATE( F( num_of_images ) )
       !
       F = 0.D0
       !
       DO image = 2, ( num_of_images - 1 )
          !
          F(image) = - ( grad_pes(:,image) .dot. tangent(:,image) )
          !
       END DO
       !
       react_coord(1) = 0.D0
       !
       DO image = 1, ( num_of_images - 1 )
          !
          d_R = pos(:,( image + 1 )) - pos(:,image)
          !
          R = norm( d_R )
          !
          react_coord(image+1) = react_coord(image) + R
          !
          ! ... cubic interpolation
          !
          a(image) = 2.D0 * ( pes(image) - pes(image+1) ) / R**(3) - &
                     ( F(image) + F(image+1) ) / R**(2)
          !
          b(image) = 3.D0 * ( pes(image+1) - pes(image) ) / R**(2) + &
                     ( 2.D0 * F(image) + F(image+1) ) / R
          !
          c(image) = - F(image)
          !
          d(image) = pes(image)
          !
       END DO          
       !
       DO image = 1, num_of_images
          !
          WRITE( UNIT = iundat, FMT = dat_fmt ) &
              ( react_coord(image) / react_coord(num_of_images) ), &
              ( pes(image) - pes(1) ) * au, error(image)
          !
       END DO
       !
       image = 1
       !
       delta_R = react_coord(num_of_images) / DBLE( max_i )
       !
       DO j = 0, max_i
          !
          R = DBLE( j ) * delta_R 
          !
          IF ( ( R > react_coord(image+1) ) .AND. &
               ( image < ( num_of_images - 1 ) ) ) image = image + 1
          !
          x = R - react_coord(image)
          !
          ener = a(image)*(x**3) + b(image)*(x**2) + c(image)*x + d(image) 
          !
          IF ( j == 0 ) ener_0 = ener
          !
          WRITE( UNIT = iunint, FMT = int_fmt ) &
              ( R / react_coord(num_of_images) ), &
              ( ener - ener_0 ) * au
          !
       END DO
       !
       DEALLOCATE( d_R )
       !
       DEALLOCATE( a )
       DEALLOCATE( b )
       DEALLOCATE( c )
       DEALLOCATE( d )
       DEALLOCATE( F )
       !
       CLOSE( UNIT = iundat )
       CLOSE( UNIT = iunint )
       !
       IF ( lcoarsegrained ) RETURN
       !
       ! ... the *.xyz file is written here
       !
       OPEN( UNIT = iunxyz, FILE = xyz_file, &
             STATUS = "UNKNOWN", ACTION = "WRITE" )
       !
       DO image = 1, num_of_images
          !
          WRITE( UNIT = iunxyz, FMT = '(I5,/)' ) nat
          !
          DO atom = 1, nat
             !
             WRITE( UNIT = iunxyz, FMT = xyz_fmt ) &
                 TRIM( atom_label( ityp( atom ) ) ), &
                 pos(3*atom-2,image)*bohr_radius_angs, &
                 pos(3*atom-1,image)*bohr_radius_angs, &
                 pos(3*atom-0,image)*bohr_radius_angs
             !
          END DO   
          !
       END DO  
       !     
       CLOSE( UNIT = iunxyz )
       !
       ! ... the *.axsf file is written here
       !
       OPEN( UNIT = iunaxsf, FILE = axsf_file, STATUS = "UNKNOWN", &
             ACTION = "WRITE" )
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
       DO image = 1, num_of_images
          !
          WRITE( UNIT = iunaxsf, FMT = '(" PRIMCOORD ",I3)' ) image
          WRITE( UNIT = iunaxsf, FMT = '(I5,"  1")' ) nat
          !
          DO atom = 1, nat
             !
             WRITE( UNIT = iunaxsf, FMT = axsf_fmt ) &
                 TRIM( atom_label(ityp(atom)) ), &
                 pos(3*atom-2,image)*bohr_radius_angs, &
                 pos(3*atom-1,image)*bohr_radius_angs, &
                 pos(3*atom-0,image)*bohr_radius_angs, &
                 - grad_pes(3*atom-2,image)/bohr_radius_angs, &
                 - grad_pes(3*atom-1,image)/bohr_radius_angs, &
                 - grad_pes(3*atom-0,image)/bohr_radius_angs
             !
          END DO   
          !
       END DO  
       !     
       CLOSE( UNIT = iunaxsf )
       !
     END SUBROUTINE write_dat_files
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_output()
       !-----------------------------------------------------------------------
       !
       USE io_files,       ONLY : iunpath
       USE path_variables, ONLY : num_of_images, error, path_length, &
                                  activation_energy, pes, pos, frozen, &
                                  first_last_opt, CI_scheme, Emax_index
       USE path_formats,   ONLY : run_info, run_output
       USE io_global,      ONLY : meta_ionode
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER        :: mode, image
       REAL (DP) :: inter_image_distance
       !
       !
       IF ( .NOT. meta_ionode ) RETURN 
       !
       WRITE( UNIT = iunpath, &
              FMT = '(/,5X,"activation energy (->) = ",F10.6," eV")' ) &
              activation_energy
       WRITE( UNIT = iunpath, &
              FMT = '(5X,"activation energy (<-) = ",F10.6," eV",/)' ) &
              activation_energy + ( pes(1) - pes(num_of_images) ) * au
       !
       WRITE( UNIT = iunpath, FMT = run_info )
       !
       path_length = 0.D0
       !
       DO image = 1, num_of_images
          !
          IF ( image > 1 ) &
             path_length = path_length + &
                           norm( pos(:,image) - pos(:,image-1) )
          !
          WRITE( UNIT = iunpath, FMT = run_output ) &
              image, pes(image) * au, error(image), frozen(image)
          !
       END DO
       !
       inter_image_distance = path_length / DBLE( num_of_images - 1 )
       !
       IF ( CI_scheme == "highest-TS" ) &
          WRITE( UNIT = iunpath, &
                 FMT = '(/,5X,"climbing image = ",I2)' ) Emax_index
       !
       WRITE( UNIT = iunpath, &
              FMT = '(/,5X,"path length",&
                     & T26," = ",F6.3," bohr")' ) path_length
       WRITE( UNIT = iunpath, &
              FMT = '(5X,"inter-image distance", &
                      & T26," = ",F6.3," bohr")' ) inter_image_distance
       !
     END SUBROUTINE write_output
     !
     !-----------------------------------------------------------------------
     SUBROUTINE new_image_init( N_in, outdir )
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine initializes the file needed for the 
       ! ... parallelization among images
       !
       USE io_files,       ONLY : iunnewimage, prefix
       USE path_variables, ONLY : tune_load_balance
       USE mp_global,      ONLY : nimage
       !
       IMPLICIT NONE
       !
       INTEGER,          INTENT(IN) :: N_in
       CHARACTER(LEN=*), INTENT(IN) :: outdir
       !
       !
       IF ( nimage == 1 .OR. .NOT. tune_load_balance ) RETURN
       !
       OPEN( UNIT = iunnewimage, FILE = TRIM( outdir ) // &
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
     SUBROUTINE get_new_image( image, outdir )
       !-----------------------------------------------------------------------
       !
       ! ... this subroutine is used to get the new image to work on
       ! ... the "prefix.BLOCK" file is needed to avoid (when present) that 
       ! ... other jobs try to read/write on file "prefix.newimage" 
       !
       USE io_files,       ONLY : iunnewimage, iunblock, prefix
       USE io_global,      ONLY : ionode
       USE path_variables, ONLY : tune_load_balance
       USE mp_global,      ONLY : nimage
       !
       IMPLICIT NONE
       !
       INTEGER,          INTENT(INOUT) :: image
       CHARACTER(LEN=*), INTENT(IN)    :: outdir
       !
       INTEGER            :: ioerr
       CHARACTER(LEN=256) :: filename
       LOGICAL            :: opened, exists
       !
       !
       IF ( .NOT. ionode ) RETURN
       !
       IF ( nimage > 1 ) THEN
          !
          IF ( tune_load_balance ) THEN
             !
             filename = TRIM( outdir ) // TRIM( prefix ) // '.BLOCK'
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
                   OPEN( UNIT = iunnewimage, FILE = TRIM( outdir ) // &
                       & TRIM( prefix ) // '.newimage' , STATUS = 'OLD' )
                   !
                   READ( iunnewimage, * ) image
                   !
                   CLOSE( UNIT = iunnewimage, STATUS = 'DELETE' )
                   !
                   OPEN( UNIT = iunnewimage, FILE = TRIM( outdir ) // &
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
END MODULE path_io_routines
