!
! Copyright (C) 2002-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! 
!
!----------------------------------------------------------------------------
MODULE io_routines
  !----------------------------------------------------------------------------
  !
  USE kinds,      ONLY :  DP
  USE constants,  ONLY :  AU, BOHR_RADIUS_ANGS
  !      
  IMPLICIT NONE
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE read_restart()
       !-----------------------------------------------------------------------
       !
       USE control_flags,    ONLY : istep, nstep
       USE io_files,         ONLY : iunneb, iunrestart, neb_file   
       USE input_parameters, ONLY : if_pos
       USE neb_variables,    ONLY : pos, vel, num_of_images, dim, PES, &
                                    PES_gradient, suspended_image,     &
                                    Emax, Emin, Emax_index,            &
                                    lquick_min , ldamped_dyn, lmol_dyn
       USE mp_global,        ONLY : mpime, my_pool_id
       USE io_global,        ONLY : ionode_id
       USE mp,               ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !    
       INTEGER :: i, j, ia
       LOGICAL :: file_exists
       !
       ! ... end of local variables
       !
       !
       IF ( mpime == 0 .AND. my_pool_id == 0 ) THEN

          WRITE( UNIT = iunneb, &
                 FMT = '(/,5X,"reading file ", A,/)') TRIM( neb_file )
          !
          OPEN( UNIT = iunrestart, FILE = neb_file, STATUS = "OLD", &
                ACTION = "READ" )
          !
          READ( UNIT = iunrestart, FMT = * )
          READ( UNIT = iunrestart, FMT = * ) istep
          READ( UNIT = iunrestart, FMT = * ) nstep
          READ( UNIT = iunrestart, FMT = * ) suspended_image 
          !
          READ( UNIT = iunrestart, FMT = * )
          !
          READ( UNIT = iunrestart, FMT = * )
          READ( UNIT = iunrestart, FMT = * ) PES(1)
          !
          ia = 0  
          !
          DO j = 1, dim, 3 
             !
             ia = ia + 1
             !
             READ( UNIT = iunrestart, FMT = * ) &
                 pos(j,1),                 &
                 pos((j+1),1),             &
                 pos((j+2),1),             &
                 PES_gradient(j,1),        &
                 PES_gradient((j+1),1),    & 
                 PES_gradient((j+2),1),    &
                 if_pos(1,ia),             &
                 if_pos(2,ia),             &
                 if_pos(3,ia) 
             !
             PES_gradient(:,1) = PES_gradient(:,1) * &
                                 REAL( RESHAPE( if_pos, (/ dim /) ) )
             !
          END DO
          !
          Emin       = PES(1)
          Emax       = PES(1)
          Emax_index = 1
          !
          DO i = 2, num_of_images
             !
             READ( UNIT = iunrestart, FMT = * )
             READ( UNIT = iunrestart, FMT = * ) PES(i)
             !
             DO j = 1, dim, 3 
                !
                READ( UNIT = iunrestart, FMT = * ) &
                    pos(j,i),                 &
                    pos((j+1),i),             &
                    pos((j+2),i),             &
                    PES_gradient(j,i),        &
                    PES_gradient((j+1),i),    &
                    PES_gradient((j+2),i)
                 !
             END DO
             !
             PES_gradient(:,i) = PES_gradient(:,i) * &
                                 REAL( RESHAPE( if_pos, (/ dim /) ) ) 
             !
             IF ( PES(i) <= Emin ) Emin = PES(i)
             !
             IF ( PES(i) >= Emax ) THEN
               !
               Emax = PES(i)
               !
               Emax_index = i
               !
             END IF             
             !  
          END DO
          !
          IF (  lquick_min .OR. ldamped_dyn .OR. lmol_dyn  ) THEN
             !
             READ( UNIT = iunrestart, FMT = * )
             ! 
             DO i = 1, num_of_images
                !
                READ( UNIT = iunrestart, FMT = * )
                !
                DO j = 1, dim, 3
                   !
                   READ( UNIT = iunrestart, FMT = * ) &
                       vel(j,i),                & 
                       vel((j+1),i),            &
                       vel((j+2),i)
                   !
                END DO
                !
                vel(:,i) = vel(:,i) * REAL( RESHAPE( if_pos, (/ dim /) ) )
                !
             END DO
             !
          END IF
          !
          CLOSE( iunrestart )
          !  
       END IF
       !
       ! ... broadcast to all nodes
       !
       CALL mp_bcast( istep, ionode_id )
       CALL mp_bcast( nstep, ionode_id )
       CALL mp_bcast( suspended_image, ionode_id )
       ! 
       CALL mp_bcast( pos, ionode_id )  
       CALL mp_bcast( if_pos, ionode_id )  
       CALL mp_bcast( PES, ionode_id )
       CALL mp_bcast( PES_gradient, ionode_id )
       !
       IF (  lquick_min .OR. ldamped_dyn .OR. lmol_dyn  ) &
          CALL mp_bcast( vel, ionode_id )
       !   
       CALL mp_bcast( Emax, ionode_id )  
       CALL mp_bcast( Emin, ionode_id )
       CALL mp_bcast( Emax_index, ionode_id )
       !
     END SUBROUTINE read_restart
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_restart()
       !-----------------------------------------------------------------------
       !
       USE control_flags,    ONLY : istep, nstep
       USE input_parameters, ONLY : if_pos       
       USE io_files,         ONLY : iunrestart, neb_file 
       USE neb_variables,    ONLY : pos, vel, num_of_images, PES, &
                                    PES_gradient, dim, suspended_image, &
                                    lquick_min , ldamped_dyn, lmol_dyn
       USE formats,          ONLY : energy, restart_first, restart_others, &
                                    velocities
       USE mp_global,        ONLY : mpime, my_pool_id
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: i, j, ia
       !
       ! ... end of local variables
       !
       !
       IF ( mpime == 0 .AND. my_pool_id == 0 ) THEN

          OPEN( UNIT = iunrestart, FILE = neb_file, STATUS = "UNKNOWN", &
                ACTION = "WRITE" )
          !
          WRITE( UNIT = iunrestart, FMT = '("RESTART INFORMATIONS")' )
          !
          WRITE( UNIT = iunrestart, FMT = '(I4)' ) istep
          WRITE( UNIT = iunrestart, FMT = '(I4)' ) nstep
          WRITE( UNIT = iunrestart, FMT = '(I4)' ) suspended_image
          !
          WRITE( UNIT = iunrestart, &
                 FMT = '("ENERGY, POSITIONS AND GRADIENTS")' )
          !
          DO i = 1, num_of_images
             !
             WRITE( UNIT = iunrestart, FMT = '("Image: ",I4)' ) i
             WRITE( UNIT = iunrestart, FMT = energy ) PES(i)
             !
             ia = 0
             !
             DO j = 1, dim, 3
                !
                ia = ia + 1
                !
                IF ( i == 1 ) THEN
                   !
                   WRITE( UNIT = iunrestart, FMT = restart_first ) &
                       pos(j,i),                             &
                       pos((j+1),i),                         &
                       pos((j+2),i),                         &
                       PES_gradient(j,i),                    &
                       PES_gradient((j+1),i),                & 
                       PES_gradient((j+2),i),                &
                       if_pos(1,ia),                         &
                       if_pos(2,ia),                         &
                       if_pos(3,ia) 
                   !
                ELSE
                   !
                   WRITE( UNIT = iunrestart, FMT = restart_others ) &
                       pos(j,i),                              &
                       pos((j+1),i),                          &
                       pos((j+2),i),                          &
                       PES_gradient(j,i),                     &
                       PES_gradient((j+1),i),                 & 
                       PES_gradient((j+2),i)
                   !
                END IF
                !
             END DO
             !
          END DO
          !
          IF (  lquick_min .OR. ldamped_dyn .OR. lmol_dyn  ) THEN
             !
             WRITE( UNIT = iunrestart, &
                    FMT = '("VELOCITIES")' )
             !
             DO i = 1, num_of_images
                !
                WRITE( UNIT = iunrestart, FMT = '("Image: ",I4)' ) i
                !
                DO j = 1, dim, 3
                   !
                   WRITE( UNIT = iunrestart, FMT = velocities ) &
                       vel(j,i),                          & 
                       vel((j+1),i),                      &
                       vel((j+2),i)
                   !
                END DO
                !
             END DO
             ! 
          END IF
          !
          CLOSE( iunrestart )
          !
       END IF
       !
     END SUBROUTINE write_restart
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_dat_files()
       !-----------------------------------------------------------------------
       !
       USE input_parameters,       ONLY : atom_label
       USE cell_base,              ONLY : alat, at
       USE ions_base,              ONLY : ityp, nat
       USE formats,                ONLY : dat_fmt, int_fmt, xyz_fmt, axsf_fmt
       USE basic_algebra_routines, ONLY : norm
       USE supercell,              ONLY : pbc
       USE neb_variables,          ONLY : dim, PES, PES_gradient, pos, &
                                          tangent, num_of_images, error, grad
       USE io_files,               ONLY : iundat, iunint, iunxyz, iunaxsf, &
                                          dat_file, int_file, xyz_file,    &
                                          axsf_file
       USE mp_global,              ONLY : mpime, my_pool_id
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER                     :: j, atom, image
       REAL (KIND=DP)              :: R, delta_R, x
       REAL (KIND=DP), ALLOCATABLE :: d_R(:)
       REAL (KIND=DP), ALLOCATABLE :: a(:), b(:), c(:), d(:), F(:)
       REAL (KIND=DP), ALLOCATABLE :: react_coord(:)
       REAL (KIND=DP)              :: E, E_0
       INTEGER, PARAMETER          :: max_i = 1000 
       !
       ! ... end of local variables
       !
       !
       IF ( mpime == 0 .AND. my_pool_id == 0 ) THEN
          ALLOCATE( d_R( dim ) )
          !
          ALLOCATE( a( num_of_images - 1 ) )
          ALLOCATE( b( num_of_images - 1 ) )
          ALLOCATE( c( num_of_images - 1 ) )
          ALLOCATE( d( num_of_images - 1 ) )
          ALLOCATE( F( num_of_images ) )
          ALLOCATE( react_coord( num_of_images ) )
          !
          F = 0.D0
          !
          DO image = 2, ( num_of_images - 1 )
             !
             F(image) = DOT_PRODUCT( - PES_gradient(:,image) , &
                                       tangent(:,image) )
             !
          END DO
          !
          react_coord(1) = 0.D0
          !
          DO image = 1, ( num_of_images - 1 )
             !
             d_R = pbc( pos(:,( image + 1 )) - pos(:,image) ) 
             !
             R = norm( d_R )
             !
             react_coord(image+1) = react_coord(image) + R
             !
             ! ... cubic interpolation
             !
             a(image) = 2.D0 * ( PES(image) - PES(image+1) ) / R**(3) - &
                        ( F(image) + F(image+1) ) / R**(2)
             !
             b(image) = 3.D0 * ( PES(image+1) - PES(image) ) / R**(2) + &
                        ( 2.D0 * F(image) + F(image+1) ) / R
             !
             c(image) = - F(image)
             !
             d(image) = PES(image)
             !
          END DO
          !
          OPEN( UNIT = iundat, FILE = dat_file, STATUS = "UNKNOWN", &
                ACTION = "WRITE" )
          !  
          DO image = 1, num_of_images
             !
             WRITE( UNIT = iundat, FMT = dat_fmt ) &
                 ( react_coord(image) / react_coord(num_of_images) ), &
                 ( PES(image) - PES(1) ) * AU, &
                 error(image) * ( AU / BOHR_RADIUS_ANGS )
             !
          END DO
          !
          CLOSE( UNIT = iundat )
          !
          OPEN( UNIT = iunint, FILE = int_file, STATUS = "UNKNOWN", &
                ACTION = "WRITE" )
          !
          image = 1
          !
          delta_R = react_coord(num_of_images) / REAL(max_i)
          !
          DO j = 0, max_i
             !
             R = REAL(j) * delta_R 
             !
             IF ( ( R > react_coord(image+1) ) .AND. &
                  ( image < ( num_of_images - 1 ) ) ) image = image + 1
             !
             x = R - react_coord(image)
             !
             E = a(image)*(x**3) + b(image)*(x**2) + c(image)*x + d(image) 
             !
             IF ( j == 0 ) E_0 = E
             !
             WRITE( UNIT = iunint, FMT = int_fmt ) &
                 ( R / react_coord(num_of_images) ), &
                 ( E - E_0 ) * AU
             !
          END DO
          !
          CLOSE( UNIT = iunint )
          !
          DEALLOCATE( d_R )
          !
          DEALLOCATE( a )
          DEALLOCATE( b )
          DEALLOCATE( c )
          DEALLOCATE( d )
          DEALLOCATE( F )
          DEALLOCATE( react_coord )
          !
          OPEN( UNIT = iunxyz, FILE = xyz_file, STATUS = "UNKNOWN", &
                ACTION = "WRITE" )
          !
          DO image = 1, num_of_images
             !
             WRITE( UNIT = iunxyz, FMT = '(I5,/)' ) nat
             !
             DO atom = 1, nat
                !
                WRITE( UNIT = iunxyz, FMT = xyz_fmt ) &
                    TRIM( atom_label(ityp(atom)) ), &
                    pos((3*atom-2),image) * BOHR_RADIUS_ANGS, &
                    pos((3*atom-1),image) * BOHR_RADIUS_ANGS, &
                    pos((3*atom),image) * BOHR_RADIUS_ANGS
                !
             END DO   
             !
          END DO  
          !     
          CLOSE( UNIT = iunxyz )
          !
          OPEN( UNIT = iunaxsf, FILE = axsf_file, STATUS = "UNKNOWN", &
                ACTION = "WRITE" )
          !
          WRITE( UNIT = iunaxsf, FMT = '(" ANIMSTEPS ",I3)' ) num_of_images
          WRITE( UNIT = iunaxsf, FMT = '(" CRYSTAL ")' )
          WRITE( UNIT = iunaxsf, FMT = '(" PRIMVEC ")' )
          WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
              at(1,1) * alat * BOHR_RADIUS_ANGS, &
              at(2,1) * alat * BOHR_RADIUS_ANGS, &
              at(3,1) * alat * BOHR_RADIUS_ANGS  
          WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
              at(1,2) * alat * BOHR_RADIUS_ANGS, &
              at(2,2) * alat * BOHR_RADIUS_ANGS, &
              at(3,2) * alat * BOHR_RADIUS_ANGS
          WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
              at(1,3) * alat * BOHR_RADIUS_ANGS, &
              at(2,3) * alat * BOHR_RADIUS_ANGS, &
              at(3,3) * alat * BOHR_RADIUS_ANGS
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
                    pos((3*atom-2),image) * BOHR_RADIUS_ANGS,  &
                    pos((3*atom-1),image) * BOHR_RADIUS_ANGS,  &
                    pos((3*atom),image) * BOHR_RADIUS_ANGS,    &
                    - grad((3*atom-2),image) / BOHR_RADIUS_ANGS, &
                    - grad((3*atom-1),image) / BOHR_RADIUS_ANGS, &
                    - grad((3*atom),image) / BOHR_RADIUS_ANGS
                !
             END DO   
             !
          END DO  
          !     
          CLOSE( UNIT = iunaxsf )
          !
       END IF
       !
     END SUBROUTINE write_dat_files
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_output()
       !-----------------------------------------------------------------------
       !
       USE io_files,       ONLY : iunneb
       USE neb_variables,  ONLY : num_of_images, PES, error
       USE formats,        ONLY : final_output
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: image
       !
       ! ... end of local variables
       !
       DO image = 1, num_of_images
          !
          WRITE( UNIT = iunneb, FMT = final_output ) &
              image, &
              PES(image) * AU, &
              error(image) * ( AU / BOHR_RADIUS_ANGS )
          !
       END DO
       !
       WRITE( UNIT = iunneb, FMT = '(/,5X,75("-")/)' )       
       !
     END SUBROUTINE write_output
     !
END MODULE io_routines
