!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE int_global_variables
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER                      :: N, dim
  INTEGER                      :: old_num_of_images, new_num_of_images
  INTEGER                      :: first_image, last_image
  REAL (DP), ALLOCATABLE  :: old_pos(:,:), new_pos(:,:)
  REAL (DP), ALLOCATABLE  :: old_PES_gradient(:,:), new_PES_gradient(:,:)
  INTEGER,        ALLOCATABLE  :: fix_atom(:)
  REAL (DP), ALLOCATABLE  :: old_V(:), new_V(:) 
  REAL (DP), ALLOCATABLE  :: d_R(:)
  REAL (DP), ALLOCATABLE  :: a(:), b(:), c(:), d(:), F(:)
  REAL (DP), ALLOCATABLE  :: old_mesh(:), new_mesh(:)
  REAL (DP), ALLOCATABLE  :: tangent(:)
  CHARACTER(LEN=256)           :: old_restart_file, new_restart_file
  !
END MODULE int_global_variables
!
!
PROGRAM images_interpolator
  !
  USE kinds,                  ONLY : DP
  USE constants,              ONLY : eps16
  USE path_formats
  USE basic_algebra_routines, ONLY : norm
  USE cell_base,              ONLY : at, alat
  USE int_global_variables
  USE splinelib
  !
  IMPLICIT NONE
  !
  INTEGER             :: i, j, ia
  INTEGER             :: istep, nstep, suspended_image
  INTEGER             :: ierr
  REAL (DP)      :: R, delta_R, x
  LOGICAL             :: no_interpolation
  LOGICAL, EXTERNAL   :: matches 
  CHARACTER (LEN=20)  :: cell_parameters
  CHARACTER (LEN=256) :: input_line
  !
  INTEGER :: iunrestart
  !
  iunrestart = 28
  !
  !
  ! ... the input file is read
  !
  READ(*,*) N
  READ(*,*) old_num_of_images
  READ(*,*) new_num_of_images
  READ(*,*) first_image
  READ(*,*) last_image
  READ(*,'(A)') old_restart_file
  READ(*,'(A)') new_restart_file
  READ(*,*) alat
  !
  READ( UNIT = *, FMT = '(A)', IOSTAT = ierr ) cell_parameters
  !    
  IF ( ierr < 0 ) THEN
     !  
     WRITE(*,'(T2,"read_input: the card CELL_PARAMETERS not found")') 
     STOP
     ! 
  ELSE IF ( ierr > 0 ) THEN
     !
     WRITE(*,'(T2,"read_input: an error occured reading CELL_PARAMETERS")')
     STOP
     !
  END IF
  !
  IF ( .NOT. matches( "CELL_PARAMETERS" , TRIM( cell_parameters ) ) ) THEN
     !
     WRITE(*,'(T2,"read_input: ",A," is not a valid card")') cell_parameters
     STOP
     !
  END IF
  !
  READ( UNIT = *, FMT = *, IOSTAT = ierr ) at
  !
  IF ( ierr < 0 ) THEN
     !   
     WRITE(*,'(T2,"read_input: lattice vectors not found")') 
     STOP
     ! 
  ELSE IF ( ierr > 0 ) THEN
     !
     WRITE(*,'(T2,"read_input:  an error occured reading laccice vectors")') 
     WRITE(*, FMT = lattice_vectors ) at    
     STOP
     !   
  END IF
  !
  dim = 3 * N
  !
  ALLOCATE( old_pos( dim, old_num_of_images ) )
  ALLOCATE( old_PES_gradient( dim, old_num_of_images ) ) 
  ALLOCATE( old_V( old_num_of_images ) )      
  ALLOCATE( new_pos( dim, new_num_of_images ) )
  ALLOCATE( new_PES_gradient( dim, new_num_of_images ) )
  ALLOCATE( new_V( new_num_of_images ) )    
  ALLOCATE( fix_atom( dim ) )   
  ALLOCATE( d_R( dim ) )
  ALLOCATE( tangent( dim ) )   
  ALLOCATE( a( old_num_of_images - 1 ) )
  ALLOCATE( b( old_num_of_images - 1 ) )
  ALLOCATE( c( old_num_of_images - 1 ) )
  ALLOCATE( d( old_num_of_images - 1 ) )
  ALLOCATE( F( old_num_of_images ) )  
  ALLOCATE( old_mesh( old_num_of_images ) )
  ALLOCATE( new_mesh( new_num_of_images ) )
  !
  ! ... the old restart file is read
  !
  OPEN( UNIT = iunrestart, FILE = old_restart_file, STATUS = "OLD", &
       ACTION = "READ" )
  !
  no_interpolation = .FALSE.
  !
  READ( UNIT = iunrestart, FMT = * ) ! RESTART INFORMATION
  READ( UNIT = iunrestart, FMT = * ) istep
  READ( UNIT = iunrestart, FMT = * ) nstep
  READ( UNIT = iunrestart, FMT = * ) suspended_image 
  READ( UNIT = iunrestart, FMT = * ) ! conv_elec
  !
  ! ... read either "ENERGIES, POSITIONS AND GRADIENTS"
  !
  repeat_loop: DO 
     !  
     READ( UNIT = iunrestart, FMT = '(256A)', IOSTAT = ierr ) input_line
     !
     IF ( matches( input_line, &
                   "ENERGIES, POSITIONS AND GRADIENTS" ) ) EXIT repeat_loop
     !
     IF ( ierr /= 0 ) THEN
        !
        WRITE( *, '(/,5X,"an error occured reading",A)' ) old_restart_file
        !
        STOP
        !
     END IF
     !
  END DO repeat_loop  
  !
  READ( UNIT = iunrestart, FMT = * ) ! Image:
  READ( UNIT = iunrestart, FMT = * ) old_V(1)
  !
  ia = 0  
  !
  DO j = 1, dim, 3 
     !
     ia = ia + 1
     !
     READ( UNIT = iunrestart, FMT = * ) &
          old_pos(j,1),                 &
          old_pos((j+1),1),             &
          old_pos((j+2),1),             &
          old_PES_gradient(j,1),        &
          old_PES_gradient((j+1),1),    & 
          old_PES_gradient((j+2),1),    &
          fix_atom(j),                  &
          fix_atom((j+1)),              &
          fix_atom((j+2))
     !
     old_PES_gradient(:,1) = old_PES_gradient(:,1) * fix_atom
     !
  END DO
  !
  DO i = 2, old_num_of_images
     !
     READ( UNIT = iunrestart, FMT = * ) ! Image:
     READ( UNIT = iunrestart, FMT = * ) old_V(i)
     !
     DO j = 1, dim, 3 
        !
        READ( UNIT = iunrestart, FMT = * ) &
             old_pos(j,i),                 &
             old_pos((j+1),i),             &
             old_pos((j+2),i),             &
             old_PES_gradient(j,i),        &
             old_PES_gradient((j+1),i),    &
             old_PES_gradient((j+2),i)
        !
     END DO
     !
     old_PES_gradient(:,i) = old_PES_gradient(:,i) * fix_atom 
     !
  END DO
  !
  CLOSE( UNIT = iunrestart )
  !
  F = 0.D0
  !
  DO i = 2, ( old_num_of_images - 1 )
     !
     ! ... tangent to the path ( normalized )
     !
     tangent(:) = 0.5D0 * ( ( old_pos(:,i+1) - old_pos(:,i) ) /     &
                            norm( old_pos(:,i+1) - old_pos(:,i) ) + &
                            ( old_pos(:,i) - old_pos(:,i-1) ) /     &
                            norm( old_pos(:,i) - old_pos(:,i-1) ) )
     !
     tangent = tangent / norm( tangent )
     !
     F(i) = DOT_PRODUCT( - old_PES_gradient(:,i) , tangent )
     !  
  END DO
  !
  old_mesh(1) = 0.D0
  !
  DO i = 1, ( old_num_of_images - 1 )
     !
     d_R = old_pos(:,(i+1)) - old_pos(:,i)
     !
     R = norm( d_R )
     !
     old_mesh(i+1) = old_mesh(i) + R
     !
     a(i) = 2.D0 * ( old_V(i) - old_V(i+1) ) / R**(3) - &
          ( F(i) + F(i+1) ) / R**(2)
     !   
     b(i) = 3.D0 * ( old_V(i+1) - old_V(i) ) / R**(2) + &
          ( 2.D0 * F(i) + F(i+1) ) / R
     !
     c(i) = - F(i)
     !
     d(i) = old_V(i)
     !
  END DO
  !
  i = first_image
  !
  delta_R = ( old_mesh(last_image) - &
              old_mesh(first_image) ) / DBLE( new_num_of_images - 1 )
  ! 
  DO j = 0, ( new_num_of_images - 1 )
     !
     R = old_mesh( first_image ) + DBLE(j) * delta_R 
     !
     new_mesh(j+1) = R
     !
     check_index: DO
        !
        IF ( ( R > old_mesh(i+1) ) .AND. &
             ( i < ( old_num_of_images - 1 ) ) ) THEN
           !
           i = i + 1
           !
        ELSE
           !
           EXIT check_index
           !
        END IF
        !
     END DO check_index
     !
     x = R - old_mesh( i )
     !
     new_V(j+1) = a(i)*(x**3) + b(i)*(x**2) + c(i)*x + d(i) 
     !
  END DO
  !
  new_mesh = new_mesh / old_mesh(old_num_of_images)  
  old_mesh = old_mesh / old_mesh(old_num_of_images)
  ! 
  CALL dosplineint( old_mesh , old_pos , new_mesh , new_pos )
  !
  new_PES_gradient = 0.D0
  !
  new_PES_gradient(:,1) = old_PES_gradient(:,1)
  !
  new_PES_gradient(:,new_num_of_images) = old_PES_gradient(:,old_num_of_images)
  !
  ! ... the new restart file is written
  !
  OPEN( UNIT = iunrestart, FILE = new_restart_file, STATUS = "UNKNOWN", &
       ACTION = "WRITE" )
  !
  WRITE( UNIT = iunrestart, FMT = '("RESTART INFORMATION")' )
  !
  ! ... by default istep and nstep are set to zero
  !
  WRITE( UNIT = iunrestart, FMT = '(I4)' ) 0
  WRITE( UNIT = iunrestart, FMT = '(I4)' ) 0
  WRITE( UNIT = iunrestart, FMT = '(I4)' ) 0
  WRITE( UNIT = iunrestart, FMT = '(A4)' ) 'F'
  !
  WRITE( UNIT = iunrestart, FMT = '("ENERGIES, POSITIONS AND GRADIENTS")' )
  !
  DO i = 1, new_num_of_images
     !
     WRITE( UNIT = iunrestart, FMT = '("Image: ",I4)' ) i
     WRITE( UNIT = iunrestart, FMT = energy ) new_V(i)
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
                new_pos(j,i),                             &
                new_pos((j+1),i),                         &
                new_pos((j+2),i),                         &
                new_PES_gradient(j,i),                    &
                new_PES_gradient((j+1),i),                & 
                new_PES_gradient((j+2),i),                &
                fix_atom(j),                              &
                fix_atom((j+1)),                          &
                fix_atom((j+2))
           !
        ELSE
           !
           WRITE( UNIT = iunrestart, FMT = restart_others ) &
                new_pos(j,i),                              &
                new_pos((j+1),i),                          &
                new_pos((j+2),i),                          &
                new_PES_gradient(j,i),                     &
                new_PES_gradient((j+1),i),                 & 
                new_PES_gradient((j+2),i)
           !
        END IF
        !
     END DO
     !
  END DO
  !
  WRITE( UNIT = iunrestart, FMT = '("END")' )
  !
  CLOSE( UNIT = iunrestart )   
  !
  IF ( ALLOCATED( old_pos ) )             DEALLOCATE( old_pos )
  IF ( ALLOCATED( old_PES_gradient ) )    DEALLOCATE( old_PES_gradient ) 
  IF ( ALLOCATED( old_V ) )               DEALLOCATE( old_V )      
  IF ( ALLOCATED( new_pos ) )             DEALLOCATE( new_pos )
  IF ( ALLOCATED( new_PES_gradient ) )    DEALLOCATE( new_PES_gradient )  
  IF ( ALLOCATED( new_V ) )               DEALLOCATE( new_V )          
  IF ( ALLOCATED( fix_atom ) )            DEALLOCATE( fix_atom )   
  IF ( ALLOCATED( d_R ) )                 DEALLOCATE( d_R )
  IF ( ALLOCATED( tangent ) )             DEALLOCATE( tangent ) 
  IF ( ALLOCATED( a ) )                   DEALLOCATE( a )
  IF ( ALLOCATED( b ) )                   DEALLOCATE( b )
  IF ( ALLOCATED( c ) )                   DEALLOCATE( c )
  IF ( ALLOCATED( d ) )                   DEALLOCATE( d )
  IF ( ALLOCATED( F ) )                   DEALLOCATE( F )
  IF ( ALLOCATED( old_mesh ) )            DEALLOCATE( old_mesh )
  IF ( ALLOCATED( new_mesh ) )            DEALLOCATE( new_mesh )
  !
END PROGRAM images_interpolator
