!
! Copyright (C) 2003-2004 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
MODULE path_variables
  !---------------------------------------------------------------------------
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... "general" variables :
  !
  LOGICAL :: &
       first_last_opt,           &! if .TRUE. the first and the last image
                                  ! will be optimized
       conv_path,                &! .TRUE. if NEB convergence has been
                                  ! achieved
       reset_vel                  ! .TRUE. if velocities have to be reset at
                                  ! restart time           
  INTEGER :: &
       dim,                      &! dimension of the configuration space
       num_of_images,            &! number of images
       deg_of_freedom,           &! number of degrees of freedom 
                                  ! ( dim - #( of fixed coordinates ) )
       suspended_image            ! last image for which scf has not been
                                  ! achieved
  REAL (KIND=DP) :: &
       path_thr,                 &! convergence threshold
       damp,                     &! damp coefficient
       temp,                     &! actual temperature ( average over images )
       temp_req,                 &! required temperature
       activation_energy,        &! forward activatation energy

       path_length                ! lentgth of the path
  LOGICAL :: &
       lsteep_des,               &! .TRUE. if minimization_scheme = "sd"
       lquick_min,               &! .TRUE. if minimization_scheme = "quick-min"
       ldamped_dyn,              &! .TRUE. if minimization_scheme = "damped-dyn"
       lmol_dyn                   ! .TRUE. if minimization_scheme = "mol-dyn"
  INTEGER :: &
       istep_path,               &! iteration in the optimization procedure
       nstep_path                 ! maximum number of iterations
  !
  ! ... "general" real space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       ds(:),                    &! the minimization step 
       pes(:),                   &! the potential enrgy along the path
       error(:)                   ! the error from the true MEP
  REAL (KIND=DP), ALLOCATABLE :: &
       pos(:,:),                 &!
       grad_pes(:,:),            &!
       tangent(:,:)               !
  LOGICAL, ALLOCATABLE :: &
       frozen(:)                  ! .TRUE. if the image or mode has not 
                                  !        to be optimized
  !
  ! ... "neb specific" variables :
  !
  LOGICAL, ALLOCATABLE :: &
       climbing(:)                ! .TRUE. if the image is required to climb
  CHARACTER (LEN=20) :: &
       CI_scheme                  ! Climbing Image scheme
  INTEGER :: &
       Emax_index                 ! index of the image with the highest energy
  !
  REAL (KIND=DP) :: &
       k_max,                    &!
       k_min,                    &!
       Eref,                     &!
       Emax,                     &!
       Emin                       !
  !
  ! ... real space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       elastic_grad(:),          &!
       k(:),                     &!  
       react_coord(:),           &! the reaction coordinate (in bohr)
       norm_grad(:)               !
  REAL (KIND=DP), ALLOCATABLE :: &
       vel(:,:),                 &!
       grad(:,:),                &!
       pos_old(:,:),             &!
       grad_old(:,:)              !
  LOGICAL, ALLOCATABLE :: &
       vel_zeroed(:)              ! .TRUE. if the velocity of this image has
                                  !        been reset
  !
  ! ... "smd specific" variables :
  !
  INTEGER :: &
       num_of_modes               ! number of modes
  INTEGER :: &
       Nft,                      &! number of discretization points in the
                                  ! discrete fourier transform
       Nft_smooth                 ! smooth real-space grid
  REAL (KIND=DP) :: &
       ft_coeff                   ! normalization in fourier transformation
  !
  ! ... real space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       pos_star(:,:),            &!
       pes_star(:),              &!
       grad_proj_star(:,:),      &!
       grad_proj(:,:)             !
  !
  ! ... reciprocal space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       ft_pos(:,:),              &!
       ft_pes(:)  ,              &!
       ft_grad(:,:),             &!
       norm_ft_grad(:),          &!
       ft_vel(:,:),              &!
       ft_pos_old(:,:),          &!
       ft_grad_old(:,:)           !
  REAL (KIND=DP), ALLOCATABLE :: &
       ft_error(:)                ! the error from the true MEP
  LOGICAL, ALLOCATABLE :: &
       ft_frozen(:)               ! .TRUE. if the image or mode has not 
                                  !        to be optimized
  LOGICAL, ALLOCATABLE :: &
       ft_vel_zeroed(:)           ! .TRUE. if the velocity of this mode has
                                  !        been reset
  !
  ! ... precomputed values of sin(n*pi*j/N) and n*pi*cos(n*pi*j/N)
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       ft_sin(:,:),              &!
       ft_cos(:,:)                !
  !
  CONTAINS
     !
     !----------------------------------------------------------------------
     SUBROUTINE path_allocation( method )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER (LEN=*), INTENT(IN) :: method
       !
       !
       SELECT CASE ( TRIM( method ) )
       CASE( 'neb' )
          !
          ALLOCATE( ds( num_of_images ) )
          
          !
          ALLOCATE( pos(     dim, num_of_images ) )
          ALLOCATE( pos_old( dim, num_of_images ) )
          !
          ALLOCATE( vel( dim,   num_of_images ) ) 
          ALLOCATE( vel_zeroed( num_of_images ) )    
          !
          ALLOCATE( grad(     dim, num_of_images ) )
          ALLOCATE( grad_old( dim, num_of_images ) )
          ALLOCATE( grad_pes( dim, num_of_images ) )
          ALLOCATE( tangent(  dim, num_of_images ) )
          !
          ALLOCATE( react_coord( num_of_images ) )
          ALLOCATE( norm_grad(   num_of_images ) )
          ALLOCATE( pes(         num_of_images ) )
          ALLOCATE( k(           num_of_images ) )
          ALLOCATE( error(       num_of_images ) )
          ALLOCATE( climbing(    num_of_images ) )
          ALLOCATE( frozen(      num_of_images ) )
          !
          ALLOCATE( elastic_grad( dim ) )
          !
       CASE( 'smd' )
          !
          ! ... real space arrays
          !       
          ALLOCATE( ds(          num_of_images ) )
          ALLOCATE( pes(         num_of_images ) )
          !
          ALLOCATE( pos(       dim, num_of_images ) )       
          ALLOCATE( grad_pes(  dim, num_of_images ) )
          ALLOCATE( grad_proj( dim, num_of_images ) )
          ALLOCATE( tangent(   dim, num_of_images ) )
          !
          ! ... real space "0:( Nft - 1 )" arrays
          !
          ALLOCATE( pes_star( 0:( Nft - 1 ) ) )
          !       
          ALLOCATE( pos_star(       dim, 0:( Nft - 1 ) ) )       
          ALLOCATE( grad_proj_star( dim, 0:( Nft - 1 ) ) )          
          !
          ! ... reciprocal space arrays
          !
          ALLOCATE( ft_pes(        ( Nft - 1 ) ) )
          ALLOCATE( norm_ft_grad(  ( Nft - 1 ) ) )  
          ALLOCATE( ft_error(      ( Nft - 1 ) ) )
          ALLOCATE( ft_frozen(     ( Nft - 1 ) ) )
          ALLOCATE( ft_vel_zeroed( ( Nft - 1 ) ) )
          !
          ALLOCATE( ft_pos(      dim, ( Nft - 1 ) ) )       
          ALLOCATE( ft_grad(     dim, ( Nft - 1 ) ) )
          ALLOCATE( ft_vel(      dim, ( Nft - 1 ) ) )
          ALLOCATE( ft_pos_old(  dim, ( Nft - 1 ) ) )
          ALLOCATE( ft_grad_old( dim, ( Nft - 1 ) ) )
          !
          ALLOCATE( ft_sin( ( Nft - 1 ), 0:( Nft - 1 ) ) )
          ALLOCATE( ft_cos( ( Nft - 1 ), 0:( Nft - 1 ) ) )
          !
          IF ( first_last_opt ) THEN
             ! 
             ALLOCATE( error(      num_of_images ) )
             ALLOCATE( frozen(     num_of_images ) )            
             ALLOCATE( vel_zeroed( num_of_images ) )
             ALLOCATE( norm_grad(  num_of_images ) )
             !
             ALLOCATE( vel(      dim, num_of_images ) )
             ALLOCATE( grad(     dim, num_of_images ) )
             ALLOCATE( pos_old(  dim, num_of_images ) )
             ALLOCATE( grad_old( dim, num_of_images ) )             
             !
          END IF
          !
       END SELECT
       !
     END SUBROUTINE path_allocation     
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE path_deallocation( method )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER (LEN=*), INTENT(IN) :: method
       !
       !
       SELECT CASE ( TRIM( method ) )
       CASE( 'neb' )
          !
          IF ( ALLOCATED( ds ) )             DEALLOCATE( ds )
          IF ( ALLOCATED( pos ) )            DEALLOCATE( pos )
          IF ( ALLOCATED( pos_old ) )        DEALLOCATE( pos_old )
          IF ( ALLOCATED( vel ) )            DEALLOCATE( vel )
          IF ( ALLOCATED( grad ) )           DEALLOCATE( grad )
          IF ( ALLOCATED( grad_old ) )       DEALLOCATE( grad_old )
          IF ( ALLOCATED( react_coord ) )    DEALLOCATE( react_coord )
          IF ( ALLOCATED( norm_grad ) )      DEALLOCATE( norm_grad )       
          IF ( ALLOCATED( pes ) )            DEALLOCATE( pes )       
          IF ( ALLOCATED( grad_pes ) )       DEALLOCATE( grad_pes )
          IF ( ALLOCATED( k ) )              DEALLOCATE( k )
          IF ( ALLOCATED( elastic_grad ) )   DEALLOCATE( elastic_grad )
          IF ( ALLOCATED( tangent ) )        DEALLOCATE( tangent ) 
          IF ( ALLOCATED( error ) )          DEALLOCATE( error )
          IF ( ALLOCATED( climbing ) )       DEALLOCATE( climbing )
          IF ( ALLOCATED( frozen ) )         DEALLOCATE( frozen )
          IF ( ALLOCATED( vel_zeroed ) )     DEALLOCATE( vel_zeroed )          
          !
       CASE( 'smd' )
          !
          ! ... "general" real space arrays
          !       
          IF ( ALLOCATED( ds ) )             DEALLOCATE( ds )
          IF ( ALLOCATED( pes ) )            DEALLOCATE( pes )
          !
          IF ( ALLOCATED( pos ) )            DEALLOCATE( pos )
          IF ( ALLOCATED( grad_pes ) )       DEALLOCATE( grad_pes )
          IF ( ALLOCATED( grad_proj ) )      DEALLOCATE( grad_proj )
          IF ( ALLOCATED( tangent ) )        DEALLOCATE( tangent )
          !
          ! ... real space "0:( Nft - 1 )" arrays
          !
          IF ( ALLOCATED( pes_star ) )       DEALLOCATE( pes_star )
          !
          IF ( ALLOCATED( pos_star ) )       DEALLOCATE( pos_star )
          IF ( ALLOCATED( grad_proj_star ) ) DEALLOCATE( grad_proj_star )
          !       
          ! ... reciprocal space arrays
          !
          IF ( ALLOCATED( ft_pes ) )         DEALLOCATE( ft_pes )
          IF ( ALLOCATED( norm_ft_grad ) )   DEALLOCATE( norm_ft_grad )  
          IF ( ALLOCATED( ft_error ) )       DEALLOCATE( ft_error )
          IF ( ALLOCATED( ft_frozen ) )      DEALLOCATE( ft_frozen )
          IF ( ALLOCATED( ft_vel_zeroed ) )  DEALLOCATE( ft_vel_zeroed )
          !
          IF ( ALLOCATED( ft_pos ) )         DEALLOCATE( ft_pos )       
          IF ( ALLOCATED( ft_grad ) )        DEALLOCATE( ft_grad )
          IF ( ALLOCATED( ft_vel ) )         DEALLOCATE( ft_vel )
          IF ( ALLOCATED( ft_pos_old ) )     DEALLOCATE( ft_pos_old )
          IF ( ALLOCATED( ft_grad_old ) )    DEALLOCATE( ft_grad_old )
          !
          IF ( ALLOCATED( ft_sin ) )         DEALLOCATE( ft_sin )
          IF ( ALLOCATED( ft_cos ) )         DEALLOCATE( ft_cos )
          !
          IF ( first_last_opt ) THEN
             !
             IF ( ALLOCATED( pos_old ) )        DEALLOCATE( pos_old )
             IF ( ALLOCATED( vel ) )            DEALLOCATE( vel )
             IF ( ALLOCATED( grad ) )           DEALLOCATE( grad )
             IF ( ALLOCATED( grad_old ) )       DEALLOCATE( grad_old )
             IF ( ALLOCATED( norm_grad ) )      DEALLOCATE( norm_grad )
             IF ( ALLOCATED( error ) )          DEALLOCATE( error )
             IF ( ALLOCATED( frozen ) )         DEALLOCATE( frozen )
             !
          END IF
          !
       END SELECT
       !
     END SUBROUTINE path_deallocation
     !
END MODULE path_variables
