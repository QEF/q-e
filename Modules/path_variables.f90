!
! Copyright (C) 2003-2005 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE path_variables
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all variables needed by path optimisations
  !
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER, PARAMETER :: history_ndim = 8
  !
  ! ... "general" variables :
  !
  LOGICAL :: &
       first_last_opt,           &! .TRUE. if the first and the last image
                                  !        have to be optimized
       conv_path,                &! .TRUE. if "path" convergence has been
                                  !        achieved
       reset_vel,                &! .TRUE. if velocities have to be reset at
                                  !        restart time
       use_multistep,            &! .TRUE. if multistep has to be used in smd
                                  !        optimization
       write_save,               &! .TRUE. if the save file has to be written
       free_energy,              &! .TRUE. for free-energy calculations
       fixed_tan,                &! if. TRUE. the projection is done using the
                                  !           tangent of the average path
       use_freezing,             &! if .TRUE. images are optimised according
                                  !           to their error (see frozen array)
       tune_load_balance          ! if .TRUE. the load balance for image
                                  !           parallelisation is tuned at 
                                  !           runtime
  INTEGER :: &
       dim,                      &! dimension of the configuration space
       num_of_images,            &! number of images
       init_num_of_images,       &! number of images used in the initial
                                  ! discretization (SMD only)
       deg_of_freedom,           &! number of degrees of freedom 
                                  ! ( dim - #( of fixed coordinates ) )
       suspended_image            ! last image for which scf has not been
                                  ! achieved
  REAL (KIND=DP) :: &
       ds,                       &! the optimization step
       path_thr,                 &! convergence threshold
       damp,                     &! damp coefficient
       temp_req,                 &! required temperature
       activation_energy,        &! forward activatation energy
       err_max,                  &! the largest error
       path_length,              &! lentgth of the path
       path_length_av             ! lentgth of the average path
  LOGICAL :: &
       lsteep_des  = .FALSE.,    &! .TRUE. if opt_scheme = "sd"
       lquick_min  = .FALSE.,    &! .TRUE. if opt_scheme = "quick-min"
       ldamped_dyn = .FALSE.,    &! .TRUE. if opt_scheme = "damped-dyn"
       lmol_dyn    = .FALSE.,    &! .TRUE. if opt_scheme = "mol-dyn"
       lbroyden    = .FALSE.,    &! .TRUE. if opt_scheme = "broyden"
       llangevin   = .FALSE.      ! .TRUE. if opt_scheme = "langevin"
  INTEGER :: &                   
       istep_path,               &! iteration in the optimization procedure
       nstep_path,               &! maximum number of iterations
       av_counter                 ! number of steps used to compute averages
  LOGICAL :: &
       reset_broyden = .FALSE.    ! used to reset the broyden subspace
  !
  ! ... "general" real space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       pes(:),                   &! the potential enrgy along the path
       norm_tangent(:),          &!
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
       lang(:,:),                &! langevin random force 
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
  !
  ! ... real space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       pos_in_av(:),             &!
       pos_fin_av(:),            &!
       pos_in_h(:,:),            &!
       pos_fin_h(:,:)             !
  REAL (KIND=DP), ALLOCATABLE :: &
       pos_star(:,:)              !
  !
  ! ... reciprocal space arrays
  !
  REAL (KIND=DP), ALLOCATABLE :: &
       ft_pos(:,:),              &!
       ft_pos_av(:,:),           &!
       ft_pos_h(:,:,:)            !
  !
  ! ... Y. Kanai variabiles for combined smd/cp dynamics :
  !
  INTEGER, PARAMETER :: smx = 20    ! max number of images
  INTEGER, PARAMETER :: smmi = 4    ! a parameter for  polynomial interpolation
                                    ! # of replicas used for interpolation
  LOGICAL :: &
       smd_cp,                     &! regular CP calculation
       smd_lm,                     &! String method w/ Lagrange Mult.
       smd_opt,                    &! CP for 2 replicas, initial & final
       smd_linr,                   &! linear interpolation
       smd_polm,                   &! polynomial interpolation
       smd_stcd
  INTEGER :: &
       smd_p,                      &! sm_p = 0 .. SM_P replica
       smd_kwnp,                   &! # of points used in polm
       smd_codfreq,                &!
       smd_forfreq,                &! frequency of calculating Lag. Mul
       smd_wfreq,                  &!
       smd_lmfreq,                 &
       smd_maxlm                    ! max_ite = # of such iteration allowed
  REAL(KIND=DP) :: &
       smd_tol,                    &! tolrance on const in terms of
                                    ! [alpha(k) - alpha(k-1)] - 1/sm_P
       smd_ene_ini = 1.D0,         &
       smd_ene_fin = 1.D0
  !
  TYPE smd_ptr
    !
    REAL(KIND=DP), POINTER :: d3(:,:)
    !
  END TYPE smd_ptr
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
       IF ( method == "smd" ) THEN
          !
          ! ... real space arrays
          !
          ALLOCATE( pos_in_av(  dim ) )
          ALLOCATE( pos_fin_av( dim ) )
          !
          ALLOCATE( pos_in_h(  dim, history_ndim ) )
          ALLOCATE( pos_fin_h( dim, history_ndim ) )
          !
          ALLOCATE( pos_star( dim, 0:( Nft - 1 ) ) )
          !
          ALLOCATE( ft_pos(    dim, ( Nft - 1 ) ) )
          ALLOCATE( ft_pos_av( dim, ( Nft - 1 ) ) )
          !
          ALLOCATE( ft_pos_h( dim, ( Nft - 1 ), history_ndim ) )
          !
          IF ( llangevin ) THEN
             !
             ALLOCATE( lang( dim, num_of_images ) )
             !
          END IF
          !
       END IF
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
       IF ( method == "smd" ) THEN
          !
          IF ( ALLOCATED( pos_in_av ) )   DEALLOCATE( pos_in_av )
          IF ( ALLOCATED( pos_fin_av ) )  DEALLOCATE( pos_fin_av )
          !
          IF ( ALLOCATED( pos_in_h ) )    DEALLOCATE( pos_in_h )
          IF ( ALLOCATED( pos_fin_h ) )   DEALLOCATE( pos_fin_h )
          !
          IF ( ALLOCATED( pos_star ) )    DEALLOCATE( pos_star )
          !
          IF ( ALLOCATED( ft_pos ) )      DEALLOCATE( ft_pos )
          IF ( ALLOCATED( ft_pos_av ) )   DEALLOCATE( ft_pos_av )
          !
          IF ( llangevin ) THEN
             !
             IF ( ALLOCATED( lang ) )     DEALLOCATE( lang )
             !
          END IF
          !
       END IF
       !
     END SUBROUTINE path_deallocation
     !
END MODULE path_variables
