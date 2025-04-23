!
! Copyright (C) 2003-2006 Quantum ESPRESSO group
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
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... "general" variables :
  !
  LOGICAL :: lneb, lsmd
  !
  LOGICAL :: restart
  !
  LOGICAL :: &
       conv_path                  ! .TRUE. when "path" convergence has been
                                  !        achieved
  LOGICAL :: &
       first_last_opt,           &! if .TRUE. the first and the last image
                                  !           are optimised too.
       use_masses,               &! if .TRUE. mass weighted coordinates are
                                  !           used
       fixed_tan,                &! if. TRUE. the projection is done using the
                                  !           tangent of the average path
       use_freezing,             &! if .TRUE. images are optimised according
                                  !           to their error (see frozen array)
       tune_load_balance          ! if .TRUE. the load balance for image
                                  !           parallelisation is tuned at
                                  !           runtime
  INTEGER :: &
       dim1,                      &! dimension of the configuration space
       num_of_images,            &! number of images
       deg_of_freedom,           &! number of degrees of freedom
                                  ! ( dim1 - #( of fixed coordinates ) )
       pending_image              ! last image for which scf has not been
                                  ! achieved
  REAL(DP) :: &
       ds,                       &! the optimization step
       path_thr,                 &! convergence threshold
       temp_req,                 &! required temperature
       activation_energy,        &! forward activatation energy
       err_max,                  &! the largest error
       path_length                ! length of the path
  LOGICAL :: &
       lsteep_des  = .FALSE.,    &! .TRUE. if opt_scheme = "sd"
       lquick_min  = .FALSE.,    &! .TRUE. if opt_scheme = "quick-min"
       lbroyden    = .FALSE.,    &! .TRUE. if opt_scheme = "broyden"
       lbroyden2   = .FALSE.,    &! .TRUE. if opt_scheme = "broyden2"
       llangevin   = .FALSE.      ! .TRUE. if opt_scheme = "langevin"
  INTEGER :: &
       istep_path,               &! iteration in the optimization procedure
       nstep_path                 ! maximum number of iterations
  !
  ! ... "general" real space arrays
  !
  REAL(DP), ALLOCATABLE :: &
       pes(:),                   &! the potential enrgy along the path
       error(:)                   ! the error from the true MEP
  REAL(DP), ALLOCATABLE :: &
       pos(:,:),                 &! reaction path
       grad_pes(:,:),            &! gradients acting on the path
       tangent(:,:)               ! tangent to the path
  INTEGER, ALLOCATABLE :: &
       fix_atom_pos(:,:)                ! 0 or 1, if 0 fixed atom
  LOGICAL, ALLOCATABLE :: &
       frozen(:)                  ! .TRUE. if the image or mode has not
                                  !        to be optimized
  !
  ! ... "neb specific" variables :
  !
  LOGICAL, ALLOCATABLE :: &
       climbing(:)                ! .TRUE. if the image is required to climb
  CHARACTER(LEN=20) :: &
       CI_scheme                  ! Climbing Image scheme
  INTEGER :: &
       Emax_index                 ! index of the image with the highest energy
  !
  REAL (DP) :: &
       k_max,                    &!
       k_min,                    &!
       Emax,                     &!
       Emin                       !
  !
  ! ... real space arrays
  !
  REAL(DP), ALLOCATABLE :: &
       elastic_grad(:),          &! elastic part of the gradients
       mass(:),                  &! atomic masses
       k(:)                       ! elastic constants
  REAL(DP), ALLOCATABLE :: &
       posold(:,:),              &! old positions (for the quick-min)
       grad(:,:),                &!
       lang(:,:)                  ! langevin random force
  !
  CONTAINS
     !
     !----------------------------------------------------------------------
     SUBROUTINE path_allocation()
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       ALLOCATE( pos( dim1, num_of_images ) )
       !
       ALLOCATE( posold(   dim1, num_of_images ) )
       ALLOCATE( grad(     dim1, num_of_images ) )
       ALLOCATE( grad_pes( dim1, num_of_images ) )
       ALLOCATE( tangent(  dim1, num_of_images ) )
       !
       ALLOCATE( pes(      num_of_images ) )
       ALLOCATE( k(        num_of_images ) )
       ALLOCATE( error(    num_of_images ) )
       ALLOCATE( climbing( num_of_images ) )
       ALLOCATE( frozen(   num_of_images ) )
       !
       ALLOCATE( mass(         dim1 ) )
       ALLOCATE( elastic_grad( dim1 ) )
       !
       ALLOCATE( lang( dim1, num_of_images ) )
       !
     END SUBROUTINE path_allocation
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE path_deallocation()
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       IF ( ALLOCATED( pos ) )          DEALLOCATE( pos )
       IF ( ALLOCATED( posold ) )       DEALLOCATE( posold )
       IF ( ALLOCATED( grad ) )         DEALLOCATE( grad )
       IF ( ALLOCATED( pes ) )          DEALLOCATE( pes )
       IF ( ALLOCATED( grad_pes ) )     DEALLOCATE( grad_pes )
       IF ( ALLOCATED( k ) )            DEALLOCATE( k )
       IF ( ALLOCATED( mass ) )         DEALLOCATE( mass )
       IF ( ALLOCATED( elastic_grad ) ) DEALLOCATE( elastic_grad )
       IF ( ALLOCATED( tangent ) )      DEALLOCATE( tangent )
       IF ( ALLOCATED( error ) )        DEALLOCATE( error )
       IF ( ALLOCATED( climbing ) )     DEALLOCATE( climbing )
       IF ( ALLOCATED( frozen ) )       DEALLOCATE( frozen )
       IF ( ALLOCATED( lang ) )         DEALLOCATE( lang )
       !
       IF ( ALLOCATED( fix_atom_pos ) )     DEALLOCATE( fix_atom_pos )      
       !
     END SUBROUTINE path_deallocation
     !
END MODULE path_variables
