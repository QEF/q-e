!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
MODULE neb_variables
  !---------------------------------------------------------------------------
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: &
       optimization,             &! if .TRUE. the first and the last image
                                  ! will be optimized (see free_minimization)
       conv_neb                   ! .TRUE. if NEB convergence has been
                                  ! achieved
  LOGICAL, ALLOCATABLE :: &
       climbing(:),              &!  
       free_minimization(:)       !
  CHARACTER (LEN=20) :: &
       CI_scheme,                &! Climbing Image scheme
       VEC_scheme                 ! Variable Elastic Constant scheme
  INTEGER :: &
       dim,                      &! dimension of the configuration space
       num_of_images,            &! number of images
       deg_of_freedom,           &! number of degrees of freedom 
                                  ! ( dim - #( of fixed coordinates ) )
       Emax_index,               &! index of the image with the highest energy
       suspended_image            ! last image for which scf has not been
                                  ! achieved
  REAL (KIND=DP) :: &
       neb_thr,                  &! convergence threshold for NEB
       ds,                       &! minimization step 
       damp,                     &! damp coefficient
       temp,                     &! actual temperature ( average over images )
       temp_req                   ! required temperature
  LOGICAL :: &
       lsteep_des,               &! .TRUE. if minimization_scheme = "sd"
       lquick_min,               &! .TRUE. if minimization_scheme = "quick-min"
       ldamped_dyn,              &! .TRUE. if minimization_scheme = "damped-dyn"
       lmol_dyn                   ! .TRUE. if minimization_scheme = "mol-dyn"
  REAL (KIND=DP), ALLOCATABLE :: &
       pos(:,:),                 &!
       vel(:,:),                 &!
       PES(:),                   &!       
       PES_gradient(:,:),        &! 
       grad(:,:),                &!
       norm_grad(:),             &!
       k(:),                     &!
       elastic_gradient(:),      &!
       tangent(:,:),             &!
       mass(:),                  &!
       error(:)                   !
  REAL (KIND=DP) :: &
       k_max,                    &!
       k_min,                    &!
       Eref,                     &!
       Emax,                     &!
       Emin                       !
  !
  CONTAINS
     !
     !----------------------------------------------------------------------
     SUBROUTINE neb_dyn_allocation()
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       ALLOCATE( pos( dim, num_of_images ) )
       !
       IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
          !
          ALLOCATE( vel( dim, num_of_images ) )     
          !
       END IF
       !
       ALLOCATE( grad( dim, num_of_images ) )
       ALLOCATE( norm_grad( num_of_images ) )
       ALLOCATE( PES( num_of_images ) )
       ALLOCATE( PES_gradient( dim, num_of_images ) )       
       ALLOCATE( k( num_of_images ) )
       ALLOCATE( elastic_gradient( dim ) )
       ALLOCATE( tangent( dim, num_of_images ) ) 
       ALLOCATE( mass( dim ) )              
       ALLOCATE( error( num_of_images ) )
       ALLOCATE( climbing( num_of_images ) )
       ALLOCATE( free_minimization( num_of_images ) )
       !           
     END SUBROUTINE neb_dyn_allocation     
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE neb_deallocation()
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       IF ( ALLOCATED( pos ) )               DEALLOCATE( pos )
       IF ( ALLOCATED( vel ) )               DEALLOCATE( vel )
       IF ( ALLOCATED( grad ) )              DEALLOCATE( grad )
       IF ( ALLOCATED( norm_grad ) )         DEALLOCATE( norm_grad )       
       IF ( ALLOCATED( PES ) )               DEALLOCATE( PES )       
       IF ( ALLOCATED( PES_gradient ) )      DEALLOCATE( PES_gradient )
       IF ( ALLOCATED( k ) )                 DEALLOCATE( k )
       IF ( ALLOCATED( elastic_gradient ) )  DEALLOCATE( elastic_gradient )
       IF ( ALLOCATED( tangent ) )           DEALLOCATE( tangent ) 
       IF ( ALLOCATED( mass ) )              DEALLOCATE( mass )       
       IF ( ALLOCATED( error ) )             DEALLOCATE( error )
       IF ( ALLOCATED( climbing ) )          DEALLOCATE( climbing )
       IF ( ALLOCATED( free_minimization ) ) DEALLOCATE( free_minimization )
       !
     END SUBROUTINE neb_deallocation
     !
END MODULE neb_variables
