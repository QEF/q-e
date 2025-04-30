!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!
MODULE ring_variables
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all variables taken from NEB path optimisations
  !
  !  
  ! ... Written by Carlo Sbraccia ( 2003-2006 )

  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE


  LOGICAL :: restart
 
  LOGICAL :: &
 
       tune_load_balance          ! if .TRUE. the load balance for image
                                  !           parallelisation is tuned at
                                  !           runtime
  INTEGER :: &
       dim1,                      &! dimension of the configuration space
       num_of_images,            &! number of images
       pending_image              ! last image for which scf has not been
                                  ! achieved
  REAL(DP) :: &
       path_length                ! length of the path

  INTEGER :: &
       istep_path,               &! iteration in the optimization procedure
       nstep_path                 ! maximum number of iterations

  REAL(DP), ALLOCATABLE :: &
       pos(:,:)                   ! reaction path

  
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
     END SUBROUTINE path_allocation
          
     !----------------------------------------------------------------------
     SUBROUTINE path_deallocation()
       !----------------------------------------------------------------------
       
       IMPLICIT NONE
       
       IF ( ALLOCATED( pos ) )          DEALLOCATE( pos )
 
       
     END SUBROUTINE path_deallocation
     
END MODULE ring_variables
