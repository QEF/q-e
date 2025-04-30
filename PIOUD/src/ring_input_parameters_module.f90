!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!
!=----------------------------------------------------------------------------=!
!

MODULE ring_input_parameters_module
!
!=----------------------------------------------------------------------------=!
! 
! ... Written by Carlo Sbraccia ( 2003-2006 )
!  this module contains
!  1) the definition of other input parameters taken from NEB
!  2) routines that allocate/deallocate data needed in input
!  Based upon original NEB implementation ( C.S. 17/10/2003 )

  !
  USE kinds,      ONLY : DP
  !
  
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: nat = 1
  REAL(DP) :: alat
  !
  CHARACTER(len=80) :: restart_mode
  ! specify how to start/restart the simulation
  CHARACTER(len=80) :: restart_mode_allowed(3)
  DATA restart_mode_allowed / 'from_scratch', 'restart', 'reset_counters' /
  !
  INTEGER :: nstep_path
  INTEGER :: input_images = 1
  INTEGER :: num_of_images = 0


  LOGICAL :: minimum_image  = .false.

  !
  CHARACTER(len=80) :: opt_scheme = 'quick-min'

  CHARACTER(len=80) :: opt_scheme_allowed(5)
  DATA opt_scheme_allowed / 'quick-min', 'broyden', 'broyden2', 'sd', 'langevin' /
  !
  NAMELIST / PATH / &
                    restart_mode, &
                     nstep_path, num_of_images, & 
                     opt_scheme,     &
                     minimum_image 
!
!    ATOMIC_POSITIONS
!
        REAL(DP), ALLOCATABLE :: pos(:,:)
        INTEGER, ALLOCATABLE :: typ(:)
        !
! ----------------------------------------------------------------------

CONTAINS

  SUBROUTINE allocate_path_input_ions( num_of_images )
    !
    INTEGER, INTENT(in) :: num_of_images
    !
    IF ( allocated( pos ) ) DEALLOCATE( pos )
    IF ( allocated( typ ) ) DEALLOCATE( typ )
    !
    ALLOCATE( pos( 3*nat, num_of_images ) )
    ALLOCATE( typ( nat ) )
    !
    pos(:,:) = 0.0
    !
    RETURN
    !
  END SUBROUTINE allocate_path_input_ions
  !
  SUBROUTINE deallocate_path_input_ions()
    !
    IF ( allocated( pos ) ) DEALLOCATE( pos )
    IF ( allocated( typ ) ) DEALLOCATE( typ )
    !
    RETURN
    !
  END SUBROUTINE deallocate_path_input_ions
  !
!=----------------------------------------------------------------------------=!
!
END MODULE ring_input_parameters_module
!
!=----------------------------------------------------------------------------=!
