!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE gcscf_variables
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: &
       lgcscf   = .FALSE.   ! if .TRUE. GC-SCF is used.
  !
  REAL(DP) :: &
       gcscf_mu = 0.0_DP    ! target Fermi energy,
                            ! in Hartree
  !
  REAL(DP), ALLOCATABLE :: &
       gcscf_nelec(:),     &! the numbers of electrons
       gcscf_ef(:)          ! the Fermi energies, in Hartree
  !
  CONTAINS
     !
     !--------------------------------------------------------------------------
     SUBROUTINE gcscf_allocation()
       !--------------------------------------------------------------------------
       !
       USE path_variables, ONLY : num_of_images
       !
       IMPLICIT NONE
       !
       ALLOCATE( gcscf_nelec( num_of_images ) )
       ALLOCATE( gcscf_ef(    num_of_images ) )
       !
     END SUBROUTINE gcscf_allocation
     !
     !--------------------------------------------------------------------------
     SUBROUTINE gcscf_deallocation()
       !--------------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       IF ( ALLOCATED( gcscf_nelec ) ) DEALLOCATE( gcscf_nelec )
       IF ( ALLOCATED( gcscf_ef ) )    DEALLOCATE( gcscf_ef )
       !
     END SUBROUTINE gcscf_deallocation
     !
END MODULE gcscf_variables
