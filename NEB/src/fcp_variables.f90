!
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_variables
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE mdiis, ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: &
       lfcp             = .FALSE.   ! if .TRUE. FCP is optimized too.
  !
  REAL(DP) :: &
       fcp_mu           = 0.0_DP    ! target Fermi energy,
                                    ! in Hartree
  !
  LOGICAL :: &
       lfcp_linmin      = .FALSE., &! .TRUE. if fcp_scheme = "lm"
       lfcp_newton      = .FALSE., &! .TRUE. if fcp_scheme = "newton"
       lfcp_coupled     = .FALSE.   ! .TRUE. if fcp_scheme = "coupled"
  !
  REAL(DP) :: &
       fcp_thr          = 0.01_DP, &! convergence threshold for FCP relaxation, in eV
       fcp_err_max      = 0.0_DP    ! the largest error, in eV
  !
  INTEGER :: &
       fcp_ndiis        = 4         ! size of DIIS for Newton algorithm
  !
  REAL(DP) :: &
       fcp_rdiis        = 1.0_DP    ! step of DIIS for Newton algorithm
  !
  REAL(DP), ALLOCATABLE :: &
       fcp_nelec(:),               &! the numbers of electrons
       fcp_ef(:),                  &! the Fermi energies, in Hartree
       fcp_dos(:),                 &! the DOSs on Fermi surfaces, in 1/Hartree
       fcp_error(:)                 ! the error of FCP, in eV
  !
  ! ... variables for Line-Minimization
  !
  REAL(DP), ALLOCATABLE :: &
       nelec0(:),                  &! - old #electrons
       force0(:)                    ! - old forces
  !
  LOGICAL, ALLOCATABLE :: &
       firstcall(:)                 ! - first call, or not ?
  !
  ! ... variables for DIIS (coupled with Newton-Raphson)
  !
  LOGICAL :: &
       init_mdiis       = .FALSE.   ! - DIIS is initialized, or not ?
  !
  TYPE(mdiis_type) :: &
       mdiist                       ! - data of DIIS
  !
  ! ... variables for Coupled Method with ionic positions
  !
  REAL(DP) :: &
       fcp_max_volt     = 0.05_DP   ! - maximum estimated voltage, in Hartree
  !
  CONTAINS
     !
     !--------------------------------------------------------------------------
     SUBROUTINE fcp_allocation()
       !--------------------------------------------------------------------------
       !
       USE path_variables, ONLY : num_of_images
       !
       IMPLICIT NONE
       !
       ALLOCATE( fcp_nelec( num_of_images ) )
       ALLOCATE( fcp_ef(    num_of_images ) )
       ALLOCATE( fcp_dos(   num_of_images ) )
       ALLOCATE( fcp_error( num_of_images ) )
       !
       IF ( lfcp_linmin ) THEN
          !
          ALLOCATE( nelec0(    num_of_images ) )
          ALLOCATE( force0(    num_of_images ) )
          ALLOCATE( firstcall( num_of_images ) )
          !
          firstcall = .TRUE.
          !
       END IF
       !
       IF ( lfcp_newton ) THEN
          !
          init_mdiis = .TRUE.
          !
          CALL allocate_mdiis( mdiist, fcp_ndiis, num_of_images, fcp_rdiis, 1 )
          !
       END IF
       !
     END SUBROUTINE fcp_allocation
     !
     !--------------------------------------------------------------------------
     SUBROUTINE fcp_deallocation()
       !--------------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       IF ( ALLOCATED( fcp_nelec ) ) DEALLOCATE( fcp_nelec )
       IF ( ALLOCATED( fcp_ef ) )    DEALLOCATE( fcp_ef )
       IF ( ALLOCATED( fcp_dos ) )   DEALLOCATE( fcp_dos )
       IF ( ALLOCATED( fcp_error ) ) DEALLOCATE( fcp_error )
       !
       IF ( ALLOCATED( nelec0 ) )    DEALLOCATE( nelec0 )
       IF ( ALLOCATED( force0 ) )    DEALLOCATE( force0 )
       IF ( ALLOCATED( firstcall ) ) DEALLOCATE( firstcall )
       !
       IF ( init_mdiis ) CALL deallocate_mdiis( mdiist )
       !
     END SUBROUTINE fcp_deallocation
     !
END MODULE fcp_variables
