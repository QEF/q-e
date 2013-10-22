!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!==-----------------------------------------------------------------------==!
! Wrappers for accessing facilities in the Modules subdirectory
!----------------------------------------------------------------------------
! These init subroutines have to be outside of Fortran modules so they
! can be called from C/C++ or Fortran code
!==-----------------------------------------------------------------------==!
! Configure qm/mm interface for MPI message passing, C version
SUBROUTINE c2qmmm_mpi_config ( qmmm_mode, inter_comm, verb, steps ) BIND(C)
  USE iso_c_binding
  USE qmmm, ONLY: qmmm_config
  IMPLICIT NONE
  !
  INTEGER(C_INT), VALUE, INTENT(in) :: qmmm_mode, inter_comm, verb, steps

  CALL qmmm_config( mode=qmmm_mode, comm=inter_comm, verbose=verb, step=steps )
END SUBROUTINE c2qmmm_mpi_config
!==-----------------------------------------------------------------------==!
! Configure qm/mm interface for MPI message passing, Fortran version
SUBROUTINE f2qmmm_mpi_config ( qmmm_mode, inter_comm, verb, steps )
  USE iso_c_binding
  USE qmmm, ONLY: qmmm_config
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: qmmm_mode, inter_comm, verb, steps

  CALL qmmm_config( mode=qmmm_mode, comm=inter_comm, verbose=verb, step=steps )
END SUBROUTINE f2qmmm_mpi_config
!==-----------------------------------------------------------------------==!
