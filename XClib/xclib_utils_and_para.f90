!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------
MODULE xclib_utils_and_para
!-----------------------------------------------------
!! MPI stuff and error vars.
!
#if defined (__MPI)
        !
        !     Include file for MPI
        !
#if defined (__MPI_MODULE)
        USE mpi
#else
        INCLUDE 'mpif.h'
#endif
#else
        ! dummy world and null communicator
        INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
        INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
        INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
#endif
        INTEGER, PARAMETER :: stdout = 6
        !! standard output unit
        !
        LOGICAL :: nowarning = .FALSE.
        !! switch for warning messages
        !
        INTEGER :: inside_error = 0
        !$acc declare copyin( inside_error )
        !! index to recover error type inside gpu regions (see error_msg)
        !
        CHARACTER(LEN=28) :: error_msg(7)
        DATA error_msg / 'Invalid ID for LDA exchange ', &
                         'Invalid ID for LDA corr.    ', &
                         'Invalid ID for GGA exchange ', &
                         'Invalid ID for GGA corr.    ', &
                         'Invalid ID for MGGA         ', &
                         'Bad args. in EXPINT function', &
                         'Bad args in wggax_analy_erfc'  /
        !
   CONTAINS
      !
      SUBROUTINE xc_inside_error( in_err )
        !! Recover the error type inside GPU kernel regions.
        IMPLICIT NONE
        !$acc routine seq
        INTEGER, INTENT(IN) :: in_err
        !$acc atomic write
        inside_error = in_err
        RETURN
      END SUBROUTINE
      !
END MODULE xclib_utils_and_para
