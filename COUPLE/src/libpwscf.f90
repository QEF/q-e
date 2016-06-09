!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE c2libpwscf(lib_comm,nim,npt,npl,nta,nbn,ndg,retval,infile) BIND(C)
  !----------------------------------------------------------------------------
  !
  ! ... C wrapper for library interface to the Pwscf
  USE ISO_C_BINDING
  !
  IMPLICIT NONE
  !
  INTEGER (kind=C_INT), VALUE :: lib_comm, nim, npt, npl, nta, nbn, ndg
  INTEGER (kind=C_INT), INTENT(OUT) :: retval
  CHARACTER (kind=C_CHAR), INTENT(IN) :: infile(*)
  INTEGER  :: i, lib_comm_, nim_, npt_, npl_, nta_, nbn_, ndg_, retval_
  CHARACTER(LEN=80)  :: infile_
  !
  ! ... Copy C data types to Fortran data types
  lib_comm_ = lib_comm
  nim_ = nim
  npt_ = npt
  npl_ = npl
  nta_ = nta
  nbn_ = nbn
  ndg_ = ndg
  retval = 0
  infile_ = ' '
  !
  ! ... Copying a string from C to Fortran is a bit ugly.
  DO i=1,80
      IF (infile(i) == C_NULL_CHAR) EXIT
      infile_ = TRIM(infile_) // infile(i)
  END DO
  !
  CALL f2libpwscf(lib_comm_,nim_,npt_,npl_,nta_,nbn_,ndg_,retval_,infile_)
  retval = retval_
  !
END SUBROUTINE c2libpwscf
!
!----------------------------------------------------------------------------
SUBROUTINE f2libpwscf(lib_comm,nim,npt,npl,nta,nbn,ndg,retval,infile)
  !----------------------------------------------------------------------------
  !
  ! ... Library interface to the Plane Wave Self-Consistent Field code
  !
  USE environment,       ONLY : environment_start
  USE mp_global,         ONLY : mp_startup
  USE read_input,        ONLY : read_input_file
  USE command_line_options, ONLY: set_command_line
  USE parallel_include
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: lib_comm, nim, npt, npl, nta, nbn, ndg
  INTEGER, INTENT(INOUT) :: retval
  CHARACTER(LEN=80)      :: infile
  !
#if defined(DEBUG_QECOUPLE)
  INTEGER :: me, num, ierr
  CALL MPI_COMM_SIZE(lib_comm,num,ierr)
  IF (ierr /= MPI_SUCCESS) THEN
      CALL MPI_ERROR_STRING(ierr, infile, 80, retval)
      PRINT*,'MPI Error: ', infile
      STOP 100
  END IF
  CALL MPI_COMM_RANK(lib_comm,me,ierr)
  IF (me == 0) THEN
      PRINT*, 'Calling PW library interface with these flags:'
      PRINT*, 'communicator index: ', lib_comm
      PRINT*, 'communicator size:  ', num
      PRINT*, 'nimage: ', nim
      PRINT*, 'npool:  ', npl
      PRINT*, 'ntaskg: ', nta
      PRINT*, 'nband:  ', nbn
      PRINT*, 'ndiag:  ', ndg
      PRINT*, 'input:  "',TRIM(infile),'"'
  END IF
#endif
  !
  CALL set_command_line( nimage=nim, npool=npl, ntg=nta, &
      nband=nbn, ndiag=ndg )
  CALL mp_startup ( my_world_comm=lib_comm )
  CALL environment_start ( 'PWSCF' )
  !
  CALL read_input_file ('PW', infile )
  !
  ! ... Perform actual calculation
  !
  CALL run_pwscf  ( retval )
  !
  CALL stop_run( retval )
  !
END SUBROUTINE f2libpwscf

