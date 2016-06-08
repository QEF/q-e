!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM qecouple
  !----------------------------------------------------------------------------
  !
  ! ... Test program for Q-E library interface
  !
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  !
  INTEGER :: i, exit_status, ierr, ncpu, me, key, new_comm, nargs
  INTEGER :: nimage, npots, npools, ntg, nband, ndiag, nres
  CHARACTER(LEN=80) :: input_file, arg
  !
  ! set defaults
  nimage = 1
  npots  = 1
  npools = 1 
  ntg    = 1
  nband  = 1 
  ndiag  = 1
  nres   = 0
  input_file = ' '
  !
  ! MPI setup
  CALL mpi_init(ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,ncpu,ierr)
  CALL mpi_Comm_rank(MPI_COMM_WORLD,me,ierr)
  !
  ! parse command line flags
  nargs = command_argument_count()
  i = 1
  DO
      CALL getarg(i,arg)
      IF (LEN_TRIM(arg) == 0) EXIT
      !
      i = i + 1
      IF (i > nargs) EXIT
      !
      SELECT CASE ( TRIM(arg) )
          !
      CASE ( '-i', '-in', '-inp', '-input' ) 
          CALL getarg(i, input_file)
          IF ( TRIM (input_file) == ' ') GO TO 15
          i = i + 1
      CASE ( '-ni', '-nimage', '-nimages' ) 
          CALL getarg(i, arg)
          READ ( arg, *, ERR = 15, END = 15) nimage
          i = i + 1
      CASE ( '-nk', '-npool', '-npools') 
          CALL getarg(i, arg)
          READ ( arg, *, ERR = 15, END = 15) npools
          i = i + 1
      CASE ( '-nt', '-ntg', '-ntask_groups') 
          CALL getarg(i, arg)
          READ ( arg, *, ERR = 15, END = 15) ntg
          i = i + 1
      CASE ( '-nb', '-nband', '-nbgrp', '-nband_group') 
          CALL getarg(i, arg)
          READ ( arg, *, ERR = 15, END = 15) nband
          i = i + 1
      CASE ( '-nd', '-ndiag', '-northo', '-nproc_diag', '-nproc_ortho') 
          CALL getarg(i, arg)
          READ ( arg, *, ERR = 15, END = 15) ndiag
          i = i + 1
      CASE ( '-nr', '-nres', '-nreserved') 
          CALL getarg(i, arg)
          READ ( arg, *, ERR = 15, END = 15) nres
          i = i + 1
      CASE DEFAULT
          PRINT*, 'unknown input flag: ',TRIM(arg)
          CALL mpi_abort(MPI_COMM_WORLD,-1,ierr)
      END SELECT
  END DO

15 CONTINUE
  key = MPI_UNDEFINED
  IF (me < (ncpu - nres)) key = 1

  CALL mpi_comm_split(MPI_COMM_WORLD, key, me, new_comm, ierr)

  IF (new_comm /= MPI_COMM_NULL) THEN
      CALL f2libpwscf(new_comm,nimage,npots,npools,ntg,nband,ndiag, &
            exit_status, input_file)
      PRINT *, 'Call to libpwscf finished with exit status', exit_status
  ELSE
      PRINT *, 'Reserved CPU rank:', me, " of", ncpu-1
      exit_status = 0
  END IF
  !
  CALL mpi_finalize(ierr)
  CALL do_stop( exit_status )
  !
  STOP
  !
END PROGRAM qecouple
