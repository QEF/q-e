!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE command_line_options
  !----------------------------------------------------------------------------
  !
  ! ... Utilities to read command-line variables and to set related variables:
  ! ... "get_command_line()" with no arguments: 
  ! ...                      reads the command line,
  ! ...                      interprets QE-specific variables,
  ! ...                      stores the corresponding values
  ! ...                      (nimage, npool, ntg, nband, ndiag),
  ! ...                      broadcasts them to all processors,
  ! ...                      leaves the rest of the command line 
  ! ...                      (including the code name) in "command_line"
  ! ... "get_command_line(input_command_line)" with a string argument:
  ! ...                      as above, but reading from the input string
  ! ... Variables are read on one processor and broadcast to all others
  ! ... because there is no guarantee that all processors have access to
  ! ... command-line options in parallel execution.
  ! ... "set_command_line" directly sets nimage, npool, ntg, nband, ndiag.
  ! ... Useful to initialize parallelism when QE is used as a library
  !
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : root, world_comm
  USE io_global, ONLY : meta_ionode
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... Number of arguments in command line
  INTEGER :: nargs = 0
  ! ... QE arguments read from command line
  INTEGER :: nimage_= 1, npool_= 1, ndiag_ = 0, nband_= 1, ntg_= 1
  ! ... Indicate if using library init
  LOGICAL :: library_init = .FALSE.
  ! ... input file name read from command line
  CHARACTER(LEN=256) :: input_file_ = ' '
  ! ... Command line arguments that were not identified and processed
  CHARACTER(LEN=512) :: command_line = ' '
  !
CONTAINS
  !
  SUBROUTINE get_command_line ( input_command_line )
     IMPLICIT NONE
     CHARACTER(LEN=*), OPTIONAL :: input_command_line 
     INTEGER :: narg
     LOGICAL :: read_string
     CHARACTER(LEN=256) :: arg 
     CHARACTER(LEN=6), EXTERNAL :: int_to_char
     !
     command_line = ' '
     read_string = PRESENT ( input_command_line )
     !
     ! command line parameters have already been set via set_command_line()
     IF (library_init) GO TO 20
     !
     IF (read_string) THEN
        nargs = my_iargc ( input_command_line )
     ELSE
        nargs = command_argument_count()
     ENDIF
     CALL mp_bcast ( nargs, root, world_comm )
     !
     ! ... Only the first node reads and broadcasts
     !
     IF ( .NOT. meta_ionode ) GO TO 20
     !
     arg = ' '
     narg=0
10   CONTINUE
        IF (read_string) THEN
           CALL my_getarg ( input_command_line, narg, arg )
        ELSE
           CALL get_command_argument ( narg, arg )
        ENDIF
        narg = narg + 1
        SELECT CASE ( TRIM(arg) )
           CASE ( '-i', '-in', '-inp', '-input' ) 
           IF (read_string) THEN
              CALL my_getarg ( input_command_line, narg, input_file_ )
           ELSE
              CALL get_command_argument ( narg, input_file_ )
           ENDIF
              IF ( TRIM (input_file_) == ' ' ) GO TO 15
              narg = narg + 1
           CASE ( '-ni', '-nimage', '-nimages', '-npot' ) 
              IF (read_string) THEN
                 CALL my_getarg ( input_command_line, narg, arg )
              ELSE
                 CALL get_command_argument ( narg, arg )
              ENDIF
              READ ( arg, *, ERR = 15, END = 15) nimage_
              narg = narg + 1
           CASE ( '-nk', '-npool', '-npools') 
              IF (read_string) THEN
                 CALL my_getarg ( input_command_line, narg, arg )
              ELSE
                 CALL get_command_argument ( narg, arg )
              ENDIF
              READ ( arg, *, ERR = 15, END = 15) npool_
              narg = narg + 1
           CASE ( '-nt', '-ntg', '-ntask_groups') 
              IF (read_string) THEN
                 CALL my_getarg ( input_command_line, narg, arg )
              ELSE
                 CALL get_command_argument ( narg, arg )
              ENDIF
              READ ( arg, *, ERR = 15, END = 15) ntg_
              narg = narg + 1
           CASE ( '-nb', '-nband', '-nbgrp', '-nband_group') 
              IF (read_string) THEN
                 CALL my_getarg ( input_command_line, narg, arg )
              ELSE
                 CALL get_command_argument ( narg, arg )
              ENDIF
              READ ( arg, *, ERR = 15, END = 15) nband_
              narg = narg + 1
           CASE ( '-nd', '-ndiag', '-northo', '-nproc_diag', '-nproc_ortho') 
              IF (read_string) THEN
                 CALL my_getarg ( input_command_line, narg, arg )
              ELSE
                 CALL get_command_argument ( narg, arg )
              ENDIF
              READ ( arg, *, ERR = 15, END = 15) ndiag_
              narg = narg + 1
           CASE DEFAULT
              command_line = TRIM(command_line) // ' ' // TRIM(arg)
        END SELECT
        IF ( narg > nargs ) GO TO 20
     GO TO 10
     ! ... something wrong: notify and continue
15   CALL infomsg ('get_command_line', 'unexpected argument # ' // &
                  & int_to_char(narg) // ':' //TRIM(arg))
     narg = narg + 1
     GO TO 10
     ! ... normal exit
20   CONTINUE
     CALL mp_bcast( command_line, root, world_comm ) 
     CALL mp_bcast( input_file_ , root, world_comm ) 
     CALL mp_bcast( nimage_, root, world_comm ) 
     CALL mp_bcast( npool_ , root, world_comm ) 
     CALL mp_bcast( ntg_   , root, world_comm ) 
     CALL mp_bcast( nband_ , root, world_comm ) 
     CALL mp_bcast( ndiag_ , root, world_comm ) 
     
  END SUBROUTINE get_command_line
  !
  INTEGER FUNCTION my_iargc ( input_command_line )
     IMPLICIT NONE
     CHARACTER(LEN=*), INTENT(IN) :: input_command_line 
     CHARACTER(LEN=1) :: previous, current
     INTEGER :: i

     my_iargc = 0
     previous = ' '
     DO i=1,LEN_TRIM(input_command_line)
        current = input_command_line(i:i)
        IF ( current /= ' ' .AND. previous == ' ' ) my_iargc = my_iargc+1
        previous = current
     END DO
  
  END FUNCTION my_iargc
  !
  SUBROUTINE my_getarg ( input_command_line, narg, arg )
     IMPLICIT NONE
     CHARACTER(LEN=*), INTENT(IN) :: input_command_line 
     INTEGER, INTENT(IN) :: narg 
     CHARACTER(LEN=*), INTENT(OUT) :: arg
     CHARACTER(LEN=1) :: previous, current
     INTEGER :: iarg, i, indx

     iarg = 0
     previous = ' '
     arg = ' '
     indx= 0
     DO i=1,LEN_TRIM(input_command_line)
        current = input_command_line(i:i)
        IF ( current /= ' ' .AND. previous == ' ' ) iarg = iarg+1
        IF ( iarg == narg ) THEN
           indx = indx + 1
           arg(indx:indx) = current           
           IF ( indx == LEN(arg) ) RETURN
        ELSE IF ( iarg > narg ) THEN
           RETURN
        END IF
        previous = current
     END DO

  END SUBROUTINE my_getarg 

  SUBROUTINE set_command_line ( nimage, npool, ntg, nband, ndiag)
     ! directly set command line options without going through the command line
     IMPLICIT NONE

     INTEGER, INTENT(IN), OPTIONAL :: nimage, npool, ntg, nband, ndiag
     !
     IF ( PRESENT(nimage) ) nimage_ = nimage
     IF ( PRESENT(npool)  ) npool_  = npool
     IF ( PRESENT(ntg)    ) ntg_    = ntg
     IF ( PRESENT(nband)  ) nband_  = nband
     IF ( PRESENT(ndiag)  ) ndiag_  = ndiag
     !
     library_init = .TRUE.
     !
  END SUBROUTINE set_command_line
  !
END MODULE command_line_options

