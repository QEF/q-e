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
  ! ... Read variables from command line and broadcast it to all processors
  ! ... ( this is done because there is no guarantee that all processors
  ! ...   have access to command-line options in parallel execution )
  ! ... Interpret QE-specific variables, store the corresponding values
  ! ... Leave the rest (including the code name) in "command_line"
  !
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : root, world_comm
  USE io_global, ONLY : meta_ionode
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Number of arguments in command line
  INTEGER :: nargs = 0
  ! ... QE arguments read from command line
  INTEGER :: nimage_= 1, npool_= 1, npot_= 1, ndiag_ = 0, nband_= 1, ntg_= 1
  CHARACTER(LEN=80) :: input_file_ = ' '
  ! ... Command line arguments not identified
  CHARACTER(LEN=256) :: command_line = ' '
  !
CONTAINS
  !
  SUBROUTINE get_command_line ( )

     IMPLICIT NONE
     INTEGER :: narg
     ! Do not define iargc as external: gfortran doesn't like it
     INTEGER :: iargc
     CHARACTER(LEN=80) :: arg 
     CHARACTER(LEN=6), EXTERNAL :: int_to_char
     !
     command_line = ' '
     nargs = iargc()
     CALL mp_bcast ( nargs, root, world_comm )
     !
     ! ... Only the first node reads and broadcasts
     !
     IF ( .NOT. meta_ionode ) GO TO 20
     !
     arg = ' '
     narg=0
10   CONTINUE
        CALL getarg ( narg, arg )
        narg = narg + 1
        SELECT CASE ( TRIM(arg) )
           CASE ( '-i', '-in', '-inp', '-input' ) 
              CALL getarg ( narg, input_file_ )
              IF ( TRIM (input_file_) == ' ' ) GO TO 15
              narg = narg + 1
           CASE ( '-ni', '-nimage', '-nimages' ) 
              CALL getarg ( narg, arg )
              READ ( arg, *, ERR = 15, END = 15) nimage_
              narg = narg + 1
           CASE ( '-npot', '-npots' ) 
              CALL getarg ( narg, arg )
              READ ( arg, *, ERR = 15, END = 15) npot_
              narg = narg + 1
           CASE ( '-nk', '-npool', '-npools') 
              CALL getarg ( narg, arg )
              READ ( arg, *, ERR = 15, END = 15) npool_
              narg = narg + 1
           CASE ( '-nt', '-ntg', '-ntask_groups') 
              CALL getarg ( narg, arg )
              READ ( arg, *, ERR = 15, END = 15) ntg_
              narg = narg + 1
           CASE ( '-nb', '-nband', '-nbgrp', '-nband_group') 
              CALL getarg ( narg, arg )
              READ ( arg, *, ERR = 15, END = 15) nband_
              narg = narg + 1
           CASE ( '-nd', '-ndiag', '-northo', '-nproc_diag', '-nproc_ortho') 
              CALL getarg ( narg, arg )
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
     CALL mp_bcast( npot_  , root, world_comm ) 
     CALL mp_bcast( npool_ , root, world_comm ) 
     CALL mp_bcast( ntg_   , root, world_comm ) 
     CALL mp_bcast( nband_ , root, world_comm ) 
     CALL mp_bcast( ndiag_ , root, world_comm ) 
     
  END SUBROUTINE get_command_line
  !
END MODULE command_line_options

