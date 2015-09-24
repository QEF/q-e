!
! Copyright (C) 2010-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_arguments()
  !-----------------------------------------------------------------------------
  !
  ! check for presence of command-line option "-plugin_name" or "--plugin_name"
  ! where "plugin_name" has to be set here. If such option is found, variable
  ! "use_plugin_name" is set and usage of "plugin_name" is thus enabled.
  ! Currently implemented: "plumed", "pw2casino" (both case-sensitive)
  !
  USE kinds,         ONLY : DP
  !
  USE io_global,     ONLY : stdout
  !
  USE plugin_flags
  !
  IMPLICIT NONE
  !
  INTEGER  :: iiarg, nargs, i, i0
  CHARACTER (len=1), EXTERNAL ::  lowercase
  CHARACTER (len=256) :: arg
  !
  nargs = command_argument_count()
  ! add here more plugins
  use_plumed = .false.
  use_pw2casino = .false.
  use_environ = .false.
  !
  DO iiarg = 1, nargs 
    CALL get_command_argument( iiarg, plugin_name)
    IF ( plugin_name(1:1) == '-') THEN
       i0 = 1
       IF ( plugin_name(2:2) == '-') i0 = 2
       arg = ' '
       DO i=i0+1, LEN_TRIM (plugin_name)
          arg(i-i0:i-i0) = lowercase (plugin_name(i:i))
       END DO
!       write(0,*) "plugin_name: ", trim(arg)
       ! add here more plugins
       IF ( TRIM(arg)=='plumed' ) THEN
          use_plumed = .true.
       END IF
       IF ( TRIM(arg)=='pw2casino' ) THEN
          use_pw2casino = .true.
       ENDIF
       IF ( TRIM(arg)=='environ' ) THEN
          use_environ = .true.
       ENDIF
    ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE plugin_arguments
!
!----------------------------------------------------------------------------
  SUBROUTINE plugin_arguments_bcast(root,comm)
  !----------------------------------------------------------------------------
  !
  ! broadcast plugin arguments
  !
  USE mp, ONLY : mp_bcast
  USE plugin_flags
  !
  IMPLICIT NONE
  !
  integer :: root
  integer :: comm
  !
  CALL mp_bcast(use_plumed,root,comm)
  !
  CALL mp_bcast(use_pw2casino,root,comm)
  !
  CALL mp_bcast(use_environ,root,comm)
  !
!  write(0,*) "use_plumed: ", use_plumed
  !
  RETURN
  !
END SUBROUTINE plugin_arguments_bcast

