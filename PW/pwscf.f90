!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM pwscf
  !-----------------------------------------------------------------------
  !
  ! ... Plane Wave Self-Consistent Field
  !
  USE wvfct,            ONLY : gamma_only
  USE varie,            ONLY : nstep, istep, conv_elec, conv_ions
  USE io_files,         ONLY : nd_nmbr
  USE global_version,   ONLY : version_number
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=9) :: code = 'PWSCF'
  EXTERNAL          :: date_and_tim
  !
  !
  ! use ".FALSE." to disable all clocks except the total cpu time clock
  ! use ".TRUE."  to enable clocks
  ! CALL init_clocks( .FALSE. )
  CALL init_clocks( .TRUE. )
  !
  CALL start_clock('PWSCF')
  !
  gamma_only = .FALSE.
  !
  CALL startup( nd_nmbr, code, version_number )
  CALL init_run
  !
  istep = 0
  main_loop : DO WHILE ( istep < nstep )
     istep = istep + 1
     CALL electrons
     IF ( .NOT. conv_elec ) CALL stop_pw( conv_elec )
     CALL ions
     IF ( conv_ions ) EXIT main_loop
     CALL hinit1
  END DO main_loop
  !
  CALL punch
  !
  CALL stop_pw( conv_ions )
  !
  STOP
  !
END PROGRAM pwscf


