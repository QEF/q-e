!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine is inspired by the former routine pw2casino of
! Norbert Nemec
! (C) 2010 by Norbert Nemec <Norbert@Nemec-online.de> 
!----------------------------------------------------------------------------
SUBROUTINE pw2casino()
  !----------------------------------------------------------------------------
  !
  !
  USE control_flags, ONLY : istep
  !
  USE plugin_flags, ONLY : use_pw2casino
  !
  IMPLICIT NONE
  !
  CHARACTER(len=4) :: postfix
  !
  IF ( use_pw2casino ) THEN
    write(postfix,'(i4.4)') istep
    CALL write_casino_wfn( &
             .false., & ! gather
             .true.,  & ! blip
             1.0d0,   & ! multiplicity
             .true.,  & ! binwrite
             .true.,  & ! single_precision_blips
             0,       & ! n_points_for_test
             '.'//postfix)   ! postfix
  ENDIF
  !
  !
END SUBROUTINE pw2casino
