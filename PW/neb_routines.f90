!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE neb_routines
  !-----------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the NEB implementation into the PWSCF code
  ! ... Written by Carlo Sbraccia ( 04-11-2003 )
  !
  USE neb_base,   ONLY : initialize_neb, search_mep
  !
  PRIVATE
  !
  PUBLIC :: initialize_neb, search_mep
  !   
END MODULE neb_routines
