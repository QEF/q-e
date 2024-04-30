!
! Copyright (C) 2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE starting_scf
  !
  !! The variables needed to describe how to start an scf calculation
  !
  SAVE
  !
  CHARACTER(len=30) :: starting_wfc
  !! It can be: 'random', 'atomic', 'file', 'atomic+random' (default)
  CHARACTER(len=30) :: starting_pot
  !! It can be 'atomic' or 'file'
  CHARACTER(len=30) :: startingconfig
  !! It can be 'input' or 'file'
  !
END MODULE starting_scf
