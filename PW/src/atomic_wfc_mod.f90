!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE basis
  !
  !! The variables needed to describe the atomic wavefunctions.
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  INTEGER :: natomwfc
  !! number of (starting) atomic wavefunctions
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)
  !! (starting) atomic wavefunctions
  COMPLEX(DP), ALLOCATABLE :: swfcatom(:,:)
  !! S * (starting) atomic wavefunctions
  CHARACTER(len=30) :: starting_wfc
  !! It can be: 'random', 'atomic', 'file', 'atomic+random' (default)
  CHARACTER(len=30) :: starting_pot
  !! It can be 'atomic' or 'file'
  CHARACTER(len=30) :: startingconfig
  !! It can be 'input' or 'file'
  !
END MODULE basis
