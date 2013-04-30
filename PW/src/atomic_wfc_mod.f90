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
  ! ... The variables needed to describe atomic wavefunctions
  !
  USE kinds, ONLY : dp
  SAVE
  !
  INTEGER :: &
       natomwfc            ! number of (starting) atomic wavefunctions
  COMPLEX(dp), ALLOCATABLE :: &
       swfcatom(:,:)       ! S * (starting) atomic wavefunctions
  CHARACTER(len=30) ::    &!
       starting_wfc,      &! 'random','atomic','file','atomic+random' (default)
       starting_pot,      &! 'atomic' or 'file'
       startingconfig      ! 'input' or 'file'
  !
END MODULE basis
