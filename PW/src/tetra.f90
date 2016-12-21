!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
!
MODULE ktetra
  !
  ! ... The variables for the tetrahedron method
  !
  SAVE
  !
  INTEGER :: &
       ntetra            ! number of tetrahedra
  INTEGER, ALLOCATABLE :: &
       tetra(:,:)        ! index of k-points in a given tetrahedron
                         ! shape (4,ntetra)
  !
END MODULE ktetra
