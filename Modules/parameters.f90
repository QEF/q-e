!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE parameters
  
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: &
       npk    = 40000,  &! max number of k-points               
       ntypx  = 10,     &! max number of different types of atom
       nsx    = ntypx,  &! max number of atomic species (CP)
       natx   = 50       ! max number of atoms for DFT+U+V calculations

END MODULE parameters
