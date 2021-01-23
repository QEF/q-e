!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE parameters
  !
  !! Upper limits on k-points, atoms and supercell size.
  !
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: npk = 40000
  !! max number of k-points
  INTEGER, PARAMETER :: ntypx = 10
  !! max number of different types of atom
  INTEGER, PARAMETER :: nsx = ntypx
  !! max number of atomic species (CP)
  INTEGER, PARAMETER :: natx = 50
  !! max number of atoms for DFT+U+V calculations
  INTEGER, PARAMETER :: sc_size = 1
  !! Defines the supercell in DFT+U+V as composed by the unit cells located
  !! by (n1,n2,n3) in primitive vectors base with \(-\text{sc_size} \leq ni
  !! \leq \text{sc_size}\) and \((2\text{sc_size}+1)^3\) is the number of cells.
  INTEGER, PARAMETER :: nsolx  = 10
  !! max number of solvents (RISM)

END MODULE parameters
