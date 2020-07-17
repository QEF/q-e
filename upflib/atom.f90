!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE atom
  !
  ! ... The variables needed to describe the atoms and related quantities
  !
  USE radial_grids, ONLY : radial_grid_type
  !
  SAVE
  !
  ! the information on atomic radial grids
  type(radial_grid_type), allocatable, target :: rgrid(:)
  !
  ! last grid point with r < rcut = 10 a.u. (currently set in read_pseudo)
  ! used to perform integration of atomic functions going to 0 at large r
  ! (e.g. Vloc(r)-Ze^2/r), cutting off numerical noise
  ! FIXME: should be included in "rgrid"
  INTEGER, ALLOCATABLE :: msh(:) 
  !
END MODULE atom
