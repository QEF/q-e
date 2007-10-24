!
! Copyright (C) 2004 PWSCF-CP-FPMD group
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
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npsx
  USE radial_grids, ONLY : radial_grid_type
  !
  SAVE
  !
  type(radial_grid_type) :: &
       rgrid(npsx)                ! the information on atomic radial grids.
                                  ! NB: some of the subsequent data are therefore redundant 
                                  ! and will be eliminated in due course asap
  INTEGER :: &
       msh(npsx)                  ! the point at rcut
  LOGICAL :: &
       nlcc(npsx)                 ! if .TRUE. the atom has nlcc
  !
END MODULE atom
