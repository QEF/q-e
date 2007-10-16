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
  USE parameters, ONLY : npsx, nchix
  USE radial_grids, ONLY : ndmx, radial_grid_type
  !
  SAVE
  !
  type(radial_grid_type) :: &
       rgrid(npsx)                ! the information on atomic radial grids.
                                  ! NB: some of the subsequent data are therefore redundant 
                                  ! and will be eliminated in due course asap
  REAL(DP) :: &
       jchi(nchix,npsx),         &! total angular momentum of atomic orbitals
       chi(ndmx,nchix,npsx),     &! radial atomic orbitals
       oc(nchix,npsx),           &! atomic level occupation
       rho_at(ndmx,npsx),        &! radial atomic charge density
       rho_atc(ndmx,npsx)         ! radial core charge density
  INTEGER :: &
       msh(npsx),                &! the point at rcut
       nchi(npsx),               &! number of atomic orbitals
       lchi(nchix,npsx)           ! angular momentum of atomic orbitals
  LOGICAL :: &
       nlcc(npsx)                 ! if .TRUE. the atom has nlcc
  !
END MODULE atom
