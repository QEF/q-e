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
  USE parameters, ONLY : npsx, ndmx, nchix
  !
  SAVE
  !
  REAL(DP) :: &
       zmesh(npsx),              &! the atomic charge for mesh generation
       xmin(npsx),               &! initial linear mesh point (obsolescent)
       dx(npsx),                 &! linear interval for logaritmic mesh (obs.)
       r(ndmx,npsx),             &! radial logaritmic mesh
       rab(ndmx,npsx),           &! derivative of the radial mesh
       jchi(nchix,npsx),         &! total angular momentum of atomic orbitals
       chi(ndmx,nchix,npsx),     &! radial atomic orbitals
       oc(nchix,npsx),           &! atomic level occupation
       rho_at(ndmx,npsx),        &! radial atomic charge density
       rho_atc(ndmx,npsx)        ! radial core charge density
  INTEGER :: &
       mesh(npsx),               &! number of mesh points
       msh(npsx),                &! the point at rcut
       nchi(npsx),               &! number of atomic orbitals
       lchi(nchix,npsx)           ! angular momentum of atomic orbitals
  LOGICAL :: &
       numeric(npsx),            &! .TRUE. for pot. in numeric form (ie always)
       nlcc(npsx)                 ! .TRUE. if the potential has nonlinear
                                  !        core correction
  !
END MODULE atom
