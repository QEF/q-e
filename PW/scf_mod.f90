!
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE scf
  !  
  !  This module contains variables and auxiliary routines needed for
  !  the self-consistent cycle
  !
  !  ROUTINES: allocate_scf_type
  !
  USE kinds,       ONLY : DP
  !
  USE lsda_mod,   ONLY : nspin
  USE ldaU,       ONLY : lda_plus_u
  USE gvect,      ONLY : nrxx, ngm
  !
  SAVE
  !
TYPE scf_type
  REAL(DP),    ALLOCATABLE :: of_r(:,:) ! the charge density in R-space
!  COMPLEX(DP), ALLOCATABLE :: of_g(:,:) ! the charge density in G-space
END TYPE scf_type

  type (scf_type) :: rho

  REAL(DP) :: v_of_0    ! vltot(G=0)      
  REAL(DP), ALLOCATABLE :: &
       vr(:,:),        &! the Hartree + xc potential in real space
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grig)
       rho_core(:),    &! the core charge in real space
       tauk(:,:),      &! kinetic energy density in real space (dense grid)
       kedtau(:,:),    &! position dependent kinetic energy enhancement factor
                        ! used in META-GGA in real space (smooth grid)
       kedtaur(:,:)     ! position dependent kinetic energy enhancement factor
                        ! used in META-GGA in real space (dense grid)
  COMPLEX(DP), ALLOCATABLE :: &
       rhog(:,:),      &! the charge density in reciprocal space
       rhog_core(:),   &! the core charge in reciprocal space
       taukg(:,:)       ! the kinetic energy density in reciprocal space

CONTAINS

 subroutine allocate_scf_type ( rho )
 type (scf_type) :: rho
 allocate (rho%of_r( nrxx, nspin))
 return
 end subroutine allocate_scf_type

 subroutine deallocate_scf_type ( rho )
 type (scf_type) :: rho
 if (allocated(rho%of_r)) deallocate(rho%of_r)
 return
 end subroutine deallocate_scf_type
 !
END MODULE scf
