!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Common for d3toten
!
module g0aux
  USE kinds, only: DP
  real(DP), pointer:: vlocg0(:,:) ! local potential at q+G for q=0
  complex(DP), pointer:: vkb0 (:,:) ! contains beta functions at q=0
  complex(DP), pointer:: d0rc(:,:)! contain the rhoc for q=0
end module g0aux
!
!   the units of the files and the record lengths
!
module units_d3
  integer:: iudqwf,        &! the unit with | Pc d/du(q) psi_{k+q} >
            iud0qwf,       &! the unit with | Pc d/du(0) psi_{k+q} >
            iud0rho         ! the unit where q=0 delta rho is written
end module units_d3
!
!    the name of the files
!
module d0rho
  character(len=256) :: fild0rho
end module d0rho
!
!     the variable needed to describe the patterns when q=0
!
module modesg0
  USE kinds, only: DP
  integer :: nsymg0,        &! the number of symmetries of the crystal
             nirrg0          ! the number of irreducible representation
!
  integer, pointer :: npertg0(:)  ! the number of perturbations per IR
  complex(DP), pointer :: ug0(:,:), tg0(:,:,:,:)
 ! ug0: transformation modes patterns
 ! tg0: the symmetry in the base of pattern (q=0)
end module modesg0
!
!    third order dynamical matrices (auxiliary)
!
module d3aux
  USE kinds, only: DP
  complex(DP), allocatable :: &
       d3dyn_aux1(:,:,:), d3dyn_aux2(:,:,:), d3dyn_aux3(:,:,:), &
       d3dyn_aux4(:,:,:), d3dyn_aux5(:,:,:), d3dyn_aux6(:,:,:), &
       d3dyn_aux7(:,:,:), d3dyn_aux8(:,:,:), d3dyn_aux9(:,:,:)
end module d3aux
!
!    third order dynamical matrix
!
module thirdorder
  USE kinds, only: DP
  complex(DP), allocatable ::  d3dyn(:,:,:)
  ! third order dynamical matrix
  complex(DP), allocatable :: psidqvpsi(:,:)
  !    <psi| dqV |psi>
  real(DP) :: ethr_ph ! eigenvalues convergence threshold
  real(DP), allocatable ::  ef_sh(:) ! E_Fermi shift
  integer :: istop
  logical :: wraux, recv
end module thirdorder
!
!    test variables
!
module testvar
  USE kinds, only: DP
  real(DP) :: testreal(50)
  integer :: testint(50)
  logical :: testflag(50)
end module testvar
!
!  the units of the files and the record lengths
!
module units_d3ph
  integer :: &
       iuef,           &! unit with ef_sh
       iupdqvp,        &! unit with <psi| dqV |psi>
       iupd0vp,        &! unit with <psi| d0V |psi>
       lrpdqvp,        &! length of <psi| dV  |psi>
       iudpdvp_1,      &! unit with <dqpsi| dqV |psi>
       iudpdvp_2,      &! unit with <dqpsi| d0V |psi>
       iudpdvp_3,      &! unit with <d0psi| dqV |psi>
       lrdpdvp          ! length of <dpsi |  dV |psi> records
end module units_d3ph
!
! In the parallel version of the program some loop on perturbations
! may be split betweem pools. npert_i and npert_f are the initial
! and final value for a counter on the modes to be split among pools
!
module npert_mod
  integer :: &
       npert_i,    &! starting value for the mode counter
       npert_f      ! final value for the mode counter
end module npert_mod
!
! Variables used for computing and writing only selected modes at q=0
! --the first index of the dthird matrix--
!
module q0modes
  integer :: q0mode_todo(300)  ! list of the q=0 modes to be computed
!
  logical :: &
       q0mode(300),     &! if .true. this mode is to be computed
       wrmode(300),     &! if .true. this mode is to be written
       allmodes          ! it is .true. if you are at gamma and you
                         !  want to compute all the modes
end module q0modes

module d3com
  use g0aux
  use units_d3
  use units_d3ph
  use d0rho
  use d3aux
  use thirdorder
  use testvar
  use modesg0
  use npert_mod
  use q0modes
end module d3com
