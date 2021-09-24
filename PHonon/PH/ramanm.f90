!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE ramanm
  !
  !! Module with Raman related variables.
  !
  USE kinds,  ONLY : dp
  !
  SAVE
  !
  REAL (DP) :: eth_rps, eth_ns
  !
  REAL (DP), ALLOCATABLE ::  ramtns (:,:,:,:)
  !! Raman tensor (3, 3, 3, nat)
  REAL (DP) :: eloptns (3,3,3)
  !! electro-optic tensor
  REAL (DP) ::dek
  !! delta_xk used for wavefunctions derivation
  !
  LOGICAL :: lraman, elop, done_lraman=.FALSE., done_elop=.FALSE.
  !
  INTEGER :: iuchf
  !! unit for chi-wavefuntion
  INTEGER :: iud2w
  !! unit for the second derivative of the wavefunction
  INTEGER :: iuba2
  !! unit for nscf part of the potential in scf calc.
  INTEGER :: lrchf
  !! length of chi-wavefunction record
  INTEGER :: lrd2w
  !! length of iud2w record
  INTEGER :: lrba2
  !! length of iuba2 record
  !
  INTEGER :: jab(3,3)
  !! 1-D index defining a 3x3 symmetric tensor
  INTEGER :: a1j(6)
  !! first components of the symmetric tensor defined by
  !! \(\text{jab}(3,3)\).
  INTEGER :: a2j(6)
  !! second component of the symmetric tensor defined by
  !! \(\text{jab}(3,3)\).
  !
  ! When dealing with a 3x3 symmetric tensor A2(i,j), the program considers a
  ! 6-components vector A1(k) whose components are defined by jab, a1j, a2j
  ! --common variables-- in the following way: A1(jab(i,j)) = A2(i,j);
  ! or A2(a1j(k),a2j(k)) = A2(a2j(k),a1j(k)) = A1(k)
  !
  data jab /1, 4, 6, 4, 2, 5, 6, 5, 3/
  data a1j /1, 2, 3, 1, 2, 1 /
  data a2j /1, 2, 3, 2, 3, 3 /
  !
END MODULE ramanm

