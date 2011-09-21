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
  USE kinds,  ONLY : dp
  !
  SAVE
  !
  REAL (dp) :: eth_rps, eth_ns
  REAL (DP), ALLOCATABLE ::  ramtns (:,:,:,:)
  ! raman tensor (3, 3, 3, nat)
  REAL (DP) :: eloptns (3,3,3), dek
  ! electro-optic tensor
  ! delta_xk used for wavefunctions derivation

  LOGICAL :: lraman, elop, done_lraman=.FALSE., done_elop=.FALSE.
  INTEGER :: iuchf, iud2w, iuba2, lrchf, lrd2w, lrba2
  ! unit for chi-wavefuntion
  ! unit for the second derivative of the wavefunction
  ! unit for nscf part of the potential in scf calc.
  ! length of chi-wavefunction record
  ! length of iud2w record
  ! length of iuba2 record

  INTEGER ::  jab(3,3), a1j(6), a2j(6)
  ! 1-D index defining a 3x3 symmetric tensor
  ! a1j and a2j are the two components of the
  ! symmetric tensor defined by jab(3,3)
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

