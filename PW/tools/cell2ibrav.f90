!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
PROGRAM ibrav2cell
!----------------------------------------------------------------------
  !
  USE Kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: ibrav, ios
  REAL(DP) :: celldm(6)
  !
  REAL(DP) :: a(3,3), omega, e(3,3), alat
  INTEGER, EXTERNAL :: at2ibrav
  !
  WRITE(*, '("Enter the lattice vectors (one per line)")')
  READ(*,*) a(:,1)
  READ(*,*) a(:,2)
  READ(*,*) a(:,3)
  WRITE(*, '(a)') "Enter alat (or 'x' to skip if axis are in bohr)"
  READ(*,*,iostat=ios) alat
  IF(ios/=0) alat=1._dp
  ibrav = at2ibrav (a(:,1), a(:,2), a(:,3))
  celldm = 0._dp
  CALL  at2celldm (ibrav,alat,a(:,1), a(:,2), a(:,3),celldm)
  WRITE(*,'("ibrav = ",i6)') ibrav
  WRITE(*,'(" celldm(1) = ",f15.8)') celldm(1)
  IF( celldm(2) /= 0._dp) WRITE(*,'(" celldm(2) = ",f15.8)') celldm(2)
  IF( celldm(3) /= 0._dp) WRITE(*,'(" celldm(3) = ",f15.8)') celldm(3)
  IF( celldm(4) /= 0._dp) WRITE(*,'(" celldm(4) = ",f15.8)') celldm(4)
  IF( celldm(5) /= 0._dp) WRITE(*,'(" celldm(5) = ",f15.8)') celldm(5)
  IF( celldm(6) /= 0._dp) WRITE(*,'(" celldm(6) = ",f15.8)') celldm(6)
  !
  CALL latgen( ibrav, celldm, e(:,1), e(:,2), e(:,3), omega )
  WRITE(*, '("New lattice vectors in bohr:")')
  WRITE(*,'(3f15.8)') e(:,1)
  WRITE(*,'(3f15.8)') e(:,2)
  WRITE(*,'(3f15.8)') e(:,3)
  WRITE(*, '("New lattice vectors in old alat:")')
  WRITE(*,'(3f15.8)') e(:,1)/alat
  WRITE(*,'(3f15.8)') e(:,2)/alat
  WRITE(*,'(3f15.8)') e(:,3)/alat
  WRITE(*, '("New lattice vectors in new alat:")')
  WRITE(*,'(3f15.8)') e(:,1)/celldm(1)
  WRITE(*,'(3f15.8)') e(:,2)/celldm(1)
  WRITE(*,'(3f15.8)') e(:,3)/celldm(1)
  WRITE(*,'("Discrepancy old/new in bohr = ", f12.6)') DSQRT(SUM((a*alat-e)**2))

  !
END PROGRAM ibrav2cell
!----------------------------------------------------------------------
