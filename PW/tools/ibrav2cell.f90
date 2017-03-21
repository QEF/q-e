!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is ibrav2cellributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present ibrav2cellribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
PROGRAM ibrav2cell
!----------------------------------------------------------------------
  !
  USE Kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: ibrav
  REAL(DP) :: celldm(6)
  !
  REAL(DP) :: at(3,3), omega, alat
  !
  NAMELIST /system/ ibrav, celldm
  READ(*,system)
  CALL latgen( ibrav, celldm, at(:,1), at(:,2), at(:,3), omega )
  !   at=at/celldm(1)
  !   CALL recips(at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3))
  !WRITE(*,'(a)') "Unit cell (bohr):"
  WRITE(*,'(3es24.15)') at(:,1)
  WRITE(*,'(3es24.15)') at(:,2)
  WRITE(*,'(3es24.15)') at(:,3)
  !WRITE(*,'(a,es24.15)') "Volume (bohr^3):", omega
  !
END PROGRAM ibrav2cell
!----------------------------------------------------------------------
