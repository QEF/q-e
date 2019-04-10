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
  USE constants, ONLY : pi
  !
  IMPLICIT NONE
  INTEGER :: ibrav
  REAL(DP) :: celldm(6), angle(3)
  !
  REAL(DP) :: at(3,3), omega, R(3,3)
  REAL(DP),PARAMETER :: grad_to_rad = pi/180
  !
  NAMELIST /system/ ibrav, celldm, angle
  WRITE(*,*) "This code reads a &system/ namelist from standard input"
  WRITE(*,*) "containing ibrav and celldm(1..6) as in INPUT_PW."
  WRITE(*,*) "If you specify three angles as angle(1..3) in degrees, the cell"
  WRITE(*,*) "will be rotated accordingly around the three cartesian axis"
  READ(*,system)
  CALL latgen( ibrav, celldm, at(:,1), at(:,2), at(:,3), omega )
  !
  IF (ANY(angle/=0._dp)) THEN
    R = rot(angle(1)*grad_to_rad, angle(2)*grad_to_rad, angle(3)*grad_to_rad) 
    at = matmul(R,at)
  ENDIF

  !   at=at/celldm(1)
  !   CALL recips(at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3))
  WRITE(*,'(a)') "Unit cell (bohr):"
  WRITE(*,'(3f25.15)') at(:,1)
  WRITE(*,'(3f25.15)') at(:,2)
  WRITE(*,'(3f25.15)') at(:,3)
  at=at/celldm(1)
  WRITE(*,'(a,f12.6,a)') "Unit cell in units of alat (=", celldm(1)," bohr):"
  WRITE(*,'(3f15.8)') at(:,1)
  WRITE(*,'(3f15.8)') at(:,2)
  WRITE(*,'(3f15.8)') at(:,3)
  !WRITE(*,'(a,es24.15)') "Volume (bohr^3):", omega
  !
 CONTAINS
 function rotx (theta) RESULT(R)
 IMPLICIT NONE
   REAL(DP),INTENT(in) :: theta
   REAL(DP) :: R(3,3)
   R(:,1) = (/ 1._dp,          0._dp,            0._dp/)
   R(:,2) = (/ 0._dp, cos(theta), -sin(theta) /)
   R(:,3) = (/ 0._dp, sin(theta), cos(theta)  /)
 endfunction
 function roty (theta) RESULT(R)
 IMPLICIT NONE
   REAL(DP),INTENT(in) :: theta
   REAL(DP) :: R(3,3)
   R(:,1) = (/ cos(theta), 0._dp, sin(theta) /)
   R(:,2) = (/ 0._dp,1._dp,0._dp/)
   R(:,3) = (/ -sin(theta), 0._dp, cos(theta)/)
 endfunction
 function rotz (theta) RESULT(R)
 IMPLICIT NONE
   REAL(DP),INTENT(in) :: theta
   REAL(DP) :: R(3,3)
   R(:,1) = (/cos(theta),-sin(theta), 0._dp /)
   R(:,2) = (/sin(theta), cos(theta), 0._dp /)
   R(:,3) = (/0._dp,0._dp,1._dp /)
 endfunction
 function rot(alpha,beta,gamma) RESULT(R)
   IMPLICIT NONE
   REAL(DP),INTENT(in) :: alpha,beta,gamma
   REAL(DP) :: R(3,3)
   R = matmul(matmul(rotx(alpha),roty(beta)), rotz(gamma))
 endfunction

END PROGRAM ibrav2cell
!----------------------------------------------------------------------




