!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE fcp_dyn_calcavg(iter, nelec, vel, acc, force, mass, temp)
  !----------------------------------------------------------------------------
  !
  ! ... calculate physical properties through the dynamics of FCP
  ! ...
  ! ... Variables:
  ! ...   iter:   number of iteration
  ! ...   nelec:  total number of electrons
  ! ...   vel:    velocity of FCP (in a.u.)
  ! ...   acc:    acceleration of FCP (in a.u.)
  ! ...   force:  force acting of FCP (in Ry)
  ! ...   mass:   mass of FCP (in a.u.)
  ! ...   temp:   temperature (in Kelvin)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: iter
  REAL(DP), INTENT(IN) :: nelec
  REAL(DP), INTENT(IN) :: vel
  REAL(DP), INTENT(IN) :: acc
  REAL(DP), INTENT(IN) :: force
  REAL(DP), INTENT(IN) :: mass
  REAL(DP), INTENT(IN) :: temp
  !
  !
  ! >>> calculate properties, here
  !
  !
END SUBROUTINE fcp_dyn_calcavg
