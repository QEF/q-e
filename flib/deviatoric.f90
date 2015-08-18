!
! Copyright (C) 2010 Davide Ceresoli
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE impose_deviatoric_strain ( at_old, at )
  !---------------------------------------------------------------------
  !
  !     Impose a pure deviatoric (volume-conserving) deformation
  !     Needed to enforce volume conservation in variable-cell MD/optimization
  !
  USE kinds, ONLY: dp
  IMPLICIT NONE
  REAL(dp), INTENT(in)    :: at_old(3,3)
  REAL(dp), INTENT(inout) :: at(3,3)
  REAL(dp) :: tr, omega, omega_old

  tr = (at(1,1)+at(2,2)+at(3,3))/3.d0
  tr = tr - (at_old(1,1)+at_old(2,2)+at_old(3,3))/3.d0
  ! Commented out, while waiting for better idea:
  ! it breaks the symmetry of hexagonal lattices - PG
  ! at(1,1) = at(1,1) - tr
  ! at(2,2) = at(2,2) - tr
  ! at(3,3) = at(3,3) - tr
  ! print '("difference in trace: ",e12.4)', tr

  CALL volume (1.d0, at_old(1,1), at_old(1,2), at_old(1,3), omega_old)
  CALL volume (1.d0, at(1,1), at(1,2), at(1,3), omega)
  at = at * (omega_old / omega)**(1.d0/3.d0)

END SUBROUTINE impose_deviatoric_strain
!
!---------------------------------------------------------------------
SUBROUTINE impose_deviatoric_strain_2d ( at_old, at )
  !---------------------------------------------------------------------
  !     Modif. of impose_deviatoric_strain but for
  !     Area conserving deformation (2DSHAPE) added by Richard Charles Andrew
  !     Physics Department, University if Pretoria,
  !     South Africa, august 2012
  !
  USE kinds, ONLY: dp
  IMPLICIT NONE
  REAL(dp), INTENT(in)    :: at_old(3,3)
  REAL(dp), INTENT(inout) :: at(3,3)
  REAL(dp) :: omega, omega_old
  INTEGER :: i, j

  CALL volume (1.d0, at_old(1,1), at_old(1,2), at_old(1,3), omega_old)
  CALL volume (1.d0, at(1,1), at(1,2), at(1,3), omega)
  DO i = 1,3
     DO j = 1,3
        IF (j==3) THEN
           at(i,j) = at(i,j) ! DON'T CHANGE IN z- DIRECTION IF 2DSHAPE
        ELSE
           at(i,j) = at(i,j) * (omega_old / omega)**(1.d0/2.d0)
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE impose_deviatoric_strain_2d
!
!---------------------------------------------------------------------
SUBROUTINE impose_deviatoric_stress ( sigma )
  !---------------------------------------------------------------------
  !
  !     Impose a pure deviatoric stress
  !
  USE kinds, ONLY: dp
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  REAL(dp), INTENT(inout) :: sigma(3,3)
  REAL(dp) :: tr

  tr = (sigma(1,1)+sigma(2,2)+sigma(3,3))/3.d0
  sigma(1,1) = sigma(1,1) - tr
  sigma(2,2) = sigma(2,2) - tr
  sigma(3,3) = sigma(3,3) - tr
  WRITE (stdout,'(5x,"Volume is kept fixed: isostatic pressure set to zero")')

END SUBROUTINE impose_deviatoric_stress
!
!---------------------------------------------------------------------
SUBROUTINE impose_deviatoric_stress_2d ( sigma )
  !---------------------------------------------------------------------
  !
  !     Modif. of impose_deviatoric_stress but for
  !     Area conserving deformation (2DSHAPE) added by Richard Charles Andrew
  !     Physics Department, University if Pretoria,
  !     South Africa, august 2012
  !
  USE kinds, ONLY: dp
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  REAL(dp), INTENT(inout) :: sigma(3,3)
  REAL(dp) :: tr
  tr = (sigma(1,1)+sigma(2,2))/2.d0
  sigma(1,1) = sigma(1,1) - tr
  sigma(2,2) = sigma(2,2) - tr
  WRITE (stdout,'(5x,"Area is kept fixed: isostatic in-plane pressure in xy set to zero")')
END SUBROUTINE impose_deviatoric_stress_2d
