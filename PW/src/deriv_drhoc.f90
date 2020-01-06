!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE deriv_drhoc( ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, drhocg )
  !--------------------------------------------------------------------------
  !! Calculates the Fourier transform of \(d\text{Rho}_c/dG\).
  !
  USE kinds
  USE constants,  ONLY : pi, fpi
  !
  IMPLICIT NONE
  !
  INTEGER :: ngl
  !! input: the number of g shell
  INTEGER :: mesh
  !! input: the number of radial mesh points
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! input: the number of G shells
  REAL(DP), INTENT(IN) :: r(mesh)
  !! input: the radial mesh
  REAL(DP), INTENT(IN) :: rab(mesh)
  !! input: the derivative of the radial mesh
  REAL(DP), INTENT(IN) :: rhoc(mesh)
  !! input: the radial core charge
  REAL(DP), INTENT(IN) :: omega
  !! input: the volume of the unit cell
  REAL(DP), INTENT(IN) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP), INTENT(OUT) :: drhocg(ngl)
  !! Fourier transform of d Rho_c/dG
  !
  ! ... local variables
  !
  REAL(DP) :: gx, rhocg1
  ! the modulus of g for a given shell
  ! the fourier transform
  REAL(DP), ALLOCATABLE :: aux(:)
  ! auxiliary memory for integration
  INTEGER :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !
  ! G=0 term
  !
  IF (gl(1) < 1.0d-8) THEN
     drhocg(1) = 0.0d0
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  ! G <> 0 term
  !
!$omp parallel private(aux, gx, rhocg1)
  !
  ALLOCATE( aux(mesh) )
!$omp do
  DO igl = igl0, ngl
     gx = SQRT( gl(igl) * tpiba2 )
     DO ir = 1, mesh
        aux(ir) = r(ir)*rhoc(ir)*( r(ir) * COS(gx*r(ir)) /       &
                                      gx - SIN(gx*r(ir)) / gx**2 )
     ENDDO
     CALL simpson( mesh, aux, rab, rhocg1 )
     drhocg(igl) = fpi / omega * rhocg1
  ENDDO
!$omp end do nowait
  DEALLOCATE( aux )
  !
!$omp end parallel
  !
  RETURN
  !
END SUBROUTINE deriv_drhoc

