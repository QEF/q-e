!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE drhoc( ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, rhocg )
  !-----------------------------------------------------------------------
  !! Calculates the Fourier transform of the core charge.
  !
  USE kinds
  USE constants,   ONLY: pi, fpi
  !
  IMPLICIT NONE
  !
  INTEGER :: ngl
  !! input: the number of g shell
  INTEGER :: mesh
  !! input: the number of radial mesh points
  REAL(DP) :: gl(ngl)
  !! input: the number of G shells
  REAL(DP) :: r(mesh)
  !! input: the radial mesh
  REAL(DP) :: rab(mesh)
  !! input: the derivative of the radial mesh
  REAL(DP) :: rhoc(mesh)
  !! input: the radial core charge
  REAL(DP) :: omega
  !! input: the volume of the unit cell
  REAL(DP) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP) :: rhocg(ngl)
  !! output: the Fourier transform of the core charge
  !
  ! ... local variables
  !
  REAL(DP) :: gx, rhocg1
  ! the modulus of g for a given shell
  ! the Fourier transform
  REAL(DP), ALLOCATABLE ::  aux(:)
  ! auxiliary memory for integration
  INTEGER :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !
  !
!$omp parallel private(aux, gx, rhocg1)
  !
  ALLOCATE( aux(mesh) )
  !
  ! G=0 term
  !
!$omp single
  IF ( gl(1) < 1.0d-8 ) THEN
     DO ir = 1, mesh
        aux(ir) = r(ir)**2 * rhoc(ir)
     ENDDO
     CALL simpson( mesh, aux, rab, rhocg1 )
     rhocg(1) = fpi * rhocg1 / omega
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
!$omp end single
  !
  ! G <> 0 term
  !
!$omp do
  DO igl = igl0, ngl
     gx = SQRT(gl(igl) * tpiba2)
     CALL sph_bes( mesh, r, gx, 0, aux )
     DO ir = 1, mesh
        aux(ir) = r(ir)**2 * rhoc(ir) * aux(ir)
     ENDDO
     CALL simpson( mesh, aux, rab, rhocg1 )
     rhocg(igl) = fpi * rhocg1 / omega
  ENDDO
!$omp end do nowait
  DEALLOCATE( aux )
  !
!$omp end parallel
  !
  RETURN
  !
END SUBROUTINE drhoc

