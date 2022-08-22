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
  USE constants,          ONLY: pi, fpi
  USE upf_acc_interfaces, ONLY: simpson, sph_bes_acc
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
  REAL(DP), ALLOCATABLE ::  aux(:), gxv(:)
  ! auxiliary memory for integration
  INTEGER :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !
  !--------PROVISIONAL------ gl def needs to be adjusted---------
  ALLOCATE( gxv(ngl) )
  DO igl = 1, ngl
    gxv(igl) = SQRT(gl(igl) * tpiba2)
  ENDDO
  !--------------------------------------------------------------
  !
#if !defined(_OPENACC)
!$omp parallel private(aux,gx,rhocg1)
#endif
  !
  ALLOCATE( aux(mesh) )
#if defined(_OPENACC)
  !$acc data present_or_copyin(gl,r,rab,rhoc) present_or_copyout(rhocg)
  !$acc data create(aux)
#endif
  !
  ! G=0 term
  !
#if !defined(_OPENACC)
!$omp single
#endif
  IF ( gl(1) < 1.0d-8 ) THEN
     !$acc parallel loop
     DO ir = 1, mesh
        aux(ir) = r(ir)**2 * rhoc(ir)
     ENDDO
     !$acc kernels
     CALL simpson( mesh, aux, rab, rhocg1 )
     rhocg(1) = fpi * rhocg1 / omega
     !$acc end kernels
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
#if !defined(_OPENACC)
!$omp end single
#endif
!$acc end data
  !
  ! G <> 0 term
  !
#if defined(_OPENACC)
!$acc parallel loop gang private(aux) copyin(gxv)
#else
!$omp do
#endif
  DO igl = igl0, ngl
     gx = gxv(igl)
     CALL sph_bes_acc( mesh, r, gx, 0, aux )
     !$acc loop vector
     DO ir = 1, mesh
       aux(ir) = r(ir)**2 * rhoc(ir) * aux(ir)
     ENDDO
     CALL simpson( mesh, aux, rab, rhocg1 )
     rhocg(igl) = fpi * rhocg1 / omega
  ENDDO
#if !defined(_OPENACC)
!$omp end do nowait
#endif
  !
  !$acc end data
  !
  DEALLOCATE( aux )
#if !defined(_OPENACC)
  !$omp end parallel
#endif
  !
  DEALLOCATE( gxv )
  !
  RETURN
  !
END SUBROUTINE drhoc
