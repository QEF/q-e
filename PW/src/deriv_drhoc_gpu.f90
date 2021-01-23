!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE compute_drhocg_gpu_m
#if defined(__CUDA)
  !
  USE cudafor
  !
  IMPLICIT NONE
  TYPE(dim3) :: drhocg_threads = dim3(32,8,1) 
  PUBLIC     :: drhocg_threads 
  CONTAINS
      ATTRIBUTES(global) SUBROUTINE compute_drhocg_gpu( n, tpiba2, omega, gl, r, rhoc, rab, mesh, drhocg )
      !
      USE cudafor
      USE kinds,     ONLY: DP 
      USE constants, ONLY: pi, fpi, eps14
      !
      IMPLICIT NONE 
      !
      INTEGER, VALUE :: n, mesh 
      REAL(DP),VALUE :: tpiba2, omega
      REAL(DP),DEVICE,INTENT(IN)  :: gl(n), r(mesh), rhoc(mesh), rab(mesh)
      REAL(DP),DEVICE,INTENT(OUT) :: drhocg(n) 
      !   
      INTEGER  :: tx, ty, igl, ir
      REAL(DP) :: mysum, val, gx, x
      ! 
      tx = threadIdx%x
      ty = threadIdx%y
      !
      igl = (blockIdx%x - 1) * blockDim%y + ty
      !
      IF (igl > n ) RETURN 
      ! 
      gx = SQRT(gl(igl) * tpiba2)
      mysum = 0_dp
      !
      DO ir = tx, mesh, blockDim%x
        !
        val = r(ir)*rhoc(ir)*( r(ir) * COS(gx*r(ir)) / &
                                  gx - SIN(gx*r(ir)) / gx**2 ) * rab(ir)
        !
        IF (ir == 1 .OR. ir == mesh) THEN
          mysum = mysum + val
        ELSEIF (MOD(ir,2)) THEN
          mysum = mysum + 2._dp*val
        ELSE
          mysum = mysum + 4._dp*val
        ENDIF
      ENDDO
      !
      ! Reduce by warp
      val = __shfl_down(mysum,1)
      mysum = mysum + val
      val = __shfl_down(mysum,2)
      mysum = mysum + val
      val = __shfl_down(mysum,4)
      mysum = mysum + val
      val = __shfl_down(mysum,8)
      mysum = mysum + val
      val = __shfl_down(mysum,16)
      mysum = mysum + val
      !
      IF (tx == 1) THEN
        drhocg(igl) = fpi * mysum / (3._dp * omega)
      ENDIF
      !
    END SUBROUTINE compute_drhocg_gpu
    !
#endif
END MODULE compute_drhocg_gpu_m 
!
#if defined(__CUDA)
!----------------------------------------------------------------------------
SUBROUTINE deriv_drhoc_gpu( ngl, gl_d, omega, tpiba2, mesh, r_d, rab_d, rhoc_d, &
                            drhocg_d )
  !--------------------------------------------------------------------------
  !! Calculates the Fourier transform of \(d\text{Rho}_c/dG\).
  !---cuda kernel version
  !
  USE cudafor
  USE kinds
  USE constants,  ONLY : pi, fpi
  USE compute_drhocg_gpu_m
  !
  IMPLICIT NONE
  !
  INTEGER :: ngl
  !! input: the number of g shell
  INTEGER :: mesh
  !! input: the number of radial mesh points
  REAL(DP), INTENT(IN), DEVICE :: gl_d(ngl)
  !! input: the number of G shells
  REAL(DP), INTENT(IN), DEVICE :: r_d(mesh)
  !! input: the radial mesh
  REAL(DP), INTENT(IN), DEVICE :: rab_d(mesh)
  !! input: the derivative of the radial mesh
  REAL(DP), INTENT(IN), DEVICE :: rhoc_d(mesh)
  !! input: the radial core charge
  REAL(DP), INTENT(IN) :: omega
  !! input: the volume of the unit cell
  REAL(DP), INTENT(IN) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP), INTENT(OUT), DEVICE :: drhocg_d(ngl)
  !! Fourier transform of d Rho_c/dG
  !
  ! ... local variables
  !
  REAL(DP) :: gx, gl1
  ! the modulus of g for a given shell
  INTEGER :: igl, blocks, igl0
  TYPE(dim3) :: threads
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !
  ! G=0 term
  !
  gl1 = gl_d(1)
  IF (gl1 < 1.0d-8) THEN
     drhocg_d(1) = 0.0_DP
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  ! G <> 0 term
  !
  threads = dim3(32,8,1)
  blocks = CEILING(REAL(ngl-igl0+1)/8)
  CALL compute_drhocg_gpu<<<blocks, threads>>>( ngl-igl0+1, tpiba2, omega, gl_d(igl0), r_d, &
                                 rhoc_d, rab_d, mesh, drhocg_d(igl0) )
  !
  RETURN
  !
END SUBROUTINE deriv_drhoc_gpu
!
!
#else
!
!
SUBROUTINE deriv_drhoc_gpu( ngl, gl_d, omega, tpiba2, mesh, r_d, rab_d, rhoc_d, &
                            drhocg_d )
  !--------------------------------------------------------------------------
  !! Calculates the Fourier transform of \(d\text{Rho}_c/dG\).
  !---cuf kernel loop version---
  !
  USE kinds
  USE constants,    ONLY: pi, fpi
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE simpsn_gpum,  ONLY: simpsn_gpu_dev
  USE device_fbuff_m,     ONLY: dev_buf
  !
  IMPLICIT NONE
  !
  INTEGER :: ngl
  !! input: the number of g shell
  INTEGER :: mesh
  !! input: the number of radial mesh points
  REAL(DP), INTENT(IN) :: gl_d(ngl)
  !! input: the number of G shells
  REAL(DP), INTENT(IN) :: r_d(mesh)
  !! input: the radial mesh
  REAL(DP), INTENT(IN) :: rab_d(mesh)
  !! input: the derivative of the radial mesh
  REAL(DP), INTENT(IN) :: rhoc_d(mesh)
  !! input: the radial core charge
  REAL(DP), INTENT(IN) :: omega
  !! input: the volume of the unit cell
  REAL(DP), INTENT(IN) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP), INTENT(OUT) :: drhocg_d(ngl)
  !! Fourier transform of d Rho_c/dG
  !
  ! ... local variables
  !
  REAL(DP) :: gx, gl1
  ! the modulus of g for a given shell
  REAL(DP), ALLOCATABLE :: aux_d(:,:)
  ! auxiliary memory for integration
  INTEGER :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  !
#if defined(__CUDA)
  attributes(DEVICE) :: gl_d, r_d, rab_d, rhoc_d, drhocg_d, aux_d
#endif
  !
  ! G=0 term
  !
  gl1 = gl_d(1)
  IF (gl1 < 1.0d-8) THEN
     drhocg_d(1) = 0.0_DP
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  ! G <> 0 term
  !
  ALLOCATE( aux_d(mesh,ngl) )
  !
  !$cuf kernel do (1)
  DO igl = igl0, ngl
     gx = SQRT( gl_d(igl) * tpiba2 )
     DO ir = 1, mesh
        aux_d(ir,igl) = r_d(ir)*rhoc_d(ir)*( r_d(ir) * COS(gx*r_d(ir)) /  &
                                      gx - SIN(gx*r_d(ir)) / gx**2 )
     ENDDO
     CALL simpsn_gpu_dev( mesh, aux_d(:,igl), rab_d, drhocg_d(igl) )
     drhocg_d(igl) = fpi / omega * drhocg_d(igl)
  ENDDO
  !
  DEALLOCATE( aux_d )
  !
  RETURN
  !
END SUBROUTINE deriv_drhoc_gpu
!
#endif
