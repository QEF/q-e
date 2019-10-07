!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#if defined(__CUDA) 
MODULE compute_rhocg_gpu_m
  USE cudafor
  IMPLICIT NONE 
  TYPE(dim3) :: rhocg_threads = dim3 (32,8,1) 
  PUBLIC     :: rhocg_threads 
  CONTAINS 
   ATTRIBUTES( global) SUBROUTINE compute_rhocg_gpu( n, tpiba2, omega, gl, r, rhoc, rab, mesh, rhocg)  
      !! implements in a cuda kernel the loop: 
      !! DO igl = igl0, ngl
      !!    gx = SQRT(gl(igl) * tpiba2)
      !!    CALL sph_bes( mesh, r, gx, 0, aux )
      !!    DO ir = 1, mesh
      !!       aux(ir) = r(ir)**2 * rhoc(ir) * aux(ir)
      !!    ENDDO
      !!    CALL simpson( mesh, aux, rab, rhocg1 )
      !!    rhocg(igl) = fpi * rhocg1 / omega
      !! ENDDO
      USE cudafor
      USE kinds,     ONLY: DP 
      USE constants, ONLY: pi, fpi, eps14  
      IMPLICIT NONE 
      !
      INTEGER, VALUE                  ::  n, mesh 
      REAL(DP),VALUE                  ::  tpiba2, omega
      REAL(DP),DEVICE,INTENT(IN)      ::  gl(n), r(mesh), rhoc(mesh), rab(mesh)
      REAL(DP),DEVICE,INTENT(OUT)     ::  rhocg(n) 
      !   
      INTEGER                         :: tx, ty, igl, ir
      REAL(DP)                        :: mysum, val, gx, x
      ! 
      tx = threadIdx%x
      ty = threadIdx%y 
      ! 
      igl =  (blockIdx%x - 1) * blockDim%y + ty
      ! 
      IF (igl > n ) RETURN 
      ! 
      gx = SQRT(gl(igl) * tpiba2)
      mysum = 0 
      !
      ! 
      IF (ABS(gx) < eps14) THEN 
         DO ir = tx, mesh, blockDim%x
            val = r(ir) * r(ir) * rhoc(ir) * rab(ir)
              IF (ir == 1 .OR. ir == mesh) THEN
            mysum = mysum + val
          ELSE IF (mod(ir,2)) THEN
            mysum = mysum + 2.d0*val
          ELSE
            mysum = mysum + 4.d0*val
          ENDIF
        END DO
      ELSE
        DO ir = tx, mesh, blockDim%x 
          x = gx * r(ir)
          IF  (ABS(x) > 0.05_DP) THEN  !from check on xseries = 0.05_DP
            val = SIN (x) / (x) * r(ir) * r(ir) * rhoc(ir)
          ELSE 
            val = ( 1.0_dp - x*x/6.0_dp * &
                  ( 1.0_dp - x*x/20.0_dp * &
                  ( 1.0_dp - x*x/42.0_dp * &
                  ( 1.0_dp - x*x/72.0_dp ) ) ) ) * r(ir) * r(ir) * rhoc(ir)
          ENDIF

          val = val * rab(ir)

          IF (ir == 1 .or. ir == mesh) THEN
            mysum = mysum + val
          ELSE IF (mod(ir,2)) THEN
            mysum = mysum + 2.d0*val
          ELSE
            mysum = mysum + 4.d0*val
          END IF
        END DO
      END IF
      

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

      IF (tx == 1) THEN 
        rhocg(igl) = fpi * mysum / (3.d0 * omega)
      ENDIF

    END SUBROUTINE compute_rhocg_gpu  

END MODULE compute_rhocg_gpu_m 
!-----------------------------------------------------------------------
SUBROUTINE drhoc_gpu( ngl, gl_d, omega, tpiba2, mesh, r_d, rab_d, rhoc_d, rhocg_d )
  !-----------------------------------------------------------------------
  !! Calculates the Fourier transform of the core charge.
  !
  USE kinds
  USE constants,   ONLY: pi, fpi, eps14 
  USE cudafor  
  USE compute_rhocg_gpu_m  
  !
  IMPLICIT NONE
  !
  INTEGER :: ngl
  !! input: the number of g shell
  INTEGER :: mesh
  !! input: the number of radial mesh points
  REAL(DP),DEVICE :: gl_d(ngl)
  !! input: the number of G shells
  REAL(DP),DEVICE :: r_d(mesh)
  !! input: the radial mesh
  REAL(DP),DEVICE  :: rab_d(mesh)
  !! input: the derivative of the radial mesh
  REAL(DP),DEVICE :: rhoc_d(mesh)
  !! input: the radial core charge
  REAL(DP) :: omega
  !! input: the volume of the unit cell
  REAL(DP) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP), DEVICE :: rhocg_d(ngl)
  !! output: the Fourier transform of the core charge
  !
  ! ... local variables
  !
  REAL(DP)         :: gl1, rhocg1
  REAL(DP),DEVICE  :: func_d(3) 
  ! the modulus of g for a given shell
  ! the Fourier transform
  REAL(DP), ALLOCATABLE,DEVICE  ::  aux_d(:)
  ! auxiliary memory for integration
  INTEGER :: ir, igl, igl0, ir0, blocks 
  TYPE(dim3) :: threads 
  INTEGER    :: n_
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !
  !
  CALL start_clock('drhoc') 
  !
  ALLOCATE( aux_d(mesh) )
  !
  ! G=0 term
  !
  gl1 = gl_d(1) 
  IF ( gl1 < 1.0d-8 ) THEN
     !$cuf kernel do(1) <<<*,*>>> 
     DO ir = 1, mesh
        aux_d(ir) = r_d(ir)**2 * rhoc_d(ir)
     END DO
     CALL simpson_gpu( mesh, aux_d, rab_d, rhocg1 )
     rhocg1 = fpi * rhocg1 / omega
     rhocg_d(1) = rhocg1 
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  ! G <> 0 term
  !
  blocks = ceiling(REAL((ngl-igl0+1)/8))
  threads = rhocg_threads 
  n_ = ngl-igl0+1
  !
  
  !
  CALL compute_rhocg_gpu<<<blocks, threads>>>(n_, tpiba2, omega, gl_d(igl0), r_d, rhoc_d, &
                                             rab_d, mesh, rhocg_d(igl0))  
     
  DEALLOCATE(aux_d) 
     ! 
     CALL stop_clock('drhoc') 
  !
  RETURN
  !
END SUBROUTINE drhoc_gpu 
#else 
SUBROUTINE drhoc_gpu()
END SUBROUTINE drhoc_gpu 
#endif 
