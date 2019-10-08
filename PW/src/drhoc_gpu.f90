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
subroutine drhoc_gpu (ngl, gl_d, omega, tpiba2, mesh, r_d, rab_d, rhoc_d, rhocg_d)
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : pi, fpi, eps14
  use cudafor
  use compute_rhocg_gpu_m
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ngl, mesh
  ! input: the number of g shell
  ! input: the number of radial mesh points

  real(DP) ::  omega, tpiba2
  real(DP), device :: gl_d(ngl), r_d(mesh), rab_d(mesh), rhoc_d(mesh), rhocg_d(ngl)
  ! input: the number of G shells
  ! input: the radial mesh
  ! input: the derivative of the radial mesh
  ! input: the radial core charge
  ! input: the volume of the unit cell
  ! input: 2 times pi / alat
  ! output: the fourier transform of the core charge
  !
  !     here the local variables
  !
  real(DP) :: rhocg1, gl1
  real(DP), device :: func_d(3)
  ! the modulus of g for a given shell
  ! the fourier transform
  real(DP), allocatable, device ::  aux_d (:)
  ! auxiliary memory for integration

  integer :: ir, igl, igl0, ir0, blocks
  type(dim3) :: threads
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  CALL start_clock( 'drhoc' )

  allocate (aux_d( mesh))     
  !
  ! G=0 term
  !
  gl1 = gl_d(1)
  if (gl1 < 1.0d-8) then
     !$cuf kernel do(1) <<<*, *>>>
     do ir = 1, mesh
        aux_d (ir) = r_d (ir) **2 * rhoc_d (ir)
     enddo

     call simpson_gpu (mesh, aux_d, rab_d, rhocg1)
     rhocg1 = fpi * rhocg1 / omega
     rhocg_d(1) = rhocg1 ! copy result to device

     igl0 = 2
  else
     igl0 = 1
  endif
  !
  ! G <> 0 term
  !
  threads = dim3(32, 8, 1)
  blocks = ceiling(real(ngl - igl0 + 1)/8)
  call compute_rhocg_gpu<<<blocks, threads>>>(ngl - igl0 + 1, tpiba2, omega, gl_d(igl0), r_d, rhoc_d, & 
                                              rab_d, mesh, rhocg_d(igl0))

  deallocate(aux_d)
  !
  CALL stop_clock( 'drhoc' )
  return
end subroutine drhoc_gpu

#else 
SUBROUTINE drhoc_gpu()
END SUBROUTINE drhoc_gpu 
#endif 
