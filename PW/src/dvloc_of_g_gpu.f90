!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE compute_dvloc_gpum
#if defined(__CUDA)
  USE cudafor
  IMPLICIT NONE
  TYPE(dim3) :: dvloc_threads = dim3(32,8,1)
  PUBLIC     :: dvloc_threads
  CONTAINS
  ATTRIBUTES( global ) SUBROUTINE compute_dvloc_gpu( mesh, n, tpiba2, omega, zp, gl, aux1, r, rab ,dvloc, do_cutoff_2D, do_comp_esm, esm_bc_l )
  !
  !
  USE cudafor
  USE kinds
  USE constants,    ONLY : pi, fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  INTEGER,  VALUE :: n, mesh
  REAL(DP), VALUE :: zp, tpiba2, omega
  LOGICAL,  VALUE :: do_cutoff_2D, do_comp_esm, esm_bc_l
  REAL(DP), DEVICE, INTENT(IN) :: rab(mesh), r(mesh), gl(n)
  REAL(DP), DEVICE, INTENT(IN) :: aux1(mesh)
  !
  REAL(DP), DEVICE, INTENT(OUT) ::  dvloc(n)
  ! the Fourier transform dVloc/dG
  !
  REAL(DP) :: val, mysum
  REAL(DP) :: g2a, gx
  !
  INTEGER :: tx, ty, ir, igl
  !
  tx = threadIdx%x
  ty = threadIdx%y
  ! 
  igl =  (blockIdx%x - 1) * blockDim%y + ty
  ! 
  IF (igl > n ) RETURN
  !
  gx = SQRT(gl(igl) * tpiba2)
  mysum = 0.d0
  !
  DO ir = tx, mesh, blockDim%x
    !
    val = aux1(ir) * ( r(ir) * COS(gx*r(ir)) / gx - SIN(gx*r(ir)) / gx**2 )
    val = val * rab(ir)
    !
    IF (ir == 1 .OR. ir == mesh) THEN
      mysum = mysum + val
    ELSEIF ( MOD(ir,2) ) THEN
      mysum = mysum + 2.d0*val
    ELSE
      mysum = mysum + 4.d0*val
    ENDIF
    !
  ENDDO
  !
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
    !
    mysum = mysum * fpi / omega / 2.0d0 / gx / 3.d0
    !
    IF ( (.NOT. do_comp_esm) .OR. esm_bc_l ) THEN
      IF (.NOT. do_cutoff_2D) THEN
        g2a  = gl(igl) * tpiba2 / 4.d0
        mysum = mysum + fpi / omega * zp * e2 * EXP(-g2a) * (g2a + &
                1.d0) / (gl(igl)*tpiba2)**2
      ENDIF
    ENDIF
    !
    dvloc(igl) = mysum
    !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE compute_dvloc_gpu
!
#endif
END MODULE compute_dvloc_gpum


!----------------------------------------------------------------------
SUBROUTINE dvloc_of_g_gpu( mesh, msh, rab_d, r_d, vloc_at_d, zp, tpiba2, &
                           ngl, gl_d, omega, dvloc_d ) !, igl0 )
  !----------------------------------------------------------------------
  !
  ! dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds
  USE constants,          ONLY : pi, fpi, e2, eps8
  USE Coul_cut_2D,        ONLY : do_cutoff_2D
  USE esm,                ONLY : do_comp_esm, esm_bc
#if defined(__CUDA)
  USE compute_dvloc_gpum, ONLY : compute_dvloc_gpu
#endif
  !
  IMPLICIT NONE
  !
  !    first the dummy variables
  !
  INTEGER, INTENT(IN) :: mesh, msh, ngl !, igl0
  ! the number of shell of G vectors
  ! max number of mesh points
  ! number of mesh points for radial integration
  !
  REAL(DP), INTENT(IN) :: rab_d(mesh), r_d(mesh), vloc_at_d(mesh), gl_d(ngl)
  REAL(DP), INTENT(IN) :: zp, tpiba2, omega
  !
  REAL(DP), INTENT(OUT) ::  dvloc_d(ngl)
  ! the Fourier transform dVloc/dG
  !
#if defined(__CUDA)
  attributes(DEVICE) ::  dvloc_d, rab_d, r_d , gl_d, vloc_at_d
#endif
  !
  ! valence pseudocharge
  ! the derivative of the radial grid
  ! the radial grid
  ! the pseudo on the radial grid
  ! 2 pi / alat
  ! the volume of the unit cell
  ! the moduli of g vectors for each s
  !
  LOGICAL :: esm_bc_lg
  !
  REAL(DP), ALLOCATABLE :: aux1_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux1_d
  !
  TYPE(dim3) :: threads
#endif
  !
  REAL(DP) :: gl1
  INTEGER :: i, igl, blocks, igl0
  !
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shell with g != 0
  !
  ! the  G=0 component is not computed
  gl1 = gl_d(1)
  IF (gl1 < eps8) THEN
     dvloc_d(1) = 0.0_DP
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  ! Pseudopotentials in numerical form (Vloc contains the local part)
  ! In order to perform the Fourier transform, a term erf(r)/r is
  ! subtracted in real space and added again in G space
  !
  !
  ALLOCATE( aux1_d(mesh) )
  !
  !   This is the part of the integrand function
  !   indipendent of |G| in real space
  !
#if defined(__CUDA)
  !$cuf kernel do(1)<<<*,*>>>
  DO i = 1, msh
     aux1_d(i) = r_d(i) * vloc_at_d(i) + zp * e2 * erf(r_d(i))
  ENDDO
  !
  esm_bc_lg = ( esm_bc == 'pbe' )
  !
  threads = dim3(32,8,1)
  blocks = CEILING(REAL(ngl-igl0+1)/8)
  CALL compute_dvloc_gpu<<<blocks, threads>>>( msh, ngl-igl0+1, tpiba2, omega, zp, gl_d(igl0), &
                                 aux1_d, r_d, rab_d, dvloc_d(igl0), do_cutoff_2D, &
                                 do_comp_esm, esm_bc_lg )
  !
#else
  CALL errore( 'dvloc_of_g_gpu' , 'GPU version of dvloc_of_g called but not compiled.', 0 )
#endif
  DEALLOCATE( aux1_d )
  !
  !
  RETURN
  !
END SUBROUTINE dvloc_of_g_gpu
!
!
!----------------------------------------------------------------------
SUBROUTINE dvloc_coul_gpu( zp, tpiba2, ngl, gl_d, omega, dvloc_d )
  !----------------------------------------------------------------------
  !! Fourier transform of the Coulomb potential - For all-electron
  !! calculations, in specific cases only, for testing purposes.
  !
  USE kinds
  USE constants,    ONLY : fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngl
  ! the number of shell of G vectors
  ! first shell with g != 0
  REAL(DP), INTENT(IN) :: zp, tpiba2, omega
  ! valence pseudocharge
  ! 2 pi / alat
  ! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl_d(ngl)
  ! the moduli of g vectors for each s
  REAL(DP), INTENT(OUT) :: dvloc_d(ngl)
  ! fourier transform: dvloc = D Vloc (g^2) / D g^2 = 4pi e^2/omegai /G^4
  INTEGER :: j, igl0
  REAL(DP) :: gl1
  !
#if defined(__CUDA)
  attributes(DEVICE) ::  dvloc_d, gl_d
#endif
  !
  ! the  G=0 component is 0
  gl1 = gl_d(1)
  IF (gl1 < eps8) then
     dvloc_d(1) = 0.0_DP
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  !$cuf kernel do(1) <<<*,*>>>
  DO j = igl0, ngl
     dvloc_d(j) = fpi * zp * e2 / omega / ( tpiba2 * gl_d(j) )**2
  ENDDO
  !
  RETURN
  !
END SUBROUTINE dvloc_coul_gpu
