!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE compute_dvloc_gpu_m
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
!---------------------------------------------------------------------
ATTRIBUTES(DEVICE) FUNCTION qe_erf_d( x )
  !---------------------------------------------------------------------
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  USE kinds,   ONLY: DP
  implicit none
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: x2, p1 (4), q1 (4)
  REAL(DP) :: qe_erf_d
  DATA p1 / 2.426679552305318E2_DP, 2.197926161829415E1_DP, &
            6.996383488619136_DP,  -3.560984370181538E-2_DP /
  DATA q1 / 2.150588758698612E2_DP, 9.116490540451490E1_DP, &
            1.508279763040779E1_DP, 1.000000000000000_DP /
  !
  IF (ABS(x) > 6.0_DP) THEN
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1
     !
     qe_erf_d = SIGN(1.0_DP, x)
  ELSE
     IF (ABS(x)  <= 0.47_DP) THEN
        x2 = x**2
        qe_erf_d=x *(p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 (4) ) ) ) &
                / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 (4) ) ) )
     ELSE
        qe_erf_d = 1.0_DP - qe_erfc_d(x)
     ENDIF
  ENDIF
  !
  RETURN
  !
END FUNCTION qe_erf_d
!
!-----------------------------------------------------------------------------
ATTRIBUTES(DEVICE) FUNCTION qe_erfc_d( x )
  !---------------------------------------------------------------------
  ! erfc(x) = 1-erf(x)  - See comments in erf
  !
  USE kinds,   ONLY: DP
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: qe_erfc_d
  REAL(DP) :: ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1
  !
  DATA p2 / 3.004592610201616E2_DP,  4.519189537118719E2_DP, &
            3.393208167343437E2_DP,  1.529892850469404E2_DP, &
            4.316222722205674E1_DP,  7.211758250883094_DP,   &
            5.641955174789740E-1_DP,-1.368648573827167E-7_DP /
  DATA q2 / 3.004592609569833E2_DP,  7.909509253278980E2_DP, &
            9.313540948506096E2_DP,  6.389802644656312E2_DP, &
            2.775854447439876E2_DP,  7.700015293522947E1_DP, &
            1.278272731962942E1_DP,  1.000000000000000_DP /
  DATA p3 /-2.996107077035422E-3_DP,-4.947309106232507E-2_DP, &
           -2.269565935396869E-1_DP,-2.786613086096478E-1_DP, &
           -2.231924597341847E-2_DP /
  DATA q3 / 1.062092305284679E-2_DP, 1.913089261078298E-1_DP, &
            1.051675107067932_DP,    1.987332018171353_DP,    &
            1.000000000000000_DP /

  DATA pim1 / 0.56418958354775629_DP /
  !        ( pim1= sqrt(1/pi) )
  ax = ABS(x)
  IF (ax > 26.0_DP) THEN
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     qe_erfc_d = 0.0_DP
  ELSEIF (ax > 4.0_DP) THEN
     x2 = x**2
     xm2 = (1.0_DP / ax) **2
     qe_erfc_d = (1.0_DP / ax) * EXP( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  ELSEIF (ax > 0.47_DP) THEN
     x2 = x**2
     qe_erfc_d = EXP( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3)   &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  ELSE
     qe_erfc_d = 1.0_DP - qe_erf_d(ax)
  ENDIF
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  IF (x < 0.0_DP) qe_erfc_d = 2.0_DP - qe_erfc_d
  !
  RETURN
  !
END FUNCTION qe_erfc_d
!
!
END MODULE compute_dvloc_gpu_m



!----------------------------------------------------------------------
SUBROUTINE dvloc_of_g_gpu( mesh, msh, rab_d, r_d, vloc_at_d, zp, tpiba2, &
                           ngl, gl_d, omega, dvloc_d, igl0 )
  !----------------------------------------------------------------------
  !
  ! dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  !
  USE cudafor
  USE kinds
  USE constants,          ONLY : pi, fpi, e2, eps8
  USE Coul_cut_2D,        ONLY : do_cutoff_2D
  USE esm,                ONLY : do_comp_esm, esm_bc
  USE compute_dvloc_gpu_m
  !
  IMPLICIT NONE
  !
  !    first the dummy variables
  !
  INTEGER, INTENT(IN) :: mesh, msh, ngl, igl0
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
  REAL(DP), ALLOCATABLE, DEVICE :: aux1_d(:)
  !
  INTEGER :: i, igl, blocks
  TYPE(dim3) :: threads
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shell with g != 0
  !
  ! the  G=0 component is not computed
  IF (igl0==2) dvloc_d(1) = 0.d0
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
  !$cuf kernel do(1)<<<*,*>>>
  DO i = 1, msh
     aux1_d(i) = r_d(i) * vloc_at_d(i) + zp * e2 * qe_erf_d(r_d(i))
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
  DEALLOCATE( aux1_d )
  !
  !
  RETURN
  !
END SUBROUTINE dvloc_of_g_gpu
!
!
!----------------------------------------------------------------------
SUBROUTINE dvloc_coul_gpu( zp, tpiba2, ngl, gl_d, omega, dvloc_d, igl0 )
  !----------------------------------------------------------------------
  !! Fourier transform of the Coulomb potential - For all-electron
  !! calculations, in specific cases only, for testing purposes.
  !
  USE kinds
  USE constants,    ONLY : fpi, e2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngl, igl0
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
  INTEGER :: j
  !
#if defined(__CUDA)
  attributes(DEVICE) ::  dvloc_d, gl_d
#endif
  !
  ! the  G=0 component is 0
  IF (igl0==2) dvloc_d(1) = 0._DP
  !
  !$cuf kernel do(1) <<<*,*>>>
  DO j = igl0, ngl
     dvloc_d(j) = fpi * zp * e2 / omega / ( tpiba2 * gl_d(j) )**2
  ENDDO
  !
  RETURN
  !
END SUBROUTINE dvloc_coul_gpu
