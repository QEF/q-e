!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dvloc_of_g( mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, gl, &
                       omega, dvloc )
  !----------------------------------------------------------------------
  !! This routine gives:
  !! \[ \text{dvloc} = D\text{Vloc}(g^2)/Dg^2 = (1/2g)\ D\text{Vloc}(g)/Dg
  !! \]
  !
  USE kinds
  USE constants,   ONLY: pi, fpi, e2, eps8
  USE Coul_cut_2D, ONLY: do_cutoff_2D
  USE esm,         ONLY: do_comp_esm, esm_bc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngl
  !! the number of shell of G vectors
  INTEGER, INTENT(IN) :: mesh
  !! max number of mesh points
  INTEGER, INTENT(IN) :: msh
  !! number of mesh points for radial integration
  REAL(DP), INTENT(IN) :: zp
  !! valence pseudocharge
  REAL(DP), INTENT(IN) :: rab(mesh)
  !! the derivative of the radial grid
  REAL(DP), INTENT(IN) :: r(mesh)
  !! the radial grid
  REAL(DP), INTENT(IN) :: vloc_at(mesh)
  !! the pseudo on the radial grid 
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the moduli of g vectors for each s
  REAL(DP), INTENT(OUT) ::  dvloc(ngl)
  !! the Fourier transform dVloc/dG
  !
  ! ... local variables
  !
  REAL(DP) :: vlcp, g2a, gx, vlcp_0, vlcp_1
  REAL(DP), ALLOCATABLE :: aux(:,:), aux1(:)
  INTEGER :: i, igl, igl0
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shell with g != 0
  REAL(DP), PARAMETER :: r12=1.0d0/3.0d0 
  !
  !$acc data present( dvloc, gl )
  !
  ! the  G=0 component is not computed
  IF (gl(1) < eps8) THEN
     !$acc kernels
     dvloc(1) = 0.0d0
     !$acc end kernels
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  ! Pseudopotentials in numerical form (Vloc contains the local part)
  ! In order to perform the Fourier transform, a term erf(r)/r is
  ! subtracted in real space and added again in G space
  !
  ALLOCATE( aux1(mesh) )
  !
  ! This is the part of the integrand function
  ! indipendent of |G| in real space
  !
  ALLOCATE( aux(mesh,ngl) )
  !
  !$acc data copyin(r,rab) create(aux1,aux)
  !
  !$acc parallel loop copyin(vloc_at)
  DO i = 1, msh
     aux1(i) = r(i)*vloc_at(i) + zp*e2*ERF(r(i))
  ENDDO
  !
#if defined(_OPENACC)
!$acc parallel loop gang present(aux,aux1,rab,r,dvloc)
#else
!$omp parallel private( gx, vlcp, vlcp_1, vlcp_0, g2a )
!$omp do
#endif
  DO igl = igl0, ngl
     !
     gx = SQRT(gl(igl)*tpiba2)
     !
     !    and here we perform the integral, after multiplying for the |G|
     !    dependent  part
     !
     ! DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
     !
     !$acc loop seq
     DO i = 1, msh
       aux(i,igl) = aux1(i)*(r(i)*COS(gx*r(i))/gx - SIN(gx*r(i))/gx**2)
     ENDDO
     !
     !----Simpson int.---
     vlcp_0 = 0.0d0    
     !$acc loop seq reduction(+:vlcp_0)
     DO i = 2, msh-1,  2
       vlcp_0 = vlcp_0 + ( aux(i-1,igl)*rab(i-1) + 4.0d0*aux(i,igl)*rab(i) + &
                           aux(i+1,igl)*rab(i+1) )*r12
     ENDDO
     !------
     vlcp_1 = vlcp_0 * fpi / omega / 2.0d0 / gx
     !
     ! DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
     !vlcp = fpi / omega / 2.0d0 / gx * vlcp
     !
     ! for ESM stress
     ! In ESM, vloc and dvloc have only short term.
     IF ( (( .NOT. do_comp_esm ) .OR. ( esm_bc .EQ. 'pbc' )) .AND. &
             .NOT.do_cutoff_2D ) THEN
        ! subtract the long-range term
        ! 2D cutoff: do not re-add LR part here re-added later in stres_loc)
        g2a = gl(igl) * tpiba2 / 4.d0
        vlcp = vlcp_1 + fpi / omega * zp * e2 * EXP(-g2a) * (g2a + 1.d0) / &
                          (gl(igl)*tpiba2)**2
     ELSE
        vlcp = vlcp_1
     ENDIF
     dvloc(igl) = vlcp
  ENDDO
#if !defined(_OPENACC)
!$omp end do nowait
!$omp end parallel
#else
!$acc end data
!$acc end data
#endif
  !
  DEALLOCATE( aux )
  DEALLOCATE( aux1 )
  !
  RETURN
  !
END SUBROUTINE dvloc_of_g
!
!----------------------------------------------------------------------
SUBROUTINE dvloc_coul( zp, tpiba2, ngl, gl, omega, dvloc )
  !----------------------------------------------------------------------
  !! Fourier transform of the Coulomb potential - For all-electron
  !! calculations, in specific cases only, for testing purposes.
  !
  USE kinds
  USE constants, ONLY : fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngl
  !! the number of shell of G vectors
  REAL(DP), INTENT(IN) :: zp
  !! valence pseudocharge
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the moduli of g vectors for each s
  REAL(DP), INTENT(OUT) :: dvloc(ngl)
  !! Fourier transform:  
  !! dvloc = D Vloc (g^2) / D g^2 = 4pi e^2/omegai /G^4
  !
  INTEGER :: igl0
  ! first shell with g != 0
  !
  !$acc data present( dvloc, gl )
  !
  ! the  G=0 component is 0
  IF (gl(1) < eps8) THEN
     !$acc kernels
     dvloc(1) = 0.0d0
     !$acc end kernels
     igl0 = 2
  ELSE
     igl0 = 1
  ENDIF
  !
  !$acc kernels
  dvloc(igl0:ngl) = fpi*zp*e2 / omega / (tpiba2*gl(igl0:ngl))**2
  !$acc end kernels
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE dvloc_coul

