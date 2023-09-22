!
! Copyright (C) 2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE vloc_mod
  !
  !! Variables and routines for local pseudopotential in numerical form
  !! Code moved to upflib and restructured by Paolo Giannozzi, 2023
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: vloc_of_g
  PUBLIC ::dvloc_of_g
  !
CONTAINS
!----------------------------------------------------------------------
SUBROUTINE vloc_of_g( nt, ngl, gl, tpiba2, modified_coulomb, omega, &
                vloc, ierr )
  !----------------------------------------------------------------------
  !! This routine computes the Fourier transform of the local
  !! part of an atomic pseudopotential, given in numerical form.
  !! A term erf(r)/r is subtracted in real space (thus making the
  !! function short-ranged) and added again in G space (for G<>0)
  !! The G=0 term contains \int (V_loc(r)+ Ze^2/r) 4pi r^2 dr.
  !! This is the "alpha" in the so-called "alpha Z" term of the energy.
  !! Atomic Ry units everywhere.
  !
  !
  USE atom,         ONLY : rgrid, msh
  USE uspp_param,   ONLY : upf
  USE m_gth,        ONLY : vloc_gth
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nt
  !! the index of type of pseudopotential
  INTEGER, INTENT(IN) :: ngl
  !! the number of shells of G vectors
  LOGICAL, INTENT(IN) :: modified_coulomb
  !! for ESM and 2D calculations - ierr = 1 is returned if not implemented 
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the (ordered!) moduli of g vectors for each shell 
  REAL(DP), INTENT(OUT) :: vloc(ngl)
  !! the fourier transform of the potential
  INTEGER, INTENT(OUT) :: ierr
  !! error code (ierr = 0 all good)
  !
  ! ... local variables
  !
  REAL(DP) :: fac, gx
  ! auxiliary variables
  REAL(DP), ALLOCATABLE :: aux(:), aux1(:)
  ! work space for integration
  INTEGER :: igl, igl0, ir
  ! igl : counter on g shells vectors
  ! igl0: position of first nonzero G
  ! ir  : counter on mesh points
  !
  IF ( modified_coulomb .AND. (upf(nt)%is_gth .OR. upf(nt)%tcoulombp) ) THEN
     ierr = 1
     RETURN
  END IF
  ierr = 0
  IF (gl (1) < eps8) THEN
     igl0 = 2
  ELSE
     igl0 = 1
  END IF
  IF ( upf(nt)%is_gth ) THEN
     ! special case: GTH pseudopotential
     CALL vloc_gth( nt, upf(nt)%zp, tpiba2, ngl, gl, omega, vloc )
  ELSE IF ( upf(nt)%tcoulombp ) THEN
     ! special case: pure Coulomb pseudopotential
     IF ( igl0 > 1 ) vloc(1) = 0.0_dp
     vloc (igl0:ngl) = - fpi * upf(nt)%zp*e2 / omega / tpiba2 / gl (igl0:ngl)
  ELSE
     ! normal case
     ALLOCATE (aux(msh(nt)), aux1(msh(nt)))
     !
     IF ( igl0 == 2 ) THEN
        !
        IF ( modified_coulomb ) THEN
           ! G=0 limit for the modified Coulomb potential (ESM, 2D)
           DO ir = 1, msh(nt)
              aux (ir) = rgrid(nt)%r(ir) * (rgrid(nt)%r(ir)*upf(nt)%vloc (ir) +&
                   upf(nt)%zp * e2 *  erf (rgrid(nt)%r(ir)) )
           END DO
        ELSE
           ! G=0 limit for Coulomb potential
           DO ir = 1, msh(nt)
              aux (ir) = rgrid(nt)%r(ir) * (rgrid(nt)%r(ir)*upf(nt)%vloc (ir) +&
                   upf(nt)%zp * e2 )
           END DO
        END IF
        !
        CALL simpson (msh(nt), aux, upf(nt)%rab, vloc(1))
        !
     END IF
     !   here the G<>0 terms with long range removed
     !
     ! G-independent term 
     DO ir = 1, msh(nt)
        aux (ir) = rgrid(nt)%r(ir)*upf(nt)%vloc (ir) + &
             upf(nt)%zp * e2 * erf (rgrid(nt)%r(ir))
     END DO
     !
     !    and here we perform the integral, after multiplying for the |G|
     !    dependent part
     !
     DO igl = igl0, ngl
        gx = sqrt ( gl(igl)*tpiba2 )
        DO ir = 1, msh(nt)
           aux1(ir) = aux (ir) * sin ( gx*upf(nt)%r(ir) ) / gx
        END DO
        CALL simpson (msh(nt), aux1, upf(nt)%rab, vloc(igl))
     END DO
     DEALLOCATE (aux, aux1)
     !
     IF ( .not. modified_coulomb ) THEN
        fac = upf(nt)%zp * e2 / tpiba2
        DO igl = igl0, ngl
           !
           !   here we re-add the analytic fourier transform of the erf function
           !
           vloc(igl) = vloc(igl) - fac * exp (-gl (igl)*tpiba2*0.25d0) / gl (igl)
        END DO
     END IF
     vloc (:) = vloc(:) * fpi / omega
  END IF

END SUBROUTINE vloc_of_g
!
!----------------------------------------------------------------------
SUBROUTINE dvloc_of_g( mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, gl, &
                       is_coulomb, modified_coulomb, omega, dvloc )
  !----------------------------------------------------------------------
  !! This routine computes:
  !! \[ \text{dvloc} = D\text{Vloc}(g^2)/Dg^2 = (1/2g)\ D\text{Vloc}(g)/Dg
  !! \]
  !
  INTEGER, INTENT(IN) :: ngl
  !! the number of shell of G vectors
  INTEGER, INTENT(IN) :: mesh
  !! max number of mesh points
  INTEGER, INTENT(IN) :: msh
  !! number of mesh points for radial integration
  LOGICAL, INTENT(IN) :: is_coulomb
  !! for pure Coulomb pseudopotentials
  LOGICAL, INTENT(IN) :: modified_coulomb
  !! for ESM and 2D calculations
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
  IF ( is_coulomb ) THEN
     !$acc kernels
     dvloc(igl0:ngl) = fpi*zp*e2 / omega / (tpiba2*gl(igl0:ngl))**2
     !$acc end kernels
     RETURN
  END IF
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
     IF ( .NOT. modified_coulomb ) THEN
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
END MODULE vloc_mod
