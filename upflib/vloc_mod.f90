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
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: vloc_of_g, vloc_coul
  SAVE
CONTAINS
!----------------------------------------------------------------------
SUBROUTINE vloc_of_g( nt, ngl, gl, tpiba2, modified_coulomb, omega, vloc )
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
  USE upf_const,    ONLY : fpi, e2, eps8
  USE atom,         ONLY : rgrid, msh
  USE uspp_param,   ONLY : upf
  ! USE m_gth,        ONLY : vloc_gth
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nt
  !! the index of type of pseudpotential
  INTEGER, INTENT(IN) :: ngl
  !! the number of shells of G vectors
  LOGICAL, INTENT(IN) :: modified_coulomb
  !! for ESM and 2D calculations
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the (ordered!) moduli of g vectors for each shell 
  REAL(DP), INTENT(OUT) :: vloc(ngl)
  !! the fourier transform of the potential
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
  ALLOCATE (aux(msh(nt)), aux1(msh(nt)))
  !
  IF (gl (1) < eps8) THEN
     !
     IF ( modified_coulomb ) THEN
        ! G=0 limit for the modified Coulomb potential (ESM, 2D)
        DO ir = 1, msh(nt)
           aux (ir) = rgrid(nt)%r(ir) * ( rgrid(nt)%r(ir)*upf(nt)%vloc (ir) + &
                      upf(nt)%zp * e2 *  erf (rgrid(nt)%r(ir)) )
        END DO
     ELSE
        ! G=0 limit for Coulomb potential
        DO ir = 1, msh(nt)
           aux (ir) = rgrid(nt)%r(ir) * ( rgrid(nt)%r(ir)*upf(nt)%vloc (ir) + &
                      upf(nt)%zp * e2 )
        END DO
     END IF
     !
     CALL simpson (msh(nt), aux, upf(nt)%rab, vloc(1))
     !
     igl0 = 2
     !
  ELSE
     !
     igl0 = 1
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
  DEALLOCATE (aux, aux1)

END SUBROUTINE vloc_of_g
!
!----------------------------------------------------------------------
SUBROUTINE vloc_coul( zp, tpiba2, ngl, gl, omega, vloc )
  !----------------------------------------------------------------------
  !! Fourier transform of the Coulomb potential - For all-electron
  !! calculations, in specific cases only, for testing purposes.
  !
  USE upf_const,    ONLY : fpi, e2, eps8
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngl
  !! the number of shells of G vectors
  REAL(DP), INTENT(IN) :: zp
  !! valence pseudocharge
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the moduli of g vectors for each shell
  !
  ! ... local variables
  !
  REAL(DP), INTENT(OUT) :: vloc(ngl)
  ! the fourier transform of the potential
  INTEGER :: igl0
  !
  IF (gl (1) < eps8) THEN
     igl0 = 2
     vloc(1) = 0.0_dp
  ELSE
     igl0 = 1
  END IF

  vloc (igl0:ngl) = - fpi * zp *e2 / omega / tpiba2 / gl (igl0:ngl)

END SUBROUTINE vloc_coul
!
END MODULE vloc_mod
