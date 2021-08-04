!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE rgen ( dtau, rmax, mxr, at, bg, r, r2, nrm)
  !-----------------------------------------------------------------------
  !! Generates neighbours shells (cartesian, in units of lattice parameter)
  !! with length < rmax, and returns them in order of increasing length.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: mxr
  !! maximum number of vectors
  INTEGER, INTENT(out) :: nrm
  !! the number of vectors with r^2 < rmax^2
  REAL(DP), INTENT(in) :: at(3,3)
  !! lattice vectors: \(a_1=at(:,1)\), \(a_2=at(:,2)\), \(a_3=at(:,3)\)
  REAL(DP), INTENT(in) :: bg(3,3)
  !! the reciprocal lattice vectors
  REAL(DP), INTENT(in) :: dtau(3)
  !! the difference of the atomic postions: tau\_s - tau\_s'
  REAL(DP), INTENT(in) :: rmax
  !! the maximum radius to consider
  REAL(DP), INTENT(out) :: r(3,mxr)
  !! \(r(:) = i*a_1(:) + j*a_2(:) + k*a_3(:)-\text{dtau}(:)\) where \(a_1\),
  !! \(a_2\), \(a_3\) are primitive lattice vectors.
  REAL(DP), INTENT(out) :: r2(mxr)
  !! r2 = r^2
  !
  ! ... local variables
  !
  INTEGER, ALLOCATABLE :: irr (:)
  INTEGER ::  nm1, nm2, nm3, i, j, k, ipol, ir, indsw, iswap
  real(DP) :: ds(3), dtau0(3)
  real(DP) :: t (3), tt, swap
  real(DP), EXTERNAL :: dnrm2
  !
  !
  nrm = 0
  IF (rmax==0.d0) RETURN

  ! bring dtau into the unit cell centered on the origin - prevents trouble
  ! if atomic positions are not centered around the origin but displaced
  ! far away (remember that translational invariance allows this!)
  !
  ds(:) = matmul( dtau(:), bg(:,:) )
  ds(:) = ds(:) - anint(ds(:))
  dtau0(:) = matmul( at(:,:), ds(:) )
  !
  ALLOCATE (irr( mxr))
  !
  ! these are estimates of the maximum values of needed integer indices
  !
  nm1 = int (dnrm2 (3, bg (1, 1), 1) * rmax) + 2
  nm2 = int (dnrm2 (3, bg (1, 2), 1) * rmax) + 2
  nm3 = int (dnrm2 (3, bg (1, 3), 1) * rmax) + 2
  !
  DO i = -nm1, nm1
     DO j = -nm2, nm2
        DO k = -nm3, nm3
           tt = 0.d0
           DO ipol = 1, 3
              t (ipol) = i*at (ipol, 1) + j*at (ipol, 2) + k*at (ipol, 3) &
                       - dtau0(ipol)
              tt = tt + t (ipol) * t (ipol)
           ENDDO
           IF (tt<=rmax**2.and.abs (tt) >1.d-10) THEN
              nrm = nrm + 1
              IF (nrm>mxr) CALL errore ('rgen', 'too many r-vectors', nrm)
              DO ipol = 1, 3
                 r (ipol, nrm) = t (ipol)
              ENDDO
              r2 (nrm) = tt
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  !   reorder the vectors in order of increasing magnitude
  !
  !   initialize the index inside sorting routine
  !
  irr (1) = 0
  IF (nrm>1) CALL hpsort (nrm, r2, irr)
  DO ir = 1, nrm - 1
20   indsw = irr (ir)
     IF (indsw/=ir) THEN
        DO ipol = 1, 3
           swap = r (ipol, indsw)
           r (ipol, indsw) = r (ipol, irr (indsw) )
           r (ipol, irr (indsw) ) = swap
        ENDDO
        iswap = irr (ir)
        irr (ir) = irr (indsw)
        irr (indsw) = iswap
        GOTO 20
     ENDIF

  ENDDO
  DEALLOCATE(irr)
  !
  RETURN
END SUBROUTINE rgen

