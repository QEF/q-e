!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE ktetra
  !
  ! ... Variables used by the tetrahedron method
  ! ... Three versions are implemented: Linear, Optimized, Bloechl
  ! ... Linear and Optimized tetrahedra contributed by Mitsuaki Kawamura,
  ! ... University of Tokyo
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  INTEGER:: &
       tetra_type = 0  ! 0 for Bloechl's correction
                       ! 1 for Linear tetrahedron method
  !                    ! 2 for Optimized tetrahedron method
  INTEGER :: &
       ntetra, &         ! number of tetrahedra
       nntetra           ! k-points per tetrahedron used to compute weights
                         ! 4 for linear / 20 for optimized tetrahedron method
  INTEGER, ALLOCATABLE :: &
       tetra(:,:)        ! index of k-points in a given tetrahedron
                         ! shape (nntetra,ntetra)
  !
  REAL(dp), ALLOCATABLE :: &
       wlsm(:,:)         ! Weights for the optimized tetrahedron method
  !
  PUBLIC :: tetra, ntetra, nntetra
  PUBLIC :: tetra_init, tetra_weights, tetra_weights_only, tetra_dos_t
  PUBLIC :: opt_tetra_init, opt_tetra_weights, opt_tetra_weights_only, &
  &         opt_tetra_dos_t, opt_tetra_partialdos, tetra_type, wlsm
  PUBLIC :: deallocate_tetra
  !
CONTAINS
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE tetra_init ( nsym, s, time_reversal, t_rev, at, bg, npk, &
     k1,k2,k3, nk1,nk2,nk3, nks, xk )
  !-----------------------------------------------------------------------
  !
  ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN):: nks, nsym, t_rev(48), s(3,3,48), npk, &
                        k1, k2, k3, nk1, nk2, nk3
  LOGICAL, INTENT (IN) :: time_reversal
  real(DP), INTENT(IN) :: at(3,3), bg(3,3), xk(3,npk)
  !
  real(DP) :: xkr(3), deltap(3), deltam(3)
  real(DP), PARAMETER:: eps=1.0d-5
  real(DP), ALLOCATABLE :: xkg(:,:)
  INTEGER :: nkr, i,j,k, ns, n, nk, ip1,jp1,kp1, &
       n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER, ALLOCATABLE:: equiv(:)
  !
  ntetra =6*nk1*nk2*nk3
  nntetra=4
  ALLOCATE ( tetra (nntetra, ntetra) )
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  nkr=nk1*nk2*nk3
  ALLOCATE (xkg( 3,nkr))
  ALLOCATE (equiv( nkr))
!
  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO

  !  locate k-points of the uniform grid in the list of irreducible k-points
  !  that was previously calculated

  !  bring irreducible k-points to crystal axis
  CALL cryst_to_cart (nks,xk,at,-1)
  !
  DO nk=1,nkr
     DO n=1,nks
        DO ns=1,nsym
           DO i=1,3
              xkr(i) = s(i,1,ns) * xk(1,n) + &
                       s(i,2,ns) * xk(2,n) + &
                       s(i,3,ns) * xk(3,n)
           ENDDO
           IF(t_rev(ns)==1) xkr = -xkr
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           DO i=1,3
              deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
              deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
           ENDDO
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           IF ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .or. ( time_reversal .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) < eps ) ) THEN
              !  equivalent irreducible k-point found
              equiv(nk) = n
              GOTO 15
           ENDIF
        ENDDO
     ENDDO
     !  equivalent irreducible k-point found - something wrong
     CALL errore('tetra_init','cannot locate  k point',nk)
15   CONTINUE
  ENDDO

  DO n=1,nks
     DO nk=1,nkr
        IF (equiv(nk)==n) GOTO 20
     ENDDO
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     CALL errore('tetra_init','cannot remap grid on k-point list',n)
20   CONTINUE
  ENDDO

  !  bring irreducible k-points back to cartesian axis
  CALL cryst_to_cart (nks,xk,bg, 1)

  !  construct tetrahedra

  DO i=1,nk1
     DO j=1,nk2
        DO k=1,nk3
           !  n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = mod(i,nk1)+1
           jp1 = mod(j,nk2)+1
           kp1 = mod(k,nk3)+1
           n1 = (  k-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n2 = (  k-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n3 = (  k-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n4 = (  k-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n5 = (kp1-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n6 = (kp1-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n7 = (kp1-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n8 = (kp1-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           !  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
           n  = 6 * ( (k-1) + (j-1)*nk3 + (i-1)*nk3*nk2 )

           tetra (1,n+1) = equiv(n1)
           tetra (2,n+1) = equiv(n2)
           tetra (3,n+1) = equiv(n3)
           tetra (4,n+1) = equiv(n6)

           tetra (1,n+2) = equiv(n2)
           tetra (2,n+2) = equiv(n3)
           tetra (3,n+2) = equiv(n4)
           tetra (4,n+2) = equiv(n6)

           tetra (1,n+3) = equiv(n1)
           tetra (2,n+3) = equiv(n3)
           tetra (3,n+3) = equiv(n5)
           tetra (4,n+3) = equiv(n6)

           tetra (1,n+4) = equiv(n3)
           tetra (2,n+4) = equiv(n4)
           tetra (3,n+4) = equiv(n6)
           tetra (4,n+4) = equiv(n8)

           tetra (1,n+5) = equiv(n3)
           tetra (2,n+5) = equiv(n6)
           tetra (3,n+5) = equiv(n7)
           tetra (4,n+5) = equiv(n8)

           tetra (1,n+6) = equiv(n3)
           tetra (2,n+6) = equiv(n5)
           tetra (3,n+6) = equiv(n6)
           tetra (4,n+6) = equiv(n7)
        ENDDO
     ENDDO
  ENDDO

  !  check

  DO n=1,ntetra
     DO i=1,nntetra
        IF ( tetra(i,n)<1 .or. tetra(i,n)>nks ) &
             CALL errore ('tetra_init','something wrong',n)
     ENDDO
  ENDDO

  DEALLOCATE(equiv)
  DEALLOCATE(xkg)

  RETURN
  END SUBROUTINE tetra_init
  !
  !-----------------------------------------------------------------------
  SUBROUTINE opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, npk, &
       k1, k2, k3, nk1, nk2, nk3, nks, xk, kstep)
  !-----------------------------------------------------------------------------
  !
  ! This routine set the corners and additional points for each tetrahedron
  !
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN):: &
       & nks, & ! Total # of k in irreducible BZ
       & nsym, & ! # of crystalline symmetries
       & t_rev(48), & ! time reversal flag, for noncolinear magnetism
       & s(3,3,48), & ! Crystalline symmetry operators in reciprocal fractional space
       & npk, & ! Maximum # of k (for dimension size)
       & k1, k2, k3, & ! = 0 for unshifted, = 1 for shifted grid
       & nk1, nk2, nk3, & ! k-point grid
       & kstep            ! = 1 for pw.x, = 2 for ph.x
  !
  LOGICAL,INTENT (IN) :: &
       & time_reversal !if .TRUE. the system has time reversal symmetry
  !
  REAL(dp),INTENT(IN) :: &
       & at(3,3), & ! Direct lattice vectors [Bohr]
       & bg(3,3), & ! Reciplocal lattice vectors [2 pi / a]
       & xk(3,npk)  ! k points [2 pi / a]
  !
  REAL(dp),PARAMETER :: eps = 1e-5_dp
  !
  INTEGER :: isym, i1, i2, i3, itet, itettot, ii, ik, jk, &
       &          ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3), &
       &          equiv(nk1 * nk2 * nk3) ! index equivalent k in irr-BZ
  !
  REAL(dp) :: xkr(3), l(4), bvec2(3,3), bvec3(3,4), xkg(3, nk1 * nk2 * nk3), &
  &           deltap(3), deltam(3)
  !
  ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
  !
  bvec2(1:3,1) = bg(1:3,1) / REAL(nk1, dp)
  bvec2(1:3,2) = bg(1:3,2) / REAL(nk2, dp)
  bvec2(1:3,3) = bg(1:3,3) / REAL(nk3, dp)
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  DO ii = 1, 4
     l(ii) = DOT_PRODUCT(bvec3(1:3, ii), bvec3(1:3, ii))
  END DO
  !
  ii = MINLOC(l(1:4),1)
  !
  ivvec0(1:4) = (/ 0, 0, 0, 0 /)
  !
  divvec(1:4,1) = (/ 1, 0, 0, 0 /)
  divvec(1:4,2) = (/ 0, 1, 0, 0 /)
  divvec(1:4,3) = (/ 0, 0, 1, 0 /)
  divvec(1:4,4) = (/ 0, 0, 0, 1 /)
  !
  ivvec0(ii) = 1
  divvec(ii, ii) = - 1
  !
  ! Divide a subcell into 6 tetrahedra
  !
  itet = 0
  DO i1 = 1, 3
     DO i2 = 1, 3
        IF(i2 == i1) CYCLE
        DO i3 = 1, 3
           IF(i3 == i1 .OR. i3 == i2) CYCLE
           !
           itet = itet + 1
           !
           ivvec(1:3,1,itet) = ivvec0(1:3)
           ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
           ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
           ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
           !
        END DO
     END DO
  END DO
  !
  ! Additional points surrounding the tetrahedron
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  !
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
  !
  ! Set the weight for the each tetrahedron method
  !
  ntetra  = 6*nk1*nk2*nk3
  !
  IF(tetra_type == 1) THEN
     !
     WRITE(stdout,*) "    [opt_tetra]  Linear tetrahedron method is used."
     !
     nntetra = 4
     ALLOCATE ( tetra(nntetra,ntetra), wlsm(4,nntetra) )
     wlsm(:,:) = 0.0_dp
     !
     wlsm(1,1) = 1.0_dp
     wlsm(2,2) = 1.0_dp
     wlsm(3,3) = 1.0_dp
     wlsm(4,4) = 1.0_dp
     !
  ELSE IF(tetra_type == 2) THEN
     !
     WRITE(stdout,*) "    [opt_tetra]  Optimized tetrahedron method is used."
     !
     nntetra = 20
     ALLOCATE ( tetra(nntetra,ntetra), wlsm(4,nntetra) )
     !
     wlsm(1, 1: 4) = REAL((/1440,    0,   30,    0/), dp)
     wlsm(2, 1: 4) = REAL((/   0, 1440,    0,   30/), dp)
     wlsm(3, 1: 4) = REAL((/  30,    0, 1440,    0/), dp)
     wlsm(4, 1: 4) = REAL((/   0,   30,    0, 1440/), dp)
     !
     wlsm(1, 5: 8) = REAL((/ -38,    7,   17,  -28/), dp)
     wlsm(2, 5: 8) = REAL((/ -28,  -38,    7,   17/), dp)
     wlsm(3, 5: 8) = REAL((/  17,  -28,  -38,    7/), dp)
     wlsm(4, 5: 8) = REAL((/   7,   17,  -28,  -38/), dp)
     !
     wlsm(1, 9:12) = REAL((/ -56,    9,  -46,    9/), dp)
     wlsm(2, 9:12) = REAL((/   9,  -56,    9,  -46/), dp)
     wlsm(3, 9:12) = REAL((/ -46,    9,  -56,    9/), dp)
     wlsm(4, 9:12) = REAL((/   9,  -46,    9,  -56/), dp)
     !
     wlsm(1,13:16) = REAL((/ -38,  -28,   17,    7/), dp)
     wlsm(2,13:16) = REAL((/   7,  -38,  -28,   17/), dp)
     wlsm(3,13:16) = REAL((/  17,    7,  -38,  -28/), dp)
     wlsm(4,13:16) = REAL((/ -28,   17,    7,  -38/), dp)
     !
     wlsm(1,17:20) = REAL((/ -18,  -18,   12,  -18/), dp)
     wlsm(2,17:20) = REAL((/ -18,  -18,  -18,   12/), dp)
     wlsm(3,17:20) = REAL((/  12,  -18,  -18,  -18/), dp)
     wlsm(4,17:20) = REAL((/ -18,   12,  -18,  -18/), dp)
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / 1260.0_dp
     !
  ELSE
     !
     CALL errore( 'opt_tetra_init','tetra_type is wrong !', tetra_type )
     !
  END IF
  !
  ! Generate a uniform grid of k-points xkg
  !
  DO i1 = 1, nk1
     DO i2 = 1, nk2
        DO i3 = 1, nk3
           xkg(1:3, 1 + (i3 - 1) + (i2 - 1) * nk3 + (i1 - 1) * nk2 * nk3) &
           &  = (REAL((/i1, i2, i3/) - 1, dp) + 0.5_dp * REAL((/k1, k2, k3/), dp)) &
           &          / REAL((/nk1, nk2, nk3/), dp)
        ENDDO ! i3
     ENDDO ! i2
  ENDDO ! i1
  !
  !  locate k-points of the uniform grid in the list of irreducible k-points
  !  that was previously calculated
  !
  !  bring irreducible k-points to crystal axis
  CALL cryst_to_cart (nks,xk,at,-1)
  !
  DO ik = 1, nk1 * nk2 * nk3
     DO jk = 1, nks, kstep
        DO isym = 1, nsym
           !
           xkr(1:3) = MATMUL(REAL(s(1:3, 1:3, isym), dp), xk(1:3, jk))
           IF(t_rev(isym) == 1) xkr(1:3) = - xkr(1:3)
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           deltap(1:3) = xkr(1:3) - xkg(1:3,ik) - nint(xkr(1:3) - xkg(1:3,ik))
           deltam(1:3) = xkr(1:3) + xkg(1:3,ik) - nint(xkr(1:3) + xkg(1:3,ik))
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           IF ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .or. ( time_reversal .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) < eps ) ) THEN
              !  equivalent irreducible k-point found
              equiv(ik) = jk
              GOTO 15
           ENDIF
           !
        ENDDO
     ENDDO
     !  equivalent irreducible k-point found - something wrong
     CALL errore('opt_tetra_init', 'cannot locate  k point', ik)
15   CONTINUE
     !
  ENDDO
  !
  DO jk = 1, nks, kstep
     DO ik = 1, nk1 * nk2 * nk3
        IF (equiv(ik) == jk) GOTO 20
     ENDDO
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     CALL errore('opt_tetra_init','cannot remap grid on k-point list',jk)
20   CONTINUE
  ENDDO
  !
  !  bring irreducible k-points back to cartesian axis
  !
  CALL cryst_to_cart (nks,xk,bg, 1)
  !
  ! Construct tetrahedra
  !
  itettot = 0
  DO i1 = 1, nk1
     DO i2 = 1, nk2
        DO i3 = 1, nk3
           !
           DO itet = 1, 6
              !
              itettot = itettot + 1
              !
              DO ii = 1, nntetra
                 !
                 ikv(1:3) = (/i1, i2, i3/) - 1
                 ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
                 ikv(1:3) = MODULO(ikv(1:3), (/nk1, nk2, nk3/))
                 !
                 ik = ikv(3) + nk3 * (ikv(2) + nk2 * ikv(1)) + 1
                 !
                 tetra(ii, itettot) = equiv(ik)
                 !
              END DO ! ii
              !
           END DO ! itet
           !
        END DO ! i3
     END DO ! i2
  END DO ! i1
  !
END SUBROUTINE opt_tetra_init
!
!--------------------------------------------------------------------
subroutine tetra_weights (nks, nspin, nbnd, nelec, ntetra, tetra, et, &
     ef, wg, is, isk )
  !--------------------------------------------------------------------
  !
  ! ... calculates Ef and weights with the tetrahedron method (P.E.Bloechl)
  ! ... Wrapper routine: computes first Ef, then the weights
  !
  USE kinds
  implicit none
  ! I/O variables
  integer, intent(in) :: nks, nspin, is, isk(nks), nbnd, ntetra, &
       tetra (4, ntetra)
  real(DP), intent(in) :: et (nbnd, nks), nelec
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  real(DP), intent(out) :: ef
  ! local variables
  real(DP), external :: efermit

  ! Calculate the Fermi energy ef

  ef = efermit (et, nbnd, nks, nelec, nspin, ntetra, tetra, is, isk)
  !
  ! if efermit cannot find a sensible value for Ef it returns Ef=1d10
  !
  if (abs(ef) > 1.0d8) call errore ('tetra_weights', 'bad Fermi energy ',1)
  !
  CALL tetra_weights_only (nks, nspin, is, isk, nbnd, nelec, ntetra, &
       tetra, et, ef, wg)
  !
  return
end subroutine tetra_weights

!--------------------------------------------------------------------
subroutine tetra_weights_only (nks, nspin, is, isk, nbnd, nelec, ntetra, &
     tetra, et, ef, wg)
  !--------------------------------------------------------------------
  !
  ! ... calculates weights with the tetrahedron method (P.E.Bloechl)
  ! ... Fermi energy has to be calculated in previous step
  ! ... Generalization to noncollinear case courtesy of Iurii Timrov

  USE kinds
  implicit none
  ! I/O variables
  integer, intent(in) :: nks, nspin, is, isk(nks), nbnd, ntetra, &
       tetra (4, ntetra)
  real(DP), intent(in) :: et (nbnd, nks), nelec, ef
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  ! local variables
  real(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra (4), dosef
  integer :: ik, ibnd, nt, nk, ns, i, kp1, kp2, kp3, kp4, itetra (4)
  integer :: nspin_lsda
  !
  do ik = 1, nks
     if (is /= 0) then
        if (isk(ik) .ne. is) cycle
     end if
     do ibnd = 1, nbnd
        wg (ibnd, ik) = 0.d0
     enddo
  enddo

  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF

  do ns = 1, nspin_lsda
     if (is /= 0) then
        if (ns .ne. is) cycle
     end if
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           !
           ! etetra are the energies at the vertexes of the nt-th tetrahedron
           !
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           !
           ! ...sort in ascending order: e1 < e2 < e3 < e4
           !
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           !
           ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
           !
           kp1 = tetra (itetra (1), nt) + nk
           kp2 = tetra (itetra (2), nt) + nk
           kp3 = tetra (itetra (3), nt) + nk
           kp4 = tetra (itetra (4), nt) + nk
           !
           ! calculate weights wg
           !
           if (ef.ge.e4) then
              wg (ibnd, kp1) = wg (ibnd, kp1) + 0.25d0 / ntetra
              wg (ibnd, kp2) = wg (ibnd, kp2) + 0.25d0 / ntetra
              wg (ibnd, kp3) = wg (ibnd, kp3) + 0.25d0 / ntetra
              wg (ibnd, kp4) = wg (ibnd, kp4) + 0.25d0 / ntetra
           elseif (ef.lt.e4.and.ef.ge.e3) then
              c4 = 0.25d0 / ntetra * (e4 - ef) **3 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              dosef = 3.d0 / ntetra * (e4 - ef) **2 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              wg (ibnd, kp1) = wg (ibnd, kp1) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + 0.25d0 / ntetra - c4 * &
                   (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
                   + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
                   et (ibnd, kp4) ) / 40.d0
           elseif (ef.lt.e3.and.ef.ge.e2) then
              c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
              c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
                   / (e4 - e1) / (e3 - e2) / (e3 - e1)
              c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
                   / (e3 - e2) / (e4 - e1)
              dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
                   (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
              wg (ibnd, kp1) = wg (ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
                   / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
                   * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + (c1 + c2) * (ef - e1) &
                   / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
                   / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
                   e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           elseif (ef.lt.e2.and.ef.ge.e1) then
              c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              wg (ibnd, kp1) = wg (ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
                   * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg (ibnd, kp2) = wg (ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg (ibnd, kp3) = wg (ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg (ibnd, kp4) = wg (ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           endif
        enddo
     enddo


  enddo
  ! add correct spin normalization (2 for LDA, 1 for all other cases)
  IF ( nspin == 1 ) wg (:,1:nks) = wg (:,1:nks) * 2.d0
  !
  return
end subroutine tetra_weights_only
!----------------------------------------------------------------------------
SUBROUTINE opt_tetra_weights(nks,nspin,nbnd,nelec,ntetra,tetra,et,ef,wg, is, isk)
  !----------------------------------------------------------------------------
  !
  ! Calculate Fermi energy by using bisection method
  !
  INTEGER,INTENT(IN) :: &
  & nks,    & ! The total # of k in irr-BZ
  & nspin,  & ! The # of spin components
  & nbnd,   & ! The # of bands
  & ntetra, & ! Total # of k in whole BZ
  & is,     & ! spin index (up/down)
  & isk(nks)  ! which k is up/down
  !
  INTEGER,INTENT(IN) :: &
  & tetra(nntetra, ntetra) ! index of k in each tetrahedra
  !
  REAL(dp),INTENT(IN) :: &
  & et(nbnd, nks), & ! Kohn Sham energy [Ry]
  & nelec            ! The # of electrons
  !
  REAL(dp),INTENT(INOUT) :: &
  & wg(nbnd, nks), & ! Intetration weight of each k
  & ef               ! The Fermi energy
  !
  INTEGER :: iter, maxiter = 300
  REAL(dp) :: elw, eup, sumkmid, eps = 1.0e-10_dp
  !
  ! find bounds for the Fermi energy.
  !
  elw = minval(et(1:nbnd,1:nks))
  eup = maxval(et(1:nbnd,1:nks))
  !
  ! Bisection method
  !
  DO iter = 1, maxiter
     !
     ef = (eup + elw) * 0.5_dp
     !
     ! Calc. # of electrons 
     !
     CALL opt_tetra_weights_only(nks, nspin, nbnd, ntetra, tetra, et, ef, wg, is, isk)
     !
     IF(is == 0) THEN
        sumkmid = sum(wg(1:nbnd,1:nks))
     ELSE IF(is == 1) THEN
        sumkmid = sum(wg(1:nbnd,1:nks / 2))
     ELSE IF(is == 2) THEN
        sumkmid = SUM(wg(1:nbnd,nks / 2 + 1:nks))
     END IF
     !
     ! convergence check
     !
     IF(ABS(sumkmid - nelec) < eps) THEN
        EXIT
     ELSEIF(sumkmid < nelec) THEN
        elw = ef
     ELSE
        eup = ef
     END IF
     !
  END DO ! iter
  IF(iter >= maxiter) CALL errore("opt_tetra_weights", "Not converged", iter)
  !
END SUBROUTINE opt_tetra_weights
!
!--------------------------------------------------------------------------------------
SUBROUTINE opt_tetra_weights_only(nks, nspin, nbnd, ntetra, tetra, et, ef, wg, is, isk)
  !------------------------------------------------------------------------------------
  !
  ! Calculate Occupation with given Fermi energy
  !
  INTEGER,INTENT(IN) :: nbnd, nks, ntetra, nspin, is, isk(nks)
  INTEGER,INTENT(IN) :: tetra(nntetra, ntetra)
  REAL(dp),INTENT(IN) :: ef, et(nbnd,nks)
  REAL(dp),INTENT(INOUT) :: wg(nbnd,nks)
  !
  INTEGER :: ns, ik, nt, ibnd, jbnd, kbnd, ii, itetra(4), nk
  REAL(dp) :: e(4), wg0(4), C(3), a(4,4), wg1
  integer :: nspin_lsda
  !
  ! Zero-clear only "is"th spin component
  !
  do ik = 1, nks
     if (is /= 0) then
        if (isk(ik) .ne. is) cycle
     end if
     do ibnd = 1, nbnd
        wg (ibnd, ik) = 0.d0
     enddo
  enddo

  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF

  DO ns = 1, nspin_lsda
     if (is /= 0) then
        if (ns .ne. is) cycle
     end if
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nks / 2
     END IF
     !
     DO nt = 1, ntetra
        !
        DO ibnd = 1, nbnd
           !
           e(1:4) = 0.0_dp
           DO ii = 1, nntetra
              !
              ik = tetra(ii, nt) + nk
              e(1:4) = e(1:4) + wlsm(1:4,ii) * et(ibnd,ik)
              !
           END DO
           !
           itetra(1) = 0
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( ef - e(1:4) ) / (e(ii) - e(1:4) )
           END DO
           !
           IF( e(1) <= ef .AND. ef < e(2) ) THEN
              !
              C(1) = a(2,1) * a(3,1) * a(4,1) * 0.25_dp
              wg0(1) = C(1) * (1.0_dp + a(1,2) + a(1,3) + a(1,4))
              wg0(2) = C(1) * a(2,1)
              wg0(3) = C(1) * a(3,1)
              wg0(4) = C(1) * a(4,1)
              !
           ELSE IF( e(2) <= ef .AND. ef < e(3)) THEN
              !
              C(1) = a(4,1) * a(3,1) * 0.25_dp
              C(2) = a(4,1) * a(3,2) * a(1,3) * 0.25_dp
              C(3) = a(4,2) * a(3,2) * a(1,4) * 0.25_dp
              !
              wg0(1) = C(1) + (C(1) + C(2)) * a(1,3) + (C(1) + C(2) + C(3)) * a(1,4)
              wg0(2) = C(1) + C(2) + C(3) + (C(2) + C(3)) * a(2,3) + C(3) * a(2,4)
              wg0(3) = (C(1) + C(2)) * a(3,1) + (C(2) + C(3)) * a(3,2)
              wg0(4) = (C(1) + C(2) + C(3)) * a(4,1) + C(3) * a(4,2)
              !
           ELSE IF( e(3) <= ef .AND. ef < e(4)) THEN
              !
              C(1) = a(1,4) * a(2,4) * a(3,4)
              !
              wg0(1) = 1.0_dp - C(1) * a(1,4)
              wg0(2) = 1.0_dp - C(1) * a(2,4)
              wg0(3) = 1.0_dp - C(1) * a(3,4)
              wg0(4) = 1.0_dp - C(1) * (1.0_dp + a(4,1) + a(4,2) + a(4,3))
              !
              wg0(1:4) = wg0(1:4) * 0.25_dp
              !
           ELSE IF(e(4) <= ef ) THEN
              !
              wg0(1:4) = 0.25_dp
              !
           ELSE
              !
              wg0(1:4) = 0.0_dp
              !
           END IF
           !
           wg0(1:4) = wg0(1:4) / REAL(ntetra, dp)
           !
           DO ii = 1, nntetra
              !
              ik = tetra(ii, nt) + nk
              wg(ibnd, ik) = wg(ibnd, ik) &
              &   + DOT_PRODUCT(wlsm(itetra(1:4),ii), wg0(1:4))
              !
           END DO
           !
        END DO ! ibnd
        !
     END DO ! nt
     !
  END DO
  !
  ! Average weights of degenerated states
  !
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        !
        wg1 = wg(ibnd,ik)
        !
        DO jbnd = ibnd + 1, nbnd
           !
           IF(ABS(et(ibnd,ik) - et(jbnd,ik)) < 1e-6_dp) THEN
              wg1 = wg1 + wg(jbnd,ik)
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 wg(kbnd,ik) = wg1 / real(jbnd - ibnd, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
  ! add correct spin normalization (2 for LDA, 1 for all other cases)
  IF ( nspin == 1 ) wg (1:nbnd,1:nks) = wg (1:nbnd,1:nks) * 2.0_dp
  !
END SUBROUTINE opt_tetra_weights_only
!
!--------------------------------------------------------------------
subroutine tetra_dos_t (et, nspin, nbnd, nks, e, dost)
  !------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  integer, intent(in) :: nspin, nbnd, nks

  real(DP), intent(in) :: et (nbnd, nks), e
  REAL(dp), INTENT(OUT):: dost (2)
  !
  integer :: itetra (4), nk, ns, nt, ibnd, i
  real(DP) :: etetra (4), e1, e2, e3, e4
  integer :: nspin0

  if (nspin==4) then
     nspin0=1
  else 
     nspin0=nspin
  endif
  do ns = 1, nspin0
     dost (ns) = 0.d0
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           do i = 1, 4
              etetra (i) = et (ibnd, tetra (i, nt) + nk)
           enddo
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           if (e.lt.e4.and.e.ge.e3) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
           elseif (e.lt.e2.and.e.gt.e1) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
           endif
        enddo

     enddo

     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 

     if ( nspin == 1 ) dost (ns) = dost (ns) * 2.d0

  enddo
  return
end subroutine tetra_dos_t

! Optimized tetrahedron method
!
!--------------------------------------------------------------------
SUBROUTINE opt_tetra_dos_t (et, nspin, nbnd, nks, e, dost)
  !------------------------------------------------------------------
  !
  implicit none
  integer :: nspin, nbnd, nks

  real(DP) :: et (nbnd, nks), e, dost (2)
  integer :: itetra (4), nk, ns, nt, ibnd, i

  real(DP) :: etetra (4), e1, e2, e3, e4
  integer :: nspin0

  if (nspin==4) then
     nspin0=1
  else 
     nspin0=nspin
  endif
  do ns = 1, nspin0
     dost (ns) = 0.d0
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then
        nk = 0
     else
        nk = nks / 2
     endif
     do nt = 1, ntetra
        do ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           etetra(1:4) = 0.0_dp
           do i = 1, nntetra
              etetra(1:4) = etetra(1:4) + wlsm(1:4,i) * et(ibnd,tetra(i, nt) + nk)
           end do
           itetra (1) = 0
           call hpsort (4, etetra, itetra)
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           if (e.lt.e4.and.e.ge.e3) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
           elseif (e.lt.e2.and.e.gt.e1) then
              dost (ns) = dost (ns) + 1.d0 / ntetra * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
           endif
        enddo

     enddo

     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 

     if ( nspin == 1 ) dost (ns) = dost (ns) * 2.d0

  enddo
  return
end SUBROUTINE opt_tetra_dos_t
!
! Compute partial dos
!
SUBROUTINE opt_tetra_partialdos(nspin0, kresolveddos,ne,natomwfc,nkseff,&
&                               Emin,DeltaE, proj, pdos, dostot, nspindos)
  !
  USE lsda_mod, ONLY: nspin, isk
  USE wvfct, ONLY: et, nbnd
  USE klist, ONLY: nkstot, wk
  USE spin_orb,   ONLY: lspinorb
  USE noncollin_module, ONLY : noncolin
  USE constants, ONLY: rytoev
  !
  implicit none
  !
  LOGICAL :: kresolveddos
  !
  INTEGER,INTENT(in) :: ne, natomwfc, nkseff, nspin0, nspindos
  !
  REAL(dp),INTENT(in) :: Emin, deltaE, proj(natomwfc, nbnd, nkstot)
  REAL(dp),INTENT(out) :: pdos(0:ne,natomwfc,nspin0,nkseff), &
  &                       dostot(0:ne,nspindos,nkseff)
  !
  INTEGER :: ns, nk, nt, ibnd, itetra(4), ii, ik, ie, nwfc, nspin1
  REAL(dp) :: etetra(4), e0, e1, e2, e3, e4, wgt(4), G, &
  &           f12, f13, f14, f21, f23, f24, f31, f32, f34, f41, f42
  !
  IF(nspin == 2) THEN
     nspin1 = nspin
  ELSE
     nspin1 = 1
  END IF
  !
  DO ns = 1, nspin1
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     ENDIF
     !
     DO nt = 1, ntetra
        !
        DO ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           etetra(1:4) = 0.0_dp
           DO ii = 1, nntetra
              etetra(1:4) = etetra(1:4) + wlsm(1:4,ii) * et(ibnd,tetra(ii, nt) + nk)
           ENDDO
           !
           itetra (1) = 0
           CALL hpsort (4, etetra, itetra)
           e1 = etetra (1)
           e2 = etetra (2)
           e3 = etetra (3)
           e4 = etetra (4)
           !
           DO ie = 0, ne
              !
              e0 = Emin + DeltaE * REAL(ie, dp)
              !
              IF ( e3 < e0 .and. e0 < e4) THEN
                 !
                 f14 = (e0-e4)/(e1-e4)
                 f24 = (e0-e4)/(e2-e4)
                 f34 = (e0-e4)/(e3-e4)
                 !
                 G  =  3.0_dp * f14 * f24 * f34 / (e4-e0)
                 wgt(1) =  f14 / 3.0_dp
                 wgt(2) =  f24 / 3.0_dp
                 wgt(3) =  f34 / 3.0_dp
                 wgt(4) =  (3.0_dp - f14 - f24 - f34 ) / 3.0_dp
                 !
              ELSE IF ( e2 < e0 .and. e0 < e3 ) THEN
                 !
                 f13 = (e0-e3)/(e1-e3)
                 f31 = 1.0_dp - f13
                 f14 = (e0-e4)/(e1-e4)
                 f41 = 1.0_dp-f14
                 f23 = (e0-e3)/(e2-e3)
                 f32 = 1.0_dp - f23
                 f24 = (e0-e4)/(e2-e4)
                 f42 = 1.0_dp - f24
                 !
                 G   =  3.0_dp * (f23*f31 + f32*f24)
                 wgt(1)  =  f14 / 3.0_dp + f13*f31*f23 / G
                 wgt(2)  =  f23 / 3.0_dp + f24*f24*f32 / G
                 wgt(3)  =  f32 / 3.0_dp + f31*f31*f23 / G
                 wgt(4)  =  f41 / 3.0_dp + f42*f24*f32 / G
                 G   =  G / (e4-e1)
                 !
              ELSE IF ( e1 < e0 .and. e0 < e2 ) THEN
                 !
                 f12 = (e0-e2)/(e1-e2)
                 f21 = 1.0_dp - f12
                 f13 = (e0-e3)/(e1-e3)
                 f31 = 1.0_dp - f13
                 f14 = (e0-e4)/(e1-e4)
                 f41 = 1.0_dp - f14
                 !
                 G  =  3.0_dp * f21 * f31 * f41 / (e0-e1)
                 wgt(1) =  (f12 + f13 + f14) / 3.0_dp
                 wgt(2) =  f21 / 3.0_dp
                 wgt(3) =  f31 / 3.0_dp
                 wgt(4) =  f41 / 3.0_dp
                 !
              ELSE
                 !
                 wgt(1:4) = 0.0_dp
                 G = 0.0_dp
                 CYCLE
                 !
              END IF
              !
              IF (kresolveddos) THEN
                 !
                 DO ii = 1, nntetra
                    !
                    ik = tetra(ii, nt) + nk
                    DO nwfc = 1, natomwfc
                       pdos(ie,nwfc,ns,ik - nk) = pdos(ie,nwfc,ns,ik - nk) &
                       & + proj (nwfc, ibnd, ik) &
                       & * DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G / wk(ik)
                    ENDDO
                    !
                    IF(nspindos == 1) THEN
                       dostot(ie,1,ik - nk) = dostot(ie,1,ik - nk) &
                       & + DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G / wk(ik)
                    ELSE
                       dostot(ie,ns,ik - nk) = dostot(ie,ns,ik - nk) &
                       & + DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G / wk(ik)
                    END IF
                    !
                 END DO
                 !
              ELSE
                 !
                 DO ii = 1, nntetra
                    !
                    ik = tetra(ii, nt) + nk
                    DO nwfc = 1, natomwfc
                       pdos(ie,nwfc,ns,1) = pdos(ie,nwfc,ns,1) &
                       &  + proj (nwfc, ibnd, ik) * DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G
                    ENDDO
                    !
                 END DO
                 !
                 IF(nspindos == 1) THEN
                    dostot(ie,1,1) = dostot(ie,1,1) + SUM(wgt(1:4)) * G
                 ELSE
                    dostot(ie,ns,1) = dostot(ie,ns,1) + SUM(wgt(1:4)) * G
                 END IF
                 !
              END IF
              !
           END DO
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  pdos(0:ne,1:natomwfc,1:nspin0,1:nkseff) = pdos(0:ne,1:natomwfc,1:nspin0,1:nkseff) &
  &                                       / (rytoev * REAL(ntetra,dp))
  dostot(0:ne,1:nspindos,1:nkseff) = dostot(0:ne,1:nspindos,1:nkseff) / (rytoev * REAL(ntetra,dp))
  !
END SUBROUTINE opt_tetra_partialdos
!
SUBROUTINE deallocate_tetra ( )
   IF ( ALLOCATED(tetra) ) DEALLOCATE (tetra)
   IF ( ALLOCATED(wlsm ) ) DEALLOCATE (wlsm )
END SUBROUTINE deallocate_tetra
!
END MODULE ktetra
