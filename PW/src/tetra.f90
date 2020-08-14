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
  !! Variables used by the tetrahedron method.  
  !! Three versions are implemented: Linear, Optimized, Bloechl.  
  !! Linear and Optimized tetrahedra contributed by Mitsuaki Kawamura,
  !! University of Tokyo.
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  INTEGER :: tetra_type = 0
  !! 0 for Bloechl's correction,  
  !! 1 for Linear tetrahedron method,  
  !! 2 for Optimized tetrahedron method.
  INTEGER :: ntetra = 0 
  !! number of tetrahedra
  INTEGER :: nntetra
  !! k-points per tetrahedron used to compute weights.  
  !! 4 for linear / 20 for optimized tetrahedron method
  INTEGER, ALLOCATABLE :: tetra(:,:)
  !! index of k-points in a given tetrahedron shape (nntetra,ntetra)
  REAL(DP), ALLOCATABLE :: wlsm(:,:)
  !! Weights for the optimized tetrahedron method
  !
  PUBLIC :: tetra, ntetra, nntetra
  PUBLIC :: tetra_init, tetra_weights, tetra_weights_only,  tetra_dos_t
  PUBLIC :: opt_tetra_init, opt_tetra_weights, opt_tetra_weights_only,  &
            opt_tetra_dos_t, opt_tetra_partialdos, tetra_type, wlsm
  PUBLIC :: deallocate_tetra
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE tetra_init( nsym, s, time_reversal, t_rev, at, bg, npk, &
                         k1,k2,k3, nk1,nk2,nk3, nks, xk )
  !-----------------------------------------------------------------------
  !! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994).
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nsym
  !! total number of crystal symmetries
  INTEGER, INTENT(IN) :: t_rev(48)
  !! time reversal flag, for noncolinear magnetism
  INTEGER, INTENT(IN) :: s(3,3,48)
  !! symmetry matrices, in crystal axis
  INTEGER, INTENT(IN) :: npk
  !! max number of k-points  
  INTEGER, INTENT(IN) :: k1
  !! the offset from the origin, direction 1
  INTEGER, INTENT(IN) :: k2
  !! the offset from the origin, direction 2
  INTEGER, INTENT(IN) :: k3
  !! the offset from the origin, direction 3
  INTEGER, INTENT(IN) :: nk1
  !! the special-point grid, direction 1
  INTEGER, INTENT(IN) :: nk2
  !! the special-point grid, direction 2
  INTEGER, INTENT(IN) :: nk3
  !! the special-point grid, direction 3
  LOGICAL, INTENT(IN) :: time_reversal
  !! if .TRUE. the system has time reversal symmetry
  REAL(DP), INTENT(IN) :: at(3,3)
  !! direct lattice primitive vectors.  
  !! at(:,i) are the lattice vectors of the simulation cell, a_i,
  !! in alat units: a_i(:) = at(:,i)/alat
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! reciprocal lattice primitive vectors.  
  !! bg(:,i) are the reciprocal lattice vectors, b_i,
  !! in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
  REAL(DP), INTENT(INOUT) :: xk(3,npk)
  !! coordinates of k points
  !
  ! ... local variables
  !
  REAL(DP) :: xkr(3), deltap(3), deltam(3)
  REAL(DP), PARAMETER:: eps=1.0d-5
  REAL(DP), ALLOCATABLE :: xkg(:,:)
  INTEGER :: nkr, i,j,k, ns, n, nk, ip1, jp1, kp1, &
             n1, n2, n3, n4, n5, n6, n7, n8
  INTEGER, ALLOCATABLE:: equiv(:)
  !
  ntetra  = 6*nk1*nk2*nk3
  nntetra = 4
  !
  IF (.NOT. ALLOCATED(tetra))  ALLOCATE( tetra(nntetra,ntetra) )
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  nkr = nk1*nk2*nk3
  !
  ALLOCATE( xkg(3,nkr) )
  ALLOCATE( equiv(nkr) )
  !
  DO i = 1, nk1
     DO j = 1, nk2
        DO k = 1, nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = DBLE(i-1)/nk1 + DBLE(k1)/2/nk1
           xkg(2,n) = DBLE(j-1)/nk2 + DBLE(k2)/2/nk2
           xkg(3,n) = DBLE(k-1)/nk3 + DBLE(k3)/2/nk3
        ENDDO
     ENDDO
  ENDDO
  !
  ! ... locate k-points of the uniform grid in the list of irreducible k-points
  ! that was previously calculated
  !
  ! ... bring irreducible k-points to crystal axis
  CALL cryst_to_cart( nks, xk, at, -1 )
  !
  DO nk = 1, nkr
     DO n = 1, nks
        DO ns = 1, nsym
           DO i = 1, 3
              xkr(i) = s(i,1,ns) * xk(1,n) + &
                       s(i,2,ns) * xk(2,n) + &
                       s(i,3,ns) * xk(3,n)
           ENDDO
           IF (t_rev(ns) == 1) xkr = -xkr
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           DO i = 1, 3
              deltap(i) = xkr(i)-xkg(i,nk) - NINT(xkr(i)-xkg(i,nk))
              deltam(i) = xkr(i)+xkg(i,nk) - NINT(xkr(i)+xkg(i,nk))
           ENDDO
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           IF ( SQRT( deltap(1)**2 + &
                      deltap(2)**2 + &
                      deltap(3)**2 ) < eps .OR. ( time_reversal .AND. &
                SQRT( deltam(1)**2 + &
                      deltam(2)**2 + &
                      deltam(3)**2 ) < eps ) ) THEN
              !  equivalent irreducible k-point found
              equiv(nk) = n
              GOTO 15
           ENDIF
        ENDDO
     ENDDO
     !  equivalent irreducible k-point found - something wrong
     CALL errore( 'tetra_init', 'cannot locate  k point', nk )
15   CONTINUE
  ENDDO
  !
  DO n = 1, nks
     DO nk = 1, nkr
        IF (equiv(nk) == n) GOTO 20
     ENDDO
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     CALL errore( 'tetra_init', 'cannot remap grid on k-point list', n )
20   CONTINUE
  ENDDO
  !
  ! bring irreducible k-points back to cartesian axis
  CALL cryst_to_cart( nks, xk, bg, 1 )
  !
  ! construct tetrahedra
  !
  DO i = 1, nk1
     DO j = 1, nk2
        DO k = 1, nk3
           ! n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = MOD(i,nk1) + 1
           jp1 = MOD(j,nk2) + 1
           kp1 = MOD(k,nk3) + 1
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
           !
           tetra (1,n+1) = equiv(n1)
           tetra (2,n+1) = equiv(n2)
           tetra (3,n+1) = equiv(n3)
           tetra (4,n+1) = equiv(n6)
           ! 
           tetra (1,n+2) = equiv(n2)
           tetra (2,n+2) = equiv(n3)
           tetra (3,n+2) = equiv(n4)
           tetra (4,n+2) = equiv(n6)
           !
           tetra (1,n+3) = equiv(n1)
           tetra (2,n+3) = equiv(n3)
           tetra (3,n+3) = equiv(n5)
           tetra (4,n+3) = equiv(n6)
           !
           tetra (1,n+4) = equiv(n3)
           tetra (2,n+4) = equiv(n4)
           tetra (3,n+4) = equiv(n6)
           tetra (4,n+4) = equiv(n8)
           !
           tetra (1,n+5) = equiv(n3)
           tetra (2,n+5) = equiv(n6)
           tetra (3,n+5) = equiv(n7)
           tetra (4,n+5) = equiv(n8)
           !
           tetra (1,n+6) = equiv(n3)
           tetra (2,n+6) = equiv(n5)
           tetra (3,n+6) = equiv(n6)
           tetra (4,n+6) = equiv(n7)
        ENDDO
     ENDDO
  ENDDO
  !
  !  check
  !
  DO n=1,ntetra
     DO i=1,nntetra
        IF ( tetra(i,n)<1 .or. tetra(i,n)>nks ) &
             CALL errore ('tetra_init','something wrong',n)
     ENDDO
  ENDDO
  !
  DEALLOCATE( equiv )
  DEALLOCATE( xkg )
  !
  RETURN
  !
  END SUBROUTINE tetra_init
  !
  !-----------------------------------------------------------------------
  SUBROUTINE opt_tetra_init( nsym, s, time_reversal, t_rev, at, bg, npk, &
                             k1, k2, k3, nk1, nk2, nk3, nks, xk, kstep )
  !-----------------------------------------------------------------------------
  !! This rouotine sets the corners and additional points for each tetrahedron.
  !
  USE io_global,    ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! Total # of k in irreducible BZ
  INTEGER, INTENT(IN) :: nsym
  !! # of crystalline symmetries
  INTEGER, INTENT(IN) :: t_rev(48)
  !! time reversal flag, for noncolinear magnetism
  INTEGER, INTENT(IN) :: s(3,3,48)
  !! Crystalline symmetry operators in reciprocal fractional space
  INTEGER, INTENT(IN) :: npk
  !! Maximum # of k (for dimension size)
  INTEGER, INTENT(IN) :: k1
  !! = 0 for unshifted, = 1 for shifted grid
  INTEGER, INTENT(IN) :: k2
  !! = 0 for unshifted, = 1 for shifted grid
  INTEGER, INTENT(IN) :: k3
  !! = 0 for unshifted, = 1 for shifted grid
  INTEGER, INTENT(IN) :: nk1
  !! k-point grid
  INTEGER, INTENT(IN) :: nk2
  !! k-point grid
  INTEGER, INTENT(IN) :: nk3
  !! k-point grid
  INTEGER, INTENT(IN) :: kstep
  !! = 1 for pw.x, = 2 for ph.x
  !
  LOGICAL, INTENT(IN) :: time_reversal
  !! if .TRUE. the system has time reversal symmetry
  !
  REAL(DP), INTENT(IN) :: at(3,3)
  !! Direct lattice vectors [Bohr]
  REAL(DP), INTENT(IN) :: bg(3,3)
  !! Reciplocal lattice vectors [2 pi / a]
  !
  REAL(DP), INTENT(INOUT) :: xk(3,npk)
  !! k points [2 pi / a]
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: eps = 1e-5_dp
  !
  INTEGER :: isym, i1, i2, i3, itet, itettot, ii, ik, jk,   &
             ivvec(3,20,6), divvec(4,4), ivvec0(4), ikv(3), &
             equiv(nk1*nk2*nk3) ! index equivalent k in irr-BZ
  !
  REAL(DP) :: xkr(3), l(4), bvec2(3,3), bvec3(3,4), xkg(3,nk1*nk2*nk3), &
              deltap(3), deltam(3)
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
  ENDDO
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
        ENDDO
     ENDDO
  ENDDO
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
  IF (tetra_type == 1) THEN
     !
     WRITE(stdout,*) "    [opt_tetra]  Linear tetrahedron method is used."
     !
     nntetra = 4
     IF(.NOT. ALLOCATED(tetra)) ALLOCATE ( tetra(nntetra,ntetra) )
     IF(.NOT. ALLOCATED(wlsm))  ALLOCATE ( wlsm(4,nntetra) )
     wlsm(:,:) = 0.0_dp
     !
     wlsm(1,1) = 1.0_dp
     wlsm(2,2) = 1.0_dp
     wlsm(3,3) = 1.0_dp
     wlsm(4,4) = 1.0_dp
     !
  ELSEIF (tetra_type == 2) THEN
     !
     WRITE(stdout,*) "    [opt_tetra]  Optimized tetrahedron method is used."
     !
     nntetra = 20
     IF (.NOT. ALLOCATED(tetra)) ALLOCATE( tetra(nntetra,ntetra) )
     IF (.NOT. ALLOCATED(wlsm))  ALLOCATE( wlsm(4,nntetra) )
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
  ENDIF
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
  CALL cryst_to_cart( nks, xk, at, -1 )
  !
  DO ik = 1, nk1 * nk2 * nk3
     DO jk = 1, nks, kstep
        DO isym = 1, nsym
           !
           xkr(1:3) = MATMUL(REAL(s(1:3, 1:3, isym), dp), xk(1:3, jk))
           IF (t_rev(isym) == 1) xkr(1:3) = - xkr(1:3)
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           deltap(1:3) = xkr(1:3) - xkg(1:3,ik) - NINT(xkr(1:3) - xkg(1:3,ik))
           deltam(1:3) = xkr(1:3) + xkg(1:3,ik) - NINT(xkr(1:3) + xkg(1:3,ik))
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           IF ( SQRT(  deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .OR. ( time_reversal .AND. &
                SQRT(  deltam(1)**2 +  &
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
     CALL errore( 'opt_tetra_init', 'cannot locate  k point', ik )
     !
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
     CALL errore( 'opt_tetra_init', 'cannot remap grid on k-point list', jk )
     !
20   CONTINUE
  ENDDO
  !
  !  bring irreducible k-points back to cartesian axis
  !
  CALL cryst_to_cart( nks, xk, bg, 1 )
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
              ENDDO ! ii
              !
           ENDDO ! itet
           !
        ENDDO ! i3
     ENDDO ! i2
  ENDDO ! i1
  !
END SUBROUTINE opt_tetra_init
!
!--------------------------------------------------------------------
SUBROUTINE tetra_weights( nks, nspin, nbnd, nelec, et, &
                          ef, wg, is, isk )
  !--------------------------------------------------------------------
  !! Calculates Ef and weights with the tetrahedron method (P.E.Bloechl).  
  !! Wrapper routine: computes first Ef, then the weights. 
  !! Needs to be called only after tetrahedra have been initialized, otherwise 
  !! it will stop with an error call 
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! Total # of k in irreducible BZ
  INTEGER, INTENT(IN) :: nspin
  !! polarization label
  INTEGER, INTENT(IN) :: is
  !! spin label
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! eigenvalues of the hamiltonian
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost.
  REAL(DP), INTENT(OUT) :: ef
  !! Fermi energy
  !
  ! ... local variables
  !
  REAL(DP), EXTERNAL :: efermit
  !
  ! Calculate the Fermi energy ef
  !
  IF ( ntetra == 0 ) &
     CALL errore('tetra weigths', 'called without initialization', 1) 

  ef = efermit( et, nbnd, nks, nelec, nspin, ntetra, tetra, is, isk )
  !
  ! if efermit cannot find a sensible value for Ef it returns Ef=1d10
  !
  IF (ABS(ef) > 1.0d8) CALL errore( 'tetra_weights', 'bad Fermi energy ', 1 )
  !
  CALL tetra_weights_only( nks, nspin, is, isk, nbnd, nelec,  et, ef, wg )
  !
  RETURN
  !
END SUBROUTINE tetra_weights
!
!
!--------------------------------------------------------------------
SUBROUTINE tetra_weights_only( nks, nspin, is, isk, nbnd, nelec, et, ef, wg )
  !--------------------------------------------------------------------
  !! Calculates weights with the tetrahedron method (P.E.Bloechl).  
  !! Fermi energy has to be calculated in previous step.  
  !! Generalization to noncollinear case courtesy of Iurii Timrov.
  !! @Note (P. Delugas 8/10/2019) Needs to be called only after initializations,  
  !!       stops the program with an error call otherwise. 
  !!       
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! Total # of k in irreducible BZ
  INTEGER, INTENT(IN) :: nspin
  !! polarization label
  INTEGER, INTENT(IN) :: is
  !! spin label
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! eigenvalues of the hamiltonian
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! the weight of each k point and band
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost.
  REAL(DP), INTENT(IN) :: ef
  !! Fermi energy
  !
  ! ... local variables
  !
  REAL(DP) :: e1, e2, e3, e4, c1, c2, c3, c4, etetra(4), dosef
  INTEGER :: ik, ibnd, nt, nk, ns, i, kp1, kp2, kp3, kp4, itetra(4)
  INTEGER :: nspin_lsda
  !
  IF ( ntetra == 0 ) &
     CALL errore('tetra_weights_only: ', 'called before initialization', 1)  
  DO ik = 1, nks
     IF (is /= 0) THEN
        IF (isk(ik) /= is) CYCLE
     ENDIF
     DO ibnd = 1, nbnd
        wg (ibnd, ik) = 0.d0
     ENDDO
  ENDDO
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  ENDIF
  !
  DO ns = 1, nspin_lsda
     IF (is /= 0) THEN
        IF (ns /= is) CYCLE
     ENDIF
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nks / 2
     ENDIF
     !
     DO nt = 1, ntetra
        DO ibnd = 1, nbnd
           !
           ! etetra are the energies at the vertexes of the nt-th tetrahedron
           !
           DO i = 1, 4
              etetra(i) = et (ibnd, tetra(i,nt) + nk)
           ENDDO
           itetra (1) = 0
           CALL hpsort( 4, etetra, itetra )
           !
           ! ...sort in ascending order: e1 < e2 < e3 < e4
           !
           e1 = etetra(1)
           e2 = etetra(2)
           e3 = etetra(3)
           e4 = etetra(4)
           !
           ! kp1-kp4 are the irreducible k-points corresponding to e1-e4
           !
           kp1 = tetra(itetra(1), nt) + nk
           kp2 = tetra(itetra(2), nt) + nk
           kp3 = tetra(itetra(3), nt) + nk
           kp4 = tetra(itetra(4), nt) + nk
           !
           ! calculate weights wg
           !
           IF (ef>=e4) THEN
              !
              wg(ibnd, kp1) = wg(ibnd, kp1) + 0.25d0 / ntetra
              wg(ibnd, kp2) = wg(ibnd, kp2) + 0.25d0 / ntetra
              wg(ibnd, kp3) = wg(ibnd, kp3) + 0.25d0 / ntetra
              wg(ibnd, kp4) = wg(ibnd, kp4) + 0.25d0 / ntetra
              !
           ELSEIF (ef<e4 .AND. ef>=e3) THEN
              !
              c4 = 0.25d0 / ntetra * (e4 - ef)**3 / (e4 - e1) / (e4 - e2) &
                   / (e4 - e3)
              dosef = 3.d0 / ntetra * (e4 - ef)**2 / (e4 - e1) / (e4 - e2) &
                      / (e4 - e3)
              wg(ibnd,kp1) = wg(ibnd,kp1) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e1) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp1) ) / 40.d0
              wg(ibnd,kp2) = wg(ibnd,kp2) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e2) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp2) ) / 40.d0
              wg(ibnd,kp3) = wg(ibnd,kp3) + 0.25d0 / ntetra - c4 * &
                   (e4 - ef) / (e4 - e3) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et &
                   (ibnd, kp3) ) / 40.d0
              wg(ibnd,kp4) = wg(ibnd,kp4) + 0.25d0 / ntetra - c4 * &
                   (4.d0 - (e4 - ef) * (1.d0 / (e4 - e1) + 1.d0 / (e4 - e2) &
                   + 1.d0 / (e4 - e3) ) ) + dosef * (e1 + e2 + e3 + e4 - 4.d0 * &
                   et(ibnd,kp4) ) / 40.d0
              !
           ELSEIF (ef<e3 .AND. ef>=e2) THEN
              !
              c1 = 0.25d0 / ntetra * (ef - e1) **2 / (e4 - e1) / (e3 - e1)
              c2 = 0.25d0 / ntetra * (ef - e1) * (ef - e2) * (e3 - ef) &
                   / (e4 - e1) / (e3 - e2) / (e3 - e1)
              c3 = 0.25d0 / ntetra * (ef - e2) **2 * (e4 - ef) / (e4 - e2) &
                   / (e3 - e2) / (e4 - e1)
              dosef = 1.d0 / ntetra / (e3 - e1) / (e4 - e1) * (3.d0 * &
                   (e2 - e1) + 6.d0 * (ef - e2) - 3.d0 * (e3 - e1 + e4 - e2) &
                   * (ef - e2) **2 / (e3 - e2) / (e4 - e2) )
              wg(ibnd, kp1) = wg(ibnd, kp1) + c1 + (c1 + c2) * (e3 - ef) &
                   / (e3 - e1) + (c1 + c2 + c3) * (e4 - ef) / (e4 - e1) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg(ibnd, kp2) = wg(ibnd, kp2) + c1 + c2 + c3 + (c2 + c3) &
                   * (e3 - ef) / (e3 - e2) + c3 * (e4 - ef) / (e4 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg(ibnd, kp3) = wg(ibnd, kp3) + (c1 + c2) * (ef - e1) &
                   / (e3 - e1) + (c2 + c3) * (ef - e2) / (e3 - e2) + dosef * &
                   (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg(ibnd, kp4) = wg(ibnd, kp4) + (c1 + c2 + c3) * (ef - e1) &
                   / (e4 - e1) + c3 * (ef - e2) / (e4 - e2) + dosef * (e1 + e2 + &
                   e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
              !
           ELSEIF (ef<e2 .AND. ef>=e1) THEN
              !
              c4 = 0.25d0 / ntetra * (ef - e1) **3 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              dosef = 3.d0 / ntetra * (ef - e1) **2 / (e2 - e1) / (e3 - e1) &
                   / (e4 - e1)
              wg(ibnd, kp1) = wg(ibnd, kp1) + c4 * (4.d0 - (ef - e1) &
                   * (1.d0 / (e2 - e1) + 1.d0 / (e3 - e1) + 1.d0 / (e4 - e1) ) ) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp1) ) / 40.d0
              wg(ibnd, kp2) = wg(ibnd, kp2) + c4 * (ef - e1) / (e2 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp2) ) / 40.d0
              wg(ibnd, kp3) = wg(ibnd, kp3) + c4 * (ef - e1) / (e3 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp3) ) / 40.d0
              wg(ibnd, kp4) = wg(ibnd, kp4) + c4 * (ef - e1) / (e4 - e1) &
                   + dosef * (e1 + e2 + e3 + e4 - 4.d0 * et (ibnd, kp4) ) / 40.d0
           ENDIF
           !
        ENDDO
     ENDDO
     !
  ENDDO
  ! add correct spin normalization (2 for LDA, 1 for all other cases)
  IF ( nspin == 1 ) wg(:,1:nks) = wg(:,1:nks) * 2.d0
  !
  RETURN
  !
END SUBROUTINE tetra_weights_only
!
!
!----------------------------------------------------------------------------
SUBROUTINE opt_tetra_weights( nks, nspin, nbnd, nelec, et, ef, &
                              wg, is, isk )
  !----------------------------------------------------------------------------
  !! Calculate Fermi energy by using bisection method.
  !! To be called only after tetrahedra initialization. Stops the program
  !! with an error call otherwise
  !
  INTEGER, INTENT(IN) :: nks
  !! The total # of k in irr-BZ
  INTEGER, INTENT(IN) :: nspin
  !! The # of spin components
  INTEGER, INTENT(IN) :: nbnd
  !! The # of bands
  INTEGER, INTENT(IN) :: is
  !! spin index (up/down)
  INTEGER, INTENT(IN) :: isk(nks)
  !! which k is up/down
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! Kohn Sham energy [Ry]
  REAL(DP), INTENT(IN) :: nelec
  !! The # of electrons
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! Intetration weight of each k
  REAL(DP), INTENT(INOUT) :: ef
  !! The Fermi energy
  !
  ! ... local variables
  !
  INTEGER :: iter, maxiter=300
  REAL(DP) :: elw, eup, sumkmid, eps=1.0e-10_dp
  !
  ! find bounds for the Fermi energy.
  !
  IF (ntetra == 0 ) &
     CALL errore('opt_tetra_weights:', 'called before initialization', 1) 
  elw = MINVAL(et(1:nbnd,1:nks))
  eup = MAXVAL(et(1:nbnd,1:nks))
  !
  ! Bisection method
  !
  DO iter = 1, maxiter
     !
     ef = (eup + elw) * 0.5_dp
     !
     ! Calc. # of electrons 
     !
     CALL opt_tetra_weights_only( nks, nspin, nbnd, et, ef, &
                                  wg, is, isk )
     !
     IF (is == 0) THEN
        sumkmid = SUM(wg(1:nbnd,1:nks))
     ELSEIF (is == 1) THEN
        sumkmid = SUM(wg(1:nbnd,1:nks / 2))
     ELSEIF (is == 2) THEN
        sumkmid = SUM(wg(1:nbnd,nks / 2 + 1:nks))
     ENDIF
     !
     ! convergence check
     !
     IF (ABS(sumkmid - nelec) < eps) THEN
        EXIT
     ELSEIF (sumkmid < nelec) THEN
        elw = ef
     ELSE
        eup = ef
     ENDIF
     !
  ENDDO ! iter
  !
  IF (iter >= maxiter) CALL errore( "opt_tetra_weights", "Not converged", iter )
  !
END SUBROUTINE opt_tetra_weights
!
!
!--------------------------------------------------------------------------------------
SUBROUTINE opt_tetra_weights_only( nks, nspin, nbnd, et, ef, &
                                   wg, is, isk )
  !------------------------------------------------------------------------------------
  !! Calculate Occupation with given Fermi energy. 
  !! @Note ( P. Delugas 8/10/2019) Needs to be called only after initialization, otherwise
  !!       stops the program with an error call.
  !
  INTEGER, INTENT(IN) :: nks
  !! The total # of k in irr-BZ
  INTEGER, INTENT(IN) :: nspin
  !! The # of spin components
  INTEGER, INTENT(IN) :: nbnd
  !! The # of bands
  INTEGER, INTENT(IN) :: is
  !! spin index (up/down)
  INTEGER, INTENT(IN) :: isk(nks)
  !! which k is up/down
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! Kohn Sham energy [Ry]
  REAL(DP), INTENT(INOUT) :: wg(nbnd,nks)
  !! Intetration weight of each k
  REAL(DP), INTENT(IN) :: ef
  !! The Fermi energy
  !
  ! ... local variables
  !
  INTEGER :: ns, ik, nt, ibnd, jbnd, kbnd, i, ii, itetra(4), nk
  REAL(DP) :: e(4), wg0(4), C(3), a(4,4), wg1
  INTEGER :: nspin_lsda
  !
  ! Zero-clear only "is"th spin component
  !
  DO ik = 1, nks
     IF (is /= 0) THEN
        IF (isk(ik) /= is) CYCLE
     ENDIF
     DO ibnd = 1, nbnd
        wg (ibnd, ik) = 0.d0
     ENDDO
  ENDDO
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  ENDIF
  !
  DO ns = 1, nspin_lsda
     IF (is /= 0) THEN
        IF (ns /= is) CYCLE
     ENDIF
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nks / 2
     ENDIF
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
           ENDDO
           !
           itetra(1) = 0
           CALL hpsort( 4, e, itetra )
           !
           DO ii = 1, 4
              DO i = 1, 4
                 IF ( ABS(e(i)-e(ii)) < 1.d-12 ) THEN
                     a(ii,i) = 0.0_dp
                 ELSE
                     a(ii,i) = ( ef - e(i) ) / (e(ii) - e(i) )
                 END IF
              ENDDO
           ENDDO
           !
           IF( e(1) <= ef .AND. ef < e(2) ) THEN
              !
              C(1) = a(2,1) * a(3,1) * a(4,1) * 0.25_dp
              wg0(1) = C(1) * (1.0_dp + a(1,2) + a(1,3) + a(1,4))
              wg0(2) = C(1) * a(2,1)
              wg0(3) = C(1) * a(3,1)
              wg0(4) = C(1) * a(4,1)
              !
           ELSEIF( e(2) <= ef .AND. ef < e(3)) THEN
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
           ELSEIF ( e(3) <= ef .AND. ef < e(4)) THEN
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
           ELSEIF (e(4) <= ef ) THEN
              !
              wg0(1:4) = 0.25_dp
              !
           ELSE
              !
              wg0(1:4) = 0.0_dp
              !
           ENDIF
           !
           wg0(1:4) = wg0(1:4) / REAL(ntetra, dp)
           !
           DO ii = 1, nntetra
              !
              ik = tetra(ii, nt) + nk
              wg(ibnd,ik) = wg(ibnd,ik) + DOT_PRODUCT(wlsm(itetra(1:4),ii), wg0(1:4))
              !
           ENDDO
           !
        ENDDO ! ibnd
        !
     ENDDO ! nt
     !
  ENDDO
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
           IF (ABS(et(ibnd,ik) - et(jbnd,ik)) < 1e-6_dp) THEN
              wg1 = wg1 + wg(jbnd,ik)
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 wg(kbnd,ik) = wg1 / REAL(jbnd - ibnd, dp)
              ENDDO
              !
              EXIT
           ENDIF
           !
        ENDDO
        !
     ENDDO
  ENDDO
  !
  ! add correct spin normalization (2 for LDA, 1 for all other cases)
  IF ( nspin == 1 ) wg(1:nbnd,1:nks) = wg(1:nbnd,1:nks) * 2.0_dp
  !
END SUBROUTINE opt_tetra_weights_only
!
!
!--------------------------------------------------------------------
SUBROUTINE tetra_dos_t( et, nspin, nbnd, nks, e, dost, dosint )
  !------------------------------------------------------------------
  !
  USE kinds,   ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspin
  !! The # of spin components
  INTEGER, INTENT(IN) :: nbnd
  !! The # of bands
  INTEGER, INTENT(IN) :: nks
  !! The total # of k in irr-BZ
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! Kohn Sham energy [Ry]
  REAL(DP), INTENT(IN) :: e
  REAL(DP), INTENT(OUT):: dost(2)
  REAL(DP), OPTIONAL, INTENT(OUT) :: dosint(2)
  !
  ! ... local variables
  !
  INTEGER :: itetra(4), nk, ns, nt, ibnd, i
  REAL(DP) :: etetra(4), e1, e2, e3, e4
  INTEGER :: nspin0
  !
  IF (nspin==4) THEN
     nspin0=1
  ELSE
     nspin0=nspin
  ENDIF
  !
  DO ns = 1, nspin0
     dost (ns) = 0.d0
     IF (PRESENT(dosint)) dosint(ns) = 0.d0
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns==1) THEN
        nk = 0
     ELSE
        nk = nks / 2
     ENDIF
     !
     DO nt = 1, ntetra
        DO ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           DO i = 1, 4
              etetra(i) = et(ibnd, tetra(i,nt) + nk)
           ENDDO
           itetra(1) = 0
           CALL hpsort(4, etetra, itetra)
           e1 = etetra(1)
           e2 = etetra(2)
           e3 = etetra(3)
           e4 = etetra(4)
           IF (e>=e4 ) THEN 
              IF (PRESENT(dosint)) dosint(ns) = dosint(ns) + (1.d0 / ntetra) 
           ELSEIF (e<e4 .AND. e>=e3) THEN
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
              IF (PRESENT(dosint)) dosint(ns) = dosint(ns) + &
                 (1.d0/ntetra) * ( 1.d0 - ( e4 - e)**3/((e4 - e1)*(e4 - e2)*(e4 - e3))) 
           ELSEIF (e<e3 .AND. e>=e2) THEN
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
              IF (PRESENT(dosint)) dosint(ns) = dosint(ns) + (1.d0/ntetra) / (e3 - e1) / (e4 - e1) * &
                        ( (e2 -e1)**2 + 3.d0 * (e2 - e1) *(e - e2) +3.d0*(e - e2)**2 - & 
                        (e3 - e1 +e4 -e2 ) / (e3 - e2 ) / (e4 - e2) * (e - e2) **3)
           ELSEIF (e<e2 .AND. e>e1) THEN 
              dost (ns) = dost (ns) + (1.d0 / ntetra) * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
              IF ( PRESENT(dosint)) &
                 dosint(ns) = dosint(ns) + (1.d0/ntetra) * ((e -e1)**3)/(e2 -e1)/(e3-e1)/(e4-e1) 
           ENDIF
        ENDDO

     ENDDO
     !
     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 
     !
     IF ( nspin == 1 ) dost(ns) = dost(ns) * 2.d0
     IF (PRESENT(dosint) .AND. nspin==1) dosint(ns) = dosint(ns) * 2.d0 
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE tetra_dos_t
!
!
!--------------------------------------------------------------------
SUBROUTINE opt_tetra_dos_t( et, nspin, nbnd, nks, e, dost, dosint )
  !------------------------------------------------------------------
  !! Optimized tetrahedron method.
  !
  IMPLICIT NONE
  !
  INTEGER :: nspin, nbnd, nks
  REAL(DP) :: et(nbnd,nks), e, dost(2)
  REAL(DP), OPTIONAL :: dosint(2) 
  !
  ! ... local variables
  !
  INTEGER :: itetra(4), nk, ns, nt, ibnd, i
  !
  REAL(DP) :: etetra(4), e1, e2, e3, e4
  INTEGER :: nspin0
  !
  IF (nspin==4) THEN
     nspin0=1
  ELSE 
     nspin0=nspin
  ENDIF
  !
  DO ns = 1, nspin0
     dost (ns) = 0.d0
     IF (PRESENT(dosint)) dosint(ns) = 0.d0
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns==1) THEN
        nk = 0
     ELSE
        nk = nks / 2
     ENDIF
     !
     DO nt = 1, ntetra
        DO ibnd = 1, nbnd
           ! these are the energies at the vertexes of the nt-th tetrahedron
           etetra(1:4) = 0.0_dp
           DO i = 1, nntetra
              etetra(1:4) = etetra(1:4) + wlsm(1:4,i) * et(ibnd,tetra(i, nt) + nk)
           ENDDO
           itetra (1) = 0
           CALL hpsort( 4, etetra, itetra )
           e1 = etetra(1)
           e2 = etetra(2)
           e3 = etetra(3)
           e4 = etetra(4)
           IF (e >= e4 ) THEN
              IF (PRESENT(dosint)) dosint(ns) = dosint(ns) + (1.d0 / ntetra) 
           ELSEIF (e<e4 .AND. e>=e3) THEN
              dost(ns) = dost(ns) + 1.d0 / ntetra * (3.0d0 * (e4 - e) **2 / &
                         (e4 - e1) / (e4 - e2) / (e4 - e3) )
              IF (PRESENT(dosint)) dosint(ns) = dosint(ns) + &
                 (1.d0/ntetra) * ( 1.d0 - ( e4 - e)**3/((e4 - e1)*(e4 - e2)*(e4 - e3))) 
           ELSEIF (e<e3 .AND. e>=e2) THEN
              dost(ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                         * (3.0d0 * (e2 - e1) + 6.0d0 * (e-e2) - 3.0d0 * (e3 - e1 + e4 - e2) &
                         / (e3 - e2) / (e4 - e2) * (e-e2) **2)
              IF (PRESENT(dosint)) dosint(ns) = dosint(ns) + (1.d0/ntetra) / (e3 - e1) / (e4 - e1) * &
                        ( (e2 -e1)**2 + 3.d0 * (e2 - e1) *(e - e2) +3.d0*(e - e2)**2 - & 
                        (e3 - e1 +e4 -e2 ) / (e3 - e2 ) / (e4 - e2) * (e - e2) **3)
           ELSEIF (e<e2 .AND. e>e1) THEN
              dost(ns) = dost (ns) + 1.d0 / ntetra * 3.0d0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
              IF ( PRESENT(dosint)) &
                 dosint(ns) = dosint(ns) + (1.d0/ntetra) * ((e -e1)**3)/(e2 -e1)/(e3-e1)/(e4-e1) 
           ENDIF
        ENDDO
        !
     ENDDO
     !
     ! add correct spin normalization : 2 for LDA, 1 for LSDA or
     ! noncollinear calculations 
     !
     IF ( nspin == 1 ) THEN
        dost(ns) = dost(ns) * 2.d0
        IF (PRESENT(dosint)) dosint(ns) = dosint(ns) * 2.d0 
     ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE opt_tetra_dos_t
!
!
!-----------------------------------------------------------------------------
SUBROUTINE opt_tetra_partialdos( nspin0, kresolveddos, ne, natomwfc, nkseff, &
                                 Emin, DeltaE, proj, pdos, dostot, nspindos )
  !-------------------------------------------------------------------------
  !! Compute partial dos.
  !
  USE lsda_mod,           ONLY : nspin, isk
  USE wvfct,              ONLY : et, nbnd
  USE klist,              ONLY : nkstot, wk
  USE spin_orb,           ONLY : lspinorb
  USE constants,          ONLY : rytoev
  !
  IMPLICIT NONE
  !
  LOGICAL :: kresolveddos
  !
  INTEGER,INTENT(IN) :: ne, natomwfc, nkseff, nspin0, nspindos
  !
  REAL(DP),INTENT(IN) :: Emin, deltaE, proj(natomwfc, nbnd, nkstot)
  REAL(DP),INTENT(OUT) :: pdos(0:ne,natomwfc,nspin0,nkseff), &
                          dostot(0:ne,nspindos,nkseff)
  !
  ! ... local variables
  !
  INTEGER :: ns, nk, nt, ibnd, itetra(4), ii, ik, ie, nwfc, nspin1
  REAL(DP) :: etetra(4), e0, e1, e2, e3, e4, wgt(4), G, spindeg,  &
              f12, f13, f14, f21, f23, f24, f31, f32, f34, f41, f42
  !
  IF(nspin == 2) THEN
     nspin1 = nspin
  ELSE
     nspin1 = 1
  ENDIF
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
           CALL hpsort( 4, etetra, itetra )
           e1 = etetra(1)
           e2 = etetra(2)
           e3 = etetra(3)
           e4 = etetra(4)
           !
           DO ie = 0, ne
              !
              e0 = Emin + DeltaE * REAL(ie, dp)
              !
              IF ( e3 < e0 .AND. e0 < e4) THEN
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
              ELSEIF ( e2 < e0 .AND. e0 < e3 ) THEN
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
              ELSEIF ( e1 < e0 .AND. e0 < e2 ) THEN
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
              ENDIF
              !
              IF (kresolveddos) THEN
                 !
                 DO ii = 1, nntetra
                    !
                    ik = tetra(ii, nt) + nk
                    DO nwfc = 1, natomwfc
                       pdos(ie,nwfc,ns,ik - nk) = pdos(ie,nwfc,ns,ik - nk) &
                         + proj (nwfc, ibnd, ik) &
                         * DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G / wk(ik)
                    ENDDO
                    !
                    IF(nspindos == 1) THEN
                       dostot(ie,1,ik - nk) = dostot(ie,1,ik - nk) &
                         + DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G / wk(ik)
                    ELSE
                       dostot(ie,ns,ik - nk) = dostot(ie,ns,ik - nk) &
                         + DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G / wk(ik)
                    ENDIF
                    !
                 ENDDO
                 !
              ELSE
                 !
                 DO ii = 1, nntetra
                    !
                    ik = tetra(ii, nt) + nk
                    DO nwfc = 1, natomwfc
                       pdos(ie,nwfc,ns,1) = pdos(ie,nwfc,ns,1) &
                         + proj (nwfc, ibnd, ik) * DOT_PRODUCT(wlsm(itetra(1:4),ii), wgt(1:4)) * G
                    ENDDO
                    !
                 ENDDO
                 !
                 IF (nspindos == 1) THEN
                    dostot(ie,1,1) = dostot(ie,1,1) + SUM(wgt(1:4)) * G
                 ELSE
                    dostot(ie,ns,1) = dostot(ie,ns,1) + SUM(wgt(1:4)) * G
                 ENDIF
                 !
              ENDIF
              !
           ENDDO
           !
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  spindeg= 1.0_dp
  IF (nspin == 1) spindeg= 2.0_dp
  pdos(0:ne,1:natomwfc,1:nspin0,1:nkseff) = pdos(0:ne,1:natomwfc,1:nspin0,1:nkseff) &
      * spindeg / (rytoev * REAL(ntetra,dp))
  dostot(0:ne,1:nspindos,1:nkseff) = dostot(0:ne,1:nspindos,1:nkseff) &
      * spindeg / (rytoev * REAL(ntetra,dp))
  !
END SUBROUTINE opt_tetra_partialdos
!
!-----------------------------------------------------------
SUBROUTINE deallocate_tetra( )
   !--------------------------------------------------------
   !! Deallocate tetra and wlsm
   !
   IF ( ALLOCATED(tetra) ) DEALLOCATE (tetra)
   IF ( ALLOCATED(wlsm ) ) DEALLOCATE (wlsm )
   !
END SUBROUTINE deallocate_tetra
!
!
END MODULE ktetra
