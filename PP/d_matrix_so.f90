!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!---------------------------------------------------------------
SUBROUTINE d_matrix_so (dyj12, dyj32, dyj52, dyj72)  
  !---------------------------------------------------------------
  !
  ! Provides symmetry operations in the j=1/2, j=3/2, j=5/2 and j=7/2 
  ! subspaces 
  !
  USE kinds, only: DP
  USE cell_base, ONLY : at, bg, ibrav, symm_type
  USE symme, ONLY:  nsym, s
  USE spin_orb,   ONLY : rot_ylm
  USE random_numbers, ONLY : rndm
  !
  IMPLICIT NONE
  !
  !     first the input/output variables
  !
  COMPLEX(DP) :: dyj12 (2, 2, 48), dyj32 (4, 4, 48), dyj52 (6, 6, 48), &
                      dyj72 (8, 8, 48)
  !
  ! output: symmetry matrices in the j=1/2, j=3/2, j=5/2 and j=7/2 subspace 
  !         respectively
  ! 
  !    here the local parameters
  !
  INTEGER, PARAMETER :: maxl = 3, maxm = 2*maxl+1, &
       maxlm = (maxl+1)*(maxl+1)
  INTEGER, PARAMETER :: maxmj = 2*maxl+2
  !
  ! maxl = max value of l allowed
  ! maxm = number of m components for l=maxl
  ! maxlm= number of l,m spherical harmonics for l <= maxl
  ! maxmj = number of mj components
  !
  REAL(DP), PARAMETER :: eps = 1.0d-9
  !
  !   and the local variables
  !  
  INTEGER :: m, n, m1, n1, ipol, isym
  INTEGER :: l, n2, ind, ind1, ind2
  REAL(DP) :: j, ylm(maxm, maxlm), ylms(maxm, maxlm), &
       rl(3,maxm), rrl (maxm), srl(3,maxm), scart (3, 3), capel
  REAL(DP) :: Ulall(maxl,maxmj+(maxmj-1),maxmj+(maxmj-1)), spinor, &
                   Ul1(6,6), Ul1_inv(6,6), Ul3(14,14), Ul3_inv(14,14)
  COMPLEX(DP) :: dy1 (3, 3, 48), dy2 (5, 5, 48), & 
                      dy3 (7, 7, 48), dy112 (6, 6, 48), &
                      dy212 (10, 10, 48), dy312 (14, 14, 48), &
                      d12_con(2,2), d32_con(4,4), d52_con(6,6), d72_con(8,8), &
                      s_spin(2,2,48), delta(8,8)
  COMPLEX(DP) :: ylm_compl(maxm, maxlm), ylms_compl(maxm, maxlm), &
                      yl1 (3, 3), yl2(5, 5), yl3(7,7), &
                      yl1_inv (3, 3), yl2_inv(5, 5),  yl3_inv(7, 7),  &
                      Ul1C(6,6), Ul1C_inv(6,6), Ul3C(14,14), Ul3C_inv(14,14)
  COMPLEX(DP), EXTERNAL :: ZDOTU
  !
  !    Here we find the true symmetries of the crystal 
  !
  IF ( ibrav == 4 .OR. ibrav == 5 ) THEN
     !
     ! ... here the hexagonal or trigonal bravais lattice
     !
     CALL hexspinsym( s_spin )
     !
  ELSE IF ( ibrav >=1  .AND. ibrav <= 14 ) THEN
     !
     ! ... here for the cubic bravais lattice
     !
     CALL cubicspinsym( s_spin )
     !
  ELSE IF ( ibrav == 0 ) THEN
     !
     IF ( symm_type == 'cubic' ) THEN
        !
        CALL cubicspinsym( s_spin )
        !
     ELSE IF ( symm_type == 'hexagonal' ) THEN
        !
        CALL hexspinsym( s_spin )
        !
     END IF
     !
  ELSE
     !
     CALL errore( 'd_matrix_so', 'wrong ibrav', 1 )
     !
  END IF
  !
  ! Transformation matrices from the | l m s s_z > basis to the 
  ! | j mj l s > basis in the l-subspace
  !
  Ulall (:,:,:) = 0.d0
  DO l = 1, 3
     j = DBLE(l) - 0.5d0
     DO m1= 1, 2*l
        m= m1 - l
        Ulall (l,m1,2*(m1-1)+1) = spinor (l,j,m,1)
        Ulall (l,m1,2*(m1-1)+4) = spinor (l,j,m,2)
     END DO
     j = DBLE(l) + 0.5d0 
     DO m1= 1, 2*l + 2
        m = m1 - l - 2
        IF (m1 == 1) THEN
           Ulall (l,m1+2*l,2*(m1-1)+2) = spinor (l,j,m,2)
        ELSE IF (m1==2*l+2) THEN
           Ulall (l,m1+2*l,2*(m1-1)-1) = spinor (l,j,m,1)
        ELSE 
           Ulall (l,m1+2*l,2*(m1-1)-1) = spinor (l,j,m,1)
           Ulall (l,m1+2*l,2*(m1-1)+2) = spinor (l,j,m,2)
        END IF
     END DO
  END DO

  !
  !  invert Ulall for l = 1 and l = 3 blocks
  !   (we don't need to invert l = 2 block, see below)
  !
  !  l = 1 block
  !
  DO m = 1, 6
     DO n = 1, 6
        Ul1 (m, n) = Ulall (1, m, n)
     END DO
  END DO
  CALL invmat (6, Ul1, Ul1_inv, capel)
  Ul1C (:,:) = CMPLX(Ul1 (:,:), 0.d0) 
  Ul1C_inv (:,:) = CMPLX(Ul1_inv (:,:), 0.d0)
  !
  !  l = 3 block
  !
  DO m = 1, 14
     DO n = 1, 14
        Ul3 (m, n) = Ulall (3, m, n)
     END DO
  END DO
  CALL invmat (14, Ul3, Ul3_inv, capel)
  Ul3C (:,:) = CMPLX(Ul3 (:,:), 0.d0) 
  Ul3C_inv (:,:) = CMPLX(Ul3_inv (:,:), 0.d0)
  ! 
  !  randomly distributed points on a sphere
  !
  DO m = 1, maxm
     rl (1, m) = rndm () - 0.5d0
     rl (2, m) = rndm () - 0.5d0
     rl (3, m) = rndm () - 0.5d0
     rrl (m) = rl (1,m)**2 + rl (2,m)**2 + rl (3,m)**2
  END DO
  CALL ylmr2 ( maxlm, 2*maxl+1, rl, rrl, ylm )
  !
  ! calculate complex spherical harmonics
  !
  ylm_compl = (0.d0,0.d0)
  DO l = 1, maxl
     DO m = -l, l
        ind = maxl + 1 + m
        ind1 = l**2 + l + 1 + m 
        DO n2 = 1, 2*l+1
           ind2 = l**2 + n2
           IF (ABS(rot_ylm(ind,n2)).gt.1.d-8) &
           ylm_compl(:,ind1) = ylm_compl(:,ind1) + rot_ylm(ind,n2)*ylm(:,ind2)
        END DO
     END DO
  END DO
  !
  !  invert Yl for each block of definite l (note the transpose operation)
  !
  !  l = 1 block
  !
  DO m = 1, 3
     DO n = 1, 3
        yl1 (m, n) = ylm_compl (n, 1+m)
     END DO
  END DO
  CALL invmat_complex (3, yl1, yl1_inv, capel)
  !
  !  l = 2 block
  !
  DO m = 1, 5
     DO n = 1, 5
        yl2 (m, n) = ylm_compl (n, 4+m)
     END DO
  END DO
  CALL invmat_complex (5, yl2, yl2_inv, capel)
  !
  !  l = 3 block
  !
  DO m = 1, 7
     DO n = 1, 7
        yl3 (m, n) = ylm_compl (n, 9+m)
     END DO
  END DO
  CALL invmat_complex (7, yl3, yl3_inv, capel)
  !
  ! now for each symmetry operation of the point-group ...
  !
  DO isym = 1, nsym
     !
     ! scart = symmetry operation in cartesian axis
     ! srl(:,m) = rotated rl(:,m) vectors
     !
     CALL s_axis_to_cart (s (1, 1, isym), scart, at, bg)
     srl = matmul (scart, rl)
     !
     CALL ylmr2 ( maxlm, maxm, srl, rrl, ylms )
     !
     !  find  D_S = Yl_S * Yl_inv (again, beware the transpose)
     !  and the rotation matrices for the | j mj l s > basis
     !
     ylms_compl = (0.d0,0.d0)
     DO l = 1, maxl
        DO m = -l, l
           ind = maxl + 1 + m
           ind1 = l**2 + l + 1 + m
           DO n2 = 1, 2*l+1
              ind2 = l**2 + n2
              IF (ABS(rot_ylm(ind,n2)).gt.1.d-8) &
              ylms_compl (:,ind1) = ylms_compl (:,ind1) &
                                 + rot_ylm (ind,n2) * ylms (:,ind2)
           END DO
        END DO
     END DO
     !
     !  l = 1 block   =>  j=1/2  and  j=3/2  
     !
     DO m = 1, 3
        DO n = 1, 3
           yl1 (m, n) = ylms_compl (n, 1+m)
        END DO
     END DO
     dy1 (:, :, isym) = matmul (yl1 (:,:), yl1_inv (:,:))
     DO m = 1, 3
        DO n = 1, 3
           DO m1 = 1, 2
              DO n1 = 1, 2
                 dy112 (2*(m-1)+m1, 2*(n-1)+n1, isym) = dy1 (m, n, isym) &
                                            * s_spin (m1, n1, isym)
              END DO
           END DO
        END DO
     END DO
     dy112 (:, :, isym) = matmul (dy112 (:,:,isym), Ul1C_inv (:,:))
     dy112 (:, :, isym) = matmul (Ul1C (:,:), dy112 (:,:,isym))
     DO m = 1, 2
        DO n = 1, 2
           dyj12 (m, n, isym) = CONJG(dy112 (n, m, isym))
        END DO
     END DO
     DO m = 1, 4
        DO n = 1, 4
           dyj32 (m ,n, isym) = CONJG(dy112 (2+n, 2+m, isym))
        END DO
     END DO
     !
     !  l = 3 block  =>  j=5/2  and  j=7/2 
     !
     DO m = 1, 7
        DO n = 1, 7
           yl3 (m, n) = ylms_compl (n, 9+m)
        END DO
     END DO
     dy3 (:, :, isym) = matmul (yl3 (:,:), yl3_inv (:,:))
     DO m = 1, 7
        DO n = 1, 7
           DO m1 = 1, 2
              DO n1 = 1, 2
                 dy312 (2*(m-1)+m1, 2*(n-1)+n1, isym) = dy3 (m, n, isym) &
                                            * s_spin (m1, n1, isym)
              END DO
           END DO
        END DO
     END DO
     dy312 (:, :, isym) = matmul (dy312 (:,:,isym), Ul3C_inv (:,:))
     dy312 (:, :, isym) = matmul (Ul3C (:,:), dy312 (:,:,isym))
     DO m = 1, 6
        DO n = 1, 6
           dyj52 (m, n, isym) = CONJG(dy312 (n, m, isym))
        END DO
     END DO
     DO m = 1, 8
        DO n = 1, 8
           dyj72 (m ,n, isym) = CONJG(dy312 (6+n, 6+m, isym))
        END DO
     END DO
     !
  END DO
  !
  ! check that D_S matrices are unitary as they should 
  !
  delta (:,:) = (0.d0,0.d0)
  DO m= 1, 8
     delta(m,m) = (1.d0,0.d0)
  END DO
  DO isym =1,nsym
     !
     !  j = 1/2 block
     !
     capel = 0.d0
     d12_con(:,:) = CONJG(dyj12 (:,:,isym))
     DO m = 1, 2
        DO n = 1, 2
           capel = capel +  &
           ABS( ZDOTU(2, d12_con(1,m), 1, dyj12(1,n,isym), 1) - delta(m,n) )**2
        END DO
     END DO
     IF (capel.gt.eps) CALL errore ('d_matrix_so', &
        'D_S (j=1/2) for this symmetry operation is not unitary',isym)
     !
     !  j = 3/2 block
     !
     capel = 0.d0
     d32_con(:,:) = CONJG(dyj32 (:,:,isym))
     DO m = 1, 4
        DO n = 1, 4
           capel = capel +  &
           ABS( ZDOTU(4, d32_con(1,m), 1, dyj32(1,n,isym), 1) - delta(m,n) )**2
        END DO
     END DO
     IF (capel.gt.eps) CALL errore ('d_matrix_so', &
        'D_S (j=3/2) for this symmetry operation is not unitary',isym)
     !
     !  j = 5/2 block
     !
     capel = 0.d0
     d52_con(:,:) = CONJG(dyj52 (:,:,isym))
     DO m = 1, 6
        DO n = 1, 6
           capel = capel +  &
           ABS( ZDOTU(6, d52_con(1,m), 1, dyj52(1,n,isym), 1) - delta(m,n) )**2
        END DO
     END DO
     IF (capel.gt.eps) CALL errore ('d_matrix_so', &
        'D_S (j=5/2) for this symmetry operation is not unitary',isym)
     !
     !  j = 7/2 block
     !
     capel = 0.d0
     d72_con(:,:) = CONJG(dyj72 (:,:,isym))
     DO m = 1, 8
        DO n = 1, 8
           capel = capel +  &
           ABS( ZDOTU(8, d72_con(1,m), 1, dyj72(1,n,isym), 1) - delta(m,n) )**2
        END DO
     END DO
     IF (capel.gt.eps) CALL errore ('d_matrix_so', &
        'D_S (j=7/2) for this symmetry operation is not unitary',isym)
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE d_matrix_so
