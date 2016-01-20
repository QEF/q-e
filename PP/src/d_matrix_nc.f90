!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE d_matrix_nc (dy012, dy112, dy212, dy312)
  !---------------------------------------------------------------
  !
  ! Provides symmetry operations in the (l, s) subspaces for l=0,1,2,3
  !
  USE kinds, ONLY: DP
  USE symm_base, ONLY:  nsym, sr
  USE random_numbers, ONLY : randy
  USE matrix_inversion
  !
  IMPLICIT NONE
  !
  !     first the input/output variables
  !
  COMPLEX(DP) :: dy012 (2, 2, 48), dy112 (6, 6, 48), dy212 (10, 10, 48), &
                      dy312 (14, 14, 48)
  !
  ! output: symmetry matrices in the l=0, l=1, l=2 and l=3 subspace resp.
  !
  !    here the local parameters
  !
  INTEGER, PARAMETER :: maxl = 3, maxm = 2*maxl+1, &
       maxlm = (maxl+1)*(maxl+1)
  !
  ! maxl = max value of l allowed
  ! maxm = number of m components for l=maxl
  ! maxlm= number of l,m spherical harmonics for l <= maxl
  !
  REAL(DP), PARAMETER :: eps = 1.0d-9
  !
  !   and the local variables
  !
  INTEGER :: m, n, m1, n1, ipol, isym
  INTEGER :: l, n2, ind, ind1, ind2
  REAL(DP) :: ylm(maxm, maxlm), ylms(maxm, maxlm), &
       rl(3,maxm), rrl (maxm), srl(3,maxm), capel
  REAL(DP) :: yl1 (3, 3), yl2(5, 5), yl3(7,7), yl1_inv (3, 3), &
                   yl2_inv(5, 5),  yl3_inv(7, 7), &
                   dy1 (3, 3, 48), dy2 (5, 5, 48), dy3 (7, 7, 48)
  COMPLEX(DP) :: dy012_con(2,2), dy112_con(6,6), dy212_con(10,10), &
                      dy312_con(14,14), s_spin(2,2,48), delta(14,14)
  COMPLEX(DP), EXTERNAL :: ZDOTU
  !
  !    Here we find the true symmetries of the crystal
  !
  CALL spinsym( s_spin )
  !
  !  randomly distributed points on a sphere
  !
  DO m = 1, maxm
     rl (1, m) = randy () - 0.5d0
     rl (2, m) = randy () - 0.5d0
     rl (3, m) = randy () - 0.5d0
     rrl (m) = rl (1,m)**2 + rl (2,m)**2 + rl (3,m)**2
  ENDDO
  CALL ylmr2 ( maxlm, 2*maxl+1, rl, rrl, ylm )
  !
  !  invert Yl for each block of definite l (note the transpose operation)
  !
  !  l = 1 block
  !
  DO m = 1, 3
     DO n = 1, 3
        yl1 (m, n) = ylm (n, 1+m)
     ENDDO
  ENDDO
  CALL invmat (3, yl1, yl1_inv)
  !
  !  l = 2 block
  !
  DO m = 1, 5
     DO n = 1, 5
        yl2 (m, n) = ylm (n, 4+m)
     ENDDO
  ENDDO
  CALL invmat (5, yl2, yl2_inv)
  !
  !  l = 3 block
  !
  DO m = 1, 7
     DO n = 1, 7
        yl3 (m, n) = ylm (n, 9+m)
     ENDDO
  ENDDO
  CALL invmat (7, yl3, yl3_inv )
  !
  ! now for each symmetry operation of the point-group ...
  !
  DO isym = 1, nsym
     !
     ! srl(:,m) = rotated rl(:,m) vectors
     !
     srl = matmul (sr(:,:,isym), rl)
     !
     CALL ylmr2 ( maxlm, maxm, srl, rrl, ylms )
     !
     !  find  D_S = Yl_S * Yl_inv (again, beware the transpose)
     !  and the rotation matrices for the | l m s s_z > basis
     !
     !  l = 0 block
     !
     DO m1 = 1, 2
        DO n1 = 1, 2
           dy012 (m1, n1, isym)= conjg( s_spin (n1, m1, isym) )
        ENDDO
     ENDDO
     !
     !  l = 1 block
     !
     DO m = 1, 3
        DO n = 1, 3
           yl1 (m, n) = ylms (n, 1+m)
        ENDDO
     ENDDO
     dy1 (:, :, isym) = matmul (yl1(:,:), yl1_inv(:,:))
     DO m = 1, 3
        DO n = 1, 3
           DO m1 = 1, 2
              DO n1 = 1, 2
                 dy112 (m+3*(m1-1), n+3*(n1-1), isym) = &
                          cmplx(dy1 (m, n, isym), 0.d0,kind=DP) * s_spin (m1, n1, isym)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
     !  l = 2 block
     !
     DO m = 1, 5
        DO n = 1, 5
           yl2 (m, n) = ylms (n, 4+m)
        ENDDO
     ENDDO
     dy2 (:, :, isym) = matmul (yl2(:,:), yl2_inv(:,:))
     DO m = 1, 5
        DO n = 1, 5
           DO m1 = 1, 2
              DO n1 = 1, 2
                 dy212 (m+5*(m1-1), n+5*(n1-1), isym) =   &
                     cmplx(dy2 (m, n, isym), 0.d0,kind=DP) * s_spin (m1, n1, isym)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
     !  l = 3 block
     !
     DO m = 1, 7
        DO n = 1, 7
           yl3 (m, n) = ylms (n, 9+m)
        ENDDO
     ENDDO
     dy3 (:, :, isym) = matmul (yl3(:,:), yl3_inv(:,:))
     DO m = 1, 7
        DO n = 1, 7
           DO m1 = 1, 2
              DO n1 = 1, 2
                 dy312 (m+7*(m1-1), n+7*(n1-1), isym) =  &
                    cmplx(dy3 (m, n, isym), 0.d0,kind=DP) * s_spin (m1, n1, isym)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO
  !
  ! check that D_S matrices are unitary as they should
  !
  delta (:,:) = (0.d0,0.d0)
  DO m= 1, 14
     delta(m,m) = (1.d0,0.d0)
  ENDDO
  DO isym =1,nsym
     !
     !  l = 0 block
     !
     capel = 0.d0
     dy012_con(:,:) = conjg( dy012(:,:,isym) )
     DO m = 1, 2
        DO n = 1, 2
           capel = capel +  &
           abs(ZDOTU(2,dy012_con(1,m),1,dy012(1,n,isym),1)-delta(m,n))**2
        ENDDO
     ENDDO
     IF (capel>eps) CALL errore ('d_matrix_nc', &
        'D_S (l=0) for this symmetry operation is not unitary',isym)
     !
     !  l = 1 block
     !
     capel = 0.d0
     dy112_con(:,:) = conjg( dy112(:,:,isym) )
     DO m = 1, 6
        DO n = 1, 6
           capel = capel +  &
           abs(ZDOTU(6, dy112_con(1,m), 1, dy112(1,n,isym), 1)-delta(m,n))**2
        ENDDO
     ENDDO
     IF (capel>eps) CALL errore ('d_matrix_nc', &
        'D_S (l=1) for this symmetry operation is not unitary',isym)
     !
     !  l = 2 block
     !
     capel = 0.d0
     dy212_con(:,:)=conjg(dy212(:,:,isym))
     DO m = 1, 10
        DO n = 1, 10
           capel = capel +  &
           abs(ZDOTU(10, dy212_con(1,m), 1, dy212(1,n,isym), 1)-delta(m,n))**2
        ENDDO
     ENDDO
     IF (capel>eps) CALL errore ('d_matrix_nc', &
        'D_S (l=2) for this symmetry operation is not unitary',isym)
     !
     !  l = 3 block
     !
     capel = 0.d0
     dy312_con(:,:)=conjg(dy312(:,:,isym))
     DO m = 1, 14
        DO n = 1, 14
           capel = capel +  &
           abs(ZDOTU(14, dy312_con(1,m), 1, dy312(1,n,isym), 1)-delta(m,n))**2
        ENDDO
     ENDDO
     IF (capel>eps) CALL errore ('d_matrix_nc', &
        'D_S (l=3) for this symmetry operation is not unitary',isym)
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE d_matrix_nc
