!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine d_matrix (dy1, dy2)
  !---------------------------------------------------------------
  !
  use pwcom
  implicit none
  real(kind=DP) :: dy1 (3, 3, 48), dy2 (5, 5, 48)
  !
  integer, parameter :: maxl = 2, maxm = 2*maxl+1, &
       maxlm = (maxl+1)*(maxl+1)
  ! maxl = max value of l allowed
  ! maxm = number of m components for l=maxl
  ! maxlm= number of l,m spherical harmonics for l <= maxl
  integer :: m, n, ipol, isym
  real(kind=DP) :: ylm(maxm, maxlm),  yl1 (3, 3), yl2(5, 5), &
       yl1_inv (3, 3), yl2_inv(5, 5), ylms(maxm, maxlm),     &
       rl(3,maxm), rrl (maxm), srl(3,maxm), scart (3, 3)
  real(kind=DP), external :: rndm
  !
  !  randomly distributed points on a sphere
  !
  do m = 1, maxm
     rl (1, m) = rndm () - 0.5
     rl (2, m) = rndm () - 0.5
     rl (3, m) = rndm () - 0.5
     rrl (m) = rl (1,m)**2 + rl (2,m)**2 + rl (3,m)**2
  enddo
  call ylmr2 ( maxlm, 2*maxl+1, rl, rrl, ylm )
  !
  !  invert Yl for each block of definite l (note the transpose operation)
  !
  !  l = 1 block
  !
  do m = 1, 3
     do n = 1, 3
        yl1 (m, n) = ylm (n, 1+m)
     end do
  end do
  call invmat (yl1, yl1_inv, 3)
  !
  !  l = 2 block
  !
  do m = 1, 5
     do n = 1, 5
        yl2 (m, n) = ylm (n, 4+m)
     end do
  end do
  call invmat (yl2, yl2_inv, 5)
  !
  ! now for each symmetry operation of the point-group ...
  !
  do isym = 1, nsym
     !
     ! scart = symmetry operation in cartesian axis
     ! srl(:,m) = rotated rl(:,m) vectors
     !
     call s_axis_to_cart (s (1, 1, isym), scart, at, bg)
     srl = matmul (scart, rl)
     !
     call ylmr2 ( maxlm, maxm, srl, rrl, ylms )
     !
     !  find  D_S = Yl_S * Yl_inv (again, beware the transpose)
     !
     !  l = 1
     !
     do m = 1, 3
        do n = 1, 3
           yl1 (m, n) = ylms (n, 1+m)
        end do
     end do
     dy1 (:, :, isym) = matmul (yl1(:,:), yl1_inv(:,:))
     !
     !  l = 2 block
     !
     do m = 1, 5
        do n = 1, 5
           yl2 (m, n) = ylms (n, 4+m)
        end do
     end do
     dy2 (:, :, isym) = matmul (yl2(:,:), yl2_inv(:,:))
     !
  enddo
  return

end subroutine d_matrix
