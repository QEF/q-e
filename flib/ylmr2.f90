!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define NEWYLMR
#ifdef NEWYLMR
!-----------------------------------------------------------------------
subroutine ylmr2 (lmax2, ng, g, gg, ylm)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical 
  !     Recipes but avoiding the calculation of factorials that generate 
  !     overflow for lmax > 11
  !
  USE kinds
  USE constants, ONLY : pi, fpi
  implicit none
  !
  ! Input
  !
  integer :: lmax2, ng
  real(DP) :: g (3, ng), gg (ng)
  !
  ! Output
  !
  real(DP) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  real(DP), allocatable :: cost (:), sent(:), phi (:), Q(:,:,:)
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm
  !
  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, 25
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  stop 'l > 25 or wrong number of Ylm required'
10 continue

  !
  if (lmax == 0) then
     ylm(:,1) =  sqrt (1.d0 / fpi)
     return
  end if
  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  allocate(cost(ng), sent(ng), phi(ng), Q(ng,0:lmax,0:lmax) )
  do ig = 1, ng
     gmod = sqrt (gg (ig) )
     if (gmod < eps) then
        cost(ig) = 0.d0
     else
        cost(ig) = g(3,ig)/gmod
     endif
     !
     !  beware the arc tan, it is defined modulo pi
     !
     if (g(1,ig) > eps) then
        phi (ig) = atan( g(2,ig)/g(1,ig) )
     else if (g(1,ig) < -1.d-9) then
        phi (ig) = atan( g(2,ig)/g(1,ig) ) + pi
     else
        phi (ig) = sign( pi/2.d0,g(2,ig) )
     end if
  enddo
  sent(:) = sqrt(max(0d0,1.d0-cost(:)**2))
  !
  !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
  !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  !
  lm = 0
  do l = 0, lmax
     c = sqrt (DBLE(2*l+1) / fpi)
     if ( l == 0 ) then
        Q (:,0,0) = 1.d0
     else if ( l == 1 ) then
        Q (:,1,0) = cost(:)
        Q (:,1,1) =-sent(:)/sqrt(2.d0)
     else
        !
        !  recursion on l for Q(:,l,m)
        !
        do m = 0, l - 2
           Q(:,l,m) = cost(:)*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(:,l-1,m) &
                    - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(:,l-2,m)
        end do
        Q(:,l,l-1) = cost(:) * sqrt(DBLE(2*l-1)) * Q(:,l-1,l-1)
        Q(:,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent(:)*Q(:,l-1,l-1) 
     end if
     !
     ! Y_lm, m = 0
     !
     lm = lm + 1
     ylm(:, lm) = c * Q(:,l,0)
     !
     do m = 1, l
        !
        ! Y_lm, m > 0
        !
        lm = lm + 1
        ylm(:, lm) = c * sqrt(2.d0) * Q(:,l,m) * cos (m*phi(:))
        !
        ! Y_lm, m < 0
        !
        lm = lm + 1
        ylm(:, lm) = c * sqrt(2.d0) * Q(:,l,m) * sin (m*phi(:))
     end do
  end do
  !
  deallocate(cost, sent, phi, Q)
  !
  return
end subroutine ylmr2

#else
!-----------------------------------------------------------------------
subroutine ylmr2 (lmax2, ng, g, gg, ylm)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm as given in Numerical Recipes
  !
  USE kinds
  USE constants, ONLY : pi, fpi
  implicit none
  !
  ! Input
  !
  integer :: lmax2, ng
  real(DP) :: g (3, ng), gg (ng)
  !
  ! Output
  !
  real(DP) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  real(DP), allocatable :: cost (:), phi (:), P(:,:,:)
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm
  integer, external:: fact, semifact
  !
  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, 6
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  call errore (' ylmr', 'l > 6 or wrong number of Ylm required',lmax2)
10 continue
  !
  if (lmax == 0) then
     ylm (:,1) =  sqrt (1.d0 / fpi)
     return
  end if
  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  allocate(cost(ng), phi(ng), P(ng,0:lmax,0:lmax) )
  do ig = 1, ng
     gmod = sqrt (gg (ig) )
     if (gmod < eps) then
        cost(ig) = 0.d0
     else
        cost(ig) = g(3,ig)/gmod
     endif
     !
     !  beware the arc tan, it is defined modulo pi
     !
     if (g(1,ig) > eps) then
        phi (ig) = atan( g(2,ig)/g(1,ig) )
     else if (g(1,ig) < -1.d-9) then
        phi (ig) = atan( g(2,ig)/g(1,ig) ) + pi
     else
        phi (ig) = sign( pi/2.d0,g(2,ig) )
     end if
  enddo
  !
  !  P(:,l,m), are the Legendre Polynomials (0 <= m <= l)
  !
  lm = 0
  do l = 0, lmax
     c = sqrt (DBLE(2*l+1) / fpi)
     if ( l == 0 ) then
        P (:,0,0) = 1.d0
     else if ( l == 1 ) then
        P (:,1,0) = cost(:)
        P (:,1,1) = -sqrt(max(0d0,1.d0-cost(:)**2))
     else
        !
        !  recursion on l for P(:,l,m)
        !
        do m = 0, l - 2
           P(:,l,m) = (cost(:)*(2*l-1)*P(:,l-1,m) - (l+m-1)*P(:,l-2,m))/(l-m)
        end do
        P(:,l,l-1) = cost(:) * (2*l-1) * P(:,l-1,l-1)
        P(:,l,l)   = (-1)**l * semifact(2*l-1) * (max(0.d0,1.d0-cost(:)**2))**(DBLE(l)/2)
     end if
     !
     ! Y_lm, m = 0
     !
     lm = lm + 1
     ylm (:, lm) = c * P(:,l,0)
     !
     do m = 1, l
        !
        ! Y_lm, m > 0
        !
        lm = lm + 1
        ylm (:, lm) = c * sqrt(DBLE(fact(l-m))/DBLE(fact(l+m))) * &
             sqrt(2.d0) * P(:,l,m) * cos (m*phi(:))
        !
        ! Y_lm, m < 0
        !
        lm = lm + 1
        ylm (:, lm) = c * sqrt(DBLE(fact(l-m))/DBLE(fact(l+m))) * &
             sqrt(2.d0) * P(:,l,m) * sin (m*phi(:))
     end do
  end do
  !
  deallocate(cost, phi, P)
  !
  return
end subroutine ylmr2

integer function fact(n)
  ! fact(n) = n!
  implicit none
  integer :: n, i

  fact = 1
  do i = n, 2, -1
     fact = i*fact
  end do
  return
end function fact

#endif

integer function semifact(n)
  ! semifact(n) = n!!
  implicit none
  integer :: n, i

  semifact = 1
  do i = n, 1, -2
     semifact = i*semifact
  end do
  return
end function semifact
