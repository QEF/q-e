!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
  USE kinds, ONLY : DP
  USE constants, ONLY : pi, fpi
  !
  IMPLICIT NONE
  !
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !
  real(DP), intent(out) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  real(DP), allocatable ::  Q(:,:,:)
  real(DP) :: cost , sent, phi 
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm
  !
  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, 25
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  call errore (' ylmr', 'l > 25 or wrong number of Ylm required',lmax2)
10 continue

  !
  if (lmax == 0) then
     ylm(:,1) =  sqrt (1.d0 / fpi)
     return
  end if
  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  allocate( Q(ng,0:lmax,0:lmax) )
  !
  
!$omp parallel do  default(shared), private(ig,gmod,lm,cost,sent,phi,l,c,m)
   do ig=1,ng
     gmod = sqrt (gg (ig) )
     if (gmod < eps) then
        cost = 0.d0
     else
        cost = g(3,ig)/gmod
     endif
     !
     !  beware the arc tan, it is defined modulo pi
     !
     if (g(1,ig) > eps) then
        phi  = atan( g(2,ig)/g(1,ig) )
     else if (g(1,ig) < -eps) then
        phi  = atan( g(2,ig)/g(1,ig) ) + pi
     else
        phi  = sign( pi/2.d0,g(2,ig) )
     end if
     sent = sqrt(max(0d0,1.d0-cost*cost))
  !
  !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
  !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  !
   
  Q (ig,0,0) = 1.d0
  Q (ig,1,0) = cost
  Q (ig,1,1) =-sent/sqrt(2.d0)
  c = sqrt (3.d0 / fpi)
  ylm(ig, 1) = sqrt (1.d0 / fpi)* Q(ig,0,0)
  ylm(ig, 2) = c* Q(ig,1,0)
  ylm(ig, 3) = c*sqrt (2.d0)* Q(ig,1,1) * cos (phi)
  ylm(ig, 4) = c*sqrt (2.d0)* Q(ig,1,1) * sin (phi)
  lm = 4
  do l = 2, lmax
     c = sqrt (DBLE(2*l+1) / fpi)
     !if ( l == 0 ) then
     !      Q (0,0) = 1.d0
     !else if ( l == 1 ) then
     !      Q (1,0) = cost
     !      Q (1,1) =-sent/sqrt(2.d0)
     !else
        !
        !  recursion on l for Q(:,l,m)
        !
        do m = 0, l - 2
           Q(ig,l,m) = cost*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(ig,l-1,m) &
                       - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(ig,l-2,m)
        end do
           Q(ig,l,l-1) = cost * sqrt(DBLE(2*l-1)) * Q(ig,l-1,l-1)
           Q(ig,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent*Q(ig,l-1,l-1) 
     !end if
     !
     ! Y_lm, m = 0
     !
     lm = lm + 1
        ylm(ig, lm) = c * Q(ig,l,0)
     !
     do m = 1, l
        !
        ! Y_lm, m > 0
        !
        !lm = lm + 1
        ! ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * cos (m*phi)
         ylm(ig, lm+2*m-1) = c * sqrt(2.d0) * Q(ig,l,m) * cos (m*phi)
        !
        ! Y_lm, m < 0
        !
        !lm = lm + 1
        !ylm(ig, lm) = c * sqrt(2.d0) * Q(ig,l,m) * sin (m*phi)
        ylm(ig, lm+2*m) = c * sqrt(2.d0) * Q(ig,l,m) * sin (m*phi)
     end do
     lm=lm+2*l
     !if (ig==1) print *,lm,l
  end do
  enddo
!$omp end parallel do
  !
  deallocate( Q)
  !
  return
end subroutine ylmr2
