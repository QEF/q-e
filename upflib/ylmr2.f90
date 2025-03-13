!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
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
  !     Last modified Jan. 2024, by PG: calls CUF version if __CUDA
  !
  USE upf_kinds, ONLY : DP
  USE upf_const, ONLY : pi, fpi
  !
  IMPLICIT NONE
  !
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  ! Incorrect results will ensue if gg != g(1)^2 + g(2)^2 +g(3)^2 
  !
  real(DP), intent(out) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  integer, parameter :: maxl = 20
  real(DP) :: cost , sent, phi 
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm, lm1, lm2
  !
  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, maxl
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  call upf_error (' ylmr2', 'l too large, or wrong number of Ylm required',lmax)
10 continue
  !
#if defined(__CUDA)
  !$acc data present_or_copyout(ylm) present_or_copyin(g, gg)
  !$acc host_data use_device(g,gg,ylm)
  call ylmr2_gpu(lmax2, ng, g, gg, ylm)
  !$acc end host_data
  !$acc end data
#else
  if (lmax == 0) then
     ylm(:,1) =  sqrt (1.d0 / fpi)
     return
  end if
  !
!$omp parallel do  default(shared), private(ig,gmod,lm,lm1,lm2,cost,sent,phi,l,c,m)
   do ig=1,ng
     gmod = sqrt (gg (ig) )
     if (gmod < eps) then
        cost = 0.d0
     else
        cost = g(3,ig)/gmod
     endif
     sent = sqrt(max(0.0_dp,1.0_dp-cost*cost))
     !
     !  cost = cos(theta), sent = sin(theta), with theta = polar angle
     !
     !  The Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m),
     !  where P(:,l,m) are the Legendre Polynomials (0 <= m <= l),
     !  and are currently stored into Ylm(:,lm) where lm = l**2+1+2*m
     !  (one might also store them into an array Q(l,m) for each ig)
     !
     ylm (ig,1) = 1.d0
     ylm (ig,2) = cost
     ylm (ig,4) =-sent/sqrt(2.d0)
     do l = 2, lmax
        !
        !  recursion on l for Q(:,l,m)
        !
        do m = 0, l - 2
           lm = (l  )**2 + 1 + 2*m
           lm1= (l-1)**2 + 1 + 2*m
           lm2= (l-2)**2 + 1 + 2*m
           ylm(ig,lm) = cost*(2*l-1)/sqrt(DBLE(l*l-m*m))*ylm(ig,lm1) &
                - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*ylm(ig,lm2)
        end do
        lm = (l  )**2 + 1 + 2*l
        lm1= (l  )**2 + 1 + 2*(l-1)
        lm2= (l-1)**2 + 1 + 2*(l-1)
        ylm(ig,lm1) = cost * sqrt(DBLE(2*l-1)) * ylm(ig,lm2)
        ylm(ig,lm ) = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent*ylm(ig,lm2) 
        !
     end do
     !
     ! now add cos(phi), sin(phi), and other factors to get the true Ylm
     ! beware the arc tan, it is defined modulo pi
     !
     if (g(1,ig) > eps) then
        phi  = atan( g(2,ig)/g(1,ig) )
     else if (g(1,ig) < -eps) then
        phi  = atan( g(2,ig)/g(1,ig) ) + pi
     else
        phi  = sign( pi/2.d0,g(2,ig) )
     end if
     lm = 1
     ylm(ig,1) = ylm(ig,1) / sqrt(fpi)
     !
     do l = 1, lmax
        !
        c = sqrt (DBLE(2*l+1) / fpi)
        !
        ! Y_lm, m = 0
        !
        lm = lm + 1
        ylm(ig, lm) = c * ylm(ig,lm)
        !
        do m = 1, l
           !
           ! Y_lm, m > 0
           !
           lm = lm + 2
           ylm(ig, lm-1) = c * sqrt(2.d0) * ylm(ig,lm) * cos (m*phi)
           !
           ! Y_lm, m < 0
           !
           ylm(ig, lm  ) = c * sqrt(2.d0) * ylm(ig,lm) * sin (m*phi)
           !
        end do
     end do
  enddo
  !$omp end parallel do
#endif
  !
  return
end subroutine ylmr2
