!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module ylmr2_gpum
#if defined(__CUDA)

! CUDA Kernel version

use cudafor
  INTEGER, PARAMETER :: maxl = 20
contains
attributes(global) subroutine ylmr2_gpu_kernel (lmax,lmax2, ng, g, gg, ylm)
  implicit none
  !
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  real(DP), intent(out) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  real(DP) ::  Q(0:maxl,0:maxl)  !Allocate Q for the maximum supported size

  real(DP) :: cost , sent, phi
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm

  attributes(value)::lmax,lmax2,ng
  attributes(device):: g,gg,Q,ylm

  ig= threadIdx%x+BlockDim%x*(BlockIdx%x-1)

  if (ig <= ng) then
    !
    if (lmax == 0) then
      ylm(ig,1) =  sqrt (1.d0 / fpi)
      return
    end if
    !
    !  theta and phi are polar angles, cost = cos(theta)
    !
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
    Q(0,0) = 1.d0
    Q(1,0) = cost
    Q(1,1) =-sent/sqrt(2.d0)
    c = sqrt (3.d0 / fpi)
    ylm(ig, 1) = sqrt (1.d0 / fpi)* Q(0,0)
    ylm(ig, 2) = c* Q(1,0)
    ylm(ig, 3) = c*sqrt (2.d0)* Q(1,1) * cos (phi)
    ylm(ig, 4) = c*sqrt (2.d0)* Q(1,1) * sin (phi)
    lm = 4
    do l = 2, lmax
      c = sqrt (DBLE(2*l+1) / fpi)
      !
      !  recursion on l for Q(:,l,m)
      !
      do m = 0, l - 2
         Q(l,m) = cost*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(l-1,m) &
                     - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(l-2,m)
      end do
      Q(l,l-1) = cost * sqrt(DBLE(2*l-1)) * Q(l-1,l-1)
      Q(l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent*Q(l-1,l-1)
      !
      !
      ! Y_lm, m = 0
      !
      lm = lm + 1
      ylm(ig, lm) = c * Q(l,0)
      !
      do m = 1, l
       !
       ! Y_lm, m > 0
       !
         ylm(ig, lm+2*m-1) = c * sqrt(2.d0) * Q(l,m) * cos (m*phi)
       !
       ! Y_lm, m < 0
       !
         ylm(ig, lm+2*m  ) = c * sqrt(2.d0) * Q(l,m) * sin (m*phi)
      end do
      lm=lm+2*l
    end do
  end if
  return
  end subroutine ylmr2_gpu_kernel
#endif
end module ylmr2_gpum

subroutine ylmr2_gpu(lmax2, ng, g, gg, ylm)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax, GPU version
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical 
  !     Recipes but avoiding the calculation of factorials that generate 
  !     overflow for lmax > 11
  !     Last modified Jan 2024, by PG
  !
#if defined(__CUDA)
  USE cudafor
  USE ylmr2_gpum, ONLY : ylmr2_gpu_kernel, maxl
#endif
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  real(DP), intent(out) :: ylm (ng,lmax2)
#if defined(__CUDA)
  attributes(device):: g,gg,ylm
  !
  ! CUDA Fortran Kernel implementation. Optimizes the use of Q (see below)
  !
  integer:: lmax
  type(dim3):: grid,tBlock

  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !
  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, maxl
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  call upf_error (' ylmr', 'l too large or wrong number of Ylm required',lmax)

10 continue

  tBlock = dim3(256,1,1)
  grid = dim3(ceiling(real(ng)/tBlock%x),1,1)
  call ylmr2_gpu_kernel<<<grid,tBlock>>>(lmax, lmax2, ng, g, gg, ylm)
  !
#else
  call upf_error('ylmr2_gpu','you should not be here, go away!',1)
#endif
  return
end subroutine ylmr2_gpu
