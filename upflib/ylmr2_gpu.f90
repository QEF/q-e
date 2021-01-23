!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

#define __USE_CUDA_KERNELS

module ylmr2_gpum
#if defined(__CUDA)
#if defined(__USE_CUDA_KERNELS)
use cudafor
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
  real(DP) ::  Q(0:4,0:4)  !Allocate Q for the maximum supported size 

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
#endif
end module ylmr2_gpum

subroutine ylmr2_gpu(lmax2, ng, g_d, gg_d, ylm_d)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical 
  !     Recipes but avoiding the calculation of factorials that generate 
  !     overflow for lmax > 11
  !
#if defined(__CUDA)
  USE cudafor
#if defined(__USE_CUDA_KERNELS)
  USE ylmr2_gpum, ONLY : ylmr2_gpu_kernel
#endif
#endif  
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g_d (3, ng), gg_d (ng)
  real(DP), intent(out) :: ylm_d (ng,lmax2)
#if defined(__CUDA)
  attributes(device):: g_d,gg_d,ylm_d
#endif
  !
#if defined(__USE_CUDA_KERNELS) && defined(__CUDA)
  !
  ! CUDA Fortran Kernel implementation. Optimizes the use of Q_d (see below)
  !
  integer:: lmax
  type(dim3):: grid,tBlock

  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !

  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, 25
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  call errore (' ylmr', 'l > 25 or wrong number of Ylm required',lmax2)

10 continue

  tBlock = dim3(256,1,1)
  grid = dim3(ceiling(real(ng)/tBlock%x),1,1)
  call ylmr2_gpu_kernel<<<grid,tBlock>>>(lmax, lmax2, ng, g_d, gg_d, ylm_d)
#else
  !
  ! Cuf kernel implementation.
  !
  real(DP), allocatable ::  Q_d(:,:,:)
#if defined(__CUDA)
  attributes(device):: Q_d
#endif
  real(DP) :: cost , sent, phi
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm
  real(DP), parameter :: eps = 1.0d-9
  !
  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, 25
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  call errore (' ylmr', 'l > 25 or wrong number of Ylm required',lmax2)
10 continue
  !
  if (lmax == 0) then
     !$cuf kernel do
     DO ig=1,ng
        ylm_d(ig,1) =  sqrt (1.d0 / fpi)
     END DO
     return
  end if
  ALLOCATE(Q_d(ng,0:lmax,0:lmax))
  !
  !  theta and phi are polar angles, cost = cos(theta)
  !
  ! present(g_d, gg_d, ylm_d) create(Q_d(ng,0:lmax,0:lmax))
  !$cuf kernel do
   do ig=1,ng
     gmod = sqrt (gg_d (ig) )
     if (gmod < eps) then
        cost = 0.d0
     else
        cost = g_d(3,ig)/gmod
     endif
     !
     !  beware the arc tan, it is defined modulo pi
     !
     if (g_d(1,ig) > eps) then
        phi  = atan( g_d(2,ig)/g_d(1,ig) )
     else if (g_d(1,ig) < -eps) then
        phi  = atan( g_d(2,ig)/g_d(1,ig) ) + pi
     else
        phi  = sign( pi/2.d0,g_d(2,ig) )
     end if
     sent = sqrt(max(0d0,1.d0-cost*cost))
     !
     !  Q_d(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
     !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
     !
     Q_d (ig,0,0) = 1.d0
     Q_d (ig,1,0) = cost
     Q_d (ig,1,1) =-sent/sqrt(2.d0)
     c = sqrt (3.d0 / fpi)
     ylm_d(ig, 1) = sqrt (1.d0 / fpi)* Q_d(ig,0,0)
     ylm_d(ig, 2) = c* Q_d(ig,1,0)
     ylm_d(ig, 3) = c*sqrt (2.d0)* Q_d(ig,1,1) * cos (phi)
     ylm_d(ig, 4) = c*sqrt (2.d0)* Q_d(ig,1,1) * sin (phi)
     lm = 4
     do l = 2, lmax
        c = sqrt (DBLE(2*l+1) / fpi)
        !if ( l == 0 ) then
        !      Q_d (0,0) = 1.d0
        !else if ( l == 1 ) then
        !      Q_d (1,0) = cost
        !      Q_d (1,1) =-sent/sqrt(2.d0)
        !else
           !
           !  recursion on l for Q_d(:,l,m)
           !
           do m = 0, l - 2
              Q_d(ig,l,m) = cost*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q_d(ig,l-1,m) &
                          - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q_d(ig,l-2,m)
           end do
           Q_d(ig,l,l-1) = cost * sqrt(DBLE(2*l-1)) * Q_d(ig,l-1,l-1)
           Q_d(ig,l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent*Q_d(ig,l-1,l-1)
        !end if
        !
        ! Y_lm, m = 0
        !
        lm = lm + 1
        ylm_d(ig, lm) = c * Q_d(ig,l,0)
        !
        do m = 1, l
           !
           ! Y_lm, m > 0
           !
           !lm = lm + 1
           ! ylm_d(ig, lm) = c * sqrt(2.d0) * Q_d(ig,l,m) * cos (m*phi)
           ylm_d(ig, lm+2*m-1) = c * sqrt(2.d0) * Q_d(ig,l,m) * cos (m*phi)
           !
           ! Y_lm, m < 0
           !
           !lm = lm + 1
           !ylm_d(ig, lm) = c * sqrt(2.d0) * Q_d(ig,l,m) * sin (m*phi)
           ylm_d(ig, lm+2*m) = c * sqrt(2.d0) * Q_d(ig,l,m) * sin (m*phi)
        end do
        lm=lm+2*l
        !if (ig==1) print *,lm,l
     end do
  end do
  !
  deallocate(Q_d)
  !
  return
#endif
end subroutine ylmr2_gpu
