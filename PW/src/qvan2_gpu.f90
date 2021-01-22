!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
#if defined(__CUDA)
module qvan2_gpum
contains
attributes(global) subroutine qvan2_kernel(ngy, ih, jh, np, qmod_d, qg_d, ylmk0_d, lmaxq, nbetam, nlx, dq)
   !-----------------------------------------------------------------------
   !
   !    This routine computes the fourier transform of the Q functions
   !    The interpolation table for the radial fourier trasform is stored 
   !    in qrad.
   !
   !    The formula implemented here is
   !
   !     q(g,i,j) = sum_lm (-i)^l ap(lm,i,j) yr_lm(g^) qrad(g,l,i,j)
   !
   !
   USE kinds, ONLY: DP
   USE us_gpum, ONLY: qrad_d
   USE uspp,      ONLY: lpl_d, lpx_d, ap_d
   USE uspp_gpum, ONLY : indv_d, nhtolm_d
   implicit none
   !
   ! Input variables
   !
   REAL(DP),intent(IN) :: dq
   integer, intent(IN) :: ngy, ih, jh, np, lmaxq, nbetam, nlx
   attributes(value):: ngy, ih, jh, np, lmaxq, nbetam, nlx, dq
   ! ngy   :   number of G vectors to compute
   ! ih, jh:   first and second index of Q
   ! np    :   index of pseudopotentials
   !
   real(DP),intent(IN) :: ylmk0_d (ngy, lmaxq * lmaxq), qmod_d (ngy)
   ! ylmk0 :  spherical harmonics
   ! qmod  :  moduli of the q+g vectors
   !
   ! output: the fourier transform of interest
   !
   real(DP),intent(OUT) :: qg_d (2,ngy)
   attributes(device):: ylmk0_d, qmod_d, qg_d
   !
   !     here the local variables
   !
   real (DP) :: sig
   ! the nonzero real or imaginary part of (-i)^L 
   real (DP), parameter :: sixth = 1.d0 / 6.d0
   !
   integer :: nb, mb, ijv, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3, ind
   real(DP) :: dqi, qm, px, ux, vx, wx, uvx, pwx, work, qm1
   real(DP) :: uvwx, pwvx, pwux, pxuvx
   !
   ig= threadIdx%x+BlockDim%x*(BlockIdx%x-1)
   !
   if (ig <= ngy) then
      !     compute the indices which correspond to ih,jh
      dqi = 1.0_DP / dq
      nb = indv_d (ih, np)
      mb = indv_d (jh, np)
      if (nb.ge.mb) then
         ijv = nb * (nb - 1) / 2 + mb
      else
         ijv = mb * (mb - 1) / 2 + nb
      endif
      ivl = nhtolm_d(ih, np)
      jvl = nhtolm_d(jh, np)
      
      qg_d(1,ig) = 0.d0
      qg_d(2,ig) = 0.d0
      
      qm = qmod_d (ig) * dqi
      px = qm - int (qm)
      ux = 1.d0 - px
      vx = 2.d0 - px
      wx = 3.d0 - px
      i0 = INT( qm ) + 1
      i1 = i0 + 1
      i2 = i0 + 2
      i3 = i0 + 3
      uvx = ux * vx * sixth
      pwx = px * wx * 0.5d0
      
      do lm = 1, lpx_d (ivl, jvl)
         lp = lpl_d (ivl, jvl, lm)
          if (lp == 1) then
             l = 1
             sig = 1.0d0
             ind = 1
          elseif ( lp <= 4) then
             l = 2
             sig =-1.0d0
             ind = 2
          elseif ( lp <= 9 ) then
             l = 3
             sig =-1.0d0
             ind = 1
          elseif ( lp <= 16 ) then
             l = 4
              sig = 1.0d0
             ind = 2
          elseif ( lp <= 25 ) then
             l = 5
             sig = 1.0d0
             ind = 1
          elseif ( lp <= 36 ) then
             l = 6
             sig =-1.0d0
             ind = 2
          else
             l = 7
             sig =-1.0d0
             ind = 1
          endif
          sig = sig * ap_d (lp, ivl, jvl)
               work = qrad_d (i0, ijv, l, np) * uvx * wx + &
                      qrad_d (i1, ijv, l, np) * pwx * vx - &
                      qrad_d (i2, ijv, l, np) * pwx * ux + &
                      qrad_d (i3, ijv, l, np) * px * uvx
          qg_d (ind,ig) = qg_d (ind,ig) + sig * ylmk0_d (ig, lp) * work
      end do
   end if
   !
   return
end subroutine qvan2_kernel
end module qvan2_gpum
#endif


subroutine qvan2_gpu (ngy, ih, jh, np, qmod_d, qg_d, ylmk0_d)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the fourier transform of the Q functions
  !    The interpolation table for the radial fourier transform is stored 
  !    in qrad.
  !
  !    The formula implemented here is
  !
  !     q(g,i,j) = sum_lm (-i)^l ap(lm,i,j) yr_lm(g^) qrad(g,l,i,j)
  !
  !
  USE kinds,       ONLY: DP
  USE uspp_param,  ONLY: lmaxq, nbetam
  USE uspp,        ONLY: nlx, nhtolm, indv
  USE us,          ONLY: dq
  !
  USE us_gpum,     ONLY : using_qrad_d
  USE uspp_gpum,   ONLY : using_indv_d, using_nhtolm_d, &
                          using_indv, using_nhtolm
#if defined(__CUDA)
  USE cudafor
  USE qvan2_gpum,  ONLY : qvan2_kernel
#endif
  !
  implicit none
  !
  ! Input variables
  !
  integer,intent(IN) :: ngy, ih, jh, np
  ! ngy   :   number of G vectors to compute
  ! ih, jh:   first and second index of Q
  ! np    :   index of pseudopotentials
  !
  real(DP),intent(IN) :: ylmk0_d (ngy, lmaxq * lmaxq), qmod_d (ngy)
  ! ylmk0 :  spherical harmonics
  ! qmod  :  moduli of the q+g vectors
  !
  ! output: the fourier transform of interest
  !
  real(DP),intent(OUT) :: qg_d (2,ngy)
  !
  !     here the local variables
  !
  integer :: nb, mb, ijv, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3, ind
  ! nb,mb  : atomic index corresponding to ih,jh
  ! ijv    : combined index (nb,mb)
  ! ivl,jvl: combined LM index corresponding to ih,jh
  ! ig     : counter on g vectors
  ! lp     : combined LM index
  ! l-1    is the angular momentum L
  ! lm     : all possible LM's compatible with ih,jh
  ! i0-i3  : counters for interpolation table
  ! ind    : ind=1 if the results is real (l even), ind=2 if complex (l odd)
  !
  !
  !     compute the indices which correspond to ih,jh
  !
#if defined(__CUDA)
  type(dim3):: grid, tBlock
  attributes(device):: ylmk0_d, qmod_d, qg_d
#endif
  !
  CALL using_indv(0)
  CALL using_nhtolm(0)
  !
  nb = indv (ih, np)
  mb = indv (jh, np)
  if (nb.ge.mb) then
     ijv = nb * (nb - 1) / 2 + mb
  else
     ijv = mb * (mb - 1) / 2 + nb
  endif
  ivl = nhtolm(ih, np)
  jvl = nhtolm(jh, np)
  if (nb > nbetam .OR. mb > nbetam) &
       call errore (' qvan2 ', ' wrong dimensions (1)', MAX(nb,mb))
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  !
  ! Sync (if needed) global variables used in kernel
  CALL using_qrad_d(0)
  CALL using_indv_d(0)
  CALL using_nhtolm_d(0)
  !
#if defined(__CUDA)
  tBlock = dim3(256,1,1)
  grid = dim3(ceiling(real(ngy)/tBlock%x),1,1)
  call qvan2_kernel<<<grid,tBlock>>>(ngy, ih, jh, np, qmod_d, qg_d, ylmk0_d, lmaxq, nbetam, nlx, dq)
  !
#else
  ! possibly change this to call to CPU version...
  CALL errore('qvan2_gpu', 'Trying to use device subroutine but code was not compiled with device support!', 1)
#endif
  !
end subroutine qvan2_gpu

