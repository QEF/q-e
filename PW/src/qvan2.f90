!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine qvan2 (ngy, ih, jh, np, qmod, qg, ylmk0)
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
  USE us, ONLY: dq, qrad
  USE uspp_param, ONLY: lmaxq, nbetam
  USE uspp, ONLY: nlx, lpl, lpx, ap, indv, nhtolm
  implicit none
  !
  ! Input variables
  !
  integer,intent(IN) :: ngy, ih, jh, np
  ! ngy   :   number of G vectors to compute
  ! ih, jh:   first and second index of Q
  ! np    :   index of pseudopotentials
  !
  real(DP),intent(IN) :: ylmk0 (ngy, lmaxq * lmaxq), qmod (ngy)
  ! ylmk0 :  spherical harmonics
  ! qmod  :  moduli of the q+g vectors
  !
  ! output: the fourier transform of interest
  !
  real(DP),intent(OUT) :: qg (2,ngy)
  !
  !     here the local variables
  !
  real (DP) :: sig
  ! the nonzero real or imaginary part of (-i)^L 
  real (DP), parameter :: sixth = 1.d0 / 6.d0
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
  real(DP) :: dqi, qm, px, ux, vx, wx, uvx, pwx, work, qm1
  ! 1 divided dq
  ! qmod/dq
  ! measures for interpolation table
  ! auxiliary variables for intepolation
  ! auxiliary variables
  !
  !     compute the indices which correspond to ih,jh
  !
  dqi = 1.0_DP / dq
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
  qg = 0.d0
  !
  !    and make the sum over the non zero LM
  !
  do lm = 1, lpx (ivl, jvl)
     lp = lpl (ivl, jvl, lm)
     if ( lp < 1 .or. lp > 49 ) call errore ('qvan2', ' lp wrong ', max(lp,1))
     !
     !     find angular momentum l corresponding to combined index lp
     !
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
     sig = sig * ap (lp, ivl, jvl)

     qm1 = -1.0_dp !  any number smaller than qmod(1)

!$omp parallel do default(shared), private(qm,px,ux,vx,wx,i0,i1,i2,i3,uvx,pwx,work)
     do ig = 1, ngy
        !
        ! calculate quantites depending on the module of G only when needed
        !
#if ! defined __OPENMP
        IF ( ABS( qmod(ig) - qm1 ) > 1.0D-6 ) THEN
#endif
           !
           qm = qmod (ig) * dqi
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
           work = qrad (i0, ijv, l, np) * uvx * wx + &
                  qrad (i1, ijv, l, np) * pwx * vx - &
                  qrad (i2, ijv, l, np) * pwx * ux + &
                  qrad (i3, ijv, l, np) * px * uvx
#if ! defined __OPENMP
           qm1 = qmod(ig)
        END IF
#endif
        qg (ind,ig) = qg (ind,ig) + sig * ylmk0 (ig, lp) * work
     enddo
!$omp end parallel do

  enddo

  return
end subroutine qvan2

