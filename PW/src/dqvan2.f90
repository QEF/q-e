!
! Copyright (C) 2001-2018 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dqvan2 (ih, jh, np, ipol, ngy, g, qmod, ylmk0, dylmk0, dqg )
  !-----------------------------------------------------------------------
  !
  !    This routine computes the derivatives of the fourier transform of
  !    the Q function needed in stress assuming that the radial fourier
  !    transform is already computed and stored in table qrad.
  !
  !    The formula implemented here is
  !
  !     dq(g,i,j) = sum_lm (-i)^l ap(lm,i,j) *
  !                ( yr_lm(g^) dqrad(g,l,i,j) + dyr_lm(g^) qrad(g,l,i,j))
  !
  USE kinds, ONLY: DP
  USE us, ONLY: dq, qrad
  USE uspp_param, ONLY: lmaxq, nbetam
  USE uspp, ONLY: nlx, lpl, lpx, ap, indv, nhtol, nhtolm
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ngy, ih, jh, np, ipol
  ! ngy: the number of G vectors to compute
  ! ih : first index of Q(ih,jh)
  ! jh : second index of Q(ih,jh)
  ! np : pseudopotential index
  ! ipol:the polarization of the derivative
  REAL(DP), INTENT(IN) ::  g(3,ngy), qmod (ngy), ylmk0 (ngy, lmaxq*lmaxq), &
       dylmk0 (ngy, lmaxq * lmaxq)
  ! g: G vectors
  ! qmod: moduli of q+G vectors
  ! ylmk0: spherical harmonics
  ! dylmk0: derivetives of spherical harmonics
  COMPLEX(DP), INTENT(OUT) :: dqg (ngy)
  ! the fourier transform of interest
  !
  COMPLEX(DP) :: sig
  ! (-i)^L
  INTEGER :: nb, mb, ijv, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3
  ! the atomic index corresponding to ih
  ! the atomic index corresponding to jh
  ! combined index (nb,mb)
  ! the lm corresponding to ih
  ! the lm corresponding to jh
  ! counter on g vectors
  ! the actual LM
  ! the angular momentum L
  ! the possible LM's compatible with ih,j
  ! counters for interpolation table

  REAL(DP) :: sixth, dqi, qm, px, ux, vx, wx, uvx, pwx, work, work1, qm1
  ! 1 divided by six
  ! 1 divided dq
  ! qmod/dq
  ! measures for interpolation table
  ! auxiliary variables for intepolation
  ! auxiliary variable
  ! auxiliary variable
  !
  !     compute the indices which correspond to ih,jh
  !
  sixth = 1.d0 / 6.d0
  dqi = 1 / dq
  nb = indv (ih, np)
  mb = indv (jh, np)
  IF (nb >= mb) THEN
     ijv = nb * (nb - 1) / 2 + mb
  ELSE
     ijv = mb * (mb - 1) / 2 + nb
  ENDIF
  ivl = nhtolm (ih, np)
  jvl = nhtolm (jh, np)

  IF (nb > nbetam .or. mb > nbetam) &
       CALL errore (' dqvan2 ', ' wrong dimensions (1)', max(nb,mb))
  IF (ivl > nlx .or. jvl > nlx) &
       CALL errore (' dqvan2 ', ' wrong dimensions (2)', max(ivl,jvl))

  dqg(:) = (0.d0,0.d0)
  !
  !    and make the sum over the non zero LM
  !
  DO lm = 1, lpx (ivl, jvl)
     lp = lpl (ivl, jvl, lm)
     !
     !     extraction of angular momentum l from lp:
     !
     IF (lp==1) THEN
        l = 1
     ELSEIF ( (lp>=2) .and. (lp<=4) ) THEN
        l = 2
     ELSEIF ( (lp>=5) .and. (lp<=9) ) THEN
        l = 3
     ELSEIF ( (lp>=10) .and. (lp<=16) ) THEN
        l = 4
     ELSEIF ( (lp>=17) .and. (lp<=25) ) THEN
        l = 5
     ELSEIF ( (lp>=26) .and. (lp<=36) ) THEN
        l = 6
     ELSEIF ( (lp>=37) .and. (lp<=49) ) THEN
        l = 7
     ELSE
        CALL errore (' dqvan2 ', ' lp.gt.49 ', lp)
     ENDIF

     sig = (0.d0, -1.d0) ** (l - 1)
     sig = sig * ap (lp, ivl, jvl)
     !
     qm1 = -1.0_dp !  any number smaller than qmod(1)
     !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(qm,px,ux,vx,wx,i0,i1,i2,i3,uvx,pwx,work,work1)
     DO ig = 1, ngy
        !
        ! calculate quantites depending on the module of G only when needed
        !
#if !defined(_OPENMP)
        IF ( abs( qmod(ig) - qm1 ) > 1.0D-6 ) THEN
#endif
           qm = qmod (ig) * dqi
           px = qm - int (qm)
           ux = 1.d0 - px
           vx = 2.d0 - px
           wx = 3.d0 - px
           i0 = qm + 1
           i1 = qm + 2
           i2 = qm + 3
           i3 = qm + 4
           uvx = ux * vx * sixth

           pwx = px * wx * 0.5d0

           work = qrad (i0, ijv, l, np) * uvx * wx + &
                  qrad (i1, ijv, l, np) * pwx * vx - &
                  qrad (i2, ijv, l, np) * pwx * ux + &
                  qrad (i3, ijv, l, np) * px * uvx
           work1 = - qrad(i0, ijv, l, np) * (ux*vx + vx*wx + ux*wx) * sixth &
                   + qrad(i1, ijv, l, np) * (wx*vx - px*wx - px*vx) * 0.5d0 &
                   - qrad(i2, ijv, l, np) * (wx*ux - px*wx - px*ux) * 0.5d0 &
                   + qrad(i3, ijv, l, np) * (ux*vx - px*ux - px*vx) * sixth

           work1 = work1 * dqi

#if !defined(_OPENMP)
           qm1 = qmod(ig)
        ENDIF
#endif

        dqg (ig) = dqg (ig) + sig * dylmk0 (ig, lp) * work
        IF (qmod (ig) > 1.d-9) dqg (ig) = dqg (ig) + &
            sig * ylmk0 (ig, lp) * work1 * g (ipol, ig) / qmod (ig)
     ENDDO
!$OMP END PARALLEL DO

  ENDDO
  RETURN

END SUBROUTINE dqvan2

