!
! Copyright (C) 2001-2018 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dqvan2( ih, jh, np, ipol, ngy, g, tpiba, qmod, ylmk0, dylmk0,&
                 dqg )
  !-----------------------------------------------------------------------
  !! This routine computes the derivatives of the Fourier transform of
  !! the Q function needed in stress assuming that the radial fourier
  !! transform is already computed and stored in table qrad.
  !! The formula implemented here is:

  !!     dq(g,i,j) = sum_lm (-i)^l ap(lm,i,j) *
  !!                ( yr_lm(g^) dqrad(g,l,i,j) + dyr_lm(g^) qrad(g,l,i,j))
  !
  USE kinds,       ONLY: DP
  USE us,          ONLY: dq, qrad
  USE uspp_param,  ONLY: lmaxq, nbetam
  USE uspp,        ONLY: nlx, lpl, lpx, ap, indv, nhtol, nhtolm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngy
  !! the number of G vectors to compute
  INTEGER, INTENT(IN) :: ih
  !! first index of Q(ih,jh)
  INTEGER, INTENT(IN) :: jh
  !! second index of Q(ih,jh)
  INTEGER, INTENT(IN) :: np
  !! pseudopotential index
  INTEGER, INTENT(IN) :: ipol
  !! the polarization of the derivative
  REAL(DP), INTENT(IN) ::  g(3,ngy)
  !! G vectors
  REAL(DP), INTENT(IN) ::  tpiba
  !! 2pi/a factor, multiplies G vectors
  REAL(DP), INTENT(IN) ::  qmod(ngy)
  !! moduli of q+G vectors
  REAL(DP), INTENT(IN) ::  ylmk0(ngy,lmaxq*lmaxq)
  !! spherical harmonics
  REAL(DP), INTENT(IN) ::  dylmk0(ngy,lmaxq*lmaxq)
  !! derivetives of spherical harmonics
  COMPLEX(DP), INTENT(OUT) :: dqg(ngy)
  !! the fourier transform of interest
  !
  ! ... local variables
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
  !
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
  nb = indv(ih, np)
  mb = indv(jh, np)
  IF (nb >= mb) THEN
     ijv = nb * (nb - 1) / 2 + mb
  ELSE
     ijv = mb * (mb - 1) / 2 + nb
  ENDIF
  ivl = nhtolm(ih, np)
  jvl = nhtolm(jh, np)
  !
  IF (nb > nbetam .OR. mb > nbetam) &
       CALL errore (' dqvan2 ', ' wrong dimensions (1)', MAX(nb,mb))
  IF (ivl > nlx .OR. jvl > nlx) &
       CALL errore (' dqvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  !
  dqg(:) = (0.d0,0.d0)
  !
  !    and make the sum over the non zero LM
  !
  DO lm = 1, lpx(ivl, jvl)
     lp = lpl(ivl, jvl, lm)
     !
     !     extraction of angular momentum l from lp:
     !
     IF (lp==1) THEN
        l = 1
     ELSEIF ( (lp>=2) .AND. (lp<=4) ) THEN
        l = 2
     ELSEIF ( (lp>=5) .AND. (lp<=9) ) THEN
        l = 3
     ELSEIF ( (lp>=10) .AND. (lp<=16) ) THEN
        l = 4
     ELSEIF ( (lp>=17) .AND. (lp<=25) ) THEN
        l = 5
     ELSEIF ( (lp>=26) .AND. (lp<=36) ) THEN
        l = 6
     ELSEIF ( (lp>=37) .AND. (lp<=49) ) THEN
        l = 7
     ELSE
        CALL errore (' dqvan2 ', ' lp.gt.49 ', lp)
     ENDIF
     !
     sig = (0.d0, -1.d0)**(l - 1)
     sig = sig * ap(lp, ivl, jvl)
     !
     qm1 = -1.0_dp !  any number smaller than qmod(1)
     !
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(qm,px,ux,vx,wx,i0,i1,i2,i3,uvx,pwx,work,work1)
     DO ig = 1, ngy
        !
        ! calculate quantites depending on the module of G only when needed
        !
#if !defined(_OPENMP)
        IF ( ABS( qmod(ig) - qm1 ) > 1.0D-6 ) THEN
#endif
           qm = qmod (ig) * dqi
           px = qm - INT(qm)
           ux = 1.d0 - px
           vx = 2.d0 - px
           wx = 3.d0 - px
           i0 = qm + 1
           i1 = qm + 2
           i2 = qm + 3
           i3 = qm + 4
           uvx = ux * vx * sixth

           pwx = px * wx * 0.5d0

           work = qrad(i0, ijv, l, np) * uvx * wx + &
                  qrad(i1, ijv, l, np) * pwx * vx - &
                  qrad(i2, ijv, l, np) * pwx * ux + &
                  qrad(i3, ijv, l, np) * px * uvx
           work1 = - qrad(i0, ijv, l, np) * (ux*vx + vx*wx + ux*wx) * sixth &
                   + qrad(i1, ijv, l, np) * (wx*vx - px*wx - px*vx) * 0.5d0 &
                   - qrad(i2, ijv, l, np) * (wx*ux - px*wx - px*ux) * 0.5d0 &
                   + qrad(i3, ijv, l, np) * (ux*vx - px*ux - px*vx) * sixth

           work1 = work1 * dqi

#if !defined(_OPENMP)
           qm1 = qmod(ig)
        ENDIF
#endif

        dqg(ig) = dqg(ig) + sig * dylmk0(ig, lp) * work / tpiba
        IF (qmod(ig) > 1.d-9) dqg(ig) = dqg(ig) + &
            sig * ylmk0(ig, lp) * work1 * tpiba * g(ipol, ig) / qmod(ig)
     ENDDO
!$OMP END PARALLEL DO
  !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE dqvan2

