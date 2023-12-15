!
! Copyright (C) 2001-2018 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dqvan2( ih, jh, np, ipol, ngy, g, tpiba, qmod, ylmk0, dylmk0, dqg )
  !-----------------------------------------------------------------------
  !! This routine computes the derivatives of the Fourier transform of
  !! the Q function needed in stress assuming that the radial Fourier
  !! transform is already computed and stored in table \(\text{tab_qrad}\).
  !! The implemented formula:
  !
  !! \[ \text{dq}(g,i,j) = \sum_lm (-i)^l \text{ap}(lm,i,j) *
  !!                       ( \text{yr}_{lm}(g^) \text{dqrad}(g,l,i,j) + 
  !!                          \text{dyr}_{lm}(g') \text{qrad}(g,l,i,j)) \]
  !
  USE upf_kinds,   ONLY: DP
  USE qrad_mod,    ONLY: dq, tab_qrad
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
  COMPLEX(DP) :: sig, dqg_bgr
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
  REAL(DP) :: sixth, dqi, qm, px, ux, vx, wx, uvx, pwx, work, work1
  ! 1 divided by six
  ! 1 divided dq
  ! qmod/dq
  ! measures for interpolation table
  ! auxiliary variables for intepolation
  ! auxiliary variable
  ! auxiliary variable
  !
  ! ... compute the indices which correspond to ih,jh
  !
  !$acc data present_or_copyin(g,qmod,ylmk0,dylmk0) present_or_copyout(dqg) present(tab_qrad)
  !
  sixth = 1.d0 / 6.d0
  nb = indv(ih, np)
  mb = indv(jh, np)
  IF (nb >= mb) THEN
     ijv = nb * (nb - 1) / 2 + mb
  ELSE
     ijv = mb * (mb - 1) / 2 + nb
  ENDIF
  !
  ivl = nhtolm(ih, np)
  jvl = nhtolm(jh, np)
  !
  IF (nb > nbetam .OR. mb > nbetam) &
       CALL upf_error (' dqvan2 ', ' wrong dimensions (1)', MAX(nb,mb))
  IF (ivl > nlx .OR. jvl > nlx) &
       CALL upf_error (' dqvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  !
  !$acc kernels
  dqg(:) = (0.d0,0.d0)
  !$acc end kernels
  !
  ! ... and make the sum over the non zero LM
  !
  dqi = 1 / dq
  DO lm = 1, lpx(ivl,jvl)
     lp = lpl(ivl,jvl,lm)
     !
     ! ... extraction of angular momentum l from lp:
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
        CALL upf_error (' dqvan2 ', ' lp.gt.49 ', lp)
     ENDIF
     !
     sig = (0.d0, -1.d0)**(l - 1)
     sig = sig * ap(lp,ivl,jvl)
     !
!$acc parallel loop
     DO ig = 1, ngy
        !
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
        !
        pwx = px * wx * 0.5d0
        !
        work = tab_qrad(i0, ijv, l, np) * uvx * wx + &
               tab_qrad(i1, ijv, l, np) * pwx * vx - &
               tab_qrad(i2, ijv, l, np) * pwx * ux + &
               tab_qrad(i3, ijv, l, np) * px * uvx
        work1 = (- tab_qrad(i0, ijv, l, np) * (ux*vx + vx*wx + ux*wx) * sixth &
                 + tab_qrad(i1, ijv, l, np) * (wx*vx - px*wx - px*vx) * 0.5d0 &
                 - tab_qrad(i2, ijv, l, np) * (wx*ux - px*wx - px*ux) * 0.5d0 &
                 + tab_qrad(i3, ijv, l, np) * (ux*vx - px*ux - px*vx) * sixth) * dqi
        !
        IF (qmod(ig) > 1.d-9) THEN
          dqg_bgr = sig * ylmk0(ig,lp) * work1 * tpiba * g(ipol,ig) / qmod(ig)
        ELSE
          dqg_bgr = (0.d0,0.d0)
        ENDIF
        !
        dqg(ig) = dqg(ig) + sig * dylmk0(ig,lp) * work / tpiba + dqg_bgr
        !
     ENDDO
     !
  ENDDO
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE dqvan2

