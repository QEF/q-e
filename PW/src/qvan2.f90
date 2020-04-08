!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE qvan2( ngy, ih, jh, np, qmod, qg, ylmk0 )
  !-----------------------------------------------------------------------
  !! This routine computes the Fourier transform of the Q functions.
  !
  !! The interpolation table for the radial Fourier transform is stored 
  !! in qrad.
  !
  !! The formula implemented here is:
  !! \[ q(g,i,j) = \sum_\text{lm} (-i)^l \text{ap}(\text{lm},i,j) 
  !! \text{yr}_\text{lm}(g) \text{qrad}(g,l,i,j) \]
  !
  USE kinds,       ONLY: DP
  USE us,          ONLY: dq, qrad
  USE uspp_param,  ONLY: lmaxq, nbetam
  USE uspp,        ONLY: nlx, lpl, lpx, ap, indv, nhtolm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngy
  !! number of G vectors to compute
  INTEGER, INTENT(IN) :: ih
  !! first index of Q
  INTEGER, INTENT(IN) :: jh
  !! second index of Q
  INTEGER, INTENT(IN) :: np
  !! index of pseudopotentials
  REAL(DP), INTENT(IN) :: ylmk0(ngy,lmaxq*lmaxq)
  !! spherical harmonics
  REAL(DP), INTENT(IN) :: qmod(ngy)
  !! moduli of the q+g vectors 
  REAL(DP), INTENT(OUT) :: qg(2,ngy)
  !! the Fourier transform of interest
  !
  ! ... local variables
  !
  REAL(DP) :: sig
  ! the nonzero real or imaginary part of (-i)^L
  !
  REAL(DP), PARAMETER :: sixth = 1.0_DP / 6.0_DP
  !
  INTEGER :: nb, mb, ijv, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3, ind
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
  REAL(DP) :: dqi, qm, px, ux, vx, wx, uvx, pwx, work, qm1
  ! 1 divided dq
  ! qmod/dq
  ! measures for interpolation table
  ! auxiliary variables for intepolation
  ! auxiliary variables
  !
  ! ... computes the indices which correspond to ih,jh
  !
  dqi = 1.0_DP / dq
  nb = indv(ih,np)
  mb = indv(jh,np)
  !
  IF (nb >= mb) THEN
     ijv = nb * (nb - 1) / 2 + mb
  ELSE
     ijv = mb * (mb - 1) / 2 + nb
  ENDIF
  !
  ivl = nhtolm(ih,np)
  jvl = nhtolm(jh,np)
  !
  IF (nb > nbetam .OR. mb > nbetam) &
       CALL errore( ' qvan2 ', ' wrong dimensions (1)', MAX(nb,mb) )
  IF (ivl > nlx .OR. jvl > nlx) &
       CALL errore( ' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl) )
  !
  qg = 0.0_DP
  !
  ! ... and makes the sum over the non zero LM
  !
  DO lm = 1, lpx(ivl,jvl)
     lp = lpl(ivl,jvl,lm)
     IF ( lp < 1 .OR. lp > 49 ) CALL errore( 'qvan2', ' lp wrong ', MAX(lp,1) )
     !
     ! ... finds angular momentum l corresponding to combined index lp (l is 
     !     actually l+1 because this is the way qrad is stored, check init_us_1)
     !
     IF ( lp == 1 ) THEN
        l = 1
        sig = 1.0_DP
        ind = 1
     ELSEIF ( lp <=  4 ) THEN
        l = 2
        sig =-1.0_DP
        ind = 2
     ELSEIF ( lp <=  9 ) THEN
        l = 3
        sig =-1.0_DP
        ind = 1
     ELSEIF ( lp <= 16 ) THEN
        l = 4
        sig = 1.0_DP
        ind = 2
     ELSEIF ( lp <= 25 ) THEN
        l = 5
        sig = 1.0_DP
        ind = 1
     ELSEIF ( lp <= 36 ) THEN
        l = 6
        sig =-1.0_DP
        ind = 2
     ELSE
        l = 7
        sig =-1.0_DP
        ind = 1
     ENDIF
     !
     sig = sig * ap(lp, ivl, jvl)
     !
     qm1 = -1.0_DP !  any number smaller than qmod(1)
     !
!$omp parallel do default(shared), private(qm,px,ux,vx,wx,i0,i1,i2,i3,uvx,pwx,work)
     DO ig = 1, ngy
        !
        ! ... calculates quantites depending on the module of G only when needed
        !
#if ! defined _OPENMP
        IF ( ABS( qmod(ig) - qm1 ) > 1.0D-6 ) THEN
#endif
           !
           qm = qmod (ig) * dqi
           px = qm - INT(qm)
           ux = 1.0_DP - px
           vx = 2.0_DP - px
           wx = 3.0_DP - px
           i0 = INT(qm) + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           uvx = ux * vx * sixth
           pwx = px * wx * 0.5_DP
           work = qrad(i0,ijv,l,np) * uvx * wx + &
                  qrad(i1,ijv,l,np) * pwx * vx - &
                  qrad(i2,ijv,l,np) * pwx * ux + &
                  qrad(i3,ijv,l,np) * px * uvx
#if ! defined _OPENMP
           qm1 = qmod(ig)
        END IF
#endif
        qg(ind,ig) = qg(ind,ig) + sig * ylmk0(ig,lp) * work
        !
     ENDDO
!$omp end parallel do
  !
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE qvan2

