!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Routines for linear-response calculation using optimized tetrahedron method
! Written by Mitsuaki Kawamura, U. Tokyo
!
!----------------------------------------------------------------------------
MODULE dfpt_tetra_mod
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: dfpt_tetra_beta, dfpt_tetra_ttheta, dfpt_tetra_delta, &
            dfpt_tetra_main, dfpt_tetra_linit
  !
  LOGICAL,SAVE :: dfpt_tetra_linit = .FALSE.
  !
  REAL(dp),ALLOCATABLE,SAVE :: &
  & dfpt_tetra_beta(:,:,:), & ! (nbnd, nbnd, nksq) beta(e_k, e_k+q)
  & dfpt_tetra_ttheta(:,:,:), & ! (nbnd, nbnd, nksq) t(ef-e_k)t(e_k - e_k+q)+ t(ef-e_k+q)t(e_k+q - e_k) 
  & dfpt_tetra_delta(:,:) ! (nbnd, nksq) delta(e_F - e_k)
  !
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_main()
  !----------------------------------------------------------------------------
  !
  ! This routine computes the integration weight as represented in Eq. (B 28)
  ! in PRB 64, 235118 (2001).
  !
  USE kinds,      ONLY : dp
  USE mp,         ONLY : mp_sum
  USE mp_images,  ONLY : me_image, nproc_image, intra_image_comm
  USE io_global,  ONLY : stdout
  USE wvfct,      ONLY : nbnd
  USE ktetra,     ONLY : ntetra
  USE klist,      ONLY : wk, nks, nkstot
  USE qpoint,     ONLY : nksq, xq
  USE control_lr, ONLY : lgamma, alpha_pv, nbnd_occ
  USE ener,       ONLY : ef
  USE cell_base,  ONLY : at
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE
  !
  INTEGER :: ntpp, rest, ik, ibnd
  INTEGER :: &
  & tfst, & ! the first tetrahedron in this Prosess Element
  & tlst, & ! the last tetrahedron in this PE
  & iq      ! = 0 for q = Gamma, = 1 for other
  !
  REAL(dp) :: qvec(3)
  !
  REAL(dp),ALLOCATABLE :: &
  & et_col(:,:), & ! (nbnd, nkstot) Kohn-Sham energy
  & delta(:,:), & ! (nbnd,nkstot) df / de
  & beta(:,:,:)  ! (nbnd,nbnd,nkstot) 
  !
  ALLOCATE(delta(nbnd,nkstot), beta(nbnd,nbnd,nkstot), et_col(nbnd, nkstot))
  !
  ! IF .NOT. setup_pw, initialize tetrahedra
  !
  CALL dfpt_tetra_setup(et_col)
  !
  ! Work Sharing
  !
  ntpp = ntetra / nproc_image
  rest = MOD(ntetra, nproc_image)
  IF(me_image < rest) THEN
     tfst = (ntpp + 1) *  me_image + 1
     tlst = (ntpp + 1) * (me_image + 1)
  ELSE
     tfst = ntpp *  me_image + 1  + rest
     tlst = ntpp * (me_image + 1) + rest
  END IF
  !
  alpha_pv = ef - MINVAL(et_col(1:nbnd,1:nkstot))
  !
  IF(lgamma) THEN
     iq = 0
  ELSE
     iq = 1
  END IF
  !
  ! Calculate derivative of occupation
  !
  CALL dfpt_tetra_calc_delta(tfst,tlst,et_col,delta)
  !
  ! Integration weight(1st & 2nd term)
  !
  beta(1:nbnd,1:nbnd,1:nkstot) = 0.0_dp
  CALL dfpt_tetra_calc_beta1(iq,tfst,tlst,et_col,beta)
  WRITE(stdout,*) ""
  !WRITE(stdout,*) "[dfpt_tetra]  beta1 : ", SUM(beta)
  CALL dfpt_tetra_calc_beta2(iq,tfst,tlst,et_col,beta)
  !WRITE(stdout,*) "[dfpt_tetra]  beta2 : ", SUM(beta)
  !
  CALL mp_sum(beta(1:nbnd,1:nbnd,1:nkstot), intra_image_comm)
  !
  delta(      1:nbnd,1:nkstot) = delta(       1:nbnd,1:nkstot) / REAL(ntetra, dp)
  beta(1:nbnd,1:nbnd,1:nkstot) = beta(1:nbnd, 1:nbnd,1:nkstot) / REAL(ntetra, dp)
  !
  IF(nspin == 1) THEN
     WRITE(stdout,'(a,e15.7)') "     [dfpt_tetra]  Dos(E_F)[/Ry] : ", SUM(delta) * 2.0_dp
  ELSE
     WRITE(stdout,'(a,e15.7)') "     [dfpt_tetra]  Dos(E_F)[/Ry] : ", SUM(delta)
  END IF
  !
  ! Map k-point 
  !
  IF(ALLOCATED(dfpt_tetra_ttheta)) DEALLOCATE(dfpt_tetra_ttheta) 
  IF(ALLOCATED(dfpt_tetra_delta )) DEALLOCATE(dfpt_tetra_delta )
  !
  ALLOCATE(dfpt_tetra_ttheta(nbnd,nbnd,nks), dfpt_tetra_delta(nbnd,nks))
  !
  CALL poolscatter( nbnd,      nkstot, delta, nks,dfpt_tetra_delta )
  CALL poolscatter( nbnd*nbnd, nkstot, beta,  nks,dfpt_tetra_ttheta)
  !
  ! Integration weight(3rd term)
  !
  beta(1:nbnd,1:nbnd,1:nkstot) = 0.0_dp
  CALL dfpt_tetra_calc_beta3(iq,tfst,tlst,et_col,beta)
  !WRITE(stdout,*) "[dfpt_tetra]  beta3 : ", SUM(beta)
  !
  beta(1:nbnd,1:nbnd,1:nkstot) = beta(1:nbnd, 1:nbnd,1:nkstot) / REAL(ntetra, dp)
  !
  CALL mp_sum(beta(1:nbnd,1:nbnd,1:nkstot), intra_image_comm)
  !
  ! Map k-point 
  !
  IF(ALLOCATED(dfpt_tetra_beta)) DEALLOCATE(dfpt_tetra_beta)
  ALLOCATE(dfpt_tetra_beta(nbnd,nbnd,nks))
  !
  CALL poolscatter( nbnd*nbnd, nkstot, beta, nks, dfpt_tetra_beta)
  CALL dfpt_tetra_average_beta(dfpt_tetra_ttheta)
  CALL dfpt_tetra_average_beta(dfpt_tetra_beta)
  !
  IF(nspin == 1) THEN
     dfpt_tetra_ttheta(1:nbnd,1:nbnd,1:nks) = dfpt_tetra_ttheta(1:nbnd,1:nbnd,1:nks) * 2.0_dp
     dfpt_tetra_beta(  1:nbnd,1:nbnd,1:nks) = dfpt_tetra_beta(  1:nbnd,1:nbnd,1:nks) * 2.0_dp
     dfpt_tetra_delta(        1:nbnd,1:nks) = dfpt_tetra_delta(        1:nbnd,1:nks) * 2.0_dp
  END IF
  !
  DO ik = 1, nks
     IF(ABS(wk(ik)) > 1e-10_dp) THEN
        dfpt_tetra_ttheta(1:nbnd,1:nbnd,ik) = dfpt_tetra_ttheta(1:nbnd,1:nbnd,ik) / wk(ik)
        dfpt_tetra_beta(  1:nbnd,1:nbnd,ik) = dfpt_tetra_beta(  1:nbnd,1:nbnd,ik) / wk(ik)
        dfpt_tetra_delta(        1:nbnd,ik) = dfpt_tetra_delta(        1:nbnd,ik) / wk(ik)
     END IF
     nbnd_occ(ik) = nbnd
  END DO
  !
  ! Drude term
  !
  qvec(1:3) = MATMUL(xq(1:3), at(1:3,1:3))
  qvec(1:3) = ABS(qvec(1:3) - REAL(NINT(qvec(1:3)), dp))
  !
  IF(lgamma .OR. MAXVAL(qvec(1:3)) < 1e-8_dp) THEN
     !
     WRITE(stdout,'(a)') "     [dfpt_tetra]  Add Drude term"
     !
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           dfpt_tetra_beta(ibnd,ibnd,ik) = dfpt_tetra_beta(ibnd, ibnd, ik) &
           &                   - 0.5_dp * dfpt_tetra_delta(      ibnd, ik)
        END DO
     END DO
  END IF
  !
  dfpt_tetra_beta(1:nbnd,1:nbnd,1:nks) = dfpt_tetra_ttheta(1:nbnd,1:nbnd,1:nks) &
  &                           + alpha_pv * dfpt_tetra_beta(1:nbnd,1:nbnd,1:nks)
  !
  DEALLOCATE(delta,beta,et_col)
  !
END SUBROUTINE dfpt_tetra_main
!
!----------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_setup(et_col)
  !--------------------------------------------------------------------------
  !
  ! This routine initialize tetrahedra, if it is not done in setup_nscf
  !
  USE kinds,      ONLY : dp
  USE lsda_mod,   ONLY : nspin
  USE wvfct,      ONLY : nbnd, et
  USE ktetra,     ONLY : tetra, opt_tetra_init
  USE klist,      ONLY : xk, nks, nkstot
  USE start_k,    ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base,  ONLY : s, t_rev, time_reversal
  USE cell_base,  ONLY : at, bg
  USE lr_symm_base,  ONLY : nsymq
  USE parameters,    ONLY : npk
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(OUT) :: et_col(nbnd, nkstot)
  !
  INTEGER :: nktot
  REAL(dp),ALLOCATABLE :: xk_col(:,:)
  !
  CALL poolcollect(nbnd, nks, et, nkstot, et_col)
  !
  IF (dfpt_tetra_linit) THEN
     !
     ALLOCATE(xk_col(3, nkstot))
     CALL poolcollect(3, nks, xk, nkstot, xk_col)
     !
     IF(nspin == 2) THEN
        nktot = nkstot / 2
     ELSE
        nktot = nkstot
     END IF
     !
     IF(ALLOCATED(tetra)) DEALLOCATE(tetra)
     CALL opt_tetra_init(nsymq, s, time_reversal, t_rev, at, bg, npk, &
                         k1, k2, k3, nk1, nk2, nk3, nktot, xk_col, 1)
     !
     DEALLOCATE(xk_col)
     !
  END IF
  !
END SUBROUTINE dfpt_tetra_setup
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_delta(tfst,tlst,et_col,delta)
  !--------------------------------------------------------------------------
  !
  ! This routine compute derivative of occupation
  !
  USE ener,  ONLY : ef
  USE kinds, ONLY : dp
  USE klist, ONLY : nkstot
  USE ktetra, ONLY : ntetra, tetra, wlsm
  USE wvfct, ONLY : nbnd
  USE lsda_mod,   ONLY : nspin
  USE mp,         ONLY : mp_sum
  USE mp_images,   ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: tfst, tlst
  REAL(dp),INTENT(IN) :: et_col(nbnd,nkstot)
  REAL(dp),INTENT(OUT) :: delta(nbnd,nkstot)
  !
  INTEGER :: nspin_lsda, ns, nt, nk, ii, ibnd, jbnd, kbnd, itetra(4), ik
  REAL(dp) :: ei(4,nbnd), e(4), wdos1(4), a(4,4), V, wdos0(4,nbnd), wdos2
  !
  delta(1:nbnd,1:nkstot) = 0.0_dp
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  DO ns = 1, nspin_lsda
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     END IF
     !
     DO nt = tfst, tlst
        !
        ei(1:4,1:nbnd) = 0.0_dp
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           !        
           DO ibnd = 1, nbnd
              ei(1:4, ibnd) = ei(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd,ik) - ef)
           END DO
           !
        END DO
        !
        wdos0(1:nbnd,1:4) = 0.0_dp
        !
        DO ibnd = 1, nbnd
           !
           itetra(1) = 0
           e(1:4) = ei(1:4,ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0.0_dp - e(1:4) ) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0.0_dp .AND. 0.0_dp < e(2)) THEN
              !
              V = a(2,1) * a(3,1) * a(4,1) / (0_dp - e(1) )
              wdos1(1) = a(1,2) + a(1,3) + a(1,4)
              wdos1(2:4) = a(2:4,1)
              wdos1(1:4) = wdos1(1:4) * V
              !
           ELSE IF(e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
              !
              V = a(2,3) * a(3,1) + a(3,2) * a(2,4)
              wdos1(1) = a(1,4) * V + a(1,3) * a(3,1) * a(2,3)
              wdos1(2) = a(2,3) * V + a(2,4) * a(2,4) * a(3,2)
              wdos1(3) = a(3,2) * V + a(3,1) * a(3,1) * a(2,3)
              wdos1(4) = a(4,1) * V + a(4,2) * a(2,4) * a(3,2)
              V  = 1.0_dp / ( e(4) - e(1) )
              wdos1(1:4) = wdos1(1:4) * V
              !
           ELSE IF(e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              V = a(1,4) * a(2,4) * a(3,4) / ( e(4) - 0.0_dp )
              wdos1(1:3)  = a(1:3,4)
              wdos1(4)  = a(4,1) + a(4,2) + a(4,3)
              wdos1(1:4) = wdos1(1:4) * V
              !
           ELSE
              !
              wdos1(1:4) = 0.0_dp
              !
           END IF
           !
           wdos0(itetra(1:4),ibnd) = wdos1(1:4)
           !
        END DO ! ib
        !
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           delta(1:nbnd,ik) = delta(1:nbnd,ik) &
           &            + MATMUL(wlsm(1:4,ii), wdos0(1:4,1:nbnd))
           !
        END DO
        !
     END DO ! it
     !
  END DO ! ik
  !
  ! Average weights of degenerated states
  !
  DO ik = 1, nkstot
     DO ibnd = 1, nbnd
        !
        wdos2 = delta(ibnd,ik)
        !
        DO jbnd = ibnd + 1, nbnd
           !
           IF(ABS(et_col(ibnd,ik) - et_col(jbnd,ik)) < 1e-6_dp) THEN
              wdos2 = wdos2 + delta(jbnd,ik)
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 delta(kbnd,ik) = wdos2 / real(jbnd - ibnd, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
  CALL mp_sum(delta(1:nbnd,1:nkstot), intra_image_comm)
  !
END SUBROUTINE dfpt_tetra_calc_delta
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_beta1(iq,tfst,tlst,et_col,beta)
  !--------------------------------------------------------------------------
  !
  ! This routine compute the first term of (B 28) in PRB 64, 235118 (2001).
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE ener,  ONLY : ef
  USE klist, ONLY : nkstot
  USE ktetra, ONLY : ntetra, tetra, wlsm
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: iq, tfst, tlst
  REAL(dp),INTENT(IN) :: et_col(nbnd, nkstot)
  REAL(dp),INTENT(INOUT) :: beta(nbnd,nbnd,nkstot)
  !
  INTEGER :: nspin_lsda, ns, nt, nk, ii, jj, ibnd, itetra(4), ik
  REAL(dp) :: thr = 1e-8_dp, V, tsmall(4,4), ei0(4,nbnd), ej0(4,nbnd), e(4), &
  &           ei1(4), ej1(4,nbnd), a(4,4), w0(nbnd,nbnd,4), w1(nbnd,4)
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  DO ns = 1, nspin_lsda
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     END IF
     !
     DO nt = tfst, tlst
        !
        ei0(1:4,1:nbnd) = 0.0_dp
        ej0(1:4,1:nbnd) = 0.0_dp
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           DO ibnd = 1, nbnd
              ei0(1:4,ibnd) = ei0(1:4,ibnd) + wlsm(1:4,ii) * (et_col(ibnd, ik)    - ef)
              ej0(1:4,ibnd) = ej0(1:4,ibnd) + wlsm(1:4,ii) * (et_col(ibnd, ik+iq) - ef)
           END DO
           !
        END DO
        !
        w0(1:nbnd,1:nbnd,1:4) = 0.0_dp
        !
        DO ibnd = 1, nbnd
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4,ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0.0_dp - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0.0_dp .AND. 0.0_dp < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(4, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(4) <= 0.0_dp ) THEN
              !
              ! D - 1
              !
              ei1(1:4) = e(1:4)
              ej1(1:4,1:nbnd) = ej0(itetra(1:4), 1:nbnd)
              !
              CALL dfpt_tetra2_theta(ei1,ej1,w1)
              !
              w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
              &                  + w1(1:nbnd,1:4)
              !
           END IF
           !
        END DO ! ib
        !
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           DO jj = 1, 4
              beta(1:nbnd,1:nbnd,ik) = beta(1:nbnd,1:nbnd,ik) &
              &               + wlsm(jj,ii) * w0(1:nbnd,1:nbnd,jj)
           END DO
           !
        END DO
        !
     END DO ! nt
     !
  END DO ! ns
  !
END SUBROUTINE dfpt_tetra_calc_beta1
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_average_beta(beta)
  !--------------------------------------------------------------------------
  !
  ! Average weights of degenerated states
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd, et
  USE klist, ONLY : nks
  USE qpoint,ONLY : nksq, ikks, ikqs
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(INOUT) :: beta(nbnd,nbnd,nks)
  !
  INTEGER :: ibnd, jbnd, kbnd, ik
  REAL(dp) :: beta2(nbnd)
  !
  DO ik = 1, nksq
     !
     DO ibnd = 1, nbnd
        !
        beta2(1:nbnd) = beta(1:nbnd,ibnd,ikks(ik))
        !
        DO jbnd = ibnd + 1, nbnd
           !
           IF(ABS(et(ibnd,ikks(ik)) - et(jbnd,ikks(ik))) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(1:nbnd,jbnd,ikks(ik))
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 beta(1:nbnd,kbnd,ikks(ik)) = beta2(1:nbnd) / real(jbnd - ibnd, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
        beta2(1:nbnd) = beta(ibnd,1:nbnd,ikks(ik))
        !
        DO jbnd = ibnd + 1, nbnd
           !
           IF(ABS(et(ibnd,ikqs(ik)) - et(jbnd,ikqs(ik))) < 1e-6_dp) THEN
              beta2(1:nbnd) = beta2(1:nbnd) + beta(jbnd,1:nbnd,ikks(ik))
           ELSE
              !
              DO kbnd = ibnd, jbnd - 1
                 beta(kbnd,1:nbnd,ikks(ik)) = beta2(1:nbnd) / real(jbnd - ibnd, dp)
              END DO
              !
              EXIT
           END IF
           !
        END DO
        !
     END DO
  END DO
  !
END SUBROUTINE dfpt_tetra_average_beta
!
!--------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_beta2(iq,tfst,tlst,et_col,beta)
  !--------------------------------------------------------------------------
  !
  ! This routine compute the second term of (B 28) in PRB 64, 235118 (2001).
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE ener,  ONLY : ef
  USE klist, ONLY : nkstot
  USE ktetra, ONLY : ntetra, tetra, wlsm
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: iq, tfst, tlst
  REAL(dp),INTENT(IN) :: et_col(nbnd,nkstot)
  REAL(dp),INTENT(INOUT) :: beta(nbnd,nbnd,nkstot)
  !
  INTEGER :: nspin_lsda, ns, nt, nk, ii, jj, ibnd, itetra(4), ik
  REAL(dp) :: thr = 1e-8_dp, V, tsmall(4,4), ei0(4,nbnd), ej0(4,nbnd), e(4), &
  &           ei1(4), ej1(4,nbnd), a(4,4), w0(nbnd,nbnd,4), w1(nbnd,4)
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  DO ns = 1, nspin_lsda
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     END IF
     !
     DO nt = tfst, tlst
        !
        ei0(1:4,1:nbnd) = 0.0_dp
        ej0(1:4,1:nbnd) = 0.0_dp
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           DO ibnd = 1, nbnd
              !
              ei0(1:4, ibnd) = ei0(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd, ik+iq) - ef)
              ej0(1:4, ibnd) = ej0(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd, ik)    - ef)
              !
           END DO
           !
        END DO
        !
        w0(1:nbnd,1:nbnd,1:4) = 0.0_dp
        !
        DO ibnd = 1, nbnd
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4, ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0.0_dp - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0.0_dp .AND. 0.0_dp < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(4, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_theta(ei1,ej1,w1)
                 !
                 w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
                 &           + V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(4) <= 0.0_dp ) THEN
              !
              ! D - 1
              !
              ei1(1:4) = e(1:4)
              ej1(1:4,1:nbnd) = ej0(itetra(1:4), 1:nbnd)
              !
              CALL dfpt_tetra2_theta(ei1,ej1,w1)
              !
              w0(ibnd,1:nbnd,itetra(1:4)) = w0(ibnd,1:nbnd,itetra(1:4)) &
              &                  + w1(1:nbnd,1:4)
              !
           END IF
           !
        END DO ! ib
        !
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           DO jj = 1, 4
              beta(1:nbnd,1:nbnd,ik) = beta(1:nbnd,1:nbnd,ik) &
              &               + wlsm(jj,ii) * w0(1:nbnd,1:nbnd,jj)
           END DO
           !
        END DO
        !
     END DO ! nt
     !
  END DO ! ns
  !
END SUBROUTINE dfpt_tetra_calc_beta2
!
!----------------------------------------------------------------------------
SUBROUTINE dfpt_tetra_calc_beta3(iq,tfst,tlst,et_col,beta)
  !--------------------------------------------------------------------------
  !
  ! This routine compute the third term of (B 28) in PRB 64, 235118 (2001).
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  USE ener,  ONLY : ef
  USE klist, ONLY : nkstot
  USE ktetra, ONLY : ntetra, tetra, wlsm
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: iq, tfst, tlst
  REAL(dp),INTENT(IN) :: et_col(nbnd,nkstot)
  REAL(dp),INTENT(INOUT) :: beta(nbnd,nbnd,nkstot)
  !
  INTEGER :: nspin_lsda, ns, nt, nk, ii, jj, ibnd, itetra(4), ik
  REAL(dp) :: thr = 1e-8_dp, V, tsmall(4,4), ei0(4,nbnd), ej0(4,nbnd), e(4), &
  &           ei1(4), ej1(4,nbnd), a(4,4), w0(nbnd,nbnd,4), w1(nbnd,4)
  !
  IF ( nspin == 2 ) THEN
     nspin_lsda = 2
  ELSE
     nspin_lsda = 1
  END IF
  !
  DO ns = 1, nspin_lsda
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     IF (ns == 1) THEN
        nk = 0
     ELSE
        nk = nkstot / 2
     END IF
     !
     DO nt = tfst, tlst
        !
        ei0(1:4, 1:nbnd) = 0.0_dp
        ej0(1:4, 1:nbnd) = 0.0_dp
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           DO ibnd = 1, nbnd
              ei0(1:4, ibnd) = ei0(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd, ik)    - ef)
              ej0(1:4, ibnd) = ej0(1:4, ibnd) + wlsm(1:4,ii) * (et_col(ibnd, ik+iq) - ef)
           END DO
           !
        END DO
        !
        w0(1:nbnd,1:nbnd,1:4) = 0.0_dp
        !
        DO ibnd = 1, nbnd
           !
           itetra(1) = 0
           e(1:4) = ei0(1:4, ibnd)
           call hpsort (4, e, itetra)
           !
           DO ii = 1, 4
              a(ii,1:4) = ( 0.0_dp - e(1:4)) / (e(ii) - e(1:4))
           END DO
           !
           IF(e(1) <= 0.0_dp .AND. 0.0_dp < e(2)) THEN
              !
              ! A - 1
              !
              V = a(2,1) * a(3,1) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(4, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) THEN
              !
              ! B - 1
              !
              V = a(3,1) * a(4,1) * a(2,4)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! B - 2
              !
              V = a(3,2) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! B - 3
              !
              V = a(2,3) * a(3,1) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) THEN
              !
              ! C - 1
              !
              V = a(4,3)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! C - 2
              !
              V = a(3,4) * a(4,2)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
              ! C - 3
              !
              V = a(3,4) * a(2,4) * a(4,1)
              !
              IF(V > thr) THEN
                 !
                 tsmall(1, 1:4) = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                 tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
                 tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
                 tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
                 !
                 ei1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
                 ej1(1:4,1:nbnd) = MATMUL(tsmall(1:4,1:4), ej0(itetra(1:4), 1:nbnd))
                 !
                 CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
                 !
                 w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
                 &           - V * MATMUL(w1(1:nbnd,1:4), tsmall(1:4,1:4))
                 !
              END IF
              !
           ELSE IF( e(4) <= 0.0_dp ) THEN
              !
              ! D - 1
              !
              ei1(1:4) = e(1:4)
              ej1(1:4,1:nbnd) = ej0(itetra(1:4), 1:nbnd)
              !
              CALL dfpt_tetra2_lindhard(ei1,ej1,w1)
              !
              w0(1:nbnd,ibnd,itetra(1:4)) = w0(1:nbnd,ibnd,itetra(1:4)) &
              &                  - w1(1:nbnd,1:4)
              !
           END IF
           !
        END DO ! ib
        !
        DO ii = 1, 20
           !
           ik = tetra(ii, nt) + nk
           DO jj = 1, 4
              beta(1:nbnd,1:nbnd,ik) = beta(1:nbnd,1:nbnd,ik) &
              &               + wlsm(jj,ii) * w0(1:nbnd,1:nbnd,jj)
           END DO
           !
        END DO
        !
     END DO ! it
     !
  END DO ! ik
  !
END SUBROUTINE dfpt_tetra_calc_beta3
!
!-----------------------------------------------------------
SUBROUTINE dfpt_tetra2_theta(ei0,ej0,w0)
  !
  ! This routine compute theta(ei - ej)
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei0(4), ej0(4,nbnd)
  REAL(dp),INTENT(OUT) :: w0(nbnd,4)
  !
  INTEGER :: ii, ibnd, itetra(4)
  REAL(dp) :: C(3), e(4), a(4,4), thr = 1.0e-8_dp
  !
  w0(1:nbnd,1:4) = 0.0_dp
  !
  DO ibnd = 1, nbnd
     !
     e(1:4) = ej0(1:4,ibnd) - ei0(1:4)
     itetra(1) = 0
     call hpsort (4, e, itetra)
     !
     DO ii = 1, 4
        a(ii,1:4) = (0_dp - e(1:4)) / (e(ii) - e(1:4))
     END DO
     !
     IF(ABS(e(1)) < thr .AND. ABS(e(4)) < thr) THEN
        !
        ! Theta(0) = 0.5
        !
        w0(ibnd,1:4) = 0.5_dp * 0.25_dp
        !
     ELSE IF((e(1) <= 0.0_dp .AND. 0.0_dp < e(2)) .OR. (e(1) < 0.0_dp .AND. 0.0_dp <= e(2))) THEN
        !
        C(1) = a(2,1) * a(3,1) * a(4,1) * 0.25_dp
        w0(ibnd,itetra(1)) = C(1) * (1.0_dp + a(1,2) + a(1,3) + a(1,4))
        w0(ibnd,itetra(2)) = C(1) * a(2,1)
        w0(ibnd,itetra(3)) = C(1) * a(3,1)
        w0(ibnd,itetra(4)) = C(1) * a(4,1)
        !
     ELSE IF((e(2) <= 0.0_dp .AND. 0.0_dp < e(3)) .OR. (e(2) < 0.0_dp .AND. 0.0_dp <= e(3))) THEN
        !
        C(1) = a(4,1) * a(3,1) * 0.25_dp
        C(2) = a(4,1) * a(3,2) * a(1,3) * 0.25_dp
        C(3) = a(4,2) * a(3,2) * a(1,4) * 0.25_dp
        !
        w0(ibnd,itetra(1)) = C(1) + (C(1) + C(2)) * a(1,3) + (C(1) + C(2) + C(3)) * a(1,4)
        w0(ibnd,itetra(2)) = C(1) + C(2) + C(3) + (C(2) + C(3)) * a(2,3) + C(3) * a(2,4)
        w0(ibnd,itetra(3)) = (C(1) + C(2)) * a(3,1) + (C(2) + C(3)) * a(3,2)
        w0(ibnd,itetra(4)) = (C(1) + C(2) + C(3)) * a(4,1) + C(3) * a(4,2)
        !
     ELSE IF((e(3) <= 0.0_dp .AND. 0.0_dp < e(4)) .OR. (e(3) < 0.0_dp .AND. 0.0_dp <= e(4))) THEN
        !
        C(1) = a(1,4) * a(2,4) * a(3,4)
        !
        w0(ibnd,itetra(1)) = 1.0_dp - C(1) * a(1,4)
        w0(ibnd,itetra(2)) = 1.0_dp - C(1) * a(2,4)
        w0(ibnd,itetra(3)) = 1.0_dp - C(1) * a(3,4)
        w0(ibnd,itetra(4)) = 1.0_dp - C(1) * (1.0_dp + a(4,1) + a(4,2) + a(4,3))
        !
        w0(ibnd,1:4) = w0(ibnd,1:4) * 0.25_dp
        !
     ELSE IF(e(4) <= 0.0_dp) THEN
        !
        ! D - 1
        !
        w0(ibnd,1:4) = 0.25_dp
        !
     END IF
     !
  END DO
  !
END SUBROUTINE dfpt_tetra2_theta
!
!-----------------------------------------------------------------------
SUBROUTINE dfpt_tetra2_lindhard(ei0,ej0,w0)
  !---------------------------------------------------------------------
  !
  ! This routine take the unoccupied region.
  !
  USE kinds, ONLY : dp
  USE wvfct, ONLY : nbnd
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei0(4), ej0(4,nbnd)
  REAL(dp),INTENT(OUT) :: w0(nbnd,4)
  !
  INTEGER :: ii, ibnd, itetra(4)
  REAL(dp) :: V, ei1(4), ej1(4), w1(4), thr = 1e-8_dp
  REAL(dp) :: e(4), a(4,4), tsmall(4,4)
  !
  w0(1:nbnd,1:4) = 0.0_dp
  !
  DO ibnd = 1, nbnd
     !
     e(1:4) = ej0(1:4,ibnd)
     itetra(1) = 0
     call hpsort (4, e, itetra)
     !
     DO ii = 1, 4
        a(ii,1:4) = ( 0.0_dp - e(1:4) ) / (e(ii) - e(1:4))
     END DO
     !
     IF(0_dp <= e(1)) THEN
        !
        ! A - 1
        !
        ej1(1:4) = e(1:4)
        ei1(1:4) = ei0(itetra(1:4))
        !
        CALL dfpt_tetra_lindhard(ei1,ej1,w1)
        !
        w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) + w1(1:4)
        !
     ELSE IF((e(1) < 0.0_dp .AND. 0.0_dp <= e(2)) .OR. (e(1) <= 0.0_dp .AND. 0.0_dp < e(2))) THEN
        !
        ! B - 1
        !
        V = a(1,2)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
           tsmall(2, 1:4) = (/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !       
        END IF
        !
        ! B - 2
        !
        V = a(1,3) * a(2,1)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
           tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        ! B - 3
        !
        V = a(1,4) * a(2,1) * a(3,1)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,2), a(2,1), 0.0_dp, 0.0_dp/)
           tsmall(2, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(3, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !       
        END IF
        !          
     ELSE IF((e(2) < 0.0_dp .AND. 0.0_dp <= e(3)) .OR. (e(2) <= 0.0_dp .AND. 0.0_dp < e(3))) THEN
        !          
        ! C - 1
        !
        V = a(2,4) * a(1,4) * a(3,1)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(2, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
           tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !      
        END IF
        !
        ! C - 2
        !
        V = a(1,3) * a(2,3)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !
        END IF
        !
        ! C - 3
        ! 
        V = a(1,3) * a(2,4) * a(3,2)
        !
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,3), 0.0_dp, a(3,1), 0.0_dp/)
           tsmall(2, 1:4) = (/0.0_dp, a(2,3), a(3,2), 0.0_dp/)
           tsmall(3, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !
        END IF
        !          
     ELSE IF((e(3) < 0.0_dp .AND. 0.0_dp <= e(4)) .OR. (e(3) <= 0.0_dp .AND. 0.0_dp < e(4))) THEN
        !
        ! D - 1
        !
        V = a(3,4) * a(2,4) * a(1,4) 
        !          
        IF(V > thr) THEN
           !
           tsmall(1, 1:4) = (/a(1,4), 0.0_dp, 0.0_dp, a(4,1)/)
           tsmall(2, 1:4) = (/0.0_dp, a(2,4), 0.0_dp, a(4,2)/)
           tsmall(3, 1:4) = (/0.0_dp, 0.0_dp, a(3,4), a(4,3)/)
           tsmall(4, 1:4) = (/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
           !
           ej1(1:4) = MATMUL(tsmall(1:4,1:4), e(1:4))
           ei1(1:4) = MATMUL(tsmall(1:4,1:4), ei0(itetra(1:4)))
           !
           CALL dfpt_tetra_lindhard(ei1,ej1,w1)
           !
           w0(ibnd,itetra(1:4)) = w0(ibnd,itetra(1:4)) &
           &                   + V * MATMUL(w1(1:4), tsmall(1:4,1:4))
           !        
        END IF
        !
     END IF
     !
  END DO
  !
END SUBROUTINE dfpt_tetra2_lindhard
!
!-----------------------------------------------------------------
SUBROUTINE dfpt_tetra_lindhard(ei,ej,w)
  !---------------------------------------------------------------
  !
  ! This routine compute 1 / (e_{k+q} - e_{k})
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: ei(4), ej(4)
  REAL(dp),INTENT(OUT) :: w(4)
  !
  INTEGER :: ii, itetra(4)
  REAL(dp) :: e(4), le(4), thr, thr2
  !
  w(1:4) = 0.0_dp
  !
  itetra(1) = 0
  e(1:4) = ej(1:4) - ei(1:4)
  call hpsort (4, e, itetra)
  !
  thr = MAXVAL(e(1:4)) * 1e-3_dp
  thr2 = 1e-8_dp
  !
  DO ii = 1, 4
     IF(e(ii) < thr2) THEN
        IF(ii == 3) THEN
           CALL errore("dfpt_tetra_lindhard", "Nesting occurs.", 0)
        END IF
        le(ii) = 0.0_dp
        e(ii) = 0.0_dp
     ELSE
        le(ii) = LOG(e(ii))
     END IF
  END DO
  !
  IF(ABS(e(4) - e(3)) < thr ) THEN
     IF(ABS(e(4) - e(2)) < thr ) THEN
        IF(ABS(e(4) - e(1)) < thr ) THEN
           !
           ! e(4) = e(3) = e(2) = e(1)
           !
           w(itetra(4)) = 0.25_dp / e(4)
           w(itetra(3)) = w(itetra(4))
           w(itetra(2)) = w(itetra(4))
           w(itetra(1)) = w(itetra(4))
           !
        ELSE
           !
           ! e(4) = e(3) = e(2)
           !
           w(itetra(4)) = dfpt_tetra_lindhard_1211(e(4),e(1),le(4),le(1))
           w(itetra(3)) = w(itetra(4))
           w(itetra(2)) = w(itetra(4))
           w(itetra(1)) = dfpt_tetra_lindhard_1222(e(1),e(4),le(1),le(4))
           !
           IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
              WRITE(*,'(100e15.5)') e(1:4)
              WRITE(*,'(100e15.5)') w(itetra(1:4))
              CALL errore("dfpt_tetra_lindhard", "4=3=2", 0)
           END IF
           !
        END IF
     ELSE IF(ABS(e(2) - e(1)) < thr ) THEN
        !
        ! e(4) = e(3), e(2) = e(1)
        !
        w(itetra(4)) = dfpt_tetra_lindhard_1221(e(4),e(2), le(4),le(2))
        w(itetra(3)) = w(itetra(4))
        w(itetra(2)) = dfpt_tetra_lindhard_1221(e(2),e(4), le(2),le(4))
        w(itetra(1)) = w(itetra(2))
        !
        IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w(itetra(1:4))
           CALL errore("dfpt_tetra_lindhard", "4=3 2=1", 0)
        END IF
        !
     ELSE
        !
        ! e(4) = e(3)
        !
        w(itetra(4)) = dfpt_tetra_lindhard_1231(e(4),e(1),e(2),le(4),le(1),le(2))
        w(itetra(3)) = w(itetra(4))
        w(itetra(2)) = dfpt_tetra_lindhard_1233(e(2),e(1),e(4),le(2),le(1),le(4))
        w(itetra(1)) = dfpt_tetra_lindhard_1233(e(1),e(2),e(4),le(1),le(2),le(4))
        !
        IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w(itetra(1:4))
           CALL errore("dfpt_tetra_lindhard", "4=3", 0)
        END IF
        !
     END IF
  ELSE IF(ABS(e(3) - e(2)) < thr) THEN
     IF(ABS(e(3) - e(1)) < thr) THEN
        !
        ! e(3) = e(2) = e(1)
        !
        w(itetra(4)) = dfpt_tetra_lindhard_1222(e(4),e(3), le(4),le(3))
        w(itetra(3)) = dfpt_tetra_lindhard_1211(e(3),e(4), le(3),le(4))
        w(itetra(2)) = w(itetra(3))
        w(itetra(1)) = w(itetra(3))
        !
        IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w(itetra(1:4))
           CALL errore("dfpt_tetra_lindhard", "3=2=1", 0)
        END IF
        !
     ELSE
        !
        ! e(3) = e(2)
        !
        w(itetra(4)) = dfpt_tetra_lindhard_1233(e(4),e(1),e(3),le(4),le(1),le(3))
        w(itetra(3)) = dfpt_tetra_lindhard_1231(e(3),e(1),e(4),le(3),le(1),le(4))
        w(itetra(2)) = w(itetra(3))
        w(itetra(1)) = dfpt_tetra_lindhard_1233(e(1),e(4),e(3),le(1),le(4),le(3))
        !
        IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
           WRITE(*,'(100e15.5)') e(1:4)
           WRITE(*,'(100e15.5)') w(itetra(1:4))
           CALL errore("dfpt_tetra_lindhard", "3=2", 0)
        END IF
        !
     END IF
  ELSE IF(ABS(e(2) - e(1)) < thr) THEN
     !
     ! e(2) = e(1)
     !
     w(itetra(4)) = dfpt_tetra_lindhard_1233(e(4),e(3),e(2),le(4),le(3),le(2))
     w(itetra(3)) = dfpt_tetra_lindhard_1233(e(3),e(4),e(2),le(3),le(4),le(2))
     w(itetra(2)) = dfpt_tetra_lindhard_1231(e(2),e(3),e(4),le(2),le(3),le(4))
     w(itetra(1)) = w(itetra(2))
     !
     IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') w(itetra(1:4))
        CALL errore("dfpt_tetra_lindhard", "2=1", 0)
     END IF
     !
  ELSE
     !
     ! DIFferent each other.
     !
     w(itetra(4)) = dfpt_tetra_lindhard_1234(e(4),e(1),e(2),e(3),le(4),le(1),le(2),le(3))
     w(itetra(3)) = dfpt_tetra_lindhard_1234(e(3),e(1),e(2),e(4),le(3),le(1),le(2),le(4))
     w(itetra(2)) = dfpt_tetra_lindhard_1234(e(2),e(1),e(3),e(4),le(2),le(1),le(3),le(4))
     w(itetra(1)) = dfpt_tetra_lindhard_1234(e(1),e(2),e(3),e(4),le(1),le(2),le(3),le(4))
     !      
     IF(ANY(w(itetra(1:4)) < 0.0_dp)) THEN
        WRITE(*,'(100e15.5)') e(1:4)
        WRITE(*,'(100e15.5)') w(itetra(1:4))
        CALL errore("dfpt_tetra_lindhard", "Something wrong.", 0)
     END IF
     !
  END IF
  !
END SUBROUTINE dfpt_tetra_lindhard
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1234(g1,g2,g3,g4,lng1,lng2,lng3,lng4) RESULT(w)
  !--------------------------------------------------------------------------
  !
  ! g1, g2, g3, g4 are different each other
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,g3,g4,lng1,lng2,lng3,lng4
  REAL(dp) :: w
  !
  REAL(dp) :: w2, w3, w4
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1.0_dp)*g2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1.0_dp)*g3/(g3 - g1)
  w4 = ((lng4 - lng1)/(g4 - g1)*g4 - 1.0_dp)*g4/(g4 - g1)
  w2 = ((w2 - w3)*g2)/(g2 - g3)
  w4 = ((w4 - w3)*g4)/(g4 - g3)
  w = (w4 - w2)/(g4 - g2)
  !
END FUNCTION dfpt_tetra_lindhard_1234
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1231(g1,g2,g3,lng1,lng2,lng3) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g1
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(dp) :: w
  !
  REAL(dp) :: w2, w3
  !
  w2 = ((lng2 - lng1)/(g2 - g1)*g2 - 1.0_dp)*g2**2/(g2 - g1) - g1/( &
  &   2.0_dp)
  w2 = w2/(g2 - g1)
  w3 = ((lng3 - lng1)/(g3 - g1)*g3 - 1.0_dp)*g3**2/(g3 - g1) - g1/( &
  &   2.0_dp)
  w3 = w3/(g3 - g1)
  w = (w3 - w2)/(g3 - g2)
  !
END FUNCTION dfpt_tetra_lindhard_1231
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1233(g1,g2,g3,lng1,lng2,lng3) RESULT(w)
  !--------------------------------------------------------------------------
  !
  ! g4 = g3
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,g3,lng1,lng2,lng3
  REAL(dp) :: w
  !
  REAL(dp) :: w2, w3
  !
  w2 = (lng2 - lng1)/(g2 - g1)*g2 - 1.0_dp
  w2 = (g2*w2)/(g2 - g1)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1.0_dp
  w3 = (g3*w3)/(g3 - g1)
  w2 = (w3 - w2)/(g3 - g2)
  w3 = (lng3 - lng1)/(g3 - g1)*g3 - 1.0_dp
  w3 = 1.0_dp - (2.0_dp*w3*g1)/(g3 - g1)
  w3 = w3/(g3 - g1)
  w = (g3*w3 - g2*w2)/(g3 - g2)
  !
END FUNCTION dfpt_tetra_lindhard_1233
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1221(g1,g2,lng1,lng2) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g1 and g3 = g2
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(dp) :: w
  !
  w = 1.0_dp - (lng2 - lng1)/(g2 - g1)*g1
  w = -1.0_dp + (2.0_dp*g2*w)/(g2 - g1)
  w = -1.0_dp + (3.0_dp*g2*w)/(g2 - g1)
  w = w/(2.0_dp*(g2 - g1))
  !
END FUNCTION dfpt_tetra_lindhard_1221
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1222(g1,g2,lng1,lng2) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g3 = g2
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1, g2, lng1, lng2
  REAL(dp) :: w
  !
  w = (lng2 - lng1)/(g2 - g1)*g2 - 1.0_dp
  w = (2.0_dp*g1*w)/(g2 - g1) - 1.0_dp
  w = (3.0_dp*g1*w)/(g2 - g1) + 1.0_dp
  w = w/(2.0_dp*(g2 - g1))
  !
END FUNCTION dfpt_tetra_lindhard_1222
!
!----------------------------------------------------------------------------
FUNCTION dfpt_tetra_lindhard_1211(g1,g2,lng1,lng2) RESULT(w)
  !----------------------------------------------------------------------------
  !
  ! g4 = g3 = g1
  !
  USE kinds, ONLY : dp
  !
  IMPLICIT NONE
  !
  REAL(dp),INTENT(IN) :: g1,g2,lng1,lng2
  REAL(dp) :: w
  !
  w = -1.0_dp + (lng2 - lng1)/(g2 - g1)*g2
  w = -1.0_dp + (2.0_dp*g2*w)/(g2 - g1)
  w = -1.0_dp + (3.0_dp*g2*w)/(2.0_dp*(g2 - g1))
  w = w/(3.0_dp*(g2 - g1))
  !
END FUNCTION dfpt_tetra_lindhard_1211
!
END MODULE dfpt_tetra_mod
