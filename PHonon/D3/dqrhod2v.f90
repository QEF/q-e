!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dqrhod2v (ipert, drhoscf)
  !-----------------------------------------------------------------------
  ! calculates the term containing the second variation of the potential
  ! and the first variation of the charge density with respect to a
  ! perturbation at a generic q
  !
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE pwcom
  USE uspp,                 ONLY : vkb, dvan
  USE uspp_param,           ONLY : nh
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : iunigk
  use qpoint,               ONLY : xq, npwq, nksq, igkq
  USE phcom
  USE d3com
  USE mp_global,            ONLY : my_pool_id
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: ipert
  ! index of the perturbation associated with drho
  COMPLEX (DP) :: drhoscf (dfftp%nnr)
  ! the variation of the charge density
  !
  ! local variables
  !
  INTEGER :: icart, jcart, na_icart, na_jcart, na, ng, nt, &
       ik, ikk, ikq, ig, ibnd, nu_i, nu_j, nu_k, ikb, jkb, nrec, ios
  ! counters

  REAL (DP) :: gtau, wgg
  ! the product G*\tau_s
  ! the weight of a K point

  COMPLEX (DP) :: zdotc, fac, alpha (8), work
  COMPLEX (DP), ALLOCATABLE :: d3dywrk (:,:), work0 (:), &
       work1 (:), work2 (:), work3 (:), work4 (:), work5 (:), work6 (:)
  ! work space

  ALLOCATE  (d3dywrk( 3 * nat, 3 * nat))
  ALLOCATE  (work0( dfftp%nnr))
  ALLOCATE  (work1( npwx))
  ALLOCATE  (work2( npwx))
  ALLOCATE  (work3( npwx))
  ALLOCATE  (work4( npwx))
  ALLOCATE  (work5( npwx))
  ALLOCATE  (work6( npwx))

  d3dywrk (:,:) = (0.d0, 0.d0)
  !
  ! Here the contribution deriving from the local part of the potential
  !
  !   ... computed only by the first pool (no sum over k needed)
  !
  IF ( my_pool_id == 0 ) THEN
     !
     work0 (:) = drhoscf(:)
     CALL fwfft ('Dense', work0, dfftp)
     DO na = 1, nat
        DO icart = 1, 3
           na_icart = 3 * (na - 1) + icart
           DO jcart = 1, 3
              na_jcart = 3 * (na - 1) + jcart
              DO ng = 1, ngm
                 gtau = tpi * ( (xq (1) + g (1, ng) ) * tau (1, na) + &
                      (xq (2) + g (2, ng) ) * tau (2, na) + &
                      (xq (3) + g (3, ng) ) * tau (3, na) )
                 fac = CMPLX(COS (gtau), - SIN (gtau) ,kind=DP)
                 d3dywrk (na_icart, na_jcart) = d3dywrk (na_icart, na_jcart) &
                      - tpiba2 * omega * (xq (icart) + g (icart, ng) ) * &
                      (xq (jcart) + g (jcart, ng) ) * &
                      vlocq (ng, ityp (na) ) * fac * CONJG (work0 (nl (ng) ) )
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
     CALL mp_sum(  d3dywrk, intra_pool_comm )
     !
  END IF
  !
  ! each pool contributes to next term
  !
  ! Here we compute the nonlocal (Kleinman-Bylander) contribution.
  !
  REWIND (unit = iunigk)

  DO ik = 1, nksq
     READ (iunigk, err = 200, iostat = ios) npw, igk
200  CALL errore ('dqrhod2v', 'reading igk', ABS (ios) )
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
        npwq = npw
     ELSE
        ikk = 2 * ik - 1
        ikq = 2 * ik
        READ (iunigk, err = 300, iostat = ios) npwq, igkq
300     CALL errore ('dqrhod2v', 'reading igkq', ABS (ios) )
     ENDIF
     wgg = wk (ikk)
     CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
     !
     ! In metallic case it necessary to know the wave function at k+q point
     ! so as to correct dpsi. dvpsi is used as working array
     !
     IF (degauss /= 0.d0) CALL davcio (dvpsi, lrwfc, iuwfc, ikq, -1)
     CALL init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     CALL init_us_2 (npw, igk, xk (1, ikk), vkb0)
     !
     ! Reads the first variation of the wavefunction projected on conduction
     !
     nrec = (ipert - 1) * nksq + ik
     CALL davcio (dpsi, lrdwf, iudqwf, nrec, - 1)
     !
     ! In the metallic case corrects dpsi so as that the density matrix
     ! will be:   Sum_{k,nu} 2 * | dpsi > < psi |
     !
     IF (degauss /= 0.d0) THEN
        nrec = ipert + (ik - 1) * 3 * nat
        CALL davcio (psidqvpsi, lrpdqvp, iupdqvp, nrec, - 1)
        CALL dpsi_corr (dvpsi, psidqvpsi, ikk, ikq, ipert)
     ENDIF
     !
     DO icart = 1, 3
        DO jcart = 1, 3
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 work1(ig)=evc(ig,ibnd)*tpiba*(xk(icart,ikk)+g(icart,igk(ig)))
                 work2(ig)=evc(ig,ibnd)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work5(ig)=   work1(ig)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
              ENDDO
              DO ig = 1, npwq
                 work3(ig)=dpsi(ig,ibnd)*tpiba*(xk(icart,ikq)+g(icart,igkq(ig)))
                 work4(ig)=dpsi(ig,ibnd)*tpiba*(xk(jcart,ikq)+g(jcart,igkq(ig)))
                 work6(ig)=    work3(ig)*tpiba*(xk(jcart,ikq)+g(jcart,igkq(ig)))
              ENDDO
              jkb=0
              DO nt = 1, ntyp
                 DO na = 1, nat
                    IF (ityp (na).EQ.nt) THEN
                       na_icart = 3 * (na - 1) + icart
                       na_jcart = 3 * (na - 1) + jcart
                       DO ikb = 1, nh (nt)
                          jkb = jkb+1
                          alpha(1) = zdotc(npw, work1, 1,vkb0(1,jkb), 1)
                          alpha(2) = zdotc(npwq,vkb(1,jkb), 1, work4, 1)
                          alpha(3) = zdotc(npw, work2, 1,vkb0(1,jkb), 1)
                          alpha(4) = zdotc(npwq,vkb(1,jkb), 1, work3, 1)
                          alpha(5) = zdotc(npw, work5, 1,vkb0(1,jkb), 1)
                          alpha(6) = zdotc(npwq,vkb(1,jkb),1,dpsi(1,ibnd),1)
                          alpha(7) = zdotc(npw, evc(1,ibnd),1,vkb0(1,jkb),1)
                          alpha(8) = zdotc(npwq,vkb(1,jkb),1,work6, 1)
                          !
                          CALL mp_sum(  alpha, intra_pool_comm )
                          !
                          d3dywrk(na_icart,na_jcart) = d3dywrk(na_icart,na_jcart) &
                            + CONJG(alpha(1) * alpha(2) + alpha(3) * alpha(4) - &
                                    alpha(5) * alpha(6) - alpha(7) * alpha(8) ) &
                                    * dvan (ikb, ikb, nt) * wgg * 2.0d0
                       ENDDO
                    ENDIF
                 ENDDO
              END DO
           END DO
        ENDDO
     ENDDO
  ENDDO
  !
  CALL mp_sum( d3dywrk, inter_pool_comm )
  !
  !  Rotate the dynamical matrix on the basis of patterns
  !  some indices do not need to be rotated
  !
  nu_k = ipert
  DO nu_i = 1, 3 * nat
     IF (q0mode (nu_i) ) THEN
        DO nu_j = 1, 3 * nat
           work = (0.0d0, 0.0d0)
           DO na = 1, nat
              DO icart = 1, 3
                 na_icart = 3 * (na - 1) + icart
                 DO jcart = 1, 3
                    na_jcart = 3 * (na - 1) + jcart
                    work = work + ug0 (na_icart, nu_i) * &
                         d3dywrk (na_icart,na_jcart) * u (na_jcart, nu_j)
                 ENDDO
              ENDDO
           ENDDO
           d3dyn (nu_i, nu_k, nu_j) = d3dyn (nu_i, nu_k, nu_j) + work
           d3dyn (nu_i, nu_j, nu_k) = d3dyn (nu_i, nu_j, nu_k) + CONJG(work)
        ENDDO
     ENDIF
  ENDDO

  DEALLOCATE (work6)
  DEALLOCATE (work5)
  DEALLOCATE (work4)
  DEALLOCATE (work3)
  DEALLOCATE (work2)
  DEALLOCATE (work1)
  DEALLOCATE (work0)
  DEALLOCATE (d3dywrk)

  RETURN
END SUBROUTINE dqrhod2v
