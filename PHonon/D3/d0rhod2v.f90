!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE d0rhod2v (ipert, drhoscf)
  !-----------------------------------------------------------------------
  ! calculates the term containing the second variation of the potential
  ! and the first variation of the charge density with respect to a
  ! perturbation at q=0
  !
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp, tau
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : iunigk
  USE kinds,                 ONLY : DP
  USE uspp,                  ONLY : dvan
  USE uspp_param,            ONLY : nh
  USE fft_base,              ONLY : dfftp
  USE fft_interfaces,        ONLY : fwfft
  USE pwcom
  USE wavefunctions_module,  ONLY : evc
  USE phcom
  USE d3com
  USE mp_global,             ONLY : my_pool_id, inter_pool_comm, intra_pool_comm
  USE mp,                    ONLY : mp_sum

  USE qpoint,     ONLY : igkq, nksq, npwq
  USE control_lr, ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER :: ipert              ! index of the perturbation associated with drho
  COMPLEX (DP) :: drhoscf (dfftp%nnr)! the variation of the charge density
  !
  INTEGER :: icart,           & ! counter on polarizations
             jcart,           & ! counter on polarizations
             na_icart,        & ! counter on modes
             na_jcart,        & ! counter on modes
             na,              & ! counter on atoms
             ng,              & ! counter on G vectors
             nt,              & ! counter on atomic types
             ik,              & ! counter on k points
             ikk,             & ! counter on k points
             ig,              & ! counter on G vectors
             ibnd,            & ! counter on bands
             nu_i,            & ! counter on modes
             nu_j,            & ! counter on modes
             nu_k,            & ! counter on modes
             ikb, jkb,        & ! counter on beta functions
             nrec,            & ! record position of dwfc
             ios                ! integer variable for I/O control

  REAL (DP) :: gtau,           & ! the product G*\tau_s
              wgg               ! the weight of a K point

  COMPLEX (DP) :: zdotc, d3dywrk (3*nat,3*nat), fac, alpha(8), work
  COMPLEX (DP), ALLOCATABLE :: work0 (:), work1 (:), work2 (:), &
                                      work3 (:), work4 (:), work5 (:), &
                                      work6 (:)
  ! auxiliary space

  ALLOCATE (work0(dfftp%nnr))
  ALLOCATE (work1(npwx))
  ALLOCATE (work2(npwx))
  ALLOCATE (work3(npwx))
  ALLOCATE (work4(npwx))
  ALLOCATE (work5(npwx))
  ALLOCATE (work6(npwx))

  d3dywrk (:,:) = (0.d0, 0.d0)
  !
  ! Here the contribution deriving from the local part of the potential
  !
  IF ( my_pool_id == 0 ) THEN
     !
     !   ... computed only by the first pool (no sum over k needed)
     !
     work0 (:) = drhoscf (:)
     CALL fwfft ('Dense', work0, dfftp)
     DO na = 1, nat
        DO icart = 1,3
           na_icart = 3*(na-1)+icart
           DO jcart = 1,3
              na_jcart = 3*(na-1)+jcart
              DO ng = 1, ngm
                 gtau = tpi * ( g(1,ng)*tau(1,na) + &
                      g(2,ng)*tau(2,na) + &
                      g(3,ng)*tau(3,na) )

                 fac = CMPLX(COS(gtau),SIN(gtau),kind=DP)

                 d3dywrk(na_icart,na_jcart) = &
                      d3dywrk(na_icart,na_jcart) - &
                      tpiba2 * g(icart,ng) * g(jcart,ng) * &
                      omega * vloc(igtongl(ng),ityp(na)) * &
                      fac*work0(nl(ng))
              ENDDO
           ENDDO
        ENDDO
        WRITE( stdout,*) na
        WRITE( stdout,'(3(2f10.6,2x))') &
             ((d3dywrk(3*(na-1)+icart,3*(na-1)+jcart), &
             jcart=1,3),icart=1,3)
     ENDDO

     CALL mp_sum( d3dywrk, intra_pool_comm )
     !
  END IF
  !
  ! each pool contributes to next term
  !
  ! Here we compute the nonlocal (Kleinman-Bylander) contribution.
  !
  REWIND (unit=iunigk)

  DO ik = 1, nksq
     READ (iunigk, err = 200, iostat = ios) npw, igk
200  CALL errore ('d0rhod2v', 'reading igk', ABS (ios) )
     IF (lgamma) THEN
        ikk = ik
        npwq = npw
     ELSE
        ikk = 2 * ik - 1
        READ (iunigk, err = 300, iostat = ios) npwq, igkq
300     CALL errore ('d0rhod2v', 'reading igkq', ABS (ios) )
        npwq = npw
     ENDIF
     wgg = wk (ikk)
     CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)

     CALL init_us_2 (npw, igk, xk (1, ikk), vkb0)
     !
     ! Reads the first variation of the wavefunction projected on conduction
     !
     nrec = (ipert - 1) * nksq + ik
     CALL davcio (dpsi, lrdwf, iudwf, nrec, - 1)
     !
     ! In the metallic case corrects dpsi so as that the density matrix
     ! will be:   Sum_{k,nu} 2 * | dpsi > < psi |
     !
     IF (degauss /= 0.d0) THEN
        nrec = ipert + (ik - 1) * 3 * nat
        CALL davcio (psidqvpsi, lrpdqvp, iupd0vp, nrec, - 1)
        CALL dpsi_corr (evc, psidqvpsi, ikk, ikk, ipert)
     ENDIF
     DO icart = 1, 3
        DO jcart = 1, 3
           DO ibnd = 1, nbnd
              DO ig = 1, npw
                 work1(ig)= evc(ig,ibnd)*tpiba*(xk(icart,ikk)+g(icart,igk(ig)))
                 work2(ig)= evc(ig,ibnd)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work3(ig)=dpsi(ig,ibnd)*tpiba*(xk(icart,ikk)+g(icart,igk(ig)))
                 work4(ig)=dpsi(ig,ibnd)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work5(ig)=    work1(ig)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work6(ig)=    work3(ig)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
              ENDDO
              jkb=0
              DO nt = 1, ntyp
                 DO na = 1, nat
                    IF (ityp (na) == nt) THEN
                       na_icart = 3 * (na - 1) + icart
                       na_jcart = 3 * (na - 1) + jcart
                       DO ikb = 1, nh (nt)
                          jkb=jkb+1
                          alpha (1) = zdotc (npw, work1, 1, vkb0(1,jkb), 1)
                          alpha (2) = zdotc (npw, vkb0(1,jkb), 1, work4, 1)
                          alpha (3) = zdotc (npw, work2, 1, vkb0(1,jkb), 1)
                          alpha (4) = zdotc (npw, vkb0(1,jkb), 1, work3, 1)
                          alpha (5) = zdotc (npw, work5, 1, vkb0(1,jkb), 1)
                          alpha (6) = zdotc (npw, vkb0(1,jkb), 1, dpsi (1,ibnd), 1)
                          alpha (7) = zdotc (npw,  evc (1,ibnd), 1, vkb0(1,jkb), 1)
                          alpha (8) = zdotc (npw, vkb0(1,jkb), 1, work6, 1)
#ifdef __MPI
                          CALL mp_sum(  alpha, intra_pool_comm )
#endif
                          d3dywrk (na_icart, na_jcart) = d3dywrk (na_icart, na_jcart) &
                               + (alpha(1)*alpha(2) + alpha(3)*alpha(4) &
                                - alpha(5)*alpha(6) - alpha(7)*alpha(8)) * &
                                  dvan (ikb,ikb,nt) * wgg * 2.0d0
                       ENDDO
                    END IF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  CALL mp_sum ( d3dywrk, inter_pool_comm )
  !
  !  Rotate the dynamical matrix on the basis of patterns
  !  first index does not need to be rotated
  !
  nu_k = ipert
  DO nu_i = 1, 3 * nat
     DO nu_j = 1, 3 * nat
        work = (0.0d0, 0.0d0)
        DO na = 1, nat
           DO icart = 1, 3
              na_icart = 3 * (na-1) + icart
              DO jcart = 1, 3
                 na_jcart = 3 * (na-1) + jcart
                 work = work + CONJG(u(na_icart,nu_i)) * &
                               d3dywrk(na_icart,na_jcart) * &
                               u(na_jcart,nu_j)
              ENDDO
           ENDDO
        ENDDO
        d3dyn(nu_k,nu_i,nu_j) = d3dyn(nu_k,nu_i,nu_j) + work
        IF (allmodes) THEN
           d3dyn(nu_j,nu_k,nu_i) = d3dyn(nu_j,nu_k,nu_i) + work
           d3dyn(nu_i,nu_j,nu_k) = d3dyn(nu_i,nu_j,nu_k) + work
        ENDIF
     ENDDO
  ENDDO

  DEALLOCATE (work6)
  DEALLOCATE (work5)
  DEALLOCATE (work4)
  DEALLOCATE (work3)
  DEALLOCATE (work2)
  DEALLOCATE (work1)
  DEALLOCATE (work0)

  RETURN
END SUBROUTINE d0rhod2v
