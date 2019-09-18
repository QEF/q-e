  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from PH/dvanqq (QE)
  ! 
  !----------------------------------------------------------------------
  SUBROUTINE dvanqq2()
  !----------------------------------------------------------------------
  !!
  !! This routine calculates two integrals of the Q functions and
  !! its derivatives with V_loc and V_eff which are used
  !! to compute term dV_bare/dtau * psi  in addusdvqpsi.
  !! The result is stored in int1, int2, int4, int5. The routine is called
  !! for each q in nqc. 
  !! int1 -> Eq. B20 of Ref.[1]
  !! int2 -> Eq. B21 of Ref.[1]
  !! int4 -> Eq. B23 of Ref.[1]
  !! int5 -> Eq. B24 of Ref.[1]
  !!
  !! [1] PRB 64, 235118 (2001).
  !! 
  !! RM - Nov/Dec 2014 
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !! Roxana Margine - Dec 2018: Updated based on QE 6.3
  !! SP: Sept. 2019 - Cleaning  
  !!
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE spin_orb,         ONLY : lspinorb
  USE cell_base,        ONLY : tpiba2, omega, tpiba
  USE gvect,            ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
  USE scf,              ONLY : v, vltot
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE phcom,            ONLY : int1, int2, int4, int4_nc, int5, int5_so, & 
                               vlocq
  USE qpoint,           ONLY : xq, eigqts
  USE uspp_param,       ONLY : upf, lmaxq, nh
  USE uspp,             ONLY : okvan, ijtoh
  USE mp_global,        ONLY : intra_pool_comm
  USE mp,               ONLY : mp_sum
  USE fft_base,         ONLY : dfftp
  USE fft_interfaces,   ONLY : fwfft
  USE constants_epw,    ONLY : zero, czero
  !
  IMPLICIT NONE
  !
  !   Local variables
  !
  INTEGER :: na
  !! counter on atoms
  INTEGER :: nb
  !! counter on atoms  
  INTEGER :: ntb
  !! counter on atomic types (species)
  INTEGER :: nta
  !! index of atomic type (specie)
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: ir
  !! counter on FFT mesh points
  INTEGER :: ih
  !! counter on beta functions per atomic type
  INTEGER :: jh
  !! counter on beta functions per atomic type
  INTEGER :: ijh
  !! correspondence beta indexes ih,jh -> composite index ijh
  INTEGER :: ipol
  !! counter on polarizations
  INTEGER :: jpol
  !! counter on polarizations
  INTEGER :: is
  !! counter on spin
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP), ALLOCATABLE :: qmod(:)
  !! the modulus of q+G
  REAL(KIND = DP), ALLOCATABLE :: qmodg(:)
  !! the modulus of G
  REAL(KIND = DP), ALLOCATABLE :: qpg(:, :)
  !! the q+G vectors
  REAL(KIND = DP), ALLOCATABLE :: ylmkq(:, :)
  !! the spherical harmonics at q+G
  REAL(KIND = DP), ALLOCATABLE ::  ylmk0(:, :)
  !! the spherical harmonics at G
  COMPLEX(KIND = DP) :: fact
  !! e^{-i q * \tau} * CONJG(e^{-i q * \tau}) 
  COMPLEX(KIND = DP) :: fact1
  !! -i * omega
  COMPLEX(KIND = DP), EXTERNAL :: ZDOTC
  !! the scalar product function
  COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:), aux2(:), aux3(:), aux5(:)
  !! Auxiallary array
  COMPLEX(KIND = DP), ALLOCATABLE :: sk(:)
  !! 
  COMPLEX(KIND = DP), ALLOCATABLE :: veff(:, :)
  !! effective potential
  COMPLEX(KIND = DP), ALLOCATABLE, TARGET :: qgm(:)
  !! the augmentation function at G
  COMPLEX(KIND = DP), POINTER :: qgmq(:)
  !! the augmentation function at q+G
  ! 
  IF (.NOT. okvan) RETURN
  !
  CALL start_clock('dvanqq2')
  ! 
  int1(:, :, :, :, :) = czero
  int2(:, :, :, :, :) = czero
  int4(:, :, :, :, :) = czero
  int5(:, :, :, :, :) = czero
  !
  ALLOCATE(sk(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating sk', 1)
  ALLOCATE(aux1(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux1', 1)
  ALLOCATE(aux2(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux2', 1)
  ALLOCATE(aux3(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux3', 1)
  ALLOCATE(aux5(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux5', 1)
  ALLOCATE(qmodg(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qmodg', 1)
  ALLOCATE(qmod(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qmod', 1)
  ALLOCATE(qgmq(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qgmq', 1)
  ALLOCATE(qgm(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qgm', 1)
  ALLOCATE(ylmk0(ngm, lmaxq * lmaxq), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating ylmk0', 1)
  ALLOCATE(ylmkq(ngm, lmaxq * lmaxq), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating ylmkq', 1)
  sk(:) = czero
  aux1(:) = czero
  aux2(:) = czero
  aux3(:) = czero
  aux5(:) = czero
  qmodg(:) = zero
  qmod(:) = zero
  qgmq(:) = czero
  qgm(:) = czero
  ylmk0(:, :) = zero
  ylmkq(:, :) = zero
  !
  ! compute spherical harmonics
  !
  CALL ylmr2(lmaxq * lmaxq, ngm, g, gg, ylmk0)
  !
  DO ig = 1, ngm
    qmodg(ig) = SQRT(gg(ig))
  ENDDO
  ! 
  ALLOCATE(qpg(3, ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qpg', 1)
  qpg(:, :) = zero
  !
  CALL setqmod(ngm, xq, g, qmod, qpg)
  CALL ylmr2(lmaxq * lmaxq, ngm, qpg, qmod, ylmkq)
  !
  DEALLOCATE(qpg, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating qpg', 1)
  DO ig = 1, ngm
    qmod(ig) = SQRT(qmod(ig))
  ENDDO
  !
  !   we start by computing the FT of the effective potential
  !
  ALLOCATE(veff(dfftp%nnr, nspin_mag), STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating veff', 1)
  veff(:, :) = czero
  !
  DO is = 1, nspin_mag
    IF (nspin_mag /= 4 .OR. is == 1) THEN
      DO ir = 1, dfftp%nnr
        veff(ir, is) = CMPLX(vltot(ir) + v%of_r(ir, is), zero, KIND = DP)
      ENDDO
    ELSE
      DO ir = 1, dfftp%nnr
        veff(ir, is) = CMPLX(v%of_r(ir, is), zero, KIND = DP)
      ENDDO
    ENDIF
    CALL fwfft('Rho', veff(:, is), dfftp)
  ENDDO
  !
  ! We compute here two of the three integrals needed in the phonon
  !
  fact1 = CMPLX(0.d0, - tpiba * omega, KIND = DP)
  !
  DO ntb = 1, ntyp
    IF (upf(ntb)%tvanp) THEN
      !
      DO ih = 1, nh(ntb)
        DO jh = ih, nh(ntb)
          ijh = ijtoh(ih, jh, ntb)
          !
          ! Compute the augmentation function
          !
          CALL qvan2(ngm, ih, jh, ntb, qmodg, qgm, ylmk0)
          CALL qvan2(ngm, ih, jh, ntb, qmod, qgmq, ylmkq)
          !
          ! NB: for this integral the moving atom and the atom of Q
          ! do not necessarily coincide
          !
          DO nb = 1, nat
            IF (ityp(nb) == ntb) THEN
              DO ig = 1, ngm
                aux1(ig) = qgmq(ig) * eigts1(mill(1, ig), nb) &
                                    * eigts2(mill(2, ig), nb) &
                                    * eigts3(mill(3, ig), nb)
              ENDDO
              !
              DO na = 1, nat
                fact = eigqts(na) * CONJG(eigqts(nb))
                !
                !    nb is the atom of the augmentation function
                !
                nta = ityp(na)
                DO ig = 1, ngm
                  sk(ig) = vlocq(ig, nta) * eigts1(mill(1, ig), na) &
                                          * eigts2(mill(2, ig), na) &
                                          * eigts3(mill(3, ig), na) 
                ENDDO
                !
                DO ipol = 1, 3
                  DO ig = 1, ngm
                    aux5(ig) = sk(ig) * (g(ipol, ig) + xq(ipol))
                  ENDDO
                  int2(ih, jh, ipol, na, nb) = fact * fact1 * ZDOTC(ngm, aux1, 1, aux5, 1)
                  ! 
                  DO jpol = 1, 3
                    IF (jpol >= ipol) THEN
                      DO ig = 1, ngm
                        aux3(ig) = aux5(ig) * (g(jpol, ig) + xq(jpol))
                      ENDDO
                      int5(ijh, ipol, jpol, na, nb) = CONJG(fact) * tpiba2 * omega * &
                         ZDOTC(ngm, aux3, 1, aux1, 1)
                    ELSE
                      int5(ijh, ipol, jpol, na, nb) = int5(ijh, jpol, ipol, na, nb)
                    ENDIF
                  ENDDO
                ENDDO !ipol
                !
              ENDDO !na
              !
              DO ig = 1, ngm
                aux1(ig) = qgm(ig) * eigts1(mill(1,ig),nb) &
                                   * eigts2(mill(2,ig),nb) &
                                   * eigts3(mill(3,ig),nb)
              ENDDO
              !
              DO is = 1, nspin_mag
                DO ipol = 1, 3
                  DO ig = 1, ngm
                    aux2(ig) = veff(dfftp%nl(ig), is) * g(ipol, ig)
                  ENDDO
                  int1(ih, jh, ipol, nb, is) = - fact1 * ZDOTC(ngm, aux1, 1, aux2, 1)
                  DO jpol = 1, 3
                    IF (jpol >= ipol) THEN
                      DO ig = 1, ngm
                         aux3(ig) = aux2(ig) * g(jpol,ig)
                      ENDDO
                      int4(ijh, ipol, jpol, nb, is) = - tpiba2 * omega * ZDOTC(ngm, aux3, 1, aux1, 1)
                    ELSE
                      int4(ijh, ipol, jpol, nb, is) = int4(ijh, jpol, ipol, nb, is)
                    ENDIF
                  ENDDO ! jpol
                ENDDO ! ipol
              ENDDO ! is
            ENDIF ! ityp
          ENDDO ! nb
        ENDDO ! jh
      ENDDO ! ih
      !
      DO ih = 1, nh(ntb)
        DO jh = ih + 1, nh(ntb)
          !
          !    We use the symmetry properties of the integral factor
          !
          DO nb = 1, nat
            IF (ityp(nb) == ntb) THEN
              DO ipol = 1, 3
                DO is = 1, nspin_mag
                  int1(jh, ih, ipol, nb, is) = int1(ih, jh, ipol, nb, is)
                ENDDO
                DO na = 1, nat
                  int2(jh, ih, ipol, na, nb) = int2(ih, jh, ipol, na, nb)
                ENDDO ! na
              ENDDO ! ipol
            ENDIF
          ENDDO ! nb
        ENDDO ! jh
      ENDDO ! ih
    ENDIF ! upf
  ENDDO ! ntb
  CALL mp_sum(int1, intra_pool_comm)
  CALL mp_sum(int2, intra_pool_comm)
  CALL mp_sum(int4, intra_pool_comm)
  CALL mp_sum(int5, intra_pool_comm)
  !
  IF (noncolin) THEN
    CALL set_int12_nc(0)
    int4_nc = czero
    IF (lspinorb) int5_so = czero
    DO ntb = 1, ntyp
      IF (upf(ntb)%tvanp) THEN
        DO na = 1, nat
          IF (ityp(na) == ntb) THEN
            IF (upf(ntb)%has_so)  THEN
              CALL transform_int4_so(int4, na)
              CALL transform_int5_so(int5, na)
            ELSE
              CALL transform_int4_nc(int4, na)
              IF (lspinorb) CALL transform_int5_nc(int5,na)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO ! ntb
  ENDIF
  !
  DEALLOCATE(sk, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating sk', 1)
  DEALLOCATE(aux1, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating aux1', 1)
  DEALLOCATE(aux2, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating aux2', 1)
  DEALLOCATE(aux3, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating aux3', 1)
  DEALLOCATE(aux5, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating aux5', 1)
  DEALLOCATE(qmodg, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating qmodg', 1)
  DEALLOCATE(qmod, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating qmod', 1)
  DEALLOCATE(qgmq, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating qgmq', 1)
  DEALLOCATE(qgm, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating qgm', 1)
  DEALLOCATE(ylmk0, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating ylmk0', 1)
  DEALLOCATE(ylmkq, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating ylmkq', 1)
  DEALLOCATE(veff, STAT = ierr)
  IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating veff', 1)
  !
  CALL stop_clock('dvanqq2')
  RETURN
  ! 
  !----------------------------------------------------------------------
  END SUBROUTINE dvanqq2
  !----------------------------------------------------------------------
