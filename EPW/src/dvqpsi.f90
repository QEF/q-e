  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE dvqpsi
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines to computes dV_bare/dtau * psi or its components.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !----------------------------------------------------------------------
    SUBROUTINE dvqpsi_us3(ik, uact, addnlcc, xxkq, xq0, igk, igkq, npw, npwq, evc)
    !----------------------------------------------------------------------
    !!
    !! This routine calculates dV_bare/dtau * psi for one perturbation
    !! with a given q. The displacements are described by a vector u.
    !! The result is stored in dvpsi. The routine is called for each k point
    !! and for each pattern u. It computes simultaneously all the bands.
    !! It implements Eq. B29 of PRB 64, 235118 (2001). The contribution
    !! of the local pseudopotential is calculated here, that of the nonlocal
    !! pseudopotential in dvqpsi_us_only3.
    !! Adapted from PH/dvqpsi_us (QE)
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !! RM - Jan 2019
    !! Updated based on QE 6.3
    !!
    !! SP - Feb 2020
    !! Pass the wfc at k (evc) explicitely to work with noncolin rotation.
    !!
    USE kinds,                 ONLY : DP
    USE pwcom,                 ONLY : nbnd
    USE ions_base,             ONLY : nat, ityp
    USE cell_base,             ONLY : tpiba
    USE fft_base,              ONLY : dfftp, dffts
    USE fft_interfaces,        ONLY : fwfft, invfft
    USE gvect,                 ONLY : eigts1, eigts2, eigts3, mill, g, ngm
    USE gvecs,                 ONLY : ngms
    USE lsda_mod,              ONLY : lsda
    USE scf,                   ONLY : rho, rho_core
    USE noncollin_module,      ONLY : nspin_lsda, nspin_gga, npol
    use uspp_param,            ONLY : upf
    USE wvfct,                 ONLY : npwx
    USE nlcc_ph,               ONLY : drc
    USE uspp,                  ONLY : nlcc_any
    USE eqv,                   ONLY : dvpsi, dmuxc, vlocq
    USE qpoint,                ONLY : eigqts
    USE klist_epw,             ONLY : isk_loc
    USE gc_lr,                 ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
    USE funct,                 ONLY : dft_is_nonlocc
    USE xc_lib,                ONLY : xclib_dft_is
    USE elph2,                 ONLY : lower_band, upper_band, ibndstart
    USE constants_epw,         ONLY : czero, eps12
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: addnlcc
    !! True if NLCC is present
    INTEGER, INTENT(in) :: ik
    !! Counter on k-point
    INTEGER, INTENT(in) :: npw
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER, INTENT(in) :: npwq
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER, INTENT(in) :: igk(npw)
    !! k+G mapping
    INTEGER, INTENT(in) :: igkq(npwq)
    !! k+G+q mapping
    REAL(KIND = DP), INTENT(in) :: xq0(3)
    !! Current coarse q-point coordinate
    REAL(KIND = DP), INTENT(in) :: xxkq(3)
    !! k+q point coordinate
    COMPLEX(KIND = DP), INTENT(in) :: uact(3 * nat)
    !! the pattern of displacements
    COMPLEX(KIND = DP), INTENT(in) :: evc(npwx * npol, nbnd)
    !! Wavefunction at k
    !
    ! Local variables
    INTEGER :: na
    !! counter on atoms
    INTEGER :: mu
    !! counter on modes
    INTEGER :: ig
    !! counter on G vectors
    INTEGER :: nt
    !! counter on atomic types
    INTEGER :: ibnd
    !! counter on bands
    INTEGER ::  ir
    !! counter on real mesh
    INTEGER :: is
    !! counter on spin
    INTEGER :: ip
    !! counter on polarizations
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: fac
    !! spin degeneracy factor
    COMPLEX(KIND = DP) :: gtau
    !! e^{-i G * \tau}
    COMPLEX(KIND = DP) :: u1, u2, u3
    !! components of displacement pattern u
    COMPLEX(KIND = DP) :: gu0
    !! scalar product q * u
    COMPLEX(KIND = DP) :: gu
    !! q * u + G * u
    COMPLEX(KIND = DP) :: fact
    !! e^{-i q * \tau}
    COMPLEX(KIND = DP), ALLOCATABLE, TARGET :: aux(:)
    !! Auxillary variable
    COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:), aux2(:)
    !! Auxillary variable
    COMPLEX(KIND = DP), POINTER :: auxs(:)
    !! Auxiallary pointer
    COMPLEX(KIND = DP), ALLOCATABLE :: drhoc(:)
    !! response core charge density
    !
    CALL start_clock('dvqpsi_us3')
    !
    IF (nlcc_any .AND. addnlcc) THEN
      ALLOCATE(drhoc(dfftp%nnr), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error allocating drhoc', 1)
      ALLOCATE(aux(dfftp%nnr), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error allocating aux', 1)
      ALLOCATE(auxs(dffts%nnr), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error allocating auxs', 1)
    ENDIF
    ALLOCATE(aux1(dffts%nnr), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error allocating aux1', 1)
    ALLOCATE(aux2(dffts%nnr), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error allocating aux2', 1)
    !
    !    We start by computing the contribution of the local potential.
    !    The computation of the derivative of the local potential is done in
    !    reciprocal space while the product with the wavefunction is done in real space
    !
    dvpsi(:, :) = czero
    aux1(:) = czero
    DO na = 1, nat
      fact = tpiba * (0.d0, -1.d0) * eigqts(na)
      mu = 3 * (na - 1)
      u1 = uact(mu + 1)
      u2 = uact(mu + 2)
      u3 = uact(mu + 3)
      IF (ABS(u1) + ABS(u2) + ABS(u3) > eps12) THEN
        nt = ityp(na)
        gu0 = xq0(1) * u1 + xq0(2) * u2 + xq0(3) * u3
        DO ig = 1, ngms
          gtau = eigts1(mill(1, ig), na) * &
                 eigts2(mill(2, ig), na) * &
                 eigts3(mill(3, ig), na)
          gu = gu0 + g(1, ig) * u1 + g(2, ig) * u2 + g(3, ig) * u3
          aux1(dffts%nl(ig)) = aux1(dffts%nl(ig)) + vlocq(ig, nt) * gu * fact * gtau
        ENDDO
      ENDIF
    ENDDO
    !
    ! add NLCC when present
    !
    IF (nlcc_any .AND. addnlcc) THEN
      drhoc(:) = czero
      DO na = 1, nat
        fact = tpiba * (0.d0, -1.d0) * eigqts(na)
        mu = 3 * (na - 1)
        u1 = uact(mu + 1)
        u2 = uact(mu + 2)
        u3 = uact(mu + 3)
        IF (ABS(u1) + ABS(u2) + ABS(u3) > eps12) THEN
          nt = ityp(na)
          gu0 = xq0(1) * u1 + xq0(2) * u2 + xq0(3) * u3
          IF (upf(nt)%nlcc) THEN
            DO ig = 1, ngm
              gtau = eigts1(mill(1, ig), na) * &
                     eigts2(mill(2, ig), na) * &
                     eigts3(mill(3, ig), na)
              gu = gu0 + g(1, ig) * u1 + g(2, ig) * u2 + g(3, ig) * u3
              drhoc(dfftp%nl(ig)) = drhoc(dfftp%nl(ig)) + drc(ig, nt) * gu * fact * gtau
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      !
      CALL invfft('Rho', drhoc, dfftp)
      !
      aux(:) = czero
      IF (.NOT. lsda) THEN
        DO ir = 1, dfftp%nnr
          aux(ir) = drhoc(ir) * dmuxc(ir, 1, 1)
        ENDDO
      ELSE
        is = isk_loc(ik)
        DO ir = 1, dfftp%nnr
          aux(ir) = drhoc(ir) * 0.5d0 * (dmuxc(ir, is, 1) + dmuxc(ir, is, 2))
        ENDDO
      ENDIF
      !
      fac = 1.d0 / DBLE(nspin_lsda)
      DO is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core
      ENDDO
      !
      IF (xclib_dft_is('gradient')) THEN
        CALL dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq0, drhoc, &
                       1, nspin_gga, g, aux)
      ENDIF
      !
      IF (dft_is_nonlocc()) THEN
        CALL dnonloccorr(rho%of_r, drhoc, xq0, aux)
      ENDIF
      !
      DO is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core
      ENDDO
      !
      CALL fwfft('Rho', aux, dfftp)
      !
      ! This is needed also when the smooth and the thick grids coincide to
      ! cut the potential at the cut-off
      !
      auxs(:) = czero
      DO ig = 1, ngms
        auxs(dffts%nl(ig)) = aux(dfftp%nl(ig))
      ENDDO
      aux1(:) = aux1(:) + auxs(:)
    ENDIF
    !
    ! Now we compute dV_loc/dtau in real space
    !
    CALL invfft('Rho', aux1, dffts)
    DO ibnd = lower_band, upper_band
      DO ip = 1, npol
        aux2(:) = czero
        IF (ip == 1) THEN
          DO ig = 1, npw
            aux2(dffts%nl(igk(ig))) = evc(ig, ibnd + ibndstart - 1)
          ENDDO
        ELSE
          DO ig = 1, npw
            aux2(dffts%nl(igk(ig))) = evc(ig + npwx, ibnd + ibndstart - 1)
          ENDDO
        ENDIF
        !
        !  This wavefunction is computed in real space
        !
        CALL invfft('Wave', aux2, dffts)
        DO ir = 1, dffts%nnr
          aux2(ir) = aux2(ir) * aux1(ir)
        ENDDO
        !
        ! and finally dV_loc/dtau * psi is transformed in reciprocal space
        !
        CALL fwfft('Wave', aux2, dffts)
        IF (ip == 1) THEN
          DO ig = 1, npwq
            dvpsi(ig, ibnd) = aux2(dffts%nl(igkq(ig)))
          ENDDO
        ELSE
          DO ig = 1, npwq
            dvpsi(ig + npwx, ibnd) = aux2(dffts%nl(igkq(ig)))
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    IF (nlcc_any .AND. addnlcc) THEN
      DEALLOCATE(drhoc, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error deallocating drhoc', 1)
      DEALLOCATE(aux, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error deallocating aux', 1)
      DEALLOCATE(auxs, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error deallocating auxs', 1)
    ENDIF
    DEALLOCATE(aux1, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error deallocating aux1', 1)
    DEALLOCATE(aux2, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvqpsi_us3', 'Error deallocating aux2', 1)
    !
    ! We add the contribution of the nonlocal potential in the US form
    ! First a term similar to the KB case.
    ! Then a term due to the change of the D coefficients in the perturbat
    !
    CALL dvqpsi_us_only3(ik, uact, xxkq, igkq, npwq)
    !
    CALL stop_clock('dvqpsi_us3')
    !
    RETURN
    !
    !----------------------------------------------------------------------
    END SUBROUTINE dvqpsi_us3
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE dvqpsi_us_only3(ik, uact, xxkq, igkq, npwq)
    !----------------------------------------------------------------------
    !!
    !! This routine calculates dV_bare/dtau * psi for one perturbation
    !! with a given q. The displacements are described by a vector uact.
    !! The result is stored in dvpsi. The routine is called for each k point
    !! and for each pattern u. It computes simultaneously all the bands.
    !! This routine implements Eq. B29 of PRB 64, 235118 (2001).
    !! Only the contribution of the nonlocal potential is calculated here.
    !! Adapted from PH/dvqpsi_us_only (QE)
    !!
    !-----------------------------------------------------------------------
    USE kinds,      ONLY : DP
    USE cell_base,  ONLY : tpiba
    USE gvect,      ONLY : g
    USE ions_base,  ONLY : nat, ityp, ntyp => nsp
    USE lsda_mod,   ONLY : lsda, current_spin, nspin
    USE wvfct,      ONLY : npwx, et
    USE uspp,       ONLY : okvan, nkb, vkb
    USE uspp_param, ONLY : nh, nhm
    USE phus,       ONLY : int1, int1_nc, int2, int2_so, alphap
    USE lrus,       ONLY : becp1
    USE eqv,        ONLY : dvpsi
    USE elph2,      ONLY : lower_band, upper_band, ibndstart
    USE noncollin_module, ONLY : noncolin, npol, lspinorb
    USE constants_epw,    ONLY : czero, zero, cone, eps12
    USE klist_epw,  ONLY : isk_loc
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! the k point
    INTEGER, INTENT(in) :: npwq
    !! Number of k+G-vectors inside 'ecut sphere'
    INTEGER, INTENT(in) :: igkq(npwq)
    !! k+G+q mapping
    REAL(KIND = DP), INTENT(in) :: xxkq(3)
    !! the k+q point (cartesian coordinates)
    COMPLEX(KIND = DP), INTENT(in) :: uact(3 * nat)
    !! the pattern of displacements
    !
    ! Local variables
    LOGICAL :: ok
    !!
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nb
    !! Counter on atoms
    INTEGER :: mu
    !! Counter on modes
    INTEGER :: nu
    !! Counter on modes
    INTEGER :: ig
    !! Counter on G vectors
    INTEGER :: igg
    !! Auxiliary counter on G vectors
    INTEGER :: nt
    !! Counter on atomic types
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ijkb0
    !! Auxiliary variable for counting
    INTEGER :: ikb
    !! Counter on becp functions
    INTEGER :: jkb
    !! Counter on becp functions
    INTEGER :: ipol
    !! Counter on polarizations
    INTEGER :: ih
    !! Counter on beta functions
    INTEGER :: jh
    !! Counter on beta functions
    INTEGER :: is
    !! Counter on polarization
    INTEGER :: js
    !! Counter on polarization
    INTEGER ::  ijs
    !! Counter on combined is and js polarization
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: deff(:, :, :)
    !
    COMPLEX(KIND = DP), ALLOCATABLE :: ps1(:, :)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: ps2(:, :, :)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: aux(:)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: deff_nc(:, :, :, :)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: ps1_nc(:, :, :)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE :: ps2_nc(:, :, :, :)
    !!
    !
    CALL start_clock('dvqpsi_us_on')
    IF (noncolin) THEN
      ALLOCATE(ps1_nc(nkb, npol, lower_band:upper_band), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating ps1_nc', 1)
      ALLOCATE(ps2_nc(nkb, npol, lower_band:upper_band, 3), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating ps2_nc', 1)
      ALLOCATE(deff_nc(nhm, nhm, nat, nspin), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating deff_nc', 1)
      ps1_nc(:, :, :) = czero
      ps2_nc(:, :, :, :) = czero
      deff_nc(:, :, :, :) = czero
    ELSE
      ALLOCATE(ps1(nkb, lower_band:upper_band), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating ps1', 1)
      ALLOCATE(ps2(nkb, lower_band:upper_band, 3), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating ps2', 1)
      ALLOCATE(deff(nhm, nhm, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating deff', 1)
      ps1(:, :) = czero
      ps2(:, :, :) = czero
      deff(:, :, :) = zero
    ENDIF
    ALLOCATE(aux(npwx), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error allocating aux', 1)
    aux(:) = czero
    !
    IF (lsda) current_spin = isk_loc(ik)
    !
    !   we first compute the coefficients of the vectors
    !
    DO ibnd = lower_band, upper_band
      IF (noncolin) THEN
        CALL compute_deff_nc(deff_nc, et(ibnd + ibndstart - 1, ik))
      ELSE
        CALL compute_deff(deff, et(ibnd + ibndstart - 1, ik))
      ENDIF
      !
      ijkb0 = 0
      DO nt = 1, ntyp
        DO na = 1, nat
          IF (ityp(na) == nt) THEN
            mu = 3 * (na - 1)
            DO ih = 1, nh(nt)
              ikb = ijkb0 + ih
              DO jh = 1, nh(nt)
                jkb = ijkb0 + jh
                DO ipol = 1, 3
                  IF (ABS(uact(mu + 1)) + ABS(uact(mu + 2)) + ABS(uact(mu + 3)) > eps12) THEN
                    IF (noncolin) THEN
                      ijs = 0
                      DO is = 1, npol
                        DO js = 1, npol
                          ijs = ijs + 1
                          ps1_nc(ikb, is, ibnd) = ps1_nc(ikb, is, ibnd) + deff_nc(ih, jh, na, ijs) * &
                                 alphap(ipol, ik)%nc(jkb, js, ibnd + ibndstart - 1) * uact(mu + ipol)
                          ps2_nc(ikb, is, ibnd, ipol) = ps2_nc(ikb, is, ibnd, ipol) + &
                                 deff_nc(ih, jh, na, ijs) * becp1(ik)%nc(jkb, js, ibnd + ibndstart - 1) * &
                                 (0.d0, -1.d0) * uact(mu + ipol) * tpiba
                        ENDDO
                      ENDDO
                    ELSE
                      ps1(ikb, ibnd) = ps1(ikb, ibnd) + &
                                       deff(ih, jh, na) * alphap(ipol, ik)%k(jkb, ibnd + ibndstart - 1) * uact(mu + ipol)
                      ps2(ikb, ibnd, ipol) = ps2(ikb, ibnd, ipol) + deff(ih, jh, na) * becp1(ik)%k(jkb, ibnd + ibndstart - 1) * &
                                             (0.d0, -1.d0) * uact(mu + ipol) * tpiba
                    ENDIF
                    IF (okvan) THEN
                      IF (noncolin) THEN
                        ijs = 0
                        DO is = 1, npol
                          DO js = 1, npol
                            ijs = ijs + 1
                            ps1_nc(ikb, is, ibnd) = ps1_nc(ikb, is, ibnd) + int1_nc(ih, jh, ipol, na, ijs) * &
                               becp1(ik)%nc(jkb, js, ibnd + ibndstart - 1) * uact(mu + ipol)
                          ENDDO
                        ENDDO
                      ELSE
                        ps1(ikb, ibnd) = ps1(ikb, ibnd) + int1(ih, jh, ipol, na, current_spin) * &
                            becp1(ik)%k(jkb, ibnd + ibndstart - 1) * uact(mu + ipol)
                      ENDIF
                    ENDIF ! okvan
                  ENDIF  ! uact>0
                  IF (okvan) THEN
                    DO nb = 1, nat
                      nu = 3 * (nb - 1)
                      IF (noncolin) THEN
                        IF (lspinorb) THEN
                          ijs = 0
                          DO is = 1, npol
                            DO js = 1, npol
                              ijs = ijs + 1
                              ps1_nc(ikb, is, ibnd) = ps1_nc(ikb, is, ibnd) + int2_so(ih, jh, ipol, nb, na, ijs) * &
                                 becp1(ik)%nc(jkb, js, ibnd + ibndstart - 1) * uact(nu + ipol)
                            ENDDO
                          ENDDO
                        ELSE
                          DO is = 1, npol
                            ps1_nc(ikb, is, ibnd) = ps1_nc(ikb, is, ibnd) + int2(ih, jh, ipol, nb, na) * &
                               becp1(ik)%nc(jkb, is, ibnd + ibndstart - 1) * uact(nu + ipol)
                          ENDDO
                        ENDIF
                      ELSE
                        ps1(ikb,ibnd) = ps1(ikb,ibnd) + int2(ih, jh, ipol, nb, na) * &
                            becp1(ik)%k(jkb, ibnd + ibndstart - 1) * uact(nu + ipol)
                      ENDIF
                    ENDDO
                  ENDIF  ! okvan
                ENDDO ! ipol
              ENDDO ! jh
            ENDDO ! ih
            ijkb0 = ijkb0 + nh(nt)
          ENDIF
        ENDDO  ! na
      ENDDO ! nt
    ENDDO ! nbnd
    !
    !      This term is proportional to beta(k+q+G)
    !
    IF (nkb > 0) THEN
      IF (noncolin) THEN
        CALL ZGEMM('n', 'n', npwq, (upper_band - lower_band + 1)*npol, nkb, &
                   cone, vkb, npwx, ps1_nc, nkb, cone, dvpsi, npwx)
      ELSE
        CALL ZGEMM('n', 'n', npwq, (upper_band - lower_band + 1), nkb, &
                   cone, vkb, npwx, ps1, nkb, cone, dvpsi, npwx)
      ENDIF
    ENDIF
    !
    !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
    !
    DO ikb = 1, nkb
      DO ipol = 1, 3
        ok = .FALSE.
        IF (noncolin) THEN
          DO ibnd = lower_band, upper_band
            ok = ok .OR. (ABS(ps2_nc(ikb, 1, ibnd, ipol) ) > eps12) .OR. &
                         (ABS(ps2_nc(ikb, 2, ibnd, ipol) ) > eps12)
          ENDDO
        ELSE
          DO ibnd = lower_band, upper_band
            ok = ok .OR. (ABS(ps2(ikb, ibnd, ipol)) > eps12)
          ENDDO
        ENDIF
        IF (ok) THEN
          DO ig = 1, npwq
            igg = igkq(ig)
            aux(ig) = vkb(ig, ikb) * (xxkq(ipol) + g(ipol, igg))
          ENDDO
          DO ibnd = lower_band, upper_band
            IF (noncolin) THEN
              CALL ZAXPY(npwq, ps2_nc(ikb, 1, ibnd, ipol), aux, 1, dvpsi(1, ibnd), 1)
              CALL ZAXPY(npwq, ps2_nc(ikb, 2, ibnd, ipol), aux, 1, dvpsi(1 + npwx, ibnd), 1)
            ELSE
              CALL ZAXPY(npwq, ps2(ikb, ibnd, ipol), aux, 1, dvpsi(1, ibnd), 1)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    DEALLOCATE(aux, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating aux', 1)
    IF (noncolin) THEN
      DEALLOCATE(ps1_nc, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating ps1_nc', 1)
      DEALLOCATE(ps2_nc, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating ps2_nc', 1)
      DEALLOCATE(deff_nc, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating deff_nc', 1)
    ELSE
      DEALLOCATE(ps1, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating ps1', 1)
      DEALLOCATE(ps2, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating ps2', 1)
      DEALLOCATE(deff, STAT = ierr)
      IF (ierr /= 0) CALL errore('dvqpsi_us_only3', 'Error deallocating deff', 1)
    ENDIF
    !
    CALL stop_clock('dvqpsi_us_on')
    !
    RETURN
    !
    !----------------------------------------------------------------------
    END SUBROUTINE dvqpsi_us_only3
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE dvanqq2()
    !----------------------------------------------------------------------
    !!
    !! This routine calculates two integrals of the Q functions and
    !! its derivatives with V_loc and V_eff which are used
    !! to compute term dV_bare/dtau * psi  in addusdvqpsi.
    !! The result is stored in int1, int2. The routine is called
    !! for each q in nqc.
    !! int1 -> Eq. B20 of Ref.[1]
    !! int2 -> Eq. B21 of Ref.[1]
    !!
    !! [1] PRB 64, 235118 (2001).
    !!
    !! RM - Nov/Dec 2014
    !! Imported the noncolinear case implemented by xlzhang
    !!
    !! Roxana Margine - Dec 2018: Updated based on QE 6.3
    !! SP: Sept. 2019 - Cleaning
    !!
    !! HL: Mar 2020 - Parallelization over G using intra image communicator 
    !!
    !
    USE kinds,            ONLY : DP
    USE ions_base,        ONLY : nat, ityp, ntyp => nsp
    USE cell_base,        ONLY : tpiba2, omega, tpiba
    USE gvect,            ONLY : ngm, gg, g, eigts1, eigts2, eigts3, mill
    USE scf,              ONLY : v, vltot
    USE noncollin_module, ONLY : noncolin, nspin_mag
    USE phcom,            ONLY : int1, int2, vlocq
    USE qpoint,           ONLY : xq, eigqts
    USE uspp_param,       ONLY : upf, lmaxq, nh
    USE uspp,             ONLY : okvan, ijtoh
    USE mp_global,        ONLY : intra_pool_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_interfaces,   ONLY : fwfft
    USE constants_epw,    ONLY : zero, czero
    USE mp_images,        ONLY : intra_image_comm
    USE elph2,            ONLY : veff, ig_s, ig_e
    !
    IMPLICIT NONE
    !
    ! Local variables
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
    INTEGER :: ngvec
    !! number of G in this group
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
    COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:), aux2(:), aux5(:)
    !! Auxiallary array
    COMPLEX(KIND = DP), ALLOCATABLE :: sk(:)
    !!
    COMPLEX(KIND = DP), ALLOCATABLE, TARGET :: qgm(:)
    !! the augmentation function at G
    COMPLEX(KIND = DP), POINTER :: qgmq(:)
    !! the augmentation function at q+G
    !
    IF (.NOT. okvan) RETURN
    !
    CALL start_clock('dvanqq2')
    !
    ngvec = ig_e - ig_s + 1
    !
    int1(:, :, :, :, :) = czero
    int2(:, :, :, :, :) = czero
    !
    ALLOCATE(sk(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating sk', 1)
    ALLOCATE(aux1(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux1', 1)
    ALLOCATE(aux2(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux2', 1)
    ALLOCATE(aux5(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating aux5', 1)
    ALLOCATE(qmodg(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qmodg', 1)
    ALLOCATE(qmod(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qmod', 1)
    ALLOCATE(qgmq(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qgmq', 1)
    ALLOCATE(qgm(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qgm', 1)
    ALLOCATE(ylmk0(ngvec, lmaxq * lmaxq), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating ylmk0', 1)
    ALLOCATE(ylmkq(ngvec, lmaxq * lmaxq), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating ylmkq', 1)
    sk(:) = czero
    aux1(:) = czero
    aux2(:) = czero
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
    CALL ylmr2(lmaxq * lmaxq, ngvec, g(:, ig_s), gg(ig_s), ylmk0)
    !
    DO ig = 1, ngvec
      qmodg(ig) = DSQRT(gg(ig + ig_s - 1))*tpiba
    ENDDO
    !
    ALLOCATE(qpg(3, ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error allocating qpg', 1)
    qpg(:, :) = zero
    !
    CALL setqmod(ngvec, xq, g(:, ig_s), qmod, qpg)
    CALL ylmr2(lmaxq * lmaxq, ngvec, qpg, qmod, ylmkq)
    !
    DEALLOCATE(qpg, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating qpg', 1)
    DO ig = 1, ngvec
      qmod(ig) = DSQRT(qmod(ig))*tpiba
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
            CALL qvan2(ngvec, ih, jh, ntb, qmodg, qgm, ylmk0)
            CALL qvan2(ngvec, ih, jh, ntb, qmod, qgmq, ylmkq)
            !
            ! NB: for this integral the moving atom and the atom of Q
            ! do not necessarily coincide
            !
            DO nb = 1, nat
              IF (ityp(nb) == ntb) THEN
                DO ig = 1, ngvec
                  aux1(ig) = qgmq(ig) * eigts1(mill(1, ig + ig_s - 1), nb) &
                                      * eigts2(mill(2, ig + ig_s - 1), nb) &
                                      * eigts3(mill(3, ig + ig_s - 1), nb)
                ENDDO
                !
                DO na = 1, nat
                  fact = eigqts(na) * CONJG(eigqts(nb))
                  !
                  !    nb is the atom of the augmentation function
                  !
                  nta = ityp(na)
                  DO ig = 1, ngvec
                    sk(ig) = vlocq(ig + ig_s - 1, nta) &
                             * eigts1(mill(1, ig + ig_s - 1), na) &
                             * eigts2(mill(2, ig + ig_s - 1), na) &
                             * eigts3(mill(3, ig + ig_s - 1), na)
                  ENDDO
                  !
                  DO ipol = 1, 3
                    DO ig = 1, ngvec
                      aux5(ig) = sk(ig) * (g(ipol, ig + ig_s - 1) + xq(ipol))
                    ENDDO
                    int2(ih, jh, ipol, na, nb) = fact * fact1 * ZDOTC(ngvec, aux1, 1, aux5, 1)
                    !
                  ENDDO !ipol
                  !
                ENDDO !na
                !
                DO ig = 1, ngvec
                  aux1(ig) = qgm(ig) * eigts1(mill(1,ig + ig_s - 1),nb) &
                                     * eigts2(mill(2,ig + ig_s - 1),nb) &
                                     * eigts3(mill(3,ig + ig_s - 1),nb)
                ENDDO
                !
                DO is = 1, nspin_mag
                  DO ipol = 1, 3
                    DO ig = 1, ngvec
                      aux2(ig) = veff(dfftp%nl(ig + ig_s - 1), is) * g(ipol, ig + ig_s - 1)
                    ENDDO
                    int1(ih, jh, ipol, nb, is) = - fact1 * ZDOTC(ngvec, aux1, 1, aux2, 1)
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
    CALL mp_sum(int1, intra_image_comm)
    CALL mp_sum(int2, intra_image_comm)
    !
    IF (noncolin) THEN
      CALL set_int12_nc(0)
    ENDIF
    !
    DEALLOCATE(sk, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating sk', 1)
    DEALLOCATE(aux1, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating aux1', 1)
    DEALLOCATE(aux2, STAT = ierr)
    IF (ierr /= 0) CALL errore('dvanqq2', 'Error deallocating aux2', 1)
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
    !
    CALL stop_clock('dvanqq2')
    RETURN
    !
    !----------------------------------------------------------------------
    END SUBROUTINE dvanqq2
    !----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE newdq2(dvscf, npe, xq0, timerev)
    !----------------------------------------------------------------------
    !!
    !! This routine computes the contribution of the selfconsistent
    !! change of the potential to the known part of the linear
    !! system and adds it to dvpsi.
    !! Adapted from LR_Modules/newdq.f90 (QE)
    !!
    !! Roxana Margine - Jan 2019: Updated based on QE 6.3
    !!
    !! HL: Mar 2020
    !! 1. Parallelization over G using intra image communicator 
    !! 2. PAW added and cleaning
    !!
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ityp, ntyp => nsp
    USE noncollin_module,     ONLY : noncolin, nspin_mag
    USE cell_base,            ONLY : omega, tpiba
    USE fft_base,             ONLY : dfftp
    USE fft_interfaces,       ONLY : fwfft
    USE gvect,                ONLY : g, ngm, mill, eigts1, eigts2, eigts3
    USE uspp,                 ONLY : okvan
    USE uspp_param,           ONLY : upf, lmaxq, nh
    USE paw_variables,        ONLY : okpaw
    USE mp_global,            ONLY : intra_pool_comm
    USE mp,                   ONLY : mp_sum
    USE lrus,                 ONLY : int3, int3_paw
    USE qpoint,               ONLY : eigqts
    USE constants_epw,        ONLY : czero, zero
    USE mp_images,            ONLY : intra_image_comm
    USE elph2,                ONLY : ig_s, ig_e
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: timerev
    !!  true if we are using time reversal
    INTEGER, INTENT(in) :: npe
    !! Number of perturbations for this irr representation
    REAL(KIND = DP), INTENT(in) :: xq0(3)
    !! The first q-point in the star (cartesian coords.)
    COMPLEX(KIND = DP), INTENT(in) :: dvscf(dfftp%nnr, nspin_mag, npe)
    !! Change of the selfconsistent potential
    !
    ! Local variables
    INTEGER :: na
    !! counter on atoms
    INTEGER :: ig
    !! counter on G vectors
    INTEGER :: nt
    !! counter on atomic types
    INTEGER ::  ir
    !! counter on real mesh
    INTEGER :: ipert
    !! counter on change of Vscf due to perturbations
    INTEGER :: is
    !! counter on spin
    INTEGER :: ih
    !! Counter on beta functions
    INTEGER :: jh
    !! Counter on beta functions
    INTEGER :: ierr
    !! Error status
    INTEGER :: ngvec
    !! number of G in this group
    REAL(KIND = DP), ALLOCATABLE :: qmod(:)
    !! the modulus of q+G
    REAL(KIND = DP), ALLOCATABLE :: qg(:, :)
    !! the values of q+G
    REAL(KIND = DP), ALLOCATABLE :: ylmk0(:, :)
    !! the spherical harmonics at q+G
    COMPLEX(KIND = DP), EXTERNAL :: ZDOTC
    !! the scalar product function
    COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:), aux2(:, :)
    !! Auxillary variable
    COMPLEX(KIND = DP), ALLOCATABLE :: qgm(:)
    !! the augmentation function at q+G
    COMPLEX(KIND = DP), ALLOCATABLE :: veff(:)
    !! effective potential
    !
    IF (.NOT. okvan) RETURN
    !
    CALL start_clock('newdq2')
    !
    ngvec = ig_e - ig_s + 1
    !
    int3(:, :, :, :, :) = czero
    ALLOCATE(aux1(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating aux1', 1)
    ALLOCATE(aux2(ngvec, nspin_mag), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating aux2', 1)
    ALLOCATE(veff(dfftp%nnr), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating veff', 1)
    ALLOCATE(ylmk0(ngvec, lmaxq * lmaxq), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating ylmk0', 1)
    ALLOCATE(qgm(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating qgm', 1)
    ALLOCATE(qmod(ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating qmod', 1)
    ALLOCATE(qg(3, ngvec), STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error allocating qg', 1)
    aux1(:)     = czero
    aux2(:, :)  = czero
    veff(:)     = czero
    ylmk0(:, :) = zero
    qgm(:)      = czero
    qmod(:)     = zero
    qg(:, :)    = zero
    !
    ! first compute the spherical harmonics
    !
    CALL setqmod(ngvec, xq0, g(:, ig_s), qmod, qg)
    CALL ylmr2(lmaxq * lmaxq, ngvec, qg, qmod, ylmk0)
    !
    DO ig = 1, ngvec
      qmod(ig) = DSQRT(qmod(ig))*tpiba
    ENDDO
    !
    ! and for each perturbation of this irreducible representation
    ! integrate the change of the self consistent potential and
    ! the Q functions
    !
    DO ipert = 1, npe
      DO is = 1, nspin_mag
        DO ir = 1, dfftp%nnr
          IF (timerev) THEN
            veff(ir) = CONJG(dvscf(ir, is, ipert))
          ELSE
            veff(ir) = dvscf(ir, is, ipert)
          ENDIF
        ENDDO
        CALL fwfft('Rho', veff, dfftp)
        DO ig = 1, ngvec
          aux2(ig, is) = veff(dfftp%nl(ig + ig_s - 1))
        ENDDO
      ENDDO
      !
      DO nt = 1, ntyp
        IF (upf(nt)%tvanp) THEN
          !
          DO ih = 1, nh(nt)
            DO jh = ih, nh(nt)
              !
              CALL qvan2(ngvec, ih, jh, nt, qmod, qgm, ylmk0)
              !
              DO na = 1, nat
                IF (ityp(na) == nt) THEN
                  DO ig = 1, ngvec
                    aux1(ig) = qgm(ig) * eigts1(mill(1, ig + ig_s - 1), na) * &
                                         eigts2(mill(2, ig + ig_s - 1), na) * &
                                         eigts3(mill(3, ig + ig_s - 1), na) * &
                                         eigqts(na)
                  ENDDO
                  DO is = 1, nspin_mag
                    int3(ih, jh, na, is, ipert) = omega * ZDOTC(ngvec, aux1, 1, aux2(1, is), 1)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO ! jh
          ENDDO ! ih
          !
          DO na = 1, nat
            IF (ityp(na) == nt) THEN
              !
              ! We use the symmetry properties of the ps factor
              !
              DO ih = 1, nh(nt)
                DO jh = ih, nh(nt)
                  DO is = 1, nspin_mag
                    int3(jh, ih, na, is, ipert) = int3(ih, jh, na, is, ipert)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF ! ityp
          ENDDO ! na
        ENDIF ! upf
      ENDDO ! nt
    ENDDO ! ipert
    !
    CALL mp_sum(int3, intra_image_comm)
    !
    ! Sum of the USPP and PAW terms 
    ! (see last two terms in Eq.(12) in PRB 81, 075123 (2010))
    !
    IF (okpaw) int3 = int3 + int3_paw
    !
    IF (noncolin) CALL set_int3_nc(npe)
    !
    DEALLOCATE(aux1, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating aux1', 1)
    DEALLOCATE(aux2, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating aux2', 1)
    DEALLOCATE(veff, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating veff', 1)
    DEALLOCATE(ylmk0, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating ylmk0', 1)
    DEALLOCATE(qgm, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating qgm', 1)
    DEALLOCATE(qmod, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating qmod', 1)
    DEALLOCATE(qg, STAT = ierr)
    IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating qg', 1)
    !
    CALL stop_clock('newdq2')
    !
    RETURN
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE newdq2
    !---------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE adddvscf2(ipert, ik)
    !----------------------------------------------------------------------
    !!
    !! This routine computes the contribution of the selfconsistent
    !! change of the potential to the known part of the linear
    !! system and adds it to dvpsi.
    !! It implements the second term in Eq. B30 of PRB 64, 235118 (2001).
    !! Only used in the case of USPP.
    !! Adapted from LR_Modules/adddvscf.f90 (QE)
    !!
    !! Roxana Margine - Jan 2019: Updated based on QE 6.3
    !! SP - Jan 2019: Clean
    !!
    USE kinds,      ONLY : DP
    USE uspp_param, ONLY : upf, nh
    USE uspp,       ONLY : vkb, okvan
    USE lsda_mod,   ONLY : lsda, current_spin
    USE klist_epw,  ONLY : isk_loc
    USE ions_base,  ONLY : ntyp => nsp, nat, ityp
    USE wvfct,      ONLY : npwx
    USE lrus,       ONLY : int3, int3_nc, becp1
    USE qpoint,     ONLY : npwq
    USE eqv,        ONLY : dvpsi
    USE elph2,      ONLY : lower_band, upper_band, ibndstart
    USE noncollin_module, ONLY : noncolin, npol
    USE constants_epw,    ONLY : czero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-point
    INTEGER, INTENT(in) :: ipert
    !! Counter on Vscf perturbations
    !
    !   Local variables
    !
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: nt
    !! Counter on atomic types
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ih
    !! Counter on beta functions
    INTEGER :: jh
    !! Counter on beta functions
    INTEGER :: ijkb0
    !! Auxiliary variable for counting
    INTEGER :: ikb
    !! Counter on becp functions
    INTEGER :: jkb
    !! Counter on becp functions
    INTEGER :: is
    !! Counter on polarization
    INTEGER :: js
    !! Counter on polarization
    INTEGER ::  ijs
    !! Counter on combined is and js polarization
    !
    COMPLEX(KIND = DP) :: sum_k
    !! auxiliary sum variable
    COMPLEX(KIND = DP) :: sum_nc(npol)
    !! auxiliary sum variable non-collinear case
    !
    IF (.NOT. okvan) RETURN
    !
    CALL start_clock('adddvscf2')
    !
    IF (lsda) current_spin = isk_loc(ik)
    !
    ijkb0 = 0
    DO nt = 1, ntyp
      IF (upf(nt)%tvanp) THEN
        DO na = 1, nat
          IF (ityp(na) == nt) THEN
            !
            ! We multiply the integral for the becp term and the beta_n
            !
            DO ibnd = lower_band, upper_band
              DO ih = 1, nh(nt)
                 ikb = ijkb0 + ih
                 IF (noncolin) THEN
                   sum_nc = czero
                 ELSE
                   sum_k = czero
                 ENDIF
                 DO jh = 1, nh(nt)
                   jkb = ijkb0 + jh
                   IF (noncolin) THEN
                     ijs = 0
                     DO is = 1, npol
                       DO js = 1, npol
                         ijs = ijs + 1
                         sum_nc(is) = sum_nc(is) + int3_nc(ih, jh, na, ijs, ipert) * becp1(ik)%nc(jkb, js, ibnd + ibndstart - 1)
                       ENDDO
                     ENDDO
                   ELSE
                     sum_k = sum_k + int3(ih,jh,na,current_spin,ipert) * &
                                 becp1(ik)%k(jkb,ibnd + ibndstart - 1)
                   ENDIF
                 ENDDO
                 IF (noncolin) THEN
                   CALL ZAXPY(npwq, sum_nc(1), vkb(1, ikb), 1, dvpsi(1, ibnd), 1)
                   CALL ZAXPY(npwq, sum_nc(2), vkb(1, ikb), 1, dvpsi(1 + npwx, ibnd), 1)
                 ELSE
                   CALL ZAXPY(npwq, sum_k, vkb(1, ikb), 1, dvpsi(1, ibnd), 1)
                 ENDIF
              ENDDO
            ENDDO
            ijkb0 = ijkb0 + nh(nt)
          ENDIF
        ENDDO
      ELSE
        DO na = 1, nat
          IF (ityp(na) == nt) ijkb0 = ijkb0 + nh(nt)
        ENDDO
      ENDIF
    ENDDO
    !
    CALL stop_clock('adddvscf2')
    !
    RETURN
    !
    !---------------------------------------------------------------------------
    END SUBROUTINE adddvscf2
    !---------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  END MODULE dvqpsi
  !-----------------------------------------------------------------------------

