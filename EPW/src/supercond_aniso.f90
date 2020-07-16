  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE supercond_aniso
  !----------------------------------------------------------------------
  !!
  !! This module contains all the subroutines linked with superconductivity using
  !! the isotropic or anisotropic Eliashberg formalism.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the anisotropic
    !! Eliashberg equations on the imaginary-axis.
    !!
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE control_flags,     ONLY : iverbosity
    USE epwcom,            ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                                  limag, lpade, lacon, fsthick, imag_read, wscut
    USE elph2,             ONLY : gtemp
    USE eliashbergcom,     ONLY : nsw, nsiw, adelta, adeltap, adeltai, adeltaip, &
                                  nkfs, nbndfs, ekfs, ef0
    USE supercond,         ONLY : free_energy, dos_quasiparticle, gen_freqgrid_iaxis, &
                                  deallocate_eliashberg_iaxis, deallocate_eliashberg_raxis, &
                                  deallocate_eliashberg_aniso, eliashberg_grid
    USE constants_epw,     ONLY : kelvin2eV, ci, zero
    USE io_global,         ONLY : ionode_id
    USE mp_global,         ONLY : inter_pool_comm
    USE mp,                ONLY : mp_bcast, mp_barrier
    USE mp_world,          ONLY : mpime
    USE io_eliashberg,     ONLY : eliashberg_read_aniso_iaxis
    USE utilities,         ONLY : mix_broyden
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE printing,          ONLY : prtheader_supercond
    !
    IMPLICIT NONE
    !
    ! Local variables
    LOGICAL :: conv
    !! True if calculation is converged
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iter
    !! Counter on iteration steps
    INTEGER :: N
    !! Maximum nr. frequency points in Pade approx
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: tcpu
    !! cpu time
    REAL(KIND = DP), EXTERNAL :: get_clock
    !! get the time spent
    REAL(KIND = DP), ALLOCATABLE :: rdeltain(:), rdeltaout(:)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(KIND = DP), ALLOCATABLE :: cdeltain(:), cdeltaout(:)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(KIND = DP), ALLOCATABLE :: df1(:, :, :, :), df2(:, :, :, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv1(:, :, :, :), dv2(:, :, :, :)
    !! Temporary variables for mix_broyden
    !
    CALL start_clock('aniso_iaxis')
    !
    CALL eliashberg_grid()
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      CALL prtheader_supercond(itemp, 1)
      CALL start_clock('iaxis_imag')
      CALL gen_freqgrid_iaxis(itemp)
      !
      IF ((limag .AND. .NOT. imag_read) .OR. (limag .AND. imag_read .AND. itemp /= 1)) THEN
        !
        IF (mpime == ionode_id) THEN
          ALLOCATE(df1(nbndfs, nkfs, nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df1', 1)
          ALLOCATE(dv1(nbndfs, nkfs, nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv1', 1)
          df1(:, :, :, :) = zero
          dv1(:, :, :, :) = zero
        ENDIF
        !
        iter = 1
        conv = .FALSE.
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL sum_eliashberg_aniso_iaxis(itemp, iter, conv)
          IF (mpime == ionode_id) THEN
            DO ik = 1, nkfs
              DO ibnd = 1, nbndfs
                IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                  CALL mix_broyden(nsiw(itemp), adeltai(ibnd, ik, :), adeltaip(ibnd, ik, :), &
                                   broyden_beta, iter, broyden_ndim, conv, &
                                   df1(ibnd, ik, :, :), dv1(ibnd, ik, :, :))
                ENDIF
              ENDDO
            ENDDO
          ENDIF
          CALL mp_bcast(adeltai,  ionode_id, inter_pool_comm)
          CALL mp_bcast(adeltaip, ionode_id, inter_pool_comm)
          CALL mp_barrier(inter_pool_comm)
          iter = iter + 1
        ENDDO ! iter
        !
        IF (mpime == ionode_id) THEN
          DEALLOCATE(df1, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df1', 1)
          DEALLOCATE(dv1, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv1', 1)
        ENDIF
        !
        IF (conv) THEN
          IF (iverbosity == 2) THEN
            IF (mpime == ionode_id) THEN
              CALL free_energy(itemp)
            ENDIF
            CALL mp_barrier(inter_pool_comm)
          ENDIF
          CALL stop_clock('iaxis_imag')
          CALL print_clock('iaxis_imag')
          WRITE(stdout,'(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          CALL deallocate_eliashberg_iaxis()
          CALL deallocate_eliashberg_aniso()
          CALL stop_clock('iaxis_imag')
          CALL print_clock('iaxis_imag')
          WRITE(stdout,'(a)') ' '
          RETURN
        ENDIF
      ELSEIF (limag .AND. imag_read .AND. itemp == 1) THEN
        CALL eliashberg_read_aniso_iaxis(itemp)
      ENDIF ! limag
      !
      IF (lpade) THEN
        CALL prtheader_supercond(itemp, 2)
        CALL start_clock('raxis_pade')
        N = 90 * nsiw(itemp) / 100
        IF (mod(N,2) /= 0 ) N = N + 1
        CALL pade_cont_aniso_iaxis_to_raxis(itemp, N)
        !
        IF (mpime == ionode_id) THEN
          CALL dos_quasiparticle(itemp)
        ENDIF
        CALL mp_barrier(inter_pool_comm)
        CALL stop_clock('raxis_pade')
        CALL print_clock('raxis_pade')
        WRITE(stdout,'(a)') ' '
      ENDIF ! lpade
      !
      IF (lacon) THEN
        CALL prtheader_supercond(itemp, 3)
        CALL start_clock('raxis_acon')
        !
        iter = 1
        conv = .FALSE.
        IF (mpime == ionode_id) THEN
          ALLOCATE(rdeltain(nsw), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating rdeltain', 1)
          ALLOCATE(cdeltain(nsw), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating cdeltain', 1)
          ALLOCATE(rdeltaout(nsw), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating rdeltaout', 1)
          ALLOCATE(cdeltaout(nsw), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating cdeltaout', 1)
          ALLOCATE(df1(nbndfs, nkfs, nsw, broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df1', 1)
          ALLOCATE(dv1(nbndfs, nkfs, nsw, broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv1', 1)
          ALLOCATE(df2(nbndfs, nkfs, nsw, broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df2', 1)
          ALLOCATE(dv2(nbndfs, nkfs, nsw, broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv2', 1)
          rdeltain(:)  = zero
          cdeltain(:)  = zero
          rdeltaout(:) = zero
          cdeltaout(:) = zero
          df1(:, :, :, :) = zero
          dv1(:, :, :, :) = zero
          df2(:, :, :, :) = zero
          dv2(:, :, :, :) = zero
        ENDIF
        !
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL analytic_cont_aniso_iaxis_to_raxis(itemp, iter, conv)
          IF (mpime == ionode_id) THEN
            DO ik = 1, nkfs
              DO ibnd = 1, nbndfs
                IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                  rdeltain(:)  = REAL(adeltap(ibnd, ik, :))
                  cdeltain(:)  = AIMAG(adeltap(ibnd, ik, :))
                  rdeltaout(:) = REAL(adelta(ibnd, ik, :))
                  cdeltaout(:) = AIMAG(adelta(ibnd, ik, :))
                  CALL mix_broyden(nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, &
                                   conv, df1(ibnd, ik, :, :), dv1(ibnd, ik, :, :))
                  CALL mix_broyden(nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, &
                                   conv, df2(ibnd, ik, :, :), dv2(ibnd, ik, :, :))
                  adeltap(ibnd, ik, :) = rdeltain(:) + ci * cdeltain(:)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
          CALL mp_bcast(adelta,  ionode_id, inter_pool_comm)
          CALL mp_bcast(adeltap, ionode_id, inter_pool_comm)
          CALL mp_barrier(inter_pool_comm)
          iter = iter + 1
        ENDDO ! iter
        !
        IF (mpime == ionode_id) THEN
          DEALLOCATE(rdeltain, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating rdeltain', 1)
          DEALLOCATE(cdeltain, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating cdeltain', 1)
          DEALLOCATE(rdeltaout, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating rdeltaout', 1)
          DEALLOCATE(cdeltaout, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating cdeltaout', 1)
          DEALLOCATE(df1, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df1', 1)
          DEALLOCATE(dv1, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv1', 1)
          DEALLOCATE(df2, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df2', 1)
          DEALLOCATE(dv2, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv2', 1)
        ENDIF
        !
        IF (conv) THEN
          IF (mpime == ionode_id) THEN
            CALL dos_quasiparticle(itemp)
          ENDIF
          CALL mp_barrier(inter_pool_comm)
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout,'(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          CALL deallocate_eliashberg_iaxis()
          CALL deallocate_eliashberg_raxis()
          CALL deallocate_eliashberg_aniso()
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout,'(a)') ' '
          RETURN
        ENDIF
        !
      ENDIF ! lacon
      !
      CALL deallocate_eliashberg_iaxis()
      IF (lpade .OR. lacon) CALL deallocate_eliashberg_raxis()
      !
      ! remove memory allocated for wsi, deltai, znormi, nznormi, adeltai, aznormi, naznormi
      imelt = (4 + 3 * nbndfs * nkfs) * nsiw(itemp)
      CALL mem_size_eliashberg(2, -imelt)
      !
      IF (lpade) THEN
        ! remove memory allocated for ws, delta, znorm, adelta, aznorm
        imelt = nsw + 2 * (2 + 2 * nbndfs * nkfs) * nsw
        CALL mem_size_eliashberg(2, -imelt)
      ELSEIF (lacon) THEN
        ! remove memory allocated for ws, delta, znorm, adelta, adeltap, aznorm, aznormp
        imelt = nsw + 2 * (2 + 4 * nbndfs * nkfs) * nsw
        CALL mem_size_eliashberg(2, -imelt)
      ENDIF
      !
      tcpu = get_clock('aniso_iaxis')
      WRITE(stdout, '(5x, a, i3, a, f18.2, a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
      WRITE(stdout,'(a)') ' '
      !
    ENDDO ! itemp
    !
    CALL deallocate_eliashberg_aniso()
    !
    CALL stop_clock('aniso_iaxis')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic Eliashberg equations on the imaginary-axis
    !!
    USE kinds,             ONLY : DP
    USE elph2,             ONLY : wqf, gtemp
    USE epwcom,            ONLY : nsiter, nstemp, muc, conv_thr_iaxis, fsthick
    USE eliashbergcom,     ONLY : nsiw, gap0, gap, agap, wsi, akeri, limag_fly, &
                                  naznormi, aznormi, adeltai, adeltaip, nznormi, znormi, &
                                  deltai, wsphmax, nkfs, nbndfs, dosef, ef0, ixkqf, ixqfs, &
                                  nqfs, wkfs, w0g, ekfs
    USE constants_epw,     ONLY : zero, czero
    USE constants,         ONLY : pi
    USE io_global,         ONLY : stdout, ionode_id
    USE mp_global,         ONLY : inter_pool_comm
    USE mp_world,          ONLY : mpime
    USE mp,                ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_eliashberg,     ONLY : eliashberg_write_iaxis
    USE division,          ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg, memlt_eliashberg
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    !
    ! Local variables
    INTEGER :: iw, iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: lambdam
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: lambdap
    !! K_{+}(n,n',T)
    REAL(KIND = DP) :: kernelm
    !! kernelm = lambdam - lambdap
    REAL(KIND = DP) :: kernelp
    !! kernelp = lambdam + lambdap
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP) :: esqrt
    !! Temporary variable
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP), ALLOCATABLE :: wesqrt(:, :, :), desqrt(:, :, :)
    !! Temporary variables
    REAL(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    ALLOCATE(wesqrt(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating wesqrt', 1)
    ALLOCATE(desqrt(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating desqrt', 1)
    wesqrt(:, :, :) = zero
    desqrt(:, :, :) = zero
    !
    IF (iter == 1) THEN
      !
      IF (itemp == 1) THEN
        ! get the size of required memory for  gap, agap
        imelt = (1 + nbndfs * nkfs) * nstemp
        CALL mem_size_eliashberg(2, imelt)
      ENDIF
      !
      ! get the size of required memory for
      ! wesqrt, desqrt, deltai, znormi, nznormi, adeltai, adeltaip, aznormi, naznormi, deltaold
      imelt = (4 + 6 * nbndfs * nkfs) * nsiw(itemp)
      CALL mem_size_eliashberg(2, imelt)
      !
      ALLOCATE(gap(nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating gap', 1)
      ALLOCATE(agap(nbndfs, nkfs, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating agap', 1)
      ALLOCATE(deltai(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating deltai', 1)
      ALLOCATE(znormi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating znormi', 1)
      ALLOCATE(nznormi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating nznormi', 1)
      ALLOCATE(adeltai(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating adeltai', 1)
      ALLOCATE(adeltaip(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating adeltaip', 1)
      ALLOCATE(aznormi(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating aznormi', 1)
      ALLOCATE(naznormi(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating naznormi', 1)
      gap(itemp) = zero
      agap(:,:,itemp) = zero
      adeltaip(:, :, :) = zero
      !
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            DO iw = 1, nsiw(itemp)
              IF (wsi(iw) < 2.d0 * wsphmax) THEN
                adeltaip(ibnd, ik, iw) = gap0
              ELSE
                adeltaip(ibnd, ik, iw) = zero
              ENDIF
            ENDDO
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
      CALL memlt_eliashberg(itemp, 'imag')
      IF (.NOT. limag_fly) CALL kernel_aniso_iaxis(itemp)
      !
    ENDIF
    deltai(:) = zero
    znormi(:) = zero
    nznormi(:) = zero
    adeltai(:, :, :) = zero
    aznormi(:, :, :) = zero
    naznormi(:, :, :) = zero
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) / dosef
                DO iw = 1, nsiw(itemp) ! loop over omega
                  DO iwp = 1, nsiw(itemp) ! loop over omega_prime
                    !
                    ! this step is performed at each iter step only for iw=1
                    IF (iw == 1) THEN
                      esqrt = 1.d0 / DSQRT(wsi(iwp) * wsi(iwp) + &
                                           adeltaip(jbnd, ixkqf(ik, iq0), iwp) * &
                                           adeltaip(jbnd, ixkqf(ik, iq0), iwp))
                      wesqrt(jbnd, ixkqf(ik, iq0), iwp) = wsi(iwp) * esqrt
                      desqrt(jbnd, ixkqf(ik, iq0), iwp) = adeltaip(jbnd, ixkqf(ik, iq0), iwp) * esqrt
                    ENDIF
                    IF (limag_fly) THEN
                      CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) - wsi(iwp), lambdam)
                      CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) + wsi(iwp), lambdap)
                    ELSE
                      lambdam = akeri(ik, iq, ibnd, jbnd, ABS(iw - iwp) + 1)
                      lambdap = akeri(ik, iq, ibnd, jbnd, ABS(iw + iwp))
                    ENDIF
                    ! Eq. (4.4) in Picket, PRB 26, 1186 (1982)
                    kernelm = lambdam - lambdap
                    kernelp = lambdam + lambdap
                    naznormi(ibnd, ik, iw) = naznormi(ibnd, ik, iw) + weight * kernelm
                    ! Eqs.(21)-(22) in Margine and Giustino, PRB 87, 024505 (2013)
                    ! using kernelm and kernelp the sum over |wp| < wscut in Eqs. (21)-(22)
                    ! is rewritten as a sum over iwp = 1, nsiw(itemp)
                    aznormi(ibnd, ik, iw) = aznormi(ibnd, ik, iw) + weight * wesqrt(jbnd, ixkqf(ik, iq0), iwp) &
                                          * kernelm
                    adeltai(ibnd, ik, iw) = adeltai(ibnd, ik, iw) + weight * desqrt(jbnd, ixkqf(ik, iq0), iwp) &
                                          * (kernelp - 2.d0 * muc)
                  ENDDO ! iwp
                ENDDO ! iw
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    DEALLOCATE(wesqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating wesqrt', 1)
    DEALLOCATE(desqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating desqrt', 1)
    !
    ! collect contributions from all pools
    CALL mp_sum(aznormi, inter_pool_comm)
    CALL mp_sum(naznormi, inter_pool_comm)
    CALL mp_sum(adeltai, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      IF (iter == 1) THEN
        ALLOCATE(deltaold(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating deltaold', 1)
        deltaold(:) = gap0
      ENDIF
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) / dosef
              nznormi(iw) = nznormi(iw) + weight * naznormi(ibnd, ik, iw)
              znormi(iw) = znormi(iw) + weight * aznormi(ibnd, ik, iw)
              deltai(iw) = deltai(iw) + weight * adeltai(ibnd, ik, iw)
              naznormi(ibnd, ik, iw) = 1.d0 + pi * gtemp(itemp) * naznormi(ibnd, ik, iw) / wsi(iw)
              ! Eqs.(21)-(22) in Margine and Giustino, PRB 87, 024505 (2013)
              aznormi(ibnd, ik, iw) = 1.d0 + pi * gtemp(itemp) * aznormi(ibnd, ik, iw) / wsi(iw)
              adeltai(ibnd, ik, iw) = pi * gtemp(itemp) * adeltai(ibnd, ik, iw) / aznormi(ibnd, ik, iw)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
        nznormi(iw) = 1.d0 + pi * gtemp(itemp) * nznormi(iw) / wsi(iw)
        znormi(iw) = 1.d0 + pi * gtemp(itemp) * znormi(iw) / wsi(iw)
        deltai(iw) = pi * gtemp(itemp) * deltai(iw) / znormi(iw)
        reldelta = reldelta + ABS(deltai(iw) - deltaold(iw))
        absdelta = absdelta + ABS(deltai(iw))
      ENDDO ! iw
      errdelta = reldelta / absdelta
      deltaold(:) = deltai(:)
      !
      IF (iter == 1) &
        WRITE(stdout, '(5x, a)') '   iter      ethr        znormi [eV]    deltai [eV]'
      WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, znormi(1), deltai(1)
!      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
!                    '   ethr = ', errdelta, '   znormi(1) = ', znormi(1), &
!                    '   deltai(1) = ', deltai(1)
      !
      IF (errdelta < conv_thr_iaxis) conv = .TRUE.
      IF (conv .OR. iter == nsiter) THEN
        gap(itemp) = deltai(1)
        gap0 = gap(itemp)
        CALL eliashberg_write_iaxis(itemp)
      ENDIF
      !
      IF (conv .OR. iter == nsiter) THEN
        DEALLOCATE(deltaold, STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating deltaold', 1)
      ENDIF
      IF (conv) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was reached in nsiter = ', iter
        WRITE(stdout,'(a)') ' '
      ELSEIF (.NOT. conv .AND. iter == nsiter) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was not reached in nsiter = ', iter
        WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_iaxis'
        WRITE(stdout,'(a)') ' '
      ENDIF
    ENDIF ! ionode_id
    CALL mp_bcast(deltai, ionode_id, inter_pool_comm)
    CALL mp_bcast(znormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(nznormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(aznormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(naznormi, ionode_id, inter_pool_comm)
    CALL mp_bcast(gap0, ionode_id, inter_pool_comm)
    CALL mp_bcast(gap, ionode_id, inter_pool_comm)
    CALL mp_bcast(agap, ionode_id, inter_pool_comm)
    CALL mp_bcast(conv, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (conv .OR. iter == nsiter) THEN
      !
      ! remove memory allocated for wesqrt, desqrt, adeltaip, deltaold
      imelt = (1 + 3 * nbndfs * nkfs) * nsiw(itemp)
      CALL mem_size_eliashberg(2, -imelt)
      !
      IF (.NOT. limag_fly) THEN
        !
        DEALLOCATE(akeri, STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating akeri', 1)
        !
        ! remove memory allocated for akeri
        imelt = (upper_bnd - lower_bnd + 1) * MAXVAL(nqfs(:)) * nbndfs**2 * (2 * nsiw(itemp))
        CALL mem_size_eliashberg(2, -imelt)
        !
      ENDIF
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE analytic_cont_aniso_iaxis_to_raxis(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !
    ! This routine does the analytic continuation of the anisotropic Eliashberg equations
    ! from the imaginary-axis to the real axis
    ! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
    !
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : wqf, wf, gtemp
    USE epwcom,        ONLY : nqstep, degaussq, nsiter, conv_thr_racon, fsthick, &
                              lpade, eps_acustic
    USE eliashbergcom, ONLY : nsw, dwsph, ws, wsph, gap, agap, gp, gm, adsumi, azsumi, &
                              delta, znorm, adelta, adeltap, aznorm, aznormp, g2, lacon_fly, &
                              a2fij, wkfs, dosef, ixkqf, ixqfs, nqfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE supercond,     ONLY : gamma_acont
    USE constants_epw, ONLY : ci, zero, one, czero, cone
    USE constants,     ONLY : pi
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : mem_size_eliashberg, memlt_eliashberg
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on the iteration number
    !
    ! Local variables
    CHARACTER(LEN = 256) :: cname
    !! character in file name
    !
    INTEGER :: i, iw, iwp
    !! Counter on frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: iwph
    !! Counter on frequency in a2f
    INTEGER :: imode
    !! Counter on phonon modes
    REAL(KIND = DP) :: inv_degaussq
    !! 1.0/degaussq. Defined for efficiency reasons
    INTEGER :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: rgammap
    !! - bose_einstein(w') - fermi_dirac(w + w')
    REAL(KIND = DP) :: rgammam
    !!   bose_einstein(w') + fermi_dirac(-w + w')
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP) :: a2f_
    !! Temporary variable for Eliashberg spectral function
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !
    COMPLEX(KIND = DP) :: az2, ad2, esqrt, root
    !! Temporary variables
    COMPLEX(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    a2f_ = zero
    inv_degaussq = one / degaussq
    !
    IF (iter == 1) THEN
      !
      ! get the size of required allocated memory for
      ! delta, znorm, deltaold, adelta, adeltap, aznorm, aznormp, gp, gm
      IF (lpade) THEN
        imelt = 2 * ( 1 + 2 * nbndfs * nkfs ) * nsw + 2 * nqstep * nsw
      ELSE
        imelt = 2 * ( 3 + 4 * nbndfs * nkfs ) * nsw + 2 * nqstep * nsw
      ENDIF
      CALL mem_size_eliashberg(2, imelt)
      !
      IF (.NOT. lpade) THEN
        ALLOCATE(delta(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating delta', 1)
        ALLOCATE(adelta(nbndfs, nkfs, nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating adelta', 1)
        ALLOCATE(znorm(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating znorm', 1)
        ALLOCATE(aznorm(nbndfs, nkfs, nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating aznorm', 1)
      ENDIF
      ALLOCATE(adeltap(nbndfs, nkfs, nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating adeltap', 1)
      ALLOCATE(aznormp(nbndfs, nkfs, nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating aznormp', 1)
      ALLOCATE(deltaold(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating deltaold', 1)
      adeltap(:, :, :) = czero
      aznormp(:, :, :) = cone
      deltaold(:) = czero
      IF (lpade) THEN
        adeltap(:, :, :) = adelta(:, :, :)
        deltaold(:) = delta(:)
      ELSE
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              adeltap(ibnd, ik, :) = agap(ibnd, ik, itemp)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
        deltaold(:) = gap(itemp)
      ENDIF
      !
      ALLOCATE(gp(nsw, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating gp', 1)
      ALLOCATE(gm(nsw, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error allocating gm', 1)
      !
      ! ! Eq.(28) in Margine and Giustino, PRB 87, 024505 (2013)
      DO iw = 1, nsw ! loop over omega
        DO iwp = 1, nqstep ! loop over omega_prime
          CALL gamma_acont(ws(iw), ws(iwp), gtemp(itemp), rgammap, rgammam)
          gp(iw, iwp) = rgammap
          gm(iw, iwp) = rgammam
        ENDDO
      ENDDO
      CALL kernel_aniso_iaxis_analytic_cont(itemp)
      CALL memlt_eliashberg(itemp, 'acon')
      IF (.NOT. lacon_fly) CALL evaluate_a2fij()
    ENDIF
    delta(:) = czero
    znorm(:) = czero
    adelta(:, :, :) = czero
    aznorm(:, :, :) = czero
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik,iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                IF (lacon_fly) THEN ! evaluate a2fij on the fly
                  DO imode = 1, nmodes
                    IF (wf(imode, iq0) > eps_acustic) THEN
                      DO iwph = 1, nqstep
                        weight = w0gauss((wsph(iwph) - wf(imode, iq0)) * inv_degaussq, 0) * inv_degaussq
                        a2f_ = weight * dosef * g2(ik, iq, ibnd, jbnd, imode)
                      ENDDO ! iwph
                    ENDIF ! wf
                  ENDDO ! imode
                ENDIF ! lacon_fly
                !
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) / dosef
                DO iw = 1, nsw ! loop over omega
                  DO iwp = 1, nqstep ! loop over omega_prime
                    !
                    i = iw + iwp - 1
                    IF (i <= nsw) THEN
                      az2 = aznormp(jbnd, ixkqf(ik, iq0), i) * aznormp(jbnd, ixkqf(ik, iq0), i)
                      ad2 = adeltap(jbnd, ixkqf(ik, iq0), i) * adeltap(jbnd, ixkqf(ik, iq0), i)
                      root = SQRT(az2 * (ws(i) * ws(i) - ad2))
                      IF (AIMAG(root) < zero) THEN
                        esqrt = aznormp(jbnd, ixkqf(ik, iq0), i) / CONJG(root)
                      ELSE
                        esqrt = aznormp(jbnd, ixkqf(ik, iq0), i) / root
                      ENDIF
                      IF (lacon_fly) THEN
                        esqrt = esqrt * weight * gp(iw, iwp) * a2f_
                      ELSE
                        esqrt = esqrt * weight * gp(iw, iwp) * a2fij(ik, iq, ibnd, jbnd, iwp)
                      ENDIF
                      aznorm(ibnd, ik, iw) = aznorm(ibnd, ik, iw) - ws(i) * esqrt
                      adelta(ibnd, ik, iw) = adelta(ibnd, ik, iw) - adeltap(jbnd, ixkqf(ik, iq0), i) * esqrt
                    ENDIF
                    !
                    i = ABS(iw - iwp) + 1
                    az2 = aznormp(jbnd, ixkqf(ik, iq0), i) * aznormp(jbnd, ixkqf(ik, iq0), i)
                    ad2 = adeltap(jbnd, ixkqf(ik, iq0), i) * adeltap(jbnd, ixkqf(ik, iq0), i)
                    root = SQRT(az2 * (ws(i) * ws(i) - ad2))
                    IF (AIMAG(root) < zero) THEN
                      esqrt = aznormp(jbnd, ixkqf(ik, iq0), i) / CONJG(root)
                    ELSE
                      esqrt = aznormp(jbnd, ixkqf(ik, iq0), i) / root
                    ENDIF
                    esqrt = esqrt * weight * gm(iw, iwp) * a2fij(ik, iq, ibnd, jbnd, iwp)
                    IF (iw < iwp) THEN
                      aznorm(ibnd, ik, iw) = aznorm(ibnd, ik, iw) - ws(i) * esqrt
                    ELSE
                      aznorm(ibnd, ik, iw) = aznorm(ibnd, ik, iw) + ws(i) * esqrt
                    ENDIF
                    adelta(ibnd, ik, iw) = adelta(ibnd, ik, iw) + adeltap(jbnd, ixkqf(ik, iq0), i) * esqrt
                  ENDDO ! iwp
                ENDDO ! iw
              ENDIF ! fsthick
            ENDDO ! jbnd
          ENDDO ! iq
          ! Eqs.(26)-(27) in Margine and Giustino, PRB 87, 024505 (2013)
          DO iw = 1, nsw ! loop over omega
            aznorm(ibnd, ik, iw) = - gtemp(itemp) * azsumi(ibnd, ik, iw) + ci * aznorm(ibnd, ik, iw) * dwsph
            adelta(ibnd, ik, iw) =   gtemp(itemp) * adsumi(ibnd, ik, iw) + ci * adelta(ibnd, ik, iw) * dwsph
          ENDDO ! iw
        ENDIF ! fsthick
      ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools
    CALL mp_sum(aznorm, inter_pool_comm)
    CALL mp_sum(adelta, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsw ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) / dosef
              znorm(iw) = znorm(iw) + weight * aznorm(ibnd, ik, iw)
              delta(iw) = delta(iw) + weight * adelta(ibnd, ik, iw)
              ! Eqs.(26)-(27) in Margine and Giustino, PRB 87, 024505 (2013)
              aznorm(ibnd, ik, iw) = 1.d0 + pi * aznorm(ibnd, ik, iw) / ws(iw)
              adelta(ibnd, ik, iw) = pi * adelta(ibnd, ik, iw) / aznorm(ibnd, ik, iw)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
        znorm(iw) = 1.0d0 + pi * znorm(iw) / ws(iw)
        delta(iw) = pi * delta(iw) / znorm(iw)
        reldelta = reldelta + ABS(delta(iw) - deltaold(iw))
        absdelta = absdelta + ABS(delta(iw))
      ENDDO ! iw
      errdelta = reldelta / absdelta
      deltaold(:) = delta(:)
      !
      IF (iter == 1) &
        WRITE(stdout, '(5x, a)') '   iter      ethr      Re[znorm] [eV] Re[delta] [eV]'
      WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, REAL(znorm(1)), REAL(delta(1))
!      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
!                    '   ethr = ', errdelta, '   Re[znorm(1)] = ', REAL(znorm(1)), &
!                    '   Re[delta(1)] = ', REAL(delta(1))
      !
      IF (errdelta < conv_thr_racon) conv = .TRUE.
      IF (conv .OR. iter == nsiter) THEN
        cname = 'acon'
        CALL eliashberg_write_raxis(itemp, cname)
      ENDIF
      !
      IF (conv) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was reached in nsiter = ', iter
        WRITE(stdout,'(a)') ' '
      ELSEIF (.NOT. conv .AND. iter == nsiter) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was not reached in nsiter = ', iter
        WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_racon'
        WRITE(stdout,'(a)') ' '
      ENDIF
    ENDIF
    CALL mp_bcast(delta, ionode_id, inter_pool_comm)
    CALL mp_bcast(znorm, ionode_id, inter_pool_comm)
    CALL mp_bcast(aznorm, ionode_id, inter_pool_comm)
    CALL mp_bcast(conv, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (conv .OR. iter == nsiter) THEN
      !
      DEALLOCATE(deltaold, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error deallocating deltaold', 1)
      DEALLOCATE(gp, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error deallocating gp', 1)
      DEALLOCATE(gm, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error deallocating gm', 1)
      DEALLOCATE(adsumi, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error deallocating adsumi', 1)
      DEALLOCATE(azsumi, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error deallocating azsumi', 1)
      !
      ! remove memory allocated for deltaold, gp, gm, adsumi, azsumi
      imelt = 2 * nsw + 2 * nqstep * nsw + 2 * (upper_bnd - lower_bnd + 1) * nbndfs * nsw
      CALL mem_size_eliashberg(2, -imelt)
      !
      IF (.NOT. lacon_fly) THEN
        !
        DEALLOCATE(a2fij, STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso_iaxis_to_raxis', 'Error deallocating a2fij', 1)
        !
        ! remove memory allocated for a2fij
        imelt = (upper_bnd - lower_bnd + 1) * MAXVAL(nqfs(:)) * nbndfs**2 * nqstep
        CALL mem_size_eliashberg(2, -imelt)
        !
      ENDIF
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE analytic_cont_aniso_iaxis_to_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE pade_cont_aniso_iaxis_to_raxis(itemp, N)
    !-----------------------------------------------------------------------
    !
    ! This routine uses pade approximants to continue the anisotropic Eliashberg equations
    ! from the imaginary-axis to the real-axis
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : fsthick
    USE eliashbergcom, ONLY : nsw, ws, nsiw, wsi, delta, znorm, &
                              adelta, aznorm, adeltai, aznormi, &
                              wkfs, dosef, w0g, nkfs, nbndfs, ef0, ekfs
    USE utilities,       ONLY : pade_coeff, pade_eval
    USE constants_epw, ONLY : cone, ci, zero, czero
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: N
    !! Nr. of frequency points in the Pade approx
    !
    ! Local variable
    CHARACTER(LEN = 256) :: cname
    !! character in file name
    !
    INTEGER :: iw
    !! Counter on frequency imag- and real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: weight
    !! Temporary variable
    !
    ! arrays used in pade_coeff and pade_eval
    COMPLEX(KIND = DP) :: omega
    !! frequency real-axis
    COMPLEX(KIND = DP) :: padapp
    !! znorm or delta on real-axis after pade_eval
    COMPLEX(KIND = DP), ALLOCATABLE :: a(:)
    !! a - pade coeff for deltai
    COMPLEX(KIND = DP), ALLOCATABLE :: b(:)
    !! b - pade coeff for znormi
    COMPLEX(KIND = DP), ALLOCATABLE :: z(:)
    !! z - frequency imag-axis
    COMPLEX(KIND = DP), ALLOCATABLE :: u(:)
    !! u - deltai
    COMPLEX(KIND = DP), ALLOCATABLE :: v(:)
    !! v - znormi
    !
    ! get the size of required allocated memory for
    ! a, b, z, u, v, delta, znorm, adelta, aznorm
    imelt = 2 * 5 * N + 2 * (2 + 2 * nbndfs * nkfs) * nsw
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(delta(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating delta', 1)
    ALLOCATE(znorm(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating znorm', 1)
    ALLOCATE(adelta(nbndfs, nkfs, nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating adelta', 1)
    ALLOCATE(aznorm(nbndfs, nkfs, nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating aznorm', 1)
    ALLOCATE(a(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating a', 1)
    ALLOCATE(b(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating b', 1)
    ALLOCATE(z(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating z', 1)
    ALLOCATE(u(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating u', 1)
    ALLOCATE(v(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error allocating v', 1)
    delta(:) = czero
    znorm(:) = czero
    adelta(:, :, :) = czero
    aznorm(:, :, :) = czero
    a(:) = czero
    b(:) = czero
    z(:) = czero
    u(:) = czero
    v(:) = czero
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iw = 1, N
            z(iw) = ci * wsi(iw)
            u(iw) = cone * adeltai(ibnd, ik, iw)
            v(iw) = cone * aznormi(ibnd, ik, iw)
          ENDDO
          CALL pade_coeff(N, z, u, a)
          CALL pade_coeff(N, z, v, b)
          DO iw = 1, nsw
            omega = cone * ws(iw)
            CALL pade_eval(N, z, a, omega, padapp)
            adelta(ibnd, ik, iw) = padapp
            CALL pade_eval(N, z, b, omega, padapp)
            aznorm(ibnd, ik, iw) = padapp
          ENDDO
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools
    CALL mp_sum(aznorm, inter_pool_comm)
    CALL mp_sum(adelta, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      DO iw = 1, nsw ! loop over omega
        DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) / dosef
              znorm(iw) = znorm(iw) + weight * aznorm(ibnd, ik, iw)
              delta(iw) = delta(iw) + weight * adelta(ibnd, ik, iw)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      !
      WRITE(stdout, '(5x, a)') '   pade Re[znorm] [eV] Re[delta] [eV]'
      WRITE(stdout, '(5x, i6, 2ES15.6)') N, REAL(znorm(1)), REAL(delta(1))
!      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6)') 'pade = ', N, &
!                    '   Re[znorm(1)] = ', REAL(znorm(1)), &
!                    '   Re[delta(1)] = ', REAL(delta(1))
      WRITE(stdout, '(a)') ' '
      !
      WRITE(stdout, '(5x, a, i6, a)') 'Convergence was reached for N = ', N, ' Pade approximants'
      WRITE(stdout, '(a)') ' '
      !
      cname = 'pade'
      CALL eliashberg_write_raxis(itemp, cname)
    ENDIF
    CALL mp_bcast(delta, ionode_id, inter_pool_comm)
    CALL mp_bcast(znorm, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(a, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error deallocating a', 1)
    DEALLOCATE(b, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error deallocating b', 1)
    DEALLOCATE(z, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error deallocating z', 1)
    DEALLOCATE(u, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error deallocating u', 1)
    DEALLOCATE(v, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso_iaxis_to_raxis', 'Error deallocating v', 1)
    !
    ! remove memory allocated for a, b, z, u, v
    imelt = 2 * 5 * N
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE pade_cont_aniso_iaxis_to_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis(itemp)
    !-----------------------------------------------------------------------
    !!
    !! Compute kernels K_{+}(ik, iq, ibnd, jbnd; n, n', T) and
    !! K_{-}(ik, iq, ibnd, jbnd; n, n', T) and store them in memory
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : fsthick
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, akeri, ekfs, ef0, ixkqf, ixqfs, nqfs
    USE constants_epw, ONLY : zero
    USE constants,     ONLY : pi
    USE division,      ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: iw, n
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bandst k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: omega
    !! frequency imag-axis
    REAL(KIND = DP) :: lambda_eph
    !! electron-phonon coupling
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ALLOCATE(akeri(lower_bnd:upper_bnd, MAXVAL(nqfs(:)), nbndfs, nbndfs, 2 * nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis', 'Error allocating akeri', 1)
    akeri(:, :, :, :, :) = zero
    !
    ! RM - if lambdar_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                DO iw = 1, 2*nsiw(itemp)
                  n = iw - 1
                  omega = DBLE(2 * n) * pi * gtemp(itemp)
                  CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, lambda_eph)
                  !CALL lambdar_aniso_ver2(ik, iq, ibnd, jbnd, omega, lambda_eph)
                  akeri(ik, iq, ibnd, jbnd, iw) = lambda_eph
                ENDDO ! iw
              ENDIF ! fsthick
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF ! fsthick
      ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, lambda_eph)
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik, iq, ibnd, jbnd; n-n')
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands k+q
    !
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency on imag-axis
    REAL(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; n-n')
    !
    ! Local variables
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: imode
    !! Counter on phonon modes
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik, iq)
    lambda_eph = zero
    DO imode = 1, nmodes  ! loop over frequency modes
      IF (wf(imode, iq0) > eps_acustic) THEN
        lambda_eph = lambda_eph + g2(ik, iq, ibnd, jbnd, imode) * wf(imode, iq0) &
                   / (wf(imode, iq0) * wf(imode, iq0) + omega * omega)
      ENDIF
    ENDDO
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdar_aniso_ver1
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver2(ik, iq, ibnd, jbnd, omega, lambda_eph)
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik, iq, ibnd, jbnd; n-n')
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, wsph, dwsph
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands at k+q
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency on imag-axis
    REAL(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; n-n')
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    ! Eq.(18) in Margine and Giustino, PRB 87, 024505 (2013)
    lambda_eph = zero
    DO iwph = 1, nqstep
      lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik, iq, ibnd, jbnd, iwph) &
                 / (wsph(iwph) * wsph(iwph) + omega * omega)
    ENDDO
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdar_aniso_ver2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis_analytic_cont(itemp)
    !-----------------------------------------------------------------------
    !!
    !! computes kernels K_{+}(w, iw_n, T) and K_{-}(w, iw_n, T)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : wqf
    USE epwcom,        ONLY : muc, fsthick
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, adeltai, nkfs, nbndfs, dosef, ixkqf, ixqfs, nqfs, &
                              w0g, ekfs, ef0, adsumi, azsumi
    USE constants_epw, ONLY : zero
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: esqrt
    !! 1 / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: kernelr
    !! 2 * Re[lambda(w - iw_n)]
    REAL(KIND = DP) :: kerneli
    !! 2 * Im[lambda(w - iw_n)]
    REAL(KIND = DP) :: weight
    !! factor in supercond. gap equations
    REAL(KIND = DP), ALLOCATABLE :: wesqrt(:, :, :)
    !! w / sqrt{w^2+\delta^2}
    REAL(KIND = DP), ALLOCATABLE :: desqrt(:, :, :)
    !! \delta / sqrt{w^2+\delta^2}
    !
    COMPLEX(KIND = DP) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; w-iw_n)
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get memory size required for wesqrt, desqrt, adsumi, azsumi
    imelt = 2 * nbndfs * nkfs * nsiw(itemp) + 2 * (upper_bnd - lower_bnd + 1) * nbndfs * nsw
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(wesqrt(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis_analytic_cont', 'Error allocating wesqrt', 1)
    ALLOCATE(desqrt(nbndfs, nkfs, nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis_analytic_cont', 'Error allocating desqrt', 1)
    !
    ALLOCATE(adsumi(nbndfs, lower_bnd:upper_bnd, nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis_analytic_cont', 'Error allocating adsumi', 1)
    ALLOCATE(azsumi(nbndfs, lower_bnd:upper_bnd, nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis_analytic_cont', 'Error allocating azsumi', 1)
    adsumi(:, :, :) = zero
    azsumi(:, :, :) = zero
    !
    ! RM - if lambdai_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) / dosef
                DO iw = 1, nsw ! loop over omega
                  DO iwp = 1, nsiw(itemp) ! loop over iw_n
                    CALL lambdai_aniso_ver1(ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph)
                    !CALL lambdai_aniso_ver2(ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph)
                    kernelr = 2.d0 * REAL(lambda_eph)
                    kerneli = 2.d0 * AIMAG(lambda_eph)
                    IF (iw == 1) THEN
                      esqrt = 1.d0 / DSQRT(wsi(iwp) * wsi(iwp) + &
                                           adeltai(jbnd, ixkqf(ik, iq0), iwp) * &
                                           adeltai(jbnd, ixkqf(ik, iq0), iwp))
                      wesqrt(jbnd, ixkqf(ik, iq0), iwp) =  wsi(iwp) * esqrt
                      desqrt(jbnd, ixkqf(ik, iq0), iwp) =  adeltai(jbnd, ixkqf(ik, iq0), iwp) * esqrt
                    ENDIF
                    azsumi(ibnd, ik, iw) = azsumi(ibnd, ik, iw) &
                                         + weight * wesqrt(jbnd, ixkqf(ik, iq0), iwp) * kerneli
                    adsumi(ibnd, ik, iw) = adsumi(ibnd, ik, iw) &
                                         + weight * desqrt(jbnd, ixkqf(ik, iq0), iwp) * (kernelr - 2.d0 * muc)
                  ENDDO ! iwp
                ENDDO ! iw
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    DEALLOCATE(wesqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis_analytic_cont', 'Error deallocating wesqrt', 1)
    DEALLOCATE(desqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis_analytic_cont', 'Error deallocating desqrt', 1)
    !
    ! remove memory allocated for wesqrt, desqrt
    imelt = 2 * nbndfs * nkfs * nsiw(itemp)
    CALL mem_size_eliashberg (2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_aniso_iaxis_analytic_cont
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_aniso_ver1(ik, iq, ibnd, jbnd, omega, omegap, lambda_eph)
    !-----------------------------------------------------------------------
    !!
    !! computes lambda_ij(k, k+q; w-iw_n)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : ci, czero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands at k+q
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency w at point iw on real-axis
    REAL(KIND = DP), INTENT(in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; w-iw_n)
    !
    ! Local variables
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: imode
    !! Counter on phonon modes
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik, iq)
    lambda_eph = czero
    DO imode = 1, nmodes  ! loop over frequency modes
      IF (wf(imode, iq0) > eps_acustic) THEN
        lambda_eph = lambda_eph +  g2(ik, iq, ibnd, jbnd, imode) * wf(imode, iq0) &
                   / (wf(imode, iq0) * wf(imode, iq0) - (omega - ci * omegap) * (omega - ci * omegap))
      ENDIF
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdai_aniso_ver1
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_aniso_ver2(ik, iq, ibnd, jbnd, omega, omegap, lambda_eph)
    !-----------------------------------------------------------------------
    !!
    !! computes lambda(w-iw_n)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, dwsph, wsph
    USE constants_epw, ONLY : ci, czero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands at k+q
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency w at point iw on real-axis
    REAL(KIND = DP), INTENT(in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; w-iw_n)
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    lambda_eph = czero
    DO iwph = 1, nqstep
      lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik, iq, ibnd, jbnd, iwph) &
                 / (wsph(iwph) * wsph(iwph) - (omega - ci * omegap) * (omega - ci * omegap))
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdai_aniso_ver2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2fij
    !-----------------------------------------------------------------------
    !!
    !! computes the anisotropic spectral function a2F(k, k', w)
    !!
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, a2fij, ixkqf, ixqfs, nqfs, ekfs, ef0, &
                              dosef, wsph
    USE constants_epw, ONLY : zero, one
    USE division,      ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: iwph
    !! Counter on frequency
    REAL(KIND = DP) :: inv_degaussq
    !! 1.0/degaussq. Defined for efficiency reasons
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: weight
    !! Factor in a2fij
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    ALLOCATE(a2fij(lower_bnd:upper_bnd, MAXVAL(nqfs(:)), nbndfs, nbndfs, nqstep), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2fij', 'Error deallocating a2fij', 1)
    a2fij(:, :, :, :, :) = zero
    !
    inv_degaussq = one / degaussq
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                DO imode = 1, nmodes
                  IF (wf(imode, iq0) > eps_acustic) THEN
                    DO iwph = 1, nqstep
                      weight  = w0gauss((wsph(iwph) - wf(imode, iq0) ) * inv_degaussq, 0) * inv_degaussq
                      a2fij(ik, iq, ibnd, jbnd, iwph) = a2fij(ik, iq, ibnd, jbnd, iwph) &
                                             + weight * dosef * g2(ik, iq, ibnd, jbnd, imode)
                    ENDDO ! iwph
                  ENDIF ! wf
                ENDDO ! imode
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE evaluate_a2fij
    !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END MODULE supercond_aniso
  !-----------------------------------------------------------------------
