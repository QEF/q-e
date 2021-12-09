  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE supercond_iso
  !----------------------------------------------------------------------
  !!
  !! This module contains all the subroutines linked with superconductivity
  !! using the isotropic Eliashberg formalism.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_iso_iaxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the isotropic
    !! Eliashberg equations on the imaginary-axis.
    !!
    !
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE control_flags,     ONLY : iverbosity
    USE epwcom,            ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                                  limag, lpade, lacon, npade
    USE elph2,             ONLY : gtemp
    USE eliashbergcom,     ONLY : nsw, nsiw, deltai, deltaip, delta, deltap
    USE constants_epw,     ONLY : kelvin2eV, ci, zero
    USE supercond,         ONLY : free_energy, dos_quasiparticle, gen_freqgrid_iaxis, &
                                  eliashberg_grid
    USE utilities,         ONLY : mix_broyden
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
    REAL(KIND = DP), ALLOCATABLE :: df1(:, :), df2(:, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv1(:, :), dv2(:, :)
    !! Temporary variables for mix_broyden
    !
    CALL start_clock('iso_iaxis')
    !
    CALL eliashberg_grid()
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      CALL prtheader_supercond(itemp, 1)
      CALL start_clock('iaxis_imag')
      CALL gen_freqgrid_iaxis(itemp)
      !
      IF (limag) THEN
        !
        ALLOCATE(df1(nsiw(itemp), broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating df1', 1)
        ALLOCATE(dv1(nsiw(itemp), broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating dv1', 1)
        df1(:, :) = zero
        dv1(:, :) = zero
        !
        iter = 1
        conv = .FALSE.
        !
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL sum_eliashberg_iso_iaxis(itemp, iter, conv)
          CALL mix_broyden(nsiw(itemp), deltai, deltaip, broyden_beta, &
                           iter, broyden_ndim, conv, df1, dv1)
          iter = iter + 1
        ENDDO ! iter
        !
        DEALLOCATE(df1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating df1', 1)
        DEALLOCATE(dv1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating dv1', 1)
        !
        IF (conv) THEN
          !
          IF (iverbosity == 2) THEN
            CALL free_energy(itemp)
          ENDIF
          CALL stop_clock('iaxis_imag')
          CALL print_clock('iaxis_imag')
          WRITE(stdout, '(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          CALL deallocate_iso_iaxis()
          CALL deallocate_iso()
          CALL stop_clock('iaxis_imag' )
          CALL print_clock('iaxis_imag')
          WRITE(stdout, '(a)') ' '
          RETURN
        ENDIF
      ENDIF
      !
      IF (lpade) THEN
        CALL prtheader_supercond(itemp, 2)
        CALL start_clock('raxis_pade')
        N = npade * nsiw(itemp) / 100
        IF (mod(N, 2) /= 0) N = N + 1
        CALL pade_cont_iso(itemp, N)
        !
        CALL dos_quasiparticle(itemp)
        CALL stop_clock('raxis_pade')
        CALL print_clock('raxis_pade')
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      IF (lacon) THEN
        CALL prtheader_supercond(itemp, 3)
        CALL start_clock('raxis_acon')
        !
        ALLOCATE(rdeltain(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating rdeltain', 1)
        ALLOCATE(cdeltain(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating cdeltain', 1)
        ALLOCATE(rdeltaout(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating rdeltaout', 1)
        ALLOCATE(cdeltaout(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating cdeltaout', 1)
        ALLOCATE(df1(nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating df1', 1)
        ALLOCATE(dv1(nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating dv1', 1)
        ALLOCATE(df2(nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating df2', 1)
        ALLOCATE(dv2(nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating dv2', 1)
        rdeltain(:)  = zero
        cdeltain(:)  = zero
        rdeltaout(:) = zero
        cdeltaout(:) = zero
        df1(:, :) = zero
        dv1(:, :) = zero
        df2(:, :) = zero
        dv2(:, :) = zero
        !
        iter = 1
        conv = .FALSE.
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL analytic_cont_iso(itemp, iter, conv)
          rdeltain(:)  =  REAL(deltap(:))
          cdeltain(:)  = AIMAG(deltap(:))
          rdeltaout(:) =  REAL(delta(:))
          cdeltaout(:) = AIMAG(delta(:))
          CALL mix_broyden(nsw, rdeltaout, rdeltain, broyden_beta, &
                           iter, broyden_ndim, conv, df1, dv1)
          CALL mix_broyden(nsw, cdeltaout, cdeltain, broyden_beta, &
                           iter, broyden_ndim, conv, df2, dv2)
          deltap(:) = rdeltain(:) + ci * cdeltain(:)
          iter = iter + 1
        ENDDO ! iter
        !
        DEALLOCATE(rdeltain, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating rdeltain', 1)
        DEALLOCATE(cdeltain, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating cdeltain', 1)
        DEALLOCATE(rdeltaout, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating rdeltaout', 1)
        DEALLOCATE(cdeltaout, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating cdeltaout', 1)
        DEALLOCATE(df1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating df1', 1)
        DEALLOCATE(dv1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating dv1', 1)
        DEALLOCATE(df2, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating df2', 1)
        DEALLOCATE(dv2, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating dv2', 1)
        !
        IF (conv) THEN
          CALL dos_quasiparticle(itemp)
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout, '(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          CALL deallocate_iso_iaxis()
          CALL deallocate_iso_raxis()
          CALL deallocate_iso()
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout, '(a)') ' '
          RETURN
        ENDIF
      ENDIF
      !
      CALL deallocate_iso_iaxis()
      !
      tcpu = get_clock('iso_iaxis')
      WRITE(stdout, '(5x, a, i3, a, f8.1, a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
      WRITE(stdout, '(a)') ' '
      !
      IF (lpade .OR. lacon) CALL deallocate_iso_raxis()
      !
    ENDDO ! itemp
    !
    CALL deallocate_iso()
    !
    CALL stop_clock('iso_iaxis')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_iso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_iso_iaxis(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the isotropic Eliashberg equations on the imaginary-axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsiw, gap0, gap, wsi, keri, &
                              deltai, deltaip, znormi, nznormi
    USE constants_epw, ONLY : zero
    USE constants,     ONLY : pi
    USE io_eliashberg, ONLY : eliashberg_write_iaxis
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
    REAL(KIND = DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    !! Temporary variables
    REAL(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    ALLOCATE(wesqrt(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating wesqrt', 1)
    ALLOCATE(desqrt(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating desqrt', 1)
    !
    IF (iter == 1) THEN
      ALLOCATE(gap(nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating gap', 1)
      ALLOCATE(deltai(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating deltai', 1)
      ALLOCATE(deltaip(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating deltaip', 1)
      ALLOCATE(znormi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating znormi', 1)
      ALLOCATE(nznormi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating nznormi', 1)
      ALLOCATE(deltaold(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating deltaold', 1)
      gap(itemp) = zero
      deltaip(:) = zero
      deltaip(:) = gap0
      deltaold(:) = gap0
      !
      CALL kernel_iso_iaxis(itemp)
    ENDIF
    deltai(:) = zero
    znormi(:) = zero
    nznormi(:) = zero
    !
    absdelta = zero
    reldelta = zero
    DO iw = 1, nsiw(itemp) ! loop over omega
      DO iwp = 1, nsiw(itemp) ! loop over omega_prime
        ! this step is performed at each iter step only for iw=1 since it is independ of wsi(iw)
        IF (iw == 1) THEN
          esqrt = 1.d0 / DSQRT(wsi(iwp)**2.d0 + deltaip(iwp)**2.d0)
          wesqrt(iwp) = wsi(iwp) * esqrt
          desqrt(iwp) = deltaip(iwp) * esqrt
        ENDIF
        lambdam = keri(ABS(iw - iwp) + 1)
        lambdap = keri(ABS(iw + iwp))
        ! Eq. (4.4) in Picket, PRB 26, 1186 (1982)
        kernelm = lambdam - lambdap
        kernelp = lambdam + lambdap
        nznormi(iw) = nznormi(iw) + kernelm
        ! Eqs.(34)-(35) in Margine and Giustino, PRB 87, 024505 (2013)
        ! using kernelm and kernelp the sum over |wp| < wscut in Eqs. (34)-(35)
        ! is rewritten as a sum over iwp = 1, nsiw(itemp)
        znormi(iw) = znormi(iw) + wesqrt(iwp) * kernelm
        deltai(iw) = deltai(iw) + desqrt(iwp) * (kernelp - 2.d0 * muc)
      ENDDO ! iwp
      znormi(iw) = 1.d0 + pi * gtemp(itemp) * znormi(iw) / wsi(iw)
      ! Eqs.(34)-(35) in Margine and Giustino, PRB 87, 024505 (2013)
      nznormi(iw) = 1.d0 + pi * gtemp(itemp) * nznormi(iw) / wsi(iw)
      deltai(iw) = pi * gtemp(itemp) * deltai(iw) / znormi(iw)
      reldelta = reldelta + ABS(deltai(iw) - deltaold(iw))
      absdelta = absdelta + ABS(deltai(iw))
    ENDDO ! iw
    errdelta = reldelta / absdelta
    deltaold(:) = deltai(:)
    !
    IF (iter == 1) &
      WRITE(stdout, '(5x, a)') '   iter      ethr          znormi      deltai [meV]'
    WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, znormi(1), deltai(1) * 1000.d0
!    WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
!                  '   ethr = ', errdelta, '   znormi(1) = ', znormi(1), &
!                  '   deltai(1) = ', deltai(1)
    !
    IF (errdelta < conv_thr_iaxis) conv = .TRUE.
    IF (errdelta < conv_thr_iaxis .OR. iter == nsiter) THEN
      gap(itemp) = deltai(1)
      gap0 = gap(itemp)
      CALL eliashberg_write_iaxis(itemp)
    ENDIF
    !
    IF (conv) THEN
      WRITE(stdout, '(5x, a, i6)') 'Convergence was reached in nsiter = ', iter
      WRITE(stdout, '(a)') ' '
    ELSEIF (.NOT. conv .AND. iter == nsiter) THEN
      WRITE(stdout, '(5x, a, i6)') 'Convergence was not reached in nsiter = ', iter
      WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_iaxis'
      WRITE(stdout, '(a)') ' '
    ENDIF
    !
    DEALLOCATE(wesqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error deallocating wesqrt', 1)
    DEALLOCATE(desqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error deallocating desqrt', 1)
    !
    IF (conv .OR. iter == nsiter) THEN
      DEALLOCATE(deltaold, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error deallocating deltaold', 1)
      DEALLOCATE(keri, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error deallocating keri', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_iso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE analytic_cont_iso(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! This routine does the analyic continuation of the isotropic Eliashberg equations
    !! from the imaginary-axis to the real-axis
    !! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, nsiter, conv_thr_racon, lpade, nstemp
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsw, dwsph, ws, gap, a2f_iso, dsumi, zsumi, &
                              delta, deltap, znorm, znormp, gp, gm
    USE constants_epw, ONLY : ci, zero, czero, cone
    USE constants,     ONLY : pi
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    USE supercond,     ONLY : gamma_acont
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
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: rgammap
    !! - bose_einstein(w') - fermi_dirac(w + w')
    REAL(KIND = DP) :: rgammam
    !!   bose_einstein(w') + fermi_dirac(-w + w')
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    !
    COMPLEX(KIND = DP) :: esqrt, root
    !! Temporary variables
    COMPLEX(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    IF (iter == 1) THEN
      IF (.NOT. lpade) THEN
        ALLOCATE(delta(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating delta', 1)
        ALLOCATE(znorm(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating znorm', 1)
      ENDIF
      ALLOCATE(deltap(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating deltap', 1)
      ALLOCATE(znormp(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating znormp', 1)
      ALLOCATE(deltaold(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating deltaold', 1)
      ALLOCATE(gp(nsw, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating gp', 1)
      ALLOCATE(gm(nsw, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error allocating gm', 1)
      deltap(:) = czero
      znormp(:) = cone
      deltaold(:) = czero
      IF (lpade) THEN
        deltap(:)   = delta(:)
        deltaold(:) = delta(:)
      ELSE
        deltap(:)   = gap(itemp)
        deltaold(:) = gap(itemp)
      ENDIF
      !
      !! Eq.(28) in Margine and Giustino, PRB 87, 024505 (2013)
      DO iw = 1, nsw ! loop over omega
        DO iwp = 1, nqstep ! loop over omega_prime
          CALL gamma_acont(ws(iw), ws(iwp), gtemp(itemp), rgammap, rgammam)
          gp(iw, iwp) = rgammap
          gm(iw, iwp) = rgammam
        ENDDO
      ENDDO
      CALL kernel_iso_analytic_cont(itemp)
    ENDIF ! iter
    !
    delta(:) = czero
    znorm(:) = czero
    !
    absdelta = zero
    reldelta = zero
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nqstep ! loop over omega_prime
        !
        i = iw + iwp - 1
        IF (i <= nsw) THEN
          root = SQRT(znormp(i) * znormp(i) * (ws(i) * ws(i) - deltap(i) * deltap(i)))
          IF (AIMAG(root) < zero) THEN
            esqrt = znormp(i) / CONJG(root)
          ELSE
            esqrt = znormp(i) / root
          ENDIF
          esqrt = esqrt * gp(iw, iwp) * a2f_iso(iwp)
          znorm(iw) = znorm(iw) - ws(i) * esqrt
          delta(iw) = delta(iw) - deltap(i) * esqrt
        ENDIF
        !
        i = ABS(iw - iwp) + 1
        root = SQRT(znormp(i) * znormp(i) * (ws(i) * ws(i) - deltap(i) * deltap(i)))
        IF (AIMAG(root) < zero) THEN
          esqrt = znormp(i) / CONJG(root)
        ELSE
          esqrt = znormp(i) / root
        ENDIF
        esqrt = esqrt * gm(iw, iwp) * a2f_iso(iwp)
        IF (iw < iwp) THEN
          znorm(iw) = znorm(iw) - ws(i) * esqrt
        ELSE
          znorm(iw) = znorm(iw) + ws(i) * esqrt
        ENDIF
        delta(iw) = delta(iw) + deltap(i) * esqrt
      ENDDO ! iwp
      znorm(iw) = 1.d0 + pi * (- gtemp(itemp) * zsumi(iw) + ci * znorm(iw) * dwsph) / ws(iw)
      delta(iw) = pi * (gtemp(itemp) * dsumi(iw) + ci * delta(iw) * dwsph) / znorm(iw)
      reldelta = reldelta + ABS(delta(iw) - deltaold(iw))
      absdelta = absdelta + ABS(delta(iw))
    ENDDO ! iw
    errdelta = reldelta / absdelta
    deltaold(:) = delta(:)
    !
    IF (iter == 1) &
      WRITE(stdout, '(5x, a)') '   iter      ethr      Re[znorm]      Re[delta] [meV]'
      WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, REAL(znorm(1)), REAL(delta(1)) * 1000.d0
!      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
!                    '   ethr = ', errdelta, '   Re[znorm(1)] = ', REAL(znorm(1)), &
!                    '   Re[delta(1)] = ', REAL(delta(1))
    !
    IF (errdelta < conv_thr_racon) conv = .TRUE.
    IF (errdelta < conv_thr_racon .OR. iter == nsiter) THEN
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
    !
    IF (conv .OR. iter == nsiter) THEN
      DEALLOCATE(deltaold, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_iso', 'Error deallocating deltaold', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE analytic_cont_iso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE pade_cont_iso(itemp, N)
    !-----------------------------------------------------------------------
    !
    ! This routine uses pade approximants to continue the isotropic Eliashberg equations
    ! from the imaginary-axis to the real-axis
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE eliashbergcom, ONLY : nsw, ws, wsi, gap, delta, znorm, deltai, znormi
    USE constants_epw, ONLY : cone, ci, zero, czero
    USE utilities,     ONLY : pade_coeff, pade_eval
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: N
    !! Nr. of frequency points in Pade approx
    !
    ! Local variable
    CHARACTER(LEN = 256) :: cname
    !! character in file name
    !
    INTEGER :: iw
    !! Counter on frequency imag- and real-axis
    INTEGER :: ierr
    !! Error status
    !
    ! arrays used in pade_coeff and pade_eval
    COMPLEX(KIND = DP) :: omega
    !! frequency on real-axis
    COMPLEX(KIND = DP) :: padapp
    !! znorm or delta on real-axis after pade_eval
    COMPLEX(KIND = DP) :: a(N)
    !! a - pade coeff for deltai
    COMPLEX(KIND = DP) :: b(N)
    !! b - pade coeff for znormi
    COMPLEX(KIND = DP) :: z(N)
    !! z - frequency imag-axis
    COMPLEX(KIND = DP) :: u(N)
    !! u - deltai
    COMPLEX(KIND = DP) :: v(N)
    !! v - znormi
    !
    ALLOCATE(delta(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_iso', 'Error allocating delta', 1)
    ALLOCATE(znorm(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_iso', 'Error allocating znorm', 1)
    znorm(:) = czero
    delta(:) = czero
    a(:) = czero
    b(:) = czero
    z(:) = czero
    u(:) = czero
    v(:) = czero
    !
    DO iw = 1, N
      z(iw) = ci * wsi(iw)
      u(iw) = cone * deltai(iw)
      v(iw) = cone * znormi(iw)
    ENDDO
    !
    CALL pade_coeff(N, z, u, a)
    CALL pade_coeff(N, z, v, b)
    !
    DO iw = 1, nsw
      omega = cone * ws(iw)
      CALL pade_eval(N, z, a, omega, padapp)
      delta(iw) = padapp
      CALL pade_eval(N, z, b, omega, padapp)
      znorm(iw) = padapp
    ENDDO
    !
    WRITE(stdout, '(5x, a)') '   pade    Re[znorm]   Re[delta] [meV]'
    WRITE(stdout, '(5x, i6, 2ES15.6)') N, REAL(znorm(1)), REAL(delta(1)) * 1000.d0
!    WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6)') 'pade = ', N, &
!                  '   Re[znorm(1)] = ', REAL(znorm(1)), &
!                  '   Re[delta(1)] = ', REAL(delta(1))
    WRITE(stdout, '(5x, a, i6, a)') 'Convergence was reached for N = ', N, ' Pade approximants'
    WRITE(stdout, '(a)') ' '
    !
    cname = 'pade'
    CALL eliashberg_write_raxis(itemp, cname)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE pade_cont_iso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_iso_iaxis(itemp)
    !-----------------------------------------------------------------------
    !
    ! computes kernels K_{+}(n, n', T) and K_{-}(n, n', T)
    ! reference W. E. Pickett, PRB 26, 1186 (1982)
    !
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero
    USE constants,     ONLY : pi
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsiw, keri
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw, n
    !! Counter on frequency imag-axis
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: omega
    !! frequency imag-axis
    REAL(KIND = DP) :: lambda_eph
    !! electron-phonon coupling
    !
    ALLOCATE(keri(2 * nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_iaxis', 'Error allocating keri', 1)
    keri(:) = zero
    !
    DO iw = 1, 2 * nsiw(itemp)
      n = iw - 1
      omega = DBLE(2 * n) * pi * gtemp(itemp)
      CALL lambdar_iso(omega, lambda_eph)
      keri(iw) = lambda_eph
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_iso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_iso(omega, lambda_eph)
    !-----------------------------------------------------------------------
    !
    ! computes lambda(n - n')
    ! reference W. E. Pickett, PRB 26, 1186 (1982)
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nqstep
    USE constants_epw, ONLY : zero
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency on imag-axis
    REAL(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda(n-n')
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    lambda_eph = zero
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
      lambda_eph = lambda_eph + wsph(iwph) * a2f_iso(iwph) &
                 / (wsph(iwph) * wsph(iwph) + omega * omega)
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdar_iso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_iso_analytic_cont(itemp)
    !-----------------------------------------------------------------------
    !
    ! computes kernels K_{+}(w, iw_n, T) and K_{-}(w, iw_n, T)
    ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : muc
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, deltai, dsumi, zsumi
    USE constants_epw, ONLY : zero
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
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: esqrt
    !! 1 / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: kernelr
    !! 2 * Re[lambda(w - iw_n)]
    REAL(KIND = DP) :: kerneli
    !! 2 * Im[lambda(w - iw_n)]
    REAL(KIND = DP), ALLOCATABLE :: wesqrt(:)
    !! w / sqrt{w^2+\delta^2}
    REAL(KIND = DP), ALLOCATABLE :: desqrt(:)
    !! \delta / sqrt{w^2+\delta^2}
    !
    COMPLEX(KIND = DP) :: lambda_eph
    !! electron-phonon coupling lambda(w - iw_n)
    !
    ALLOCATE(wesqrt(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error allocating wesqrt', 1)
    ALLOCATE(desqrt(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error allocating desqrt', 1)
    ALLOCATE(dsumi(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error allocating dsumi', 1)
    ALLOCATE(zsumi(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error allocating zsumi', 1)
    wesqrt(:) = zero
    desqrt(:) = zero
    dsumi(:) = zero
    zsumi(:) = zero
    !
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nsiw(itemp) ! loop over iw_n
        CALL lambdai_iso(ws(iw), wsi(iwp), lambda_eph)
        kernelr = 2.d0 * REAL(lambda_eph)
        kerneli = 2.d0 * AIMAG(lambda_eph)
        IF (iw == 1) THEN
          esqrt = 1.d0 / DSQRT(wsi(iwp) * wsi(iwp) + deltai(iwp) * deltai(iwp))
          wesqrt(iwp) =  wsi(iwp) * esqrt
          desqrt(iwp) =  deltai(iwp) * esqrt
        ENDIF
        zsumi(iw) = zsumi(iw) + kerneli * wesqrt(iwp)
        dsumi(iw) = dsumi(iw) + (kernelr - 2.d0 * muc) * desqrt(iwp)
      ENDDO
    ENDDO
    !
    DEALLOCATE(wesqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error deallocating wesqrt', 1)
    DEALLOCATE(desqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error deallocating desqrt', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_iso_analytic_cont
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_iso(omega, omegap, lambda_eph)
    !-----------------------------------------------------------------------
    !
    ! computes lambda(w-iw_n)
    ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !
    USE kinds, ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
    USE constants_epw, ONLY : ci, czero
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency w at point iw on real-axis
    REAL(KIND = DP), INTENT(in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda(w - iw_n)
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    lambda_eph = czero
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
      lambda_eph = lambda_eph + wsph(iwph) * a2f_iso(iwph) &
                 / (wsph(iwph) * wsph(iwph) - (omega - ci * omegap) * (omega - ci * omegap))
    ENDDO ! iwph
    lambda_eph = lambda_eph * 2.d0 * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdai_iso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_iso_raxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the isotropic
    !! Eliashberg equations directly on the real-axis.
    !!
    !
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE epwcom,            ONLY : nsiter, nstemp, broyden_beta, broyden_ndim
    USE elph2,             ONLY : gtemp
    USE eliashbergcom,     ONLY : nsw, delta, deltap, gap
    USE constants_epw,     ONLY : kelvin2eV, ci, zero
    USE supercond,         ONLY : gen_freqgrid_raxis, eliashberg_grid
    USE utilities,         ONLY : mix_broyden
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
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: tcpu
    !! cpu time
    REAL(KIND = DP), EXTERNAL :: get_clock
    !! get the time spent
    REAL(KIND = DP), ALLOCATABLE :: rdeltain(:), rdeltaout(:)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: cdeltain(:), cdeltaout(:)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: df1(:, :), df2(:, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv1(:, :), dv2(:, :)
    !! Temporary variables for mix_broyden

    !
    CALL start_clock('iso_raxis')
    !
    CALL eliashberg_grid()
    CALL gen_freqgrid_raxis()
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      CALL prtheader_supercond(itemp, 4)
      !
      iter = 1
      conv = .FALSE.
      ALLOCATE(rdeltain(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating rdeltain', 1)
      ALLOCATE(cdeltain(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating cdeltain', 1)
      ALLOCATE(rdeltaout(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating rdeltaout', 1)
      ALLOCATE(cdeltaout(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating cdeltaout', 1)
      ALLOCATE(df1(nsw, broyden_ndim), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating df1', 1)
      ALLOCATE(dv1(nsw, broyden_ndim), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating dv1', 1)
      ALLOCATE(df2(nsw, broyden_ndim), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating df2', 1)
      ALLOCATE(dv2(nsw, broyden_ndim), STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error allocating dv2', 1)
      rdeltain(:)  = zero
      cdeltain(:)  = zero
      rdeltaout(:) = zero
      cdeltaout(:) = zero
      df1(:, :) = zero
      dv1(:, :) = zero
      df2(:, :) = zero
      dv2(:, :) = zero
      !
      iter = 1
      conv = .FALSE.
      DO WHILE(.NOT. conv .AND. iter <= nsiter)
        CALL integrate_eliashberg_iso_raxis(itemp, iter, conv)
        rdeltain(:)  = REAL(deltap(:))
        cdeltain(:)  = AIMAG(deltap(:))
        rdeltaout(:) = REAL(delta(:))
        cdeltaout(:) = AIMAG(delta(:))
        CALL mix_broyden(nsw, rdeltaout, rdeltain, broyden_beta, &
                         iter, broyden_ndim, conv, df1, dv1)
        CALL mix_broyden(nsw, cdeltaout, cdeltain, broyden_beta, &
                         iter, broyden_ndim, conv, df2, dv2)
        deltap(:) = rdeltain(:) + ci * cdeltain(:)
        iter = iter + 1
      ENDDO ! iter
      !
      DEALLOCATE(rdeltain, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating rdeltain', 1)
      DEALLOCATE(cdeltain, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating cdeltain', 1)
      DEALLOCATE(rdeltaout, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating rdeltaout', 1)
      DEALLOCATE(cdeltaout, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating cdeltaout', 1)
      DEALLOCATE(df1, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating df1', 1)
      DEALLOCATE(dv1, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating dv1', 1)
      DEALLOCATE(df2, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating df2', 1)
      DEALLOCATE(dv2, STAT = ierr)
      IF (ierr /= 0) CALL errore('eliashberg_iso_raxis', 'Error deallocating dv2', 1)
      !
      WRITE(stdout, '(5x, a, i3, a, f8.4, a, a, i3, a, f10.6, a, a, f10.6, a)') &
                    'temp(', itemp, ') = ', gtemp(itemp) / kelvin2eV, ' K ', &
                    '  gap_edge(', itemp, ') = ', gap(itemp), ' eV ', &
                    '  Re[delta(1)] = ', REAL(delta(1)), ' eV '
      WRITE(stdout, '(a)') '    '
      tcpu = get_clock('iso_raxis')
      WRITE( stdout, '(5x, a, i3, a, f8.1, a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
      !
      IF (conv) THEN
        WRITE(stdout, '(a)') '    '
        CALL print_clock('iso_raxis')
        WRITE(stdout, '(a)') '    '
      ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
        WRITE(stdout, '(a)') '    '
        CALL stop_clock('iso_raxis')
        CALL print_clock('iso_raxis')
        CALL errore('integrate_eliashberg_iso_raxis', 'converged was not reached', 1)
        RETURN
      ENDIF
      !
    ENDDO ! itemp
    !
    CALL stop_clock('iso_raxis')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_iso_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE integrate_eliashberg_iso_raxis(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the isotropic Eliashberg equations directly on the real-axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilker
    USE io_global,     ONLY : stdout
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : nswfc, nqstep, nsiter, muc, conv_thr_raxis, &
                              kerwrite, kerread, nstemp
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsw, ws, dws, gap0, gap, bewph, fdwp, &
                              kp, km, delta, deltap, znorm, wsph
    USE constants_epw, ONLY : kelvin2eV, ci, eps6, zero, czero
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! output file name
    CHARACTER(LEN = 256) :: cname
    !! character in output file name
    !
    INTEGER :: iw, iwp
    !! Counter on frequency real-axis
    INTEGER :: iwph
    !! Counter on frequency
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: a, b, c, d
    !! Temporary variables for reading kernelp and kernelm from file
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    REAL(KIND = DP), ALLOCATABLE :: wesqrt(:)
    !! w / sqrt{w^2+\delta^2}
    REAL(KIND = DP), ALLOCATABLE :: desqrt(:)
    !! \delta / sqrt{w^2+\delta^2}
    !
    COMPLEX(KIND = DP) :: esqrt
    !! 1 / sqrt{w^2+\delta^2}
    COMPLEX(KIND = DP) :: kernelp
    !! Temporary array for kernel K_{+}(w', w, T)
    COMPLEX(KIND = DP) :: kernelm
    !! Temporary array for kernel K_{-}(w', w, T)
    COMPLEX(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    ALLOCATE(wesqrt(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating wesqrt', 1)
    ALLOCATE(desqrt(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating desqrt', 1)
    wesqrt(:) = zero
    desqrt(:) = zero
    !
    IF (iter == 1) THEN
      ALLOCATE(gap(nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating gap', 1)
      ALLOCATE(delta(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating delta', 1)
      ALLOCATE(deltap(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating deltap', 1)
      ALLOCATE(znorm(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating znorm', 1)
      ALLOCATE(fdwp(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating fdwp', 1)
      ALLOCATE(kp(nsw, nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating kp', 1)
      ALLOCATE(km(nsw, nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating km', 1)
      gap(itemp) = zero
      deltap(:)  = czero
      deltap(:)  = gap0
      kp(:, :) = czero
      km(:, :) = czero
      !
      ! Fermi Dirac distribution
      fdwp(iw) = zero
      DO iw = 1, nsw
        IF (ABS(gtemp(itemp)) >  eps6) THEN
          fdwp(iw) = wgauss(-ws(iw) / gtemp(itemp), -99)
        ENDIF
      ENDDO
      !
      IF (kerwrite) THEN
        ! Bose-Einstein distribution
        ALLOCATE(bewph(nqstep), STAT = ierr)
        IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating bewph', 1)
        bewph(:) = zero
        DO iwph = 1, nqstep  ! loop over omega (integration variable)
          IF (ABS(gtemp(itemp)) > eps6) THEN
            bewph(iwph) = wgauss(-wsph(iwph) / gtemp(itemp), -99)
            bewph(iwph) = bewph(iwph) / (1.d0 - 2.d0 * bewph(iwph))
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    delta(:) = czero
    znorm(:) = czero
    !
    temp = gtemp(itemp) / kelvin2eV
    IF (temp < 10.d0) THEN
      WRITE(name1, '(a, a7, f4.2)') TRIM(prefix), '.ker_00', temp
    ELSEIF (temp >= 10.d0) THEN
      WRITE(name1, '(a, a6, f5.2)') TRIM(prefix), '.ker_0', temp
    ELSEIF (temp >= 100.d0) THEN
      WRITE(name1, '(a, a5, f6.2)') TRIM(prefix), '.ker_', temp
    ENDIF
    OPEN(UNIT = iufilker, FILE = name1, STATUS = 'unknown', FORM = 'unformatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'error opening file' // name1, iufilker)
    !
    IF (iter == 1) THEN
      ALLOCATE(deltaold(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating deltaold', 1)
      deltaold(:) = gap0
    ENDIF
    !
    absdelta = zero
    reldelta = zero
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nsw ! loop over omega_prime
        IF (iter == 1) THEN
          !
          ! read the kernels from file if they were calculated before otherwise calculate them
          IF (kerread) THEN
            READ(iufilker) a, b, c, d
            kp(iwp, iw) = a + ci * b
            km(iwp, iw) = c + ci * d
          ENDIF
          IF (kerwrite) THEN
            CALL kernel_raxis(iw, iwp, kernelp, kernelm)
            kp(iwp, iw) = kernelp
            km(iwp, iw) = kernelm
            WRITE(iufilker) REAL(kp(iwp, iw)), AIMAG(kp(iwp, iw)), &
                            REAL(km(iwp, iw)), AIMAG(km(iwp, iw))
          ENDIF
        ENDIF
        !
        ! this step is performed at each iter step only for iw=1 since it is independent of w(iw)
        IF (iw == 1) THEN
          esqrt = 1.d0 / SQRT(ws(iwp) * ws(iwp) - deltap(iwp) * deltap(iwp))
          wesqrt(iwp) =  REAL(ws(iwp) * esqrt)
          desqrt(iwp) =  REAL(deltap(iwp) * esqrt)
        ENDIF
        !
        znorm(iw) = znorm(iw) + dws(iwp) * wesqrt(iwp) * km(iwp, iw)
        delta(iw) = delta(iw) + dws(iwp) * desqrt(iwp) &
                  * (kp(iwp, iw) - muc * (1.d0 - 2.d0 * fdwp(iwp)))
      ENDDO ! iwp
      znorm(iw) = 1.d0 - znorm(iw) / ws(iw)
      delta(iw) = delta(iw) / znorm(iw)
      reldelta = reldelta + ABS(delta(iw) - deltaold(iw))
      absdelta = absdelta + ABS(delta(iw))
    ENDDO ! iw
    CLOSE(iufilker)
    errdelta = reldelta / absdelta
    deltaold(:) = delta(:)
    !
    IF (iter == 1) &
      WRITE(stdout, '(5x, a)') '   iter      ethr         Re[znorm]   Re[delta] [eV]'
    WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, REAL(znorm(1)), REAL(delta(1))
!    WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
!                  '   ethr = ', errdelta, '   Re[znorm(1)] = ', REAL(znorm(1)), &
!                  '   Re[delta(1)] = ', REAL(delta(1))
    !
    IF (errdelta < conv_thr_raxis) conv = .TRUE.
    IF (errdelta < conv_thr_raxis .OR. iter == nsiter) THEN
      cname = 'real'
      CALL eliashberg_write_raxis(itemp, cname)
      gap0 = gap(itemp)
    ENDIF
    !
    DEALLOCATE(wesqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error deallocating wesqrt', 1)
    DEALLOCATE(desqrt, STAT = ierr)
    IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error deallocating desqrt', 1)
    !
    IF (conv .OR. iter == nsiter) THEN
      DEALLOCATE(deltaold, STAT = ierr)
      IF (ierr /= 0) CALL errore('integrate_eliashberg_iso_raxis', 'Error allocating deltaold', 1)
    ENDIF
    IF (conv) THEN
      WRITE(stdout, '(5x, a, i6)') 'Convergence was reached in nsiter = ', iter
      WRITE(stdout,'(a)') ' '
    ELSEIF (.NOT. conv .AND. iter == nsiter) THEN
      WRITE(stdout, '(5x, a, i6)') 'Convergence was not reached in nsiter = ', iter
      WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_raxis'
      WRITE(stdout,'(a)') ' '
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE integrate_eliashberg_iso_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_raxis(iw, iwp, kernelp, kernelm)
    !-----------------------------------------------------------------------
    !
    ! computes kernels K_{+}(w', w, T) and K_{-}(w', w, T)
    ! reference M. J. Holcomb, PRB 54, 6648 (1996)
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : ci, eps6, zero, czero, one
    USE constants,     ONLY : pi
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph, ws, bewph, fdwp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iw
    !! index frequency w : ws(iw)
    INTEGER, INTENT(in) :: iwp
    !! index frequency w' : ws(iwp)
    COMPLEX(KIND = DP), INTENT(out) :: kernelp
    !! phonon kernel K_{+}(w', w, T)
    COMPLEX(KIND = DP), INTENT(out) :: kernelm
    !! phonon kernel K_{-}(w', w, T)
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: degaussw0
    !! smearing
    REAL(KIND = DP) :: inv_degaussw0
    !! define inverse smearing for efficiency
    REAL(KIND = DP) :: f1, f2, f3, f4, w1, w2, w3, w4, var1, var2
    !! Temporaty variables
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !
    COMPLEX(KIND = DP) :: e1, e2, e3, e4, g1, g2, g3, g4
    !! Temporary variables
    !
    degaussw0 = dwsph
    inv_degaussw0 = one / degaussw0
    !
    f1 = zero
    f2 = zero
    f3 = zero
    f4 = zero
    var1 = zero
    var2 = zero
    kernelp = czero
    kernelm = czero
    e1 = czero
    e2 = czero
    e3 = czero
    e4 = czero
    g1 = czero
    g2 = czero
    g3 = czero
    g4 = czero
    !
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
      !
      ! a small complex number is added to denominator to move the pole away from the real-axis
      !
      ! in order to reduce the numerical noise at very small frequencies coming from
      ! the complex number added in the denominator, the contribution of the imaginary part
      ! is reestimated using delta function (RM notes)
      !
      ! subtract the imaginary part coming from e1 to e4 and add instead the imaginary part
      ! coming from f1 to f4
      !
      w1 = wsph(iwph) + ws(iwp) + ws(iw)
      w2 = wsph(iwph) + ws(iwp) - ws(iw)
      w3 = wsph(iwph) - ws(iwp) + ws(iw)
      w4 = wsph(iwph) - ws(iwp) - ws(iw)
      !
      e1 = one / (w1 + ci * degaussw0)
      e2 = one / (w2 - ci * degaussw0)
      e3 = one / (w3 + ci * degaussw0)
      e4 = one / (w4 - ci * degaussw0)
      !
      ! estimate of the imaginary part using delta function
      f1 = w0gauss(w1 * inv_degaussw0, 0) * inv_degaussw0
      f2 = w0gauss(w2 * inv_degaussw0, 0) * inv_degaussw0
      f3 = w0gauss(w3 * inv_degaussw0, 0) * inv_degaussw0
      f4 = w0gauss(w4 * inv_degaussw0, 0) * inv_degaussw0
      !
      g1 = e1 - ci * AIMAG(e1) - ci * pi * f1
      g2 = e2 - ci * AIMAG(e2) + ci * pi * f2
      g3 = e3 - ci * AIMAG(e3) - ci * pi * f3
      g4 = e4 - ci * AIMAG(e4) + ci * pi * f4
      var1 = one - fdwp(iwp) + bewph(iwph)
      var2 = fdwp(iwp) + bewph(iwph)
      kernelp = kernelp + a2f_iso(iwph) * (var1 * (g1 + g2) - var2 * (g3 + g4))
      kernelm = kernelm + a2f_iso(iwph) * (var1 * (g1 - g2) + var2 * (g3 - g4))
    ENDDO ! iwph
    kernelp = kernelp * dwsph
    kernelm = kernelm * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_raxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE crit_temp_iso()
    !-----------------------------------------------------------------------
    !!
    !! SH: Routine to find the Tc by solving the linearized Eliashberg
    !!        (isotropic) equation for the specified range of temperatures
    !! HP: updated 5/11/2021
    !
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE epwcom,            ONLY : limag, nstemp, muc, tc_linear, tc_linear_solver
    USE elph2,             ONLY : gtemp
    USE eliashbergcom,     ONLY : nsiw, wsi
    USE constants_epw,     ONLY : Kelvin2eV, zero, pi
    USE supercond,         ONLY : gen_freqgrid_iaxis, eliashberg_grid, &
                                  crit_temp_solver
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iw, iwp, iws
    !! Counter for loops
    INTEGER :: ierr
    !! Error status
    INTEGER :: N
    !! Dimension of the main matrix = 2 * nsiw(itemp) - 1
    INTEGER :: numiter
    !! Number of iterations
    !
    REAL(KIND = DP) :: lambda_eph
    !! Electron-phonon coupling
    REAL(KIND = DP), ALLOCATABLE :: S(:, :)
    !! Main matrix for finding eigenvalues
    INTEGER, ALLOCATABLE :: freqn(:)
    !! Index of Matsubara frequencies
    REAL(KIND = DP), ALLOCATABLE :: freqv(:)
    !! Value of Matsubara frequencies
    REAL(KIND = DP), ALLOCATABLE :: freqk(:, :)
    !! Isotropic lambda for (iw,iwp)
    REAL(KIND = DP) :: maxeigen
    !! Maximum eigenvalue for each itemp
    REAL(KIND = DP) :: selement
    !! Used for summation of matrix elements
    !
    WRITE(stdout, '(/5x, a, a6/)') &
      'Start: Solving (isotropic) linearized Eliashberg equation with solver = ', tc_linear_solver
    !
    CALL eliashberg_grid()
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, a, f7.2, a)') 'For the first Temp. ', gtemp(itemp) / Kelvin2eV, '  K'
        WRITE(stdout, '(7x, a, i6, a, i6)') 'Total number of frequency points nsiw(', itemp, ') = ', nsiw(itemp)
        WRITE(stdout, '(7x, a, f10.4/)') 'Cutoff frequency wscut = ', (2.d0 * nsiw(itemp) + 1) * pi * gtemp(itemp)
        WRITE(stdout, '(5x, a)') &
           'Superconducting transition temp. Tc is the one which has Max. eigenvalue close to 1'
      ENDIF
      CALL gen_freqgrid_iaxis(itemp)
      !
      ! n = dimension of matrix: # of symmetric Matsubara frequencies
      n = 2 * nsiw(itemp) - 1
      !
      ALLOCATE(s(n, n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error allocating s', 1)
      ALLOCATE(freqn(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error allocating freqn', 1)
      ALLOCATE(freqv(n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error allocating freqv', 1)
      ALLOCATE(freqk(n, n), STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error allocating freqk', 1)
      !
      selement    = zero
      s(:, :)     = zero
      freqn(:)    = 0
      freqv(:)    = zero
      freqk(:, :) = zero
      !
      ! Re-construct full range of frequency indices
      iwp = -nsiw(itemp)
      DO iw = 1, n
        freqn(iw) = iwp + iw
      ENDDO
      ! Re-construct full range of frequency vlaues
      DO iw = 1, n
        freqv(iw) = DBLE(2 * freqn(iw) + 1) * pi * gtemp(itemp)
      ENDDO
      !
      ! Re-construct full isotropic kernel
      DO iw = 1, n
        DO iwp = 1, n
          CALL lambdar_iso(freqv(iw) - freqv(iwp), lambda_eph)
          freqk(iw, iwp) = lambda_eph
        ENDDO
      ENDDO
      !
      ! Find the input S matrix components (Allen et al., Eqns. 11.10-11)
      DO iw = 1, n
        DO iwp = 1, n
          selement = freqk(iw, iwp) - muc
          IF (iw == iwp) THEN
            DO iws = 1, n
              selement = selement - (freqk(iw, iws) * (freqv(iw) * freqv(iws)) &
                / ABS(freqv(iw) * freqv(iws)))
            ENDDO
          ENDIF
          s(iw, iwp) = ABS( 1.d0 / DBLE(2 * freqn(iwp) + 1)) * selement
        ENDDO
      ENDDO
      !
      ! Main solver
      CALL crit_temp_solver(n, s, maxeigen, numiter)
      !
      ! Print the output: largest eigenvalue
      IF (itemp == 1) WRITE(stdout, '(6x, a)') REPEAT('-', 65)
      IF (itemp == 1) WRITE(stdout, 102) 'Temp.', 'Max.', 'nsiw', 'wscut', 'Nr. of iters'
      IF (itemp == 1) WRITE(stdout, 103) '(K)', 'eigenvalue', '(itemp)', '(eV)', 'to Converge'
      IF (itemp == 1) WRITE(stdout, '(6x, a)') REPEAT('-', 65)
      WRITE(stdout, 104) gtemp(itemp) / Kelvin2eV, maxeigen, nsiw(itemp), &
                         (2.d0 * nsiw(itemp) + 1) * pi * gtemp(itemp), numiter
      102 FORMAT (7x, a11, a11, a13, a10, a16)
      103 FORMAT (7x, a9, a17, a11, a7, a17)
      104 FORMAT(5x, f12.2, f15.7, i9, f12.4, i9)
      IF (itemp == nstemp) WRITE(stdout, '(6x, a)') REPEAT('-', 65)
      !
      DEALLOCATE(s, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error deallocating A', 1)
      DEALLOCATE(freqn, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error deallocating freqn', 1)
      DEALLOCATE(freqv, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error deallocating freqv', 1)
      DEALLOCATE(freqk, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error deallocating freqk', 1)
      DEALLOCATE(wsi, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error deallocating wsi', 1)
      !
    ENDDO ! itemp
    !
    CALL deallocate_iso()
    !
    WRITE(stdout, '(5x, a/)') 'Finish: Solving (isotropic) linearized Eliashberg equation'
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE crit_temp_iso
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_iso_iaxis()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated for imag-axis solutions
    !!
    !----------------------------------------------------------------------
    !
    USE eliashbergcom, ONLY : wsi, gap, deltai, deltaip, znormi, nznormi
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! gen_freqgrid_iaxis
    DEALLOCATE(wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating wsi', 1)
    ! sum_eliashberg_iso_iaxis
    DEALLOCATE(gap, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating gap', 1)
    DEALLOCATE(deltai, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating deltai', 1)
    DEALLOCATE(deltaip, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating deltaip', 1)
    DEALLOCATE(znormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating znormi', 1)
    DEALLOCATE(nznormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating nznormi', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_iso_iaxis
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_iso_raxis()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated for real-axis solutions
    !!
    USE epwcom,        ONLY : limag, lacon, lreal
    USE eliashbergcom, ONLY : ws, delta, deltap, znorm, znormp, &
                              gp, gm, dsumi, zsumi, &
                              dws, fdwp, bewph, kp, km
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! gen_freqgrid_iaxis or gen_freqgrid_raxis
    DEALLOCATE(ws, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating ws', 1)
    !
    ! pade_cont_iso or analytic_cont_iso or integrate_eliashberg_iso_raxis
    DEALLOCATE(delta, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating delta', 1)
    DEALLOCATE(znorm, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating znorm', 1)
    !
    ! analytic_cont_iso   
    IF (limag .AND. lacon) THEN
      DEALLOCATE(deltap, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating deltap', 1)
      DEALLOCATE(znormp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating znormp', 1)
      DEALLOCATE(gp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating gp', 1)
      DEALLOCATE(gm, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating gm', 1)
      ! kernel_iso_analytic_cont
      DEALLOCATE(dsumi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating dsumi', 1)
      DEALLOCATE(zsumi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating zsumi', 1)
    ENDIF
    !
    ! integrate_eliashberg_iso_raxis
    IF (lreal) THEN
      DEALLOCATE(dws, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating dws', 1)
      DEALLOCATE(fdwp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating fdwp', 1)
      DEALLOCATE(bewph, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating bewph', 1)
      DEALLOCATE(kp, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating kp', 1)
      DEALLOCATE(km, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating km', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_iso_raxis
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_iso()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by eliashberg_init, read_a2f, 
    !!  eliashberg_grid
    !!
    USE epwcom,        ONLY : limag
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : a2f_iso, wsph, nsiw
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! eliashberg_init
    DEALLOCATE(gtemp, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating gtemp', 1)
    ! read_a2f
    DEALLOCATE(wsph, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating wsph', 1)
    DEALLOCATE(a2f_iso, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating a2f_iso', 1)
    ! eliashberg_grid
    IF (limag) THEN
      DEALLOCATE(nsiw, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating nsiw', 1)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_iso
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE supercond_iso
  !-----------------------------------------------------------------------
