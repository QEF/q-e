  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
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
    !! SH: Modified to allow for fbw calculations (Nov 2021).
    !
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE input,             ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                                  limag, lpade, lacon, npade, fbw
    USE supercond_common,  ONLY : nsw, nsiw, deltai, deltaip, delta, deltap, &
                                  znormi, znormip, shifti, shiftip
    USE ep_constants,      ONLY : kelvin2eV, ci, zero
    USE supercond,         ONLY : dos_quasiparticle, gen_freqgrid_iaxis, &
                                  eliashberg_grid
    USE utilities,         ONLY : mix_wrap
    USE printing,          ONLY : prtheader_supercond
    USE io_supercond,      ONLY : read_dos
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
    REAL(KIND = DP), ALLOCATABLE :: df1(:, :), df2(:, :), df3(:, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv1(:, :), dv2(:, :), dv3(:, :)
    !! Temporary variables for mix_broyden
    !
    CALL start_clock('iso_iaxis')
    !
    CALL eliashberg_grid()
    !
    ! SH: read electronic dos (only for fbw calculations)
    IF (fbw) CALL read_dos()
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      CALL start_clock('iaxis_imag')
      CALL gen_freqgrid_iaxis(itemp)
      CALL prtheader_supercond(itemp, 1)
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
        IF (fbw) THEN
          ALLOCATE(df2(nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating df2', 1)
          ALLOCATE(dv2(nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating dv2', 1)
          ALLOCATE(df3(nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating df3', 1)
          ALLOCATE(dv3(nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error allocating dv3', 1)
          df2(:, :) = zero
          df3(:, :) = zero
          dv2(:, :) = zero
          dv3(:, :) = zero
        ENDIF
        !
        iter = 1
        conv = .FALSE.
        !
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL sum_eliashberg_iso_iaxis(itemp, iter, conv)
          ! SH: mix_wrap is used to allow for linear mixing
          CALL mix_wrap(nsiw(itemp), deltai, deltaip, broyden_beta, &
                        iter, broyden_ndim, conv, df1, dv1)
          ! SH: mixing for fbw runs;
          IF (fbw) THEN
            CALL mix_wrap(nsiw(itemp), znormi, znormip, broyden_beta, &
                          iter, broyden_ndim, conv, df2, dv2)
            CALL mix_wrap(nsiw(itemp), shifti, shiftip, broyden_beta, &
                          iter, broyden_ndim, conv, df3, dv3)
          ENDIF
          iter = iter + 1
        ENDDO ! iter
        !
        DEALLOCATE(df1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating df1', 1)
        DEALLOCATE(dv1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating dv1', 1)
        ! SH: deallocaating fbw-related arrays
        IF (fbw) THEN
          DEALLOCATE(df2, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating df2', 1)
          DEALLOCATE(dv2, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating dv2', 1)
          DEALLOCATE(df3, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating df3', 1)
          DEALLOCATE(dv3, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_iso_iaxis', 'Error deallocating dv3', 1)
        ENDIF
        !
        IF (conv) THEN
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
      ! SH: acon is not implemented for fbw runs
      IF (.NOT. fbw .AND. lacon) THEN
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
          ! SH: mix_wrap is used to allow for linear mixing
          CALL mix_wrap(nsw, rdeltaout, rdeltain, broyden_beta, &
                        iter, broyden_ndim, conv, df1, dv1)
          CALL mix_wrap(nsw, cdeltaout, cdeltain, broyden_beta, &
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
    !! SH: Modified to allow for fbw calculations (Nov 2021).
    !! SM: Added paralleilization of freq. among images (April 2025).
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE input,         ONLY : nsiter, muc, conv_thr_iaxis, &
                              fbw, dos_del, muchem
    USE global_var,    ONLY : gtemp
    USE mp_global,     ONLY : inter_image_comm, my_image_id
    USE mp,            ONLY : mp_sum
    USE supercond_common,     ONLY : nsiw, gap0, wsi, wsn, keri, muintr, &
                              deltai, deltaip, znormi, nznormi, &
                              znormip, shifti, shiftip, ef0, &
                              en, dosen, ndos, dosef
    USE ep_constants,  ONLY : kelvin2eV, zero, one
    USE ep_constants,  ONLY : pi
    USE io_supercond,  ONLY : eliashberg_write_iaxis
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
    INTEGER :: ie
    !! Counter on energy grid
    INTEGER :: startiw, lastiw
    !! Lower/upper bound index for iw-points
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
    REAL(KIND = DP) :: esqrt, zesqrt, desqrt, sesqrt
    !! Temporary variables
    REAL(KIND = DP), SAVE :: nel, nstate
    !! mu_inter parameters
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: omega
    !! sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: dFE
    !! free energy difference between supercond and normal states
    REAL(KIND = DP), ALLOCATABLE, SAVE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    REAL(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    ! SM: added image parallelization for distributing freq. points into images
    CALL divide(inter_image_comm, nsiw(itemp), startiw, lastiw) 
    !
    IF (iter == 1) THEN
      ! Print for user information
      WRITE(stdout, '(5x, "   startiw = ", i0, ", lastiw = ", i0, ", nsiw(itemp) = ", i0)') &
      startiw, lastiw, nsiw(itemp)
      !
      ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating inv_wsi', 1)
      !
      inv_wsi(:) = zero
      DO iw = 1, nsiw(itemp)
        inv_wsi(iw) = one / wsi(iw)
      ENDDO
      !
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
      deltaip(:)  = gap0
      deltaold(:) = gap0
      !
      IF (fbw) THEN
        ALLOCATE(znormip(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating znormip', 1)
        ALLOCATE(shifti(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating shifti', 1)
        ALLOCATE(shiftip(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error allocating shiftip', 1)
        znormip(:)  = one
        shiftip(:)  = zero
      ENDIF
      !
      CALL kernel_iso_iaxis(itemp)
      !
      IF (fbw .AND. itemp == 1) THEN
        nel = zero
        nstate = zero
        DO ie = 1, ndos
          ! HP: one of the 2 factor comes form the spin in DOS
          nstate = nstate + dos_del * dosen(ie)
          IF ((en(ie) - ef0) < zero) nel = nel + 2.d0 * dos_del * dosen(ie)
        ENDDO
      ENDIF
      !
      muintr = ef0
    ENDIF !iter
    !
    IF (fbw) THEN
      ! SH: update the chemical potential from the inital guess
      IF (muchem) CALL mu_inter_iso(itemp, muintr, nel, nstate)
      !
      inv_dos = one / dosef
      !
      nznormi(:) = zero
      deltai(:)  = zero
      znormi(:)  = zero
      shifti(:)  = zero
      DO iwp = 1, nsiw(itemp) ! loop over omega_prime
        zesqrt = zero
        desqrt = zero
        sesqrt = zero
        DO ie = 1, ndos
          esqrt = dos_del / (znormip(iwp)**2.d0 * (wsi(iwp)**2.d0 + deltaip(iwp)**2.d0) + &
            (en(ie) - muintr + shiftip(iwp))**2.d0)
          zesqrt = zesqrt + dosen(ie) * esqrt * wsi(iwp) * znormip(iwp)
          desqrt = desqrt + dosen(ie) * esqrt * deltaip(iwp) * znormip(iwp)
          sesqrt = sesqrt + dosen(ie) * esqrt * (en(ie)- muintr + shiftip(iwp))
        ENDDO
        zesqrt = zesqrt * inv_dos
        desqrt = desqrt * inv_dos
        sesqrt = sesqrt * inv_dos
        !
        !DO iw = 1, nsiw(itemp) ! loop over omega
        ! SM: distribute into images
        DO iw = startiw, lastiw
          ! SH: For general case (including sparse sampling)
          !       "actual" matsubara indices n1/n2 are needed instead of iw/iwp
          lambdam = keri(ABS(wsn(iw) - wsn(iwp)) + 1)
          lambdap = keri(ABS(wsn(iw) + wsn(iwp) + 1) + 1)
          ! Eq. (4.4) in Picket, PRB 26, 1186 (1982)
          kernelm = lambdam - lambdap
          kernelp = lambdam + lambdap
          nznormi(iw) = nznormi(iw) + kernelm
          ! Eqs. (4.1-4.3) in Picket, PRB 26, 1186 (1982) for FBW
          ! using kernelm and kernelp the sum over |wp| < wscut
          ! is rewritten as a sum over iwp = 1, nsiw(itemp)
          znormi(iw) = znormi(iw) + zesqrt * kernelm
          deltai(iw) = deltai(iw) + desqrt * (kernelp - 2.d0 * muc)
          shifti(iw) = shifti(iw) + sesqrt * kernelp
        ENDDO ! iw
      ENDDO ! iwp
      !
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsiw(itemp) ! loop over omega
        nznormi(iw) = 1.d0 + gtemp(itemp) * nznormi(iw) * inv_wsi(iw)
        ! Eqs.(34)-(35) in Margine and Giustino, PRB 87, 024505 (2013)
        znormi(iw) = 1.d0 + gtemp(itemp) * znormi(iw) * inv_wsi(iw)
        deltai(iw) =   gtemp(itemp) * deltai(iw) / znormi(iw)
        shifti(iw) = - gtemp(itemp) * shifti(iw)
        reldelta   = reldelta + ABS(deltai(iw) - deltaold(iw))
        absdelta   = absdelta + ABS(deltai(iw))
      ENDDO ! iw
      ! Collect from images
      CALL mp_sum(nznormi, inter_image_comm)
      CALL mp_sum(znormi, inter_image_comm)
      CALL mp_sum(deltai, inter_image_comm)
      CALL mp_sum(shifti, inter_image_comm)
      CALL mp_sum(reldelta, inter_image_comm)
      CALL mp_sum(absdelta, inter_image_comm)
    ELSE !FSR
      deltai(:)  = zero
      znormi(:)  = zero
      nznormi(:) = zero
      DO iwp = 1, nsiw(itemp) ! loop over omega_prime
        esqrt = 1.d0 / DSQRT(wsi(iwp)**2.d0 + deltaip(iwp)**2.d0)
        zesqrt = wsi(iwp) * esqrt
        desqrt = deltaip(iwp) * esqrt
        !
        DO iw = 1, nsiw(itemp) ! loop over omega
          ! SH: For general case (including sparse sampling)
          !       "actual" matsubara indices n1/n2 are needed instead of iw/iwp
          lambdam = keri(ABS(wsn(iw) - wsn(iwp)) + 1)
          lambdap = keri(ABS(wsn(iw) + wsn(iwp) + 1) + 1)
          ! Eq. (4.4) in Picket, PRB 26, 1186 (1982)
          kernelm = lambdam - lambdap
          kernelp = lambdam + lambdap
          nznormi(iw) = nznormi(iw) + kernelm
          ! Eqs.(34)-(35) in Margine and Giustino, PRB 87, 024505 (2013)
          ! using kernelm and kernelp the sum over |wp| < wscut in Eqs. (34)-(35)
          ! is rewritten as a sum over iwp = 1, nsiw(itemp)
          znormi(iw) = znormi(iw) + zesqrt * kernelm
          deltai(iw) = deltai(iw) + desqrt * (kernelp - 2.d0 * muc)
        ENDDO ! iw
      ENDDO ! iwp
      !
      absdelta = zero
      reldelta = zero
      !DO iw = 1, nsiw(itemp) ! loop over omega
      !SM: distribute into images
      DO iw = startiw, lastiw
        znormi(iw)  = 1.d0 + pi * gtemp(itemp) * znormi(iw) * inv_wsi(iw)
        ! Eqs.(34)-(35) in Margine and Giustino, PRB 87, 024505 (2013)
        nznormi(iw) = 1.d0 + pi * gtemp(itemp) * nznormi(iw) * inv_wsi(iw)
        deltai(iw)  = pi * gtemp(itemp) * deltai(iw) / znormi(iw)
        reldelta = reldelta + ABS(deltai(iw) - deltaold(iw))
        absdelta = absdelta + ABS(deltai(iw))
      ENDDO ! iw
      ! Collect from images
      CALL mp_sum(nznormi, inter_image_comm)
      CALL mp_sum(znormi, inter_image_comm)
      CALL mp_sum(deltai, inter_image_comm)
      CALL mp_sum(reldelta, inter_image_comm)
      CALL mp_sum(absdelta, inter_image_comm)
    ENDIF ! fbw
    errdelta = reldelta / absdelta
    deltaold(:) = deltai(:)
    !
    IF (iter == 1 .AND. fbw) WRITE(stdout, '(5x, a)') &
      '   iter      ethr          znormi      deltai [meV]   shifti [meV]    mu (eV)'
    IF (iter == 1 .AND. .NOT. fbw) &
      WRITE(stdout, '(5x, a)') '   iter      ethr          znormi      deltai [meV]'
    IF (fbw) THEN
      WRITE(stdout, '(5x, i6, 5ES15.6)') &
        iter, errdelta, znormi(1), deltai(1) * 1000.d0, shifti(1) * 1000.d0, muintr
    ELSE
      WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, znormi(1), deltai(1) * 1000.d0
      !WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
      !            '   ethr = ', errdelta, '   znormi(1) = ', znormi(1), &
      !            '   deltai(1) = ', deltai(1)
    ENDIF
    !
    IF (errdelta < conv_thr_iaxis) conv = .TRUE.
    IF (errdelta < conv_thr_iaxis .OR. iter == nsiter) THEN
      gap0 = deltai(1)
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
    IF (conv .OR. iter == nsiter) THEN
      IF (fbw) THEN
        WRITE(stdout, '(5x, a, i3, a, ES20.10, a)') &
                    'Chemical potential (itemp = ', itemp, ') = ', muintr, ' eV'
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      ! Compute the free energy difference between the superconducting and normal states
      dFE = zero
      DO iw = 1, nsiw(itemp)
        omega = DSQRT(wsi(iw) * wsi(iw) + deltai(iw) * deltai(iw))
        dFE = dFE - (omega - wsi(iw)) &
            * (znormi(iw) - nznormi(iw) * wsi(iw) / omega)
      ENDDO
      dFE = dFE * pi * gtemp(itemp)
      ! HM: We multiply it by two to account for the contribution 
      ! from the negative Matsubara frequencies.
      dFE = dFE * 2.0d0
      !
      WRITE(stdout, '(5x, a, i3, a, f8.3, a, a, f12.6, a)') &
              'Temp (itemp = ', itemp, ') = ', gtemp(itemp) / kelvin2eV, ' K', &
              '  Free energy = ', dFE * 1000.0, ' meV'
      WRITE(stdout, '(a)') ' '
      !
      DEALLOCATE(inv_wsi, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_iso_iaxis', 'Error deallocating inv_wsi', 1)
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
    USE input,         ONLY : nqstep, nsiter, conv_thr_racon, lpade
    USE global_var,    ONLY : gtemp
    USE supercond_common,     ONLY : nsw, dwsph, ws, gap0, a2f_tmp, dsumi, zsumi, &
                              delta, deltap, znorm, znormp, gp, gm
    USE ep_constants,  ONLY : ci, zero, czero, cone
    USE ep_constants,  ONLY : pi
    USE io_supercond,  ONLY : eliashberg_write_raxis
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
    REAL(KIND = DP) :: root_im
    !! Temporary variable
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
        deltap(:)   = gap0
        deltaold(:) = gap0
      ENDIF
      !
      !! Eq.(28) in Margine and Giustino, PRB 87, 024505 (2013)
      DO iwp = 1, nqstep ! loop over omega_prime
        DO iw = 1, nsw ! loop over omega
          CALL gamma_acont(ws(iw), ws(iwp), gtemp(itemp), rgammap, rgammam)
          gp(iw, iwp) = rgammap
          gm(iw, iwp) = rgammam
        ENDDO
      ENDDO
      CALL kernel_iso_analytic_cont(itemp)
    ENDIF ! iter
    !
    ! RM - calculate esqrt (esqrt is stored as znormp to avoid allocation of a new array)
    ! Z(w \pm w') / sqrt{[Z(w \pm w')]^2 * [(w \pm w')^2] - [D(w \pm w')]^2]}
    ! nvfortran bug workaround https://gitlab.com/QEF/q-e/-/issues/593
#if defined(__NVCOMPILER) && ( __NVCOMPILER_MAJOR__ > 23 || ( __NVCOMPILER_MAJOR__ == 23 &&  __NVCOMPILER_MINOR__ >= 3) )
    !pgi$l novector
#endif
    DO iw = 1, nsw
      root = SQRT(znormp(iw) * znormp(iw) * (ws(iw) * ws(iw) - deltap(iw) * deltap(iw)))
      root_im=AIMAG(root)
      IF (root_im < zero) &
        root = CONJG(root)
      esqrt = znormp(iw) / root
      znormp(iw) = esqrt
    ENDDO
    !
    delta(:) = czero
    znorm(:) = czero
    DO iwp = 1, nqstep ! loop over omega_prime
      DO iw = 1, nsw ! loop over omega
        !
        i = iw + iwp
        IF (i <= nsw) THEN
          esqrt = gp(iw, iwp) * a2f_tmp(iwp) * znormp(i)
          znorm(iw) = znorm(iw) - esqrt * ws(i)
          delta(iw) = delta(iw) - esqrt * deltap(i)
        ENDIF
        !
        i = ABS(iw - iwp)
        IF (i > 0) THEN
          esqrt = gm(iw, iwp) * a2f_tmp(iwp) * znormp(i)
          znorm(iw) = znorm(iw) + esqrt * ws(i) * SIGN(1, iw - iwp)
          delta(iw) = delta(iw) + esqrt * deltap(i)
        ENDIF
      ENDDO ! iw
    ENDDO ! iwp
    !
    absdelta = zero
    reldelta = zero
    DO iw = 1, nsw ! loop over omega
      znorm(iw) = 1.d0 + pi * (- gtemp(itemp) * zsumi(iw) + ci * znorm(iw) * dwsph) / ws(iw)
      delta(iw) = pi * (gtemp(itemp) * dsumi(iw) + ci * delta(iw) * dwsph) / znorm(iw)
      reldelta = reldelta + ABS(delta(iw) - deltaold(iw))
      absdelta = absdelta + ABS(delta(iw))
    ENDDO ! iw
    errdelta = reldelta / absdelta
    deltaold(:) = delta(:)
    !
    !RM - update znormp for next iteration
    znormp(:) = znorm(:)
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
    !!
    !! This routine uses pade approximants to continue the isotropic Eliashberg equations
    !! from the imaginary-axis to the real-axis
    !!
    !! SH: Modified to allow for fbw calculations (Nov 2021).
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE input,         ONLY : fbw
    USE supercond_common,     ONLY : nsw, ws, wsi, delta, znorm, deltai, znormi, shift, shifti
    USE ep_constants,  ONLY : cone, ci, zero, czero
    USE utilities,     ONLY : pade_coeff, pade_eval
    USE io_supercond,  ONLY : eliashberg_write_raxis
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
    COMPLEX(KIND = DP) :: c(N)
    !! b - pade coeff for shifti
    COMPLEX(KIND = DP) :: z(N)
    !! z - frequency imag-axis
    COMPLEX(KIND = DP) :: u(N)
    !! u - deltai
    COMPLEX(KIND = DP) :: v(N)
    !! v - znormi
    COMPLEX(KIND = DP) :: w(N)
    !! w - shifti
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
    IF (fbw) THEN
      ALLOCATE(shift(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_iso', 'Error allocating shift', 1)
      shift(:) = czero
      c(:) = czero
      w(:) = czero
    ENDIF
    !
    IF (fbw) THEN
      ! SH: for fbw case of the Pade approximation
      DO iw = 1, N
        z(iw) = ci * wsi(iw)
        u(iw) = cone * deltai(iw)
        v(iw) = cone * znormi(iw)
        w(iw) = cone * shifti(iw)
      ENDDO
      !
      CALL pade_coeff(N, z, u, a)
      CALL pade_coeff(N, z, v, b)
      CALL pade_coeff(N, z, w, c)
      !
      DO iw = 1, nsw
        omega = cone * ws(iw)
        CALL pade_eval(N, z, a, omega, padapp)
        delta(iw) = padapp
        CALL pade_eval(N, z, b, omega, padapp)
        znorm(iw) = padapp
        CALL pade_eval(N, z, c, omega, padapp)
        shift(iw) = padapp
      ENDDO
    ELSE
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
    ENDIF ! fbw
    !
    IF (fbw) THEN
      WRITE(stdout, '(5x, a)') '   pade    Re[znorm]   Re[delta] [meV]   Re[shift] [meV]'
      WRITE(stdout, '(5x, i6, 3ES15.6)') N, REAL(znorm(1)), REAL(delta(1)) * 1000.d0, REAL(shift(1)) * 1000.d0
    ELSE
      WRITE(stdout, '(5x, a)') '   pade    Re[znorm]   Re[delta] [meV]'
      WRITE(stdout, '(5x, i6, 2ES15.6)') N, REAL(znorm(1)), REAL(delta(1)) * 1000.d0
    ENDIF
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
    !!
    !! computes kernels K_{+}(n, n', T) and K_{-}(n, n', T)
    !! reference W. E. Pickett, PRB 26, 1186 (1982)
    !!
    !! SH: Modified to allow for sparse sampling of Matsubara freq. (Nov 2021).
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : zero
    USE ep_constants,  ONLY : pi
    USE global_var,    ONLY : gtemp
    USE supercond_common,     ONLY : nsiw, keri, wsn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: n
    !! Nr. of Matsubara frequencies
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: omega
    !! frequency imag-axis
    REAL(KIND = DP) :: lambda_eph
    !! electron-phonon coupling
    !
    ! SH: nsiw is replaced with "1+largest Matsubara index"
    !       for compatibility with the general case of sparse sampling
    ! RM: n = nsiw(itemp) for uniform sampling
    !
    n = wsn(nsiw(itemp)) + 1
    !
    ALLOCATE(keri(2 * n), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_iaxis', 'Error allocating keri', 1)
    keri(:) = zero
    !
    DO iw = 1, 2 * n
      omega = DBLE(2 * (iw - 1)) * pi * gtemp(itemp)
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
    !!
    !! computes lambda(n - n')
    !! reference W. E. Pickett, PRB 26, 1186 (1982)
    !!
    USE kinds,            ONLY : DP
    USE input,            ONLY : nqstep
    USE ep_constants,     ONLY : zero
    USE supercond_common, ONLY : a2f_tmp, wsph, dwsph
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
      lambda_eph = lambda_eph + wsph(iwph) * a2f_tmp(iwph) &
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
    !!
    !! computes kernels K_{+}(w, iw_n, T) and K_{-}(w, iw_n, T)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE input,         ONLY : muc
    USE supercond_common,     ONLY : nsw, nsiw, ws, wsi, deltai, dsumi, zsumi
    USE ep_constants,  ONLY : zero
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
    REAL(KIND = DP) :: zesqrt
    !! w / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: desqrt
    !! \delta / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: kernelr
    !! 2 * Re[lambda(w - iw_n)]
    REAL(KIND = DP) :: kerneli
    !! 2 * Im[lambda(w - iw_n)]
    !
    COMPLEX(KIND = DP) :: lambda_eph
    !! electron-phonon coupling lambda(w - iw_n)
    !
    ALLOCATE(dsumi(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error allocating dsumi', 1)
    ALLOCATE(zsumi(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_iso_analytic_cont', 'Error allocating zsumi', 1)
    dsumi(:) = zero
    zsumi(:) = zero
    !
    DO iwp = 1, nsiw(itemp) ! loop over iw_n
      esqrt = 1.d0 / DSQRT(wsi(iwp) * wsi(iwp) + deltai(iwp) * deltai(iwp))
      zesqrt =  wsi(iwp) * esqrt
      desqrt =  deltai(iwp) * esqrt
      DO iw = 1, nsw ! loop over omega
        CALL lambdai_iso(ws(iw), wsi(iwp), lambda_eph)
        kernelr = 2.d0 * REAL(lambda_eph)
        kerneli = 2.d0 * AIMAG(lambda_eph)
        zsumi(iw) = zsumi(iw) + zesqrt * kerneli
        dsumi(iw) = dsumi(iw) + desqrt * (kernelr - 2.d0 * muc)
      ENDDO
    ENDDO
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
    !!
    !! computes lambda(w-iw_n)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds, ONLY : DP
    USE input,         ONLY : nqstep
    USE supercond_common,     ONLY : a2f_tmp, wsph, dwsph
    USE ep_constants,  ONLY : ci, czero
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
      lambda_eph = lambda_eph + wsph(iwph) * a2f_tmp(iwph) &
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
    USE input,             ONLY : nsiter, nstemp, broyden_beta, broyden_ndim
    USE global_var,        ONLY : gtemp
    USE supercond_common,         ONLY : nsw, delta, deltap, gap0
    USE ep_constants,      ONLY : kelvin2eV, ci, zero
    USE supercond,         ONLY : gen_freqgrid_raxis, eliashberg_grid
    USE utilities,         ONLY : mix_wrap
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
        ! SH: mix_wrap is used to allow for linear mixing
        CALL mix_wrap(nsw, rdeltaout, rdeltain, broyden_beta, &
                      iter, broyden_ndim, conv, df1, dv1)
        CALL mix_wrap(nsw, cdeltaout, cdeltain, broyden_beta, &
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
                    '  gap_edge(', itemp, ') = ', gap0, ' eV ', &
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
    USE input,         ONLY : nqstep, nsiter, muc, conv_thr_raxis, &
                              kerwrite, kerread
    USE global_var,    ONLY : gtemp
    USE supercond_common,     ONLY : nsw, ws, dwsph, gap0, bewph, fdwp, &
                              kp, km, delta, deltap, znorm, wsph
    USE ep_constants,  ONLY : kelvin2eV, ci, eps6, zero, one, czero
    USE io_supercond,  ONLY : eliashberg_write_raxis
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
    REAL(KIND = DP) :: inv_etemp
    !! Invese temperature inv_etemp = 1/temp. Defined for efficiency reason
    REAL(KIND = DP) :: a, b, c, d
    !! Temporary variables for reading kernelp and kernelm from file
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP) :: zesqrt
    !! w / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: desqrt
    !! \delta / sqrt{w^2+\delta^2}
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
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
    IF (iter == 1) THEN
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
      deltap(:) = czero
      deltap(:) = gap0
      kp(:, :) = czero
      km(:, :) = czero
      !
      inv_etemp = one / gtemp(itemp)
      !
      ! Fermi Dirac distribution
      fdwp(iw) = zero
      DO iw = 1, nsw
        IF (ABS(gtemp(itemp)) >  eps6) THEN
          fdwp(iw) = wgauss(-ws(iw) * inv_etemp, -99)
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
            bewph(iwph) = wgauss(-wsph(iwph) * inv_etemp, -99)
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
    IF (iter == 1) THEN
      DO iwp = 1, nsw ! loop over omega_prime
        DO iw = 1, nsw ! loop over omega
          !
          ! read the kernels from file if they were calculated before otherwise calculate them
          IF (kerread) THEN
            READ(iufilker) a, b, c, d
            kp(iw, iwp) = a + ci * b
            km(iw, iwp) = c + ci * d
          ENDIF
          IF (kerwrite) THEN
            CALL kernel_raxis(iw, iwp, kernelp, kernelm)
            kp(iw, iwp) = kernelp
            km(iw, iwp) = kernelm
            WRITE(iufilker) REAL(kp(iw, iwp)), AIMAG(kp(iw, iwp)), &
                            REAL(km(iw, iwp)), AIMAG(km(iw, iwp))
          ENDIF
        ENDDO ! iw
      ENDDO ! iwp
    ENDIF
    !
    DO iwp = 1, nsw ! loop over omega_prime
      esqrt = 1.d0 / SQRT(ws(iwp) * ws(iwp) - deltap(iwp) * deltap(iwp))
      zesqrt =  REAL(ws(iwp) * esqrt)
      desqrt =  REAL(deltap(iwp) * esqrt)
      DO iw = 1, nsw ! loop over omega
        znorm(iw) = znorm(iw) + zesqrt * km(iw, iwp)
        delta(iw) = delta(iw) + desqrt &
                  * (kp(iw, iwp) - muc * (1.d0 - 2.d0 * fdwp(iwp)))
      ENDDO ! iw
    ENDDO ! iwp
    !
    absdelta = zero
    reldelta = zero
    DO iw = 1, nsw ! loop over omega
      znorm(iw) = znorm(iw) * dwsph
      znorm(iw) = 1.d0 - znorm(iw) / ws(iw)
      delta(iw) = delta(iw) * dwsph
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
    ENDIF
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
    !!
    !! computes kernels K_{+}(w', w, T) and K_{-}(w', w, T)
    !! reference M. J. Holcomb, PRB 54, 6648 (1996)
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : ci, eps6, zero, czero, one
    USE ep_constants,  ONLY : pi
    USE input,         ONLY : nqstep
    USE supercond_common,     ONLY : a2f_tmp, wsph, dwsph, ws, bewph, fdwp
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
      kernelp = kernelp + a2f_tmp(iwph) * (var1 * (g1 + g2) - var2 * (g3 + g4))
      kernelm = kernelm + a2f_tmp(iwph) * (var1 * (g1 - g2) + var2 * (g3 - g4))
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
    !!       (isotropic) equation for the specified range of temperatures
    !! HP: updated 5/11/2021
    !
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE input,             ONLY : nstemp, muc, tc_linear_solver
    USE global_var,        ONLY : gtemp
    USE supercond_common,         ONLY : nsiw, wsi, wsn
    USE ep_constants,      ONLY : Kelvin2eV, zero, pi
    USE supercond,         ONLY : gen_freqgrid_iaxis, eliashberg_grid, crit_temp_solver
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
      DEALLOCATE(wsn, STAT = ierr)
      IF (ierr /= 0) CALL errore('crit_temp_iso', 'Error deallocating wsn', 1)
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
    !-----------------------------------------------------------------------
    SUBROUTINE mu_inter_iso(itemp, muintr, nel, nstate)
    !-----------------------------------------------------------------------
    !!
    !! SH: To find the superconducting state chemical potential
    !!       using the Newton-Raphson method (Nov 2021).
    !!
    !! HP: updated for the isotropic calculation (Jan 2022)
    !!
    USE kinds,          ONLY : DP
    USE supercond_common,      ONLY : wsi, nsiw, deltaip, znormip, shiftip, &
                               en, ndos, dosen, dosef
    USE global_var,     ONLY : gtemp
    USE input,          ONLY : dos_del, broyden_beta, nsiter
    USE ep_constants,   ONLY : zero, one, eps6
    !
    IMPLICIT NONE
    !
    INTEGER         :: itemp
    !! Temperature
    REAL(KIND = DP) :: muintr
    !! Interacting chemical potential: initial value
    REAL(KIND = DP) :: nel
    !! Without SOC: Nr. of electrons within Fermi window
    !! With SOC: 0.5 * Nr. of electrons within Fermi window
    REAL(KIND = DP) :: nstate
    !! Nr. of states within Fermi window
    !
    ! Local variables
    LOGICAL :: conv
    !! Convergence parameter
    INTEGER :: iw
    !! Index for Matsubara frequencies
    INTEGER :: ie
    !! Counter on energy grid
    INTEGER :: iter
    !! Counter for iterations
    REAL(KIND = DP) :: delta
    !! Temporary variable to store energy difference
    REAL(KIND = DP) :: theta
    !! Theta in Eliashberg equations
    REAL(KIND = DP) :: muout, muin
    !! Variables for Newton's minimization
    REAL(KIND = DP) :: fmu
    !! To store f(mu) value
    REAL(KIND = DP) :: dmu
    !! To store d_f(mu)/d_mu value
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    !
    inv_dos = one / dosef
    muin = muintr
    !
    iter = 1
    conv = .FALSE.
    !
    DO WHILE (.NOT. conv .AND. iter <= nsiter)
      !
      fmu = zero ! f(mu)
      dmu = zero ! d_f(mu)/d_mu
      !
      !
      ! SH: Here, we have an explicit summation over full range of frequencies
      !      for even functions, and multiply that by two (loop over iw).
      !
      DO ie = 1, ndos
        DO iw = 1, nsiw(itemp)
          delta = en(ie) - muin + shiftip(iw)
          theta = (wsi(iw) * znormip(iw))**2.d0 + &
                  (en(ie) - muin + shiftip(iw))**2.d0 + &
                  (znormip(iw) * deltaip(iw))**2.d0
          fmu = fmu + 2.d0 * dos_del * dosen(ie) * delta / theta
          dmu = dmu + 2.d0 * dos_del * dosen(ie) * (2.d0 * delta**2.d0 - theta) / (theta**2.d0)
        ENDDO ! iw
      ENDDO ! ie
      !
      fmu = nstate - 2.d0 * gtemp(itemp) * inv_dos * fmu - nel
      dmu = - 2.d0 * gtemp(itemp) * inv_dos * dmu
      !
      muout = muin - fmu / dmu
      !
      ! HP: linear mixing
      muout = (1.d0 - ABS(broyden_beta)) * muin + ABS(broyden_beta) * muout
      !
      IF (ABS((muout - muin) / muin) <= eps6) conv = .TRUE.
      !
      muin = muout
      iter = iter + 1
    END DO
    !
    IF (.NOT. conv .AND. (iter - 1) == nsiter) &
      CALL errore('mu_inter_iso', 'Error failed to find the mu_inter_iso value',1)
    !
    muintr = muout
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mu_inter_iso
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
    USE input,            ONLY : fbw
    USE supercond_common, ONLY : wsi, wsn, deltai, deltaip, znormi, nznormi, &
                                 znormip, shifti, shiftip
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! gen_freqgrid_iaxis
    DEALLOCATE(wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating wsi', 1)
    DEALLOCATE(wsn, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating wsn', 1)
    ! sum_eliashberg_iso_iaxis
    DEALLOCATE(deltai, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating deltai', 1)
    DEALLOCATE(deltaip, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating deltaip', 1)
    DEALLOCATE(znormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating znormi', 1)
    DEALLOCATE(nznormi, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating nznormi', 1)
    ! SH: variables for the fbw calculations
    IF (fbw) THEN
      DEALLOCATE(znormip, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating znormip', 1)
      DEALLOCATE(shifti, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating shifti', 1)
      DEALLOCATE(shiftip, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_iaxis', 'Error deallocating shiftip', 1)
    ENDIF
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
    USE input,         ONLY : limag, lacon, lreal, fbw
    USE supercond_common,     ONLY : ws, delta, deltap, znorm, znormp, shift, &
                              gp, gm, dsumi, zsumi, fdwp, bewph, kp, km
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ! gen_freqgrid_iaxis and gen_freqgrid_raxis
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
    ! SH: to deallocate fbw-related arrays
    IF (fbw) THEN
      DEALLOCATE(shift, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso_raxis', 'Error deallocating shift', 1)
    ENDIF
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
    USE input,         ONLY : limag, fbw
    USE global_var,    ONLY : gtemp
    USE supercond_common, ONLY : a2f_tmp, wsph, nsiw, en, dosen
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
    DEALLOCATE(a2f_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating a2f_tmp', 1)
    ! eliashberg_grid
    IF (limag) THEN
      DEALLOCATE(nsiw, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating nsiw', 1)
    ENDIF
    ! SH: deallocate fbw-related array en, dosen
    IF (fbw) THEN
      ! read_dos
      DEALLOCATE(en, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating en', 1)
      DEALLOCATE(dosen, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_iso', 'Error deallocating dosen', 1)
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
