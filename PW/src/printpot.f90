!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE printpot_module
  !--------------------------------------------------------------------------
  !!
  !! Planar-averaged potential output for arbitrary directions.
  !!
CONTAINS

  SUBROUTINE printpot(rhog, rhor, idir)
    !------------------------------------------------------------------------
    !!
    !! Compute planar average of charge density, Hartree potential, local
    !! potential (from vltot), and their sum including the dipole correction
    !! sawtooth (via add_efield) along direction idir. Mirrors pp.x plot_num=11,
    !! and then average.x. Output is written to file <prefix>.ef1.
    !!
    !! For LSDA (nspin=2), spin-up and spin-down charge density columns are
    !! added between total charge and the potentials.
    !!
    USE kinds,             ONLY : DP
    USE constants,         ONLY : e2, fpi, AUTOEV, BOHR_RADIUS_ANGS
    USE cell_base,         ONLY : omega, alat, at, tpiba2
    USE gvect,             ONLY : ngm, gg
    USE scf,               ONLY : vltot
    USE lsda_mod,          ONLY : nspin
    USE extfield,          ONLY : tefield, dipfield
    USE fft_base,          ONLY : dfftp
    USE generate_function, ONLY : planar_average
    USE io_files,          ONLY : prefix, tmp_dir
    USE io_global,         ONLY : ionode
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(in) :: rhog(ngm)
    !! Total charge density in G-space (used for V_H via Poisson).
    REAL(DP),    INTENT(in) :: rhor(dfftp%nnr, nspin)
    !! Charge density in real space, all spin components.
    !! rhor(:,1) is passed to add_efield for the dipole computation.
    INTEGER,     INTENT(in) :: idir
    !! Averaging direction: 1, 2, or 3.
    !
    INTEGER :: ix, naxis
    INTEGER :: ig
    REAL(DP) :: L, S, dummy
    CHARACTER(len=1)   :: axis_label
    CHARACTER(len=256) :: filename
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    REAL(DP),    ALLOCATABLE :: pot3d(:), work(:)
    REAL(DP),    ALLOCATABLE :: rho1d(:), spup1d(:), spdn1d(:)
    REAL(DP),    ALLOCATABLE :: vh1d(:), vloc1d(:), vtot1d(:)
    !
    SELECT CASE (idir)
    CASE (1); naxis = dfftp%nr1; axis_label = 'x'
    CASE (2); naxis = dfftp%nr2; axis_label = 'y'
    CASE (3); naxis = dfftp%nr3; axis_label = 'z'
    CASE DEFAULT
      CALL errore('printpot', 'idir must be 1, 2, or 3', 1)
    END SELECT
    !
    L = alat * SQRT(at(1,idir)**2 + at(2,idir)**2 + at(3,idir)**2)
    S = omega / L
    !
    ALLOCATE(aux(dfftp%nnr), pot3d(dfftp%nnr), work(dfftp%nnr))
    ALLOCATE(rho1d(naxis), spup1d(naxis), spdn1d(naxis))
    ALLOCATE(vh1d(naxis), vloc1d(naxis), vtot1d(naxis))
    !
    ! Total charge density: planar average directly from real space
    !
    work(:) = rhor(:, 1)
    CALL planar_average(dfftp%nnr, naxis, idir, 0, .FALSE., work, rho1d)
    !
    ! Spin-up and spin-down for LSDA (nspin=2)
    !
    IF (nspin == 2) THEN
      work(:) = (rhor(:, 1) + rhor(:, 2)) * 0.5d0
      CALL planar_average(dfftp%nnr, naxis, idir, 0, .FALSE., work, spup1d)
      work(:) = (rhor(:, 1) - rhor(:, 2)) * 0.5d0
      CALL planar_average(dfftp%nnr, naxis, idir, 0, .FALSE., work, spdn1d)
    END IF
    !
    ! Hartree potential V_H(G) = e2*4pi/G^2 * rho(G), G -> R, then planar average
    !
    aux(:) = (0.d0, 0.d0)
    DO ig = 1, ngm
      IF (gg(ig) .GT. 0.d0) aux(dfftp%nl(ig)) = e2 * fpi / (gg(ig) * tpiba2) * rhog(ig)
    END DO
    CALL g_avg(aux, pot3d, idir, vh1d)
    !
    ! Local potential from vltot, planar average
    !
    CALL planar_average(dfftp%nnr, naxis, idir, 0, .FALSE., vltot, vloc1d)
    !
    ! Total: V_H (in pot3d) + vltot, add sawtooth via add_efield, then planar average
    !
    pot3d(:) = pot3d(:) + vltot(:)
    IF (tefield .AND. dipfield) CALL add_efield(pot3d, dummy, rhor(:, 1), .TRUE.)
    CALL planar_average(dfftp%nnr, naxis, idir, 0, .FALSE., pot3d, vtot1d)
    !
    IF (ionode) THEN
      filename = TRIM(tmp_dir)//TRIM(prefix)//".ef1"
      OPEN(UNIT=4, FILE=filename, STATUS="UNKNOWN", ACTION="WRITE")
      IF (nspin == 2) THEN
        WRITE(UNIT=4, FMT=9060) axis_label
        DO ix = 1, naxis
          WRITE(UNIT=4, FMT=9061) &
            DBLE(ix - 1) / DBLE(naxis) * L * BOHR_RADIUS_ANGS, &
            rho1d(ix)  * S / BOHR_RADIUS_ANGS, &
            spup1d(ix) * S / BOHR_RADIUS_ANGS, &
            spdn1d(ix) * S / BOHR_RADIUS_ANGS, &
            vh1d(ix)   * AUTOEV / e2, &
            vloc1d(ix) * AUTOEV / e2, &
            vtot1d(ix) * AUTOEV / e2
        END DO
      ELSE
        WRITE(UNIT=4, FMT=9050) axis_label
        DO ix = 1, naxis
          WRITE(UNIT=4, FMT=9051) &
            DBLE(ix - 1) / DBLE(naxis) * L * BOHR_RADIUS_ANGS, &
            rho1d(ix)  * S / BOHR_RADIUS_ANGS, &
            vh1d(ix)   * AUTOEV / e2, &
            vloc1d(ix) * AUTOEV / e2, &
            vtot1d(ix) * AUTOEV / e2
        END DO
      END IF
      CLOSE(UNIT=4)
    END IF
    !
    DEALLOCATE(aux, pot3d, work, rho1d, spup1d, spdn1d, vh1d, vloc1d, vtot1d)
    !
9050 FORMAT('#', A1, ' (A)', 2X, 'Avg charge (e/A)', 2X, 'Avg V_hartree (eV)', 2X, &
            'Avg V_local (eV)', 2X, 'Avg V_hart+V_loc (eV)')
9051 FORMAT(F8.4, F20.7, F20.7, F18.7, F18.7)
9060 FORMAT('#', A1, ' (A)', 2X, 'Avg charge (e/A)', 2X, 'Avg spin-up (e/A)', 2X, &
            'Avg spin-dn (e/A)', 2X, 'Avg V_hartree (eV)', 2X, &
            'Avg V_local (eV)', 2X, 'Avg V_hart+V_loc (eV)')
9061 FORMAT(F8.4, F20.7, F20.7, F20.7, F20.7, F18.7, F18.7)
  END SUBROUTINE printpot

  SUBROUTINE g_avg(aux, work, idir, f1d)
    !!
    !! Apply gamma_only symmetry, FFT G -> R, then planar average along idir.
    !!
    USE kinds,             ONLY : DP
    USE control_flags,     ONLY : gamma_only
    USE fft_base,          ONLY : dfftp
    USE fft_interfaces,    ONLY : invfft
    USE generate_function, ONLY : planar_average
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)
    !! G-space array (modified in-place by invfft).
    REAL(DP),    INTENT(out)   :: work(dfftp%nnr)
    !! Real-space workspace; holds V(r) on exit.
    INTEGER,     INTENT(in)    :: idir
    REAL(DP),    INTENT(out)   :: f1d(:)
    !
    IF (gamma_only) aux(dfftp%nlm(:)) = CONJG(aux(dfftp%nl(:)))
    CALL invfft('Rho', aux, dfftp)
    work(:) = REAL(aux(:), KIND=DP)
    CALL planar_average(dfftp%nnr, SIZE(f1d), idir, 0, .FALSE., work, f1d)
  END SUBROUTINE g_avg

END MODULE printpot_module
