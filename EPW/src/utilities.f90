  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE utilities
  !----------------------------------------------------------------------
  !!
  !! This module contains the routines associated with Broyden's method,
  !! Pade' approximants, DOS and Fermi level determination
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !!!!!
    !-----------------------------------------------------------------------
    SUBROUTINE mix_wrap(ndim, deltaout, deltain, alphamix, iter, n_iter, conv, df, dv)
    !-----------------------------------------------------------------------
    !!
    !!  SH: Wrapper for the linear/broyden mixings (Nov 2021).
    !!        Note: the linear mixing option is implemented for 
    !!        benchmarking/development purposes and can be invoked by
    !!        setting the broyden_beta parameter to a negative value.
    !!        
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: conv
    !! If true convergence reached
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(KIND = DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: deltaout(ndim)
    !! output delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: deltain(ndim)
    !! delta at previous iteration
    REAL(KIND = DP), INTENT(inout) :: df(ndim, n_iter)
    !! arrays containing info from previous iterations
    REAL(KIND = DP), INTENT(inout) :: dv(ndim, n_iter)
    !! arrays containing info from previous iterations
    !
    IF (alphamix < zero ) THEN
      IF ((iter <= 3) .AND. (alphamix < -0.2d0)) THEN
        CALL mix_linear(ndim, deltaout, deltain, -0.2d0)
      ELSE
        CALL mix_linear(ndim, deltaout, deltain, alphamix)
      ENDIF
    ELSE
      IF ((iter <= 3) .AND. (alphamix > 0.2d0)) THEN
        CALL mix_broyden(ndim, deltaout, deltain, 0.2d0, iter, n_iter, conv, df, dv)
      ELSE
        CALL mix_broyden(ndim, deltaout, deltain, alphamix, iter, n_iter, conv, df, dv)
      ENDIF
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_wrap
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mix_linear(ndim, arout, arin, mixf)
    !-----------------------------------------------------------------------
    !!
    !!  SH: Simple linear mixing for gap, normalization, shift, etc (Nov 2021).
    !!
    !
    USE kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    !
    REAL(KIND = DP), INTENT(in) :: mixf
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: arout(ndim)
    !! output delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: arin(ndim)
    !! delta at previous iteration
    !
    ! Local variables
    INTEGER :: i
    !
    DO i = 1, ndim
      arin(i) = (1.d0 - DABS(mixf)) * arin(i) + DABS(mixf) * arout(i)
    ENDDO
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_linear
    !-----------------------------------------------------------------------
    !
    !!!!!
    !-----------------------------------------------------------------------
    SUBROUTINE mix_broyden(ndim, deltaout, deltain, alphamix, iter, n_iter, conv, df, dv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nsiter
    USE eliashbergcom, ONLY : nkfs, nbndfs
    USE constants_epw, ONLY : eps2, zero, one, two
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: conv
    !! If true convergence reached
    !
    INTEGER, INTENT(in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT(in) :: iter
    !! Current iteration number
    INTEGER, INTENT(in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL(KIND = DP), INTENT(in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL(KIND = DP), INTENT(inout) :: deltaout(ndim)
    !! output delta at current iteration
    REAL(KIND = DP), INTENT(inout) :: deltain(ndim)
    !! delta at previous iteration
    REAL(KIND = DP), INTENT(inout) :: df(ndim, n_iter)
    !! arrays containing info from previous iterations
    REAL(KIND = DP), INTENT(inout) :: dv(ndim, n_iter)
    !! arrays containing info from previous iterations
    !
    ! Local variables
    INTEGER, PARAMETER :: maxter = 8
    !! max number of iterations used in mixing: n_iter must be <= maxter
    INTEGER :: n
    !! Counter on deltain/deltaout dimension (1 to nmin)
    INTEGER :: i, j
    !! Counter on iterations (1 to iter_used)
    INTEGER :: iter_used
    !! number of iterations used in mixing
    INTEGER :: ipos
    !! position at which results from the present iteraction are stored
    INTEGER :: inext
    !! position at which results for the next iteration are stored
    INTEGER :: ierr
    !! Error status
    INTEGER :: info
    !! Exit info DSYTRF and DSYTRI lapack subroutines
    INTEGER :: iwork(maxter)
    !! Workspace array in DSYTRF and DSYTRI lapack subroutines
    !
    REAL(KIND = DP) :: gammamix
    !! Mixing parameter
    REAL(KIND = DP) :: norm
    !! norm of df
    REAL(KIND = DP) :: inv_norm
    !! 1.0/norm. Defined for efficiency reasons
    REAL(KIND = DP) :: wg0
    !!
    REAL(KIND = DP) :: work(maxter)
    !!
    REAL(KIND = DP) :: wg(maxter)
    !!
    REAL(KIND = DP), ALLOCATABLE :: deltainsave(:)
    !! Array to store deltain from previous iteration
    REAL(KIND = DP) :: beta(maxter, maxter)
    !!
    REAL(KIND = DP), EXTERNAL :: DDOT
    !! Inner product of two vectors
    REAL(KIND = DP), EXTERNAL :: DNRM2
    !! Norm of a vector
    !!
    !
    IF (conv) RETURN
    !
    ! adjustable parameters as suggested in the original paper
    wg0 = eps2
    wg  = maxter * one
    !
    IF (iter < 1) CALL errore('mix_broyden', 'iter is smaller than 1', 1)
    IF (n_iter > maxter) CALL errore('mix_broyden', 'n_iter is too big', 1)
    IF (ndim <= 0) CALL errore('mix_broyden', 'ndim <= 0', 1)
    !
    ALLOCATE(deltainsave(ndim), STAT = ierr)
    IF (ierr /= 0) CALL errore('mix_broyden', 'Error allocating deltainsave', 1)
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = MIN(iter - 1, n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos = iter - 1 until ipos = n_iter, then back to 1, 2,...
    !
    ipos = iter - 1 - ((iter - 2) / n_iter) * n_iter
    !
    DO n = 1, ndim
      deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF (iter > 1) THEN
      DO n = 1, ndim
        df(n, ipos) = deltaout(n) - df(n, ipos)
        dv(n, ipos) = deltain(n)  - dv(n, ipos)
      ENDDO
      norm = (DNRM2(ndim, df(1, ipos), 1))**two
      norm = DSQRT(norm)
      inv_norm = one / norm
      ! DSCAL scales df and dv by inv_norm
      CALL DSCAL(ndim, inv_norm, df(1, ipos), 1)
      CALL DSCAL(ndim, inv_norm, dv(1, ipos), 1)
    ENDIF
    !
    DO i = 1, iter_used
      DO j = i + 1, iter_used
        beta(i, j) = wg(i) * wg(j) * DDOT(ndim, df(1, j), 1, df(1, i), 1)
      ENDDO
      beta(i, i) = wg0**two + wg(i)**two
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix
    !
    CALL DSYTRF('u', iter_used, beta, maxter, iwork, work, maxter, info)
    CALL errore('mix_broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('u', iter_used, beta, maxter, iwork, work, info)
    CALL errore('mix_broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
      DO j = i + 1, iter_used
        beta(j, i) = beta(i, j)
      ENDDO
    ENDDO
    !
    DO i = 1, iter_used
      work(i) = DDOT(ndim, df(1, i), 1, deltaout, 1)
    ENDDO
    !
    DO n = 1, ndim
      deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
      gammamix = zero
      DO j = 1, iter_used
        gammamix = gammamix + beta(j, i) * wg(j) * work(j)
      ENDDO
      !
      DO n = 1, ndim
        deltain(n) = deltain(n) - wg(i) * gammamix * (alphamix * df(n, i) + dv(n, i))
      ENDDO
    ENDDO
    !
    inext = iter - ((iter - 1) / n_iter) * n_iter
    df(:, inext) = deltaout(:)
    dv(:, inext) = deltainsave(:)
    !
    DEALLOCATE(deltainsave, STAT = ierr)
    IF (ierr /= 0) CALL errore('mix_broyden', 'Error deallocating deltainsave', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mix_broyden
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    ! General note on PADE:
    ! Lebegue, Arnaud, Alouani, and Blochel [PRB 67, 155208 (2003)]
    ! state that when they use Pade of order N = 12 (resulting in
    ! numerator or order (N-2)/2 = 5 and denominator N/2 = 6),
    ! they obtain extremely stable fits and the quasiparticle energies
    ! are essentially identical to those obtained using the contour
    ! deformation method.
    !
    ! using this sub:
    !
    ! INTEGER :: N
    ! COMPLEX(KIND = DP) :: z(N), u(N), a(N), w, padapp
    !
    ! CALL pade_coeff(N, z, u, a)
    ! CALL pade_eval(N, z, a, w, padapp)
    !
    !-----------------------------------------------------------
    SUBROUTINE pade_coeff(N, z, u, a)
    !-----------------------------------------------------------
    ! N-point Pade' approximant - find the Pade' coefficients
    !
    ! This routine uses the recursive algorithm described in
    ! HJ Vidberg and JW Serene, "Solving the Eliashberg equations
    ! by means of N-point Pade' approximants", J Low Temp Phys
    ! 29, 179 (1977). The notation adopted here is the same as
    ! in the above manuscript.
    !
    !-----------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE, INTRINSIC :: IEEE_ARITHMETIC
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: N
    !! order of the Pade' approximant
    COMPLEX(KIND = DP), INTENT(in) :: z(N)
    !! points at which the original function is known
    COMPLEX(KIND = DP), INTENT(in) :: u(N)
    !! values of the function at the z points
    COMPLEX(KIND = DP), INTENT(out) :: a(N)
    !! coefficients of the continued fraction
    !
    ! Local variables
    INTEGER :: i, p
    !! Counter over Pade' approximants
    REAL(KIND = DP) :: ar, ai
    !! Real/complex part of a
    COMPLEX(KIND = DP) :: tmp1, tmp2
    !! Temporary variables
    COMPLEX(KIND = DP) :: g(N, N)
    ! g(p, i) = g_p(z_i) in the notation of Vidberg and Serene
    !
    DO p = 1, N
      IF (p == 1) THEN
        DO i = 1, N
           g(p, i) = u(i)
        ENDDO
      ELSE
        DO i = p, N
          !  g(p, i) = (g(p - 1, p - 1) - g(p-1, i)) &
          !          / ((z(i) - z(p - 1)) * g(p - 1, i))
          !
          ! this seems necessary to avoid nasty NaN when
          ! still don't quite understand why the procedure
          ! becomes unstable - certainly it happens only
          ! when u(:) is very small
          !
          !IF (ABS(g(p - 1, i)) == 0) THEN
          !  WRITE(stdout, '(4x, "fitting parameter too small. g(p-1,i)= ", 2f9.5)') g(p - 1, i)
          !  STOP
          !ENDIF
          !
          tmp1 = g(p - 1, p - 1) / g(p - 1, i)
          tmp2 = g(p - 1, i)     / g(p - 1, i)
          g(p, i) = (tmp1 - tmp2) / (z(i) - z(p - 1))
          !
        ENDDO
      ENDIF
      a(p) = g(p, p)
      !
      ! check whether a(p) is not NaN
      !
      ar = REAL(a(p))
      ai = AIMAG(a(p))
      !
      IF (IEEE_IS_NAN(ar) .OR. IEEE_IS_NAN(ai)) THEN
        !WRITE(stdout, *) (z(i), i = 1, N)
        !WRITE(stdout, *) (u(i), i = 1, N)
        !WRITE(stdout, *) (a(i), i = 1, N)
        WRITE(stdout, *) 'One or more Pade coefficients are NaN'
        !CALL errore('pade_coeff', 'one or more coefficients are NaN', 1)
      ENDIF
      !
    ENDDO
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE pade_coeff
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------
    SUBROUTINE pade_eval(N, z, a, w, padapp)
    !-----------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : zero, one
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: N
    !! order of the Pade' approximant
    COMPLEX(KIND = DP), INTENT(in) :: w
    !! point at which we need the approximant
    COMPLEX(KIND = DP), INTENT(out) :: padapp
    !! value of the approximant at the point w
    COMPLEX(KIND = DP), INTENT(in) :: z(N)
    !! points at which the original function is known
    COMPLEX(KIND = DP), INTENT(in) :: a(N)
    !! coefficients of the continued fraction
    !
    ! Local variables
    INTEGER :: i
    !! Counter over Pade' approximants
    COMPLEX(KIND = DP) :: acap(0:N)
    !!
    COMPLEX(KIND = DP) :: bcap(0:N)
    !!
    !
    acap(0) = zero
    acap(1) = a(1)
    bcap(0) = one
    bcap(1) = one
    !
    DO i = 2, N
      acap(i) = acap(i - 1) + (w - z(i - 1)) * a(i) * acap(i - 2)
      bcap(i) = bcap(i - 1) + (w - z(i - 1)) * a(i) * bcap(i - 2)
    ENDDO
    padapp = acap(N) / bcap(N)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE pade_eval
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_dos(itemp, ef0, dos)
    !-----------------------------------------------------------------------
    !!
    !! This routine computes the density of states at a given fermi level.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : two, eps16, ryd2mev
    USE epwcom,        ONLY : ngaussw, nstemp, nbndsub, degaussw
    USE elph2,         ONLY : etf, nkqf, wkf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(inout) :: dos(nstemp)
    !! DOS to compute for the temperature itemp.
    !
    ! Local variables
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !
    ! divide by two to have DOS/spin
    IF (ABS(degaussw) < eps16) THEN
      ! use 1 meV instead
      dos(itemp) = dos_ef(ngaussw, 1.0d0 / ryd2mev, ef0(itemp), etf, wkf, nkqf, nbndsub) / two
    ELSE
      dos(itemp) = dos_ef(ngaussw, degaussw, ef0(itemp), etf, wkf, nkqf, nbndsub) / two
    ENDIF
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_dos
    !-----------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE broadening(ik, ikk, ikq, w, vmefp, eta)
    !--------------------------------------------------------------------------
    !!
    !! This routine computes the adaptative broadening
    !! It requires electronic and phononic velocities
    !! The implemented equation is Eq. 18 of Computer Physics Communications 185, 1747 (2014)
    !! 2019: Samuel Ponce & Francesco Macheda
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : alat, bg
    USE elph2,         ONLY : nbndfst, nkf, dmef, vmef, ibndmin, etf
    USE epwcom,        ONLY : vme, nqf1, nqf2, nqf3
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : eps40, ryd2mev, twopi, zero, two, eps6, eps8, eps4
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Current k-point on that core
    INTEGER, INTENT(in) :: ikk
    !! Current k point on that core (ikk = 2 * ik + 1)
    INTEGER, INTENT(in) :: ikq
    !! k+q point on that core
    REAL(KIND = DP), INTENT(in) :: w(nmodes)
    !! Phonon frequencies
    REAL(KIND = DP), INTENT(out) :: eta(nmodes, nbndfst, nkf)
    !! Adaptative smearing value
    COMPLEX(KIND = DP), INTENT(in) :: vmefp(3, nmodes, nmodes)
    !! Phonon velocity
    !
    ! Local variables
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: imode, jmode
    !! Mode index
    INTEGER :: n_av
    !! To average eta_av
    REAL(KIND = DP) :: vel_diff(3)
    !! Velocity difference when computed adaptative broadening
    REAL(KIND = DP) :: eta_tmp(3)
    !! Temporary adaptative broadening
    REAL(KIND = DP) :: eta_deg(nmodes, nbndfst)
    !! Average eta over degenerate states
    REAL(KIND = DP) :: e_1
    !! Eigenvalue 1 for deg. testing
    REAL(KIND = DP) :: e_2
    !! Eigenvalue 2 for deg. testing
    REAL(KIND = DP) :: w_1
    !! Phonon frequency for degeneracy checking
    REAL(KIND = DP) :: w_2
    !! Phonon frequency for degeneracy checking
    REAL(KIND = DP) :: vmeq(3, nmodes)
    !! Local phonon velocity
    REAL(KIND = DP) :: vmek(3, nbndfst)
    !! Local electron velocity
    REAL(KIND = DP) :: vmeq_av(3)
    !! Average phonon velocity
    REAL(KIND = DP) :: vmek_av(3)
    !! Average phonon velocity
    !
    eta_deg(:, :) = zero
    vmeq(:, :) = zero
    vmek(:, :) = zero
    !
    ! First average the phonon velocities
    DO imode = 1, nmodes
      w_1 = w(imode)
      vmeq_av(:) = zero
      n_av = 0
      DO jmode = 1, nmodes
        w_2 = w(jmode)
        IF (ABS(w_2 - w_1) < eps6) THEN
          n_av   = n_av + 1
          vmeq_av(:) = vmeq_av(:) + REAL(vmefp(:, jmode, jmode), KIND = DP)
        ENDIF
      ENDDO
      vmeq(:, imode) = vmeq_av(:) / FLOAT(n_av)
    ENDDO
    !
    ! Average electron velocity
    DO ibnd = 1, nbndfst
      e_1 = etf(ibndmin - 1 + ibnd, ikk)
      vmek_av(:) = zero
      n_av   = 0
      DO jbnd = 1, nbndfst
        e_2 = etf(ibndmin - 1 + jbnd, ikk)
        IF (ABS(e_2 - e_1) < eps4) THEN
          n_av = n_av + 1
          IF (vme == 'wannier') THEN
            vmek_av(:) = vmek_av(:) + REAL(vmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq), KIND = DP)
          ELSE
            vmek_av(:) = vmek_av(:) + REAL(dmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq), KIND = DP)
          ENDIF
        ENDIF
      ENDDO
      vmek(:, ibnd) = vmek_av(:) / FLOAT(n_av)
    ENDDO
    !
    ! vmefp and vmef are obtained using irvec, which are without alat; therefore I multiply them to bg without alat
    DO ibnd = 1, nbndfst
      DO imode = 1, nmodes
        IF (w(imode) > 0) THEN
          vel_diff(:) = vmeq(:, imode) / (2d0 * w(imode)) - vmek(:, ibnd)
          !IF (vme == 'wannier') THEN
          !  vel_diff(:) = REAL(vmefp(:, imode, imode) / &
          !                    (2d0 * w(imode)) - vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikq))
          !ELSE
          !  vel_diff(:) = REAL(vmefp(:, imode ,imode) / &
          !                    (2d0 * w(imode)) - dmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikq))
          !ENDIF
          IF (SQRT(DOT_PRODUCT(vel_diff, vel_diff)) < eps40) THEN
            eta(imode, ibnd, ik) = 1.0d0 / ryd2mev
          ELSE
            eta_tmp(1) = (twopi / alat) * ABS(DOT_PRODUCT(vel_diff(:), bg(:, 1)) / DBLE(nqf1))
            eta_tmp(2) = (twopi / alat) * ABS(DOT_PRODUCT(vel_diff(:), bg(:, 2)) / DBLE(nqf2))
            eta_tmp(3) = (twopi / alat) * ABS(DOT_PRODUCT(vel_diff(:), bg(:, 3)) / DBLE(nqf3))
            !eta(imode, ibnd, ik) = MAXVAL(eta_tmp) !Eq. (24) of PRB 97 075405 (2015)
            !eta(imode, ibnd, ik) = DSQRT(eta_tmp(1)**2+eta_tmp(2)**2+eta_tmp(3)**2)/DSQRT(12d0)
            ! Eq. (18) of Computer Physics Communications 185 (2014), 1747-1758
            ! The prefactor 0.5 is arbitrary and is to speedup convergence
            eta(imode, ibnd, ik) = 0.5d0 * DSQRT(eta_tmp(1)**two + eta_tmp(2)**two + eta_tmp(3)**two) / SQRT(12.0d0)
          ENDIF
          !
          ! If the smearing is too small, set 1 meV. Too small value are numerically unstable.
          IF (eta(imode, ibnd, ik) * ryd2mev < 1.0d0) THEN
            eta(imode, ibnd, ik) = 1.0d0 / ryd2mev
          ENDIF
        ELSE
          ! Fixed value 1 meV
          eta(imode, ibnd, ik) = 1.0d0 / ryd2mev
        ENDIF
      ENDDO
    ENDDO
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE broadening
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE broadening_imp(ik, ikk, ikq, eta)
    !--------------------------------------------------------------------------
    !!
    !! This routine computes the adaptative broadening
    !! It requires electronic and phononic velocities
    !! The implemented equation is Eq. 18 of Computer Physics Communications 185, 1747 (2014)
    !! 2019: Samuel Ponce & Francesco Macheda
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : alat, bg
    USE elph2,         ONLY : nbndfst, nkf, dmef, vmef, ibndmin, etf
    USE epwcom,        ONLY : vme, nqf1, nqf2, nqf3, nkf1, nkf2, nkf3
    USE modes,         ONLY : nmodes
    USE constants_epw, ONLY : eps40, ryd2mev, twopi, zero, two, eps6, eps8, eps4
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Current k-point on that core
    INTEGER, INTENT(in) :: ikk
    !! Current k point on that core (ikk = 2 * ik + 1)
    INTEGER, INTENT(in) :: ikq
    !! k+q point on that core
    REAL(KIND = DP), INTENT(out) :: eta(nbndfst, nkf)
    !! Adaptative smearing value
    !
    ! Local variables
    INTEGER :: ibnd
    !! Band index
    INTEGER :: jbnd
    !! Band index
    INTEGER :: n_av
    !! To average eta_av
    REAL(KIND = DP) :: vel_diff(3)
    !! Velocity difference when computed adaptative broadening
    REAL(KIND = DP) :: eta_tmp(3)
    !! Temporary adaptative broadening
    REAL(KIND = DP) :: e_1
    !! Eigenvalue 1 for deg. testing
    REAL(KIND = DP) :: e_2
    !! Eigenvalue 2 for deg. testing
    REAL(KIND = DP) :: vmek(3, nbndfst)
    !! Local electron velocity
    REAL(KIND = DP) :: vmek_av(3)
    !! Average phonon velocity
    !
    vmek(:, :) = zero
    !
    ! Average electron velocity
    DO ibnd = 1, nbndfst
      e_1 = etf(ibndmin - 1 + ibnd, ikk)
      vmek_av(:) = zero
      n_av   = 0
      DO jbnd = 1, nbndfst
        e_2 = etf(ibndmin - 1 + jbnd, ikk)
        IF (ABS(e_2 - e_1) < eps4) THEN
          n_av = n_av + 1
          IF (vme == 'wannier') THEN
            vmek_av(:) = vmek_av(:) + REAL(vmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq), KIND = DP)
          ELSE
            vmek_av(:) = vmek_av(:) + REAL(dmef(:, ibndmin - 1 + jbnd, ibndmin - 1 + jbnd, ikq), KIND = DP)
          ENDIF
        ENDIF
      ENDDO
      vmek(:, ibnd) = vmek_av(:) / FLOAT(n_av)
    ENDDO
    !
    DO ibnd = 1, nbndfst
      IF (SQRT(DOT_PRODUCT(vmek(:, ibnd),vmek(:,ibnd))) < eps40) THEN
        eta(ibnd, ik) = 1.0d0 / ryd2mev
      ELSE
        !eta(ibnd, ik) = 1.0d0 / ryd2mev
        eta_tmp(1) = (twopi / alat) * ABS(DOT_PRODUCT(vmek(:,ibnd),bg(:, 1)) / DBLE(nkf1))
        eta_tmp(2) = (twopi / alat) * ABS(DOT_PRODUCT(vmek(:,ibnd),bg(:, 2)) / DBLE(nkf2))
        eta_tmp(3) = (twopi / alat) * ABS(DOT_PRODUCT(vmek(:,ibnd),bg(:, 3)) / DBLE(nkf3))
        eta(ibnd, ik) = (1.0d0/SQRT(12.0d0))*(1.0d0/3.0d0)*(eta_tmp(1)+eta_tmp(2)+eta_tmp(3))
      ENDIF
      IF (eta(ibnd, ik)*ryd2mev < 1.0d0) THEN
        eta(ibnd, ik) = 1.0d0 / ryd2mev
      ENDIF
    ENDDO
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE broadening_imp
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fermicarrier(itemp, etemp, ef0, efcb, ctype)
    !-----------------------------------------------------------------------
    !!
    !!  This routine computes the Fermi energy associated with a given
    !!  carrier concentration using bissection for insulators or
    !!  semi-conductors.
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : omega, alat, at
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : etf, nkf, wkf, efnew, nkqf, partion
    USE constants_epw, ONLY : ryd2ev, bohr2ang, ang2cm, eps5, kelvin2eV, &
                              zero, eps80
    USE noncollin_module, ONLY : noncolin
    USE pwcom,     ONLY : nelec
    USE epwcom,    ONLY : int_mob, nbndsub, ncarrier, nstemp, fermi_energy, &
                          system_2d, carrier, efermi_read, assume_metal, ngaussw
    USE klist_epw, ONLY : isk_dummy
    USE mp,        ONLY : mp_barrier, mp_sum, mp_max, mp_min
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    INTEGER, INTENT(out) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in kBT [Ry] unit.
    REAL(KIND = DP), INTENT(inout) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(inout) :: efcb(nstemp)
    !! Second fermi level for the temperature itemp
    REAL(KIND = DP), EXTERNAL :: efermig
    !! External function to calculate the fermi energy
    !
    ! Local variables
    INTEGER :: i
    !! Index for the bisection iteration
    INTEGER :: ik
    !! k-point index per pool
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ivbm
    !! Index of the VBM
    INTEGER :: icbm
    !! Index of the CBM
    INTEGER, PARAMETER :: maxiter = 500 ! 300
    !! Maximum interation
    REAL(KIND = DP) :: fermi
    !! Fermi level returned
    REAL(KIND = DP) :: fermicb
    !! Fermi level returned for second Fermi level
    REAL(KIND = DP) :: fnk
    !! Fermi-Diract occupation
    REAL(KIND = DP) :: ks_exp(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT
    REAL(KIND = DP) :: ks_expcb(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT for CB
    REAL(KIND = DP) :: fermi_exp
    !! Fermi level in exponential format
    REAL(KIND = DP) :: rel_err
    !! Relative error
    REAL(KIND = DP) :: factor
    !! Factor that goes from number of carrier per unit cell to number of
    !! carrier per cm^-3
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: evbm
    !! Energy of the VBM
    REAL(KIND = DP) :: ecbm
    !! Energy of the CBM
    REAL(KIND = DP) :: ef_tmp
    !! Energy of the current Fermi level for the bisection method
    REAL(KIND = DP) :: elw
    !! Energy lower bound for the bisection method
    REAL(KIND = DP) :: eup
    !! Energy upper bound for the bisection method
    REAL(KIND = DP) :: hole_density
    !! Hole carrier density
    REAL(KIND = DP) :: electron_density
    !! Electron carrier density
    REAL(KIND = DP), PARAMETER :: maxarg = 200.d0
    !! Maximum value for the argument of the exponential
    REAL(KIND = DP) :: ncarrierp
    !! ncarrier*fraction of ionized impurities
    !
    ncarrierp = ncarrier*partion(itemp)
    !
    IF (assume_metal) THEN
      !! set conduction band chemical potential to 0 since it is irrelevent
      ctype = -1  ! act like it's for holes
      efcb(itemp) = zero
      ef0(itemp) = efermig(etf, nbndsub, nkqf, nelec, wkf, etemp, ngaussw, 0, isk_dummy)
      RETURN
    ENDIF
    ef_tmp  = zero
    fermi   = zero
    fermicb = zero
    inv_cell = 1.0d0 / omega
    !
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF (system_2d) inv_cell = inv_cell * at(3, 3) * alat
    ! vbm index
    IF (noncolin) THEN
      ivbm = FLOOR(nelec / 1.0d0)
    ELSE
      ivbm = FLOOR(nelec / 2.0d0)
    ENDIF
    icbm = ivbm + 1 ! Nb of bands
    !
    ! If we only Wannierze valence bands.
    IF (icbm > nbndsub) icbm = 0
    !
    ! Initialization value. Should be large enough ...
    evbm = -10000d0
    ecbm = 10000d0 ! In Ry
    !
    ! We Wannerize both the CB and VB
    IF (ivbm > 0 .AND. icbm > 0) THEN
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          IF (ibnd < ivbm + 1) THEN
            IF (etf(ibnd, ikk) > evbm) THEN
              evbm = etf(ibnd, ikk)
            ENDIF
          ENDIF
          ! Find cbm index
          IF (ibnd > ivbm) THEN
            IF (etf(ibnd, ikk) < ecbm) THEN
              ecbm = etf(ibnd, ikk)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      ! Find max and min across pools
      CALL mp_max(evbm, inter_pool_comm)
      CALL mp_min(ecbm, inter_pool_comm)
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, "Valence band maximum    = ", f10.6, " eV")') evbm * ryd2ev
        WRITE(stdout, '(5x, "Conduction band minimum = ", f10.6, " eV")') ecbm * ryd2ev
      ENDIF
    ENDIF ! ivbm > 0 .AND. icbm > 0
    !
    ! We only Wannierze the valence bands
    IF (icbm == 0) THEN
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          IF (etf(ibnd, ikk) > evbm) THEN
            evbm = etf(ibnd, ikk)
          ENDIF
        ENDDO
      ENDDO
      ! Find max across pools
      CALL mp_max(evbm, inter_pool_comm)
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, "Valence band maximum    = ", f10.6, " eV")') evbm * ryd2ev
      ENDIF
    ENDIF ! icbm == 0
    !
    ! If we only Wannierized the conduction bands
    IF (ivbm == 0) THEN
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          IF (etf(ibnd, ikk) < ecbm) THEN
            ecbm = etf(ibnd, ikk)
          ENDIF
        ENDDO
      ENDDO
      ! Find min across pools
      CALL mp_min(ecbm, inter_pool_comm)
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, "Conduction band minimum = ", f10.6, " eV")') ecbm * ryd2ev
      ENDIF
    ENDIF ! ivbm == 0
    !
    ! Store e^(e_nk/kbT) on each core
    IF (ivbm > 0) THEN
      DO ik = 1, nkf
        DO ibnd = 1, nbndsub
          ikk = 2 * ik - 1
          ! Because the number are so large. It does lead to instabilities
          ! Therefore we rescale everything to the VBM
          IF (ABS(etemp) < eps80) THEN
            CALL errore('fermicarrier', 'etemp cannot be 0', 1)
          ELSE
            arg = (etf(ibnd, ikk) - evbm) / etemp
          ENDIF
          !
          IF (arg < - maxarg) THEN
            ks_exp(ibnd, ik) = zero
          ELSE
            ks_exp(ibnd, ik) = EXP(arg)
          ENDIF
        ENDDO
      ENDDO
    ENDIF ! ivbm > 0
    !
    ! Store e^(e_nk/kbT) on each core for the electrons (CBM only)
    IF (icbm > 0) THEN
      DO ik = 1, nkf
        DO ibnd = 1, nbndsub
          ikk = 2 * ik - 1
          ! Because the number are so large. It does lead to instabilities
          ! Therefore we rescale everything to the CBM
          arg = (etf(ibnd, ikk) - ecbm) / etemp
          !
          IF (arg > maxarg) THEN
            ks_expcb(ibnd, ik) = 1.0d200
          ELSE
            ks_expcb(ibnd, ik) = EXP(arg)
          ENDIF
        ENDDO
      ENDDO
    ENDIF ! icbm > 0
    !
    ! Case 1 : Intrinsic mobilities (electron and hole concentration are the same)
    ! Starting bounds energy for the biscection method. The energies are rescaled to the VBM
    elw = 1.0d0  ! This is e^0 = 1.0
    eup = 1d-160 ! This is e^(-large) = 0.0 (small)
    IF (int_mob .AND. .NOT. carrier) THEN
      ! Use bisection method
      DO i = 1, maxiter
        !
        !WRITE(stdout,*),'Iteration ',i
        ! We want ef_tmp = (eup + elw) / 2.d0 but the variables are exp therefore:
        ef_tmp = DSQRT(eup) * DSQRT(elw)
        !
        !WRITE(stdout,*),'ef_tmp ', - log (ef_tmp) * etemp * ryd2ev
        hole_density = zero
        electron_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ! Compute hole carrier concentration
          DO ibnd = 1, ivbm
            ! Discard very large numbers
            IF (ks_exp(ibnd, ik) * ef_tmp > 1d60) THEN
              fnk = zero
            ELSE
              fnk = 1.0d0 / (ks_exp(ibnd, ik) * ef_tmp  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk)
          ENDDO
          ! Compute electron carrier concentration
          DO ibnd = icbm, nbndsub
            ! Discard very large numbers
            IF (ks_exp(ibnd, ik) * ef_tmp > 1d60) THEN
              fnk = zero
            ELSE
              fnk = 1.0d0 / (ks_exp(ibnd, ik) * ef_tmp  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            electron_density = electron_density + wkf(ikk) * fnk
          ENDDO
          !
        ENDDO
        !
        CALL mp_sum(hole_density, inter_pool_comm)
        CALL mp_sum(electron_density, inter_pool_comm)
        !
        ! WRITE(stdout,*),'hole_density ',hole_density * (1.0d0/omega) * ( bohr2ang * ang2cm  )**(-3)
        ! WRITE(stdout,*),'electron_density ',electron_density * (1.0d0/omega) * (bohr2ang * ang2cm  )**(-3)
        ! CALL FLUSH(stdout)
        IF (ABS(hole_density) < eps80) THEN
          rel_err = -1000d0
        ELSE
          rel_err = (hole_density - electron_density) / hole_density
        ENDIF
        !
        IF (ABS(rel_err) < eps5) THEN
          fermi_exp = ef_tmp
          fermi = evbm - (LOG(fermi_exp) * etemp)
          EXIT
        ELSEIF ((rel_err) > eps5) THEN
          elw = ef_tmp
        ELSE
          eup = ef_tmp
        ENDIF
      ENDDO ! iteration
    ENDIF
    !
    ! Case 2 :
    ! Hole doped mobilities (Carrier concentration should be larger than 1E5 cm^-3)
    IF (system_2d) THEN
      factor = inv_cell * (bohr2ang * ang2cm)**(-2.d0)
    ELSE
      factor = inv_cell * (bohr2ang * ang2cm)**(-3.d0)
    ENDIF
    eup = 1d-160 ! e^(-large) = 0.0 (small)
    elw = 1.0d40 ! e^0 = 1
    IF (ncarrierp < -1E5 .OR. (int_mob .AND. carrier)) THEN
      IF (int_mob .AND. carrier) ncarrierp = -ABS(ncarrierp)
      ! Use bisection method
      DO i = 1, maxiter
        ! We want ef_tmp = (eup + elw) / 2.d0 but the variables are exp therefore:
        ef_tmp = DSQRT(eup) * DSQRT(elw)
        !
        hole_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ! Compute hole carrier concentration
          DO ibnd = 1, ivbm
            ! Discard very large numbers
            IF (ks_exp(ibnd, ik) * ef_tmp > 1d60) THEN
              fnk = 0.0d0
            ELSE
              fnk = 1.0d0 / (ks_exp(ibnd, ik) * ef_tmp  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
          ENDDO
          !
        ENDDO
        !
        CALL mp_sum(hole_density, inter_pool_comm)
        !
        ! WRITE(stdout,*),'hole_density ',hole_density * (1.0d0/omega) * ( bohr2ang * ang2cm  )**(-3)
        ! CALL FLUSH(stdout)
        IF (ABS(hole_density) < eps80) THEN
          rel_err = -1000.0d0
        ELSE
          ! In this case ncarrier is a negative number
          rel_err = (hole_density - ABS(ncarrierp)) / hole_density
        ENDIF
        !
        IF (ABS(rel_err) < eps5) THEN
          fermi_exp = ef_tmp
          fermi = evbm - (LOG(fermi_exp) * etemp)
          EXIT
        ELSEIF ((rel_err) > eps5) THEN
          elw = ef_tmp
        ELSE
          eup = ef_tmp
        ENDIF
      ENDDO ! iteration
    ENDIF
    !
    ! Case 3 : Electron doped mobilities (Carrier concentration should be larger than 1E5 cm^-3)
    eup = 1.0d-40 ! e^(0) =1
    elw = 1.0d80 ! e^large yields fnk = 1
    IF (ncarrierp > 1E5 .OR. (int_mob .AND. carrier)) THEN
      IF (int_mob .AND. carrier) ncarrierp = ABS(ncarrierp)
      ! Use bisection method
      DO i = 1, maxiter
        ! We want ef_tmp = (eup + elw) / 2.d0 but the variables are exp therefore:
        ef_tmp = DSQRT(eup) * DSQRT(elw)
        !
        electron_density = zero
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ! Compute electron carrier concentration
          DO ibnd = icbm, nbndsub
            ! Discard very large numbers
            IF (ks_expcb(ibnd, ik) * ef_tmp > 1d60) THEN
              fnk = zero
            ELSE
              fnk = 1.0d0 / (ks_expcb(ibnd, ik) * ef_tmp  + 1.0d0)
            ENDIF
            ! The wkf(ikk) already include a factor 2
            electron_density = electron_density + wkf(ikk) * fnk * factor
          ENDDO
          !
        ENDDO
        !
        CALL mp_sum(electron_density, inter_pool_comm)
        !
        IF (ABS(electron_density) < eps80) THEN
          rel_err = 1000.0d0
        ELSE
          ! In this case ncarrier is a negative number
          rel_err = (electron_density - ncarrierp) / electron_density
        ENDIF
        !
        IF (ABS(rel_err) < eps5) THEN
          fermi_exp = ef_tmp
          fermicb = ecbm - (LOG(fermi_exp) * etemp)
          EXIT
        ELSEIF ((rel_err) > eps5) THEN
          eup = ef_tmp
        ELSE
          elw = ef_tmp
        ENDIF
      ENDDO ! iteration
    ENDIF
    !
    IF (i == maxiter) THEN
      WRITE(stdout, '(5x, "Warning: too many iterations in bisection"/ &
                    5x, "ef_tmp = ", f10.6)' ) fermi * ryd2ev
    ENDIF
    !
    ! Print results
    !
    WRITE(stdout, '(/5x, "Temperature ", f8.3, " K")' ) etemp * ryd2ev / kelvin2eV
    !
    ! Small gap semiconductor. Computes intrinsic mobility by placing
    ! the Fermi level such that carrier density is equal for electron and holes
    IF (int_mob .AND. .NOT. carrier) THEN
      !
      ef0(itemp) = fermi
      WRITE(stdout, '(5x, "Mobility Fermi level = ", f10.6, " eV")' )  ef0(itemp) * ryd2ev
      ! We only compute 1 Fermi level so we do not need the other
      efcb(itemp) = zero
      ctype = -1
      !
    ENDIF
    !
    ! Large bandgap semiconductor. Place the gap at the value ncarrier.
    ! The user want both VB and CB mobilities.
    IF (int_mob .AND. carrier) THEN
      !
      ef0(itemp) = fermi
      WRITE(stdout, '(5x, "Mobility VB Fermi level = ", f10.6, " eV")' )  ef0(itemp) * ryd2ev
      !
      efcb(itemp) = fermicb
      WRITE(stdout, '(5x, "Mobility CB Fermi level = ", f10.6, " eV")' )  efcb(itemp) * ryd2ev
      ctype = 0
      !
    ENDIF
    !
    ! User decide the carrier concentration and choose to only look at VB or CB
    IF (.NOT. int_mob .AND. carrier) THEN
      !
      ! VB only
      IF (ncarrierp < 0.0d0) THEN
        ef0(itemp) = fermi
        WRITE(stdout, '(5x, "Mobility VB Fermi level = ", f10.6, " eV")' )  ef0(itemp) * ryd2ev
        ! We only compute 1 Fermi level so we do not need the other
        efcb(itemp) = zero
        ctype = -1
      ELSE ! CB
        efcb(itemp) = fermicb
        WRITE(stdout, '(5x, "Mobility CB Fermi level = ", f10.6, " eV")' )  efcb(itemp) * ryd2ev
        ! We only compute 1 Fermi level so we do not need the other
        ef0(itemp) = 0
        ctype = 1
      ENDIF
    ENDIF
    !
    ! In the case were we do not want mobility (just scattering rates)
    IF (.NOT. int_mob .AND. .NOT. carrier) THEN
      IF (efermi_read) THEN
        !
        ef0(itemp) = fermi_energy
        !
      ELSE !SP: This is added for efficiency reason because the efermig routine is slow
        ef0(itemp) = efnew
      ENDIF
      ! We only compute 1 Fermi level so we do not need the other
      efcb(itemp) = zero
      ctype = -1
      !
    ENDIF
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE fermicarrier
    !-----------------------------------------------------------------------
    !
    !!!!
    !-----------------------------------------------------------------------
    SUBROUTINE calcpartion(itemp, etemp, ctype)
    !-----------------------------------------------------------------------
    !!
    !!  This routine computes the Fermi energy associated with a given
    !!  carrier concentration using bissection for insulators or
    !!  semi-conductors.
    !!
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : omega, alat, at
    USE io_global, ONLY : stdout
    USE elph2,     ONLY : etf, nkf, wkf, efnew, nkqf, partion
    USE constants_epw, ONLY : ryd2ev, bohr2ang, ang2cm, eps5, kelvin2eV, &
                              zero, eps80, cc2cb
    USE noncollin_module, ONLY : noncolin
    USE pwcom,     ONLY : nelec
    USE epwcom,    ONLY : int_mob, nbndsub, ncarrier, nstemp, fermi_energy, &
                          system_2d, carrier, efermi_read, assume_metal, ngaussw, &
                          ii_eda, ii_n
    USE klist_epw, ONLY : isk_dummy
    USE mp,        ONLY : mp_barrier, mp_sum, mp_max, mp_min
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    INTEGER, INTENT(out) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in kBT [Ry] unit.
    !
    ! Local variables
    INTEGER :: i
    !! Index for the bisection iteration
    INTEGER :: ik
    !! k-point index per pool
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ivbm
    !! Index of the VBM
    INTEGER :: icbm
    !! Index of the CBM
    INTEGER, PARAMETER :: maxiter = 500 ! 300
    !! Maximum interation
    REAL(KIND = DP) :: fermi
    !! Fermi level returned
    REAL(KIND = DP) :: fermicb
    !! Fermi level returned for second Fermi level
    REAL(KIND = DP) :: fnk
    !! Fermi-Diract occupation
    REAL(KIND = DP) :: ks_exp(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT
    REAL(KIND = DP) :: ks_expcb(nbndsub, nkf)
    !! Exponential of the eigenvalues divided by kBT for CB
    REAL(KIND = DP) :: fermi_exp
    !! Fermi level in exponential format
    REAL(KIND = DP) :: rel_err
    !! Relative error
    REAL(KIND = DP) :: factor
    !! Factor that goes from number of carrier per unit cell to number of
    !! carrier per cm^-3
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    REAL(KIND = DP) :: inv_cell
    !! Inverse of the volume in [Bohr^{-3}]
    REAL(KIND = DP) :: evbm
    !! Energy of the VBM
    REAL(KIND = DP) :: ecbm
    !! Energy of the CBM
    REAL(KIND = DP) :: ef_tmp
    !! Energy of the current Fermi level for the bisection method
    REAL(KIND = DP) :: elw
    !! Energy lower bound for the bisection method
    REAL(KIND = DP) :: eup
    !! Energy upper bound for the bisection method
    REAL(KIND = DP) :: hole_density
    !! Hole carrier density
    REAL(KIND = DP) :: electron_density
    !! Electron carrier density
    REAL(KIND = DP), PARAMETER :: maxarg = 200.d0
    !! Maximum value for the argument of the exponential
    REAL(KIND = DP) :: eimp
    !! Ionization energy of the impurity state, relative to band energies
    REAL(KIND = DP) :: dummy1
    !! Dummy var, intermediate for sums
    REAL(KIND = DP) :: impurity_partion
    !! density of ionized impurities
    REAL(KIND = DP) :: impurity_density
    !! density of dopants in per cubic Bohr
    !
    ! Convert impurity density to per cubic Bohr
    impurity_density = ii_n * omega / cc2cb
    !
    ! Decide ctype
    ctype = 0
    IF (ncarrier > 0.0d0) THEN
      ctype = 1
    ELSEIF (ncarrier < 0.0d0) THEN
      ctype = -1
    ENDIF
    IF (int_mob) THEN
      ctype = 0
    ENDIF
    !
    IF (ctype == 0) THEN
      WRITE(stdout, '(5x, "Warning! partial ionization not implemented for ctype=0")')
      RETURN
    ENDIF
    IF (assume_metal) THEN
      !! set conduction band chemical potential to 0 since it is irrelevent
      !!ctype = -1  ! act like it's for holes
      !!efcb(itemp) = zero
      !!ef0(itemp) = efermig(etf, nbndsub, nkqf, nelec, wkf, etemp, ngaussw, 0, isk_dummy)
      WRITE(stdout, '(5x, "Warning! No partial ionization for metallic systems, exiting")')
      RETURN
    ENDIF
    !ef_tmp  = zero
    !fermi   = zero
    !fermicb = zero
    inv_cell = 1.0d0 / omega
    !
    ! for 2d system need to divide by area (vacuum in z-direction)
    IF (system_2d) inv_cell = inv_cell * at(3, 3) * alat
    ! vbm index
    IF (noncolin) THEN
      ivbm = FLOOR(nelec / 1.0d0)
    ELSE
      ivbm = FLOOR(nelec / 2.0d0)
    ENDIF
    icbm = ivbm + 1 ! Nb of bands
    !
    ! If we only Wannierze valence bands.
    IF (icbm > nbndsub) icbm = 0
    !
    ! Initialization value. Should be large enough ...
    evbm = -10000d0
    ecbm = 10000d0 ! In Ry
    !
    ! We Wannerize both the CB and VB
    IF (ivbm > 0 .AND. icbm > 0) THEN
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          IF (ibnd < ivbm + 1) THEN
            IF (etf(ibnd, ikk) > evbm) THEN
              evbm = etf(ibnd, ikk)
            ENDIF
          ENDIF
          ! Find cbm index
          IF (ibnd > ivbm) THEN
            IF (etf(ibnd, ikk) < ecbm) THEN
              ecbm = etf(ibnd, ikk)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      ! Find max and min across pools
      CALL mp_max(evbm, inter_pool_comm)
      CALL mp_min(ecbm, inter_pool_comm)
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, "Valence band maximum    = ", f10.6, " eV")') evbm * ryd2ev
        WRITE(stdout, '(5x, "Conduction band minimum = ", f10.6, " eV")') ecbm * ryd2ev
      ENDIF
    ENDIF ! ivbm > 0 .AND. icbm > 0
    !
    ! We only Wannierze the valence bands
    IF (icbm == 0) THEN
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          IF (etf(ibnd, ikk) > evbm) THEN
            evbm = etf(ibnd, ikk)
          ENDIF
        ENDDO
      ENDDO
      ! Find max across pools
      CALL mp_max(evbm, inter_pool_comm)
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, "Valence band maximum    = ", f10.6, " eV")') evbm * ryd2ev
      ENDIF
    ENDIF ! icbm == 0
    !
    ! If we only Wannierized the conduction bands
    IF (ivbm == 0) THEN
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          IF (etf(ibnd, ikk) < ecbm) THEN
            ecbm = etf(ibnd, ikk)
          ENDIF
        ENDDO
      ENDDO
      ! Find min across pools
      CALL mp_min(ecbm, inter_pool_comm)
      IF (itemp == 1) THEN
        WRITE(stdout, '(5x, "Conduction band minimum = ", f10.6, " eV")') ecbm * ryd2ev
      ENDIF
    ENDIF ! ivbm == 0
    !
    IF (ctype == -1) THEN
      dummy1 = 0.0d0
      ! Calc impuirty energy relative to bands
      eimp = evbm + ii_eda
      ! Now, loop through temp, ik, and ibnd
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          dummy1 = dummy1 + (1.0d0 - (1.0d0 / (EXP((etf(ibnd, ikk)-eimp)/(etemp))+1.0d0)))*wkf(ikk)
        ENDDO
      ENDDO
      CALL mp_sum(dummy1, inter_pool_comm)
      !! now calculate fraction of ionized impurities in partion(itemp)
      partion(itemp) = SQRT(impurity_density*dummy1) / impurity_density
      !! equation will give values greater than 1, 
      !! enforce fraction == 1.0 for full ionization
      IF ( partion(itemp) > 1.0d0 ) THEN
        partion(itemp) = 1.0d0
      ENDIF
      !WRITE(stdout, '(/5x, "Temperature ", f8.3, " K")' ) etemp * ryd2ev / kelvin2eV
      WRITE(stdout, '(/5x, "Fraction of ionized impurities at ", f8.3," K: ", f16.12)' ) etemp * ryd2ev / kelvin2eV, partion(itemp)
    ENDIF
    !
    IF (ctype == 1) THEN
      dummy1 = 0.0d0
      ! Calc impuirty energy relative to bands
      eimp = ecbm - ii_eda
      ! Now, loop through temp, ik, and ibnd
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndsub
          dummy1 = dummy1 + (1.0d0 / (EXP((etf(ibnd, ikk)-eimp)/(etemp))+1.0d0))*wkf(ikk)
        ENDDO
      ENDDO
      CALL mp_sum(dummy1, inter_pool_comm)
      !! now calculate fraction of ionized impurities in partion(itemp)
      partion(itemp) = SQRT(impurity_density*dummy1) / impurity_density
      !! equation will give values greater than 1, 
      !! enforce fraction == 1.0 for full ionization
      IF ( partion(itemp) > 1.0d0 ) THEN
        partion(itemp) = 1.0d0
      ENDIF
      !WRITE(stdout, '(/5x, "Temperature ", f8.3, " K")' ) etemp * ryd2ev / kelvin2eV
      WRITE(stdout, '(/5x, "Fraction of ionized impurities at ", f8.3," K: ", f16.12)' ) etemp * ryd2ev / kelvin2eV, partion(itemp)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE calcpartion
    !-----------------------------------------------------------------------
    !!!!!
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fermiwindow()
    !-----------------------------------------------------------------------
    !!
    !! Find the band indices of the first
    !! and last state falling within the window e_fermi+-efermithickness
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf
    USE epwcom,        ONLY : fsthick, nbndsub
    USE pwcom,         ONLY : ef
    USE mp,            ONLY : mp_max, mp_min
    USE mp_global,     ONLY : inter_pool_comm
    USE constants_epw, ONLY : ryd2ev
    !
    IMPLICIT NONE
    !
    INTEGER :: ik
    !! Counter on k-points in the pool
    INTEGER :: ibnd
    !! Counter on bands
    REAL(KIND = DP) :: ebnd
    !! Eigenvalue at etf(ibnd, ik)
    REAL(KIND = DP) :: ebndmin
    !! Minimum eigenvalue
    REAL(KIND = DP) :: ebndmax
    !! Maximum eigenvalue
    REAL(KIND = DP) :: tmp
    !
    !
    ibndmin = 100000
    ibndmax = 0
    ebndmin =  1.d8
    ebndmax = -1.d8
    !
    DO ik = 1, nkqf
      DO ibnd = 1, nbndsub
        ebnd = etf(ibnd, ik)
        !
        IF (ABS(ebnd - ef) < fsthick) THEN
          ibndmin = MIN(ibnd, ibndmin)
          ibndmax = MAX(ibnd, ibndmax)
          ebndmin = MIN(ebnd, ebndmin)
          ebndmax = MAX(ebnd, ebndmax)
        ENDIF
        !
      ENDDO
    ENDDO
    !
    tmp = DBLE(ibndmin)
    CALL mp_min(tmp, inter_pool_comm)
    ibndmin = NINT(tmp)
    CALL mp_min(ebndmin, inter_pool_comm)
    !
    tmp = DBLE(ibndmax)
    CALL mp_max(tmp, inter_pool_comm)
    ibndmax = NINT(tmp)
    CALL mp_max(ebndmax, inter_pool_comm)
    !
    WRITE(stdout,'(/14x,a,i5,2x,a,f9.3,a)') 'ibndmin = ', ibndmin, 'ebndmin = ', ebndmin * ryd2ev, ' eV'
    WRITE(stdout,'(14x,a,i5,2x,a,f9.3,a/)') 'ibndmax = ', ibndmax, 'ebndmax = ', ebndmax * ryd2ev, ' eV'
    !
    !----------------------------------------------------------------------
    END SUBROUTINE fermiwindow
    !---------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  END MODULE utilities
  !-----------------------------------------------------------------------------
