  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .         
  !                                                                            
  ! Adapted from QE  
  !
  !
  !----------------------------------------------------------------------
  MODULE broyden
  !----------------------------------------------------------------------
  !! 
  !! This module contains the routines associated with Broyden's method 
  !! and Pade' approximants
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE mix_broyden(ndim, deltaout, deltain, alphamix, iter, n_iter, conv, df, dv)
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nsiter
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
#if ! defined(__GFORTRAN__) || (__GNUC__ > 4 )
    USE, INTRINSIC :: IEEE_ARITHMETIC
#endif
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
#if ! defined(__GFORTRAN__) || (__GNUC__ > 4 )
      IF (IEEE_IS_NAN(ar) .OR. IEEE_IS_NAN(ai)) THEN
        !WRITE(stdout, *) (z(i), i = 1, N)
        !WRITE(stdout, *) (u(i), i = 1, N)
        !WRITE(stdout, *) (a(i), i = 1, N)
        WRITE(stdout, *) 'One or more Pade coefficients are NaN'
        !CALL errore('pade_coeff', 'one or more coefficients are NaN', 1)
      ENDIF
#endif
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
  !-----------------------------------------------------------------------------
  END MODULE broyden
  !-----------------------------------------------------------------------------
