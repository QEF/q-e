  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2001-2008 Quantum-Espresso group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE rigid_epw
  !----------------------------------------------------------------------
  !!
  !! This module contains routine linked with the calculation of the rigid-ion
  !! (long-range) term for q.
  !!
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  SAVE

  INTEGER :: igmin(3), igmin_qG(3)
  !!
  REAL(KIND = DP) :: qqcut
  !!
  REAL(KIND = DP), PARAMETER :: alph = 1.d0
  !! Ewald parameter, units (2pi/alat)^{2}
  REAL(KIND = DP), PARAMETER :: gmax = 14.d0
  !! Cutoff criteria for G-sum:
  !! e^{ - (q+G) \cdot \epsilon \cdot (q+G) / (4 \alpha) } < e^{ - gmax } terms are neglected.
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------
    COMPLEX(KIND = DP) FUNCTION H_eps(z)
    !--------------------------------------------------------------------------
    !!
    !! Function used in the Lindhard function. See Eq.(56) of Hedin (1965)
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : eps10
    !
    IMPLICIT NONE
    !
    COMPLEX(KIND = DP), INTENT(in) :: z
    !! Argument of the Lindhard function
    !
    IF (ABS(z - 1.d0) > eps10) THEN
      IF (ABS((z + 1.d0) / (z - 1.d0)) > eps10) THEN
        H_eps = 2.d0 * z + (1.d0 - z**2) * LOG((z + 1.d0) / (z - 1.d0))
      ENDIF
    ENDIF
    !
    RETURN
    !
    !--------------------------------------------------------------------------
    END FUNCTION H_eps
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_min_qG(q)
    !-----------------------------------------------------------------------
    !!
    !! Find igmin_qG = min_{G}|q+G|
    !! to shift G-sum in long-range g
    !! in order to ensure periodicity g(q+G')=g(q)
    !!
    USE constants_epw, ONLY : eps8
    USE cell_base,     ONLY : bg, at
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid, in crystal coords.
    !
    ! Local variables
    INTEGER         :: m1
    !
    INTEGER         :: m2
    !
    INTEGER         :: m3
    !! Loop over G-vectors
    REAL(KIND = DP) :: g1
    !!
    REAL(KIND = DP) :: g2
    !!
    REAL(KIND = DP) :: g3
    !! G-vector in cartesian coords.
    REAL(KIND = DP) :: qq
    !! |q+G|
    REAL(KIND = DP) :: qqmin
    !! min|q+G|
    REAL(KIND = DP) :: qtmp(3)
    !!
    !
    ! Move q to 1BZ
    qtmp(:) = q(:) - INT(q(:))
    ! Transform to cartesian coords
    CALL cryst_to_cart(1, qtmp, bg, 1)
    !
    igmin_qG = 0
    qqmin = 1E10
    DO m1 = -2, 2
      DO m2 = -2, 2
        DO m3 = -2, 2
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + qtmp(1)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + qtmp(2)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + qtmp(3)
          qq = g1 * g1 + g2 * g2 + g3 * g3
          IF (qqmin > qq) THEN
            qqmin = qq
            igmin_qG(1) = m1 - INT(q(1))
            igmin_qG(2) = m2 - INT(q(2))
            igmin_qG(3) = m3 - INT(q(3))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE find_min_qG
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_gmin(q)
    !-----------------------------------------------------------------------
    !!
    !! Find the index of G vectors for the minimum distance and
    !! the cut-off for the complete first shell.
    !!
    USE constants_epw, ONLY : eps8
    USE cell_base,     ONLY : bg
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    !
    ! Local variables
    INTEGER         :: m1
    !! Loop over q-points
    INTEGER         :: m2
    !! Loop over q-points
    INTEGER         :: m3
    !! Loop over q-points
    REAL(KIND = DP) :: g1
    !!
    REAL(KIND = DP) :: g2
    !!
    REAL(KIND = DP) :: g3
    !!
    REAL(KIND = DP) :: qq
    !!
    REAL(KIND = DP) :: qqmin
    !!
    REAL(KIND = DP) :: qtmp(3)
    !!
    !
    qtmp(:) = q(:) - INT(q(:))
    CALL cryst_to_cart(1, qtmp, bg, 1)
    qqmin = 1E10
    DO m1 = -2, 2
      DO m2 = -2, 2
        DO m3 = -2, 2
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) - qtmp(1)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) - qtmp(2)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) - qtmp(3)
          qq = g1 * g1 + g2 * g2 + g3 * g3
          IF (qqmin > qq) THEN
            qqmin = qq
            igmin(1) = m1 + INT(q(1))
            igmin(2) = m2 + INT(q(2))
            igmin(3) = m3 + INT(q(3))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    qqcut = -1E10
    DO m1 = 1, 3
      qqcut = MAX(qqcut, SUM(bg(:, m1)**2))
    ENDDO
    qqcut = qqcut + eps8
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE find_gmin
    !--------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rgd_blk(nqc1, nqc2, nqc3, nat, dyn, q, tau, epsil, zeu, signe)
    !-----------------------------------------------------------------------
    !! This is adapted from QE PH/rigid.f90
    !!
    !! compute the rigid-ion (long-range) term for q
    !! The long-range term used here, to be added to or subtracted from the
    !! dynamical matrices, is exactly the same of the formula introduced
    !! in:
    !! X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is
    !! implemented: the Ewald parameter alpha must be large enough to
    !! have negligible r-space contribution
    !!
    !! This implements Eq. 98 of Rev. Mod. Phys., 73, 515 (2001)
    !! SP: 04/2019 - Using nrx1 is overkill.
    !! SP: 11/2019 - Addition of system_2d (we assume z is the vacuum direction).
    !! SP: 08/2020 - Restoration of nrx and parallelization.
    !! SP: 04/2021 - Addition of quadrupoles
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi, fpi, e2
    USE cell_base,     ONLY : bg, omega, alat, at
    USE constants_epw, ONLY : eps6, ci, zero, czero, twopi, eps8
    USE io_global,     ONLY : ionode_id
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : world_comm
    USE division,      ONLY : para_bounds
    USE mp,            ONLY : mp_bcast, mp_sum
    USE epwcom,        ONLY : lpolar, system_2d
    USE elph2,         ONLY : area, Qmat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL(KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! effective charges tensor
    REAL(KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract rigid-ion term
    REAL(KIND = DP), INTENT(in) :: tau(3, nat)
    !! Atomic positions
    COMPLEX(KIND = DP), INTENT(inout) :: dyn(3 * nat, 3 * nat)
    !! Dynamical matrix
    !
    ! Local variables
    INTEGER :: na
    !! Atom index 1
    INTEGER :: nb
    !! Atom index 2
    INTEGER :: i
    !! Cartesian direction 1
    INTEGER :: j
    !! Cartesian direction 1
    INTEGER :: ipol, jpol, kpol, lpol
    !! Polarization direction
    INTEGER :: mm, m1, m2, m3
    !! Loop over q-points
    INTEGER :: mmax
    !! Max = nr1x * nr2x * nr3x
    INTEGER :: nr1x, nr2x, nr3x
    !! Minimum supercell size to include all vector such that G^2 < geg
    INTEGER :: mm_start
    !! Starting ir for this cores
    INTEGER :: mm_stop
    !! Ending ir for this pool
    INTEGER :: diff
    !! Difference between starting and ending on master core
    INTEGER :: add
    !! Additional element
    REAL(KIND = DP):: metric
    !! (2*pi/a)^2
    REAL(KIND = DP):: geg
    !! <q+G| epsil | q+G>
    !REAL(KIND = DP) :: alph
    !!! Ewald parameter
    REAL(KIND = DP) :: fac
    !! General prefactor
    REAL(KIND = DP) :: gg(3)
    !! G-vectors
    REAL(KIND = DP) :: facgd
    !! fac * EXP(-geg / (alph * 4.0d0)) / geg
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    !REAL(KIND = DP) :: gmax
    !!! Maximum G
    REAL(KIND = DP) :: zag(3)
    !! Z * G
    REAL(KIND = DP) :: qag(3)
    !! Q * G
    REAL(KIND = DP) :: zcg(3)
    !! Z * G
    REAL(KIND = DP) :: qcg(3)
    !! Q * G
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    REAL(KIND = DP) :: qtmp(3)
    !! Temporary q vector to find min_{G}|q+G|
    COMPLEX(KIND = DP) :: fnat(3)
    !! Z with \delta_kk' summed
    COMPLEX(KIND = DP) :: qnat(3)
    !! Q with \delta_kk' summed
    REAL(KIND = DP) :: reff(2, 2)
    !! Effective screening length for 2D materials
    REAL(KIND = DP) :: grg
    !! G-vector * reff * G-vector
    COMPLEX(KIND = DP) :: Qqq
    !! Quadrupole-quadrupole term
    COMPLEX(KIND = DP) :: Qdq
    !! Quadrupole-dipole term
    COMPLEX(KIND = DP) :: Qdd
    !! Dipole-dipole term
    COMPLEX(KIND = DP) :: facg
    !! Atomic position exponential
    COMPLEX(KIND = DP) :: dyn_tmp(3 * nat, 3 * nat)
    !! Temporary dyn. matrice
    !
    ! In case you have a non polar materials but you want quadrupole effects
    ! Enforce 0 Born effective charges
    IF (.NOT. lpolar) THEN
      zeu(:, :, :) = zero
    ENDIF
    Qqq = zero
    Qdq = zero
    Qdd = zero
    metric = (twopi / alat)**2
    !
    ! alph is the Ewald parameter, geg is an estimate of G^2
    ! such that the G-space sum is convergent for that alph
    ! very rough estimate: geg/4/alph > gmax = 14
    ! (exp (-14) = 10^-6)
    !
    IF (ABS(ABS(signe) - 1.0) > eps6) CALL errore('rgd_blk', ' wrong value for signe ', 1)
    !
    IF (system_2d) THEN
      ! Vacuum size in Bohr unit
      c = alat / bg(3, 3)
      ! (e^2 * 2\pi) / Area
      fac = (signe * e2 * twopi) / area
      ! Effective screening length
      ! reff = (epsil - 1) * c/2
      reff(:, :) = zero
      reff(:, :) = epsil(1:2, 1:2) * 0.5d0 * c ! eps * c/2
      reff(1, 1) = reff(1, 1) - 0.5d0 * c ! (-1) * c/2
      reff(2, 2) = reff(2, 2) - 0.5d0 * c ! (-1) * c/2
    ELSE
      ! (e^2 * 4\pi) / Volume
      fac = (signe * e2 * fpi) / omega
    ENDIF
    !
    !gmax = 14.d0
    !alph = 1.0d0
    geg = gmax * alph * 4.0d0
    !
    ! Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
    ! Only for dimensions where periodicity is present, e.g. if nr1=1
    ! and nr2=1, then the G-vectors run along nr3 only.
    ! (useful if system is in vacuum, e.g. 1D or 2D)
    IF (nqc1 == 1) THEN
      nr1x = 0
    ELSE
      nr1x = INT(SQRT(geg) / SQRT(bg(1, 1)**2 + bg(2, 1)**2 + bg(3, 1)**2)) + 1
    ENDIF
    IF (nqc2 == 1) THEN
      nr2x = 0
    ELSE
      nr2x = INT(SQRT(geg) / SQRT(bg(1, 2)**2 + bg(2, 2)**2 + bg(3, 2)**2)) + 1
    ENDIF
    IF (nqc3 == 1) THEN
      nr3x = 0
    ELSE
      nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
    ENDIF
    !
    mmax = (2 * nr1x + 1) * (2 * nr2x + 1) * (2 * nr3x + 1)
    !
    ! Distribute the cpu
    CALL para_bounds(mm_start, mm_stop, mmax)
    !
    IF (mpime == ionode_id) THEN
      diff = mm_stop - mm_start
    ENDIF
    CALL mp_bcast(diff, ionode_id, world_comm)
    !
    ! If you are the last cpu with less element
    IF (mm_stop - mm_start /= diff) THEN
      add = 1
    ELSE
      add = 0
    ENDIF
    !
    !JLB: Find igmin_qG = min_{G}|q+G|
    qtmp=q
    CALL cryst_to_cart(1, qtmp, at, -1)
    CALL find_min_qG(qtmp)
    !
    dyn_tmp(:, :) = czero
    !
    ! DO mm = 1, mmax
    DO mm = mm_start, mm_stop + add
      IF (add == 1 .AND. mm == mm_stop + add) CYCLE
      !
      m1 = -nr1x + FLOOR(1.0d0 * (mm - 1) / ((2 * nr3x + 1) * (2 * nr2x + 1)))
      m2 = -nr2x + MOD(FLOOR(1.0d0 * (mm - 1) / (2 * nr3x + 1)), (2 * nr2x + 1))
      m3 = -nr3x + MOD(1.0d0 * (mm - 1), 1.0d0 * (2 * nr3x + 1))
      !
      ! Special case of q = 0
      gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1,3)) * (twopi / alat)
      gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2,3)) * (twopi / alat)
      gg(3) = (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3,3)) * (twopi / alat)
      !
      IF (system_2d) THEN
        geg = gg(1)**2 + gg(2)**2 + gg(3)**2
        grg = zero
        IF (gg(1)**2 + gg(2)**2 > eps8) THEN
          grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
          grg = grg / (gg(1)**2 + gg(2)**2)
        ENDIF
      ELSE
        !
        geg = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
               gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
               gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
      ENDIF
      !
      IF (geg > 0.0d0 .AND. geg / (metric * alph * 4.0d0) < gmax) THEN
        !
        IF (system_2d) THEN
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / SQRT(geg) / (1.0 + grg * SQRT(geg))
        ELSE
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / geg
        ENDIF
        !
        DO na = 1, nat
          zag(:) = zero
          DO i = 1, 3
            DO ipol = 1, 3
              zag(i) = zag(i) + gg(ipol) * zeu(ipol, i, na)
            ENDDO
          ENDDO
          !
          qag(:) = zero
          DO i = 1, 3
            DO ipol = 1, 3
              DO jpol = 1, 3
                qag(i) = qag(i) + gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
              ENDDO
            ENDDO
          ENDDO
          !
          fnat(:) = zero
          qnat(:) = zero
          DO nb = 1, nat
            arg = alat * (gg(1) * (tau(1, na) - tau(1, nb)) + &
                          gg(2) * (tau(2, na) - tau(2, nb)) + &
                          gg(3) * (tau(3, na) - tau(3, nb)))
            ! Dipole-dipole
            zcg(:) = zero
            DO j = 1, 3 ! Cartesian direction alpha
              DO jpol = 1, 3 !
                zcg(j) = zcg(j) + gg(jpol) * zeu(jpol, j, nb)
              ENDDO
              fnat(j) = fnat(j) + zcg(j) * CMPLX(COS(arg), SIN(arg), KIND = DP)
            ENDDO
            !
            qcg(:) = zero
            DO j = 1, 3 ! Cartesian direction alpha
              DO jpol = 1, 3
                DO kpol = 1, 3
                  qcg(j) = qcg(j) + gg(jpol) * gg(kpol) *  Qmat(nb, j, jpol, kpol)
                ENDDO
              ENDDO
              qnat(j) = qnat(j) + qcg(j) * CMPLX(COS(arg), SIN(arg), KIND = DP)
            ENDDO
            !
          ENDDO ! nb
          DO j = 1, 3
            DO i = 1, 3
              ! Dipole-dipole
              Qdd = zag(i) * fnat(j)
              ! Dipole-quad
              !Qdq = 0.5d0 * (zag(i) * qnat(j) - fnat(i) * qag(i))
              Qdq = 0.5d0 * (qag(i) * fnat(j) - zag(i) * qnat(j))
              ! Quad-quad
              Qqq = 0.25d0 * qag(i) * qnat(j)
              !
              dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) &
                                           - facgd * (Qdd + ci * Qdq + Qqq)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! nat
      ENDIF ! geg
      !
      ! Case q =/ 0
      !gg(1) = gg(1) + q(1) * (twopi / alat)
      !gg(2) = gg(2) + q(2) * (twopi / alat)
      !gg(3) = gg(3) + q(3) * (twopi / alat)
      !
      !JLB: shift G-sum and center around min_{G}|q+G|, to ensure periodicity
      gg(1) = ( (m1+igmin_qG(1)) * bg(1, 1) + (m2+igmin_qG(2)) * bg(1, 2) + (m3+igmin_qG(3)) * bg(1,3) + q(1)) * (twopi / alat)
      gg(2) = ( (m1+igmin_qG(1)) * bg(2, 1) + (m2+igmin_qG(2)) * bg(2, 2) + (m3+igmin_qG(3)) * bg(2,3) + q(2)) * (twopi / alat)
      gg(3) = ( (m1+igmin_qG(1)) * bg(3, 1) + (m2+igmin_qG(2)) * bg(3, 2) + (m3+igmin_qG(3)) * bg(3,3) + q(3)) * (twopi / alat)
      !
      !
      IF (system_2d) THEN
        geg = gg(1)**2 + gg(2)**2 + gg(3)**2
        grg = zero
        IF (gg(1)**2 + gg(2)**2 > eps8) THEN
          grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
          grg = grg / (gg(1)**2 + gg(2)**2)
        ENDIF
      ELSE
        geg = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
               gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
               gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
      ENDIF
      !
      IF (geg > 0.0d0 .AND. geg / (metric * alph * 4.0d0) < gmax) THEN
        !
        IF (system_2d) THEN
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / (SQRT(geg) * (1.0 + grg * SQRT(geg)))
        ELSE
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / geg
        ENDIF
        !
        DO nb = 1, nat ! kappa
          DO na = 1, nat ! kappa'
            arg = alat * (gg(1) * (tau(1, na) - tau(1 ,nb)) + &
                          gg(2) * (tau(2, na) - tau(2, nb)) + &
                          gg(3) * (tau(3, na) - tau(3, nb)) )
            !
            facg = facgd * CMPLX(COS(arg), SIN(arg), DP)
            !
            DO j = 1, 3 ! Cartesian direction alpha
              DO i = 1, 3 ! Carestian direction beta
                ! Dipole - dipole term
                Qdd = zero
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    Qdd = Qdd + gg(ipol) * zeu(ipol, j, nb) * gg(jpol) * zeu(jpol, i, na)
                  ENDDO
                ENDDO
                !
                ! Dipole - quadrupole term
                Qdq = zero
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    DO kpol = 1, 3
                      Qdq = Qdq + 0.5 * (gg(ipol) * zeu(ipol, j, nb) * gg(jpol) * gg(kpol) * Qmat(na, i, jpol, kpol) &
                                       - gg(ipol) * gg(jpol) * Qmat(nb, j, ipol, jpol) * gg(kpol) * zeu(kpol, i, na))
                    ENDDO
                  ENDDO
                ENDDO
                !
                ! Quadrupole - quadrupole term
                Qqq = zero
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    DO kpol = 1, 3
                      DO lpol = 1, 3
                        Qqq = Qqq + 0.25 *  gg(ipol) * gg(jpol) * Qmat(nb, j, ipol, jpol) * &
                                            gg(kpol) * gg(lpol) * Qmat(na, i, kpol, lpol)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                !
                dyn_tmp((na - 1) * 3 + i, (nb - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (nb - 1) * 3 + j) &
                                             + facg * (Qdd + ci * Qdq + Qqq)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! na
        ENDDO ! nb
      ENDIF
    ENDDO ! mm
    !
    CALL mp_sum(dyn_tmp, world_comm)
    dyn(:, :) = dyn_tmp(:, :) + dyn(:, :)
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE rgd_blk
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE rgd_blk_epw(nqc1, nqc2, nqc3, q, uq, epmat, nmodes, epsil, zeu, bmat, signe)
    !-------------------------------------------------------------------------------
    !!
    !! Compute the long range term for the e-ph vertex
    !! to be added or subtracted from the vertex
    !!
    !! The long-range part can be computed using Eq. (4) of PRL 115, 176401 (2015).
    !! The sum over G is converged using the Ewald summation technique (see for example
    !! F.2, p.500 in Martin Electronic structure book) where the Ewald factor is ((q+G)**2)/alph/4.0d0.
    !!
    !! Technical note: From the solution of the Poisson equation, there is an additional factor
    !! e^{-i(q+G)\tau_\kappa} with respect to Eq. (4) of PRL 115, 176401 (2015).
    !! The full equation can be found in Eq. (S4) of the supplemental materials of PRL 115, 176401 (2015).
    !! In practical calculations the G-vector sum is restricted to small |q + G| via the cutoff
    !! function e^{-({\bf q}+{\bf G})^2/4\alpha}. See footnote 7 of RMP 89, 015003 (2017).
    !!
    !!
    !! The final implemented formula is:
    !!
    !! $$ g_{mn\nu}^{\mathcal L}({\bf k},{\bf q) = i\frac{4\pi e^2}{\Omega} \sum_{\kappa}
    !!   \left(\frac{\hbar}{2 {M_\kappa \omega_{{\bf q}\nu}}}\right)^{\!\!\frac{1}{2}}
    !!   \sum_{{\bf G}\ne -{\bf q}} e^{-({\bf q}+{\bf G})^2/4\alpha}
    !! \frac{ ({\bf q}+{\bf G})\cdot{\bf Z}^*_\kappa \cdot {\bf e}_{\kappa\nu}({\bf q}) }
    !!  {({\bf q}+{\bf G})\cdot\bm\epsilon^\infty\!\cdot({\bf q}+{\bf G})}\,
    !!   \left[ U_{{\bf k}+{\bf q}}\:U_{{\bf k}}^{\dagger} \right]_{mn} $$
    !!
    USE kinds,         ONLY : dp
    USE cell_base,     ONLY : bg, omega, alat, at
    USE ions_base,     ONLY : tau, nat
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12, zero, eps8
    USE epwcom,        ONLY : shortrange, lpolar, system_2d
    USE elph2,         ONLY : area, Qmat
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nmodes
    !! Max number of modes
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL(KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! effective charges tensor
    REAL(KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract long range term
    COMPLEX(KIND = DP), INTENT(in) :: uq(nmodes, nmodes)
    !! phonon eigenvec associated with q
    COMPLEX(KIND = DP), INTENT(inout) :: epmat(nmodes)
    !! e-ph matrix elements
    COMPLEX(KIND = DP), INTENT(in) :: bmat
    !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
    !
    ! Local variables
    INTEGER :: na
    !! Atom index 1
    INTEGER :: i
    !! Cartesian index direction
    INTEGER :: ipol, jpol
    !! Polarison direction
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: nr1x, nr2x, nr3x
    !! Minimum supercell size to include all vector such that G^2 < geg
    INTEGER :: mmin(3), mmax(3)
    !! Shifted G-loop to be centered around min_{G}|q+G|
    INTEGER :: nGtest
    !! Number of G-vectors within cutoff (for testing purposes only)
    REAL(KIND = DP):: metric
    !! (2*pi/a)^2
    REAL(KIND = DP) :: qeq
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: arg
    !! Argument of the cos,sin for the Euler formula
    REAL(KIND = DP) :: zaq
    !! Z^* \cdot (q+g)
    REAL(KIND = DP) :: gg(3)
    !! G-vector
    !REAL(KIND = DP) :: gmax
    !!! Max G-vector
    !REAL(KIND = DP) :: alph
    !!! Ewald factor (arbitrary, here chosen to be 1)
    REAL(KIND = DP) :: geg
    !!  <q+G| epsil | q+G>
    REAL(KIND = DP) :: Qqq
    !! In the case of Si, its a single value
    REAL(KIND = DP) :: reff(2, 2)
    !! Effective screening length for 2D materials
    REAL(KIND = DP) :: grg
    !! G-vector * reff * G-vector
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    REAL(KIND = DP) :: qtmp(3)
    !! Temporary q vector to find min_{G}|q+G|
    COMPLEX(KIND = DP) :: fac
    !! General prefactor
    COMPLEX(KIND = DP) :: facqd
    !! Exp function
    COMPLEX(KIND = DP) :: facq
    !! Atomic position exponential
    COMPLEX(KIND = DP) :: epmatl(nmodes)
    !! Long-range part of the el-ph matrix elements
    !
    ! If non-polar materials with quadrupoles
    IF (.NOT. lpolar) THEN
      zeu(:, :, :) = zero
    ENDIF
    !
    IF(ABS(ABS(signe) - 1.0) > eps12) CALL errore('rgd_blk_epw', 'Wrong value for signe ', 1)
    !
    IF (system_2d) THEN
      ! Vacuum size in Bohr unit
      c = alat / bg(3, 3)
      ! (e^2 * 2\pi * ci) / Area
      fac = (signe * e2 * twopi * ci) / area
      ! Effective screening length
      ! reff = (epsil - 1) * c/2
      reff(:, :) = zero
      reff(:, :) = epsil(1:2, 1:2) * 0.5d0 * c ! eps * c/2
      reff(1, 1) = reff(1, 1) - 0.5d0 * c ! (-1) * c/2
      reff(2, 2) = reff(2, 2) - 0.5d0 * c ! (-1) * c/2
    ELSE
      ! (e^2 * 4\pi * i) / Volume
      fac = (signe * e2 * fpi * ci) / omega
    ENDIF
    !
    !gmax = 14.d0
    !alph = 1.0d0
    metric = (twopi / alat)**2
    geg = gmax * alph * 4.0d0
    !
    ! Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
    ! Only for dimensions where periodicity is present.
    IF (nqc1 == 1) THEN
      nr1x = 0
    ELSE
      nr1x = INT(SQRT(geg) / SQRT(bg(1, 1)**2 + bg(2, 1)**2 + bg(3, 1)**2)) + 1
    ENDIF
    IF (nqc2 == 1) THEN
      nr2x = 0
    ELSE
      nr2x = INT(SQRT(geg) / SQRT(bg(1, 2)**2 + bg(2, 2)**2 + bg(3, 2)**2)) + 1
    ENDIF
    IF (nqc3 == 1) THEN
      nr3x = 0
    ELSE
      nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
    ENDIF
    !
    !JLB: Find igmin_qG = min_{G}|q+G|
    qtmp=q
    CALL cryst_to_cart(1, qtmp, at, -1)
    CALL find_min_qG(qtmp)
    ! shift G-sum and center around min_{G}|q+G|, to ensure periodicity
    mmin(1) = -nr1x + igmin_qG(1)
    mmax(1) =  nr1x + igmin_qG(1)
    mmin(2) = -nr2x + igmin_qG(2)
    mmax(2) =  nr2x + igmin_qG(2)
    mmin(3) = -nr3x + igmin_qG(3)
    mmax(3) =  nr3x + igmin_qG(3)
    !
    epmatl(:) = czero
    !DO m1 = -nr1x, nr1x
    !  DO m2 = -nr2x, nr2x
    !    DO m3 = -nr3x, nr3x
    !
    DO m1 = mmin(1), mmax(1)
      DO m2 = mmin(2), mmax(2)
        DO m3 = mmin(3), mmax(3)
          !
          gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)) * (twopi / alat)
          gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)) * (twopi / alat)
          gg(3) = (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)) * (twopi / alat)
          !
          IF (system_2d) THEN
            qeq = gg(1)**2 + gg(2)**2 + gg(3)**2
            grg = zero
            IF (gg(1)**2 + gg(2)**2 > eps8) THEN
              grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + &
                    gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
              grg = grg / (gg(1)**2 + gg(2)**2)
            ENDIF
          ELSE
            qeq = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
                   gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
                   gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
          ENDIF
          IF (qeq > 0.0d0 .AND. qeq / (metric * alph * 4.0d0) < gmax) THEN
            !
            IF (system_2d) THEN
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / (SQRT(qeq) * (1.0 + grg * SQRT(qeq)))
            ELSE
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / qeq  !<-- this is correct
              !facqd = fac * EXP(-qeq * DSQRT(metric) / (metric * alph * 4.0d0)) / qeq ! <-- this is to keep as previous
              !facqd = fac / qeq ! no Ewald
            ENDIF
            !
            DO na = 1, nat
              arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na) + gg(3) * tau(3, na))
              facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
              ! Cartesian direction
              DO i = 1, 3
                zaq = zero
                DO ipol = 1, 3
                  zaq = zaq + gg(ipol) * zeu(ipol, i, na)
                ENDDO
                !
                Qqq = zero
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    Qqq = Qqq + 0.5 * gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
                  ENDDO
                ENDDO
                !
                epmat = epmat + facq * (zaq - ci * Qqq) * uq(3 * (na - 1) + i, :) * bmat
                epmatl = epmatl + facq * (zaq - ci * Qqq) * uq(3 * (na - 1) + i, :) * bmat
                !
              ENDDO !ipol
            ENDDO !nat
          ENDIF
          !
        ENDDO ! m3
      ENDDO ! m2
    ENDDO ! m1
    !
    ! In case we want only the short-range we do
    ! g_s = DSQRT(g*g - g_l*g_l)
    !
    ! Important notice: It is possible that (g*g - g_l*g_l) < 0, in which
    ! case the sqrt will give an pure imaginary number. If it is positive we
    ! will get a pure real number.
    ! In any case, when g_s will be squared both will become real numbers.
    IF (shortrange) THEN
      !epmat = ZSQRT(epmat*CONJG(epmat) - epmatl*CONJG(epmatl))
      epmat = SQRT(epmat * CONJG(epmat) - epmatl * CONJG(epmatl))
    ENDIF
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE rgd_blk_epw
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE rgd_blk_epw_fine(nqc1, nqc2, nqc3, q, uq, epmat, nmodes, epsil, zeu, bmat, signe)
    !-------------------------------------------------------------------------------
    !!
    !! Compute the long range term for the e-ph vertex
    !! to be added or subtracted from the vertex
    !!
    !! The long-range part can be computed using Eq. (4) of PRL 115, 176401 (2015).
    !! The sum over G is converged using the Ewald summation technique (see for example
    !! F.2, p.500 in Martin Electronic structure book) where the Ewald factor is ((q+G)**2)/alph/4.0d0.
    !!
    !! Technical note: From the solution of the Poisson equation, there is an additional factor
    !! e^{-i(q+G)\tau_\kappa} with respect to Eq. (4) of PRL 115, 176401 (2015).
    !! The full equation can be found in Eq. (S4) of the supplemental materials of PRL 115, 176401 (2015).
    !! In practical calculations the G-vector sum is restricted to small |q + G| via the cutoff
    !! function e^{-({\bf q}+{\bf G})^2/4\alpha}. See footnote 7 of RMP 89, 015003 (2017).
    !!
    !! The final implemented formula is:
    !!
    !! $$ g_{mn\nu}^{\mathcal L}({\bf k},{\bf q) = i\frac{4\pi e^2}{\Omega} \sum_{\kappa}
    !!   \left(\frac{\hbar}{2 {M_\kappa \omega_{{\bf q}\nu}}}\right)^{\!\!\frac{1}{2}}
    !!   \sum_{{\bf G}\ne -{\bf q}} e^{-({\bf q}+{\bf G})^2/4\alpha}
    !! \frac{ ({\bf q}+{\bf G})\cdot{\bf Z}^*_\kappa \cdot {\bf e}_{\kappa\nu}({\bf q}) }
    !!  {({\bf q}+{\bf G})\cdot\bm\epsilon^\infty\!\cdot({\bf q}+{\bf G})}\,
    !!   \left[ U_{{\bf k}+{\bf q}}\:U_{{\bf k}}^{\dagger} \right]_{mn} $$
    !!
    !! 10/2016 - SP: Optimization
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg, omega, alat, at
    USE ions_base,     ONLY : tau, nat
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12, zero, eps8
    USE epwcom,        ONLY : shortrange, nbndsub, lpolar, system_2d
    USE elph2,         ONLY : area, Qmat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nmodes
    !! Max number of modes
    REAL (KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL (KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL (KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! effective charges tensor
    REAL (KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract long range term
    COMPLEX (KIND = DP), INTENT(in) :: uq(nmodes, nmodes)
    !! phonon eigenvec associated with q
    COMPLEX (KIND = DP), INTENT(inout) :: epmat(nbndsub, nbndsub, nmodes)
    !! e-ph matrix elements
    COMPLEX (KIND = DP), INTENT(in) :: bmat(nbndsub, nbndsub)
    !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
    !
    ! Local variables
    INTEGER :: na
    !! Atom index 1
    INTEGER :: i
    !! Cartesian index direction
    INTEGER :: ipol, jpol
    !! Polarison
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: imode
    !! Mode index
    INTEGER :: nr1x, nr2x, nr3x
    !! Minimum supercell size to include all vector such that G^2 < geg
    INTEGER :: mmin(3), mmax(3)
    !! Shifted G-loop to be centered around min_{G}|q+G|
    INTEGER :: nGtest
    !! Number of G-vectors within cutoff (for testing purposes only)
    REAL(KIND = DP):: metric
    !! (2*pi/a)^2
    REAL(KIND = DP) :: qeq
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: arg
    !! Argument of the cos,sin for the Euler formula
    REAL(KIND = DP) :: zaq
    !! Z^* \cdot (q+g)
    REAL(KIND = DP) :: gg(3)
    !! G-vector
    !REAL(KIND = DP) :: gmax
    !!!  Max G-vector
    !REAL(KIND = DP) :: alph
    !!! Ewald factor (arbitrary, here chosen to be 1)
    REAL(KIND = DP) :: geg
    !! <G| epsil | G>
    REAL(KIND = DP) :: reff(2, 2)
    !! Effective screening length for 2D materials
    REAL(KIND = DP) :: grg
    !! G-vector * reff * G-vector
    REAL(KIND = DP) :: Qqq
    !! In the case of Si, its a single value
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    REAL(KIND = DP) :: qtmp(3)
    !! Temporary q vector to find min_{G}|q+G|
    COMPLEX(KIND = DP) :: fac
    !! General prefactor
    COMPLEX(KIND = DP) :: facqd
    !! Ewald filtering
    COMPLEX(KIND = DP) :: facq
    !! Atomic position exponential
    COMPLEX(KIND = DP) :: epmatl(nbndsub, nbndsub, nmodes)
    !! Long-range part of the matrix element
    !
    ! Impose zero Born effective charge in case of non-polar materials with quadrupoles
    IF (.NOT. lpolar) THEN
      zeu(:, :, :) = zero
    ENDIF
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore ('rgd_blk_epw_fine', 'Wrong value for signe ', 1)
    !
    IF (system_2d) THEN
      ! Vacuum size in Bohr unit
      c = alat / bg(3, 3)
      ! (e^2 * 2\pi * ci) / Area
      fac = (signe * e2 * twopi * ci) / area
      ! Effective screening length
      ! reff = (epsil - 1) * c/2
      reff(:, :) = zero
      reff(:, :) = epsil(1:2, 1:2) * 0.5d0 * c ! eps * c/2 in 2pi/a units
      reff(1, 1) = reff(1, 1) - 0.5d0 * c ! (-1) * c/2 in 2pi/a units
      reff(2, 2) = reff(2, 2) - 0.5d0 * c ! (-1) * c/2 in 2pi/a units
    ELSE
      ! (e^2 * 4\pi * i) / Volume
      fac = (signe * e2 * fpi * ci) / omega
    ENDIF
    !
    !gmax = 14.d0
    !alph = 1.0d0
    metric = (twopi / alat)**2
    geg = gmax * alph * 4.0d0
    !
    ! Estimate of nr1x, nr2x, nr3x generating all vectors up to G^2 < geg
    ! Only for dimensions where periodicity is present.
    IF (nqc1 == 1) THEN
      nr1x = 0
    ELSE
      nr1x = INT(SQRT(geg) / SQRT(bg(1, 1)**2 + bg(2, 1)**2 + bg(3, 1)**2)) + 1
    ENDIF
    IF (nqc2 == 1) THEN
      nr2x = 0
    ELSE
      nr2x = INT(SQRT(geg) / SQRT(bg(1, 2)**2 + bg(2, 2)**2 + bg(3, 2)**2)) + 1
    ENDIF
    IF (nqc3 == 1) THEN
      nr3x = 0
    ELSE
      nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
    ENDIF
    !
    !JLB: Find igmin_qG = min_{G}|q+G|
    qtmp=q
    CALL cryst_to_cart(1, qtmp, at, -1)
    CALL find_min_qG(qtmp)
    ! shift G-sum and center around min_{G}|q+G|, to ensure periodicity
    mmin(1) = -nr1x + igmin_qG(1)
    mmax(1) =  nr1x + igmin_qG(1)
    mmin(2) = -nr2x + igmin_qG(2)
    mmax(2) =  nr2x + igmin_qG(2)
    mmin(3) = -nr3x + igmin_qG(3)
    mmax(3) =  nr3x + igmin_qG(3)
    !
    !nGtest = 0
    epmatl(:, :, :) = czero
    !DO m1 = -nr1x, nr1x
    !  DO m2 = -nr2x, nr2x
    !    DO m3 = -nr3x, nr3x
    DO m1 = mmin(1), mmax(1)
      DO m2 = mmin(2), mmax(2)
        DO m3 = mmin(3), mmax(3)
          !
          gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)) * (twopi / alat)
          gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)) * (twopi / alat)
          gg(3) = (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)) * (twopi / alat)
          !
          IF (system_2d) THEN
            qeq = gg(1)**2 + gg(2)**2 + gg(3)**2
            grg = zero
            IF (gg(1)**2 + gg(2)**2 > eps8) THEN
              grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + &
                    gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
              grg = grg / (gg(1)**2 + gg(2)**2)
            ENDIF
          ELSE
            !
            qeq = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
                   gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
                   gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
          ENDIF
          !
          IF (qeq > 0.0d0 .AND. qeq / (metric * alph * 4.0d0) < gmax) THEN
            !
            !nGtest = nGtest + 1
            !
            IF (system_2d) THEN
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / (SQRT(qeq) * (1.0 + grg * SQRT(qeq)))
            ELSE
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / qeq  !<-- this is correct
              !facqd = fac * EXP(-qeq * DSQRT(metric) / (metric * alph * 4.0d0)) / qeq ! <-- this is to keep as previous
              !facqd = fac / qeq ! no Ewald
            ENDIF
            !
            DO na = 1, nat
              arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na) + gg(3) * tau(3, na))
              facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
              ! Cartesian index direction
              DO i = 1, 3
                zaq = zero
                DO ipol = 1, 3
                  zaq = zaq + gg(ipol) * zeu(ipol, i, na)
                ENDDO
                !
                Qqq = zero
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    Qqq = Qqq + 0.5 * gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
                  ENDDO
                ENDDO
                !
                DO imode = 1, nmodes
                  CALL ZAXPY(nbndsub**2, facq * (zaq - ci * Qqq) * uq(3 * (na - 1) + i, imode),&
                              bmat(:, :), 1, epmat(:, :, imode), 1)
                  CALL ZAXPY(nbndsub**2, facq * (zaq - ci * Qqq) * uq(3 * (na - 1) + i, imode), &
                              bmat(:, :), 1, epmatl(:, :, imode), 1)
                ENDDO
                !
              ENDDO !ipol
            ENDDO !nat
          ENDIF
          !
        ENDDO ! m3
      ENDDO ! m2
    ENDDO ! m1
    !
    !write(*,*) "Number of G-vectors within cutoff:", nGtest
    !
    ! In case we want only the short-range we do
    ! g_s = DSQRT(g*g - g_l*g_l)
    !
    ! Important notice: It is possible that (g*g - g_l*g_l) < 0, in which
    ! case the sqrt will give an pure imaginary number. If it is positive we
    ! will get a pure real number.
    ! In any case, when g_s will be squared both will become real numbers.
    IF (shortrange) THEN
      !epmat = ZSQRT(epmat*CONJG(epmat) - epmatl*CONJG(epmatl))
      epmat = SQRT(epmat * CONJG(epmat) - epmatl * CONJG(epmatl))
    ENDIF
    !
    !-----------------------------------------------------------------------------
    END SUBROUTINE rgd_blk_epw_fine
    !-----------------------------------------------------------------------------
    !
    !!!!!
    !-------------------------------------------------------------------------------
    SUBROUTINE rgd_imp_epw_fine(nqc1, nqc2, nqc3, q, epmat, epsil, bmat, signe)
    !-------------------------------------------------------------------------------
    !!
    !! JL 01/2022
    !! Compute the unscreened carrier-ionized impurity matrix element
    !! The element is a Coulomb potential
    !! (4*pi*Z*e^2)/(eps_inf*q^2)*U_{n,k}*U_{m,k+q}
    !! Calculated on the fine grid only
    !! Valid for the long range interaction, as U matricies for q=0 only
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg, omega, alat
    USE ions_base,     ONLY : tau, nat
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12, cc2cb
    USE epwcom,        ONLY : shortrange, nbndsub, ii_n, ii_charge, ii_eps0
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    !INTEGER, INTENT(in) :: nmodes
    !! Max number of modes
    REAL (KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL (KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL (KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract long range term
    COMPLEX (KIND = DP), INTENT(inout) :: epmat(nbndsub, nbndsub)
    !! e-ph matrix elements
    COMPLEX (KIND = DP), INTENT(in) :: bmat(nbndsub, nbndsub)
    !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
    !
    ! Local variables
    INTEGER :: na
    !! Atom index 1
    INTEGER :: ipol
    !! Polarison
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: imode
    !! Mode index
    REAL(KIND = DP) :: qeq
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: arg
    !!
    REAL(KIND = DP) :: zaq
    !!
    REAL(KIND = DP) :: g1, g2, g3
    !!
    REAL(KIND = DP) :: gmax
    !!
    REAL(KIND = DP) :: alph
    !!
    REAL(KIND = DP) :: geg
    !!
    REAL(KIND = DP) :: impurity_density
    !! Impurity density per unit cell (1/N factor, where N is the number of 
    !! cells to contain one imp)
    COMPLEX(KIND = DP) :: fac
    !!
    COMPLEX(KIND = DP) :: facqd
    !!
    COMPLEX(KIND = DP) :: facq
    !!
    REAL(KIND = DP) :: m1f
    !! fixed value of m1 to minimize |q - G|
    REAL(KIND = DP) :: m2f
    !! fixed value of m2 to minimize |q - G|
    REAL(KIND = DP) :: m3f
    !! fixed value of m3 to minimize |q - G|
    REAL(KIND = DP) :: qmG
    !! value of |q - G| to be minimized with respect to choices mNf
    REAL(KIND = DP) :: val
    !! comparison value
    REAL(KIND = DP) :: eps_loc(3,3)
    !! local epsilon
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore ('rgd_blk_epw_fine', 'Wrong value for signe ', 1)
    !
    eps_loc(:, :) = 0.0d0
    !
    IF (ii_eps0 == 0.0d0) THEN
      eps_loc(1,1) = epsil(1,1)
      eps_loc(2,2) = epsil(2,2)
      eps_loc(3,3) = epsil(3,3)
      eps_loc(1,2) = epsil(1,2)
      eps_loc(2,1) = epsil(2,1)
      eps_loc(1,3) = epsil(1,3)
      eps_loc(3,1) = epsil(3,1)
      eps_loc(2,3) = epsil(2,3)
      eps_loc(3,2) = epsil(3,2)
    ELSEIF (ii_eps0 > 0.0d0) THEN
      eps_loc(1,1) = ii_eps0
      eps_loc(2,2) = ii_eps0
      eps_loc(3,3) = ii_eps0
    ENDIF
    !
    gmax = 14.d0
    alph = 1.0d0
    geg = gmax * alph * 4.0d0
    fac = signe * e2 * fpi * ii_charge / omega * ci
    !
    ! JL : look for G0 that minimizes | q + G | with G = 0
    ! INITIAL VALUES
    val = 1.0d24
    m1f = 0
    m2f = 0
    m3f = 0
    !
    DO m1 = -5, 5
      DO m2 = -5, 5
        DO m3 = -5, 5
          g1 = q(1) + (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3))
          g2 = q(2) + (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3))
          g3 = q(3) + (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3))
          qmG = SQRT(g1**2+g2**2+g3**2)
          IF (qmG < val) THEN
            val = qmG
            m1f = m1
            m2f = m2
            m3f = m3
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    ! JL: Now, sum Gi on first shell displaced by G0, that is q + Gi + G0 
    !
    DO m1 = -1, 1 !-nqc1, nqc1
      DO m2 = -1, 1 !-nqc2, nqc2
        DO m3 = -1, 1 !-nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1) + m1f * bg(1, 1) + &
                  m2f * bg(1, 2) + m3f * bg(1, 3)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2) + m1f * bg(2, 1) + &
                  m2f * bg(2, 2) + m3f * bg(2, 3)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3) + m1f * bg(3, 1) + &
                  m2f * bg(3, 2) + m3f * bg(3, 3)
          !
          qeq = (g1 * (eps_loc(1, 1) * g1 + eps_loc(1, 2) * g2 + eps_loc(1, 3) * g3) + &
                 g2 * (eps_loc(2, 1) * g1 + eps_loc(2, 2) * g2 + eps_loc(2, 3) * g3) + &
                 g3 * (eps_loc(3, 1) * g1 + eps_loc(3, 2) * g2 + eps_loc(3, 3) * g3)) !*twopi/alat
          !
          IF (qeq > 0.0d0) THEN !  .AND. qeq / (alph * 4.0d0) < gmax) THEN
            !
            qeq = qeq * (twopi / alat)**2.0d0
            facqd = fac / qeq !fac * EXP(-qeq / (alph * 4.0d0)) / qeq !/(two*wq)
            !
            facq = facqd * 1.0d0 ! For now add impurty at position 0, 0, 0, JL
            CALL ZAXPY(nbndsub**2, facq , bmat(:, :), 1, epmat(:, :), 1)
          ENDIF
          !
        ENDDO ! m3
      ENDDO ! m2
    ENDDO ! m1
    !
    !-----------------------------------------------------------------------------
    END SUBROUTINE rgd_imp_epw_fine
    !-----------------------------------------------------------------------------
    !
    !!!!!
    !-----------------------------------------------------------------------------
    SUBROUTINE rpa_epsilon(q, w, nmodes, epsil, eps_rpa)
    !-----------------------------------------------------------------------------
    !!
    !!  Compute the Lindhard dielectric function for the homogeneous electron gas
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, omega, alat
    USE constants_epw, ONLY : twopi, ha2ev, cone, ci, eps5, eps10
    USE constants,     ONLY : pi
    USE epwcom,        ONLY : meff, fermi_diff, nel, smear_rpa, system_2d
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Number of phonon modes
    REAL (KIND = DP), INTENT(inout) :: q(3)
    !! q vector (in crystal coordinates
    REAL (KIND = DP), INTENT(inout) :: w(nmodes)
    !! phonon frequencies associated with q
    REAL (KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    COMPLEX(KIND = DP), INTENT(out) :: eps_rpa(nmodes)
    !! electronic screening
    !
    ! Local variable
    LOGICAL, SAVE :: first_call = .TRUE.
    !! Logical for first_call the routine
    INTEGER :: im
    !! Mode counter
    REAL(KIND = DP) :: n
    !! Electron density in atomic units
    REAL(KIND = DP) :: rs
    !! Prefactor for the dielectric screening
    REAL(KIND = DP) :: EF
    !! Fermi-level in eV
    REAL(KIND = DP) :: kF
    !! Fermi wavevector
    REAL(KIND = DP) :: pref
    !! Prefactor for the dielectric function
    REAL(KIND = DP) :: eta
    !! Broadening for the dielectric function
    REAL(KIND = DP) :: q2
    !! q-point square
    REAL(KIND = DP) :: qm
    !! Internal units for Hedin's formula
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    COMPLEX(KIND = DP) :: u
    !! Complex frequency argument
    !
    IF (system_2d) THEN
      n = (at(3, 3) * alat) * nel / omega
    ELSE
      n = nel / omega
    ENDIF
    EF = fermi_diff / ha2ev
    kF = (3.d0 * pi**2 * n)**(1.d0 / 3.d0)
    eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    rs  = (3.d0 / ( 4.d0 * pi * n ) )**(1.d0 / 3.d0) * meff / eps_ave
    w = w * 0.5d0 / EF / 4.d0 !Ha&internal units for Hedin's formula
    pref = (4.d0 / 9.d0 / pi )**(1.d0 / 3.0) * (rs / 8.d0 / pi)
    eta = smear_rpa / ha2ev / EF / 4.d0
    !
    IF (first_call) THEN
      first_call = .FALSE.
      WRITE(stdout, '(5x,"Calculation of Lindhard screening: use with care")')
      WRITE(stdout, '(5x,"Warning: current implementation for doubly degenerate band, one valley")')
      WRITE(stdout, '(5x,a,f12.8,a,f12.8,a,f12.8)') 'Nel = ', nel, ', n = ', n, ' au^-3, meff = ', meff
      WRITE(stdout, '(5x,a,f12.8,a,f12.8,a,f12.8)') 'EF = ', EF * ha2ev, ' eV, kF = ', kF, ' au^-1, rs = ', rs
      IF (eps_ave < eps5) WRITE(stdout, '(5x,"Warning: dielectric constant not found; set to 1")')
    ENDIF
    IF (eps_ave < eps5) eps_ave = 1.d0
    !
    CALL cryst_to_cart(1, q, bg, 1)
    q2 = q(1)**2 + q(2)**2 + q(3)**2
    qm = DSQRT(q2) * (twopi / alat) / kF / 2.d0 ! internal units for Hedin's formula
    !
    IF (ABS(qm) > eps10) THEN
      DO im = 1, nmodes
        u = w(im) + SIGN(eta, w(im)) * ci
        eps_rpa(im) = 1.d0 + pref * (H_eps(qm + u / qm) + H_eps(qm - u / qm)) / qm**3
      ENDDO
    ELSE
      eps_rpa = cone
    ENDIF
    !
    w = w / (0.5d0 / EF / 4.d0)
    CALL cryst_to_cart(1, q, at, -1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE rpa_epsilon
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE tf_epsilon(q, nmodes, epsil, eps_tf)
    !--------------------------------------------------------------------------
    !!
    !!  Compute the Thomas-Fermi dielectric screening
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, omega, alat
    USE constants_epw, ONLY : twopi, ha2ev, cone, eps5, eps10
    USE constants,     ONLY : pi
    USE epwcom,        ONLY : fermi_diff, nel
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! Number of phonon modes
    REAL(KIND = DP), INTENT(inout) :: q(3)
    !! q vector (in crystal coordinates)
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    COMPLEX(KIND = DP), INTENT(out) :: eps_tf(nmodes)
    !! electronic screening
    !
    ! Local variable
    LOGICAL, SAVE :: first_call = .TRUE.
    !! Logical for first_call the routine
    REAL(KIND = DP) :: n
    !! Electron density in atomic units
    REAL(KIND = DP) :: EF
    !! Fermi-level in eV
    REAL(KIND = DP) :: q2
    !! q-point square
    REAL(KIND = DP) :: qtf
    !! Thomas-Fermi wavector
    REAL(KIND = DP) :: qtfc
    !! Thomas-Fermi wavector in unit of 2pi/a
    REAL(KIND = DP) :: qm
    !! Modulus of q
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    !
    n = nel / omega
    EF = fermi_diff / ha2ev
    eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    qtf = (6.d0 * pi * n / EF / eps_ave )**(1.d0 / 2.d0)
    qtfc = qtf / (twopi / alat)
    !
    IF (first_call) THEN
      first_call = .FALSE.
      WRITE(stdout, '(5x,"Calculation of Thomas-Fermi screening: use with care")')
      WRITE(stdout, '(5x,"Warning: current implementation for doubly degenerate band, one valley")')
      WRITE(stdout, '(5x,a,f12.8,a,f12.8,a,f12.8)') 'Nel = ', nel, ', n = ', n, ' au^-3, EF (eV) = ', EF*ha2ev
      WRITE(stdout, '(5x,a,f12.8,a,f12.8)') 'q_tf (au-1) = ', qtf, ', q_tf (tpiba) = ', qtfc
      IF (eps_ave < eps5) WRITE(stdout, '(5x,"Warning: dielectric constant not found; set to 1")')
    ENDIF
    IF (eps_ave < eps5) eps_ave = 1.d0
    !
    CALL cryst_to_cart(1, q, bg, 1)
    q2 = q(1)**2 + q(2)**2 + q(3)**2
    qm = DSQRT(q2) ! in tpiba
    IF (ABS(qm) > eps10) THEN
      eps_tf = 1.d0 + qtfc**2 / q2
    ELSE
      eps_tf = cone
    ENDIF
    !
    CALL cryst_to_cart(1, q, at, -1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE tf_epsilon
    !--------------------------------------------------------------------------
    !
    !!!!!
    !--------------------------------------------------------------------------
    SUBROUTINE calc_qtf2_therm(itemp, etemp, ef0, efcb, ctype, epsil)
    !--------------------------------------------------------------------------
    !!
    !! JL 01/2022
    !! Compute thermal Thomas fermi wave vector for doped Semiconductores
    !! ctype = +/1 1 only!! assume_metal = .FALSE. only (for now...) 
    !! q_{tf}^{2}(T) = \frac{4\pi}{\omega*\varepsilon_\infty}\sum_{n,\mathbf{k}}
    !! w_{n\mathbf k}\frac{\partial f_{n,\mathbf k}}{\partial \epsilon_{n\mathbf
    !! k}}
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : omega, at, alat
    !
    USE epwcom,        ONLY : nstemp, nbndsub, ii_eps0
    USE elph2,         ONLY : ibndmin, etf, nkf, wkf, &
                              nbndfst, qtf2_therm
    USE constants,     ONLY : pi
    USE io_global,     ONLY : stdout
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! temp index
    INTEGER, INTENT(in) :: ctype
    !! -1 for holes, 1 for electrons, 0 for both (forbidden)
    REAL(KIND = DP), INTENT(in)  :: etemp
    !! energy of temperatre at index itemp
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! fermi energy of p-type or metal material
    REAL(KIND = DP), INTENT(in) :: efcb(nstemp)
    !! fermi_energy of n-type materials
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric tensor
    !
    ! local variables
    INTEGER :: dummy1
    !! dummy
    INTEGER :: ik
    !! counter on k ind
    INTEGER :: ikk
    !! counter over k on fine grid
    INTEGER :: ibnd
    !! counter over bands
    REAL(KIND=DP) :: ekk
    !! Energy on fine grid at ekk in fermi window
    REAL(KIND=DP) :: dfnk
    !! derivative of the fermi-dirac distribution
    REAL(KIND=DP) :: inv_etemp
    !! beta = 1/(kB*T), inverse thermal temperature in Ryd
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    !
    IF (ii_eps0 == 0.0) THEN
      eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    ELSEIF (ii_eps0 > 0.0) THEN
      eps_ave = ii_eps0
    ENDIF
    !eps_ave = 1.0d0
    !
    IF (ctype==0) THEN
      WRITE(stdout, '(5x,"ctype=0 not supported, keeping epf_tf_therm(:) == 1.0d0")')
      RETURN
    ENDIF
    !
    !IF (my_pool_id == 0) THEN
    !  IF (ctype == -1) THEN
    !    WRITE(stdout, '(5x,"Using Fermi energy of: ")')
    !    WRITE(stdout,*) ef0(itemp)
    !  ELSEIF (ctype == 1) THEN
    !    WRITE(stdout, '(5x,"Using Fermi energy of: ")')
    !    WRITE(stdout,*) efcb(itemp)
    !  ENDIF
    !ENDIF
    !
    inv_etemp = 1.0d0 / etemp
    !
    ! calc qtf^2(T) for holes in case of ctype=-1
    IF (ctype==-1) THEN
      DO ik = 1, nkf
        ikk = 2*ik-1
        DO ibnd = 1, nbndsub
          ekk  = etf(ibnd, ikk)
          dfnk = w0gauss((ekk-ef0(itemp))*inv_etemp,-99)*inv_etemp
          qtf2_therm(itemp) = qtf2_therm(itemp) + (1.0d0/(eps_ave*omega))*8.0d0*pi*( wkf(ikk)*dfnk)
        ENDDO
      ENDDO
    ENDIF
    !
    ! calc qtf^2(T) for holes in case of ctype=-1
    IF (ctype==1) THEN
      DO ik = 1, nkf
        ikk = 2*ik-1
        DO ibnd = 1, nbndsub
          ekk  = etf(ibnd, ikk)
          dfnk = w0gauss((ekk-efcb(itemp))*inv_etemp,-99)*inv_etemp
          qtf2_therm(itemp) = qtf2_therm(itemp) + (1.0d0/(omega*eps_ave))*8.0d0*pi*( wkf(ikk)*dfnk  )
        ENDDO
      ENDDO
    ENDIF
    !
    CALL mp_sum(qtf2_therm(itemp), inter_pool_comm)
    !
    !IF (my_pool_id == 0) THEN
    !  !IF (iqq == 0 ) THEN
    !  WRITE(stdout, '(5x,"ctype, itemp, q_tf_therm^2: ")')
    !  WRITE(stdout,*) ctype, itemp, qtf2_therm(itemp)
    !  !ENDIF
    !ENDIF
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_qtf2_therm
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE calc_epstf_therm(q, nstemp, epsil)
    !--------------------------------------------------------------------------
    !!
    !! JL: 01/2022
    !! Calculate the thermal thomas fermi dielectric function (1 + qtf^2(T)/q^2)
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, omega, alat
    USE constants_epw, ONLY : twopi, ha2ev, cone, eps5, eps10, one
    USE elph2,         ONLY : qtf2_therm, epstf_therm
    !
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nstemp
    !! number of temperature points
    REAL(KIND = DP), INTENT(inout) :: q(3)
    !! q vector (in crystal coordinates)
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    !
    ! Local variable
    LOGICAL, SAVE :: first_call = .TRUE.
    !! Logical for first_call the routine
    INTEGER :: itemp
    !! counter of temperature index
    INTEGER :: m1, m2, m3
    !! counter for G vecs
    INTEGER :: m1f, m2f, m3f
    !! index of G that minimizes q+G
    REAL(KIND = DP) :: eps_ave
    !! Average dielectric function (semiconductor/insulator)
    REAL(KIND = DP) :: q2
    !! squared phonon wavevector
    REAL(KIND = DP) :: q2inv
    !! inverse of q2, which is summed over G
    REAL(KIND = DP) :: qm
    !! sqrt of q2
    REAL(KIND = DP) :: qtfc2
    !! square Thomas fermi wv in units of (two*pi/alat)
    !!!!!!qtfc = qtf / (twopi / alat)
    REAL(KIND = DP) :: val
    !! test value for finding G that minimizes q+G
    REAL(KIND = DP) :: g1, g2, g3, qeq, qmG 
    !!
    REAL(KIND = DP) :: eps_loc(3,3)
    !
    eps_ave = (epsil(1, 1) + epsil(2, 2) + epsil(3, 3)) / 3.d0
    eps_loc(1,1) = 1.0d0
    eps_loc(2,2) = 1.0d0
    eps_loc(3,3) = 1.0d0
    eps_loc(1,2) = 0.0d0
    eps_loc(2,1) = 0.0d0
    eps_loc(1,3) = 0.0d0
    eps_loc(3,1) = 0.0d0
    eps_loc(2,3) = 0.0d0
    eps_loc(3,2) = 0.0d0
    !
    IF (first_call) THEN
      first_call = .FALSE.
      WRITE(stdout, '(5x,"Calculation of thermal Thomas-Fermi screening: use with care")')
      DO itemp = 1, nstemp
        WRITE(stdout, '(5x,a,i5)') 'itemp=', itemp
        WRITE(stdout, '(5x,a,f22.16)') 'q_tf (au^-1) = ', SQRT(qtf2_therm(itemp))
      ENDDO
      IF (eps_ave < eps5) WRITE(stdout, '(5x,"Warning: dielectric constant not found; set to 1")')
    ENDIF
    !
    CALL cryst_to_cart(1, q, bg, 1)
    !
    ! JL: Look for G that minimizes q+G --> enforce periodicity
    val = 1.0d24
    m1f = 0
    m2f = 0
    m3f = 0
    !
    DO m1 = -5, 5
      DO m2 = -5, 5
        DO m3 = -5, 5
          g1 = q(1) + (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3))
          g2 = q(2) + (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3))
          g3 = q(3) + (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3))
          qmG = SQRT(g1**2+g2**2+g3**2)
          IF (qmG < val) THEN
            val = qmG
            m1f = m1
            m2f = m2
            m3f = m3
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    q2inv = 0.0d0
    !
    DO m1 = -1, 1 !-nqc1, nqc1
      DO m2 = -1, 1 !-nqc2, nqc2
        DO m3 = -1, 1 !-nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1) + m1f * bg(1, 1) + &
                  m2f * bg(1, 2) + m3f * bg(1, 3)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2) + m1f * bg(2, 1) + &
                  m2f * bg(2, 2) + m3f * bg(2, 3)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3) + m1f * bg(3, 1) + &
                  m2f * bg(3, 2) + m3f * bg(3, 3)
          !
          !qeq = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
          !       g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
          !       g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3)) !*twopi/alat
          qeq = (g1 * (eps_loc(1, 1) * g1 + eps_loc(1, 2) * g2 + eps_loc(1, 3) * g3) + &
                 g2 * (eps_loc(2, 1) * g1 + eps_loc(2, 2) * g2 + eps_loc(2, 3) * g3) + &
                 g3 * (eps_loc(3, 1) * g1 + eps_loc(3, 2) * g2 + eps_loc(3, 3) * g3)) !*twopi/alat
          IF (qeq > 0.0d0) THEN 
            q2inv = q2inv + (1.0d0 / qeq)
          ELSE
            q2inv = q2inv + 0.0d0 
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    ! 
    !! JL now calculate contributions to epstf_therm(itmep)
    !
    IF (q2inv > 0.0d0) THEN
      q2 = 1.0d0 / q2inv
    ELSE
      q2 = 0.0
    ENDIF 
    !q2 = q(1)**2 + q(2)**2 + q(3)**2
    qm = DSQRT(q2) ! in tpiba
    DO itemp = 1, nstemp
      IF (ABS(qm) > eps10) THEN
        qtfc2 = qtf2_therm(itemp) / ((twopi / alat)**2.0d0)
        epstf_therm(itemp) = 1.d0 + (qtfc2 / q2)
      ELSE
        epstf_therm(itemp) = one
      ENDIF
    ENDDO
    !
    CALL cryst_to_cart(1, q, at, -1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE calc_epstf_therm
    !--------------------------------------------------------------------------
    !
    !!!!!
    !-----------------------------------------------------------------------
    SUBROUTINE compute_umn_f(nbnd, cufkk, cufkq, bmatf)
    !-----------------------------------------------------------------------
    !!
    !! Calculates $$ U_{k+q} U_k^\dagger = <\Psi_{mk+q}|e^{i{q+G}r}|\Psi_{nk}> $$
    !! in the approximation q+G->0 on the fine grids.
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, cone
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands (possibly in the optimal subspace)
    COMPLEX(KIND = DP), INTENT(in) :: cufkk(nbnd, nbnd)
    !! rotation matrix U(k)^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(in) :: cufkq(nbnd, nbnd)
    !! rotation matrix U(k+q)^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(out) :: bmatf(nbnd, nbnd)
    ! overlap wfcs in Bloch representation, fine grid
    !
    ! Every pool works with its own subset of k points on the fine grid
    bmatf = czero
    !
    !  U(k'+q')^\dagger * U(k')
    !
    CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, cufkq, nbnd, cufkk, nbnd, czero, bmatf, nbnd)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_umn_f
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_umn_c(nbnd, nbndsub, nks, cuk, cukq, bmat)
    !-----------------------------------------------------------------------
    !!
    !! Calculates $$ U_{k+q} U_k^\dagger = <\Psi_{mk+q}|e^{i(q+G)r}|\Psi_{nk}> $$
    !! in the approximation q+G->0 on the coarse grids.
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : czero, cone
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands
    INTEGER, INTENT(in) :: nks
    !! number of kpoint blocks, per pool
    INTEGER, INTENT(in) :: nbndsub
    !! Number of band on the subspace of Wannier
    COMPLEX(KIND = DP), INTENT(in) :: cuk(nbnd, nbndsub, nks)
    !! rotation matrix U(k), coarse mesh
    COMPLEX(KIND = DP), INTENT(in) :: cukq(nbnd, nbndsub, nks)
    !! rotation matrix U(k+q), coarse mesh
    COMPLEX(KIND = DP), INTENT(out) :: bmat(nbnd, nbnd, nks)
    !! overlap wfcs in Bloch representation, fine grid
    !
    ! Local variables
    INTEGER :: ik
    !! k-point index
    !
    ! Every pool works with its own subset of k points on the fine grid
    bmat = czero
    !
    !  U(k+q) * U(k)^\dagger
    !
    DO ik = 1, nks
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbndsub, cone, cukq(:, :, ik), &
                  nbnd, cuk(:, :, ik), nbnd, czero, bmat(:, :, ik), nbnd)
    ENDDO
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_umn_c
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rgd_blk_der(nqc1, nqc2, nqc3, nat, dyn_der, q, tau, epsil, zeu, signe)
    !-----------------------------------------------------------------------
    !!
    !! Compute the rigid-ion (long-range) derivative term for q
    !! Note 1: the derivative is only made on the dipole (no quadrupole).
    !! Note 2: would need to be updated for 2D but not important for now as only used for smearing.
    !! 2019 - Samuel Ponce & Francesco Macheda
    !! 2021 - SP: update to add nr1x and (G+q) in unit of 2pi/alat
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi, e2, ci, twopi
    USE constants,     ONLY : pi
    USE cell_base,     ONLY : bg, omega, alat
    USE constants_epw, ONLY : eps6
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nat
    !! Number of atoms
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL(KIND = DP), INTENT(in) :: zeu(3, 3, nat)
    !! effective charges tensor
    REAL(KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract rigid-ion term
    REAL(KIND = DP), INTENT(in) :: tau(3, nat)
    !! Atomic positions
    COMPLEX(KIND = DP), INTENT(inout) :: dyn_der(3, 3 * nat, 3 * nat)
    !! Dynamical matrix
    COMPLEX(KIND = DP) :: dyn_der_part(7, 3, 3 * nat, 3 * nat)
    !! Dynamical matrix, partial summations
    !
    ! Local variables
    INTEGER :: na
    !! Atom index 1
    INTEGER :: nb
    !! Atom index 2
    INTEGER :: i
    !! Cartesian direction 1
    INTEGER :: j
    !! Cartesian direction 1
    INTEGER :: isum
    !! Index to sum the different component of the derivative
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: nr1x, nr2x, nr3x
    !! Minimum supercell size to include all vector such that G^2 < geg
    REAL(KIND = DP):: metric
    !! (2*pi/a)^2
    REAL(KIND = DP):: geg
    !! <q+G| epsil | q+G>
    !REAL(KIND = DP) :: alph
    !!! Ewald parameter
    REAL(KIND = DP) :: fac
    !! General prefactor
    REAL(KIND = DP) :: gg(3)
    !! G-vectors
    REAL(KIND = DP) :: facgd
    !! fac * EXP(-geg / (alph * 4.0d0)) / geg
    REAL(KIND = DP) :: arg
    !! Argument of the exponential
    !REAL(KIND = DP) :: gmax
    !!! Maximum G
    REAL(KIND = DP) :: arg_no_g(3)
    !! Difference of atomic position
    REAL(KIND = DP) :: zag(3)
    !! Z * G
    REAL(KIND = DP) :: zbg(3)
    !! Z * G
    REAL(KIND = DP) :: zbg_der(3, 3)
    !! Z derivative
    REAL(KIND = DP) :: zag_der(3, 3)
    !! Z derivative
    COMPLEX(KIND = DP) :: facg
    !! Atomic position exponential
    !
    ! alph is the Ewald parameter, geg is an estimate of G^2
    ! such that the G-space sum is convergent for that alph
    ! very rough estimate: geg/4/alph > gmax = 14
    ! (exp (-14) = 10^-6)
    !
    !gmax = 14.d0
    !alph = 1.0d0
    metric = (twopi / alat)**2
    geg = gmax * alph * 4.0d0
    !
    IF (ABS(ABS(signe) - 1.0) > eps6) CALL errore('rgd_blk_der', ' wrong value for signe ', 1)
    !
    !gmax = 14.d0
    !alph = 1.0d0
    geg = gmax * alph * 4.0d0
    fac = signe * e2 * fpi / omega
    !
    IF (nqc1 == 1) THEN
      nr1x = 0
    ELSE
      nr1x = INT(SQRT(geg) / SQRT(bg(1, 1)**2 + bg(2, 1)**2 + bg(3, 1)**2)) + 1
    ENDIF
    IF (nqc2 == 1) THEN
      nr2x = 0
    ELSE
      nr2x = INT(SQRT(geg) / SQRT(bg(1, 2)**2 + bg(2, 2)**2 + bg(3, 2)**2)) + 1
    ENDIF
    IF (nqc3 == 1) THEN
      nr3x = 0
    ELSE
      nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
    ENDIF
    !
    DO m1 = -nr1x, nr1x
      DO m2 = -nr2x, nr2x
        DO m3 = -nr3x, nr3x
          !
          gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)) * (twopi / alat)
          gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)) * (twopi / alat)
          gg(3) = (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)) * (twopi / alat)
          !
          geg = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
                 gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
                 gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
          !
          IF (geg > 0.0d0 .AND. geg / (metric * alph * 4.0d0) < gmax) THEN
            !
            facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / geg
            !
            DO nb = 1, nat
              zbg(:) = gg(1) * zeu(1, :, nb) + gg(2) * zeu(2, :, nb) + gg(3) * zeu(3, :, nb)
              zbg_der(:, :) = zeu(:, :, nb)
              DO na = 1, nat
                zag(:) = gg(1) * zeu(1, :, na) + gg(2) * zeu(2, :, na) + gg(3) * zeu(3, :, na)
                zag_der(:, :) = zeu(:, :, na)
                arg = alat * (gg(1) * (tau(1, na) - tau(1, nb)) + gg(2) * (tau(2, na) - tau(2, nb)) + &
                              gg(3) * (tau(3, na) - tau(3, nb)))
                arg_no_g(:) = alat * (tau(:,na) - tau(:,nb))
                !
                facg = facgd * CMPLX(COS(arg), SIN(arg), KIND=DP)
                DO j = 1, 3
                  DO i = 1, 3
                    dyn_der_part(1, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) = facg * zag_der(:, i) * zbg(j)
                    dyn_der_part(2, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) = facg * zag(i) * zbg_der(:, j)
                    dyn_der_part(3, :, (na - 1) * 3 + i,( nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(:, 1) * gg(1) + epsil(:, 2) * gg(2) + epsil(:, 3) * gg(3)) / geg
                    dyn_der_part(4, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(1, :) * gg(1) + epsil(2, :) * gg(2) + epsil(3, :) * gg(3)) / geg
                    dyn_der_part(5, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) = facg * zag(i) * zbg(j) * ci * arg_no_g(:)
                    dyn_der_part(6, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(1, :) * gg(1) + epsil(2, :) * gg(2) + epsil(3, :) * gg(3)) / (4d0 * alph)
                    dyn_der_part(7, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(:, 1) * gg(1) + epsil(:, 2) * gg(2) + epsil(:, 3) * gg(3)) / (4d0 * alph)
                    DO isum = 1, 7
                      dyn_der(:, (na - 1) * 3 + i, (nb - 1) * 3 + j) = dyn_der(:, (na - 1) * 3 + i, (nb - 1) * 3 + j) &
                        + dyn_der_part(isum, :, (na - 1) * 3 + i, (nb - 1) * 3 + j)
                    ENDDO ! isum
                  ENDDO ! i
                ENDDO ! j
              ENDDO ! na
            ENDDO ! nb
          ENDIF ! geg >
        ENDDO ! m3
      ENDDO ! m2
    ENDDO ! m1
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE rgd_blk_der
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE rgd_blk_epw_fine_mem(imode, nqc1, nqc2, nqc3, q, uq, epmat, nmodes, epsil, zeu, bmat, signe)
    !-------------------------------------------------------------------------------
    !!
    !! Compute the long range term for the e-ph vertex
    !! to be added or subtracted from the vertex
    !!
    !! The long-range part can be computed using Eq. (4) of PRL 115, 176401 (2015).
    !! The sum over G is converged using the Ewald summation technique (see for example
    !! F.2, p.500 in Martin Electronic structure book) where the Ewald factor is ((q+G)**2)/alph/4.0_DP.
    !!
    !! Technical note: From the solution of the Poisson equation, there is an additional factor
    !! e^{-i(q+G)\tau_\kappa} with respect to Eq. (4) of PRL 115, 176401 (2015).
    !! The full equation can be found in Eq. (S4) of the supplemental materials of PRL 115, 176401 (2015).
    !!
    !! The final implemented formula is:
    !!
    !! $$ g_{mn\nu}^{\mathcal L}({\bf k},{\bf q) = i\frac{4\pi e^2}{\Omega} \sum_{\kappa}
    !!   \left(\frac{\hbar}{2 {M_\kappa \omega_{{\bf q}\nu}}}\right)^{\!\!\frac{1}{2}}
    !!   \sum_{{\bf G}\ne -{\bf q}} e^{-({\bf q}+{\bf G})^2/4\alpha}
    !! \frac{ ({\bf q}+{\bf G})\cdot{\bf Z}^*_\kappa \cdot {\bf e}_{\kappa\nu}({\bf q}) }
    !!  {({\bf q}+{\bf G})\cdot\bm\epsilon^\infty\!\cdot({\bf q}+{\bf G})}\,
    !!   \left[ U_{{\bf k}+{\bf q}}\:U_{{\bf k}}^{\dagger} \right]_{mn} $$
    !!
    !! 10/2016 - SP: Optimization
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg, omega, alat
    USE ions_base,     ONLY : tau, nat
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12, zero, eps8
    USE epwcom,        ONLY : shortrange, nbndsub, lpolar, system_2d
    USE elph2,         ONLY : area, Qmat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: imode
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nmodes
    !! Max number of modes
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL(KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! effective charges tensor
    REAL(KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract long range term
    COMPLEX(KIND = DP), INTENT(in) :: uq(nmodes, nmodes)
    !! phonon eigenvec associated with q
    COMPLEX(KIND = DP), INTENT(inout) :: epmat(nbndsub, nbndsub)
    !! e-ph matrix elements
    COMPLEX(KIND = DP), INTENT(in) :: bmat(nbndsub, nbndsub)
    !! Overlap matrix elements $$<U_{mk+q}|U_{nk}>$$
    !
    ! Local variables
    INTEGER :: na
    !! Atom index 1
    INTEGER :: i
    !! Cartesian index direction
    INTEGER :: ipol, jpol
    !! Polarison
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: nr1x, nr2x, nr3x
    !! Minimum supercell size to include all vector such that G^2 < geg
    REAL(KIND = DP):: metric
    !! (2*pi/a)^2
    REAL(KIND = DP) :: qeq
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: arg
    !! Argument of the cos,sin for the Euler formula
    REAL(KIND = DP) :: zaq
    !! Z^* \cdot (q+g)
    REAL(KIND = DP) :: gg(3)
    !! G-vector
    !REAL(KIND = DP) :: gmax
    !!!  Max G-vector
    !REAL(KIND = DP) :: alph
    !!! Ewald factor (arbitrary, here chosen to be 1)
    REAL(KIND = DP) :: geg
    !! <G| epsil | G>
    REAL(KIND = DP) :: reff(2, 2)
    !! Effective screening length for 2D materials
    REAL(KIND = DP) :: grg
    !! G-vector * reff * G-vector
    REAL(KIND = DP) :: Qqq
    !! In the case of Si, its a single value
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    COMPLEX(KIND = DP) :: fac
    !! General prefactor
    COMPLEX(KIND = DP) :: facqd
    !! Ewald filtering
    COMPLEX(KIND = DP) :: facq
    !! Atomic position exponential
    COMPLEX(KIND = DP) :: epmatl(nbndsub, nbndsub)
    !! Long-range part of the matrix element
    !
    ! Impose zero Born effective charge in case of non-polar materials with quadrupoles
    IF (.NOT. lpolar) THEN
      zeu(:, :, :) = zero
    ENDIF
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore('rgd_blk_epw_fine_mem', ' wrong value for signe ', 1)
    !
    IF (system_2d) THEN
      ! Vacuum size in Bohr unit
      c = alat / bg(3, 3)
      ! (e^2 * 2\pi * ci) / Area
      fac = (signe * e2 * twopi * ci) / area
      ! Effective screening length
      ! reff = (epsil - 1) * c/2
      reff(:, :) = zero
      reff(:, :) = epsil(1:2, 1:2) * 0.5d0 * c ! eps * c/2 in 2pi/a units
      reff(1, 1) = reff(1, 1) - 0.5d0 * c ! (-1) * c/2 in 2pi/a units
      reff(2, 2) = reff(2, 2) - 0.5d0 * c ! (-1) * c/2 in 2pi/a units
    ELSE
      ! (e^2 * 4\pi * i) / Volume
      fac = (signe * e2 * fpi * ci) / omega
    ENDIF
    !
    !gmax = 14.d0
    !alph = 1.0d0
    metric = (twopi / alat)**2
    geg  = gmax * alph * 4.0d0
    !
    ! Estimate of nr1x, nr2x, nr3x generating all vectors up to G^2 < geg
    ! Only for dimensions where periodicity is present.
    IF (nqc1 == 1) THEN
      nr1x = 0
    ELSE
      nr1x = INT(SQRT(geg) / SQRT(bg(1, 1)**2 + bg(2, 1)**2 + bg(3, 1)**2)) + 1
    ENDIF
    IF (nqc2 == 1) THEN
      nr2x = 0
    ELSE
      nr2x = INT(SQRT(geg) / SQRT(bg(1, 2)**2 + bg(2, 2)**2 + bg(3, 2)**2)) + 1
    ENDIF
    IF (nqc3 == 1) THEN
      nr3x = 0
    ELSE
      nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
    ENDIF
    !
    epmatl(:, :) = czero
    DO m1 = -nr1x, nr1x
      DO m2 = -nr2x, nr2x
        DO m3 = -nr3x, nr3x
          !
          gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)) * (twopi / alat)
          gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)) * (twopi / alat)
          gg(3) = (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)) * (twopi / alat)
          !
          IF (system_2d) THEN
            qeq = gg(1)**2 + gg(2)**2 + gg(3)**2
            grg = zero
            IF (gg(1)**2 + gg(2)**2 > eps8) THEN
              grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + &
                    gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
              grg = grg / (gg(1)**2 + gg(2)**2)
            ENDIF
          ELSE
            !
            qeq = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
                   gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
                   gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
          ENDIF
          !
          IF (qeq > 0.0d0 .AND. qeq / (metric * alph * 4.0d0) < gmax) THEN
            !
            IF (system_2d) THEN
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / (SQRT(qeq) * (1.0 + grg * SQRT(qeq)))
            ELSE
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / qeq  !<-- this is correct
              !facqd = fac * EXP(-qeq * DSQRT(metric) / (metric * alph * 4.0d0)) / qeq ! <-- this is to keep as previous
            ENDIF
            !
            DO na = 1, nat
              arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na) + gg(3) * tau(3, na))
              facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
              ! Cartesian index direction
              DO i = 1, 3
                zaq = zero
                DO ipol = 1, 3
                  zaq = zaq + gg(ipol) * zeu(ipol, i, na)
                ENDDO
                !
                Qqq = zero
                DO ipol = 1, 3
                  DO jpol = 1, 3
                    Qqq = Qqq + 0.5 * gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
                  ENDDO
                ENDDO
                CALL ZAXPY(nbndsub**2, facq * (zaq - ci * Qqq) * uq(3 * (na - 1) + i, imode), &
                           bmat(:, :), 1, epmat(:, :), 1)
                CALL ZAXPY(nbndsub**2, facq * (zaq - ci * Qqq) * uq(3 * (na - 1) + i, imode), &
                           bmat(:, :), 1, epmatl(:, :), 1)
                !
              ENDDO !ipol
            ENDDO !nat
          ENDIF ! (qeq > 0.0_DP .AND. qeq / alph / 4.0_DP < gmax)
        ENDDO ! m3
      ENDDO ! m2
    ENDDO ! m1
    !
    ! In case we want only the short-range we do
    ! g_s = DSQRT(g*g - g_l*g_l)
    !
    ! Important notice: It is possible that (g*g - g_l*g_l) < 0, in which
    ! case the sqrt will give an pure imaginary number. If it is positive we
    ! will get a pure real number.
    ! In any case, when g_s will be squared both will become real numbers.
    IF (shortrange) THEN
      !epmat = ZSQRT(epmat*CONJG(epmat) - epmatl*CONJG(epmatl))
      epmat = SQRT(epmat * CONJG(epmat) - epmatl * CONJG(epmatl))
    ENDIF
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE rgd_blk_epw_fine_mem
    !-------------------------------------------------------------------------------
  !---------------------------------------------------------------------------------
  END MODULE rigid_epw
  !---------------------------------------------------------------------------------
