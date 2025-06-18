  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2001-2008 Quantum-Espresso group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE longrange
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
  REAL(KIND = DP), PARAMETER :: gmax = 14.d0
  !! Cutoff criteria for G-sum:
  !! e^{ - (q+G) \cdot \epsilon \cdot (q+G) / (4 \alpha) } < e^{ - gmax } terms are neglected.
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_min_qG(q)
    !-----------------------------------------------------------------------
    !!
    !! Find igmin_qG = min_{G}|q+G|
    !! to shift G-sum in long-range g
    !! in order to ensure periodicity g(q+G')=g(q)
    !!
    USE ep_constants,  ONLY : eps8
    USE cell_base,     ONLY : bg
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
    USE ep_constants,  ONLY : eps8
    USE cell_base,     ONLY : bg
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
    SUBROUTINE rgd_blk(L, nqc1, nqc2, nqc3, nat, dyn, q, tau, epsil, zeu, signe)
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
    !! SP: 05/2021 - Addition of out-of-plane 2d electrostatic with L
    !! SP: 11/2023 - Addition of dipole_sh. For now dipole_sh = dipole_sp but will
    !!               be implemented in the future.
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg, omega, alat, at
    USE ep_constants,  ONLY : eps6, ci, zero, czero, twopi, eps8, pi, fpi, e2
    USE io_global,     ONLY : ionode_id, ionode
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : world_comm, inter_pool_comm
    USE parallelism,   ONLY : para_bounds
    USE mp,            ONLY : mp_bcast, mp_sum
    USE input,         ONLY : lpolar, system_2d
    USE global_var,    ONLY : Qmat, alph
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
    REAL(KIND = DP), INTENT(in) :: L
    !! Range separation length in 2D
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
    LOGICAL :: criteria
    !! Criteria to neglect components
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
    REAL(KIND = DP) :: reff(2, 2)
    !! Effective screening length for 2D materials
    REAL(KIND = DP) :: grg
    !! G-vector * reff * G-vector
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    REAL(KIND = DP) :: qtmp(3)
    !! Temporary q vector to find min_{G}|q+G|
    REAL(KIND = DP) :: f
    !! Range-separation function
    REAL(KIND = DP) :: qnorm
    !! Norm of (G+q)
    REAL(KIND = DP) :: alpha_para
    !! Polarizability in plane
    REAL(KIND = DP) :: alpha_perp
    !! Polarizability out of plane
    REAL(KIND = DP) :: epsilon_para
    !! Dielectic function in plane
    REAL(KIND = DP) :: epsilon_perp
    !! Dielectric function out of plane
    COMPLEX(KIND = DP) :: zag_para(3)
    !! Z * G in-plane
    COMPLEX(KIND = DP) :: zag_perp(3)
    !! Z * G out of plane
    COMPLEX(KIND = DP) :: zbg_para(3)
    !! Z * G in-plane
    COMPLEX(KIND = DP) :: zbg_perp(3)
    !! Z * G out of plane
    COMPLEX(KIND = DP) :: zcg_para(3)
    !! Z * G in-plane
    COMPLEX(KIND = DP) :: zcg_perp(3)
    !! Z * G out of plane
    COMPLEX(KIND = DP) :: fnat_para(3)
    !! Z * G * cos(arg) in plane
    COMPLEX(KIND = DP) :: fnat_perp(3)
    !! Z * G * cos(arg) out of plane
    COMPLEX(KIND = DP) :: fnat(3)
    !! Z with \delta_kk' summed
    COMPLEX(KIND = DP) :: qnat(3)
    !! Q with \delta_kk' summed
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
    IF (.NOT. lpolar) zeu(:, :, :) = zero
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
    ! Vacuum size in Bohr unit
    c = alat / bg(3, 3)
    IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
      !
      ! (e^2 * 2\pi) / Area
      fac = (signe * e2 * twopi) * (c / omega)
      !
    ELSEIF (system_2d == 'gaussian') THEN
      ! (e^2 * 2\pi) / Area
      fac = (signe * e2 * twopi) * (c / omega)
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
    IF (ionode) THEN
      diff = mm_stop - mm_start
    ENDIF
    CALL mp_bcast(diff, ionode_id, inter_pool_comm)
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
      IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
        geg = gg(1)**2 + gg(2)**2
        qnorm = SQRT(geg)
        f = 1.0 - TANH(qnorm * L * 0.5)
        alpha_para = (gg(1) * (gg(1) * (epsil(1, 1) - 1.0) + gg(2) * epsil(2, 1)) + &
                      gg(2) * (gg(1) * epsil(1, 2) + gg(2) * (epsil(2, 2) - 1.0))) * c / fpi
        alpha_perp = (epsil(3, 3) - 1.0) * c / fpi
        criteria = geg > zero
        IF (criteria) THEN
          epsilon_para = 1.0 + twopi * f * alpha_para / qnorm
          epsilon_perp = 1.0 - twopi * qnorm * f * alpha_perp
        ENDIF
      ELSEIF (system_2d == 'gaussian') THEN
        geg = gg(1)**2 + gg(2)**2 + gg(3)**2
        grg = zero
        IF (gg(1)**2 + gg(2)**2 > eps8) THEN
          grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
          grg = grg / (gg(1)**2 + gg(2)**2)
        ENDIF
        criteria = geg > 0.0d0 .AND. geg / (metric * alph * 4.0d0) < gmax
      ELSE
        !
        geg = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
               gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
               gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
        criteria = geg > 0.0d0 .AND. geg / (metric * alph * 4.0d0) < gmax
      ENDIF
      !
      IF (criteria) THEN
        !
        IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
          facgd = fac * f / qnorm
        ELSEIF (system_2d == 'gaussian') THEN
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / (SQRT(geg) * (1.0 + grg * SQRT(geg)))
        ELSE
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / geg
        ENDIF
        !
        IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
          DO na = 1, nat
            zag_para(:) = zero
            zag_perp(:) = zero
            DO i = 1, 3
              ! Dipole
              DO ipol = 1, 2
                zag_para(i) = zag_para(i) + gg(ipol) * zeu(ipol, i, na)
              ENDDO
              zag_perp(i) = qnorm * zeu(3, i, na)
              !
              ! quadrupole
              DO ipol = 1, 2
                DO jpol = 1, 2
                  zag_para(i) = zag_para(i) - 0.5 * ci *  gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
                ENDDO
              ENDDO
              zag_para(i) = zag_para(i) + 0.5 * ci * Qmat(na, i, 3, 3) * qnorm * qnorm
              DO ipol = 1, 2
                zag_perp(i) = zag_perp(i) - ci * qnorm * gg(ipol) * Qmat(na, i, 3, ipol)
              ENDDO
            ENDDO
            fnat_para(:) = czero
            fnat_perp(:) = czero
            DO nb = 1, nat
              ! alat is needed to express tau in Bohr
              arg = alat * (gg(1) * (tau(1, na) - tau(1, nb)) + &
                            gg(2) * (tau(2, na) - tau(2, nb)))
              zcg_para(:) = czero
              zcg_perp(:) = czero
              DO i = 1, 3
                ! Dipole
                DO ipol = 1, 2
                  zcg_para(i) = zcg_para(i) + gg(ipol) * zeu(ipol, i, nb)
                ENDDO
                zcg_perp(i) = qnorm * zeu(3, i, nb)
                !
                ! quadrupole
                DO ipol = 1, 2
                  DO jpol = 1, 2
                    zcg_para(i) = zcg_para(i) - 0.5 * ci * gg(ipol) * gg(jpol) * Qmat(nb, i, ipol, jpol)
                  ENDDO
                ENDDO
                zcg_para(i) = zcg_para(i) + 0.5 * ci * Qmat(nb, i, 3, 3) * qnorm * qnorm
                DO ipol = 1, 2
                  zcg_perp(i) = zcg_perp(i) - ci * qnorm * gg(ipol) * Qmat(nb, i, 3, ipol)
                ENDDO
              ENDDO ! i
              fnat_para(:) = fnat_para(:) + zcg_para(:) * CMPLX(COS(arg), SIN(arg), KIND=DP)
              fnat_perp(:) = fnat_perp(:) + zcg_perp(:) * CMPLX(COS(arg), SIN(arg), KIND=DP)
            ENDDO !nb
            DO j = 1, 3
              DO i = 1, 3
                dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) - facgd * &
                   ((CONJG(zag_para(i)) * fnat_para(j)) / epsilon_para - (CONJG(zag_perp(i)) * fnat_perp(j)) / epsilon_perp)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! nat
        ELSE
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
                ! Dipole-quadrupole
                !Qdq = 0.5d0 * (zag(i) * qnat(j) - fnat(i) * qag(i))
                Qdq = 0.5d0 * (qag(i) * fnat(j) - zag(i) * qnat(j))
                ! Quadrupole-quadrupole
                Qqq = 0.25d0 * qag(i) * qnat(j)
                !
                dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) &
                                             - facgd * (Qdd + ci * Qdq + Qqq)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! nat
        ENDIF ! system_2d .AND. L > 0.001
      ENDIF ! geg
      !
      ! Case q =/ 0
      !gg(1) = gg(1) + q(1) * (twopi / alat)
      !gg(2) = gg(2) + q(2) * (twopi / alat)
      !gg(3) = gg(3) + q(3) * (twopi / alat)
      !
      !JLB: shift G-sum and center around min_{G}|q+G|, to ensure periodicity
      gg(1) = ((m1 + igmin_qG(1)) * bg(1, 1) + (m2 + igmin_qG(2)) * bg(1, 2) + (m3 + igmin_qG(3)) * bg(1,3) + q(1)) * (twopi / alat)
      gg(2) = ((m1 + igmin_qG(1)) * bg(2, 1) + (m2 + igmin_qG(2)) * bg(2, 2) + (m3 + igmin_qG(3)) * bg(2,3) + q(2)) * (twopi / alat)
      gg(3) = ((m1 + igmin_qG(1)) * bg(3, 1) + (m2 + igmin_qG(2)) * bg(3, 2) + (m3 + igmin_qG(3)) * bg(3,3) + q(3)) * (twopi / alat)
      !
      IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
        geg = gg(1)**2 + gg(2)**2
        qnorm = SQRT(geg)
        f = 1.0 - TANH(qnorm * L * 0.5)
        alpha_para = (gg(1) * (gg(1) * (epsil(1, 1) - 1.0) + gg(2) * epsil(2, 1))  + &
                      gg(2) * (gg(1) * epsil(1, 2) + gg(2) * (epsil(2, 2) - 1.0))) * c / fpi
        alpha_perp = (epsil(3, 3) - 1.0) * c / fpi
        criteria = geg > zero
        IF (criteria) THEN
          epsilon_para = 1.0 + twopi * f * alpha_para / qnorm
          epsilon_perp = 1.0 - twopi * qnorm * f * alpha_perp
        ENDIF
      ELSEIF (system_2d == 'gaussian') THEN
        geg = gg(1)**2 + gg(2)**2 + gg(3)**2
        grg = zero
        IF (gg(1)**2 + gg(2)**2 > eps8) THEN
          grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
          grg = grg / (gg(1)**2 + gg(2)**2)
        ENDIF
        criteria = geg > zero .AND. geg / (metric * alph * 4.0d0) < gmax
      ELSE
        geg = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
               gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
               gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
        criteria = geg > zero .AND. geg / (metric * alph * 4.0d0) < gmax
      ENDIF
      !
      IF (criteria) THEN
        !
        IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
          facgd = fac * f / qnorm
        ELSEIF (system_2d == 'gaussian') THEN
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / (SQRT(geg) * (1.0 + grg * SQRT(geg)))
        ELSE
          facgd = fac * EXP(-geg / (metric * alph * 4.0d0)) / geg
        ENDIF
        !
        IF (system_2d == 'dipole_sp' .OR. system_2d == 'dipole_sh' .OR. system_2d == 'quadrupole') THEN
          DO nb = 1, nat
            zbg_para(:) = czero
            zbg_perp(:) = czero
            DO i = 1, 3
              ! Dipole
              DO ipol = 1, 2
                zbg_para(i) = zbg_para(i) + gg(ipol) * zeu(ipol, i, nb)
              ENDDO
              zbg_perp(i) = qnorm * zeu(3, i, nb)
              !
              ! Quadrupole
              DO ipol = 1, 2
                DO jpol = 1, 2
                  zbg_para(i) = zbg_para(i) - 0.5 * ci *  gg(ipol) * gg(jpol) * Qmat(nb, i, ipol, jpol)
                ENDDO
              ENDDO
              zbg_para(i) = zbg_para(i) + 0.5 * ci * Qmat(nb, i, 3, 3) * qnorm * qnorm
              !
              DO ipol = 1, 2
                zbg_perp(i) = zbg_perp(i) - ci * qnorm * gg(ipol) * Qmat(nb, i, 3, ipol)
              ENDDO
            ENDDO ! i
            !
            DO na = 1, nat
              ! alat is needed to express tau in Bohr unit
              arg = alat * (gg(1) * (tau(1, na) - tau(1, nb)) + &
                            gg(2) * (tau(2, na) - tau(2, nb)))
              zag_para(:) = czero
              zag_perp(:) = czero
              DO i = 1, 3
                ! Dipole
                DO ipol = 1, 2
                  zag_para(i) = zag_para(i) + gg(ipol) * zeu(ipol, i, na)
                ENDDO
                zag_perp(i) = qnorm * zeu(3, i, na)
                !
                ! Quadrupole
                DO ipol = 1, 2
                  DO jpol = 1, 2
                    zag_para(i) = zag_para(i) - 0.5 * ci * gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
                  ENDDO
                ENDDO
                zag_para(i) = zag_para(i) + 0.5 * ci * Qmat(na, i, 3, 3) * qnorm * qnorm
                !
                DO ipol = 1, 2
                  zag_perp(i) = zag_perp(i) - ci * qnorm * gg(ipol) * Qmat(na, i, 3, ipol)
                ENDDO
              ENDDO ! i
              !
              facg = facgd * CMPLX(COS(arg), SIN(arg), KIND=DP)
              DO j = 1, 3
                DO i = 1, 3
                  dyn_tmp((na - 1) * 3 + i, (nb - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (nb - 1) * 3 + j) + facg * &
                     ((CONJG(zag_para(i)) * zbg_para(j)) / epsilon_para - (CONJG(zag_perp(i)) * zbg_perp(j)) / epsilon_perp)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! na
          ENDDO ! nb
        ELSE ! system_2d .AND. L > 0.001
          !
          DO nb = 1, nat ! kappa
            DO na = 1, nat ! kappa'
              arg = alat * (gg(1) * (tau(1, na) - tau(1 ,nb)) + &
                            gg(2) * (tau(2, na) - tau(2, nb)) + &
                            gg(3) * (tau(3, na) - tau(3, nb)))
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
        ENDIF ! system_2d .AND. L > 0.001
      ENDIF ! criteria
    ENDDO ! mm
    !
    CALL mp_sum(dyn_tmp, inter_pool_comm)
    dyn(:, :) = dyn_tmp(:, :) + dyn(:, :)
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE rgd_blk
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE rgd_blk_epw(nqc1, nqc2, nqc3, q, epmat, nbndsub, nmodes, epsil, zeu, A, signe)
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
    !! In the case of 2D materials
    !!   for dipole_sp and quadrupole, the implemented formula are Eqs. 43-47 of PRB 107, 155424 (2023).
    !!   for dipole_sh, the implemented formula is Eq. 34 of PRB 105, 115414 (2022).
    !!
    !! 10/2016 - SP: Optimization
    !! 03/2024 - JML: Optimization, merge with rgd_blk_epw_fine
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : bg, omega, alat, at
    USE ions_base,     ONLY : tau, nat
    USE ep_constants,  ONLY : twopi, fpi, e2, ci, czero, eps12, zero, eps8, cone, &
                              zero, one, half, two, four, eps5
    USE input,         ONLY : shortrange, lpolar, system_2d
    USE global_var,    ONLY : Qmat, L, do_cutoff_2D_epw, alph, qrpl
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
    INTEGER, INTENT(in) :: nbndsub
    !! Number of bands in optimal space
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid in Cartesian coordinates
    REAL(KIND = DP), INTENT(in) :: epsil(3, 3)
    !! dielectric constant tensor
    REAL(KIND = DP), INTENT(inout) :: zeu(3, 3, nat)
    !! effective charges tensor
    REAL(KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract long range term
    COMPLEX(KIND = DP), INTENT(inout) :: epmat(nbndsub, nbndsub, nmodes)
    !! e-ph matrix elements
    COMPLEX(KIND = DP), INTENT(inout) :: A(3, nbndsub, nbndsub)
    !! Berry connection (position matrix elements)
    !
    ! Local variables
    LOGICAL :: criteria
    !! Criteria to neglect components
    INTEGER :: na
    !! Atom index 1
    INTEGER :: i
    !! Cartesian index direction
    INTEGER :: ipol, jpol
    !! Polarization direction
    INTEGER :: imode
    !! Mode index
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: nr1x, nr2x, nr3x
    !! Minimum supercell size to include all vector such that G^2 < geg
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: mmin(3), mmax(3)
    !! Shifted G-loop to be centered around min_{G}|q+G|
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
    REAL(KIND = DP) :: geg
    !!  <q+G| epsil | q+G>
    REAL(KIND = DP) :: reff(2, 2)
    !! Effective screening length for 2D materials
    REAL(KIND = DP) :: grg
    !! G-vector * reff * G-vector
    REAL(KIND = DP) :: Qqq
    !! Quadrupole for given q
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    REAL(KIND = DP) :: qtmp(3)
    !! Temporary q vector to find min_{G}|q+G|
    REAL(KIND = DP) :: f
    !! 2D range separation function (envelope function)
    REAL(KIND = DP) :: qnorm
    !! Norm of (G+q)
    REAL(KIND = DP) :: alpha_para
    !! Polarizability in plane
    REAL(KIND = DP) :: alpha_perp
    !! Polarizability out of plane
    REAL(KIND = DP) :: epsilon_para
    !! Dielectic function in plane
    REAL(KIND = DP) :: epsilon_perp
    !! Dielectric function out of plane
    REAL(KIND = DP) :: eps_inf_env
    !! Epsilon inf in vacuum  = 1
    REAL(KIND = DP) :: gpsq
    !! In-plane norm of G
    COMPLEX(KIND = DP) :: zaqc
    !! Complex Z^* \cdot (q+g) for 2D
    COMPLEX(KIND = DP) :: zbg_para(3)
    !! Z * G in-plane
    COMPLEX(KIND = DP) :: zbg_perp(3)
    !! Z * G out of plane
    COMPLEX(KIND = DP) :: zbg_dip(3)
    !! Dipole inplane only
    COMPLEX(KIND = DP) :: fac
    !! General prefactor
    COMPLEX(KIND = DP) :: facqd
    !! Ewald filtering
    COMPLEX(KIND = DP) :: facq
    !! Atomic position exponential
    COMPLEX(KIND = DP) :: epmatl(nbndsub, nbndsub, nmodes)
    !! Long-range part of the el-ph matrix elements
    COMPLEX(KIND = DP) :: delta(nbndsub, nbndsub)
    !! Kronecker delta
    COMPLEX(KIND = DP) :: sk
    !! 2D Kernel
    COMPLEX(KIND = DP) :: sdk
    !! Derivative of 2D Kernel
    COMPLEX(KIND = DP) :: coeff_delta(nmodes)
    !! Coefficient of the delta(iw, jw) term
    COMPLEX(KIND = DP) :: coeff_r(3, nmodes)
    !! Coefficient of the rrk(3, iw, jw) term
    !
    CALL start_clock('rgd_blk_epw')
    !
    ! Impose zero Born effective charge in case of non-polar materials with quadrupoles
    IF (.NOT. lpolar) zeu(:, :, :) = zero
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore('rgd_blk_epw', 'Wrong value for signe ', 1)
    !
    ! Vacuum size in Bohr unit
    c = alat / bg(3, 3)
    IF (system_2d == 'dipole_sp' .OR. system_2d == 'quadrupole') THEN
      ! (e^2 * 2\pi * i) / Area
      fac = (e2 * twopi) * (c / omega)
    ELSEIF (system_2d == 'dipole_sh') THEN
      eps_inf_env = one
      fac = (e2 * fpi / omega) * half / eps_inf_env
    ELSEIF (system_2d == 'gaussian') THEN
      ! (e^2 * 2\pi * i) / Area
      fac = (e2 * twopi) * (c / omega)
      ! Effective screening length
      ! reff = (epsil - 1) * c/2
      reff(:, :) = zero
      reff(:, :) = epsil(1:2, 1:2) * 0.5d0 * c ! eps * c/2 in 2pi/a units
      reff(1, 1) = reff(1, 1) - 0.5d0 * c ! (-1) * c/2 in 2pi/a units
      reff(2, 2) = reff(2, 2) - 0.5d0 * c ! (-1) * c/2 in 2pi/a units
    ELSE
      ! (e^2 * 4\pi * i) / Volume
      fac = (e2 * fpi) / omega
    ENDIF
    !
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
    IF (system_2d == 'dipole_sh') THEN
      nr3x = INT(SQRT(geg) / SQRT(bg(1, 3)**2 + bg(2, 3)**2 + bg(3, 3)**2)) + 1
    ENDIF
    !
    !JLB: Find igmin_qG = min_{G}|q+G|
    qtmp = q
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
    IF (do_cutoff_2D_epw .AND. system_2d == 'dipole_sh') THEN
      mmin(3) = 0
      mmax(3) = 0
      IF (ABS(q(3)) > eps12) CALL errore('rgd_blk_epw', 'Error q_z must be zero in 2D limit', 1)
    ENDIF
    !
    coeff_delta(:) = czero
    coeff_r(:, :) = czero
    !
    DO m1 = mmin(1), mmax(1)
      DO m2 = mmin(2), mmax(2)
        DO m3 = mmin(3), mmax(3)
          !
          gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)) * (twopi / alat)
          gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)) * (twopi / alat)
          gg(3) = (m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)) * (twopi / alat)
          !
          IF (system_2d == 'dipole_sp' .OR. system_2d == 'quadrupole') THEN
            qeq = gg(1)**2 + gg(2)**2
            qnorm = SQRT(qeq)
            f = 1.0 - TANH(qnorm * L * 0.5)
            alpha_para = (gg(1) * (gg(1) * (epsil(1, 1) - 1.0) + gg(2) * epsil(2, 1)) + &
                          gg(2) * (gg(1) * epsil(1, 2) + gg(2) * (epsil(2, 2) - 1.0))) * c / fpi
            alpha_perp = (epsil(3, 3) - 1.0) * c / fpi
            criteria = qeq > zero
            IF (criteria) THEN
              epsilon_para = 1.0 + twopi * f * alpha_para / qnorm
              epsilon_perp = 1.0 - twopi * qnorm * f * alpha_perp
            ENDIF
          ELSEIF (system_2d == 'dipole_sh') THEN
            gg(1) = (m1 * bg(1, 1) + m2 * bg(1, 2) + q(1)) * (twopi / alat)
            gg(2) = (m1 * bg(2, 1) + m2 * bg(2, 2) + q(2)) * (twopi / alat)
            gg(3) = (m3 * bg(3, 3) + q(3)) * (twopi / alat)
            qeq = gg(1)**2 + gg(2)**2 + gg(3)**2
            criteria = qeq > zero .AND. qeq / (metric * alph * four) < gmax
          ELSEIF (system_2d == 'gaussian') THEN
            qeq = gg(1)**2 + gg(2)**2 + gg(3)**2
            grg = zero
            IF (gg(1)**2 + gg(2)**2 > eps8) THEN
              grg = gg(1) * reff(1, 1) * gg(1) + gg(1) * reff(1, 2) * gg(2) + &
                    gg(2) * reff(2, 1) * gg(1) + gg(2) * reff(2, 2) * gg(2)
              grg = grg / (gg(1)**2 + gg(2)**2)
            ENDIF
            criteria = qeq > zero .AND. qeq / (metric * alph * four) < gmax
          ELSE
            qeq = (gg(1) * (epsil(1, 1) * gg(1) + epsil(1, 2) * gg(2) + epsil(1, 3) * gg(3)) + &
                   gg(2) * (epsil(2, 1) * gg(1) + epsil(2, 2) * gg(2) + epsil(2, 3) * gg(3)) + &
                   gg(3) * (epsil(3, 1) * gg(1) + epsil(3, 2) * gg(2) + epsil(3, 3) * gg(3)))
            criteria = qeq > zero .AND. qeq / (metric * alph * four) < gmax
          ENDIF
          !
          IF (criteria) THEN
            !
            ! nGtest = nGtest + 1
            !
            IF (system_2d == 'dipole_sp' .OR. system_2d == 'quadrupole') THEN
              facqd = fac * f / qnorm
            ELSEIF (system_2d == 'dipole_sh') THEN
              gpsq = NORM2((/gg(1), gg(2), zero/))
              IF (gpsq < eps12) gpsq = eps5
              facqd = fac * EXP(-qeq / (metric * alph * four)) / gpsq
            ELSEIF (system_2d == 'gaussian') THEN
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / (SQRT(qeq) * (1.0 + grg * SQRT(qeq)))
            ELSE
              facqd = fac * EXP(-qeq / (metric * alph * 4.0d0)) / qeq  !<-- this is correct
              !facqd = fac * EXP(-qeq * DSQRT(metric) / (metric * alph * 4.0d0)) / qeq ! <-- this is to keep as previous
              !facqd = fac / qeq ! no Ewald
            ENDIF
            !
            IF (system_2d == 'dipole_sp') THEN
              DO na = 1, nat
                arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na))
                facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
                zbg_para(:) = czero
                zbg_perp(:) = czero
                DO i = 1, 3
                  imode = 3 * (na - 1) + i
                  !
                  ! Dipole
                  DO ipol = 1, 2
                    zbg_para(i) = zbg_para(i) + ci * gg(ipol) * zeu(ipol, i, na)
                  ENDDO
                  !
                  coeff_delta(imode) = coeff_delta(imode) + facq * (zbg_para(i) / epsilon_para)
                  !
                ENDDO ! i
              ENDDO ! na
            ELSEIF (system_2d == 'quadrupole') THEN
              ! Here quadrupole means that we are adding dipole (first order) + quadrupole + Berry connection (second order)
              DO na = 1, nat
                arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na))
                facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
                zbg_para(:) = czero
                zbg_perp(:) = czero
                zbg_dip(:) = czero
                DO i = 1, 3
                  imode = 3 * (na - 1) + i
                  !
                  ! Dipole
                  DO ipol = 1, 2
                    zbg_para(i) = zbg_para(i) + ci * gg(ipol) * zeu(ipol, i, na)
                    zbg_dip(i)  = zbg_dip(i) + gg(ipol) * zeu(ipol, i, na)
                  ENDDO
                  zbg_perp(i) = qnorm * zeu(3, i, na)
                  !
                  ! Quadrupole
                  DO ipol = 1, 2
                    DO jpol = 1, 2
                      zbg_para(i) = zbg_para(i) + 0.5 *  gg(ipol) * gg(jpol) * Qmat(na, i, ipol, jpol)
                    ENDDO
                  ENDDO
                  zbg_para(i) = zbg_para(i) - 0.5 * Qmat(na, i, 3, 3) * qnorm * qnorm
                  !
                  coeff_delta(imode) = coeff_delta(imode) + facq * (zbg_para(i) / epsilon_para)
                  coeff_r(1, imode) = coeff_r(1, imode) - facq * (zbg_dip(i) / epsilon_para) * gg(1)
                  coeff_r(2, imode) = coeff_r(2, imode) - facq * (zbg_dip(i) / epsilon_para) * gg(2)
                  coeff_r(3, imode) = coeff_r(3, imode) + facq * qnorm * (zbg_perp(i) / epsilon_perp)
                  !
                ENDDO ! i
              ENDDO ! na
            ELSEIF (system_2d == 'dipole_sh') THEN
              DO na = 1, nat
                arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na))
                facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
                IF (do_cutoff_2D_epw) THEN
                  CALL kernel_2d((/gg(1), gg(2), 0.d0/), tau(3, na), sk, sdk)
                ELSE
                  CALL kernel_slab((/gg(1), gg(2), gg(3)/), q(3), m3, na, sk, sdk)
                ENDIF
                ! Cartesian direction
                DO i = 1, 3
                  imode = 3 * (na - 1) + i
                  zaqc = (gg(1) * zeu(1, i, na) + gg(2) * zeu(2, i, na)) * sk * ci - zeu(3, i, na) * sdk
                  ! SP - FIXME - add support for quadrupoles
                  !
                  coeff_delta(imode) = coeff_delta(imode) + facq * zaqc
                ENDDO !ipol
                !
              ENDDO !nat
            ELSE
              DO na = 1, nat
                arg = - alat * (gg(1) * tau(1, na) + gg(2) * tau(2, na) + gg(3) * tau(3, na))
                facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
                ! Cartesian direction
                DO i = 1, 3
                  imode = 3 * (na - 1) + i
                  !
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
                  coeff_delta(imode) = coeff_delta(imode) + facq * (ci * zaq + Qqq)
                  IF (qrpl) coeff_r(:, imode) = coeff_r(:, imode) - facq * zaq * gg(:)
                  !
                ENDDO !ipol
              ENDDO !nat
            ENDIF ! system_2d
          ENDIF ! criteria
        ENDDO ! m3
      ENDDO ! m2
    ENDDO ! m1
    !
    epmatl(:, :, :) = czero
    !
    ! Add the \delta_{iw, jw} term.
    delta(:, :) = czero
    DO ibnd = 1, nbndsub
      delta(ibnd, ibnd) = cone  ! Kronecker delta
    ENDDO
    DO imode = 1, nmodes
      CALL ZAXPY(nbndsub**2, coeff_delta(imode), delta, 1, epmatl(:, :, imode), 1)
    ENDDO
    !
    ! Add the A_{ipol, iw, jw} term.
    ! (3, nbndsub^2)^T * (3, nmodes) -> (nbndsub^2, nmodes)
    CALL ZGEMM('T', 'N', nbndsub**2, nmodes, 3, cone, A, 3, coeff_r, 3, cone, epmatl, nbndsub**2)
    !
    ! Add or subtract long-range part to epmat
    !
    epmat = epmat + signe * epmatl
    !
    ! WRITE(*,*) "Number of G-vectors within cutoff:", nGtest
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
    CALL stop_clock('rgd_blk_epw')
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE rgd_blk_epw
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE rgd_imp_epw_fine(nqc1, nqc2, nqc3, q, epmat, cufkk, cufkq, signe)
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
    USE ep_constants,  ONLY : twopi, fpi, e2, ci, czero, eps12, cc2cb, cone
    USE input,         ONLY : nbndsub, ii_charge
    USE global_var,    ONLY : alph, epsi_s
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc1
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc2
    !! Coarse q-point grid
    INTEGER, INTENT(in) :: nqc3
    !! Coarse q-point grid
    REAL (KIND = DP), INTENT(in) :: q(3)
    !! q-vector from the full coarse or fine grid.
    REAL (KIND = DP), INTENT(in) :: signe
    !! signe=+/-1.0 ==> add/subtract long range term
    COMPLEX (KIND = DP), INTENT(inout) :: epmat(nbndsub, nbndsub)
    !! e-ph matrix elements
    COMPLEX(KIND = DP), INTENT(in) :: cufkk(nbndsub, nbndsub)
    !! rotation matrix U(k)^\dagger, fine mesh
    COMPLEX(KIND = DP), INTENT(in) :: cufkq(nbndsub, nbndsub)
    !! rotation matrix U(k+q)^\dagger, fine mesh
    !
    ! Local variables
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    REAL(KIND = DP) :: qeq
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: g1, g2, g3
    !!
    REAL(KIND = DP) :: gmax
    !!
    REAL(KIND = DP) :: geg
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
    COMPLEX(KIND = DP) :: fac
    !! Pre-factor
    COMPLEX(KIND = DP) :: facqd
    !! Exp factor
    COMPLEX(KIND = DP) :: facq
    !! Exp factor
    COMPLEX(KIND = DP) :: bmatf(nbndsub, nbndsub)
    ! overlap wfcs in Bloch representation, fine grid
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore ('rgd_imp_epw_fine', 'Wrong value for signe ', 1)
    !
    ! U(k'+q') * U(k')^\dagger
    bmatf = czero
    CALL ZGEMM('n', 'c', nbndsub, nbndsub, nbndsub, cone, cufkq, nbndsub,&
                    cufkk, nbndsub, czero, bmatf, nbndsub)
    !
    gmax = 14.d0
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
          qeq = (g1 * (epsi_s(1, 1) * g1 + epsi_s(1, 2) * g2 + epsi_s(1, 3) * g3) + &
                 g2 * (epsi_s(2, 1) * g1 + epsi_s(2, 2) * g2 + epsi_s(2, 3) * g3) + &
                 g3 * (epsi_s(3, 1) * g1 + epsi_s(3, 2) * g2 + epsi_s(3, 3) * g3)) !*twopi/alat
          !
          IF (qeq > 0.0d0) THEN !  .AND. qeq / (alph * 4.0d0) < gmax) THEN
            !
            qeq = qeq * (twopi / alat)**2.0d0
            facqd = fac / qeq !fac * EXP(-qeq / (alph * 4.0d0)) / qeq !/(two*wq)
            !
            facq = facqd * 1.0d0 ! For now add impurty at position 0, 0, 0, JL
            CALL ZAXPY(nbndsub**2, facq , bmatf(:, :), 1, epmat(:, :), 1)
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
    USE ep_constants,  ONLY : fpi, e2, ci, twopi
    USE ep_constants,  ONLY : pi
    USE cell_base,     ONLY : bg, omega, alat
    USE ep_constants,  ONLY : eps6
    USE global_var,    ONLY : alph
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
    ! S. Tiwari: alph has changed for IFC read and should be read directly from
    ! elph2
    !gmax = 14.d0
    metric = (twopi / alat)**2
    geg = gmax * alph * 4.0d0
    !
    IF (ABS(ABS(signe) - 1.0) > eps6) CALL errore('rgd_blk_der', ' wrong value for signe ', 1)
    !
    !gmax = 14.d0
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
    SUBROUTINE epsi_thickn_2d()
    !-------------------------------------------------------------------------------
    !!
    !! Compute effective epsilon_inf and thickness for 2D system
    !! Using Eq. (64-66) of [PRB 105, 115414 (2022)].
    !! Added by V.-A. Ha (2023)
    !! Clean SP (2023)
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, alat, bg
    USE ions_base,     ONLY : tau, nat
    USE global_var,    ONLY : epsi, thickn_2d, epsi_2d, tz_ref, emib3tz, eib3d, do_cutoff_2D_epw
    USE ep_constants,  ONLY : ci, twopi, eps6, eps8, one, zero, two, four
    USE Coul_cut_2D,   ONLY : do_cutoff_2D
    USE input,         ONLY : epbwrite, epwread
    USE io_global,     ONLY : stdout
    USE mp,            ONLY : mp_bcast
    !
    IMPLICIT NONE
    !
    INTEGER :: INFO
    !! FIXME
    INTEGER :: LWORK
    !! FIXME
    INTEGER :: I
    !! FIXME
    INTEGER :: J(1)
    !! FIXME
    INTEGER :: ierr
    !! Error status
    INTEGER :: na
    !! Atom index
    INTEGER, PARAMETER :: LWMAX = 500
    !! FIXME
    REAL(KIND = DP) :: VL(3,3)
    !! FIXME
    REAL(KIND = DP) :: VR(3,3)
    !! FIXME
    REAL(KIND = DP) :: WR(3)
    !! FIXME
    REAL(KIND = DP) :: WI(3)
    !! FIXME
    REAL(KIND = DP) :: WORK(LWMAX)
    !! FIXME
    REAL(KIND = DP) :: UV(3,3)
    !! FIXME
    REAL(KIND = DP) :: AA(3,3)
    !! FIXME
    REAL(KIND = DP) :: DDOT
    !! FIXME
    REAL(KIND = DP) :: epsi_xyz(3)
    !! FIXME
    REAL(KIND = DP) :: epsi_xy
    !! FIXME
    REAL(KIND = DP) :: epsi_vec(3, 3)
    !! FIXME
    REAL(KIND = DP) :: sdot(3)
    !! FIXME
    REAL(KIND = DP) :: tt
    !! Effective epsilon
    REAL(KIND = DP) :: kk
    !! FIXME
    REAL(KIND = DP) :: bb
    !! FIXME
    REAL(KIND = DP) :: delta
    !! FIXME
    REAL(KIND = DP), PARAMETER :: epsi_xy_env = one
    !! FIXME
    REAL(KIND = DP), PARAMETER ::  epsi_z_env = one
    !! FIXME
    !
    IF ((.NOT. epwread) .AND. epbwrite) do_cutoff_2D_epw = do_cutoff_2D
    ! Get eigenvalues of epsi matrix
    UV = RESHAPE((/one, zero, zero, zero, one, zero, zero, zero, one/), (/3,3/))
    AA = epsi
    LWORK = -1
    CALL DGEEV( 'N', 'V', 3, AA, 3, WR, WI, VL, 3, VR, 3, WORK, LWORK, INFO)
    LWORK = MIN( LWMAX, INT(WORK(1)) )
    CALL DGEEV( 'N', 'V', 3, AA, 3, WR, WI, VL, 3, VR, 3, WORK, LWORK, INFO)
    IF (INFO > zero) CALL errore('epsi_thickn_2d', 'Error diagonalizing epsilon_inf matrix', 1)
    IF (ABS(SUM(WI)) > eps6) CALL errore('epsi_thickn_2d', 'Eigenvalues of epsilon_inf are complex', 1)
    !
    DO I = 1, 3
      sdot(1) = DDOT(3, VR(I,:), 1, UV(1,:), 1)
      sdot(2) = DDOT(3, VR(I,:), 1, UV(2,:), 1)
      sdot(3) = DDOT(3, VR(I,:), 1, UV(3,:), 1)
      J = MAXLOC(ABS(sdot))
      epsi_xyz(J(1)) = WR(I)
      epsi_vec(J(1),:) = VR(I,:)
    ENDDO
    WRITE(stdout, '(5x,a)') ' '
    WRITE(stdout, '(5x,a)') 'eigenvalues and eigenvectors of epsilon_inf'
    WRITE(stdout, '(10x,4f16.10)') epsi_xyz(1), epsi_vec(1,:)
    WRITE(stdout, '(10x,4f16.10)') epsi_xyz(2), epsi_vec(2,:)
    WRITE(stdout, '(10x,4f16.10)') epsi_xyz(3), epsi_vec(3,:)
    !
    epsi_xy = (epsi_xyz(1) + epsi_xyz(2)) / two
    ! Compute effective epsilon_inf and thickness
    tt = epsi_xy - epsi_xy_env
    kk = one / epsi_xyz(3) - one / epsi_z_env
    bb = tt - kk * epsi_xy_env * epsi_z_env
    delta = bb * bb + four * kk * tt * epsi_z_env * epsi_z_env
    IF (delta < zero) CALL errore('epsi_thickn_2d', 'No root for effective 2D epsilon_inf', 1)
    epsi_2d = (-bb - DSQRT(delta)) / (two * kk * epsi_z_env)
    thickn_2d = alat * at(3,3) * (epsi_xy - epsi_xy_env) / (epsi_2d - epsi_xy_env)
    WRITE(stdout, '(5x, a, f16.10)') 'effective 2D epsilon_inf ', epsi_2d
    WRITE(stdout, '(5x, a, f16.10, a)') 'effective 2D thickness ', thickn_2d, ' [Bohr]'
    IF ((epsi_2d < one) .OR. (thickn_2d < eps8)) THEN
      CALL errore('epsi_thickn_2d', 'Invalid effective 2D parameters', 1)
    ENDIF
    tz_ref = (MAXVAL(tau(3, :)) + MINVAL(tau(3,:)) + thickn_2d / alat) / two
    !
    IF (.NOT. do_cutoff_2D_epw) THEN
      ALLOCATE(emib3tz(nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('epsi_thickn_2d', 'Error allocating emib3tz', 1)
      DO na = 1, nat
        emib3tz(na) = EXP(-ci * twopi * bg(3, 3) * (tau(3, na) - tz_ref))
      ENDDO
      eib3d = EXP(ci * twopi * bg(3, 3) * thickn_2d / alat)
    ENDIF
    !
    RETURN
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE epsi_thickn_2d
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE kernel_2d(q, tauz, kern, dkern)
    !------------------------------------------------------------------------------
    !! Compute kernel K in eq. (34) and its derivative in 2D limit using
    !! eq. (37-39) [PRB 105, 115414 (2022)]. (assume_isolated = '2D')
    !! Coded by W.-H. Sio (2022) then modified and merged by V.-A. Ha (2023)
    !! Clean/adapted by SP (2023)
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : pi, twopi, zero, half, one, two, four, eps5, eps12
    USE cell_base,     ONLY : at, alat
    USE global_var,    ONLY : thickn_2d, epsi_2d, tz_ref, alph
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: tauz
    !! Atom position
    REAL(KIND = DP), INTENT(in) :: q(3)
    !! Q-point
    COMPLEX(KIND = DP), INTENT(out) :: kern
    !! 2D Kernel
    COMPLEX(KIND = DP), INTENT(out) :: dkern
    !! Derivative of 2D kernel
    !
    ! local variables
    REAL(KIND = DP) :: qp, ap, sap, metric, tz, Alp, B, h0, h1, h2, hdr, e(11)
    REAL(KIND = DP) :: eps_inf_env
    !
    eps_inf_env = one
    Alp = half * (one + eps_inf_env / epsi_2d)
    B = one - Alp
    qp = norm2(q(1:2) )
    IF (qp < eps12) qp = eps5
    !
    tz = (tauz - tz_ref) * alat        ! in Bohr
    metric = (twopi / alat)**2       ! in Bohr^-2
    ap = one / (metric * alph * four)  ! in Bohr^2
    sap = SQRT(ap)                   ! in Bohr
    !
    e(1) = EXP(qp * thickn_2d)
    e(2) = one / e(1)             ! exp(-qp*d)
    e(3) = e(2) * e(2)            ! exp(-2qp*d)
    e(4) = EXP(qp * tz)
    e(5) = one / e(4)             ! exp(-qp*tz)
    e(6) = exp(ap * qp * qp)
    e(7) = ERF(sap * qp)
    e(8) = ERF(-tz / two / sap - sap * qp)
    e(9) = ERF(-tz / two / sap + sap * qp)
    e(10) = ERF(thickn_2d / two / sap - sap * qp)
    e(11) = ERF(thickn_2d / two / sap + sap * qp)
    ! tau_z = 0
    h0 = pi * e(6) * (one - e(7))
    ! tau_z != 0
    h1 = half * pi * e(6) *(e(4) * (e(8) + one) + e(5) * (one - e(9)))
    ! tau_z = -d
    h2 = half * pi * e(6) * (e(2) * (e(10) + one) + e(1) * (one - e(11)))
    !
    hdr = half * pi * e(6) * (e(4) * (e(8) + one) - e(5) * (one - e(9)))
    ! c = at(3,3) * alat
    kern = at(3, 3) * alat * (Alp - B) / pi * (h1 + B / (Alp * Alp - B * B * e(3)) * &
           ((Alp * e(4) + B * e(5) * e(3)) * h0 + (B * e(4) * e(2) + Alp * e(5) * e(2) ) * h2))
    dkern = at(3, 3) * alat * (Alp - B) / pi * qp * (hdr + B / (Alp * Alp - B * B * e(3) ) * &
           ((Alp * e(4) - B * e(5) * e(3)) * h0 + (B * e(4) * e(2) - Alp * e(5) * e(2)) * h2))
    !
    RETURN
    !
    !-------------------------------------------------------------------------------
    END SUBROUTINE kernel_2d
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    SUBROUTINE kernel_slab(q, q3, m3, na, kern, dkern)
    !-------------------------------------------------------------------------------
    !! Compute kernel K in eq. (34) and its derivative in the form for slab
    !! using eq. (A1) [PRB 105, 115414 (2022)]. (assume_isolated == 'none')
    !! Coded by W.-H. Sio (2022) then modified and merged by V.-A. Ha (2023)
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : ci, cone, twopi, zero, half, one, two, thre, four, eps5, eps12
    USE global_var,    ONLY : thickn_2d, epsi_2d, tz_ref, emib3tz, eib3d
    USE cell_base,     ONLY : at, alat
    USE ions_base,     ONLY : tau
    !
    INTEGER, INTENT(in) :: m3
    !! FIXME
    INTEGER, INTENT(in) :: na
    !! atom index
    REAL(KIND = DP), INTENT(in)     :: q(3)
    !! Q-point
    REAL(KIND = DP), INTENT(in)     :: q3
    !! FIXME
    COMPLEX(KIND = DP), INTENT(out) :: kern
    !! Output 2D kernel
    COMPLEX(KIND = DP), INTENT(out) :: dkern
    !! Output 2D kernel
    !
    ! Local variables
    REAL(KIND = DP)     :: tz
    !! FIXME
    REAL(KIND = DP)     :: Alp
    !! FIXME
    REAL(KIND = DP)     :: B
    !! FIXME
    REAL(KIND = DP)     :: cd
    !! FIXME
    REAL(KIND = DP)     :: c
    !! FIXME
    REAL(KIND = DP)     :: et
    !! FIXME
    REAL(KIND = DP)     :: QSQ
    !! FIXME
    REAL(KIND = DP)     :: qp
    !! FIXME
    REAL(KIND = DP)     :: qz
    !! FIXME
    REAL(KIND = DP)     :: limit
    !! FIXME
    REAL(KIND = DP)     :: fac
    !! FIXME
    REAL(KIND = DP)     :: e(29)
    !! FIXME
    REAL(KIND = DP)     :: eps_inf_env
    !! FIXME
    REAL(KIND = DP)     :: d
    !! FIXME
    COMPLEX(KIND = DP)  :: emiq3tz
    !! FIXME
    COMPLEX(KIND = DP)  :: eiq3c
    !! FIXME
    COMPLEX(KIND = DP)  :: eiq3d
    !! FIXME
    COMPLEX(KIND = DP)  :: emiq3cd
    !! FIXME
    COMPLEX(KIND = DP)  :: I1
    !! FIXME
    COMPLEX(KIND = DP)  :: I2
    !! FIXME
    COMPLEX(KIND = DP)  :: T(6)
    !! FIXME
    COMPLEX(KIND = DP)  :: eiqzc
    !! FIXME
    COMPLEX(KIND = DP)  :: eiqzd
    !! FIXME
    COMPLEX(KIND = DP)  :: emiqzcd
    !! FIXME
    COMPLEX(KIND = DP)  :: eib3dm3
    !! FIXME
    COMPLEX(KIND = DP)  :: Id1
    !! FIXME
    COMPLEX(KIND = DP)  :: Id2
    !! FIXME
    COMPLEX(KIND = DP)  :: Td(6)
    !! FIXME
    !
    eps_inf_env = one
    d = thickn_2d
    !
    IF (abs(q3) < eps12) THEN
      emiq3tz  = cone
      eiq3c    = cone
      eiq3d    = cone
      emiq3cd  = cone
    ELSE
      emiq3tz = EXP(-ci * twopi * q3 * tau(3, na))
      eiq3c = EXP(ci * twopi * q3 * at(3, 3))
      eiq3d = EXP(ci * twopi * q3 * thickn_2d / alat)
      emiq3cd = eiq3d / eiq3c
    ENDIF
    !
    limit = 30.0_DP
    c = at(3,3) * alat                        ! c in Bohr
    cd = c - d                                ! capital D in Bohr, d in Bohr
    Alp = half *( one + eps_inf_env / epsi_2d )  ! \alpha (Eq. 20)
    B = one - Alp                             ! \beta  (Eq. 21)
    qp = norm2( q(1:2) )                      ! Q paralell in Bohr^-1
    IF (qp < eps12) qp = eps5
    QSQ = DOT_PRODUCT(q,q)                    ! Q^2 in Bohr^-2
    qz = ABS(q(3))                            ! Qz  in Bohr^-1
    tz = (tau(3, na) - tz_ref) * alat         ! tz in Bohr
    !
    et = ACOSH(COSH(qp * (cd - d)) + two * Alp * Alp / (two * Alp - one ) * SINH(qp * cd) * SINH(qp * d)) ! \eta (Eq. 23)
    !
    I1 = (four * Alp - two) * qp / QSQ * emib3tz(na)**m3 * emiq3tz
    Id1 = -ci * qz * I1
    !
    eiqzc = eiq3c
    eib3dm3 = eib3d**m3
    eiqzd = eib3dm3 * eiq3d
    emiqzcd = eib3dm3 * emiq3cd
    !
    IF (qp * cd < limit) THEN
      e(1) = COSH(et)
      e(2) = COSH(qp * tz)
      e(3) = COSH(qp * cd + qp * tz)
      e(4) = SINH(qp * tz)
      e(5) = SINH(qp * cd)
      e(6) = COSH(qp * d + qp * tz)
      e(7) = SINH(qp * cd - qp * d)
      e(8) = SINH(qp * d)
      e(9) = COSH(qp * cd)
      e(10) = SINH(qp * cd + qp * tz)
      e(11) = SINH(qp * d + qp * tz)
      e(12) = COSH(qp * cd - qp * d)
      e(13) = COSH(qp * d)
      !
      ! SP: WARNING - is this correct ? eiqzc is complex and fac is real.
      ! FIXME
      fac = qp / (QSQ * two * (e(1) - eiqzc / two - one / (eiqzc * two)))
      !
      T(1) = fac * four * (Alp - one) * e(2) * eiqzc
      Td(1) = fac * qp * four * (Alp - one) * e(4) * eiqzc
      !
      T(2) = fac * eiqzd * (two*e(3) - four * Alp * e(4) * e(5))
      Td(2) = fac * qp * eiqzd * (two * e(10) - four * Alp * e(2) * e(5))
      !
      T(3) = fac * emiqzcd * (four * Alp - two) * e(6)
      Td(3) = fac * qp * emiqzcd * (four * Alp - two) * e(11)
      !
      fac = fac /( e(7) + two * Alp * e(8) * e(9))
      T(4) = fac * two * (Alp - one) * ((four * Alp * e(4) * e(9) - two * e(10)) * e(9) &
           + two * (two - four * Alp) * e(9) * e(11) * e(1) - (two - four * Alp) * e(12) * e(11))
      Td(4) = fac * qp * two * (Alp - one) * ((four * Alp * e(2) * e(9) - two * e(3)) * e(9) &
            + two * (two - four * Alp) * e(9) * e(6) * e(1) - (two - four * Alp) * e(12) * e(6))
      !
      T(5) = fac*emiqzcd*( two*(two-four*Alp)*e(11)*e(1) &
           + (four*Alp*e(4)*e(9) - two*e(10)) - (two-four*Alp)*(e(12) + two*Alp*e(8)*e(5))*e(11) )
      Td(5) = fac*qp*emiqzcd*( two*(two-four*Alp)*e(6)*e(1) &
            + (four*Alp*e(2)*e(9) - two*e(3)) - (two-four*Alp)*(e(12) + two*Alp*e(8)*e(5))*e(6) )
      !
      T(6) = fac*eiqzd*( (e(12) - two*Alp*e(13)*e(9))*(four*Alp*e(4)*e(9) - two*e(10)) &
           - two*(two*Alp-one)**2 *e(11) + (four*Alp-two)*(four*Alp*e(4)*e(9) - two*e(10))*e(1) )
      Td(6) = fac*qp*eiqzd*( (e(12) - two*Alp*e(13)*e(9))*(four*Alp*e(2)*e(9) - two*e(3)) &
            - two*(two*Alp-one)**2 *e(6) + (four*Alp-two)*(four*Alp*e(2)*e(9) - two*e(3))*e(1) )
      !
    ELSE IF ( qp*cd > limit ) THEN
      e(1) = EXP(qp*tz-et)
      e(2) = EXP(-qp*tz-et)
      e(3) = EXP(qp*cd+qp*tz-et)
      e(4) = EXP(-qp*cd-qp*tz-et)
      e(5) = EXP(qp*cd-qp*tz-et)
      e(6) = EXP(qp*d+qp*tz-et)
      e(7) = EXP(-qp*d-qp*tz-et)
      e(8) = EXP(qp*cd-qp*d+qp*tz-et)
      e(9) = EXP(qp*cd-qp*d-qp*tz-et)
      e(10) = EXP(-qp*cd-qp*d-qp*tz-et)
      e(11) = EXP(qp*tz)
      e(12) = EXP(-two*qp*d-qp*tz)
      e(13) = EXP(-two*qp*d)
      e(14) = EXP(-qp*cd+qp*d+qp*tz-et)
      e(15) = EXP(-qp*cd+qp*tz)
      e(16) = EXP(-qp*cd-two*qp*d-qp*tz)
      e(17) = EXP(-qp*d+qp*tz-et)
      e(18) = EXP(-two*qp*cd-qp*d-qp*tz-et)
      e(19) = EXP(-thre*qp*cd-qp*tz-et)
      e(20) = EXP(-two*qp*cd+qp*d+qp*tz-et)
      e(21) = EXP(qp*cd-two*d+qp*tz-et)
      e(22) = EXP(-qp*cd+qp*tz-et)
      e(23) = EXP(qp*cd-two*qp*d-qp*tz-et)
      e(24) = EXP(-qp*cd-two*qp*d-qp*tz-et)
      e(25) = EXP(-thre*qp*cd-qp*tz-et)
      e(26) = EXP(-qp*d+qp*tz)
      e(27) = EXP(-qp*d-qp*tz)
      e(28) = EXP(-two*qp*cd-qp*d-qp*tz)
      e(29) = EXP(-two*qp*cd)
      !
      fac = qp / QSQ
      !
      T(1) = fac*two*(Alp-one)*(e(1) + e(2))*eiqzc
      Td(1) = fac*two*(Alp-one)*qp*(e(1) - e(2))*eiqzc
      !
      T(2) = fac*eiqzd*( (one-Alp)*e(3) + e(4) + Alp*e(5))
      Td(2) = fac*eiqzd*qp*( (one-Alp)*e(3) - e(4) - Alp*e(5))
      !
      T(3) = fac*emiqzcd*(two*Alp-one)*(e(6) + e(7))
      Td(3) = fac*emiqzcd*qp*(two*Alp-one)*(e(6) - e(7))
      !
      fac = fac/(one - e(29))
      T(4) = fac*(Alp-one)*( two*( (Alp-one)*e(8) - Alp*e(9) + e(10) ) &
           + two*(one-two*Alp)*( e(11) - e(12) - (one+e(13))*(e(14) - e(10)) ) )
      Td(4) = fac*(Alp-one)*qp*( two*( (Alp-one)*e(8) + Alp*e(9) - e(10) ) &
            + two*(one-two*Alp)*( e(11) + e(12) - (one+e(13))*(e(14) + e(10)) ) )
      !
      T(5) = fac*emiqzcd*( two*(one-two*Alp)*(e(15) - e(16)) + two*(Alp-one)*e(17) &
           - two*Alp*e(7) + two*e(18) - (one-two*Alp)*( (one-Alp)*e(17) &
           + (Alp-one)*e(19) + Alp*e(6) - Alp*e(7) + e(20) - e(18) ) )
      Td(5) = fac*emiqzcd*qp*( two*(one-two*Alp)*(e(15) + e(16)) + two*(Alp-one)*e(17) &
            + two*Alp*e(7) - two*e(18) - (one-two*Alp)*( (one-Alp)*e(17) &
            - (Alp-one)*e(19) + Alp*e(6) + Alp*e(7) + e(20) + e(18) ) )
      !
      T(6) = fac*eiqzd*( (-(one-Alp)**2*e(21) - Alp*(Alp-one)*e(3) + (Alp-one)*e(22) &
           - Alp*(one-Alp)*e(23) + Alp**2*e(5) - Alp*e(4) + (one-Alp)*e(24) - Alp*e(4) + e(25) ) &
           - two*(two*Alp-one)**2*(e(22) - e(24)) + two*(two*Alp-one)*( &
           (Alp-one)*e(26) - Alp*e(27) + e(28) ) )
      Td(6) = fac*eiqzd*qp*( (-(one-Alp)**2*e(21) - Alp*(Alp-one)*e(3) + (Alp-one)*e(22) &
            + Alp*(one-Alp)*e(23) - Alp**2*e(5) + Alp*e(4) - (one-Alp)*e(24) + Alp*e(4) - e(25) ) &
            - two*(two*Alp-one)**2*(e(22) + e(24)) + two*(two*Alp-one)*( &
            (Alp-one)*e(26) + Alp*e(27) - e(28) ) )
      !
    ENDIF
    !
    I2 = SUM(T)
    Id2 = SUM(Td)
    kern = I1 + I2
    dkern = Id1 + Id2
    !
    RETURN
    !-------------------------------------------------------------------------------
    END SUBROUTINE kernel_slab
    !-------------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_umn_f(nbnd, cufkk, cufkq, bmatf)
    !-----------------------------------------------------------------------
    !!
    !! Calculates $$ U_{k+q} U_k^\dagger = <\Psi_{mk+q}|e^{i{q+G}r}|\Psi_{nk}> $$
    !! in the approximation q+G->0 on the fine grids.
    !!
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero, cone
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
    !  U(k'+q') * U(k')^\dagger
    !
    CALL ZGEMM('n', 'c', nbnd, nbnd, nbnd, cone, cufkq, nbnd, cufkk, nbnd, czero, bmatf, nbnd)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE compute_umn_f
    !-----------------------------------------------------------------------
  !---------------------------------------------------------------------------------
  END MODULE longrange
  !---------------------------------------------------------------------------------
