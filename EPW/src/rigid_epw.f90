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
  IMPLICIT NONE
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
    !!
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : fpi, e2
    USE constants,     ONLY : pi
    USE cell_base,     ONLY : bg, omega
    USE constants_epw, ONLY : eps6, czero
    USE epwcom,        ONLY : system_2d
    USE io_global,     ONLY : ionode_id
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : world_comm
    USE division,      ONLY : para_bounds
    USE mp,            ONLY : mp_bcast, mp_sum
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
    REAL(KIND = DP):: geg
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: alph
    !! Ewald parameter
    REAL(KIND = DP) :: fac
    !! Missing definition
    REAL(KIND = DP) :: g1, g2, g3
    !! Missing definition
    REAL(KIND = DP) :: facgd
    !! fac * EXP(-geg / (alph * 4.0d0)) / geg
    REAL(KIND = DP) :: arg
    !! Missing definition
    REAL(KIND = DP) :: gmax
    !! Maximum G
    REAL(KIND = DP) :: zag(3)
    !! Z * G
    REAL(KIND = DP) :: zbg(3)
    !! Z * G
    REAL(KIND = DP) :: zcg(3)
    !! Z * G
    REAL(KIND = DP) :: fnat(3)
    !! Missing definition
    COMPLEX(KIND = DP) :: facg
    !! Missing definition
    COMPLEX(KIND = DP) :: dyn_tmp(3 * nat, 3 * nat)
    !! Temporary dyn. matrice
    !
    ! alph is the Ewald parameter, geg is an estimate of G^2
    ! such that the G-space sum is convergent for that alph
    ! very rough estimate: geg/4/alph > gmax = 14
    ! (exp (-14) = 10^-6)
    !
    IF (ABS(ABS(signe) - 1.0) > eps6) CALL errore('rgd_blk', ' wrong value for signe ', 1)
    ! 
    IF (system_2d) THEN
      !fac = (signe * e2 * fpi) / omega*0.5d0*alat/bg(3,3)
      !reff=0.0d0
      !DO i=1,2
      !   DO j=1,2
      !      reff(i,j)=epsil(i,j)*0.5d0*tpi/bg(3,3) ! (eps)*c/2 in 2pi/a units
      !   ENDDO
      !ENDDO
      !DO i=1,2
      !   reff(i,i)=reff(i,i)-0.5d0*tpi/bg(3,3) ! (-1)*c/2 in 2pi/a units
      !ENDDO
    ELSE
      ! (e^2 * 4\pi) / Volume
      fac = (signe * e2 * fpi) / omega
    ENDIF
    !
    gmax = 14.d0
    alph = 1.0d0
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
      g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1,3)
      g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2,3)
      g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3,3)
      !
      geg = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
             g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
             g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3))
      !
      IF (geg > 0.0d0 .AND. geg / (alph * 4.0d0) < gmax) THEN
        !
        facgd = fac * EXP(-geg / (alph * 4.0d0)) / geg
        !
        DO na = 1, nat
          zag(:) = g1 * zeu(1, :, na) + g2 * zeu(2, :, na) + g3 * zeu(3, :, na)
          fnat(:) = 0.d0
          DO nb = 1, nat
            arg = 2.d0 * pi * (g1 * (tau(1, na) - tau(1, nb)) + &
                               g2 * (tau(2, na) - tau(2, nb)) + &
                               g3 * (tau(3, na) - tau(3, nb)))
            zcg(:)  = g1 * zeu(1, :, nb) + g2 * zeu(2, :, nb) + g3 * zeu(3, :, nb)
            fnat(:) = fnat(:) + zcg(:) * COS(arg)
          ENDDO
          DO j = 1, 3
            DO i = 1, 3
              dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (na - 1) * 3 + j) &
                                           - facgd * zag(i) * fnat(j)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! nat
      ENDIF ! geg
      !
      g1 = g1 + q(1)
      g2 = g2 + q(2)
      g3 = g3 + q(3)
      !
      geg = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
             g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
             g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3))
      !
      IF (geg > 0.0d0 .AND. geg / (alph * 4.0d0) < gmax) THEN
        !
        facgd = fac * EXP(-geg / (alph * 4.0d0)) / geg
        !
        DO nb = 1, nat
          zbg(:) = g1 * zeu(1, :, nb) + g2 * zeu(2, :, nb) + g3 * zeu(3, :, nb)
          DO na = 1, nat
            zag(:) = g1 * zeu(1, :, na) + g2 * zeu(2, :, na) + g3 * zeu(3, :, na)
            arg = 2.d0 * pi * (g1 * (tau(1, na) - tau(1 ,nb)) + &
                            g2 * (tau(2, na) - tau(2, nb)) + &
                            g3 * (tau(3, na) - tau(3, nb)) )
            !
            facg = facgd * CMPLX(COS(arg), SIN(arg), DP)
            DO j = 1, 3
              DO i = 1, 3
                dyn_tmp((na - 1) * 3 + i, (nb - 1) * 3 + j) = dyn_tmp((na - 1) * 3 + i, (nb - 1) * 3 + j) &
                                             + facg * zag(i) * zbg(j)
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
    USE cell_base,     ONLY : bg, omega, alat
    USE ions_base,     ONLY : tau, nat
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12
    USE epwcom,        ONLY : shortrange
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
    REAL(KIND = DP), INTENT(in) :: zeu(3, 3, nat)
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
    INTEGER :: ipol
    !! Polarison
    INTEGER :: m1, m2, m3
    !! Loop over q-points
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
    COMPLEX(KIND = DP) :: fac
    !!
    COMPLEX(KIND = DP) :: facqd
    !!
    COMPLEX(KIND = DP) :: facq
    !!
    COMPLEX(KIND = DP) :: epmatl(nmodes)
    !! Long-range part of the el-ph matrix elements
    !
    IF(ABS(ABS(signe) - 1.0) > eps12) CALL errore('rgd_blk_epw', 'Erong value for signe ', 1)
    !
    gmax = 14.d0
    alph = 1.0d0
    geg = gmax * alph * 4.0d0
    fac = signe * e2 * fpi / omega * ci
    !
    !
    epmatl(:) = czero
    !DO m1 = -nrx1, nrx1
    DO m1 = -nqc1, nqc1
      DO m2 = -nqc2, nqc2
        DO m3 = -nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)
          !
          qeq = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
                 g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
                 g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3)) !*twopi/alat
          !
          IF (qeq > 0.0d0 .AND. qeq / (alph * 4.0d0) < gmax) THEN
            !
            qeq = qeq * twopi / alat
            facqd = fac * EXP(-qeq / (alph * 4.0d0)) / qeq !/(two*wq)
            !
            DO na = 1, nat
              arg = - twopi * (g1 * tau(1, na) + g2 * tau(2, na) + g3 * tau(3, na))
              facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
              DO ipol = 1, 3
                zaq = g1 * zeu(1, ipol, na) + g2 * zeu(2, ipol, na) + g3 * zeu(3, ipol, na)
                !
                epmat = epmat + facq * zaq * uq(3 * (na - 1) + ipol, :) * bmat
                epmatl = epmatl + facq * zaq * uq(3 * (na - 1) + ipol, :) * bmat
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
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12
    USE epwcom,        ONLY : shortrange, nbndsub
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
    REAL (KIND = DP), INTENT(in) :: zeu(3, 3, nat)
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
    COMPLEX(KIND = DP) :: fac
    !!
    COMPLEX(KIND = DP) :: facqd
    !!
    COMPLEX(KIND = DP) :: facq
    !!
    COMPLEX(KIND = DP) :: epmatl(nbndsub, nbndsub, nmodes)
    !! Long-range part of the matrix element
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore ('rgd_blk_epw_fine', 'Wrong value for signe ', 1)
    !
    gmax = 14.d0
    alph = 1.0d0
    geg = gmax * alph * 4.0d0
    fac = signe * e2 * fpi / omega * ci
    !
    epmatl(:, :, :) = czero
    !
    DO m1 = -nqc1, nqc1
      DO m2 = -nqc2, nqc2
        DO m3 = -nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)
          !
          qeq = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
                 g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
                 g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3)) !*twopi/alat
          !
          IF (qeq > 0.0d0 .AND. qeq / (alph * 4.0d0) < gmax) THEN
            !
            qeq = qeq * twopi / alat
            facqd = fac * EXP(-qeq / (alph * 4.0d0)) / qeq !/(two*wq)
            !
            DO na = 1, nat
              arg = -twopi * (g1 * tau(1, na) + g2 * tau(2, na) + g3 * tau(3, na))
              facq = facqd * CMPLX(COS(arg), SIN(arg), kind=DP)
              DO ipol = 1, 3
                zaq = g1 * zeu(1, ipol, na) + g2 * zeu(2, ipol, na) + g3 * zeu(3, ipol, na)
                !
                DO imode = 1, nmodes
                  CALL ZAXPY(nbndsub**2, facq * zaq * uq(3 * (na - 1) + ipol, imode), bmat(:, :), 1, epmat(:, :, imode), 1)
                  CALL ZAXPY(nbndsub**2, facq * zaq * uq(3 * (na - 1) + ipol, imode), bmat(:, :), 1, epmatl(:, :, imode), 1)
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
    USE epwcom,        ONLY : meff, fermi_diff, nel, smear_rpa
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
    n = nel / omega
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
    !! compute the rigid-ion (long-range) derivative term for q
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
    INTEGER :: m1, m2, m3
    !! Loop over q-points
    INTEGER :: isum
    !! Index to sum the different component of the derivative
    REAL(KIND = DP):: geg
    !! <q+G| epsil | q+G>
    REAL(KIND = DP) :: alph
    !! Missing definition
    REAL(KIND = DP) :: fac
    !! Missing definition
    REAL(KIND = DP) :: g1, g2, g3
    !! Missing definition
    REAL(KIND = DP) :: facgd
    !! Missing definition
    REAL(KIND = DP) :: arg
    !! Missing definition
    REAL(KIND = DP) :: gmax
    !! Missing definition
    REAL(KIND = DP) :: arg_no_g(3)
    !! Missing definition
    REAL(KIND = DP) :: zag(3)
    !! Missing definition
    REAL(KIND = DP) :: zbg(3)
    !! Missing definition
    REAL(KIND = DP) :: zbg_der(3, 3)
    !! Missing definition
    REAL(KIND = DP) :: zag_der(3, 3)
    !! Missing definition
    COMPLEX(KIND = DP) :: facg
    !! Missing definition
    !
    ! alph is the Ewald parameter, geg is an estimate of G^2
    ! such that the G-space sum is convergent for that alph
    ! very rough estimate: geg/4/alph > gmax = 14
    ! (exp (-14) = 10^-6)
    !
    gmax = 14.d0
    alph = 1.0d0
    geg = gmax * alph * 4.0d0
    !
    IF (ABS(ABS(signe) - 1.0) > eps6) CALL errore('rgd_blk_der', ' wrong value for signe ', 1)
    !
    gmax = 14.d0
    alph = 1.0d0
    geg = gmax * alph * 4.0d0
    fac = signe * e2 * fpi / omega
    !
    DO m1 = -nqc1, nqc1
      DO m2 = -nqc2, nqc2
        DO m3 = -nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3)
          !
          g1 = g1 + q(1)
          g2 = g2 + q(2)
          g3 = g3 + q(3)
          !
          geg = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3) + &
                 g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3) + &
                 g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3))
          !
          IF (geg > 0.0d0 .AND. geg / (alph * 4.0d0) < gmax) THEN
            !
            facgd = fac * EXP(-geg / (alph * 4.0d0) ) / geg * (alat / twopi)
            !
            DO nb = 1, nat
              zbg(:) = g1 * zeu(1, :, nb) + g2 * zeu(2, :, nb) + g3 * zeu(3, :, nb)
              zbg_der(:, :) = zeu(:, :, nb)
              DO na = 1, nat
                zag(:) = g1 * zeu(1, :, na) + g2 * zeu(2, :, na) + g3 * zeu(3, :, na)
                zag_der(:, :) = zeu(:, :, na)
                arg = 2.d0 * pi * (g1 * (tau(1, na) - tau(1, nb)) + &
                                   g2 * (tau(2, na) - tau(2, nb)) + &
                                   g3 * (tau(3, na) - tau(3, nb)))
                arg_no_g(:) = 2.d0 * pi * (tau(:,na) - tau(:,nb))
                !
                facg = facgd * CMPLX(COS(arg), SIN(arg), DP)
                DO j = 1, 3
                  DO i = 1, 3
                    dyn_der_part(1, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) = facg * zag_der(:, i) * zbg(j)
                    dyn_der_part(2, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) = facg * zag(i) * zbg_der(:, j)
                    dyn_der_part(3, :, (na - 1) * 3 + i,( nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(:, 1) * g1 + epsil(:, 2) * g2 + epsil(:, 3) * g3) / geg
                    dyn_der_part(4, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(1, :) * g1 + epsil(2, :) * g2 + epsil(3, :) * g3) / geg
                    dyn_der_part(5, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) = facg * zag(i) * zbg(j) * ci * arg_no_g(:)
                    dyn_der_part(6, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(1, :) * g1 + epsil(2, :) * g2 + epsil(3, :) * g3) / (4d0 * alph)
                    dyn_der_part(7, :, (na - 1) * 3 + i, (nb - 1) * 3 + j) =-facg * zag(i) * zbg(j) &
                      * (epsil(:, 1) * g1 + epsil(:, 2) * g2 + epsil(:, 3) * g3) / (4d0 * alph)
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
    USE constants_epw, ONLY : twopi, fpi, e2, ci, czero, eps12
    USE epwcom,        ONLY : shortrange, nbndsub
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
    REAL(KIND = DP), INTENT(in) :: zeu(3, 3, nat)
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
    INTEGER :: ipol
    !! Polarison
    INTEGER :: m1, m2, m3
    !! Loop over q-points
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
    COMPLEX(KIND = DP) :: fac
    !!
    COMPLEX(KIND = DP) :: facqd
    !!
    COMPLEX(KIND = DP) :: facq
    !!
    COMPLEX(KIND = DP) :: epmatl(nbndsub, nbndsub)
    !! Long-range part of the matrix element
    !
    IF (ABS(ABS(signe) - 1.0) > eps12) CALL errore('rgd_blk_epw_fine_mem', ' wrong value for signe ', 1)
    !
    gmax = 14.d0
    alph = 1.0d0
    geg  = gmax * alph * 4.0d0
    fac  = signe * e2 * fpi / omega * ci
    !
    epmatl(:, :) = czero
    !
    DO m1 = -nqc1, nqc1
      DO m2 = -nqc2, nqc2
        DO m3 = -nqc3, nqc3
          !
          g1 = m1 * bg(1, 1) + m2 * bg(1, 2) + m3 * bg(1, 3) + q(1)
          g2 = m1 * bg(2, 1) + m2 * bg(2, 2) + m3 * bg(2, 3) + q(2)
          g3 = m1 * bg(3, 1) + m2 * bg(3, 2) + m3 * bg(3, 3) + q(3)
          !
          qeq = (g1 * (epsil(1, 1) * g1 + epsil(1, 2) * g2 + epsil(1, 3) * g3 ) + &
                 g2 * (epsil(2, 1) * g1 + epsil(2, 2) * g2 + epsil(2, 3) * g3 ) + &
                 g3 * (epsil(3, 1) * g1 + epsil(3, 2) * g2 + epsil(3, 3) * g3 )) !*twopi/alat
          !
          IF (qeq > 0.0_DP .AND. qeq / alph / 4.0_DP < gmax) THEN
            !
            qeq = qeq * twopi / alat
            facqd = fac * EXP(-qeq / alph / 4.0d0) / qeq !/(two*wq)
            !
            DO na = 1, nat
              arg = -twopi * (g1 * tau(1, na) + g2 * tau(2, na) + g3 * tau(3, na))
              facq = facqd * CMPLX(COS(arg), SIN(arg), KIND = DP)
              DO ipol = 1, 3
                zaq = g1 * zeu(1, ipol, na) + g2 * zeu(2, ipol, na) + g3 * zeu(3, ipol, na)
                !
                CALL ZAXPY(nbndsub**2, facq * zaq * uq(3 * (na - 1) + ipol, imode), bmat(:, :), 1, epmat(:, :), 1)
                CALL ZAXPY(nbndsub**2, facq * zaq * uq(3 * (na - 1) + ipol, imode), bmat(:, :), 1, epmatl(:, :), 1)
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
