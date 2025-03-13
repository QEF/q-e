MODULE esm_stres_mod

  USE kinds,    ONLY : DP
  USE esm_common_mod, ONLY : esm_w, esm_bc, &
                             mill_2d, imill_2d, ngm_2d, exp_erfc
  IMPLICIT NONE

  ! Workaround for Cray bug - note that exp, cosh, sinh with complex argument
  ! are in the F2008 standard so qe_exp, qe_cosh, qe_sinh are no longer needed

#if defined(__CRAY)
#define QE_EXP  exp
#define QE_COSH cosh
#define QE_SINH sinh
#endif  
CONTAINS

  !-----------------------------------------------------------------------
  !--------------ESM STRESS SUBROUTINE------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE esm_stres_har(sigmahar, rhog)
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    IMPLICIT NONE
    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)

    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_har must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_har_bc1(sigmahar, rhog)
    CASE ('bc2')
      CALL esm_stres_har_bc2(sigmahar, rhog)
    CASE ('bc3')
      CALL esm_stres_har_bc3(sigmahar, rhog)
    CASE ('bc4')
      STOP 'esm_stres_har has not yet implemented for esm_bc = bc4'
    END SELECT

    RETURN
  END SUBROUTINE esm_stres_har

  SUBROUTINE esm_stres_ewa(sigmaewa)
    !-----------------------------------------------------------------------
    !
    ! Calculates Ewald stresswith both G- and R-space terms.
    ! Determines optimal alpha. Should hopefully work for any structure.
    !
    USE kinds,     ONLY : DP
    USE constants, ONLY : tpi
    USE cell_base, ONLY : tpiba2
    USE ions_base, ONLY : zv, nat, ityp
    USE gvect,     ONLY : gcutm
    IMPLICIT NONE
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    ! output: the ewald stress
    !
    !    here the local variables
    !
    INTEGER :: ia
    ! counter on atoms

    REAL(DP) :: charge, alpha, upperbound
    ! total ionic charge in the cell
    ! alpha term in ewald sum
    ! the maximum radius to consider real space sum
    REAL(DP) :: sigmaewg(3, 3), sigmaewr(3, 3)
    ! ewald stress computed in reciprocal space
    ! ewald stress computed in real space

    charge = sum(zv(ityp(:)))

    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    alpha = 2.9d0
    DO
      alpha = alpha - 0.1d0
      IF (alpha .le. 0.d0) CALL errore('esm_stres_ewa', 'optimal alpha not found', 1)
      upperbound = 2.d0*charge**2*sqrt(2.d0*alpha/tpi)* &
                   erfc(sqrt(tpiba2*gcutm/4.d0/alpha))
      IF (upperbound < 1.0d-7) EXIT
    END DO

    ! G-space sum here.
    ! Determine if this processor contains G=0 and set the constant term
    CALL esm_stres_ewg(alpha, sigmaewg)

    ! R-space sum here (only for the processor that contains G=0)
    CALL esm_stres_ewr(alpha, sigmaewr)

    sigmaewa(:, :) = sigmaewg(:, :) + sigmaewr(:, :)

    RETURN
  END SUBROUTINE esm_stres_ewa

  SUBROUTINE esm_stres_loclong(sigmaloclong, rhog)
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    IMPLICIT NONE
    REAL(DP), INTENT(out)   :: sigmaloclong(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)

    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_loclong must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_loclong_bc1(sigmaloclong, rhog)
    CASE ('bc2')
      CALL esm_stres_loclong_bc2(sigmaloclong, rhog)
    CASE ('bc3')
      CALL esm_stres_loclong_bc3(sigmaloclong, rhog)
    CASE ('bc4')
      STOP 'esm_stres_loclong has not yet implemented for esm_bc = bc4'
    END SELECT

    RETURN
  END SUBROUTINE esm_stres_loclong

  SUBROUTINE esm_stres_har_bc1(sigmahar, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)

    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, jz
    REAL(DP) :: L, S, z0, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum1c, sum2c, sum2p, sum2m
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dgp2_deps(2, 2)  !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps

    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0

    ! initialize
    sigmahar(:, :) = 0.0d0

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))

    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)

    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3

      ! expand function symmetrically to gz<0
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      END IF
    END DO ! ig

    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))

      IF (gp == 0.0d0) CYCLE ! skip gp=0

      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dgp2_deps(la, mu) = -g(la)*g(mu)*2.0d0
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO

      ! summations over gz
      sum1p = (0.d0, 0.d0)
      sum1m = (0.d0, 0.d0)
      sum2p = (0.d0, 0.d0)
      sum2m = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, igp)
        sum1p = sum1p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)
        sum1m = sum1m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)
        sum2p = sum2p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)**2
        sum2m = sum2m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)**2
      END DO ! igz

      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        dVr_deps(iz, :, :) = &
          -(dgp_deps(:, :)*tpi/gp**2*(gp*(z - z0) - 1.0d0) &
          - delta(:, :)*tpi/gp) &
          *EXP(+gp*(z - z0))*sum1p &
          + dgp_deps(:, :)*tpi/gp*EXP(+gp*(z - z0))*sum2p &
          + (dgp_deps(:, :)*tpi/gp**2*(gp*(z + z0) + 1.0d0) &
          + delta(:, :)*tpi/gp) &
          *EXP(-gp*(z + z0))*sum1m &
          + dgp_deps(:, :)*tpi/gp*EXP(-gp*(z + z0))*sum2m
      END DO ! iz

      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO

      ! add bare coulomn terms to dV(gz)/deps
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, igp)

        dVg_deps(iz, :, :) = dVg_deps(iz, :, :) &
                              - delta(:, :)*fpi*rg3/(gp**2 + gz**2) &
                              - dgp2_deps(:, :)*fpi*rg3/(gp**2 + gz**2)**2
      END DO ! igz

      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, :, :))
      END DO ! igz
    END DO ! igp

    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! summations over gz
      sum1c = (0.d0, 0.d0)
      sum2c = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, imill_2d(0, 0))
        sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
        sum2c = sum2c + rg3*cos(gz*z0)/gz**2
      END DO ! igz

      ! calculate V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        rg3 = rhog3(1, imill_2d(0, 0))
        Vr(iz) = &
          -tpi*z**2*rg3 &
          - tpi*z0**2*rg3 &
          - fpi*z*sum1c &
          - fpi*sum2c
      END DO ! iz

      ! separation by polynomial
      z_l = -z0
      z_r = +z0
      f1 = -tpi*z_r**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_r*sum1c &
           - fpi*sum2c
      f2 = -tpi*z_l**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_l*sum1c &
           - fpi*sum2c
      f3 = -fpi*z_r*rg3 &
           - fpi*sum1c
      f4 = -fpi*z_l*rg3 &
           - fpi*sum1c
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3

      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)

      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0

      ! add bare coulomn terms to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, imill_2d(0, 0))
        Vg(iz) = Vg(iz) + fpi*rg3/gz**2
      END DO ! igz

      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz

    END IF ! imill_2d(0,0) > 0

    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:, :) = sigmahar(:, :)*(-0.5d0*e2)

    CALL mp_sum(sigmahar, intra_bgrp_comm)

    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)

    RETURN
  END SUBROUTINE esm_stres_har_bc1

  SUBROUTINE esm_stres_har_bc2(sigmahar, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)

    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, jz
    REAL(DP) :: L, S, z0, z1, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum1c, sum2c, sum2p, sum2m
    COMPLEX(DP) :: sum1sp, sum1sm, sum1cp, sum1cm, sum2sp, sum2sm
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dgp2_deps(2, 2)  !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps

    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w

    ! initialize
    sigmahar(:, :) = 0.0d0

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))

    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)

    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3

      ! expand function symmetrically to gz<0
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      END IF
    END DO ! ig

    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))

      IF (gp == 0.0d0) CYCLE ! skip gp=0

      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dgp2_deps(la, mu) = -g(la)*g(mu)*2.0d0
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO

      ! summations over gz
      sum1p = (0.d0, 0.d0)
      sum1m = (0.d0, 0.d0)
      sum2p = (0.d0, 0.d0)
      sum2m = (0.d0, 0.d0)
      sum1sp = (0.d0, 0.d0)
      sum1sm = (0.d0, 0.d0)
      sum1cp = (0.d0, 0.d0)
      sum1cm = (0.d0, 0.d0)
      sum2sp = (0.d0, 0.d0)
      sum2sm = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, igp)

        sum1p = sum1p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)
        sum1m = sum1m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)
        sum2p = sum2p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)**2
        sum2m = sum2m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)**2

        sum1sp = sum1sp + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)
        sum1sm = sum1sm + rg3*QE_SINH(gp*z0 - ci*gz*z0)/(gp - ci*gz)

        sum1cp = sum1cp + rg3*QE_COSH(gp*z0 + ci*gz*z0)/(gp + ci*gz)*z0
        sum1cm = sum1cm + rg3*QE_COSH(gp*z0 - ci*gz*z0)/(gp - ci*gz)*z0

        sum2sp = sum2sp + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)**2
        sum2sm = sum2sm + rg3*QE_SINH(gp*z0 - ci*gz*z0)/(gp - ci*gz)**2
      END DO ! igz

      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        !! BC1 terms
        dVr_deps(iz, :, :) = &
          -(dgp_deps(:, :)*tpi/gp**2*(gp*(z - z0) - 1.0d0) &
          - delta(:, :)*tpi/gp) &
          *EXP(+gp*(z - z0))*sum1p &
          + dgp_deps(:, :)*tpi/gp*EXP(+gp*(z - z0))*sum2p &
          + (dgp_deps(:, :)*tpi/gp**2*(gp*(z + z0) + 1.0d0) &
          + delta(:, :)*tpi/gp) &
          *EXP(-gp*(z + z0))*sum1m &
          + dgp_deps(:, :)*tpi/gp*EXP(-gp*(z + z0))*sum2m

        !! BC2 terms
        dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                             + dgp_deps(:, :)*( &
                             - tpi/gp**2*(EXP(-gp*(z + 2*z1)) - EXP(+gp*z))/sinh(2*gp*z1) &
                             + tpi/gp*((-z - 2*z1)*EXP(-gp*(z + 2*z1)) - z*EXP(+gp*z))/sinh(2*gp*z1) &
                             - tpi/gp*(EXP(-gp*(z + 2*z1)) - EXP(+gp*z))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)) &
                             * sum1sp &
                             + tpi/gp*(EXP(-gp*(z + 2*z1)) - EXP(+gp*z))/sinh(2*gp*z1) &
                             * (-delta(:, :)*sum1sp + dgp_deps(:, :)*(sum1cp - sum2sp)) &
                             + dgp_deps(:, :)*( &
                             - tpi/gp**2*(EXP(+gp*(z - 2*z1)) - EXP(-gp*z))/sinh(2*gp*z1) &
                             + tpi/gp*((+z - 2*z1)*EXP(+gp*(z - 2*z1)) + z*EXP(-gp*z))/sinh(2*gp*z1) &
                             - tpi/gp*(EXP(+gp*(z - 2*z1)) - EXP(-gp*z))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)) &
                             * sum1sm &
                             + tpi/gp*(EXP(+gp*(z - 2*z1)) - EXP(-gp*z))/sinh(2*gp*z1) &
                             * (-delta(:, :)*sum1sm + dgp_deps(:, :)*(sum1cm - sum2sm))

      END DO ! iz

      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO

      ! add bare couloum terms to dV(gz)/deps
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, igp)

        dVg_deps(iz, :, :) = dVg_deps(iz, :, :) &
                              - delta(:, :)*fpi*rg3/(gp**2 + gz**2) &
                              - dgp2_deps(:, :)*fpi*rg3/(gp**2 + gz**2)**2
      END DO ! igz

      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, :, :))
      END DO ! igz
    END DO ! igp

    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! summations over gz
      sum1c = (0.d0, 0.d0)
      sum2c = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, imill_2d(0, 0))
        sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
        sum2c = sum2c + rg3*cos(gz*z0)/gz**2
      END DO ! igz

      ! calculate V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        rg3 = rhog3(1, imill_2d(0, 0))

        !! BC1 terms
        Vr(iz) = &
          - tpi*z**2*rg3 &
          - tpi*z0**2*rg3 &
          - fpi*z*sum1c &
          - fpi*sum2c

        !! BC2 terms
        Vr(iz) = Vr(iz) &
                 + tpi*z1*2*z0*rg3 - tpi*(-z/z1)*2*z0*sum1c
      END DO ! iz

      ! separation by polynomial
      z_l = -z0
      z_r = +z0
      f1 = -tpi*z_r**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_r*sum1c &
           - fpi*sum2c
      f1 = f1 &
           + tpi*z1*2*z0*rg3 - tpi*(-z_r/z1)*2*z0*sum1c

      f2 = -tpi*z_l**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_l*sum1c &
           - fpi*sum2c
      f2 = f2 &
           + tpi*z1*2*z0*rg3 - tpi*(-z_l/z1)*2*z0*sum1c

      f3 = -fpi*z_r*rg3 &
           - fpi*sum1c
      f3 = f3 &
           - tpi*(-1.0d0/z1)*2*z0*sum1c

      f4 = -fpi*z_l*rg3 &
           - fpi*sum1c
      f4 = f4 &
           - tpi*(-1.0d0/z1)*2*z0*sum1c

      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3

      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)

      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0

      ! add bare coulomn terms to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, imill_2d(0, 0))
        Vg(iz) = Vg(iz) + fpi*rg3/gz**2
      END DO ! igz

      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END IF ! imill_2d(0,0) > 0

    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:, :) = sigmahar(:, :)*(-0.5d0*e2)

    CALL mp_sum(sigmahar, intra_bgrp_comm)

    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)

    RETURN
  END SUBROUTINE esm_stres_har_bc2

  SUBROUTINE esm_stres_har_bc3(sigmahar, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)

    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, jz
    REAL(DP) :: L, S, z0, z1, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum2p, sum2m, sum1c, sum2c
    COMPLEX(DP) :: sum1sh, sum1ch, sum2sh
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dgp2_deps(2, 2)  !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps

    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)

    REAL(DP) :: sigmahar_bc1(3, 3)

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w

    ! initialize
    sigmahar(:, :) = 0.0d0

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))

    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)

    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3

      ! expand function symmetrically to gz<0
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      END IF
    END DO ! ig

    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))

      IF (gp == 0.0d0) CYCLE ! skip gp=0

      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dgp2_deps(la, mu) = -g(la)*g(mu)*2.0d0
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO

      ! summations over gz
      sum1p = (0.d0, 0.d0)
      sum1m = (0.d0, 0.d0)
      sum2p = (0.d0, 0.d0)
      sum2m = (0.d0, 0.d0)
      sum1sh = (0.d0, 0.d0)
      sum1ch = (0.d0, 0.d0)
      sum2sh = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, igp)
        sum1p = sum1p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)
        sum1m = sum1m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)
        sum2p = sum2p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)**2
        sum2m = sum2m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)**2
        sum1sh = sum1sh + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)
        sum1ch = sum1ch + rg3*QE_COSH(gp*z0 + ci*gz*z0)/(gp + ci*gz)*z0
        sum2sh = sum2sh + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)**2
      END DO ! igz

      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        !! BC1 terms
        dVr_deps(iz, :, :) = &
          -(dgp_deps(:, :)*tpi/gp**2*(gp*(z - z0) - 1.0d0) &
            - delta(:, :)*tpi/gp) &
          *EXP(+gp*(z - z0))*sum1p &
          + dgp_deps(:, :)*tpi/gp*EXP(+gp*(z - z0))*sum2p &
          + (dgp_deps(:, :)*tpi/gp**2*(gp*(z + z0) + 1.0d0) &
             + delta(:, :)*tpi/gp) &
          *EXP(-gp*(z + z0))*sum1m &
          + dgp_deps(:, :)*tpi/gp*EXP(-gp*(z + z0))*sum2m

        !! BC3 termn
        dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                             - dgp_deps(:, :)*( &
                             -fpi/gp**2*EXP(-gp*(-z + 2*z1)) &
                             - fpi/gp*(-z + 2*z1)*EXP(-gp*(-z + 2*z1)) &
                             )*sum1sh &
                             - fpi/gp*EXP(-gp*(-z + 2*z1))*( &
                             -delta(:, :)*sum1sh &
                             + dgp_deps(:, :)*(sum1ch - sum2sh))
      END DO ! iz

      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO

      ! add bare coulomn terms to dV(gz)/deps
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, igp)

        dVg_deps(iz, :, :) = dVg_deps(iz, :, :) &
                              - delta(:, :)*fpi*rg3/(gp**2 + gz**2) &
                              - dgp2_deps(:, :)*fpi*rg3/(gp**2 + gz**2)**2
      END DO ! igz

      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, :, :))
      END DO ! igz
    END DO ! igp

    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! summations over gz
      sum1c = (0.d0, 0.d0)
      sum2c = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, imill_2d(0, 0))
        sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
        sum2c = sum2c + rg3*cos(gz*z0)/gz**2
      END DO ! igz

      ! calculate V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        rg3 = rhog3(1, imill_2d(0, 0))
        !! BC1 terms
        Vr(iz) = &
          - tpi*z**2*rg3 &
          - tpi*z0**2*rg3 &
          - fpi*z*sum1c &
          - fpi*sum2c

        !! BC3 terms
        Vr(iz) = Vr(iz) - tpi*(z - 2*z1)*2*z0*rg3 + fpi*z0*sum1c
      END DO ! iz

      ! separation by polynomial
      z_l = -z0
      z_r = +z0
      f1 = -tpi*z_r**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_r*sum1c &
           - fpi*sum2c
      f1 = f1 &
           - tpi*(z_r - 2*z1)*2*z0*rg3 + fpi*z0*sum1c

      f2 = -tpi*z_l**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_l*sum1c &
           - fpi*sum2c
      f2 = f2 &
           - tpi*(z_l - 2*z1)*2*z0*rg3 + fpi*z0*sum1c

      f3 = -fpi*z_r*rg3 &
           - fpi*sum1c
      f3 = f3 &
           - tpi*(1.0d0)*2*z0*rg3

      f4 = -fpi*z_l*rg3 &
           - fpi*sum1c
      f4 = f4 &
           - tpi*(1.0d0)*2*z0*rg3

      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3

      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)

      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0

      ! add bare coulomn terms to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        rg3 = rhog3(iz, imill_2d(0, 0))
        Vg(iz) = Vg(iz) + fpi*rg3/gz**2
      END DO ! igz

      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))

        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END IF ! imill_2d(0,0) > 0

    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:, :) = sigmahar(:, :)*(-0.5d0*e2)

    CALL mp_sum(sigmahar, intra_bgrp_comm)

    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)

    RETURN
  END SUBROUTINE esm_stres_har_bc3

  SUBROUTINE esm_stres_ewr(alpha, sigmaewa)
    USE kinds, ONLY : DP
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_ewa must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_ewr_pbc(alpha, sigmaewa)
    CASE ('bc2')
      CALL esm_stres_ewr_pbc(alpha, sigmaewa)
    CASE ('bc3')
      CALL esm_stres_ewr_pbc(alpha, sigmaewa)
    CASE ('bc4')
      STOP 'esm_stres_ewa has not yet implemented for esm_bc = bc4'
    END SELECT

    RETURN
  END SUBROUTINE esm_stres_ewr

  SUBROUTINE esm_stres_ewr_pbc(alpha, sigmaewa)
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base, ONLY : omega, alat, at, tpiba, bg
    USE ions_base, ONLY : zv, nat, tau, ityp
    USE gvect,     ONLY : gstart
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    INTEGER, PARAMETER :: mxr = 50
    ! the maximum number of R vectors included in r sum
    INTEGER  :: ia, ib, nr, nrm, la, mu
    REAL(DP) :: Qa, Qb, dtau(3), rmax
    REAL(DP) :: salp, r(3, mxr), r2(mxr), rr, fac

    salp = sqrt(alpha)

    ! initialize
    sigmaewa(:, :) = 0.d0

    !
    ! R-space sum here (only for the processor that contains G=0)
    !
    IF (gstart == 2) THEN
      rmax = 4.0d0/salp/alat
      !
      ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
      !
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          !
          !     generates nearest-neighbors shells r(i)=R(i)-dtau(i)
          !
          dtau(:) = tau(:, ib) - tau(:, ia)
          CALL rgen(dtau, rmax, mxr, at, bg, r, r2, nrm)

          DO nr = 1, nrm
            rr = sqrt(r2(nr))*alat
            r(:, nr) = r(:, nr)*alat

            fac = Qb*Qa/rr**3 &
                  *(erfc(salp*rr) &
                  + rr*2.0d0*salp*sqrtpm1*EXP(-alpha*rr**2))
            DO la = 1, 3
              DO mu = 1, 3
                sigmaewa(la, mu) = sigmaewa(la, mu) + fac*r(la, nr)*r(mu, nr)
              END DO ! mu
            END DO ! la
          END DO ! nr
        END DO ! ia
      END DO ! ib
    END IF

    sigmaewa(:, :) = sigmaewa(:, :)*(e2/2.0d0/omega)

    CALL mp_sum(sigmaewa, intra_bgrp_comm)

    RETURN
  END SUBROUTINE esm_stres_ewr_pbc

  SUBROUTINE esm_stres_ewg(alpha, sigmaewa)
    USE kinds, ONLY : DP
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_ewa must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_ewg_bc1(alpha, sigmaewa)
    CASE ('bc2')
      CALL esm_stres_ewg_bc2(alpha, sigmaewa)
    CASE ('bc3')
      CALL esm_stres_ewg_bc3(alpha, sigmaewa)
    CASE ('bc4')
      STOP 'esm_stres_ewa must not be called for esm_bc = bc4'
    END SELECT

    RETURN
  END SUBROUTINE esm_stres_ewg

  SUBROUTINE esm_stres_ewg_bc1(alpha, sigmaewa)
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    INTEGER  :: ia, ib, igp, iga, igb, la, mu, iz
    REAL(DP) :: L, S, salp
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP) :: zbza, isalp, gpzbza, gp2a, mgazz, pgazz, fact

    REAL(DP) :: dE_deps(2, 2)
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    salp  = sqrt(alpha)
    isalp = 1.0_DP/salp
    fact  = 0.5_DP * isalp

    ! initialize
    sigmaewa(:, :) = 0.0d0

    !****For gp!=0 case ********************
    DO ib = 1, nat
      Qb = (-1.0d0)*zv(ityp(ib))
      rb(1:2) = tau(1:2, ib)*alat
      zb = tau(3, ib)*alat
      IF (zb > L*0.5d0) THEN
        zb = zb - L
      END IF

      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        ! distance between atoms is defined
        zbza = zb - za

        ! summations over gp
        dE_deps(:, :) = 0.0d0
        DO igp = 1, ngm_2d
          iga = mill_2d(1, igp)
          igb = mill_2d(2, igp)
          g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
          gp = sqrt(g(1)*g(1) + g(2)*g(2))

          IF (gp == 0.0d0) CYCLE ! skip gp=0

          ! define exp phases
          gpzbza = gp * zbza
          gp2a   = gp * 0.5_DP * isalp
          mgazz  = gp2a - salp * zbza
          pgazz  = gp2a + salp * zbza

          ! derivatives by strain tensor
          DO la = 1, 2
            DO mu = 1, 2
              dgp_deps(la, mu) = -g(la)*g(mu)/gp
              dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
            END DO
          END DO

          ! coefficients
          cosgpr = cos(g(1)*(rb(1) - ra(1)) + g(2)*(rb(2) - ra(2)))
          experfcm = exp_erfc(-gpzbza, mgazz)
          experfcp = exp_erfc(+gpzbza, pgazz)
          dexperfcm_dgp = -zbza*experfcm &
                          -exp_gauss( -gpzbza, mgazz ) * fact
          dexperfcp_dgp = +zbza*experfcp &
                          -exp_gauss( +gpzbza, pgazz ) * fact
          !
          ! Old code is not safe, because diverged terms are included.
          ! However, this code is a faithful for original formula.
          ! For this reason, we leave following old codes as comment.
          !
          ! experfcm = exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za))
          ! experfcp = exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za))
          ! dexperfcm_dgp = -(zb - za)*exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za)) &
          !                 - EXP(-gp*(zb - za))*qe_gauss(gp/2.d0/salp - salp*(zb - za))/2.d0/salp
          ! dexperfcp_dgp = +(zb - za)*exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za)) &
          !                 - EXP(+gp*(zb - za))*qe_gauss(gp/2.d0/salp + salp*(zb - za))/2.d0/salp
          !
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcm &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcm &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcm_dgp &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcp &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcp &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcp_dgp
        END DO ! igp

        ! modifications
        IF (gamma_only) THEN
          dE_deps(:, :) = dE_deps(:, :)*2.0d0
        END IF

        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega

      END DO ! ia
    END DO ! ib

    !****For gp=0 case ********************
    IF (gstart == 2) THEN
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        rb(1:2) = tau(1:2, ib)*alat
        zb = tau(3, ib)*alat
        IF (zb > L*0.5d0) THEN
          zb = zb - L
        END IF

        Vr = 0.0d0
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF

          Vr = Vr - tpi*Qa/S &
               *((zb - za)*erf(salp*(zb - za)) &
                 + EXP(-alpha*(zb - za)**2)*sqrtpm1/salp)
        END DO ! ia

        dE_deps(1:2, 1:2) = -delta(1:2, 1:2)*Vr*Qb

        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
      END DO ! ib
    END IF

    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:, :) = sigmaewa(:, :)*(0.5d0*e2)

    CALL mp_sum(sigmaewa, intra_bgrp_comm)

    RETURN
  END SUBROUTINE esm_stres_ewg_bc1

  SUBROUTINE esm_stres_ewg_bc2(alpha, sigmaewa)
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    INTEGER  :: ia, ib, igp, iga, igb, la, mu
    REAL(DP) :: L, S, salp, z0, z1
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP) :: exph1, exph2, exph3
    REAL(DP) :: zbza, isalp, gpzbza, gp2a, mgazz, pgazz, fact

    REAL(DP) :: dE_deps(2, 2)
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    salp = sqrt(alpha)
    isalp = 1.0_DP/salp
    fact  = 0.5_DP * isalp

    ! initialize
    sigmaewa(:, :) = 0.0d0

    !****For gp!=0 case ********************
    DO ib = 1, nat
      Qb = (-1.0d0)*zv(ityp(ib))
      rb(1:2) = tau(1:2, ib)*alat
      zb = tau(3, ib)*alat
      IF (zb > L*0.5d0) THEN
        zb = zb - L
      END IF

      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        ! distance between atoms is defined
        zbza = zb - za

        ! summations over gp
        dE_deps(:, :) = 0.0d0
        DO igp = 1, ngm_2d
          iga = mill_2d(1, igp)
          igb = mill_2d(2, igp)
          g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
          gp = sqrt(g(1)*g(1) + g(2)*g(2))

          IF (gp == 0.0d0) CYCLE ! skip gp=0

          ! define exp phases
          gpzbza = gp * zbza
          gp2a   = gp * 0.5_DP * isalp
          mgazz  = gp2a - salp * zbza
          pgazz  = gp2a + salp * zbza

          ! derivatives by strain tensor
          DO la = 1, 2
            DO mu = 1, 2
              dgp_deps(la, mu) = -g(la)*g(mu)/gp
              dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
            END DO
          END DO

          ! coefficients
          cosgpr = cos(g(1)*(rb(1) - ra(1)) + g(2)*(rb(2) - ra(2)))
          experfcm = exp_erfc(-gpzbza, mgazz)
          experfcp = exp_erfc(+gpzbza, pgazz)
          dexperfcm_dgp = -zbza*experfcm &
                          -exp_gauss( -gpzbza, mgazz ) * fact
          dexperfcp_dgp = +zbza*experfcp &
                          -exp_gauss( +gpzbza, pgazz ) * fact
          !
          ! Old code is not safe, because diverged terms are included.
          ! However, this code is a faithful for original formula.
          ! For this reason, we leave following old codes as comment.
          !
          ! experfcm = exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za))
          ! experfcp = exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za))
          ! dexperfcm_dgp = -(zb - za)*exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za)) &
          !                 - EXP(-gp*(zb - za))*qe_gauss(gp/2.d0/salp - salp*(zb - za))/2.d0/salp
          ! dexperfcp_dgp = +(zb - za)*exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za)) &
          !                 - EXP(+gp*(zb - za))*qe_gauss(gp/2.d0/salp + salp*(zb - za))/2.d0/salp
          !
          exph1 = (cosh(gp*(zb - za))*EXP(-2*gp*z1) - cosh(gp*(zb + za)))/sinh(2*gp*z1)
          exph2 = ((zb - za)*sinh(gp*(zb - za))*EXP(-2*gp*z1) &
                  - 2*z1*cosh(gp*(zb - za))*EXP(-2*gp*z1) &
                  - (zb + za)*sinh(gp*(zb + za)))/sinh(2*gp*z1)
          exph3 = -(cosh(gp*(zb - za))*EXP(-2*gp*z1) - cosh(gp*(zb + za)))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)

          !! BC1 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcm &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcm &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcm_dgp &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcp &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcp &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcp_dgp

          !! BC2 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*tpi/gp*Qb*Qa/S*cosgpr*exph1 &
                          - tpi/gp*delta(:, :)*Qb*Qa/S*cosgpr*exph1 &
                          + tpi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*(exph2 + exph3)
        END DO ! igp

        ! modifications
        IF (gamma_only) THEN
          dE_deps(:, :) = dE_deps(:, :)*2.0d0
        END IF

        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega

      END DO ! ia
    END DO ! ib

    !****For gp=0 case ********************
    IF (gstart == 2) THEN
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        rb(1:2) = tau(1:2, ib)*alat
        zb = tau(3, ib)*alat
        IF (zb > L*0.5d0) THEN
          zb = zb - L
        END IF

        ! [note] this Vr does not contain a term due to efield z*efield
        ! because it vanishes in the differentiation with respect to strain.
        Vr = 0.0d0
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF

          !! BC1 terms
          Vr = Vr - tpi*Qa/S &
               *((zb - za)*erf(salp*(zb - za)) &
                 + EXP(-alpha*(zb - za)**2)*sqrtpm1/salp)

          !! BC2 terms
          Vr = Vr + tpi*Qa/S*(-zb*za + z1*z1)/z1
        END DO ! ia

        dE_deps(1:2, 1:2) = -delta(1:2, 1:2)*Vr*Qb

        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
      END DO ! ib
    END IF

    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:, :) = sigmaewa(:, :)*(0.5d0*e2)

    CALL mp_sum(sigmaewa, intra_bgrp_comm)

    RETURN
  END SUBROUTINE esm_stres_ewg_bc2

  SUBROUTINE esm_stres_ewg_bc3(alpha, sigmaewa)
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)

    INTEGER  :: ia, ib, igp, iga, igb, la, mu
    REAL(DP) :: L, S, salp, z0, z1
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp, expm
    REAL(DP) :: zbza, isalp, gpzbza, gp2a, mgazz, pgazz, fact

    REAL(DP) :: dE_deps(2, 2)
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    salp = sqrt(alpha)
    isalp = 1.0_DP/salp
    fact  = 0.5_DP * isalp

    ! initialize
    sigmaewa(:, :) = 0.0d0

    !****For gp!=0 case ********************
    DO ib = 1, nat
      Qb = (-1.0d0)*zv(ityp(ib))
      rb(1:2) = tau(1:2, ib)*alat
      zb = tau(3, ib)*alat
      IF (zb > L*0.5d0) THEN
        zb = zb - L
      END IF

      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        ! distance between atoms is defined
        zbza = zb - za

        ! summations over gp
        dE_deps(:, :) = 0.0d0
        DO igp = 1, ngm_2d
          iga = mill_2d(1, igp)
          igb = mill_2d(2, igp)
          g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
          gp = sqrt(g(1)*g(1) + g(2)*g(2))

          IF (gp == 0.0d0) CYCLE ! skip gp=0

          ! define exp phases
          gpzbza = gp * zbza
          gp2a   = gp * 0.5_DP * isalp
          mgazz  = gp2a - salp * zbza
          pgazz  = gp2a + salp * zbza

          ! derivatives by strain tensor
          DO la = 1, 2
            DO mu = 1, 2
              dgp_deps(la, mu) = -g(la)*g(mu)/gp
              dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
            END DO
          END DO

          ! coefficients
          cosgpr = cos(g(1)*(rb(1) - ra(1)) + g(2)*(rb(2) - ra(2)))
          experfcm = exp_erfc(-gpzbza, mgazz)
          experfcp = exp_erfc(+gpzbza, pgazz)
          dexperfcm_dgp = -zbza*experfcm &
                          -exp_gauss( -gpzbza, mgazz ) * fact
          dexperfcp_dgp = +zbza*experfcp &
                          -exp_gauss( +gpzbza, pgazz ) * fact
          !
          ! Old code is not safe, because diverged terms are included.
          ! However, this code is a faithful for original formula.
          ! For this reason, we leave following old codes as comment.
          !
          ! experfcm = exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za))
          ! experfcp = exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za))
          ! dexperfcm_dgp = -(zb - za)*exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za)) &
          !                 - EXP(-gp*(zb - za))*qe_gauss(gp/2.d0/salp - salp*(zb - za))/2.d0/salp
          ! dexperfcp_dgp = +(zb - za)*exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za)) &
          !                 - EXP(+gp*(zb - za))*qe_gauss(gp/2.d0/salp + salp*(zb - za))/2.d0/salp
          expm = EXP(-gp*(-zb + 2*z1 - za))
          !
          !! BC1 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcm &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcm &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcm_dgp &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcp &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcp &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcp_dgp

          !! BC3 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          - gp*dinvgp_deps(:, :)*tpi/gp*Qb*Qa/S*cosgpr*expm &
                          + tpi/gp*delta(:, :)*Qb*Qa/S*cosgpr*expm &
                          + tpi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*(-zb + 2*z1 - za)*expm
        END DO ! igp

        ! modifications
        IF (gamma_only) THEN
          dE_deps(:, :) = dE_deps(:, :)*2.0d0
        END IF

        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega

      END DO ! ia
    END DO ! ib

    !****For gp=0 case ********************
    IF (gstart == 2) THEN
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        rb(1:2) = tau(1:2, ib)*alat
        zb = tau(3, ib)*alat
        IF (zb > L*0.5d0) THEN
          zb = zb - L
        END IF

        Vr = 0.0d0
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF

          !! BC1 terms
          Vr = Vr - tpi*Qa/S &
               *((zb - za)*erf(salp*(zb - za)) &
               + EXP(-alpha*(zb - za)**2)*sqrtpm1/salp)

          !! BC3 terms
          Vr = Vr + tpi*Qa/S*(-zb + 2*z1 - za)
        END DO ! ia

        dE_deps(1:2, 1:2) = -delta(1:2, 1:2)*Vr*Qb

        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
      END DO ! ib
    END IF

    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:, :) = sigmaewa(:, :)*(0.5d0*e2)

    CALL mp_sum(sigmaewa, intra_bgrp_comm)

    RETURN
  END SUBROUTINE esm_stres_ewg_bc3

  SUBROUTINE esm_stres_loclong_bc1(sigmaloclong, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(out) :: sigmaloclong(3, 3)
    COMPLEX(DP) :: rhog(ngm)   !  n(G)

    INTEGER  :: ig, iga, igb, igz, igp, la, mu, iz, jz, ia
    REAL(DP) :: L, S, z0, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP)    :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
    REAL(DP) :: isalp, fact
    REAL(DP) :: zza, mgza, pgza, g2a_maza, g2a_paza

    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))

    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)

    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3

      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      ENDIF
    END DO ! ig

    ! cell settings
    L = at(3, 3)*alat
    S     = omega/L
    z0    = L/2.d0
    alpha = 1.0d0
    salp  = sqrt(alpha)
    ! useful values are setted
    isalp = 1.d0/sqrt(alpha)
    fact  = 1.d0/2.d0/salp
    ! initialize
    sigmaloclong(:, :) = 0.0d0

    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))

      IF (gp == 0.0d0) CYCLE ! skip gp=0

      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO

      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        ! summations over all atoms
        dVr_deps(iz, :, :) = (0.0d0, 0.0d0)
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
! --------------------------------------------------------------------------------------------------
!         Following code is old version. However, this code explicitly shows the formulation of 
!         stress tensor within ESM scheme. For this reason, we left the code as comment.
!          expimgpr = qe_exp(-ci*(g(1)*ra(1) + g(2)*ra(2)))
!          experfcm = exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za))
!          experfcp = exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za))
!          dexperfcm_dgp = -(z - za)*exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za)) &
!                          - EXP(-gp*(z - za))*qe_gauss(gp/2.d0/salp - salp*(z - za))/2.d0/salp
!          dexperfcp_dgp = +(z - za)*exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za)) &
!                          - EXP(+gp*(z - za))*qe_gauss(gp/2.d0/salp + salp*(z - za))/2.d0/salp
!---------------------------------------------------------------------------------------------------
          !
          ! ... Set useful values
          zza = z - za
          mgza = -gp*zza
          pgza = +gp*zza
          g2a_maza = gp*0.5d0*isalp - salp * zza
          g2a_paza = gp*0.5d0*isalp + salp * zza
          !
          expimgpr = QE_EXP(-ci*(g(1)*ra(1) + g(2)*ra(2)))
          experfcm = exp_erfc(mgza, g2a_maza)
          experfcp = exp_erfc(pgza, g2a_paza)
          dexperfcm_dgp = -zza * experfcm &
                          -exp_gauss( mgza, g2a_maza ) * fact
          dexperfcp_dgp = +zza * experfcp &
                          -exp_gauss( pgza, g2a_paza ) * fact
          !
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcm &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcm &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcm_dgp

          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcp &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcp &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcp_dgp
        END DO ! ia
      END DO ! iz

      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO

      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END DO ! igp

    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! calculate V(z)
      Vr(:) = 0.0d0
      ! separation by polynomial
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      z_l = -z0
      z_r = +z0
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        DO iz = 1, dfftp%nr3
          jz = iz - 1
          IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            jz = jz - dfftp%nr3
          END IF
          z = DBLE(jz) / DBLE(dfftp%nr3) * L

          Vr(iz) = Vr(iz) - tpi*Qa/S &
                   *((z - za)*erf(salp*(z - za)) &
                     + EXP(-alpha*(z - za)**2)*sqrtpm1/salp)
        END DO ! iz

        f1 = f1 - tpi*Qa/S &
             *((z_r - za)*erf(salp*(z_r - za)) &
               + EXP(-alpha*(z_r - za)**2)*sqrtpm1/salp)
        f2 = f2 - tpi*Qa/S &
             *((z_l - za)*erf(salp*(z_l - za)) &
               + EXP(-alpha*(z_l - za)**2)*sqrtpm1/salp)
        f3 = f3 - tpi*Qa/S &
             *erf(salp*(z_r - za))
        f4 = f4 - tpi*Qa/S &
             *erf(salp*(z_l - za))
      END DO ! ia

      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3

      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)

      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0

      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    ENDIF ! imill_2d(0,0) > 0

    ! e2 means hartree -> Ry.
    sigmaloclong(:, :) = sigmaloclong(:, :)*(e2)

    CALL mp_sum(sigmaloclong, intra_bgrp_comm)

    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)

    RETURN
  END SUBROUTINE esm_stres_loclong_bc1

  SUBROUTINE esm_stres_loclong_bc2(sigmaloclong, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(out) :: sigmaloclong(3, 3)
    COMPLEX(DP) :: rhog(ngm)   !  n(G)

    INTEGER  :: ig, iga, igb, igz, igp, la, mu, iz, jz, ia
    REAL(DP) :: L, S, z0, z1, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    COMPLEX(DP) :: exph1, exph2, exph3
    REAL(DP)    :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
    REAL(DP) :: isalp, fact
    REAL(DP) :: zza, mgza, pgza, g2a_maza, g2a_paza

    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))

    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)

    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3

      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      ENDIF
    END DO ! ig

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    alpha = 1.0d0
    salp = sqrt(alpha)
    ! useful values are setted
    isalp = 1.d0/sqrt(alpha)
    fact  = 1.d0/2.d0/salp
    ! initialize
    sigmaloclong(:, :) = 0.0d0

    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))

      IF (gp == 0.0d0) CYCLE ! skip gp=0

      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO

      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        ! summations over all atoms
        dVr_deps(iz, :, :) = (0.0d0, 0.0d0)
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
! --------------------------------------------------------------------------------------------------
!         Following code is old version. However, this code explicitly shows the formulation of 
!         stress tensor within ESM scheme. For this reason, we left the code as comment.
!          expimgpr = qe_exp(-ci*(g(1)*ra(1) + g(2)*ra(2)))
!          experfcm = exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za))
!          experfcp = exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za))
!          dexperfcm_dgp = -(z - za)*exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za)) &
!                          - EXP(-gp*(z - za))*qe_gauss(gp/2.d0/salp - salp*(z - za))/2.d0/salp
!          dexperfcp_dgp = +(z - za)*exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za)) &
!                          - EXP(+gp*(z - za))*qe_gauss(gp/2.d0/salp + salp*(z - za))/2.d0/salp
!---------------------------------------------------------------------------------------------------
          !
          ! ... Set useful values
          zza = z - za
          mgza = -gp*zza
          pgza = +gp*zza
          g2a_maza = gp*0.5d0*isalp - salp * zza
          g2a_paza = gp*0.5d0*isalp + salp * zza
          !
          expimgpr = QE_EXP(-ci*(g(1)*ra(1) + g(2)*ra(2)))
          experfcm = exp_erfc(mgza, g2a_maza)
          experfcp = exp_erfc(pgza, g2a_paza)
          dexperfcm_dgp = -zza * experfcm &
                          -exp_gauss( mgza, g2a_maza ) * fact
          dexperfcp_dgp = +zza * experfcp &
                          -exp_gauss( pgza, g2a_paza ) * fact
          !
          exph1 = (cosh(gp*(z - za))*EXP(-2*gp*z1) - cosh(gp*(z + za)))/sinh(2*gp*z1)
          exph2 = ((z - za)*sinh(gp*(z - za))*EXP(-2*gp*z1) &
                   - 2*z1*cosh(gp*(z - za))*EXP(-2*gp*z1) &
                   - (z + za)*sinh(gp*(z + za)))/sinh(2*gp*z1)
          exph3 = -(cosh(gp*(z - za))*EXP(-2*gp*z1) - cosh(gp*(z + za)))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)

          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcm &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcm &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcm_dgp

          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcp &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcp &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcp_dgp

          !! BC2 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*tpi/gp*Qa/S*expimgpr*exph1 &
                               - tpi/gp*delta(:, :)*Qa/S*expimgpr*exph1 &
                               + tpi/gp*Qa/S*expimgpr*dgp_deps(:, :)*(exph2 + exph3)
        END DO ! ia
      END DO ! iz

      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO

      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END DO ! igp

    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! calculate V(z)
      ! [note] this Vr does not contain a term due to efield z*efield
      ! because it vanishes in the differentiation with respect to strain.
      Vr(:) = 0.0d0
      ! separation by polynomial
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      z_l = -z0
      z_r = +z0
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        DO iz = 1, dfftp%nr3
          jz = iz - 1
          IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            jz = jz - dfftp%nr3
          END IF
          z = DBLE(jz) / DBLE(dfftp%nr3) * L

          !! BC1 terms
          Vr(iz) = Vr(iz) - tpi*Qa/S &
                   *((z - za)*erf(salp*(z - za)) &
                     + EXP(-alpha*(z - za)**2)*sqrtpm1/salp)

          !! BC2 terms
          Vr(iz) = Vr(iz) + tpi*Qa/S*(-z*za + z1*z1)/z1
        END DO ! iz

        f1 = f1 - tpi*Qa/S &
             *((z_r - za)*erf(salp*(z_r - za)) &
               + EXP(-alpha*(z_r - za)**2)*sqrtpm1/salp)
        f1 = f1 + tpi*Qa/S*(-z_r*za + z1*z1)/z1

        f2 = f2 - tpi*Qa/S &
             *((z_l - za)*erf(salp*(z_l - za)) &
               + EXP(-alpha*(z_l - za)**2)*sqrtpm1/salp)
        f2 = f2 + tpi*Qa/S*(-z_l*za + z1*z1)/z1

        f3 = f3 - tpi*Qa/S &
             *erf(salp*(z_r - za))
        f3 = f3 + tpi*Qa/S*(-za)/z1

        f4 = f4 - tpi*Qa/S &
             *erf(salp*(z_l - za))
        f4 = f4 + tpi*Qa/S*(-za)/z1

      END DO ! ia

      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3

      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)

      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0

      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    ENDIF ! imill_2d(0,0) > 0

    ! e2 means hartree -> Ry.
    sigmaloclong(:, :) = sigmaloclong(:, :)*(e2)

    CALL mp_sum(sigmaloclong, intra_bgrp_comm)

    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)

    RETURN
  END SUBROUTINE esm_stres_loclong_bc2

  SUBROUTINE esm_stres_loclong_bc3(sigmaloclong, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE

    REAL(DP), INTENT(out) :: sigmaloclong(3, 3)
    COMPLEX(DP) :: rhog(ngm)   !  n(G)

    INTEGER  :: ig, iga, igb, igz, igp, la, mu, iz, jz, ia
    REAL(DP) :: L, S, z0, z1, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    COMPLEX(DP) :: expm
    REAL(DP)    :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
    REAL(DP) :: isalp, fact
    REAL(DP) :: zza, mgza, pgza, g2a_maza, g2a_paza, tz1

    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)

    REAL(DP) :: sigmaloclong_bc1(3, 3)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))

    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)

    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3

      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      ENDIF
    END DO ! ig

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    alpha = 1.0d0
    salp = sqrt(alpha)
    ! useful values are setted
    isalp = 1.d0/sqrt(alpha)
    fact  = 1.d0/2.d0/salp
    tz1 = 2.d0 * z1
    !
    ! initialize
    sigmaloclong(:, :) = 0.0d0

    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))

      IF (gp == 0.0d0) CYCLE ! skip gp=0

      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO

      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L

        ! summations over all atoms
        dVr_deps(iz, :, :) = (0.0d0, 0.0d0)
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
! --------------------------------------------------------------------------------------------------
!         Following code is old version. However, this code explicitly shows the formulation of 
!         stress tensor within ESM scheme. For this reason, we left the code as comment.
!          expimgpr = qe_exp(-ci*(g(1)*ra(1) + g(2)*ra(2)))
!          experfcm = exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za))
!          experfcp = exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za))
!          dexperfcm_dgp = -(z - za)*exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za)) &
!                          - EXP(-gp*(z - za))*qe_gauss(gp/2.d0/salp - salp*(z - za))/2.d0/salp
!          dexperfcp_dgp = +(z - za)*exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za)) &
!                          - EXP(+gp*(z - za))*qe_gauss(gp/2.d0/salp + salp*(z - za))/2.d0/salp
!---------------------------------------------------------------------------------------------------
          !
          ! ... Set useful values
          zza = z - za
          mgza = -gp*zza
          pgza = +gp*zza
          g2a_maza = gp*0.5d0*isalp - salp * zza
          g2a_paza = gp*0.5d0*isalp + salp * zza
          !
          expimgpr = QE_EXP(-ci*(g(1)*ra(1) + g(2)*ra(2)))
          experfcm = exp_erfc(mgza, g2a_maza)
          experfcp = exp_erfc(pgza, g2a_paza)
          dexperfcm_dgp = -zza * experfcm &
                          -exp_gauss( mgza, g2a_maza ) * fact
          dexperfcp_dgp = +zza * experfcp &
                          -exp_gauss( pgza, g2a_paza ) * fact
          !
          expm = EXP(-gp*(-z + tz1 - za))

          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcm &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcm &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcm_dgp

          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcp &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcp &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcp_dgp

          !! BC3 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               - gp*dinvgp_deps(:, :)*tpi/gp*Qa/S*expimgpr*expm &
                               + tpi/gp*delta(:, :)*Qa/S*expimgpr*expm &
                               + tpi/gp*Qa/S*expimgpr*dgp_deps(:, :)*(-z + 2*z1 - za)*expm
        END DO ! ia
      END DO ! iz

      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO

      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END DO ! igp

    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! calculate V(z)
      Vr(:) = 0.0d0
      ! separation by polynomial
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      z_l = -z0
      z_r = +z0
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        DO iz = 1, dfftp%nr3
          jz = iz - 1
          IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            jz = jz - dfftp%nr3
          END IF
          z = DBLE(jz) / DBLE(dfftp%nr3) * L

          !! BC1 terms
          Vr(iz) = Vr(iz) - tpi*Qa/S &
                   *((z - za)*erf(salp*(z - za)) &
                   + EXP(-alpha*(z - za)**2)*sqrtpm1/salp)

          !! BC3 terms
          Vr(iz) = Vr(iz) + tpi*Qa/S*(-z + 2*z1 - za)
        END DO ! iz

        f1 = f1 - tpi*Qa/S &
             *((z_r - za)*erf(salp*(z_r - za)) &
             + EXP(-alpha*(z_r - za)**2)*sqrtpm1/salp)
        f1 = f1 + tpi*Qa/S*(-z_r + 2*z1 - za)

        f2 = f2 - tpi*Qa/S &
             *((z_l - za)*erf(salp*(z_l - za)) &
             + EXP(-alpha*(z_l - za)**2)*sqrtpm1/salp)
        f2 = f2 + tpi*Qa/S*(-z_l + 2*z1 - za)

        f3 = f3 - tpi*Qa/S &
             *erf(salp*(z_r - za))
        f3 = f3 + tpi*Qa/S*(-1.0d0)

        f4 = f4 - tpi*Qa/S &
             *erf(salp*(z_l - za))
        f4 = f4 + tpi*Qa/S*(-1.0d0)
      END DO ! ia

      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3

      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)

      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L

        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0

      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz

      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    ENDIF ! imill_2d(0,0) > 0

    ! e2 means hartree -> Ry.
    sigmaloclong(:, :) = sigmaloclong(:, :)*(e2)

    CALL mp_sum(sigmaloclong, intra_bgrp_comm)

    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)

    RETURN
  END SUBROUTINE esm_stres_loclong_bc3

  COMPLEX(DP) FUNCTION qe_exp(x)
    COMPLEX(DP), INTENT(in) :: x
    REAL(DP) :: r, i, c, s

    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)

    qe_exp = EXP(r)*cmplx(c, s, kind=DP)

  END FUNCTION qe_exp

  COMPLEX(DP) FUNCTION qe_sinh(x)
    COMPLEX(DP), INTENT(in) :: x
    REAL(DP) :: r, i, c, s

    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)

    qe_sinh = 0.5d0*(EXP(r)*cmplx(c, s, kind=DP) - EXP(-r)*cmplx(c, -s, kind=DP))

  END FUNCTION qe_sinh

  COMPLEX(DP) FUNCTION qe_cosh(x)
    COMPLEX(DP), INTENT(in) :: x
    REAL(DP) :: r, i, c, s

    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)

    qe_cosh = 0.5d0*(EXP(r)*cmplx(c, s, kind=DP) + EXP(-r)*cmplx(c, -s, kind=DP))

  END FUNCTION qe_cosh

  FUNCTION qe_gauss(x) result(gauss)
    USE kinds,     ONLY : DP
    USE constants, ONLY : sqrtpm1  ! 1/sqrt(pi)
    IMPLICIT NONE

    REAL(DP), INTENT(in) :: x
    REAL(DP) :: gauss

    gauss = 2.0d0*sqrtpm1*EXP(-x*x)

  END FUNCTION qe_gauss

  FUNCTION exp_gauss( x, y )
    USE kinds,     ONLY : DP
    USE constants, ONLY : sqrtpm1 !1/sqrt(pi)

    REAL(DP), INTENT(IN) :: x, y
    REAL(DP) :: exp_gauss

    exp_gauss = 2._DP*sqrtpm1*EXP( x - y*y )

  END FUNCTION exp_gauss
END MODULE esm_stres_mod
