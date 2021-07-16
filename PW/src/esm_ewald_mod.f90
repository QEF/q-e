MODULE esm_ewald_mod

  USE kinds, ONLY : DP
  USE esm_common_mod, ONLY : esm_efield, esm_w, esm_a, esm_bc, mill_2d, ngm_2d, &
                             vl11j0, vl12j0, vl21j0, vl22j0, vl11, vl22, &
                             qromb, exp_erfc, esm_rgen_2d
  IMPLICIT NONE

CONTAINS

  FUNCTION esm_ewald()
    !-----------------------------------------------------------------------
    !
    ! Calculates Ewald energy with both G- and R-space terms.
    ! Determines optimal alpha. Should hopefully work for any structure.
    !
    USE kinds,     ONLY : DP
    USE constants, ONLY : tpi, e2
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    USE cell_base, ONLY : tpiba2
    USE ions_base, ONLY : zv, nat, ityp
    USE gvect,     ONLY : gcutm
    IMPLICIT NONE

    REAL(DP) :: esm_ewald
    ! output: the ewald energy
    !
    !    here the local variables
    !
    INTEGER :: na
    ! counter on atoms

    REAL(DP) :: charge, ewaldg, ewaldr, alpha, upperbound
    ! total ionic charge in the cell
    ! ewald energy computed in reciprocal space
    ! ewald energy computed in real space
    ! alpha term in ewald sum
    ! the maximum radius to consider real space sum

    charge = 0.d0
    DO na = 1, nat
      charge = charge + zv(ityp(na))
    ENDDO

    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    alpha = 2.9d0
    DO
      alpha = alpha - 0.1d0
      IF (alpha .le. 0.d0) CALL errore('esm_ewald', 'optimal alpha not found', 1)
      upperbound = 2.d0*charge**2*sqrt(2.d0*alpha/tpi)* &
                   erfc(sqrt(tpiba2*gcutm/4.d0/alpha))
      IF (upperbound < 1.0d-7) EXIT
    END DO

    ! G-space sum here.
    ! Determine if this processor contains G=0 and set the constant term
    CALL esm_ewaldg(alpha, ewaldg)

    ! R-space sum here (only for the processor that contains G=0)
    CALL esm_ewaldr(alpha, ewaldr)

    esm_ewald = 0.5d0*e2*(ewaldg + ewaldr)

    CALL mp_sum(esm_ewald, intra_bgrp_comm)
    !write( *,'(5x,"alpha used in ewald term: ",f5.2 )')alpha

    RETURN
  END FUNCTION esm_ewald

  SUBROUTINE esm_ewaldr(alpha_g, ewr)
    USE kinds, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewr

    IF (esm_bc == 'pbc') THEN
      CALL esm_ewaldr_pbc(alpha_g, ewr)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_ewaldr_pbc(alpha_g, ewr)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_ewaldr_pbc(alpha_g, ewr)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_ewaldr_pbc(alpha_g, ewr)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_ewaldr_bc4(alpha_g, ewr)
    END IF

  END SUBROUTINE esm_ewaldr

  SUBROUTINE esm_ewaldg(alpha_g, ewg)
    USE kinds, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewg

    IF (esm_bc == 'pbc') THEN
      CALL esm_ewaldg_pbc(alpha_g, ewg)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_ewaldg_bc1(alpha_g, ewg)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_ewaldg_bc2(alpha_g, ewg)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_ewaldg_bc3(alpha_g, ewg)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_ewaldg_bc4(alpha_g, ewg)
    END IF

  END SUBROUTINE esm_ewaldg

!-----------------------------------------------------------------------
!--------------ESM EWALD RSUM SUBROUTINE--------------------------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_ewaldr_pbc(alpha_g, ewr)

    USE io_global,        ONLY : stdout
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : iverbosity
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewr
    !
    !    here the local variables
    !
    INTEGER            :: na, nb, nr, nrm, np, ip
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER :: mxr = 500
    ! the maximum number of R vectors included in r
    REAL(DP)           :: dtau(3), r(3, mxr), r2(mxr)
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    ! buffer variable
    ! buffer variable
    !
    ! ESM variables
    !
    REAL(DP)            :: tmp, fac, ss, ew, rmax0, rr
    !
    ewr = 0.d0

    tmp = sqrt(alpha_g)
    rmax0 = 4.d0/tmp/alat

    ip = mp_rank(intra_bgrp_comm)
    np = mp_size(intra_bgrp_comm)

    ew = 0.d0
    DO na = ip + 1, nat, np
      DO nb = 1, nat
        dtau(:) = tau(:, na) - tau(:, nb)
        fac = zv(ityp(nb))*zv(ityp(na))
        !
        ! generates nearest-neighbors shells
        !
        CALL rgen(dtau, rmax0, mxr, at, bg, r, r2, nrm)
        !
        ! and sum to the real space part
        !
        DO nr = 1, nrm
          rr = sqrt(r2(nr))*alat
          ew = ew + fac*erfc(tmp*rr)/rr
        ENDDO
      ENDDO
      ! Here add the other constant term
      ew = ew - zv(ityp(na))**2*tmp/sqrt(pi)*2.d0 ! 2.d0: fit to original code
    ENDDO
    ewr = ewr + ew

  END SUBROUTINE esm_ewaldr_pbc

  SUBROUTINE esm_ewaldr_bc4(alpha_g, ewr)

    USE io_global,        ONLY : stdout
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : iverbosity
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewr
    !
    !    here the local variables
    !
    INTEGER            :: na, nb, nr, nrm, np, ip
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER :: mxr = 500
    ! the maximum number of R vectors included in r
    REAL(DP)           :: dtau(3), r(3, mxr), r2(mxr), rxy, rxyz
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    ! buffer variable
    ! buffer variable
    !
    ! ESM variables
    !
    REAL(DP)            :: L, z, zp, z0, z1, aaa, tmp, fac, ss, ew, err, ss0, &
                           gpmax, rmax0, rmax, zbuff, znrm, rr
    ! gpmax: upper bound of g_parallel integral
    ! rmax: the maximum radius to consider real space sum
    ! zbuff: smearing width to avoid the singularity of the Force
    ! znrm: threashold value for normal RSUM and Smooth-ESM's RSUM
    REAL(DP), PARAMETER :: eps = 1.d-11, epsneib = 1.d-6
    !
    ewr = 0.d0
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    aaa = esm_a
    tmp = sqrt(alpha_g)
    zbuff = 1.d0
    !
    ! Define upperbound for g_parallel integral
    err = 1.d0; ss0 = 0.d0; gpmax = 1.d0
    DO
      gpmax = gpmax + 1.d0
      IF (gpmax .gt. 1000.d0) &
        CALL errore('esm_ewaldr', 'optimal gpmax not found', 1)
      CALL qromb(vl11, aaa, tmp, z1, z1 - zbuff, z1 - zbuff, 0.0_DP, gpmax, ss)
      err = ABS(ss - ss0); ss0 = ss
      IF (err .lt. eps) EXIT
    ENDDO
    ! Define znrm using the deviation from the constant term in RSUM
    znrm = z1
    DO
      znrm = znrm - 0.01d0
      IF (znrm .le. -z0) &
        CALL errore('esm_ewaldr', 'optimal znrm not found', 1)
      CALL qromb(vl11, aaa, tmp, z1, znrm, znrm, 0.0_DP, gpmax, ss)
      err = -2.d0*tmp/sqrt(pi) - ss*2.d0
      IF (ABS(err) .lt. eps) EXIT
    ENDDO
    ! Define rmax for real space sum
    rmax = 1.d0
    DO
      rmax = rmax + 1.d0
      IF (rmax .gt. 200.d0) &
        CALL errore('esm_ewaldr', 'optimal rmax not found', 1)
      CALL qromb(vl11j0, aaa, tmp, z1, z1 - zbuff, z1 - zbuff, rmax, gpmax, ss)
      err = 1.d0/rmax + ss*2.d0
      IF (ABS(err) .lt. epsneib) EXIT
    ENDDO
    rmax = rmax/alat
    IF (iverbosity > 0) THEN
      write (stdout, '(5x,"=== Smooth-ESM RSUM parameters (Energy) ===")')
      write (stdout, '(5x,A,F10.2,A)') &
        'Upper bound of g_parallel integral:      ', gpmax, ' (1/a.u.)'
      write (stdout, '(5x,A,F10.2,A)') &
        'Boundary for normal RSUM|Smooth-ESM RSUM:', z1 - znrm, ' (a.u.)'
      write (stdout, '(5x,A,F10.2,A)') &
        'Upper bound of real-space summation:     ', rmax*alat, ' (a.u.)'
      write (stdout, '(5x,"===========================================")')
    ENDIF

    ip = mp_rank(intra_bgrp_comm)
    np = mp_size(intra_bgrp_comm)

    ew = 0.d0
    DO na = ip + 1, nat, np
      z = tau(3, na)
      IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
      z = z*alat
      DO nb = 1, nat
        zp = tau(3, nb)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        dtau(1:2) = tau(1:2, na) - tau(1:2, nb)
        dtau(3) = (z - zp)/alat
        fac = zv(ityp(nb))*zv(ityp(na))
        IF (z < znrm) THEN
          IF (zp < znrm) THEN ! z in I, zp in I (normal RSUM)
            rmax0 = 4.d0/tmp/alat
            !
            ! generates nearest-neighbors shells
            !
            CALL rgen(dtau, rmax0, mxr, at, bg, r, r2, nrm)
            !
            ! and sum to the real space part
            !
            DO nr = 1, nrm
              rr = sqrt(r2(nr))*alat
              ew = ew + fac*erfc(tmp*rr)/rr
            ENDDO
          ELSEIF (zp < z1) THEN ! z in I, zp in I
            CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              rxyz = sqrt(r2(nr) + dtau(3)**2)*alat
              CALL qromb(vl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              ew = ew + fac*(1.d0/rxyz + ss*2.d0)
            ENDDO
          ELSE ! z in I, zp in II
            CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              CALL qromb(vl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              ew = ew + fac*ss*2.d0
            ENDDO
          ENDIF ! IF for zp
        ELSEIF (z < z1) THEN ! znrm < z < z1
          CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
          IF (zp < z1) THEN ! z in I, zp in I
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              rxyz = sqrt(r2(nr) + dtau(3)**2)*alat
              CALL qromb(vl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              ew = ew + fac*(1.d0/rxyz + ss*2.d0)
            ENDDO
          ELSE ! z in I, zp in II
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              CALL qromb(vl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              ew = ew + fac*ss*2.d0
            ENDDO
          ENDIF ! IF for zp
        ELSE ! z1 < z
          CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
          IF (zp < z1) THEN ! z in II, zp in I
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              CALL qromb(vl21j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              ew = ew + fac*ss*2.d0
            ENDDO
          ELSE ! z in II, zp in II
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              rxyz = sqrt(r2(nr) + dtau(3)**2)*alat
              CALL qromb(vl22j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              ew = ew + fac*(EXP(-aaa*(rxyz + z + zp - 2.d0*z1))/rxyz + ss*2.d0)
            ENDDO
          ENDIF
        ENDIF ! IF for z
      ENDDO
      IF (z < znrm) THEN
        ss = -tmp/sqrt(pi)
      ELSEIF (z < z1) THEN
        CALL qromb(vl11, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss)
      ELSE
        CALL qromb(vl22, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss)
      ENDIF
      ew = ew + zv(ityp(na))**2*ss*2.d0 ! 2.0: fit to original code
    ENDDO
    ewr = ewr + ew

  END SUBROUTINE esm_ewaldr_bc4

!-----------------------------------------------------------------------
!--------------ESM EWALD GSUM SUBROUTINE--------------------------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_ewaldg_pbc(alpha_g, ewg)

    USE constants,        ONLY : tpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, tpiba2
    USE ions_base,        ONLY : zv, nat, nsp, ityp
    USE control_flags,    ONLY : gamma_only
    USE gvect,            ONLY : ngm, gg
    USE vlocal,           ONLY : strf

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewg

    INTEGER     :: ng
    REAL(DP)    :: charge, fact
    COMPLEX(DP) :: rhon

    charge = sum(zv(ityp(1:nat)))

    ! same of the GSUM part in ewald.f90
    IF (gstart == 2) THEN
      ewg = -charge**2/alpha_g/4.0d0
    ELSE
      ewg = 0.0d0
    ENDIF
    IF (gamma_only) THEN
      fact = 2.d0
    ELSE
      fact = 1.d0
    END IF
    DO ng = gstart, ngm
      rhon = sum(zv(1:nsp)*CONJG(strf(ng, 1:nsp)))
      ewg = ewg + fact*ABS(rhon)**2* &
            EXP(-gg(ng)*tpiba2/alpha_g/4.d0)/gg(ng)/tpiba2
    ENDDO
    ewg = 2.d0*tpi/omega*ewg

  END SUBROUTINE esm_ewaldg_pbc

  SUBROUTINE esm_ewaldg_bc1(alpha_g, ewg)

    USE constants,        ONLY : pi, tpi, fpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewg
    !
    !    here the local variables
    !
    INTEGER  :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
                arg101, arg102

    ewg = 0.d0
    L = at(3, 3)*alat
    z0 = L/2.d0
    tmp = sqrt(alpha_g)
    sa = omega/L
    ew = 0d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc1
        arg001 = -tmp**2*(z - zp)**2
        arg101 = tmp*(z - zp)
        kk1 = 0.5d0*(-(z - zp)*erf(arg101) - EXP(arg001)/tmp/sqrt(pi))
        kk2 = 0.d0

        cc1 = 0.d0
        cc2 = 0.d0
        DO ng_2d = 1, ngm_2d
          k1 = mill_2d(1, ng_2d)
          k2 = mill_2d(2, ng_2d)
          IF (k1 == 0 .and. k2 == 0) CYCLE
          t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
          gp2 = sum(t(:)*t(:))*tpiba2
          gp = sqrt(gp2)
          ff = ((k1*bg(1, 1) + k2*bg(1, 2))*(tau(1, it1) - tau(1, it2)) &
                + (k1*bg(2, 1) + k2*bg(2, 2))*(tau(2, it1) - tau(2, it2)))*tpi
          ! bc1
          arg001 = -gp*(z - zp)
          arg002 = gp*(z - zp)
          arg101 = gp/2.d0/tmp - tmp*(z - zp)
          arg102 = gp/2.d0/tmp + tmp*(z - zp)
          t1 = exp_erfc(arg001, arg101)
          t2 = exp_erfc(arg002, arg102)
          cc1 = cc1 + cos(ff)*(t1 + t2)/4.d0/gp
          cc2 = 0.d0
        ENDDO

        IF (gamma_only) THEN
          cc1 = cc1*2d0
          cc2 = cc2*2d0
        ENDIF
        ew = ew + tt*(cc1 + cc2)
        IF (gstart == 2) ew = ew + tt*(kk1 + kk2)
      ENDDO
    ENDDO
    ewg = ewg + ew

    RETURN
  END SUBROUTINE esm_ewaldg_bc1

  SUBROUTINE esm_ewaldg_bc2(alpha_g, ewg)

    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewg
    !
    !    here the local variables
    !
    INTEGER  :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
                arg003, arg004, arg005, arg006, arg007, arg101, &
                arg102

    ewg = 0.d0
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    tmp = sqrt(alpha_g)
    sa = omega/L
    ew = 0d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa

        IF (gstart == 2) THEN
          IF (it1 == it2) THEN
            !! add coulomb energy of ions under efield
            ew = ew - zv(ityp(it1))*(z1 - z)*esm_efield/e2*2.0
          END IF
        END IF

        ! bc2
        arg001 = -tmp**2*(z - zp)**2
        arg101 = tmp*(z - zp)
        kk1 = 0.5d0*(-(z - zp)*erf(arg101) - EXP(arg001)/tmp/sqrt(pi))
        kk2 = 0.5d0*(z1 - z*zp/z1)

        cc1 = 0.d0
        cc2 = 0.d0
        DO ng_2d = 1, ngm_2d
          k1 = mill_2d(1, ng_2d)
          k2 = mill_2d(2, ng_2d)
          IF (k1 == 0 .and. k2 == 0) CYCLE
          t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
          gp2 = sum(t(:)*t(:))*tpiba2
          gp = sqrt(gp2)
          ff = ((k1*bg(1, 1) + k2*bg(1, 2))*(tau(1, it1) - tau(1, it2)) &
                + (k1*bg(2, 1) + k2*bg(2, 2))*(tau(2, it1) - tau(2, it2)))*tpi
          ! bc2
          arg001 = -gp*(z - zp)
          arg002 = gp*(z - zp)
          arg003 = -gp*(z + zp + 2.d0*z1)
          arg004 = gp*(z + zp - 2.d0*z1)
          arg005 = -gp*(z - zp + 4.d0*z1)
          arg006 = gp*(z - zp - 4.d0*z1)
          arg007 = -4.d0*gp*z1
          arg101 = gp/2.d0/tmp - tmp*(z - zp)
          arg102 = gp/2.d0/tmp + tmp*(z - zp)
          t1 = exp_erfc(arg001, arg101)
          t2 = exp_erfc(arg002, arg102)
          cc1 = cc1 + cos(ff)*(t1 + t2)/4.d0/gp
          cc2 = cc2 + cos(ff)*(EXP(arg006) + EXP(arg005) &
                               - EXP(arg004) - EXP(arg003)) &
                /(1.d0 - EXP(arg007))/2.d0/gp
        ENDDO

        IF (gamma_only) THEN
          cc1 = cc1*2d0
          cc2 = cc2*2d0
        ENDIF
        ew = ew + tt*(cc1 + cc2)
        IF (gstart == 2) ew = ew + tt*(kk1 + kk2)
      ENDDO
    ENDDO
    ewg = ewg + ew

    RETURN
  END SUBROUTINE esm_ewaldg_bc2

  SUBROUTINE esm_ewaldg_bc3(alpha_g, ewg)

    USE constants,        ONLY : pi, tpi, fpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewg
    !
    !    here the local variables
    !
    INTEGER  :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
                arg003, arg101, arg102

    ewg = 0.d0
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    tmp = sqrt(alpha_g)
    sa = omega/L
    ew = 0d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc3
        arg001 = -tmp**2*(z - zp)**2
        arg101 = tmp*(z - zp)
        kk1 = 0.5d0*(-(z - zp)*erf(arg101) - EXP(arg001)/tmp/sqrt(pi))
        kk2 = 0.5d0*(2.d0*z1 - z - zp)

        cc1 = 0.d0
        cc2 = 0.d0
        DO ng_2d = 1, ngm_2d
          k1 = mill_2d(1, ng_2d)
          k2 = mill_2d(2, ng_2d)
          IF (k1 == 0 .and. k2 == 0) CYCLE
          t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
          gp2 = sum(t(:)*t(:))*tpiba2
          gp = sqrt(gp2)
          ff = ((k1*bg(1, 1) + k2*bg(1, 2))*(tau(1, it1) - tau(1, it2)) &
                + (k1*bg(2, 1) + k2*bg(2, 2))*(tau(2, it1) - tau(2, it2)))*tpi
          ! bc3
          arg001 = -gp*(z - zp)
          arg002 = gp*(z - zp)
          arg003 = gp*(z + zp - 2.d0*z1)
          arg101 = gp/2.d0/tmp - tmp*(z - zp)
          arg102 = gp/2.d0/tmp + tmp*(z - zp)
          t1 = exp_erfc(arg001, arg101)
          t2 = exp_erfc(arg002, arg102)
          cc1 = cc1 + cos(ff)*(t1 + t2)/4.d0/gp
          cc2 = cc2 + cos(ff)*(-EXP(arg003))/2.d0/gp
        ENDDO

        IF (gamma_only) THEN
          cc1 = cc1*2d0
          cc2 = cc2*2d0
        ENDIF
        ew = ew + tt*(cc1 + cc2)
        IF (gstart == 2) ew = ew + tt*(kk1 + kk2)
      ENDDO
    ENDDO
    ewg = ewg + ew

    RETURN
  END SUBROUTINE esm_ewaldg_bc3

  SUBROUTINE esm_ewaldg_bc4(alpha_g, ewg)

    USE constants,        ONLY : pi, tpi, fpi
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, nsp, ityp, tau
    USE control_flags,    ONLY : gamma_only

    IMPLICIT NONE
    REAL(DP), INTENT(in)  :: alpha_g
    REAL(DP), INTENT(out) :: ewg
    !
    !    here the local variables
    !
    INTEGER  :: k1, k2, it1, it2, ng_2d
    REAL(DP) :: gp2, t(2), gp, sa, z, zp, z1, z0, L, t1, t2, tt, &
                tmp, cc1, cc2, kk1, kk2, ff, ew, arg001, arg002, &
                arg003, arg005, arg006, arg007, arg008, &
                arg009, arg011, arg101, arg102, arg103, arg104, &
                arg106, arg107, arg109, arg111, arg113, aaa, t3, &
                alpha, beta, kappa, lambda, xi, chi

    ewg = 0.d0
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    aaa = esm_a
    tmp = sqrt(alpha_g)
    sa = omega/L
    ew = 0d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        tt = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc4
        arg001 = -tmp**2*(z - zp)**2
        arg002 = -tmp**2*(z1 - zp)**2
        arg005 = -2.d0*aaa*(z - z1)
        arg006 = aaa**2/tmp**2 + 2.d0*aaa*(z1 - zp)
        arg101 = tmp*(z - zp)
        arg102 = tmp*(z1 - zp)
        arg104 = aaa/tmp + tmp*(z - zp)
        arg106 = aaa/tmp + tmp*(z1 - zp)
        IF (z < z1) THEN
          t1 = -(z - zp)*erf(arg101) + (0.5d0/aaa + z1 - zp)*erf(arg102)
          t2 = 0.5d0/aaa*exp_erfc(arg006, arg106)
          t3 = 0.5d0/aaa - (z - z1) + EXP(arg002)/tmp/sqrt(pi) &
               - EXP(arg001)/tmp/sqrt(pi)
          kk1 = (t1 + t2)/2.d0
          kk2 = t3/2.d0
        ELSE
          t1 = -exp_erfc(arg005, arg101)/aaa
          t2 = exp_erfc(arg006, arg104)/aaa
          t3 = EXP(arg005)/aaa
          kk1 = (t1 + t2)/4.d0
          kk2 = t3/2.d0
        ENDIF

        cc1 = 0.d0
        cc2 = 0.d0
        DO ng_2d = 1, ngm_2d
          k1 = mill_2d(1, ng_2d)
          k2 = mill_2d(2, ng_2d)
          IF (k1 == 0 .and. k2 == 0) CYCLE
          t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
          gp2 = sum(t(:)*t(:))*tpiba2
          gp = sqrt(gp2)
          ff = ((k1*bg(1, 1) + k2*bg(1, 2))*(tau(1, it1) - tau(1, it2)) &
                + (k1*bg(2, 1) + k2*bg(2, 2))*(tau(2, it1) - tau(2, it2)))*tpi
          ! bc4
          alpha = aaa + gp + sqrt(aaa**2 + gp**2)
          beta = aaa + gp - sqrt(aaa**2 + gp**2)
          kappa = aaa - gp + sqrt(aaa**2 + gp**2)
          xi = aaa + sqrt(aaa**2 + gp**2)
          chi = aaa - sqrt(aaa**2 + gp**2)
          lambda = sqrt(aaa**2 + gp**2)
          arg001 = gp*(z - zp)
          arg002 = -gp*(z - zp)
          arg003 = gp*(z + zp - 2.d0*z1)
          arg005 = -gp*(z1 - zp) - xi*(z - z1)
          arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
          arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
          arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
          arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
          arg101 = gp/2.d0/tmp + tmp*(z - zp)
          arg102 = gp/2.d0/tmp - tmp*(z - zp)
          arg103 = gp/2.d0/tmp + tmp*(z1 - zp)
          arg104 = gp/2.d0/tmp - tmp*(z1 - zp)
          arg107 = xi/2.d0/tmp + tmp*(z - zp)
          arg109 = xi/2.d0/tmp + tmp*(z1 - zp)
          arg111 = chi/2.d0/tmp + tmp*(z - zp)
          arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
          IF (z < z1) THEN
            t1 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103)
            t2 = exp_erfc(arg002, arg102) &
                 - kappa/alpha*exp_erfc(arg003, arg104)
            t3 = exp_erfc(arg006, arg109)/alpha
            cc1 = cc1 + cos(ff)*(t1 + t2)/4.d0/gp
            cc2 = cc2 + cos(ff)*t3/2.d0
          ELSE
            t1 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111)
            t2 = exp_erfc(arg008, arg107) &
                 - beta/alpha*exp_erfc(arg009, arg109)
            t3 = exp_erfc(arg005, arg104)/alpha
            cc1 = cc1 + cos(ff)*(t1 + t2)/4.d0/lambda
            cc2 = cc2 + cos(ff)*t3/2.d0
          ENDIF
        ENDDO

        IF (gamma_only) THEN
          cc1 = cc1*2d0
          cc2 = cc2*2d0
        ENDIF
        ew = ew + tt*(cc1 + cc2)
        IF (gstart == 2) ew = ew + tt*(kk1 + kk2)
      ENDDO
    ENDDO
    ewg = ewg + ew

    RETURN
  END SUBROUTINE esm_ewaldg_bc4

END MODULE esm_ewald_mod
