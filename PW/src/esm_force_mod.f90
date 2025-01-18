MODULE esm_force_mod

  USE kinds,    ONLY : DP
  USE esm_common_mod, ONLY : esm_efield, esm_w, esm_a, esm_bc, &
                             mill_2d, imill_2d, ngm_2d, &
                             vl11j0, vl11j1, vl12j0, vl12j1, vl21j1, vl22j1, vl11, vl22, dvl11, dvl22, &
                             dvl11j0, dvl12j0, dvl21j0, dvl22j0, &
                             qromb, exp_erfc, esm_rgen_2d

  IMPLICIT NONE

CONTAINS

  SUBROUTINE esm_force_ew(forceion)
    !-----------------------------------------------------------------------
    !
    !  This routine computes the Ewald contribution to the forces,
    !  both the real- and reciprocal-space terms are present
    !
    USE kinds
    USE constants, ONLY : tpi, e2
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    USE ions_base, ONLY : zv, nat, ityp
    USE gvect,     ONLY : gcutm
    USE cell_base, ONLY : tpiba2
    IMPLICIT NONE

    REAL(DP), INTENT(out) :: forceion(3, nat)
    ! output: the ewald part of the forces
    !
    REAL(DP) :: alpha, charge, upperbound
    ! the alpha parameter
    ! the total charge
    ! used to determine alpha

    forceion(:, :) = 0.d0
    charge = sum(zv(ityp(1:nat)))
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error ON THE ENERGY
    !
    alpha = 2.9d0
    DO
      alpha = alpha - 0.1d0
      IF (alpha == 0.d0) THEN
        CALL errore('esm_force_ew', 'optimal alpha not found', 1)
      END IF
      upperbound = e2*charge**2*sqrt(2.d0*alpha/tpi)* &
                   erfc(sqrt(tpiba2*gcutm/4.d0/alpha))
      IF (upperbound < 1.0d-7) EXIT
    END DO
    !write(*,'(5X,A,F5.2)')'alpha used in esm ewald force :',alpha

    CALL esm_force_ewg(alpha, forceion)

    CALL esm_force_ewr(alpha, forceion)

    CALL mp_sum(forceion, intra_bgrp_comm)

    RETURN
  END SUBROUTINE esm_force_ew

  !-----------------------------------------------------------------------
  !--------------ESM FORCE SUBROUTINE-------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE esm_force_ewr(alpha_g, forceion)
    USE kinds,     ONLY : DP
    USE ions_base, ONLY : nat
    IMPLICIT NONE
    REAL(DP),INTENT(in)    :: alpha_g
    REAL(DP),INTENT(inout) :: forceion(3,nat)

    IF (esm_bc == 'pbc') THEN
      CALL esm_force_ewr_pbc(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_force_ewr_pbc(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_force_ewr_pbc(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_force_ewr_pbc(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_force_ewr_bc4(alpha_g, forceion)
    END IF

  END SUBROUTINE esm_force_ewr

  SUBROUTINE esm_force_ewg(alpha_g, forceion)
    USE kinds,     ONLY : DP
    USE ions_base, ONLY : nat
    IMPLICIT NONE
    REAL(DP), INTENT(in)    :: alpha_g
    REAL(DP), INTENT(out)   :: forceion(3, nat)

    IF (esm_bc == 'pbc') THEN
      CALL esm_force_ewg_pbc(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_force_ewg_bc1(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_force_ewg_bc2(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_force_ewg_bc3(alpha_g, forceion)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_force_ewg_bc4(alpha_g, forceion)
    END IF

  END SUBROUTINE esm_force_ewg

  SUBROUTINE esm_force_lc(aux, forcelc)
    USE kinds, ONLY: DP
    USE ions_base, ONLY: nat
    USE fft_base, ONLY: dfftp
    IMPLICIT NONE
    COMPLEX(DP), INTENT(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
    REAL(DP), INTENT(inout) :: forcelc(3, nat)

    IF (esm_bc == 'pbc') THEN
      CALL esm_force_lc_pbc(aux, forcelc)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_force_lc_bc1(aux, forcelc)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_force_lc_bc2(aux, forcelc)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_force_lc_bc3(aux, forcelc)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_force_lc_bc4(aux, forcelc)
    END IF

  END SUBROUTINE esm_force_lc

!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE (RSUM) SUBROUTINE ---------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_force_ewr_pbc(alpha_g, forceion)
    USE constants,        ONLY : pi, e2
    USE cell_base,        ONLY : alat, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm

    IMPLICIT NONE
    INTEGER               :: na, nb, nr, nrm, ip, np
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER    :: mxr = 1000
    ! the maximum number of R vectors included in r
    REAL(DP)              :: dtau(3), r(3, mxr), r2(mxr)
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    REAL(DP), INTENT(in)   :: alpha_g
    REAL(DP), INTENT(inout):: forceion(3, nat)
    !
    ! ESM variables
    !
    REAL(DP)              :: tmp, fac, rmax0, rr
    ! rmax0: the maximum radius to consider real space sum
    REAL(DP), ALLOCATABLE :: force(:, :)

    tmp = sqrt(alpha_g)
    rmax0 = 5.d0/tmp/alat

    ip = mp_rank(intra_bgrp_comm)
    np = mp_size(intra_bgrp_comm)

    ALLOCATE (force(3, nat))
    force(:, :) = 0.d0
    DO na = ip + 1, nat, np
      DO nb = 1, nat
        IF (nb .eq. na) CYCLE
        dtau(:) = tau(:, na) - tau(:, nb)
        fac = zv(ityp(na))*zv(ityp(nb))*e2
        !
        ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
        !
        CALL rgen(dtau, rmax0, mxr, at, bg, r, r2, nrm)
        !
        ! and sum to the real space part
        !
        DO nr = 1, nrm
          rr = sqrt(r2(nr))*alat
          force(:, na) = force(:, na) &
                         - fac/rr**2*(erfc(tmp*rr)/rr + 2.d0*tmp/sqrt(pi) &
                                      *EXP(-tmp**2*rr**2))*r(:, nr)*alat
        ENDDO
      ENDDO
    ENDDO
    forceion(:, :) = forceion(:, :) + force(:, :)
    DEALLOCATE (force)

  END SUBROUTINE esm_force_ewr_pbc

  SUBROUTINE esm_force_ewr_bc4(alpha_g, forceion)
    USE io_global,        ONLY : stdout
    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : gstart
    USE cell_base,        ONLY : alat, tpiba2, at, bg
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE control_flags,    ONLY : iverbosity
    USE mp,               ONLY : mp_rank, mp_size
    USE mp_bands,         ONLY : intra_bgrp_comm

    IMPLICIT NONE
    INTEGER               :: na, nb, nr, nrm, ipol, ip, np
    ! counter on atoms
    ! counter on atoms
    ! counter over direct vectors
    ! number of R vectors included in r sum
    INTEGER, PARAMETER    :: mxr = 1000
    ! the maximum number of R vectors included in r
    REAL(DP)              :: dtau(3), r(3, mxr), r2(mxr), rxy, rxyz
    ! the difference tau_s - tau_s'
    ! neighbering shell vector
    ! the square modulus of R_j-tau_s-tau_s'
    ! buffer variable
    ! buffer variable
    REAL(DP), INTENT(in)   :: alpha_g
    REAL(DP), INTENT(inout):: forceion(3, nat)
    !
    ! ESM variables
    !
    REAL(DP)              :: L, z, zp, z0, z1, aaa, tmp, ss, fac, err, ss0, &
                             gpmax, rmax0, rmax, zbuff, znrm, rr
    ! gpmax: upper bound of g_parallel integral
    ! rmax: the maximum radius to consider real space sum
    ! zbuff: smearing width to avoid the singularity of the Force
    ! znrm: threashold value for normal RSUM and Smooth-ESM's RSUM
    REAL(DP), PARAMETER :: eps = 1.d-11, epsneib = 1.d-6
    REAL(DP), ALLOCATABLE :: force(:, :)

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
        CALL errore('esm_force_ewr', 'optimal gpmax not found', 1)
      CALL qromb(vl11, aaa, tmp, z1, z1 - zbuff, z1 - zbuff, 0.0_DP, gpmax, ss)
      err = ABS(ss - ss0); ss0 = ss
      IF (err .lt. eps) EXIT
    ENDDO
    ! Define znrm using the deviation from the constant term in RSUM
    znrm = z1
    DO
      znrm = znrm - 0.01d0
      IF (znrm .le. -z0) &
        CALL errore('esm_force_ewr', 'optimal znrm not found', 1)
      CALL qromb(vl11, aaa, tmp, z1, znrm, znrm, 0.0_DP, gpmax, ss)
      err = -2.d0*tmp/sqrt(pi) - ss*2.d0
      IF (ABS(err) .lt. eps) EXIT
    ENDDO
    ! Define rmax for real space sum
    rmax = 1.d0
    DO
      rmax = rmax + 1.d0
      IF (rmax .gt. 200.d0) &
        CALL errore('esm_force_ewr', 'optimal rmax not found', 1)
      CALL qromb(dvl11j0, aaa, tmp, z1, z1 - zbuff, z1 - zbuff, rmax, gpmax, ss)
      err = ss
      IF (ABS(err) .lt. epsneib) EXIT
    ENDDO
    rmax = rmax/alat
    IF (iverbosity > 0) THEN
      write (stdout, '(5x,"=== Smooth-ESM RSUM parameters (Force) ===")')
      write (stdout, '(5x,A,F10.2,A)') &
        'Upper bound of g_parallel integral:      ', gpmax, ' (1/a.u.)'
      write (stdout, '(5x,A,F10.2,A)') &
        'Boundary for normal RSUM|Smooth-ESM RSUM:', z1 - znrm, ' (a.u.)'
      write (stdout, '(5x,A,F10.2,A)') &
        'Upper bound of real-space summation:     ', rmax*alat, ' (a.u.)'
      write (stdout, '(5x,"==========================================")')
    ENDIF
    !
    ip = mp_rank(intra_bgrp_comm)
    np = mp_size(intra_bgrp_comm)

    ALLOCATE (force(3, nat))
    force(:, :) = 0.d0
    DO na = ip + 1, nat, np
      z = tau(3, na)
      IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
      z = z*alat
      DO nb = 1, nat
        IF (nb .eq. na) CYCLE
        zp = tau(3, nb)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        dtau(1:2) = tau(1:2, na) - tau(1:2, nb)
        dtau(3) = (z - zp)/alat
        fac = zv(ityp(na))*zv(ityp(nb))*e2
        IF (z < znrm) THEN
          IF (zp < znrm) THEN ! z in I, zp in I (normal RSUM)
            rmax0 = 5.d0/tmp/alat
            !
            ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
            !
            CALL rgen(dtau, rmax0, mxr, at, bg, r, r2, nrm)
            !
            ! and sum to the real space part
            !
            DO nr = 1, nrm
              rr = sqrt(r2(nr))*alat
              DO ipol = 1, 3
                force(ipol, na) = force(ipol, na) &
                                  - fac/rr**2*(erfc(tmp*rr)/rr + 2.d0*tmp/sqrt(pi) &
                                  *EXP(-tmp**2*rr**2))*r(ipol, nr)*alat
              ENDDO
            ENDDO
          ELSEIF (zp < z1) THEN ! z in I, zp in I
            CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              rxyz = sqrt(r2(nr) + dtau(3)**2)*alat
              CALL qromb(vl11j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(1:2, nb) = force(1:2, nb) &
                               - fac*(1.d0/rxyz**3 + 1.d0/rxy*ss)*r(1:2, nr)*alat
              CALL qromb(dvl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(3, nb) = force(3, nb) - fac*((z - zp)/rxyz**3 + ss)
            ENDDO
          ELSE ! z in I, zp in II
            CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              CALL qromb(vl12j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(1:2, nb) = force(1:2, nb) &
                               - fac*ss/rxy*r(1:2, nr)*alat
              CALL qromb(dvl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(3, nb) = force(3, nb) - fac*ss
            ENDDO
          ENDIF ! IF for zp
        ELSEIF (z < z1) THEN ! znrm < z < z1
          CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
          IF (zp < z1) THEN ! z in I, zp in I
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              rxyz = sqrt(r2(nr) + dtau(3)**2)*alat
              CALL qromb(vl11j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(1:2, nb) = force(1:2, nb) &
                               - fac*(1.d0/rxyz**3 + 1.d0/rxy*ss)*r(1:2, nr)*alat
              CALL qromb(dvl11j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(3, nb) = force(3, nb) - fac*((z - zp)/rxyz**3 + ss)
            ENDDO
          ELSE ! z in I, zp in II
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              CALL qromb(vl12j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(1:2, nb) = force(1:2, nb) &
                               - fac*ss/rxy*r(1:2, nr)*alat
              CALL qromb(dvl12j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(3, nb) = force(3, nb) - fac*ss
            ENDDO
          ENDIF ! IF for zp
        ELSE ! z1 < z
          CALL esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
          IF (zp < z1) THEN ! z in II, zp in I
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              CALL qromb(vl21j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(1:2, nb) = force(1:2, nb) &
                               - fac*ss/rxy*r(1:2, nr)*alat
              CALL qromb(dvl21j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(3, nb) = force(3, nb) - fac*ss
            ENDDO
          ELSE ! z in II, zp in II
            DO nr = 1, nrm
              rxy = sqrt(r2(nr))*alat
              rxyz = sqrt(r2(nr) + dtau(3)**2)*alat
              CALL qromb(vl22j1, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(1:2, nb) = force(1:2, nb) &
                               - (EXP(-aaa*(rxyz + z + zp - 2.d0*z1))*(aaa + 1.d0/rxyz)/rxyz**2 &
                               + ss/rxy)*fac*r(1:2, nr)*alat
              CALL qromb(dvl22j0, aaa, tmp, z1, z, zp, rxy, gpmax, ss)
              force(3, nb) = force(3, nb) &
                             - (EXP(-aaa*(rxyz + z + zp - 2.d0*z1))*(aaa + 1.d0/rxyz)/rxyz**2 &
                             *(z - zp) - aaa*EXP(-aaa*(rxyz + z + zp - 2.d0*z1))/rxyz + ss)*fac
            ENDDO
          ENDIF ! IF for zp
        ENDIF
      ENDDO
      IF (z < znrm) THEN
        ss = 0.d0
      ELSEIF (z < z1) THEN
        CALL qromb(dvl11, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss)
      ELSE
        CALL qromb(dvl22, aaa, tmp, z1, z, z, 0.0_DP, gpmax, ss)
      ENDIF
      ! factor e2: hartree -> Ry.
      force(3, na) = force(3, na) - zv(ityp(na))**2*e2*ss
    ENDDO
    forceion(:, :) = forceion(:, :) + force(:, :)
    DEALLOCATE (force)

  END SUBROUTINE esm_force_ewr_bc4

!-----------------------------------------------------------------------
!--------------ESM EWALD-DERIVED FORCE (GSUM) SUBROUTINE ---------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_force_ewg_pbc(alpha_g, forceion)

    USE constants,        ONLY : tpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, gg, g
    USE vlocal,           ONLY : strf

    IMPLICIT NONE
    REAL(DP), INTENT(in)     :: alpha_g
    REAL(DP), INTENT(out)    :: forceion(3, nat)
    INTEGER                  :: nt, ig, na, ipol
    REAL(DP)                 :: fact, arg, sumnb
    COMPLEX(DP), ALLOCATABLE :: aux(:)

    forceion(:, :) = 0.d0

    ! same of the GSUM part in force_ew.f90
    ALLOCATE (aux(ngm))
    aux(:) = (0.d0, 0.d0)

    DO nt = 1, nsp
      DO ig = gstart, ngm
        aux(ig) = aux(ig) + zv(nt)*CONJG(strf(ig, nt))
      ENDDO
    ENDDO
    DO ig = gstart, ngm
      aux(ig) = aux(ig)*EXP(-gg(ig)*tpiba2/alpha_g/4.d0) &
                /(gg(ig)*tpiba2)
    ENDDO

    IF (gamma_only) THEN
      fact = 4.d0
    ELSE
      fact = 2.d0
    END IF
    DO na = 1, nat
      DO ig = gstart, ngm
        arg = tpi*(g(1, ig)*tau(1, na) + g(2, ig)*tau(2, na) &
                   + g(3, ig)*tau(3, na))
        sumnb = cos(arg)*AIMAG(aux(ig)) - sin(arg)*DBLE(aux(ig))
        forceion(1, na) = forceion(1, na) + g(1, ig)*sumnb
        forceion(2, na) = forceion(2, na) + g(2, ig)*sumnb
        forceion(3, na) = forceion(3, na) + g(3, ig)*sumnb
      ENDDO
      DO ipol = 1, 3
        forceion(ipol, na) = -zv(ityp(na))*fact*e2*tpi**2/ &
                             omega/alat*forceion(ipol, na)
      ENDDO
    ENDDO
    DEALLOCATE (aux)

    RETURN
  END SUBROUTINE esm_force_ewg_pbc

  SUBROUTINE esm_force_ewg_bc1(alpha_g, forceion)

    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g

    IMPLICIT NONE
    REAL(DP), INTENT(in)    :: alpha_g
    REAL(DP), INTENT(out)   :: forceion(3, nat)
    !
    !    here the local variables
    !
    INTEGER                 :: it1, it2, k1, k2, ng_2d
    REAL(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                               c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                               ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                               arg001, arg002, arg101, arg102

    forceion(:, :) = 0.d0
    for_g(:, :) = 0.d0
    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    tmp = sqrt(alpha_g)

    for = 0.d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        IF (gamma_only) THEN
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
        ELSE
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ENDIF
        t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc1
        kk1_for = 0.5d0*erf(tmp*(z - zp))
        kk2_for = 0.d0

        c1_for(:) = 0.d0; c2_for(:) = 0.d0
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
          c1_for(1) = c1_for(1) + sin(ff)*(t1 + t2)/4.d0/gp*k1
          c1_for(2) = c1_for(2) + sin(ff)*(t1 + t2)/4.d0/gp*k2
          c1_for(3) = c1_for(3) + cos(ff)*(t1 - t2)/4.d0
        ENDDO
        for(:, it2) = for(:, it2) + t1_for*(c1_for(:) + c2_for(:))
        IF (gstart == 2) THEN
          for(3, it2) = for(3, it2) + t2_for*(kk1_for + kk2_for)
        ENDIF

      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)

    for_g(:, :) = for_g(:, :)*e2 ! factor e2: hartree -> Ry.

    DO it1 = 1, nat
      forceion(1, it1) = -sum(for_g(1:2, it1)*bg(1, 1:2))*sqrt(tpiba2)
      forceion(2, it1) = -sum(for_g(1:2, it1)*bg(2, 1:2))*sqrt(tpiba2)
      forceion(3, it1) = -for_g(3, it1)
    ENDDO

    RETURN
  END SUBROUTINE esm_force_ewg_bc1

  SUBROUTINE esm_force_ewg_bc2(alpha_g, forceion)

    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g

    IMPLICIT NONE
    REAL(DP), INTENT(in)    :: alpha_g
    REAL(DP), INTENT(out)   :: forceion(3, nat)
    !
    !    here the local variables
    !
    INTEGER                 :: it1, it2, k1, k2, ng_2d
    REAL(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                               c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                               ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                               arg001, arg002, arg003, arg004, arg005, &
                               arg006, arg007, arg101, arg102

    forceion(:, :) = 0.d0
    for_g(:, :) = 0.d0
    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    tmp = sqrt(alpha_g)

    for = 0.d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        IF (gamma_only) THEN
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
        ELSE
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ENDIF
        t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc2
        kk1_for = 0.5d0*erf(tmp*(z - zp))
        kk2_for = -0.5d0*(z/z1)

        c1_for(:) = 0.d0; c2_for(:) = 0.d0
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
          c1_for(1) = c1_for(1) + sin(ff)*(t1 + t2)/4.d0/gp*k1
          c1_for(2) = c1_for(2) + sin(ff)*(t1 + t2)/4.d0/gp*k2
          c1_for(3) = c1_for(3) + cos(ff)*(t1 - t2)/4.d0
          c2_for(1) = c2_for(1) + sin(ff)*(EXP(arg006) + EXP(arg005) &
                      - EXP(arg004) - EXP(arg003))/(1.d0 - EXP(arg007))/2.d0/gp*k1
          c2_for(2) = c2_for(2) + sin(ff)*(EXP(arg006) + EXP(arg005) &
                      - EXP(arg004) - EXP(arg003))/(1.d0 - EXP(arg007))/2.d0/gp*k2
          c2_for(3) = c2_for(3) - cos(ff)*(EXP(arg006) - EXP(arg005) &
                      + EXP(arg004) - EXP(arg003))/(1.d0 - EXP(arg007))/2.d0
        ENDDO
        for(:, it2) = for(:, it2) + t1_for*(c1_for(:) + c2_for(:))
        IF (gstart == 2) THEN
          for(3, it2) = for(3, it2) + t2_for*(kk1_for + kk2_for)
        ENDIF

      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)

    for_g(:, :) = for_g(:, :)*e2 ! factor e2: hartree -> Ry.

    DO it1 = 1, nat
      forceion(1, it1) = -sum(for_g(1:2, it1)*bg(1, 1:2))*sqrt(tpiba2)
      forceion(2, it1) = -sum(for_g(1:2, it1)*bg(2, 1:2))*sqrt(tpiba2)
      forceion(3, it1) = -for_g(3, it1)
      IF (gstart == 2) THEN
        !! add coulomb fource of ions under efield
        forceion(3, it1) = forceion(3, it1) - zv(ityp(it1))*esm_efield
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE esm_force_ewg_bc2

  SUBROUTINE esm_force_ewg_bc3(alpha_g, forceion)

    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g

    IMPLICIT NONE
    REAL(DP), INTENT(in)    :: alpha_g
    REAL(DP), INTENT(out)   :: forceion(3, nat)
    !
    !    here the local variables
    !
    INTEGER                 :: it1, it2, k1, k2, ng_2d
    REAL(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                               c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                               ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                               arg001, arg002, arg003, arg101, arg102

    forceion(:, :) = 0.d0
    for_g(:, :) = 0.d0
    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    tmp = sqrt(alpha_g)

    for = 0.d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        IF (gamma_only) THEN
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
        ELSE
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ENDIF
        t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc3
        kk1_for = 0.5d0*erf(tmp*(z - zp))
        kk2_for = -0.5d0

        c1_for(:) = 0.d0; c2_for(:) = 0.d0
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
          c1_for(1) = c1_for(1) + sin(ff)*(t1 + t2)/4.d0/gp*k1
          c1_for(2) = c1_for(2) + sin(ff)*(t1 + t2)/4.d0/gp*k2
          c1_for(3) = c1_for(3) + cos(ff)*(t1 - t2)/4.d0
          c2_for(1) = c2_for(1) + sin(ff)*(-EXP(arg003))/2.d0/gp*k1
          c2_for(2) = c2_for(2) + sin(ff)*(-EXP(arg003))/2.d0/gp*k2
          c2_for(3) = c2_for(3) + cos(ff)*(-EXP(arg003))/2.d0
        ENDDO
        for(:, it2) = for(:, it2) + t1_for*(c1_for(:) + c2_for(:))
        IF (gstart == 2) THEN
          for(3, it2) = for(3, it2) + t2_for*(kk1_for + kk2_for)
        ENDIF

      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)

    for_g(:, :) = for_g(:, :)*e2 ! factor e2: hartree -> Ry.

    DO it1 = 1, nat
      forceion(1, it1) = -sum(for_g(1:2, it1)*bg(1, 1:2))*sqrt(tpiba2)
      forceion(2, it1) = -sum(for_g(1:2, it1)*bg(2, 1:2))*sqrt(tpiba2)
      forceion(3, it1) = -for_g(3, it1)
    ENDDO

    RETURN
  END SUBROUTINE esm_force_ewg_bc3

  SUBROUTINE esm_force_ewg_bc4(alpha_g, forceion)

    USE constants,        ONLY : tpi, fpi, e2
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, nsp, tau, ityp
    USE gvect,            ONLY : gstart, ngm, g

    IMPLICIT NONE
    REAL(DP), INTENT(in)    :: alpha_g
    REAL(DP), INTENT(out)   :: forceion(3, nat)
    !
    !    here the local variables
    !
    INTEGER                 :: it1, it2, k1, k2, ng_2d
    REAL(DP)                :: for(3, nat), for_g(3, nat), t1_for, t2_for, &
                               c1_for(3), c2_for(3), kk1_for, kk2_for, t1, t2, &
                               ff, z0, z1, z, zp, tmp, gp2, gp, t(2), L, sa, &
                               arg001, arg002, arg003, arg004, arg005, &
                               arg006, arg007, arg008, arg009, arg010, &
                               arg011, arg012, arg101, arg102, arg103, &
                               arg104, arg105, arg106, arg107, arg108, &
                               arg109, arg110, arg111, arg112, arg113, &
                               arg114, aaa, t3, alpha, beta, kappa, lambda, &
                               xi, chi
    ! auxiliary space

    forceion(:, :) = 0.d0
    for_g(:, :) = 0.d0
    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    aaa = esm_a
    z1 = z0 + esm_w
    tmp = sqrt(alpha_g)

    for = 0.d0
    DO it1 = 1, nat
      DO it2 = 1, nat
        z = tau(3, it1)
        IF (z .gt. at(3, 3)*0.5) z = z - at(3, 3)
        z = z*alat
        zp = tau(3, it2)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        IF (gamma_only) THEN
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa*2.d0
        ELSE
          t1_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ENDIF
        t2_for = zv(ityp(it1))*zv(ityp(it2))*fpi/sa
        ! bc4
        arg004 = -2.d0*aaa*(zp - z1)
        arg006 = aaa**2/tmp**2 + 2.d0*aaa*(z1 - zp)
        arg101 = tmp*(z - zp)
        arg102 = tmp*(z1 - zp)
        arg104 = aaa/tmp + tmp*(z - zp)
        arg106 = aaa/tmp + tmp*(z1 - zp)
        IF (z < z1) THEN  ! factor 1/2 <- non-reciprocality
          IF (zp < z1) THEN
            kk1_for = 0.5d0*(erf(arg101) - erf(arg102))/2.d0 &
                      - 0.5d0*exp_erfc(arg006, arg106)/2.d0
            kk2_for = -0.5d0*erfc(arg101)/2.d0
          ELSE
            kk1_for = 0.5d0*(erf(arg101) - erf(arg102))/2.d0 &
                      - 0.5d0*exp_erfc(arg006, arg106)/2.d0
            kk2_for = -0.5d0*exp_erfc(arg004, arg101)/2.d0
          ENDIF
        ELSE
          IF (zp < z1) THEN
            kk1_for = -0.5d0*exp_erfc(arg006, arg104)/2.d0
            kk2_for = -0.5d0*erfc(arg101)/2.d0
          ELSE
            kk1_for = -0.5d0*exp_erfc(arg006, arg104)/2.d0
            kk2_for = -0.5d0*exp_erfc(arg004, arg101)/2.d0
          ENDIF
        ENDIF
        c1_for(:) = 0.d0; c2_for(:) = 0.d0
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
          arg004 = gp*(z - z1) + xi*(z1 - zp)
          arg005 = -gp*(z1 - zp) - xi*(z - z1)
          arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
          arg007 = aaa/2.d0/tmp**2*xi - gp*(z1 - zp) - xi*(z - z1)
          arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
          arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
          arg010 = aaa/2.d0/tmp**2*xi + chi*(z1 - zp) - xi*(z - z1)
          arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
          arg012 = aaa/2.d0/tmp**2*chi + xi*(z1 - zp) - chi*(z - z1)
          arg101 = gp/2.d0/tmp + tmp*(z - zp)
          arg102 = gp/2.d0/tmp - tmp*(z - zp)
          arg103 = gp/2.d0/tmp + tmp*(z1 - zp)
          arg104 = gp/2.d0/tmp - tmp*(z1 - zp)
          arg105 = gp/2.d0/tmp + tmp*(z - z1)
          arg106 = gp/2.d0/tmp - tmp*(z - z1)
          arg107 = xi/2.d0/tmp + tmp*(z - zp)
          arg108 = xi/2.d0/tmp - tmp*(z - zp)
          arg109 = xi/2.d0/tmp + tmp*(z1 - zp)
          arg110 = xi/2.d0/tmp - tmp*(z - z1)
          arg111 = chi/2.d0/tmp + tmp*(z - zp)
          arg112 = chi/2.d0/tmp - tmp*(z - zp)
          arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
          arg114 = chi/2.d0/tmp - tmp*(z - z1)
          IF (z < z1) THEN ! factor 1/2 <- non-reciprocality
            IF (zp < z1) THEN
              t1 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103)
              t2 = exp_erfc(arg002, arg102) &
                   - kappa/alpha*exp_erfc(arg003, arg104)
              t3 = exp_erfc(arg006, arg109)/alpha
              c1_for(1) = c1_for(1) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k1/2.d0
              c1_for(2) = c1_for(2) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k2/2.d0
              t1 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106)
              t2 = exp_erfc(arg001, arg101) &
                   - kappa/alpha*exp_erfc(arg003, arg105)
              t3 = exp_erfc(arg007, arg110)/alpha
              c2_for(1) = c2_for(1) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k1/2.d0
              c2_for(2) = c2_for(2) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k2/2.d0
              t1 = exp_erfc(arg001, arg103) - exp_erfc(arg001, arg101)
              t2 = exp_erfc(arg002, arg102) &
                   - kappa/alpha*exp_erfc(arg003, arg104)
              t3 = -xi/alpha*exp_erfc(arg006, arg109)
              c1_for(3) = c1_for(3) + cos(ff)*((t1 + t2)/4.d0 + t3/2.d0)/2.d0
              t1 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106)
              t2 = -exp_erfc(arg001, arg101) &
                   - kappa/alpha*exp_erfc(arg003, arg105)
              t3 = gp/alpha*exp_erfc(arg007, arg110)
              c2_for(3) = c2_for(3) + cos(ff)*((t1 + t2)/4.d0 + t3/2.d0)/2.d0
            ELSE
              t1 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103)
              t2 = exp_erfc(arg002, arg102) &
                   - kappa/alpha*exp_erfc(arg003, arg104)
              t3 = exp_erfc(arg006, arg109)/alpha
              c1_for(1) = c1_for(1) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k1/2.d0
              c1_for(2) = c1_for(2) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k2/2.d0
              t1 = exp_erfc(arg012, arg114) - exp_erfc(arg012, arg112)
              t2 = exp_erfc(arg010, arg108) &
                   - beta/alpha*exp_erfc(arg009, arg110)
              t3 = exp_erfc(arg004, arg105)/alpha
              c2_for(1) = c2_for(1) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k1/2.d0
              c2_for(2) = c2_for(2) + sin(ff)*((t1 + t2)/4.d0/gp + t3/2.d0)*k2/2.d0
              t1 = exp_erfc(arg001, arg103) - exp_erfc(arg001, arg101)
              t2 = exp_erfc(arg002, arg102) &
                   - kappa/alpha*exp_erfc(arg003, arg104)
              t3 = -xi/alpha*exp_erfc(arg006, arg109)
              c1_for(3) = c1_for(3) + cos(ff)*((t1 + t2)/4.d0 + t3/2.d0)/2.d0
              t1 = xi*(exp_erfc(arg012, arg112) - exp_erfc(arg012, arg114))
              t2 = -chi*exp_erfc(arg010, arg108) &
                   + xi*beta/alpha*exp_erfc(arg009, arg110)
              t3 = -xi/alpha*exp_erfc(arg004, arg105)
              c2_for(3) = c2_for(3) + cos(ff)*((t1 + t2)/4.d0 + t3/2.d0)/2.d0
            ENDIF
          ELSE
            IF (zp < z1) THEN
              t1 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111)
              t2 = exp_erfc(arg008, arg107) &
                   - beta/alpha*exp_erfc(arg009, arg109)
              t3 = exp_erfc(arg005, arg104)/alpha
              c1_for(1) = c1_for(1) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k1/2.d0
              c1_for(2) = c1_for(2) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k2/2.d0
              t1 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106)
              t2 = exp_erfc(arg001, arg101) &
                   - kappa/alpha*exp_erfc(arg003, arg105)
              t3 = exp_erfc(arg007, arg110)/alpha
              c2_for(1) = c2_for(1) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k1/2.d0
              c2_for(2) = c2_for(2) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k2/2.d0
              t1 = chi*(exp_erfc(arg011, arg111) - exp_erfc(arg011, arg113))
              t2 = -xi*exp_erfc(arg008, arg107) &
                   + xi*beta/alpha*exp_erfc(arg009, arg109)
              t3 = gp/alpha*exp_erfc(arg005, arg104)
              c1_for(3) = c1_for(3) + cos(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)/2.d0
              t1 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106)
              t2 = -exp_erfc(arg001, arg101) &
                   - kappa/alpha*exp_erfc(arg003, arg105)
              t3 = gp/alpha*exp_erfc(arg007, arg110)
              c2_for(3) = c2_for(3) + cos(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)/2.d0
            ELSE
              t1 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111)
              t2 = exp_erfc(arg008, arg107) &
                   - beta/alpha*exp_erfc(arg009, arg109)
              t3 = exp_erfc(arg005, arg104)/alpha
              c1_for(1) = c1_for(1) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k1/2.d0
              c1_for(2) = c1_for(2) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k2/2.d0
              t1 = exp_erfc(arg012, arg114) - exp_erfc(arg012, arg112)
              t2 = exp_erfc(arg010, arg108) &
                   - beta/alpha*exp_erfc(arg009, arg110)
              t3 = exp_erfc(arg004, arg105)/alpha
              c2_for(1) = c2_for(1) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k1/2.d0
              c2_for(2) = c2_for(2) + sin(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)*k2/2.d0
              t1 = chi*(exp_erfc(arg011, arg111) - exp_erfc(arg011, arg113))
              t2 = -xi*exp_erfc(arg008, arg107) &
                   + xi*beta/alpha*exp_erfc(arg009, arg109)
              t3 = gp/alpha*exp_erfc(arg005, arg104)
              c1_for(3) = c1_for(3) + cos(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)/2.d0
              t1 = xi*(exp_erfc(arg012, arg112) - exp_erfc(arg012, arg114))
              t2 = -chi*exp_erfc(arg010, arg108) &
                   + xi*beta/alpha*exp_erfc(arg009, arg110)
              t3 = -xi/alpha*exp_erfc(arg004, arg105)
              c2_for(3) = c2_for(3) + cos(ff)*((t1 + t2)/4.d0/lambda + t3/2.d0)/2.d0
            ENDIF
          ENDIF
        ENDDO
        for(:, it2) = for(:, it2) + t1_for*(c1_for(:) + c2_for(:))
        IF (gstart == 2) THEN
          for(3, it2) = for(3, it2) + t2_for*(kk1_for + kk2_for)
        ENDIF

      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)

    for_g(:, :) = for_g(:, :)*e2 ! factor e2: hartree -> Ry.

    DO it1 = 1, nat
      forceion(1, it1) = -sum(for_g(1:2, it1)*bg(1, 1:2))*sqrt(tpiba2)
      forceion(2, it1) = -sum(for_g(1:2, it1)*bg(2, 1:2))*sqrt(tpiba2)
      forceion(3, it1) = -for_g(3, it1)
    ENDDO

    RETURN
  END SUBROUTINE esm_force_ewg_bc4

!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL-DERIVED FORCE SUBROUTINE-------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_force_lc_pbc(aux, forcelc)
    USE ions_base, ONLY: nat
    USE fft_base,  ONLY: dfftp
    IMPLICIT NONE
    COMPLEX(DP), INTENT(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
    REAL(DP),    INTENT(inout) :: forcelc(3, nat)

    STOP 'esm_force_lc must not be called for esm_bc = pbc'

  END SUBROUTINE esm_force_lc_pbc

  SUBROUTINE esm_force_lc_bc1(aux, forcelc)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z

    IMPLICIT NONE
    COMPLEX(DP), INTENT(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
    REAL(DP),    INTENT(inout) :: forcelc(3, nat)
    !
    !    here are the local variables
    !
    INTEGER                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE   :: for(:, :), for_g(:, :)
    REAL(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                               t2, z, zp, L, tmp, r1, r2, f1(3), f2(3), &
                               arg001, arg002, arg101, arg102
    COMPLEX(DP), ALLOCATABLE :: vg_f(:, :), vg_f_r(:, :), rhog3(:, :)
    COMPLEX(DP)              :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rhog3(n3, ng_2d) = aux(ng)
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(aux(ng))
      ENDIF
    ENDDO

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0 + esm_w

    ALLOCATE (for_g(3, nat))
    for_g(:, :) = 0.d0

!**** for gp!=0 *********
    ALLOCATE (for(3, nat), vg_f(dfftp%nr3x, 3), vg_f_r(dfftp%nr3x, 3))
    for(:, :) = 0.d0
    vg_f_r(:, :) = (0.d0, 0.d0)
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE

      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)

      DO it = 1, nat
        IF (gamma_only) THEN
          tt = -fpi*zv(ityp(it))/sa*2.d0
        ELSE
          tt = -fpi*zv(ityp(it))/sa
        ENDIF
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                   + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc1
          arg001 = gp*(z - zp)
          arg002 = -gp*(z - zp)
          arg101 = gp/2.d0/tmp + tmp*(z - zp)
          arg102 = gp/2.d0/tmp - tmp*(z - zp)
          t1 = exp_erfc(arg002, arg102)
          t2 = exp_erfc(arg001, arg101)
          c1(1) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k1
          c1(2) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k2
          c1(3) = CMPLX(cc, ss, kind=DP)*(t1 - t2)/4.d0
          c2(:) = (0.d0, 0.d0)
          vg_f_r(iz, :) = tt*(c1(:) + c2(:))
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        CALL cft_1z(vg_f_r(:, 2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 2))
        CALL cft_1z(vg_f_r(:, 3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 3))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(:) = dble(vg_f(iz, :))
          f2(:) = aimag(vg_f(iz, :))
          for(:, it) = for(:, it) - r1*f1(:) - r2*f2(:)
        ENDDO
      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)
    DEALLOCATE (for, vg_f, vg_f_r)

!***** for gp==0********
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg_f(dfftp%nr3x, 1), vg_f_r(dfftp%nr3x, 1))
      vg_f_r(:, 1) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc1
          cc1 = 0.5d0*erf(tmp*(z - zp))
          cc2 = (0.d0, 0.d0)

          vg_f_r(iz, 1) = tt*(cc1 + cc2)
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(3) = dble(vg_f(iz, 1))
          f2(3) = aimag(vg_f(iz, 1))
          for_g(3, it) = for_g(3, it) - r1*f1(3) - r2*f2(3)
        ENDDO
      ENDDO

      DEALLOCATE (vg_f, vg_f_r)
    ENDIF ! IF( ng_2d > 0 )

!***** sum short_range part and long_range part in local potential force
!***** at cartecian coordinate

    DO it = 1, nat
                                             ! factor e2: hartree -> Ry.
      forcelc(1, it) = forcelc(1, it) &
                       + sum(for_g(1:2, it)*bg(1, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(2, it) = forcelc(2, it) &
                       + sum(for_g(1:2, it)*bg(2, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(3, it) = forcelc(3, it) + for_g(3, it)*omega*e2
    ENDDO

    DEALLOCATE (for_g)
    DEALLOCATE (rhog3)

    RETURN
  END SUBROUTINE esm_force_lc_bc1

  SUBROUTINE esm_force_lc_bc2(aux, forcelc)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z

    IMPLICIT NONE
    COMPLEX(DP), INTENT(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
    REAL(DP),    INTENT(inout) :: forcelc(3, nat)
    !
    !    here are the local variables
    !
    INTEGER                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE   :: for(:, :), for_g(:, :)
    REAL(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                               t2, z, zp, L, tmp, r1, r2, f1(3), f2(3), &
                               arg001, arg002, arg003, arg005, &
                               arg006, arg008, arg009, arg101, arg102
    COMPLEX(DP), ALLOCATABLE :: vg_f(:, :), vg_f_r(:, :), rhog3(:, :)
    COMPLEX(DP)              :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rhog3(n3, ng_2d) = aux(ng)
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(aux(ng))
      ENDIF
    ENDDO

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0 + esm_w

    ALLOCATE (for_g(3, nat))
    for_g(:, :) = 0.d0

!**** for gp!=0 *********
    ALLOCATE (for(3, nat), vg_f(dfftp%nr3x, 3), vg_f_r(dfftp%nr3x, 3))
    for(:, :) = 0.d0
    vg_f_r(:, :) = (0.d0, 0.d0)
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE

      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)

      DO it = 1, nat
        IF (gamma_only) THEN
          tt = -fpi*zv(ityp(it))/sa*2.d0
        ELSE
          tt = -fpi*zv(ityp(it))/sa
        ENDIF
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                   + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc2
          arg001 = gp*(z - zp)
          arg002 = -gp*(z - zp)
          arg003 = -gp*(z + zp + 2.d0*z1)
          arg005 = gp*(z + zp - 2.d0*z1)
          arg006 = -gp*(z - zp + 4.d0*z1)
          arg008 = gp*(z - zp - 4.d0*z1)
          arg009 = -4.d0*gp*z1
          arg101 = gp/2.d0/tmp + tmp*(z - zp)
          arg102 = gp/2.d0/tmp - tmp*(z - zp)
          t1 = exp_erfc(arg002, arg102)
          t2 = exp_erfc(arg001, arg101)
          c1(1) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k1
          c1(2) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k2
          c1(3) = CMPLX(cc, ss, kind=DP)*(t1 - t2)/4.d0
          c2(1) = CMPLX(ss, -cc, kind=DP)*(EXP(arg008) + EXP(arg006) &
                  - EXP(arg005) - EXP(arg003))/(1.d0 - EXP(arg009))/2.d0/gp*k1
          c2(2) = CMPLX(ss, -cc, kind=DP)*(EXP(arg008) + EXP(arg006) &
                  - EXP(arg005) - EXP(arg003))/(1.d0 - EXP(arg009))/2.d0/gp*k2
          c2(3) = CMPLX(cc, ss, kind=DP)*(-EXP(arg008) + EXP(arg006) &
                  - EXP(arg005) + EXP(arg003))/(1.d0 - EXP(arg009))/2.d0
          vg_f_r(iz, :) = tt*(c1(:) + c2(:))
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        CALL cft_1z(vg_f_r(:, 2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 2))
        CALL cft_1z(vg_f_r(:, 3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 3))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(:) = dble(vg_f(iz, :))
          f2(:) = aimag(vg_f(iz, :))
          for(:, it) = for(:, it) - r1*f1(:) - r2*f2(:)
        ENDDO
      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)
    DEALLOCATE (for, vg_f, vg_f_r)

!***** for gp==0********
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg_f(dfftp%nr3x, 1), vg_f_r(dfftp%nr3x, 1))
      vg_f_r(:, 1) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc2
          cc1 = 0.5d0*erf(tmp*(z - zp))
          cc2 = -0.5d0*(z/z1)
          vg_f_r(iz, 1) = tt*(cc1 + cc2)
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(3) = dble(vg_f(iz, 1))
          f2(3) = aimag(vg_f(iz, 1))
          for_g(3, it) = for_g(3, it) - r1*f1(3) - r2*f2(3)
        ENDDO
      ENDDO
      DEALLOCATE (vg_f, vg_f_r)
    ENDIF ! IF( ng_2d > 0 )

!***** sum short_range part and long_range part in local potential force
!***** at cartecian coordinate

    DO it = 1, nat
      ! factor e2: hartree -> Ry.
      forcelc(1, it) = forcelc(1, it) &
                       + sum(for_g(1:2, it)*bg(1, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(2, it) = forcelc(2, it) &
                       + sum(for_g(1:2, it)*bg(2, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(3, it) = forcelc(3, it) + for_g(3, it)*omega*e2
    ENDDO

    DEALLOCATE (for_g)
    DEALLOCATE (rhog3)

    RETURN
  END SUBROUTINE esm_force_lc_bc2

  SUBROUTINE esm_force_lc_bc3(aux, forcelc)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z

    IMPLICIT NONE
    COMPLEX(DP), INTENT(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
    REAL(DP),    INTENT(inout) :: forcelc(3, nat)
    !
    !    here are the local variables
    !
    INTEGER                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE   :: for(:, :), for_g(:, :)
    REAL(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                               t2, z, zp, L, tmp, r1, r2, f1(3), f2(3), &
                               arg001, arg002, arg003, arg101, arg102
    COMPLEX(DP), ALLOCATABLE :: vg_f(:, :), vg_f_r(:, :), rhog3(:, :)
    COMPLEX(DP)              :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rhog3(n3, ng_2d) = aux(ng)
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(aux(ng))
      ENDIF
    ENDDO

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0 + esm_w

    ALLOCATE (for_g(3, nat))
    for_g(:, :) = 0.d0

!**** for gp!=0 *********
    ALLOCATE (for(3, nat), vg_f(dfftp%nr3x, 3), vg_f_r(dfftp%nr3x, 3))
    for(:, :) = 0.d0
    vg_f_r(:, :) = (0.d0, 0.d0)
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE

      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)

      DO it = 1, nat
        IF (gamma_only) THEN
          tt = -fpi*zv(ityp(it))/sa*2.d0
        ELSE
          tt = -fpi*zv(ityp(it))/sa
        ENDIF
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                   + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc3
          arg001 = gp*(z - zp)
          arg002 = -gp*(z - zp)
          arg003 = gp*(z + zp - 2.d0*z1)
          arg101 = gp/2.d0/tmp + tmp*(z - zp)
          arg102 = gp/2.d0/tmp - tmp*(z - zp)
          t1 = exp_erfc(arg002, arg102)
          t2 = exp_erfc(arg001, arg101)
          c1(1) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k1
          c1(2) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k2
          c1(3) = CMPLX(cc, ss, kind=DP)*(t1 - t2)/4.d0
          c2(1) = CMPLX(ss, -cc, kind=DP)*(-EXP(arg003))/2.d0/gp*k1
          c2(2) = CMPLX(ss, -cc, kind=DP)*(-EXP(arg003))/2.d0/gp*k2
          c2(3) = CMPLX(cc, ss, kind=DP)*(-EXP(arg003))/2.d0

          vg_f_r(iz, :) = tt*(c1(:) + c2(:))
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        CALL cft_1z(vg_f_r(:, 2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 2))
        CALL cft_1z(vg_f_r(:, 3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 3))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(:) = dble(vg_f(iz, :))
          f2(:) = aimag(vg_f(iz, :))
          for(:, it) = for(:, it) - r1*f1(:) - r2*f2(:)
        ENDDO
      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)
    DEALLOCATE (for, vg_f, vg_f_r)

!***** for gp==0********
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg_f(dfftp%nr3x, 1), vg_f_r(dfftp%nr3x, 1))
      vg_f_r(:, 1) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc3
          cc1 = 0.5d0*erf(tmp*(z - zp))
          cc2 = -0.5d0
          vg_f_r(iz, 1) = tt*(cc1 + cc2)
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(3) = dble(vg_f(iz, 1))
          f2(3) = aimag(vg_f(iz, 1))
          for_g(3, it) = for_g(3, it) - r1*f1(3) - r2*f2(3)
        ENDDO
      ENDDO

      DEALLOCATE (vg_f, vg_f_r)
    ENDIF ! IF( ng_2d > 0 )

!***** sum short_range part and long_range part in local potential force
!***** at cartecian coordinate

    DO it = 1, nat
                                             ! factor e2: hartree -> Ry.
      forcelc(1, it) = forcelc(1, it) &
                       + sum(for_g(1:2, it)*bg(1, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(2, it) = forcelc(2, it) &
                       + sum(for_g(1:2, it)*bg(2, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(3, it) = forcelc(3, it) + for_g(3, it)*omega*e2
    ENDDO

    DEALLOCATE (for_g)
    DEALLOCATE (rhog3)

    RETURN
  END SUBROUTINE esm_force_lc_bc3

  SUBROUTINE esm_force_lc_bc4(aux, forcelc)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z

    IMPLICIT NONE
    COMPLEX(DP), INTENT(in)    :: aux(dfftp%nnr) ! aux contains n(G) (input)
    REAL(DP),    INTENT(inout) :: forcelc(3, nat)
    !
    !    here are the local variables
    !
    INTEGER                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, ng_2d
    REAL(DP), ALLOCATABLE   :: for(:, :), for_g(:, :)
    REAL(DP)                :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, t1, &
                               t2, z, zp, L, tmp, r1, r2, f1(3), &
                               f2(3), arg001, arg002, arg003, arg005, &
                               arg006, arg008, arg009, arg011, arg101, &
                               arg102, arg103, arg104, arg106, arg107, &
                               arg109, arg111, arg113, aaa, t3, alpha, beta, &
                               kappa, lambda, xi, chi
    COMPLEX(DP), ALLOCATABLE :: vg_f(:, :), vg_f_r(:, :), rhog3(:, :)
    COMPLEX(DP)              :: c1(3), c2(3), cc1, cc2

! Map to FFT mesh
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rhog3(n3, ng_2d) = aux(ng)
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(aux(ng))
      ENDIF
    ENDDO

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0
    z1 = z0 + esm_w
    aaa = esm_a

    ALLOCATE (for_g(3, nat))
    for_g(:, :) = 0.d0

!**** for gp!=0 *********
    ALLOCATE (for(3, nat), vg_f(dfftp%nr3x, 3), vg_f_r(dfftp%nr3x, 3))
    for(:, :) = 0.d0
    vg_f_r(:, :) = (0.d0, 0.d0)
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE

      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)

      DO it = 1, nat
        IF (gamma_only) THEN
          tt = -fpi*zv(ityp(it))/sa*2.d0
        ELSE
          tt = -fpi*zv(ityp(it))/sa
        ENDIF
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                   + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
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
            c1(1) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k1
            c1(2) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/gp*k2
            c2(1) = CMPLX(ss, -cc, kind=DP)*t3/2.d0*k1
            c2(2) = CMPLX(ss, -cc, kind=DP)*t3/2.d0*k2
            t1 = exp_erfc(arg001, arg103) - exp_erfc(arg001, arg101)
            t2 = exp_erfc(arg002, arg102) &
                 - kappa/alpha*exp_erfc(arg003, arg104)
            t3 = -xi/alpha*exp_erfc(arg006, arg109)
            c1(3) = CMPLX(cc, ss, kind=DP)*(t1 + t2)/4.d0
            c2(3) = CMPLX(cc, ss, kind=DP)*t3/2.d0
          ELSE
            t1 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111)
            t2 = exp_erfc(arg008, arg107) &
                 - beta/alpha*exp_erfc(arg009, arg109)
            t3 = exp_erfc(arg005, arg104)/alpha
            c1(1) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/lambda*k1
            c1(2) = CMPLX(ss, -cc, kind=DP)*(t1 + t2)/4.d0/lambda*k2
            c2(1) = CMPLX(ss, -cc, kind=DP)*t3/2.d0*k1
            c2(2) = CMPLX(ss, -cc, kind=DP)*t3/2.d0*k2
            t1 = chi*(exp_erfc(arg011, arg111) - exp_erfc(arg011, arg113))
            t2 = -xi*(exp_erfc(arg008, arg107) &
                      + beta/alpha*exp_erfc(arg009, arg109))
            t3 = gp/alpha*exp_erfc(arg005, arg104)
            c1(3) = CMPLX(cc, ss, kind=DP)*(t1 + t2)/4.d0/lambda
            c2(3) = CMPLX(cc, ss, kind=DP)*t3/2.d0
          ENDIF
          vg_f_r(iz, :) = tt*(c1(:) + c2(:))
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        CALL cft_1z(vg_f_r(:, 2), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 2))
        CALL cft_1z(vg_f_r(:, 3), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 3))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(:) = dble(vg_f(iz, :))
          f2(:) = aimag(vg_f(iz, :))
          for(:, it) = for(:, it) - r1*f1(:) - r2*f2(:)
        ENDDO
      ENDDO
    ENDDO
    for_g(:, :) = for_g(:, :) + for(:, :)
    DEALLOCATE (for, vg_f, vg_f_r)

!***** for gp==0********
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg_f(dfftp%nr3x, 1), vg_f_r(dfftp%nr3x, 1))
      vg_f_r(:, 1) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        zp = tau(3, it)
        IF (zp .gt. at(3, 3)*0.5) zp = zp - at(3, 3)
        zp = zp*alat
        DO iz = 1, dfftp%nr3
          k3 = iz - 1
          IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            k3 = k3 - dfftp%nr3
          END IF
          z = DBLE(k3) / DBLE(dfftp%nr3) * L
          ! bc4
          arg006 = aaa**2/tmp**2 + 2.d0*aaa*(z1 - zp)
          arg101 = tmp*(z - zp)
          arg102 = tmp*(z1 - zp)
          arg104 = aaa/tmp + tmp*(z - zp)
          arg106 = aaa/tmp + tmp*(z1 - zp)
          IF (z < z1) THEN
            cc1 = 0.5d0*(erf(arg101) - erf(arg102))
            cc2 = -0.5d0*exp_erfc(arg006, arg106)
          ELSE
            cc1 = 0.d0
            cc2 = -0.5d0*exp_erfc(arg006, arg104)
          ENDIF
          vg_f_r(iz, 1) = tt*(cc1 + cc2)
        ENDDO
        CALL cft_1z(vg_f_r(:, 1), 1, dfftp%nr3, dfftp%nr3, -1, vg_f(:, 1))
        DO iz = 1, dfftp%nr3
          r1 = dble(rhog3(iz, ng_2d))
          r2 = aimag(rhog3(iz, ng_2d))
          f1(3) = dble(vg_f(iz, 1))
          f2(3) = aimag(vg_f(iz, 1))
          for_g(3, it) = for_g(3, it) - r1*f1(3) - r2*f2(3)
        ENDDO
      ENDDO
      DEALLOCATE (vg_f, vg_f_r)
    ENDIF ! IF( ng_2d > 0 )

!***** sum short_range part and long_range part in local potential force
!***** at cartecian coordinate

    DO it = 1, nat
                                             ! factor e2: hartree -> Ry.
      forcelc(1, it) = forcelc(1, it) &
                       + sum(for_g(1:2, it)*bg(1, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(2, it) = forcelc(2, it) &
                       + sum(for_g(1:2, it)*bg(2, 1:2))*sqrt(tpiba2)*omega*e2
      forcelc(3, it) = forcelc(3, it) + for_g(3, it)*omega*e2
    ENDDO

    DEALLOCATE (for_g)
    DEALLOCATE (rhog3)

    RETURN
  END SUBROUTINE esm_force_lc_bc4
END MODULE esm_force_mod
