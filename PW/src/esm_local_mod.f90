MODULE esm_local_mod

  USE kinds, ONLY : DP
  USE esm_common_mod, ONLY : esm_nfit, esm_efield, esm_w, esm_a, esm_bc, &
                           & mill_2d, imill_2d, ngm_2d, exp_erfc
  IMPLICIT NONE

CONTAINS
  SUBROUTINE esm_local(aux)
    USE kinds,    ONLY: DP
    USE fft_base, ONLY: dfftp
    IMPLICIT NONE
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)

    IF (esm_bc == 'pbc') THEN
      CALL esm_local_pbc(aux)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_local_bc1(aux)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_local_bc2(aux)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_local_bc3(aux)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_local_bc4(aux)
    END IF

  END SUBROUTINE esm_local

!-----------------------------------------------------------------------
!--------------ESM LOCAL POTENTIAL SUBROUTINE---------------------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_local_pbc(aux)
    USE fft_base, ONLY : dfftp
    IMPLICIT NONE
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)

    STOP 'esm_local must not be called for esm_bc = pbc'

  END SUBROUTINE esm_local_pbc

  SUBROUTINE esm_local_bc1(aux)

    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)
    !
    !    here the local variables
    !
    REAL(DP)                :: t(2), tt, gp, gp2, sa, z0, pp, cc, ss, &
                               t1, t2, z, zp, tmp, L, z_l, z_r, &
                               arg001, arg002, arg101, arg102
    INTEGER                 :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                               nz_r, ng_2d
    COMPLEX(DP)             :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:, :), vg(:), vg_r(:)

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    ALLOCATE (vloc3(dfftp%nr3, ngm_2d))

! for gp!=0
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      vg_r(:) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                   + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        cs = CMPLX(cc, ss, kind=DP)
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
          cc1 = cs*(t1 + t2)/4.d0/gp
          cc2 = (0.d0, 0.d0)
          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO
    ENDDO
    DEALLOCATE (vg, vg_r)

    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      vg_r(:) = (0.d0, 0.d0)
! for smoothing
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = dble(nz_l - 1)*L/dble(dfftp%nr3) - L
      z_r = dble(nz_r - 1)*L/dble(dfftp%nr3)
! for gp=0
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
          arg001 = -tmp**2*(z - zp)**2
          arg101 = tmp*(z - zp)
          cc1 = 0.5d0*(-(z - zp)*erf(arg101) - exp(arg001)/tmp/sqrt(pi))
          cc2 = (0.d0, 0.d0)

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
        ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc1
        f1 = f1 + tt*0.5d0*(-(z_r - zp)*erf(tmp*(z_r - zp)) &
             - exp(-tmp**2*(z_r - zp)**2)/tmp/sqrt(pi))
        f2 = f2 + tt*0.5d0*(-(z_l - zp)*erf(tmp*(z_l - zp)) &
             - exp(-tmp**2*(z_l - zp)**2)/tmp/sqrt(pi))
        f3 = f3 - tt*0.5d0*erf(tmp*(z_r - zp))
        f4 = f4 - tt*0.5d0*erf(tmp*(z_l - zp))
      ENDDO
      ! for smoothing
      ! factor e2: hartree -> Ry.
      f1 = f1*e2; f2 = f2*e2; f3 = f3*e2; f4 = f4*e2
      z_r = z_r
      z_l = z_l + L
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
      DO iz = nz_r, nz_l
        z = dble(iz - 1)/dble(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = aux(ng) + vloc3(n3,ng_2d)    
    ENDDO

    DEALLOCATE (vloc3)

    RETURN
  END SUBROUTINE esm_local_bc1

  SUBROUTINE esm_local_bc2(aux)

    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)
    !
    !    here the local variables
    !
    REAL(DP)                 :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                                t1, t2, z, zp, v0, tmp, L, z_l, z_r, &
                                arg001, arg002, arg003, arg004, arg005, &
                                arg006, arg007, arg101, arg102
    INTEGER                  :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                                nz_r, ng_2d
    COMPLEX(DP)              :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:, :), vg(:), vg_r(:)

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    z1 = z0 + esm_w
    v0 = esm_efield*z1*2.d0/e2 ! factor 1/e2: unit Ry. -> hartree
    ALLOCATE (vloc3(dfftp%nr3, ngm_2d))

! for gp!=0
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      vg_r(:) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
             + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        cs = CMPLX(cc, ss, kind=DP)
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
          arg004 = gp*(z + zp - 2.d0*z1)
          arg005 = -gp*(z - zp + 4.d0*z1)
          arg006 = gp*(z - zp - 4.d0*z1)
          arg007 = -4.d0*gp*z1
          arg101 = gp/2.d0/tmp + tmp*(z - zp)
          arg102 = gp/2.d0/tmp - tmp*(z - zp)
          t1 = exp_erfc(arg002, arg102)
          t2 = exp_erfc(arg001, arg101)
          cc1 = cs*(t1 + t2)/4.d0/gp
          cc2 = cs*(exp(arg006) + exp(arg005) - exp(arg004) - exp(arg003)) &
                /(1.d0 - exp(arg007))/2.d0/gp

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO
    ENDDO
    DEALLOCATE (vg, vg_r)

    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      vg_r(:) = (0.d0, 0.d0)
! for smoothing
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = dble(nz_l - 1)*L/dble(dfftp%nr3) - L
      z_r = dble(nz_r - 1)*L/dble(dfftp%nr3)
! add constant potential (capacitor term)
      ! bc2
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        vg_r(iz) = -0.5d0*v0*(z - z1)/z1*e2 ! factor e2: hartree -> Ry.
      ENDDO
      f1 = -0.5d0*v0*(z_r - z1)/z1 ! unit: hartree
      f2 = -0.5d0*v0*(z_l - z1)/z1 ! unit: hartree
      f3 = -0.5d0*v0/z1 ! unit: hartree/a.u.
      f4 = -0.5d0*v0/z1 ! unit: harteee/a.u.

! for gp=0
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
          arg001 = -tmp**2*(z - zp)**2
          arg101 = tmp*(z - zp)
          cc1 = 0.5d0*(-(z - zp)*erf(arg101) - exp(arg001)/tmp/sqrt(pi))
          cc2 = 0.5d0*(z1 - z*zp/z1)

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
        ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc2
        f1 = f1 + tt*0.5d0*(-(z_r - zp)*erf(tmp*(z_r - zp)) &
                            - exp(-tmp**2*(z_r - zp)**2)/tmp/sqrt(pi))
        f2 = f2 + tt*0.5d0*(-(z_l - zp)*erf(tmp*(z_l - zp)) &
                            - exp(-tmp**2*(z_l - zp)**2)/tmp/sqrt(pi))
        f3 = f3 - tt*0.5d0*erf(tmp*(z_r - zp))
        f4 = f4 - tt*0.5d0*erf(tmp*(z_l - zp))
        f1 = f1 + tt*0.5d0*(z1 - z_r*zp/z1)
        f2 = f2 + tt*0.5d0*(z1 - z_l*zp/z1)
        f3 = f3 + tt*(-0.5d0*(zp/z1))
        f4 = f4 + tt*(-0.5d0*(zp/z1))
      ENDDO
      ! for smoothing
      ! factor e2: hartree -> Ry.
      f1 = f1*e2; f2 = f2*e2; f3 = f3*e2; f4 = f4*e2
      z_r = z_r
      z_l = z_l + L
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
      DO iz = nz_r, nz_l
        z = dble(iz - 1)/dble(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = aux(ng) + vloc3(n3,ng_2d)
    ENDDO

    DEALLOCATE (vloc3)

    RETURN
  END SUBROUTINE esm_local_bc2

  SUBROUTINE esm_local_bc3(aux)

    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)
    !
    !    here the local variables
    !
    REAL(DP)                 :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                                t1, t2, z, zp, tmp, L, z_l, z_r, &
                                arg001, arg002, arg003, arg101, arg102
    INTEGER                  :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                                nz_r, ng_2d
    COMPLEX(DP)              :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:, :), vg(:), vg_r(:)

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    z1 = z0 + esm_w
    ALLOCATE (vloc3(dfftp%nr3, ngm_2d))

! for gp!=0
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      vg_r(:) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                   + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        cs = CMPLX(cc, ss, kind=DP)
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
          cc1 = cs*(t1 + t2)/4.d0/gp
          cc2 = cs*(-exp(arg003))/2.d0/gp

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO
    ENDDO
    DEALLOCATE (vg, vg_r)

    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      vg_r(:) = (0.d0, 0.d0)
! for smoothing
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = dble(nz_l - 1)*L/dble(dfftp%nr3) - L
      z_r = dble(nz_r - 1)*L/dble(dfftp%nr3)

! for gp=0
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
          arg001 = -tmp**2*(z - zp)**2
          arg101 = tmp*(z - zp)
          cc1 = 0.5d0*(-(z - zp)*erf(arg101) - exp(arg001)/tmp/sqrt(pi))
          cc2 = 0.5d0*(2.d0*z1 - z - zp)

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
        ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc3
        f1 = f1 + tt*0.5d0*(-(z_r - zp)*erf(tmp*(z_r - zp)) &
             - exp(-tmp**2*(z_r - zp)**2)/tmp/sqrt(pi))
        f2 = f2 + tt*0.5d0*(-(z_l - zp)*erf(tmp*(z_l - zp)) &
             - exp(-tmp**2*(z_l - zp)**2)/tmp/sqrt(pi))
        f3 = f3 - tt*0.5d0*erf(tmp*(z_r - zp))
        f4 = f4 - tt*0.5d0*erf(tmp*(z_l - zp))
        f1 = f1 + tt*0.5d0*(2.d0*z1 - z_r - zp)
        f2 = f2 + tt*0.5d0*(2.d0*z1 - z_l - zp)
        f3 = f3 - tt*0.5d0
        f4 = f4 - tt*0.5d0
      ENDDO
      ! for smoothing
      ! factor e2: hartree -> Ry.
      f1 = f1*e2; f2 = f2*e2; f3 = f3*e2; f4 = f4*e2
      z_r = z_r
      z_l = z_l + L
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
      DO iz = nz_r, nz_l
        z = dble(iz - 1)/dble(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = aux(ng) + vloc3(n3,ng_2d)
    ENDDO

    DEALLOCATE (vloc3)

    RETURN
  END SUBROUTINE esm_local_bc3

  SUBROUTINE esm_local_bc4(aux)

    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE control_flags,    ONLY : gamma_only
    USE cell_base,        ONLY : at, bg, alat, tpiba2, omega
    USE ions_base,        ONLY : zv, nat, tau, ityp
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    ! aux contains v_loc_short(G) (input) and v_loc(G) (output)
    COMPLEX(DP), INTENT(inout) :: aux(dfftp%nnr)
    !
    !    here the local variables
    !
    REAL(DP)                 :: t(2), tt, gp, gp2, sa, z1, z0, pp, cc, ss, &
                                t1, t2, z, zp, tmp, L, z_l, z_r, &
                                arg001, arg002, arg003, arg005, &
                                arg006, arg008, arg009, arg011, &
                                arg101, arg102, arg103, arg104, arg106, &
                                arg107, arg109, arg111, arg113, aaa, t3, &
                                alpha, beta, kappa, lambda, xi, chi
    INTEGER                  :: iz, it, k1, k2, k3, ng, n1, n2, n3, nz_l, &
                                nz_r, ng_2d
    COMPLEX(DP)              :: cs, cc1, cc2, a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP), ALLOCATABLE :: vloc3(:, :), vg(:), vg_r(:)

    L = at(3, 3)*alat
    sa = omega/L
    z0 = L/2.d0
    tmp = 1.d0   ! Gaussian width
    z1 = z0 + esm_w
    aaa = esm_a
    ALLOCATE (vloc3(dfftp%nr3, ngm_2d))

! for gp!=0
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      vg_r(:) = (0.d0, 0.d0)
      DO it = 1, nat
        tt = -fpi*zv(ityp(it))/sa
        pp = -tpi*(tau(1, it)*(k1*bg(1, 1) + k2*bg(1, 2)) &
                 + tau(2, it)*(k1*bg(2, 1) + k2*bg(2, 2)))
        cc = cos(pp)
        ss = sin(pp)
        cs = CMPLX(cc, ss, kind=DP)
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
          beta  = aaa + gp - sqrt(aaa**2 + gp**2)
          kappa = aaa - gp + sqrt(aaa**2 + gp**2)
          xi    = aaa      + sqrt(aaa**2 + gp**2)
          chi   = aaa      - sqrt(aaa**2 + gp**2)
          lambda = sqrt(aaa**2 + gp**2)
          arg001 =  gp*(z - zp)
          arg002 = -gp*(z - zp)
          arg003 =  gp*(z + zp - 2.d0*z1)
          arg005 = -gp*(z1 - zp) - xi*(z - z1)
          arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
          arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
          arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
          arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
          arg101 =  gp/2.d0/tmp + tmp*(z - zp)
          arg102 =  gp/2.d0/tmp - tmp*(z - zp)
          arg103 =  gp/2.d0/tmp + tmp*(z1 - zp)
          arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
          arg107 =  xi/2.d0/tmp + tmp*(z - zp)
          arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
          arg111 = chi/2.d0/tmp + tmp*(z - zp)
          arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
          IF (z < z1) THEN
            t1 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103)
            t2 = exp_erfc(arg002, arg102) &
                 - kappa/alpha*exp_erfc(arg003, arg104)
            t3 = exp_erfc(arg006, arg109)/alpha
            cc1 = cs*(t1 + t2)/4.d0/gp
            cc2 = cs*t3/2.d0
          ELSE
            t1 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111)
            t2 = exp_erfc(arg008, arg107) &
                 - beta/alpha*exp_erfc(arg009, arg109)
            t3 = exp_erfc(arg005, arg104)/alpha
            cc1 = cs*(t1 + t2)/4.d0/lambda
            cc2 = cs*t3/2.d0
          ENDIF

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO
    ENDDO
    DEALLOCATE (vg, vg_r)

    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      vg_r(:) = (0.d0, 0.d0)
! for smoothing
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = dble(nz_l - 1)*L/dble(dfftp%nr3) - L
      z_r = dble(nz_r - 1)*L/dble(dfftp%nr3)
! for gp=0
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
            t3 = 0.5d0/aaa - (z - z1) + exp(arg002)/tmp/sqrt(pi) &
                 - exp(arg001)/tmp/sqrt(pi)
            cc1 = (t1 + t2)/2.d0
            cc2 = t3/2.d0
          ELSE
            t1 = -exp_erfc(arg005, arg101)/aaa
            t2 =  exp_erfc(arg006, arg104)/aaa
            t3 =  exp(arg005)/aaa
            cc1 = (t1 + t2)/4.d0
            cc2 = t3/2.d0
          ENDIF

          vg_r(iz) = vg_r(iz) + tt*(cc1 + cc2)*e2 ! factor e2: hartree -> Ry.
        ENDDO
        ! smoothing cell edge potential (avoiding unphysical oscillation)
        ! bc4
        arg002 = -tmp**2*(z1 - zp)**2
        arg006 = aaa**2/tmp**2 + 2.d0*aaa*(z1 - zp)
        arg102 = tmp*(z1 - zp)
        arg106 = aaa/tmp + tmp*(z1 - zp)
        !-right only
        arg005 = -2.d0*aaa*(z_r - z1)
        arg101 = tmp*(z_r - zp)
        arg104 = aaa/tmp + tmp*(z_r - zp)
        !--
        t1 = -exp_erfc(arg005, arg101)/aaa
        t2 = exp_erfc(arg006, arg104)/aaa
        t3 = exp(arg005)/aaa
        f1 = f1 + tt*((t1 + t2)/2.d0 + t3)/2.d0
        f3 = f3 - tt*0.5d0*exp(arg005)*(1.d0 + erf(arg101))
        !-left only
        arg001 = -tmp**2*(z_l - zp)**2
        arg101 = tmp*(z_l - zp)
        !--
        t1 = -(z_l - zp)*erf(arg101) + (0.5d0/aaa + z1 - zp)*erf(arg102)
        t2 = 0.5d0/aaa*exp_erfc(arg006, arg106)
        t3 = 0.5d0/aaa - (z_l - z1) + exp(arg002)/tmp/sqrt(pi) &
             - exp(arg001)/tmp/sqrt(pi)
        f2 = f2 + tt*(t1 + t2 + t3)/2.d0
        f4 = f4 - tt*0.5d0*(1.d0 + erf(arg101))
      ENDDO
      ! for smoothing
      ! factor e2: hartree -> Ry.
      f1 = f1*e2; f2 = f2*e2; f3 = f3*e2; f4 = f4*e2
      z_r = z_r
      z_l = z_l + L
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
      DO iz = nz_r, nz_l
        z = dble(iz - 1)/dble(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      DO iz = 1, dfftp%nr3
        vloc3(iz, ng_2d) = vg(iz)
      ENDDO

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Map to FFT mesh (dfftp%nrx)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = aux(ng) + vloc3(n3,ng_2d)
    ENDDO

    DEALLOCATE (vloc3)

    RETURN
  END SUBROUTINE esm_local_bc4
END MODULE esm_local_mod
