MODULE esm_hartree_mod

  USE kinds, ONLY : DP
  USE esm_common_mod, ONLY : esm_nfit, esm_w, esm_a, esm_bc, &
                           & mill_2d, imill_2d, ngm_2d
  IMPLICIT NONE

CONTAINS
  !-----------------------------------------------------------------------
  !--------------ESM ENERGY AND POTENTIAL SUBROUTINE----------------------
  !-----------------------------------------------------------------------
  SUBROUTINE esm_hartree(rhog, ehart, aux)
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    USE fft_base, ONLY : dfftp
    IMPLICIT NONE
    REAL(DP)    :: ehart             !  Hartree energy
    COMPLEX(DP) :: rhog(ngm)         !  n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)    !  v_h(G)

    IF (esm_bc == 'pbc') THEN
      CALL esm_hartree_pbc(rhog, ehart, aux)
    ELSE IF (esm_bc == 'bc1') THEN
      CALL esm_hartree_bc1(rhog, ehart, aux)
    ELSE IF (esm_bc == 'bc2') THEN
      CALL esm_hartree_bc2(rhog, ehart, aux)
    ELSE IF (esm_bc == 'bc3') THEN
      CALL esm_hartree_bc3(rhog, ehart, aux)
    ELSE IF (esm_bc == 'bc4') THEN
      CALL esm_hartree_bc4(rhog, ehart, aux)
    END IF

  END SUBROUTINE esm_hartree

!
!-----------------------------------------------------------------------
!--------------ESM HARTREE SUBROUTINE-----------------------------------
!-----------------------------------------------------------------------
  SUBROUTINE esm_hartree_pbc(rhog, ehart, aux)
    USE gvect,    ONLY : ngm
    USE fft_base, ONLY : dfftp
    IMPLICIT NONE
    REAL(DP)    :: ehart             !  Hartree energy
    COMPLEX(DP) :: rhog(ngm)         !  n(G)
    COMPLEX(DP) :: aux(dfftp%nnr)    !  v_h(G)

    STOP 'esm_hartree must not be called for esm_bc = pbc'

  END SUBROUTINE esm_hartree_pbc

  SUBROUTINE esm_hartree_bc1(rhog, ehart, aux)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z

    !
    IMPLICIT NONE
    !
    REAL(DP)                 :: ehart             !  Hartree energy
    COMPLEX(DP)              :: rhog(ngm)         !  n(G)
    COMPLEX(DP)              :: aux(dfftp%nnr)    !  v_h(G)
    !
    !    here the local variables
    !
    INTEGER                  :: k1, k2, k3, iz, igz, ng, n1, n2, n3, nz_r, nz_l, &
                                ng_2d
    REAL(DP)                 :: t(2), z, z0, gp, gp2, kn, cc0, ss0, L, &
                                z_l, z_r, eh, arg1, arg2
    COMPLEX(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                                s_r, s_l, rg3, tmp1, tmp2, tmp3
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :), vg(:), vg_r(:), vg3(:, :)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d), vg3(dfftp%nr3, ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rg3 = rhog(ng)
      rhog3(n3, ng_2d) = rg3
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(rg3)
      ENDIF
    ENDDO
! End mapping
!
    vg3(:, :) = (0.d0, 0.d0)
    L = at(3, 3)*alat
    z0 = L/2.d0
    ci = (0.d0, 1.d0)

!****For gp!=0 case ********************
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0);
      vg(:) = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        rg3 = rhog3(iz, ng_2d)
        ! bc1
        vg(iz) = fpi*rg3/(gp**2 + kn**2)
        tmp1 = tmp1 + rg3*(cc0 + ci*ss0)/(gp - ci*kn)
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/(gp + ci*kn)
      ENDDO
      vg3(:, ng_2d) = vg(:)

      ! real part
      vg_r(:) = (0.d0, 0.d0)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        ! bc1
        arg1 = gp*(z - z0)
        arg2 = -gp*(z + z0)
        vg_r(iz) = -tpi/gp*(EXP(arg1)*tmp1 + EXP(arg2)*tmp2)
      ENDDO

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = (vg3(:, ng_2d) + vg(:))*e2 ! factor e2: hartree -> Ry.
    ENDDO
    DEALLOCATE (vg, vg_r)

!****For gp=0 case ********************
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0); tmp3 = (0.d0, 0.d0);
      vg(:) = (0.d0, 0.d0);
      rg3 = rhog3(1, ng_2d)
      vg(1) = -tpi*z0**2*rg3
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        rg3 = rhog3(iz, ng_2d)
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        tmp1 = tmp1 + rg3*ci*(cc0 + ci*ss0)/kn
        tmp2 = tmp2 + rg3*ci*(cc0 - ci*ss0)/kn
        tmp3 = tmp3 + rg3*cc0/kn**2
        vg(iz) = fpi*rg3/(kn**2)
      ENDDO
      vg3(:, ng_2d) = vg(:)

      ! real part
      vg_r(:) = (0.d0, 0.d0)
      rg3 = rhog3(1, ng_2d)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        ! bc1
        vg_r(iz) = -tpi*z**2*rg3 &
                   - tpi*(z - z0)*tmp1 &
                   - tpi*(z + z0)*tmp2 &
                   - fpi*tmp3
      ENDDO

      ! start smoothing
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = DBLE(nz_l - 1)*L/DBLE(dfftp%nr3) - L
      z_r = DBLE(nz_r - 1)*L/DBLE(dfftp%nr3)
      f1 = -tpi*z_r**2*rg3 &
           - tpi*(z_r - z0)*tmp1 &
           - tpi*(z_r + z0)*tmp2 &
           - fpi*tmp3
      f2 = -tpi*z_l**2*rg3 &
           - tpi*(z_l - z0)*tmp1 &
           - tpi*(z_l + z0)*tmp2 &
           - fpi*tmp3
      f3 = -fpi*z_r*rg3 &
           - tpi*tmp1 &
           - tpi*tmp2
      f4 = -fpi*z_l*rg3 &
           - tpi*tmp1 &
           - tpi*tmp2
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
        z = DBLE(iz - 1)/DBLE(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
      ! end smoothing

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = (vg3(:, ng_2d) + vg(:))*e2 ! factor e2: hartree -> Ry.

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      eh = eh + sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
    ENDDO
    ehart = ehart + eh
    IF (gamma_only) THEN
      ehart = ehart*2d0
      ng_2d = imill_2d(0, 0)
      IF (ng_2d > 0) THEN
        ehart = ehart - sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
      ENDIF
    ENDIF
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum(ehart, intra_bgrp_comm)
    !
! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = vg3(n3,ng_2d)
    ENDDO

    DEALLOCATE (rhog3, vg3)

    RETURN
  END SUBROUTINE esm_hartree_bc1

  SUBROUTINE esm_hartree_bc2(rhog, ehart, aux)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP)                 :: ehart             !  Hartree energy
    COMPLEX(DP)              :: rhog(ngm)         !  n(G)
    COMPLEX(DP)              :: aux(dfftp%nnr)    !  v_h(G)
    !
    !    here the local variables
    !
    INTEGER                  :: k1, k2, k3, iz, igz, ng, n1, n2, n3, nz_r, nz_l, &
                                ng_2d
    REAL(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                                z_l, z_r, eh, arg1, arg2, arg3, arg4, arg5
    COMPLEX(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                                s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :), vg(:), vg_r(:), vg3(:, :)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d), vg3(dfftp%nr3, ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rg3 = rhog(ng)
      rhog3(n3, ng_2d) = rg3
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(rg3)
      ENDIF
    ENDDO
! End mapping
!
    vg3(:, :) = (0.d0, 0.d0)
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    ci = (0.d0, 1.d0)

!****For gp!=0 case ********************
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0);
      vg(:) = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        rg3 = rhog3(iz, ng_2d)
        ! bc2
        arg1 = gp*(z1 - z0)
        arg2 = -gp*(z1 - z0)
        vg(iz) = fpi*rg3/(gp**2 + kn**2)
        tmp = ((gp + ci*kn)*EXP(arg1) + (gp - ci*kn)*EXP(arg2))/(2.d0*gp)
        tmp1 = tmp1 + rg3*(cc0 + ci*ss0)/(gp**2 + kn**2)*tmp
        tmp = ((gp - ci*kn)*EXP(arg1) + (gp + ci*kn)*EXP(arg2))/(2.d0*gp)
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/(gp**2 + kn**2)*tmp
      ENDDO
      vg3(:, ng_2d) = vg(:)

      ! real part
      vg_r(:) = (0.d0, 0.d0)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        ! bc2
        arg1 = gp*(z - z1)
        arg2 = -gp*(z + z1)
        arg3 = gp*(z - 3.d0*z1)
        arg4 = -gp*(z + 3.d0*z1)
        arg5 = -4.d0*gp*z1
        vg_r(iz) = -fpi*(EXP(arg1) - EXP(arg4))*tmp1/(1.d0 - EXP(arg5)) &
                   + fpi*(EXP(arg3) - EXP(arg2))*tmp2/(1.d0 - EXP(arg5))
      ENDDO
      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = (vg3(:, ng_2d) + vg(:))*e2 ! factor e2: hartree -> Ry.
    ENDDO
    DEALLOCATE (vg, vg_r)

!****For gp=0 case ********************
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0); tmp3 = (0.d0, 0.d0); tmp4 = (0.d0, 0.d0)
      vg(:) = (0.d0, 0.d0);
      rg3 = rhog3(1, ng_2d)
      vg(1) = tpi*(2.d0*z1 - z0)*z0*rg3
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        rg3 = rhog3(iz, ng_2d)
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        tmp1 = tmp1 + rg3*(cc0 + ci*ss0)/kn**2
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/kn**2
        tmp3 = tmp3 + rg3*ci*cc0/kn
        tmp4 = tmp4 + rg3*ss0/kn
        vg(iz) = fpi*rg3/(kn**2)
      ENDDO
      vg3(:, ng_2d) = vg(:)

      ! real part
      vg_r(:) = (0.d0, 0.d0)
      rg3 = rhog3(1, ng_2d)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        vg_r(iz) = -tpi*z**2*rg3 &
                   - tpi*(z + z1)*tmp1/z1 &
                   + tpi*(z - z1)*tmp2/z1 &
                   - fpi*z*(z1 - z0)/z1*tmp3 &
                   + fpi*(z1 - z0)*tmp4
      ENDDO

      ! start smoothing
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = DBLE(nz_l - 1)*L/DBLE(dfftp%nr3) - L
      z_r = DBLE(nz_r - 1)*L/DBLE(dfftp%nr3)
      f1 = -tpi*z_r**2*rg3 &
           - tpi*(z_r + z1)*tmp1/z1 &
           + tpi*(z_r - z1)*tmp2/z1 &
           - fpi*z_r*(z1 - z0)/z1*tmp3 &
           + fpi*(z1 - z0)*tmp4
      f2 = -tpi*z_l**2*rg3 &
           - tpi*(z_l + z1)*tmp1/z1 &
           + tpi*(z_l - z1)*tmp2/z1 &
           - fpi*z_l*(z1 - z0)/z1*tmp3 &
           + fpi*(z1 - z0)*tmp4
      f3 = -fpi*z_r*rg3 &
           - tpi*tmp1/z1 &
           + tpi*tmp2/z1 &
           - fpi*(z1 - z0)/z1*tmp3
      f4 = -fpi*z_l*rg3 &
           - tpi*tmp1/z1 &
           + tpi*tmp2/z1 &
           - fpi*(z1 - z0)/z1*tmp3
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
        z = DBLE(iz - 1)/DBLE(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
      ! end smoothing

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = (vg3(:, ng_2d) + vg(:))*e2 ! factor e2: hartree -> Ry.

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      eh = eh + sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
    ENDDO
    ehart = ehart + eh
    IF (gamma_only) THEN
      ehart = ehart*2d0
      ng_2d = imill_2d(0, 0)
      IF (ng_2d > 0) THEN
        ehart = ehart - sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
      ENDIF
    ENDIF
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum(ehart, intra_bgrp_comm)
    !
! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = vg3(n3,ng_2d)
    ENDDO

    DEALLOCATE (rhog3, vg3)

    RETURN
  END SUBROUTINE esm_hartree_bc2

  SUBROUTINE esm_hartree_bc3(rhog, ehart, aux)

    USE constants,        ONLY : tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP)                 :: ehart             !  Hartree energy
    COMPLEX(DP)              :: rhog(ngm)         !  n(G)
    COMPLEX(DP)              :: aux(dfftp%nnr)    !  v_h(G)
    !
    !    here the local variables
    !
    INTEGER                  :: k1, k2, k3, iz, igz, ng, n1, n2, n3, nz_r, nz_l, &
                                ng_2d
    REAL(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                                z_l, z_r, eh, arg1, arg2, arg3
    COMPLEX(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                                s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :), vg(:), vg_r(:), vg3(:, :)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d), vg3(dfftp%nr3, ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rg3 = rhog(ng)
      rhog3(n3, ng_2d) = rg3
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(rg3)
      ENDIF
    ENDDO
! End mapping
!
    vg3(:, :) = (0.d0, 0.d0)
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    ci = (0.d0, 1.d0)

!****For gp!=0 case ********************
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0);
      vg(:) = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        rg3 = rhog3(iz, ng_2d)
        ! bc3
        arg1 = gp*(z1 - z0)
        arg2 = -gp*(z1 - z0)
        vg(iz) = fpi*rg3/(gp**2 + kn**2)
        tmp = ((gp + ci*kn)*EXP(arg1) + (gp - ci*kn)*EXP(arg2))/(2.d0*gp)
        tmp1 = tmp1 + rg3*(cc0 + ci*ss0)/(gp**2 + kn**2)*tmp
        tmp = (gp - ci*kn)/gp
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/(gp**2 + kn**2)*tmp
      ENDDO
      vg3(:, ng_2d) = vg(:)

      ! real part
      vg_r(:) = (0.d0, 0.d0)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        ! bc3
        arg1 = gp*(z - z1)
        arg2 = -gp*(z + z0)
        arg3 = gp*(z - z0 - 2.d0*z1)
        vg_r(iz) = -fpi*EXP(arg1)*tmp1 + tpi*(EXP(arg3) - EXP(arg2))*tmp2
      ENDDO

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = (vg3(:, ng_2d) + vg(:))*e2 ! factor e2: hartree -> Ry.
    ENDDO
    DEALLOCATE (vg, vg_r)

!****For gp=0 case ********************
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3))
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0); tmp3 = (0.d0, 0.d0);
      vg(:) = (0.d0, 0.d0);
      rg3 = rhog3(1, ng_2d)
      vg(1) = tpi*(4.d0*z1 - z0)*z0*rg3
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        rg3 = rhog3(iz, ng_2d)
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        tmp1 = tmp1 + rg3*(cc0 + ci*ss0)/kn**2
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/kn
        tmp3 = tmp3 + rg3*(cc0 + ci*ss0)/kn
        vg(iz) = fpi*rg3/(kn**2)
      ENDDO
      vg3(:, ng_2d) = vg(:)

      ! real part
      vg_r(:) = (0.d0, 0.d0)
      rg3 = rhog3(1, ng_2d)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        vg_r(iz) = -tpi*(z**2 + 2.d0*z*z0)*rg3 &
                   - fpi*tmp1 &
                   - fpi*ci*(z - z1)*tmp2 &
                   - fpi*ci*(z1 - z0)*tmp3
      ENDDO

      ! start smoothing
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = DBLE(nz_l - 1)*L/DBLE(dfftp%nr3) - L
      z_r = DBLE(nz_r - 1)*L/DBLE(dfftp%nr3)
      f1 = -tpi*(z_r**2 + 2.d0*z_r*z0)*rg3 &
           - fpi*tmp1 &
           - fpi*ci*(z_r - z1)*tmp2 &
           - fpi*ci*(z1 - z0)*tmp3
      f2 = -tpi*(z_l**2 + 2.d0*z_l*z0)*rg3 &
           - fpi*tmp1 &
           - fpi*ci*(z_l - z1)*tmp2 &
           - fpi*ci*(z1 - z0)*tmp3
      f3 = -tpi*(2.d0*z_r + 2.d0*z0)*rg3 &
           - fpi*ci*tmp2
      f4 = -tpi*(2.d0*z_l + 2.d0*z0)*rg3 &
           - fpi*ci*tmp2
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
        z = DBLE(iz - 1)/DBLE(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
      ! end smoothing

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = (vg3(:, ng_2d) + vg(:))*e2 ! factor e2: hartree -> Ry.

      DEALLOCATE (vg, vg_r)
    ENDIF ! IF( ng_2d > 0 )

! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      eh = eh + sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
    ENDDO
    ehart = ehart + eh
    IF (gamma_only) THEN
      ehart = ehart*2d0
      ng_2d = imill_2d(0, 0)
      IF (ng_2d > 0) THEN
        ehart = ehart - sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
      ENDIF
    ENDIF
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum(ehart, intra_bgrp_comm)
    !
! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = vg3(n3,ng_2d)
    ENDDO

    DEALLOCATE (rhog3, vg3)

    RETURN
  END SUBROUTINE esm_hartree_bc3

  SUBROUTINE esm_hartree_bc4(rhog, ehart, aux)

    USE constants,        ONLY : pi, tpi, fpi, e2
    USE gvect,            ONLY : ngm, mill
    USE cell_base,        ONLY : omega, alat, tpiba2, at, bg
    USE control_flags,    ONLY : gamma_only
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE fft_base,         ONLY : dfftp
    USE fft_scalar,       ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    REAL(DP)                 :: ehart             !  Hartree energy
    COMPLEX(DP)              :: rhog(ngm)         !  n(G)
    COMPLEX(DP)              :: aux(dfftp%nnr)    !  v_h(G)
    !
    !    here the local variables
    !
    INTEGER                  :: k1, k2, k3, iz, igz, ng, n1, n2, n3, nz_r, nz_l, &
                                ng_2d
    REAL(DP)                 :: t(2), z, z0, z1, gp, gp2, kn, cc0, ss0, L, &
                                z_l, z_r, eh, aaa, cc1, ss1, alpha, beta, &
                                chi, xi, kappa, lambda, arg1, arg2, arg3, &
                                arg4, argr1, argr2, argr3, argr4, argr5
    COMPLEX(DP)              :: ci, f1, f2, f3, f4, a0, a1, a2, a3, c_r, c_l, &
                                s_r, s_l, rg3, tmp, tmp1, tmp2, tmp3, tmp4, &
                                tmpr1, tmpr2, tmpr3, tmpr4
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :), vg(:), vg_r(:), vg3(:, :), vr(:), &
                                vr_r(:)

    ALLOCATE (rhog3(dfftp%nr3, ngm_2d), vg3(dfftp%nr3, ngm_2d))
!
! Map to FFT mesh (dfftp%nr3,ngm_2d)
    rhog3(:, :) = (0.d0, 0.d0)
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      rg3 = rhog(ng)
      rhog3(n3, ng_2d) = rg3
      IF (gamma_only .and. n1 == 0 .and. n2 == 0) THEN
        n3 = -mill(3, ng) + 1
        IF (n3 < 1) n3 = n3 + dfftp%nr3
        rhog3(n3, ng_2d) = CONJG(rg3)
      ENDIF
    ENDDO
! End mapping
!
    vg3(:, :) = (0.d0, 0.d0)
    L = at(3, 3)*alat
    z0 = L/2.d0
    z1 = z0 + esm_w
    aaa = esm_a
    ci = (0.d0, 1.d0)

!****For gp!=0 case ********************
    ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3), vr(dfftp%nr3), vr_r(dfftp%nr3))
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      IF (k1 == 0 .and. k2 == 0) CYCLE
      t(1:2) = k1*bg(1:2, 1) + k2*bg(1:2, 2)
      gp2 = sum(t(:)*t(:))*tpiba2
      gp = sqrt(gp2)
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0); tmp3 = (0.d0, 0.d0)
      tmp4 = (0.d0, 0.d0); tmpr1 = (0.d0, 0.d0); tmpr2 = (0.d0, 0.d0)
      tmpr3 = (0.d0, 0.d0); tmpr4 = (0.d0, 0.d0)
      vr(:) = (0.d0, 0.d0); vg(:) = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        rg3 = rhog3(iz, ng_2d)
        ! bc4
        vg(iz) = fpi*rg3/(gp**2 + kn**2)
        vr(iz) = fpi*rg3/(gp**2 + kn**2 + ci*aaa*kn)
        cc1 = cos(kn*z1)
        ss1 = sin(kn*z1)
        alpha = aaa + gp + sqrt(aaa**2 + gp**2)
        beta = aaa + gp - sqrt(aaa**2 + gp**2)
        kappa = aaa - gp + sqrt(aaa**2 + gp**2)
        xi = aaa + sqrt(aaa**2 + gp**2)
        chi = aaa - sqrt(aaa**2 + gp**2)
        lambda = sqrt(aaa**2 + gp**2)
        tmp1 = tmp1 + rg3*(cc0 + ci*ss0)/(xi - ci*kn)/alpha
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/(gp + ci*kn)/gp
        tmp3 = tmp3 + rg3*kappa/alpha*(cc0 - ci*ss0)/(gp + ci*kn)/gp
        tmp4 = tmp4 + rg3*kappa*(cc1 + ci*ss1)/(xi - ci*kn)/(gp**2 + kn**2)
        tmpr1 = tmpr1 + rg3*(cc0 - ci*ss0)/(gp + ci*kn)/alpha
        tmpr2 = tmpr2 + rg3*(cc0 + ci*ss0)/(xi - ci*kn)/lambda
        tmpr3 = tmpr3 + rg3*beta/alpha*(cc0 + ci*ss0)/(xi - ci*kn)/lambda
        tmpr4 = tmpr4 + rg3*beta*(cc1 + ci*ss1)/(gp + ci*kn) &
                /(gp**2 + kn**2 + ci*2.d0*aaa*kn)
      ENDDO

      CALL cft_1z(vg, 1, dfftp%nr3, dfftp%nr3, 1, vg_r)
      ! bc4
      CALL cft_1z(vr, 1, dfftp%nr3, dfftp%nr3, 1, vr_r)

      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        ! bc4
        arg1 = gp*(z - z1) - xi*(z0 - z1)
        arg2 = -gp*(z + z0)
        arg3 = -gp*(z0 + z1) + gp*(z - z1)
        arg4 = gp*(z - z1)
        argr1 = -gp*(z0 + z1) - xi*(z - z1)
        argr2 = -xi*(z0 - z1) - chi*(z - z1)
        argr3 = -xi*(z - z1) - xi*(z0 - z1)
        argr4 = -xi*(z - z1)
        argr5 = -2.d0*aaa*(z - z1)
        IF (z < z1) THEN
          vg_r(iz) = vg_r(iz) - fpi*EXP(arg1)*tmp1 - tpi*EXP(arg2)*tmp2 &
                     + tpi*EXP(arg3)*tmp3 - fpi*EXP(arg4)*tmp4
        ELSE
          vg_r(iz) = vr_r(iz)*EXP(argr5) &
                     - fpi*EXP(argr1)*tmpr1 - tpi*EXP(argr2)*tmpr2 &
                     + tpi*EXP(argr3)*tmpr3 - fpi*EXP(argr4)*tmpr4
        ENDIF
      ENDDO

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)
      vg3(:, ng_2d) = vg(:)*e2 ! factor e2: hartree -> Ry.
    ENDDO
    DEALLOCATE (vg, vg_r, vr, vr_r)

!****For gp=0 case ********************
    ng_2d = imill_2d(0, 0)
    IF (ng_2d > 0) THEN
      ALLOCATE (vg(dfftp%nr3), vg_r(dfftp%nr3), vr(dfftp%nr3), vr_r(dfftp%nr3))
      tmp1 = (0.d0, 0.d0); tmp2 = (0.d0, 0.d0); tmp3 = (0.d0, 0.d0); tmp4 = (0.d0, 0.d0)
      vg(:) = (0.d0, 0.d0); vr(:) = (0.d0, 0.d0)
      !for smoothing
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      nz_l = dfftp%nr3/2 + 1 + esm_nfit
      nz_r = dfftp%nr3/2 + 1 - esm_nfit
      z_l = DBLE(nz_l - 1)*L/DBLE(dfftp%nr3) - L
      z_r = DBLE(nz_r - 1)*L/DBLE(dfftp%nr3)
      !
      rg3 = rhog3(1, ng_2d)
      ! bc4
      arg1 = -2.d0*aaa*(z0 - z1)
      vg(1) = tpi*((z0 + z1)/aaa + 2.d0*z0*z1 + z1**2)*rg3 &
              - pi*(EXP(arg1) - 1.d0)/aaa**2*rg3
      vr(1) = tpi*(z0 + 0.5d0/aaa)/aaa*rg3

      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        kn = DBLE(igz)*tpi/L
        rg3 = rhog3(iz, ng_2d)
        ! bc4
        cc0 = cos(kn*z0)
        ss0 = sin(kn*z0)
        cc1 = cos(kn*z1)
        ss1 = sin(kn*z1)
        tmp1 = tmp1 + rg3*(cc1 + ci*ss1)/(2.d0*aaa - ci*kn)/kn**2
        tmp2 = tmp2 + rg3*(cc0 - ci*ss0)/kn
        tmp3 = tmp3 + rg3*(cc0 + ci*ss0)/(2.d0*aaa - ci*kn)
        tmp4 = tmp4 + (0.d0, 0.d0)

        vg(iz) = fpi*rg3/(kn**2)
        ! bc4
        vr(iz) = fpi*rg3/(kn**2 + ci*2.d0*aaa*kn)

        !for smoothing
        c_r = cos(kn*z_r)
        s_r = sin(kn*z_r)
        c_l = cos(kn*z_l)
        s_l = sin(kn*z_l)
        ! bc4
        f1 = f1 + fpi*   rg3*(c_r + ci*s_r)/(kn**2 + ci*2.d0*aaa*kn)
        f2 = f2 + fpi*   rg3*(c_l + ci*s_l)/kn**2
        f3 = f3 + fpi*ci*rg3*(c_r + ci*s_r)/kn
        f4 = f4 + fpi*ci*rg3*(c_l + ci*s_l)/kn
        !
      ENDDO

      CALL cft_1z(vg, 1, dfftp%nr3, dfftp%nr3, 1, vg_r)
      ! bc4
      CALL cft_1z(vr, 1, dfftp%nr3, dfftp%nr3, 1, vr_r)

      rg3 = rhog3(1, ng_2d)
      DO iz = 1, dfftp%nr3
        k3 = iz - 1
        IF (k3 >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          k3 = k3 - dfftp%nr3
        END IF
        z = DBLE(k3) / DBLE(dfftp%nr3) * L
        ! bc4
        arg1 = -2.d0*aaa*(z0 - z1)
        arg2 = -2.d0*aaa*(z - z1)
        IF (z < z1) THEN
          vg_r(iz) = vg_r(iz) &
                     - fpi*2.d0*aaa*tmp1 &
                     - tpi*ci*(2.d0*(z - z1) - 1.d0/aaa)*tmp2 &
                     - tpi*EXP(arg1)/aaa*tmp3 &
                     - tpi*z*(z + 2.d0*z0)*rg3
        ELSE
          vg_r(iz) = vr_r(iz)*EXP(arg2) &
                     + tpi*ci*EXP(arg2)/aaa*tmp2 &
                     - tpi*EXP(arg1)/aaa*tmp3 &
                     + tpi*EXP(arg2)*z/aaa*rg3 &
                     - pi*EXP(arg1)/aaa**2*rg3
        ENDIF
      ENDDO

      !for smoothing
      ! bc4
      arg1 = -2.d0*aaa*(z0 - z1)
      arg2 = -2.d0*aaa*(z_r - z1)
      f1 = f1 + tpi*(z0 + 0.5d0/aaa)/aaa*rg3
      f1 = f1*EXP(arg2) &
           + tpi*ci*EXP(arg2)/aaa*tmp2 &
           - tpi*EXP(arg1)/aaa*tmp3 &
           + tpi*EXP(arg2)*z_r/aaa*rg3 &
           - pi*EXP(arg1)/aaa**2*rg3
      f2 = f2 + tpi*((z0 + z1)/aaa + 2.d0*z0*z1 + z1**2)*rg3 &
           - pi*(EXP(arg1) - 1.d0)/aaa**2*rg3
      f2 = f2 &
           - fpi*2.d0*aaa*tmp1 &
           - tpi*ci*(2.d0*(z_l - z1) - 1.d0/aaa)*tmp2 &
           - tpi*EXP(arg1)/aaa*tmp3 &
           - tpi*z_l*(z_l + 2.d0*z0)*rg3
      f3 = f3*EXP(arg2) &
           - fpi*ci*EXP(arg2)*tmp2 &
           - fpi*(z_r + z0)*EXP(arg2)*rg3
      f4 = f4 - fpi*ci*tmp2 - fpi*(z_l + z0)*rg3
      ! for smoothing
      !factor e2 will be multiplied later (at vg3 <= vg)
      !f1=f1*e2; f2=f2*e2; f3=f3*e2; f4=f4*e2
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
        z = DBLE(iz - 1)/DBLE(dfftp%nr3)*L
        vg_r(iz) = (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO

      CALL cft_1z(vg_r, 1, dfftp%nr3, dfftp%nr3, -1, vg)

      vg3(:, ng_2d) = vg(:)*e2 ! factor e2: hartree -> Ry.

      DEALLOCATE (vg, vg_r, vr, vr_r)
    ENDIF ! IF( ng_2d > 0 )

! Hartree Energy
    ehart = 0.d0
    eh = 0d0
    DO ng_2d = 1, ngm_2d
      k1 = mill_2d(1, ng_2d)
      k2 = mill_2d(2, ng_2d)
      eh = eh + sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
    ENDDO
    ehart = ehart + eh
    IF (gamma_only) THEN
      ehart = ehart*2d0
      ng_2d = imill_2d(0, 0)
      IF (ng_2d > 0) THEN
        ehart = ehart - sum(vg3(:, ng_2d)*conjg(rhog3(:, ng_2d)))
      ENDIF
    ENDIF
    ehart = ehart*omega*0.5d0
    !
    CALL mp_sum(ehart, intra_bgrp_comm)
    !
! Map to FFT mesh (dfftp%nrx)
    aux = 0.0d0
    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      ng_2d = imill_2d(n1, n2)
      n3 = mill(3, ng) + 1
      IF (n3 < 1) n3 = n3 + dfftp%nr3
      aux(ng) = vg3(n3,ng_2d)
    ENDDO

    DEALLOCATE (rhog3, vg3)

    RETURN
  END SUBROUTINE esm_hartree_bc4

END MODULE esm_hartree_mod
