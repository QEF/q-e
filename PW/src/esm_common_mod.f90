MODULE esm_common_mod

  USE kinds, ONLY : DP
  IMPLICIT NONE

  INTEGER                  :: esm_nfit
  REAL(DP)                 :: esm_efield, esm_w, esm_a
  CHARACTER(LEN=3)         :: esm_bc
  INTEGER, ALLOCATABLE     :: mill_2d(:, :), imill_2d(:, :)
  INTEGER                  :: ngm_2d = 0
  LOGICAL                  :: do_comp_esm = .FALSE.

CONTAINS

  SUBROUTINE esm_init(reject_charged_bc1)
    USE fft_base,  ONLY : dfftp
    USE constants, ONLY : tpi
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN), OPTIONAL :: reject_charged_bc1
    !
    INTEGER  :: iz, igz
    REAL(DP) :: phi
    !
    IF (PRESENT(reject_charged_bc1)) THEN
      CALL esm_check(reject_charged_bc1)
    ELSE
      CALL esm_check(.FALSE.)
    END IF
    !
    CALL esm_ggen_2d()
    !
  END SUBROUTINE esm_init

  SUBROUTINE esm_rgen_2d(dtau, rmax, mxr, at, bg, r, r2, nrm)
    !-----------------------------------------------------------------------
    !
    !   generates neighbours shells (cartesian, in units of lattice parameter)
    !   with length < rmax,and returns them in order of increasing length:
    !      r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:),   r2 = r^2
    !   where a1, a2, a3 are primitive lattice vectors. Other input variables:
    !     mxr = maximum number of vectors
    !     at  = lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
    !     bg  = reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
    !   Other output variables:
    !     nrm = the number of vectors with r^2 < rmax^2
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    INTEGER, INTENT(in)   :: mxr
    INTEGER, INTENT(out)  :: nrm
    REAL(DP), INTENT(in)  :: at(3, 3), bg(3, 3), dtau(3), rmax
    REAL(DP), INTENT(out) :: r(3, mxr), r2(mxr)
    !
    !    and here the local variables
    !
    INTEGER, ALLOCATABLE :: irr(:)
    INTEGER  ::  nm1, nm2, i, j, ipol, ir, indsw, iswap
    REAL(DP) :: ds(3), dtau0(3)
    REAL(DP) :: t(3), tt, swap
    REAL(DP), EXTERNAL :: dnrm2
    !
    !
    nrm = 0
    IF (rmax == 0.d0) RETURN

    ! bring dtau into the unit cell centered on the origin - prevents trouble
    ! if atomic positions are not centered around the origin but displaced
    ! far away (remember that translational invariance allows this!)
    !
    ds(:) = matmul(dtau(:), bg(:, :))
    ds(:) = ds(:) - anint(ds(:))
    dtau0(:) = matmul(at(:, :), ds(:))
    !
    ALLOCATE (irr(mxr))
    !
    ! these are estimates of the maximum values of needed integer indices
    !
    nm1 = int(dnrm2(3, bg(1, 1), 1)*rmax) + 2
    nm2 = int(dnrm2(3, bg(1, 2), 1)*rmax) + 2
    !
    DO i = -nm1, nm1
      DO j = -nm2, nm2
        tt = 0.d0
        DO ipol = 1, 3
          t(ipol) = i*at(ipol, 1) + j*at(ipol, 2) - dtau0(ipol)
          tt = tt + t(ipol)*t(ipol)
        ENDDO
        IF (tt <= rmax**2 .and. ABS(tt) > 1.d-10) THEN
          nrm = nrm + 1
          IF (nrm > mxr) CALL errore('esm_rgen_2d', 'too many r-vectors', nrm)
          DO ipol = 1, 3
            r(ipol, nrm) = t(ipol)
          ENDDO
          r2(nrm) = tt
        ENDIF
      ENDDO
    ENDDO
    !
    !   reorder the vectors in order of increasing magnitude
    !
    !   initialize the index inside sorting routine
    !
    irr(1) = 0
    IF (nrm > 1) CALL hpsort(nrm, r2, irr)
    DO ir = 1, nrm - 1
20    indsw = irr (ir)
      IF (indsw /= ir) THEN
        DO ipol = 1, 3
          swap = r(ipol, indsw)
          r(ipol, indsw) = r(ipol, irr(indsw))
          r(ipol, irr(indsw)) = swap
        ENDDO
        iswap = irr(ir)
        irr(ir) = irr(indsw)
        irr(indsw) = iswap
        GOTO 20
      ENDIF
    ENDDO
    DEALLOCATE (irr)
    !
    RETURN
  END SUBROUTINE esm_rgen_2d

  SUBROUTINE esm_ggen_2d()
    USE fft_base, ONLY: dfftp
    USE gvect, ONLY: ngm, mill
    !
    IMPLICIT NONE
    !
    INTEGER              :: n1xh, n2xh, ng, n1, n2, ng_2d
    LOGICAL, ALLOCATABLE  :: do_mill_2d(:, :)
    !
    !     Make g parallel array
    !
    n1xh = dfftp%nr1x/2
    n2xh = dfftp%nr2x/2
    ALLOCATE (do_mill_2d(-n1xh:n1xh, -n2xh:n2xh))
    do_mill_2d(:, :) = .false.

    DO ng = 1, ngm
      n1 = mill(1, ng)
      n2 = mill(2, ng)
      do_mill_2d(n1, n2) = .true.
    ENDDO
    ngm_2d = COUNT(do_mill_2d)
!*** do_mill_2d(h,k) = .true. means there is an h,k vector on this proc
!*** ngm_2d = total number of vectors (h,k) on this proc, excluding
!*** duplicates with different l values

    IF (ALLOCATED(mill_2d))   DEALLOCATE(mill_2d)
    IF (ALLOCATED(imill_2d))  DEALLOCATE(imill_2d)

    ALLOCATE (mill_2d(2, ngm_2d))
    ALLOCATE (imill_2d(-n1xh:n1xh, -n2xh:n2xh))

    mill_2d(:, :) = 0
    imill_2d(:, :) = 0
    ng_2d = 1
    DO n1 = -n1xh, n1xh
    DO n2 = -n2xh, n2xh
      IF (do_mill_2d(n1, n2)) THEN
        mill_2d(1, ng_2d) = n1
        mill_2d(2, ng_2d) = n2
        imill_2d(n1, n2) = ng_2d
        ng_2d = ng_2d+1
      ENDIF
    ENDDO
    ENDDO
    DEALLOCATE (do_mill_2d)
!**** mill_2d(:,ig) = h,k indices of vector ig
!**** imill_2d(h,k) = 2d index of vector with h,k indices
!**** ng_2d = total number of 2d g vectors on this proc

    RETURN
  END SUBROUTINE esm_ggen_2d

!-----------------------------------------------------------------------
!--------------ESM FINAL PRINTOUT SUBROUTINE----------------------------
!-----------------------------------------------------------------------
!
! Prints out vlocal and vhartree to stdout once electrons are converged
! Format: z, rho(r), v_hartree, v_local, (v_hartree + v_local)
!
  SUBROUTINE esm_printpot(rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill, igtongl
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE constants,     ONLY : AUTOEV, BOHR_RADIUS_ANGS
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp, ntyp => nsp
    USE vlocal,        ONLY : strf, vloc
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE io_files,      ONLY : prefix, tmp_dir
    IMPLICIT NONE

    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)

    INTEGER :: ig, iga, igb, igz, iz, jz, ia
    REAL(DP) :: L, S, z0, z1, z, gz, alpha, salp
    REAL(DP) :: Qa, ra(2), za
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3, vg3, sum1c, sum2c
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp

    COMPLEX(DP), ALLOCATABLE :: rho0r(:), rho0g(:), rho0g_tmp(:)
    COMPLEX(DP), ALLOCATABLE :: Vhar0r(:), Vhar0g(:)
    COMPLEX(DP), ALLOCATABLE :: Vloc0r(:), Vloc0g(:)
    CHARACTER(len=256)     :: esm1_file = 'os.esm1'

    IF (imill_2d(0, 0) == 0) RETURN

    !****For gp=0 case ********************

    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    alpha = 1.0d0
    salp = sqrt(alpha)

    ALLOCATE (rho0r(dfftp%nr3), rho0g(dfftp%nr3), rho0g_tmp(dfftp%nr3))
    ALLOCATE (Vhar0r(dfftp%nr3), Vhar0g(dfftp%nr3))
    ALLOCATE (Vloc0r(dfftp%nr3), Vloc0g(dfftp%nr3))

    rho0g(:) = (0.d0, 0.d0)

    !!---- calculate density potential
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1

      IF (.not. (iga == 0 .and. igb == 0)) CYCLE

      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      rg3 = rhog(ig)
      rho0g(igz) = rg3

      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rho0g(igz) = CONJG(rg3)
      END IF
    END DO ! ig

    rho0g_tmp = rho0g
    CALL cft_1z(rho0g_tmp, 1, dfftp%nr3, dfftp%nr3, +1, rho0r)

    !!---- calculate hartree potential
    Vhar0g(:) = 0.0d0
    DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
      IF (igz == 0) CYCLE
      iz = igz + 1
      IF (iz < 1) iz = iz + dfftp%nr3
      gz = dble(igz)*tpi/L

      rg3 = rho0g(iz)
      Vhar0g(iz) = fpi*rg3/gz**2
    END DO ! igz

    CALL cft_1z(Vhar0g, 1, dfftp%nr3, dfftp%nr3, +1, Vhar0r)

    ! summations over gz
    sum1c = (0.d0, 0.d0)
    sum2c = (0.d0, 0.d0)
    DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
      IF (igz == 0) CYCLE
      iz = igz + 1
      IF (iz < 1) iz = iz + dfftp%nr3
      gz = dble(igz)*tpi/L

      rg3 = rho0g(iz)
      sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
      sum2c = sum2c + rg3*cos(gz*z0)/gz**2
    END DO ! igz

    rg3 = rho0g(1)
    DO iz = 1, dfftp%nr3
      jz = iz - 1
      IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
        jz = jz - dfftp%nr3
      END IF
      z = DBLE(jz) / DBLE(dfftp%nr3) * L

      IF (do_comp_esm .AND. ( esm_bc /= 'pbc' )) THEN
        !! BC1 terms
        Vhar0r(iz) = Vhar0r(iz) &
                     - tpi*z**2*rg3 &
                     - tpi*z0**2*rg3 &
                     - fpi*z*sum1c &
                     - fpi*sum2c
        IF (esm_bc == 'bc2') THEN
          !! BC2 terms
          Vhar0r(iz) = Vhar0r(iz) &
                       + tpi*z1*2*z0*rg3 - tpi*(-z/z1)*2*z0*sum1c
        ELSE IF (esm_bc == 'bc3') THEN
          !! BC3 terms
          Vhar0r(iz) = Vhar0r(iz) &
                       - tpi*(z - 2*z1)*2*z0*rg3 + fpi*z0*sum1c
        END IF
      END IF
    END DO ! iz

    !!---- calculate local potential
    ! short range
    Vloc0g(:) = 0.0d0
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)

      IF (.not. (iga == 0 .and. igb == 0)) CYCLE

      igz = mill(3, ig) + 1
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF

      vg3 = 0.0d0
      DO ia = 1, ntyp
        vg3 = vg3 + vloc(igtongl(ig), ia)*strf(ig, ia)/e2
      END DO
      Vloc0g(igz) = vg3

      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF

        Vloc0g(igz) = CONJG(vg3)
      END IF
    END DO

    CALL cft_1z(Vloc0g, 1, dfftp%nr3, dfftp%nr3, +1, Vloc0r)

    ! long range
    DO iz = 1, dfftp%nr3
      jz = iz - 1
      IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
        jz = jz - dfftp%nr3
      END IF
      z = DBLE(jz) / DBLE(dfftp%nr3) * L

      IF (esm_bc == 'bc2') THEN
        Vloc0r(iz) = Vloc0r(iz) + (z1 - z)*esm_efield/e2
      END IF

      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF

        IF (do_comp_esm .AND. ( esm_bc /= 'pbc' )) THEN
          !! BC1 terms
          Vloc0r(iz) = Vloc0r(iz) - tpi*Qa/S &
                       *((z - za)*erf(salp*(z - za)) &
                       + exp(-alpha*(z - za)**2)*sqrtpm1/salp)
          IF (esm_bc == 'bc2') THEN
            !! BC2 terms
            Vloc0r(iz) = Vloc0r(iz) &
                         + tpi*Qa/S*(-z*za + z1*z1)/z1
          ELSE IF (esm_bc == 'bc3') THEN
            !! BC3 terms
            Vloc0r(iz) = Vloc0r(iz) &
                         + tpi*Qa/S*(-z + 2*z1 - za)
          END IF
        END IF

      END DO ! ia
    END DO ! iz

    !!---- output potentials
    esm1_file = TRIM(tmp_dir)//TRIM(prefix)//".esm1"
    open (UNIT=4, FILE=esm1_file, STATUS="UNKNOWN", &
          ACTION="WRITE")
    write (UNIT=4, FMT=9050)

    DO iz = (dfftp%nr3 - dfftp%nr3/2 + 1), dfftp%nr3
      z = DBLE(iz - 1 - dfftp%nr3) / DBLE(dfftp%nr3) * L

      write (UNIT=4, FMT=9051) &
        z*BOHR_RADIUS_ANGS, &
        REAL(rho0r(iz))*S/BOHR_RADIUS_ANGS, &
        REAL(Vhar0r(iz))*AUTOEV, &
        REAL(Vloc0r(iz))*AUTOEV, &
        REAL(Vhar0r(iz))*AUTOEV + REAL(Vloc0r(iz))*AUTOEV
    END DO

    DO iz = 1, (dfftp%nr3 - dfftp%nr3/2)
      z = DBLE(iz - 1) / DBLE(dfftp%nr3) * L

      write (UNIT=4, FMT=9051) &
        z*BOHR_RADIUS_ANGS, &
        REAL(rho0r(iz))*S/BOHR_RADIUS_ANGS, &
        REAL(Vhar0r(iz))*AUTOEV, &
        REAL(Vloc0r(iz))*AUTOEV, &
        REAL(Vhar0r(iz))*AUTOEV + REAL(Vloc0r(iz))*AUTOEV
    END DO
    close (UNIT=4)

    DEALLOCATE (rho0r, rho0g, rho0g_tmp)
    DEALLOCATE (Vhar0r, Vhar0g)
    DEALLOCATE (Vloc0r, Vloc0g)

    RETURN

9050 FORMAT( '#z (A)',2X,'Tot chg (e/A)',2X,'Avg v_hartree (eV)',2X,&
            &'Avg v_local (eV)',2x,'Avg v_hart+v_loc (eV)' )
9051 FORMAT( F6.2,F20.7,F20.7,F18.7,F18.7 )
  END SUBROUTINE esm_printpot
!
!-----------------------------------------------------------------------
!--------------ESM SUMMARY PRINTOUT SUBROUTINE--------------------------
!-----------------------------------------------------------------------
!
! Prints summary of ESM parameters to stdout
!
  SUBROUTINE esm_summary()
    !
    USE io_global, ONLY : stdout, ionode
    USE constants, ONLY : rytoev, BOHR_RADIUS_ANGS
    USE klist,     ONLY : tot_charge
    !
    IMPLICIT NONE
    !
    IF (.NOT. ionode) RETURN
    !
    WRITE (UNIT = stdout, &
           FMT  = '(/,5x, "Effective Screening Medium Method",     &
           &/,5x, "=================================")' )
    !
    SELECT CASE(TRIM(esm_bc))
    !
    CASE ('pbc')
      WRITE (UNIT = stdout, &
             FMT = '(5x, "Ordinary Periodic Boundary Conditions")')
    CASE ('bc1')
      WRITE (UNIT = stdout, &
             FMT = '(5x, "Boundary Conditions: Vacuum-Slab-Vacuum")')
    CASE ('bc2')
      WRITE (UNIT = stdout, &
             FMT = '(5x, "Boundary Conditions: Metal-Slab-Metal")')
    CASE ('bc3')
      WRITE (UNIT = stdout, &
             FMT = '(5x, "Boundary Conditions: Vacuum-Slab-Metal")')
    CASE ('bc4')
      WRITE (UNIT = stdout, &
             FMT = '(5x, "Boundary Conditions: Vacuum-Slab-smooth ESM)")')
    END SELECT
    !
    WRITE(UNIT = stdout, FMT = 9055) tot_charge
    IF (esm_efield /= 0.0_DP) THEN
      WRITE (UNIT = stdout, FMT = 9051) esm_efield
    END IF
    !
    IF (esm_w /= 0.0_DP) THEN
      WRITE (UNIT = stdout, FMT = 9052) esm_w*BOHR_RADIUS_ANGS, esm_w
    END IF
    !
    IF (esm_bc .EQ. 'bc4') THEN
      WRITE (UNIT = stdout, FMT = 9054) esm_a
    ENDIF
    !
    WRITE (UNIT = stdout, FMT = 9053) esm_nfit
    !
    WRITE (stdout, *)
    !
9051 FORMAT( '     field strength                   = ', F8.4,' Ry/a.u.')
9052 FORMAT( '     ESM offset from cell edge        = ', F8.2,' A' &
            /'                                      = ', F8.2,' a.u.')
9053 FORMAT( '     grid points for fit at edges     = ', I8,' ')
9054 FORMAT( '     smoothness parameter             = ', F8.2,' 1/a.u.' )
9055 FORMAT( '     total charge in unit cell        = ', F8.4)

  END SUBROUTINE esm_summary
!
!-----------------------------------------------------------------------
!--------------ESM CHECK CONDITIONS SUBROUTINE--------------------------
!-----------------------------------------------------------------------
!
! Checks conditions to perform ESM
!
  SUBROUTINE esm_check(reject_charged_bc1)
    !
    USE cell_base, ONLY : at, iforceh
    USE cellmd,    ONLY : lmovecell
    USE constants, ONLY : eps14
    USE exx_base,  ONLY : x_gamma_extrapolation
    USE xc_lib,    ONLY : exx_is_active
    USE ions_base, ONLY : nat, tau
    USE klist,     ONLY : nkstot, xk, tot_charge
    USE lsda_mod,  ONLY : lsda
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: reject_charged_bc1
    !
    INTEGER :: ia
    INTEGER :: ik, nk_
    !
    ! ... correct cell shape ?
    IF (ABS(at(1, 3)) > eps14 .OR. ABS(at(3, 1)) > eps14 .OR. &
      & ABS(at(2, 3)) > eps14 .OR. ABS(at(3, 2)) > eps14) THEN
      CALL errore('esm_check', 'incorrect unit cell for ESM', 1)
    END IF
    !
    ! ... correct atomic positions ?
    DO ia = 1, nat
      IF (tau(3, ia) <= (-0.5_DP * at(3, 3)) .OR. (0.5_DP * at(3, 3)) <= tau(3, ia)) THEN
        CALL errore('esm_check', 'incorrect atomic position for ESM', ia)
      END IF
    END DO
    !
    ! ... correct k-points ?
    IF (lsda) THEN
      nk_ = nkstot / 2
    ELSE
      nk_ = nkstot
    END IF
    !
    DO ik = 1, nk_
      IF (ABS(xk(3, ik)) > eps14) THEN
        CALL errore('esm_check', 'incorrect k-point for ESM', ik)
      END IF
    END DO
    !
    ! ... correct Vexx(G=0) ?
    IF (exx_is_active() .AND. (.NOT. x_gamma_extrapolation)) THEN
      CALL errore('esm_check', 'ESM requires Vexx(G=0)', 1)
    END IF
    !
    IF (reject_charged_bc1) THEN
      ! ... cannot have charge, if BC1
      IF (TRIM(esm_bc) == 'bc1' .AND. ABS(tot_charge) > eps14) THEN
        CALL errore('esm_check', 'cannot have charge, when ESM-BC1', 1)
      END IF
    END IF
    !
    ! ... vc-relax, only supports 2Dxy
    IF (lmovecell) THEN
      IF (iforceh(3, 1) /= 0 .OR. iforceh(3, 2) /= 0 .OR. iforceh(3, 3) /= 0 .OR. &
        & iforceh(1, 3) /= 0 .OR. iforceh(2, 3) /= 0) THEN
        CALL errore('esm_check', 'ESM only supports cell_dofree = "2Dxy"', 1)
      END IF
    END IF
    !
  END SUBROUTINE esm_check

! exp(x) * erfc(y)
! This function is to avoid INFINITY * ZERO for large positive x and y.
  REAL(8) FUNCTION exp_erfc(x, y)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x, y
    !REAL(8)             :: ym, ym2, nume, deno
    !REAL(8), PARAMETER  :: rtpim = 0.564189583547756279d0 ! 1 / sqrt(PI)
    !REAL(8), PARAMETER  :: r(0:4) = (/ &
    !                       -2.99610707703542174d-3, -4.94730910623250734d-2, &
    !                       -2.26956593539686930d-1, -2.78661308609647788d-1, &
    !                       -2.23192459734184686d-2/)
    !REAL(8), PARAMETER  :: s(0:4) = (/ &
    !                       1.06209230528467918d-2, 1.91308926107829841d-1, &
    !                       1.05167510706793207d0, 1.98733201817135256d0, &
    !                       1.00000000000000000d0/)
    !
    !IF (x < 709.0d0 .or. y < 4.0d0) THEN
    !  exp_erfc = exp(x)*erfc(y)
    !ELSE
    !  ym = 1d0/y
    !  ym2 = ym**2
    !  nume = (((r(4)*ym2 + r(3))*ym2 + r(2))*ym2 + r(1))*ym2 + r(0)
    !  deno = (((s(4)*ym2 + s(3))*ym2 + s(2))*ym2 + s(1))*ym2 + s(0)
    !  exp_erfc = exp(-y**2 + x)*ym*(rtpim + ym2*nume/deno)
    !END IF

    exp_erfc = exp(x + log(erfc(y)))

    RETURN
  END FUNCTION exp_erfc

  SUBROUTINE qromb(func, aaa, tmp, z1, z, zp, rxy, b, ss)
    USE kinds, ONLY : DP
    INTEGER, PARAMETER    ::jmax = 20, jmaxp = jmax + 1, k = 5, km = k - 1
    INTEGER               ::j
    REAL(DP), INTENT(in)  ::aaa, tmp, z1, z, zp, rxy, b
    REAL(DP), INTENT(out) ::ss
    REAL(DP), PARAMETER   ::a = 0.0_DP, eps = 1.e-12
    REAL(DP), EXTERNAL    ::func
    REAL(DP)              ::dss, h(jmaxp), s(jmaxp)
    !
    ! ss=int_a^b func(gp,aaa,tmp,z1,z,zp,rxy) dgp
    !
    h(1) = 1.0_DP
    DO j = 1, jmax
      CALL trapzd(func, aaa, tmp, z1, z, zp, rxy, a, b, s(j), j)
      IF (j .ge. k) THEN
        CALL polint(h(j - km), s(j - km), k, 0.0_DP, ss, dss)
        IF (ABS(ss) .le. 1.e-8) RETURN
        IF (ABS(dss) .le. eps*ABS(ss)) RETURN
      ENDIF
      s(j + 1) = s(j)
      h(j + 1) = 0.25*h(j)
    ENDDO
    STOP 'too many steps in qromb'
  END SUBROUTINE qromb

  SUBROUTINE trapzd(func, aaa, tmp, z1, z, zp, rxy, a, b, s, n)
    USE kinds, ONLY : DP
    INTEGER, INTENT(in)::n
    REAL(DP), INTENT(in)::aaa, tmp, z1, z, zp, rxy, a, b
    REAL(DP), INTENT(inout)::s
    REAL(DP), EXTERNAL ::func
    INTEGER::it, j
    REAL(DP)::del, sum, tnm, x

    IF (n .eq. 1) THEN
      s = 0.5*(b - a)*(func(a, aaa, tmp, z1, z, zp, rxy) + func(b, aaa, tmp, z1, z, zp, rxy))
    ELSE
      it = 2**(n - 2)
      tnm = it
      del = (b - a)/tnm
      x = a + 0.5*del
      sum = 0.
      DO j = 1, it
        sum = sum + func(x, aaa, tmp, z1, z, zp, rxy)
        x = x + del
      ENDDO
      s = 0.5*(s + (b - a)*sum/tnm)
    ENDIF
    RETURN
  END SUBROUTINE trapzd

  FUNCTION vl11j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY: DP
    IMPLICIT NONE
    REAL(DP)             :: vl11j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                           arg003, arg006, arg101, arg102, arg103, arg104, &
                           arg109
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi = aaa + sqrt(aaa**2 + gp**2)
    arg001 = gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 = gp*(z + zp - 2.d0*z1)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg101 = gp/2.d0/tmp + tmp*(z - zp)
    arg102 = gp/2.d0/tmp - tmp*(z - zp)
    arg103 = gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 = gp/2.d0/tmp - tmp*(z1 - zp)
    arg109 = xi/2.d0/tmp + tmp*(z1 - zp)
    t1 = -exp(arg003)*kappa/alpha
    t2 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103) &
         + exp_erfc(arg002, arg102) - kappa/alpha*exp_erfc(arg003, arg104)
    t3 = exp_erfc(arg006, arg109)/alpha
    vl11j0 = (t1/2.d0 - (t2/4.d0 + gp*t3/2.d0))*dbesj0(gp*rxy)
    !
    RETURN
  END FUNCTION vl11j0

  FUNCTION vl11j1(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl11j1
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                           arg003, arg006, arg007, arg101, arg102, arg103, arg104, &
                           arg105, arg106, arg109, arg110
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi = aaa + sqrt(aaa**2 + gp**2)
    arg001 = gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 = gp*(z + zp - 2.d0*z1)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg007 = aaa/2.d0/tmp**2*xi - gp*(z1 - zp) - xi*(z - z1)
    arg101 = gp/2.d0/tmp + tmp*(z - zp)
    arg102 = gp/2.d0/tmp - tmp*(z - zp)
    arg103 = gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 = gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 = gp/2.d0/tmp + tmp*(z - z1)
    arg106 = gp/2.d0/tmp - tmp*(z - z1)
    arg109 = xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 = xi/2.d0/tmp - tmp*(z - z1)
    t1 = -exp(arg003)*kappa/alpha
    t2 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103) &
         + exp_erfc(arg002, arg102) - kappa/alpha*exp_erfc(arg003, arg104) &
         + exp_erfc(arg006, arg109)*2.d0*gp/alpha
    t3 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106) &
         + exp_erfc(arg001, arg101) - kappa/alpha*exp_erfc(arg003, arg105) &
         + exp_erfc(arg007, arg110)*2.d0*gp/alpha
    vl11j1 = gp*(t1 - (t2 + t3)/4.d0)*dbesj1(gp*rxy)
    !
    RETURN
  END FUNCTION vl11j1

  FUNCTION vl12j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl12j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                           arg003, arg004, arg006, arg101, arg102, arg103, &
                           arg104, arg109
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi = aaa + sqrt(aaa**2 + gp**2)
    arg001 = gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 = gp*(z + zp - 2.d0*z1)
    arg004 = gp*(z - z1) + xi*(z1 - zp)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg101 = gp/2.d0/tmp + tmp*(z - zp)
    arg102 = gp/2.d0/tmp - tmp*(z - zp)
    arg103 = gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 = gp/2.d0/tmp - tmp*(z1 - zp)
    arg109 = xi/2.d0/tmp + tmp*(z1 - zp)
    t1 = exp(arg004)/alpha
    t2 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103) &
         + exp_erfc(arg002, arg102) - kappa/alpha*exp_erfc(arg003, arg104)
    t3 = exp_erfc(arg006, arg109)/alpha
    vl12j0 = (gp*t1 - (t2/4.d0 + gp*t3/2.d0))*dbesj0(gp*rxy)
  END FUNCTION vl12j0

  FUNCTION vl12j1(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl12j1
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, kappa, xi, chi, lambda, t1, t2, t3, &
                           arg001, arg002, arg003, arg004, arg006, arg009, &
                           arg010, arg012, arg101, arg102, arg103, arg104, &
                           arg105, arg108, arg109, arg110, arg112, arg114
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg004 =  gp*(z - z1) + xi*(z1 - zp)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg010 = aaa/2.d0/tmp**2*xi + chi*(z1 - zp) - xi*(z - z1)
    arg012 = aaa/2.d0/tmp**2*chi + xi*(z1 - zp) - chi*(z - z1)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg103 =  gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg108 =  xi/2.d0/tmp - tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg112 = chi/2.d0/tmp - tmp*(z - zp)
    arg114 = chi/2.d0/tmp - tmp*(z - z1)
    t1 = exp(arg004)/alpha
    t2 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103) &
         + exp_erfc(arg002, arg102) - kappa/alpha*exp_erfc(arg003, arg104) &
         + exp_erfc(arg006, arg109)*2.d0*gp/alpha
    t3 = exp_erfc(arg012, arg114) - exp_erfc(arg012, arg112) &
         + exp_erfc(arg010, arg108) - beta/alpha*exp_erfc(arg009, arg110) &
         + exp_erfc(arg004, arg105)*2.d0*lambda/alpha
    vl12j1 = (2.d0*t1*gp**2 - (gp*t2 + gp**2*t3/lambda)/4.d0)*dbesj1(gp*rxy)
  END FUNCTION vl12j1

  FUNCTION vl21j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl21j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg005, &
                           arg008, arg009, arg011, arg104, arg107, arg109, &
                           arg111, arg113
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    t1 = exp(arg005)/alpha
    t2 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111) &
         + exp_erfc(arg008, arg107) - beta/alpha*exp_erfc(arg009, arg109)
    t3 = exp_erfc(arg005, arg104)/alpha
    vl21j0 = gp*(t1 - (t2/4.d0/lambda + t3/2.d0))*dbesj0(gp*rxy)
  END FUNCTION vl21j0

  FUNCTION vl21j1(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl21j1
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, kappa, xi, chi, lambda, t1, t2, t3, &
                           arg001, arg002, arg003, arg005, arg007, arg008, &
                           arg009, arg011, arg101, arg102, arg104, arg105, &
                           arg106, arg107, arg109, arg110, arg111, arg113
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg007 = aaa/2.d0/tmp**2*xi - gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg106 =  gp/2.d0/tmp - tmp*(z - z1)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    t1 = exp(arg005)/alpha
    t2 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111) &
         + exp_erfc(arg008, arg107) - beta/alpha*exp_erfc(arg009, arg109) &
         + exp_erfc(arg005, arg104)*2.d0*lambda/alpha
    t3 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106) &
         + exp_erfc(arg001, arg101) - kappa/alpha*exp_erfc(arg003, arg105) &
         + exp_erfc(arg007, arg110)*2.d0*gp/alpha
    vl21j1 = (2.d0*t1*gp**2 - (gp**2*t2/lambda + gp*t3)/4.d0)*dbesj1(gp*rxy)
  END FUNCTION vl21j1

  FUNCTION vl22j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl22j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                           arg005, arg008, arg009, arg011, arg104, arg107, &
                           arg109, arg111, arg113
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg000 = -xi*(z + zp - 2.d0*z1)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    t1 = -exp(arg000)*beta/alpha
    t2 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111) &
         + exp_erfc(arg008, arg107) - beta/alpha*exp_erfc(arg009, arg109)
    t3 = exp_erfc(arg005, arg104)/alpha
    vl22j0 = gp*(t1/2.d0/lambda - (t2/4.d0/lambda + t3/2.d0))*dbesj0(gp*rxy)
    !
    RETURN
  END FUNCTION vl22j0

  FUNCTION vl22j1(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl22j1
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                           arg004, arg005, arg008, arg009, arg010, arg011, &
                           arg012, arg104, arg105, arg107, arg108, arg109, &
                           arg110, arg111, arg112, arg113, arg114
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg000 = -xi*(z + zp - 2.d0*z1)
    arg004 =  gp*(z - z1) + xi*(z1 - zp)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg010 = aaa/2.d0/tmp**2*xi + chi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg012 = aaa/2.d0/tmp**2*chi + xi*(z1 - zp) - chi*(z - z1)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg108 =  xi/2.d0/tmp - tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg112 = chi/2.d0/tmp - tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    arg114 = chi/2.d0/tmp - tmp*(z - z1)
    t1 = -exp(arg000)*beta/alpha
    t2 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111) &
         + exp_erfc(arg008, arg107) - beta/alpha*exp_erfc(arg009, arg109) &
         + exp_erfc(arg005, arg104)*2.d0*lambda/alpha
    t3 = exp_erfc(arg012, arg114) - exp_erfc(arg012, arg112) &
         + exp_erfc(arg010, arg108) - beta/alpha*exp_erfc(arg009, arg110) &
         + exp_erfc(arg004, arg105)*2.d0*lambda/alpha
    vl22j1 = gp**2*(t1 - (t2 + t3)/4.d0)*dbesj1(gp*rxy)/lambda
    !
    RETURN
  END FUNCTION vl22j1

  FUNCTION vl11(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl11
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                           arg003, arg006, arg101, arg102, arg103, arg104, &
                           arg109
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi = aaa + sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg103 =  gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    t1 = -exp(arg003)*kappa/alpha
    t2 = exp_erfc(arg001, arg101) - exp_erfc(arg001, arg103) &
         + exp_erfc(arg002, arg102) - kappa/alpha*exp_erfc(arg003, arg104)
    t3 = exp_erfc(arg006, arg109)/alpha
    vl11 = t1/2.d0 - (t2/4.d0 + gp*t3/2.d0)
    !
    RETURN
  END FUNCTION vl11

  FUNCTION vl22(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: vl22
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, xi, chi, lambda, t1, t2, t3, arg000, &
                           arg005, arg008, arg009, arg011, arg104, arg107, &
                           arg109, arg111, arg113
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg000 = -xi*(z + zp - 2.d0*z1)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    t1 = -exp(arg000)*beta/alpha
    t2 = exp_erfc(arg011, arg113) - exp_erfc(arg011, arg111) &
         + exp_erfc(arg008, arg107) - beta/alpha*exp_erfc(arg009, arg109)
    t3 = exp_erfc(arg005, arg104)/alpha
    vl22 = gp*t1/2.d0/lambda - gp*(t2/4.d0/lambda + t3/2.d0)
    !
    RETURN
  END FUNCTION vl22

  FUNCTION dvl11(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: dvl11
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                           arg003, arg006, arg007, arg101, arg102, arg103, &
                           arg104, arg105, arg106, arg109, arg110
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg007 = aaa/2.d0/tmp**2*xi - gp*(z1 - zp) - xi*(z - z1)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg103 =  gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg106 =  gp/2.d0/tmp - tmp*(z - z1)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    t1 = -exp(arg003)*kappa/alpha
    t2 = exp_erfc(arg001, arg103) - exp_erfc(arg001, arg101) &
         + exp_erfc(arg002, arg102) - exp_erfc(arg003, arg104)*kappa/alpha &
         - exp_erfc(arg006, arg109)*xi/alpha*2.d0
    t3 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106) &
         - exp_erfc(arg001, arg101) - exp_erfc(arg003, arg105)*kappa/alpha &
         + exp_erfc(arg007, arg110)*gp/alpha*2.d0
    dvl11 = gp*(t1 - (t2 + t3)/4.d0)
    !
    RETURN
  END FUNCTION dvl11

  FUNCTION dvl22(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: dvl22
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, kappa, xi, chi, lambda, arg000, &
                           arg004, arg005, arg008, arg009, arg010, arg011, &
                           arg012, arg104, arg105, arg107, arg108, arg109, &
                           arg110, arg111, arg112, arg113, arg114, t1, t2, t3
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg000 = -xi*(z + zp - 2.d0*z1)
    arg004 =  gp*(z - z1) + xi*(z1 - zp)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg010 = aaa/2.d0/tmp**2*xi + chi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg012 = aaa/2.d0/tmp**2*chi + xi*(z1 - zp) - chi*(z - z1)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg108 =  xi/2.d0/tmp - tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg112 = chi/2.d0/tmp - tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    arg114 = chi/2.d0/tmp - tmp*(z - z1)
    t1 = exp(arg000)*beta*xi/alpha/lambda
    t2 = (exp_erfc(arg011, arg111) - exp_erfc(arg011, arg113))*chi/lambda &
         - exp_erfc(arg008, arg107)*xi/lambda &
         + exp_erfc(arg009, arg109)*xi*beta/alpha/lambda &
         + exp_erfc(arg005, arg104)*gp/alpha*2.d0
    t3 = (exp_erfc(arg012, arg112) - exp_erfc(arg012, arg114))*xi/lambda &
         - exp_erfc(arg010, arg108)*chi/lambda &
         + exp_erfc(arg009, arg110)*xi*beta/alpha/lambda &
         - exp_erfc(arg004, arg105)*xi/alpha*2.d0
    dvl22 = gp*(t1 - (t2 + t3)/4.d0)
    !
    RETURN
  END FUNCTION dvl22

  FUNCTION dvl11j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: dvl11j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, kappa, xi, t1, t2, t3, arg001, arg002, &
                           arg003, arg006, arg007, arg101, arg102, arg103, &
                           arg104, arg105, arg106, arg109, arg110
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg007 = aaa/2.d0/tmp**2*xi - gp*(z1 - zp) - xi*(z - z1)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg103 =  gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg106 =  gp/2.d0/tmp - tmp*(z - z1)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    t1 = -exp(arg003)*kappa/alpha
    t2 = exp_erfc(arg001, arg103) - exp_erfc(arg001, arg101) &
         + exp_erfc(arg002, arg102) - exp_erfc(arg003, arg104)*kappa/alpha &
         - exp_erfc(arg006, arg109)*xi/alpha*2.d0
    t3 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106) &
         - exp_erfc(arg001, arg101) - exp_erfc(arg003, arg105)*kappa/alpha &
         + exp_erfc(arg007, arg110)*gp/alpha*2.d0
    dvl11j0 = gp*(t1 - (t2 + t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
  END FUNCTION dvl11j0

  FUNCTION dvl12j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: dvl12j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, kappa, xi, chi, lambda, arg001, &
                           arg002, arg003, arg004, arg006, arg009, arg010, &
                           arg012, arg101, arg102, arg103, arg104, arg105, &
                           arg108, arg109, arg110, arg112, arg114, t1, t2, t3
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg004 =  gp*(z - z1) + xi*(z1 - zp)
    arg006 = aaa/2.d0/tmp**2*xi + gp*(z - z1) + xi*(z1 - zp)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg010 = aaa/2.d0/tmp**2*xi + chi*(z1 - zp) - xi*(z - z1)
    arg012 = aaa/2.d0/tmp**2*chi + xi*(z1 - zp) - chi*(z - z1)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg103 =  gp/2.d0/tmp + tmp*(z1 - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg108 =  xi/2.d0/tmp - tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg112 = chi/2.d0/tmp - tmp*(z - zp)
    arg114 = chi/2.d0/tmp - tmp*(z - z1)
    t1 = -exp(arg004)*xi/alpha
    t2 = exp_erfc(arg001, arg103) - exp_erfc(arg001, arg101) &
         + exp_erfc(arg002, arg102) - exp_erfc(arg003, arg104)*kappa/alpha &
         - exp_erfc(arg006, arg109)*xi/alpha*2.d0
    t3 = (exp_erfc(arg012, arg112) - exp_erfc(arg012, arg114))*xi/lambda &
         - exp_erfc(arg010, arg108)*chi/lambda &
         + exp_erfc(arg009, arg110)*xi*beta/alpha/lambda &
         - exp_erfc(arg004, arg105)*xi/alpha*2.d0
    dvl12j0 = gp*(2.d0*t1 - (t2 + t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
  END FUNCTION dvl12j0

  FUNCTION dvl21j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: dvl21j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, kappa, xi, chi, lambda, arg001, &
                           arg002, arg003, arg005, arg007, arg008, arg009, &
                           arg011, arg101, arg102, arg104, arg105, arg106, &
                           arg107, arg109, arg110, arg111, arg113, t1, t2, t3
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg001 =  gp*(z - zp)
    arg002 = -gp*(z - zp)
    arg003 =  gp*(z + zp - 2.d0*z1)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg007 = aaa/2.d0/tmp**2*xi - gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg101 =  gp/2.d0/tmp + tmp*(z - zp)
    arg102 =  gp/2.d0/tmp - tmp*(z - zp)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg106 =  gp/2.d0/tmp - tmp*(z - z1)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    t1 = exp(arg005)*gp/alpha
    t2 = (exp_erfc(arg011, arg111) - exp_erfc(arg011, arg113))*chi/lambda &
         - exp_erfc(arg008, arg107)*xi/lambda &
         + exp_erfc(arg009, arg109)*xi*beta/alpha/lambda &
         + exp_erfc(arg005, arg104)*gp/alpha*2.d0
    t3 = exp_erfc(arg002, arg102) - exp_erfc(arg002, arg106) &
         - exp_erfc(arg001, arg101) - exp_erfc(arg003, arg105)*kappa/alpha &
         + exp_erfc(arg007, arg110)*gp/alpha*2.d0
    dvl21j0 = gp*(2.d0*t1 - (t2 + t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
  END FUNCTION dvl21j0

  FUNCTION dvl22j0(gp, aaa, tmp, z1, z, zp, rxy)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP)             :: dvl22j0
    REAL(DP), INTENT(in) :: gp, aaa, tmp, z1, z, zp, rxy
    !local
    REAL(DP)             ::alpha, beta, kappa, xi, chi, lambda, arg000, &
                           arg004, arg005, arg008, arg009, arg010, arg011, &
                           arg012, arg104, arg105, arg107, arg108, arg109, &
                           arg110, arg111, arg112, arg113, arg114, t1, t2, t3
    !
    alpha = aaa + gp + sqrt(aaa**2 + gp**2)
    beta  = aaa + gp - sqrt(aaa**2 + gp**2)
    kappa = aaa - gp + sqrt(aaa**2 + gp**2)
    xi    = aaa + sqrt(aaa**2 + gp**2)
    chi   = aaa - sqrt(aaa**2 + gp**2)
    lambda = sqrt(aaa**2 + gp**2)
    arg000 = -xi*(z + zp - 2.d0*z1)
    arg004 =  gp*(z - z1) + xi*(z1 - zp)
    arg005 = -gp*(z1 - zp) - xi*(z - z1)
    arg008 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - chi*(z - z1)
    arg009 = aaa/2.d0/tmp**2*xi + xi*(z1 - zp) - xi*(z - z1)
    arg010 = aaa/2.d0/tmp**2*xi + chi*(z1 - zp) - xi*(z - z1)
    arg011 = aaa/2.d0/tmp**2*chi + chi*(z1 - zp) - xi*(z - z1)
    arg012 = aaa/2.d0/tmp**2*chi + xi*(z1 - zp) - chi*(z - z1)
    arg104 =  gp/2.d0/tmp - tmp*(z1 - zp)
    arg105 =  gp/2.d0/tmp + tmp*(z - z1)
    arg107 =  xi/2.d0/tmp + tmp*(z - zp)
    arg108 =  xi/2.d0/tmp - tmp*(z - zp)
    arg109 =  xi/2.d0/tmp + tmp*(z1 - zp)
    arg110 =  xi/2.d0/tmp - tmp*(z - z1)
    arg111 = chi/2.d0/tmp + tmp*(z - zp)
    arg112 = chi/2.d0/tmp - tmp*(z - zp)
    arg113 = chi/2.d0/tmp + tmp*(z1 - zp)
    arg114 = chi/2.d0/tmp - tmp*(z - z1)
    t1 = exp(arg000)*beta*xi/alpha/lambda
    t2 = (exp_erfc(arg011, arg111) - exp_erfc(arg011, arg113))*chi/lambda &
         - exp_erfc(arg008, arg107)*xi/lambda &
         + exp_erfc(arg009, arg109)*xi*beta/alpha/lambda &
         + exp_erfc(arg005, arg104)*gp/alpha*2.d0
    t3 = (exp_erfc(arg012, arg112) - exp_erfc(arg012, arg114))*xi/lambda &
         - exp_erfc(arg010, arg108)*chi/lambda &
         + exp_erfc(arg009, arg110)*xi*beta/alpha/lambda &
         - exp_erfc(arg004, arg105)*xi/alpha*2.d0
    dvl22j0 = gp*(t1 - (t2 + t3)/4.d0)*dbesj0(gp*rxy)
    !
    RETURN
  END FUNCTION dvl22j0

!
! Bessel J_0(x) function in double precision
!
  REAL(8) FUNCTION dbesj0(x)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x
    REAL(8), PARAMETER  :: pi4 = 0.78539816339744830962d0
    REAL(8), PARAMETER  :: a(0:7) = (/ &
       -0.0000000000023655394d0, 0.0000000004708898680d0, &
       -0.0000000678167892231d0, 0.0000067816840038636d0, &
       -0.0004340277777716935d0, 0.0156249999999992397d0, &
       -0.2499999999999999638d0, 0.9999999999999999997d0/)
    REAL(8), PARAMETER  :: b(0:12, 0:4) = reshape((/ &
        0.0000000000626681117d0, -0.0000000022270614428d0, &
        0.0000000662981656302d0, -0.0000016268486502196d0, &
        0.0000321978384111685d0, -0.0005005237733315830d0, &
        0.0059060313537449816d0, -0.0505265323740109701d0, &
        0.2936432097610503985d0, -1.0482565081091638637d0, &
        1.9181123286040428113d0, -1.1319199475221700100d0, &
       -0.1965480952704682000d0, &
        0.0000000000457457332d0, -0.0000000015814772025d0, &
        0.0000000455487446311d0, -0.0000010735201286233d0, &
        0.0000202015179970014d0, -0.0002942392368203808d0, &
        0.0031801987726150648d0, -0.0239875209742846362d0, &
        0.1141447698973777641d0, -0.2766726722823530233d0, &
        0.1088620480970941648d0,  0.5136514645381999197d0, &
       -0.2100594022073706033d0, &
        0.0000000000331366618d0, -0.0000000011119090229d0, &
        0.0000000308823040363d0, -0.0000006956602653104d0, &
        0.0000123499947481762d0, -0.0001662951945396180d0, &
        0.0016048663165678412d0, -0.0100785479932760966d0, &
        0.0328996815223415274d0, -0.0056168761733860688d0, &
       -0.2341096400274429386d0,  0.2551729256776404262d0, &
        0.2288438186148935667d0, &
        0.0000000000238007203d0, -0.0000000007731046439d0, &
        0.0000000206237001152d0, -0.0000004412291442285d0, &
        0.0000073107766249655d0, -0.0000891749801028666d0, &
        0.0007341654513841350d0, -0.0033303085445352071d0, &
        0.0015425853045205717d0,  0.0521100583113136379d0, &
       -0.1334447768979217815d0, -0.1401330292364750968d0, &
        0.2685616168804818919d0, &
        0.0000000000169355950d0, -0.0000000005308092192d0, &
        0.0000000135323005576d0, -0.0000002726650587978d0, &
        0.0000041513240141760d0, -0.0000443353052220157d0, &
        0.0002815740758993879d0, -0.0004393235121629007d0, &
       -0.0067573531105799347d0,  0.0369141914660130814d0, &
        0.0081673361942996237d0, -0.2573381285898881860d0, &
        0.0459580257102978932d0/), (/13, 5/))
    REAL(8), PARAMETER  :: c(0:13, 0:4) = reshape((/ &
       -0.00000000003009451757d0, -0.00000000014958003844d0, &
        0.00000000506854544776d0,  0.00000001863564222012d0, &
       -0.00000060304249068078d0, -0.00000147686259937403d0, &
        0.00004714331342682714d0,  0.00006286305481740818d0, &
       -0.00214137170594124344d0, -0.00089157336676889788d0, &
        0.04508258728666024989d0, -0.00490362805828762224d0, &
       -0.27312196367405374426d0,  0.04193925184293450356d0,  &
       -0.00000000000712453560d0, -0.00000000041170814825d0, &
        0.00000000138012624364d0,  0.00000005704447670683d0, &
       -0.00000019026363528842d0, -0.00000533925032409729d0, &
        0.00001736064885538091d0,  0.00030692619152608375d0, &
       -0.00092598938200644367d0, -0.00917934265960017663d0, &
        0.02287952522866389076d0,  0.10545197546252853195d0, &
       -0.16126443075752985095d0, -0.19392874768742235538d0,  &
        0.00000000002128344556d0, -0.00000000031053910272d0, &
       -0.00000000334979293158d0,  0.00000004507232895050d0, &
        0.00000036437959146427d0, -0.00000446421436266678d0, &
       -0.00002523429344576552d0,  0.00027519882931758163d0, &
        0.00097185076358599358d0, -0.00898326746345390692d0, &
       -0.01665959196063987584d0,  0.11456933464891967814d0, &
        0.07885001422733148815d0, -0.23664819446234712621d0,  &
        0.00000000003035295055d0,  0.00000000005486066835d0, &
       -0.00000000501026824811d0, -0.00000000501246847860d0, &
        0.00000058012340163034d0,  0.00000016788922416169d0, &
       -0.00004373270270147275d0,  0.00001183898532719802d0, &
        0.00189863342862291449d0, -0.00113759249561636130d0, &
       -0.03846797195329871681d0,  0.02389746880951420335d0, &
        0.22837862066532347461d0, -0.06765394811166522844d0,  &
        0.00000000001279875977d0,  0.00000000035925958103d0, &
       -0.00000000228037105967d0, -0.00000004852770517176d0, &
        0.00000028696428000189d0,  0.00000440131125178642d0, &
       -0.00002366617753349105d0, -0.00024412456252884129d0, &
        0.00113028178539430542d0,  0.00708470513919789080d0, &
       -0.02526914792327618386d0, -0.08006137953480093426d0, &
        0.16548380461475971846d0,  0.14688405470042110229d0/), (/14, 5/))
    REAL(8), PARAMETER  :: d(0:12, 0:3) = reshape((/ &
        1.059601355592185731d-14, -2.71150591218550377d-13, &
        8.6514809056201638d-12,   -4.6264028554286627d-10, &
        5.0815403835647104d-8,    -1.76722552048141208d-5, &
        0.16286750396763997378d0,  2.949651820598278873d-13, &
       -8.818215611676125741d-12,  3.571119876162253451d-10, &
       -2.631924120993717060d-8,   4.709502795656698909d-6, &
       -5.208333333333283282d-3, &
        7.18344107717531977d-15,  -2.51623725588410308d-13, &
        8.6017784918920604d-12,   -4.6256876614290359d-10, &
        5.0815343220437937d-8,    -1.76722551764941970d-5, &
        0.16286750396763433767d0,  2.2327570859680094777d-13, &
       -8.464594853517051292d-12,  3.563766464349055183d-10, &
       -2.631843986737892965d-8,   4.709502342288659410d-6, &
       -5.2083333332278466225d-3, &
        5.15413392842889366d-15,  -2.27740238380640162d-13, &
        8.4827767197609014d-12,   -4.6224753682737618d-10, &
        5.0814848128929134d-8,    -1.76722547638767480d-5, &
        0.16286750396748926663d0,  1.7316195320192170887d-13, &
       -7.971122772293919646d-12,  3.544039469911895749d-10, &
       -2.631443902081701081d-8,   4.709498228695400603d-6, &
       -5.2083333315143653610d-3, &
        3.84653681453798517d-15,  -2.04464520778789011d-13, &
        8.3089298605177838d-12,   -4.6155016158412096d-10, &
        5.0813263696466650d-8,    -1.76722528311426167d-5, &
        0.16286750396650065930d0,  1.3797879972460878797d-13, &
       -7.448089381011684812d-12,  3.512733797106959780d-10, &
       -2.630500895563592722d-8,   4.709483934775839193d-6, &
       -5.2083333227940760113d-3/), (/13, 4/))
    REAL(8) :: w, t, y, v, theta
    INTEGER :: k, i

    w = ABS(x)
    IF (w < 1.0d0) THEN
      t = w*w
      y = a(0)
      DO i = 1, 7
        y = y*t + a(i)
      END DO
    ELSE IF (w < 8.5d0) THEN
      t = w*w*0.0625d0
      k = int(t)
      t = t - (k + 0.5d0)
      y = b(0, k)
      DO i = 1, 12
        y = y*t + b(i, k)
      END DO
    ELSE IF (w < 12.5d0) THEN
      k = int(w)
      t = w - (k + 0.5d0)
      k = k - 8
      y = c(0, k)
      DO i = 1, 13
        y = y*t + c(i, k)
      END DO
    ELSE
      v = 24.0d0/w
      t = v*v
      k = int(t)
      y = d(0, k)
      DO i = 1, 6
        y = y*t + d(i, k)
      END DO
      y = y*sqrt(v)
      theta = d(7, k)
      DO i = 8, 12
        theta = theta*t + d(i, k)
      END DO
      theta = theta*v - pi4
      y = y*cos(w + theta)
    END IF
    dbesj0 = y
  END FUNCTION dbesj0

  REAL(8) FUNCTION dbesj1(x)
    IMPLICIT NONE
    REAL(8), INTENT(in) :: x
    REAL(8), PARAMETER :: pi4 = 0.78539816339744830962d0
    REAL(8), PARAMETER :: a(0:7) = (/ &
       -0.00000000000014810349d0,  0.00000000003363594618d0, &
       -0.00000000565140051697d0,  0.00000067816840144764d0, &
       -0.00005425347222188379d0,  0.00260416666666662438d0, &
       -0.06249999999999999799d0,  0.49999999999999999998d0/)
    REAL(8), PARAMETER :: b(0:12, 0:4) = reshape((/ &
        0.00000000000243721316d0, -0.00000000009400554763d0, &
        0.00000000306053389980d0, -0.00000008287270492518d0, &
        0.00000183020515991344d0, -0.00003219783841164382d0, &
        0.00043795830161515318d0, -0.00442952351530868999d0, &
        0.03157908273375945955d0, -0.14682160488052520107d0, &
        0.39309619054093640008d0, -0.47952808215101070280d0, &
        0.14148999344027125140d0, &
        0.00000000000182119257d0, -0.00000000006862117678d0, &
        0.00000000217327908360d0, -0.00000005693592917820d0, &
        0.00000120771046483277d0, -0.00002020151799736374d0, &
        0.00025745933218048448d0, -0.00238514907946126334d0, &
        0.01499220060892984289d0, -0.05707238494868888345d0, &
        0.10375225210588234727d0, -0.02721551202427354117d0, &
       -0.06420643306727498985d0, &
        0.000000000001352611196d0, -0.000000000049706947875d0, &
        0.000000001527944986332d0, -0.000000038602878823401d0, &
        0.000000782618036237845d0, -0.000012349994748451100d0, &
        0.000145508295194426686d0, -0.001203649737425854162d0, &
        0.006299092495799005109d0, -0.016449840761170764763d0, &
        0.002106328565019748701d0,  0.058527410006860734650d0, &
       -0.031896615709705053191d0, &
        0.000000000000997982124d0, -0.000000000035702556073d0, &
        0.000000001062332772617d0, -0.000000025779624221725d0, &
        0.000000496382962683556d0, -0.000007310776625173004d0, &
        0.000078028107569541842d0, -0.000550624088538081113d0, &
        0.002081442840335570371d0, -0.000771292652260286633d0, &
       -0.019541271866742634199d0,  0.033361194224480445382d0, &
        0.017516628654559387164d0, &
        0.000000000000731050661d0, -0.000000000025404499912d0, &
        0.000000000729360079088d0, -0.000000016915375004937d0, &
        0.000000306748319652546d0, -0.000004151324014331739d0, &
        0.000038793392054271497d0, -0.000211180556924525773d0, &
        0.000274577195102593786d0,  0.003378676555289966782d0, &
       -0.013842821799754920148d0, -0.002041834048574905921d0, &
        0.032167266073736023299d0/), (/13, 5/))
    REAL(8), PARAMETER :: c(0:13, 0:4) = reshape((/ &
        -0.00000000001185964494d0,  0.00000000039110295657d0, &
         0.00000000180385519493d0, -0.00000005575391345723d0, &
        -0.00000018635897017174d0,  0.00000542738239401869d0, &
         0.00001181490114244279d0, -0.00033000319398521070d0, &
        -0.00037717832892725053d0,  0.01070685852970608288d0, &
         0.00356629346707622489d0, -0.13524776185998074716d0, &
         0.00980725611657523952d0,  0.27312196367405374425d0,  &
        -0.00000000003029591097d0,  0.00000000009259293559d0, &
         0.00000000496321971223d0, -0.00000001518137078639d0, &
        -0.00000057045127595547d0,  0.00000171237271302072d0, &
         0.00004271400348035384d0, -0.00012152454198713258d0, &
        -0.00184155714921474963d0,  0.00462994691003219055d0, &
         0.03671737063840232452d0, -0.06863857568599167175d0, &
        -0.21090395092505707655d0,  0.16126443075752985095d0,  &
        -0.00000000002197602080d0, -0.00000000027659100729d0, &
         0.00000000374295124827d0,  0.00000003684765777023d0, &
        -0.00000045072801091574d0, -0.00000327941630669276d0, &
         0.00003571371554516300d0,  0.00017664005411843533d0, &
        -0.00165119297594774104d0, -0.00485925381792986774d0, &
         0.03593306985381680131d0,  0.04997877588191962563d0, &
        -0.22913866929783936544d0, -0.07885001422733148814d0,  &
         0.00000000000516292316d0, -0.00000000039445956763d0, &
        -0.00000000066220021263d0,  0.00000005511286218639d0, &
         0.00000005012579400780d0, -0.00000522111059203425d0, &
        -0.00000134311394455105d0,  0.00030612891890766805d0, &
        -0.00007103391195326182d0, -0.00949316714311443491d0, &
         0.00455036998246516948d0,  0.11540391585989614784d0, &
        -0.04779493761902840455d0, -0.22837862066532347460d0,  &
         0.00000000002697817493d0, -0.00000000016633326949d0, &
        -0.00000000433134860350d0,  0.00000002508404686362d0, &
         0.00000048528284780984d0, -0.00000258267851112118d0, &
        -0.00003521049080466759d0,  0.00016566324273339952d0, &
         0.00146474737522491617d0, -0.00565140892697147306d0, &
        -0.02833882055679300400d0,  0.07580744376982855057d0, &
         0.16012275906960187978d0, -0.16548380461475971845d0/), (/14, 5/))
    REAL(8), PARAMETER :: d(0:12, 0:3) = reshape((/ &
        -1.272346002224188092d-14,  3.370464692346669075d-13, &
        -1.144940314335484869d-11,  6.863141561083429745d-10, &
        -9.491933932960924159d-8,   5.301676561445687562d-5,  &
         0.1628675039676399740d0,  -3.652982212914147794d-13, &
         1.151126750560028914d-11, -5.165585095674343486d-10, &
         4.657991250060549892d-8,  -1.186794704692706504d-5,  &
         1.562499999999994026d-2, &
        -8.713069680903981555d-15,  3.140780373478474935d-13, &
        -1.139089186076256597d-11,  6.862299023338785566d-10, &
        -9.491926788274594674d-8,   5.301676558106268323d-5,  &
         0.1628675039676466220d0,  -2.792555727162752006d-13, &
         1.108650207651756807d-11, -5.156745588549830981d-10, &
         4.657894859077370979d-8,  -1.186794650130550256d-5,  &
         1.562499999987299901d-2, &
        -6.304859171204770696d-15,  2.857249044208791652d-13, &
        -1.124956921556753188d-11,  6.858482894906716661d-10, &
        -9.491867953516898460d-8,   5.301676509057781574d-5,  &
         0.1628675039678191167d0,  -2.185193490132496053d-13, &
         1.048820673697426074d-11, -5.132819367467680132d-10, &
         4.657409437372994220d-8,  -1.186794150862988921d-5,  &
         1.562499999779270706d-2, &
        -4.740417209792009850d-15,  2.578715253644144182d-13, &
        -1.104148898414138857d-11,  6.850134201626289183d-10, &
        -9.491678234174919640d-8,   5.301676277588728159d-5,  &
         0.1628675039690033136d0,  -1.755122057493842290d-13, &
         9.848723331445182397d-12, -5.094535425482245697d-10, &
         4.656255982268609304d-8,  -1.186792402114394891d-5,  &
         1.562499998712198636d-2/), (/13, 4/))
    REAL(8) :: w, t, y, v, theta
    INTEGER :: k, i

    w = ABS(x)
    IF (w < 1.0d0) THEN
      t = w*w
      y = a(0)
      DO i = 1, 7
        y = y*t + a(i)
      END DO
      y = y*w
    ELSE IF (w < 8.5d0) THEN
      t = w*w*0.0625d0
      k = int(t)
      t = t - (k + 0.5d0)
      y = b(0, k)
      DO i = 1, 12
        y = y*t + b(i, k)
      END DO
      y = y*w
    ELSE IF (w < 12.5d0) THEN
      k = int(w)
      t = w - (k + 0.5d0)
      k = k - 8
      y = c(0, k)
      DO i = 1, 13
        y = y*t + c(i, k)
      END DO
    ELSE
      v = 24.0d0/w
      t = v*v
      k = int(t)
      y = d(0, k)
      DO i = 1, 6
        y = y*t + d(i, k)
      END DO
      y = y*sqrt(v)
      theta = d(7, k)
      DO i = 8, 12
        theta = theta*t + d(i, k)
      END DO
      theta = theta*v - pi4
      y = y*sin(w + theta)
    END IF
    IF (x < 0.0d0) y = -y
    dbesj1 = y
  END FUNCTION dbesj1

  SUBROUTINE polint(xa, ya, n, x, y, dy)
    USE kinds, ONLY : DP
    INTEGER, INTENT(in)  ::n
    INTEGER, PARAMETER   ::nmax = 10
    INTEGER              ::i, m, ns
    REAL(DP), INTENT(in) ::x, xa(n), ya(n)
    REAL(DP), INTENT(out)::y, dy
    REAL(DP)             ::den, dif, dift, ho, hp, w, c(nmax), d(nmax)

    ns = 1
    dif = ABS(x - xa(1))
    DO i = 1, n
      dift = ABS(x - xa(i))
      IF (dift .lt. dif) THEN
        ns = i
        dif = dift
      ENDIF
      c(i) = ya(i)
      d(i) = ya(i)
    ENDDO
    y = ya(ns)
    ns = ns - 1
    DO m = 1, n - 1
      DO i = 1, n - m
        ho = xa(i) - x
        hp = xa(i + m) - x
        w = c(i + 1) - d(i)
        den = ho - hp
        IF (den .eq. 0.) STOP 'failure in polint'
        den = w/den
        d(i) = hp*den
        c(i) = ho*den
      ENDDO
      IF (2*ns .lt. n - m) THEN
        dy = c(ns + 1)
      ELSE
        dy = d(ns)
        ns = ns - 1
      ENDIF
      y = y + dy
    ENDDO
    RETURN
  END SUBROUTINE polint
  !
  ! Checks inversion symmetry along z-axis
  !
  LOGICAL FUNCTION esm_z_inv(lrism)
    !
    USE constants, ONLY : eps14
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lrism
    !
    esm_z_inv = .TRUE.
    !
    IF (do_comp_esm) THEN
      IF (TRIM(esm_bc) == 'bc1') THEN
        esm_z_inv = (.NOT. lrism)
      ELSE IF (TRIM(esm_bc) == 'bc2') THEN
        esm_z_inv = (ABS(esm_efield) < eps14)
      ELSE IF (TRIM(esm_bc) == 'bc3') THEN
        esm_z_inv = .FALSE.
      ELSE IF (TRIM(esm_bc) == 'bc4') THEN
        esm_z_inv = .FALSE.
      END IF
    END IF
    !
  END FUNCTION esm_z_inv

END MODULE esm_common_mod
