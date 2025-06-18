  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE bloch2wannier
  !----------------------------------------------------------------------
  !!
  !! Contains all the routines that transforms quantities from Bloch space
  !! to Wannier space.
  !!
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------------
    SUBROUTINE hambloch2wan(nbnd, nbndsub, nks, nkstot, et, xk, cu, &
       lwin, exband, nrr, irvec, wslen, chw)
    !--------------------------------------------------------------------------
    !!
    !!  From the Hamiltonian in Bloch representationi (coarse mesh),
    !!  find the corresponding Hamiltonian in Wannier representation
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, alat
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero, zero, ryd2ev
    USE io_global, ONLY : ionode_id
    USE io_var,    ONLY : iudecayH
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE global_var,ONLY : nbndep
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lwin(nbndep, nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands
    INTEGER, INTENT(in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nks
    !! number of kpoints in this pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of kpoints
    INTEGER, INTENT(in) ::  nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! coordinates of WS points
    !
    REAL(KIND = DP), INTENT(in) :: et(nbnd, nks)
    !! hamiltonian eigenvalues, coarse mesh
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbndep, nbndsub, nks)
    !! rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(out) :: chw(nbndsub, nbndsub, nrr)
    !! Hamiltonian in Wannier basis
    !
    ! Local variables
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: ibnd1
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: mbnd
    !! Counter on band index
    INTEGER :: i
    !! Counter on band index
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}
    REAL(KIND = DP) :: tmp
    !! Maximum value of the real space Hamiltonian
    REAL(KIND = DP) :: et_opt(nbndep, nks)
    !! hamiltonian eigenvalues within the outer window in the first ndimwin(ik) entries
    !
    COMPLEX(KIND = DP) :: cfac
    !! $$e^{-i\mathbf{r}\cdot\mathbf{k}}$$
    COMPLEX(KIND = DP) :: ctmp
    !! Temporary variable to store the Hamiltonian
    COMPLEX(KIND = DP), ALLOCATABLE :: chs(:, :, :)
    !! Hamiltonian in Bloch basis, coarse k-mesh
    !
    !----------------------------------------------------------
    !    STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    CALL start_clock ('Ham: step 1')
    !
    ALLOCATE(chs(nbndsub, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('hambloch2wan', 'Error allocating chs', 1)
    !
    ! slim down et to contain states within the outer window
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF
    ENDDO
    !
    et_opt = zero
    IF (nexband_tmp > 0) THEN
      DO ik = 1, nks
        ibnd = 0
        ibnd1 = 0
        DO i = 1, nbnd
          IF (exband(i)) CYCLE
          ibnd1 = ibnd1 + 1
          IF (lwin(ibnd1, ik)) THEN
            ibnd = ibnd + 1
            et_opt(ibnd, ik) = et(i, ik)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbndep
          IF (lwin(i, ik)) THEN
            ibnd = ibnd + 1
            et_opt(ibnd, ik) = et(i, ik)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    ! [Eqn. 26 of PRB 76, 165108 (2007)]
    ! H~(k) = U(k)^\dagger * H(k) * U(k)
    ! H~(k) is chs( nbndsub, nbndsub, ik )
    ! Note H~(k) is H^(W)(k) in PRB 74, 195118 (2006) notations
    !
    chs(:, :, :) = czero
    DO ik = 1, nks
      !
      DO jbnd = 1, nbndsub
        DO ibnd = 1, jbnd
          !
          ctmp = czero
          DO mbnd = 1, nbndep
            ctmp = ctmp + CONJG(cu(mbnd, ibnd, ik)) * et_opt(mbnd, ik) * cu(mbnd, jbnd, ik)
          ENDDO
          !
          chs(ibnd, jbnd, ik) = ctmp
          chs(jbnd, ibnd, ik) = CONJG(ctmp)
          !
        ENDDO
      ENDDO
      !
    ENDDO
    !
    CALL stop_clock('Ham: step 1')
    !
    !----------------------------------------------------------
    !  STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 26 of PRB 76, 165108 (2007)]
    ! H(R) = (1/nk) sum_k e^{-ikR} H~(k)
    ! H(R) is chw(nbndsub, nbndsub, nrr)
    ! Note H(R) is H^(W)(R) in PRB 74, 195118 (2006) notations
    !
    CALL start_clock ('Ham: step 2')
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nks, xk, at, -1)
    !
    chw( :, :, :) = czero
    !
    DO ir = 1, nrr
      DO ik = 1, nks
         !
         rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
         cfac = EXP(-ci * rdotk ) / DBLE(nkstot)
         chw(:, :, ir) = chw(:, :, ir ) + cfac * chs(:, :, ik)
         !
      ENDDO
    ENDDO
    CALL mp_sum(chw,inter_pool_comm)
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nks, xk, bg, 1)
    !
    ! check spatial decay of Hamiltonian in Wannier basis
    !
    ! Hamiltonian matrix is in Ryd units and spatial dimensions in bohr
    ! [mind when comparing with wannier code (eV and angstrom units) with
    ! write_hr=.TRUE.]
    !
    IF (mpime == ionode_id) THEN
      !
      OPEN(UNIT = iudecayH,FILE = 'decay.H')
      WRITE(iudecayH, '(a)') '# Spatial decay of Hamiltonian in Wannier basis'
      WRITE(iudecayH, '(a)') '# R_e [Ang]      max_{n,m} |H(n,m)| [Ry] '
      DO ir = 1, nrr
        !
        tmp = MAXVAL(ABS(chw(:, :, ir)))
        WRITE(iudecayH,*) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      ! RMDB
      !DO ir = 1, nrr
      !  DO jbnd = 1, nbndsub
      !    DO ibnd = 1, nbndsub
      !      WRITE(iudecayH, '(5I5,6F12.6)') irvec(:, ir), ibnd, jbnd, chw(ibnd, jbnd, ir) * ryd2ev
      !    ENDDO
      !  ENDDO
      !ENDDO
      !
      CLOSE(iudecayH)
      !
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(chs, STAT = ierr)
    IF (ierr /= 0) CALL errore('hambloch2wan', 'Error deallocating chs', 1)
    !
    CALL stop_clock('Ham: step 2')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE hambloch2wan
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dmebloch2wan(nbnd, nbndsub, nks, nkstot, dmec, xk, cu, &
                            nrr, irvec, wslen, lwin, exband, cdmew)
    !--------------------------------------------------------------------------
    !!
    !!  From the Dipole in Bloch representationi (coarse mesh),
    !!  find the corresponding Dipole in Wannier representation
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE global_var,    ONLY : nbndep
    USE io_var,        ONLY : iudecayP
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero, cone
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! number of bands
    INTEGER, INTENT(in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nks
    !! number of kpoints in this pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of kpoints
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! Coordinate of Wannier space points
    !
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(KIND = DP), INTENT(in) :: dmec(3, nbndep, nbndep, nks)
    !! Dipole matrix elements on coarse mesh
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbndep, nbndsub, nks)
    !! rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(inout) :: cdmew(3, nbndsub, nbndsub, nrr)
    !! Dipole matrix elements in Wannier basis
    !
    LOGICAL, INTENT(in) :: lwin(nbndep, nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    !
    ! Local variables
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: i
    !! Counter on total band index
    INTEGER :: j
    !! Counter on total band index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: ibnd1
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: jbnd1
    !! Counter on band index
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: dmec_utmp(nbndep, nbndsub)
    !! dmec after multiplication with the Wannier rotation matrix cu.
    COMPLEX(KIND = DP), ALLOCATABLE :: cps(:, :, :, :)
    !! Dipole in smooth Bloch basis, coarse mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: dmec_opt(:, :, :, :)
    !! dmec computed in pmn, rescaled down if skipping lwin.
    !
    CALL start_clock ('Dipole: step 1')
    !
    ALLOCATE(cps(3, nbndsub, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('dmebloch2wan', 'Error allocating cps', 1)
    ALLOCATE(dmec_opt(3, nbndep, nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('dmebloch2wan', 'Error allocating dmec_opt', 1)
    !
    !--------------------------------------------------------------
    !    STEP 0: Rescale the optical matrix on the coarse grid down
    !            This is if you skip band during the Wannierization
    !--------------------------------------------------------------
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF
    ENDDO
    !
    dmec_opt = czero
    IF (nexband_tmp > 0) THEN
      DO ik = 1,nks
        jbnd = 0
        jbnd1 = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          jbnd1 = jbnd1 + 1
          IF (lwin(jbnd1,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            ibnd1 = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              ibnd1 = ibnd1 + 1
              IF (lwin(ibnd1,ik)) THEN
                ibnd = ibnd + 1
                dmec_opt(:, ibnd, jbnd, ik) = dmec(:, ibnd1, jbnd1, ik)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbndep
          IF (lwin(j, ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbndep
              IF (lwin(i, ik)) THEN
                ibnd = ibnd + 1
                dmec_opt(:, ibnd, jbnd, ik) = dmec(:, i, j, ik)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    !----------------------------------------------------------
    !    STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    ! p~(k) = U(k)^\dagger * p(k) * U(k)
    ! p~(k) is cps( ipol, nbndsub, nbndsub, ik )
    ! Note p~(k) is p^(W)(k) in PRB 74, 195118 (2006) notations
    !
    dmec_utmp(:, :) = czero
    cps(:, :, :, :) = czero
    !
    DO ik = 1, nks
      DO ipol = 1, 3
        !
        ! dmec_utmp(:, :) = MATMUL( dmec_opt(ipol,:,:,ik), cu(:,:,ik) )
        ! cps(ipol,:,:,ik) = MATMUL( CONJG(transpose( cu(:,:,ik))), dmec_utmp(:, :) )
        !
        CALL ZGEMM('n', 'n', nbndep, nbndsub, nbndep, cone, dmec_opt(ipol, :, :, ik), &
                   nbndep, cu(:, :, ik), nbndep, czero, dmec_utmp(:, :), nbndep)
        CALL ZGEMM('c', 'n', nbndsub, nbndsub, nbndep, cone, cu(:, :, ik), &
                   nbndep, dmec_utmp(:, :), nbndep, czero, cps(ipol, :, :, ik), nbndsub)
        !
      ENDDO
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! p(R) = (1/nk) sum_k e^{-ikR} p~(k)
    ! p(R) is cdmew(ipol, nbndsub, nbndsub, nrr)
    ! Note p(R) is p^(W)(R) in PRB 74, 195118 (2006) notations
    !
    CALL start_clock('Dipole: step 2')
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nks, xk, at, -1)
    !
    cdmew( :, :, :, :) = czero
    !
    DO ir = 1, nrr
      DO ik = 1, nks
         !
         rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
         cfac = EXP(-ci * rdotk) / DBLE(nkstot)
         cdmew(:, :, :, ir) = cdmew(:, :, :, ir) + cfac * cps(:, :, :, ik)
         !
      ENDDO
    ENDDO
    CALL mp_sum(cdmew, inter_pool_comm)
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nks, xk, bg, 1)
    !
    ! Check spatial decay of Dipole in Wannier basis
    ! the unit in r-space is angstrom
    !
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iudecayP, FILE = 'decay.P')
      WRITE(iudecayP, '(a)') '# Spatial decay of dipole in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  MAXVAL(ABS(cdmew(:, :,:,ir)))
        WRITE(iudecayP,*) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iudecayP)
    ENDIF
    !CALL mp_barrier(inter_pool_comm)
    !
    DEALLOCATE(cps, STAT = ierr)
    IF (ierr /= 0) CALL errore('dmebloch2wan', 'Error deallocating cps', 1)
    DEALLOCATE(dmec_opt, STAT = ierr)
    IF (ierr /= 0) CALL errore('dmebloch2wan', 'Error deallocating dmec_opt', 1)
    !
    CALL stop_clock('Dipole: step 2')
    !
    !------------------------------------------------------------------------
    END SUBROUTINE dmebloch2wan
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    SUBROUTINE dynbloch2wan(nmodes, nq, xk, dynq, nrr, irvec, wslen )
    !------------------------------------------------------------------------
    !!
    !!  From the Dynamical Matrix in Bloch representation (coarse mesh),
    !!  find the corresponding matrix in Wannier representation
    !!
    !!  NOTA BENE: it seems to be very important that the matrix is kept real
    !!  physically these are truely the interatomic force constants.
    !!  If you use a complex matrix instead, you may get some spurious
    !!  oscillations when you interpolate the phonon dispersions.
    !!
    !!  Note also that the acoustic sum rule for the q=0 case has been imposed
    !!  already in readmat_shuffle.f90
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE ions_base,     ONLY : nat, tau
    USE global_var,    ONLY : rdw, epsi, zstar, qrpl, L
    USE input,         ONLY : lpolar, nqc1, nqc2, nqc3, system_2d
    USE io_var,        ONLY : iudecaydyn
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero, fpi, eps6
    USE io_global,     ONLY : ionode_id, stdout
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier
    USE mp_global,     ONLY : inter_pool_comm
    USE longrange,     ONLY : rgd_blk
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq
    !! number of qpoints
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! coordinates of WS points
    REAL(KIND = DP), INTENT(in) :: xk(3, nq)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    COMPLEX(KIND = DP), INTENT(inout) :: dynq(nmodes, nmodes, nq)
    !! dynamical matrix in bloch representation (Cartesian coordinates)
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER, PARAMETER :: maxiter = 30
    !! Maximum interation
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    REAL(KIND = DP) :: Lmin, Lmax, L1, L2
    !! Min and max value of L
    REAL(KIND = DP) :: ifc1, ifc2, ifc
    !! sum of IFC for a given L
    REAL(KIND = DP) :: c
    !! vacuum size (supercell length along the z direction) in case of 2D
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    !
    ! In the case of 2D, one need to compute the optimal range separation length L [Eq. 63 of PRB 107, 155424 (2023)]
    IF (system_2d == 'dipole_sp' .OR. system_2d == 'quadrupole') THEN
      WRITE(stdout, '(5x, a)') 'Find optimal range separation length L'
      WRITE(stdout, '(5x, a)') ' '
      ! We find the L that minimize the sum of the real-space IFC using trisection [bisection useless for parabolas].
      Lmin = 1
      Lmax = 80 ! Bohr
      !
      DO WHILE (Lmax - Lmin > eps6)
        L1 = Lmin + (Lmax - Lmin) / 3.0d0
        L2 = Lmin + (Lmax - Lmin) *2.0d0 / 3.0d0
        CALL findL(L1, nmodes, nq, xk, dynq, nrr, irvec, ifc1)
        CALL findL(L2, nmodes, nq, xk, dynq, nrr, irvec, ifc2)
        ! Trisection
        IF (ifc1 < ifc2) THEN
          Lmax = L2
        ELSE
          Lmin = L1
        ENDIF
        L = (Lmin + Lmax) / 2.0_dp
        CALL findL(L, nmodes, nq, xk, dynq, nrr, irvec, ifc)
        WRITE(stdout, '(5x, a, f12.5, a, f20.12)') 'L ', L, ' Bohr with IFC = ', ifc
      ENDDO
      c = alat / bg(3, 3)
      WRITE(stdout, *) '     '
      WRITE(stdout, *) '     In-plane polarizability in cartesian axis'
      WRITE(stdout, *) '     ', (epsi(1, 1) - 1.0) * c / fpi, '  ',epsi(1, 2) * c / fpi
      WRITE(stdout, *) '     ', epsi(2, 1) * c / fpi, '  ',(epsi(2, 2) - 1.0) * c / fpi
      WRITE(stdout, *) '     Out-of-plane polarizability in cartesian axis ',&
                                (epsi(3, 3) - 1.0) * c / fpi
      WRITE(stdout, *) ' '
    ENDIF
    !
    ! Subtract the long-range term from D(q)
    IF (lpolar .OR. qrpl) THEN
      DO ik = 1, nq
        !xk has to be in cart. coord.
        CALL rgd_blk(L, nqc1, nqc2, nqc3, nat, dynq(:, :, ik), xk(:, ik), tau, epsi, zstar, -1.d0)
        !
      ENDDO
    ENDIF
    !
    !----------------------------------------------------------
    !  Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 29 of PRB 76, 165108 (2007)]
    ! D(R) = (1/nk) sum_k e^{-ikR} D~(k)
    ! D(R) is rdw(nmodes, nmodes, ir)
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nq, xk, at, -1)
    !
    rdw( :, :, :) = czero
    !
    DO ir = 1, nrr
      !
      DO ik = 1, nq
        !
        rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
        cfac = EXP(-ci * rdotk)/ DBLE(nq)
        !DBSP - real was commented
        rdw(:, :, ir) = rdw(:, :, ir) +  cfac * dynq(:, :, ik)
        !rdw( :, :, ir ) = rdw( :, :, ir ) + real ( cfac * dynq( :, :, ik ) )
        !                                    ^^^^
        !                                 note this
      ENDDO
      !
    ENDDO
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nq, xk, bg, 1)
    !
    ! check spatial decay of dynamical matrix in Wannier basis
    ! the unit in r-space is angstrom, and I am plotting the matrix for the first mode only
    !
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iudecaydyn, FILE = 'decay.dynmat')
      WRITE(iudecaydyn, '(a)') '# Spatial decay of Dynamical matrix in Wannier basis'
      WRITE(iudecaydyn, '(a)') '# R_p [Ang]      max_{mu,nu} |D(mu,nu)| [Ry] '
      DO ir = 1, nrr
        !
        tmp =  MAXVAL(ABS(rdw(:, :, ir)))
        WRITE(iudecaydyn, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      CLOSE(iudecaydyn)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE dynbloch2wan
    !--------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    SUBROUTINE findL(L, nmodes, nq, xk, dynq, nrr, irvec, ifc)
    !------------------------------------------------------------------------
    !!
    !! Compute the sum of the IFC in real space for a given separation length L
    !! See Eq. 63 of PRB 107, 155424 (2023).
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE ions_base,     ONLY : nat, tau
    USE global_var,    ONLY : epsi, zstar
    USE input,         ONLY : nqc1, nqc2, nqc3
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero
    USE mp,            ONLY : mp_barrier
    USE longrange,     ONLY : rgd_blk
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq
    !! number of qpoints
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! coordinates of WS points
    REAL(KIND = DP), INTENT(in) :: L
    !! Range separation length in Bohr
    REAL(KIND = DP), INTENT(in) :: xk(3, nq)
    !! kpoint coordinates (cartesian in units of 2piba)
    COMPLEX(KIND = DP), INTENT(in) :: dynq(nmodes, nmodes, nq)
    !! dynamical matrix in bloch representation (Cartesian coordinates)
    REAL(KIND = DP), INTENT(inout) :: ifc
    !! Sum of real space IFC
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: iat, jat, i, j
    !! Counter mode index
    INTEGEr :: ierr
    !! Error status
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP), ALLOCATABLE :: dynq_tmp(:, :, :)
    !! dynamical matrix in bloch representation (Cartesian coordinates)
    COMPLEX(KIND = DP), ALLOCATABLE :: rdw(:, :, :)
    !! Real-space IFC
    !
    ALLOCATE(dynq_tmp(nmodes, nmodes, nq), STAT = ierr)
    IF (ierr /= 0) CALL errore('findL', 'Error allocating dynq_tmp', 1)
    ALLOCATE(rdw(nmodes, nmodes, nrr), STAT = ierr)
    IF (ierr /= 0) CALL errore('findL', 'Error allocating rdw', 1)
    !
    dynq_tmp(:, :, :) = dynq(:, :, :)
    !
    ! Subtract the long-range term from D(q)
    DO ik = 1, nq
      !xk has to be in cart. coord.
      CALL rgd_blk(L, nqc1, nqc2, nqc3, nat, dynq_tmp(:, :, ik), xk(:, ik), tau, epsi, zstar, -1.d0)
      !
    ENDDO
    !
    CALL cryst_to_cart(nq, xk, at, -1)
    rdw(:, :, :) = czero
    DO ir = 1, nrr
      DO ik = 1, nq
        rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
        cfac = EXP(-ci * rdotk)/ DBLE(nq)
        rdw(:, :, ir) = rdw(:, :, ir) +  cfac * dynq_tmp(:, :, ik)
      ENDDO
    ENDDO
    CALL cryst_to_cart(nq, xk, bg, 1)
    !
    ifc = 0.d0
    DO ir = 1, nrr
      DO iat = 1, nat
        DO i = 1, 3
          DO jat = 1, nat
            DO j = 1, 3
              ! We have to exclude the l=0, \kappa=\kappa' term.
              IF (iat == jat .AND. ir == INT((nrr + 1) / 2) ) THEN
                continue
              ELSE
                ifc = ifc + ABS(DBLE(rdw(3 * (iat - 1) + i, 3 * (jat - 1) + j, ir))) / 2.0d0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ifc = ifc / nrr
    !
    DEALLOCATE(dynq_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('findL', 'Error deallocating dynq_tmp', 1)
    DEALLOCATE(rdw, STAT = ierr)
    IF (ierr /= 0) CALL errore('findL', 'Error deallocating rdw', 1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE findL
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE vmebloch2wan(dims, nbnd, nbndsub, nks, nkstot, xk, cu, &
                            nrr, irvec, wslen, lwin, exband, A, w_centers, ndegen_k)
    !--------------------------------------------------------------------------
    !!
    !! Calculate the velocity matrix elements in the Wannier basis
    !! at no point do we actually have the coarse mesh v-ME.
    !!
    !! RM 03/2018: debugged and updated
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE global_var,    ONLY : cvmew, nbndep, crrw, qrpl
    USE ep_constants,  ONLY : twopi, one, zero, ci, czero, cone, bohr2ang
    USE input,         ONLY : use_ws
    USE io_var,        ONLY : iummn, iubvec, iudecayv
    USE io_files,      ONLY : prefix
    USE io_global,     ONLY : ionode_id, stdout
    USE mp_global,     ONLY : inter_pool_comm, world_comm
    USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_world,      ONLY : mpime
    USE parallelism,   ONLY : fkbounds
    USE kfold,         ONLY : ktokpmq
    !
    IMPLICIT NONE
    !
    ! Input variables
    LOGICAL, INTENT(in) :: lwin(nbndep,nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    INTEGER, INTENT(in) :: dims
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: nbnd
    !! number of bands
    INTEGER, INTENT(in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nks
    !! number of kpoints in this pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of kpoints
    INTEGER, INTENT(in) :: nrr
    !! number of WS points
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! Coordinate of Wannier space points
    INTEGER, INTENT(in) :: ndegen_k(nrr, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    REAL(KIND = DP), INTENT(in) :: w_centers(3, nbndsub)
    !! Wannier centers
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbndep, nbndsub, nks)
    !! rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(inout) :: A(3, nbndsub, nbndsub, nks)
    !! Berry connection in Wannier representation
    !
    ! Local variables
    CHARACTER(LEN = 256) :: tempfile
    !! Temporary file
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ikstart
    !! Starting ik for this pool
    INTEGER :: ikstop
    !! Ending ik for this pool
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: ibnd1
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: jbnd1
    !! Counter on band index
    INTEGER :: i
    !! Counter on band index
    INTEGER :: j
    !! Counter on band index
    INTEGER :: nnb
    !! total number of neighbours for each k-point
    INTEGER :: ib
    !! Counter on b-vectors
    INTEGER :: ipool
    !! Current pool
    INTEGER :: nkb
    !! in-pool index of k+b-point
    INTEGER :: nkb_abs
    !! absolute index of k+b-point
    INTEGER :: nkk
    !! in-pool index of k-point
    INTEGER :: nkk_abs
    !! absolute index of k-point
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    INTEGER :: nkstot_tmp
    !! Number of k-point to test
    INTEGER :: ierr
    !! Error status
    INTEGER :: iw
    !! Counter on bands when use_ws == .TRUE.
    INTEGER :: iw2
    !! Counter on bands when use_ws == .TRUE.
    REAL(KIND = DP) :: bvec_crys(3)
    !! b vector in crystal coordinate
    REAL(KIND = DP), ALLOCATABLE :: bvec(:, :, :)
    !! b-vectors connecting each k-point to its nearest neighbors
    REAL(KIND = DP), ALLOCATABLE :: wb(:)
    !! weight of the nnb-th b-vector (ordered along shells)
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: cfac2
    !! $$ e^{-ib(rn-rm - R)/2}
    REAL(KIND = DP) ::  b_tmp(3)
    !! temporary b-vectors
    REAL(KIND = DP) :: zero_vect(3)
    !! temporary zero vector
    REAL(KIND = DP) :: delta
    !! \delta_nm = 1 if n == m and 0 if n /= m
    REAL(KIND = DP) :: xk_cart(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP) :: xk_crys(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP) :: w_centers_crys(3, nbndsub)
    !! Wannier centers
    COMPLEX(KIND = DP) :: M_mn_utmp(nbndep, nbndsub)
    !! M_mn after multiplication with the Wannier rotation matrix cu.
    COMPLEX(KIND = DP), ALLOCATABLE :: cu_big(:, :, :)
    !! rotation matrix from wannier code
    COMPLEX(KIND = DP), ALLOCATABLE :: Apos(:, :, :, :)
    !! A^W_{mn,\alpha}(k)
    COMPLEX(KIND = DP), ALLOCATABLE :: M_mn(:, :, :, :)
    !! M_mn(k,b)
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mat_opt(:, :, :, :)
    !! M_mn(k,b) computed in mmn, rescaled down if skipping lwin.
    COMPLEX(KIND = DP), ALLOCATABLE :: cvs(:, :, :, :)
    !! M_mn in smooth Bloch basis, coarse k-mesh
    !
    ALLOCATE(cu_big(nbndep, nbndsub, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating cu_big', 1)
    ALLOCATE(Apos(3, nbndsub, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating Apos', 1)
    !
    xk_cart(:, :) = xk(:, :)
    xk_crys(:, :) = xk(:, :)
    CALL cryst_to_cart(nks, xk_crys, at, -1)
    !
    w_centers_crys(:, :) = w_centers(:, :)
    CALL cryst_to_cart(nbndsub, w_centers_crys, bg, -1)
    !
    ! setup rotation matrix - we need access to all for the k+b
    cu_big = czero
    CALL fkbounds(nkstot, ikstart, ikstop)
    cu_big(:, :, ikstart:ikstop) = cu(:, :, :)
    CALL mp_sum(cu_big, inter_pool_comm)
    !
    !--------------------------------------------------------------
    !  STEP 0: Read in b-vectors bvec and their weights wb, and
    !          M_mn matrix
    !--------------------------------------------------------------
    !
    ! RM - bvec are writen on file if write_bvec = .true. is specified
    !      in the wannier input
    !
    IF (mpime == ionode_id) THEN
      tempfile = TRIM(prefix) // '.bvec'
      OPEN(iubvec, FILE = tempfile, ACTION = 'read', IOSTAT = ios)
      IF (ios /= 0) THEN
        !
        CALL errore ('vmebloch2wan','You selected vme =.true. but error opening' // tempfile, 1)
      ELSE
        READ(iubvec,*) tempfile
        READ(iubvec,*) nkstot_tmp, nnb
        IF (nkstot_tmp /= nkstot) CALL errore('vmebloch2wan', 'Unexpected number of k-points in .bvec file', 1)
        ALLOCATE(bvec(3, nnb, nkstot), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating bvec', 1)
        ALLOCATE(wb(nnb), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating wb', 1)
        DO ik = 1, nkstot
          DO ib = 1, nnb
            READ(iubvec,*) bvec(:, ib, ik), wb(ib)
          ENDDO
        ENDDO
        CLOSE(iubvec)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(nnb, ionode_id, world_comm)
    !
    ! All other cpu than master
    IF (mpime /= ionode_id) THEN
      ALLOCATE(bvec(3, nnb, nkstot), STAT = ierr)
      IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating bvec', 1)
      ALLOCATE(wb(nnb), STAT = ierr)
      IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating wb', 1)
    ENDIF
    !
    CALL mp_bcast(bvec, ionode_id, world_comm)
    CALL mp_bcast(wb, ionode_id, world_comm)
    !
    ! b-vectors in units of Ang^-1 and wb-weights in units of Ang^2
    ! when read from file
    bvec = bvec * bohr2ang
    wb = wb / bohr2ang**2
    !
    !  read Mmn for velocity calculation
    !
    ALLOCATE(M_mn(nbndep, nbndep, nnb, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating M_mn', 1)
    M_mn = czero
    !
    IF (mpime == ionode_id) THEN
      tempfile = TRIM(prefix)//'.mmn'
      OPEN(iummn, FILE = tempfile, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
      !
      IF (ios /= 0) THEN
        ! if it doesn't exist, then we just set the mmn to zero
        CALL errore ('vmebloch2wan','error opening' // tempfile, 1)
      ELSE
        !
        DO ik = 1, nkstot
          DO ib = 1, nnb
            DO jbnd = 1, nbndep
              DO ibnd = 1, nbndep
                READ(iummn,*) M_mn(ibnd, jbnd, ib, ik)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iummn)
        !
      ENDIF
    ENDIF
    CALL mp_bcast(M_mn, ionode_id, world_comm)
    !
    CALL start_clock('Velocity: step 1')
    !
    ! slim down M_mn matrix to contain states within the outer window
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF
    ENDDO
    !
    ALLOCATE(m_mat_opt(nbndep, nbndep, nnb, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating m_mat_opt', 1)
    m_mat_opt(:, :, :, :) = czero
    zero_vect(:) = zero
    !
    IF (nexband_tmp > 0) THEN
      DO ik = 1, nks
        CALL ktokpmq(xk_cart(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        jbnd = 0
        jbnd1 = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          jbnd1 = jbnd1 + 1
          IF (lwin(jbnd1,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            ibnd1 = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              ibnd1 = ibnd1 + 1
              IF (lwin(ibnd1,ik)) THEN
                ibnd = ibnd + 1
                m_mat_opt(ibnd, jbnd, :, ik) = M_mn(ibnd1, jbnd1, :, nkk_abs)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        CALL ktokpmq(xk_cart(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        jbnd = 0
        DO j = 1, nbndep
          IF (lwin(j, ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbndep
              IF (lwin(i, ik)) THEN
                ibnd = ibnd + 1
                m_mat_opt(ibnd, jbnd, :, ik) = M_mn(i, j, :, nkk_abs)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    DEALLOCATE(M_mn, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating M_mn', 1)
    !
    !----------------------------------------------------------
    ! STEP 1: Calculate A^(W)_{mn,\alpha}(k) [Eqn. 44 of PRB 74, 195118 (2006)]
    ! A^(W)_{mn,\alpha}(k) = i \sum_b wb b_{\alpha} (<u^(W)_nk | u^(W)_m(k+b)> - delta_mn)
    !----------------------------------------------------------
    !
    ! <u^(W)_nk | u^(W)_m(k+b)> = U(k)^\dagger * M_mn(k) * U(k+b)
    !
    ! <u^(W)_nk | u^(W)_m(k+b)> is cvs(nbndsub,nbndsub,nnb,nks)
    !
    ALLOCATE(cvs(nbndsub, nbndsub, nnb, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating cvs', 1)
    !
    cvs(:, :, :, :) = czero
    M_mn_utmp(:, :) = czero
    !
    DO ik = 1, nks
      !
      CALL ktokpmq(xk_cart(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
      !
      DO ib = 1, nnb
        !
        ! bring bvec to units of 2piba since xk is cartesian units of 2piba
        b_tmp(:) = alat / (twopi) * bvec(:, ib, nkk_abs)
        CALL ktokpmq(xk_cart(:, ik), b_tmp(:), +1, ipool, nkb, nkb_abs)
        !
        ! M_mn_utmp(:, :) = MATMUL( m_mat_opt(:,:,ib,ik), cu_big(:,:,nkb_abs) )
        ! cvs(:,:,ib,ik) = MATMUL( CONJG(transpose(cu(:,:,ik))), M_mn_utmp(:, :) )
        !
        CALL ZGEMM('n', 'n', nbndep, nbndsub, nbndep, cone, m_mat_opt(:, :, ib, ik), &
                   nbndep, cu_big(:, :, nkb_abs), nbndep, czero, M_mn_utmp(:, :), nbndep)
        CALL ZGEMM('c', 'n', nbndsub, nbndsub, nbndep, cone, cu(:, :, ik), &
                   nbndep, M_mn_utmp(:, :), nbndep, czero, cvs(:, :, ib, ik), nbndsub)
      ENDDO
      !
    ENDDO
    !
    DEALLOCATE(m_mat_opt, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating m_mat_opt', 1)
    !
    ! A^(W)_{mn,\alpha}(k) is Apos(3,nbndsub,nbndsub,nks)
    !
    Apos = czero
    !
    DO ik = 1, nks
      !
      CALL ktokpmq(xk_cart(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
      !
      DO ib = 1, nnb
        !
        DO jbnd = 1, nbndsub
          DO ibnd = 1, nbndsub
            !
            delta = zero
            IF (jbnd == ibnd) delta = one
            DO ipol = 1, 3
              Apos(ipol, ibnd, jbnd, ik) = Apos(ipol, ibnd, jbnd, ik) + &
                 ci * wb(ib) * bvec(ipol, ib, nkk_abs) * (cvs(ibnd, jbnd, ib, ik) - delta)
            ENDDO
            !
          ENDDO
        ENDDO
        !
      ENDDO
      !
    ENDDO
    !
    !
    CALL stop_clock('Velocity: step 1')
    !
    WRITE(stdout,'(/5x,a)') 'Inside velocity step 1'
    WRITE(stdout,'(a)') ' '
    !
    !----------------------------------------------------------
    ! STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! Calculate position matrix [Eqn. 43 of PRB 74, 195118 (2006)]
    ! r_{\alpha}(R) = (1/nk) \sum_k e^(-ikR) A^(W)_{mn,\alpha}(k)
    !
    CALL start_clock('Velocity: step 2')
    !
    ! bring xk in crystal coordinates
    !
    !CALL cryst_to_cart(nks, xk, at, -1)
    !
    ! r_{\alpha}(R) is cvmew(3,nbndsub,nbndsub,nrr)
    !
    cvmew(:, :, :, :) = czero
    !
    DO ir = 1, nrr
      DO ik = 1, nks
        !
        rdotk = twopi * DOT_PRODUCT(xk_crys(:, ik), DBLE(irvec(:, ir)))
        cfac = EXP(-ci * rdotk) / DBLE(nkstot)
        cvmew(:, :, :, ir) = cvmew(:, :, :, ir) + cfac * Apos(:, :, :, ik)
        !
      ENDDO
    ENDDO
    !
    CALL mp_sum(cvmew, inter_pool_comm)
    !
    ! Now compute the position operator r_ijR
    !print*,'(INT(nrr/2) + MOD(nrr,2)) ',(INT(nrr/2) + MOD(nrr,2))
    !print*,'irvec(:, ir) ',irvec(:, (INT(nrr/2) + MOD(nrr,2)))
    !
    crrw(:,:,:,:) = czero
    A(:, :, :, :) = czero
    ! The Berry connection is only computed in case of quadrupole
    IF (qrpl) THEN
      DO ir = 1, nrr
        DO ik = 1, nks
          CALL ktokpmq(xk_cart(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
          rdotk = twopi * DOT_PRODUCT(xk_crys(:, ik), DBLE(irvec(:, ir)))
          cfac = EXP(-ci * rdotk) / DBLE(nkstot)
          DO ib = 1, nnb
            bvec_crys(:) = bvec(:, ib, nkk_abs) * alat / (twopi)
            CALL cryst_to_cart(1, bvec_crys, at, -1)

            DO jbnd = 1, nbndsub
              DO ibnd = 1, nbndsub
                DO ipol = 1, 3
                  cfac2 = EXP(twopi * ci * DOT_PRODUCT(bvec_crys(:), &
                             (w_centers_crys(:, ibnd) + w_centers_crys(:, jbnd) - DBLE(irvec(:, ir)))/2.0d0))
                  crrw(ipol, ibnd, jbnd, ir) = crrw(ipol, ibnd, jbnd, ir) + ci * cfac * wb(ib) * bvec(ipol, ib, nkk_abs) &
                                           * cfac2 * cvs(ibnd, jbnd, ib, ik)
                ENDDO ! ipol
              ENDDO ! ibnd
            ENDDO ! jbnd
          ENDDO ! ib
        ENDDO ! ik
      ENDDO ! ir
      !
      CALL mp_sum(crrw, inter_pool_comm)
      !
      ! For the diagonal part at R=\Gamma, which is half the nrr
      DO ibnd = 1, nbndsub
        DO ipol = 1, 3
          crrw(ipol, ibnd, ibnd, (INT(nrr/2) + MOD(nrr,2))) = w_centers(ipol, ibnd) * alat
        ENDDO ! ipol
      ENDDO ! ibnd
      !
      !print*,'crrw(3, 2, 2, R=0) ',crrw(3, 2, 2, (INT(nrr/2) + MOD(nrr,2)))
      !write(100,'(a)') 'ir   ibnd  jbnd  ipol   Re r_ij   Im r_ij '
      !DO ir=1, nrr
      !  DO ibnd = 1, nbndsub
      !    DO jbnd = 1, nbndsub
      !      DO ipol = 1, 3
      !        write(100,'(4i9, 2f20.10)') ir, ibnd, jbnd, ipol, REAL(crrw(ipol, ibnd, ibnd, ir)), AIMAG(crrw(ipol, ibnd, ibnd, ir))
      !      ENDDO
      !    ENDDO
      !  ENDDO
      !ENDDO
      !
      ! Fourier transform (Wannier to Bloch) crrw -> A_nmk^(W)
      DO ik = 1, nks
        IF (use_ws) THEN
          DO iw = 1, dims
            DO iw2 = 1, dims
              DO ir = 1, nrr
                IF (ndegen_k(ir, iw2, iw) > 0) THEN
                  rdotk = twopi * DOT_PRODUCT(xk_crys(:, ik), DBLE(irvec(:, ir)))
                  cfac = EXP(ci * rdotk) / ndegen_k(ir, iw2, iw)
                  A(:, iw, iw2, ik) = A(:, iw, iw2, ik) + cfac * crrw(:, iw, iw2, ir)
                ENDIF
              ENDDO ! ir
            ENDDO
          ENDDO ! iw
        ELSE ! use_ws
          DO ir = 1, nrr
            rdotk = twopi * DOT_PRODUCT(xk_crys(:, ik), DBLE(irvec(:, ir)))
            cfac = EXP(ci * rdotk) / ndegen_k(ir, 1, 1)
            A(:, :, :, ik) = A(:, :, :, ik) + cfac * crrw(:, :, :, ir)
          ENDDO
        ENDIF ! use_ws
      ENDDO ! ik
    ENDIF ! qrpl
    !
    !
    DEALLOCATE(bvec, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating bvec', 1)
    DEALLOCATE(wb, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating wb', 1)
    !
    ! check spatial decay of position matrix elements in Wannier basis
    !
    ! position matrix cvmew and spatial dimensions are in units of bohr
    ! [mind when comparing with wannier code (angstrom units) with write_rmn=.TRUE.]
    !
    IF (mpime == ionode_id) then
      OPEN(UNIT = iudecayv, FILE = 'decay.v')
      WRITE(iudecayv, '(a)') '# Spatial decay of Velocity matrix element in Wannier basis'
      WRITE(iudecayv, '(a)') '# R_e [Ang]       max_{a,m,n} |v(a,m,n)| [Bohr/2.418E-17s]    '
      DO ir = 1, nrr
        !
        tmp =  MAXVAL(ABS(cvmew(:, :, :, ir)))
        WRITE(iudecayv, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iudecayv)
    ENDIF
    IF (mpime == ionode_id) then
      OPEN(UNIT = iudecayv, FILE = 'decay.r')
      WRITE(iudecayv, '(a)') '# Spatial decay of Position matrix element in Wannier basis'
      WRITE(iudecayv, '(a)') '# R_e [Ang]       max_{a,m,n} |v(a,m,n)| [Bohr/2.418E-17s]    '
      DO ir = 1, nrr
        !
        tmp =  MAXVAL(ABS(crrw(:, :, :, ir)))
        WRITE(iudecayv, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iudecayv)
    ENDIF
    !
    CALL stop_clock('Velocity: step 2')
    !
    WRITE(stdout,'(/5x,a)') 'Velocity matrix elements calculated'
    WRITE(stdout,'(a)') ' '
    !
    DEALLOCATE(cu_big, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating cu_big', 1)
    DEALLOCATE(Apos, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating Apos', 1)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE vmebloch2wan
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE ephbloch2wane(iq, xq, nbnd, nbndsub, nmodes, nks, nkstot, xk, &
         cu, cuq, epmatk, nrr, irvec, wslen, epmatw, A, is_dw)
    !-----------------------------------------------------------------------
    !!
    !!  From the EP matrix elements in Bloch representation (coarse
    !!  mesh), find the corresponding matrix elements in electron-Wannier
    !!  representation and phonon-Bloch representation
    !!
    !-----------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero, cone, eps8, one
    USE io_global,     ONLY : ionode_id
    USE io_var,        ONLY : iundwdecay
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_sum
    USE mp_world,      ONLY : mpime
    USE global_var,    ONLY : zstar, epsi, qrpl
    USE input,         ONLY : lpolar
    USE longrange,     ONLY : rgd_blk_epw
    USE input,         ONLY : nqc1, nqc2, nqc3
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: iq
    !! Q-point index.
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands
    INTEGER, INTENT(in) :: nbndsub
    !! Number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nmodes
    !! Number of modes
    INTEGER, INTENT(in) :: nks
    !! Number of kpoints in this pool
    INTEGER, INTENT(in) :: nrr
    !! Number of WS points
    INTEGER, INTENT(in) :: nkstot
    !! Total number of kpoints
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! Coordinates of WS points
    !
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! Current coarse q-point
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(KIND = DP), INTENT(inout) :: cu(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(inout) :: cuq(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(inout) :: epmatk(nbnd, nbnd, nks, nmodes)
    !! e-p matrix in bloch representation, coarse mesh
    COMPLEX(KIND = DP), INTENT(inout), OPTIONAL :: A(3, nbndsub, nbndsub, nks)
    !! Berry connection in Wannier representation
    !
    LOGICAL, INTENT(in), OPTIONAL :: is_dw
    !! If .true., Fourier transforming the Debye-Waller matrix element
    !
    ! output variables
    !
    COMPLEX(KIND = DP), INTENT(out) :: epmatw(nbndsub, nbndsub, nrr, nmodes)
    ! EP vertex (Wannier el and Bloch ph)
    !
    ! Work variables
    !
    LOGICAL :: is_dw_
    !! If .true., the input is Debye-Waller matrix element. Default is false.
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: imode
    !! Mode index
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Max electron-phonon value
    !
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: eptmp(nbndsub, nbnd)
    !!  e-p matrix, temporary
    COMPLEX(KIND = DP), ALLOCATABLE :: epmats(:, :, :, :)
    !!  e-p matrix  in smooth Bloch basis, coarse mesh
    !
    ALLOCATE(epmats(nbndsub, nbndsub, nks, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wane', 'Error allocating epmats', 1)
    !
    is_dw_ = .FALSE.
    IF (PRESENT(is_dw)) is_dw_ = is_dw
    !
    !----------------------------------------------------------
    !  STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g~(k,q) = U(k+q)^\dagger * g(k,q) * U(k)
    !
    ! g(k,q)   is epmatk (ibnd, jbnd, ik)
    ! g~(k,q)  is epmats (ibnd, jbnd, ik)
    !
    CALL start_clock ('ep: step 1')
    !
    epmats(:, :, :, :) = czero
    DO imode = 1, nmodes
      DO ik = 1, nks
        !
        ! the two zgemm calls perform the following ops:
        ! epmats  = [ cu(ikq)^\dagger * epmatk ] * cu(ikk)
        ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub]
        !
        CALL ZGEMM('c', 'n', nbndsub, nbnd, nbnd, cone, cuq(:, :, ik),  &
                  nbnd, epmatk(:, :, ik, imode), nbnd, czero, eptmp, nbndsub)
        CALL ZGEMM('n', 'n', nbndsub, nbndsub, nbnd, cone, eptmp,     &
                  nbndsub, cu(:, :, ik), nbnd, czero, epmats(:, :, ik, imode), nbndsub)
        !
      ENDDO
    ENDDO
    !
    CALL stop_clock ('ep: step 1')
    !
    !----------------------------------------------------------
    ! STEP 2: First remove long range.
    !----------------------------------------------------------
    !
    IF ((.NOT. is_dw_) .AND. (lpolar .OR. qrpl)) THEN
      DO ik = 1, nks
        CALL rgd_blk_epw(nqc1, nqc2, nqc3, xq, epmats(:, :, ik, :), nbndsub, nmodes, epsi, zstar, A(:, :, :, ik), -one)
      ENDDO
    ENDIF
    !
    !----------------------------------------------------------------------
    !  STEP 3: Fourier transform to obtain matrix elements in electron wannier basis
    !----------------------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,q) = (1/nkc) sum_k e^{-ikR_e} g~(k,q)
    ! g(R_e,q) is epmatw (nbndsub,nbndsub,ir)
    !
    CALL start_clock('ep: step 2')
    !
    epmatw(:, :, :, :) = czero
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nks, xk, at, -1)
    !
    DO imode = 1, nmodes
      DO ir = 1, nrr
        DO ik = 1, nks
          !
          rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
          cfac = EXP(-ci * rdotk) / DBLE(nkstot)
          epmatw( :, :, ir, imode) = epmatw( :, :, ir, imode) + cfac * epmats(:, :, ik, imode)
          !
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum(epmatw, inter_pool_comm)
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nks, xk, bg, 1)
    !
    !
    !  Check spatial decay of EP matrix elements in electron-Wannier basis
    !  the unit in r-space is angstrom, and I am plotting
    !  the matrix for the first mode only
    !
    IF ((mpime == ionode_id) .AND. is_dw_) THEN
      OPEN(UNIT = iundwdecay, FILE = 'decay.dwwane')
      WRITE(iundwdecay, '(a)') '# Spatial decay of Debye-Waller matrix elements in electron Wannier basis'
      DO ir = 1, nrr
        !
        tmp = MAXVAL(ABS(epmatw(:, :, ir, :)))
        WRITE(iundwdecay, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iundwdecay)
    ENDIF
    !
    DEALLOCATE(epmats, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wane', 'Error deallocating epmats', 1)
    !
    CALL stop_clock ('ep: step 2')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE ephbloch2wane
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE ephbloch2wanp(nbnd, nmodes, xk, nq, irvec_k, irvec_g, nrr_k, nrr_g, epmatwe)
    !--------------------------------------------------------------------------
    !!
    !!  From the EP matrix in electron-Wannier representation and
    !!  phonon-Bloch representation (coarse mesh), find the corresponding matrix
    !!  electron-Wannier representation and phonon-Wannier representation
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE global_var,    ONLY : epmatwp
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero
    USE io_var,        ONLY : iuwanep, iuwane
    USE mp,            ONLY : mp_barrier
    !
    IMPLICIT NONE
    !
    !  Input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of electronic bands
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq
    !! number of qpoints
    INTEGER, INTENT(in) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nrr_g
    !! number of el-ph WS points
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Coordinates of real space vector for electrons
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of real space vector for electron-phonon
    !
    REAL(KIND = DP), INTENT(in) :: xk(3, nq)
    !! Kpoint coordinates (cartesian in units of 2piba)
    !
    COMPLEX(KIND = DP), INTENT(in) :: epmatwe(nbnd, nbnd, nrr_k, nmodes, nq)
    !! EP matrix in electron-Wannier representation and phonon-Bloch representation
    !!   (Cartesian coordinates)
    !
    ! Work variables
    !
    INTEGER :: iq
    !! Counter on q-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ire
    !! Counter on WS points
    !
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp(nrr_k, nrr_g)
    !! Max electron-phonon value
    REAL(KIND = DP) :: rvec1(3)
    !! WS vectors
    REAL(KIND = DP) :: rvec2(3)
    !! WS vectors
    REAL(KIND = DP) :: len1(nrr_k)
    !! Electron distance when printing the decay files
    REAL(KIND = DP) :: len2(nrr_g)
    !! Phonon distance when printing the decay files
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    !
    !----------------------------------------------------------
    !  Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,R_p) = (1/nq) sum_q e^{-iqR_p} g(R_e,q)
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nq, xk, at, -1)
    !
    epmatwp = czero
    !
    DO ir = 1, nrr_g
      !
      DO iq = 1, nq
        !
        rdotk = twopi * DOT_PRODUCT(xk(:, iq), DBLE(irvec_g(:, ir)))
        cfac = EXP(-ci * rdotk) / DBLE(nq)
        epmatwp(:, :, :, :, ir) = epmatwp(:, :, :, :, ir) + cfac * epmatwe(:, :, :, :, iq)
        !
      ENDDO
      !
      ! Check spatial decay of EP matrix elements in wannier basis - electrons + phonons
      !
      rvec2 = DBLE(irvec_g(1, ir)) * at(:, 1) + DBLE(irvec_g(2, ir)) * at(:, 2) + DBLE(irvec_g(3, ir)) * at(:, 3)
      ! phonon - electron0 distance
      len2(ir) = DSQRT(rvec2(1)**2.d0 + rvec2(2)**2.d0 + rvec2(3)**2.d0)
      DO ire = 1, nrr_k
        tmp(ire, ir)  = MAXVAL(ABS(epmatwp(:, :, ire, :, ir)))
      ENDDO
      !
    ENDDO ! ir
    !
    !  Check spatial decay of EP matrix elements in wannier basis - electrons + phonons
    !  We plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
    ! electron-electron0 distance
    DO ire = 1, nrr_k
      rvec1 = DBLE(irvec_k(1, ire)) * at(:, 1) + DBLE(irvec_k(2, ire)) * at(:, 2) + DBLE(irvec_k(3, ire)) * at(:, 3)
      len1(ire) = DSQRT(rvec1(1)**2.d0 + rvec1(2)**2.d0 + rvec1(3)**2.d0)
    ENDDO
    !
    OPEN(UNIT = iuwane, FILE = 'decay.epmate', STATUS = 'unknown')
    WRITE(iuwane, '(a)') '#   R_e [Ang]    max_{m,n,nu} |g(m, n, nu, R_e, :)| [Ry] '
    DO ire = 1, nrr_k
      WRITE(iuwane, '(2f15.10, 1E20.10)') len1(ire) * alat * bohr2ang, MAXVAL(tmp(ire, :))
    ENDDO
    CLOSE(iuwane)
    !
    OPEN(UNIT = iuwanep, FILE = 'decay.epmatp', STATUS = 'unknown')
    WRITE(iuwanep, '(a)') '#   R_p [Ang]    max_{m,n,nu} |g(m, n, nu, :, R_p)| [Ry] '
    DO ir = 1, nrr_g
      WRITE(iuwanep, '(2f15.10, 1E20.10)') len2(ir) * alat * bohr2ang, MAXVAL(tmp(:, ir))
    ENDDO
    CLOSE(iuwanep)
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nq, xk, bg, 1)
    !
    ! -----------------------------------------------------------
    END SUBROUTINE ephbloch2wanp
    ! -----------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE ephbloch2wanp_mem(nbnd, nmodes, xk, nq, irvec_k, irvec_g, nrr_k, nrr_g)
    !--------------------------------------------------------------------------
    !
    !!  From the EP matrix in electron-Wannier representation and
    !!  phonon-Bloch representation (coarse mesh), find the corresponding matrix
    !!  electron-Wannier representation and phonon-Wannier representation
    !!
    !!  SP - June 2020 - Update decays printing.
    !!  SP - Apr 2020 - MPI-IO parallelization.
    !
    USE kinds,            ONLY : DP
    USE cell_base,        ONLY : at, bg, alat
    USE ep_constants,     ONLY : bohr2ang, twopi, ci, czero, zero
    USE io_var,           ONLY : iunepmatwe, iunepmatwp, iuwanep, iuwane
    USE io_global,        ONLY : ionode_id, stdout
    USE mp_global,        ONLY : world_comm
    USE mp,               ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_world,         ONLY : mpime
    USE io,               ONLY : rwepmatw
    USE parallelism,      ONLY : para_bounds
    USE io_files,         ONLY : prefix, diropn, tmp_dir
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, MPI_MODE_RDONLY, &
                                 MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, &
                                 MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_INFO_NULL, &
                                 MPI_MODE_DELETE_ON_CLOSE
#endif
    !
    IMPLICIT NONE
    !
    ! Input variables - note irvec is dimensioned with nrr_k (which is assumed to be larger than nrr_q)
    INTEGER, INTENT(in) :: nbnd
    !! Number of electronic bands
    INTEGER, INTENT(in) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nrr_g
    !! number of el-h WS points
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq
    !! number of qpoints
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Coordinates of real space vector
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of real space vector
    REAL(KIND = DP), INTENT(in) :: xk(3, nq)
    !! K-point coordinates (cartesian in units of 2piba)
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    INTEGER :: iq
    !! Counter on q-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ire
    !! Counter on WS points
    INTEGER :: ir_start
    !! Starting ir for this cores
    INTEGER :: ir_stop
    !! Ending ir for this pool
    INTEGER :: ierr
    !! Error status
    INTEGER :: diff
    !! Difference between starting and ending on master core
    INTEGER :: add
    !! Additional element
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER(KIND = MPI_OFFSET_KIND) :: lsize
    !! Offset to tell where to start reading the file
#else
    INTEGER(KIND = 8) :: lrepmatw
    !! Offset to tell where to start reading the file
    INTEGER(KIND = 4) :: lsize
    !! Offset to tell where to start reading the file
#endif
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp(nrr_k, nrr_g)
    !! Max electron-phonon value
    REAL(KIND = DP) :: rvec1(3)
    !! WS vectors
    REAL(KIND = DP) :: rvec2(3)
    !! WS vectors
    REAL(KIND = DP) :: len1(nrr_k)
    !! Electron distance when printing the decay files
    REAL(KIND = DP) :: len2(nrr_g)
    !! Phonon distance when printing the decay files
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwe(:, :, :, :)
    !! EP matrix in electron-Wannier representation and phonon-Bloch representation
    !! (Cartesian coordinates)
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwp_mem(:, :, :, :)
    !!  e-p matrix in Wannier basis
    LOGICAL :: exst
    !! Check if backup files exist
    !
    ALLOCATE(epmatwe(nbnd, nbnd, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'Error allocating epmatwe', 1)
    ALLOCATE(epmatwp_mem(nbnd, nbnd, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'Error allocating epmatwp_mem', 1)
    !
    len1(:) = zero
    len2(:) = zero
    tmp(:, :) = zero
    !----------------------------------------------------------
    !  Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,R_p) = (1/nq) sum_q e^{-iqR_p} g(R_e,q)
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nq, xk, at, -1)
    !
    ! Distribute the cpu
    CALL para_bounds(ir_start, ir_stop, nrr_g)
    !
    IF (mpime == ionode_id) THEN
      diff = ir_stop - ir_start
    ENDIF
    CALL mp_bcast(diff, ionode_id, world_comm)
    !
    ! If you are the last cpu with less element
    IF (ir_stop - ir_start /= diff) THEN
      add = 1
    ELSE
      add = 0
    ENDIF
    !
#if defined(__MPI)
    ! Size of the read array
    lsize = 2_MPI_OFFSET_KIND * INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                INT(nrr_k, KIND = MPI_OFFSET_KIND) * &
                                INT(nmodes, KIND = MPI_OFFSET_KIND)
    !
    ! Open the epmatwe file
    filint = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwe1'
    CALL MPI_FILE_OPEN(world_comm, filint, MPI_MODE_RDONLY + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, iunepmatwe, ierr)
    !CALL MPI_FILE_OPEN(world_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatwe, ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'error in MPI_FILE_OPEN epmatwe', 1)
    !
    ! Open the epmatwp file
    filint = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
    CALL MPI_FILE_OPEN(world_comm, filint, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iunepmatwp, ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'error in MPI_FILE_OPEN epmatwp', 1)
#else
    ! Size of the read array
    lsize = INT(2 * nbnd * nbnd * nrr_k * nmodes, KIND = 4)
    filint   = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwe'
    CALL diropn(iunepmatwe, 'epmatwe', lsize, exst)
    IF (.NOT. exst) CALL errore('ephbloch2wanp_mem', 'file ' // TRIM(filint) // ' not found', 1)
    !
    CALL diropn(iunepmatwp, 'epmatwp', lsize, exst)
#endif
    !
    DO ir = ir_start, ir_stop + add
      WRITE(stdout, '(a,i10,a,i10)' ) '     Bloch2wanp: ',ir - ir_start + 1,' / ', ir_stop + add - ir_start + 1
      !
#if defined(__MPI)
      IF (add == 1 .AND. ir == ir_stop + add) lsize = 0_MPI_OFFSET_KIND
#endif
      !
      epmatwp_mem = czero
      !
      DO iq = 1, nq
        !
#if defined(__MPI)
        lrepmatw = 2_MPI_OFFSET_KIND * 8_MPI_OFFSET_KIND * &
                                     INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                     INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                     INT(nrr_k, KIND = MPI_OFFSET_KIND) * &
                                     INT(nmodes, KIND = MPI_OFFSET_KIND) * &
                                    (INT(iq, KIND = MPI_OFFSET_KIND) - 1_MPI_OFFSET_KIND)
        !
        ! Parallel read using MPI-IO of epmatwe for this iq
        epmatwe = czero
        CALL MPI_FILE_READ_AT(iunepmatwe, lrepmatw, epmatwe, lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
        IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'error in MPI_FILE_READ_AT', 1)
        IF (add == 1 .AND. ir == ir_stop + add) CYCLE
#else
        epmatwe = czero
        CALL rwepmatw(epmatwe, nbnd, nrr_k, nmodes, iq, iunepmatwe, -1)
#endif
        !
        !
        rdotk = twopi * DOT_PRODUCT(xk(:, iq), DBLE(irvec_g(:, ir)))
        cfac = EXP(-ci * rdotk) / DBLE(nq)
        epmatwp_mem = epmatwp_mem + cfac * epmatwe
        !
      ENDDO
      !
#if defined(__MPI)
      IF (add == 1 .AND. ir == ir_stop + add) CYCLE
      lrepmatw = 2_MPI_OFFSET_KIND * 8_MPI_OFFSET_KIND * &
                                   INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                   INT(nbnd , KIND = MPI_OFFSET_KIND) * &
                                   INT(nrr_k, KIND = MPI_OFFSET_KIND) * &
                                   INT(nmodes, KIND = MPI_OFFSET_KIND) * &
                                  (INT(ir, KIND = MPI_OFFSET_KIND) - 1_MPI_OFFSET_KIND)

      CALL MPI_FILE_WRITE_AT(iunepmatwp, lrepmatw, epmatwp_mem, lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'error in MPI_FILE_WRITE_AT', 1)
#else
      ! direct write of epmatwp_mem for this ir
      CALL rwepmatw(epmatwp_mem, nbnd, nrr_k, nmodes, ir, iunepmatwp, +1)
#endif
      !
      ! Check spatial decay of EP matrix elements in wannier basis - electrons + phonons
      !
      rvec2 = DBLE(irvec_g(1, ir)) * at(:, 1) + DBLE(irvec_g(2, ir)) * at(:, 2) + DBLE(irvec_g(3, ir)) * at(:, 3)
      ! phonon - electron0 distance
      len2(ir) = DSQRT(rvec2(1)**2.d0 + rvec2(2)**2.d0 + rvec2(3)**2.d0)
      DO ire = 1, nrr_k
        tmp(ire, ir)  = MAXVAL(ABS(epmatwp_mem(:, :, ire, :)))
      ENDDO
      !
    ENDDO ! ir
    !
    !  Check spatial decay of EP matrix elements in wannier basis - electrons + phonons
    !  We plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
    CALL mp_sum(len2, world_comm)
    CALL mp_sum(tmp, world_comm)
    IF (mpime == ionode_id) THEN
      ! electron-electron0 distance
      DO ire = 1, nrr_k
        rvec1 = DBLE(irvec_k(1, ire)) * at(:, 1) + DBLE(irvec_k(2, ire)) * at(:, 2) + DBLE(irvec_k(3, ire)) * at(:, 3)
        len1(ire) = DSQRT(rvec1(1)**2.d0 + rvec1(2)**2.d0 + rvec1(3)**2.d0)
      ENDDO
      !
      OPEN(UNIT = iuwane, FILE = 'decay.epmate', STATUS = 'unknown')
      WRITE(iuwane, '(a)') '#   R_e [Ang]    max_{m,n,nu} |g(m, n, nu, R_e, :)| [Ry] '
      DO ire = 1, nrr_k
        WRITE(iuwane, '(2f15.10, 1E20.10)') len1(ire) * alat * bohr2ang, MAXVAL(tmp(ire, :))
      ENDDO
      CLOSE(iuwane)
      !
      OPEN(UNIT = iuwanep, FILE = 'decay.epmatp', STATUS = 'unknown')
      WRITE(iuwanep, '(a)') '#   R_p [Ang]    max_{m,n,nu} |g(m, n, nu, :, R_p)| [Ry] '
      DO ir = 1, nrr_g
        WRITE(iuwanep, '(2f15.10, 1E20.10)') len2(ir) * alat * bohr2ang, MAXVAL(tmp(:, ir))
      ENDDO
      CLOSE(iuwanep)
    ENDIF
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nq, xk, bg, 1)
    !
#if defined(__MPI)
    CALL MPI_FILE_CLOSE(iunepmatwe, ierr)
    CALL MPI_FILE_CLOSE(iunepmatwp, ierr)
#else
    CLOSE(iunepmatwe, STATUS = 'delete')
    CLOSE(iunepmatwp, STATUS = 'keep')
#endif
    !
    DEALLOCATE(epmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'Error deallocating epmatwe', 1)
    DEALLOCATE(epmatwp_mem, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'Error deallocating epmatwp_mem', 1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE ephbloch2wanp_mem
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dgbloch2wane(nbnd, nbndsub, nks, nkstot, etk, etq, &
      ahc_win_min, ahc_win_max, xk, xq, cu, cuq, epmatk, nrr, &
      irvec, wslen, dgmatw)
    !--------------------------------------------------------------------------
    !!
    !! From the electron-phonon matrix elements in Bloch representation
    !! (coarse mesh), find the deltaH elements in Wannier representation
    !!
    !! dg_pq(k, q) = U(k+q)^dagger_pm * E_{m,k+q} * Q(k+q)_mp
    !!             * g(k,q)_pn / (E_{n,k} - E_{p,k+q}) * U(k)_nq
    !!
    !! Sum over n is done only for bands inside the frozen window
    !!
    !! Only bands inside the outer window are considered. Bands outside the
    !! outer window need not be included. They are already removed from the
    !! input matrices.
    !!
    !! dg_pq(R_e, q) = (1/nkc) sum_k e^{-ikR_e} dg_pq(k, q)
    !!
    !! Q(k+q) is a projection out of the Wannier subspace:
    !! Q(k+q)_mn = delta_mn - sum_p U(k+q)_np U(k+q)^dagger_pm
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero, cone, rytoev
    USE mp,            ONLY : mp_sum
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : inter_pool_comm
    USE io_global,     ONLY : ionode_id
    USE io_var,        ONLY : iuwane
    USE cell_base,     ONLY : at, alat
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands (nbndep)
    INTEGER, INTENT(in) :: nbndsub
    !! Number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nks
    !! Number of kpoints in this pool
    INTEGER, INTENT(in) :: nrr
    !! Number of WS points
    INTEGER, INTENT(in) :: nkstot
    !! Total number of kpoints
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! Coordinates of WS points
    !
    REAL(KIND = DP), INTENT(in) :: etk(nbnd, nks)
    !! Hamiltonian eigenvalues at k within the outer window
    REAL(KIND = DP), INTENT(in) :: etq(nbnd, nks)
    !! Hamiltonian eigenvalues at k+q within the outer window
    REAL(KIND = DP), INTENT(in) :: ahc_win_min
    !! AHC window of Wannierization (in eV)
    REAL(KIND = DP), INTENT(in) :: ahc_win_max
    !! AHC window of Wannierization (in eV)
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! qpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(in) :: cuq(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(in) :: epmatk(nbnd, nbnd, nks)
    !! e-p matrix in bloch representation, coarse mesh
    !
    ! output variables
    !
    COMPLEX(KIND = DP), INTENT(out) :: dgmatw(nbndsub, nbndsub, nrr)
    !! delta g matrix in Wannier basis
    !
    ! Work variables
    !
    LOGICAL :: inside_froz_etq(nbnd, nks)
    !! True if etq(ibnd, ik) is inside the frozen window
    LOGICAL :: inside_froz_etk(nbnd, nks)
    !! True if et(ibnd, ik) is inside the frozen window
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: kbnd
    !! Counter on bands
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xk_cry(3, nks)
    !! k point coordniates in crystal coordinate
    REAL(KIND = DP) :: xq_cry(3)
    !! q point coordniates in crystal coordinate
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(KIND = DP) :: qmat(nbnd, nbnd)
    !! I - U(k+q) * U(k+q)^\dagger
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: eptmp(nbndsub, nbnd)
    !!  e-p matrix, temporary
    COMPLEX(KIND = DP) :: eptmp1(nbnd, nbnd)
    !!  e-p matrix, temporary
    COMPLEX(KIND = DP) :: eptmp2(nbnd, nbnd)
    !!  e-p matrix, temporary
    COMPLEX(KIND = DP), ALLOCATABLE :: dgmats(:, :, :)
    !!  delta g matrix in smooth Bloch basis, coarse mesh
    !
    !----------------------------------------------------------
    !  STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    ! dg_pq(k, q) = U(k+q)^dagger_pm * E_{m,k+q} * Q(k+q)_mp
    !             * g(k,q)_pn / (E_{n,k} - E_{p,k+q}) * U(k)_nq
    !
    ! Sum over n is done only for bands inside the frozen window
    ! Sum over m and p can be done only for bands outside the frozen window
    !
    ! g(k,q)      is epmatk (ibnd, jbnd, ik)
    ! dg_pq(k,q)  is dgmats (ibnd, jbnd, ik)
    !
    CALL start_clock ('dg: step 1')
    !
    ALLOCATE(dgmats(nbndsub, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('dgbloch2wane', 'Error allocating dgmats', 1)
    !
    dgmats = czero
    !
    ! .true. if energy inside the frozen window
    inside_froz_etq = (etq > ahc_win_min) .AND. (etq < ahc_win_max)
    inside_froz_etk = (etk > ahc_win_min) .AND. (etk < ahc_win_max)
    !
    DO ik = 1, nks
      !
      ! qmat = I - U(k+q) * U(k+q)^\dagger
      !
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbndsub, &
        cone, cuq(:, :, ik), nbnd, cuq(:, :, ik), nbnd, czero, qmat, nbnd)
      !
      qmat = -qmat
      DO ibnd = 1, nbnd
        qmat(ibnd, ibnd) = qmat(ibnd, ibnd) + cone
      ENDDO
      !
      ! Compute Q(k+q)_jk * g(k,q)_ki / (E_{i,k} - E_{k,k+q})
      !
      ! i : inside frozen window
      ! j, k: outside frozen window
      !
      eptmp1 = czero
      !
      ! Compute and use eptmp1(k, i) = g(k,q)_ki / (E_{i,k} - E_{k,k+q})
      !
      DO ibnd = 1, nbnd
        IF (.NOT. inside_froz_etk(ibnd, ik)) CYCLE
        !
        DO kbnd = 1, nbnd
          IF (inside_froz_etq(kbnd, ik)) CYCLE
          !
          eptmp1(kbnd, ibnd) = epmatk(kbnd, ibnd, ik) / (etk(ibnd, ik) - etq(kbnd, ik))
          !
        ENDDO ! kbnd
      ENDDO ! ibnd
      !
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, cone, qmat, nbnd, &
          eptmp1, nbnd, czero, eptmp2, nbnd)
      !
      DO ibnd = 1, nbnd
        eptmp2(ibnd, :) = eptmp2(ibnd, :) * etq(ibnd, ik)
      ENDDO
      !
      ! the two zgemm calls perform the following ops:
      ! dgmats  = [ cu(ikq)^\dagger * eptmp2 ] * cu(ikk)
      ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub]
      !
      CALL ZGEMM('c', 'n', nbndsub, nbnd, nbnd, cone, cuq(:, :, ik), &
                nbnd, eptmp2, nbnd, czero, eptmp, nbndsub)
      CALL ZGEMM('n', 'n', nbndsub, nbndsub, nbnd, cone, eptmp, &
                nbndsub, cu(:, :, ik), nbnd, czero, dgmats(:, :, ik), nbndsub)
      !
    ENDDO
    !
    CALL stop_clock ('dg: step 1')
    !
    !----------------------------------------------------------------------
    !  STEP 2: Fourier transform to obtain matrix elements in wannier basis
    !----------------------------------------------------------------------
    !
    ! dg_pq(R_e, q) = (1/nkc) sum_k e^{-ikR_e} dg_pq(k, q)
    !               + (1/nkc) sum_k e^{-i(k+q)R_e} [dg_pq(k, q)]^*
    ! dg_pq(R_e, q) is dgmatw (nbndsub, nbndsub, ir)
    ! dg_pq(k, q)   is dgmats (nbndsub, nbndsub, nks)
    !
    CALL start_clock('dg: step 2')
    !
    dgmatw(:, :, :) = czero
    !
    ! bring xk and xq in crystal coordinates
    !
    xk_cry = xk
    xq_cry = xq
    CALL cryst_to_cart(nks, xk_cry, at, -1)
    CALL cryst_to_cart(1, xq_cry, at, -1)
    !
    DO ir = 1, nrr
      DO ik = 1, nks
        !
        ! First term (1/nkc) sum_k e^{-ikR_e} dg_pq(k, q)
        rdotk = twopi * DOT_PRODUCT(xk_cry(:, ik), REAL(irvec(:, ir), DP))
        cfac = EXP(-ci * rdotk) / DBLE(nkstot)
        dgmatw(:, :, ir) = dgmatw(:, :, ir) + cfac * dgmats(:, :, ik)
        !
        ! Second term (1/nkc) sum_k e^{-i(k+q)R_e} [dg_pq(k, q)]^*
        rdotk = twopi * DOT_PRODUCT((xk_cry(:, ik) + xq_cry(:)), REAL(irvec(:, ir), DP))
        cfac = EXP(-ci * rdotk) / DBLE(nkstot)
        dgmatw(:, :, ir) = dgmatw(:, :, ir) + cfac * CONJG(dgmats(:, :, ik))
        !
      ENDDO
    ENDDO
    !
    CALL mp_sum(dgmatw, inter_pool_comm)
    !
    !  Check spatial decay of matrix elements in Wannier basis
    !  the unit in r-space is angstrom, and I am plotting
    !  the matrix for the first mode only
    !
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iuwane, FILE = 'decay.dgwane')
      WRITE(iuwane, '(a)') '# Spatial decay of delta g matrix elements in Wannier basis'
      WRITE(iuwane, '(a)') '# Contains only the Sternheimer contribution, not the epmat part.'
      DO ir = 1, nrr
        !
        tmp = MAXVAL(ABS(dgmatw(:, :, ir)))
        WRITE(iuwane, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iuwane)
    ENDIF
    !
    DEALLOCATE(dgmats, STAT = ierr)
    IF (ierr /= 0) CALL errore('dgbloch2wane', 'Error deallocating dgmats', 1)
    !
    CALL stop_clock ('dg: step 2')
    !
    ! -------------------------------------------------------------------------
    END SUBROUTINE dgbloch2wane
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE sthbloch2wane(nbnd, nbndsub, nks, nkstot, etk, etq, &
      ahc_win_min, ahc_win_max, xk, cu, cuq, epmatk1, epmatk2, sthmatk, &
      nrr, irvec, wslen, sthmatw)
    !--------------------------------------------------------------------------
    !!
    !! From the Sternheimer matrix elements in Bloch representation
    !! (coarse mesh), find the Sternheimer matrix in Wannier representation
    !!
    !! s_pq(k, q) = sum_{m in F, n}U(k+q)^dagger_pm * s_mn(k, q) * U(k)_mq
    !! s_mn(k, q) = sthmatk_mn(k,q)
    !!            + sum_{a,b} [g1_am(k,q)]^* * Q_ab(k+q) * g2_bn(k,q)
    !!              / (E_{m,k} - E_{a,k+q})
    !!
    !! (This equation is Eq.(33) of [1], with w^{(1)} as defined in Eq.(8).)
    !!
    !! Sum over m is done only for bands inside the frozen window
    !! Sum over a and b is done only for bands outside the frozen window
    !!
    !! Only bands inside the outer window are considered. Bands outside the
    !! outer window need not be included. They are already removed from the
    !! input matrices.
    !!
    !! s_pq(R_e, q) = (1/nkc) sum_k e^{-ikR_e} s_pq(k, q)
    !!
    !! Q(k+q) is a projection out of the Wannier subspace:
    !! Q(k+q)_mn = delta_mn - sum_p U(k+q)_np U(k+q)^dagger_pm
    !!
    !! References
    !! [1] Lihm and Park, PRX 11, 041053 (2021)
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : bohr2ang, twopi, ci, czero, cone, rytoev
    USE mp,            ONLY : mp_sum
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : inter_pool_comm
    USE io_var,        ONLY : iuwane
    USE io_global,     ONLY : ionode_id
    USE cell_base,     ONLY : at, alat
    !
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands (nbndep)
    INTEGER, INTENT(in) :: nbndsub
    !! Number of bands in the optimal subspace
    INTEGER, INTENT(in) :: nks
    !! Number of kpoints in this pool
    INTEGER, INTENT(in) :: nrr
    !! Number of WS points
    INTEGER, INTENT(in) :: nkstot
    !! Total number of kpoints
    INTEGER, INTENT(in) :: irvec(3, nrr)
    !! Coordinates of WS points
    !
    REAL(KIND = DP), INTENT(in) :: etk(nbnd, nks)
    !! Hamiltonian eigenvalues at k within the outer window
    REAL(KIND = DP), INTENT(in) :: etq(nbnd, nks)
    !! Hamiltonian eigenvalues at k+q within the outer window
    REAL(KIND = DP), INTENT(in) :: ahc_win_min
    !! AHC window of Wannierization (in eV)
    REAL(KIND = DP), INTENT(in) :: ahc_win_max
    !! AHC window of Wannierization (in eV)
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(in) :: cuq(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(in) :: epmatk1(nbnd, nbnd, nks)
    !! e-p matrix in bloch representation, coarse mesh
    COMPLEX(KIND = DP), INTENT(in) :: epmatk2(nbnd, nbnd, nks)
    !! e-p matrix in bloch representation, coarse mesh
    COMPLEX(KIND = DP), INTENT(in) :: sthmatk(nbnd, nbnd, nks)
    !! Sternheimer matrix in bloch representation, coarse mesh
    !
    ! output variables
    !
    COMPLEX(KIND = DP), INTENT(out) :: sthmatw(nbndsub, nbndsub, nrr)
    !! Sternheimer matrix in Wannier basis
    !
    ! Work variables
    !
    LOGICAL :: inside_froz_etq(nbnd, nks)
    !! True if etq(ibnd, ik) is inside the frozen window
    LOGICAL :: inside_froz_etk(nbnd, nks)
    !! True if et(ibnd, ik) is inside the frozen window
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xk_cry(3, nks)
    !! k point coordniates in crystal coordinate
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(KIND = DP) :: qmat(nbnd, nbnd)
    !! I - U(k+q) * U(k+q)^\dagger
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: eptmp(nbndsub, nbnd)
    !! e-p matrix, temporary
    COMPLEX(KIND = DP) :: eptmp1(nbnd, nbnd)
    !! e-p matrix, temporary
    COMPLEX(KIND = DP) :: eptmp2(nbnd, nbnd)
    !! e-p matrix, temporary
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmats(:, :, :)
    !! Sternheimer matrix in smooth Bloch basis, coarse mesh
    !
    !----------------------------------------------------------
    !  STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    !  s_pq(k, q) = sum_{m in F, n} U(k+q)^dagger_pm *
    !  [ sthmatk_mn(k,q)
    !   + sum_{a,b} [g1_am(k,q)]^* * Q_ab(k+q) * g2_bn(k,q) / (E_{m,k} - E_{a,k+q})
    !  ] * U(k)_mq
    !
    ! Sum over n is done only for bands inside the frozen window
    ! Sum over a and b can be done only for bands outside the frozen window
    !
    ! g1(k,q)      is epmatk1 (ibnd, jbnd, ik)
    ! g2(k,q)      is epmatk2 (ibnd, jbnd, ik)
    ! sthmatk(k,q) is sthmatk (ibnd, jbnd, ik)
    ! s_pq(k,q)    is sthmats (ibnd, jbnd, ik)
    !
    CALL start_clock ('sth: step 1')
    !
    ALLOCATE(sthmats(nbndsub, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('sthbloch2wane', 'Error allocating sthmats', 1)
    !
    sthmats = czero
    !
    ! True if energy inside the frozen window
    inside_froz_etq = (etq > ahc_win_min) .AND. (etq < ahc_win_max)
    inside_froz_etk = (etk > ahc_win_min) .AND. (etk < ahc_win_max)
    !
    DO ik = 1, nks
      !
      ! qmat = I - U(k+q) * U(k+q)^\dagger
      !
      CALL ZGEMM('n', 'c', nbnd, nbnd, nbndsub, cone, cuq(:, :, ik), &
                nbnd, cuq(:, :, ik), nbnd, czero, qmat, nbnd)
      !
      qmat = -qmat
      DO ibnd = 1, nbnd
        qmat(ibnd, ibnd) = qmat(ibnd, ibnd) + cone
      ENDDO
      !
      eptmp1 = czero
      !
      ! Compute and use eptmp1(i, j) = [g1_ji]^* / (E_{i,k} - E_{j,k+q})
      !
      DO ibnd = 1, nbnd
        IF (.NOT. inside_froz_etk(ibnd, ik)) CYCLE
        !
        DO jbnd = 1, nbnd
          IF (inside_froz_etq(jbnd, ik)) CYCLE
          !
          eptmp1(ibnd, jbnd) = CONJG(epmatk1(jbnd, ibnd, ik)) / (etk(ibnd, ik) - etq(jbnd, ik))
          !
        ENDDO ! jbnd
      ENDDO ! ibnd
      !
      ! Compute eptmp1_ma * Q_ab(k+q) * g2_bn(k,q)
      !
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, &
        cone, eptmp1, nbnd, qmat, nbnd, czero, eptmp2, nbnd)
      CALL ZGEMM('n', 'n', nbnd, nbnd, nbnd, &
        cone, eptmp2, nbnd, epmatk2(:, :, ik), nbnd, czero, eptmp1, nbnd)
      !
      ! Add sthmatk part to eptmp1
      !
      eptmp1 = eptmp1 + sthmatk(:, :, ik)
      !
      ! the two zgemm calls perform the following ops:
      ! sthmats  = [ cu(ikq)^\dagger * eptmp1 ] * cu(ikk)
      ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub]
      !
      CALL ZGEMM('c', 'n', nbndsub, nbnd, nbnd, cone, cu(:, :, ik), nbnd, &
          eptmp1, nbnd, czero, eptmp, nbndsub)
      CALL ZGEMM('n', 'n', nbndsub, nbndsub, nbnd, cone, eptmp, nbndsub, &
          cu(:, :, ik), nbnd, czero, sthmats(:, :, ik), nbndsub)
      !
    ENDDO
    !
    CALL stop_clock ('sth: step 1')
    !
    !----------------------------------------------------------------------
    !  STEP 2: Fourier transform to obtain matrix elements in wannier basis
    !----------------------------------------------------------------------
    !
    ! sth_pq(R_e, q) = (1/nkc) sum_k e^{-ikR_e} sth_pq(k, q)
    !
    ! sth_pq(R_e, q) is sthmatw(1:nbndsub, 1:nbndsub, ir)
    ! sth_pq(k, q)   is sthmats(1:nbndsub, 1:nbndsub, ik)
    !
    CALL start_clock('sth: step 2')
    !
    sthmatw(:, :, :) = czero
    !
    ! bring xk and xq in crystal coordinates
    !
    xk_cry = xk
    CALL cryst_to_cart(nks, xk_cry, at, -1)
    !
    DO ir = 1, nrr
      DO ik = 1, nks
        !
        rdotk = twopi * DOT_PRODUCT(xk_cry(:, ik), REAL(irvec(:, ir), DP))
        cfac = EXP(-ci * rdotk) / DBLE(nkstot)
        !
        sthmatw(:, :, ir) = sthmatw(:, :, ir) + cfac * sthmats(:, :, ik)
        !
      ENDDO
    ENDDO
    !
    CALL mp_sum(sthmatw, inter_pool_comm)
    !
    !  Check spatial decay of matrix elements in Wannier basis
    !  the unit in r-space is angstrom, and I am plotting
    !  the matrix for the first mode only
    !
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iuwane, FILE = 'decay.sthwane')
      WRITE(iuwane, '(a)') '# Spatial decay of Sternheimer matrix elements in Wannier basis'
      DO ir = 1, nrr
        !
        tmp = MAXVAL(ABS(sthmatw(:, :, ir)))
        WRITE(iuwane, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iuwane)
    ENDIF
    !
    DEALLOCATE(sthmats, STAT = ierr)
    IF (ierr /= 0) CALL errore('sthbloch2wane', 'Error deallocating sthmats', 1)
    !
    CALL stop_clock ('sth: step 2')
    !
    ! -------------------------------------------------------------------------
    END SUBROUTINE sthbloch2wane
    !--------------------------------------------------------------------------
    !
  !----------------------------------------------------------------------------
  END MODULE bloch2wannier
  !----------------------------------------------------------------------------
