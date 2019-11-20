  !
  !k Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE bloch2wan
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
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, zero, ryd2ev
    USE io_global, ONLY : ionode_id
    USE io_var,    ONLY : iudecayH
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: lwin(nbnd, nks)
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
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    COMPLEX(KIND = DP), INTENT(out) :: chw(nbndsub, nbndsub, nrr)
    !! Hamiltonian in Wannier basis
    !
    ! Local variables
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
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
    !
    REAL(KIND = DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}
    REAL(KIND = DP) :: tmp
    !! Maximum value of the real space Hamiltonian
    REAL(KIND = DP) :: et_opt(nbnd, nks)
    !! hamiltonian eigenvalues within the outer window in the first ndimwin(ik) entries
    REAL(KIND = DP) :: et_tmp(nbnd, nks)
    !! temporary array for hamiltonian eigenvalues 
    !
    COMPLEX(KIND = DP) :: chs(nbndsub, nbndsub, nks)
    !! Hamiltonian in Bloch basis, coarse k-mesh
    COMPLEX(KIND = DP) :: cfac
    !! $$e^{-i\mathbf{r}\cdot\mathbf{k}}$$
    COMPLEX(KIND = DP) :: ctmp
    !! Temporary variable to store the Hamiltonian
    !
    !----------------------------------------------------------
    !    STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    CALL start_clock ('Ham: step 1')
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
    et_tmp = zero
    et_opt = zero
    IF (nexband_tmp > 0) THEN 
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbnd
          IF (exband(i)) CYCLE
          IF (lwin(i, ik)) THEN
            ibnd = ibnd + 1
            et_tmp(ibnd, ik) = et(i, ik)
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbnd
          IF (exband(i)) CYCLE
          ibnd = ibnd + 1
          et_opt(i, ik) = et_tmp(ibnd, ik)
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbnd
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
          DO mbnd = 1, nbnd
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
      WRITE(iudecayH, '(/3x,a/)') '#Spatial decay of Hamiltonian in Wannier basis'
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
    CALL stop_clock('Ham: step 2')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE hambloch2wan 
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dmebloch2wan(nbnd, nbndsub, nks, nkstot, dmec, xk, cu, &
                            nrr, irvec, wslen, lwin, exband)
    !--------------------------------------------------------------------------
    !!
    !!  From the Dipole in Bloch representationi (coarse mesh), 
    !!  find the corresponding Dipole in Wannier representation 
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE elph2,         ONLY : cdmew
    USE io_var,        ONLY : iudecayP
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
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
    COMPLEX(KIND = DP), INTENT(in) :: dmec(3, nbnd, nbnd, nks)
    !! Dipole matrix elements on coarse mesh
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    !
    LOGICAL, INTENT(in) :: lwin(nbnd, nks)
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
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    !
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(KIND = DP) :: cps(3, nbndsub, nbndsub, nks)
    !! Dipole in smooth Bloch basis, coarse mesh
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: dmec_utmp(nbnd, nbndsub)
    !! dmec after multiplication with the Wannier rotation matrix cu.
    COMPLEX(KIND = DP) :: dmec_opt(3, nbnd, nbnd, nks)
    !! dmec computed in pmn, rescaled down if skipping lwin.
    COMPLEX(KIND = DP) :: dmec_tmp(3, nbnd, nbnd, nks)
    !! temporary dmec matrix
    !
    CALL start_clock ('Dipole: step 1')
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
    dmec_tmp = czero
    IF (nexband_tmp > 0) THEN
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          IF (lwin(j,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              IF (lwin(i,ik)) THEN
                ibnd = ibnd + 1
                dmec_tmp(:, ibnd, jbnd, ik) = dmec(:, i, j, ik)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          jbnd = jbnd + 1
          ibnd = 0
          DO i = 1, nbnd
            IF (exband(i)) CYCLE
            ibnd = ibnd + 1
            dmec_opt(:, i, j, ik) = dmec_tmp(:, ibnd, jbnd, ik)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (lwin(j, ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
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
        CALL ZGEMM('n', 'n', nbnd, nbndsub, nbnd, cone, dmec_opt(ipol, :, :, ik), &
                   nbnd, cu(:, :, ik), nbnd, czero, dmec_utmp(:, :), nbnd)
        CALL ZGEMM('c', 'n', nbndsub, nbndsub, nbnd, cone, cu(:, :, ik), &
                   nbnd, dmec_utmp(:, :), nbnd, czero, cps(ipol, :, :, ik), nbndsub)
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
      WRITE(iudecayP, '(/3x,a/)') '#Spatial decay of dipole in Wannier basis'
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
    USE elph2,         ONLY : rdw, epsi, zstar
    USE epwcom,        ONLY : lpolar, nqc1, nqc2, nqc3
    USE io_var,        ONLY : iudecaydyn
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_global,     ONLY : ionode_id
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier
    USE mp_global,     ONLY : inter_pool_comm
    USE rigid_epw,     ONLY : rgd_blk
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
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    !
    !  subtract the long-range term from D(q)
    IF (lpolar) THEN
      DO ik = 1, nq
        !xk has to be in cart. coord.
        CALL rgd_blk(nqc1, nqc2, nqc3, nat, dynq(:, :, ik), xk(:, ik), tau, epsi, zstar, -1.d0)
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
    rdw( :, :, :) = 0.d0
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
      WRITE(iudecaydyn, '(/3x,a/)') '#Spatial decay of Dynamical matrix in Wannier basis'
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
    !--------------------------------------------------------------------------
    SUBROUTINE vmebloch2wan(nbnd, nbndsub, nks, nkstot, xk, cu, &
                            nrr, irvec, wslen, lwin, exband)
    !--------------------------------------------------------------------------
    !!
    !!  Calculate the velocity matrix elements in the Wannier basis
    !!  at no point do we actually have the coarse mesh v-ME. 
    !!
    !! RM 03/2018: debugged and updated
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, alat
    USE elph2,     ONLY : cvmew
    USE constants_epw, ONLY : twopi, one, zero, ci, czero, cone, bohr2ang
    USE io_var,    ONLY : iummn, iubvec, iudecayv
    USE io_files,  ONLY : prefix
    USE io_global, ONLY : ionode_id, stdout
    USE mp_global, ONLY : inter_pool_comm, world_comm
    USE mp,        ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
    USE kfold,     ONLY : ktokpmq  
    !
    IMPLICIT NONE
    !
    ! Input variables
    LOGICAL, INTENT(in) :: lwin(nbnd,nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
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
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(KIND = DP), INTENT(in) :: wslen(nrr)
    !! WS vectors length (alat units)
    COMPLEX(KIND = DP), INTENT(in) :: cu(nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
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
    INTEGER :: jbnd
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
    COMPLEX(KIND = DP) :: cu_big(nbnd, nbndsub, nkstot)
    !! rotation matrix from wannier code
    REAL(KIND = DP) ::  b_tmp(3)
    !! temporary b-vectors
    REAL(KIND = DP) :: zero_vect(3)
    !! temporary zero vector
    REAL(KIND = DP) :: delta
    !! \delta_nm = 1 if n == m and 0 if n /= m
    COMPLEX(KIND = DP) :: Apos(3, nbndsub, nbndsub, nks)
    !! A^W_{mn,\alpha}(k)
    COMPLEX(KIND = DP), ALLOCATABLE :: M_mn(:, :, :, :)
    !! M_mn(k,b)
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mat_opt(:, :, :, :)
    !! M_mn(k,b) computed in mmn, rescaled down if skipping lwin.
    COMPLEX(KIND = DP), ALLOCATABLE :: m_mat_tmp(:, :, :, :)
    !! temporary M_mn matrix
    COMPLEX(KIND = DP), ALLOCATABLE :: cvs(:, :, :, :)
    !! M_mn in smooth Bloch basis, coarse k-mesh
    COMPLEX(KIND = DP) :: M_mn_utmp(nbnd, nbndsub)
    !! M_mn after multiplication with the Wannier rotation matrix cu.
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
        ! if it doesn't exist, then we just set the bvec and wb to zero
        !
        nnb = 1
        ALLOCATE(bvec(3, nnb, nkstot), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating bvec', 1)
        ALLOCATE(wb(nnb), STAT = ierr)
        IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating wb', 1)
        bvec = zero
        wb   = zero
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
    ! RM - M_mn matrix is writen on file in pw2wan90epw.f90/compute_mmn_para
    !    - dimensions of M_mn are M_mn(nbnd, nbnd, nnb, nkstot)
    !
    ALLOCATE(M_mn(nbnd, nbnd, nnb, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating M_mn', 1)
    M_mn = czero
    !
    IF (mpime == ionode_id) THEN
      tempfile = TRIM(prefix)//'.mmn'
      OPEN(iummn, FILE = tempfile, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
      !
      IF (ios /= 0) THEN
        ! if it doesn't exist, then we just set the mmn to zero
        CALL errore ('vmebloch2wan','error opening' // tempfile, 0)
      ELSE
        !
        DO ik = 1, nkstot
          DO ib = 1, nnb
            DO jbnd = 1, nbnd
              DO ibnd = 1, nbnd
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
    ALLOCATE(m_mat_opt(nbnd, nbnd, nnb, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating m_mat_opt', 1)
    ALLOCATE(m_mat_tmp(nbnd, nbnd, nnb, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error allocating m_mat_tmp', 1)
    m_mat_opt(:, :, :, :) = czero
    m_mat_tmp(:, :, :, :) = czero
    zero_vect(:) = zero
    !
    IF (nexband_tmp > 0) THEN
      DO ik = 1, nks
        CALL ktokpmq(xk(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          IF (lwin(j,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              IF (lwin(i,ik)) THEN
                ibnd = ibnd + 1
                m_mat_tmp(ibnd, jbnd, :, ik) = M_mn(i, j, :, nkk_abs)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          jbnd = jbnd + 1
          ibnd = 0
          DO i = 1, nbnd
            IF (exband(i)) CYCLE
            ibnd = ibnd + 1
            m_mat_opt(i, j, :, ik) = m_mat_tmp(ibnd, jbnd, :, ik)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        CALL ktokpmq(xk(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        jbnd = 0
        DO j = 1, nbnd
          IF (lwin(j, ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
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
    DEALLOCATE(m_mat_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating m_mat_tmp', 1)
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
      CALL ktokpmq(xk(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
      !
      DO ib = 1, nnb
        !
        ! bring bvec to units of 2piba since xk is cartesian units of 2piba
        b_tmp(:) = alat / (twopi) * bvec(:, ib, nkk_abs)
        CALL ktokpmq(xk(:, ik), b_tmp(:), +1, ipool, nkb, nkb_abs)
        !
        ! M_mn_utmp(:, :) = MATMUL( m_mat_opt(:,:,ib,ik), cu_big(:,:,nkb_abs) )
        ! cvs(:,:,ib,ik) = MATMUL( CONJG(transpose(cu(:,:,ik))), M_mn_utmp(:, :) )
        !
        CALL ZGEMM('n', 'n', nbnd, nbndsub, nbnd, cone, m_mat_opt(:, :, ib, ik), &
                   nbnd, cu_big(:, :, nkb_abs), nbnd, czero, M_mn_utmp(:, :), nbnd)
        CALL ZGEMM('c', 'n', nbndsub, nbndsub, nbnd, cone, cu(:, :, ik), &
                   nbnd, M_mn_utmp(:, :), nbnd, czero, cvs(:, :, ib, ik), nbndsub)
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
      CALL ktokpmq(xk(:, ik), zero_vect, +1, ipool, nkk, nkk_abs)
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
    DEALLOCATE(bvec, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating bvec', 1)
    DEALLOCATE(wb, STAT = ierr)
    IF (ierr /= 0) CALL errore('vmebloch2wan', 'Error deallocating wb', 1)
    !
    CALL stop_clock('Velocity: step 1')
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
    CALL cryst_to_cart(nks, xk, at, -1)
    !
    ! r_{\alpha}(R) is cvmew(3,nbndsub,nbndsub,nrr)
    !
    cvmew(:, :, :, :) = czero
    ! 
    DO ir = 1, nrr
      DO ik = 1, nks
        !
        rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
        cfac = EXP(-ci * rdotk) / DBLE(nkstot)
        cvmew(:, :, :, ir) = cvmew(:, :, :, ir) + cfac * Apos(:, :, :, ik)
        !
      ENDDO
    ENDDO
    !
    CALL mp_sum(cvmew, inter_pool_comm) 
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nks, xk, bg, 1)
    !
    ! check spatial decay of position matrix elements in Wannier basis
    !
    ! position matrix cvmew and spatial dimensions are in units of bohr
    ! [mind when comparing with wannier code (angstrom units) with write_rmn=.TRUE.]
    !
    IF (mpime == ionode_id) then
      OPEN(UNIT = iudecayv, FILE = 'decay.v')
      WRITE(iudecayv, '(/3x,a/)') '#Spatial decay of Velocity matrix element in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  MAXVAL(ABS(cvmew(:, :, :, ir)))
        WRITE(iudecayv, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iudecayv)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    CALL stop_clock('Velocity: step 2')
    !
    WRITE(stdout,'(/5x,a)') 'Velocity matrix elements calculated'
    WRITE(stdout,'(a)') ' '
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE vmebloch2wan
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE ephbloch2wane(nbnd, nbndsub, nks, nkstot, xk, &
         cu, cuq, epmatk, nrr, irvec, wslen, epmatw)
    !-----------------------------------------------------------------------
    !!
    !!  From the EP matrix elements in Bloch representation (coarse 
    !!  mesh), find the corresponding matrix elements in electron-Wannier 
    !!  representation and phonon-Bloch representation
    !!
    !-----------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, alat
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
    USE io_var,    ONLY : iuwane
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp       , ONLY : mp_sum 
    USE mp_world,  ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    !  input variables
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of bands
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
    REAL(KIND = DP), INTENT(in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
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
    COMPLEX(KIND = DP), INTENT(out) :: epmatw(nbndsub, nbndsub, nrr)
    ! EP vertex (Wannier el and Bloch ph)
    !
    ! Work variables
    !
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP) :: epmats(nbndsub, nbndsub, nks)
    !!  e-p matrix  in smooth Bloch basis, coarse mesh
    COMPLEX(KIND = DP) :: eptmp(nbndsub, nbnd)
    !!  e-p matrix, temporary
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
    DO ik = 1, nks
      !
      ! the two zgemm calls perform the following ops:
      ! epmats  = [ cu(ikq)^\dagger * epmatk ] * cu(ikk)
      ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub] 
      !
      CALL ZGEMM('c', 'n', nbndsub, nbnd, nbnd, cone, cuq(:, :, ik),  &
                nbnd, epmatk(:, :, ik), nbnd, czero, eptmp, nbndsub)
      CALL ZGEMM('n', 'n', nbndsub, nbndsub, nbnd, cone, eptmp,     &
                nbndsub, cu(:, :, ik), nbnd, czero, epmats(:, :, ik), nbndsub)
      !
    ENDDO
    !
    CALL stop_clock ('ep: step 1')
    !
    !----------------------------------------------------------------------
    !  STEP 2: Fourier transform to obtain matrix elements in electron wannier basis
    !----------------------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,q) = (1/nkc) sum_k e^{-ikR_e} g~(k,q)
    ! g(R_e,q) is epmatw (nbndsub,nbndsub,ir)
    !
    CALL start_clock('ep: step 2')
    !
    epmatw(:, :, :) = czero
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart(nks, xk, at, -1)
    !
    DO ir = 1, nrr
      DO ik = 1, nks
        !
        rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
        cfac = EXP(-ci * rdotk) / DBLE(nkstot)
        epmatw( :, :, ir) = epmatw( :, :, ir) + cfac * epmats(:, :, ik)
        !
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
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iuwane, FILE = 'decay.epwane')
      WRITE(iuwane, '(a)') '# Spatial decay of e-p matrix elements in electron Wannier basis'
      DO ir = 1, nrr
        ! 
        tmp = MAXVAL(ABS(epmatw(:, :, ir))) 
        WRITE(iuwane, *) wslen(ir) * alat * bohr2ang, tmp
        !
      ENDDO
      !
      CLOSE(iuwane)
    ENDIF
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
    USE elph2,         ONLY : epmatwp
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_var,        ONLY : iuwanep
    USE io_global,     ONLY : ionode_id
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
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
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    REAL(KIND = DP) :: rvec1(3)
    !!
    REAL(KIND = DP) :: rvec2(3)
    !!
    REAL(KIND = DP) :: len1
    !!
    REAL(KIND = DP) :: len2
    !!
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
      !  check spatial decay of EP matrix elements in wannier basis - electrons + phonons
      !
      !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
      !
      IF (mpime == ionode_id) THEN
        IF (ir == 1) OPEN(UNIT = iuwanep, FILE = 'decay.epmat_wanep', STATUS = 'unknown')
        IF (ir == 1) WRITE(iuwanep, '(a)') '#  R_e,    R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)| '
        DO ire = 1, nrr_k
          !
          rvec1 = DBLE(irvec_k(1, ire)) * at(:, 1) + &
                  DBLE(irvec_k(2, ire)) * at(:, 2) + &
                  DBLE(irvec_k(3, ire)) * at(:, 3)
          rvec2 = DBLE(irvec_g(1, ir)) * at(: ,1) + &
                  DBLE(irvec_g(2, ir)) * at(: ,2) + &
                  DBLE(irvec_g(3, ir)) * at(: ,3)
          len1 = DSQRT(rvec1(1)**2.d0 + rvec1(2)**2.d0 + rvec1(3)**2.d0)
          len2 = DSQRT(rvec2(1)**2.d0 + rvec2(2)**2.d0 + rvec2(3)**2.d0)
          tmp =  MAXVAL(ABS(epmatwp(:, :, ire, :, ir)))
          !
          ! rvec1 : electron-electron0 distance
          ! rvec2 : phonon - electron0 distance
          !
          WRITE(iuwanep, '(5f15.10)') len1 * alat * bohr2ang, len2 * alat * bohr2ang, tmp
        ENDDO
        IF (ir == nrr_g) CLOSE(iuwanep)
      ENDIF
      !
    ENDDO
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
    SUBROUTINE ephbloch2wanp_mem(nbnd, nmodes, xk, nq, irvec_k, irvec_g, nrr_k, nrr_g, epmatwe)
    !--------------------------------------------------------------------------
    !
    !!  From the EP matrix in electron-Wannier representation and 
    !!  phonon-Bloch representation (coarse mesh), find the corresponding matrix 
    !!  electron-Wannier representation and phonon-Wannier representation 
    !
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, alat
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_var,        ONLY : iunepmatwe, iunepmatwp, iuwanep
    USE io_global,     ONLY : ionode_id
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_epw, ONLY : rwepmatw
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
    COMPLEX(KIND = DP), INTENT(inout) :: epmatwe(nbnd, nbnd, nrr_k, nmodes)
    !! EP matrix in electron-Wannier representation and phonon-Bloch representation
    !!   (Cartesian coordinates)
    !
    ! work variables
    !
    INTEGER :: iq
    !! Counter on q-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ire
    !! Counter on WS points
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(KIND = DP) :: tmp
    !! Temporary variables
    REAL(KIND = DP) :: rvec1(3)
    !!
    REAL(KIND = DP) :: rvec2(3)
    !!
    REAL(KIND = DP) :: len1
    !!
    REAL(KIND = DP) :: len2
    !!
    COMPLEX(KIND = DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwp_mem(:, :, :, :)
    !!  e-p matrix in Wannier basis
    !
    ALLOCATE(epmatwp_mem(nbnd, nbnd, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'Error allocating epmatwp_mem', 1)
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
    DO ir = 1, nrr_g
      !
      epmatwp_mem = czero
      ! 
      DO iq = 1, nq
        !
        ! direct read of epmatwe for this iq 
        CALL rwepmatw(epmatwe, nbnd, nrr_k, nmodes, iq, iunepmatwe, -1)
        !
        rdotk = twopi * DOT_PRODUCT(xk(:, iq), DBLE(irvec_g(:, ir)))
        cfac = EXP(-ci * rdotk) / DBLE(nq)
        epmatwp_mem = epmatwp_mem + cfac * epmatwe
        !
      ENDDO
      !
      ! direct write of epmatwp_mem for this ir 
      CALL rwepmatw(epmatwp_mem, nbnd, nrr_k, nmodes, ir, iunepmatwp, +1)
      !  check spatial decay of EP matrix elements in wannier basis - electrons + phonons
      !
      !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
      !
      IF (mpime == ionode_id) THEN
        IF (ir == 1) OPEN(UNIT = iuwanep, FILE = 'decay.epmat_wanep', STATUS = 'unknown')
        IF (ir == 1) WRITE(iuwanep, '(a)') '#  R_e,    R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)| '
        DO ire = 1, nrr_k
          !
          rvec1 = DBLE(irvec_k(1, ire)) * at(:, 1) + &
                  DBLE(irvec_k(2, ire)) * at(:, 2) + &
                  DBLE(irvec_k(3, ire)) * at(:, 3)
          rvec2 = DBLE(irvec_g(1, ir)) * at(:, 1) + &
                  DBLE(irvec_g(2, ir)) * at(:, 2) + &
                  DBLE(irvec_g(3, ir)) * at(:, 3)
          len1 = DSQRT(rvec1(1)**2.d0 + rvec1(2)**2.d0 + rvec1(3)**2.d0)
          len2 = DSQRT(rvec2(1)**2.d0 + rvec2(2)**2.d0 + rvec2(3)**2.d0)
          tmp =  MAXVAL(ABS(epmatwp_mem(:, :, ire, :)))
          !
          ! rvec1 : electron-electron0 distance
          ! rvec2 : phonon - electron0 distance
          !
          WRITE(iuwanep, '(5f15.10)') len1 * alat * bohr2ang, len2 * alat * bohr2ang, tmp
        ENDDO
        IF (ir == nrr_g) CLOSE(iuwanep)
      ENDIF
      !
    ENDDO
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart(nq, xk, bg, 1)
    !
    DEALLOCATE(epmatwp_mem, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephbloch2wanp_mem', 'Error deallocating epmatwp_mem', 1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE ephbloch2wanp_mem
    !--------------------------------------------------------------------------
    ! 
  !----------------------------------------------------------------------------
  END MODULE bloch2wan
  !----------------------------------------------------------------------------
