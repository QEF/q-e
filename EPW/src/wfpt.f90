  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !------------------------------------------------------------------------------
  MODULE wfpt
  !------------------------------------------------------------------------------
  !!
  !! This module contains variables and subroutines of polaron
  !! Authored by Jae-Mo Lihm
  !!  [1] J.-M. Lihm and C.-H. Park, PRX 11, 041053 (2021)
  !!
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: dwmatf_trunc(:, :, :, :, :)
  !! Debye-Waller matrix in smooth Bloch basis, fine mesh, calculated with
  !! truncation
  !
  CONTAINS
    !
    !------------------------------------------------------------------------------
    SUBROUTINE wfpt_bloch2wan_setup(nrr_k)
    !------------------------------------------------------------------------------
    !!
    !! Initialize variables used for Bloch to Wannier transformation in WFPT
    !!
    !------------------------------------------------------------------------------
    !
    USE io_global,     ONLY : ionode
    USE io_files,      ONLY : prefix, diropn
    USE symm_base,     ONLY : nsym, time_reversal
    USE pwcom,         ONLY : nks, nkstot
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : czero
    USE input,         ONLY : nbndsub
    USE global_var,    ONLY : nbndep, dw_mat, dwmatwe, cpmew, sthmatq, dgmatwe, &
                              sthmatwe
    USE io_var,        ONLY : iusthwe, iudgwe
    USE symmetry,      ONLY : read_sym_ktok, sthmatq_save
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nrr_k
    !! Number of WS points for electrons
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    INTEGER :: lrsthmatw
    !! record length for sthmatwe
    INTEGER :: lrdgmatw
    !! record length for dgmatwe
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(dw_mat(nbndep, nbndep, nks, 3, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating dw_mat', 1)
    dw_mat(:, :, :, :, :) = czero
    !
    ALLOCATE(dwmatwe(nbndsub, nbndsub, nrr_k, 3, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating dwmatwe', 1)
    dwmatwe(:, :, :, :, :) = czero
    !
    ALLOCATE(cpmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating cpmew', 1)
    cpmew(:, :, :, :) = czero
    !
    ALLOCATE(sthmatq(nbndep, nbndep, nks, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating sthmatq', 1)
    sthmatq(:, :, :, :, :) = czero
    !
    ALLOCATE(sthmatq_save(nbndep, nbndep, nks, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating sthmatq_save', 1)
    sthmatq_save(:, :, :, :, :) = czero
    !
    ALLOCATE(dgmatwe(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating dgmatwe', 1)
    dgmatwe(:, :, :, :) = czero
    !
    ALLOCATE(sthmatwe(nbndsub, nbndsub, nrr_k, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_setup', 'Error allocating sthmatwe', 1)
    sthmatwe(:, :, :, :, :) = czero
    !
    CALL read_sym_ktok(time_reversal, nsym, nkstot)
    !
    ! Open the prefix.sthmatwe file
    IF (ionode) THEN
      lrsthmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint = TRIM(prefix)//'.sthmatwe'
      CALL diropn(iusthwe, 'sthmatwe', lrsthmatw, exst)
      !
      lrdgmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint = TRIM(prefix)//'.dgmatwe'
      CALL diropn(iudgwe, 'dgmatwe', lrdgmatw, exst)
    ENDIF
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE wfpt_bloch2wan_setup
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE wfpt_bloch2wan_iq(iq, nqc, xqc, nrr_k, irvec_k, wslen_k, epmatwe_iq)
    !------------------------------------------------------------------------------
    !!
    !! Wannier interpolation of electron-phonon vertex
    !!
    !! The routine is divided in two parts.
    !! 1) Full coarse grid to real-space Wannier transformation
    !! 2) Real-space Wannier to fine grid Bloch space interpolation
    !------------------------------------------------------------------------------
    !
    USE io_global,     ONLY : ionode
    USE pwcom,         ONLY : nbnd, nks, nkstot
    USE modes,         ONLY : nmodes
    USE ep_constants,  ONLY : eps8, zero
    USE input,         ONLY : et_loc, xk_loc, et_all, nbndsub, ahc_win_min, ahc_win_max
    USE global_var,    ONLY : pmec, exband, nbndep, cu, cuq, lwin, lwinq, dwmatwe, &
                              dw_mat, cpmew, epmatq, dgmatwe, sthmatwe, sthmatq
    USE io_var,        ONLY : iusthwe, iudgwe
    USE bloch2wannier, ONLY : dmebloch2wan, ephbloch2wane, dgbloch2wane, sthbloch2wane
    USE io_ahc,        ONLY : read_dwmat
    USE kfold,         ONLY : ktokpmq
    USE symmetry,      ONLY : unfold_sthmat
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iq
    !! q-point index
    INTEGER, INTENT(in) :: nqc
    !! number of qpoints in the coarse grid
    REAL(KIND = DP), INTENT(IN) :: xqc(3, nqc)
    !! qpoint list, coarse mesh
    INTEGER, INTENT(IN) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(IN) :: irvec_k(3, nrr_k)
    !! Integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors for electrons
    REAL(KIND = DP), INTENT(IN) :: wslen_k(3, nrr_k)
    !! Wigner-Size supercell vectors, store in real instead of integer
    COMPLEX(KIND = DP), INTENT(IN) :: epmatwe_iq(nbndsub, nbndsub, nrr_k, nmodes)
    !! e-p matrix in wannier basis - electrons, for a single q point
    !
    INTEGER :: ik
    !! Counter on coarse k-point grid
    INTEGER :: ibnd, jbnd, ibnd1
    !! Counter on band
    INTEGER :: imode, jmode
    !! Counter on mode
    INTEGER :: nkq
    !! Index of k+q-point in the pool
    INTEGER :: nkq_abs
    !! Absolute index of k+q-point
    INTEGER :: ipool
    !! Cpu index.
    INTEGER :: lrsthmatw
    !! record length for sthmatwe
    INTEGER :: lrdgmatw
    !! record length for dgmatwe
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xq(3)
    !! q-point vector
    REAL(KIND = DP), ALLOCATABLE :: etk_opt(:, :)
    !! Hamiltonian eigenvalues at k within the outer window
    REAL(KIND = DP), ALLOCATABLE :: etq_opt(:, :)
    !! Hamiltonian eigenvalues at k+q within the outer window
    !
    ALLOCATE(etk_opt(nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_iq', 'Error allocating etk_opt', 1)
    ALLOCATE(etq_opt(nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_iq', 'Error allocating etq_opt', 1)
    !
    xq = xqc(:, iq)
    lrsthmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
    lrdgmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
    !
    ! If xq = Gamma, compute Debye-Waller matrix
    !
    IF (ALL(ABS(xq) < eps8)) THEN
      !
      ! Compute momentum matrix elements
      CALL dmebloch2wan(nbnd, nbndsub, nks, nkstot, pmec, xk_loc, cu, nrr_k, irvec_k, wslen_k, lwin, exband, cpmew)
      !
      CALL read_dwmat(lwin)
      !
      CALL ephbloch2wane(iq, xq, nbndep, nbndsub, 3 * nmodes, nks, nkstot, &
        xk_loc, cu, cuq, dw_mat, nrr_k, irvec_k, wslen_k, dwmatwe, is_dw = .TRUE.)
      !
    ENDIF
    !
    ! Slim down etk and etq to contain states only inside the outer window
    !
    etk_opt = zero
    etq_opt = zero
    !
    DO ik = 1, nks
      !
      ! slim down etk
      !
      ibnd = 0
      ibnd1 = 0
      DO jbnd = 1, nbnd
        IF (exband(jbnd)) CYCLE
        !
        ibnd1 = ibnd1 + 1
        IF (lwin(ibnd1, ik)) THEN
          ibnd = ibnd + 1
          etk_opt(ibnd, ik) = et_loc(jbnd, ik)
        ENDIF
        !
      ENDDO
      !
      ! slim down etq
      !
      CALL ktokpmq(xk_loc(:, ik), xq, +1, ipool, nkq, nkq_abs)
      !
      ibnd = 0
      ibnd1 = 0
      DO jbnd = 1, nbnd
        IF (exband(jbnd)) CYCLE
        !
        ibnd1 = ibnd1 + 1
        IF (lwinq(ibnd1, ik)) THEN
          ibnd = ibnd + 1
          etq_opt(ibnd, ik) = et_all(jbnd, nkq_abs)
        ENDIF
        !
      ENDDO
      !
    ENDDO ! ik
    !
    ! Compute the correction to the change of the hopping parameter
    !
    DO imode = 1, nmodes
      CALL dgbloch2wane(nbndep, nbndsub, nks, nkstot, etk_opt, etq_opt, &
          ahc_win_min, ahc_win_max, xk_loc, xq, cu, cuq, &
          epmatq(:, :, :, imode, iq), nrr_k, irvec_k, wslen_k, &
          dgmatwe(:, :, :, imode))
    ENDDO
    !
    ! Add epmat conribution to dgmat
    !
    dgmatwe = dgmatwe + epmatwe_iq
    !
    ! Get sthmatq at iq.
    ! At an irreducible q-point, read directly from file.
    ! At other q-points, rotate using overlap matrix dmat.
    !
    CALL unfold_sthmat(iq, nqc, xqc, lwin, lwinq)
    !
    ! Compute the upper Fan matrix element in electron Wannier representation
    !
    DO imode = 1, nmodes
      DO jmode = 1, nmodes
        CALL sthbloch2wane(nbndep, nbndsub, nks, nkstot, etk_opt, etq_opt, &
            ahc_win_min, ahc_win_max, xk_loc, cu, cuq, &
            epmatq(:, :, :, imode, iq), epmatq(:, :, :, jmode, iq), &
            sthmatq(:, :, :, imode, jmode), &
            nrr_k, irvec_k, wslen_k, sthmatwe(:, :, :, imode, jmode))
      ENDDO
    ENDDO
    !
    ! Write sthmatwe and dgmatwe to file
    !
    IF (ionode) THEN
      DO imode = 1, nmodes
        CALL davcio(sthmatwe(:, :, :, :, imode), lrsthmatw, iusthwe, &
          imode + (iq - 1) * nmodes, +1)
      ENDDO
      !
      CALL davcio(dgmatwe, lrdgmatw, iudgwe, iq, +1)
    ENDIF
    !
    DEALLOCATE(etk_opt, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_iq', 'Error deallocating etk_opt', 1)
    DEALLOCATE(etq_opt, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_iq', 'Error deallocating etq_opt', 1)
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE wfpt_bloch2wan_iq
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE wfpt_bloch2wan_finalize()
    !------------------------------------------------------------------------------
    !!
    !! Finalize (deallocate) variables used for Bloch to Wannier transformation in WFPT
    !!
    !------------------------------------------------------------------------------
    !
    USE io_global, ONLY : ionode
    USE global_var,ONLY : pmec, dw_mat, sthmatq, dgmatwe, sthmatwe
    USE io_var,    ONLY : iusthwe, iudgwe
    USE symmetry,  ONLY : sthmatq_save, sym_ktok
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (ionode) THEN
      CLOSE(iusthwe, STATUS = 'keep')
      CLOSE(iudgwe, STATUS = 'keep')
    ENDIF
    !
    DEALLOCATE(dw_mat, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating dw_mat', 1)
    DEALLOCATE(pmec, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating pmec', 1)
    DEALLOCATE(sthmatq, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating sthmatq', 1)
    DEALLOCATE(sthmatq_save, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating sthmatq_save', 1)
    DEALLOCATE(sym_ktok, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating sym_ktok', 1)
    DEALLOCATE(sthmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating sthmatwe', 1)
    DEALLOCATE(dgmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_bloch2wan_finalize', 'Error deallocating dgmatwe', 1)
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE wfpt_bloch2wan_finalize
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE wfpt_wan2bloch_setup(nrr_k, nrr_g, irvec_k, irvec_g, ndegen_k, ndegen_g, dims, dims2)
    !------------------------------------------------------------------------------
    !!
    !! Initialize variables used for Wannier to Bloch transformation in WFPT
    !!
    !------------------------------------------------------------------------------
    !
    USE modes,            ONLY : nmodes
    USE ions_base,        ONLY : nat
    USE ep_constants,     ONLY : czero, zero
    USE input,            ONLY : elecselfen, nstemp
    USE global_var,       ONLY : nbndfst, nkf, nktotf, dwf17, sigma_ahc_hdw, &
                                 sigma_ahc_uf, sigmar_dw_all
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier function if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: dims2
    !! Is equal to nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: nrr_k
    !! Number of real-space grid points for electrons
    INTEGER, INTENT(in) :: nrr_g
    !! Number of real-space grid points for electron-phonon
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Integer components of the ir-th Wigner-Seitz grid point for electrons
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Integer components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz degeneracy for electrons
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
    !! Wigner-Seitz degeneracy for electron-phonon
    !
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(dwmatf_trunc(nbndfst, nbndfst, 3, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating dwmatf_trunc', 1)
    ALLOCATE(dwf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating dwf17', 1)
    dwmatf_trunc = czero
    dwf17 = czero
    !
    IF (elecselfen) THEN
      ALLOCATE(sigmar_dw_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating sigmar_dw_all', 1)
      ALLOCATE(sigma_ahc_hdw(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating sigma_ahc_hdw', 1)
      ALLOCATE(sigma_ahc_uf(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating sigma_ahc_uf', 1)
      sigmar_dw_all = zero
      sigma_ahc_hdw = zero
      sigma_ahc_uf = zero
    ENDIF
    !
    CALL compute_dw_truncated(nrr_k, nrr_g, irvec_k, irvec_g, ndegen_k, ndegen_g, dims, dims2)
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE wfpt_wan2bloch_setup
    !------------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE compute_dw_truncated(nrr_k, nrr_g, irvec_k, irvec_g, ndegen_k, ndegen_g, dims, dims2)
    !--------------------------------------------------------------------------
    !! Compute DW term truncated to the active space at ik.
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : czero, ryd2ev, ci, eps8, twopi, zero
    USE modes,         ONLY : nmodes
    USE ions_base,     ONLY : nat
    USE input,         ONLY : nbndsub, ahc_win_min, ahc_win_max
    USE global_var,    ONLY : ibndmin, nbndfst, chw, nkf, xkf
    USE wannier2bloch, ONLY : dmewan2bloch, ephwan2blochp, ephwan2bloch, hamwan2bloch
    USE wigner,        ONLY : wigner_divide_ndegen
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: dims
    !! Is equal to the number of Wannier function if use_ws == .TRUE. Is equal to 1 otherwise.
    INTEGER, INTENT(in) :: dims2
    !! Is equal to nat if use_ws or 1 if not
    INTEGER, INTENT(in) :: nrr_k
    !! Number of real-space grid points for electrons
    INTEGER, INTENT(in) :: nrr_g
    !! Number of real-space grid points for electron-phonon
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Integer components of the ir-th Wigner-Seitz grid point for electrons
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Integer components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz degeneracy for electrons
    INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
    !! Wigner-Seitz degeneracy for electron-phonon
    !
    ! Local variables
    INTEGER :: ik
    !! kpoint index
    INTEGER :: ib
    !! Band index
    INTEGER :: jb
    !! Band index
    INTEGER :: pb
    !! Band index
    INTEGER :: idir
    !! Cartesian direction
    INTEGER :: imode
    !! Mode index
    INTEGER :: ib_full
    !! Index on the full band space
    INTEGER :: jb_full
    !! Index on the full band space
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xxk(3)
    !! k vector
    REAL(KIND = DP) :: xxq(3)
    !! q vector
    REAL(KIND = DP) :: etf(nbndsub)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP) :: etf_ks(nbndsub)
    !! Eigenenergies on the fine grid
    REAL(KIND = DP) :: rdotk(nrr_k)
    !! Exponential factor
    COMPLEX(KIND = DP) :: rrf(3, nbndsub, nbndsub)
    !! position matrix elements
    COMPLEX(KIND = DP) :: cuf(nbndsub, nbndsub)
    !! Rotation matrix U^\dagger, fine mesh
    COMPLEX(KIND = DP) :: cfac(nrr_k)
    !! Exponential factor
    COMPLEX(KIND = DP) :: epmatf(nbndsub, nbndsub, nmodes)
    !! e-p matrix in smooth Bloch basis, fine mesh, phonon Catresian basis
    COMPLEX(KIND = DP) :: pmef(3, nbndsub, nbndsub)
    !! momentum matrix elements on the fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :, :)
    !! e-p matrix in Wannier basis, fine mesh, phonon Catresian basis
    !
    ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_dw_truncated', 'Error allocating epmatwef', 1)
    !
    ! Debye-Waller matrix element in phonon Cartesian basis is q-independent
    ! and is defined in terms of the e-ph matrix elements at q=0.
    ! So, we compute the matrix at only once.
    !
    dwmatf_trunc(:, :, :, :, :) = czero
    rrf(:, :, :) = czero
    !
    xxq = (/ zero, zero, zero /)
    !
    CALL ephwan2blochp(nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatwef, nbndsub, nrr_k, dims, dims2)
    CALL wigner_divide_ndegen(epmatwef, 1, nbndsub, nrr_k, nmodes, ndegen_k, dims)
    !
    DO ik = 1, nkf
      !
      xxk = xkf(:, 2*ik-1)
      !
      CALL DGEMV('t', 3, nrr_k, twopi, REAL(irvec_k, KIND=DP), 3, xxk, 1, 0.0_DP, rdotk, 1)
      cfac(:) = EXP(ci * rdotk(:))
      !
      CALL hamwan2bloch(nbndsub, nrr_k, cuf, etf, chw, cfac, .FALSE.)
      etf_ks(:) = zero
      !
      CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cuf, cuf, epmatf, nmodes, cfac, xxq, rrf)
      CALL dmewan2bloch(nbndsub, nrr_k, cuf, pmef, etf, etf_ks, cfac, use_momentum = .TRUE.)
      !
      DO pb = 1, nbndsub
        !
        ! Skip states outside the AHC window
        IF (etf(pb) < ahc_win_min .OR. etf(pb) > ahc_win_max) CYCLE
        !
        DO ib = 1, nbndfst
          DO jb = 1, nbndfst
            ib_full = ib + ibndmin - 1
            jb_full = jb + ibndmin - 1
            !
            DO idir = 1, 3
              DO imode = 1, nmodes
                dwmatf_trunc(ib, jb, idir, imode, ik) &
                = dwmatf_trunc(ib, jb, idir, imode, ik) &
                + ci * CONJG(epmatf(pb, ib_full, imode)) * pmef(idir, pb, jb) &
                - ci * CONJG(pmef(idir, pb, ib)) * epmatf(pb, jb_full, imode)
              ENDDO
            ENDDO
            !
          ENDDO ! jb
        ENDDO ! ib
      ENDDO ! pb
    ENDDO ! ik
    !
    DEALLOCATE(epmatwef, STAT = ierr)
    IF (ierr /= 0) CALL errore('compute_dw_truncated', 'Error deallocating epmatwef', 1)
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE compute_dw_truncated
    !--------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    SUBROUTINE wfpt_wan2bloch_finalize()
    !------------------------------------------------------------------------------
    !!
    !! Finalize (deallocate) variables used for Wannier to Bloch transformation in WFPT
    !!
    !------------------------------------------------------------------------------
    !
    USE input,         ONLY : elecselfen
    USE global_var,    ONLY : dwf17, sigma_ahc_hdw, sigma_ahc_uf, sigmar_dw_all
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    DEALLOCATE(dwmatf_trunc, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_finalize', 'Error deallocating dwmatf_trunc', 1)
    DEALLOCATE(dwf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_finalize', 'Error deallocating dwf17', 1)
    !
    IF (elecselfen) THEN
      DEALLOCATE(sigmar_dw_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('wfpt_wan2bloch_finalize', 'Error deallocating sigmar_dw_all', 1)
      DEALLOCATE(sigma_ahc_hdw, STAT = ierr)
      IF (ierr /= 0) CALL errore('wfpt_wan2bloch_finalize', 'Error deallocating sigma_ahc_hdw', 1)
      DEALLOCATE(sigma_ahc_uf, STAT = ierr)
      IF (ierr /= 0) CALL errore('wfpt_wan2bloch_finalize', 'Error deallocating sigma_ahc_uf', 1)
    ENDIF
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE wfpt_wan2bloch_finalize
    !------------------------------------------------------------------------------
    !
    !
    !--------------------------------------------------------------------------
    SUBROUTINE selfen_elec_ahc_static(iq, nqc, wf_coarse)
    !--------------------------------------------------------------------------
    !!
    !!  Compute the electron self-energy due to electron-phonon interaction
    !!  using the Allen-Heine-Cardona (AHC) theory.
    !!
    !!  Use e-ph matrix elements, Debye-Waller matrix, upper Fan matrix,
    !!  electronic eigenvalues, and phonon frequencies from Wannier
    !!  interpolation.
    !!
    !!  This routine computes the contribution from the coarse phonon iq
    !!  point to all k-points. The outer loop in ephwann_shuffle.f90 over the
    !!  coarse q-grid will loop over all iq points.
    !!  The contribution from each iq is summed at the end of this subroutine
    !!  to recover the per-ik electron self energy
    !!
    !!  In contrast to selfen_elec_q, here this routine computes both the real
    !!  and imaginary part of the self-energy. This can be done because the
    !!  Debye-Waller and upper Fan matrices are interpolated.
    !!  Also, the outer loop loops only over the coarse q points. This is done
    !!  because the sum of the upper Fan and Debye-Waller self-energy converge
    !!  much faster than the lowe Fan self-energy. Accordingly, the phonon
    !!  momentum interpolation of upper Fan matrix is not implemented.
    !!
    !!
    !!  Implemented by Jae-Mo Lihm.
    !!
    !!  [1] J.-M. Lihm and C.-H. Park, PRB 101, 121102(R) (2020)
    !!  [2] J.-M. Lihm and C.-H. Park, PRX 11, 041053 (2021)
    !!
    !--------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE ep_constants,  ONLY : one, two, zero, czero, ryd2ev
    USE mp,            ONLY : mp_barrier, mp_sum
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE input,         ONLY : nbndsub, nstemp, eps_acoustic, ahc_win_min, &
                              ahc_win_max, lwfpt
    USE global_var,    ONLY : etf, ibndmin, nbndfst, nkf, epf17, gtemp, &
                              sigma_ahc_hdw, sigma_ahc_uf, lower_bnd, sthf17, dwf17
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Q-point index from full grid
    INTEGER, INTENT(in) :: nqc
    !! Total number of coarse q-points
    REAL(KIND = DP), INTENT(in) :: wf_coarse(nmodes)
    !! Phonon frequency at coarse q-points
    !
    ! Local variables
    LOGICAL :: skip_mode(nmodes)
    !! If the phonon frequency is too small or negative, skip mode.
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ik_gl
    !! Global k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: itemp
    !! Counter on temperatures
    REAL(KIND = DP) :: coeff
    !! Self-energy coefficient
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: ekq
    !! Eigen energy at k+q on the fine grid relative to the Fermi level
    REAL(KIND = DP) :: inv_eptemp
    !! Inverse of temperature define for efficiency reasons
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Fermi-Dirac distribution function (when -99)
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! This function computes the derivative of the Fermi-Dirac function
    !! It is therefore an approximation for a delta function
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency on the fine grid
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
    REAL(KIND = DP) :: wgq(nmodes)
    !! Bose occupation factor $n_{q\nu}(T)$
    COMPLEX(KIND = DP) :: sthmat(nbndfst, nmodes)
    !! Sternheimer matrix at each k-point
    !
    IF (.NOT. lwfpt) CALL errore("selfen_elec_ahc_static", &
      "To compute static AHC self-energies, wfpt must be true", 1)
    !
    IF (iq == 1) THEN
      !
      WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
      WRITE(stdout, '(5x, "Electron Self-Energy using the Allen-Heine-Cardona Theory")')
      WRITE(stdout, '(5x, " - Static part")')
      WRITE(stdout, '(5x, a)') REPEAT('=', 67)
      !
      ! Write ahc window
      WRITE(stdout, '(5x, a, f10.3, a, f10.3, a)' ) &
        'AHC active state window: ', ahc_win_min * ryd2ev, ' eV to ', ahc_win_max * ryd2ev, ' eV'
      WRITE(stdout, *)
      !
    ENDIF
    !
    ! loop over all k points of the fine mesh
    !
    DO ik = 1, nkf
      !
      ik_gl = ik + lower_bnd - 1
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      sthmat(:, :) = czero
      !
      ! Sternheimer contribution to the Upper Fan matrix.
      !
      DO imode = 1, nmodes
        DO ibnd = 1, nbndfst
          sthmat(ibnd, imode) = sthf17(ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, imode, ik)
        ENDDO
      ENDDO
      !
      ! Add inactive Wannier band contribution to upper Fan matrix.
      !
      DO imode = 1, nmodes
        DO ibnd = 1, nbndfst
          !
          ! the energy of the electron at k
          ekk = etf(ibndmin - 1 + ibnd, ikk)
          !
          ! JML : If ekk etf(ibndmin - 1 + ibnd, ikk) outside
          !
          DO jbnd = 1, nbndsub
            !
            ! the energy of the electron at k+q
            ekq = etf(jbnd, ikq)
            !
            ! Skip active states inside the ahc window
            IF (ahc_win_min <= ekq .AND. ekq <= ahc_win_max) CYCLE
            !
            sthmat(ibnd, imode) = sthmat(ibnd, imode) &
            + CONJG(epf17(jbnd, ibnd, imode, ik)) * epf17(jbnd, ibnd, imode, ik) / (ekk - ekq)
            !
          ENDDO ! jbnd
        ENDDO ! ibnd
      ENDDO ! imode
      !
      ! Calculate self-energy
      !
      DO itemp = 1, nstemp ! loop over temperatures
        !
        inv_eptemp = one / gtemp(itemp)
        !
        wq(:) = zero
        skip_mode(:) = .FALSE.
        !
        DO imode = 1, nmodes
          wq(imode) = wf_coarse(imode)
          IF (wq(imode) > eps_acoustic) THEN
            skip_mode(imode) = .FALSE.
            wgq(imode)    = wgauss(-wq(imode) * inv_eptemp, -99)
            wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
            inv_wq(imode) = one / (two * wq(imode))
          ELSE
            skip_mode(imode) = .TRUE.
            wgq(imode)    = zero
            inv_wq(imode) = zero
          ENDIF
          !
          IF (skip_mode(imode)) CYCLE
          !
          ! coeff = (n + 0.5) / (2 * omega) / nqc
          coeff = (wgq(imode) + one / two) * inv_wq(imode) / REAL(nqc, DP)
          !
          DO ibnd = 1, nbndfst
            !
            ekk = etf(ibndmin - 1 + ibnd, ikk)
            !
            ! Upper Fan
            sigma_ahc_uf(ibnd, ik_gl, itemp) = sigma_ahc_uf(ibnd, ik_gl, itemp) &
              + coeff * REAL((sthmat(ibnd, imode) + CONJG(sthmat(ibnd, imode))))
            !
            ! Debye-Waller
            sigma_ahc_hdw(ibnd, ik_gl, itemp) = sigma_ahc_hdw(ibnd, ik_gl, itemp) &
            + coeff * REAL(dwf17(ibnd, ibnd, imode, ik))
            !
          ENDDO ! ibnd
          !
        ENDDO ! imode
      ENDDO ! itemp
      !
    ENDDO ! ik
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE selfen_elec_ahc_static
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    SUBROUTINE ahc_run_static_wfpt(nrr_k, nrr_q, nrr_g, irvec_k, irvec_q, irvec_g, &
                                   ndegen_k, ndegen_q, ndegen_g, dims, dims2, rws, nrws)
    !--------------------------------------------------------------------------
    !!
    !! Allen-Heine-Cardona theory for the real part of electron self-energy
    !!
    !! The Debye-Waller and upper Fan terms are interpolated only for k-points,
    !! not for q-points. So, we calculate them here, after the q-point loop.
    !!
    ! -------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE mp_world,         ONLY : world_comm
    USE mp,               ONLY : mp_bcast
    USE io_global,        ONLY : ionode_id, ionode
    USE io_files,         ONLY : prefix, diropn
    USE cell_base,        ONLY : at, bg
    USE ions_base,        ONLY : amass, ityp
    USE modes,            ONLY : nmodes
    USE ep_constants,     ONLY : eps8, czero, twopi, ci, cone, zero
    USE input,            ONLY : nbndsub, lifc, nqc1, nqc2, nqc3, eig_read, nw_specfun, &
                                 nstemp, specfun_el, ahc_win_min, ahc_win_max
    USE io_var,           ONLY : iusthwe, iudgwe, iuxqc
    USE global_var,       ONLY : dwmatwe, dgmatwe, sthmatwe, dwf17, sthf17, dgf17,   &
                                 epf17, nkf, nbndfst, xkf,etf, etf_ks, chw_ks, chw,  &
                                 ibndmin, ibndmax, esigmar_all, sigma_ahc_hdw, sigma_ahc_uf, &
                                 lower_bnd
    USE wigner,           ONLY : wigner_divide_ndegen
    USE wannier2bloch,    ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch,           &
                                 ephwan2blochp, ephwan2bloch, dynifc2blochf,         &
                                 rrwan2bloch, dwwan2blochp, sthwan2blochp,           &
                                 dgwan2blochp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nrr_k
    !! Number of WS points for electrons
    INTEGER, INTENT(IN) :: nrr_q
    !! Number of WS points for phonons
    INTEGER, INTENT(IN) :: nrr_g
    !! Number of WS points for electron-phonons
    INTEGER, INTENT(IN) :: nrws
    !! Number of real-space Wigner-Seitz
    INTEGER, INTENT(IN) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(IN) :: dims2
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(IN) :: irvec_k(3, nrr_k)
    !! INTEGER components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, INTENT(IN) :: irvec_q(3, nrr_q)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, INTENT(IN) :: irvec_g(3, nrr_g)
    !! INTEGER components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, INTENT(IN) :: ndegen_k(nrr_k, dims, dims)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, INTENT(IN) :: ndegen_q(nrr_q, dims2, dims2)
    !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, INTENT(IN) :: ndegen_g(dims, nrr_g, dims2)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    !! atomic positions $R - \tau(na)$
    REAL(KIND = DP), INTENT(IN) :: rws(0:3, nrws)
    !! Real-space wigner-Seitz vectors
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: imode, nu, mu
    !! Counter on modes
    INTEGER :: na
    !! Counter on atoms
    INTEGER :: ik, ikk, ikq
    !! Counter on k points
    INTEGER :: iw
    !! Counter on Wannier functions
    INTEGER :: ibnd, jbnd
    !! Counter on bands
    INTEGER :: iq
    !! Counter on coarse q-point grid
    INTEGER :: itemp
    !! Counter on temperatures
    INTEGER :: nqc_ahc
    !! Coarse q-point grid size
    INTEGER :: lrsthmatw
    !! record length while reading file
    INTEGER :: lrdgmatw
    !! record length while reading file
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: ekk
    !! Eigen energy at k on the fine grid
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: xkk(3)
    !! Current k-point on the fine grid
    REAL(KIND = DP) :: xkq2(3)
    !! Current k+q point on the fine grid
    REAL(KIND = DP), ALLOCATABLE :: wf_coarse(:, :)
    !! Phonon frequency at coarse q-points
    REAL(KIND = DP), ALLOCATABLE :: w2(:)
    !! Interpolated phonon frequency
    REAL(KIND = DP), ALLOCATABLE :: xqc_loc(:, :)
    !! Coarse q-point grid
    REAL(KIND = DP), ALLOCATABLE :: irvec_r(:, :)
    !! Wigner-Size supercell vectors, store in real instead of integer
    REAL(KIND = DP), ALLOCATABLE :: rdotk(:)
    !! $r\cdot k$
    REAL(KIND = DP), ALLOCATABLE :: rdotk2(:)
    !! $r\cdot k$
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
    !! Rotation matrix, fine mesh, points k
    COMPLEX(KIND = DP), ALLOCATABLE :: cufkq(:, :)
    !! Rotation matrix, fine mesh, points k+q
    COMPLEX(KIND = DP), ALLOCATABLE :: uf(:, :)
    !! Rotation matrix for phonons
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :, :)
    !! e-p matrix  in el wannier - fine Bloch phonon grid
    COMPLEX(KIND = DP), ALLOCATABLE :: dwmatf(:, :, :, :)
    !! Debye-Waller matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmatf(:, :, :, :)
    !! Sternheimer matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: dgmatf(:, :, :)
    !! Hopping correction matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: sthmatf2(:, :, :)
    !! Sternheimer matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: dgmatf2(:, :, :)
    !! Hopping correction matrix in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
    !! Used to store $e^{2\pi r \cdot k}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: cfacq(:)
    !! Used to store $e^{2\pi r \cdot k+q}$ exponential
    COMPLEX(KIND = DP), ALLOCATABLE :: rrf(:, :, :)
    !! Interpolation position matrix elements on the fine mesh (ipol, nbnd, nbnd)
    COMPLEX(KIND = DP), ALLOCATABLE :: tmp(:, :, :)
    !! Overlap k and k+q
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatf(:, :, :)
    !! e-p matrix  in smooth Bloch basis, fine mesh
    COMPLEX(KIND = DP), ALLOCATABLE :: dwmatf2(:, :, :)
    !! Debye-Waller matrix in smooth Bloch basis, fine mesh
    !
    CALL start_clock('ep-int-ahc')
    !
    nqc_ahc = nqc1 * nqc2 * nqc3
    !
    ALLOCATE(rrf(3, nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating rrf', 1)
    rrf(:, :, :) = czero
    ALLOCATE(cfac(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cfac', 1)
    ALLOCATE(cfacq(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cfacq', 1)
    ALLOCATE(rdotk(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating rdotk', 1)
    ALLOCATE(rdotk2(nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating rdotk2', 1)
    ! This is simply because dgemv take only real number (not integer)
    ALLOCATE(irvec_r(3, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating irvec_r', 1)
    irvec_r = REAL(irvec_k, KIND = DP)
    ALLOCATE(tmp(nbndfst, nbndfst, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating tmp', 1)
    ALLOCATE(epmatf(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating epmatf', 1)
    !
    ! Zeroing everything - initialization is important !
    cfac(:)  = czero
    cfacq(:) = czero
    rdotk(:)     = zero
    rdotk2(:)    = zero
    !
    ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating epmatwef', 1)
    ALLOCATE(w2(nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating w2', 1)
    ALLOCATE(cufkk(nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cufkk', 1)
    ALLOCATE(cufkq(nbndsub, nbndsub), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating cufkq', 1)
    ALLOCATE(uf(nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating uf', 1)
    ALLOCATE(xqc_loc(3, nqc_ahc), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating xqc_loc', 1)
    ALLOCATE(wf_coarse(nmodes, nqc_ahc), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating wf_coarse', 1)
    ALLOCATE(dwmatf(nbndsub, nbndsub, 3, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dwmatf', 1)
    ALLOCATE(sthmatwe(nbndsub, nbndsub, nrr_k, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthmatwe', 1)
    ALLOCATE(dgmatwe(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgmatwe', 1)
    ALLOCATE(sthmatf(nbndsub, nbndsub, nmodes, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthmatf', 1)
    ALLOCATE(dgmatf(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgmatf', 1)
    ALLOCATE(sthmatf2(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthmatf2', 1)
    ALLOCATE(dgmatf2(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgmatf2', 1)
    ALLOCATE(sthf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating sthf17', 1)
    ALLOCATE(dgf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error allocating dgf17', 1)
    ALLOCATE(dwmatf2(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating dwmatf2', 1)
    !
    ! Read xqc from file
    !
    IF (ionode) THEN
      CALL diropn(iuxqc, 'xqc', 3 * nqc_ahc, exst)
      CALL davcio(xqc_loc, 3 * nqc_ahc, iuxqc, 1, -1)
      CLOSE(iuxqc)
      ! Transform to crystal coordinates
      CALL cryst_to_cart(nqc_ahc, xqc_loc, at, -1)
    ENDIF
    CALL mp_bcast(xqc_loc, ionode_id, world_comm)
    !
    ! Open the prefix.sthmatwe file
    IF (ionode) THEN
      lrsthmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint = TRIM(prefix)//'.sthmatwe'
      CALL diropn(iusthwe, 'sthmatwe', lrsthmatw, exst)
      !
      lrdgmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint = TRIM(prefix)//'.dgmatwe'
      CALL diropn(iudgwe, 'dgmatwe', lrdgmatw, exst)
    ENDIF
    !
    CALL wigner_divide_ndegen(dwmatwe, 1, nbndsub, nrr_k, 3 * nmodes, ndegen_k, dims)
    !
    ! Loop over the coarse q-points.
    !
    DO iq = 1, nqc_ahc
      !
      xxq = xqc_loc(:, iq)
      !
      epf17(:, :, :, :) = czero
      dwf17(:, :, :, :) = czero
      sthf17(:, :, :, :) = czero
      dgf17(:, :, :, :) = czero
      !
      ! ------------------------------------------------------
      ! dynamical matrix : Wannier -> Bloch
      ! ------------------------------------------------------
      !
      IF (.NOT. lifc) THEN
        CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq, uf, w2)
      ELSE
        CALL dynifc2blochf(nmodes, rws, nrws, xxq, uf, w2)
      ENDIF
      !
      ! ...then take into account the mass factors and square-root the frequencies...
      !
      DO nu = 1, nmodes
        !
        ! wf are the interpolated eigenfrequencies
        ! (omega on fine grid)
        !
        IF (w2(nu) > -eps8) THEN
          wf_coarse(nu, iq) =  DSQRT(ABS(w2(nu)))
        ELSE
          wf_coarse(nu, iq) = -DSQRT(ABS(w2(nu)))
        ENDIF
        !
        DO mu = 1, nmodes
          na = (mu - 1) / 3 + 1
          uf(mu, nu) = uf(mu, nu) / DSQRT(amass(ityp(na)))
        ENDDO
      ENDDO ! nu
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
      ! --------------------------------------------------------------
      !
      CALL ephwan2blochp(nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatwef, nbndsub, nrr_k, dims, dims2)
      CALL wigner_divide_ndegen(epmatwef, 1, nbndsub, nrr_k, nmodes, ndegen_k, dims)
      !
      ! Read sthmatwe and dgmatwe from file
      !
      sthmatwe = czero
      !
      IF (ionode) THEN
        DO imode = 1, nmodes
          CALL davcio(sthmatwe(:, :, :, :, imode), lrsthmatw, iusthwe, &
            imode + (iq - 1) * nmodes, -1)
        ENDDO
        CALL davcio(dgmatwe, lrdgmatw, iudgwe, iq, -1)
      ENDIF
      !
      CALL mp_bcast(dgmatwe, ionode_id, world_comm)
      CALL mp_bcast(sthmatwe, ionode_id, world_comm)
      !
      CALL wigner_divide_ndegen(sthmatwe, 1, nbndsub, nrr_k, nmodes**2, ndegen_k, dims)
      CALL wigner_divide_ndegen(dgmatwe, 1, nbndsub, nrr_k, nmodes, ndegen_k, dims)
      !
      ! This is a loop over k blocks in the pool (size of the local k-set)
      !
      DO ik = 1, nkf
        !
        ! xkf is assumed to be in crys coord
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        xkk = xkf(:, ikk)
        xkq2 = xkk + xxq
        !
        cufkk(:, :) = czero
        cufkq(:, :) = czero
        !
        CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
        CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq2, 1, 0.0_DP, rdotk2, 1)
        cfac(:)  = EXP(ci * rdotk(:))
        cfacq(:) = EXP(ci * rdotk2(:))
        !
        ! ------------------------------------------------------
        ! hamiltonian : Wannier -> Bloch
        ! ------------------------------------------------------
        !
        ! Kohn-Sham first, then get the rotation matricies for following interp.
        IF (eig_read) THEN
          CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_ks(:, ikk), chw_ks, cfac)
          CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_ks(:, ikq), chw_ks, cfacq)
        ENDIF
        !
        CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ikk), chw, cfac)
        CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf(:, ikq), chw, cfacq)
        !
        ! --------------------------------------------------------------
        ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
        ! --------------------------------------------------------------
        !
        CALL rrwan2bloch(nbndsub, nrr_k, cfac, rrf)
        !
        CALL cryst_to_cart(1, xxq, bg, 1)
        CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, xxq, rrf)
        CALL cryst_to_cart(1, xxq, at, -1)
        !
        ! Store epmatf in memory
        !
        DO jbnd = ibndmin, ibndmax
          DO ibnd = ibndmin, ibndmax
            epf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = epmatf(ibnd, jbnd, :)
          ENDDO
        ENDDO
        !
        ! Now do the eigenvector rotation: epmatf(j) = sum_i eptmp(i) * uf(i,j)
        tmp(:, :, :) = epf17(:, :, :, ik)
        CALL ZGEMM('n', 'n', nbndfst * nbndfst, nmodes, nmodes, cone, tmp(:,:, :), &
                   nbndfst * nbndfst, uf, nmodes, czero, epf17(:, :, :,ik), nbndfst * nbndfst)
        !
        ! -------------------------------------------------------------------
        ! dwmat : Wannier el, phonon Cartesian -> Bloch el, phonon eigenmode
        ! -------------------------------------------------------------------
        ! We calculate Debye-Waller only at q = Gamma.
        !
        dwmatf2(:, :, :) = czero
        !
        ! dwmatf (Debye-Waller in electron Bloch, phonon Cartesian basis)
        ! is independent of q.
        ! To reuse it, we need to store it for all k-points. (Not implemented)
        !
        ! Rotate dwmatwe from electron Wannier to Bloch basis
        ! We need to use cufkk and cufkk, not cufkk and cufkq
        ! nmodes argument is set to 3 * nmodes to loop over both idir and imode of dwmatwe.
        CALL ephwan2bloch(nbndsub, nrr_k, dwmatwe, cufkk, cufkk, dwmatf, &
            3*nmodes, cfac, xxq, rrf, long = .FALSE.)
        !
        ! Subtract the active space term
        dwmatf = dwmatf - dwmatf_trunc(:, :, :, :, ik)
        !
        ! Rotate dwmatf (phonon Cartesian) to dwmatf2 (phonon eigenmode)
        CALL dwwan2blochp(nbndsub, dwmatf, uf, dwmatf2, nmodes)
        !
        ! Slim down to states inside the fsthick window.
        !
        DO jbnd = ibndmin, ibndmax
          DO ibnd = ibndmin, ibndmax
            dwf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = dwmatf2(ibnd, jbnd, :)
          ENDDO
        ENDDO
        !
        ! -------------------------------------------------------------------
        ! sthmat : el Wannier, phonon Cartesian -> el Bloch, phonon eigenmode
        ! -------------------------------------------------------------------
        !
        ! Rotate electron Wannier -> Bloch
        ! We need to use cufkk and cufkk, not cufkk and cufkq
        !
        sthmatf = czero
        CALL ephwan2bloch(nbndsub, nrr_k, sthmatwe, cufkk, cufkk, sthmatf, &
            nmodes**2, cfac, xxq, rrf, long = .FALSE.)
        !
        ! Rotate sthmatf (phonon Cartesian) -> sthmatf2 (phonon eigenmode)
        !
        CALL sthwan2blochp(nbndsub, sthmatf, uf, sthmatf2, nmodes)
        !
        ! -------------------------------------------------------------------
        ! dgmat : el Wannier, phonon Cartesian -> el Bloch, phonon eigenmode
        ! -------------------------------------------------------------------
        !
        dgmatf = czero
        CALL ephwan2bloch(nbndsub, nrr_k, dgmatwe, cufkk, cufkq, dgmatf, &
            nmodes, cfac, xxq, rrf, long = .FALSE.)
        !
        ! Rotate dgmatf (phonon Cartesian) -> dgf (phonon eigenmode)
        CALL dgwan2blochp(nbndsub, dgmatf, uf, dgmatf2, nmodes)
        !
        !
        ! Slim down to states inside the fsthick window.
        !
        DO jbnd = ibndmin, ibndmax
          DO ibnd = ibndmin, ibndmax
            sthf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = sthmatf2(ibnd, jbnd, :)
            dgf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = dgmatf2(ibnd, jbnd, :)
          ENDDO
        ENDDO
        !
      ENDDO ! ik
      !
      CALL selfen_elec_ahc_static(iq, nqc_ahc, wf_coarse(:, iq))
      !
    ENDDO ! iq
    !
    ! Upper Fan self-energy is valid only if the energy is inside the AHC window.
    ! Set it to zero for states outside the window.
    !
    DO itemp = 1, nstemp
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        !
        DO ibnd = 1, nbndfst
          ekk = etf(ibndmin - 1 + ibnd, ikk)
          IF (ekk < ahc_win_min .OR. ahc_win_max < ekk) THEN
            sigma_ahc_uf(ibnd, ik + lower_bnd - 1, itemp) = zero
            sigma_ahc_hdw(ibnd, ik + lower_bnd - 1, itemp) = zero
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
    !
    ! Add WFPT contributions (upper Fan and high-energy DW) to the frequency-dependent
    ! self-energy esigmar_all.
    !
    IF (specfun_el) THEN
      DO itemp = 1, nstemp
        DO iw = 1, nw_specfun
          esigmar_all(:, :, iw, itemp) = esigmar_all(:, :, iw, itemp) + sigma_ahc_uf(:, :, itemp)
          esigmar_all(:, :, iw, itemp) = esigmar_all(:, :, iw, itemp) + sigma_ahc_hdw(:, :, itemp)
        ENDDO
      ENDDO
    ENDIF
    !
    IF (ionode) THEN
      CLOSE(iusthwe, STATUS = 'keep')
      CLOSE(iudgwe, STATUS = 'keep')
    ENDIF
    !
    DEALLOCATE(wf_coarse, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating wf_coarse', 1)
    DEALLOCATE(xqc_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating xqc_loc', 1)
    DEALLOCATE(dwmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dwmatf', 1)
    DEALLOCATE(dwmatf2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dwmatf2', 1)
    DEALLOCATE(dgmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgmatf', 1)
    DEALLOCATE(sthmatf, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthmatf', 1)
    DEALLOCATE(dgmatf2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgmatf2', 1)
    DEALLOCATE(sthmatf2, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthmatf2', 1)
    DEALLOCATE(dwmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dwmatwe', 1)
    DEALLOCATE(dgmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgmatwe', 1)
    DEALLOCATE(sthmatwe, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthmatwe', 1)
    DEALLOCATE(sthf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating sthf17', 1)
    DEALLOCATE(dgf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('ahc_run_static_wfpt', 'Error deallocating dgf17', 1)
    !
    CALL stop_clock('ep-int-ahc')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE ahc_run_static_wfpt
    !--------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE wfpt
  !------------------------------------------------------------------------------
