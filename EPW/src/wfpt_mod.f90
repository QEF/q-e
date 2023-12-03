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
  MODULE wfpt_mod
  !!
  !! This module contains variables and subroutines of polaron
  !! Authored by Jae-Mo Lihm
  !!  [1] J.-M. Lihm and C.-H. Park, PRX 11, 041053 (2021)
  !!
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(KIND = DP), ALLOCATABLE :: dwmatf2(:, :, :)
  !! Debye-Waller matrix in smooth Bloch basis, fine mesh
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
    USE io_global,        ONLY : ionode
    USE io_files,         ONLY : prefix, diropn
    USE symm_base,        ONLY : nsym, time_reversal
    USE pwcom,            ONLY : nks, nkstot
    USE modes,            ONLY : nmodes
    USE constants_epw,    ONLY : czero
    USE epwcom,           ONLY : nbndsub
    USE elph2,            ONLY : nbndep, dw_mat, dwmatwe, cpmew, sthmatq, dgmatwe, &
                                 sthmatwe
    USE io_var,           ONLY : iusthwe, iudgwe
    USE rotate_wavefunction, ONLY : read_sym_ktok, sthmatq_save
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
    USE io_global,        ONLY : ionode
    USE pwcom,            ONLY : nbnd, nks, nkstot
    USE modes,            ONLY : nmodes
    USE constants_epw,    ONLY : eps8, zero
    USE klist_epw,        ONLY : et_loc, xk_loc, et_all
    USE epwcom,           ONLY : nbndsub, ahc_win_min, ahc_win_max
    USE elph2,            ONLY : pmec, exband, nbndep, cu, cuq, lwin, lwinq, dwmatwe, &
                                 dw_mat, cpmew, epmatq, dgmatwe, sthmatwe, sthmatq
    USE io_var,           ONLY : iusthwe, iudgwe
    USE bloch2wan,        ONLY : dmebloch2wan, ephbloch2wane, dgbloch2wane, sthbloch2wane
    USE read_ahc_files,   ONLY : read_dwmat
    USE kfold,            ONLY : ktokpmq
    USE rotate_wavefunction, ONLY : unfold_sthmat
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
    USE io_global,        ONLY : ionode
    USE elph2,            ONLY : pmec, dw_mat, sthmatq, dgmatwe, sthmatwe
    USE io_var,           ONLY : iusthwe, iudgwe
    USE rotate_wavefunction, ONLY : sthmatq_save, sym_ktok
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
    SUBROUTINE wfpt_wan2bloch_setup()
    !------------------------------------------------------------------------------
    !!
    !! Initialize variables used for Wannier to Bloch transformation in WFPT
    !!
    !------------------------------------------------------------------------------
    !
    USE modes,            ONLY : nmodes
    USE constants_epw,    ONLY : czero
    USE epwcom,           ONLY : nbndsub
    USE elph2,            ONLY : nbndfst, nkf, dwf17
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(dwmatf2(nbndsub, nbndsub, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating dwmatf2', 1)
    ALLOCATE(dwmatf_trunc(nbndfst, nbndfst, 3, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating dwmatf_trunc', 1)
    ALLOCATE(dwf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('wfpt_wan2bloch_setup', 'Error allocating dwf17', 1)
    dwmatf2 = czero
    dwmatf_trunc = czero
    dwf17 = czero
    !
    !------------------------------------------------------------------------------
    END SUBROUTINE wfpt_wan2bloch_setup
    !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE wfpt_mod
  !------------------------------------------------------------------------------
