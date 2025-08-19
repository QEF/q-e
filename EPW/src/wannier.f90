  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE wannier
  !----------------------------------------------------------------------
  !!
  !! This module contains routine to build the Wannier el-ph
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !----------------------------------------------------------------------
    SUBROUTINE build_wannier(nqc, xqc, w_centers, dims, dims2, nrr_k, irvec_k, &
                             ndegen_k, wslen_k, nrr_q, irvec_q, ndegen_q, wslen_q, &
                             nrr_g, irvec_g, ndegen_g, wslen_g)
    !---------------------------------------------------------------------
    !!
    !! Building the quantities using Wannier function in real-space.
    !! Full coarse grid to real-space Wannier transformation
    !-----------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE pwcom,            ONLY : nbnd, nks, nkstot, xk, et
    USE input,            ONLY : et_loc, xk_loc
    USE ions_base,        ONLY : nat, tau
    USE modes,            ONLY : nmodes
    USE input,            ONLY : use_ws, system_2d, nkc1, nkc2, nkc3, lpolar,        &
                                 etf_mem, epwread, epwwrite, epbread, lifc,          &
                                 eig_read, nbndsub, nqc1, nqc2, nqc3, lwfpt,         &
                                 calc_nelec_wann, epw_memdist
    USE ep_constants,     ONLY : zero, czero
    USE io_files,         ONLY : diropn
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE io_var,           ONLY : iunepmatwe
    USE global_var,       ONLY : cu, cuq, lwin, lwinq, chw, chw_ks, cvmew, cdmew,    &
                                 rdw, epmatwp, epmatq, dynq, et_ks, dmec, exband,    &
                                 xkq, nbndep, qrpl, crrw, cpmew, nbndskip,           &
                                 epmatwp_dist
    USE bloch2wannier,    ONLY : hambloch2wan, dmebloch2wan, dynbloch2wan,           &
                                 vmebloch2wan, ephbloch2wane, ephbloch2wanp,         &
                                 ephbloch2wanp_mem
    USE wigner,           ONLY : wigner_seitz_wrap, wigner_divide_ndegen,            &
                                 wigner_divide_ndegen_epmat,                         &
                                 wigner_divide_ndegen_epmat_dist
    USE io,               ONLY : rwepmatw, epw_read, epw_write, loadumat,            &
                                 epw_read_ws_data, epw_write_ws_data
    USE mp,               ONLY : mp_bcast
    USE mp_world,         ONLY : world_comm
    USE low_lvl,          ONLY : system_mem_usage
    USE utilities,        ONLY : get_nbndskip, epmatwp_redistribution
    USE longrange,        ONLY : epsi_thickn_2d
    USE wfpt,             ONLY : wfpt_bloch2wan_setup, wfpt_bloch2wan_iq,            &
                                 wfpt_bloch2wan_finalize
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nqc
    !! number of qpoints in the coarse grid
    INTEGER, INTENT(inout) :: dims
    !! Dims is either nbndsub if use_ws or 1 if not
    INTEGER, INTENT(inout) :: dims2
    !! Dims is either nat if use_ws or 1 if not
    INTEGER, INTENT(inout) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(inout) :: nrr_q
    !! number of phonon WS points
    INTEGER, INTENT(inout) :: nrr_g
    !! number of el-ph WS points
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_k(:, :)
    !! Coordinates of real space vector for electrons
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_q(:, :)
    !! Coordinates of real space vector for phonons
    INTEGER, ALLOCATABLE, INTENT(inout) :: irvec_g(:, :)
    !! Coordinates of real space vector for electron-phonon
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_k(:, :, :)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_q(:, :, :)
    !! degeneracy of WS points for phonon
    INTEGER, ALLOCATABLE, INTENT(inout) :: ndegen_g(:, :, :)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(KIND = DP), ALLOCATABLE, INTENT(inout) :: wslen_g(:)
    !! real-space length for electron-phonons, in units of alat
    REAL(KIND = DP), INTENT(in) :: xqc(3, nqc)
    !! qpoint list, coarse mesh
    REAL(KIND = DP), INTENT(inout) :: w_centers(3, nbndsub)
    !! Wannier centers
    !
    ! Local variables
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: iq
    !! Counter on coarse q-point grid
    INTEGER :: lrepmatw
    !! record length while reading file
    INTEGER(KIND = 8) :: valueRSS(2)
    !! Return virtual and resisdent memory from system
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: xxq(3)
    !! Current q-point
    REAL(KIND = DP) :: dummy(3)
    !! Dummy variable
    REAL(KIND = DP), ALLOCATABLE :: w2(:)
    !! Interpolated phonon frequency
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwe(:, :, :, :, :)
    !! e-p matrix  in wannier basis - electrons
    COMPLEX(KIND = DP), ALLOCATABLE :: epmatwe_mem(:, :, :, :)
    !! e-p matrix  in wannier basis - electrons (written on disk)
    COMPLEX(KIND = DP), ALLOCATABLE :: A(:, :, :, :)
    !! Berry connection in Wannier representation
    !
    CALL start_clock('build_wannier')
    !
    IF (nbndsub /= nbndep) WRITE(stdout, '(/,5x,a,i4)' ) 'Band disentanglement is used: nbndsub = ', nbndsub
    !
    IF (.NOT. (epwread .AND. .NOT. epbread)) THEN
      ALLOCATE(cu(nbndep, nbndsub, nks), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating cu', 1)
      ALLOCATE(cuq(nbndep, nbndsub, nks), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating cuq', 1)
      ALLOCATE(lwin(nbndep, nks), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating lwin', 1)
      ALLOCATE(lwinq(nbndep, nks), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating lwinq', 1)
      ALLOCATE(exband(nbnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating exband', 1)
      ALLOCATE(A(3, nbndsub, nbndsub, nks), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating A', 1)
      cu(:, :, :)  = czero
      cuq(:, :, :) = czero
      lwin(:, :)   = .FALSE.
      lwinq(:, :)  = .FALSE.
      exband(:)    = .FALSE.
      A(:, :, :, :) = czero
    ENDIF
    !
    ALLOCATE(w2(3 * nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating w2', 1)
    w2(:) = zero
    !
    IF (lpolar) THEN
      WRITE(stdout, '(/,5x,a)' ) 'Computes the analytic long-range interaction for polar materials [lpolar]'
      WRITE(stdout, '(5x,a)' )   ' '
    ENDIF
    IF (.NOT. lpolar .AND. qrpl) THEN
      WRITE(stdout, '(/,5x,a)' ) 'Computes the analytic quadrupole long-range interaction for non-polar materials [Q1,Q2]'
      WRITE(stdout, '(5x,a)' )   ' '
    ENDIF
    !
    ! At this point, we will interpolate the Wannier rep to the Bloch rep
    IF (epwread) THEN
      !
      ! Read all quantities in Wannier representation from file
      ! in parallel case all pools read the same file
      ! Inside we allocate irvec_k, irvec_q, irvec_g, ndegen_k, ndegen_q, ndegen_g,
      !                    wslen_k,  wslen_q,  wslen_g
      IF (use_ws) THEN
        ! Use Wannier-centers to contstruct the WS for electonic part and el-ph part
        ! Use atomic position to contstruct the WS for the phonon part
        dims  = nbndsub
        dims2 = nat
        CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q, irvec_g, &
                               ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                               w_centers, dims, tau, dims2)
      ELSE
        ! Center the WS at Gamma for electonic part, the phonon part and el-ph part
        dims  = 1
        dims2 = 1
        dummy(:) = (/0.0, 0.0, 0.0/)
        CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q, irvec_g, &
                               ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                               dummy, dims, dummy, dims2)
      ENDIF
      !
      ! Determine the size of the respective WS sets based on the length of the matrices
      nrr_k = SIZE(irvec_k, 2)
      nrr_q = SIZE(irvec_q, 2)
      nrr_g = SIZE(irvec_g, 2)
      IF (use_ws) THEN
        WRITE(stdout, '(5x,a)' )    'Construct the Wigner-Seitz cell using Wannier centers and atomic positions '
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ',nrr_k
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ',nrr_q
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electron-phonon ',nrr_g
        WRITE(stdout, '(5x,a,i8)' ) 'Maximum number of cores for efficient parallelization ',nrr_g * nat
      ELSE
        WRITE(stdout, '(5x,a)' )    'Use zone-centred Wigner-Seitz cells '
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ',nrr_k
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ',nrr_q
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electron-phonon ',nrr_g
        WRITE(stdout, '(5x,a,i8)' ) 'Maximum number of cores for efficient parallelization ',nrr_g * nmodes
        WRITE(stdout, '(5x,a)' )    'Results may improve by using use_ws == .TRUE. '
      ENDIF
      !
      CALL epw_read(nrr_k, nrr_q, nrr_g)
      !
      !CALL epw_read_ws_data(dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, &
       ! nrr_q, irvec_q, ndegen_q, wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
      !
      ! Calculate effective dielectric constant and thickness for 2D system (V.-A.H)
      IF (lpolar .AND. (system_2d == 'dipole_sh')) CALL epsi_thickn_2d()
      !
    ELSE !if not epwread (i.e. need to calculate fmt file)
      !
      ! Determine Wigner-Seitz points
      ! For this we need the Wannier centers and w_centers is allocated inside loadumat
      xxq = 0.d0
      ALLOCATE(xkq(3, nkstot), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating xkq', 1)
      CALL loadumat(nbndep, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband, w_centers)
      DEALLOCATE(xkq, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating xkq', 1)
      !
      ! Inside we allocate irvec_k, irvec_q, irvec_g, ndegen_k, ndegen_q, ndegen_g,
      !                    wslen_k,  wslen_q,  wslen_g
      IF (use_ws) THEN
        ! Use Wannier-centers to contstruct the WS for electonic part and el-ph part
        ! Use atomic position to contstruct the WS for the phonon part
        dims  = nbndsub
        dims2 = nat
        CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q, irvec_g, &
                               ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                               w_centers, dims, tau, dims2)
      ELSE
        ! Center the WS at Gamma for electonic part, the phonon part and el-ph part
        dims  = 1
        dims2 = 1
        dummy(:) = (/0.0, 0.0, 0.0/)
        CALL wigner_seitz_wrap(nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, irvec_k, irvec_q, irvec_g, &
                               ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                               dummy, dims, dummy, dims2)
      ENDIF
      !
      ! Determine the size of the respective WS sets based on the length of the matrices
      nrr_k = SIZE(irvec_k, 2)
      nrr_q = SIZE(irvec_q, 2)
      nrr_g = SIZE(irvec_g, 2)
      IF (use_ws) THEN
        WRITE(stdout, '(5x,a)' )    'Construct the Wigner-Seitz cell using Wannier centers and atomic positions '
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ',nrr_k
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ',nrr_q
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electron-phonon ',nrr_g
        WRITE(stdout, '(5x,a,i8)' ) 'Maximum number of cores for efficient parallelization ',nrr_g * nat
      ELSE
        WRITE(stdout, '(5x,a)' )    'Use zone-centred Wigner-Seitz cells '
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electrons ',nrr_k
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for phonons ',nrr_q
        WRITE(stdout, '(5x,a,i8)' ) 'Number of WS vectors for electron-phonon ',nrr_g
        WRITE(stdout, '(5x,a,i8)' ) 'Maximum number of cores for efficient parallelization ',nrr_g * nmodes
        WRITE(stdout, '(5x,a)' )    'Results may improve by using use_ws == .TRUE. '
      ENDIF
      !
      IF (lpolar) THEN
        IF (system_2d == 'dipole_sh') CALL epsi_thickn_2d()
      ENDIF
      ! ------------------------------------------------------
      !   Bloch to Wannier transform
      ! ------------------------------------------------------
      !
      ALLOCATE(chw(nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating chw', 1)
      ALLOCATE(chw_ks(nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating chw_ks', 1)
      ALLOCATE(rdw(nmodes, nmodes, nrr_q), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating rdw', 1)
      ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating cvmew', 1)
      ALLOCATE(crrw(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating crrw', 1)
      ALLOCATE(cdmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating cdmew', 1)
      cvmew(:, :, :, :) = czero
      chw(:, :, :) = czero
      chw_ks(:, :, :) = czero
      rdw(:, :, :) = czero
      crrw(:, :, :, :) = czero
      cdmew(:, :, :, :) = czero
      !
      ! SP : Let the user chose. If false use files on disk
      IF (etf_mem == 0) THEN
        ALLOCATE(epmatwe(nbndsub, nbndsub, nrr_k, nmodes, nqc), STAT = ierr)
        IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating epmatwe', 1)
        ALLOCATE(epmatwp(nbndsub, nbndsub, nrr_k, nmodes, nrr_g), STAT = ierr)
        IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating epmatwp', 1)
        epmatwe(:, :, :, :, :) = czero
        epmatwp(:, :, : ,: ,:) = czero
      ENDIF
      ALLOCATE(epmatwe_mem(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating epmatwe_mem', 1)
      epmatwe_mem(:, :, :, :) = czero
      !
      IF (lwfpt) CALL wfpt_bloch2wan_setup(nrr_k)
      !
      ! Hamiltonian
      !
      CALL hambloch2wan(nbnd, nbndsub, nks, nkstot, et_loc, xk_loc, cu, lwin, exband, nrr_k, irvec_k, wslen_k, chw)
      !
      ! Kohn-Sham eigenvalues
      !
      IF (eig_read) THEN
        WRITE (stdout,'(5x,a)') "Interpolating MB and KS eigenvalues"
        CALL hambloch2wan(nbnd, nbndsub, nks, nkstot, et_ks, xk_loc, cu, lwin, exband, nrr_k, irvec_k, wslen_k, chw_ks)
      ENDIF
      !
      A(:, :, :, :) = czero
      ! Transform of position matrix elements - PRB 74 195118  (2006)
      CALL vmebloch2wan(dims, nbnd, nbndsub, nks, nkstot, xk_loc, cu, nrr_k, irvec_k, wslen_k, lwin, exband, &
                        A(:, :, :, :), w_centers, ndegen_k)
      ! Transform of velocity matrix elements (dH/dk)
      CALL dmebloch2wan(nbnd, nbndsub, nks, nkstot, dmec, xk_loc, cu, nrr_k, irvec_k, wslen_k, lwin, exband, cdmew)
      !
      ! Dynamical Matrix
      !
      IF (.NOT. lifc) CALL dynbloch2wan(nmodes, nqc, xqc, dynq, nrr_q, irvec_q, wslen_q)
      !
      !
      ! Electron-Phonon vertex (Bloch el and Bloch ph -> Wannier el and Bloch ph)
      !
      ! Open the prefix.epmatwe file
      IF ((etf_mem == 1) .AND. ionode) THEN
        lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
        CALL diropn(iunepmatwe, 'epmatwe', lrepmatw, exst)
      ENDIF
      !
      WRITE(stdout, '(a)' ) ' '
      DO iq = 1, nqc
        WRITE(stdout, '(a,i10,a,i10)' ) '     Bloch2wane: ',iq,' / ', nqc1 * nqc2 * nqc3
        !
        xxq = xqc(:, iq)
        !
        ! we need the cu again for the k+q points, we generate the map here
        !
        cuq(:, :, :) = czero
        ALLOCATE(xkq(3, nkstot), STAT = ierr)
        IF (ierr /= 0) CALL errore('build_wannier', 'Error allocating xkq', 1)
        CALL loadumat(nbndep, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband, w_centers)
        DEALLOCATE(xkq, STAT = ierr)
        IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating xkq', 1)
        !
        IF (etf_mem == 0) THEN
          CALL ephbloch2wane(iq, xxq, nbndep, nbndsub, nmodes, nks, nkstot, xk_loc, cu, cuq, &
            epmatq(:, :, :, :, iq), nrr_k, irvec_k, wslen_k, epmatwe(:, :, :, :, iq), A(:, :, :, :))
        ELSE
          CALL ephbloch2wane(iq, xxq, nbndep, nbndsub, nmodes, nks, nkstot, xk_loc, cu, cuq, &
            epmatq(:, :, :, :, iq), nrr_k, irvec_k, wslen_k, epmatwe_mem(:, :, :, :), A(:, :, :, :))
          !
        ENDIF
        ! Only the master node writes
        IF ((etf_mem == 1) .AND. ionode) THEN
          ! direct write of epmatwe for this iq
          CALL rwepmatw(epmatwe_mem, nbndsub, nrr_k, nmodes, iq, iunepmatwe, +1)
          !
        ENDIF
        !
        IF (lwfpt) THEN
          IF (etf_mem == 0) THEN
            epmatwe_mem = epmatwe(:, :, :, :, iq)
          ENDIF
          CALL wfpt_bloch2wan_iq(iq, nqc, xqc, nrr_k, irvec_k, wslen_k, epmatwe_mem)
        ENDIF ! wfpt
        !
      ENDDO ! iq
      !
      IF ((etf_mem == 1) .AND. ionode) CLOSE(iunepmatwe, STATUS = 'keep')
      !
      ! Electron-Phonon vertex (Wannier el and Bloch ph -> Wannier el and Wannier ph)
      !
      IF (etf_mem == 0) THEN
        IF (ionode) CALL ephbloch2wanp(nbndsub, nmodes, xqc, nqc, irvec_k, irvec_g, nrr_k, nrr_g, epmatwe)
        CALL mp_bcast(epmatwp, ionode_id, world_comm)
      ENDIF
      IF (etf_mem > 0) CALL ephbloch2wanp_mem(nbndsub, nmodes, xqc, nqc, irvec_k, irvec_g, nrr_k, nrr_g)
      !
      ! Determine the number of bands to skip when calculating Fermi level
      !
      IF (calc_nelec_wann) nbndskip = get_nbndskip(nkstot, nbnd, xk, et, nrr_k, irvec_k, ndegen_k, dims)
      !
      IF (epwwrite) THEN
        CALL epw_write(nrr_k, nrr_q, nrr_g, w_centers)
        !
        CALL epw_write_ws_data(dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, &
          nrr_q, irvec_q, ndegen_q, wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
      ENDIF
      !
      IF (etf_mem == 0 .AND. epw_memdist) CALL epmatwp_redistribution(nrr_k, nrr_q, nrr_g)
      !
      DEALLOCATE(epmatq, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating epmatq', 1)
      DEALLOCATE(dynq, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating dynq', 1)
      DEALLOCATE(dmec, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating dmec', 1)
      IF (etf_mem == 0) THEN
        DEALLOCATE(epmatwe, STAT = ierr)
        IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating epmatwe', 1)
      ELSE
        DEALLOCATE(epmatwe_mem, STAT = ierr)
        IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating epmatwe_mem', 1)
      ENDIF
      DEALLOCATE(cu, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating cu', 1)
      DEALLOCATE(cuq, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating cuq', 1)
      DEALLOCATE(lwin, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating lwin', 1)
      DEALLOCATE(lwinq, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating lwinq', 1)
      DEALLOCATE(exband, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating exband', 1)
      DEALLOCATE(A, STAT = ierr)
      IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating A', 1)
      !
      IF (lwfpt) CALL wfpt_bloch2wan_finalize()
    ENDIF ! (epwread .AND. .NOT. epbread)
    !
    DEALLOCATE(w2, STAT = ierr)
    IF (ierr /= 0) CALL errore('build_wannier', 'Error deallocating w2', 1)
    !
    ! Divide the degeneracy factor for electron Wannier R-points
    !
    CALL wigner_divide_ndegen(chw, 1, nbndsub, nrr_k, 1, ndegen_k, dims)
    CALL wigner_divide_ndegen(cvmew, 3, nbndsub, nrr_k, 1, ndegen_k, dims)
    CALL wigner_divide_ndegen(crrw, 3, nbndsub, nrr_k, 1, ndegen_k, dims)
    CALL wigner_divide_ndegen(cdmew, 3, nbndsub, nrr_k, 1, ndegen_k, dims)
    IF (lwfpt) CALL wigner_divide_ndegen(cpmew, 3, nbndsub, nrr_k, 1, ndegen_k, dims)
    IF (eig_read) CALL wigner_divide_ndegen(chw_ks, 1, nbndsub, nrr_k, 1, ndegen_k, dims)
    IF (etf_mem == 0) THEN
      IF (.NOT. epw_memdist) THEN
        CALL wigner_divide_ndegen_epmat(epmatwp, nbndsub, nrr_k, nrr_g, nmodes, ndegen_g, dims, dims2)
      ELSE
        CALL wigner_divide_ndegen_epmat_dist(epmatwp_dist, nbndsub, nrr_k, nrr_g, nmodes, ndegen_g, dims, dims2)
      ENDIF
    ENDIF
    !
    ! Check Memory usage
    CALL system_mem_usage(valueRSS)
    !
    WRITE(stdout, '(a)' )             '     ==================================================================='
    WRITE(stdout, '(a,i10,a)' ) '     Memory usage:  VmHWM =',valueRSS(2)/1024,'Mb'
    WRITE(stdout, '(a,i10,a)' ) '                   VmPeak =',valueRSS(1)/1024,'Mb'
    WRITE(stdout, '(a)' )             '     ==================================================================='
    WRITE(stdout, '(a)' )             '     '
    !
    CALL stop_clock('build_wannier')
    !
    !--------------------------------------------------------------------------
    END SUBROUTINE build_wannier
    !--------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------
  END MODULE wannier
  !-------------------------------------------------------------------------------------------
