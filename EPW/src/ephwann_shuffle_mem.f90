  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  SUBROUTINE ephwann_shuffle_mem(nqc, xqc, w_centers)
  !---------------------------------------------------------------------
  !!
  !! Wannier interpolation of electron-phonon vertex
  !!
  !! The routine is divided in two parts.
  !! 1) Full coarse grid to real-space Wannier transformation
  !! 2) Real-space Wannier to fine grid Bloch space interpolation
  !!
  !!  Here we do an additional mode parallelism to save memory
  !!  (only fine interpolation part)
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP, i4b, i8b
  USE pwcom,         ONLY : nbnd, nks, nkstot, ef, nelec
  USE klist_epw,     ONLY : et_loc, xk_loc, isk_dummy
  USE cell_base,     ONLY : at, bg
  USE ions_base,     ONLY : nat, amass, ityp, tau
  USE modes,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange,               &
                            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
                            scr_typ, nw_specfun, specfun_pl, lindabs, use_ws,   &
                            elecselfen, phonselfen, nest_fn, a2f, specfun_ph,   &
                            vme, eig_read, ephwrite, nkf1, nkf2, nkf3,          &
                            efermi_read, fermi_energy, specfun_el, band_plot,   &
                            scattering, nstemp, int_mob, scissor, carrier,      &
                            iterative_bte, longrange, scatread, nqf1, prtgkk,   &
                            nqf2, nqf3, mp_mesh_k, restart, plselfen, epbread,  &
                            epmatkqread, selecqread, restart_step, nsmear,      &
                            nkc1, nkc2, nkc3, nqc1, nqc2, nqc3, assume_metal,   &
                            cumulant, eliashberg, fermi_plot,                   &
                            nomega, omegamin, omegamax, omegastep, neta
  USE control_flags, ONLY : iverbosity
  USE noncollin_module, ONLY : noncolin
  USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, zero, czero, cone,       &
                            twopi, ci, kelvin2eV, eps6, eps8, eps16, byte2Mb
  USE io_files,      ONLY : prefix, diropn, tmp_dir
  USE io_global,     ONLY : stdout, ionode
  USE io_var,        ONLY : lambda_phself, linewidth_phself, iunepmatwe,        &
                            iunepmatwp2, iunrestart, iuntau, iuntaucb
  USE elph2,         ONLY : cu, cuq, lwin, lwinq, map_rebal, map_rebal_inv,     &
                            chw, chw_ks, cvmew, cdmew, rdw, adapt_smearing,     &
                            epmatq, wf, etf, etf_ks, xqf, xkf,                  &
                            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks,    &
                            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef,     &
                            sigmai_all, sigmai_mode, gamma_all, epsi, zstar,    &
                            efnew, sigmar_all, zi_all, nkqtotf, eps_rpa,        &
                            sigmar_all, zi_allvb, inv_tau_all, eta, nbndfst,    &
                            inv_tau_allcb, zi_allcb, exband, nbndskip,          &
                            gamma_v_all, esigmar_all, esigmai_all,              &
                            a_all, a_all_ph, wscache, lambda_v_all, threshold,  &
                            nktotf,  gtemp, xkq, lower_bnd, upper_bnd, dos,&
                            nbndep, ef0_fca, epsilon2_abs, epsilon2_abs_lorenz, &
                            epsilon2_abs_all, epsilon2_abs_lorenz_all
  USE wan2bloch,     ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch,           &
                            ephwan2blochp, ephwan2bloch, vmewan2bloch,          &
                            dynifc2blochf, ephwan2blochp_mem, ephwan2bloch_mem
  USE bloch2wan,     ONLY : hambloch2wan, dmebloch2wan, dynbloch2wan,           &
                            vmebloch2wan, ephbloch2wane, ephbloch2wanp,         &
                            ephbloch2wanp_mem
  USE wigner,        ONLY : wigner_seitz_wrap
  USE io_eliashberg, ONLY : write_ephmat, count_kpoints, kmesh_fine,kqmap_fine, &
                            check_restart_ephwrite
  USE transport,     ONLY : transport_coeffs, scattering_rate_q
  USE grid,          ONLY : qwindow
  USE printing,      ONLY : print_gkk, plot_band, plot_fermisurface
  USE io_epw,        ONLY : rwepmatw, epw_read, epw_write
  USE io_transport,  ONLY : tau_read, iter_open, print_ibte, iter_merge
  USE io_selfen,     ONLY : selfen_el_read, spectral_read
  USE transport_iter,ONLY : iter_restart
  USE close_epw,     ONLY : iter_close
  USE division,      ONLY : fkbounds
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : npool
  USE mp_world,      ONLY : mpime, world_comm
  USE low_lvl,       ONLY : system_mem_usage, mem_size
  USE utilities,     ONLY : compute_dos, fermicarrier, fermiwindow
  USE grid,          ONLY : loadqmesh_serial, loadkmesh_para, load_rebal
  USE selfen,        ONLY : selfen_phon_q, selfen_elec_q, selfen_pl_q,          &
                            nesting_fn_q
  USE spectral_func, ONLY : spectral_func_el_q, spectral_func_ph_q, a2f_main,   &
                            spectral_func_pl_q
  USE rigid_epw,     ONLY : rpa_epsilon, tf_epsilon, compute_umn_f, rgd_blk_epw_fine_mem !, &
!                            find_gmin ! Temporarily commented by H. Lee
  USE indabs,        ONLY : indabs_main, renorm_eig, fermi_carrier_indabs
  USE io_indabs,     ONLY : indabs_read
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_MODE_RDONLY, MPI_INFO_NULL, MPI_OFFSET_KIND, &
                               MPI_OFFSET
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nqc
  !! number of qpoints in the coarse grid
  REAL(KIND = DP), INTENT(in) :: xqc(3, nqc)
  !! qpoint list, coarse mesh
  REAL(KIND = DP), INTENT(in) :: w_centers(3, nbndsub)
  !! Wannier centers
  !
  ! Local  variables
  CHARACTER(LEN = 20) :: tp
  !! string for temperature
  CHARACTER(LEN = 256) :: filephselfen
  !! file name of phonon selfenergy
  CHARACTER(LEN = 256) :: filephlinewid
  !! file name of phonon linewidth
  CHARACTER(LEN = 256) :: filint
  !! Name of the file to write/read
  CHARACTER(LEN = 30)  :: myfmt
  !! Variable used for formatting output
  LOGICAL :: already_skipped
  !! Skipping band during the Wannierization
  LOGICAL :: exst
  !! If the file exist
  LOGICAL :: first_cycle
  !! Check wheter this is the first cycle after a restart.
  LOGICAL :: first_time
  !! Check wheter this is the first timeafter a restart.
  LOGICAL :: homogeneous
  !! Check if the k and q grids are homogenous and commensurate.
  INTEGER :: ios
  !! INTEGER variable for I/O control
  INTEGER :: iq
  !! Counter on coarse q-point grid
  INTEGER :: iqq
  !! Counter on coarse q-point grid
  INTEGER :: iq_restart
  !! Counter on coarse q-point grid
  INTEGER :: ik
  !! Counter on coarse k-point grid
  INTEGER :: ikk
  !! Counter on k-point when you have paired k and q
  INTEGER :: ikq
  !! Paired counter so that q is adjacent to its k
  INTEGER :: ibnd
  !! Counter on band
  INTEGER :: jbnd
  !! Counter on band
  INTEGER :: imode
  !! Counter on mode
  INTEGER :: na
  !! Counter on atom
  INTEGER :: mu
  !! counter on mode
  INTEGER :: nu
  !! counter on mode
  INTEGER :: fermicount
  !! Number of states at the Fermi level
  INTEGER :: lrepmatw
  !! record length while reading file
  INTEGER :: ir
  !! Counter for WS loop
  INTEGER :: nrws
  !! Number of real-space Wigner-Seitz
  INTEGER :: valueRSS(2)
  !! Return virtual and resisdent memory from system
  INTEGER :: ierr
  !! Error status
  INTEGER :: nrr_k
  !! Number of WS points for electrons
  INTEGER :: nrr_q
  !! Number of WS points for phonons
  INTEGER :: nrr_g
  !! Number of WS points for electron-phonons
  INTEGER :: dims
  !! Dims is either nbndsub if use_ws or 1 if not
  INTEGER :: dims2
  !! Dims is either nat if use_ws or 1 if not
  INTEGER :: iw
  !! Counter on bands when use_ws == .TRUE.
  INTEGER :: iw2
  !! Counter on bands when use_ws == .TRUE.
  INTEGER :: itemp
  !! Temperature index
  INTEGER :: itempphen
  !! Temperature counter for writing phonon selfen
  INTEGER :: icbm
  !! Index of the CBM
  INTEGER :: totq
  !! Total number of q-points within the fsthick window.
  INTEGER :: ipool
  !! Cpu index.
  INTEGER :: npool_tmp
  !! Temporary number of pools
  INTEGER, ALLOCATABLE :: irvec_k(:, :)
  !! INTEGER components of the ir-th Wigner-Seitz grid point in the basis
  !! of the lattice vectors for electrons
  INTEGER, ALLOCATABLE :: irvec_q(:, :)
  !! INTEGER components of the ir-th Wigner-Seitz grid point for phonons
  INTEGER, ALLOCATABLE :: irvec_g(:, :)
  !! INTEGER components of the ir-th Wigner-Seitz grid point for electron-phonon
  INTEGER, ALLOCATABLE :: ndegen_k(:, :, :)
  !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
  INTEGER, ALLOCATABLE :: ndegen_q(:, :, :)
  !! Wigner-Seitz weights for the phonon grid that depend on atomic positions $R + \tau(nb) - \tau(na)$
  INTEGER, ALLOCATABLE :: ndegen_g(:, :, :)
  !! Wigner-Seitz weights for the electron-phonon grid that depend on
  !! atomic positions $R - \tau(na)$
  INTEGER, ALLOCATABLE :: selecq(:)
  !! Selected q-points within the fsthick window
  INTEGER, PARAMETER :: nrwsx = 200
  !! Maximum number of real-space Wigner-Seitz
  INTEGER :: lrepmatw2_restart(npool)
  !! To restart opening files
  INTEGER :: lrepmatw5_restart(npool)
  !! To restart opening files
  INTEGER :: ctype
  !! Calculation type: -1 = hole, +1 = electron and 0 = both.
#if defined(__MPI)
  INTEGER(KIND = MPI_OFFSET_KIND) :: ind_tot
  !! Total number of points store on file
  INTEGER(KIND = MPI_OFFSET_KIND) :: ind_totcb
  !! Total number of points store on file (CB)
#else
  INTEGER(KIND = 8) :: ind_tot
  !! Total number of points store on file
  INTEGER(KIND = 8) :: ind_totcb
  !! Total number of points store on file (CB)
#endif
  REAL(KIND = DP) :: xxq(3)
  !! Current q-point
  REAL(KIND = DP) :: xxk(3)
  !! Current k-point on the fine grid
  REAL(KIND = DP) :: xkk(3)
  !! Current k-point on the fine grid
  REAL(KIND = DP) :: xkq2(3)
  !! Current k+q point on the fine grid
  REAL(KIND = DP) :: rws(0:3, nrwsx)
  !! Real-space wigner-Seitz vectors
  REAL(KIND = DP) :: atws(3, 3)
  !! Maximum vector: at*nq
  REAL(KIND = DP) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND = DP) :: etemp_fca
  !! Temperature for free carrier absorption
  REAL(KIND = DP) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND = DP) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp
  REAL(KIND = DP) :: dummy(3)
  !! Dummy variable
  REAL(KIND = DP), EXTERNAL :: efermig
  !! External function to calculate the fermi energy
  REAL(KIND = DP), ALLOCATABLE :: etf_all(:, :)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  REAL(KIND = DP), ALLOCATABLE :: w2(:)
  !! Interpolated phonon frequency
  REAL(KIND = DP), ALLOCATABLE :: irvec_r(:, :)
  !! Wigner-Size supercell vectors, store in real instead of integer
  REAL(KIND = DP), ALLOCATABLE :: rdotk(:)
  !! $r\cdot k$
  REAL(KIND = DP), ALLOCATABLE :: rdotk2(:)
  !! $r\cdot k$
  REAL(KIND = DP), ALLOCATABLE :: wslen_k(:)
  !! real-space length for electrons, in units of alat
  REAL(KIND = DP), ALLOCATABLE :: wslen_q(:)
  !! real-space length for phonons, in units of alat
  REAL(KIND = DP), ALLOCATABLE :: wslen_g(:)
  !! real-space length for electron-phonons, in units of alat
  REAL(KIND = DP), ALLOCATABLE :: vkk_all(:, :, :)
  !! velocity from all the k-point
  REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
  !! k-point weights for all the k-points
  COMPLEX(KIND = DP), ALLOCATABLE :: epmatwe_mem(:, :, :, :)
  !! e-p matrix  in wannier basis - electrons (written on disk)
  COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :)
  !! e-p matrix  in el wannier - fine Bloch phonon grid
  COMPLEX(KIND = DP), ALLOCATABLE :: epmatf(:, :)
  !! e-p matrix  in smooth Bloch basis, fine mesh
  COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
  !! Rotation matrix, fine mesh, points k
  COMPLEX(KIND = DP), ALLOCATABLE :: cufkq(:, :)
  !! the same, for points k+q
  COMPLEX(KIND = DP), ALLOCATABLE :: uf(:, :)
  !! Rotation matrix for phonons
  COMPLEX(KIND = DP), ALLOCATABLE :: bmatf(:, :)
  !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:, :, :)
  !! Used to store $e^{2\pi r \cdot k}$ exponential
  COMPLEX(KIND = DP), ALLOCATABLE :: cfacq(:, :, :)
  !! Used to store $e^{2\pi r \cdot k+q}$ exponential
  COMPLEX(KIND = DP), ALLOCATABLE :: eptmp(:, :, :, :)
  !! Temporary el-ph matrices.
  COMPLEX(KIND = DP), ALLOCATABLE :: epmatlrT(:, :, :, :)
  !! Long-range temp. save
  COMPLEX(KIND = DP), ALLOCATABLE :: vmefp(:, :, :)
  !! Phonon velocity
  !
  CALL start_clock('ephwann')
  !
  IF (nbndsub /= nbndep) WRITE(stdout, '(/,5x,a,i4)' ) 'Band disentanglement is used: nbndsub = ', nbndsub
  !
  IF (.NOT. (epwread .AND. .NOT. epbread)) THEN
    ALLOCATE(cu(nbndep, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cu', 1)
    ALLOCATE(cuq(nbndep, nbndsub, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cuq', 1)
    ALLOCATE(lwin(nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating lwin', 1)
    ALLOCATE(lwinq(nbndep, nks), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating lwinq', 1)
    ALLOCATE(exband(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating exband', 1)
    cu(:, :, :)  = czero
    cuq(:, :, :) = czero
    lwin(:, :)   = .FALSE.
    lwinq(:, :)  = .FALSE.
    exband(:)    = .FALSE.
  ENDIF
  !
  ALLOCATE(w2(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating w2', 1)
  w2(:) = zero
  !
  IF (lpolar) THEN
    WRITE(stdout, '(/,5x,a)') 'Computes the analytic long-range interaction for polar materials [lpolar]'
    WRITE(stdout, '(5x,a)')   ' '
  ENDIF
  !
  ! Determine Wigner-Seitz points
  !
  ! For this we need the Wannier centers
  ! w_centers is allocated inside loadumat
  IF (.NOT. epwread) THEN
    xxq = 0.d0
    ALLOCATE(xkq(3, nkstot), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating xkq', 1)
    CALL loadumat(nbndep, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband, w_centers)
    DEALLOCATE(xkq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating xkq', 1)
  ENDIF
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
  nrr_k = SIZE(irvec_k(1, :))
  nrr_q = SIZE(irvec_q(1, :))
  nrr_g = SIZE(irvec_g(1, :))
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
  ! At this point, we will interpolate the Wannier rep to the Bloch rep
  IF (epwread .AND. .NOT. epbread) THEN
    !
    ! Read all quantities in Wannier representation from file
    ! in parallel case all pools read the same file
    CALL epw_read(nrr_k, nrr_q, nrr_g)
    !
  ELSE !if not epwread (i.e. need to calculate fmt file)
    !
    ! ------------------------------------------------------
    !   Bloch to Wannier transform
    ! ------------------------------------------------------
    !
    ALLOCATE(chw(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating chw', 1)
    ALLOCATE(chw_ks(nbndsub, nbndsub, nrr_k), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating chw_ks', 1)
    ALLOCATE(rdw(nmodes, nmodes, nrr_q), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating rdw', 1)
    IF (vme == 'wannier') THEN
      ALLOCATE(cvmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cvmew', 1)
      cvmew(:, :, :, :) = czero
    ELSE
      ALLOCATE(cdmew(3, nbndsub, nbndsub, nrr_k), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cdmew', 1)
      cdmew(:, :, :, :) = czero
    ENDIF
    !
    ! SP : Let the user chose. If false use files on disk
    ALLOCATE(epmatwe_mem(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epmatwe_mem', 1)
    epmatwe_mem(:, :, :, :) = czero
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
    IF (vme == 'wannier') THEN
      ! Transform of position matrix elements
      ! PRB 74 195118  (2006)
      CALL vmebloch2wan(nbnd, nbndsub, nks, nkstot, xk_loc, cu, nrr_k, irvec_k, wslen_k, lwin, exband)
    ELSE
      ! Dipole
      CALL dmebloch2wan(nbnd, nbndsub, nks, nkstot, dmec, xk_loc, cu, nrr_k, irvec_k, wslen_k, lwin, exband)
    ENDIF
    !
    ! Dynamical Matrix
    !
    IF (.NOT. lifc) CALL dynbloch2wan(nmodes, nqc, xqc, dynq, nrr_q, irvec_q, wslen_q)
    !
    !
    ! Electron-Phonon vertex (Bloch el and Bloch ph -> Wannier el and Bloch ph)
    !
    ! Open the prefix.epmatwe file
    IF (ionode) THEN
      lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint   = TRIM(prefix)//'.epmatwe'
      CALL diropn(iunepmatwe, 'epmatwe', lrepmatw, exst)
    ENDIF
    !
    DO iq = 1, nqc
      !
      xxq = xqc(:, iq)
      !
      ! we need the cu again for the k+q points, we generate the map here
      !
      ALLOCATE(xkq(3, nkstot), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating xkq', 1)
      CALL loadumat(nbndep, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband, w_centers)
      DEALLOCATE(xkq, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating xkq', 1)
      !
      DO imode = 1, nmodes
        !
        CALL ephbloch2wane(nbndep, nbndsub, nks, nkstot, xk_loc, cu, cuq, &
          epmatq(:, :, :, imode, iq), nrr_k, irvec_k, wslen_k, epmatwe_mem(:, :, :, imode))
        !
      ENDDO
      ! Only the master node writes
      IF (ionode) THEN
        ! direct write of epmatwe for this iq
        CALL rwepmatw(epmatwe_mem, nbndsub, nrr_k, nmodes, iq, iunepmatwe, +1)
        !
      ENDIF
      !
    ENDDO
    !
    IF (ionode) CLOSE(iunepmatwe, STATUS = 'keep')
    !
    ! Electron-Phonon vertex (Wannier el and Bloch ph -> Wannier el and Wannier ph)
    CALL ephbloch2wanp_mem(nbndsub, nmodes, xqc, nqc, irvec_k, irvec_g, nrr_k, nrr_g)
    !
    IF (epwwrite) THEN
      CALL epw_write(nrr_k, nrr_q, nrr_g, w_centers)
    ENDIF
    !
    DEALLOCATE(epmatq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating epmatq', 1)
    DEALLOCATE(dynq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating dynq', 1)
    IF (vme == 'dipole') DEALLOCATE(dmec)
    DEALLOCATE(epmatwe_mem, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating epmatwe_mem', 1)
    DEALLOCATE(cu, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cu', 1)
    DEALLOCATE(cuq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cuq', 1)
    DEALLOCATE(lwin, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating lwin', 1)
    DEALLOCATE(lwinq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating lwinq', 1)
    DEALLOCATE(exband, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating exband', 1)
  ENDIF ! (epwread .AND. .NOT. epbread)
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
  ! At this point, we will interpolate the Wannier rep to the Bloch rep
  ! for electrons, phonons and the ep-matrix
  !
  ! Load the fine-grid q and k grids.
  ! nkqtotf is computed inside
  CALL loadqmesh_serial
  CALL loadkmesh_para
  ! Defines the total number of k-points
  nktotf = nkqtotf / 2
  !
  ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epmatwef', 1)
  ALLOCATE(wf(nmodes, nqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating wf', 1)
  ALLOCATE(etf(nbndsub, nkqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating etf', 1)
  ALLOCATE(etf_ks(nbndsub, nkqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating etf_ks', 1)
  ALLOCATE(epmatf(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epmatf', 1)
  ALLOCATE(cufkk(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cufkk', 1)
  ALLOCATE(cufkq(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cufkq', 1)
  ALLOCATE(uf(nmodes, nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating uf', 1)
  ALLOCATE(bmatf(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating bmatf', 1)
  ALLOCATE(eps_rpa(nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating eps_rpa', 1)
  ALLOCATE(isk_dummy(nkqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating isk_dummy', 1)
  epmatwef(:, :, :) = czero
  wf(:, :)          = zero
  etf(:, :)         = zero
  etf_ks(:, :)      = zero
  epmatf(:, :)      = czero
  cufkk(:, :)       = czero
  cufkq(:, :)       = czero
  uf(:, :)          = czero
  bmatf(:, :)       = czero
  eps_rpa(:)        = czero
  isk_dummy(:)      = 0
  !
  ! Allocate velocity and dipole matrix elements after getting grid size
  !
  IF (vme == 'wannier') THEN
    ALLOCATE(vmef(3, nbndsub, nbndsub, 2 * nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating vmef', 1)
    vmef(:, :, :, :) = czero
  ELSE
    ALLOCATE(dmef(3, nbndsub, nbndsub, 2 * nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating dmef', 1)
    dmef(:, :, :, :) = czero
  ENDIF
  !
  ALLOCATE(cfac(nrr_k, dims, dims), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cfac', 1)
  ALLOCATE(cfacq(nrr_k, dims, dims), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating cfacq', 1)
  ALLOCATE(rdotk(nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating rdotk', 1)
  ALLOCATE(rdotk2(nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating rdotk2', 1)
  ! This is simply because dgemv take only real number (not integer)
  ALLOCATE(irvec_r(3, nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating irvec_r', 1)
  irvec_r = REAL(irvec_k, KIND = DP)
  !
  ! Zeroing everything - initialization is important !
  cfac(:, :, :)  = czero
  cfacq(:, :, :) = czero
  rdotk(:)       = zero
  rdotk2(:)      = zero
  !
  ! ------------------------------------------------------
  ! Hamiltonian : Wannier -> Bloch (preliminary)
  ! ------------------------------------------------------
  !
  ! We here perform a preliminary interpolation of the hamiltonian
  ! in order to determine the fermi window ibndmin:ibndmax for later use.
  ! We will interpolate again afterwards, for each k and k+q separately
  !
  xxq = 0.d0
  !
  ! nkqf is the number of kpoints in the pool
  DO ik = 1, nkqf
    !
    xxk = xkf(:, ik)
    !
    IF (2 * (ik / 2) == ik) THEN
      !
      !  this is a k+q point : redefine as xkf (:, ik-1) + xxq
      !
      CALL cryst_to_cart(1, xxq, at, -1)
      xxk = xkf(:, ik - 1) + xxq
      CALL cryst_to_cart(1, xxq, bg, 1)
      !
    ENDIF
    !
    ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
    ! + optimize the 2\pi r\cdot k with Blas
    CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xxk, 1, 0.0_DP, rdotk, 1)
    !
    DO iw = 1, dims
      DO iw2 = 1, dims
        DO ir = 1, nrr_k
          IF (ndegen_k(ir, iw2, iw) > 0) cfac(ir, iw2, iw) = EXP(ci * rdotk(ir)) / ndegen_k(ir, iw2, iw)
        ENDDO
      ENDDO
    ENDDO
    !
    CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ik), chw, cfac, dims)
  ENDDO
  !
  WRITE(stdout,'(/5x,a,f10.6,a)') 'Fermi energy coarse grid = ', ef * ryd2ev, ' eV'
  !
  IF (efermi_read) THEN
    !
    ef = fermi_energy
    WRITE(stdout,'(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(/5x,a,f10.6,a)') 'Fermi energy is read from the input file: Ef = ', ef * ryd2ev, ' eV'
    WRITE(stdout,'(/5x,a)') REPEAT('=',67)
    !
    ! SP: even when reading from input the number of electron needs to be correct
    already_skipped = .FALSE.
    IF (nbndskip > 0) THEN
      IF (.NOT. already_skipped) THEN
        IF (noncolin) THEN
          nelec = nelec - one * nbndskip
        ELSE
          nelec = nelec - two * nbndskip
        ENDIF
        already_skipped = .TRUE.
        WRITE(stdout, '(/5x,"Skipping the first ", i4, " bands:")') nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ", f9.5, " electrons")') nelec
      ENDIF
    ENDIF
    !
  ELSEIF (band_plot) THEN
    !
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(/5x,"Fermi energy corresponds to the coarse k-mesh")')
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    !
  ELSE
    ! here we take into account that we may skip bands when we wannierize
    ! (spin-unpolarized)
    ! RM - add the noncolin case
    already_skipped = .FALSE.
    IF (nbndskip > 0) THEN
      IF (.NOT. already_skipped) THEN
        IF (noncolin) THEN
          nelec = nelec - one * nbndskip
        ELSE
          nelec = nelec - two * nbndskip
        ENDIF
        already_skipped = .TRUE.
        WRITE(stdout, '(/5x,"Skipping the first ", i4, " bands:")') nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ", f9.5, " electrons")') nelec
      ENDIF
    ENDIF
    !
    ! Fermi energy
    !
    ! Since wkf(:,ikq) = 0 these bands do not bring any contribution to Fermi level
    IF (ABS(degaussw) < eps16) THEN
      ! Use 1 meV instead
      efnew = efermig(etf, nbndsub, nkqf, nelec, wkf, 1.0d0 / ryd2mev, ngaussw, 0, isk_dummy)
    ELSE
      efnew = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
    ENDIF
    !
    WRITE(stdout, '(/5x,a,f10.6,a)') &
        'Fermi energy is calculated from the fine k-mesh: Ef = ', efnew * ryd2ev, ' eV'
    !
    ! if 'fine' Fermi level differs by more than 250 meV, there is probably something wrong
    ! with the wannier functions, or 'coarse' Fermi level is inaccurate
    IF (ABS(efnew - ef) * ryd2eV > 0.250d0 .AND. (.NOT. eig_read)) &
      WRITE(stdout,'(/5x,a)') 'Warning: check if difference with Fermi level fine grid makes sense'
    WRITE(stdout,'(/5x,a)') REPEAT('=',67)
    !
    ef = efnew
  ENDIF
  ! ------------------------------------------------------------
  ! Apply a possible shift to eigenenergies (applied later)
  icbm = 1
  IF (ABS(scissor) > eps6) THEN
    IF (assume_metal) THEN
      CALL errore("ephwann_shuffle", "A scissor shift is applied but the material is a metal...", 1)
    ENDIF
    IF (noncolin) THEN
      icbm = FLOOR(nelec / 1.0d0) + 1
    ELSE
      icbm = FLOOR(nelec / 2.0d0) + 1
    ENDIF
    etf(icbm:nbndsub, :) = etf(icbm:nbndsub, :) + scissor
    !
    WRITE(stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the CB ",i6)' ) scissor * ryd2ev, icbm
  ENDIF
  !
  ! Identify the bands within fsthick from the Fermi level
  ! Return ibndmin and ibndmax
  CALL fermiwindow()
  nbndfst = ibndmax - ibndmin + 1
  !
  ! Define it only once for the full run.
  CALL fkbounds(nktotf, lower_bnd, upper_bnd)
  !
  ! Re-order the k-point according to weather they are in or out of the fshick
  ! windows
  IF ((iterative_bte .OR. ephwrite) .AND. mp_mesh_k) THEN
    CALL load_rebal
  ENDIF
  !
  ! In the case of crystal ASR
  IF (lifc) THEN
    !
    ! build the WS cell corresponding to the force constant grid
    !
    atws(:, 1) = at(:, 1) * DBLE(nqc1)
    atws(:, 2) = at(:, 2) * DBLE(nqc2)
    atws(:, 3) = at(:, 3) * DBLE(nqc3)
    rws(:, :)  = zero
    nrws       = 0
    ! initialize WS r-vectors
    CALL wsinit(rws, nrwsx, nrws, atws)
  ELSE
    atws(:, :) = zero
    rws(:, :)  = zero
    nrws       = 0
  ENDIF
  !
  ! Open the ephmatwp file here
#if defined(__MPI)
  ! Check for directory given by "outdir"
  !
  filint = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
  CALL MPI_FILE_OPEN(world_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatwp2, ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'error in MPI_FILE_OPEN', 1)
#endif
  !
  ! get the size of the matrix elements stored in each pool
  ! for informational purposes.  Not necessary
  !
  CALL mem_size(nmodes, nkf)
  !
  ALLOCATE(etf_all(nbndfst, nktotf), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating etf_all', 1)
  etf_all(:, :) = zero
  !
  ! ------------------------------------------------
  ! The IBTE implement works in two steps
  ! 1) compute the dominant scattering rates and store them to file
  ! 2) read them from file and solve the IBTE where all important element are in memory
  ! ------------------------------------------------
  !
  ! Initialization and restart when doing IBTE
  IF (iterative_bte) THEN
    ALLOCATE(inv_tau_all(nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating inv_tau_all', 1)
    ALLOCATE(inv_tau_allcb(nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating inv_tau_allcb', 1)
    inv_tau_all(:, :, :)   = zero
    inv_tau_allcb(:, :, :) = zero
    lrepmatw2_restart(:)   = 0
    lrepmatw5_restart(:)   = 0
    ! We save matrix elements that are smaller than machine precision (1d-16).
    ! The sum of all the elements must be smaller than that
    ! nkf1 * nkf2 * nkf3 * nqf1 * nqf2 * nqf3 * (nbndfst) * (nbndfst)
    ! must be smaller than 1d-16
    ! To avoid overflow we need to use DP
    threshold = 1d-16 / (INT(nkf1, KIND = 8) * INT(nkf2, KIND = 8) * INT(nkf3, KIND = 8) * &
                         INT(nqf1, KIND = 8) * INT(nqf2, KIND = 8) * INT(nqf3, KIND = 8) * &
                         INT((nbndfst), KIND = 8) * INT((nbndfst), KIND = 8))
  ENDIF
  !
  IF (iterative_bte .AND. epmatkqread) THEN
    ALLOCATE(vkk_all(3, nbndfst, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating vkk_all', 1)
    ALLOCATE(wkf_all(nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating wkf_all', 1)
    !
    CALL iter_restart(etf_all, wkf_all, vkk_all, ind_tot, ind_totcb, ef0, efcb)
    !
    DEALLOCATE(vkk_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating vkk_all', 1)
    DEALLOCATE(wkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wkf_all', 1)
    DEALLOCATE(inv_tau_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating inv_tau_all', 1)
    DEALLOCATE(inv_tau_allcb, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating inv_tau_allcb', 1)
    !
  ELSE ! (iterative_bte .AND. epmatkqread)
    IF (iterative_bte) THEN
      ! Open the required files
      CALL iter_open(ind_tot, ind_totcb, lrepmatw2_restart, lrepmatw5_restart)
    ENDIF
    !
    IF (lifc) THEN
      ALLOCATE(wscache(-2 * nqc3:2 * nqc3, -2 * nqc2:2 * nqc2, -2 * nqc1:2 * nqc1, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating wscache', 1)
      wscache(:, :, :, :, :) = zero
    ENDIF
    !
    ! -----------------------------------------------------------------------
    ! Determines which q-points falls within the fsthick windows
    ! Store the result in the selecq.fmt file
    ! If the file exists, automatically restart from the file
    ! -----------------------------------------------------------------------
    !
    ! Check if the grids are homogeneous and commensurate
    homogeneous = .FALSE.
    IF ((nkf1 /= 0) .AND. (nkf2 /= 0) .AND. (nkf3 /= 0) .AND. &
        (nqf1 /= 0) .AND. (nqf2 /= 0) .AND. (nqf3 /= 0) .AND. &
        (MOD(nkf1, nqf1) == 0) .AND. (MOD(nkf2, nqf2) == 0) .AND. (MOD(nkf3, nqf3) == 0)) THEN
      homogeneous = .TRUE.
    ELSE
      homogeneous = .FALSE.
    ENDIF
    !
    totq = 0
    !
    ! HP: This is implemented for the Superconductivity
    !IF (ephwrite) THEN
    !  !
    !  totq = nqf
    !  ALLOCATE(selecq(nqf), STAT = ierr)
    !  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating selecq', 1)
    !  DO iq = 1, nqf
    !    selecq(iq) = iq
    !  ENDDO
    !  !
    !ELSE ! ephwrite
      ! Check if the file has been pre-computed
      IF (mpime == ionode_id) THEN
        INQUIRE(FILE = 'selecq.fmt', EXIST = exst)
      ENDIF
      CALL mp_bcast(exst, ionode_id, world_comm)
      !
      IF (exst) THEN
        IF (selecqread) THEN
          WRITE(stdout, '(5x,a)')' '
          WRITE(stdout, '(5x,a)')'Reading selecq.fmt file. '
          CALL qwindow(exst, nrr_k, dims, totq, selecq, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
        ELSE
          WRITE(stdout, '(5x,a)')' '
          WRITE(stdout, '(5x,a)')'A selecq.fmt file was found but re-created because selecqread == .FALSE. '
          CALL qwindow(.FALSE., nrr_k, dims, totq, selecq, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
        ENDIF
      ELSE ! exst
        IF (selecqread) THEN
          CALL errore( 'ephwann_shuffle', 'Variable selecqread == .TRUE. but file selecq.fmt not found.',1 )
        ELSE
          CALL qwindow(exst, nrr_k, dims, totq, selecq, irvec_r, ndegen_k, cufkk, cufkq, homogeneous)
        ENDIF
      ENDIF
      !
      WRITE(stdout, '(5x,a,i8,a)')'We only need to compute ', totq, ' q-points'
      WRITE(stdout, '(5x,a)')' '
      !
    !ENDIF ! ephwrite
    !
    ! -----------------------------------------------------------------------
    ! Possible restart during step 1)
    ! -----------------------------------------------------------------------
    iq_restart = 1
    first_cycle = .FALSE.
    first_time = .TRUE.
    !
    ! Fine mesh set of g-matrices.  It is large for memory storage
    ALLOCATE(epf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epf17', 1)
    ALLOCATE(eptmp(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating eptmp', 1)
    ALLOCATE(epmatlrT(nbndsub, nbndsub, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epmatlrT', 1)
    IF (phonselfen) THEN
      ALLOCATE(lambda_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating lambda_all', 1)
      ALLOCATE(lambda_v_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating lambda_v_all', 1)
      ALLOCATE(gamma_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating gamma_all', 1)
      ALLOCATE(gamma_v_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating gamma_v_all', 1)
      lambda_all(:, :, :, :)   = zero
      lambda_v_all(:, :, :, :) = zero
      gamma_all(:, :, :, :)    = zero
      gamma_v_all(:, :, :, :)  = zero
    ENDIF
    IF (specfun_el .OR. specfun_pl) THEN
      ALLOCATE(esigmar_all(nbndfst, nktotf, nw_specfun, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating esigmar_all', 1)
      ALLOCATE(esigmai_all(nbndfst, nktotf, nw_specfun, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating esigmai_all', 1)
      ALLOCATE(a_all(nw_specfun, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating a_all', 1)
      esigmar_all(:, :, :, :) = zero
      esigmai_all(:, :, :, :) = zero
      a_all(:, :, :) = zero
    ENDIF
    IF (specfun_ph) THEN
      ALLOCATE(a_all_ph(nw_specfun, totq, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating a_all_ph', 1)
      a_all_ph(:, :, :) = zero
    ENDIF
    IF (scattering .AND. .NOT. iterative_bte) THEN
      ALLOCATE(inv_tau_all(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating inv_tau_all', 1)
      ALLOCATE(zi_allvb(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating zi_allvb', 1)
      inv_tau_all(:, :, :) = zero
      zi_allvb(:, :, :)    = zero
    ENDIF
    IF (int_mob .AND. carrier) THEN
      ALLOCATE(inv_tau_allcb(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating inv_tau_allcb', 1)
      ALLOCATE(zi_allcb(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating zi_allcb', 1)
      inv_tau_allcb(:, :, :) = zero
      zi_allcb(:, :, :)      = zero
    ENDIF
    IF (elecselfen .OR. plselfen) THEN
      ALLOCATE(sigmar_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating sigmar_all', 1)
      ALLOCATE(sigmai_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating sigmai_all', 1)
      ALLOCATE(zi_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating zi_all', 1)
      sigmar_all(:, :, :) = zero
      sigmai_all(:, :, :) = zero
      zi_all(:, :, :)     = zero
      IF (iverbosity == 3) THEN
        ALLOCATE(sigmai_mode(nbndfst, nmodes, nktotf, nstemp), STAT = ierr)
        IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating sigmai_mode', 1)
        sigmai_mode(:, :, :, :) = zero
      ENDIF
    ENDIF ! elecselfen
    IF (lindabs) THEN
      ! Calculate the number of frequency points
      nomega = INT((omegamax - omegamin) / omegastep) + 1
      ALLOCATE(epsilon2_abs(3, nomega, neta, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epsilon2_abs', 1)
      ALLOCATE(epsilon2_abs_lorenz(3, nomega, neta, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating epsilon2_abs_lorenz', 1)
      ALLOCATE(epsilon2_abs_all(3, nomega, neta, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error allocating epsilon2_abs_all', 1)
      ALLOCATE(epsilon2_abs_lorenz_all(3, nomega, neta, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shiffle', 'Error allocating epsilon2_abs_lorenz_all', 1)
    ENDIF ! indabs
    !
    ! Restart in SERTA case or self-energy (electron or plasmon) case
    IF (restart) THEN
      IF (elecselfen .OR. plselfen) THEN
        CALL selfen_el_read(iq_restart, totq, nktotf, sigmar_all, sigmai_all, zi_all)
      ENDIF
      IF (specfun_el .OR. specfun_pl) THEN
        CALL spectral_read(iq_restart, totq, nktotf, esigmar_all, esigmai_all)
      ENDIF
      IF (scattering) THEN
        IF (int_mob .AND. carrier) THEN
          ! Here inv_tau_all and inv_tau_allcb gets updated
          CALL tau_read(iq_restart, totq, nktotf, .TRUE.)
        ELSE
          ! Here inv_tau_all gets updated
          CALL tau_read(iq_restart, totq, nktotf, .FALSE.)
        ENDIF
      ENDIF
      !
      ! Potential restart in indirect optics
      IF (lindabs) THEN
        CALL indabs_read(iq_restart, totq, epsilon2_abs_all, epsilon2_abs_lorenz_all)
      ENDIF
      !
      ! If you restart from reading a file. This prevent
      ! the case were you restart but the file does not exist
      IF (iq_restart > 1) first_cycle = .TRUE.
      !
    ENDIF ! restart
    !
    ! Scatread assumes that you alread have done the full q-integration
    ! We just do one loop to get interpolated eigenenergies.
    IF(scatread) iq_restart = totq - 1
    !
    ! Restart in IBTE and Superconductivity cases
    IF (iterative_bte .OR. ephwrite) THEN
      IF (mpime == ionode_id) THEN
        INQUIRE(FILE = 'restart.fmt', EXIST = exst)
      ENDIF
      CALL mp_bcast(exst, ionode_id, world_comm)
      !
      IF (exst) THEN
        IF (mpime == ionode_id) THEN
          OPEN(UNIT = iunrestart, FILE = 'restart.fmt', STATUS = 'old', IOSTAT = ios)
          READ(iunrestart, *) iq_restart
          READ(iunrestart, *) ind_tot
          READ(iunrestart, *) ind_totcb
          READ(iunrestart, *) npool_tmp
          DO ipool = 1, npool
            READ(iunrestart, *) lrepmatw2_restart(ipool)
          ENDDO
          DO ipool = 1, npool
            READ(iunrestart, *) lrepmatw5_restart(ipool)
          ENDDO
          CLOSE(iunrestart)
        ENDIF
        CALL mp_bcast(iq_restart, ionode_id, world_comm)
        CALL mp_bcast(npool_tmp, ionode_id, world_comm)
        CALL mp_bcast(lrepmatw2_restart, ionode_id, world_comm)
        CALL mp_bcast(lrepmatw5_restart, ionode_id, world_comm)
        IF (npool /= npool_tmp) CALL errore('ephwann_shuffle','Number of cores is different',1)
        !
        IF (iterative_bte) THEN
          IF (mpime == ionode_id) THEN
            OPEN(UNIT = iuntau, FORM = 'unformatted', FILE = 'inv_tau_tmp', STATUS = 'old')
            READ(iuntau) inv_tau_all
            CLOSE(iuntau)
            !
            OPEN(UNIT = iuntaucb, FORM = 'unformatted', FILE = 'inv_taucb_tmp', STATUS = 'old')
            READ(iuntaucb) inv_tau_allcb
            CLOSE(iuntaucb)
          ENDIF
          CALL mp_bcast(inv_tau_all, ionode_id, world_comm)
          CALL mp_bcast(inv_tau_allcb, ionode_id, world_comm)
          IF (lower_bnd - 1 >= 1) THEN
            inv_tau_all(:, 1:lower_bnd - 1, :) = 0d0
            inv_tau_allcb(:, 1:lower_bnd - 1, :) = 0d0
          ENDIF
          IF (upper_bnd + 1 <= nktotf) THEN
            inv_tau_all(:, upper_bnd + 1:nktotf, :) = 0d0
            inv_tau_allcb(:, upper_bnd + 1:nktotf, :) = 0d0
          ENDIF
        ENDIF
        !
#if defined(__MPI)
        CALL MPI_BCAST(ind_tot,   1, MPI_OFFSET, ionode_id, world_comm, ierr)
        CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
#endif
        IF (ierr /= 0) CALL errore('ephwann_shuffle', 'error in MPI_BCAST', 1)
        !
        IF(ephwrite .AND. iq_restart > 1) THEN
          first_cycle = .TRUE.
          CALL check_restart_ephwrite
        ENDIF
        !
        ! Now, the iq_restart point has been done, so we need to do the next
        iq_restart = iq_restart + 1
        !
        IF (iq_restart < totq) THEN
          WRITE(stdout, '(5x,a,i8,a)')'We restart from ', iq_restart, ' q-points'
        ELSE
          WRITE(stdout, '(5x,a)')'All q-points are done, no need to restart !!'
        ENDIF
        !
      ENDIF ! exst
    ENDIF
    ! Adaptative smearing when degauss = 0
    adapt_smearing = .FALSE.
    IF (ABS(degaussw) < eps16) THEN
      ALLOCATE(eta(nmodes, nbndfst, nkf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating eta', 1)
      ALLOCATE(vmefp(3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating vmefp', 1)
      eta(:, :, :)   = zero
      vmefp(:, :, :) = czero
      adapt_smearing = .TRUE.
    ENDIF
    !
    DO iqq = iq_restart, totq
      CALL start_clock ('ep-interp')
      !
      epf17(:, :, :, :) = czero
      eptmp(:, :, :, :)    = czero
      epmatlrT(:, :, :, :) = czero
      cufkk(:, :) = czero
      cufkq(:, :) = czero
      !
      iq = selecq(iqq)
      !
      !
      ! In case of big calculation, show progression of iq (especially usefull when
      ! elecselfen = true as nothing happen during the calculation otherwise.
      !
      IF (.NOT. phonselfen) THEN
        IF (MOD(iqq, restart_step) == 0) THEN
          WRITE(stdout, '(5x,a,i10,a,i10)' ) 'Progression iq (fine) = ', iqq, '/', totq
        ENDIF
      ENDIF
      !
      xxq = xqf(:, iq)
      ! Temporarily commented by H. Lee
!      CALL find_gmin(xxq)
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
          wf(nu, iq) =  DSQRT(ABS(w2(nu)))
        ELSE
          wf(nu, iq) = -DSQRT(ABS(w2(nu)))
        ENDIF
        !
        DO mu = 1, nmodes
          na = (mu - 1) / 3 + 1
          uf(mu, nu) = uf(mu, nu) / DSQRT(amass(ityp(na)))
        ENDDO
      ENDDO
      !
      ! --------------------------------------------------------------
      ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
      ! --------------------------------------------------------------
      !
      DO imode = 1, nmodes
        epmatwef(:, :, :) = czero
        !
        IF (.NOT. longrange) THEN
          CALL ephwan2blochp_mem(imode, nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatwef, nbndsub, nrr_k, dims, dims2)
        ENDIF
        !
        !
        !  number of k points with a band on the Fermi surface
        fermicount = 0
        !
        IF (lscreen) THEN
          IF (scr_typ == 0) CALL rpa_epsilon (xxq, wf(imode, iq), nmodes, epsi, eps_rpa)
          IF (scr_typ == 1) CALL tf_epsilon (xxq, nmodes, epsi, eps_rpa)
        ENDIF
        !
        ! -------------------------------------------------------------
        ! Adaptative smearing: computes phonon velocity
        ! -------------------------------------------------------------
        ! Needs to be adapted to work per mode
        !IF (adapt_smearing) THEN
        !  CALL vmewan2blochp(xxq, nmodes, nrr_q, irvec_q, ndegen_q, uf, vmefp(:, :, :), wf(imode, iq), rws, nrws)
        !ENDIF
        !
        ! this is a loop over k blocks in the pool
        ! (size of the local k-set)
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
          CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
          CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq2, 1, 0.0_DP, rdotk2, 1)
          !
          IF (use_ws) THEN
            DO iw = 1, dims
              DO iw2 = 1, dims
                DO ir = 1, nrr_k
                  IF (ndegen_k(ir, iw2, iw) > 0) THEN
                    cfac(ir, iw2, iw)  = EXP(ci * rdotk(ir))  / ndegen_k(ir, iw2, iw)
                    cfacq(ir, iw2, iw) = EXP(ci * rdotk2(ir)) / ndegen_k(ir, iw2, iw)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ELSE
            cfac(:, 1, 1)  = EXP(ci * rdotk(:))  / ndegen_k(:, 1, 1)
            cfacq(:, 1, 1) = EXP(ci * rdotk2(:)) / ndegen_k(:, 1, 1)
          ENDIF
          !
          ! ------------------------------------------------------
          ! hamiltonian : Wannier -> Bloch
          ! ------------------------------------------------------
          !
          ! Kohn-Sham first, then get the rotation matricies for following interp.
          IF (eig_read) THEN
             CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_ks(:, ikk), chw_ks, cfac, dims)
             CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_ks(:, ikq), chw_ks, cfacq, dims)
          ENDIF
          !
          CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ikk), chw, cfac, dims)
          CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf(:, ikq), chw, cfacq, dims)
          !
          ! Apply a possible scissor shift
          etf(icbm:nbndsub, ikk) = etf(icbm:nbndsub, ikk) + scissor
          etf(icbm:nbndsub, ikq) = etf(icbm:nbndsub, ikq) + scissor
          !
          IF (.NOT. scattering) THEN
            IF (vme == 'wannier') THEN
              !
              ! ------------------------------------------------------
              !  velocity: Wannier -> Bloch
              ! ------------------------------------------------------
              !
              IF (eig_read) THEN
                ! Renormalize the eigenvalues and vmef with the read eigenvalues
                CALL renorm_eig(ikk, ikq, nrr_k, dims, ndegen_k, irvec_k, irvec_r, cufkk, cufkq, cfac, cfacq)
                !
              ELSE ! eig_read
                CALL vmewan2bloch &
                     (nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), chw, cfac, dims)
                CALL vmewan2bloch &
                     (nbndsub, nrr_k, irvec_k, cufkq, vmef(:, :, :, ikq), etf(:, ikq), etf_ks(:, ikq), chw, cfacq, dims)
              ENDIF
            ELSE
              !
              ! ------------------------------------------------------
              !  dipole: Wannier -> Bloch
              ! ------------------------------------------------------
              !
              CALL dmewan2bloch(nbndsub, nrr_k, cufkk, dmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), cfac, dims)
              CALL dmewan2bloch(nbndsub, nrr_k, cufkq, dmef(:, :, :, ikq), etf(:, ikq), etf_ks(:, ikq), cfacq, dims)
              !
            ENDIF
          ENDIF ! not scattering
          !
          IF (.NOT. scatread) THEN
            ! interpolate only when (k,k+q) both have at least one band
            ! within a Fermi shell of size fsthick
            !
            IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. &
                (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
              !
              ! Compute velocities
              !
              IF (scattering) THEN
                IF (vme == 'wannier') THEN
                  CALL vmewan2bloch &
                       (nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), chw, cfac, dims)
                  CALL vmewan2bloch &
                       (nbndsub, nrr_k, irvec_k, cufkq, vmef(:, :, :, ikq), etf(:, ikq), etf_ks(:, ikq), chw, cfacq, dims)
                ELSE
                  CALL dmewan2bloch(nbndsub, nrr_k, cufkk, dmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), cfac, dims)
                  CALL dmewan2bloch(nbndsub, nrr_k, cufkq, dmef(:, :, :, ikq), etf(:, ikq), etf_ks(:, ikq), cfacq, dims)
                ENDIF
              ENDIF
              !
              ! --------------------------------------------------------------
              ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
              ! --------------------------------------------------------------
              !
              ! SP: Note: In case of polar materials, computing the long-range and short-range term
              !     separately might help speed up the convergence. Indeed the long-range term should be
              !     much faster to compute. Note however that the short-range term still contains a linear
              !     long-range part and therefore could still be a bit more difficult to converge than
              !     non-polar materials.
              !
              IF (longrange) THEN
                !
                epmatf(:, :) = czero
                !
              ELSE
                !
                epmatf(:, :) = czero
                CALL ephwan2bloch_mem(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, cfac, dims)
                !
              ENDIF
              !
              IF (lpolar) THEN
                !
                CALL compute_umn_f(nbndsub, cufkk, cufkq, bmatf)
                !
                IF ((ABS(xxq(1)) > eps8) .OR. (ABS(xxq(2)) > eps8) .OR. (ABS(xxq(3)) > eps8)) THEN
                  !
                  CALL cryst_to_cart(1, xxq, bg, 1)
                  CALL rgd_blk_epw_fine_mem(imode, nqc1, nqc2, nqc3, xxq, uf, epmatlrT(:, :, imode, ik), &
                                        nmodes, epsi, zstar, bmatf, one)
                  CALL cryst_to_cart(1, xxq, at, -1)
                  !
                ENDIF
                !
              ENDIF
              !
              ! Store epmatf in memory
              !
              DO jbnd=ibndmin, ibndmax
                DO ibnd=ibndmin, ibndmax
                  !
                  IF (lscreen) THEN
                    eptmp(ibnd - ibndmin + 1, jbnd - ibndmin + 1, imode ,ik) = epmatf(ibnd, jbnd) / eps_rpa(imode)
                  ELSE
                    eptmp(ibnd - ibndmin + 1, jbnd - ibndmin + 1, imode, ik) = epmatf(ibnd, jbnd)
                  ENDIF
                  !
                ENDDO
              ENDDO
              !
            ENDIF
          ENDIF ! scatread
        ENDDO  ! end loop over k points
      ENDDO ! modes
      !
      ! Now do the eigenvector rotation:
      ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
      !
      DO ik = 1, nkf
        CALL ZGEMM('n', 'n', nbndfst * nbndfst, nmodes, nmodes, cone, eptmp(:, :, :, ik), &
                   nbndfst * nbndfst, uf, nmodes, czero, epf17(:, :, :, ik), nbndfst * nbndfst)
        !
      ENDDO
      !
      ! After the rotation, add the long-range that is already rotated
      DO jbnd = ibndmin, ibndmax
        DO ibnd = ibndmin, ibndmax
          epf17(ibnd - ibndmin+1, jbnd - ibndmin + 1, :, :) = epf17(ibnd-ibndmin+1, jbnd - ibndmin + 1, :, :) &
                                                            + epmatlrT(ibnd, jbnd, :, :)
        ENDDO
      ENDDO
      !
      !
      !
      IF (prtgkk    ) CALL print_gkk(iq)
      IF (phonselfen) CALL selfen_phon_q(iqq, iq, totq)
      IF (elecselfen) CALL selfen_elec_q(iqq, iq, totq, first_cycle)
      IF (plselfen .AND. vme == 'dipole') CALL selfen_pl_q(iqq, iq, totq, first_cycle)
      IF (nest_fn   ) CALL nesting_fn_q(iqq, iq)
      IF (specfun_el) CALL spectral_func_el_q(iqq, iq, totq, first_cycle)
      IF (specfun_ph) CALL spectral_func_ph_q(iqq, iq, totq)
      IF (specfun_pl .AND. vme == 'dipole') CALL spectral_func_pl_q(iqq, iq, totq, first_cycle)
      IF (ephwrite) THEN
        IF (first_cycle .OR. iq == 1) THEN
          CALL kmesh_fine
          CALL kqmap_fine
          CALL count_kpoints
          first_cycle = .FALSE.
        ENDIF
        CALL write_ephmat(iqq, iq, totq)
      ENDIF
      !
      IF (.NOT. scatread) THEN
        !
        ! Indirect absorption ---------------------------------------------------------
        ! If Indirect absortpion, keep unshifted values:
        IF (lindabs .AND. .NOT. scattering) THEN
           etf_ks(:, :) = etf(:, :)
           ! We remove the scissor
           etf_ks(icbm:nbndsub, :) = etf_ks(icbm:nbndsub, :) - scissor
        ENDIF
        !
        ! Indirect absorption
        IF (lindabs .AND. .NOT. scattering) THEN
          IF (carrier .and. (iq == iq_restart)) THEN
            ALLOCATE(ef0_fca(nstemp), STAT = ierr)
            IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error allocating ef0_fca', 1)
            DO itemp = 1, nstemp
              etemp_fca = gtemp(itemp)
              CALL fermi_carrier_indabs(itemp, etemp_fca, ef0_fca)
            ENDDO
          ENDIF
          CALL indabs_main(iq, totq, first_cycle, iq_restart)
        ENDIF
        !
        ! Conductivity ---------------------------------------------------------
        IF (scattering) THEN
          !
          ! If we want to compute intrinsic mobilities, call fermicarrier to  correctly positionned the ef0 level.
          ! This is only done once for the first iq. Also compute the dos at the same time
          IF (iqq == iq_restart) THEN
            DO itemp = 1, nstemp
              etemp = gtemp(itemp)
              CALL fermicarrier(itemp, etemp, ef0, efcb, ctype)
              ! compute dos for metals
              IF (assume_metal) THEN
                CALL compute_dos(itemp, ef0, dos)
              ENDIF
            ENDDO
          ENDIF ! iqq=0
          !
          IF (.NOT. iterative_bte) THEN
            CALL scattering_rate_q(iqq, iq, totq, ef0, efcb, first_cycle)
            ! Computes the SERTA mobility
            IF (iqq == totq) CALL transport_coeffs(ef0, efcb)
          ENDIF
          !
          IF (iterative_bte) THEN
            CALL start_clock('print_ibte')
            CALL print_ibte(iqq, iq, totq, xxq, ef0, efcb, first_cycle, ind_tot, ind_totcb, &
                            lrepmatw2_restart, lrepmatw5_restart, ctype)
            CALL stop_clock('print_ibte')
            !
            ! Finished, now compute SERTA and IBTE mobilities
            IF (iqq == totq) THEN
              WRITE(stdout, '(5x,a)')' '
              WRITE(stdout, '(5x,"epmatkqread automatically changed to .TRUE. as all scattering have been computed.")')
              WRITE(stdout, '(5x,a)')' '
              ! Close files
              CALL iter_close()
              ! Merge files
              CALL iter_merge
              !
            ENDIF
          ENDIF
          !
        ENDIF ! scattering
        ! --------------------------------------
        !
        CALL stop_clock('ep-interp')
        !
      ENDIF ! scatread
    ENDDO  ! end loop over q points
    !
    ! Check Memory usage
    CALL system_mem_usage(valueRSS)
    !
    WRITE(stdout, '(a)' )             '     ==================================================================='
    WRITE(stdout, '(a,i10,a)' ) '     Memory usage:  VmHWM =',valueRSS(2)/1024,'Mb'
    WRITE(stdout, '(a,i10,a)' ) '                   VmPeak =',valueRSS(1)/1024,'Mb'
    WRITE(stdout, '(a)' )             '     ==================================================================='
    WRITE(stdout, '(a)' )
    !
    ! ---------------------------------------------------------------------------------------
    ! ---------------------------------------------------------------------------------------
    !
    ! SP: Added lambda and phonon lifetime writing to file.
    !
    IF (mpime == ionode_id) THEN
      !
      IF (phonselfen) THEN
        DO itempphen = 1, nstemp
          WRITE(tp, "(f8.3)") gtemp(itempphen) * ryd2ev / kelvin2eV
          filephselfen = 'lambda.phself.' // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = lambda_phself, FILE = filephselfen)
          WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
          WRITE(lambda_phself, *) '#Modes     ',(imode, imode = 1, nmodes)
          DO iqq = 1, nqtotf
              !
              !myfmt = "(*(3x,E15.5))"  This does not work with PGI
            myfmt = "(1000(3x,E15.5))"
            WRITE(lambda_phself,'(i9,4x)', ADVANCE = 'no') iqq
            WRITE(lambda_phself, FMT = myfmt) (REAL(lambda_all(imode, iqq, 1, itempphen)), imode = 1, nmodes)
              !
          ENDDO
          CLOSE(lambda_phself)
          !
          ! SP - 03/2019
          ! \Gamma = 1/\tau = phonon lifetime
          ! \Gamma = - 2 * Im \Pi^R where \Pi^R is the retarted phonon self-energy.
          ! Im \Pi^R = pi*k-point weight*[f(E_k+q) - f(E_k)]*delta[E_k+q - E_k - w_q]
          ! Since gamma_all = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q] we have
          ! \Gamma = 2 * gamma_all
          filephlinewid = 'linewidth.phself.' // trim(adjustl(tp)) // 'K'
          OPEN(UNIT = linewidth_phself, FILE = filephlinewid)
          WRITE(linewidth_phself, '(a)') '# Phonon frequency and phonon lifetime in meV '
          WRITE(linewidth_phself, '(a)') '# Q-point  Mode   Phonon freq (meV)   Phonon linewidth (meV)'
          DO iqq = 1, nqtotf
            !
            DO imode = 1, nmodes
              WRITE(linewidth_phself, '(i9,i6,E20.8,E22.10)') iqq, imode, &
                                     ryd2mev * wf(imode, iqq), 2.0d0 * ryd2mev * REAL(gamma_all(imode, iqq, 1, itempphen))
            ENDDO
            !
          ENDDO
          CLOSE(linewidth_phself)
        ENDDO ! itempphen
      ENDIF
    ENDIF
    IF (band_plot) CALL plot_band()
    !
    IF (fermi_plot) CALL plot_fermisurface()
    !
    IF (a2f) CALL a2f_main()
    !
    ! if scattering is read then Fermi level and scissor have not been computed.
    IF (scatread) THEN
      IF (ABS(scissor) > 0.000001) THEN
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          ikq = ikk + 1
          DO ibnd = icbm, nbndsub
            etf(ibnd, ikk) = etf(ibnd, ikk) + scissor
            etf(ibnd, ikq) = etf(ibnd, ikq) + scissor
          ENDDO
        ENDDO
        WRITE( stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the conduction states")' ) scissor * ryd2ev
      ENDIF
      DO itemp = 1, nstemp
        etemp = gtemp(itemp)
        IF (int_mob .OR. carrier) THEN
          ! SP: Determination of the Fermi level for intrinsic or doped carrier
          !     One also need to apply scissor before calling it.
          CALL fermicarrier(itemp, etemp, ef0, efcb, ctype)
          ! only compute dos for metals
          IF (assume_metal) THEN
            CALL compute_dos(itemp, ef0, dos)
          ENDIF
        ELSE
          IF (efermi_read) THEN
            ef0(itemp) = fermi_energy
          ELSE !SP: This is added for efficiency reason because the efermig routine is slow
            ef0(itemp) = efnew
          ENDIF
        ENDIF
      ENDDO ! itemp
      IF (.NOT. iterative_bte) CALL transport_coeffs(ef0, efcb)
    ENDIF ! if scattering
    !
    ! Now deallocate
    DEALLOCATE(epf17, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating epf17', 1)
    DEALLOCATE(selecq, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating selecq', 1)
    IF (scattering .AND. .NOT. iterative_bte) THEN
      DEALLOCATE(inv_tau_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating inv_tau_all', 1)
      DEALLOCATE(zi_allvb, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating zi_allvb', 1)
    ENDIF
    IF (int_mob .AND. carrier) THEN
      DEALLOCATE(inv_tau_allcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating inv_tau_allcb', 1)
      DEALLOCATE(zi_allcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating zi_allcb', 1)
    ENDIF
    IF (elecselfen .OR. plselfen) THEN
      DEALLOCATE(sigmar_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating sigmar_all', 1)
      DEALLOCATE(sigmai_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating sigmai_all', 1)
      DEALLOCATE(zi_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating zi_all', 1)
      IF (iverbosity == 3) THEN
        DEALLOCATE(sigmai_mode, STAT = ierr)
        IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating sigmai_mode', 1)
      ENDIF
    ENDIF
    IF (phonselfen) THEN
      DEALLOCATE(lambda_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating lambda_all', 1)
      DEALLOCATE(lambda_v_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating lambda_v_all', 1)
      DEALLOCATE(gamma_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating gamma_all', 1)
      DEALLOCATE(gamma_v_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating gamma_v_all', 1)
    ENDIF
    IF (specfun_el .OR. specfun_pl) THEN
      DEALLOCATE(esigmar_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating esigmar_all', 1)
      DEALLOCATE(esigmai_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating esigmai_all', 1)
      DEALLOCATE(a_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating a_all', 1)
    ENDIF
    IF (specfun_ph) THEN
      DEALLOCATE(a_all_ph, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating a_all_ph', 1)
    ENDIF
    IF (lifc) THEN
      DEALLOCATE(wscache, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wscache', 1)
    ENDIF
    !
    ! Now do the second step of mobility
    IF (iterative_bte) THEN
      ALLOCATE(vkk_all(3, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating vkk_all', 1)
      ALLOCATE(wkf_all(nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error allocating wkf_all', 1)
      vkk_all(:, :, :) = zero
      wkf_all(:) = zero
      !
      CALL iter_restart(etf_all, wkf_all, vkk_all, ind_tot, ind_totcb, ef0, efcb)
      !
      DEALLOCATE(vkk_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating vkk_all', 1)
      DEALLOCATE(wkf_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wkf_all', 1)
      DEALLOCATE(inv_tau_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating inv_tau_all', 1)
      DEALLOCATE(inv_tau_allcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating inv_tau_allcb', 1)
    ENDIF
    !
  ENDIF ! (iterative_bte .AND. epmatkqread)
  !
  IF ((iterative_bte .OR. ephwrite) .AND. mp_mesh_k) THEN
    DEALLOCATE(map_rebal, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating map_rebal', 1)
    DEALLOCATE(map_rebal_inv, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating map_rebal_inv', 1)
  ENDIF
  IF (vme == 'wannier') THEN
    DEALLOCATE(vmef, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating vmef', 1)
    DEALLOCATE(cvmew, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cvmew', 1)
  ELSE
    DEALLOCATE(cdmew, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cdmew', 1)
    DEALLOCATE(dmef, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating dmef', 1)
  ENDIF
  !
  DEALLOCATE(chw, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating chw', 1)
  DEALLOCATE(chw_ks, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating chw_ks', 1)
  DEALLOCATE(rdw, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating rdw', 1)
  DEALLOCATE(epsi, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating epsi', 1)
  DEALLOCATE(zstar, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating zstar', 1)
  DEALLOCATE(epmatwef, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating epmatwef', 1)
  DEALLOCATE(wf, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wf', 1)
  DEALLOCATE(etf, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating etf', 1)
  DEALLOCATE(etf_ks, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating etf_ks', 1)
  DEALLOCATE(epmatf, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating epmatf', 1)
  DEALLOCATE(cufkk, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cufkk', 1)
  DEALLOCATE(cufkq, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cufkq', 1)
  DEALLOCATE(uf, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating uf', 1)
  DEALLOCATE(isk_dummy, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating isk_dummy', 1)
  DEALLOCATE(eps_rpa, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating eps_rpa', 1)
  DEALLOCATE(bmatf, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating bmatf', 1)
  DEALLOCATE(w2, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating w2', 1)
  DEALLOCATE(cfac, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cfac', 1)
  DEALLOCATE(cfacq, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating cfacq', 1)
  DEALLOCATE(rdotk, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating rdotk', 1)
  DEALLOCATE(rdotk2, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating rdotk2', 1)
  DEALLOCATE(irvec_r, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating irvec_r', 1)
  DEALLOCATE(irvec_k, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating irvec_k', 1)
  DEALLOCATE(irvec_q, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating irvec_q', 1)
  DEALLOCATE(irvec_g, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating irvec_g', 1)
  DEALLOCATE(ndegen_k, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating ndegen_k', 1)
  DEALLOCATE(ndegen_q, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating ndegen_q', 1)
  DEALLOCATE(ndegen_g, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating ndegen_g', 1)
  DEALLOCATE(wslen_k, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wslen_k', 1)
  DEALLOCATE(wslen_q, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wslen_q', 1)
  DEALLOCATE(wslen_g, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating wslen_g', 1)
  DEALLOCATE(etf_all, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating etf_all', 1)
  ! Deallocate temperature when no cumulant or supercond
  IF ((.NOT. cumulant) .AND. (.NOT. eliashberg)) THEN
    DEALLOCATE(gtemp, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating gtemp', 1)
  ENDIF
  !
  DEALLOCATE(et_ks, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwann_shuffle_mem', 'Error deallocating et_ks', 1)
  IF (assume_metal) THEN
    DEALLOCATE(dos, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating dos', 1)
  ENDIF
  !
  IF (lindabs .AND. (.NOT. scattering)) THEN
    IF (carrier) THEN
      DEALLOCATE(ef0_fca, STAT = ierr)
      IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating ef0_fca', 1)
    ENDIF
    DEALLOCATE(epsilon2_abs, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating epsilon2_abs', 1)
    DEALLOCATE(epsilon2_abs_lorenz, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating epsilon2_abs_lorenz', 1)
    DEALLOCATE(epsilon2_abs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating epsilon2_abs_all', 1)
    DEALLOCATE(epsilon2_abs_lorenz_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwann_shuffle', 'Error deallocating epsilon2_abs_lorenz_all', 1)
  ENDIF
  !
  CALL stop_clock('ephwann')
  !
  !---------------------------------------------------------------------------------------
  END SUBROUTINE ephwann_shuffle_mem
  !---------------------------------------------------------------------------------------
