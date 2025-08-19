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
  SUBROUTINE use_wannier(dims, dims2, nrr_k, irvec_k, ndegen_k, wslen_k, nrr_q, &
                         irvec_q, ndegen_q, wslen_q, nrr_g, irvec_g, ndegen_g, wslen_g)
  !---------------------------------------------------------------------
  !!
  !! Wannier interpolation of electron-phonon vertex
  !!
  !! Real-space Wannier to fine grid Bloch space interpolation
  !-----------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE pwcom,            ONLY : ef, nelec
  USE cell_base,        ONLY : at, bg
  USE ions_base,        ONLY : nat, amass, ityp, tau
  USE modes,            ONLY : nmodes
  USE stop,             ONLY : ephf_deallocate
  USE input,            ONLY : nbndsub, fsthick, epwread, system_2d, plrn, ngaussw,&
                               degaussw, lifc, lscreen, etf_mem, scr_typ, loptabs, &
                               nw_specfun, elecselfen, phonselfen, nest_fn, a2f,   &
                               specfun_ph, vme, eig_read, ephwrite, efermi_read,   &
                               fermi_energy, specfun_el, band_plot, scattering,    &
                               nstemp, int_mob, scissor, carrier, iterative_bte,   &
                               longrange_only, scatread, nqf1, nqf2, nqf3,         &
                               mp_mesh_k, restart, plselfen, prtgkk, specfun_pl,   &
                               lindabs, fermi_plot, epmatkqread, restart_step,     &
                               nsmear, ii_partion, nqc1, nqc2, nqc3, assume_metal, &
                               eliashberg, meshnum, time_rev_U_plrn, lwfpt, ii_g,  &
                               start_mesh, ii_lscreen, isk_dummy, lopt_w2b,        &
                               explrn, gb_scattering, lfast_kmesh, epw_memdist,    &
                               plot_explrn_e, plot_explrn_h, dos_tetra, a2f_iso
  USE control_flags,    ONLY : iverbosity
  USE noncollin_module, ONLY : noncolin
  USE ep_constants,     ONLY : ryd2ev, ryd2mev, one, two, zero, czero, cone,       &
                               twopi, ci, eps6, eps16
  USE io_files,         ONLY : prefix, tmp_dir
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE io_var,           ONLY : iuntaucb, iunepmatwp, iunepmatwp2, iunrestart,      &
                               iuntau
  USE global_var,       ONLY : ifc, do_cutoff_2D_epw, chw, chw_ks, adapt_smearing, &
                               wf, etf, etf_ks, xqf, xkf, wkf, nqtotf, nkqf, epf17,&
                               nkf, nqf, ibndmin, ibndmax, lambda_all, vmef,       &
                               sigmai_all, sigmai_mode, gamma_all, epsi,           &
                               efnew, sigmar_all, zi_all, nkqtotf, eps_rpa,        &
                               zi_allvb, inv_tau_all, eta, nbndfst,                &
                               inv_tau_all_mode, inv_tau_allcb_mode,               &
                               inv_tau_all_freq, inv_tau_allcb_freq,               &
                               inv_tau_allcb, zi_allcb, gamma_v_all, esigmar_all,  &
                               esigmai_all, lower_bnd, upper_bnd, a_all, a_all_ph, &
                               wscache, lambda_v_all, nktotf, gtemp, dos,          &
                               nbndskip, ef0_fca, selecq, emib3tz,                 &
                               epsilon2_abs_all, epsilon2_abs_lorenz_all, eimpf17, &
                               epstf_therm, qtf2_therm, partion, eta_imp, sum_E,   &
                               E_grid, maxdim, dwf17, totq, sigmar_dw_all,         & 
                               image_arr, startq, lastq, pi_0, ctype,              &
                               pir_all, gammai_all, fermi_energies_t
  USE wannier2bloch,    ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch,           &
                               ephwan2blochp, ephwan2bloch, vmewan2bloch,          &
                               dynifc2blochf, vmewan2blochp, rrwan2bloch,          &
                               dwwan2blochp
  USE wigner,           ONLY : wigner_divide_ndegen
  USE io_supercond,     ONLY : write_ephmat, count_kpoints, kmesh_fine, kqmap_fine,&
                               check_restart_ephwrite, write_dos, write_phdos,     &
                               file_open_ephmat, write_dos_tetra
  USE transport_legacy, ONLY : transport_coeffs, scattering_rate_q
  USE bzgrid,           ONLY : qwindow_wrap, loadkmesh_fst, xqf_otf,               &
                               loadqmesh_serial, loadkmesh_para, load_rebal,       &
                               load_mesh
  USE printing,         ONLY : print_gkk, plot_band, plot_fermisurface
  USE io_transport,     ONLY : transport_read, tau_read, iter_open, print_ibte,    &
                               iter_merge
  USE io_selfen,        ONLY : selfen_el_read, spectral_read, selfen_ph_write,     &
                               selfen_el_read_wfpt
  USE transport,        ONLY : transport_setup, transport_main, transport_prepare
  USE close,            ONLY : iter_close
  USE parallelism,      ONLY : fkbounds, image_division
  USE mp,               ONLY : mp_bcast, mp_sum, mp_barrier
  USE io_global,        ONLY : ionode_id, stdout
  USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id, my_image_id
  USE mp_world,         ONLY : mpime, world_comm
  USE mp_images,        ONLY : inter_image_comm, nimage
  USE low_lvl,          ONLY : system_mem_usage, mem_size
  USE utilities,        ONLY : compute_dos, broadening, fermicarrier, fermiwindow, &
                               calcpartion, broadening_imp, fast_fermi
  USE selfen,           ONLY : selfen_phon_q, selfen_elec_q, selfen_pl_q,          &
                               nesting_fn_q, selfen_elec_print
  USE spectral,         ONLY : spectral_func_el_print, spectral_func_ph_q,         &
                               a2f_main, spectral_func_pl_q
  USE longrange,        ONLY : rgd_imp_epw_fine
  USE screening,        ONLY : rpa_epsilon, tf_epsilon, calc_qtf2_therm,           &
                               calc_epstf_therm, calc_eps_static
  USE indabs,           ONLY : indabs_main, renorm_eig, fermi_carrier_indabs,      &
                               prepare_indabs
  USE io_indabs,        ONLY : indabs_read
  USE polaron,          ONLY : plrn_flow_select, plrn_prepare, plrn_save_g_to_file,&
                               is_mirror_q, is_mirror_k, kpg_map, ikq_all
  USE qdabs,            ONLY : qdabs_main, prepare_qdabs
  USE wfpt,             ONLY : wfpt_wan2bloch_setup, dwmatf_trunc,                 &
                               wfpt_wan2bloch_finalize, ahc_run_static_wfpt
  USE wannier2bloch_opt,ONLY : wan2bloch_opt_finalize
  USE supercond_vertex, ONLY : prepare_a2f
#if defined(__HDF5)
  USE expolaron,        ONLY : run_explrn
#endif   
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_MODE_RDONLY, MPI_INFO_NULL, MPI_OFFSET_KIND,    &
                               MPI_OFFSET
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: dims
  !! Dims is either nbndsub if use_ws or 1 if not
  INTEGER, INTENT(in) :: dims2
  !! Dims is either nat if use_ws or 1 if not
  INTEGER, INTENT(in) :: nrr_k
  !! number of electronic WS points
  INTEGER, INTENT(in) :: nrr_q
  !! number of phonon WS points
  INTEGER, INTENT(in) :: nrr_g
  !! number of el-ph WS points
  INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
  !! Coordinates of real space vector for electrons
  INTEGER, INTENT(in) :: irvec_q(3, nrr_k)
  !! Coordinates of real space vector for phonons
  INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
  !! Coordinates of real space vector for electron-phonon
  INTEGER, INTENT(in) :: ndegen_k(nrr_k, dims, dims)
  !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
  INTEGER, INTENT(in) :: ndegen_q(nrr_q, dims2, dims2)
  !! degeneracy of WS points for phonon
  INTEGER, INTENT(in) :: ndegen_g(dims, nrr_g, nat)
  !! Wigner-Seitz weights for the electron-phonon grid that depend on
  REAL(KIND = DP), INTENT(in) :: wslen_k(nrr_k)
  !! real-space length for electron, in units of alat
  REAL(KIND = DP), INTENT(in) :: wslen_q(nrr_q)
  !! real-space length for phonons, in units of alat
  REAL(KIND = DP), INTENT(in) :: wslen_g(nrr_g)
  !! real-space length for electron-phonons, in units of alat
  !
  ! Local variables
  CHARACTER(LEN = 256) :: filint
  !! Name of the file to write/read
  LOGICAL :: already_skipped
  !! Skipping band during the Wannierization
  LOGICAL :: exst
  !! If the file exist
  LOGICAL :: first_cycle
  !! Check wheter this is the first cycle after a restart.
  LOGICAL :: mirror_k
  !! Mirror symmetry for k-points
  LOGICAL :: mirror_q
  !! Mirror symmetry for q-points
  LOGICAL :: mirror_kpq
  !! Mirror symmetry for (k+q)-points
  INTEGER :: ios
  !! INTEGER variable for I/O control
  INTEGER :: iq
  !! Counter on coarse q-point grid
  INTEGER :: iqq
  !! Counter on coarse q-point grid
  INTEGER :: iq_restart
  !! Counter on coarse q-point grid
  INTEGER :: qind
  !! Index of q-point either on full q-grid or within qwindow
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
  INTEGER :: nrws
  !! Number of real-space Wigner-Seitz
  INTEGER(KIND = 8) :: valueRSS(2)
  !! Return virtual and resisdent memory from system
  INTEGER :: ierr
  !! Error status
  INTEGER :: itemp
  !! Temperature index
  INTEGER :: icbm
  !! Index of the CBM
  INTEGER :: ipool
  !! Cpu index.
  INTEGER :: npool_tmp
  !! Temporary number of pools
  INTEGER :: nimage_tmp
  !! Temporary number of images
  INTEGER :: QD_mesh
  !! mesh index for quasidegenerate perturbation
  INTEGER :: mesh_dum, mesh_leng
  !! mesh dummy index for quasidegenerate perturbation
  INTEGER :: cycle_quad
  !! building Hamiltonian or performing calculation
  INTEGER :: nc = 1
  !! Index of the outer interpolation loop
  INTEGER :: lrepmatw2_restart(npool)
  !! To restart opening files
  INTEGER :: lrepmatw5_restart(npool)
  !! To restart opening files
  INTEGER :: lrepmatw
  !! record length while reading file
  INTEGER*8 :: unf_recl
  !! Record length
  INTEGER :: direct_io_factor
  !! Type of IOlength
  INTEGER :: size_image
  !! Size of image 
  INTEGER, PARAMETER :: nrwsx = 200
  !! Maximum number of real-space Wigner-Seitz
  INTEGER :: size_image_arr(nimage)
  !! size of the image array in all images
  INTEGER, ALLOCATABLE :: image_array(:)
  !! Meshgrids for QD which are in this image
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
  REAL(KIND = DP) :: dummy(3)
  !! Dummy variable
  REAL(KIND = DP) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND = DP) :: etemp_fca
  !! Temperature for free carrier absorption
  REAL(KIND = DP) :: ef0(nstemp)
  !! Fermi level for the temperature itemp
  REAL(KIND = DP) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp
  REAL(KIND = DP) :: valmin(npool)
  !! Temporary broadening min value
  REAL(KIND = DP) :: valmax(npool)
  !! Temporary broadening max value
  REAL(KIND = DP) :: xxq_r(3)
  !! Q-points
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
  COMPLEX(KIND = DP), ALLOCATABLE :: epmatwef(:, :, :, :)
  !! e-p matrix  in el wannier - fine Bloch phonon grid
  COMPLEX(KIND = DP), ALLOCATABLE :: epmatf(:, :, :)
  !! e-p matrix  in smooth Bloch basis, fine mesh
  COMPLEX(KIND = DP), ALLOCATABLE :: cufkk(:, :)
  !! Rotation matrix, fine mesh, points k
  COMPLEX(KIND = DP), ALLOCATABLE :: cufkq(:, :)
  !! the same, for points k+q
  COMPLEX(KIND = DP), ALLOCATABLE :: uf(:, :)
  !! Rotation matrix for phonons
  COMPLEX(KIND = DP), ALLOCATABLE :: cfac(:)
  !! Used to store $e^{2\pi r \cdot k}$ exponential
  COMPLEX(KIND = DP), ALLOCATABLE :: cfacq(:)
  !! Used to store $e^{2\pi r \cdot k+q}$ exponential
  COMPLEX(KIND = DP), ALLOCATABLE :: vmefp(:, :, :)
  !! Phonon velocity
  COMPLEX(KIND = DP), ALLOCATABLE :: tmp(:, :, :)
  !! Overlap k and k+q
  COMPLEX(KIND = DP), ALLOCATABLE :: rrf(:, :, :)
  !! Interpolation position matrix elements on the fine mesh (ipol, nbnd, nbnd)
  COMPLEX(KIND = DP), ALLOCATABLE :: eimpmatf(:, :)
  !! carrier-ionized impurity matrix in smooth Bloch basis
  CHARACTER(LEN = 256) :: my_image_id_ch
  !! image id
  !
  CALL start_clock('use_wannier')
  !
  IF (scattering) CALL transport_prepare()
  !
  ! Load the fine-grid q and k grids (nkqtotf is computed inside).
  ! The k-point mapping between BZ and IBZ is allocated and computed here as well bztoibz and s_bztoibz
  ! lfast_kmesh is a special optimization level to deal with ultra dense fine homogeneous grids.
  !
  ! Compute the Fermi energy as a function of temperature if required when using lfast_kmesh option
  IF (lfast_kmesh) CALL fast_fermi(nrr_k, irvec_k, efcb)
  !
  CALL load_mesh(nrr_k, irvec_k)
  !
  IF (.NOT. lfast_kmesh) THEN
    ALLOCATE(wf(nmodes, nqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating wf', 1)
    wf(:, :) = zero
  ENDIF
  !
  ALLOCATE(epmatwef(nbndsub, nbndsub, nrr_k, nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating epmatwef', 1)
  ALLOCATE(etf(nbndsub, nkqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating etf', 1)
  ALLOCATE(etf_ks(nbndsub, nkqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating etf_ks', 1)
  ALLOCATE(epmatf(nbndsub, nbndsub, nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating epmatf', 1)
  ALLOCATE(eimpmatf(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating eimpmatf', 1)
  ALLOCATE(cufkk(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating cufkk', 1)
  ALLOCATE(cufkq(nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating cufkq', 1)
  ALLOCATE(uf(nmodes, nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating uf', 1)
  ALLOCATE(eps_rpa(nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating eps_rpa', 1)
  ALLOCATE(epstf_therm(nstemp), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating epstf_therm', 1)
  ALLOCATE(isk_dummy(nkqf), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating isk_dummy', 1)
  ALLOCATE(w2(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating w2', 1)
  epmatwef(:, :, :, :) = czero
  etf(:, :)            = zero
  etf_ks(:, :)         = zero
  epmatf(:, :, :)      = czero
  eimpmatf(:, :)       = czero
  cufkk(:, :)          = czero
  cufkq(:, :)          = czero
  uf(:, :)             = czero
  eps_rpa(:)           = czero
  epstf_therm(:)       = zero
  w2(:)                = zero
  isk_dummy(:)         = 0
  !
  ! Allocate velocity and dipole matrix elements after getting grid size
  ALLOCATE(vmef(3, nbndsub, nbndsub, 2 * nkf), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating vmef', 1)
  vmef(:, :, :, :) = czero
  ALLOCATE(rrf(3, nbndsub, nbndsub), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating rrf', 1)
  rrf(:, :, :) = czero
  ALLOCATE(cfac(nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating cfac', 1)
  ALLOCATE(cfacq(nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating cfac', 1)
  ALLOCATE(rdotk(nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating rdotk', 1)
  ALLOCATE(rdotk2(nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating rdotk2', 1)
  ! This is simply because dgemv take only real number (not integer)
  ALLOCATE(irvec_r(3, nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating irvec_r', 1)
  irvec_r = REAL(irvec_k, KIND = DP)
  !
  ! Zeroing everything - initialization is important !
  cfac(:)   = czero
  cfacq(:)  = czero
  rdotk(:)  = zero
  rdotk2(:) = zero
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
    cfac(:) = EXP(ci * rdotk(:))
    !
    CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ik), chw, cfac)
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
        WRITE(stdout, '(/5x,"Skipping ", i4, " occupied bands:")') nbndskip
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
        WRITE(stdout, '(/5x,"Skipping ", i4, " occupied bands:")') nbndskip
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
      CALL errore("use_wannier", "A scissor shift is applied but the material is a metal...", 1)
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
  ! Re-order the k-point according to weather they are in or out of the fshick windows
  ! In the case of lfast_kmesh, we only have k-points within the fsthick (no need for re-ordering).
  ! SM: using a2f_iso we calculate isotropic a2f without writing to file
  IF ((iterative_bte .OR. ephwrite .OR. a2f_iso) .AND. mp_mesh_k .AND. (.NOT. lfast_kmesh)) THEN
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
  IF (etf_mem == 1) then
    ! Check for directory given by "outdir"
    !
    filint = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
    CALL MPI_FILE_OPEN(world_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatwp2, ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'error in MPI_FILE_OPEN', 1)
  ENDIF
#else
  lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
  filint = TRIM(tmp_dir) // TRIM(prefix)//'.epmatwp'
  INQUIRE(IOLENGTH = direct_io_factor) dummy(1)
  unf_recl = direct_io_factor * INT(lrepmatw, KIND = KIND(unf_recl))
  IF (unf_recl <= 0) CALL errore('use_wannier', 'wrong record length', 3)
  OPEN(iunepmatwp, FILE = TRIM(ADJUSTL(filint)), IOSTAT = ierr, FORM='unformatted', &
       STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
  IF (ierr /= 0) CALL errore('use_wannier', 'error opening ' // TRIM(filint), 1)
#endif
  !
  ! get the size of the matrix elements stored in each pool
  ! for informational purposes.  Not necessary
  !
  CALL mem_size(nmodes, nkf)
  !
  ALLOCATE(etf_all(nbndfst, nktotf), STAT = ierr)
  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating etf_all', 1)
  etf_all(:, :) = zero
  !
  ! Allocate dos we do metals
  IF (assume_metal) THEN
    ALLOCATE(dos(nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore("use_wannier", "Error allocating dos", 1)
  ENDIF
  !
  ! Charge transport using the iteratibe BTE
  IF (iterative_bte) THEN
    !The transport works in two steps:
    ! 1) compute the dominant scattering rates and store them to file
    ! 2) read them from file and solve the IBTE where all important element are in memory
    CALL transport_setup(lrepmatw2_restart, lrepmatw5_restart)
    !
    IF (epmatkqread) THEN
      CALL transport_main(etf_all, ind_tot, ind_totcb, ef0, efcb)
    ELSE
      ! Open the required files
      CALL iter_open(ind_tot, ind_totcb, lrepmatw2_restart, lrepmatw5_restart)
    ENDIF ! epmatkqread
  ENDIF ! iterative_bte
  !
  IF (.NOT. epmatkqread) THEN
    !
    IF (lifc) THEN
      ALLOCATE(wscache(-2 * nqc3:2 * nqc3, -2 * nqc2:2 * nqc2, -2 * nqc1:2 * nqc1, nat, nat), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating wscache', 1)
      wscache(:, :, :, :, :) = zero
    ENDIF
    !
    ! -----------------------------------------------------------------------
    ! Determines which q-points falls within the fsthick windows
    ! Store the result in the selecq.fmt file
    ! If the file exists, automatically restart from the file
    ! -----------------------------------------------------------------------
    CALL qwindow_wrap(totq, nrr_k, irvec_r, cufkk, cufkq)
    !
    ! -----------------------------------------------------------------------
    ! Possible restart during step 1)
    ! -----------------------------------------------------------------------
    iq_restart = 1
    first_cycle = .FALSE.
    !
    ! Fine mesh set of g-matrices.  It is large for memory storage
    ALLOCATE(epf17(nbndfst, nbndfst, nmodes, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating epf17', 1)
    ALLOCATE(tmp(nbndfst, nbndfst, nmodes), STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating tmp', 1)
    ALLOCATE(eimpf17(nbndfst, nbndfst, nkf), STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating eimpf17', 1)
    ! We allocate the phonon frequency on the q-points within the window
    IF (lfast_kmesh) THEN
      ALLOCATE(wf(nmodes, totq), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating wf', 1)
      wf(:, :) = zero
    ENDIF
    IF (phonselfen .OR. specfun_ph) THEN
      ALLOCATE(pi_0(nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating pi_0', 1)
      ALLOCATE(lambda_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating lambda_all', 1)
      ALLOCATE(lambda_v_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating lambda_v_all', 1)
      ALLOCATE(gamma_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating gamma_all', 1)
      ALLOCATE(gamma_v_all(nmodes, totq, nsmear, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating gamma_v_all', 1)
      pi_0(:)                  = zero
      lambda_all(:, :, :, :)   = zero
      lambda_v_all(:, :, :, :) = zero
      gamma_all(:, :, :, :)    = zero
      gamma_v_all(:, :, :, :)  = zero
    ENDIF
    IF (specfun_el .OR. specfun_pl) THEN
      ALLOCATE(esigmar_all(nbndfst, nktotf, nw_specfun, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating esigmar_all', 1)
      ALLOCATE(esigmai_all(nbndfst, nktotf, nw_specfun, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating esigmai_all', 1)
      ALLOCATE(a_all(nw_specfun, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating a_all', 1)
      esigmar_all(:, :, :, :) = zero
      esigmai_all(:, :, :, :) = zero
      a_all(:, :, :) = zero
    ENDIF
    IF (specfun_ph) THEN
      ALLOCATE(pir_all(nw_specfun, nmodes, nstemp, nsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating pir_all', 1)
      ALLOCATE(gammai_all(nw_specfun, nmodes, nstemp, nsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating gammai_all', 1)
      ALLOCATE(a_all_ph(nw_specfun, totq, nstemp, nsmear), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating a_all_ph', 1)
      pir_all(:, :, :, :)     = zero
      gammai_all(:, :, :, :)  = zero
      a_all_ph(:, :, :, :) = zero
    ENDIF
    IF (scattering .AND. .NOT. iterative_bte) THEN
      ALLOCATE(inv_tau_all(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating inv_tau_all', 1)
      ALLOCATE(zi_allvb(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating zi_allvb', 1)
      inv_tau_all(:, :, :) = zero
      zi_allvb(:, :, :)    = zero
    ENDIF
    IF (int_mob .AND. carrier .AND. .NOT. iterative_bte) THEN
      ALLOCATE(inv_tau_allcb(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating inv_tau_allcb', 1)
      ALLOCATE(zi_allcb(nstemp, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating zi_allcb', 1)
      inv_tau_allcb(:, :, :) = zero
      zi_allcb(:, :, :)      = zero
    ENDIF
    IF (elecselfen .OR. plselfen) THEN
      ALLOCATE(sigmar_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating sigmar_all', 1)
      ALLOCATE(sigmai_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating sigmai_all', 1)
      ALLOCATE(zi_all(nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating zi_all', 1)
      sigmar_all(:, :, :) = zero
      sigmai_all(:, :, :) = zero
      zi_all(:, :, :)     = zero
      IF (iverbosity == 3) THEN
        ALLOCATE(sigmai_mode(nbndfst, nmodes, nktotf, nstemp), STAT = ierr)
        IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating sigmai_mode', 1)
        sigmai_mode(:, :, :, :) = zero
      ENDIF
    ENDIF ! elecselfen
    !
    IF (loptabs) THEN
      nc = 2
      IF(lindabs) THEN
        WRITE(stdout, '(/5x,a)') 'lindabs not allowed with loptabs= .true.'
        lindabs = .False.
      ENDIF
    ENDIF
    !
    IF (lindabs) CALL prepare_indabs()
    !
    IF (lwfpt) CALL wfpt_wan2bloch_setup(nrr_k, nrr_g, irvec_k, irvec_g, ndegen_k, ndegen_g, dims, dims2)
    !
    ! Polaron calculations: iq_restart may be set to totq in restart/interpolation mode.
    IF (plrn) call plrn_prepare(totq, iq_restart)
    ! Exciton polaron calculations: iq_restart = totq + 1 to skip the interpolation in case of plotting
    IF (explrn) THEN
      IF (plot_explrn_e .OR. plot_explrn_h) iq_restart = totq + 1
    ENDIF
    !
    ! Restart in SERTA case or self-energy (electron or plasmon) case
    IF (restart) THEN
      IF (elecselfen .OR. plselfen) THEN
        IF (lwfpt) THEN
          CALL selfen_el_read_wfpt(iq_restart, totq, nktotf, sigmar_all, sigmai_all, zi_all, sigmar_dw_all)
        ELSE
          CALL selfen_el_read(iq_restart, totq, nktotf, sigmar_all, sigmai_all, zi_all)
        ENDIF
      ENDIF
      IF (specfun_el .OR. specfun_pl) THEN
        CALL spectral_read(iq_restart, totq, nktotf, esigmar_all, esigmai_all)
      ENDIF
      IF (scattering) THEN
        IF (iterative_bte) CALL transport_read(iq_restart, totq, lrepmatw2_restart, lrepmatw5_restart, &
                                               ind_tot, ind_totcb, first_cycle)
        IF (int_mob .AND. carrier) THEN
          ! Here inv_tau_all and inv_tau_allcb gets updated
          CALL tau_read(iq_restart, totq, nktotf, .TRUE.)
        ELSE
          ! Here inv_tau_all gets updated
          CALL tau_read(iq_restart, totq, nktotf, .FALSE.)
        ENDIF
      ENDIF
      !
      ! Restart in indirect optics
      IF (lindabs) THEN
        CALL indabs_read(iq_restart, totq, epsilon2_abs_all, epsilon2_abs_lorenz_all)
      ENDIF
      ! If you restart from reading a file. This prevent
      ! the case were you restart but the file does not exist
      IF (iq_restart > 1) first_cycle = .TRUE.
      !
    ENDIF ! restart
    !
    ! Scatread assumes that you alread have done the full q-integration
    ! We just do one loop to get interpolated eigenenergies.
    IF (scatread) iq_restart = totq - 1
    !
    ! Restart in Superconductivity case
    IF (ephwrite) THEN
      WRITE(my_image_id_ch, "(I0)") my_image_id
      IF (ionode) THEN
        INQUIRE(FILE = 'restart' // '_' // TRIM(my_image_id_ch) // '.fmt', EXIST = exst)
      ENDIF
      CALL mp_bcast(exst, ionode_id, inter_pool_comm)
      !
      IF (exst) THEN
        IF (ionode) THEN
          OPEN(UNIT = iunrestart, FILE = 'restart' // '_' // TRIM(my_image_id_ch) // '.fmt', STATUS = 'old', IOSTAT = ios)
          READ(iunrestart, *) iq_restart
          READ(iunrestart, *) ind_tot
          READ(iunrestart, *) ind_totcb
          READ(iunrestart, *) npool_tmp
          READ(iunrestart, *) nimage_tmp
          DO ipool = 1, npool
            READ(iunrestart, *) lrepmatw2_restart(ipool)
          ENDDO
          DO ipool = 1, npool
            READ(iunrestart, *) lrepmatw5_restart(ipool)
          ENDDO
          CLOSE(iunrestart)
        ENDIF
        CALL mp_bcast(iq_restart, ionode_id, inter_pool_comm)
        CALL mp_bcast(npool_tmp, ionode_id, inter_pool_comm)
        CALL mp_bcast(lrepmatw2_restart, ionode_id, inter_pool_comm)
        CALL mp_bcast(lrepmatw5_restart, ionode_id, inter_pool_comm)
        IF (npool /= npool_tmp) CALL errore('use_wannier','Number of pools is different',1)
        IF (nimage /= nimage_tmp) CALL errore('use_wannier','Number of nimages is different',1)
        !
#if defined(__MPI)
        CALL MPI_BCAST(ind_tot,   1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
        CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
#endif
        IF (ierr /= 0) CALL errore('use_wannier', 'error in MPI_BCAST', 1)
        !
        IF(iq_restart > 1) THEN
          first_cycle = .TRUE.
          IF (ephwrite .AND. iq_restart + 1 <= totq) THEN
            CALL check_restart_ephwrite(iq_restart)
          ENDIF
        ENDIF
        !
        ! Now, the iq_restart point has been done, so we need to do the next
        iq_restart = iq_restart + 1
        !
        IF (iq_restart <= totq) THEN
          WRITE(stdout, '(5x,a,i8,a)')'We restart from ', iq_restart, ' q-points'
        ELSE
          WRITE(stdout, '(5x,a)')'All q-points are done, no need to restart !!'
        ENDIF
        !
      ENDIF ! exst
    ENDIF
    !
    ! Restart in Superconductivity case
    ! SM: To do: Image paralleilization for anisotropic superconducvtvity calculations
    IF (ephwrite) THEN
      IF (ionode) THEN
        INQUIRE(FILE = 'restart.fmt', EXIST = exst)
      ENDIF
      CALL mp_bcast(exst, ionode_id, inter_pool_comm)
      !
      IF (exst) THEN
        IF (my_pool_id == ionode_id) THEN
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
        IF (npool /= npool_tmp) CALL errore('use_wannier','Number of pools is different',1)
        !
#if defined(__MPI)
        CALL MPI_BCAST(ind_tot,   1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
        CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
#endif
        IF (ierr /= 0) CALL errore('use_wannier', 'error in MPI_BCAST', 1)
        !
        IF(iq_restart > 1) THEN
          first_cycle = .TRUE.
          IF (ephwrite .AND. iq_restart + 1 <= totq) THEN
            CALL check_restart_ephwrite(iq_restart)
          ENDIF
        ENDIF
        !
        ! Now, the iq_restart point has been done, so we need to do the next
        iq_restart = iq_restart + 1
        !
        IF (iq_restart <= totq) THEN
          WRITE(stdout, '(5x,a,i8,a)')'We restart from ', iq_restart, ' q-points'
        ELSE
          WRITE(stdout, '(5x,a)')'All q-points are done, no need to restart !!'
        ENDIF
        !
      ENDIF ! exst
    ENDIF
    !
    ! Restart in Superconductivity case
    ! SM: To do: Image paralleilization for anisotropic superconducvtvity calculations
    IF (ephwrite) THEN
      IF (ionode) THEN
        INQUIRE(FILE = 'restart.fmt', EXIST = exst)
      ENDIF
      CALL mp_bcast(exst, ionode_id, world_comm)
      !
      IF (exst) THEN
        IF (my_pool_id == ionode_id) THEN
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
        IF (npool /= npool_tmp) CALL errore('use_wannier','Number of cores is different',1)
        !
        !
#if defined(__MPI)
        CALL MPI_BCAST(ind_tot,   1, MPI_OFFSET, ionode_id, world_comm, ierr)
        CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
#endif
        IF (ierr /= 0) CALL errore('use_wannier', 'error in MPI_BCAST', 1)
        !
        IF(iq_restart > 1) THEN
          first_cycle = .TRUE.
          IF (ephwrite .AND. iq_restart + 1 <= totq) THEN
            CALL check_restart_ephwrite(iq_restart)
          ENDIF
        ENDIF
        !
        ! Now, the iq_restart point has been done, so we need to do the next
        iq_restart = iq_restart + 1
        !
        IF (iq_restart <= totq) THEN
          WRITE(stdout, '(5x,a,i8,a)')'We restart from ', iq_restart, ' q-points'
        ELSE
          WRITE(stdout, '(5x,a)')'All q-points are done, no need to restart !!'
        ENDIF
        !
      ENDIF ! exst
    ENDIF
    !------------------------------------------------------------------------
    !
    ! Adaptative smearing when degauss = 0
    adapt_smearing = .FALSE.
    IF (ABS(degaussw) < eps16) THEN
      ALLOCATE(eta(nmodes, nbndfst, nkf), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating eta', 1)
      ALLOCATE(vmefp(3, nmodes, nmodes), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating vmefp', 1)
      ALLOCATE(eta_imp(nbndfst, nkf), STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating eta_imp', 1)
      eta(:, :, :)   = zero
      eta_imp(:, :)  = zero
      vmefp(:, :, :) = czero
      adapt_smearing = .TRUE.
    ENDIF
    !
    !SM : For isotropic a2f with image paralleilization, the a2f is calculated on-the-fly
    ! without storing g2 matrices to file. 
    IF (a2f_iso) THEN
      CALL prepare_a2f(nrr_k, nrr_q, nrr_g, irvec_r, irvec_q, irvec_g, &
                        ndegen_k, ndegen_q, ndegen_g, nrws, rws, dims, dims2)
      WRITE (stdout, '(5x, a)') 'We computed a2f using image paralleilization !'
      ! TYK: A workaround to skip the main interpolation loop without changing much
      !      the main interpolation loop. iq_restart set this way is always greater
      !      than totq, which is to be reset after image_division is called is called
      iq_restart = totq + 1
    ENDIF
    !
    ! V.-A. Ha: Low-frequency dielectric tensor for ionized-impurity scattering
    IF (ii_g) CALL calc_eps_static(nrr_q, irvec_q, ndegen_q, rws, nrws)
    !----------------------------------------------------------------------------!
    ! S. Tiwari: QD loop, There are two extra loops that appear when we do QDPT
    ! Loop 1 runs over the many-body energy slices which we call mesh and loop 2
    ! runs over nc=2 for two cycles, cycle 1 for building QDPT hamiltonian and
    ! cycle 2 for summing imaginary dielectric constant. In case of
    ! loptabs=False, nc=1 and mesh_leng=1, and hence, it does not affect other
    ! calculations.
    IF (loptabs) CALL prepare_qdabs(nrr_q, irvec_q, ndegen_q, rws, nrws)
    QD_mesh = 1
    IF (start_mesh == 0) THEN
      mesh_leng = meshnum
    ELSE
      mesh_leng = meshnum - start_mesh + 3
    ENDIF
    ! S.T: Below we divide the q-grid into images 
    !
    CALL image_division(totq, selecq, image_array, size_image_arr, size_image)
    ! 
    totq = size_image
    lastq = image_array(totq)
    startq = image_array(1)
    !
    DO mesh_dum = 1, mesh_leng
      QD_mesh = mesh_dum
      DO cycle_quad = 1, nc
        IF (cycle_quad == 1) THEN
          maxdim = 1
        ENDIF
        IF (mesh_dum > 1) THEN
          IF (mesh_dum < meshnum) THEN
            IF (start_mesh == 0) THEN
              QD_mesh = mesh_dum + start_mesh
            ELSE
              QD_mesh = mesh_dum + start_mesh - 2
            ENDIF
          ENDIF
        ENDIF
        sum_E = 0.d0
        IF (QD_mesh > 1) THEN
          ! E_grid is defined in prepare_qdabs and only non-zero after first QD_mesh
          IF ((E_grid(QD_mesh) == 0) .AND. (QD_mesh /= meshnum)) CYCLE
          ! Do not enter if no quasidegenerate state
        ENDIF
        !------------------------------------------------------------------------------!
        DO iqq = iq_restart, totq
          CALL start_clock ('ep-interp')
          !
          epf17(:, :, :, :) = czero
          tmp(:, :, :)      = czero
          eimpf17(:, :, :)  = czero
          cufkk(:, :) = czero
          cufkq(:, :) = czero
          !
          iq = image_array(iqq)
          !
          IF (iq == -1) CYCLE
          !
          qind = iq
          !
          ! In case of big calculation, show progression of iq (especially usefull when
          ! elecselfen = true as nothing happen during the calculation otherwise.
          !
          IF (.NOT. phonselfen) THEN
            IF (MOD(iqq, restart_step) == 0) THEN
              WRITE(stdout, '(5x, a, i10, a, i10)' ) 'Progression iq (fine) = ', iqq, '/', totq
            ENDIF
          ENDIF
          !
          ! xqf has to be in crystal coordinate
          IF (lfast_kmesh) THEN
            ! The q-point coordinate is generate on the fly for each q-point
            qind = iqq
            CALL xqf_otf(iq, xxq)
          ELSE
            xxq = xqf(:, iq)
          ENDIF
          !
          ! Added by Chao Lian for enforcing the time-rev symmetry of e_q
          ! for polaron calculations, xxq_r is coordinates of the mirror point of xxq
          ! for other calculations, xxq_r is xxq.
          IF (plrn) THEN
            IF (is_mirror_q (iq)) THEN
              xxq_r = xqf(:, kpg_map(iq))
              mirror_q = .TRUE.
            ELSE
              xxq_r = xxq
              mirror_q = .FALSE.
            ENDIF
          ELSE
            xxq_r = xxq
            mirror_q = .FALSE.
          ENDIF
          !
          ! ------------------------------------------------------
          ! dynamical matrix : Wannier -> Bloch
          ! ------------------------------------------------------
          !
          IF (.NOT. lifc) THEN
            CALL dynwan2bloch(nmodes, nrr_q, irvec_q, ndegen_q, xxq_r, uf, w2, mirror_q)
          ELSE
            !TODO: apply degeneracy lift in dynifc2blochf
            CALL dynifc2blochf(nmodes, rws, nrws, xxq_r, uf, w2, mirror_q)
          ENDIF
          !
          ! ...then take into account the mass factors and square-root the frequencies...
          !
          DO nu = 1, nmodes
            !
            ! wf are the interpolated eigenfrequencies
            ! (omega on fine grid)
            !
            IF (w2(nu) > 0.0d0) THEN
              wf(nu, qind) =  DSQRT(ABS(w2(nu)))
            ELSE
              wf(nu, qind) = -DSQRT(ABS(w2(nu)))
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
          IF (.NOT. longrange_only) THEN
            CALL ephwan2blochp(nmodes, xxq, irvec_g, ndegen_g, nrr_g, epmatwef, nbndsub, nrr_k, dims, dims2)
          ENDIF
          CALL wigner_divide_ndegen(epmatwef, 1, nbndsub, nrr_k, nmodes, ndegen_k, dims)
          !
          ! Number of k points with a band on the Fermi surface
          fermicount = 0
          !
          IF (lscreen) THEN
            IF (scr_typ == 0) CALL rpa_epsilon(xxq, wf(:, qind), nmodes, epsi, eps_rpa)
            IF (scr_typ == 1) CALL tf_epsilon(xxq, nmodes, epsi, eps_rpa)
          ENDIF
          !
          ! -------------------------------------------------------------
          ! Adaptative smearing: computes phonon velocity
          ! -------------------------------------------------------------
          !
          IF (adapt_smearing) THEN
            CALL vmewan2blochp(xxq, nmodes, nrr_q, irvec_q, ndegen_q, uf, vmefp(:, :, :), wf(:, qind), rws, nrws)
          ENDIF
          !
          ! This is a loop over k blocks in the pool (size of the local k-set)
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
            IF (plrn .AND. time_rev_U_plrn) THEN
              mirror_k = is_mirror_k(ik)
            ELSE
              mirror_k = .FALSE.
            ENDIF
            !
            ! note that ikq_all return the global index of k+q
            ! Since ikq2 is global index, we need is_mirror_q instead of is_mirror_k
            ! is_mirror_q uses global k/q index.
            ! this is different from is_mirror_k which needs local index k
            IF (plrn .AND. time_rev_U_plrn) THEN
              ! kpg_map return global index, thus xkf_all is needed
              mirror_kpq = is_mirror_q(ikq_all(ik, iq))
            ELSE
              mirror_kpq = .FALSE.
            ENDIF
            !
            CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1)
            CALL DGEMV('t', 3, nrr_k, twopi, irvec_r, 3, xkq2, 1, 0.0_DP, rdotk2, 1)
            cfac(:) = EXP(ci * rdotk(:))
            cfacq(:) = EXP(ci * rdotk2(:))
            !
            ! ------------------------------------------------------
            ! hamiltonian : Wannier -> Bloch
            ! ------------------------------------------------------
            !
            ! Kohn-Sham first, then get the rotation matricies for following interp.
            IF (eig_read) THEN
              !
              CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf_ks(:, ikk), chw_ks, cfac, mirror_k)
              CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf_ks(:, ikq), chw_ks, cfacq, mirror_kpq)
              !
            ENDIF
            !
            CALL hamwan2bloch(nbndsub, nrr_k, cufkk, etf(:, ikk), chw, cfac, mirror_k)
            CALL hamwan2bloch(nbndsub, nrr_k, cufkq, etf(:, ikq), chw, cfacq, mirror_kpq)
            !
            ! Apply a possible scissor shift
            etf(icbm:nbndsub, ikk) = etf(icbm:nbndsub, ikk) + scissor
            etf(icbm:nbndsub, ikq) = etf(icbm:nbndsub, ikq) + scissor
            !
            IF (.NOT. scatread) THEN
              !
              ! Interpolate only when (k,k+q) both have at least one band
              ! within a Fermi shell of size fsthick
              IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
                !
                ! Compute velocities
                !
                IF (vme == 'wannier') THEN
                  !
                  ! ------------------------------------------------------
                  !  velocity: Wannier -> Bloch
                  ! ------------------------------------------------------
                  !
                  IF (eig_read) THEN
                    !
                    ! Renormalize the eigenvalues and vmef with the read eigenvalues
                    CALL renorm_eig(ikk, ikq, nrr_k, irvec_k, irvec_r, cufkk, cufkq, cfac, cfacq)
                    !
                  ELSE ! eig_read
                    CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :, ikk), &
                            etf(:, ikk), etf_ks(:, ikk), chw, cfac)
                    CALL vmewan2bloch(nbndsub, nrr_k, irvec_k, cufkq, vmef(:, :, :, ikq), &
                            etf(:, ikq), etf_ks(:, ikq), chw, cfacq)
                  ENDIF
                  !
                ELSE
                  !
                  ! ------------------------------------------------------
                  !  dipole: Wannier -> Bloch
                  ! ------------------------------------------------------
                  !
                  CALL dmewan2bloch(nbndsub, nrr_k, cufkk, vmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), cfac)
                  CALL dmewan2bloch(nbndsub, nrr_k, cufkq, vmef(:, :, :, ikq), etf(:, ikq), etf_ks(:, ikq), cfacq)
                  !
                ENDIF
                !
                ! Compute Berry connection
                !
                CALL rrwan2bloch(nbndsub, nrr_k, cfac, rrf)
                !
                ! Computes adaptative smearing
                !
                IF (adapt_smearing) THEN
                  ! Return the value of the adaptative broadening eta
                  CALL broadening(ik, ikk, ikq, wf(:, qind), vmefp, eta)
                  IF (ii_g) CALL broadening_imp(ik, ikk, ikq, eta_imp)
                  !
                ENDIF ! adapt_smearing
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
                !IF (longrange_only) THEN
                !  !
                !  epmatf(:, :, :) = czero
                !  !
                !ENDIF
                !
                CALL cryst_to_cart(1, xxq, bg, 1)
                epmatf(:, :, :) = czero
                CALL ephwan2bloch(nbndsub, nrr_k, epmatwef, cufkk, cufkq, epmatf, nmodes, cfac, xxq, rrf(:, :, :))
                IF (ii_g) THEN
                  eimpmatf(:, :) = czero
                  CALL rgd_imp_epw_fine(nqc1, nqc2, nqc3, xxq, eimpmatf, cufkk, cufkq, one)
                ENDIF
                CALL cryst_to_cart(1, xxq, at, -1)
                !
                ! Store epmatf in memory
                !
                DO jbnd = ibndmin, ibndmax
                  DO ibnd = ibndmin, ibndmax
                    IF (lscreen) THEN
                      epf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = epmatf(ibnd, jbnd, :) / eps_rpa(:)
                    ELSE
                      epf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, :, ik) = epmatf(ibnd, jbnd, :)
                    ENDIF
                    IF (ii_g) eimpf17(ibnd - ibndmin + 1, jbnd - ibndmin + 1, ik) = eimpmatf(ibnd, jbnd)
                  ENDDO
                ENDDO
                ! Now do the eigenvector rotation: epmatf(j) = sum_i eptmp(i) * uf(i,j)
                tmp(:, :, :) = epf17(:, :, :, ik)
                CALL ZGEMM('n', 'n', nbndfst * nbndfst, nmodes, nmodes, cone, tmp(:,:, :), &
                          nbndfst * nbndfst, uf, nmodes, czero, epf17(:, :, :,ik), nbndfst * nbndfst)
                !
              ELSE
                ! Velocity of all k-points is needed for gb_scattering
                IF (scattering .AND. (iqq == 1 .OR. first_cycle) .AND. gb_scattering) THEN
                  IF (vme == 'wannier') THEN
                    CALL vmewan2bloch &
                        (nbndsub, nrr_k, irvec_k, cufkk, vmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), chw, cfac)
                  ELSE
                    CALL dmewan2bloch &
                        (nbndsub, nrr_k, cufkk, vmef(:, :, :, ikk), etf(:, ikk), etf_ks(:, ikk), cfac)
                  ENDIF
                ENDIF
                !
              ENDIF ! in fsthick
            ENDIF ! scatread
            !
            IF (lwfpt) THEN
              CALL dwwan2blochp(nbndfst, dwmatf_trunc(:, :, :, :, ik), uf, dwf17(:, :, :, ik), nmodes)
            ENDIF ! wfpt
            !
          ENDDO  ! end loop over k points
          !
          IF (MOD(iqq, restart_step) == 0 .AND. adapt_smearing) THEN
            ! Min non-zero value
            valmin(:) = zero
            valmin(my_pool_id + 1) = 100.0d0
            valmax(:) = zero
            DO ik = 1, nkf
              DO ibnd = 1, nbndfst
                DO imode = 1, nmodes
                  IF (eta(imode, ibnd, ik) < valmin(my_pool_id + 1) .AND. ABS(eta(imode, ibnd, ik)) > eps16) THEN
                    valmin(my_pool_id + 1) = eta(imode, ibnd, ik)
                  ENDIF
                  IF (eta(imode, ibnd, ik) > valmax(my_pool_id + 1)) THEN
                    valmax(my_pool_id + 1) = eta(imode, ibnd, ik)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            CALL mp_sum(valmin, inter_pool_comm)
            CALL mp_sum(valmax, inter_pool_comm)
            WRITE(stdout, '(7x, a, f12.6, a)' ) 'Adaptative smearing el-ph = Min: ', DSQRT(2.d0) * MINVAL(valmin) * ryd2mev,' meV'
            WRITE(stdout, '(7x, a, f12.6, a)' ) '                            Max: ', DSQRT(2.d0) * MAXVAL(valmax) * ryd2mev,' meV'
            !
            ! SP - the impurity parts should be called only when ii_g is .true.
            IF (ii_g) THEN
              valmin(:) = zero
              valmin(my_pool_id + 1) = 100.0d0
              valmax(:) = zero
              DO ik = 1, nkf
                DO ibnd = 1, nbndfst
                  IF (eta_imp(ibnd, ik) < valmin(my_pool_id + 1) .AND. ABS(eta_imp(ibnd, ik)) > eps16) THEN
                    valmin(my_pool_id + 1) = eta_imp(ibnd, ik)
                  ENDIF
                  IF (eta_imp(ibnd, ik) > valmax(my_pool_id + 1)) THEN
                    valmax(my_pool_id + 1) = eta_imp(ibnd, ik)
                  ENDIF
                ENDDO
              ENDDO
              CALL mp_sum(valmin, inter_pool_comm)
              CALL mp_sum(valmax, inter_pool_comm)
              WRITE(stdout, '(7x, a, f12.6, a)' ) 'Adaptative smearing el-imp = Min: ', &
                                DSQRT(2.d0) * MINVAL(valmin) * ryd2mev,' meV'
              WRITE(stdout, '(7x, a, f12.6, a)' ) '                             Max: ', &
                                DSQRT(2.d0) * MAXVAL(valmax) * ryd2mev,' meV'
            ENDIF
          ENDIF
          !
          IF (plrn      ) CALL plrn_save_g_to_file(iq, epf17, wf)
          IF (prtgkk    ) CALL print_gkk(iq)
          IF (phonselfen) CALL selfen_phon_q(iqq, iq, totq)
          IF (elecselfen .OR. specfun_el) CALL selfen_elec_q(iqq, iq, totq, first_cycle)
          IF (plselfen .AND. vme == 'dipole') CALL selfen_pl_q(iqq, iq, totq, first_cycle)
          IF (nest_fn   ) CALL nesting_fn_q(iqq, iq)
          IF (specfun_ph) CALL spectral_func_ph_q(iqq, iq, totq)
          IF (specfun_pl .AND. vme == 'dipole') CALL spectral_func_pl_q(iqq, iq, totq, first_cycle)
          IF (ephwrite) THEN
            IF (first_cycle .OR. iqq == 1) THEN
              CALL kmesh_fine
              CALL kqmap_fine
              CALL count_kpoints
              CALL file_open_ephmat(lrepmatw2_restart)
              CALL write_ephmat(iqq, iq, lrepmatw2_restart, first_cycle)
              first_cycle = .FALSE.
            ELSE
              CALL write_ephmat(iqq, iq, lrepmatw2_restart, first_cycle)
            ENDIF
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
              !
              IF (carrier) THEN
                IF (iqq == iq_restart) THEN
                  !
                  ALLOCATE(ef0_fca(nstemp), STAT = ierr)
                  IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating ef0_fca', 1)
                  !
                  IF (ii_g) THEN
                    ALLOCATE(partion(nstemp), STAT = ierr)
                    IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating partion', 1)
                    ALLOCATE(qtf2_therm(nstemp), STAT = ierr)
                    IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating qtf2_therm', 1)
                    qtf2_therm(:) = zero
                  ENDIF
                  !
                  DO itemp = 1, nstemp
                    etemp_fca = gtemp(itemp)
                    IF (.NOT. lfast_kmesh) THEN 
                      CALL fermi_carrier_indabs(itemp, etemp_fca, ef0_fca, ctype)
                    ELSE
                      ef0_fca(itemp) = fermi_energies_t(itemp)
                    ENDIF
                    !
                    IF (ii_g) THEN
                      !
                      IF (ii_partion) THEN
                        CALL calcpartion(itemp, etemp_fca, ctype)
                      ELSE
                        partion(:) = 1.0d0
                      ENDIF
                      !
                      IF (ii_lscreen) THEN
                        CALL calc_qtf2_therm(itemp, etemp_fca, ef0_fca, ef0_fca, ctype)
                      ENDIF
                      !
                    ENDIF ! ii_g
                  ENDDO ! itemp
                  !
                ENDIF ! iq_restart
                !
                IF (ii_g .AND. ii_lscreen) THEN
                  epstf_therm(:) = zero
                  CALL calc_epstf_therm(xxq, nstemp, epsi)
                ENDIF
                !
              ENDIF ! carrier
              !
              CALL indabs_main(iqq, iq, totq, first_cycle, iq_restart)
            ENDIF
            IF ((loptabs) .AND. (.NOT. scattering)) THEN
              CALL qdabs_main(iq,QD_mesh,cycle_quad)
            ENDIF
            !
            ! Conductivity ---------------------------------------------------------
            IF (scattering) THEN
              !
              ! If we want to compute intrinsic mobilities, call fermicarrier to  correctly positionned the ef0 level.
              ! This is only done once for the first iq. Also compute the dos at the same time
              IF (iqq == iq_restart) THEN
                ALLOCATE(partion(nstemp), STAT = ierr)
                IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating partion', 1)
                ALLOCATE(qtf2_therm(nstemp), STAT = ierr)
                IF (ierr /= 0) CALL errore('use_wannier', 'Error allocating qtf2_therm', 1)
                qtf2_therm(:) = zero
                DO itemp = 1, nstemp
                  etemp = gtemp(itemp)
                  IF (ii_partion) THEN
                    CALL calcpartion(itemp, etemp, ctype)
                  ELSE
                    partion(:) = 1.0d0
                  ENDIF 
                  IF (.NOT. lfast_kmesh) THEN 
                    CALL fermicarrier(itemp, etemp, ef0, efcb, ctype)
                  ELSE
                    ef0(itemp) = fermi_energies_t(itemp)
                  ENDIF
                  IF (ii_g .AND. ii_lscreen) THEN
                    CALL calc_qtf2_therm(itemp, etemp, ef0, efcb, ctype)
                  ENDIF
                  ! compute dos for metals
                  IF (assume_metal) THEN
                    CALL compute_dos(itemp, ef0, dos)
                  ENDIF
                ENDDO
              ENDIF
              IF (ii_g .AND. ii_lscreen) THEN
                epstf_therm(:) = zero
                CALL calc_epstf_therm(xxq, nstemp, epsi)
              ELSE
                epstf_therm(:) = one
              ENDIF
              !
              IF (.NOT. iterative_bte .AND. (.NOT. lfast_kmesh)) THEN
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
                  CALL iter_merge()
                  !
                ENDIF
              ENDIF
              !
            ENDIF ! scattering
            CALL stop_clock('ep-interp')
          ENDIF ! scatread
        ENDDO  ! end loop over q points
      ENDDO  ! cycle loop QDPT
    ENDDO  ! loop over QDPT mesh
    !---------------------------------------------------------------------------------!
    !
    IF (plrn) CALL plrn_flow_select(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims)
#if defined(__HDF5)
    !ZD-TODO: skip the wannier interpolation for epw3 calculation
    IF (explrn) CALL run_explrn(nrr_k, ndegen_k, irvec_r, nrr_q, ndegen_q, irvec_q, rws, nrws, dims) 
#endif
    IF (elecselfen .AND. lwfpt) THEN
      CALL ahc_run_static_wfpt(nrr_k, nrr_q, nrr_g, irvec_k, irvec_q, irvec_g, &
                               ndegen_k, ndegen_q, ndegen_g, dims, dims2, rws, nrws)
    ENDIF
    !
    IF (elecselfen) CALL selfen_elec_print()
    IF (specfun_el) CALL spectral_func_el_print()
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
    ! SH: Write the electronic and phonon dos files
    IF ((.NOT. band_plot) .AND. (eliashberg .OR. a2f_iso)) THEN
      ! SM: added dos calculated using tetrahedron method
      IF (dos_tetra) THEN
        CALL write_dos_tetra(ef, nelec)
      ELSE 
        CALL write_dos(ef, nelec)
      ENDIF
    ENDIF
    IF (((.NOT. band_plot) .AND. eliashberg) .OR. a2f .OR. a2f_iso) THEN
      ! HM: In write_phdos, the frequencies will be calculated
      !     for all the q points, including points excluded in qwindow.
      CALL write_phdos(nrr_q, irvec_q, ndegen_q, nrws, rws)
    ENDIF
    !
    IF (phonselfen) CALL selfen_ph_write()
    IF (band_plot)  CALL plot_band()
    IF (fermi_plot) CALL plot_fermisurface()
    IF (a2f)        CALL a2f_main()
    !
    ! if scattering is read then Fermi level and scissor have not been computed.
    IF (scatread) THEN
      IF (ABS(scissor) > 0.000001) THEN
        IF (assume_metal) THEN
          CALL errore("use_wannier", "Cannot apply scissor shift for metals.", 1)
        ENDIF
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
          ! SP: Determination of the Fermi level and dos for intrinsic or doped carrier
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
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating epf17', 1)
    DEALLOCATE(tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating tmp', 1)
    DEALLOCATE(eimpf17, STAT= ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating eimpf17', 1)
    DEALLOCATE(selecq, STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating selecq', 1)
    IF (scattering .AND. .NOT. iterative_bte) THEN
      DEALLOCATE(inv_tau_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating inv_tau_all', 1)
      DEALLOCATE(zi_allvb, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating zi_allvb', 1)
    ENDIF
    IF (int_mob .AND. carrier .AND. .NOT. iterative_bte) THEN
      DEALLOCATE(inv_tau_allcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating inv_tau_allcb', 1)
      DEALLOCATE(zi_allcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating zi_allcb', 1)
    ENDIF
    IF (elecselfen .OR. plselfen) THEN
      DEALLOCATE(sigmar_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating sigmar_all', 1)
      DEALLOCATE(sigmai_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating sigmai_all', 1)
      DEALLOCATE(zi_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating zi_all', 1)
      IF (iverbosity == 3) THEN
        DEALLOCATE(sigmai_mode, STAT = ierr)
        IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating sigmai_mode', 1)
      ENDIF
    ENDIF
    IF (phonselfen .OR. specfun_ph) THEN
      DEALLOCATE(pi_0, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating pi_0', 1)
      DEALLOCATE(lambda_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating lambda_all', 1)
      DEALLOCATE(lambda_v_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating lambda_v_all', 1)
      DEALLOCATE(gamma_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating gamma_all', 1)
      DEALLOCATE(gamma_v_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating gamma_v_all', 1)
    ENDIF
    IF (specfun_el .OR. specfun_pl) THEN
      DEALLOCATE(esigmar_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating esigmar_all', 1)
      DEALLOCATE(esigmai_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating esigmai_all', 1)
      DEALLOCATE(a_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating a_all', 1)
    ENDIF
    IF (specfun_ph) THEN
      DEALLOCATE(pir_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating pir_all', 1)
      DEALLOCATE(gammai_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating gammai_all', 1)
      DEALLOCATE(a_all_ph, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating a_all_ph', 1)
    ENDIF
    IF (lifc) THEN
      DEALLOCATE(wscache, STAT = ierr)
      IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating wscache', 1)
    ENDIF
    IF (lwfpt) CALL wfpt_wan2bloch_finalize()
    !
    IF (lopt_w2b) CALL wan2bloch_opt_finalize()
    !
    ! Now do the second step of mobility
    IF (iterative_bte) CALL transport_main(etf_all, ind_tot, ind_totcb, ef0, efcb)
    !
  ENDIF ! (.NOT. epmatkqread)
  !
  IF (lifc) THEN
    DEALLOCATE(ifc, STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating ifc', 1)
  ENDIF
  !
  IF (system_2d == 'dipole_sh' .AND. .NOT. do_cutoff_2D_epw) THEN
    DEALLOCATE(emib3tz, STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating emib3tz', 1)
  ENDIF
  !
  IF (epwread) THEN
    DEALLOCATE(tau, STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating tau', 1)
    DEALLOCATE(ityp, STAT = ierr)
    IF (ierr /= 0) CALL errore('use_wannier', 'Error deallocating ityp', 1)
  ENDIF ! epwread
  !
  CALL ephf_deallocate(epmatwef, cufkk, cufkq, uf, w2, cfac, cfacq, rdotk, rdotk2, &
                       irvec_r, etf_all, epmatf, eimpmatf)
  !
  CALL stop_clock('use_wannier')
  !
  !--------------------------------------------------------------------------
  END SUBROUTINE use_wannier
  !--------------------------------------------------------------------------
