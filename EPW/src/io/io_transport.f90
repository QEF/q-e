  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !
  !----------------------------------------------------------------------
  MODULE io_transport
  !----------------------------------------------------------------------
  !!
  !! This module contains various writing or reading routines related to transport.
  !! Most of them are for restart purposes.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_ibte(iqq, iq, totq, xxq, ef0, efcb, first_cycle, ind_tot, ind_totcb, &
                          lrepmatw2_restart, lrepmatw5_restart, ctype)
    !-----------------------------------------------------------------------
    !!
    !! This subroutine computes the transition probability and the scattering rates.
    !! Only the elements larger than threshold are saved on file.
    !!
    USE kinds,         ONLY : DP, i4b, i8b
    USE cell_base,     ONLY : omega, at, alat
    USE io_global,     ONLY : stdout
    USE modes,         ONLY : nmodes
    USE input,         ONLY : fsthick, eps_acoustic, degaussw, nstemp, ncarrier,      &
                              assume_metal, lfast_kmesh, nqf1, nqf2, nqf3, system_2d, &
                              mob_maxfreq, mob_nfreq, ii_g, ii_scattering, ii_n,      &
                              ii_only, gb_scattering, gb_only, restart_step
    USE pwcom,         ONLY : ef
    USE global_var,    ONLY : ibndmin, etf, nkf, vmef, wf, wqf,                       &
                              epf17, inv_tau_all, inv_tau_allcb, adapt_smearing,      &
                              wkf, eta, gtemp, lower_bnd, dos,                        &
                              nbndfst, nktotf, vkk_all, carrier_density,              &
                              inv_tau_all_mode, inv_tau_allcb_mode,                   &
                              inv_tau_all_freq, inv_tau_allcb_freq, eimpf17,          &
                              inv_tau_all_MPI, inv_tau_allcb_MPI,                     &
                              inv_tau_all_mode_MPI, inv_tau_allcb_mode_MPI,           &
                              inv_tau_all_freq_MPI, inv_tau_allcb_freq_MPI,           &
                              epstf_therm, partion, eta_imp, inv_tau_gb, evbm, ecbm,  &
                              startq, lastq
    USE ep_constants,  ONLY : zero, one, two, pi, ryd2mev, kelvin2eV, ryd2ev, eps4, eps8,   &
                              eps6, eps20, bohr2ang, ang2cm, hbarJ, eps160, cc2cb,          &
                              electronvolt_si
    USE io_files,      ONLY : diropn
    USE control_flags, ONLY : iverbosity
    USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,     ONLY : world_comm, my_pool_id, npool, my_image_id,             &
                              inter_image_comm, inter_pool_comm, nimage
    USE io_global,     ONLY : ionode_id, ionode
    USE io_var,        ONLY : iunepmat, iunepmatcb, iufilibtev_sup, iunrestart, iuntau,   &
                              iunsparseq, iunsparseqcb, iuntaucb, iufilmu_q
    USE global_var,    ONLY : lrepmatw2_merge, lrepmatw5_merge, threshold
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_SEEK_SET, MPI_INTEGER8, &
                                 MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_INTEGER4, &
                                 MPI_MODE_CREATE, MPI_INFO_NULL, MPI_MODE_WRONLY, MPI_OFFSET
#endif
    USE clib_wrappers,  ONLY : f_copy
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: first_cycle
    !! Use to determine weather this is the first cycle after restart
    INTEGER, INTENT(in) :: iqq
    !! Q-point index from selecq
    INTEGER, INTENT(in) :: iq
    !! Q-point index
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points in selecq
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! Current position inside the file during writing
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! Current position inside the file during writing (electron)
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of element written to file
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of element written to file
#else
    INTEGER(KIND = 8), INTENT(inout) :: ind_tot
    !! Total number of element written to file
    INTEGER(KIND = 8), INTENT(inout) :: ind_totcb
    !! Total number of element written to file
#endif
    INTEGER, INTENT(in) :: ctype
    !! Calculation type: -1 = hole, +1 = electron and 0 = both.
    REAL(KIND = DP), INTENT(in) :: xxq(3)
    !! Current q-point in crystal coordinate.
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: efcb(nstemp)
    !! Second Fermi level for the temperature itemp. Could be unused (0).
    !
    ! Local variables
    INTEGER :: n
    !! Integer for the degenerate average over eigenstates
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: imode
    !! Local mode index
    INTEGER :: itemp
    !! Index over temperature range
    INTEGER :: ifil
    !! Index over the files
    INTEGER :: nu
    !! Index for modes
    INTEGER :: pbnd
    !! Index for bands
    INTEGER :: ierr
    !! Error
    INTEGER :: ipool
    !! Pool index
    INTEGER :: i,j
    !! Cartesian index
    INTEGER :: ifreq
    !! Index on frequency
    INTEGER :: ind(npool)
    !! Nb of Matrix elements that are non-zero
    INTEGER :: indcb(npool)
    !! Nb of Matrix elements that are non-zero in the cb
    INTEGER(KIND = i4b) :: sparse_q(nbndfst * nbndfst * nstemp * nkf)
    !! Index of q-points for mapping
    INTEGER(KIND = i4b) :: sparse_k(nbndfst * nbndfst * nstemp * nkf)
    !! Index of k-points for mapping
    INTEGER(KIND = i4b) :: sparse_i(nbndfst * nbndfst * nstemp * nkf)
    !! Index of i-bands for mapping
    INTEGER(KIND = i4b) :: sparse_j(nbndfst * nbndfst * nstemp * nkf)
    !! Index of j-bands for mapping
    INTEGER(KIND = i4b) :: sparse_t(nbndfst * nbndfst * nstemp * nkf)
    !! Index of temperature for mapping
    INTEGER(KIND = i4b) :: sparsecb_q(nbndfst * nbndfst * nstemp * nkf)
    !! Index of q-points for cb for mapping
    INTEGER(KIND = i4b) :: sparsecb_k(nbndfst * nbndfst * nstemp * nkf)
    !! Index of k-points for cb for mapping
    INTEGER(KIND = i4b) :: sparsecb_i(nbndfst * nbndfst * nstemp * nkf)
    !! Index of i-band for cb for mapping
    INTEGER(KIND = i4b) :: sparsecb_j(nbndfst * nbndfst * nstemp * nkf)
    !! Index of j-band for cb for mapping
    INTEGER(KIND = i4b) :: sparsecb_t(nbndfst * nbndfst * nstemp * nkf)
    !! Index of temeprature for cb for mapping
    REAL(KIND = DP) :: tmp
    !! Temporary variable
    REAL(KIND = DP) :: tmp2
    !! Temporary variable
    REAL(KIND = DP) :: dfnk
    !! Derivative of f_nk with respect to \varepsilon_nk
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: ekq
    !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
    REAL(KIND = DP) :: g2
    !! Electron-phonon matrix elements squared (g2 is Ry^2)
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: w0g1
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q}] $$
    REAL(KIND = DP) :: w0g2
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q}] $$
    REAL(KIND = DP) :: inv_wq(nmodes)
    !! Inverse phonon frequency. Defined for efficiency reasons.
    REAL(KIND = DP) :: inv_etemp
    !! Invese temperature inv_etemp = 1/etemp. Defined for efficiency reasons.
    REAL(KIND = DP) :: g2_tmp(nmodes)
    !! Used to set component to 0 if the phonon freq. is too low. This is defined
    !! for efficiency reasons as if statement should be avoided in inner-most loops.
    REAL(KIND = DP) :: wq(nmodes)
    !! Phonon frequency $$\omega_{q\nu}$$ on the fine grid.
    REAL(KIND = DP) :: wgq(nmodes)
    !! Bose-Einstein occupation function $$n_{q\nu}$$
    REAL(KIND = DP) :: fmkq
    !! Fermi-Dirac occupation function $$f_{m\mathbf{k+q}}$$
    REAL(KIND = DP) :: trans_prob(nbndfst * nbndfst * nstemp * nkf)
    !! Temporary array to store the scattering rates
    REAL(KIND = DP) :: trans_probcb(nbndfst * nbndfst * nstemp * nkf)
    !! Temporary array to store the scattering rates
    REAL(KIND = DP) :: wkf_all(nktotf)
    !! Weights from all the cores
    REAL(KIND = DP) :: inv_eta(nmodes, nbndfst, nktotf)
    !! Inverse of the eta for speed purposes
    REAL(KIND = DP) :: inv_eta_imp(nbndfst, nktotf)
    !! Inverse of the eta impurity for speed purposes
    REAL(KIND = DP) :: etf_all(nbndfst, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP) :: epf2_deg(nbndfst, nbndfst, nmodes)
    !! Epc in degeneracies
    REAL(KIND = DP) :: w_1
    !! Temporary electronic energy
    REAL(KIND = DP) :: w_2
    !! Temporary electronic energy
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac
    REAL(KIND = DP) :: tmpq(nmodes)
    !! Temporary file for mode resolved scattering rates
    REAL(KIND = DP) :: mobilityq(3, 3, nmodes, nstemp)
    !! Mode-resolved mobility
    REAL(KIND = DP) :: inv_cell
    !! cell volume
    REAL(KIND = DP) :: mob(3, 3, nmodes)
    !! Temporary inverse mobility
    REAL(KIND = DP) :: wqf_loc
    !! Local q-point weight
    REAL(KIND = DP) :: step
    !! Energy step in Ry for the spectral decomposition
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Function that returns the Fermi energy
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP) :: w0gimp
    !! $$ \delta[\varepsilon_{nk} - \varepsilon_{mk+q}] $$
    REAL(KIND = DP) :: impurity_density
    !! prefactor, density of impurities per cell
    REAL(KIND = DP) :: invepst2
    !! squared inverse of the thermal thomas-fermi WV at a given itemp
    REAL(KIND = DP) :: eimpf2_deg(nbndfst, nbndfst)
    !! Eimpc in degeneracies
    REAL(KIND = DP) :: tmpimp
    !! Temporary variable
    LOGICAL :: exst
    !! Check if backup files exist
    LOGICAL :: is_restart_loop
    !! Check if this is a restart loop
    INTEGER :: ios
    !! IO status for copying backup files
    CHARACTER(LEN = 256) :: my_image_id_ch
    !! Image id 
    !
    WRITE(my_image_id_ch, "(I0)") my_image_id
    !
    IF (system_2d == 'no') THEN
      inv_cell = 1.0d0 / omega
    ELSE
      ! for 2d system need to divide by area (vacuum in z-direction)
      inv_cell = ( 1.0d0 / omega ) * at(3, 3) * alat
    ENDIF
    !
    ! Weight of the q-points
    IF (lfast_kmesh) THEN
      wqf_loc = 1.0d0 / REAL(nqf1 * nqf2 * nqf3, KIND = DP)
    ELSE
      wqf_loc = wqf(iq)
    ENDIF
    !
    IF (iverbosity == 3) THEN
      ! Energy steps for spectral decomposition
      step = mob_maxfreq / mob_nfreq
      IF (step < eps20) THEN
        CALL errore('print_ibte', 'Too small energy step for spectral decomposition', 1)
      ENDIF
    ENDIF
    !
    IF (iqq == 1 .OR. first_cycle) THEN
      !
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Scattering rate for IBTE")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      IF (ii_scattering) THEN
        WRITE(stdout, '(/5x,a)') REPEAT('=',67)
        WRITE(stdout, '(5x,"Including ionized impurity scattering")')
        WRITE(stdout, '(5x,a,e15.8)' ) 'Using an ionized impurity density of', ii_n
        WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      ENDIF
      IF (ii_only) THEN
        WRITE(stdout, '(/5x,a)') REPEAT('=',67)
        WRITE(stdout, '(5x,"Detected ii_only=.true., omitting carrier-phonon scattering")')
        WRITE(stdout, '(5x,a/)') REPEAT('=',67)
      ENDIF
      WRITE(stdout, '(5x,"No intermediate mobility will be shown.")')
      !
      IF (fsthick < 1.d3) THEN
        WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a)' ) 'This is computed with respect to the fine Fermi level ',ef * ryd2ev, ' eV'
        WRITE(stdout, '(5x,a,f10.6,a,f10.6,a)' ) 'Only states between ', (ef - fsthick) * ryd2ev, ' eV and ', &
                (ef + fsthick) * ryd2ev, ' eV will be included'
        WRITE(stdout, '(5x,a/)')
      ENDIF
      !
      ! We save matrix elements larger than threshold defined in ephwann_shuffle
      WRITE(stdout,'(5x,a,1E20.12)') 'Save matrix elements larger than threshold: ', threshold
      WRITE(stdout,'(5x," ")')
      !
      IF (iverbosity == 3) THEN
        ALLOCATE(vkk_all(3, nbndfst, nktotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error allocating vkk_all', 1)
        wkf_all(:) = zero
        vkk_all(:, :, :) = zero
        etf_all(:, :) = zero
        ! Computes the k-velocity
        DO ik = 1, nkf
          ikk = 2 * ik - 1
          wkf_all(ik + lower_bnd -1 ) = wkf(ikk)
          DO ibnd = 1, nbndfst
            vkk_all(:, ibnd, ik + lower_bnd - 1) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
            etf_all(ibnd, ik + lower_bnd - 1) = etf(ibndmin - 1 + ibnd, ikk)
          ENDDO
        ENDDO
        CALL mp_sum(vkk_all, inter_pool_comm)
        CALL mp_sum(etf_all, inter_pool_comm)
        CALL mp_sum(wkf_all, inter_pool_comm)
        !
        IF (.NOT. assume_metal) THEN
          ALLOCATE(carrier_density(nstemp), STAT = ierr)
          IF (ierr /= 0) CALL errore('print_ibte', 'Error allocating carrier_density', 1)
          carrier_density(:) = zero
          IF (ncarrier < 0.0) THEN ! VB
            CALL carr_density(carrier_density, etf_all, wkf_all, ef0)
          ELSE ! CB
            CALL carr_density(carrier_density, etf_all, wkf_all, efcb)
          ENDIF ! ncarrier
        ENDIF ! assume_metal
        !
        IF (my_pool_id == 0) THEN
          OPEN(UNIT = iufilmu_q, FILE = 'mobility_nuq.fmt')
          WRITE(iufilmu_q, '(a)') '# Mode-resolved contribution in limiting the carrier mobility (Vs/(cm^2))'
          WRITE(iufilmu_q, '(a)') '#     \mu(alpha,beta) = 1.0 / (sum_{\nu q} T_{\nu q}(alpha,beta))'
        ENDIF
      ENDIF ! iverbosity
      !
      IF (gb_scattering) THEN
        ALLOCATE(inv_tau_gb(nbndfst, nktotf), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error allocating inv_tau_gb', 1)
        CALL grain_boundary_scatt(inv_tau_gb)
      ENDIF
      !
    ENDIF ! iqq == 1 .OR. first_cycle
    !
    IF (iverbosity == 3) mobilityq(:, :, :, :) = zero
    !
    IF (first_cycle) THEN
      first_cycle = .FALSE.
    ENDIF
    !
    IF (.TRUE.) THEN
      !
      ! To avoid if branching in the loop
      inv_eta(:, :, :) = zero
      inv_eta_imp(:, :) = zero
      IF (adapt_smearing) THEN
        DO ik = 1, nkf
          DO ibnd = 1, nbndfst
            IF (ii_g) inv_eta_imp(ibnd, ik) = 1.0d0 / (DSQRT(2.0d0) * eta_imp(ibnd, ik))
            DO imode = 1, nmodes
              inv_eta(imode, ibnd, ik) = 1.0d0 / (DSQRT(2d0) * eta(imode, ibnd, ik))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ik = 1, nkf
          DO ibnd = 1, nbndfst
            inv_eta_imp(ibnd, ik) = 1.0d0 / degaussw
            DO imode = 1, nmodes
              inv_eta(imode, ibnd, ik) = 1.0d0 / degaussw
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !
      ! Average the el-ph matrix elements on degenerate bands and phonon modes.
      ! This is important to ensure that the mobility tensor perfectly respects crystal symmetry.
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        wkf_all(ik + lower_bnd - 1) = wkf(ikk)
        !
        ! Average over the k electrons
        DO nu = 1, nmodes
          DO jbnd = 1, nbndfst
            DO ibnd = 1, nbndfst
              w_1 = etf(ibndmin - 1 + ibnd, ikk)
              g2  = zero
              n   = 0
              DO pbnd = 1, nbndfst
                w_2 = etf(ibndmin - 1 + pbnd, ikk)
                IF (ABS(w_2 - w_1) < eps6) THEN
                  n = n + 1
                  g2 = g2 + ABS(epf17(jbnd, pbnd, nu, ik))**two
                ENDIF
              ENDDO
              epf2_deg(jbnd, ibnd, nu) = DSQRT(g2 / FLOAT(n))
            ENDDO
          ENDDO
        ENDDO
        epf17(:, :, :, ik) = epf2_deg(:, :, :)
        !
        ! Average over the k+q electrons
        DO nu = 1, nmodes
          DO jbnd = 1, nbndfst
            DO ibnd = 1, nbndfst
              w_1 = etf(ibndmin - 1 + jbnd, ikq)
              g2 = 0.d0
              n  = 0
              DO pbnd = 1, nbndfst
                w_2 = etf(ibndmin - 1 + pbnd, ikq)
                IF (ABS(w_2 - w_1) < eps6) THEN
                  n = n + 1
                  g2 = g2 + ABS(epf17(pbnd, ibnd, nu, ik))**two
                ENDIF
              ENDDO
              epf2_deg(jbnd, ibnd, nu) = g2 / FLOAT(n)
            ENDDO
          ENDDO
        ENDDO
        !
        ! Note that we already took the square above
        epf17(:, :, :, ik) = epf2_deg(:, :, :)
        !
        ! average impurity matrix elements
        !
        IF (ii_g .AND. ii_scattering) THEN
          DO jbnd = 1, nbndfst
            DO ibnd = 1, nbndfst
              w_1 = etf(ibndmin - 1 + ibnd, ikk)
              g2  = zero
              n   = 0
              DO pbnd = 1, nbndfst
                w_2 = etf(ibndmin - 1 + pbnd, ikk)
                IF (ABS(w_2-w_1) < eps6) THEN
                  n = n + 1
                  g2 = g2 + ABS(eimpf17(jbnd, pbnd, ik))**two
                ENDIF
              ENDDO
              eimpf2_deg(jbnd, ibnd) = DSQRT(g2 / FLOAT(n))
            ENDDO
          ENDDO
          !
          eimpf17(:, :, ik) = eimpf2_deg(:, :)
          !
          ! Average over the k+q electrons
          DO jbnd = 1, nbndfst
            DO ibnd = 1, nbndfst
              w_1 = etf(ibndmin - 1 + jbnd, ikq)
              g2 = 0.d0
              n  = 0
              DO pbnd = 1, nbndfst
                w_2 = etf(ibndmin - 1 + pbnd, ikq)
                IF (ABS(w_2 - w_1) < eps6) THEN
                  n = n + 1
                  g2 = g2 + ABS(eimpf17(pbnd, ibnd, ik))**two
                ENDIF
              ENDDO
              eimpf2_deg(jbnd, ibnd) = g2 / FLOAT(n)
            ENDDO
          ENDDO
          !
          ! Note that we already took the square above
          eimpf17(:, :, ik) = eimpf2_deg(:, :)
        ENDIF
        !
      ENDDO ! ik
      !
      trans_prob(:)    = zero
      sparse_q(:)      = zero
      sparse_k(:)      = zero
      sparse_i(:)      = zero
      sparse_j(:)      = zero
      sparse_t(:)      = zero
      trans_probcb(:)  = zero
      sparsecb_q(:)    = zero
      sparsecb_k(:)    = zero
      sparsecb_i(:)    = zero
      sparsecb_j(:)    = zero
      sparsecb_t(:)    = zero
      etf_all(:, :)    = zero
      ind(:)           = 0
      indcb(:)         = 0
      !
      ! compute impurities per unit cell
      ! impurity_density = ii_n * omega / 6.74822779181357d24
      !
      ! loop over temperatures
      DO itemp = 1, nstemp
        impurity_density = partion(itemp)*ii_n*omega / cc2cb
        ! Define the inverse so that we can efficiently multiply instead of dividing
        etemp = gtemp(itemp)
        inv_etemp = 1.0 / etemp
        invepst2 = (1.0 / epstf_therm(itemp))**2.0d0
        !
        ! Now pre-treat phonon modes for efficiency for this specific current q-point.
        ! Treat phonon frequency and Bose occupation
        wq(:) = zero
        DO imode = 1, nmodes
          IF (lfast_kmesh) THEN
            wq(imode) = wf(imode, iqq)
          ELSE
            wq(imode) = wf(imode, iq)
          ENDIF
          IF (wq(imode) > eps_acoustic) THEN
            g2_tmp(imode) = 1.0d0
            wgq(imode)    = wgauss(-wq(imode) * inv_etemp, -99)
            wgq(imode)    = wgq(imode) / (one - two * wgq(imode))
            inv_wq(imode) =  1.0d0 / (two * wq(imode))
          ELSE
            g2_tmp(imode) = 0.0
            wgq(imode)    = 0.0
            inv_wq(imode) = 0.0
          ENDIF
        ENDDO
        !
        DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF ((MINVAL(ABS(etf(:, ikk) - ef)) < fsthick) .AND. (MINVAL(ABS(etf(:, ikq) - ef)) < fsthick)) THEN
            IF (ctype == 0 .OR. ctype == -1) THEN ! hole
              DO ibnd = 1, nbndfst
                ! Energy at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                ! Fermi-Dirac occupation $f_{nk}$
                fnk = wgauss(-ekk * inv_etemp, -99)
                !
                ! This is to know if we need to store the data
                ! derivative Fermi distribution -df_nk/dE_nk = (f_nk)*(1-f_nk)/ (k_B T)
                dfnk = w0gauss(ekk * inv_etemp, -99 ) * inv_etemp
                !
                DO jbnd = 1, nbndfst
                  !
                  ! Energy and fermi occupation at k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0(itemp)
                  ! Fermi-Dirac occupation $f_{mk+q}$
                  fmkq = wgauss(-ekq * inv_etemp, -99)
                  !
                  ! We perform a sum over the modes
                  tmp  = zero
                  tmp2 = zero
                  !
                  tmpimp = zero
                  IF (ii_scattering .AND. ii_g .AND. (.NOT. gb_only)) THEN
                    w0gimp = w0gauss((ekk - ekq) * inv_eta_imp(ibnd, ik), 0) * inv_eta_imp(ibnd, ik)
                    tmpimp = tmpimp + two * pi * wqf_loc * impurity_density * invepst2 * REAL(eimpf17(jbnd, ibnd, ik)) * w0gimp
                  ENDIF
                  !
                  IF ((.NOT. ii_only) .AND. (.NOT. gb_only)) THEN
                    DO imode = 1, nmodes
                      !
                      ! Here we take into account the zero-point
                      ! DSQRT(hbar/2M\omega)
                      ! with hbar = 1 and M already contained in the eigenmodes
                      ! g2 is Ry^2, wkf must already account for the spin factor
                      ! Note that epf17 has already been squared above during
                      ! averaging.
                      g2 = REAL(epf17(jbnd, ibnd, imode, ik)) * inv_wq(imode) * g2_tmp(imode)
                      !
                      ! delta[E_k - E_k+q + w_q]
                      w0g1 = w0gauss((ekk - ekq + wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      ! delta[E_k - E_k+q - w_q]
                      w0g2 = w0gauss((ekk - ekq - wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      !
                      ! Transition probability - See Eq. 41 of arXiv:1908.01733
                      ! (2019).
                      ! (2 pi/hbar) * (k+q-point weight) * g2 *
                      ! { [f(E_k) + n(w_q)] * delta[E_k - E_k+q - w_q] +
                      !   [1 - f(E_k) + n(w_q)] * delta[E_k - E_k+q + w_q] }
                      !
                      ! This is summed over modes
                      tmp = tmp + two * pi * wqf_loc * g2 * ((fnk + wgq(imode)) * w0g2 + (one - fnk + wgq(imode)) * w0g1)
                      !
                      ! Here we compute the scattering rate - Eq. 64 of
                      ! arXiv:1908.01733 (2019).
                      ! inv_tau = (2 pi/hbar) * (k+q-point weight) * g2 *
                      ! { [f(E_k+q) + n(w_q)] * delta[E_k - E_k+q + w_q] + [1 -
                      ! f(E_k+q) + n(w_q)] * delta[E_k - E_k+q - w_q] }
                      tmp2 = tmp2 + two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      !
                    ENDDO !imode
                  ENDIF
                  !
                  inv_tau_all(ibnd, ik + lower_bnd - 1, itemp) = inv_tau_all(ibnd, ik + lower_bnd - 1, itemp) + tmp2 + tmpimp
                  !
                  ! Only save the elements that really contribute
                  ! The check is made on the SERTA mobility - See Eq. 44 of arXiv:1908.01733 (2019).
                  !
                  IF (ABS((tmp2+tmpimp) * dfnk) > threshold .OR. gb_only) THEN
                    !
                    ind(my_pool_id + 1) = ind(my_pool_id + 1) + 1
                    trans_prob(ind(my_pool_id + 1)) = tmp + tmpimp
                    sparse_q(ind(my_pool_id + 1)) = iq
                    sparse_k(ind(my_pool_id + 1)) = ik + lower_bnd - 1
                    sparse_i(ind(my_pool_id + 1)) = ibnd
                    sparse_j(ind(my_pool_id + 1)) = jbnd
                    sparse_t(ind(my_pool_id + 1)) = itemp
                    !
                  ENDIF
                ENDDO !jbnd
                ! Add grain boundary scattering
                IF (gb_scattering .AND. (iqq == totq) .AND. (.NOT. ii_only)) THEN
                  inv_tau_all(ibnd, ik + lower_bnd - 1, itemp) = inv_tau_all(ibnd, ik + lower_bnd - 1, itemp) &
                                                               + inv_tau_gb(ibnd, ik + lower_bnd - 1)
                ENDIF
                !
              ENDDO ! ibnd
              !
              ! Compute the mode-resolved scattering
              IF (iverbosity == 3) THEN
                DO ibnd = 1, nbndfst
                  ekk = etf(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
                  fnk = wgauss(-ekk * inv_etemp, -99)
                  dfnk = w0gauss(ekk * inv_etemp, -99 ) * inv_etemp
                  tmpq(:) = zero
                  DO jbnd = 1, nbndfst
                    ekq = etf(ibndmin - 1 + jbnd, ikq) - ef0(itemp)
                    fmkq = wgauss(-ekq * inv_etemp, -99)
                    DO imode = 1, nmodes
                      g2 = REAL(epf17(jbnd, ibnd, imode, ik)) * inv_wq(imode) * g2_tmp(imode)
                      w0g1 = w0gauss((ekk - ekq + wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      w0g2 = w0gauss((ekk - ekq - wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      tmpq(imode) = tmpq(imode) &
                          + two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      inv_tau_all_mode(imode, ibnd, ik + lower_bnd - 1, itemp) = &
                      inv_tau_all_mode(imode, ibnd, ik + lower_bnd - 1, itemp) + &
                            two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      !
                      ! Spectral decomposition histogram
                      !
                      ifreq = NINT(wq(imode) / step) + 1
                      !
                      IF(ifreq <= mob_nfreq) THEN
                        inv_tau_all_freq(ifreq, ibnd, ik + lower_bnd - 1) = inv_tau_all_freq(ifreq, ibnd, ik + lower_bnd - 1) &
                            + two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      ENDIF
                    ENDDO !imode
                  ENDDO ! jbnd
                  DO j = 1, 3
                    DO i = 1, 3
                      tmp = dfnk * vkk_all(i, ibnd, ik + lower_bnd - 1) * vkk_all(j, ibnd, ik + lower_bnd - 1)
                      IF (ABS(tmp) > eps20) THEN
                        DO imode = 1, nmodes
                          mobilityq(i, j, imode, itemp) = mobilityq(i, j, imode, itemp) + &
                                                          wkf_all(ik + lower_bnd - 1) * tmpq(imode) / tmp
                        ENDDO
                      ENDIF
                    ENDDO ! i
                  ENDDO ! j
                ENDDO ! ibnd
              ENDIF ! iverbosity
            ENDIF ! ctype
            !
            ! In this case we are also computing the scattering rate for another Fermi level position
            ! This is used to compute both the electron and hole mobility at the same time.
            IF (ctype == 0 .OR. ctype == 1) THEN
              !
              DO ibnd = 1, nbndfst
                ! Energy at k (relative to Ef)
                ekk = etf(ibndmin - 1 + ibnd, ikk) - efcb(itemp)
                ! Fermi-Diract distribution $f_{nk}$
                fnk = wgauss(-ekk * inv_etemp, -99)
                !
                ! Derivative Fermi distribution (-df_nk/dE_nk) = (f_nk)*(1-f_nk)/ (k_B T)
                dfnk = w0gauss(ekk * inv_etemp, -99) * inv_etemp
                !
                DO jbnd = 1, nbndfst
                  !
                  !  energy and fermi occupation at k+q
                  ekq = etf(ibndmin - 1 + jbnd, ikq) - efcb(itemp)
                  fmkq = wgauss(-ekq * inv_etemp, -99)
                  !
                  tmp  = zero
                  tmp2 = zero
                  tmpimp = zero
                  IF (ii_scattering .AND. ii_g .AND. (.NOT. gb_only)) THEN
                    w0gimp = w0gauss((ekk - ekq) * inv_eta_imp(ibnd, ik), 0) * inv_eta_imp(ibnd, ik)
                    tmpimp = tmpimp + two * pi * wqf_loc * impurity_density * invepst2 * REAL(eimpf17(jbnd, ibnd, ik)) * w0gimp
                  ENDIF
                  !
                  IF ((.NOT. ii_only) .AND. (.NOT. gb_only)) THEN
                    DO imode = 1, nmodes
                      ! Same as above but for conduction bands (electrons)
                      g2 = REAL(epf17(jbnd, ibnd, imode, ik)) * inv_wq(imode) * g2_tmp(imode)
                      w0g1 = w0gauss((ekk - ekq + wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      w0g2 = w0gauss((ekk - ekq - wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      tmp = tmp  + two * pi * wqf_loc * g2 * ((fnk + wgq(imode)) * w0g2 + (one - fnk + wgq(imode)) * w0g1)
                      tmp2 = tmp2 + two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                    ENDDO ! imode
                  ENDIF
                  !
                  inv_tau_allcb(ibnd, ik + lower_bnd - 1, itemp) = inv_tau_allcb(ibnd, ik + lower_bnd - 1, itemp) + tmp2 + tmpimp
                  !
                  IF (ABS((tmp2+tmpimp) * dfnk) > threshold .OR. gb_only) THEN
                    indcb (my_pool_id + 1) = indcb(my_pool_id + 1) + 1
                    trans_probcb(indcb(my_pool_id + 1)) = tmp + tmpimp
                    sparsecb_q(indcb(my_pool_id + 1)) = iq
                    sparsecb_k(indcb(my_pool_id + 1)) = ik + lower_bnd - 1
                    sparsecb_i(indcb(my_pool_id + 1)) = ibnd
                    sparsecb_j(indcb(my_pool_id + 1)) = jbnd
                    sparsecb_t(indcb(my_pool_id + 1)) = itemp
                  ENDIF
                ENDDO !jbnd
                ! Add grain boundary scattering
                IF (gb_scattering .AND. (iqq == totq) .AND. (.NOT. ii_only)) THEN
                  inv_tau_allcb(ibnd, ik + lower_bnd - 1, itemp) = inv_tau_allcb(ibnd, ik + lower_bnd - 1, itemp) &
                                                                 + inv_tau_gb(ibnd, ik + lower_bnd - 1)
                ENDIF
                !
              ENDDO !ibnd
              ! Compute the mode-resolved scattering
              IF (iverbosity == 3) THEN
                DO ibnd = 1, nbndfst
                  ekk = etf(ibndmin - 1 + ibnd, ikk) - efcb(itemp)
                  fnk = wgauss(-ekk * inv_etemp, -99)
                  dfnk = w0gauss(ekk * inv_etemp, -99 ) * inv_etemp
                  tmpq(:) = zero
                  DO jbnd = 1, nbndfst
                    ekq = etf(ibndmin - 1 + jbnd, ikq) - efcb(itemp)
                    fmkq = wgauss(-ekq * inv_etemp, -99)
                    DO imode = 1, nmodes
                      g2 = REAL(epf17(jbnd, ibnd, imode, ik)) * inv_wq(imode) * g2_tmp(imode)
                      w0g1 = w0gauss((ekk - ekq + wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      w0g2 = w0gauss((ekk - ekq - wq(imode)) * inv_eta(imode, ibnd, ik), 0) * inv_eta(imode, ibnd, ik)
                      tmpq(imode) = tmpq(imode) &
                          + two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      inv_tau_allcb_mode(imode, ibnd, ik + lower_bnd - 1, itemp) = &
                      inv_tau_allcb_mode(imode, ibnd, ik + lower_bnd - 1, itemp) + &
                            two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      !
                      ! Spectral decomposition histogram
                      !
                      ifreq = NINT(wq(imode) / step) + 1
                      !
                      IF(ifreq <= mob_nfreq) THEN
                        inv_tau_allcb_freq(ifreq, ibnd, ik + lower_bnd - 1) = inv_tau_allcb_freq(ifreq, ibnd, ik + lower_bnd - 1) &
                          + two * pi * wqf_loc * g2 * ((fmkq + wgq(imode)) * w0g1 + (one - fmkq + wgq(imode)) * w0g2)
                      ENDIF
                    ENDDO !imode
                  ENDDO ! jbnd
                  DO j = 1, 3
                    DO i = 1, 3
                      tmp = dfnk * vkk_all(i, ibnd, ik + lower_bnd - 1) * vkk_all(j, ibnd, ik + lower_bnd - 1)
                      IF (ABS(tmp) > eps20) THEN
                        DO imode = 1, nmodes
                          mobilityq(i, j, imode, itemp) = mobilityq(i, j, imode, itemp) &
                                                          + wkf_all(ik + lower_bnd - 1) * tmpq(imode) / tmp
                        ENDDO
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO ! ibnd
              ENDIF ! iverbosity
            ENDIF ! ctype
          ENDIF ! endif fsthick
        ENDDO ! end loop on k
      ENDDO ! itemp
      ! If the q-point is taken, write on file
      CALL mp_sum(ind, inter_pool_comm)
      CALL mp_sum(indcb, inter_pool_comm)
      !
      IF (assume_metal .AND. iverbosity == 3 .AND. iqq == totq) THEN 
        DEALLOCATE(vkk_all, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating vkk_all', 1)
      ENDIF
      !
      IF (iverbosity == 3 .AND. .NOT. assume_metal) THEN
        CALL mp_sum(mobilityq, inter_pool_comm)
        IF (my_pool_id == 0) THEN
          WRITE(iufilmu_q, '(a, 3f12.8)') '# q-point (crystal)', xxq(:)
          WRITE(iufilmu_q, '(a)') &
              '#  Temp (K)  Mode   Phonon freq (meV)       T_nuq(alpha, beta) (Vs/cm^2)'
          mob(:, :, :) = zero
          DO itemp = 1, nstemp
            etemp = gtemp(itemp)
            mob(:, :, :) = mobilityq(:, :, :, itemp) * (hbarJ * carrier_density(itemp)) &
                           / (electronvolt_si * (bohr2ang * ang2cm) ** 2)
            DO imode = 1, nmodes
              WRITE(iufilmu_q, '(1f12.6, i6, 1f14.8, 3E18.8)')  etemp * ryd2ev / kelvin2eV, imode, wq(imode) * ryd2mev, &
                                                mob(1, 1, imode), mob(1, 2, imode), mob(1, 3, imode)
              WRITE(iufilmu_q, '(32x, 3E18.8)') mob(2, 1, imode), mob(2, 2, imode), mob(2, 3, imode)
              WRITE(iufilmu_q, '(32x, 3E18.8)') mob(3, 1, imode), mob(3, 2, imode), mob(3, 3, imode)
            ENDDO
          ENDDO
          IF (iqq == totq) CLOSE(iufilmu_q)
        ENDIF ! my_pool_id
        IF (iqq == totq) THEN
          DEALLOCATE(vkk_all, STAT = ierr)
          IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating vkk_all', 1)
          DEALLOCATE(carrier_density, STAT = ierr)
          IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating carrier_density', 1)
        ENDIF  ! iqq
      ENDIF ! iverbosity
      ! SP - IBTE only with if EPW compiled with MPI
      IF (SUM(ind) > 0) THEN
        !
        IF (ionode) ind_tot = ind_tot + SUM(ind)
#if defined(__MPI)
        CALL MPI_BCAST(ind_tot, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
#endif
!       WRITE(stdout,'(a,i9,E22.8)') '     Total number of element written ',ind_tot
        IF (ind(my_pool_id + 1) > 0) THEN
          WRITE(iunepmat) trans_prob(1:ind(my_pool_id + 1))
          ! The flush is crucial otherwise restart wont work correctly.
          FLUSH(iunepmat)
          DO ifil = 1, ind(my_pool_id + 1)
            WRITE(iunsparseq) sparse_q(ifil)
            WRITE(iunsparseq) sparse_k(ifil)
            WRITE(iunsparseq) sparse_i(ifil)
            WRITE(iunsparseq) sparse_j(ifil)
            WRITE(iunsparseq) sparse_t(ifil)
          ENDDO
          FLUSH(iunsparseq)
        ENDIF
        !
        ! Offset for the next q iteration
        lrepmatw2_merge = lrepmatw2_merge + ind(my_pool_id + 1)
      ENDIF
      IF (SUM(indcb) > 0) THEN
        !
        IF (ionode) ind_totcb = ind_totcb + SUM(indcb)
#if defined(__MPI)
        CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
#endif
        !
        IF (indcb(my_pool_id + 1) > 0) THEN
          WRITE(iunepmatcb) trans_probcb(1:indcb(my_pool_id + 1))
          FLUSH(iunepmatcb)
          DO ifil = 1, indcb(my_pool_id + 1)
            WRITE(iunsparseqcb) sparsecb_q(ifil)
            WRITE(iunsparseqcb) sparsecb_k(ifil)
            WRITE(iunsparseqcb) sparsecb_i(ifil)
            WRITE(iunsparseqcb) sparsecb_j(ifil)
            WRITE(iunsparseqcb) sparsecb_t(ifil)
          ENDDO
          FLUSH(iunsparseqcb)
        ENDIF
        !
        ! Offset for the next q iteration
        lrepmatw5_merge = lrepmatw5_merge + indcb(my_pool_id + 1)
        !
      ENDIF ! indcb
      !
      is_restart_loop = (MOD(iqq, restart_step) == 0)  .OR. (iqq == totq)
      !
      IF (is_restart_loop) THEN
        ! Save to file restart information in formatted way for possible restart
        lrepmatw2_restart(:) = 0
        lrepmatw5_restart(:) = 0
        lrepmatw2_restart(my_pool_id + 1) = lrepmatw2_merge
        lrepmatw5_restart(my_pool_id + 1) = lrepmatw5_merge
        CALL mp_sum(lrepmatw2_restart, inter_pool_comm)
        CALL mp_sum(lrepmatw5_restart, inter_pool_comm)
        !
        ! Scattering rates
        inv_tau_all_MPI = inv_tau_all
        inv_tau_allcb_MPI = inv_tau_allcb
        CALL mp_sum(inv_tau_all_MPI, inter_pool_comm)
        CALL mp_sum(inv_tau_allcb_MPI, inter_pool_comm)
        IF (iverbosity == 3) THEN
          ! Scattering rates (modes)
          inv_tau_all_mode_MPI = inv_tau_all_mode
          inv_tau_allcb_mode_MPI = inv_tau_allcb_mode
          CALL mp_sum(inv_tau_all_mode_MPI, inter_pool_comm)
          CALL mp_sum(inv_tau_allcb_mode_MPI, inter_pool_comm)
          ! Scattering rates (freq)
          inv_tau_all_freq_MPI = inv_tau_all_freq
          inv_tau_allcb_freq_MPI = inv_tau_allcb_freq
          CALL mp_sum(inv_tau_all_freq_MPI, inter_pool_comm)
          CALL mp_sum(inv_tau_allcb_freq_MPI, inter_pool_comm)
        ENDIF ! iverbosity
        !
        IF (ionode) THEN
          INQUIRE(FILE = 'restart' // '_' // TRIM(my_image_id_ch) // '.fmt', EXIST = exst)
          IF (exst) ios = f_copy('restart' // '_' // TRIM(my_image_id_ch) // '.fmt', &
            'restart' // '_' // TRIM(my_image_id_ch) // '.fmt.bak')
          !
          INQUIRE(FILE = 'inv_tau_tmp' // '_' // TRIM(my_image_id_ch), EXIST = exst)
          IF (exst) ios = f_copy('inv_tau_tmp' // '_' // TRIM(my_image_id_ch), &
            'inv_tau_tmp' // '_' // TRIM(my_image_id_ch) // '.bak')
          !
          INQUIRE(FILE = 'inv_taucb_tmp' // '_' // TRIM(my_image_id_ch), EXIST = exst)
          IF (exst) ios = f_copy('inv_taucb_tmp' // '_' // TRIM(my_image_id_ch), &
            'inv_taucb_tmp' // '_' // TRIM(my_image_id_ch) // '.bak')
          !
          OPEN(UNIT = iunrestart, FILE = 'restart' // '_' // TRIM(my_image_id_ch) // '.fmt')
          WRITE(iunrestart, *) iqq
          WRITE(iunrestart, *) ind_tot
          WRITE(iunrestart, *) ind_totcb
          WRITE(iunrestart, *) npool
          WRITE(iunrestart, *) nimage
          DO ipool = 1, npool
            WRITE(iunrestart, *) lrepmatw2_restart(ipool)
          ENDDO
          DO ipool = 1, npool
            WRITE(iunrestart, *) lrepmatw5_restart(ipool)
          ENDDO
          CLOSE(iunrestart)
          !
          ! Scattering rates
          OPEN(UNIT = iuntau, FORM = 'unformatted', FILE = 'inv_tau_tmp' // '_' // TRIM(my_image_id_ch))
          WRITE(iuntau) inv_tau_all_MPI
          CLOSE(iuntau)
          !
          OPEN(UNIT = iuntaucb, FORM = 'unformatted', FILE = 'inv_taucb_tmp' // '_' // TRIM(my_image_id_ch))
          WRITE(iuntaucb) inv_tau_allcb_MPI
          CLOSE(iuntaucb)
          IF (iverbosity == 3) THEN
            ! Scattering rates (modes)
            OPEN(UNIT = iuntau, FORM = 'unformatted', FILE = 'inv_tau_mode_tmp' // '_' // TRIM(my_image_id_ch))
            WRITE(iuntau) inv_tau_all_mode_MPI
            CLOSE(iuntau)
            !
            OPEN(UNIT = iuntaucb, FORM = 'unformatted', FILE = 'inv_taucb_mode_tmp' // '_' // TRIM(my_image_id_ch))
            WRITE(iuntaucb) inv_tau_allcb_mode_MPI
            CLOSE(iuntaucb)
            ! Scattering rates (freq)
            OPEN(UNIT = iuntau, FORM = 'unformatted', FILE = 'inv_tau_freq_tmp' // '_' // TRIM(my_image_id_ch))
            WRITE(iuntau) inv_tau_all_freq_MPI
            CLOSE(iuntau)
            !
            OPEN(UNIT = iuntaucb, FORM = 'unformatted', FILE = 'inv_taucb_freq_tmp' // '_' // TRIM(my_image_id_ch))
            WRITE(iuntaucb) inv_tau_allcb_freq_MPI
            CLOSE(iuntaucb)
          ENDIF ! iverbosity
        ENDIF ! ionode
        !
      ENDIF ! is_restart_loop
    ENDIF ! first_cycle
    !
    IF (iqq == totq) THEN
      IF (gb_scattering) THEN
        DEALLOCATE(inv_tau_gb, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating vkk_all', 1)
      ENDIF
      IF (gb_scattering .AND. iverbosity /= 3) THEN
        DEALLOCATE(vkk_all, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating vkk_all', 1)
      ENDIF
      ALLOCATE(vkk_all(3, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('print_ibte', 'Error allocating vkk_all', 1)
      !
      vkk_all(:, :, :) = zero
      wkf_all(:) = zero
      ! Computes the k-velocity
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        !
        wkf_all(ik + lower_bnd -1 ) = wkf(ikk)
        !
        DO ibnd = 1, nbndfst
          vkk_all(:, ibnd, ik + lower_bnd - 1) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
          etf_all(ibnd, ik + lower_bnd - 1) = etf(ibndmin - 1 + ibnd, ikk)
        ENDDO
      ENDDO
      CALL mp_sum(vkk_all, inter_pool_comm)
      CALL mp_sum(etf_all, inter_pool_comm)
      CALL mp_sum(wkf_all, inter_pool_comm)
      CALL mp_sum(inv_tau_all, inter_pool_comm)
      CALL mp_sum(inv_tau_allcb, inter_pool_comm)
      CALL mp_sum(inv_tau_all, inter_image_comm)
      CALL mp_sum(inv_tau_allcb, inter_image_comm)
      IF (iverbosity == 3) THEN
        CALL mp_sum(inv_tau_all_mode, inter_pool_comm)
        CALL mp_sum(inv_tau_allcb_mode, inter_pool_comm)
        CALL mp_sum(inv_tau_all_freq, inter_pool_comm)
        CALL mp_sum(inv_tau_allcb_freq, inter_pool_comm)
        ! Image sum after pool sum
        CALL mp_sum(inv_tau_all_mode, inter_image_comm)
        CALL mp_sum(inv_tau_allcb_mode, inter_image_comm)
        CALL mp_sum(inv_tau_all_freq, inter_image_comm)
        CALL mp_sum(inv_tau_allcb_freq, inter_image_comm)
      ENDIF
      ! Now write total number of q-point inside and k-velocity
      IF (ionode) THEN
        ! Now write total number of q-point inside and k-velocity
        OPEN(iufilibtev_sup, FILE = 'IBTEvel_sup' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
        WRITE(iufilibtev_sup, '(a)') '# Number of elements in hole and electrons'
        WRITE(iufilibtev_sup, '(2i16)') ind_tot, ind_totcb
        WRITE(iufilibtev_sup, '(a)') '# evbm     ecbm'
        WRITE(iufilibtev_sup, '(2E22.12)') evbm, ecbm
        WRITE(iufilibtev_sup, '(a)') '# itemp    ef0    efcb'
        DO itemp = 1, nstemp
          WRITE(iufilibtev_sup, '(i8,2E22.12)') itemp, ef0(itemp), efcb(itemp)
        ENDDO
        WRITE(iufilibtev_sup, '(a)') '# ik     ibnd                      velocity (x,y,z)                eig      weight'
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            WRITE(iufilibtev_sup, '(i8,i6,5E22.12)') ik, ibnd, vkk_all(:, ibnd, ik), etf_all(ibnd, ik), wkf_all(ik)
          ENDDO
        ENDDO
        IF (assume_metal) THEN
          DO itemp = 1, nstemp
            WRITE(iufilibtev_sup, '(i8,1E22.12)') itemp, dos(itemp)
          ENDDO
        ENDIF
        CLOSE(iufilibtev_sup)
        !
        ! Save the inv_tau and inv_tau_all on file (formatted)
        OPEN(iufilibtev_sup, FILE = 'inv_tau' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
        WRITE(iufilibtev_sup, '(a)') '# Hole relaxation time [Multiply the relaxation time by 20670.6944033 to get 1/ps]'
        WRITE(iufilibtev_sup, '(a)') '# itemp    kpt      ibnd    energy [Ry]   relaxation time [Ry]'
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              WRITE(iufilibtev_sup, '(i5,i8,i6,2E22.12)') itemp, ik, ibnd, etf_all(ibnd, ik), inv_tau_all(ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilibtev_sup)
        !
        ! Save the mode resolvedinv_tau and inv_tau_all_mode on file (formatted)
        IF (iverbosity == 3) THEN
          WRITE(stdout,'(5x," ")')
          WRITE(stdout, '(5x,a,f10.6,a)') "NOTE: Spectral decomposition of scattering rates will be computed only up to " &
              , mob_maxfreq * ryd2mev, ' meV'
          WRITE(stdout,'(5x," ")')
          OPEN(iufilibtev_sup, FILE = 'inv_tau_mode' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
          WRITE(iufilibtev_sup, '(a)') '# Hole relaxation time [Multiply the relaxation time by 20670.6944033 to get 1/ps]'
          WRITE(iufilibtev_sup, '(a)') '# itemp    kpt      ibnd     imode   energy [Ry]   relaxation time [Ry]'
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                DO imode = 1, nmodes
                  WRITE(iufilibtev_sup, '(i5,i8,i6,i6,2E22.12)') itemp, ik, ibnd, imode, etf_all(ibnd, ik), &
                                                                 inv_tau_all_mode(imode, ibnd, ik, itemp)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilibtev_sup)
        ENDIF
        !
        ! Save the inv_tau and inv_tau_all on file (formatted)
        OPEN(iufilibtev_sup, FILE = 'inv_taucb' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
        WRITE(iufilibtev_sup, '(a)') '# Electron relaxation time [Multiply the relaxation time by 20670.6944033 to get 1/ps]'
        WRITE(iufilibtev_sup, '(a)') '# itemp    kpt      ibnd    energy [Ry]   relaxation time [Ry]'
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              WRITE(iufilibtev_sup, '(i5,i8,i6,2E22.12)') itemp, ik, ibnd, etf_all(ibnd, ik), inv_tau_allcb(ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufilibtev_sup)
        !
        IF (iverbosity == 3) THEN
          OPEN(iufilibtev_sup, FILE = 'inv_taucb_mode' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
          WRITE(iufilibtev_sup, '(a)') '# Electron relaxation time [Multiply the relaxation time by 20670.6944033 to get 1/ps]'
          WRITE(iufilibtev_sup, '(a)') '# itemp    kpt      ibnd    imode    energy [Ry]   relaxation time [Ry]'
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                DO imode = 1, nmodes
                  WRITE(iufilibtev_sup, '(i5,i8,i6,i6,2E22.12)') itemp, ik, ibnd, imode, etf_all(ibnd, ik), &
                                                                 inv_tau_allcb_mode(imode, ibnd, ik, itemp)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilibtev_sup)
          !
          ! Compute and save spectral decomposition
          OPEN(iufilibtev_sup, FILE = 'inv_tau_freq' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
          WRITE(iufilibtev_sup, '(a)') '# Electron relaxation time [Multiply the relaxation time by 20670.6944033 to get 1/ps]'
          WRITE(iufilibtev_sup, '(a)') '#  kpt      ibnd    energy [Ry]  freq (meV)    relaxation time [Ry]'
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              DO ifreq = 1, mob_nfreq
                IF (inv_tau_all_freq(ifreq, ibnd, ik) > eps160) THEN
                  WRITE(iufilibtev_sup, '(i8,i6,3E22.12)') ik, ibnd, etf_all(ibnd, ik), ifreq * step * ryd2mev, &
                                                               inv_tau_all_freq(ifreq, ibnd, ik)
                ENDIF
              ENDDO ! ifreq
            ENDDO ! ibnd
          ENDDO ! ik
          CLOSE(iufilibtev_sup)
          OPEN(iufilibtev_sup, FILE = 'inv_taucb_freq' // '_' // TRIM(my_image_id_ch) // '.fmt', FORM = 'formatted')
          WRITE(iufilibtev_sup, '(a)') '# Electron relaxation time [Multiply the relaxation time by 20670.6944033 to get 1/ps]'
          WRITE(iufilibtev_sup, '(a)') '#  kpt      ibnd    energy [Ry]  freq (meV)    relaxation time [Ry]'
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              DO ifreq = 1, mob_nfreq
                IF (inv_tau_allcb_freq(ifreq, ibnd, ik) > eps160) THEN
                  WRITE(iufilibtev_sup, '(i8,i6,3E22.12)') ik, ibnd, etf_all(ibnd, ik), ifreq * step * ryd2mev, &
                                                               inv_tau_allcb_freq(ifreq, ibnd, ik)
                ENDIF
              ENDDO ! ifreq
            ENDDO ! ibnd
          ENDDO ! ik
          CLOSE(iufilibtev_sup)
        ENDIF ! iverbosity
        !
      ENDIF ! master
      !
      ! Now print the carrier density for checking (for non-metals)
      IF (.NOT. assume_metal) THEN
        ALLOCATE(carrier_density(nstemp), STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error allocating carrier_density', 1)
        carrier_density(:) = zero
        IF (ncarrier < 0.0) THEN ! VB
          CALL carr_density(carrier_density, etf_all, wkf_all, ef0)
          IF (system_2d == 'no') THEN
            carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3.0d0)
          ELSE
            carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-2.0d0)
          ENDIF
          DO itemp = 1, nstemp
            etemp = gtemp(itemp)
            WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6)') etemp * ryd2ev / kelvin2eV, ef0(itemp) * ryd2ev,  carrier_density(itemp)
          ENDDO
        ELSE ! CB
          CALL carr_density(carrier_density, etf_all, wkf_all, efcb)
          IF (system_2d == 'no') THEN
            carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-3.0d0)
          ELSE
            carrier_density = carrier_density * inv_cell * (bohr2ang * ang2cm)**(-2.0d0)
          ENDIF
          DO itemp = 1, nstemp
            etemp = gtemp(itemp)
            WRITE(stdout, '(5x, 1f8.3, 1f12.4, 1E19.6)') etemp * ryd2ev / kelvin2eV, efcb(itemp) * ryd2ev,  carrier_density(itemp)
          ENDDO
        ENDIF ! ncarrier
        DEALLOCATE(carrier_density, STAT = ierr)
        IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating carrier_density', 1)
      ENDIF ! assume_metal
      DEALLOCATE(vkk_all, STAT = ierr)
      IF (ierr /= 0) CALL errore('print_ibte', 'Error deallocating vkk_all', 1)
    ENDIF ! iqq
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE print_ibte
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE carr_density(carrier_density, etf_all, wkf_all, ef0)
    !-----------------------------------------------------------------------
    !!
    !! This routine computes the carrier density (in a.u.)
    !!
    USE kinds,         ONLY : DP
    USE input,         ONLY : ncarrier, nstemp
    USE global_var,    ONLY : nbndfst, gtemp, nktotf, lower_bnd, nkf, evbm, ecbm
    USE ep_constants,  ONLY : eps10
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : world_comm, inter_pool_comm
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(inout) :: carrier_density(nstemp)
    !! Carrier density
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weights from all the cores
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    !
    ! Local variables
    INTEGER :: itemp
    !! temperature index
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    REAL(KIND = DP) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: fnk
    !! Fermi-Dirac occupation function
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Here computes Fermi-Dirac
    !
    carrier_density(:) = 0.0
    DO itemp = 1, nstemp
      etemp = gtemp(itemp)
      IF (ncarrier < 0.0) THEN ! VB
        DO ik = 1, nkf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for hole conduction
            IF (etf_all(ibnd, ik + lower_bnd - 1) < (evbm + eps10)) THEN
              ! Energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik + lower_bnd - 1) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density(itemp) = carrier_density(itemp) + wkf_all(ik + lower_bnd - 1) * (1.0d0 - fnk)
            ENDIF
          ENDDO
        ENDDO
      ELSE ! CB
        DO ik = 1, nkf
          DO ibnd = 1, nbndfst
            ! This selects only valence bands for hole conduction
            IF (etf_all(ibnd, ik + lower_bnd - 1) > (ecbm - eps10)) THEN
              !  energy at k (relative to Ef)
              ekk = etf_all(ibnd, ik + lower_bnd - 1) - ef0(itemp)
              fnk = wgauss(-ekk / etemp, -99)
              ! The wkf(ikk) already include a factor 2
              carrier_density(itemp) = carrier_density(itemp) + wkf_all(ik + lower_bnd - 1) *  fnk
            ENDIF
          ENDDO
        ENDDO
      ENDIF ! ncarrier
    ENDDO
    CALL mp_sum(carrier_density, inter_pool_comm)
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE carr_density
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE fin_write(iter, F_in, av_mob_old, elec)
    !----------------------------------------------------------------------------
    !!
    !! Writes the F without magnetic field for restart
    !!
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilFi_all
    USE io_files,      ONLY : diropn
    USE input,         ONLY : nstemp
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE global_var,    ONLY : nktotf, nbndfst
    USE ep_constants,      ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iter
    !! Iteration number
    REAL(KIND = DP), INTENT(in) :: f_in(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i
    REAL(KIND = DP), INTENT(in) :: av_mob_old(nstemp)
    !! Error in the hole mobility
    LOGICAL, INTENT(in) :: elec
    !! IF true we do electron mobility, if false the hole one.
    !
    ! Local variable
    LOGICAL :: exst
    !! File exist
    INTEGER :: i
    !! Running index for the vector
    INTEGER :: lfi_all
    !! Length of the vector
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: idir
    !! Direction index
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
    exst = .FALSE.
    aux(:) = zero
    IF (mpime == ionode_id) THEN
      !
      lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
      ! First element is the iteration number
      aux(1) = iter
      !
      i = 1
      DO itemp = 1, nstemp
        i = i + 1
        ! Value of the previous h mobility (used for error evaluation)
        aux(i) = av_mob_old(itemp)
      ENDDO
      !
      i = 1 + nstemp
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            DO idir = 1, 3
              i = i +1
              aux(i) = f_in(idir, ibnd, ik, itemp)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      ! Electron mobility
      IF (elec) THEN
        CALL diropn(iufilFi_all, 'Fin_restartcb', lfi_all, exst)
        CALL davcio(aux, lfi_all, iufilFi_all, 1, +1)
      ELSE
        CALL diropn(iufilFi_all, 'Fin_restart', lfi_all, exst)
        CALL davcio( aux, lfi_all, iufilFi_all, 1, +1)
      ENDIF
      CLOSE(iufilFi_all)
      !
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE fin_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE fin_read(iter, F_in, av_mob_old, elec)
    !----------------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilFi_all
    USE input,     ONLY : nstemp
    USE ep_constants,      ONLY : zero
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    USE global_var,ONLY : nbndfst, nktotf
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iter
    !! Iteration number
    REAL(KIND = DP), INTENT(inout) :: f_in(3, nbndfst, nktotf, nstemp)
    !! In solution for iteration i
    REAL(KIND = DP), INTENT(inout) :: av_mob_old(nstemp)
    !! Error in the hole mobility
    LOGICAL, INTENT(in) :: elec
    !! IF true we do electron mobility, if false the hole one.
    !
    ! Local variable
    CHARACTER(LEN = 256) :: name1
    !! Variable name
    LOGICAL :: exst
    !! Does the file exist
    INTEGER :: i
    !! Running index for the vector
    INTEGER :: lfi_all
    !! Length of the vector
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    INTEGER :: idir
    !! Direction index
    INTEGER :: itemp
    !! Temperature index
    REAL(KIND = DP) :: aux(3 * nbndfst * nktotf * nstemp + nstemp + 1)
    !! Vector to store the array
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
      IF (elec) THEN
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restartcb1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restartcb'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        !
        IF (exst) THEN ! read the file
          !
          lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
          CALL diropn(iufilFi_all, 'Fin_restartcb', lfi_all, exst)
          CALL davcio(aux, lfi_all, iufilFi_all, 1, -1)
          !
          ! First element is the iteration number
          iter = INT(aux(1))
          !
          i = 1
          DO itemp = 1, nstemp
            i = i + 1
            ! Last value of hole mobility
            av_mob_old(itemp) = aux(i)
          ENDDO
          !
          i = 1 + nstemp
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                DO idir = 1, 3
                  i = i + 1
                  f_in(idir, ibnd, ik, itemp) = aux(i)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilFi_all)
        ENDIF
      ELSE ! hole
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restart1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.Fin_restart'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        !
        IF (exst) THEN ! read the file
          !
          lfi_all = 3 * nbndfst * nktotf * nstemp + nstemp + 1
          CALL diropn(iufilFi_all, 'Fin_restart', lfi_all, exst)
          CALL davcio(aux, lfi_all, iufilFi_all, 1, -1)
          !
          ! First element is the iteration number
          iter = INT(aux(1))
          !
          i = 1
          DO itemp = 1, nstemp
            i = i + 1
            ! Last value of hole mobility
            av_mob_old(itemp) = aux(i)
          ENDDO
          !
          i = 1 + nstemp
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, (nbndfst)
                DO idir = 1, 3
                  i = i + 1
                  f_in(idir, ibnd, ik, itemp) = aux(i)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufilFi_all)
        ENDIF
      ENDIF
    ENDIF ! mpime
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iter,       ionode_id, world_comm)
      CALL mp_bcast(F_in,       ionode_id, world_comm)
      CALL mp_bcast(av_mob_old, ionode_id, world_comm)
      WRITE(stdout, '(a,i10)' ) '     Restart from iter: ', iter
    ENDIF ! exists
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE fin_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_merge()
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE io_var,           ONLY : iunepmat_merge, iunepmat, iunepmatcb_merge,              &
                                 iunepmatcb, iunsparseq_merge, iunsparsek_merge,          &
                                 iunsparsej_merge,iunsparset_merge, iunepmatcb_merge,     &
                                 iunsparseqcb_merge, iunsparsekcb_merge, iunsparsei_merge,&
                                 iunsparseicb_merge, iunsparsejcb_merge, iunsparsetcb_merge
    USE mp_global,        ONLY : my_pool_id, npool, world_comm, my_image_id, inter_pool_comm
    USE io_files,         ONLY : tmp_dir, prefix
    USE mp,               ONLY : mp_sum, mp_barrier
    USE global_var,       ONLY : lrepmatw2_merge, lrepmatw5_merge, ctype
    USE input,            ONLY : int_mob, carrier, ncarrier, assume_metal
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE,MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION,          &
                                 MPI_STATUS_IGNORE, MPI_INTEGER
#endif
    !
    IMPLICIT NONE
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! Pool number, character
    CHARACTER(LEN = 256) :: dirname(2)
    !! Name of the directory to hold files
    CHARACTER(LEN = 256) :: filename(6)
    !! Name of the files to merge files
    CHARACTER(LEN = 256) :: path_to_files(2)
    !! Name of the path to files
    INTEGER :: i2
    !! Indexes to loop over file sizes
    INTEGER :: lrepmatw2_tot(npool)
    !! Lenght of each file
    INTEGER :: lrepmatw5_tot(npool)
    !! Lenght of each file
    INTEGER :: ich
    !! Loop over directories
    INTEGER :: ifil
    !! Index over the files
    INTEGER :: io_u(6)
    !! Input output units
    INTEGER :: ierr
    !! Error status
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: lsize
    !! Size of what we write
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw
    !! Offset while writing scattering to files
    INTEGER(KIND = MPI_OFFSET_KIND) :: tmp_sum
    !! Temporary sum
#else
    INTEGER(KIND = 8) :: lsize
    !! Size of what we write
    INTEGER(KIND = 8) :: lrepmatw
    !! Offset while writing scattering to files
   INTEGER(KIND = 8) :: tmp_sum
    !! Temporary sum
#endif
    INTEGER, ALLOCATABLE :: sparse(:, :)
    !! Vaariable for reading and writing the files
    INTEGER, ALLOCATABLE :: sparsecb(:, :)
    !! Vaariable for reading and writing the files
    REAL(KIND = DP), ALLOCATABLE :: trans_prob(:)
    !! Variable for reading and writing trans_prob
    REAL(KIND = DP), ALLOCATABLE :: trans_probcb(:)
    !! Variable for reading and writing trans_prob
    CHARACTER(LEN = 256) :: my_image_id_ch
    !! image id for writing 
    !
    WRITE(my_image_id_ch, "(I0)") my_image_id
    !
    ! for metals merge like it's for holes
    IF (ctype == 0 .OR. ctype == -1) THEN
      !
      ALLOCATE(trans_prob(lrepmatw2_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error allocating trans_prob', 1)
      ALLOCATE(sparse(5, lrepmatw2_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error allocating sparse', 1)
      !
      io_u(1) = iunepmat_merge
      io_u(2) = iunsparseq_merge
      io_u(3) = iunsparsek_merge
      io_u(4) = iunsparsei_merge
      io_u(5) = iunsparsej_merge
      io_u(6) = iunsparset_merge
      !
      dirname(1) = 'Fepmatkq1' // '_' // TRIM(my_image_id_ch)
      dirname(2) = 'Fsparse' // '_' // TRIM(my_image_id_ch)
      !
      filename(1) = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_image_id_ch)
      filename(2) = 'sparseq' // '_' // TRIM(my_image_id_ch)
      filename(3) = 'sparsek' // '_' // TRIM(my_image_id_ch) 
      filename(4) = 'sparsei' // '_' // TRIM(my_image_id_ch)
      filename(5) = 'sparsej' // '_' // TRIM(my_image_id_ch)
      filename(6) = 'sparset' // '_' // TRIM(my_image_id_ch)
      !
      path_to_files(1) = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkq1' // '_'
      path_to_files(2) = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_'
      !
      lrepmatw2_tot = 0
      lrepmatw2_tot(my_pool_id + 1) = lrepmatw2_merge
      CALL mp_sum(lrepmatw2_tot, inter_pool_comm)
#if defined(__MPI)
      DO ich = 1, 6
        CALL MPI_FILE_OPEN(inter_pool_comm, filename(ich), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, io_u(ich), ierr)
      ENDDO
#else
      OPEN(UNIT = io_u(1), FILE = filename(1), IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = lrepmatw2_merge * 8)
      DO ich = 2, 6
        OPEN(UNIT = io_u(ich), FILE = filename(ich), IOSTAT = ierr, FORM = 'unformatted', &
             STATUS = 'unknown', ACCESS = 'direct', RECL = lrepmatw2_merge * 4)
      ENDDO
#endif
      IF (ierr /= 0) CALL errore('iter_merge', 'Error in opening .epmatkq1 or .sparseX file', 1)

      !
      DO ich = 1, 2
        ! Read files per processor
        WRITE(my_pool_id_ch, "(I0)") my_pool_id
        filint = TRIM(path_to_files(ich)) // TRIM(my_pool_id_ch)
        OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'old', FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
        IF (ich == 1) THEN
          DO i2 = 1, lrepmatw2_merge
            READ(iunepmat) trans_prob(i2)
          ENDDO
        ELSE
          DO i2 = 1, lrepmatw2_merge
            DO ifil = 1, 5
              READ(iunepmat) sparse(ifil, i2)
            ENDDO
          ENDDO
        ENDIF
        CLOSE(iunepmat, STATUS = 'keep')
        IF (ich == 1) THEN
#if defined(__MPI)
          tmp_sum = 0
          DO i2 = 1, my_pool_id + 1
            tmp_sum = tmp_sum + lrepmatw2_tot(i2)
          ENDDO
          lrepmatw = INT(tmp_sum - lrepmatw2_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND
          lsize = INT(lrepmatw2_merge, KIND = MPI_OFFSET_KIND)
          CALL MPI_FILE_WRITE_AT(io_u(1), lrepmatw, trans_prob(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
#else
          WRITE(UNIT = io_u(1), REC = 1, IOSTAT = ierr) trans_prob
#endif
          IF (ierr /= 0) CALL errore('iter_merge', 'Error in writing .epmatkq1 file', 1)
        ELSE
          DO ifil = 1, 5
#if defined(__MPI)
            tmp_sum = 0
            DO i2 = 1, my_pool_id + 1
              tmp_sum = tmp_sum + lrepmatw2_tot(i2)
            ENDDO
            lrepmatw = INT(tmp_sum - lrepmatw2_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND
            lsize = INT(lrepmatw2_merge, KIND = MPI_OFFSET_KIND)
            CALL MPI_FILE_WRITE_AT(io_u(ifil + 1), lrepmatw, sparse(ifil, :), lsize, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#else
            WRITE(UNIT = io_u(ifil + 1), REC = 1, IOSTAT = ierr) sparse(ifil, :)
#endif
            IF (ierr /= 0) CALL errore('iter_merge', 'Error in writing .sparseX file', 1)
          ENDDO
        ENDIF
      ENDDO
      !
      DO ich = 1, 6
#if defined(__MPI)
        CALL MPI_FILE_CLOSE(io_u(ich), ierr)
#else
        CLOSE(io_u(ich), STATUS = 'keep')
#endif
      ENDDO
      !
      DEALLOCATE(trans_prob, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error deallocating trans_prob', 1)
      DEALLOCATE(sparse, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error deallocating sparse', 1)
      !
    ENDIF
    IF (ctype == 0 .OR. ctype == 1) THEN
      !
      ALLOCATE(trans_probcb(lrepmatw5_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error allocating trans_probcb', 1)
      ALLOCATE(sparsecb(5, lrepmatw5_merge), STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error allocating sparsecb', 1)
      !
      io_u(1) = iunepmatcb_merge
      io_u(2) = iunsparseqcb_merge
      io_u(3) = iunsparsekcb_merge
      io_u(4) = iunsparseicb_merge
      io_u(5) = iunsparsejcb_merge
      io_u(6) = iunsparsetcb_merge
      !
      dirname(1) = 'Fepmatkqcb1' // '_' // TRIM(my_image_id_ch)
      dirname(2) = 'Fsparsecb' // '_' // TRIM(my_image_id_ch)
      !
      filename(1) = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkqcb1' // '_' // TRIM(my_image_id_ch)
      filename(2) = 'sparseqcb' // '_' // TRIM(my_image_id_ch)
      filename(3) = 'sparsekcb' // '_' // TRIM(my_image_id_ch)
      filename(4) = 'sparseicb' // '_' // TRIM(my_image_id_ch)
      filename(5) = 'sparsejcb' // '_' // TRIM(my_image_id_ch)
      filename(6) = 'sparsetcb' // '_' // TRIM(my_image_id_ch)
      !
      path_to_files(1) = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkqcb1' // '_'
      path_to_files(2) = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparsecb' // '_'
      !
      lrepmatw5_tot = 0
      lrepmatw5_tot(my_pool_id + 1) = lrepmatw5_merge
      CALL mp_sum(lrepmatw5_tot, inter_pool_comm)
#if defined(__MPI)
      DO ich = 1, 6
        CALL MPI_FILE_OPEN(inter_pool_comm, filename(ich), MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, io_u(ich), ierr)
      ENDDO
#else
      OPEN(UNIT = io_u(1), FILE = filename(1), IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = lrepmatw5_merge * 8)
      DO ich = 2, 6
        OPEN(UNIT = io_u(ich), FILE = filename(ich), IOSTAT = ierr, FORM = 'unformatted', &
             STATUS = 'unknown', ACCESS = 'direct', RECL = lrepmatw5_merge * 4)
      ENDDO
#endif
      !
      DO ich = 1, 2
        ! Read files per processor
        WRITE(my_pool_id_ch, "(I0)") my_pool_id
        filint = TRIM(path_to_files(ich)) // TRIM(my_pool_id_ch)
        OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
        IF (ich == 1) THEN
          DO i2 = 1, lrepmatw5_merge
            READ(iunepmatcb) trans_probcb(i2)
          ENDDO
        ELSE
          DO i2 = 1, lrepmatw5_merge
            DO ifil = 1 ,5
              READ(iunepmatcb) sparsecb(ifil, i2)
            ENDDO
          ENDDO
        ENDIF
        CLOSE(iunepmatcb, STATUS = 'keep')
        IF (ich == 1) THEN
#if defined(__MPI)
          tmp_sum = 0
          DO i2 = 1, my_pool_id + 1
            tmp_sum = tmp_sum + lrepmatw5_tot(i2)
          ENDDO
          lrepmatw = INT(tmp_sum - lrepmatw5_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND
          lsize = INT(lrepmatw5_merge, KIND = MPI_OFFSET_KIND)
          CALL MPI_FILE_WRITE_AT(io_u(1), lrepmatw, trans_probcb(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
#else
          WRITE(UNIT = io_u(1), REC = 1, IOSTAT = ierr) trans_probcb
#endif
          IF (ierr /= 0) CALL errore('iter_merge', 'Error in writing .epmatkqcb1 file', 1)
        ELSE
          DO ifil = 1, 5
#if defined(__MPI)
            tmp_sum = 0
            DO i2 = 1, my_pool_id + 1
              tmp_sum = tmp_sum + lrepmatw5_tot(i2)
            ENDDO
            lrepmatw = INT(tmp_sum - lrepmatw5_tot(my_pool_id + 1), KIND = MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND
            lsize = INT(lrepmatw5_merge, KIND = MPI_OFFSET_KIND)
            CALL MPI_FILE_WRITE_AT(io_u(ifil + 1), lrepmatw, sparsecb(ifil, :), lsize, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#else
            WRITE(UNIT = io_u(ifil + 1), REC = 1, IOSTAT = ierr) sparsecb(ifil, :)
#endif
            IF (ierr /= 0) CALL errore('iter_merge', 'Error in writing .sparsecbX file', 1)
          ENDDO
        ENDIF
      ENDDO
      !
      DO ich = 1, 6
#if defined(__MPI)
        CALL MPI_FILE_CLOSE(io_u(ich), ierr)
#else
        CLOSE(io_u(ich), STATUS = 'keep')
#endif
      ENDDO
      !
      DEALLOCATE(trans_probcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error deallocating trans_probcb', 1)
      DEALLOCATE(sparsecb, STAT = ierr)
      IF (ierr /= 0) CALL errore('iter_merge', 'Error deallocating sparsecb', 1)
      !
    ENDIF ! in all other cases it is still to decide which files to open
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_merge
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_open(ind_tot, ind_totcb, lrepmatw2_restart, lrepmatw5_restart)
    !----------------------------------------------------------------------------
    !
    ! This routine opens all the files needed to save scattering rates for the IBTE.
    !
    USE kinds,            ONLY : DP, i8b
    USE io_files,         ONLY : prefix, create_directory, delete_if_present
    USE io_var,           ONLY : iunepmat, iunsparseq,              &
                                 iunsparseqcb, iunepmatcb, iunrestart
    USE mp_global,        ONLY : world_comm, my_pool_id, npool, inter_pool_comm
    USE mp_images,        ONLY : my_image_id
    USE mp,               ONLY : mp_barrier, mp_bcast
    USE global_var,       ONLY : lrepmatw2_merge, lrepmatw5_merge, ctype
    USE input,            ONLY : int_mob, carrier, ncarrier, assume_metal
    USE io_global,        ONLY : ionode_id, ionode
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_MODE_WRONLY, MPI_MODE_CREATE, MPI_INFO_NULL, &
                                 MPI_OFFSET_KIND
#endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! To restart opening files
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! To restart opening files
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of component for the conduction band
#else
    INTEGER(KIND = 8), INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER(KIND = 8), INTENT(inout) :: ind_totcb
    !! Total number of component for conduction band
#endif
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    CHARACTER(LEN = 256) :: my_pool_id_ch
    !! my_pool_id in character type
    CHARACTER(LEN = 256) :: dirname(2), dirnamecb(2)
    !! Name of the directory to hold files
    LOGICAL :: exst
    !! Logical for existence of files
    LOGICAL :: exst2
    !! Logical for existence of files
    INTEGER :: dummy_int
    !! Dummy INTEGER for reading
    INTEGER :: ipool
    !! Pool index
    INTEGER(KIND = 8) :: position_byte
    !! Position in the file in byte
    REAL(KIND = DP) :: dummy_real
    !! Dummy variable for reading
    CHARACTER(LEN = 256) :: my_image_id_ch
    !! Image id
    !
    WRITE(my_pool_id_ch, "(I0)") my_pool_id
    !
    WRITE(my_image_id_ch, "(I0)") my_image_id
    !
    dirname(1)   = 'Fepmatkq1' // '_' // TRIM(my_image_id_ch)
    dirname(2)   = 'Fsparse' // '_' // TRIM(my_image_id_ch)
    dirnamecb(1) = 'Fepmatkqcb1' // '_' // TRIM(my_image_id_ch)
    dirnamecb(2) = 'Fsparsecb' // '_' // TRIM(my_image_id_ch)
    !
    INQUIRE(FILE = 'restart' // '_' // TRIM(my_image_id_ch) // '.fmt', EXIST = exst)
    !
    IF (ionode) THEN
      IF (exst) THEN
        OPEN(UNIT = iunrestart, FILE = 'restart' // '_' // TRIM(my_image_id_ch) // '.fmt', STATUS = 'old')
        READ(iunrestart, *)
        READ(iunrestart, *)
        READ(iunrestart, *)
        READ(iunrestart, *)
        READ(iunrestart, *)
        DO ipool = 1, npool
          READ(iunrestart, *) lrepmatw2_restart(ipool)
        ENDDO
        DO ipool = 1, npool
          READ(iunrestart, *) lrepmatw5_restart(ipool)
        ENDDO
        CLOSE(iunrestart)
      ENDIF
    ENDIF
    CALL mp_bcast(exst, ionode_id, inter_pool_comm )
    CALL mp_bcast(lrepmatw2_restart, ionode_id, inter_pool_comm )
    CALL mp_bcast(lrepmatw5_restart, ionode_id, inter_pool_comm )
    !
    ! The restart.fmt exist - we try to restart
    IF (exst) THEN
      ! Hole (or metals)
      IF (ctype == 0 .OR. ctype == -1) THEN
        !
        filint = './' // ADJUSTL(TRIM(dirname(1))) // '/'//TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_pool_id_ch)
        !
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw2_restart(my_pool_id + 1) - 1) * 8 + 1
            READ(iunepmat, POS=position_byte) dummy_real
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart.fmt is present but not the Fepmatkq1 folder', 1)
        ENDIF
        !
        filint = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          IF (lrepmatw2_restart(my_pool_id + 1) > 0) THEN
            position_byte = (5 * lrepmatw2_restart(my_pool_id + 1) - 1) * 4 + 1
            READ(iunsparseq, POS = position_byte) dummy_int
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart.fmt is present but not the Fsparse folder', 1)
        ENDIF
        !
      ENDIF ! Hole
      ! Electron
      IF (ctype == 0 .OR. ctype == 1) THEN
        !
        filint = './' // ADJUSTL(TRIM(dirnamecb(1))) // '/' // TRIM(prefix) // '.epmatkqcb1'//'_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          ! This is done to move the pointer to the right position after a restart (the position is in byte)
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (lrepmatw5_restart(my_pool_id + 1) - 1) * 8 + 1
            READ(iunepmatcb, POS = position_byte) dummy_real
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart.fmt is present but not the Fepmatkqcb1 folder', 1)
        ENDIF
        !
        filint = './' // ADJUSTL(TRIM(dirnamecb(2))) // '/' // 'sparsecb' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'old', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'rewind', ACTION = 'readwrite')
          IF (lrepmatw5_restart(my_pool_id + 1) > 0) THEN
            position_byte = (5 * lrepmatw5_restart(my_pool_id + 1) - 1) * 4 + 1
            READ(iunsparseqcb, POS = position_byte) dummy_int
          ENDIF
        ELSE
          CALL errore('iter_open', 'A restart.fmt is present but not the Fsparse folder', 1)
        ENDIF
        !
      ENDIF ! electron
      lrepmatw2_merge = lrepmatw2_restart(my_pool_id + 1)
      lrepmatw5_merge = lrepmatw5_restart(my_pool_id + 1)
      !
    ELSE ! no restart file present
      ! Hole or metals
      IF (ctype == 0 .OR. ctype == -1) THEN
        !
        CALL create_directory(ADJUSTL(TRIM(dirname(1))))
        CALL create_directory(ADJUSTL(TRIM(dirname(2))))
        !
        filint = './' // ADJUSTL(TRIM(dirname(1))) // '/' // TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint, .TRUE.)
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunepmat, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        !
        filint = './' // ADJUSTL(TRIM(dirname(2))) // '/' // 'sparse' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint, .TRUE.)
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunsparseq, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        !
      ENDIF ! Hole
      ! Electron
      IF (ctype == 0 .OR. ctype == 1) THEN
        !
        CALL create_directory(ADJUSTL(TRIM(dirnamecb(1))))
        CALL create_directory(ADJUSTL(TRIM(dirnamecb(2))))
        !
        filint = './' // ADJUSTL(TRIM(dirnamecb(1))) // '/' // TRIM(prefix) // '.epmatkqcb1' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint, .TRUE.)
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunepmatcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        !
        filint = './' // ADJUSTL(TRIM(dirnamecb(2))) // '/' // 'sparsecb' // '_' // TRIM(my_pool_id_ch)
        INQUIRE(FILE = filint, EXIST = exst2)
        !
        IF (exst2) THEN
          ! The file should not exist, we remove it
          CALL delete_if_present(filint, .TRUE.)
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ELSE
          OPEN(UNIT = iunsparseqcb, FILE = filint, STATUS = 'new', FORM = 'unformatted', &
               ACCESS = 'stream', POSITION = 'append', ACTION = 'write')
        ENDIF
        !
      ENDIF !electron
      lrepmatw2_merge = 0
      lrepmatw5_merge = 0
      !
    ENDIF ! restart
    !
    ind_tot   = 0
    ind_totcb = 0
    lrepmatw2_restart(:) = 0
    lrepmatw5_restart(:) = 0
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_open
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE scattering_write(itemp, etemp, ef0, etf_all)
    !----------------------------------------------------------------------------
    !!
    !! Write scattering rates
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilscatt_rate
    USE global_var,ONLY : ibndmin, nkqtotf, inv_tau_all, nbndfst, nktotf
    USE input,     ONLY : nbndsub, nstemp
    USE ep_constants,      ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                              meV2invps, eps4
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Temperature index
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndsub, nkqtotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file.
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: temp
    !! Temporary file name used to write scattering rate to file.
    !
    WRITE(stdout, '(/5x,"Writing scattering rate to file"/)')
    !
    IF (mpime == ionode_id) THEN
      !
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1,'(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1,'(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1,'(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate, FILE = name1, FORM = 'formatted')
      WRITE(iufilscatt_rate, '(a)') '# Inverse scattering time (ps)'
      WRITE(iufilscatt_rate, '(a)') '#      ik       ibnd                 E(ibnd)    scattering rate(1/ps)'
      !
      DO ik = 1, nktotf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        DO ibnd = 1, nbndfst
          !
          ! note that ekk does not depend on q
          ekk = etf_all(ibndmin - 1 + ibnd, ikk) - ef0(itemp)
          !
          WRITE(iufilscatt_rate, '(i9,2x)', ADVANCE = 'no') ik
          WRITE(iufilscatt_rate, '(i9,2x)', ADVANCE = 'no') ibndmin - 1 + ibnd
          WRITE(iufilscatt_rate, '(E22.14)', ADVANCE = 'no') ryd2ev * ekk
          WRITE(iufilscatt_rate, '(E26.16E3)') ryd2mev * meV2invps * inv_tau_all(itemp, ibnd, ik)
          !
        ENDDO
        !
      ENDDO
      !
      CLOSE(iufilscatt_rate)
    ENDIF
    !CALL mp_barrier(inter_pool_comm)
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE scattering_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE scattering_read(etemp, ef0, etf_all, inv_tau_all)
    !----------------------------------------------------------------------------
    !!
    !! Read scattering files
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE io_var,    ONLY : iufilscatt_rate
    USE global_var,ONLY : ibndmin, nktotf, nbndfst
    USE input,     ONLY : nbndsub, nstemp
    USE ep_constants,      ONLY : ryd2mev, kelvin2eV, ryd2ev, &
                              meV2invps, eps4
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: etemp
    !! Temperature in Ry (this includes division by kb)
    REAL(KIND = DP), INTENT(in) :: ef0
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(out) :: etf_all(nbndsub, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), INTENT(out) :: inv_tau_all(nstemp, nbndfst, nktotf)
    !! Inverse scattering rates
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name used to write scattering rates to file.
    CHARACTER(LEN = 256) :: dummy1
    !! Dummy variable to store the text of the scattering_rate file
    INTEGER :: ik
    !! K-point index
    INTEGER :: ik_tmp
    !! K-point index read from file
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ibnd_tmp
    !! Local band index read from file
    INTEGER :: ios
    !! Status of reading file
    REAL(KIND = DP) :: temp
    !! Temporary file name used to write scattering rate to file.
    !
    WRITE(stdout,'(/5x,"Reading scattering rate from file"/)')
    !
    IF (mpime == ionode_id) THEN
      ! Write to file
      temp = etemp * ryd2ev / kelvin2eV
      IF (temp < 10.d0 - eps4) THEN
        WRITE(name1, '(a18,f4.2)') 'scattering_rate_00', temp
      ELSEIF (temp >= 10.d0 - eps4 .AND. temp < 100.d0 -eps4) THEN
        WRITE(name1, '(a17,f5.2)') 'scattering_rate_0', temp
      ELSEIF (temp >= 100.d0 -eps4) THEN
        WRITE(name1, '(a16,f6.2)') 'scattering_rate_', temp
      ENDIF
      OPEN(iufilscatt_rate, FILE = name1, STATUS = 'old', IOSTAT = ios)
      WRITE(stdout,'(a16,a22)') '     Open file: ',name1
      ! There are two comment line at the beginning of the file
      READ(iufilscatt_rate, *) dummy1
      READ(iufilscatt_rate, *) dummy1
      !
      DO ik = 1, nktotf
        !
        DO ibnd = 1, nbndfst
          !
          READ(iufilscatt_rate, *) ik_tmp, ibnd_tmp, etf_all(ibndmin - 1 + ibnd, ik), inv_tau_all(1, ibnd, ik)
          inv_tau_all(1, ibnd, ik) = inv_tau_all(1, ibnd, ik) / (ryd2mev * meV2invps)
          !
          ! Check that the file corresponds to the run we are making
          IF (ABS(ibnd_tmp - ibndmin - ibnd + 1) > 0)  CALL errore('scattering_read', &
            'Band read from the scattering_rate file do not match current calculation ', 1)
          !
        ENDDO
        ! Check that the file corresponds to the run we are making
        IF (ABS(ik_tmp - ik) > 0)  CALL errore('scattering_read', &
          'k-point read from the scattering_rate file do not match current calculation ', 1)
        !
      ENDDO
      !
      etf_all = etf_all / ryd2ev
      etf_all = etf_all + ef0
      !
      CLOSE(iufilscatt_rate)
    ENDIF
    CALL mp_bcast(etf_all, ionode_id, world_comm)
    CALL mp_bcast(inv_tau_all, ionode_id, world_comm)
    !
    WRITE(stdout,'(/5x,"Scattering rate read from file"/)')
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE scattering_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE tau_write(iqq, totq, nktotf, second)
    !----------------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE input,     ONLY : nstemp
    USE io_global, ONLY : meta_ionode_id
    USE global_var,ONLY : inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : diropn
    USE mp,        ONLY : mp_barrier
    USE mp_world,  ONLY : mpime
    USE ep_constants,      ONLY : zero
    !!
    !! Write scattering rates
    !!
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iqq
    !! q-point from the selected ones within the fstick window.
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    LOGICAL, INTENT(in) :: second
    !! IF we have two Fermi level
    !
    ! Local variable
    LOGICAL :: exst
    !! Does the file exists
    INTEGER :: i
    !! Running index for the vector
    INTEGER :: itemp
    !! Running index for the temperature
    INTEGER :: ltau_all
    !! Length of the vector
    INTEGER :: ik
    !! k-point index
    INTEGER :: ibnd
    !! band index
    REAL(KIND = DP) :: aux(2 * nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array inv_tau_all and zi_all
    !
    IF (mpime == meta_ionode_id) THEN
      !
      ltau_all = 2 * nstemp * (nbndfst) * nktotf + 2
      ! First element is the iteration number
      aux(1) = REAL(iqq - 1, KIND = DP)   ! -1 because we will start at the next one.
      aux(2) = REAL(totq, KIND = DP)
      i = 2
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i + 1
            aux(i) = inv_tau_all(itemp, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            i = i +1
            aux(i) = zi_allvb(itemp, ibnd, ik)
          ENDDO
        ENDDO
      ENDDO
      CALL diropn(iufiltau_all, 'tau_restart', ltau_all, exst)
      CALL davcio(aux, ltau_all, iufiltau_all, 1, +1 )
      CLOSE(iufiltau_all)
      !
      IF (second) THEN
        ! First element is the iteration number
        aux(1) = iqq - 1   ! -1 because we will start at the next one.
        aux(2) = totq
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              aux(i) = inv_tau_allcb(itemp, ibnd, ik)
            ENDDO
          ENDDO
        ENDDO
        !
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              aux(i) = zi_allcb(itemp, ibnd, ik)
            ENDDO
          ENDDO
        ENDDO
        !
        CALL diropn(iufiltau_all, 'tau_restart_CB', ltau_all, exst)
        CALL davcio(aux, ltau_all, iufiltau_all, 1, +1)
        CLOSE(iufiltau_all)
      ENDIF
      !
    ENDIF
    !
    ! Make everythin 0 except the range of k-points we are working on
    IF (lower_bnd > 1) inv_tau_all(:, :, 1:lower_bnd - 1) = zero
    IF (upper_bnd < nktotf) inv_tau_all(:, :, upper_bnd + 1:nktotf) = zero
    IF (second) THEN
      IF (lower_bnd > 1) inv_tau_allcb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) inv_tau_allcb(:, :, upper_bnd + 1:nktotf) = zero
    ENDIF
    ! Same for the Znk factor
    IF (lower_bnd > 1) zi_allvb(:, :, 1:lower_bnd - 1) = zero
    IF (upper_bnd < nktotf) zi_allvb(:, :, upper_bnd + 1:nktotf) = zero
    IF (second) THEN
      IF (lower_bnd > 1) zi_allcb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) zi_allcb(:, :, upper_bnd + 1:nktotf) = zero
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE tau_write
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE tau_read(iqq, totq, nktotf, second)
    !----------------------------------------------------------------------------
    !!
    !! Scattering read
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, meta_ionode_id
    USE global_var,ONLY : inv_tau_all, inv_tau_allcb, zi_allvb, zi_allcb, &
                          lower_bnd, upper_bnd, nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : prefix, tmp_dir, diropn
    USE input,     ONLY : nstemp
    USE ep_constants,      ONLY : zero
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_global, ONLY : world_comm
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: iqq
    !! Current q-point from selecq.fmt
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    LOGICAL, INTENT(in) :: second
    !! IF we have two Fermi level
    !
    ! Local variables
    LOGICAL :: exst
    !! Does the file exist
    CHARACTER(LEN = 256) :: name1
    !! Name of the file
    INTEGER :: i
    !! Iterative index
    INTEGER :: itemp
    !! Iterative temperature
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ltau_all
    !! Length of the vector
    INTEGER :: nqtotf_read
    !! Total number of q-point read
    REAL(KIND = DP) :: aux(2 * nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array
    !
    IF (mpime == meta_ionode_id) THEN
      !
      ! First inquire if the file exists
#if defined(__MPI)
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart1'
#else
      name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart'
#endif
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        ltau_all = 2 * nstemp * nbndfst * nktotf + 2
        CALL diropn(iufiltau_all, 'tau_restart', ltau_all, exst)
        CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
        !
        ! First element is the iteration number
        iqq = INT(aux(1))
        iqq = iqq + 1 ! we need to start at the next q
        nqtotf_read = INT(aux(2))
        IF (nqtotf_read /= totq) CALL errore('tau_read',&
          &'Error: The current total number of q-point is not the same as the read one. ', 1)
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              inv_tau_all(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        !
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i + 1
              zi_allvb(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufiltau_all)
      ENDIF
      !
      IF (second) THEN
        ! First inquire if the file exists
#if defined(__MPI)
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart_CB1'
#else
        name1 = TRIM(tmp_dir) // TRIM(prefix) // '.tau_restart_CB'
#endif
        INQUIRE(FILE = name1, EXIST = exst)
        !
        IF (exst) THEN ! read the file
          !
          ltau_all = nstemp * nbndfst * nktotf + 2
          CALL diropn(iufiltau_all, 'tau_restart_CB', ltau_all, exst)
          CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
          !
          ! First element is the iteration number
          iqq = INT(aux(1))
          iqq = iqq + 1 ! we need to start at the next q
          nqtotf_read = INT(aux(2))
          IF (nqtotf_read /= totq) CALL errore('tau_read',&
            &'Error: The current total number of q-point is not the same as the read one. ', 1)
          !
          i = 2
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                i = i + 1
                inv_tau_allcb(itemp, ibnd, ik) = aux(i)
              ENDDO
            ENDDO
          ENDDO
          !
          DO itemp = 1, nstemp
            DO ik = 1, nktotf
              DO ibnd = 1, nbndfst
                i = i + 1
                zi_allcb(itemp, ibnd, ik) = aux(i)
              ENDDO
            ENDDO
          ENDDO
          CLOSE(iufiltau_all)
          WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau_CB: ', iqq, '/', totq
        ENDIF
      ENDIF ! second
    ENDIF
    !
    CALL mp_bcast(exst, meta_ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(iqq, meta_ionode_id, world_comm)
      CALL mp_bcast(inv_tau_all, meta_ionode_id, world_comm)
      CALL mp_bcast(zi_allvb, meta_ionode_id, world_comm)
      IF (second) CALL mp_bcast(inv_tau_allcb, meta_ionode_id, world_comm)
      IF (second) CALL mp_bcast(zi_allcb, meta_ionode_id, world_comm)
      !
      ! Make everythin 0 except the range of k-points we are working on
      IF (lower_bnd > 1)      inv_tau_all(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) inv_tau_all(:, :, upper_bnd + 1:nktotf) = zero
      IF (lower_bnd > 1)      zi_allvb(:, :, 1:lower_bnd - 1) = zero
      IF (upper_bnd < nktotf) zi_allvb(:, :, upper_bnd + 1:nktotf) = zero
      !
      IF (second) THEN
        ! Make everythin 0 except the range of k-points we are working on
        IF (lower_bnd > 1)      inv_tau_allcb(:, :, 1:lower_bnd - 1) = zero
        IF (upper_bnd < nktotf) inv_tau_allcb(:, :, upper_bnd + 1:nktotf) = zero
        IF (lower_bnd > 1)      zi_allcb(:, :, 1:lower_bnd - 1) = zero
        IF (upper_bnd < nktotf) zi_allcb(:, :, upper_bnd + 1:nktotf) = zero
      ENDIF
      !
      WRITE(stdout, '(a,i10,a,i10)' ) '     Restart from tau: ', iqq, '/', totq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE tau_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE merge_read(nktotf, nqtotf_new, inv_tau_all_new)
    !----------------------------------------------------------------------------
    !!
    !! File merging
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE global_var,ONLY : nbndfst
    USE io_var,    ONLY : iufiltau_all
    USE io_files,  ONLY : tmp_dir, diropn
    USE input,     ONLY : nstemp, restart_filq
    USE mp,        ONLY : mp_barrier, mp_bcast
    USE mp_world,  ONLY : mpime, world_comm
    USE io_global, ONLY : ionode_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nktotf
    !! Total number of k-points
    INTEGER, INTENT(out) :: nqtotf_new
    !! Total number of q-points
    REAL(KIND = DP), INTENT(inout) :: inv_tau_all_new(nstemp, nbndfst, nktotf)
    !! Scattering rate read from file restart_filq
    !
    ! Local variables
    CHARACTER(LEN = 256) :: name1
    !! Name of the file
    LOGICAL :: exst
    !! Does the variable exist
    INTEGER :: i, iq, ios
    !! Iterative index
    INTEGER :: itemp
    !! Iterative temperature
    INTEGER :: ik
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ltau_all
    !! Length of the vector
    INTEGER :: unf_recl
    !! Record length unit
    REAL(KIND = DP) :: aux(nstemp * nbndfst * nktotf + 2)
    !! Vector to store the array
    REAL(KIND = DP) :: dummy
    !! Test what the record length is
    !
    IF (mpime == ionode_id) THEN
      !
      ! First inquire if the file exists
      name1 = TRIM(tmp_dir) // TRIM(restart_filq)
      INQUIRE(FILE = name1, EXIST = exst)
      !
      IF (exst) THEN ! read the file
        !
        ltau_all = nstemp * nbndfst * nktotf + 2
        !CALL diropn (iufiltau_all, 'tau_restart', ltau_all, exst)
        !
        INQUIRE(IOLENGTH = unf_recl) dummy
        unf_recl = unf_recl * INT(ltau_all, KIND = KIND(unf_recl))
        OPEN(UNIT = iufiltau_all, FILE = restart_filq, IOSTAT = ios, FORM ='unformatted', &
             STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
        !
        CALL davcio(aux, ltau_all, iufiltau_all, 1, -1)
        !
        ! First element is the iteration number
        iq = INT(aux(1))
        iq = iq + 1 ! we need to start at the next q
        nqtotf_new = INT(aux(2))
        !
        i = 2
        DO itemp = 1, nstemp
          DO ik = 1, nktotf
            DO ibnd = 1, nbndfst
              i = i +1
              inv_tau_all_new(itemp, ibnd, ik) = aux(i)
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iufiltau_all)
      ENDIF
    ENDIF
    !
    CALL mp_bcast(exst, ionode_id, world_comm)
    !
    IF (exst) THEN
      CALL mp_bcast(nqtotf_new, ionode_id, world_comm)
      CALL mp_bcast(inv_tau_all_new, ionode_id, world_comm)
      !
      WRITE(stdout, '(a,a)' ) '     Correctly read file ', restart_filq
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE merge_read
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE grain_boundary_scatt(inv_tau_gb)
    !----------------------------------------------------------------------------
    !!
    !! Calculate grain-boundary scattering rate \tau_nk = |v_nk|/L, with L: grain-size
    !! V.-A. Ha added 2024
    USE kinds,         ONLY : DP
    USE input,         ONLY : gb_only, gb_size
    USE global_var,    ONLY : vmef, nkf, vkk_all, nbndfst, nktotf, ibndmin, lower_bnd
    USE ep_constants,  ONLY : zero, bohr2nm
    USE control_flags, ONLY : iverbosity
    USE mp,            ONLY : mp_sum
    USE mp_global,     ONLY : world_comm, inter_pool_comm
    USE io_global,     ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(inout) :: inv_tau_gb(nbndfst, nktotf)
    !! grain-boundary scattering
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: ierr
    !! Error
    !
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(5x,a,f10.2,a)') 'Including grain boundary scattering, grain size ', gb_size*bohr2nm, ' nm'
    WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    ! 
    IF (gb_only) THEN
      WRITE(stdout, '(/5x,a)') REPEAT('=',67)
      WRITE(stdout, '(5x,"Detected ii_gb=.true., only grain boundary scattering included")')
      WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    ENDIF
    !
    ! Calculate velocity if not yet
    IF (iverbosity /= 3) THEN
      ALLOCATE(vkk_all(3, nbndfst, nktotf), STAT = ierr)
      vkk_all(:,:,:) = zero
      IF (ierr /= 0) CALL errore('print_ibte', 'Error allocating vkk_all',1)
      ! Computes the k-velocity in Ry*Bohr/hbar (hbar=1 and Bohr = 1 in Atomic Rydberg Units)
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        DO ibnd = 1, nbndfst
          vkk_all(:, ibnd, ik + lower_bnd - 1) = REAL(vmef(:, ibndmin - 1 + ibnd, ibndmin - 1 + ibnd, ikk))
        ENDDO
      ENDDO
      CALL mp_sum(vkk_all, inter_pool_comm)
    ENDIF
    ! grain boundary scattering rate in Ry/hbar (hbar=1 in Atomic Rydberg Units)
    inv_tau_gb(:, :) = zero
    DO ik = 1, nkf
      DO ibnd = 1, nbndfst
        inv_tau_gb(ibnd, ik + lower_bnd - 1) = SQRT(SUM(vkk_all(:, ibnd, ik + lower_bnd - 1)**2)) / gb_size
      ENDDO
    ENDDO
    CALL mp_sum(inv_tau_gb, inter_pool_comm)
    !
    RETURN
    !----------------------------------------------------------------------------
    END SUBROUTINE grain_boundary_scatt
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE transport_read(iq_restart, totq, lrepmatw2_restart, lrepmatw5_restart, ind_tot, &
                              ind_totcb, first_cycle)
    !----------------------------------------------------------------------------
    !
    USE control_flags,    ONLY : iverbosity    
    USE mp_global,        ONLY : inter_pool_comm, npool, my_pool_id, my_image_id, nimage
    USE mp,               ONLY : mp_bcast, mp_sum, mp_barrier
    USE io_global,        ONLY : ionode_id, stdout, ionode
    USE io_var,           ONLY : iuntaucb, iunrestart, iuntau
    USE global_var,       ONLY : nktotf, inv_tau_all, inv_tau_all_mode, inv_tau_allcb_mode, &
                                 inv_tau_all_freq, inv_tau_allcb_freq, inv_tau_allcb,       &
                                 lower_bnd, upper_bnd
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND, MPI_OFFSET
#endif
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: iq_restart
    !! Current position inside the file during writing
    INTEGER, INTENT(in) ::  totq                  
    !! total number of q-points within the fsthick window.
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! Current position inside the file during writing
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! Current position inside the file during writing (electron)
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(out) :: ind_tot
    !! Total number of points store on file
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(out) :: ind_totcb
    !! Total number of points store on file (CB)
#else
    INTEGER(KIND = 8), INTENT(inout) :: ind_tot
    !! Total number of element written to file
    INTEGER(KIND = 8), INTENT(inout) :: ind_totcb
    !! Total number of element written to file
#endif
    CHARACTER(LEN = 256) :: my_image_id_ch
    !! image id
    LOGICAL :: exst
    !! If the file exist
    INTEGER :: ipool
    !! Pool index
    INTEGER :: npool_tmp
    !! Temporary number of pools
    INTEGER :: nimage_tmp
    !! Temporary number of images
    LOGICAL :: first_cycle
    !! Check wheter this is the first cycle after a restart.
    INTEGER :: ierr
    !! Error status
    INTEGER :: ios
    !! Status of reading file
    !
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
      CALL mp_bcast(nimage_tmp, ionode_id, inter_pool_comm)
      CALL mp_bcast(lrepmatw2_restart, ionode_id, inter_pool_comm)
      CALL mp_bcast(lrepmatw5_restart, ionode_id, inter_pool_comm)
      IF (npool /= npool_tmp) CALL errore('transport_read','Number of pools is different',1)
      IF (nimage /= nimage_tmp) CALL errore('transport_read','Number of nimages is different',1)
      !
      IF (ionode) THEN
        ! scattering rate
        OPEN(UNIT = iuntau, FORM = 'unformatted', &
                            FILE = 'inv_tau_tmp' // '_' // TRIM(my_image_id_ch), STATUS = 'old')
        READ(iuntau) inv_tau_all
        CLOSE(iuntau)
        !
        OPEN(UNIT = iuntaucb, FORM = 'unformatted', &
                              FILE = 'inv_taucb_tmp' // '_' // TRIM(my_image_id_ch), STATUS = 'old')
        READ(iuntaucb) inv_tau_allcb
        CLOSE(iuntaucb)
        !
        IF (iverbosity == 3) THEN
          ! mode
          OPEN(UNIT = iuntau, FORM = 'unformatted', &
                              FILE = 'inv_tau_mode_tmp' // '_' // TRIM(my_image_id_ch), STATUS = 'old')
          READ(iuntau) inv_tau_all_mode
          CLOSE(iuntau)
          !
          OPEN(UNIT = iuntaucb, FORM = 'unformatted', &
                                FILE = 'inv_taucb_mode_tmp' // '_' // TRIM(my_image_id_ch), STATUS = 'old')
          READ(iuntaucb) inv_tau_allcb_mode
          CLOSE(iuntaucb)
          ! freq
          OPEN(UNIT = iuntau, FORM = 'unformatted', &
                              FILE = 'inv_tau_freq_tmp' // '_' // TRIM(my_image_id_ch), STATUS = 'old')
          READ(iuntau) inv_tau_all_freq
          CLOSE(iuntau)
          !
          OPEN(UNIT = iuntaucb, FORM = 'unformatted', &
                                FILE = 'inv_taucb_freq_tmp' // '_' // TRIM(my_image_id_ch), STATUS = 'old')
          READ(iuntaucb) inv_tau_allcb_freq
          CLOSE(iuntaucb)
        ENDIF ! iverbosity
      ENDIF
      !
      ! scattering rate
      CALL mp_bcast(inv_tau_all, ionode_id, inter_pool_comm)
      CALL mp_bcast(inv_tau_allcb, ionode_id, inter_pool_comm)
      IF (lower_bnd - 1 >= 1) THEN
        inv_tau_all(:, 1:lower_bnd - 1, :) = 0d0
        inv_tau_allcb(:, 1:lower_bnd - 1, :) = 0d0
      ENDIF
      IF (upper_bnd + 1 <= nktotf) THEN
        inv_tau_all(:, upper_bnd + 1:nktotf, :) = 0d0
        inv_tau_allcb(:, upper_bnd + 1:nktotf, :) = 0d0
      ENDIF
      !
      IF (iverbosity == 3) THEN
        ! scattering rate (mode)
        CALL mp_bcast(inv_tau_all_mode, ionode_id, inter_pool_comm)
        CALL mp_bcast(inv_tau_allcb_mode, ionode_id, inter_pool_comm)
        IF (lower_bnd - 1 >= 1) THEN
          inv_tau_all_mode(:, :, 1:lower_bnd - 1, :) = 0d0
          inv_tau_allcb_mode(:, :, 1:lower_bnd - 1, :) = 0d0
        ENDIF
        IF (upper_bnd + 1 <= nktotf) THEN
          inv_tau_all_mode(:, :, upper_bnd + 1:nktotf, :) = 0d0
          inv_tau_allcb_mode(:, :, upper_bnd + 1:nktotf, :) = 0d0
        ENDIF
        ! scattering rate (freq)
        CALL mp_bcast(inv_tau_all_freq, ionode_id, inter_pool_comm)
        CALL mp_bcast(inv_tau_allcb_freq, ionode_id, inter_pool_comm)
        IF (lower_bnd - 1 >= 1) THEN
          inv_tau_all_freq(:, :, 1:lower_bnd - 1) = 0d0
          inv_tau_allcb_freq(:, :, 1:lower_bnd - 1) = 0d0
        ENDIF
        IF (upper_bnd + 1 <= nktotf) THEN
          inv_tau_all_freq(:, :, upper_bnd + 1:nktotf) = 0d0
          inv_tau_allcb_freq(:, :, upper_bnd + 1:nktotf) = 0d0
        ENDIF
      ENDIF ! iverbosity
      !
#if defined(__MPI)
      CALL MPI_BCAST(ind_tot,   1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
      CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
#endif
      IF (ierr /= 0) CALL errore('use_wannier', 'error in MPI_BCAST', 1)
      !
      IF(iq_restart > 1) first_cycle = .TRUE.
      !
      ! Now, the iq_restart point has been done, so we need to do the next
      iq_restart = iq_restart + 1
      !
      IF (iq_restart <= totq) THEN
        WRITE(stdout, '(5x,a,i8,a)')'We restart from ', iq_restart, ' q-points'
      ELSE
        WRITE(stdout, '(5x,a)') 'All q-points are done, no need to restart!!'
        WRITE(stdout, '(5x,a)') 'Consider use epmatkqread = .TRUE. to restart from iteration !!'
      ENDIF
      !
    ENDIF ! exst
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE transport_read
    !----------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  END MODULE io_transport
  !------------------------------------------------------------------------------
