  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE transport
  !----------------------------------------------------------------------
  !!
  !! This module contains all the routines linked with self-consistent electronic transport
  !!   - transport_prepare
  !!   - transport_setup
  !!   - transport_main
  !!   - ibte
  !!   - transport_ibte
  !!
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(kind = DP) :: nelec_to_ncarrier          
  !! factor for the conversion from number of electrons per cell to input unit
  !
  CONTAINS
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE transport_prepare()
    !-----------------------------------------------------------------------
    !! 
    !!  This subroutine is used to prepare the parameters and flags used in the 
    !!  transport module. 
    !! 
    !-----------------------------------------------------------------------
    USE cell_base,       ONLY : omega, alat, at
    USE input,           ONLY : system_2d
    USE ep_constants,    ONLY : ryd2ev, bohr2ang, ang2cm
    USE global_var,      ONLY : ctype, inv_cell
    USE input,           ONLY : assume_metal, int_mob, carrier, ncarrier
    !
    IMPLICIT NONE
    !
    REAL(kind = DP) :: neperuc_to_cmm2 
    !! factor for the conversion from number of electrons per cell to cm^-2
    REAL(kind = DP) :: neperuc_to_cmm3 
    !! factor for the conversion from number of electrons per cell to cm^-3
    !
    IF (system_2d == 'no') THEN
      inv_cell = 1.0d0 / omega
    ELSE
      ! for 2d system need to divide by area (vacuum in z-direction)
      inv_cell = (1.0d0 / omega) * at(3, 3) * alat
    ENDIF
    !
    neperuc_to_cmm2 = inv_cell * (bohr2ang * ang2cm)**(-2.d0)
    neperuc_to_cmm3 = inv_cell * (bohr2ang * ang2cm)**(-3.d0)
    !
    IF (system_2d == 'no') THEN
      nelec_to_ncarrier = neperuc_to_cmm3
    ELSE
      ! for 2d system need to divide by area (vacuum in z-direction)
      nelec_to_ncarrier = neperuc_to_cmm2
    ENDIF
    !
    IF (assume_metal) THEN
      ctype = -1  ! act like it's for holes
    ELSEIF (int_mob) THEN
      IF (carrier) THEN
        ctype = 0 
      ELSE
        ctype = -1
      ENDIF
    ELSEIF (carrier) THEN
      IF (ncarrier > 1.E5) THEN
        ctype = 1
      ELSEIF (ncarrier < -1.E5) THEN
        ctype = -1
      ENDIF
    ELSE
      ctype = -1
    ENDIF
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE transport_prepare
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE transport_setup(lrepmatw2_restart, lrepmatw5_restart)
    !-----------------------------------------------------------------------
    !!
    !! Allocate and set initial arrays used for the transport module
    !!
    USE input,           ONLY : nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, mob_nfreq, nstemp
    USE global_var,      ONLY : inv_tau_all, inv_tau_allcb, inv_tau_all_mode,             &
                                inv_tau_allcb_mode, inv_tau_all_freq, inv_tau_allcb_freq, &
                                inv_tau_all_MPI, inv_tau_allcb_MPI, inv_tau_all_mode_MPI, &
                                inv_tau_allcb_mode_MPI, inv_tau_all_freq_MPI,             &
                                inv_tau_allcb_freq_MPI, threshold, nbndfst, nktotf
    USE ep_constants,    ONLY : zero
    USE modes,           ONLY : nmodes
    USE control_flags,   ONLY : iverbosity
    USE mp_global,       ONLY : npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(inout) :: lrepmatw2_restart(npool)
    !! To restart opening files
    INTEGER, INTENT(inout) :: lrepmatw5_restart(npool)
    !! To restart opening files
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    !
    ALLOCATE(inv_tau_all(nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_all', 1)
    ALLOCATE(inv_tau_allcb(nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_allcb', 1)
    inv_tau_all(:, :, :)   = zero
    inv_tau_allcb(:, :, :) = zero
    ! Dynamically allocate MPI auxilaries
    ALLOCATE(inv_tau_all_MPI(nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_all_MPI', 1)
    ALLOCATE(inv_tau_allcb_MPI(nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_allcb_MPI', 1)
    inv_tau_all_MPI(:, :, :)   = zero
    inv_tau_allcb_MPI(:, :, :) = zero
    !
    lrepmatw2_restart(:)   = 0
    lrepmatw5_restart(:)   = 0
    IF (iverbosity == 3) THEN
      ALLOCATE(inv_tau_all_mode(nmodes, nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_all_mode', 1)
      ALLOCATE(inv_tau_allcb_mode(nmodes, nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_allcb_mode', 1)
      ALLOCATE(inv_tau_all_freq(mob_nfreq, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_all_freq', 1)
      ALLOCATE(inv_tau_allcb_freq(mob_nfreq, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_allcb_freq', 1)
      inv_tau_all_mode(:, :, :, :)   = zero
      inv_tau_allcb_mode(:, :, :, :) = zero
      inv_tau_all_freq(:, :, :)      = zero
      inv_tau_allcb_freq(:, :, :)    = zero
      ! Dynamically allocate MPI auxilaries
      ALLOCATE(inv_tau_all_mode_MPI(nmodes, nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_all_mode_MPI', 1)
      ALLOCATE(inv_tau_allcb_mode_MPI(nmodes, nbndfst, nktotf, nstemp), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_allcb_mode_MPI', 1)
      ALLOCATE(inv_tau_all_freq_MPI(mob_nfreq, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_all_freq_MPI', 1)
      ALLOCATE(inv_tau_allcb_freq_MPI(mob_nfreq, nbndfst, nktotf), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_setup', 'Error allocating inv_tau_allcb_freq_MPI', 1)
      inv_tau_all_mode_MPI(:, :, :, :)   = zero
      inv_tau_allcb_mode_MPI(:, :, :, :) = zero
      inv_tau_all_freq_MPI(:, :, :)      = zero
      inv_tau_allcb_freq_MPI(:, :, :)    = zero
    ENDIF
    ! We save matrix elements that are smaller than machine precision (1d-16).
    ! The sum of all the elements must be smaller than that
    ! nkf1 * nkf2 * nkf3 * nqf1 * nqf2 * nqf3 * (nbndfst) * (nbndfst)
    ! must be smaller than 1d-16
    ! To avoid overflow we need to use DP
    threshold = 1d-16 / (INT(nkf1, KIND = 8) * INT(nkf2, KIND = 8) * INT(nkf3, KIND = 8) * &
                         INT(nqf1, KIND = 8) * INT(nqf2, KIND = 8) * INT(nqf3, KIND = 8) * &
                         INT((nbndfst), KIND = 8) * INT((nbndfst), KIND = 8))
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE transport_setup
    !----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE transport_main(etf_all, ind_tot, ind_totcb, ef0, efcb)
    !-----------------------------------------------------------------------
    !!
    !! Transport calculation
    !!
    USE input,            ONLY : nstemp
    USE global_var,       ONLY : nbndfst, nktotf, inv_tau_all, inv_tau_allcb, inv_tau_all_mode, &
                                 inv_tau_all_freq, inv_tau_allcb_freq, inv_tau_allcb_mode,      &
                                 inv_tau_all_MPI, inv_tau_allcb_MPI, inv_tau_all_mode_MPI,      &
                                 inv_tau_allcb_mode_MPI, inv_tau_all_freq_MPI,                  &
                                 inv_tau_allcb_freq_MPI
                                 
    USE ep_constants,     ONLY : zero
    USE control_flags,    ONLY : iverbosity
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND
#endif
    !
    IMPLICIT NONE
    !
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of component for the conduction band
#else
    INTEGER(KIND = 8), INTENT(inout) :: ind_tot
    !! Tota number of component for valence band
    INTEGER(KIND = 8), INTENT(inout) :: ind_totcb
    !! Total number of component for conduction band
#endif
    REAL(KIND = DP), INTENT(inout) :: etf_all(nbndfst, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel cas
    REAL(KIND = DP), INTENT(inout) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(inout) :: efcb(nstemp)
    !! Fermi level for the temperature itemp for cb band
    !
    ! Local variable
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP), ALLOCATABLE :: vkk_all(:, :, :)
    !! velocity from all the k-point
    REAL(KIND = DP), ALLOCATABLE :: wkf_all(:)
    !! k-point weights for all the k-points
    !
    ALLOCATE(vkk_all(3, nbndfst, nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error allocating vkk_all', 1)
    ALLOCATE(wkf_all(nktotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error allocating wkf_all', 1)
    vkk_all(:, :, :) = zero
    wkf_all(:) = zero
    !
    CALL transport_ibte(etf_all, wkf_all, vkk_all, ind_tot, ind_totcb, ef0, efcb)
    !
    DEALLOCATE(vkk_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating vkk_all', 1)
    DEALLOCATE(wkf_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating wkf_all', 1)
    DEALLOCATE(inv_tau_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_all', 1)
    DEALLOCATE(inv_tau_allcb, STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_allcb', 1)
    ! Deallocate MPI auxilaries
    DEALLOCATE(inv_tau_all_MPI, STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_all_MPI', 1)
    DEALLOCATE(inv_tau_allcb_MPI, STAT = ierr)
    IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_allcb_MPI', 1)
    !
    IF (iverbosity == 3) THEN
      DEALLOCATE(inv_tau_all_mode, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_all_mode', 1)
      DEALLOCATE(inv_tau_allcb_mode, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_allcb_mode', 1)
      DEALLOCATE(inv_tau_all_freq, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_all_freq', 1)
      DEALLOCATE(inv_tau_allcb_freq, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_allcb_freq', 1)
      ! Deallocate the MPI auxilaries
      DEALLOCATE(inv_tau_all_mode_MPI, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_all_mode_MPI', 1)
      DEALLOCATE(inv_tau_allcb_mode_MPI, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_allcb_mode_MPI', 1)
      DEALLOCATE(inv_tau_all_freq_MPI, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_all_freq_MPI', 1)
      DEALLOCATE(inv_tau_allcb_freq_MPI, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_main', 'Error deallocating inv_tau_allcb_freq_MPI', 1)
    ENDIF

    !----------------------------------------------------------------------------
    END SUBROUTINE transport_main
    !----------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE ibte(nind, etf_all, vkk_all, wkf_all, trans_prob, ef0, &
                    sparse_q, sparse_k, sparse_i, sparse_j, sparse_t, &
                    inv_tau)
    !-----------------------------------------------------------------------
    !!
    !! This routine computes the transport properties with the iterative BTE.
    !! The fine k-point and q-point grid have to be commensurate.
    !! The k-point grid uses crystal symmetry to decrease computational cost.
    !!
    USE kinds,            ONLY : sgl
    USE io_global,        ONLY : stdout
    USE cell_base,        ONLY : at, bg
    USE input,            ONLY : mob_maxiter, nstemp, broyden_beta,            &
                                 mp_mesh_k, nkf1, nkf2, nkf3, bfieldx, bfieldy,&
                                 bfieldz, lfast_kmesh, ltrans_crta, sr_crta,   &
                                 assume_metal
    USE global_var,       ONLY : nkqf, xkf, nkqtotf, nbndfst, &
                                 nktotf, xqf, gtemp, ixkqf_tr,&
                                 s_bztoibz_full, f_in_b, f_serta_b, df_in_b,   &
                                 f_out_b, kpt_bz2bztau, &
                                 kpt_bztau2bz, vkk_all_b, etf_all_b,     &
                                 wkf_all_b, inv_tau_b, nkpt_bztau_max,     &
                                 ind_map, nsym_sp
    USE ep_constants,     ONLY : zero, one, two, pi, kelvin2eV, ryd2ev, eps10, &
                                 bohr2ang, ang2cm, hbarJ, eps6, eps8, byte2Mb, &
                                 eps2, eps4, eps20, eps80, eps160, hbar, cm2m, &
                                 electronvolt_si, ry2thz_sr
    USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,        ONLY : inter_pool_comm, inter_image_comm
    USE symmetry,         ONLY : s => s_k, nsym => nsym_k
    USE printing,         ONLY : print_mob, print_mob_sym, print_mob_b, print_hall, &
                                 print_meff, print_meff_sym 
    USE io_transport,     ONLY : fin_write, fin_read
    USE io_files,         ONLY : diropn
    USE control_flags,    ONLY : iverbosity
    USE bzgrid,           ONLY : k_avg, special_points, kpoint_grid_epw, symm_mapping, &
                                 loadkmesh_fullBZ, kpmq_map
    USE parallelism,      ONLY : poolgather2
    USE transport_mag,    ONLY : select_k, unfold_k, unfold_all, create_all,   &
                                 create_indkq, k_derivative_cart, neighbk,     &
                                 k_derivative_crys, size_indkq
    !
    IMPLICIT NONE
    !
    INTEGER(KIND = 8), INTENT(in) :: nind
    !! Total number of elements per cpu
    INTEGER, INTENT(in) :: sparse_q(nind)
    !! Q-point mapping index
    INTEGER, INTENT(in) :: sparse_k(nind)
    !! K-point mapping index
    INTEGER, INTENT(in) :: sparse_i(nind)
    !! Band mapping index
    INTEGER, INTENT(in) :: sparse_j(nind)
    !! Band mapping index
    INTEGER, INTENT(in) :: sparse_t(nind)
    !! Temperature mapping index
    REAL(KIND = DP), INTENT(in) :: etf_all(nbndfst, nktotf)
    !! Eigenenergies
    REAL(KIND = DP), INTENT(inout) :: vkk_all(3, nbndfst, nktotf)
    !! Velocity of k
    REAL(KIND = DP), INTENT(in) :: wkf_all(nktotf)
    !! Weight of k
    REAL(KIND = DP), INTENT(in) :: trans_prob(nind)
    !! Transition probability
    REAL(KIND = DP), INTENT(in) :: ef0(nstemp)
    !! The Fermi level
    REAL(KIND = DP), INTENT(in) :: inv_tau(nbndfst, nktotf, nstemp)
    !! inv_tau (either inv_tau_all or inv_tau_allcb)
    !
    ! Local variables
    LOGICAL :: cart_der
    !! Can we perform Cartesian derivative
    LOGICAL :: lsp(nind)
    !! At exit, is .true. if the k-point corresponding to the index "ind" is a special point
    INTEGER(KIND = 8) :: ind
    !! Index for sparse matrix
    INTEGER :: iter
    !! Innter IBTE loop
    INTEGER :: iq
    !! q-point
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! K-point index
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: nkq_abs
    !! Index of the k+q point from the full grid.
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: bztoibz_mat(nsym, nktotf)
    !! For a given k-point in the IBZ gives the k-point index of all the
    !! k-point in the full BZ that are connected to the current one by symmetry.
    !! nsym (including TR) is the max number of symmetry
    INTEGER :: nb_sp
    !! Number of special points
    INTEGER :: ierr
    !! Error status
    INTEGER :: nrws
    !! Maximum number of WS vectors
    INTEGER :: ind1
    !! Index on the IBZ of k
    INTEGER :: ind2
    !! Index on the IBZ of k+q
    INTEGER :: nb
    !! Symmetry index
    INTEGER :: nkpt_max
    !! Maximum nb of kpoint which have non-zero inv_tau on full BZ corresponding to the ind managed by the current cpu.
    INTEGER :: nsymk(nind)
    !! For a given k-point corresponding to the index "ind", give the number of equivalent kpt by symmetry
    INTEGER, PARAMETER :: nrwsx = 200
    !! Variable for WS folding
    INTEGER, ALLOCATABLE :: xkf_sp(:, :)
    !! Special k-points
    INTEGER, ALLOCATABLE :: map_neigh(:, :, :)
    !! Map of the 6 k-point neighbor for k-derivative by finite difference.
    REAL(KIND = DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND = DP) :: carrier_density(nstemp)
    !! Carrier density [nb of carrier per unit cell]
    REAL(KIND = DP) :: f_rot(3)
    !! Rotated Fi_in by the symmetry operation
    REAL(KIND = DP) :: error(nstemp)
    !! Error in the Hall mobility
    REAL(KIND = DP) :: av_mob_old(nstemp)
    !! Average hole mobility from previous iteration
    REAL(KIND = DP) :: max_mob(nstemp)
    !! Maximum mobility use for error calculations
    REAL(KIND = DP) :: dfnk
    !! Local variable
    REAL(KIND = DP) :: etemp
    !! Local variable
    REAL(KIND = DP) :: xkk(3)
    !! K-point index for printing
    REAL(KIND = DP) :: rws(4, nrwsx)
    !! Real WS vectors
    REAL(KIND = DP) :: b_abs
    !! Absolute magnetic field (Tesla)
    REAL(KIND = DP) :: vec_min(3)
    !! Norm of the smallest vector in the x, y and z direction
    REAL(KIND = DP) :: vb(3)
    !! v_nk(i) cross B(j)
    REAL(KIND = DP) :: sigma_serta(3, 3, nstemp)
    !! SERTA conductivity tensor without Magnetic field.
    REAL(KIND = DP) :: sigma_bte(3, 3, nstemp)
    !! BTE conductivity tensor without Magnetic field.
    REAL(KIND = DP) :: sigmab_serta(3, 3, nstemp)
    !! SERTA conductivity tensor with Magnetic field.
    REAL(KIND = DP) :: sigmab_bte(3, 3, nstemp)
    !! BTE conductivity tensor with Magnetic field.
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    REAL(KIND = DP), ALLOCATABLE :: xkf_all(:, :)
    !! Collect k-point coordinate (and k+q) from all pools in parallel case
    REAL(KIND = DP), ALLOCATABLE :: f_crta(:, :, :, :)
    !! SERTA solution
    REAL(KIND = DP), ALLOCATABLE :: f_serta(:, :, :, :)
    !! SERTA solution
    REAL(KIND = DP), ALLOCATABLE :: f_in(:, :, :, :)
    !! In solution for iteration i
    REAL(KIND = DP), ALLOCATABLE :: f_out(:, :, :, :)
    !! In solution for iteration i
    REAL(KIND = DP) :: inv_tau_crta
    !! scattering rate for CRTA
    !
    inv_tau_crta = sr_crta
    !
    !! S.T: Allocations are needed to avoid heap-array problem
    !
    ALLOCATE(xkf_all(3, nkqtotf), STAT = ierr)
    IF (ierr /= 0) CALL errore('ibte', 'Error deallocating xkf_all', 1)  
    ALLOCATE(f_serta(3, nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_serta', 1)  
    ALLOCATE(f_crta(3, nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_crta', 1) 
    ALLOCATE(f_in(3, nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_in', 1) 
    ALLOCATE(f_out(3, nbndfst, nktotf, nstemp), STAT = ierr)
    IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_out', 1) 
    !
    carrier_density(:) = zero
    xkf_all(:, :) = zero
    !
#if defined(__MPI)
    CALL poolgather2(3, nkqtotf, nkqf, xkf, xkf_all)
#else
    xkf_all = xkf
#endif
    !
    IF (iverbosity == 4) THEN
      ! Array size reporting
      WRITE(stdout, '(5x,a)') 'Big array size reporting [Mb]'
      WRITE(stdout, '(5x,a)') '-- ibte --'
      WRITE(stdout, '(5x,a,f12.6)') 'bztoibz : ',  nkf1 * nkf2 * nkf3 * byte2Mb / 2 !INTEGER(4)
      WRITE(stdout, '(5x,a,f12.6)') 's_bztoibz : ',  nkf1 * nkf2 * nkf3 * byte2Mb / 2
      WRITE(stdout, '(5x,a,f12.6)') 'f_crta : ',  (3 * nbndfst * nktotf * nstemp) * byte2Mb!REAL(8)
      WRITE(stdout, '(5x,a,f12.6)') 'f_serta : ',  (3 * nbndfst * nktotf * nstemp) * byte2Mb!REAL(8)
      WRITE(stdout, '(5x,a,f12.6)') 'bztoibz_mat : ',  nsym * nktotf * byte2Mb / 2 !INTEGER(4)
      rws(:, :) = zero
      CALL wsinit(rws, nrwsx, nrws, bg)
    ENDIF
    !
    ! Deal with symmetries
    IF (mp_mesh_k) THEN
      ALLOCATE(ixkqf_tr(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error allocating ixkqf_tr', 1)
      ALLOCATE(s_bztoibz_full(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error allocating s_bztoibz_full', 1)
      ixkqf_tr(:)  = 0
      s_bztoibz_full(:) = 0
      bztoibz_mat(:, :) = 0
      !
      ! ixkqf_tr and s_bztoibz_full gets defined inside
      CALL symm_mapping(nind, bztoibz_mat, xkf_all(:, 1:nkqtotf:2), sparse_q, sparse_k)
      !
      ! Determines the special k-points are k-points that are sent to themselves via a non-identity  symmetry operation.
      CALL special_points(nb_sp, xkf_all(:, 1:nkqtotf:2), xkf_sp)
      !
    ENDIF ! mp_mesh_k
    !
    f_crta(:, :, :, :) = zero
    f_serta(:, :, :, :) = zero
    !
    IF (iverbosity == 4) THEN
      WRITE(stdout, '(2X, "   temp   k-index  ibnd", 15X, "k-point", 18X, "eig[Ry]", 19X, "f_serta")')    
    ENDIF
    DO itemp = 1, nstemp
      etemp = gtemp(itemp)
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          IF (ABS(inv_tau(ibnd, ik, itemp)) > eps160) THEN
            ekk = etf_all(ibnd, ik) - ef0(itemp)
            dfnk = w0gauss(ekk / etemp, -99) / etemp
            ! (-) sign is because w0gauss is - df/de
            f_crta(:, ibnd, ik, itemp) = - dfnk * vkk_all(:, ibnd, ik) / (inv_tau_crta)
            f_serta(:, ibnd, ik, itemp) = - dfnk * vkk_all(:, ibnd, ik) / (inv_tau(ibnd, ik, itemp))
            !
            IF (iverbosity == 4) THEN
              IF (SUM(ABS(f_serta(:, ibnd, ik, itemp))) > eps160) THEN
                xkk = xkf_all(:, 2 * ik - 1)
                CALL cryst_to_cart(1, xkk, bg, 1)
                WRITE(stdout, '(3i8,4f12.6,3E14.5)') itemp, ik, ibnd, xkk, ekk, f_serta(:, ibnd, ik, itemp)
              ENDIF
            ENDIF ! iverbosity 4
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    ! We now do CRTA with and without k-point symmetries
    IF (ltrans_crta) THEN
      WRITE(stdout, '(5x,a)') ' '
      WRITE(stdout, '(5x,a)') REPEAT('=',93)
      WRITE(stdout, '(5x,"BTE in the constant relaxation time approximation (CRTA)")')
      WRITE(stdout, '(5x,"Scattering rates for all states are set as ", F15.6, " THz.")') sr_crta*ry2thz_sr
      WRITE(stdout, '(5x,a)') REPEAT('=',93)
      ! Store conductivity without B-field
      sigma_serta(:, :, :) = zero
      max_mob(:) = zero
      ! K-point symmetry.
      IF (mp_mesh_k) THEN
        ! Averages points which leaves the k-point unchanged by symmetry in F and v.
        ! e.g. k=[1,1,1] and q=[1,0,0] with the symmetry that change x and y gives k=[1,1,1] and q=[0,1,0].
        CALL k_avg(f_crta, vkk_all, nb_sp, xkf_sp)
        CALL print_mob_sym(f_crta, bztoibz_mat, vkk_all, etf_all, wkf_all, ef0, sigma_serta)
      ELSE
        CALL print_mob(f_crta, vkk_all, etf_all, wkf_all, ef0, sigma_serta)
      ENDIF
    ENDIF
    !
    ! We now do SERTA with and without k-point symmetries
    WRITE(stdout, '(5x,a)') ' '
    WRITE(stdout, '(5x,a)') REPEAT('=',93)
    WRITE(stdout, '(5x,"BTE in the self-energy relaxation time approximation (SERTA)")')
    WRITE(stdout, '(5x,a)') REPEAT('=',93)
    ! Store conductivity without B-field
    sigma_serta(:, :, :) = zero
    max_mob(:) = zero
    f_crta(:, :, :, :) = f_crta(:, :, :, :) * inv_tau_crta
    ! K-point symmetry.
    IF (mp_mesh_k) THEN
      ! Averages points which leaves the k-point unchanged by symmetry in F and v.
      ! e.g. k=[1,1,1] and q=[1,0,0] with the symmetry that change x and y gives k=[1,1,1] and q=[0,1,0].
      CALL k_avg(f_crta, vkk_all, nb_sp, xkf_sp)
      CALL k_avg(f_serta, vkk_all, nb_sp, xkf_sp)
      CALL print_meff_sym(f_crta, bztoibz_mat, vkk_all, etf_all, wkf_all, ef0, sigma_serta)
      CALL print_mob_sym(f_serta, bztoibz_mat, vkk_all, etf_all, wkf_all, ef0, sigma_serta)
    ELSE
      CALL print_meff(f_crta, vkk_all, etf_all, wkf_all, ef0, sigma_serta)
      CALL print_mob(f_serta, vkk_all, etf_all, wkf_all, ef0, sigma_serta)
    ENDIF
    !
    ! Possibily read from file
    iter = 1
    f_in(:, :, :, :) = zero
    !IF (ncarrier > 1E5) THEN
    !  CALL fin_read(iter, F_in, av_mob_old, .TRUE.)
    !ENDIF
    !!
    !IF (ncarrier < -1E5) THEN
    !  CALL fin_read(iter, F_in, av_mob_old, .FALSE.)
    !ENDIF
    ! If it is the first time, put to SERTA
    IF (iter == 1) THEN
      f_in(:, :, :, :) = f_serta(:, :, :, :)
      av_mob_old(:) = zero
    ENDIF
    !
    f_out(:, :, :, :) = zero
    sigma_bte(:, :, :) = zero
    error(:) = 1000
    ! Now compute the Iterative solution for electron or hole
    WRITE(stdout, '(5x,a)') ' '
    WRITE(stdout, '(5x,a)') REPEAT('=',93)
    WRITE(stdout, '(5x,"Start solving iterative Boltzmann Transport Equation")')
    WRITE(stdout, '(5x,a/)') REPEAT('=',93)
    !
    DO WHILE(MAXVAL(error) > eps6)
      WRITE(stdout, '(5x,"Iteration number:", i10)') iter
      !
      IF (iter > mob_maxiter) THEN
        WRITE(stdout, '(5x,a)') REPEAT('=',93)
        WRITE(stdout, '(5x,"The iteration reached the maximum but did not converge.")')
        WRITE(stdout, '(5x,a/)') REPEAT('=',93)
        EXIT
      ENDIF
      !
      ! SP: The reason for the "-" sign is because f_out and F_in are actually
      !     equal to $-\partial_E f_{nk}$ which is used for the mobility.
      !     However in the BTE, you need F_in to be $\partial_E f_{mk+q}$
      IF (mp_mesh_k) THEN
        ! Use k-point symmetry
        DO ind = 1, nind
          !
          f_rot(:) = zero
          iq    = sparse_q(ind)
          ik    = sparse_k(ind)
          ibnd  = sparse_i(ind)
          jbnd  = sparse_j(ind)
          itemp = sparse_t(ind)
          !
          CALL cryst_to_cart(1, f_in(:, jbnd, ixkqf_tr(ind), itemp), at, -1)
          CALL DGEMV('n', 3, 3, 1.d0, REAL(s(:, :, s_bztoibz_full(ind)), KIND = DP), &
                     3, f_in(:, jbnd, ixkqf_tr(ind), itemp), 1, 0.d0, f_rot(:), 1)
          CALL cryst_to_cart(1, f_in(:, jbnd, ixkqf_tr(ind), itemp), bg, 1)
          CALL cryst_to_cart(1, F_rot, bg, 1)
          !
          f_out(:, ibnd, ik, itemp) = f_out(:, ibnd, ik, itemp) + trans_prob(ind) * f_rot(:)
          !
        ENDDO
      ELSE
        DO ind = 1, nind
          !
          iq    = sparse_q(ind)
          ik    = sparse_k(ind)
          ibnd  = sparse_i(ind)
          jbnd  = sparse_j(ind)
          itemp = sparse_t(ind)
          ! We need F_in at k+q point
          CALL kpmq_map(xkf_all(:, 2 * ik - 1), xqf(:, iq), +1, nkq_abs)
          !
          f_out(:, ibnd, ik, itemp) = f_out(:, ibnd, ik, itemp) + trans_prob(ind) * f_in(:, jbnd, nkq_abs, itemp)
          !
        ENDDO
      ENDIF
      !
      CALL mp_sum(f_out, inter_pool_comm)
      CALL mp_sum(f_out, inter_image_comm)
      !
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            IF (ABS(inv_tau(ibnd, ik, itemp)) > eps160) THEN
              f_out(:, ibnd, ik, itemp) = f_serta(:, ibnd, ik, itemp) + &
                               f_out(:, ibnd, ik, itemp) / (inv_tau(ibnd, ik, itemp))
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !
      IF (mp_mesh_k) THEN
        CALL k_avg(f_out, vkk_all, nb_sp, xkf_sp)
        CALL print_mob_sym(f_out, bztoibz_mat, vkk_all, etf_all, wkf_all, ef0, sigma_bte, max_mob)
      ELSE
        CALL print_mob(f_out, vkk_all, etf_all, wkf_all, ef0, sigma_bte, max_mob)
      ENDIF
      !
      ! Computes the error
      DO itemp = 1, nstemp
        error(itemp) = ABS(max_mob(itemp) - av_mob_old(itemp))
      ENDDO
      av_mob_old = max_mob
      WRITE(stdout, '(a)')
      WRITE(stdout, '(50x, 1E16.6, a)') MAXVAL(error), '    Max error'
      !
      ! Save F_in
      ! Linear mixing
      f_in = (1.0 - broyden_beta) * f_in + broyden_beta * f_out
      f_out = zero
      !
      iter = iter + 1
      !
      ! Save F_in to file:
      !IF (ncarrier > 1E5) THEN
      !  CALL fin_write(iter, F_in, av_mob_old, .TRUE.)
      !ENDIF
      !!
      !IF (ncarrier < -1E5) THEN
      !  CALL fin_write(iter, F_in, av_mob_old, .FALSE.)
      !ENDIF
      !
    ENDDO ! end of while loop
    !
    ! Print nk-resolved mobility $\mu_{nk}^{\alpha\beta}$
    IF (iverbosity == 3  .AND. .NOT. assume_metal) THEN
      IF (mp_mesh_k) THEN
        CALL print_mob_sym(f_in, bztoibz_mat, vkk_all, etf_all, wkf_all, ef0, sigma_bte, max_mob, xkf_all(:, 1:nkqtotf:2))
      ENDIF
    ENDIF
    !
    ! ------------------
    ! - Magnetic field -
    ! ------------------
    ! Possibly add the contribution from an external B-field.
    ! We need to use the full BZ because the magnetic field breaks some crystal symmetries.
    ! Therefore we consider no crystal symmetry and reconstruct all.
    !
    b_abs = ABS(bfieldx) + ABS(bfieldy) + ABS(bfieldz)
    IF (b_abs > eps80) THEN
      !
      WRITE(stdout, '(5x,a)') ' '
      WRITE(stdout, '(5x,a)') REPEAT('=',93)
      WRITE(stdout, '(5x,"BTE in the SERTA with B-field")')
      WRITE(stdout, '(5x,a)') REPEAT('=',93)
      !
      ! Generate a grid of k-points on the full BZ without symmetries
      ! The variable xkf_bz gets allocated here
      IF (.NOT. lfast_kmesh) CALL loadkmesh_fullBZ()
      !
      ! Select the k-point that have non-zero inv_tau scattering rate
      ! In this routine the nkpt_ibztau and kpt_ibztau2ibz are allocated
      CALL select_k(inv_tau)
      !
      IF (mp_mesh_k) THEN
        ! Since the external B-field breaks some crystal symmetries, we need to unfold
        ! every quantities in the full BZ. We cannot use k-point symmetries.
        ! In this routine the kpt_bz2bztau and kpt_bztau2bz are allocated
        ! and the nkpt_bztau_max is defined.
        CALL unfold_k(bztoibz_mat)
        !
        ! Unfold the etf, vkk, xkf, inv_tau and f_serta to full BZ composed of non-zero inv_tau points
        ! This routines allocate etf_all_b, vkk_all_b, wkf_all_b, f_serta_b, f_in_b, f_out_b, inv_tau_b
        CALL unfold_all(etf_all, vkk_all, inv_tau, f_serta)
        !
        ! Precompute the maximum size of the ind1 and ind2 pairs
        ! This routines computes nkpt_max.
        CALL size_indkq(nind, nb_sp, xkf_sp, bztoibz_mat, sparse_q, sparse_k, sparse_i, sparse_j, sparse_t, &
                        inv_tau, nkpt_max)
        !
        ! For a given ind, save the index k and k+q (ind1 and ind2) that are related by crystal symmetry.
        ! In the case of special k-points, count them in nb_sp_map(nkpt_bzfstd) where nkpt_bzfstd is the number of
        ! k-points with non-zero inv_tau for the ind range consider by this specific core.
        ! This routines allocates nsym_sp and ind_map.
        lsp(:) = .FALSE.
        nsymk(:) = 0
        CALL create_indkq(nind, nkpt_max, nb_sp, xkf_sp, bztoibz_mat, sparse_q, sparse_k, sparse_i, sparse_j, &
                          sparse_t, inv_tau, nsymk, lsp)
        !
      ELSE ! Homogeneous mesh without k-point symmetry
        !
        ! This routines allocate etf_all_b, vkk_all_b, wkf_all_b, f_serta_b, f_in_b, f_out_b, inv_tau_b
        CALL create_all(etf_all, vkk_all, inv_tau, f_serta)
        !
      ENDIF ! mp_mesh_k
      !
      ! f_in_b is equal to f_serta_b for the first iteration
      f_in_b = f_serta_b
      !
      iter = 0
      error(:) = 1000
      ! The tolerence is more strict for B-field IBTE
      DO WHILE(MAXVAL(error) > eps20)
        IF (iter == 1) THEN
          WRITE(stdout, '(5x,a)') ' '
          WRITE(stdout, '(5x,a)') REPEAT('=',93)
          WRITE(stdout, '(5x,"Start solving iterative Boltzmann Transport Equation with B-field")')
          WRITE(stdout, '(5x,a/)') REPEAT('=',93)
          WRITE(stdout, '(5x,"Iteration number:", i10)') iter
        ENDIF
        IF (iter > 1) THEN
          WRITE(stdout, '(5x,a)') ' '
          WRITE(stdout, '(5x,"Iteration number:", i10)') iter
        ENDIF ! iter >0
        IF (iter > mob_maxiter) THEN
          WRITE(stdout, '(5x,a)') REPEAT('=',93)
          WRITE(stdout, '(5x,"The iteration reached the maximum but did not converge.")')
          WRITE(stdout, '(5x,a/)') REPEAT('=',93)
          EXIT
        ENDIF
        !
        f_out_b(:, :, :, :) = zero
        !
        ! iter == 0 corresponds to SERTA with B-field
        IF (iter > 0) THEN
          ikk = 0
          DO ind = 1, nind
            iq    = sparse_q(ind)
            ik    = sparse_k(ind)
            ibnd  = sparse_i(ind)
            jbnd  = sparse_j(ind)
            itemp = sparse_t(ind)
            !
            IF (mp_mesh_k) THEN
              IF (lsp(ind)) THEN
                ! The k-point is a special k-point for which some symmetries send it to himself
                DO nb = 1, nsymk(ind)
                  ikk = ikk + 1
                  ind1 = ind_map(1, ikk)
                  ind2 = ind_map(2, ikk)
                  f_out_b(:, ibnd, ind1, itemp) = f_out_b(:, ibnd, ind1, itemp)  &
                      + trans_prob(ind) * f_in_b(:, jbnd, ind2, itemp) / DBLE(nsym_sp(ikk))

                ENDDO ! nb
              ELSE ! not a special point
                DO nb = 1, nsymk(ind)
                  ikk = ikk + 1
                  ind1 = ind_map(1, ikk)
                  ind2 = ind_map(2, ikk)
                  f_out_b(:, ibnd, ind1, itemp) = f_out_b(:, ibnd, ind1, itemp)  &
                      + trans_prob(ind) * f_in_b(:, jbnd, ind2, itemp)

                ENDDO ! nb
              ENDIF ! Logical special point
            ELSE ! mp_mesh_k
              !
              CALL kpmq_map(xkf_all(:, 2 * ik - 1), xqf(:, iq), +1, nkq_abs)
              f_out_b(:, ibnd, ik, itemp) = f_out_b(:, ibnd, ik, itemp)  &
                        + trans_prob(ind) * f_in_b(:, jbnd, nkq_abs, itemp)
            ENDIF ! mp_mesh_k
          ENDDO ! ind
          CALL mp_sum(F_out_b, inter_pool_comm)
          CALL mp_sum(F_out_b, inter_image_comm)
        ENDIF ! iter > 0
        !
        ! Multiply with the scattering rate inv_tau_b
        DO itemp = 1, nstemp
          DO ik = 1, nkpt_bztau_max
            DO ibnd = 1, nbndfst
              IF (ABS(inv_tau_b(ibnd, ik, itemp)) > eps160) THEN
                f_out_b(:, ibnd, ik, itemp) = f_serta_b(:, ibnd, ik, itemp) + &
                            f_out_b(:, ibnd, ik, itemp) / (inv_tau_b(ibnd, ik, itemp))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        !
        IF (iter == 0) THEN
          ! Computing the derivative of the population takes time. Hence compute the mapping between neighbour k-points
          ! only once and store it here.
          ! We also check if we can do the derivative in Cartesian unit (more accurate) or if we have to do in crystal unit.
          ! In some tilded system with too coarse grids, increment points along the pure Cartesian direction
          ! (x,0,0), (0,y,0) and (0,0,z) cannot be found and one has to use crystal derivative
          ALLOCATE(map_neigh(6, nkpt_bztau_max, nstemp), STAT = ierr)
          IF (ierr /= 0) CALL errore('ibte', 'Error allocating map_neigh', 1)
          map_neigh(:, :, :) = 0
          cart_der = .TRUE.
          CALL neighbk(map_neigh, cart_der, vec_min)
        ENDIF
        !
        ! Computes df_in/dk using finite difference
        IF (cart_der) THEN
          CALL k_derivative_cart(map_neigh, f_in_b, vec_min, df_in_b)
        ELSE
          CALL k_derivative_crys(map_neigh, f_in_b, df_in_b)
        ENDIF
        !
        ! Now add the finite magnetic field part
        DO itemp = 1, nstemp
          etemp = gtemp(itemp)
          DO ikbz = 1, nkpt_bztau_max
            DO ibnd = 1, nbndfst
              IF (ABS(inv_tau_b(ibnd, ikbz, itemp)) > eps160) THEN
                vb(1) = vkk_all_b(2, ibnd, ikbz, itemp) * bfieldz - vkk_all_b(3, ibnd, ikbz, itemp) * bfieldy
                vb(2) = vkk_all_b(3, ibnd, ikbz, itemp) * bfieldx - vkk_all_b(1, ibnd, ikbz, itemp) * bfieldz
                vb(3) = vkk_all_b(1, ibnd, ikbz, itemp) * bfieldy - vkk_all_b(2, ibnd, ikbz, itemp) * bfieldx
                f_out_b(1, ibnd, ikbz, itemp) =  DOT_PRODUCT(vb(:), df_in_b(:, 1, ibnd, ikbz, itemp)) &
                                / (inv_tau_b(ibnd, ikbz, itemp)) + f_out_b(1, ibnd, ikbz, itemp)
                f_out_b(2, ibnd, ikbz, itemp) =  DOT_PRODUCT(vb(:), df_in_b(:, 2, ibnd, ikbz, itemp)) &
                                / (inv_tau_b(ibnd, ikbz, itemp)) + f_out_b(2, ibnd, ikbz, itemp)
                f_out_b(3, ibnd, ikbz, itemp) =  DOT_PRODUCT(vb(:), df_in_b(:, 3, ibnd, ikbz, itemp)) &
                                / (inv_tau_b(ibnd, ikbz, itemp)) + f_out_b(3, ibnd, ikbz, itemp)
              ENDIF ! inv_tau_b
            ENDDO ! ibnd
          ENDDO ! ik
        ENDDO ! itemp
        !
        IF (iter == 0) THEN
          CALL print_mob_b(f_out_b, vkk_all_b, etf_all_b, wkf_all_b, ef0, &
                           carrier_density, sigmab_serta, max_mob)
        ELSE
          CALL print_mob_b(f_out_b, vkk_all_b, etf_all_b, wkf_all_b, ef0, &
                           carrier_density, sigmab_bte, max_mob)
        ENDIF
        !
        ! Computes the error
        DO itemp = 1, nstemp
          error(itemp) = ABS(max_mob(itemp) - av_mob_old(itemp))
        ENDDO
        av_mob_old = max_mob
        WRITE(stdout, '(a)')
        WRITE(stdout, '(50x, 1E16.6, a)') MAXVAL(error), '    Max error'
        !
        ! Save F_in
        ! Linear mixing
        f_in_b = (1.0 - broyden_beta) * f_in_b + broyden_beta * f_out_b
        f_out_b = zero
        !
        iter = iter + 1
        !
      ENDDO ! end of while loop
      !
      ! Once we are converged print Hall factor
      WRITE(stdout, '(/5x, a)') REPEAT('=', 93)
      WRITE(stdout, '(5x, a)') 'Summary and Hall factor'
      WRITE(stdout, '(5x, a)') REPEAT('=', 93)
      !
      CALL print_hall(carrier_density, sigma_serta, sigma_bte, sigmab_serta, sigmab_bte)
      !
      IF (mp_mesh_k) THEN
        DEALLOCATE(kpt_bz2bztau, STAT = ierr)
        IF (ierr /= 0) CALL errore('ibte', 'Error deallocating kpt_bz2bztau', 1)
        DEALLOCATE(kpt_bztau2bz, STAT = ierr)
        IF (ierr /= 0) CALL errore('ibte', 'Error deallocating kpt_bztau2bz', 1)
      ENDIF ! mp_mesh_k
      DEALLOCATE(vkk_all_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating vkk_all_b', 1)
      DEALLOCATE(etf_all_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating etf_all_b', 1)
      DEALLOCATE(wkf_all_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating wkf_all_b', 1)
      DEALLOCATE(f_serta_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_serta_b', 1)
      DEALLOCATE(f_in_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_in_b', 1)
      DEALLOCATE(f_out_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating f_out_b', 1)
      DEALLOCATE(inv_tau_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating inv_tau_b', 1)
      DEALLOCATE(df_in_b, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating df_in_b', 1)
    ENDIF ! B_abs
    !
    ! Deallocate
    IF (mp_mesh_k) THEN
      DEALLOCATE(xkf_sp, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating xkf_sp', 1)
      DEALLOCATE(ixkqf_tr, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating ixkqf_tr', 1)
      DEALLOCATE(s_bztoibz_full, STAT = ierr)
      IF (ierr /= 0) CALL errore('ibte', 'Error deallocating s_bztoibz_full', 1)
    ENDIF
    !
    RETURN
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE ibte
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE transport_ibte(etf_all, wkf_all, vkk_all, ind_tot, ind_totcb, ef0, efcb)
    !----------------------------------------------------------------------------
    !!
    !! This routines opens all the required files to restart an IBTE calculation
    !! then call the ibte SUBROUTINE to perform the iterations.
    !! This routine requires that the scattering rates have been computed previously.
    !!
    ! ----------------------------------------------------------------------------
    USE kinds,            ONLY : i4b
    USE global_var,       ONLY : inv_tau_all, inv_tau_allcb, nbndfst, nktotf, dos, &
                                 evbm, ecbm, ctype
    USE mp_world,         ONLY : mpime, world_comm
    USE mp_images,        ONLY : my_image_id
    USE mp_global,        ONLY : my_pool_id, inter_pool_comm
    USE io_global,        ONLY : ionode_id, ionode, stdout
    USE io_files,         ONLY : tmp_dir, prefix
    USE input,            ONLY : nstemp, ncarrier, assume_metal
    USE ep_constants,     ONLY : zero
    USE io_var,           ONLY : iufilibtev_sup, iunepmat, iunsparseq, iunsparsek, &
                                 iunsparsei, iunsparsej, iunsparset, iunsparseqcb, &
                                 iunsparsekcb, iunsparseicb, iunsparsejcb,         &
                                 iunsparsetcb, iunepmatcb
    USE mp,               ONLY : mp_bcast
    USE parallelism,      ONLY : fkbounds2
    USE symmetry,         ONLY : nsym => nsym_k
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET, MPI_MODE_RDONLY, MPI_INFO_NULL, &
                                 MPI_SEEK_SET, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, &
                                 MPI_OFFSET_KIND, MPI_INTEGER4
#endif
    !
    IMPLICIT NONE
    !
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_tot
    !! Total number of component for valence band
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(inout) :: ind_totcb
    !! Total number of component for the conduction band
#else
    INTEGER(KIND = 8), INTENT(inout) :: ind_tot
    !! Tota number of component for valence band
    INTEGER(KIND = 8), INTENT(inout) :: ind_totcb
    !! Total number of component for conduction band
#endif
    !
    REAL(KIND = DP), INTENT(inout) :: etf_all(nbndfst, nktotf)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(KIND = DP), INTENT(inout) :: wkf_all(nktotf)
    !! k-point weights from all the cpu
    REAL(KIND = DP), INTENT(inout) :: vkk_all(3, nbndfst, nktotf)
    !! velocity from all the k-points
    REAL(KIND = DP), INTENT(inout) :: ef0(nstemp)
    !! Fermi level for the temperature itemp
    REAL(KIND = DP), INTENT(inout) :: efcb(nstemp)
    !! Fermi level for the temperature itemp for cb band
    !
    ! Local variables
    CHARACTER(LEN = 256) :: filint
    !! Name of the file to write/read
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! K-point
    INTEGER :: ios
    !! IO error message
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ibnd
    !! Counter on bandA
    INTEGER :: iktmp
    !! Dummy counter for k-points
    INTEGER :: ibtmp
    !! Dummy counter for bands
    INTEGER :: direct_io_factor
    !! Direct IO
    INTEGER(KIND = i4b) :: dum_int
    !! Dummy integer
    INTEGER(KIND = 8) :: nind
    !! Number of local elements per cores.
    INTEGER(KIND = 8) :: unf_recl
    !! double precision to prevent integer overflow
    INTEGER(KIND = i4b), ALLOCATABLE :: sparse_q(:)
    !! Index mapping for q-points
    INTEGER(KIND = i4b), ALLOCATABLE :: sparse_k(:)
    !! Index mapping for k-points
    INTEGER(KIND = i4b), ALLOCATABLE :: sparse_i(:)
    !! Index mapping for i bands
    INTEGER(KIND = i4b), ALLOCATABLE :: sparse_j(:)
    !! Index mapping for j bands
    INTEGER(KIND = i4b), ALLOCATABLE :: sparse_t(:)
    !! Index mapping for temperature
    INTEGER(KIND = i4b), ALLOCATABLE :: sparsecb_q(:)
    !! Index mapping for q-points for cb
    INTEGER(KIND = i4b), ALLOCATABLE :: sparsecb_k(:)
    !! Index mapping for k-points for cb
    INTEGER(KIND = i4b), ALLOCATABLE :: sparsecb_i(:)
    !! Index mapping for i bands for cb
    INTEGER(KIND = i4b), ALLOCATABLE :: sparsecb_j(:)
    !! Index mapping for j bands for cb
    INTEGER(KIND = i4b), ALLOCATABLE :: sparsecb_t(:)
    !! Index mapping for temperature for cb
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw2
    !! Local core offset for reading
    INTEGER(KIND = MPI_OFFSET_KIND) :: lrepmatw4
    !! Local core offset for reading
    INTEGER(KIND = MPI_OFFSET_KIND) :: lsize
    !! Offset to tell where to start reading the file
    INTEGER(KIND = MPI_OFFSET_KIND) :: lower_bnd
    !! start for current CPU
    INTEGER(KIND = MPI_OFFSET_KIND) :: upper_bnd
    !! end for current CPU
#else
    INTEGER(KIND = 8) :: lrepmatw2
    !! Local core offset for reading
    INTEGER(KIND = i4b) :: lrepmatw4
    !! Local core offset for reading
    INTEGER(KIND = 8) :: lsize
    !! Offset to tell where to start reading the file
    INTEGER(KIND = 8) :: lower_bnd
    !! start for current CPU
    INTEGER(KIND = 8) :: upper_bnd
    !! end for current CPU
#endif
    REAL(KIND = DP) :: dum1
    !! Dummy variable
    REAL(KIND = DP), ALLOCATABLE :: trans_prob(:)
    !! Transition probabilities
    REAL(KIND = DP), ALLOCATABLE :: trans_probcb(:)
    !! Transition probabilities for cb
    CHARACTER(LEN = 256) :: my_image_id_ch
    !! Image ID
    WRITE(my_image_id_ch, "(I0)") my_image_id
    !
    etf_all(:, :)    = zero
    wkf_all(:)       = zero
    vkk_all(:, :, :) = zero
    !
    ! SP - The implementation only works with MPI so far
    ! Read velocities
    IF (ionode) THEN
      !
      OPEN(UNIT = iufilibtev_sup, FILE = 'IBTEvel_sup' // '_' // TRIM(my_image_id_ch) // '.fmt', STATUS = 'old', IOSTAT = ios)
      READ(iufilibtev_sup, '(a)')
      READ(iufilibtev_sup, *) ind_tot, ind_totcb
      READ(iufilibtev_sup, '(a)')
      READ(iufilibtev_sup, *) evbm, ecbm
      READ(iufilibtev_sup, '(a)')
      DO itemp = 1, nstemp
        READ(iufilibtev_sup, *) dum1, ef0(itemp), efcb(itemp)
      ENDDO
      READ(iufilibtev_sup, '(a)')
      !
      DO ik = 1, nktotf
        DO ibnd = 1, nbndfst
          READ(iufilibtev_sup, *) iktmp, ibtmp, vkk_all(:, ibnd, ik), etf_all(ibnd, ik), wkf_all(ik)
        ENDDO
      ENDDO
      !
      IF (assume_metal) THEN
        DO itemp = 1, nstemp
          READ(iufilibtev_sup, *) dum1, dos(itemp)
        ENDDO
      ENDIF
      !
      CLOSE(iufilibtev_sup)
      !
      inv_tau_all(:, :, :) = zero
      OPEN(UNIT = iufilibtev_sup, FILE = 'inv_tau' // '_' // TRIM(my_image_id_ch) // '.fmt', STATUS = 'old', IOSTAT = ios)
      READ(iufilibtev_sup, '(a)')
      READ(iufilibtev_sup, '(a)')
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            READ(iufilibtev_sup, *) dum1, dum1, dum1, dum1, inv_tau_all(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iufilibtev_sup)
      !
      inv_tau_allcb(:, :, :) = zero
      OPEN(UNIT = iufilibtev_sup, FILE = 'inv_taucb' // '_' // TRIM(my_image_id_ch) // '.fmt', STATUS = 'old', IOSTAT = ios)
      READ(iufilibtev_sup, '(a)')
      READ(iufilibtev_sup, '(a)')
      DO itemp = 1, nstemp
        DO ik = 1, nktotf
          DO ibnd = 1, nbndfst
            READ(iufilibtev_sup, *) dum1, dum1, dum1, dum1, inv_tau_allcb(ibnd, ik, itemp)
          ENDDO
        ENDDO
      ENDDO
      CLOSE(iufilibtev_sup)
    ENDIF
    !
#if defined(__MPI)
    CALL MPI_BCAST(ind_tot, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
    CALL MPI_BCAST(ind_totcb, 1, MPI_OFFSET, ionode_id, inter_pool_comm, ierr)
#endif
    CALL mp_bcast(ef0, ionode_id, inter_pool_comm)
    CALL mp_bcast(evbm, ionode_id, inter_pool_comm)
    CALL mp_bcast(ecbm, ionode_id, inter_pool_comm)
    CALL mp_bcast(efcb, ionode_id, inter_pool_comm)
    CALL mp_bcast(vkk_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(wkf_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(etf_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(nsym, ionode_id, inter_pool_comm)
    CALL mp_bcast(inv_tau_all, ionode_id, inter_pool_comm)
    CALL mp_bcast(inv_tau_allcb, ionode_id, inter_pool_comm)
    !
    ! Now choose hole or electron (the implementation does not support both)
    ! hole (or metals)
    IF (ctype < 0) THEN
      !
      ! Split all the matrix elements across all cores.
      CALL fkbounds2(ind_tot, lower_bnd, upper_bnd)
      !
      ! Allocate the local size
      nind = upper_bnd - lower_bnd + 1
      WRITE(stdout, '(5x,a,i10)') 'Number of elements per core ', nind
      ALLOCATE(trans_prob(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating trans_prob', 1)
      trans_prob(:) = 0.0d0
      !
      ! Open file containing trans_prob
      filint = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkq1' // '_' // TRIM(my_image_id_ch)
#if defined(__MPI)
      CALL MPI_FILE_OPEN(inter_pool_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmat, ierr)
#else
      ! Note : For unformatted RECL, the size must be expressed as an even multiple of four
      INQUIRE(IOLENGTH = direct_io_factor) dum1
      unf_recl = direct_io_factor * INT(nind, KIND = KIND(unf_recl))
      !INQUIRE(FILE = 'si.epmatkq1', SIZE = unf_recl)
      !print*,'The read record length is ',unf_recl
      OPEN(UNIT = iunepmat, FILE = filint, IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'old', ACCESS = 'direct', RECL = unf_recl)
#endif
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN X.epmatkq1', 1)
      !
#if defined(__MPI)
      ! Offset depending on CPU
      lrepmatw2 = INT(lower_bnd - 1_MPI_OFFSET_KIND, KIND = MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND
      !
      ! Size of what we read
      lsize = INT(nind, KIND = MPI_OFFSET_KIND)
      !
      CALL MPI_FILE_SEEK(iunepmat, lrepmatw2, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK', 1)
      CALL MPI_FILE_READ(iunepmat, trans_prob(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ', 1)
      !
      ! Now open the sparse matrix mapping
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparseq' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparseq, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparseq', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparsek' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsek, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparsek', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparsei' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsei, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparsei', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparsej' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsej, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparsej', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparset' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparset, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparset', 1)
#else
      READ(UNIT = iunepmat, REC = 1, IOSTAT = ierr) trans_prob
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading X.epmatkq1', 1)
      !
      ! Now open the sparse matrix mapping
      INQUIRE(IOLENGTH = direct_io_factor) dum_int
      unf_recl = direct_io_factor * INT(nind, KIND = KIND(unf_recl))
      OPEN(UNIT = iunsparseq, FILE = 'sparseq', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparseq', 1)
      OPEN(UNIT = iunsparsek, FILE = 'sparsek', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparsek', 1)
      OPEN(UNIT = iunsparsei, FILE = 'sparsei', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparsei', 1)
      OPEN(UNIT = iunsparsej, FILE = 'sparsej', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparsej', 1)
      OPEN(UNIT = iunsparset, FILE = 'sparset', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparset', 1)
#endif
      !
      ALLOCATE(sparse_q(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparse_q', 1)
      ALLOCATE(sparse_k(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparse_k', 1)
      ALLOCATE(sparse_i(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparse_i', 1)
      ALLOCATE(sparse_j(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparse_j', 1)
      ALLOCATE(sparse_t(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparse_t', 1)
      sparse_q(:) = 0.0d0
      sparse_k(:) = 0.0d0
      sparse_i(:) = 0.0d0
      sparse_j(:) = 0.0d0
      sparse_t(:) = 0.0d0
      !
#if defined(__MPI)
      lrepmatw4 = INT(lower_bnd - 1_MPI_OFFSET_KIND, KIND = MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND
      !
      CALL MPI_FILE_SEEK(iunsparseq, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK', 1)
      CALL MPI_FILE_READ(iunsparseq, sparse_q(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ', 1)
      CALL MPI_FILE_SEEK(iunsparsek, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK', 1)
      CALL MPI_FILE_READ(iunsparsek, sparse_k(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ', 1)
      CALL MPI_FILE_SEEK(iunsparsei, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK', 1)
      CALL MPI_FILE_READ(iunsparsei, sparse_i(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ', 1)
      CALL MPI_FILE_SEEK(iunsparsej, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK', 1)
      CALL MPI_FILE_READ(iunsparsej, sparse_j(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ', 1)
      CALL MPI_FILE_SEEK(iunsparset, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK', 1)
      CALL MPI_FILE_READ(iunsparset, sparse_t(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ', 1)
#else
      READ(iunsparseq, REC = 1, IOSTAT = ierr) sparse_q
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparse_q', 1)
      READ(iunsparsek, REC = 1, IOSTAT = ierr) sparse_k
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparse_k', 1)
      READ(iunsparsei, REC = 1, IOSTAT = ierr) sparse_i
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparse_i', 1)
      READ(iunsparsej, REC = 1, IOSTAT = ierr) sparse_j
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparse_j', 1)
      READ(iunsparset, REC = 1, IOSTAT = ierr) sparse_t
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparse_t', 1)
#endif
      !
      ! Now call the ibte to solve the BTE iteratively until convergence
      CALL ibte(nind, etf_all, vkk_all, wkf_all, trans_prob, ef0, sparse_q, sparse_k, &
                sparse_i, sparse_j, sparse_t, inv_tau_all)
      !
#if defined(__MPI)
      CALL MPI_FILE_CLOSE(iunepmat, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparseq, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparsek, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparsei, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparsej, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparset, ierr)
#else
      CLOSE(iunepmat, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing X.epmatkq1', 1)
      CLOSE(iunsparseq, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparseq', 1)
      CLOSE(iunsparsek, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparsek', 1)
      CLOSE(iunsparsei, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparsei', 1)
      CLOSE(iunsparsej, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparsej', 1)
      CLOSE(iunsparset, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparset', 1)
#endif
      !
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      DEALLOCATE(trans_prob, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating trans_prob', 1)
      DEALLOCATE(sparse_q, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparse_q', 1)
      DEALLOCATE(sparse_k, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparse_k', 1)
      DEALLOCATE(sparse_i, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparse_i', 1)
      DEALLOCATE(sparse_j, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparse_j', 1)
      DEALLOCATE(sparse_t, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparse_t', 1)
      !
    ENDIF
    !
    ! Electrons
    IF (ctype >= 0) THEN
      !
      CALL fkbounds2(ind_totcb, lower_bnd, upper_bnd)
      ! Allocate the local size
      nind = upper_bnd - lower_bnd + 1
      WRITE(stdout, '(5x,a,i10)') 'Number of elements per core ', nind
      ALLOCATE(trans_probcb(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating trans_probcb', 1)
      trans_probcb(:) = 0.0d0
      !
      ! Open file containing trans_prob
      filint = TRIM(tmp_dir) // TRIM(prefix) // '.epmatkqcb1' // '_' // TRIM(my_image_id_ch)
#if defined(__MPI)
      CALL MPI_FILE_OPEN(inter_pool_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatcb, ierr)
#else
      INQUIRE(IOLENGTH = direct_io_factor) dum1
      unf_recl = direct_io_factor * INT(nind, KIND = KIND(unf_recl))
      OPEN(UNIT = iunepmatcb, FILE = filint, IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'old', ACCESS = 'direct', RECL = unf_recl)
#endif
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in opening X.epmatkqcb1', 1)
      !
#if defined(__MPI)
      ! Offset depending on CPU
      lrepmatw2 = INT(lower_bnd - 1_MPI_OFFSET_KIND, KIND = MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND
      !
      ! Size of what we read
      lsize = INT(nind, KIND = MPI_OFFSET_KIND)
      !
      CALL MPI_FILE_SEEK(iunepmatcb, lrepmatw2, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK iunepmatcb', 1)
      CALL MPI_FILE_READ(iunepmatcb, trans_probcb(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ iunepmatcb', 1)
      !
      ! Now read the sparse matrix mapping
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparseqcb' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparseqcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparseqcb', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparsekcb' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsekcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparsekcb', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparseicb' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparseicb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparseicb', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparsejcb' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsejcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparsejcb', 1)
      CALL MPI_FILE_OPEN(inter_pool_comm, 'sparsetcb' // '_' // TRIM(my_image_id_ch), &
        MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsetcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_OPEN sparsetcb', 1)
#else
      READ(UNIT = iunepmatcb, REC = 1, IOSTAT = ierr) trans_probcb
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading X.epmatkq1', 1)
      !
      ! Now open the sparse matrix mapping
      INQUIRE(IOLENGTH = direct_io_factor) dum_int
      unf_recl = direct_io_factor * INT(nind, KIND = KIND(unf_recl))
      OPEN(UNIT = iunsparseqcb, FILE = 'sparseqcb', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparseqcb', 1)
      OPEN(UNIT = iunsparsekcb, FILE = 'sparsekcb', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparsekcb', 1)
      OPEN(UNIT = iunsparseicb, FILE = 'sparseicb', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparseicb', 1)
      OPEN(UNIT = iunsparsejcb, FILE = 'sparsejcb', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparsejcb', 1)
      OPEN(UNIT = iunsparsetcb, FILE = 'sparsetcb', IOSTAT = ierr, FORM = 'unformatted', &
           STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error in reading sparsetcb', 1)
#endif
      !
      ALLOCATE(sparsecb_q(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparsecb_q', 1)
      ALLOCATE(sparsecb_k(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparsecb_k', 1)
      ALLOCATE(sparsecb_i(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparsecb_i', 1)
      ALLOCATE(sparsecb_j(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparsecb_j', 1)
      ALLOCATE(sparsecb_t(nind), STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error allocating sparsecb_t', 1)
      sparsecb_q(:) = 0.0d0
      sparsecb_k(:) = 0.0d0
      sparsecb_i(:) = 0.0d0
      sparsecb_j(:) = 0.0d0
      sparsecb_t(:) = 0.0d0
      !
#if defined(__MPI)
      lrepmatw4 = INT(lower_bnd - 1_MPI_OFFSET_KIND, KIND = MPI_OFFSET_KIND) * 4_MPI_OFFSET_KIND
      !
      CALL MPI_FILE_SEEK(iunsparseqcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK iunsparseqcb', 1)
      CALL MPI_FILE_READ(iunsparseqcb, sparsecb_q(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ iunsparseqcb', 1)
      CALL MPI_FILE_SEEK(iunsparsekcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK iunsparsekcb', 1)
      CALL MPI_FILE_READ(iunsparsekcb, sparsecb_k(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ iunsparsekcb', 1)
      CALL MPI_FILE_SEEK(iunsparseicb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK iunsparseicb', 1)
      CALL MPI_FILE_READ(iunsparseicb, sparsecb_i(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ iunsparseicb', 1)
      CALL MPI_FILE_SEEK(iunsparsejcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK iunsparsejcb', 1)
      CALL MPI_FILE_READ(iunsparsejcb, sparsecb_j(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ iunsparsejcb', 1)
      CALL MPI_FILE_SEEK(iunsparsetcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_SEEK iunsparsetcb', 1)
      CALL MPI_FILE_READ(iunsparsetcb, sparsecb_t(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_READ iunsparsetcb', 1)
#else
      READ(iunsparseqcb, REC = 1, IOSTAT = ierr) sparsecb_q
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparsecb_q', 1)
      READ(iunsparsekcb, REC = 1, IOSTAT = ierr) sparsecb_k
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparsecb_k', 1)
      READ(iunsparseicb, REC = 1, IOSTAT = ierr) sparsecb_i
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparsecb_i', 1)
      READ(iunsparsejcb, REC = 1, IOSTAT = ierr) sparsecb_j
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparsecb_j', 1)
      READ(iunsparsetcb, REC = 1, IOSTAT = ierr) sparsecb_t
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in reading sparsecb_t', 1)
#endif
      !
      CALL ibte(nind, etf_all, vkk_all, wkf_all, trans_probcb, efcb, &
                sparsecb_q, sparsecb_k, sparsecb_i, sparsecb_j, sparsecb_t, inv_tau_allcb)
      !
#if defined(__MPI)
      CALL MPI_FILE_CLOSE(iunepmatcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparseqcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparsekcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparseicb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparsejcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
      CALL MPI_FILE_CLOSE(iunsparsetcb, ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in MPI_FILE_CLOSE', 1)
#else
      CLOSE(iunepmatcb, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing X.epmatkqcb1', 1)
      CLOSE(iunsparseqcb, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparseqcb', 1)
      CLOSE(iunsparsekcb, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparsekcb', 1)
      CLOSE(iunsparseicb, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparseicb', 1)
      CLOSE(iunsparsejcb, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparsejcb', 1)
      CLOSE(iunsparsetcb, STATUS = 'keep', IOSTAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'error in closing sparsetcb', 1)
#endif
      DEALLOCATE(trans_probcb, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating trans_probcb', 1)
      DEALLOCATE(sparsecb_q, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparsecb_q', 1)
      DEALLOCATE(sparsecb_k, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparsecb_k', 1)
      DEALLOCATE(sparsecb_i, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparsecb_i', 1)
      DEALLOCATE(sparsecb_j, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparsecb_j', 1)
      DEALLOCATE(sparsecb_t, STAT = ierr)
      IF (ierr /= 0) CALL errore('transport_ibte', 'Error deallocating sparsecb_t', 1)
      !
    ENDIF
    !
    !----------------------------------------------------------------------------
    END SUBROUTINE transport_ibte
    !----------------------------------------------------------------------------
    !
  !------------------------------------------------------------------------------
  END MODULE transport
  !------------------------------------------------------------------------------
