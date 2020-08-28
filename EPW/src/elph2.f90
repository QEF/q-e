  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !--------------------------------------------------------------------------
  MODULE elph2
  !--------------------------------------------------------------------------
  !!
  !! This module contains all global variables used in EPW.
  !!
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL ::                 &!
    elph,                    &!  Perform electron-phonon interpolation
    adapt_smearing            !  Adaptative smearing
  LOGICAL, ALLOCATABLE ::    &!
    lwin(:, :),              &!  identify bands within outer energy windows (when disentanglement is used)
    lwinq(:, :),             &!  Excluded bands
    exband(:),               &!  k-point independent list of bands excluded from the calculation of overlap and projection matrices in W90
    done_elph(:)
  INTEGER ::                 &!
    nkqf,                    &!  number of k+q points per pool (fine grid)
    nkf,                     &!  number of k points per pool (fine grid)
    nqf,                     &!  number of q points per pool (fine grid)
    nkqtotf,                 &!  total number of k+q points (fine grid)
    nktotf,                  &!  total number of k points (fine grid)
    nqtotf,                  &!  total number of q points (fine grid)
    nrr,                     &!  number of wigner-seitz points (elec interp only)
    ibndmin,                 &!  Lower band bound for slimming down electron-phonon matrix
    ibndmax,                 &!  Upper band bound for slimming down electron-phonon matrix
    lower_band,              &!  Lower band index for image (band) parallelization
    upper_band,              &!  Upper band index for image (band) parallelization
    lower_bnd,               &!  Lower bound for the k-depend index among the mpi pools
    upper_bnd,               &!  Upper bound for the k-depend index among the mpi pools
    lrepmatw2_merge,         &!  File merging dimensions
    lrepmatw5_merge,         &!  File merging dimensions
    nbndfst,                 &!  Number of bands within the fsthick window.
    ibndstart,               &!  First band index after excluding bands in Wannierization step
    ibndend,                 &!  Last band index after excluding bands in Wannierization step
    nbndep,                  &!  Number of remaining bands after excluding bands in Wannierizatin step
    nbndskip,                &!  Number of bands to be skipped in Wannierization step, leading to
                              !  the exclusion from the original Hamiltonian
    ngxx,                    &!  Maximum number of G-vectors over all pools
    ngxxf,                   &!  Maximum number of G-vectors over all pools for k+q folding
    ig_s,                    &!  First G index within each core in case of G parallelization
    ig_e,                    &!  Last G index within each core in case of G parallelization
    num_wannier_plot,        &!  Number of Wannier functions to plot
    ng0vec                    ! number of inequivalent such translations (125)
  INTEGER, ALLOCATABLE ::    &!
    igk(:),                  &!  Index for k+G vector
    igkq(:),                 &!  Index for k+q+G vector
    igk_k_all(:, :),         &!  Global index (in case of parallel)
    ngk_all(:),              &!  Global number of plane wave for each global k-point
    map_rebal(:),            &!  Map between the k-point and their load rebalanced one
    map_rebal_inv(:),        &!  Map between the k-point and their load rebalanced one
    shift(:),                &!  for every k+q, index of the G0 which folds k+q into k+q+G0 of the first BZ
    gmap(:),                 &!  the map G -> G-G_0 in the large (density) G vectors set, for every G_0
    ixkqf_tr(:),             &!  Mapping matrix from k+q (where q is full BZ) to IBZ
    s_bztoibz_full(:),       &!  Rotation that brink that k-point from BZ to IBZ
    mapg(:,:,:),             &!  Map between FFT grid and G-vector index
    wanplotlist(:)            !  List of Wannier functions to plot
  REAL(KIND = DP) ::         &!
    efnew,                   &!  Fermi level on the fine grid. Added globaly for efficiency reason
    deltaq,                  &!  Displacement of fine-mesh k-points for velocity corrections
    threshold,               &!  Threshold below which the transition probabilities are not written to file in transport.
    area,                    &!  Area of the 2D unit cell.
    g0vec_all_r(3, 125)       ! G-vectors needed to fold the k+q grid into the k grid, cartesian coord.
  REAL(KIND = DP), ALLOCATABLE ::&
    a_all(:, :, :),          &!  electronic spectral function du to electron-phonon interaction
    a_all_ph(:, :, :),       &!  phononic spectral function du to electron-phonon interaction
    dos(:),                  &!  Density of states at the chemical potential.
    et_ks(:, :),             &!  lda eigenvalues
    xkq(:, :),               &!  local k+q grid, coarse (3, nks)
    etq(:, :),               &!  eigenvalues of k+q wavefunctions
    xkf(:, :),               &!  fine k point grid (3, nkqf)
    xkfd(:, :, :),           &!  fine k point grid, displaced along each cartesian direction +/delta_q (3, nkqf, 6)
    wkf(:),                  &!  weights on the fine grid (nkqf)
    xqf(:, :),               &!  fine q point grid
    wqf(:),                  &!  weights on the fine q grid
    etf(:, :),               &!  interpolated eigenvalues (nbnd, nkqf)
    etf_k(:, :),             &!  Saved interpolated KS eigenenergies for later used in q-parallelization (nbnd, nkqf)
    etf_ks(:, :),            &!  interpolated eigenvalues (nbnd, nkqf) KS eigenvalues in the case of eig_read
    wf(:, :),                &!  interpolated eigenfrequencies
    gamma_all(:, :, :, :),      &!  Gamma
    gamma_nest(:, :),        &!  Nesting function in the case of q-parallelization
    gamma_v_all(:, :, :, :),    &!  Gamma
    lambda_all(:, :, :, :),     &!  Electron-phonon coupling parameter
    lambda_v_all(:, :, :, :),   &!  Electron-phonon coupling parameter (transport)
    sigmar_all(:, :, :),        &!  Real part of the electron-phonon self-energy
    sigmai_all(:, :, :),        &!  Imaginary part of the electron-phonon self-energy
    sigmai_mode(:, :, :, :),    &!  Mode resolved imaginary electron self-energy
    zi_all(:, :, :),            &!  Z renormalization factor
    eta(:, :, :),            &!  Adaptative smearing
    esigmar_all(:, :, :, :), &!  energy of the real self-energy
    esigmai_all(:, :, :, :), &!  energy of the imaginary self-energy
    jdos(:),                 &!  j-DOS
    spectra(:, :, :, :, :, :), &!  dipole absorption spectra, polarizations, nomega, nsmear, dme/vme, absorption/emission
    zstar(:, :, :),          &!  Born effective charges
    epsi(:, :),              &!  dielectric tensor
    inv_tau_all(:, :, :),    &!  scattering rate
    inv_tau_allcb(:, :, :),  &!  Second scattering rate (for both)
    zi_allvb(:, :, :),       &!  Z-factor in scattering rate
    zi_allcb(:, :, :),       &!  Second Z-factor in scattering rate (for both VB and CB calculations)
    ifc(:,:,:,:,:,:,:),      &!  Interatomic force constant in real space
    omegap(:),               &!  Photon energy for phonon-assisted absorption
    epsilon2_abs(:, :, :, :),   &!  Imaginary part of dielectric function for phonon-assisted absorption, vs omega, vs broadening
    wscache(:, :, :, :, :),  &!  Use as cache when doing IFC when lifc = .TRUE.
    epsilon2_abs_lorenz(:, :, :, :), &! Imaginary part of dielectric function for phonon-assisted absorption, vs omega, vs broadening
    gtemp(:),                &!  Temperature used globally (units of Ry)
    mobilityh_save(:),       &!  Error in the hole mobility
    mobilityel_save(:)        !  Error in the electron mobility
  COMPLEX(KIND = DP), ALLOCATABLE :: &
    el_ph_mat(:, :, :, :),    &!  e-p matrix  (nbnd, nbnd, nks, 3*nat)
    cu(:, :, :),              &!  rot matrix for wannier interpolation of k point, coarse mesh (nbnd*nbnd*nkstot)
    cuq(:, :, :),             &!  rot matrix for wannier interpolation of k+q point, coarse mesh (nbnd*nbnd*nkstot)
    chw(:, :, :),             &!  Hamiltonian in wannier basis
    chw_ks(:, :, :),          &!  Hamiltonian in wannier basis  (Kohn-Sham eigenvalues if many-body eigenvalues are read in)
    cdmew(:, :, :, :),        &!  Dipole matrix in wannier basis
    cvmew(:, :, :, :),        &!  Velocity matrix in wannier basis
    rdw(:, :, :),             &!  dynamical matrix in wannier basis (real)
    epmatwp(:, :, :, :, :),   &!  e-p matrix  in wannier basis - electrons and phonons
    umat(:, :, :),            &!  the rotation matrix for the unique setting of the wfs gauge -- on the local pool
    umatq(:, :, :),           &!  the rotation matrix for the unique setting of the wfs gauge for the k + q-- on the local pool
    umat_all(:, :, :),        &!  the rotation matrix for the unique setting of the wfs gauge -- for all k points
    umatq_all(:, :, :),       &!  the rotation matrix for the unique setting of the wfs gauge -- for all k+q points
    dynq(:, :, :),            &!  dynamical matrix for every q (nmode, nmodes, nqtot)
    epmatq(:, :, :, :, :),    &!  e-p matrix for every q (nbnd, nbnd, nks, nmodes, nqtot)
    epf17(:, :, :, :),        &!  full ep matrix in bloch rep stored in mem (nkqtotf, nbnd, nbnd, nmodes)-nbnd inside wndw
    dmec(:, :, :, :),         &!  dipole matrix elements on the coarse mesh (ipol, nbnd, nbnd, nks)
    dmef(:, :, :, :),         &!  dipole matrix elements on the fine   mesh (ipol, nbnd, nbnd, nks)
    vmef(:, :, :, :),         &!  velocity matrix elements on the fine mesh (ipol, nbnd, nbnd, nks)
    bmat(:, :, :, :),         &!  overlap U_k+q U_k^\dagger on the coarse mesh (nbnd, nbnd, nks, nqtot)
    eps_rpa(:),               &!  rpa model screening
    veff(:, :)                 !  effective potential
  !
  ! -------------------------------------------------------------------------
  ! Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
  ! Shell implementation for future use
  INTEGER ::                  &! Add descriptions
     N1_dim,                  &
     hh,                      &
     Np
   REAL(KIND = DP), ALLOCATABLE ::&
     Ec(:),                   &
     ekf(:),                  &
     etfq(:,:),               &
     etf_qs(:,:)
   COMPLEX(KIND = DP), ALLOCATABLE :: &
     g2_4  (:, :, :, :),     &
     g2_all(:,:,:,:,:),      &
     Ac(:),                  &
     Hkk(:,:),               &
     H_copy(:,:),            &
     gq(:),                  &
     dtau(:,:,:)
  ! -------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  END MODULE elph2
  !--------------------------------------------------------------------------
