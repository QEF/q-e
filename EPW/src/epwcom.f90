  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from the code PH/phcom - Quantum-ESPRESSO group
  !-----------------------------------------------------------------------
  MODULE control_epw
  !-----------------------------------------------------------------------
  !!
  !! Common variables for the epw program
  !!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npk
  !
  SAVE
  !
  CHARACTER(LEN = 100) :: dvscf_dir ='./'
  !! directory for .dvscf and .dyn files (wannier interpolation)
  CHARACTER(LEN = 80) :: fileig
  !! output file for the electron-phonon coefficients
  CHARACTER(LEN = 256) :: proj(200)
  !! projections for W90
  CHARACTER(LEN = 80) :: scdm_entanglement
  !! disentanglement type in the SCDM algorithm
  CHARACTER(LEN = 256) :: bands_skipped
  !! k-point independent list of bands excluded from the calculation
  !! of overlap and projection matrices in W90
  CHARACTER(LEN = 256) :: wdata(200)
  !! any extra info for W90
  CHARACTER(LEN = 75) :: title
  !! ...  title of the simulation
  CHARACTER(LEN = 10)  :: asr_typ
  !! type of ASR if lifc=.TRUE.
  CHARACTER(LEN = 10) :: vme
  !! if 'dipole' then computes the velocity as dipole+commutator = <\psi_mk|p+i[V_NL,r]|\psi_nk>
  !! if 'wannier' then computes the velocity as dH_nmk/dk - i(e_nk-e_mk)A_nmk where A is the Berry connection
  !!!!!
  CHARACTER(LEN = 10) :: tc_linear_solver
  !! algorithm to solve T_c eigenvalue problem
  !!!!!
  !
  LOGICAL :: elecselfen
  !! if .TRUE. calculate electron selfenergy due to e-p interaction
  LOGICAL :: phonselfen
  !! if .TRUE. calculate phonon selfenergy due to e-p interaction
  LOGICAL :: plselfen
  !! if .TRUE. calculate the electron-plason self-energy
  LOGICAL :: epbread
  !! if .TRUE. read epmatq from files .epb
  LOGICAL :: epbwrite
  !! if .TRUE. write epmatq to files .epb
  LOGICAL :: epwread
  !! if .TRUE. read all quantities in Wannier representation from file epwdata.fmt
  LOGICAL :: epwwrite
  !! if .TRUE. write all quantities in Wannier representation to file epwdata.fmt
  LOGICAL :: restart
  !! if .TRUE. restart a calculation stopped during the interpolation phase from XXX.restart file.
  LOGICAL :: specfun_el
  !! if .TRUE. calculate spectral electron function due to e-p interaction
  LOGICAL :: specfun_ph
  !! if .TRUE. calculate spectral phonon function due to e-p interaction
  LOGICAL :: specfun_pl
  !! if .TRUE. calculate plasmon spectral function
  LOGICAL :: wannierize
  !! if .TRUE. run the wannier90 code
  LOGICAL :: auto_projections
  !! to automatically generate initial projections for W90
  LOGICAL :: scdm_proj
  !! wannier functions via the SCDM algorithm
  LOGICAL :: a2f
  !! if .TRUE. calculate Eliashberg spectral electron function from selfen_phon
  LOGICAL :: write_wfn
  !! if .TRUE. write UNK files in wannier90
  LOGICAL :: kmaps
  !! if .TRUE. read kmap and kgmap from disk.
  LOGICAL :: nest_fn
  !! if .TRUE. calculate the electronic nesting function (metals only)
  LOGICAL :: rand_q
  !! if .TRUE. use random points for the fine q-mesh
  LOGICAL :: rand_k
  !! if .TRUE. use random points for the fine k-mesh
  LOGICAL :: mp_mesh_q
  !! if .TRUE. use points in the irreducible wedge for the uniform fine q-mesh
  LOGICAL :: mp_mesh_k
  !! if .TRUE. use points in the irreducible wedge for the uniform fine k-mesh
  LOGICAL :: eig_read
  !! if .TRUE. read a set of electronic eigenvalues in eV to replace the calcualted ones
  LOGICAL :: wepexst
  !! if .TRUE. prefix.epmatwe files are already on disk. don't recalculate. debugging param
  LOGICAL :: epexst
  !! if .TRUE. prefix.epmatwp files are already on disk. don't recalculate. debugging param
  LOGICAL :: band_plot
  !! if .TRUE. write files to plot band structure and phonon dispersion
  LOGICAL :: fermi_plot
  !! if .TRUE. write files to plot Fermi surface
  LOGICAL :: lpolar
  !! if .TRUE. enable the correct Wannier interpolation in the case of polar material.
  LOGICAL :: lscreen
  !! if .TRUE. the e-ph matrix elements are screened by the RPA or TF dielectric function
  LOGICAL :: lifc
  !! if .TRUE. read interatomic force constants produced by q2r.x for phonon interpolation
  LOGICAL :: cumulant
  !! if .TRUE. calculate the electron spectral function using the cumulant expansion method
  LOGICAL :: delta_approx
  !! if .TRUE. the double delta approximation is used for the phonon self energy
  LOGICAL :: ep_coupling
  !! if .TRUE. run e-p coupling calculation
  LOGICAL :: efermi_read
  !! if .TRUE. fermi energy is read from the input file
  LOGICAL :: system_2d
  !! if .TRUE. the system is 2 dimensional (vaccum is in z-direction)
  LOGICAL :: prtgkk
  !! if .TRUE. print the |g| vertex in [meV].
  LOGICAL :: lphase
  !! if .TRUE. fix the gauge when diagonalizing the interpolated dynamical matrix and electronic Hamiltonian.
  LOGICAL :: lrot
  !! if .TRUE. fix the gauge when diagonalizing the interpolated dynamical matrix and electronic Hamiltonian.
  LOGICAL :: lindabs
  !! if .TRUE. perform phonon-assisted absorption calculations
  LOGICAL :: use_ws
  !! if .TRUE. use Wannier-centers to compute the Wigner-Seitz cell.
  LOGICAL :: epmatkqread
  !! if .TRUE. restart and IBTE calculation from the scattering rates written to files.
  LOGICAL :: selecqread
  !! if .TRUE. restart from the selecq.fmt file
  LOGICAL :: assume_metal
  !! if .TRUE. we are dealing with a metal.
  !
  ! Superconductivity
  LOGICAL :: ephwrite
  !! if .TRUE. write e-ph matrix elements on the fine mesh to file
  LOGICAL :: lreal
  !! if .TRUE. solve real-axis Eliashberg eqautions
  LOGICAL :: limag
  !! if .TRUE. solve imag-axis Eliashberg eqautions
  LOGICAL :: lpade
  !! if .TRUE. use pade approximants to continue imag-axis Eliashberg equtions
  !to real-axis
  LOGICAL :: lacon
  !! if .TRUE. use analytic continuation to continue imag-axis Eliashberg
  !equtions to real-axis
  LOGICAL :: liso
  !! if .TRUE. solve isotropic case
  LOGICAL :: laniso
  !! if .TRUE. solve anisotropic case
  !!!!!
  !  LOGICAL :: lunif
  !  !! if .TRUE. a uniform grid is defined between wsfc and wc for real-axis calculations
  !!!!!
  LOGICAL :: kerwrite
  !! if .TRUE. write kp and km to files .ker for real-axis calculations
  LOGICAL :: kerread
  !! if .TRUE. read kp and km from files .ker for real-axis calculations
  LOGICAL :: imag_read
  !! if .TRUE. read from file Delta and Znorm on the imaginary-axis
  LOGICAL :: eliashberg
  !! if .TRUE. solve the Eliashberg equations
  LOGICAL :: tc_linear
  !! if .TRUE. linearized Eliashberg eqn. for T_c will be solved
  !!!!!
  ! CHARACTER(LEN = 10) :: tc_linear_solver
  ! !! algorithm to solve T_c eigenvalue problem
  LOGICAL :: fbw
  !! if .TRUE. full-bandwidth calculations will be performed
  LOGICAL :: muchem
  !! if .TURE. the chem. pot. is updated in fbw calculations
  !!!!!
  !
  ! Conductivity
  LOGICAL :: scattering
  !! if .TRUE. scattering rates are calculated
  LOGICAL :: scattering_serta
  !! if .TRUE. scattering rates are calculated using self-energy relaxation-time-approx
  LOGICAL :: scatread
  !! if .TRUE. the scattering rates are read from file.
  LOGICAL :: scattering_0rta
  !! if .TRUE. scattering rates are calculated using 0th order relaxation-time-approx
  LOGICAL :: int_mob
  !! if .TRUE. compute the intrinsic mobilities. This means that the electron and hole carrier density is equal.
  LOGICAL :: iterative_bte
  !! if .TRUE. the iterative solution for BTE is computed. A first run with scattering_serta = .TRUE. is required.
  LOGICAL :: carrier
  !! if .TRUE. compute the doped electronic mobilities.
  LOGICAL :: longrange
  !! if .TRUE. compute the long range interaction of el-ph. Can only be .TRUE. if lpolar is also true.
  LOGICAL :: shortrange
  !! if .TRUE. compute the long range interaction of el-ph. Can only be .TRUE. if lpolar is also true.
  LOGICAL :: fixsym
  !! if .TRUE. try to fix the symmetry-related issues
  LOGICAL :: epw_no_t_rev
  !! if .TRUE. set t_rev = 0
  LOGICAL :: epw_tr
  !! set time_reversal to epw_tr
  LOGICAL :: epw_nosym
  !! If .TRUE. set nsym=1 (nrot=1) and invsym=.FALSE.
  LOGICAL :: epw_noinv
  !! If .TRUE. set imq to non-zero value
  LOGICAL :: epw_crysym
  !! If .TRUE. calculate the symmetry of the crystal
  !
  INTEGER :: ngaussw
  !! smearing type for Fermi surface average in e-ph coupling after wann. interp.
  INTEGER :: nw
  !! nr. of bins for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  INTEGER :: nbndsub
  !! nr. of bands in the optimal subspace (when disentanglement is used)
  INTEGER :: num_iter
  !! nr. of steps used in creating the maximally localized Wannier functions
  INTEGER :: iprint
  !! verbosity of the wannier90 code
  INTEGER :: nsmear
  !! nr. of temperature-like smearings calculated in selfen_phon
  INTEGER :: rand_nq
  !! use random points for the fine q-mesh
  INTEGER :: rand_nk
  !! use random points for the fine k-mesh
  INTEGER :: nqc1, nqc2, nqc3
  !! qx,qy,qz sizes of the uniform phonon coarse mesh to be used
  INTEGER :: nkc1, nkc2, nkc3
  !! kx,ky,kz sizes of the uniform electron coarse mesh to be used
  INTEGER :: nqf1, nqf2, nqf3
  !! qx,qy,qz sizes of the uniform phonon fine mesh to be used
  INTEGER :: nkf1, nkf2, nkf3
  !! kx,ky,kz sizes of the uniform electron fine mesh to be used
  INTEGER :: nqsmear
  !! nr. of smearings used to calculate a2f
  INTEGER :: nqstep
  !! nr. of steps used to calculate a2f
  INTEGER :: iswitch
  !! Switch for symmetry operations
  INTEGER :: nwanxx = 200
  !! parameter used in writing prefix.win file.
  INTEGER :: ntempxx = 50
  !! Maximum number of temperatures
  INTEGER :: etf_mem
  !! If 0, all in memory. If 1, less is stored in memory (read files).
  INTEGER :: scr_typ
  !! If 0 calculates the Lindhard screening, if 1 the Thomas-Fermi screening
  INTEGER :: bnd_cum
  !! band index for which the cumulant calculation is done
  INTEGER :: mob_maxiter
  !! Maximum number of iteration for the IBTE
  INTEGER :: nstemp
  !! nr. of temperature points for temperature dependent caclulations
  REAL(KIND = DP) :: tempsmin
  !! min. temperature in Eliashberg equations - deprecated as an input parameter
  REAL(KIND = DP) :: tempsmax
  !! max. temperature - deprecated as an input parameter
  REAL(KIND = DP) :: temps(50)
  !! input temperature array (units of Kelvin)
  !
  ! Superconductivity
  INTEGER :: nswi
  !! nr. of grid points for Eliashberg equations of imaginary axis
  INTEGER :: nsiter
  !! nr. of iterations for self-consistency
  INTEGER :: npade
  !! percentange of Matsubara points used in Pade continuation
  INTEGER :: broyden_ndim
  !! nr. of iterations used in broyden mixing scheme
  INTEGER :: nw_specfun
  !! nr. of bins for frequency in electron spectral function due to e-p interaction
  INTEGER :: restart_step
  !! Create a restart point during the interpolation part every restart_step q/k-points.
  !!!!!
  INTEGER :: gridsamp
  !! Type of the Matsubara freq. sampling (-1= read from file; 0= uniform; 1= sparse)
  !!!!!
  !
  REAL(KIND = DP) :: degaussw
  !! smearing width for Fermi surface average in e-ph coupling after wann interp
  REAL(KIND = DP) :: fsthick
  !! thickness of the Fermi shell for averaging the e-ph matrix element
  REAL(KIND = DP) :: wmin
  !! min frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  REAL(KIND = DP) :: wmax
  !! max frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  REAL(KIND = DP) :: delta_smear
  !! change in energy for each additional smearing in the selfen_phon
  !
  ! Wannierization
  CHARACTER(LEN = 255) :: wannier_plot_list
  !! Field read for parsing Wannier function list
  LOGICAL :: wannier_plot
  !! if .TRUE. plot Wannier functions
  LOGICAL :: reduce_unk
  !! if .TRUE. plot Wannier functions on reduced grids
  INTEGER :: wannier_plot_supercell(3)
  !! Size of supercell for plotting Wannier functions
  REAL(KIND = DP) :: dis_win_min
  !! min energy of the Wannier disentanglement window
  REAL(KIND = DP) :: dis_win_max
  !! max energy of the Wannier disentanglement window
  REAL(KIND = DP) :: dis_froz_min
  !! min energy of the frozen Wannier disentanglement window
  REAL(KIND = DP) :: dis_froz_max
  !! max energy of the frozen Wannier disentanglement window
  REAL(KIND = DP) :: scdm_mu
  !! parameter for Wannier functions via SCDM algorithm
  REAL(KIND = DP) :: scdm_sigma
  !! parameter for Wannier functions via SCDM algorithm
  REAL(KIND = DP) :: wannier_plot_scale
  !! Scaling parameter for cube files
  REAL(KIND = DP) :: wannier_plot_radius
  !! Cut-off radius for plotting Wannier functions
  !
  ! Superconductivity
  REAL(KIND = DP) :: eps_acustic
  !! min. phonon frequency for e-p and a2f calculations
  REAL(KIND = DP) :: degaussq
  !! smearing for sum over q in e-ph coupling
  REAL(KIND = DP) :: delta_qsmear
  !! change in energy for each additional smearing in the a2f
  REAL(KIND = DP) :: muc
  !! effective Coulomb potential in Eliashberg equations
  !!!!!
  !  REAL(KIND = DP) :: wsfc
  !  !! intermediate freqeuncy between (0,wscut)
  !  REAL(KIND = DP) :: pwc
  !  !! power used to define a non-uniform grid between wsfc and wscut
  !!!!!
  REAL(KIND = DP) :: wscut
  !! upper limit cutoff frequency in Eliashberg equations (at least 5 times wsphmax)
  REAL(KIND = DP) :: broyden_beta
  !! mixing factor for broyden mixing
  REAL(KIND = DP) :: conv_thr_raxis
  !! convergence threshold for iterative solution of real-axis Eliashberg equations
  REAL(KIND = DP) :: conv_thr_iaxis
  !! convergence threshold for iterative solution of imag-axis Eliashberg equations
  REAL(KIND = DP) :: conv_thr_racon
  !! convergence threshold for iterative solution of analytic continuation of
  !! Eliashberg equations from imag- to real-axis
  REAL(KIND = DP) :: gap_edge
  !! initial guess of the superconducting gap
  REAL(KIND = DP) :: max_memlt
  !! maximum memory that can be allocated per pool
  REAL(KIND = DP) :: fermi_energy
  !! fermi energy is given in the input file
  REAL(KIND = DP) :: wmin_specfun
  !! min frequency in electron spectral function due to e-p interaction
  REAL(KIND = DP) :: wmax_specfun
  !! max frequency in electron spectral function due to e-p `interaction
  !!!!!
  REAL(KIND = DP) :: dos_del
  !! Delta_E in electronic dos for Fermi window (in eV)
  REAL(KIND = DP) :: griddens
  !! Measure of sparsity of the grid
  !!!!!
  !
  ! Conductivity
  INTEGER :: mob_nfreq
  !! Number of frequency for the spectral decomposition of mobility
  REAL(KIND = DP) :: scissor
  !! Value of the scissor shift in eV.
  REAL(KIND = DP) :: ncarrier
  !! Amount of carrier concentration in cm^-3 when doping a semiconductors
  REAL(KIND = DP) :: nc
  !! Number of carrier per unit cell that participate to the conduction in the Ziman resistivity formula
  REAL(KIND = DP) :: bfieldx
  !! Magnetic field along the x-direction
  REAL(KIND = DP) :: bfieldy
  !! Magnetic field along the y-direction
  REAL(KIND = DP) :: bfieldz
  !! Magnetic field along the z-direction
  REAL(KIND = DP) :: mob_maxfreq
  !! Maximum frequency for the spectral decomposition of mobility. Typically that freq. is the highest phonon freq.
  !!!!!
  !
  ! Conductivity - ionized impurity scattering
  LOGICAL :: ii_g
  !! If .TRUE. calculate the ionized impurity matrix elements
  LOGICAL :: ii_scattering
  !! If .TRUE. calculate the carrier-impurity scattering rate
  LOGICAL :: ii_only
  !! If .TRUE., calculate only the ionized-impurity-limitted mobility
  LOGICAL :: ii_lscreen
  !! If .TRUE. (DEFAULT), calculate and use free-carrier
  !! screening of ii_g matrix elements, if .FALSE., calculations will diverge with 
  !! increasing density of k-point due to 1/(q^4) divergence for ionized impurity 
  !! scattering rate
  LOGICAL :: ii_partion
  !! If .TRUE., account for partial ionization 
  REAL(KIND = DP) :: ii_charge
  !! charge of the ionized impurities, units of electron charge for input
  REAL(KIND = DP) :: ii_n
  !! density of impurities, input in units of cm^-3
  REAL(KIND = DP) :: ii_eda
  !! ionization energy of the donor or accpetor impurity, for partial ionization calcs
  REAL(KIND = DP) :: ii_eps0
  !! low frequency dielectric constant
  !!!!! 
  !
  ! Plasmon
  REAL(KIND = DP) :: nel
  !! fractional number of electrons in the unit cell
  REAL(KIND = DP) :: meff
  !! Density-of-state effective mass (in unit of the electron mass)
  REAL(KIND = DP) :: epsiheg
  !! Dielectric constant at zero doping.
  REAL(KIND = DP) :: fermi_diff
  !! difference between Fermi energy and band edge (in eV)
  REAL(KIND = DP) :: smear_rpa
  !! smearing for the calculation of the Lindhard function (in eV)
  !
  ! Phonon-assisted absorption
  INTEGER :: neta = 9
  !! Number of broadening parameters
  INTEGER :: nomega
  !! Number of frequency (photon energy) points
  REAL(KIND = DP) :: omegamin
  !! Photon energy minimum (in eV)
  REAL(KIND = DP) :: omegamax
  !! Photon energy maximum (in eV)
  REAL(KIND = DP) :: omegastep
  !! Photon energy step (in eV)
  REAL(KIND = DP) :: sigma_ref
  !! Reference conductivity for Drude dielectric function
  !
  ! ----------------------------------------------------------------------------------
  ! Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
  ! Shell implementation for future use.
  INTEGER :: start_band_plrn, end_band_plrn
  !! Start and end band index in matrix element
  INTEGER :: nstate_plrn
  !! Number of polaron states calculated
  INTEGER :: ndos_plrn
  !! Number of grid in polaron DOS calculation
  INTEGER :: type_plrn
  !! polaron type (electron/hole)
  INTEGER :: init_plrn
  !! initial polaron wavefuntion with 1:Gaussian package and 2:Random number
  INTEGER :: niter_plrn
  !! Maximum number of polaron SCF loops.
  INTEGER :: seed_plrn
  !! The seed number to generate the random initial polaron wavefunction
  INTEGER :: nhblock_plrn
  !! The seed number to generate the random initial polaron wavefunction
  INTEGER :: g_start_band_plrn, g_end_band_plrn
  !! Start and end band in saving g matrix
  INTEGER :: step_wf_grid_plrn
  !! number of grid to skip in output of real space wavefunction
  INTEGER :: io_lvl_plrn
  !! number of grid to skip in output of real space wavefunction
  INTEGER :: scell_mat(3,3)
  !! Supercell transformation matrix for polaron
  INTEGER :: init_ntau_plrn
  !! Number of displacements to be read from file to calculate polaron energy landscape  
  INTEGER :: nethrdg_plrn
  !! Number of steps in the diagonalization if adaptive_ethrdg_plrn=.true.
  REAL(KIND = DP) :: conv_thr_plrn
  !!  convergent threshold for polaron calculation
  REAL(KIND = DP) :: edos_min_plrn, edos_max_plrn, edos_sigma_plrn
  !! Electron Energy range in polaron DOS calculation
  REAL(KIND = DP) :: pdos_min_plrn, pdos_max_plrn, pdos_sigma_plrn
  !! Phonon Energy range in polaron DOS calculation
  REAL(KIND = DP) :: init_sigma_plrn
  !! decay radius of polaron wavefunction in Gaussian initialization
  REAL(KIND = DP) :: init_k0_plrn(3)
  !! center k of Gaussian initialization
  REAL(KIND = DP) :: mixing_plrn
  !! mixing weight in plrn iteration
  REAL(KIND = DP) :: ethrdg_plrn
  !! threshold in diagonalization of the polaron Hamiltonian
  REAL(KIND = DP) :: beta_plrn(3)
  !! Mixing weight in Self-consistency
  REAL(KIND = DP) :: r0_plrn(3)
  !! Mixing weight in Self-consistency
  REAL(KIND = DP) :: g_start_energy_plrn
  !! Mixing weight in Self-consistency
  REAL(KIND = DP) :: g_end_energy_plrn
  !! Mixing weight in Self-consistency
  REAL(KIND = DP) :: kappa_plrn
  !! Mixing weight in Self-consistency
  REAL(KIND = DP) :: g_power_order_plrn
  !! kappa
  REAL(KIND = DP) :: omega_LO_plrn
  !! Higher optical phonon frequency
  REAL(KIND = DP) :: m_eff_plrn
  !! effective mass in LP model
  REAL(KIND = DP) :: g_tol_plrn
  !!
  REAL(KIND = DP) :: as(3,3)
  !! Supercell lattice vectors transformed by Smat, as=S*at
  REAL(KIND = DP) :: bs(3,3)
  !! Supercell reciprocal lattice vectors transformed by Sbar=transpose(inverse(Smat), bs=Sbar*bg  
  REAL(KIND = DP) :: init_ethrdg_plrn
  !! Initial coarse threshold to be used if adaptive_ethrdg_plrn=.true. 
  !!
  LOGICAL :: plrn
  !! if .true. calculates perturbated part of the wavefunction
  LOGICAL :: model_vertex_plrn
  !! if .true. calculates vertex model
  LOGICAL :: model_enband_plrn
  !! if .true. calculates vertex model
  LOGICAL :: model_phfreq_plrn
  !! if .true. calculates vertex model
  LOGICAL :: restart_plrn
  !! if .true. Using Ack from written outputs
  LOGICAL :: full_diagon_plrn
  !! if .true. diagonalizing the polaron Hamiltonian with direct diagonalization
  LOGICAL :: cal_psir_plrn
  !! if .true. Generating a 3D-plot for polaron wavefunction
  LOGICAL :: interp_Ank_plrn
  !! if .true. interpolating polaron A(k) from A(Re)
  LOGICAL :: interp_Bqu_plrn
  !! if .true. interpolating polaron phonon component bq from A(Re) and An(k)
  LOGICAL :: debug_plrn
  !! if .true. interpolating polaron phonon component bq from A(Re) and An(k)
  LOGICAL :: time_rev_A_plrn
  !! if .true. time_rev_A_plrn
  LOGICAL :: time_rev_U_plrn
  !! if .true. time_rev_A_plrn
  LOGICAL :: Mmn_plrn
  !! if .true. time_rev_A_plrn
  LOGICAL :: recal_Mmn_plrn
  !! if .true. time_rev_A_plrn
  LOGICAL :: scell_mat_plrn
  !! if .true. activate supercell transformation for polaron
  LOGICAL :: adapt_ethrdg_plrn
  !! if .true. activate adaptive threshold in polaron davidson diagonalization
  ! ----------------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  END MODULE control_epw
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  MODULE klist_epw
  !-----------------------------------------------------------------------
  !!
  !! The variable for the k-points
  !!
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npk
  !
  SAVE
  !
  INTEGER :: kmap(npk)
  !! map of k+q grid into k grid
  INTEGER, ALLOCATABLE :: isk_all(:)
  !! Spin index of each k-point (used in LSDA calculations only)
  INTEGER, ALLOCATABLE :: isk_loc(:)
  !! Spin index of local k-point (used in LSDA calculations only)
  INTEGER, ALLOCATABLE :: isk_dummy(:)
  !! Spin index on the fine grid - dummy at the moment
  REAL(KIND = DP), ALLOCATABLE :: xk_loc(:, :)
  !! List of local (each cores) kpoints in cartesian coordinates
  REAL(KIND = DP), ALLOCATABLE :: xk_all(:, :)
  !! List of all kpoints in cartesian coordinates
  REAL(KIND = DP), ALLOCATABLE :: xk_cryst(:, :)
  !! List of all kpoints in crystal coordinates
  REAL(KIND = DP), ALLOCATABLE :: et_all(:, :)
  !! Eigenenergies on the full coarse k-grid
  REAL(KIND = DP), ALLOCATABLE :: et_loc(:, :)
  !! Eigenenergies on the local (each core) coarse k-grid
  !-----------------------------------------------------------------------
  END MODULE klist_epw
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  MODULE output_epw
  !-----------------------------------------------------------------------
  !!
  !! ... the name of the files
  !!
  SAVE
  !
  CHARACTER(LEN = 80) :: filqf
  !! input  file for the fine q mesh
  CHARACTER(LEN = 80) :: filkf
  !! input  file for the fine k mesh
  CHARACTER(LEN = 80) :: filukk
  !! input  file for the rotation matrix U(k)
  CHARACTER(LEN = 80) :: fildvscf0
  !! output file for the deltavscf used as a fake perturbation to set phases
  CHARACTER(LEN = 80) :: fila2f
  !! input file containing eliashberg spectral function
  CHARACTER(LEN = 80) :: restart_filq
  !! input  file to restart from an exisiting q-file
  !
  !-----------------------------------------------------------------------
  END MODULE output_epw
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  MODULE epwcom
  !-----------------------------------------------------------------------
  USE control_epw
  USE output_epw
  USE klist_epw
  !-----------------------------------------------------------------------
  END MODULE epwcom
  !-----------------------------------------------------------------------
