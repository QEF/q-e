  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from the code PH/phcom - Quantum-ESPRESSO group
  !-----------------------------------------------------------------------
  MODULE input
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
  CHARACTER(LEN = 12) :: elecselfen_type
  !! if 'adiabatic' then compute the adiabatic version of electron self-energy
  !! if 'nonadiabatic' then compute the non-adiabatic version of electron self-energy
  CHARACTER(LEN = 10) :: system_2d
  !! if 'no' then 3D bulk materials [default]
  !! if 'gaussian' then the long-range terms include dipoles only and the range separation
  !! function is approximated by a Gaussian following PRB 94, 085415 (2016)
  !! if 'dipole_sp' then the long-range terms include dipoles only following PRB 107, 155424 (2023)
  !! if 'quadrupole' then the long-range terms include dipoles and quadrupoles terms
  !! following PRL 130, 166301 (2023) and requires the presence of a "quadrupole.fmt" file.
  !! if 'dipole_sh' then the long-range terms include dipoles only following [PRB 105, 115414 (2022)]
  CHARACTER(LEN = 10) :: tc_linear_solver
  !! algorithm to solve T_c eigenvalue problem
  CHARACTER(LEN = 255) :: wannier_plot_list
  !! Field read for parsing Wannier function list
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
  CHARACTER(LEN = 80) :: filirobj
  !! input file containing the objects of IR basis
  CHARACTER(LEN = 80) :: filnscf_coul
  !! input file containing electronic eigenvalues for Coulomb contribution
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
  LOGICAL :: lcumulant
  !! if .TRUE. calculate the electron spectral function using the cumulant expansion method
  LOGICAL :: delta_approx
  !! if .TRUE. the double delta approximation is used for the phonon self energy
  LOGICAL :: ep_coupling
  !! if .TRUE. run e-p coupling calculation
  LOGICAL :: efermi_read
  !! if .TRUE. fermi energy is read from the input file
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
  !! if .TRUE. use analytic continuation to continue imag-axis Eliashberg eqs to real-axis
  LOGICAL :: liso
  !! if .TRUE. solve isotropic case
  LOGICAL :: laniso
  !! if .TRUE. solve anisotropic case
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
  LOGICAL :: fbw
  !! if .TRUE. full-bandwidth calculations will be performed
  LOGICAL :: muchem
  !! if .TURE. the chem. pot. is updated in fbw calculations
  LOGICAL :: positive_matsu
  !! if .TRUE. negative Matsubara frequencies will not be used for sampling.
  LOGICAL :: scattering
  !! if .TRUE. scattering rates are calculated
  LOGICAL :: ltrans_crta
  !! if .TRUE. transport properties are calculated with constant relaxation time approximation (CRTA)
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
  LOGICAL :: longrange_only
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
  LOGICAL :: wannier_plot
  !! if .TRUE. plot Wannier functions
  LOGICAL :: reduce_unk
  !! if .TRUE. plot Wannier functions on reduced grids
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
  LOGICAL :: gb_scattering
  !! If .TRUE. calculate grain boundary scattering rate
  LOGICAL :: gb_only
  !! If .TRUE. calculate only the grain-boundary-limitted mobility
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
  LOGICAL :: cal_acous_plrn
  !! if .true. Generating the acoustic contribution for polaron dtau
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
  LOGICAL :: loptabs
  ! If Quasi absorption or indirect absorption
  LOGICAL :: do_CHBB
  ! Do CHBB calculation within qdabs module
  LOGICAL :: lwfpt
  !! if .TRUE. enable Wannier function perturbation theory calculations
  LOGICAL :: compute_dmat
  !! if .TRUE. compute dmat, the overlap matrix between psi(Sk) and S*psi(k)
  LOGICAL :: calc_nelec_wann
  !! if .TRUE. compute number of electrons in wannierized band
  LOGICAL :: lopt_w2b
  !! if .TRUE. use optimized version of Wannier-to-Bloch Fourier transformation
  LOGICAL :: lfast_kmesh
  !! If .TRUE. use on-the-fly generation of k point mesh within fsthick.
  LOGICAL :: epw_memdist
  !! if .TRUE. distributed storage of epmatwp in MPI processes, only works with etf_mem = 0
  LOGICAL :: dos_tetra
  !! if .true. calculate DOS using tetrahedron method
  LOGICAL :: fd
  !! if .true. ifc file came from finite displacement
  LOGICAL :: a2f_iso
  !! if .TRUE. Eliashberg a2f is calculated without storing g2 matrix elements to files
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
  INTEGER :: gridsamp
  !! Type of the Matsubara freq. sampling (-1= read from file; 0= uniform; 1= sparse; 2= sparse-ir, 3= uniform (FFT))
  INTEGER :: icoulomb
  !! this tag controls how to calculate the coulomb kernel in the Eliashberg subroutines.
  INTEGER :: wannier_plot_supercell(3)
  !! Size of supercell for plotting Wannier functions
  INTEGER :: mob_nfreq
  !! Number of frequency for the spectral decomposition of mobility
  INTEGER :: neta = 9
  !! Number of broadening parameters
  INTEGER :: nomega
  !! Number of frequency (photon energy) points
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
  REAL(KIND = DP) :: g_scale_plrn
  !! scaling factor of g matrix
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
  INTEGER :: david_ndim_plrn
  !! Dimension of workspace for Davidson diagonalization
  INTEGER :: len_mesh
  ! The number of meshgrids for QD
  INTEGER :: meshnum
  ! Last meshgrid for QD calculation (restart)
  INTEGER :: start_mesh
  ! Starting mesh point for restart calculations
  INTEGER :: DW
  ! Second order term of QDPT
  INTEGER :: nq_init
  ! Temperature integration method
  INTEGER :: wf_quasi
  ! Write Eigenvectors for QD
  INTEGER :: mode_res
  ! Number of modes to resolve
  INTEGER :: ahc_nbnd
  !! Number of bands included in ph.x Allen-Heine-Cardona calculation
  INTEGER :: ahc_nbndskip
  !! Number of low-lying bands exclude in ph.x Allen-Heine-Cardona calculation
  INTEGER :: kmap(npk)
  !! map of k+q grid into k grid
  INTEGER, ALLOCATABLE :: isk_all(:)
  !! Spin index of each k-point (used in LSDA calculations only)
  INTEGER, ALLOCATABLE :: isk_loc(:)
  !! Spin index of local k-point (used in LSDA calculations only)
  INTEGER, ALLOCATABLE :: isk_dummy(:)
  !! Spin index on the fine grid - dummy at the moment
  !
  REAL(KIND = DP) :: tempsmin
  !! min. temperature in Eliashberg equations - deprecated as an input parameter
  REAL(KIND = DP) :: tempsmax
  !! max. temperature - deprecated as an input parameter
  REAL(KIND = DP) :: temps(50)
  !! input temperature array (units of Kelvin)
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
  REAL(KIND = DP) :: eps_acoustic
  !! min. phonon frequency for e-p and a2f calculations
  REAL(KIND = DP) :: acoustic_plrn
  !! max. phonon frequency for polaron dtau
  REAL(KIND = DP) :: dtau_max_plrn
  !! max. polaron displacements in a physical solution
  REAL(KIND = DP) :: degaussq
  !! smearing for sum over q in e-ph coupling
  REAL(KIND = DP) :: delta_qsmear
  !! change in energy for each additional smearing in the a2f
  REAL(KIND = DP) :: muc
  !! effective Coulomb potential in Eliashberg equations
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
  REAL(KIND = DP) :: dos_del
  !! Delta_E in electronic dos for Fermi window (in eV)
  REAL(KIND = DP) :: griddens
  !! Measure of sparsity of the grid
  REAL(KIND = DP) :: scissor
  !! Value of the scissor shift in eV.
  REAL(KIND = DP) :: ncarrier
  !! Amount of carrier concentration in cm^-3 (cm^-2 in 2D) when doping a semiconductors
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
  REAL(KIND = DP) :: sr_crta
  !! Input scattering rate (in THz) for constant relaxation time approximation (CRTA)
  REAL(KIND = DP) :: ii_charge
  !! charge of the ionized impurities, units of electron charge for input
  REAL(KIND = DP) :: ii_n
  !! density of impurities, input in units of cm^-3 (cm^-2 in 2D)
  REAL(KIND = DP) :: ii_eda
  !! ionization energy of the donor or accpetor impurity, for partial ionization calcs
  REAL(KIND = DP) :: gb_size
  !! average size of grain in polycrystallite in units of nm
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
  REAL(KIND = DP) :: omegamin
  !! Photon energy minimum (in eV)
  REAL(KIND = DP) :: omegamax
  !! Photon energy maximum (in eV)
  REAL(KIND = DP) :: omegastep
  !! Photon energy step (in eV)
  REAL(KIND = DP) :: sigma_ref
  !! Reference conductivity for Drude dielectric function
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
  !! FIXME
  REAL(KIND = DP) :: as(3,3)
  !! Supercell lattice vectors transformed by Smat, as=S*at
  REAL(KIND = DP) :: bs(3,3)
  !! Supercell reciprocal lattice vectors transformed by Sbar=transpose(inverse(Smat), bs=Sbar*bg
  REAL(KIND = DP) :: init_ethrdg_plrn
  !! Initial coarse threshold to be used if adaptive_ethrdg_plrn=.true.
  REAL(KIND=DP) :: QD_min
  ! scaling factor for g-matrix scaling
  REAL(KIND=DP) :: QD_bin
  ! scaling factor for g-matrix scaling
  REAL(KIND = DP) :: ahc_win_min
  !! Lower bound of AHC window for the lower Fan term.
  REAL(KIND = DP) :: ahc_win_max
  !! Upper bound of AHC window for the lower Fan term.
  REAL(KIND = DP) :: a_gap0
  !! This determines the shape of initial guess of gap function
  REAL(KIND = DP) :: emax_coulomb, emin_coulomb
  !! Electron Energy range for Coulomb kernel in Eliashberg calculation (used only when icoulomb > 0)
  REAL(KIND = DP) :: eps_cut_ir
  !! This is a threshold to ignore negligibly small IR coefficients.
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
  ! ZD: For ex-plrn calculation
  LOGICAL :: exciton
  !! Calculate exciton-phonon coupling and optical transition
  LOGICAL :: explrn
  !! if .true. activates the calculation of excitonic polaron. Must go with epwread=true
  LOGICAL :: plot_explrn_e
  !! if .true. plot the e density of excitonic polaron. Will skip the explrn calculations
  LOGICAL :: plot_explrn_h
  !! if .true. plot the h density of excitonic polaron. Will skip the explrn calculations
  LOGICAL :: only_c_explrn
  !! if .true. only use conduction band to construct ex-ph matrix
  LOGICAL :: only_v_explrn
  !! if .true. only use valence band to construct ex-ph matrix
  LOGICAL :: only_pos_modes_explrn
  !! Whether include positive-frequency modes only 
  INTEGER :: nbndc_explrn
  !! Number of conduction bands for exciton calculation
  INTEGER :: nbndv_explrn
  !! Number of valence bands for exciton calculation
  INTEGER :: negnv_explrn
  ! Number of exciton states to be included in the ex-polaron calculation
  integer :: step_k1_explrn, step_k2_explrn, step_k3_explrn
  !! used to model smaller excitonic polaron
  ! ZD 
  !! ymPan: Time dependent Boltzmann Equation
  LOGICAL :: do_tdbe
  !! if .true. start time propagating BTE
  REAL(KIND = DP):: dt_tdbe       
  !! time step in fs 
  INTEGER :: nt_tdbe
  !! nb of time step, should be defined in epwcom
  INTEGER :: twrite_tdbe 
  REAL(KIND = DP) :: temp_el_tdbe
  !! initial electronic temperature 
  REAL(KIND = DP) :: temp_ph_tdbe
  !! initial phononic temperature
  CHARACTER(LEN = 75) :: init_type_tdbe
  !! type of initial carrier distribution
  REAL(KIND = DP) :: init_sigma_tdbe
  !! gaussian smearing parameter for initial carrier
  !! distribution
  CHARACTER(LEN = 75) :: solver_tdbe
  !! time propagation methods
  REAL(KIND = DP) :: ef_c_tdbe
  !! initial energy of hot electron or chemical potential
  !! of electron
  REAL(KIND = DP) :: ef_v_tdbe
  !! initial energy of hot hole or chemical potential 
  !! of hole
  INTEGER :: carr_dyn_tdbe
  !! if = 1, real time dynamic of electron
  !! if = 2, real time dynamic of hole
  !! if =3 , real time dynamic for both 
  LOGICAL :: dg_tdbe
  !! use two different grid for energy and e-ph matrix
  INTEGER :: nqf1d, nqf2d, nqf3d
  !! double grid for phonon q points
  INTEGER :: nkf1d, nkf2d, nkf3d
  !! double grid for electron k points
  INTEGER :: dph_tdbe
  !! the factor between time step for e-ph interaction and that
  !! for ph-ph interaction
  LOGICAL :: restart_tdbe
  !! restart TDBE from an previous interrupted calculation 
  LOGICAL :: phph_tdbe
  !! include phonon-phonon scattering in TDBE simulations
  REAL(KIND=DP) :: phwmin_tdbe
  !! min. phonon frequency for ph-ph calculations
  CHARACTER(LEN = 100) :: ephmat_dir ='./'
  !! path to the e-ph matrix elements prefix.ephmat
  LOGICAL :: lscreen_tdbe
  !! Apply screening to e-ph matrix elements in TDBE simulations
  !
  !-----------------------------------------------------------------------
  END MODULE input
  ! -----------------------------------------------------------------------
