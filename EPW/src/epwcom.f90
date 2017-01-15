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
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npk
  !!
  !! ... the variable controlling the EPW run
  !!
  SAVE
  ! 
  INTEGER :: ngaussw
  !! smearing type for Fermi surface average in e-ph coupling after wann. interp.
  INTEGER :: nw
  !! nr. of bins for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  INTEGER :: nbndsub
  !! nr. of bands in the optimal subspace (when disentanglement is used)
  INTEGER :: nbndskip
  !! nr. of bands to be skipped when we use only a subspace (this is nfirstwin-1 in Marzari's notation)
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
  INTEGER :: ntempxx = 25
  !! Maximum number of wannier functions
  !
  ! Superconductivity
  INTEGER :: nswfc
  !! nr. of grid points between (0,wsfc)
  INTEGER :: nswc
  !! nr. of grid points between (wsfc,wscut)
  INTEGER :: nswi
  !! nr. of grid points for Eliashberg equations of imaginary axis
  INTEGER :: nstemp
  !! nr. of temperature points for Eliashberg equations
  INTEGER :: nsiter
  !! nr. of iterations for self-consistency
  INTEGER :: broyden_ndim
  !! nr. of iterations used in broyden mixing scheme
  INTEGER :: nw_specfun
  !! nr. of bins for frequency in electron spectral function due to e-p interaction 
  INTEGER :: restart_freq
  !! Create a restart point during the interpolation part every restart_freq q/k-points. 
  !
  REAL (KIND=DP) :: degaussw
  !! smearing width for Fermi surface average in e-ph coupling after wann interp
  REAL (KIND=DP) :: fsthick
  !! thickness of the Fermi shell for averaging the e-ph matrix element
  REAL (KIND=DP) :: eptemp
  ! temperature for the electronic Fermi occupations in the e-p calculation 
  REAL (KIND=DP) :: wmin
  !! min frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  REAL (KIND=DP) :: wmax
  !! max frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  REAL (KIND=DP) :: dis_win_min
  !! min energy of the Wannier disentanglement window
  REAL (KIND=DP) :: dis_win_max
  !! max energy of the Wannier disentanglement window
  REAL (KIND=DP) :: dis_froz_min
  !! min energy of the frozen Wannier disentanglement window
  REAL (KIND=DP) :: dis_froz_max
  !! max energy of the frozen Wannier disentanglement window
  REAL (KIND=DP) :: delta_smear
  !! change in energy for each additional smearing in the selfen_phon
  ! 
  ! Superconductivity
  REAL (KIND=DP) :: eps_acustic
  !! min. phonon frequency for e-p and a2f calculations
  REAL (KIND=DP) :: degaussq
  !! smearing for sum over q in e-ph coupling  
  REAL (KIND=DP) :: delta_qsmear
  !! change in energy for each additional smearing in the a2f
  REAL (KIND=DP) :: muc
  !! effective Coulomb potential in Eliashberg equations
  REAL (KIND=DP) :: wsfc
  !! intermediate freqeuncy between (0,wscut)
  REAL (KIND=DP) :: pwc
  !! power used to define a non-uniform grid between wsfc and wscut
  REAL (KIND=DP) :: wscut
  !! upper limit cutoff frequency in Eliashberg equations (at least 5 times wsphmax)
  REAL (KIND=DP) :: tempsmin
  !! min. temperature in Eliashberg equations
  REAL (KIND=DP) :: tempsmax
  !! max. temperature
  REAL (KIND=DP) :: broyden_beta
  !! mixing factor for broyden mixing
  REAL (KIND=DP) :: conv_thr_raxis
  !! convergence threshold for iterative solution of real-axis Eliashberg equations
  REAL (KIND=DP) :: conv_thr_iaxis
  !! convergence threshold for iterative solution of imag-axis Eliashberg equations
  REAL (KIND=DP) :: conv_thr_racon
  !! convergence threshold for iterative solution of analytic continuation of
  !! Eliashberg equations from imag- to real-axis
  REAL (KIND=DP) :: gap_edge
  !! initial guess of the superconducting gap
  REAL (KIND=DP) :: max_memlt
  !! maximum memory that can be allocated per pool
  REAL (KIND=DP) :: fermi_energy
  !! fermi energy is given in the input file 
  REAL (KIND=DP) :: wmin_specfun
  !! min frequency in electron spectral function due to e-p interaction 
  REAL (KIND=DP) :: wmax_specfun
  !! max frequency in electron spectral function due to e-p `interaction
  REAL (kind=DP), dimension(50) :: temps 
  !! temperature entering in the Eliashberg equtions (units of Kelvin)
  !
  ! Conductivity
  REAL (KIND=DP) :: scissor
  !! Value of the scissor shift in eV.
  REAL (KIND=DP) :: ncarrier
  !! Amount of carrier concentration in cm^-3 when doping a semiconductors
  !
  !LOGICAL :: tphases
  !! tphases:  if .TRUE. set absolute reference for unitary gauge of the eigenvectors
  LOGICAL :: elecselfen
  !! if .TRUE. calculate electron selfenergy due to e-p interaction
  LOGICAL :: phonselfen
  !! if .TRUE. calculate phonon selfenergy due to e-p interaction
  LOGICAL :: epbread
  !! if .TRUE. read epmatq from files .epb
  LOGICAL :: epbwrite
  !! if .TRUE. write epmatq to files .epb
  LOGICAL :: epwread
  !! if .TRUE. read all quantities in Wannier representation from file epwdata.fmt
  LOGICAL :: epwwrite
  !! if .TRUE. write all quantities in Wannier representation to file epwdata.fmt
  LOGICAL :: restart
  !! if .TRUE. restart a calculation stopped during the interpolation phase from reading 
  !! the XXX.restart file. 
  LOGICAL :: specfun
  !! if .TRUE. calculate spectral electron function due to e-p interaction
  LOGICAL :: wannierize
  !! if .TRUE. run the wannier90 code
  LOGICAL :: parallel_k
  !! if .TRUE. scatter the electron k-points on the fine mesh among pools (not q)
  LOGICAL :: parallel_q
  !! if .TRUE. scatter the phonon q-points on the fine mesh among pools (not k)
  LOGICAL :: a2f
  !! if .TRUE. calculate Eliashberg spectral electron function from selfen_phon
  LOGICAL :: etf_mem
  !! If .true., the fine Bloch-space e-ph matrix elements are stored in memory
  LOGICAL :: write_wfn
  !! if .TRUE. write out UNK files in wannier90
  LOGICAL :: kmaps
  !! if .TRUE. read kmap and kgmap from disk.  Do not calculate
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
  !! if .true. then readin a set of electronic eigenvalues in eV to replace the calcualted ones
  LOGICAL :: wepexst
  !! if .TRUE. prefix.epmatwe files are already on disk.  don't recalculate. debugging param
  LOGICAL :: epexst
  !! if .TRUE. prefix.epmatwp files are already on disk.  don't recalculate  debugging param
  LOGICAL :: vme
  !! if .TRUE. calculate velocity matrix elements
  LOGICAL :: band_plot
  ! band_plot : if .true. write filrs to plot band structure and phonon dispersion
  LOGICAL :: lpolar 
  !! if .true. enable the correct Wannier interpolation in the case of polar material.  
  LOGICAL :: lifc
  !! if .true. reads interatomic force constants produced by q2r.x for phonon interpolation
  LOGICAL :: delta_approx
  !! if .true. the double delta approximation is used for the phonon self energy
  LOGICAL :: ep_coupling
  !! if .true. run e-p coupling calculation
  LOGICAL :: efermi_read
  !! if .true. fermi energy is read from the input file
  LOGICAL :: system_2d
  !! if .true. the system is 2 dimensional (vaccum is in z-direction)
  !
  ! Superconductivity
  LOGICAL :: ephwrite
  !! if .true. write el-ph matrix elements on the fine mesh to file
  LOGICAL :: lreal
  !! if .true. solve real-axis Eliashberg eqautions
  LOGICAL :: limag
  !! if .true. solve imag-axis Eliashberg eqautions
  LOGICAL :: lpade
  !! if .true. use pade approximants to continue imag-axis Eliashberg equtions to real-axis
  LOGICAL :: lacon
  !! if .true. use analytic continuation to continue imag-axis Eliashberg equtions to real-axis
  LOGICAL :: liso
  !! if .true. solve isotropic case
  LOGICAL :: laniso
  !! if .true. solve anisotropic case
  LOGICAL :: lunif
  !! if .true. a uniform grid is defined between wsfc and wc for real-axis calculations
  LOGICAL :: kerwrite
  !! if .true. write Kp and Km to files .ker for real-axis calculations
  LOGICAL :: kerread
  !! if .true. read Kp and Km from files .ker for real-axis calculations
  LOGICAL :: imag_read
  !! if .true. read from file Delta and Znorm on the imaginary-axis
  LOGICAL :: eliashberg
  !! if .true. solve the Eliashberg equations 
  !
  ! Conductivity
  LOGICAL :: scattering
  !! if .true. scattering rates are calculated
  LOGICAL :: scattering_serta
  !! if .true. scattering rates are calculated using self-energy relaxation-time-approx
  LOGICAL :: scattering_0rta
  !! if .true. scattering rates are calculated using 0th order relaxation-time-approx
  LOGICAL :: int_mob
  !! if .true. computes the intrinsic mobilities. This means that the electron and hole carrier density is equal.
  LOGICAL :: iterative_bte
  !! if .true. the iterative solution for BTE is compute. A first run with scattering_serta = .true. is required. 
  LOGICAL :: carrier
  !! if .true. computes the doped electronic mobilities.
  LOGICAL :: longrange
  !! if .true. computes the long range interaction of el-ph. Can only be .true. if lpolar is also true.
  LOGICAL :: shortrange
  !! if .true. computes the long range interaction of el-ph. Can only be .true. if lpolar is also true.  
  !
  CHARACTER(len=100) :: dvscf_dir ='./'
  !! directory for .dvscf and .dyn files (wannier interpolation)
  CHARACTER(len=80) :: fileig 
  !! output file for the electron-phonon coefficients
  CHARACTER(len=256), dimension(200) :: proj, wdata 
  !! projections and any extra info for W90 
  CHARACTER(LEN=75) :: title 
  !! ...  title of the simulation  
  CHARACTER(LEN=10)  :: asr_typ
  !! type of ASR if lifc=.true.
  !
END MODULE control_epw
!
!
MODULE klist_epw
  !!        
  !! The variable for the k-points 
  !! 
  USE kinds, ONLY: DP
  USE parameters, ONLY :npk
  !
  SAVE
  !
  INTEGER :: kmap(npk)  ! map of k+q grid into k grid 
  REAL(DP) :: xk_cryst(3,npk) ! List of all kpoints in crystal coordinates
  ! 
END MODULE klist_epw
!
MODULE units_epw
  !!
  !! ... the units of the files and the record lengths
  !!
  SAVE
  !
  INTEGER :: lrepmatf
  !INTEGER :: iudvscf0
  !
  ! iudvscf0: the unit where the delta Vscf is read to generate the fake perturbation 
  ! lrepmatf: the length of the record for the electron-phonon matrix
  !
  logical, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_epw
!
!
MODULE output_epw
  !!
  !! ... the name of the files
  !!
  SAVE
  !
  CHARACTER (LEN=80) :: filqf
  !! input  file for the fine q mesh
  CHARACTER (LEN=80) :: filkf
  !! input  file for the fine k mesh
  CHARACTER (LEN=80) :: filukk
  !! input  file for the rotation matrix U(k)
  CHARACTER (LEN=80) :: filukq
  !! input  file for the rotation matrix U(k+q)
  CHARACTER (LEN=80) :: fildvscf0
  !! output file for the deltavscf used as a fake perturbation to set phases
  CHARACTER (LEN=80) :: fila2f
  !! input file containing eliashberg spectral function
  !
END MODULE output_epw
!
MODULE epwcom
  USE control_epw
  USE units_epw
  USE output_epw
  USE klist_epw
END MODULE epwcom
