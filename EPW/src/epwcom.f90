  !                                                                            
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
  !
  ! Common variables for the epw program
  !  
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY :  DP
  USE parameters, ONLY: npk
  !
  ! ... the variable controlling the phonon run
  !
  SAVE
  ! 
  INTEGER :: ngaussw, nw, selfen_type, nbndsub, nbndskip, num_iter, iprint, &
       nsmear, rand_nq, rand_nk, nqf1, nqf2, nqf3, nkf1, nkf2, nkf3, &
       nqsmear, nqstep, neptemp, nswfc, nswc, nswi, nstemp, nsiter, broyden_ndim, &
       nw_specfun
  ! ngaussw: smearing type for Fermi surface average in e-ph coupling after wann. interp.
  ! nw: nr. of bins for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  ! selfen_type: choice of real/imag part of phonon self-energy when epstrict = .true.
  ! nbndsub: nr. of bands in the optimal subspace (when disentanglement is used)
  ! nbndskip:  nr. of bands to be skipped when we use only a subspace (this is nfirstwin-1 in Marzari's notation)
  ! @JN
  ! num_iter: nr. of steps used in creating the maximally localized Wannier functions
  ! iprint: verbosity of the wannier90 code
  ! nsmear: nr. of temperature-like smearings calculated in selfen_phon, indabs
  ! rand_nq: use random points for the fine q-mesh
  ! rand_nk: use random points for the fine k-mesh
  ! nqf1, nqf2, nqf3: qx,qy,qz sizes of the uniform phonon fine mesh to be used
  ! nkf1,nkf2,nkf3: kx,ky,kz sizes of the uniform electron fine mesh to be used
  ! nqsmear: nr. of smearings used to calculate a2f 
  ! nqstep: nr. of steps used to calculate a2f 
  ! neptemp: nr. temperature values for the electronic Fermi occupations in the e-p calculation (not used for now)
  ! @RM
  ! nswfc    : nr. of grid points between (0,wsfc)
  ! nswc     : nr. of grid points between (wsfc,wscut)
  ! nswi     : nr. of grid points for Eliashberg equations of imaginary axis
  ! nstemp   : nr. of temperature points for Eliashberg equations
  ! nsiter   : nr. of iterations for self-consistency
  ! broyden_ndim : nr. of iterations used in broyden mixing scheme
  ! nw_specfun   : nr. of bins for frequency in electron spectral function due to e-p interaction 
  !
  INTEGER :: nwanxx = 200, ntempxx = 25
  ! parameter used in writing prefix.win file.  Maximum number of wannier functions
  REAL (KIND=DP) :: degaussw, fsthick, wmin, wmax, ecutse, dis_win_min, dis_win_max, &
       dis_froz_min, dis_froz_max, delta_smear, eminabs, emaxabs, deltaeabs, eps_acustic, &
       degaussq, delta_qsmear, muc, wsfc, pwc, wscut, tempsmin, tempsmax, & 
       broyden_beta, conv_thr_raxis, conv_thr_iaxis, conv_thr_racon, gap_edge, & 
       max_memlt, fermi_energy, wmin_specfun, wmax_specfun 
  ! degaussw: smearing width for Fermi surface average in e-ph coupling after wann interp
  ! fsthick: thickness of the Fermi shell for averaging the e-ph matrix element
  ! temperature for the electronic Fermi occupations in the e-p calculation 
  ! wmin/wmax: min/max frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  ! dis_win_min: min energy of the Wannier disentanglement window
  ! dis_win_max: max energy of the Wannier disentanglement window
  ! dis_froz_min: min energy of the frozen Wannier disentanglement window
  ! dis_froz_max: max energy of the frozen Wannier disentanglement window
  ! delta_smear: change in energy for each additional smearing in the selfen_phon and indabs
  ! @RM
  ! eps_acustic: min. phonon frequency for e-p and a2f calculations
  ! degaussq: smearing for sum over q in e-ph coupling  
  ! delta_qsmear: change in energy for each additional smearing in the a2f
  ! muc     : effective Coulomb potential in Eliashberg equations
  ! wsfc    : intermediate freqeuncy between (0,wscut)
  ! pwc     : power used to define a non-uniform grid between wsfc and wscut
  ! wscut   : upper limit cutoff frequency in Eliashberg equations (at least 5 times wsphmax)
  ! tempsmin : min. temperature in Eliashberg equations
  ! tempsmax : max. temperature
  ! broyden_beta : mixing factor for broyden mixing
  ! conv_thr_raxis : convergence threshold for iterative solution of real-axis Eliashberg equations
  ! conv_thr_iaxis : convergence threshold for iterative solution of imag-axis Eliashberg equations
  ! conv_thr_racon : convergence threshold for iterative solution of analytic continuation of
  !                  Eliashberg equations from imag- to real-axis
  ! gap_edge : initial guess of the superconducting gap
  ! max_memlt : maximum memory that can be allocated per pool
  ! fermi_energy : fermi energy is given in the input file 
  ! wmin_specfun : min frequency in electron spectral function due to e-p interaction 
  ! wmax_specfun : max frequency in electron spectral function due to e-p `interaction
  !
  logical :: phinterp , elinterp, tphases, epstrict, tshuffle2, elecselfen, phonselfen, &
       epbread, epbwrite, epwread, epwwrite, specfun, wannierize, spinors, parallel_k,      &
       parallel_q, a2f, etf_mem, write_wfn, kmaps, nest_fn, rand_q, rand_k, &
       mp_mesh_q, mp_mesh_k, indabs, eig_read, wepexst, epexst, vme, twophoton, & 
       ephwrite, band_plot, lreal, limag, lpade, lacon, liso, laniso, lpolar, lunif, & 
       kerwrite, kerread, imag_read, eliashberg, ep_coupling, & 
       efermi_read, system_2d, delta_approx
  ! phinterp: if .TRUE. perform phonon interpolation of e-p matrix
  ! elinterp: if .TRUE. perform electron interpolation of e-p matrix
  ! tphases:  if .TRUE. set absolute reference for unitary gauge of the eigenvectors
  ! epstrict: if .TRUE. use strict selection rule in phonon linewidth calculation
  ! tshuffle2: shuffle mode for electrons + load all phonons at once
  ! elecselfen: if .TRUE. calculate electron selfenergy due to e-p interaction
  ! phonselfen: if .TRUE. calculate phonon selfenergy due to e-p interaction
  ! epbread: if .TRUE. read epmatq from files .epb
  ! epbwrite: if .TRUE. write epmatq to files .epb
  ! epwread: if .TRUE. read all quantities in Wannier representation from file epwdata.fmt
  ! epwwrite: if .TRUE. write all quantities in Wannier representation to file epwdata.fmt
  ! specfun: if .TRUE. calculate spectral electron function due to e-p interaction
  ! wannierize: if .TRUE. run the wannier90 code
  ! spinors: if .TRUE. run the wannier90 with spinors
  ! parallel_k: if .TRUE. scatter the electron k-points on the fine mesh among pools (not q)
  ! parallel_q: if .TRUE. scatter the phonon q-points on the fine mesh among pools (not k)
  ! a2f: if .TRUE. calculate Eliashberg spectral electron function from selfen_phon
  ! etf_mem: If .true., the fine Bloch-space e-ph matrix elements are stored in memory
  ! write_wfn: if .TRUE. write out UNK files in wannier90
  ! kmaps: if .TRUE. read kmap and kgmap from disk.  Do not calculate
  ! nest_fn: if .TRUE. calculate the electronic nesting function (metals only)
  ! rand_q: if .TRUE. use random points for the fine q-mesh
  ! rand_k: if .TRUE. use random points for the fine k-mesh
  ! mp_mesh_q: if .TRUE. use points in the irreducible wedge for the uniform fine q-mesh
  ! mp_mesh_k: if .TRUE. use points in the irreducible wedge for the uniform fine k-mesh
  ! indabs: if .TRUE. calculated indirect absorption spectrum
  ! eig_read: if .true. then readin a set of electronic eigenvalues in eV to replace the calcualted ones
  ! wepexst: if .TRUE. prefix.epmatwe files are already on disk.  don't recalculate. debugging param
  ! epexst: if .TRUE. prefix.epmatwp files are already on disk.  don't recalculate  debugging param
  ! vme: if .TRUE. calculate velocity matrix elements
  ! twophoton: 
  ! @RM
  ! ephwrite : if .true. write el-ph matrix elements on the fine mesh to file
  ! band_plot : if .true. write filrs to plot band structure and phonon dispersion
  ! lreal    : if .true. solve real-axis Eliashberg eqautions
  ! limag    : if .true. solve imag-axis Eliashberg eqautions
  ! lpade    : if .true. use pade approximants to continue imag-axis 
  !            Eliashberg equtions to real-axis
  ! lacon    : if .true. use analytic continuation to continue imag-axis 
  !            Eliashberg equtions to real-axis
  ! liso     : if .true. solve isotropic case
  ! laniso   : if .true. solve anisotropic case
  ! lpolar   : if .true. enable the correct Wannier interpolation in the case of polar material.  
  ! lunif    : if .true. a uniform grid is defined between wsfc and wc for real-axis calculations
  ! kerwrite : if .true. write Kp and Km to files .ker for real-axis calculations
  ! kerread  : if .true. read Kp and Km from files .ker for real-axis calculations
  ! imag_read : if .true. read from file Delta and Znorm on the imaginary-axis
  ! eliashberg : if .true. solve the Eliashberg equations 
  ! ep_coupling : if .true. run e-p coupling calculation
  ! efermi_read : if .true. fermi energy is read from the input file
  ! system_2d : if .true. the system is 2 dimensional (vaccum is in z-direction)
  ! delta_approx : if .true. the double delta approximation is used for the phonon self energy
  CHARACTER(len=100) :: dvscf_dir ='./' ! directory for .dvscf and .dyn files (wannier interpolation)
  CHARACTER(len=80) :: filelph, fileig ! output file for the electron-phonon coefficients
  CHARACTER(len=256), dimension(200) :: proj, wdata ! projections and any extra info for W90 
  CHARACTER(LEN=75) :: title ! ...  title of the simulation  
  REAL (kind=DP), dimension(25) :: eptemp 
  REAL (kind=DP), dimension(50) :: temps ! temperature entering in the Eliashberg equtions (units of Kelvin)
  integer :: iswitch

END MODULE control_epw
!
!
MODULE klist_epw
! The variable for the k-points 
USE kinds, ONLY: DP
USE parameters, ONLY :npk
!
SAVE
!
INTEGER :: kmap(npk)  ! map of k+q grid into k grid (only when tshuffle=.true.)
LOGICAL :: tshuffle   ! if .TRUE. refold k+q grid into k grid
REAL(DP) :: xk_cryst(3,npk) ! List of all kpoints in crystal coordinates

END MODULE klist_epw

!
MODULE units_epw
  !
  ! ... the units of the files and the record lengths
  !
  SAVE
  !
  INTEGER :: lrcuf, lrepmatf, lretf
  !INTEGER :: iudvscf0
  !
  ! iudvscf0: the unit where the delta Vscf is read to generate the fake perturbation 
  ! lrcuf: the length of the record for the rotation matrix
  ! lrepmatf: the length of the record for the electron-phonon matrix
  ! lretf: the length of the record for the interpolated hamiltonian eigenvalues
  !
  logical, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_epw
!
!
MODULE output_epw
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=80) :: filqf, filkf, filukk, filukq, fildvscf0, &
                        fila2f
  !
  ! filqf: input  file for the fine q mesh
  ! filkq: input  file for the fine k mesh
  ! filukk: input  file for the rotation matrix U(k)
  ! filukq: input  file for the rotation matrix U(k+q)
  ! fildvscf0: output file for the deltavscf used as a fake perturbation to set phases
  ! @RM
  ! fila2f  : input file containing eliashberg spectral function
  !
END MODULE output_epw
!
MODULE epwcom
  USE control_epw
  USE units_epw
  USE output_epw
  USE klist_epw
  END MODULE epwcom
