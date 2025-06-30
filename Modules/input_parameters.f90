!
! Copyright (C) 2002-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------
! TB
! included gate related stuff, search 'TB'
!---------------------------------------------
!
!=----------------------------------------------------------------------------=!
MODULE input_parameters
!=----------------------------------------------------------------------------=!
!!  This module contains:
!
!!  * the definitions of all input parameters (both those read from namelists 
!!    and those read from cards);
!!  * the definitions of all namelists;
!!  * routines that allocate data needed in input.
!
!!  Note that all values are initialized, but the default values should be
!!  set in the appropriate routines contained in module \(\texttt{read_namelists}\).  
!!  The documentation of input variables can be found in Doc/INPUT_PW.*
!!  (for pw.x) or in Doc/INPUT_CP (for cp.x).
!
!!  Originally written by Carlo Cavazzoni for FPMD.
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nsx, natx, sc_size, nsolx
  USE wannier_new,ONLY : wannier_data
  USE upf_params, ONLY : lqmax
  !
  IMPLICIT NONE
  !
  SAVE
  !
!=----------------------------------------------------------------------------=!
! BEGIN manual
!
!
! * DESCRIPTION OF THE INPUT FILE
!  (to be given as standard input)
!
!  The input file has the following layout:
!
!     &CONTROL
!       control_parameter_1,
!       control_parameter_2,
!       .......
!       control_parameter_Lastone
!     /
!     &SYSTEM
!       sistem_parameter_1,
!       sistem_parameter_2,
!       .......
!       sistem_parameter_Lastone
!     /
!     &ELECTRONS
!       electrons_parameter_1,
!       electrons_parameter_2,
!       .......
!       electrons_parameter_Lastone
!     /
!     &IONS
!       ions_parameter_1,
!       ions_parameter_2,
!       .......
!       ions_parameter_Lastone
!     /
!     &CELL
!       cell_parameter_1,
!       cell_parameter_2,
!       .......
!       cell_parameter_Lastone
!     /
!     ATOMIC_SPECIES
!      slabel_1 mass_1 pseudo_file_1
!      slabel_2 mass_2 pseudo_file_2
!      .....
!     ATOMIC_POSITIONS
!      alabel_1  px_1 py_1 pz_1
!      alabel_2  px_2 py_2 pz_2
!      .....
!     CARD_3
!     ....
!     CARD_N
!
!  -- end of input file --
!
!=----------------------------------------------------------------------------=!
!  CONTROL Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        CHARACTER(len=80) :: title = ' '
        !! a string describing the current job

        CHARACTER(len=80) :: calculation = 'none'
        !! specify the type of the simulation

        CHARACTER(len=80) :: calculation_allowed(15)
        !! allowed calculations (names)
        DATA calculation_allowed / 'scf', 'nscf', 'relax', 'md', 'cp', &
          'vc-relax', 'vc-md', 'vc-cp', 'bands', 'neb', 'smd', 'cp-wf', &
          'vc-cp-wf', 'cp-wf-nscf', 'ensemble'/
        CHARACTER(len=80) :: verbosity = 'default'
        !! define the verbosity of the code output
        CHARACTER(len=80) :: verbosity_allowed(6)
        !! allowed verbosity options
        DATA verbosity_allowed / 'debug', 'high', 'medium', 'default', &
                                 'low', 'minimal' /

        CHARACTER(len=80) :: restart_mode = 'restart'
        !! specify how to start/restart the simulation
        CHARACTER(len=80) :: restart_mode_allowed(3)
        !! allowed restart options
        DATA restart_mode_allowed / 'from_scratch', 'restart', 'reset_counters' /

        INTEGER :: nstep = 10
        !! number of simulation steps, see "restart_mode"

        INTEGER :: iprint = 10
        !! number of steps/scf iterations between successive writings
        !! of relevant physical quantities to standard output
        INTEGER :: max_xml_steps = 0 
        !! max number of steps between successive appending of an xml step 
        !! in the xml data file, default 0 means all steps are printed
        INTEGER :: isave = 100
        !! number of steps between successive savings of
        !! information needed to restart the run (see "ndr", "ndw").
        !! Used only in CP

        LOGICAL :: tstress = .true.
        !! if TRUE calculate the stress tensor

        LOGICAL :: tprnfor = .true.
        !! if TRUE calculate the atomic forces

        REAL(DP) :: dt = 1.0_DP
        !! time step for molecular dynamics simulation, in atomic units.  
        !! CP: 1 a.u. of time = 2.4189 * 10^-17 s, PW: twice that much.  
        !! Typical values for CP simulations are between 1 and 10 a.u.  
        !! For Born-Oppenheimer simulations, larger values can be used,
        !! since it mostly depends only upon the mass of ions.

        INTEGER :: ndr = 50
        !! Fortran unit from which the code reads the restart file

        INTEGER :: ndw = 50
        !! Fortran unit to which the code writes the restart file

        CHARACTER(len=256) :: outdir = './'
        !! specify the directory where the code opens output and restart
        !! files. When possible put this directory in the fastest available
        !! filesystem ( not NFS! )

        CHARACTER(len=256) :: prefix = 'prefix'
        !! specify the prefix for the output file, if not specified the
        !! files are opened as standard fortran units.

        CHARACTER(len=256) :: pseudo_dir = './'
        !! specify the directory containing the pseudopotentials

        REAL(DP) :: refg = 0.05_DP
        !! Accurancy of the interpolation table, interval between
        !! table values in Rydberg

        CHARACTER(len=256) :: wfcdir = 'undefined'
        !! scratch directory that is hopefully local to the node
        !! to store large, usually temporary files.

        REAL(DP) :: max_seconds = 1.0E+7_DP
        !! smoothly terminate program after the specified number of seconds
        !! this parameter is typically used to prevent an hard kill from
        !! the queuing system.

        REAL(DP) :: ekin_conv_thr = 1.0E-5_DP
        !! convergence criterion for electron minimization. 
        !! This criterion is met when "ekin < ekin_conv_thr",
        !! convergence is achieved when all criteria are met.

        REAL(DP) :: etot_conv_thr = 1.0E-4_DP
        !! convergence criterion for ion minimization.
        !! This criterion is met when \(\text{etot}(n+1)-\text{etot}(n)
        !! < \text{etot_conv_thr}\), where "n" is the step index,
        !! \(\text{etot}\) the DFT energy. Convergence is achieved when
        !! all criteria are met.

        REAL(DP) :: forc_conv_thr = 1.0E-3_DP
        !! convergence criterion for ion minimization-
        !! This criterion is met when \(\text{MAXVAL(fion)} < 
        !! \text{forc_conv_thr}\), where \(\text{fion}\) are the atomic
        !! forces. Convergence is achieved when all criteria are met.

        CHARACTER(len=80) :: disk_io = 'default'
        !! Specify the amount of I/O activities.

        LOGICAL :: twochem = .FALSE.
        !!if TRUE, the system is simulated using two chemical potentials,
        !!one for the electrons and one for the holes (photoexcited system)

        LOGICAL :: symmetry_with_labels = .FALSE. 
        !! if .TRUE. the  symmetries are checked comparing the first two letters of the 
        !! atomic species labels and the collinear or vector magnetization  

        LOGICAL :: use_spinflip = .FALSE. 
        !! if .TRUE. symmetries in the collinear magnetic case
        !! can  include the spinflip operation. 

        LOGICAL :: tefield  = .false.
        !! if TRUE a sawtooth potential simulating a finite electric field
        !! is added to the local potential - only used in PW

          ! TB - added gate also below to the namelist
        LOGICAL :: gate = .FALSE.
        !! if TRUE a charged plate in charged systems is added with a
        !! total charge which is opposite to the charge of the system

          LOGICAL :: tefield2  = .false.
          !! if TRUE a second finite electric field is added to the local potential
          !! only used in CP

        LOGICAL :: lelfield = .false.
        !! if TRUE a static homogeneous electric field is present
        !! via the modern theory of polarizability - differs from tefield!

        LOGICAL :: lorbm = .false.
        !! if TRUE an orbital magnetization is computed (Kubo terms)

        LOGICAL :: dipfield = .false.
        !! if TRUE the dipole field is subtracted only used in PW for 
        !! surface calculations

        LOGICAL :: lberry = .false.
        !! if TRUE, use modern theory of the polarization

        INTEGER :: gdir = 0
        !! G-vector for polarization calculation (related to \(\(\text{lberry}\)).
        !! Used in PW only

        INTEGER :: nppstr = 0
        !! number of k-points (parallel vector) (related to \(\text{lberry}\)).
        !! Used in PW only.

        INTEGER  :: nberrycyc = 1
        !! number of covergence cycles on electric field

        LOGICAL :: wf_collect = .true.
        !! This flag controls the way wavefunctions are stored to disk:  
        !! TRUE:  collect wavefunctions from all processors, store them
        !!        into a single restart file on a single processors;  
        !! FALSE: do not collect wavefunctions, store them into distributed
        !!        files - NO LONGER IMPLEMENTED SINCE v.6.3 .

        LOGICAL :: saverho = .true.
        !! This flag controls the saving of charge density in CP codes:  
        !! TRUE:  save charge density to restart dir;  
        !! FALSE: do not save charge density.

        LOGICAL :: tabps = .false.
        !! for ab-initio pressure and/or surface calculations

        LOGICAL :: use_wannier = .false.
        !! use or not Wannier functions

        LOGICAL :: lecrpa = .FALSE.
        !! if TRUE symmetry in scf run is neglected for RPA Ec calculation
        ! 
        LOGICAL :: lfcp = .FALSE.
        !! FCP calculation. enabled if TRUE

        LOGICAL :: tqmmm = .FALSE.
        !! QM/MM coupling. enabled if TRUE

        LOGICAL :: trism = .FALSE.    ! 3D-RISM/KS-DFT coupling. enabled if .true.

        CHARACTER(len=256) :: vdw_table_name = ' '

        CHARACTER(len=10) :: point_label_type='SC'

        CHARACTER(len=80) :: memory = 'default' 
        !! controls memory usage 
        CHARACTER(len=80) :: memory_allowed(3)
        !! \(\text{small}\): QE tries to use (when implemented) algorithms using less memory,
        !! even if they are slower than the default;  
        !! \(\text{large}\): QE tries to use (when implemented) algorithms using more memory
        !! to enhance performance.
        DATA memory_allowed / 'small', 'default', 'large' /

        !
        CHARACTER(len=256) :: input_xml_schema_file = ' '
        !! location of xml input according to xsd schema

        NAMELIST / control / title, calculation, verbosity, restart_mode, &
          nstep, iprint, isave, tstress, tprnfor, dt, ndr, ndw, outdir,   &
          prefix, wfcdir, max_seconds, ekin_conv_thr, etot_conv_thr,      &
          forc_conv_thr, pseudo_dir, disk_io, tefield, dipfield, lberry,  &
          gdir, nppstr, wf_collect, lelfield, nberrycyc, refg,            &
          tefield2, saverho, tabps, use_wannier, lecrpa,                  &
          lfcp, tqmmm, vdw_table_name, lorbm, memory, point_label_type,   &
          input_xml_schema_file, gate, trism, twochem, use_spinflip,      &
          symmetry_with_labels 
!
!=----------------------------------------------------------------------------=!
!  SYSTEM Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        INTEGER :: ibrav = 14
        !! index of the the Bravais lattice.  
        !! Note: in variable cell CP molecular dynamics, usually one does
        !!       not want to put constraints on the cell symmetries, thus
        !!       ibrav = 14 is used

        REAL(DP) :: celldm(6) = 0.0_DP
        !! dimensions of the cell (lattice parameters and angles)
        
        ! Alternate definition of the cell - use either this or celldm
        REAL(DP) :: a = 0.0_DP
        !! lattice parameters (a,b,c)
        REAL(DP) :: c = 0.0_DP
        REAL(DP) :: b = 0.0_DP
        REAL(DP) :: cosab = 0.0_DP
        !! lattice angles (cosab, cosac, cosbc)
        REAL(DP) :: cosac = 0.0_DP
        REAL(DP) :: cosbc = 0.0_DP

        REAL(DP) :: ref_alat = 0.0_DP
        !! reference cell alat in a.u. (see REF\_CELL\_PARAMETERS card) 

        INTEGER :: nat = 0
        !! total number of atoms

        INTEGER :: ntyp = 0
        !! number of atomic species

        INTEGER :: nbnd = 0
        !! number of electronic states, this parameter is MANDATORY in CP

        REAL(DP):: tot_charge = 0.0_DP
        !! total system charge

        REAL(DP) :: tot_magnetization = -10000.0_DP
        !! majority - minority spin. A value = -10000 means unspecified

        REAL(DP) :: ecutwfc = 0.0_DP
        !! energy cutoff for wave functions in k-space ( in Rydberg ).
        !! this parameter is MANDATORY

        REAL(DP) :: ecutrho = 0.0_DP
        !! energy cutoff for charge density in k-space (in Rydberg).
        !! By default its value is \(4 \text{ecutwfc}\)

        INTEGER :: nr1 = 0
        !! dimensions of the real space grid for charge and potentials.
        !! Presently NOT used in CP.
        INTEGER :: nr2 = 0
        INTEGER :: nr3 = 0

        INTEGER :: nr1s = 0
        !! dimensions of the real space grid for wavefunctions.
        !! Presently NOT used in CP
        INTEGER :: nr2s = 0
        INTEGER :: nr3s = 0

        INTEGER :: nr1b = 0
        !! dimensions of the "box" grid for Ultrasoft pseudopotentials
        INTEGER :: nr2b = 0
        INTEGER :: nr3b = 0

        CHARACTER(len=80) :: occupations = 'fixed'
        !! select the way electronic states are filled.
        !! See card 'OCCUPATIONS' if \(\text{ocupations}=\text{from_input}\)

        CHARACTER(len=80) :: smearing = 'gaussian'
        !! select the way electronic states are filled for metalic systems

        REAL(DP) :: degauss = 0.0_DP
        !! parameter for the smearing functions - NOT used in CP

        INTEGER :: nspin = 1
        !! number of spinors:  
        !! \(\text{nspin}=1\) for LDA simulations;  
        !! \(\text{nspin}=2\) for LSD simulations;  
        !! \(\text{nspin}=4\) for NON COLLINEAR simulations.

        LOGICAL :: nosym = .true.
        !! (do not) use symmetry
        LOGICAL :: noinv = .false.
        !! (do not) use q => -q symmetry in k-point generation
        LOGICAL :: nosym_evc = .false.
        !! if TRUE use symmetry only to symmetrize k points

        LOGICAL :: force_symmorphic = .false.
        !! if TRUE disable fractionary translations (nonsymmorphic groups)
        LOGICAL :: use_all_frac = .false.
        !! if TRUE enable usage of all fractionary translations, 
        !! disabling check if they are commensurate with FFT grid

        REAL(DP) :: ecfixed = 0.0_DP, qcutz = 0.0_DP, q2sigma = 0.0_DP
        !! parameters for modified kinetic energy functional to be used
        !! in variable-cell constant cut-off simulations

        CHARACTER(len=80) :: input_dft = 'none'
        !! Variable used to overwrite dft definition contained in
        !! pseudopotential files; 'none' means DFT is read from pseudos.
        !! Used in PW only - allowed values: any legal DFT value

        REAL(DP) :: starting_charge( nsx ) = 0.0_DP
        !! PW ONLY

        REAL(DP) :: starting_magnetization( nsx ) = 0.0_DP
        !! PW ONLY

        !  PARAMETERS FOR TWO-CHEM-CALCULATIONS
        REAL(DP) :: degauss_cond = 0.0_DP 
        !broadening for conduction band
        INTEGER ::  nbnd_cond = 0 
        ! n_bands in conduction
        REAL(DP) :: nelec_cond =0.0_DP 
        !number of electrons in the conduction bands

        ! DFT+Hubbard
        ! Old input parameters in the SYSTEM naqmelist (removed since v7.1):
        CHARACTER(len=80) :: U_projection_type = ''  ! obsolete
        CHARACTER(len=80) :: Hubbard_parameters = '' ! obsolete
        REAL(DP) :: Hubbard_U_back(nsx)  = 0.0_DP    ! obsolete
        !
        ! the following are the parameters for DFT+Hubbard
        LOGICAL :: lda_plus_u = .false.              
        INTEGER :: lda_plus_u_kind = -1            
        INTEGER, PARAMETER :: nspinx=2 ! lqmax is taken from upf_params
        REAL(DP) :: starting_ns_eigenvalue(lqmax,nspinx,nsx) = -1.0_DP
        INTEGER  :: Hubbard_l(nsx)  = -1
        INTEGER  :: Hubbard_n(nsx)  = -1
        INTEGER  :: Hubbard_l2(nsx) = -1
        INTEGER  :: Hubbard_n2(nsx) = -1
        INTEGER  :: Hubbard_l3(nsx) = -1
        INTEGER  :: Hubbard_n3(nsx) = -1
        REAL(DP) :: Hubbard_U(nsx)  = 0.0_DP
        REAL(DP) :: Hubbard_U2(nsx) = 0.0_DP
        REAL(DP) :: Hubbard_Um(lqmax,nspinx,nsx) = 0.0_DP
        REAL(DP) :: Hubbard_Um_nc(2*lqmax,nsx) = 0.0_DP
        REAL(DP) :: Hubbard_V(natx,natx*(2*sc_size+1)**3,4) = 0.0_DP 
        REAL(DP) :: Hubbard_J0(nsx) = 0.0_DP
        REAL(DP) :: Hubbard_J(3,nsx) = 0.0_DP
        REAL(DP) :: Hubbard_alpha(nsx) = 0.0_DP
        REAL(DP) :: Hubbard_alpha_back(nsx) = 0.0_DP
        REAL(DP) :: Hubbard_alpha_m(lqmax,nspinx,nsx) = 0.0_DP
        REAL(DP) :: Hubbard_alpha_m_nc(2*lqmax,nsx) = 0.0_DP
        REAL(DP) :: Hubbard_beta(nsx) = 0.0_DP
        REAL(DP) :: Hubbard_occ(nsx,3) = -1.0_DP
        CHARACTER(len=80) :: Hubbard_projectors = ''
        LOGICAL :: reserv(nsx) = .FALSE.
        LOGICAL :: reserv_back(nsx) = .FALSE.
        LOGICAL :: hub_pot_fix = .FALSE.
        LOGICAL :: orbital_resolved = .FALSE.
        LOGICAL :: backall(nsx) = .FALSE.

          ! For linking to DMFT calculations
        LOGICAL :: dmft = .FALSE.
        CHARACTER(len=256) :: dmft_prefix = 'dmft_prefix'

        LOGICAL :: la2F = .false.
          ! For electron-phonon calculations
          !
        LOGICAL :: step_pen=.false.
        REAL(DP) :: A_pen(10,nspinx) = 0.0_DP
        REAL(DP) :: sigma_pen(10) = 0.01_DP
        REAL(DP) :: alpha_pen(10) = 0.0_DP

          ! next group of variables PWSCF ONLY
          ! 
          !
        REAL(DP) :: exx_fraction = -1.0_DP
        !! exact exchange fraction. If negative, use defaults
        REAL(DP) :: screening_parameter = -1.0_DP
        INTEGER  :: nqx1 = 0
        !! use the same values as \(\text{nk1, nk2, nk3}\)
        INTEGER  :: nqx2 = 0, nqx3=0  
        !
        REAL(DP) :: gau_parameter = -1.0_DP
        !! gau-pbe
          !
        CHARACTER(len=80) :: exxdiv_treatment = 'gygi-baldereschi'
        !! define how ro cure the Coulomb divergence in EXX.
        CHARACTER(len=80) :: exxdiv_treatment_allowed(6)
        !! Allowed values for \(\text{exxdiv_treatment}\)
        DATA exxdiv_treatment_allowed / 'gygi-baldereschi', 'gygi-bald', 'g-b',&
             'vcut_ws', 'vcut_spherical', 'none' /
          !
        LOGICAL  :: x_gamma_extrapolation = .true.
        REAL(DP) :: yukawa = 0.0_DP
        REAL(DP) :: ecutvcut = 0.0_DP
          ! auxiliary variables to define exxdiv treatment
        LOGICAL  :: adaptive_thr = .FALSE.
        REAL(DP) :: conv_thr_init = 0.001_DP
        REAL(DP) :: conv_thr_multi = 0.1_DP
        REAL(DP) :: ecutfock = -1.d0
          ! variables used in Lin Lin's ACE and SCDM
        REAL(DP) :: localization_thr = 0.0_dp, scdmden=1.0d0, scdmgrd=1.0d0
        INTEGER  :: nscdm = 1        
        INTEGER  :: n_proj  = 0
        LOGICAL  :: scdm=.FALSE.
        LOGICAL  :: ace=.TRUE.

          ! parameters for external electric field
        INTEGER  :: edir = 0
        REAL(DP) :: emaxpos = 0.0_DP
        REAL(DP) :: eopreg = 0.0_DP
        REAL(DP) :: eamp = 0.0_DP

          ! TB parameters for charged plate representing the gate
          ! and a possible potential barrier
          ! added also below to the namelist
        REAL(DP) :: zgate = 0.5
        LOGICAL  :: relaxz = .false.
        LOGICAL  :: block = .false.
        REAL(DP) :: block_1 = 0.45
        REAL(DP) :: block_2 = 0.55
        REAL(DP) :: block_height = 0.1

          ! Various parameters for noncollinear calculations
        LOGICAL  :: noncolin = .false.
        LOGICAL  :: lspinorb = .false.
        LOGICAL  :: lforcet=.FALSE.
        LOGICAL  :: starting_spin_angle=.FALSE.
        REAL(DP) :: lambda = 1.0_DP
        REAL(DP) :: fixed_magnetization(3) = 0.0_DP
        REAL(DP) :: angle1(nsx) = 0.0_DP
        REAL(DP) :: angle2(nsx) = 0.0_DP
        INTEGER  :: report =-1
        LOGICAL  :: no_t_rev = .FALSE.

        CHARACTER(len=80) :: constrained_magnetization = 'none'
        REAL(DP) :: B_field(3) = 0.0_DP
        !! A fixed magnetic field defined by the vector B_field is added
        !! to the exchange and correlation magnetic field.

        CHARACTER(len=80) :: sic = 'none'
        !! CP only - SIC correction (D'avezac Mauri)
        !
        ! Parameters for SIC calculation
        REAL(DP) :: sic_epsilon = 0.0_DP
        REAL(DP) :: sic_alpha   = 0.0_DP
        LOGICAL   :: force_pairing = .false.
        CHARACTER(len=80) :: pol_type = 'none'
        REAL(DP) :: sic_gamma = 0.0_DP
        LOGICAL  :: sic_energy = .false.
        REAL(DP) :: sci_vb = 0.0_DP
        REAL(DP) :: sci_cb = 0.0_DP

        LOGICAL :: one_atom_occupations=.false.

        CHARACTER(len=80) :: assume_isolated = 'none'
        !! Possible corrections for isolated systems:  
        !!  'none', 'makov-payne', 'martyna-tuckerman', 'esm'  
        !! Plus ENVIRON-specific:  
        !!  'slabx', 'slaby', 'slabz', 'pcc'

        CHARACTER(len=80) :: vdw_corr = 'none'
        !! semi-empirical van der Waals corrections (not to be confused with 
        !! nonlocal functionals, specified in input_dft!). Default is 'none',
        !! allowed values:  
        !! 'dft-d' or 'grimme-d2' [S.Grimme, J.Comp.Chem. 27, 1787 (2006)]  
        !! 'ts', 'ts-vdW', 'tkatchenko-scheffler'
        !!  (Tkatchenko & Scheffler, Phys. Rev. Lett. 102, 073005 (2009))  
        !!  'xdm' (Otero de la Roza and Johnson, J. Chem. Phys. 136 (2012) 174109)

        LOGICAL   :: london = .false.
        !! OBSOLESCENT: same as vdw_corr='grimme-d2'  
        !
        ! Other DFT-D parameters ( see Modules/mm_dispersion.f90 )  
        !
        REAL(DP) :: london_s6   =   0.75_DP
        !! default global scaling parameter for PBE
        REAL(DP) :: london_rcut = 200.00_DP
        REAL(DP) :: london_c6( nsx ) = -1.0_DP
        !! user specified atomic C6 coefficients
        REAL(DP) :: london_rvdw( nsx ) = -1.0_DP
        !! user specified atomic vdw radii

          ! Grimme-D3 (DFT-D3) dispersion correction.
          ! version=3 is Grimme-D3, version=2 is Grimme-D2,
          ! version=3 dftd3_threebody = .false. similar to grimme d2 but with
          ! the two_body parameters of d3
        integer  ::    dftd3_version = 3
        logical ::     dftd3_threebody = .true.

        LOGICAL   :: ts_vdw = .false.
        !! OBSOLESCENT: same as vdw_corr='Tkatchenko-Scheffler'
        LOGICAL   :: mbd_vdw = .false.
        !! added for consistency with \(\text{ts_vdw}\)
        LOGICAL :: ts_vdw_isolated = .FALSE.
        !! if TRUE, TS-vdW correction for isolated system;  
        !! if FALSE, TS-vdW correction for periodic system.
        REAL(DP) :: ts_vdw_econv_thr = 1.0E-6_DP
        !! convergence criterion for TS-vdW energy for periodic system 
        !

        LOGICAL :: xdm = .FALSE.
        !! OBSOLESCENT: same as vdw_corr='xdm'
        REAL(DP) :: xdm_a1 = 0.0_DP
        REAL(DP) :: xdm_a2 = 0.0_DP
        !
        CHARACTER(LEN=3) :: esm_bc = 'pbc'
        !! 'pbc': ordinary calculation with periodic boundary conditions;  
        !! 'bc1': vacuum-slab-vacuum;  
        !! 'bc2': metal-slab-metal;  
        !! 'bc3': vacuum-slab-metal.

        REAL(DP) :: esm_efield = 0.0_DP
        !! applied electronic field [Ryd/a.u.] (used only for esm_bc='bc2')

        REAL(DP) :: esm_w = 0.0_DP
        !! position of effective screening medium from z0=L_z/2 [a.u.]
        !! note: z1 is given by \(z1=z0+\text{abs}(\text{esm_w})\)

        REAL(DP) :: esm_a = 0.0_DP
        !! smoothness parameter for smooth-ESM \(\exp{2a(z-z1)}\)

        REAL(DP) :: esm_zb = 0.0_DP
        !! smearing width for Ewald summation (RSUM) in smooth-ESM

        INTEGER :: esm_nfit = 4
        !! number of z-grid points for polynomial fitting at cell edge

        LOGICAL :: esm_debug = .FALSE.
        !! used to enable debug mode (output \(\text{v_hartree}\) and
        !! \(\text{v_local}\))

        INTEGER :: esm_debug_gpmax = 0
        !! if esm_debug is TRUE, calculate \(\text{v_hartree}\) and \(\text{v_local}\)
        !! for \(\text{abs}(gp)\leq\text{esm_debug_gpmax}\) (gp is integer and has 
        !! \(\text{tpiba}\) unit)

        LOGICAL :: lgcscf = .FALSE.
        !! if TRUE, GC-SCF is used

        LOGICAL :: gcscf_ignore_mun = .FALSE.
        !! if TRUE, ignore the term of -mu * N

        REAL(DP) :: gcscf_mu = 0.0_DP
        !! target Fermi energy of GC-SCF (in eV)

        REAL(DP) :: gcscf_conv_thr = 1.0E-2_DP
        !! convergence threshold of GC-SCF (in eV)

        REAL(DP) :: gcscf_gk = 0.4_DP
        !! wavenumber shift for Kerker operator (in 1/bohr)

        REAL(DP) :: gcscf_gh = 1.5_DP
        !! wavenumber shift for Hartree metric (in 1/bohr)

        REAL(DP) :: gcscf_beta = 0.05_DP
        !! mixing rate of Fermi energy

        INTEGER :: space_group = 0
        !! space group number for coordinates given in crystallographic form
        !
        LOGICAL :: uniqueb=.FALSE.
        !! if TRUE for monoclinic lattice choose the b unique primitive 
        !! vectors
        !
        INTEGER :: origin_choice = 1 
        !! for space groups that have more than one origin choice, choose
        !! the origin (can be 1 or 2)
        !
        LOGICAL :: rhombohedral = .TRUE.
        !! if TRUE for rhombohedral space groups give the coordinates 
        !! in rhombohedral axes. If FALSE in hexagonal axes, that are
        !! converted internally in rhombohedral axes.  
        !
        INTEGER :: nextffield = 0 
        !! Number of activated external force fields 
        !



        NAMELIST / system / ibrav, celldm, a, b, c, cosab, cosac, cosbc, nat, &
             ntyp, nbnd, ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s,         &
             nr3s, nr1b, nr2b, nr3b, nosym, nosym_evc, noinv, use_all_frac,   &
             force_symmorphic, starting_charge, starting_magnetization,       &
             occupations, degauss, nspin, ecfixed, qcutz, q2sigma,            &
             degauss_cond,nbnd_cond,nelec_cond,                       &
             lda_plus_u, lda_plus_u_kind, U_projection_type, Hubbard_parameters, & ! obsolete
             Hubbard_U, Hubbard_J0, Hubbard_J, Hubbard_V, Hubbard_U_back,     & ! moved to HUBBARD card 
             Hubbard_alpha, Hubbard_alpha_back, Hubbard_beta, Hubbard_occ,    &
             hub_pot_fix, orbital_resolved, reserv, reserv_back, dmft,        &
             dmft_prefix, edir, emaxpos, eopreg, eamp, smearing,              &
             starting_ns_eigenvalue, input_dft, la2F, assume_isolated,        &
             nqx1, nqx2, nqx3, ecutfock, localization_thr, scdm, ace,         &
             scdmden, scdmgrd, nscdm, n_proj,                                 &
             exxdiv_treatment, x_gamma_extrapolation, yukawa, ecutvcut,       &
             exx_fraction, screening_parameter, ref_alat,                     &
             noncolin, lspinorb, starting_spin_angle, lambda, angle1, angle2, &
             report, lforcet,                                                 &
             constrained_magnetization, B_field, fixed_magnetization,         &
             sic, sic_epsilon, force_pairing, sic_alpha,                      &
             pol_type, sic_gamma, sic_energy, sci_vb, sci_cb,                 &
             tot_charge, tot_magnetization, one_atom_occupations,             &
             vdw_corr, london, london_s6, london_rcut, london_c6, london_rvdw,&
             dftd3_version, dftd3_threebody,                                  &
             ts_vdw, ts_vdw_isolated, ts_vdw_econv_thr,                       &
             xdm, xdm_a1, xdm_a2,                                             &
             mbd_vdw,                                                         &
             step_pen, A_pen, sigma_pen, alpha_pen, no_t_rev,                 &
             esm_bc, esm_efield, esm_w, esm_nfit, esm_debug, esm_debug_gpmax, &
             esm_a, esm_zb,                                                   &
             lgcscf, gcscf_ignore_mun, gcscf_mu, gcscf_conv_thr,              &
             gcscf_gk, gcscf_gh, gcscf_beta,                                  &
             space_group, uniqueb, origin_choice, rhombohedral,               &
             zgate, relaxz, block, block_1, block_2, block_height,            &
             nextffield

!=----------------------------------------------------------------------------=!
!  ELECTRONS Namelist Input Parameters
!=----------------------------------------------------------------------------=!

        REAL(DP) :: emass = 0.0_DP
        !! effective electron mass in the CP Lagrangian,
        !! in atomic units (\(1 \text{a.u. of mass} = 1/1822.9 \text{a.m.u.} = 
        !! 9.10939\cdot 10^{-31} kg\)).  
        !! Typical values in CP simulation are between 100. and 1000.

        REAL(DP) :: emass_cutoff = 0.0_DP
        !! mass cut-off (in Rydbergs) for the Fourier acceleration.  
        !! Effective mass is rescaled for "G" vector components with kinetic
        !! energy above \(\text{emass_cutoff}\)
        !! Use a value greater than \(\text{ecutwfc}\) to disable Fourier 
        !! acceleration.

        CHARACTER(len=80) :: orthogonalization = 'ortho'
        !! selects the orthonormalization method for electronic wave functions:  
        !! -'Gram-Schmidt'  uses Gram-Schmidt algorithm;  
        !! -'ortho'         uses iterative algorithm.

        REAL(DP) :: ortho_eps = 1.E-9_DP
        !! meaningful only if orthogonalization = 'ortho'.
        !! Tolerance for iterative orthonormalization,
        !! a value of 1.d-8 is usually sufficent.

        INTEGER   :: ortho_max = 50
        !! meaningful only if orthogonalization = 'ortho'.
        !! Maximum number of iterations for orthonormalization
        !! usually between 20 and 300.

        INTEGER :: exx_maxstep = 1000
        !! maximum number of steps in the outer loop of electronic minimization
        !! when exx is active (hybrid functionals).
        INTEGER :: electron_maxstep = 1000
        !! maximum number of steps in electronic minimization.
        !! This parameter applies only when using 'cg' electronic or
        !! ionic dynamics and electron\_dynamics = 'CP-BO'
        LOGICAL :: scf_must_converge = .true.
        !! stop or continue if SCF does not converge

        CHARACTER(len=80) :: electron_dynamics = 'none'
        !! set how electrons should be moved
        CHARACTER(len=80) :: electron_dynamics_allowed(7)
        DATA electron_dynamics_allowed &
          / 'default', 'sd', 'cg', 'damp', 'verlet', 'none', 'cp-bo' /

        REAL(DP) :: electron_damping = 0.0_DP
        !! meaningful only if electron\_dynamics = 'damp'.
        !! Damping frequency times delta t, optimal values could be
        !! calculated with the formula:  
        !!        \(\sqrt{0.5\log{(E1-E2)/(E2-E3)}}\)  
        !! where E1 E2 E3 are successive values of the DFT total energy
        !! in a steepest descent simulations

        CHARACTER(len=80) :: electron_velocities = 'default'
        !! electron velocities:  
        !! -'zero'    restart setting electronic velocities to zero;  
        !! -'default' restart using electronic velocities of the previous run.

        CHARACTER(len=80) :: electron_temperature = 'not_controlled'
        !! electron temperature:  
        !! -'nose'           control electronic temperature using Nose thermostat
        !!                  see parameter "fnosee" and "ekincw";  
        !! -'rescaling'      control electronic temperature via velocities rescaling;  
        !! -'not_controlled' electronic temperature is not controlled-

        REAL(DP) :: ekincw = 0.0_DP
        !! meaningful only with electron\_temperature different from 'not\_controlled'.
        !! Value of the average kinetic energy (in atomic units) forced by the 
        !! temperature control

        REAL(DP) :: fnosee = 0.0_DP
        !! meaningful only with "electron\_temperature = 'nose' ".
        !! Oscillation frequency of the nose thermostat (in terahertz)

        CHARACTER(len=80) :: startingwfc = 'random'
        !! define how the code should initialize the wave function:  
        !! -'atomic'   start from superposition of atomic wave functions;  
        !! -'atomic+random' as above, plus randomization;  
        !! -'random'   start from random wave functions;  
        !! -'file'     read wavefunctions from file.

        REAL(DP) :: ampre = 0.0_DP
        !! meaningful only if startingwfc='random'. Amplitude of the
        !! randomization ( allowed values: 0.0 - 1.0 )

        REAL(DP) :: grease = 0.0_DP
        !! a number smaller or equal to 1, very close to 1: the damping in
        !! electronic damped dynamics is multiplied at each time step by 
        !! "grease" (avoids overdamping close to convergence: obsolete?).
        !! grease = 1 : normal damped dynamics. Used only in CP

        INTEGER :: diis_size = 0
        !! meaningful only with electron\_dynamics='diis'.
        !! Size of the matrix used for the inversion in the iterative subspace
        !! default is 4, allowed value 1-5

        INTEGER :: diis_nreset = 0
        !! meaningful only with electron\_dynamics='diis'.
        !! Number of steepest descendent step after a reset of the diis
        !! iteration, default value is 3

        REAL(DP) :: diis_hcut = 0.0_DP
        !! meaningful only with electron\_dynamics = 'diis'.
        !! Energy cutoff (a.u.), above which an approximate diagonal
        !! hamiltonian is used in finding the direction to the minimum
        !! default is 1.0

        REAL(DP) :: diis_wthr = 1.E-4_DP
        !! meaningful only with electron\_dynamics='diis'. Convergence threshold
        !! for wave function. This criterion is satisfied when the maximum change
        !! in the wave functions component between two diis steps is less than 
        !! this threshold. Default value is \(\text{ekin_conv_thr}\).

        REAL(DP) :: diis_delt = 1.0_DP
        !! meaningful only with electron\_dynamics='diis'.
        !! Electronic time step used in the steepest descendent step.
        !! Default is \(\text{dt}\)

        INTEGER :: diis_maxstep = 100
        !! meaningful only with electron_dynamics='diis'.
        !! Maximum number of iteration in the diis minimization.
        !! Default is \(\text{electron_maxstep}\).

        LOGICAL :: diis_rot = .false.
        !! meaningful only with electron\_dynamics='diis'.
        !! If TRUE enable diis with charge mixing and rotations.
        !! Default value is FALSE.

        REAL(DP) :: diis_fthr = 1.E-3_DP
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! Convergence threshold for ionic force. This criterion is satisfied
        !! when the maximum change in the atomic force between two diis steps
        !! is less than this threshold. Default value is 0.0.

        REAL(DP) :: diis_temp = 0.0_DP
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        ! Electronic temperature, significant only if ???

        REAL(DP) :: diis_achmix  = 0.0_DP
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! "A" parameter in the charge mixing formula:  
        !! \(\text{chmix}=A\cdot G^2/(G^2 + G0^2)\), where G represents reciprocal
        !! lattice vectors.

        REAL(DP) :: diis_g0chmix  = 0.0_DP
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! \(G0^2\) is the parameter in the charge mixing formula.

        INTEGER :: diis_nchmix = 0
        !! meaningful only with electron\_dynamics='diis' and "diis\_rot=TRUE.
        !! Dimension of the charge mixing.

        REAL(DP) :: diis_g1chmix = 0.0_DP
        !! meaningful only with electron\_dynamics='diis' and "diis\_rot=TRUE.
        !! \(G1^2\) parameter in the charge mixing formula.  
        !! \(\text{metric}=(G^2+G1^2)/G^2\), G represents reciprocal lattice vectors

        INTEGER :: diis_nrot(3) = 0
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! Start upgrading the charge density every diis\_nrot(1) steps,
        !! then every diis\_nrot(2), and at the end every diis\_nrot(3),
        !! depending on diis\_rothr

        REAL(DP) :: diis_rothr(3) = 1.E-4_DP
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! Threshold on the charge difference between two diis step.
        !! When max charge difference is less than diis\_rothr(1), switch
        !! between the diis\_nrot(1) upgrade frequency to diis\_nrot(2),
        !! then when the max charge difference is less than diis\_rothr(2),
        !! switch between diis\_nrot(2) and diis\_nrot(3), upgrade frequency,
        !! finally when the max charge difference is less than diis\_nrot(3)
        !! convergence is achieved.

        REAL(DP) :: diis_ethr = 1.E-4_DP
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! Convergence threshold for energy. This criterion is satisfied when 
        !! the change in the energy between two diis steps is less than this 
        !! threshold. Default value is \(\text{etot_conv_thr}\)

        LOGICAL :: diis_chguess = .false.
        !! meaningful only with electron\_dynamics='diis' and diis\_rot=TRUE.
        !! If diis\_chguess=TRUE enable charge density guess between two diis step,
        !! defaut value is FALSE.

        CHARACTER(len=80) :: mixing_mode = 'default'
        !! type of mixing algorithm for charge self-consistency.
        !! Used in PWscf only.

        REAL(DP) :: mixing_beta = 0.0_DP
        !! parameter for mixing algorithm. Used in PWscf only.

        INTEGER :: mixing_ndim = 0
        !! dimension of mixing subspace. Used in PWscf only.

        CHARACTER(len=80) :: diagonalization = 'david'
        !! diagonalization = 'david', 'cg', 'paro' or 'rmm'
        !! algorithm used by PWscf for iterative diagonalization

        REAL(DP) :: diago_thr_init = 0.0_DP
        !! convergence threshold for the first iterative diagonalization.
        !! Used in PWscf only.

        INTEGER :: diago_cg_maxiter = 100
        !! max number of iterations for the first iterative diagonalization.
        !! Using conjugate-gradient algorithm - used in PWscf only.

        INTEGER :: diago_david_ndim = 4
        !! dimension of the subspace used in Davidson diagonalization
        !! used in PWscf only.

        INTEGER :: diago_rmm_ndim = 4
          ! dimension of the subspace used in RMM-DIIS diagonalization
          ! used only in PWscf

        LOGICAL :: diago_rmm_conv = .false.
          ! if .TRUE., RMM-DIIS is performed up to converge
          ! if .FALSE., RMM-DIIS is performed only once
          ! used only in PWscf

        INTEGER :: diago_gs_nblock = 16
          ! blocking size in Gram-Schmidt orthogonalization
          ! used only in PWscf

        LOGICAL :: diago_full_acc = .false.

        REAL(DP) :: conv_thr = 1.E-6_DP
        !! convergence threshold in electronic ONLY minimizations.
        !! Used in PWscf only.
        !! Used with electron\_dynamics = 'CP-BO' in CP and PWscf

        INTEGER :: mixing_fixed_ns  = 0
        !! For DFT+U calculations, PWscf only.

        CHARACTER(len=80) :: startingpot = 'potfile'
        !! specify the file containing the DFT potential of the system.
        !! Used in PWscf only.

        INTEGER :: n_inner = 2
        !! number of inner loop per CG iteration. Used in CP only.

        INTEGER :: niter_cold_restart = 1
        !! frequency of full cold smearing inner cycle (in iterations).

        REAL(DP) :: lambda_cold
        !! step for not complete cold smearing inner cycle.

        LOGICAL :: tgrand = .false.
        !! whether to do grand-canonical calculations.

        REAL(DP) :: fermi_energy = 0.0_DP
        !! chemical potential of the grand-canonical ensemble.

        CHARACTER(len=80) :: rotation_dynamics = "line-minimization"
        !! evolution of the rotational degrees of freedom.

        CHARACTER(len=80) :: occupation_dynamics = "line-minimization"
        !! evolution of the occupational degrees of freedom.

        REAL(DP) :: rotmass = 0
        !! mass for the rotational degrees of freedom.

        REAL(DP) :: occmass = 0
        !! mass for the occupational degrees of freedom.

        REAL(DP) :: occupation_damping = 0
        !! damping for the rotational degrees of freedom.

        REAL(DP) :: rotation_damping = 0
        !! damping for the occupational degrees of freedom.

        LOGICAL :: tcg = .true.
        !! if TRUE perform in cpv conjugate gradient minimization of electron energy

        LOGICAL :: pre_state = .false.
        !! if TRUE, in CP's conjugate gradient routine, precondition each band
        !! with its kinetic energy (see CPV/src/cg_sub.f90)

        INTEGER :: maxiter = 100
        !! max number of conjugate gradient iterations

        REAL(DP)  :: etresh =1.0E-7_DP
        !! treshhold on energy

        REAL(DP) :: passop =0.3_DP
        !! small step for parabolic interpolation

        INTEGER :: niter_cg_restart
        !! frequency of restart for the conjugate gradient algorithm in iterations.

        INTEGER  :: epol = 3
        !! electric field direction.

        REAL(DP) :: efield =0.0_DP
        !! electric field intensity in atomic units.

          LOGICAL :: real_space = .false.
          !! real_space routines for US pps.

        REAL(DP) :: efield_cart(3)
        !! electric field vector in cartesian system of reference.
        
       CHARACTER(len=80) :: efield_phase='none'
       !! for Berry's phase electric field selection of string phases.

       INTEGER  :: epol2 = 3
       !! electric field direction

        REAL(DP) :: efield2 =0.0_DP
        !! electric field intensity in atomic units.

        LOGICAL :: tqr = .false.
        !! US contributions are added in real space.

        LOGICAL :: tq_smoothing = .false.
        !! US augmentation charge is smoothed before use.

        LOGICAL :: tbeta_smoothing = .false.
        !! beta function are smoothed before use.

        LOGICAL :: occupation_constraints = .false.
        !! If TRUE perform CP dynamics with constrained occupations.
        !! To be used together with penalty functional ...

        !
        ! ... CP-BO ...
        LOGICAL :: tcpbo = .FALSE.
        !! if TRUE perform CP-BO minimization of electron energy

        REAL(DP) :: emass_emin = 0.0_DP
        !! meaningful only if electron_dynamics = 'CP-BO'.
        !! Effective electron mass used in CP-BO electron minimization in
        !! atomic units (\( 1 \text{a.u. of mass} = 1/1822.9 \text{a.m.u.} =
        !! 9.10939\cdot 10^{-31} kg\))

        REAL(DP) :: emass_cutoff_emin = 0.0_DP
        !! meaningful only if electron_dynamics = 'CP-BO'.
        !! Mass cut-off (in Rydbergs) for the Fourier acceleration in CP-BO
        !! electron minimization.

        REAL(DP) :: electron_damping_emin = 0.0_DP
        !! meaningful only if electron_dynamics = 'CP-BO'.
        !! Damping parameter utilized in CP-BO electron minimization.

        REAL(DP) :: dt_emin = 0.0_DP
        !! meaningful only if electron_dynamics = 'CP-BO'.
        !! Time step for CP-BO electron minimization dynamics, in atomic units.
        !! CP: \(1 \text{a.u. of time} = 2.4189\cdot 10^{-17} s\), PW: twice that much.

        NAMELIST / electrons / emass, emass_cutoff, orthogonalization, &
          exx_maxstep, electron_maxstep, scf_must_converge, ortho_eps, ortho_max, electron_dynamics,   &
          electron_damping, electron_velocities, electron_temperature, &
          ekincw, fnosee, ampre, grease,                               &
          diis_size, diis_nreset, diis_hcut,                           &
          diis_wthr, diis_delt, diis_maxstep, diis_rot, diis_fthr,     &
          diis_temp, diis_achmix, diis_g0chmix, diis_g1chmix,          &
          diis_nchmix, diis_nrot, diis_rothr, diis_ethr, diis_chguess, &
          mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns,      &
          tqr, tq_smoothing, tbeta_smoothing,                          &
          diago_cg_maxiter, diago_david_ndim, diago_rmm_ndim,          &
          diago_rmm_conv, diago_gs_nblock, diagonalization,            &
          startingpot, startingwfc , conv_thr,                         &
          adaptive_thr, conv_thr_init, conv_thr_multi,                 &
          diago_thr_init, n_inner, fermi_energy, rotmass, occmass,     &
          rotation_damping, occupation_damping, rotation_dynamics,     &
          occupation_dynamics, tcg, maxiter, etresh, passop, epol,     &
          efield, epol2, efield2, diago_full_acc,                      &
          occupation_constraints, niter_cg_restart,                    &
          niter_cold_restart, lambda_cold, efield_cart, real_space,    &
          tcpbo,emass_emin, emass_cutoff_emin, electron_damping_emin,  &
          dt_emin, efield_phase, pre_state

!
!=----------------------------------------------------------------------------=!
!  IONS Namelist Input Parameters
!=----------------------------------------------------------------------------=!


        CHARACTER(len=80) :: ion_dynamics = 'none'
        !! set how ions should be moved
        CHARACTER(len=80) :: ion_dynamics_allowed(12)
        !! allowed options for ion\_dynamics.
        DATA ion_dynamics_allowed / 'none', 'sd', 'cg', 'langevin', &
                                    'damp', 'verlet', 'velocity-verlet', 'bfgs', 'beeman',& 
                                    'langevin-smc', 'ipi', 'fire' /

        REAL(DP) :: ion_radius(nsx) = 0.5_DP
        !! pseudo-atomic radius of the i-th atomic species (CP only).
        !! For Ewald summation: typical values range between 0.5 and 2.0 
       INTEGER :: iesr = 1
       !! perform Ewald summation on \(\text{iesr}^3\) cells - CP only

        REAL(DP) :: ion_damping = 0.2_DP
        !! meaningful only if ion\_dynamics='damp'. Damping frequency times
        !! delta t, optimal values could be calculated with the formula:  
        !! \( \sqrt(0.5\cdot\log{(E1-E2)/(E2-E3)}) \)  
        !! where \(E1\ E2\ E3\) are successive values of the DFT total energy
        !! in a ionic steepest descent simulation.

        CHARACTER(len=80) :: ion_positions = 'default'
        !! allowed options:  
        !! -'default'    restart the simulation with atomic positions read
        !!               from the restart file;  
        !! -'from_input' restart the simulation with atomic positions read
        !!               from standard input (see the card 'ATOMIC_POSITIONS').

        CHARACTER(len=80) :: ion_velocities = 'default'
        !! allowed options:  
        !! -'default'    restart the simulation with atomic velocities read
        !!               from the restart file;  
        !! -'random'     start the simulation with random atomic velocities;  
        !! -'from_input' restart the simulation with atomic velocities read
        !!               from standard input (see the card 'ATOMIC_VELOCITIES');  
        !! -'zero'       restart the simulation with atomic velocities set to zero.

        CHARACTER(len=80) :: ion_temperature = 'not_controlled'
        !! allowed options:  
        !! -'nose'           control ionic temperature using Nose thermostat,
        !!                   see parameters "fnosep" and "tempw";  
        !! -'rescaling'      control ionic temperature via velocity rescaling,
        !!                   see parameters "tempw" and "tolp";  
        !! -'rescale-v'      control ionic temperature via velocity rescaling,
        !!                   see parameters "tempw" and "nraise";  
        !! -'rescale-T'      control ionic temperature via velocity rescaling,
        !!                   see parameter "delta_t";  
        !! -'reduce-T'       reduce ionic temperature, see parameters "nraise", delta_t";  
        !! -'berendsen'      control ionic temperature using "soft" velocity
        !!                   rescaling, see parameters "tempw" and "nraise";  
        !! -'andersen'       control ionic temperature using Andersen thermostat,
        !!                   see parameters "tempw" and "nraise";  
        !! -'not_controlled' ionic temperature is not controlled.

        REAL(DP) :: tempw = 300.0_DP
        !! meaningful only when ion\_temperature different from 'not\_controlled'.
        !! Value of the ionic temperature (in Kelvin) forced by the temperature control.

        INTEGER, PARAMETER :: nhclm   = 4
        !! max length for the chain; it can be easily increased
        !! since the restart file should be able to handle it,
        !! perhaps better to align \(\text{nhclm}\) by 4.
        REAL(DP) :: fnosep( nhclm )  = 50.0_DP
        !! meaningful only with ion\_temperature='nose'. Oscillation frequency of
        !! the nose thermostat (in terahertz).
        
        INTEGER   ::  nhpcl = 0
        !! non-zero only with ion\_temperature = 'nose'.
        !! This defines the length of the Nose-Hoover chain.

        INTEGER   :: nhptyp = 0
        !! this parameter set the nose hoover thermostat to more than one

        INTEGER   ::  nhgrp(nsx)=0
        !! this is the array to assign thermostats to atomic types.
        !! Allows to use various thermostat setups.

        INTEGER   ::  ndega = 0
        !! this is the parameter to control active degrees of freedom.
        !! Used for temperature control and the Nose-Hoover chains.

        REAL(DP) :: tolp = 50.0_DP
        !! meaningful only with ion\_temperature='rescaling'.
        !! Tolerance (in Kelvin) of the rescaling. When ionic temperature
        !! differs from "tempw" more than "tolp" apply rescaling.

        REAL(DP)  ::  fnhscl(nsx)=-1.0_DP
        !! this is to scale the target energy, in case there are constraints
        !! the dimension is the same as \(\text{nhgrp}\), meaning that atomic
        !! type \(i\) with a group \(\text{nhgrp}(i)\) is scaled by \(\text{fnhscl}(i)\)

        LOGICAL   :: tranp(nsx) = .false.
        !! \(\text{tranp}(i)\) controls the randomization of the i-th atomic specie:  
        !! -TRUE   randomize ionic positions ( see "amprp" );  
        !! -FALSE  do nothing.

        REAL(DP) :: amprp(nsx) = 0.0_DP
        !! \(\text{amprp}(i)\) meaningful only if \(\text{tranp}(i)=\text{TRUE}\), amplitude
        !! of the randomization (allowed values: 0.0 - 1.0) for the i-th atomic specie.  
        !! Add to the positions a random displacements vector (in Bohr radius) defined as:
        !! \(\text{amprp}(i)\cdot (X,Y,Z)\), where X, Y, Z are pseudo random number in the interval
        !! \([-0.5, 0.5]\).

        REAL(DP) :: greasp = 0.0_DP
        !! same as "grease", for ionic damped dynamics.
        !! NOT used in FPMD

        INTEGER   :: ion_nstepe = 1
        !! number of electronic steps for each ionic step

        INTEGER   :: ion_maxstep = 1000
        !! maximum number of step in ionic minimization

        REAL(DP) :: upscale = 100.0_DP
        !! Max reduction allowed in scf threshold during optimization

        CHARACTER(len=80) :: pot_extrapolation = 'default', &
                             wfc_extrapolation = 'default'
          !  These variables are used only by PWSCF

        LOGICAL :: refold_pos
        LOGICAL :: remove_rigid_rot = .false.

        !
        ! ... delta_T, nraise, tolp are used to change temperature in PWscf
        !

        REAL(DP) :: delta_t = 1.0_DP
        !! used to change temperature in PWscf
        INTEGER :: nraise = 1
        !! used to change temperature in PWscf
        !
        ! ... variables added for new BFGS algorithm
        !

        INTEGER ::  bfgs_ndim = 1
        LOGICAL ::  tgdiis_step = .TRUE. 

        REAL(DP)  :: trust_radius_max = 0.8_DP
        REAL(DP)  :: trust_radius_min = 1.E-3_DP
        REAL(DP)  :: trust_radius_ini = 0.5_DP

        REAL(DP)  :: w_1 = 0.5E-1_DP
        REAL(DP)  :: w_2 = 0.5_DP

        !
        ! Parameters for minimization with the FIRE algorithm   
        !
        INTEGER  :: fire_nmin = 5 ! minimum number of steps for time step increase 
        REAL(DP) :: fire_f_inc = 1.1_DP ! factor for time step increase  
        REAL(DP) :: fire_f_dec = 0.5_DP ! factor for time step decrease
        REAL(DP) :: fire_alpha_init = 0.2_DP ! initial value of mixing factor
        REAL(DP) :: fire_falpha = 0.99_DP ! modify the mixing factor
        REAL(DP) :: fire_dtmax = 10.0_DP ! maximum time step; calculated as dtmax = fire_dtmax*dt 
        !

        !
        NAMELIST / ions / ion_dynamics, iesr, ion_radius, ion_damping,         &
                          ion_positions, ion_velocities, ion_temperature,      &
                          tempw, fnosep, nhgrp, fnhscl, nhpcl, nhptyp, ndega, tranp,   &
                          amprp, greasp, tolp, ion_nstepe, ion_maxstep,        &
                          refold_pos, upscale, delta_t, pot_extrapolation,     &
                          wfc_extrapolation, nraise, remove_rigid_rot,         &
                          trust_radius_max, trust_radius_min,                  &
                          trust_radius_ini, w_1, w_2, bfgs_ndim,tgdiis_step,   &
                          fire_nmin, fire_f_inc, fire_f_dec, fire_alpha_init,  &
                          fire_falpha, fire_dtmax 



!=----------------------------------------------------------------------------=!
!  CELL Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        CHARACTER(len=80) :: cell_parameters = 'default'
        !! allowed options:  
        !! -'default'    restart the simulation with cell parameters read
        !!               from the restart file or "celldm" if restart='from_scratch';  
        !! -'from_input' restart the simulation with cell parameters
        !!               from standard input (see the card 'CELL_PARAMETERS').

        CHARACTER(len=80) :: cell_dynamics  = 'none'
        !! set how the cell should be moved
        CHARACTER(len=80) :: cell_dynamics_allowed(8)
        DATA cell_dynamics_allowed / 'sd', 'pr', 'none', 'w', 'damp-pr', &
                                     'damp-w', 'bfgs', 'ipi'  /

        CHARACTER(len=80) :: cell_velocities = 'default'
        !! allowed options:  
        !! -'zero'    restart setting cell velocitiy to zero;  
        !! -'default' restart using cell velocity of the previous run

        REAL(DP) :: press = 0.0_DP
        !! external pressure (in GPa, remember \(1\text{kbar}=10^8\text{Pa}\))

        REAL(DP) :: wmass = 0.0_DP
        !! effective cell mass in the Parrinello-Rahman Lagrangian (in atomic units).
        !! Of the order of magnitude of the total atomic mass
        !! (sum of the mass of the atoms) within the simulation cell.
        !! If you do not specify this parameters, the code will compute
        !! its value based on some physical consideration.

        CHARACTER(len=80) :: cell_temperature  = 'not_controlled'
        !! Allowed options:  
        !! -'nose'           control cell temperature using Nose thermostat,
        !!                   see parameters "fnoseh" and "temph";  
        !! -'rescaling'      control cell temperature via velocities rescaling;  
        !! -'not_controlled' cell temperature is not controlled.  
        !! NOT used in FPMD

        REAL(DP) :: temph = 0.0_DP
        !! meaningful only with cell\_temperature different from 'not_controlled'.
        !! Value of the cell temperature (in Kelvin) forced by the temperature control.

        REAL(DP) :: fnoseh = 1.0_DP
        !! meaningful only with cell\_temperature = 'nose'.
        !! Oscillation frequency of the nose thermostat (in terahertz)

        REAL(DP) :: greash = 0.0_DP
        !! same as "grease", for cell damped dynamics

        CHARACTER(len=80) :: cell_dofree = 'all'
        !! selects which of the cell parameters should be moved. Allowed options:  
        !! -'all':    all axis and angles are propagated (default);  
        !! -'volume': the cell is simply rescaled, without changing the shape;  
        !! -'x':      only the "x" axis is moved;  
        !! -'y':      only the "y" axis is moved;  
        !! -'z':      only the "z" axis is moved;  
        !! -'xy':     only the "x" and "y" axis are moved, angles are unchanged;  
        !! -'xz':     only the "x" and "z" axis are moved, angles are unchanged;  
        !! -'yz':     only the "y" and "z" axis are moved, angles are unchanged;  
        !! -'xyz':    "x", "y" and "z" axis are moved, angles are unchanged.

        REAL(DP) :: cell_factor = 0.0_DP
        !! NOT used in FPMD

        INTEGER   :: cell_nstepe = 1
        !! number of electronic steps for each cell step

        REAL(DP) :: cell_damping = 0.1_DP
        !! meaningful only if cell\_dynamics='damp'.
        !! Damping frequency times delta t, optimal values could be
        !! calculated with the formula:  
        !! \( \sqrt(0.5\cdot\log{(E1-E2)/(E2-E3))}\),  
        !! where E1 E2 E3 are successive values of the DFT total energy
        !! in a ionic steepest descent simulation.

        REAL(DP) :: press_conv_thr = 0.5_DP
        
        LOGICAL :: treinit_gvecs  = .FALSE. 
        !! if TRUE all the quantities related to fft g vectors are updated at 
        !! step of variable cell structural optimization 

        NAMELIST / cell / cell_parameters, cell_dynamics, cell_velocities, &
                          press, wmass, cell_temperature, temph, fnoseh,   &
                          cell_dofree, greash, cell_factor, cell_nstepe,   &
                          cell_damping, press_conv_thr, treinit_gvecs 

!
!=----------------------------------------------------------------------------=!!
! PRESS_AI Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
!
      LOGICAL  :: abivol = .false.
      LOGICAL  :: abisur = .false.
      LOGICAL  :: pvar   = .false.
      LOGICAL  :: fill_vac=.false.
      LOGICAL  :: scale_at=.false.
      LOGICAL  :: t_gauss =.false.
      LOGICAL  :: jellium= .false.
      LOGICAL  :: cntr(nsx)=.false.
      REAL(DP) :: P_ext = 0.0_DP
      REAL(DP) :: P_in  = 0.0_DP
      REAL(DP) :: P_fin = 0.0_DP
      REAL(DP) :: rho_thr = 0.0_DP
      REAL(DP) :: step_rad(nsx) = 0.0_DP
      REAL(DP) :: Surf_t = 0.0_DP
      REAL(DP) :: dthr = 0.0_DP
      REAL(DP) :: R_j = 0.0_DP
      REAL(DP) :: h_j = 0.0_DP
      REAL(DP) :: delta_eps = 0.0_DP
      REAL(DP) :: delta_sigma=0.0_DP
      INTEGER  :: n_cntr = 0
      INTEGER  :: axis = 0

      NAMELIST / press_ai / abivol, P_ext, pvar, P_in, P_fin, rho_thr,  &
     &                      step_rad, delta_eps, delta_sigma, n_cntr,   &
     &                      fill_vac, scale_at, t_gauss, abisur,        &
     &                      Surf_t, dthr, cntr, axis, jellium, R_j, h_j

!=----------------------------------------------------------------------------=!
!  WANNIER Namelist Input Parameters
!=----------------------------------------------------------------------------=!

          LOGICAL :: wf_efield
          LOGICAL :: wf_switch
          !
          INTEGER :: sw_len
          !
          REAL(DP) :: efx0, efy0, efz0
          REAL(DP) :: efx1, efy1, efz1
          !
          INTEGER  :: wfsd
          !
          REAL(DP) :: wfdt
          REAL(DP) :: maxwfdt
          REAL(DP) :: wf_q
          REAL(DP) :: wf_friction
!=======================================================================
!exx_wf related
          INTEGER  :: vnbsp
          INTEGER  :: exx_neigh
          REAL(DP) :: exx_poisson_eps
          REAL(DP) :: exx_dis_cutoff
          REAL(DP) :: exx_ps_rcut_self
          REAL(DP) :: exx_ps_rcut_pair
          REAL(DP) :: exx_me_rcut_self
          REAL(DP) :: exx_me_rcut_pair
          LOGICAL  :: exx_use_cube_domain
!=======================================================================

          INTEGER :: nit
          INTEGER :: nsd
          INTEGER :: nsteps
          !
          REAL(DP) :: tolw
          !
          LOGICAL :: adapt
          !
          INTEGER :: calwf
          INTEGER :: nwf
          INTEGER :: wffort
          !
          LOGICAL :: writev
!==============================================================================
!exx_wf related
          NAMELIST / wannier / wf_efield, wf_switch, sw_len, efx0, efy0, efz0,&
                               efx1, efy1, efz1, wfsd, wfdt,exx_neigh,exx_poisson_eps,&
                               exx_dis_cutoff,exx_ps_rcut_self, exx_me_rcut_self,   &
                               exx_use_cube_domain,                      &
                               exx_ps_rcut_pair, exx_me_rcut_pair, vnbsp,&
                               maxwfdt, wf_q, wf_friction, nit, nsd, nsteps,  & 
                               tolw, adapt, calwf, nwf, wffort, writev
!===============================================================================
!  END manual
! ----------------------------------------------------------------------

!=----------------------------------------------------------------------------=!
!  WANNIER_NEW Namelist Input Parameters
!=----------------------------------------------------------------------------=!

          LOGICAL :: plot_wannier = .false.
          !! if TRUE wannier number plot_wan\_num is plotted
          LOGICAL :: use_energy_int = .false.
          !! if TRUE energy interval is used to generate wannier
          LOGICAL :: print_wannier_coeff = .false.
          ! if .TRUE.
          INTEGER, PARAMETER :: nwanx = 50
          !! max number of wannier functions
          INTEGER :: nwan
          !! number of wannier functions
          INTEGER :: plot_wan_num = 0
          !! number of wannier for plotting
          INTEGER :: plot_wan_spin = 1
          !! spin of wannier for plotting
          REAL(DP) :: constrain_pot(nwanx,2)
          !! constrained potential for wannier
          NAMELIST / wannier_ac / plot_wannier, use_energy_int, nwan, &
                                   plot_wan_num, plot_wan_spin, constrain_pot, print_wannier_coeff

!  END manual
! ----------------------------------------------------------------------
!
!=----------------------------------------------------------------------------=!
!  FCP Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        REAL(DP) :: fcp_mu = 0.0_DP
        !! target Fermi energy (in eV)

        CHARACTER(LEN=16) :: fcp_dynamics = 'none'
        !! available options:  
        !! -'none':            Not specified;  
        !! -'lm':              Line-Minimization;  
        !! -'newton':          Newton-Raphson algorithm (with DIIS);  
        !! -'bfgs':            BFGS algorithm (coupling with ions);  
        !! -'damp':            Damped dynamics (quick-min Verlet);  
        !! -'verlet':          Verlet dynamics;  
        !! -'velocity-verlet': Velocity-Verlet dynamics.

        CHARACTER(LEN=16) :: fcp_dynamics_allowed(7)
        DATA fcp_dynamics_allowed / 'none', 'lm', 'newton', 'bfgs', &
                                    'damp', 'verlet', 'velocity-verlet' /

        REAL(DP) :: fcp_conv_thr = 1.0E-2_DP
        !! convergence threshold for FCP relaxation (in eV)

        INTEGER :: fcp_ndiis = 4
        !! size of DIIS for Newton-Raphson algorithm

        REAL(DP) :: fcp_rdiis = 1.0_DP
        !! step of DIIS for Newton-Raphson algorithm

        REAL(DP) :: fcp_mass = -1.0_DP
        !! mass for the FCP

        REAL(DP) :: fcp_velocity = 0.0_DP
        !! initial velocity for the FCP

        CHARACTER(LEN=80) :: fcp_temperature = 'not_controlled'
        !! Allowed options:  
        !! -'rescaling':      control FCP's temperature via velocity rescaling,
        !!                    see parameters "fcp\_tempw" and "fcp\_tolp";  
        !! -'rescale-v':      control FCP's temperature via velocity rescaling,
        !!                    see parameters "fcp\_tempw" and "fcp\_nraise";  
        !! -'rescale-T':      control FCP's temperature via velocity rescaling,
        !!                    see parameter "fcp\_delta\_t";  
        !! -'reduce-T':       reduce FCP's temperature,
        !!                    see parameters "fcp\_nraise", "fcp\_delta\_t";  
        !! -'berendsen':      control FCP's temperature using "soft" velocity
        !!                    rescaling - see parameters "fcp\_tempw" and "fcp\_nraise"
        !! -'andersen':       control FCP's temperature using Andersen thermostat,
        !!                    see parameters "fcp\_tempw" and "fcp\_nraise";  
        !! -'initial':        initialize ion velocities to temperature fcp\_tempw
        !!                    and leave uncontrolled further on;  
        !! -'not_controlled': FCP's temperature is not controlled.

        REAL(DP) :: fcp_tempw = 300.0_DP
        !! meaningful only with fcp\_temperature different from 'not\_controlled'.
        !! Value of the FCP's temperature (in Kelvin) forced by the temperature control

        REAL(DP) :: fcp_tolp = 100.0_DP
        !! parameter to control temperature

        REAL(DP) :: fcp_delta_t = 1.0_DP
        !! parameter to control temperature

        INTEGER :: fcp_nraise = 1
        !! parameter to control temperature

        LOGICAL :: freeze_all_atoms = .FALSE.
        !! freeze (or fix) all atoms.
        !! To perform relaxation or dynamics only with FCP.

        NAMELIST / fcp / fcp_mu, fcp_dynamics, fcp_conv_thr, fcp_ndiis, fcp_rdiis, &
                         fcp_mass, fcp_velocity, fcp_temperature, &
                         fcp_tempw, fcp_tolp, fcp_delta_t, fcp_nraise, &
                         freeze_all_atoms

!  END manual
! ----------------------------------------------------------------------
!
!=----------------------------------------------------------------------------=!
!  RISM Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        INTEGER :: nsolv = 0
          ! number of solvents

        CHARACTER(len=80) :: closure = 'kh'
          ! closure = 'hnc' | 'kh'*
          ! select type of closure equation
          ! 'hnc'  HyperNetted-Chain model
          ! 'kh'   Kovalenko and Hirata's model
        CHARACTER(len=80) :: closure_allowed(2)
        DATA closure_allowed / 'hnc', 'kh' /

        REAL(DP) :: tempv = 300.0_DP
          ! value of the solvent temperature (in Kelvin)
          ! during 1D- and 3D-RISM calculations.

        REAL(DP) :: ecutsolv = 0.0_DP
          ! energy cutoff for 3D-RISM in k-space (in Rydberg)
          ! by default its value is "4 * ecutwfc"

        CHARACTER(len=80) :: solute_lj(nsx) = 'uff'
          ! solute_lj = 'none' | 'uff'* | 'clayff' | 'opls-aa'
          ! select type of Lennard-Jones force fields for solutes
          ! 'uff'      Universal Force Field
          ! 'clayff'   Clay's Force Field
          ! 'opls-aa'  OPLS-AA (generic parameters for QM/MM)
        CHARACTER(len=80) :: solute_lj_allowed(4)
        DATA solute_lj_allowed / 'none', 'uff', 'clayff', 'opls-aa' /

        REAL(DP) :: solute_epsilon(nsx) = -1.0_DP
          ! Lennard-Jones parameters `epsilon' for solutes (in kcal/mol)

        REAL(DP) :: solute_sigma(nsx) = -1.0_DP
          ! Lennard-Jones parameters `sigma' for solutes (in angstrom)

        REAL(DP) :: rmax_lj = 5.0_DP
          ! maximum radius of Lennard-Jones for 3D-RISM (in sigma)

        REAL(DP) :: rmax1d = 1000.0_DP
          ! maximum inter-site radius for 1D-RISM (in bohr)

        CHARACTER(len=80) :: starting1d = 'zero'
          ! starting1d = 'zero'* | 'file' | 'fix'
          ! define how the code should initialize the 1D-RISM's correlation function
          ! 'zero'  start from 0
          ! 'file'  read from file
          ! 'fix'   read from file, and fix correlation function
        CHARACTER(len=80) :: starting1d_allowed(3)
        DATA starting1d_allowed / 'zero', 'file', 'fix' /

        CHARACTER(len=80) :: starting3d = 'zero'
          ! starting3d = 'zero'* | 'file'
          ! define how the code should initialize the 3D-RISM's correlation function
          ! 'zero'  start from 0
          ! 'file'  read from file
        CHARACTER(len=80) :: starting3d_allowed(2)
        DATA starting3d_allowed / 'zero', 'file' /

        REAL(DP) :: smear1d = 2.0_DP
          ! smearing radius for 1D-RISM (in bohr)

        REAL(DP) :: smear3d = 2.0_DP
          ! smearing radius for 3D-RISM (in bohr)

        INTEGER :: rism1d_maxstep = 50000
          ! maximum number of steps in 1D-RISM calculation

        INTEGER :: rism3d_maxstep = 5000
          ! maximum number of steps in 3D-RISM calculation

        REAL(DP) :: rism1d_conv_thr = 1.0E-8_DP
          ! convergence threshold for 1D-RISM calculation
          ! convergence is achieved when RMS of residual vector < rism1d_conv_thr

        REAL(DP) :: rism3d_conv_thr = 1.0E-5_DP
          ! convergence threshold for 3D-RISM calculation
          ! convergence is achieved when RMS of residual vector < rism3d_conv_thr

        INTEGER :: mdiis1d_size = 20
          ! size of MDIIS algorithm in 1D-RISM calculation

        INTEGER :: mdiis3d_size = 10
          ! size of MDIIS algorithm in 3D-RISM calculation

        REAL(DP) :: mdiis1d_step = -1.0_DP
          ! step width of MDIIS algorithm in 1D-RISM calculation

        REAL(DP) :: mdiis3d_step = -1.0_DP
          ! step width of MDIIS algorithm in 3D-RISM calculation

        REAL(DP) :: rism1d_bond_width = 0.0_DP
          ! gaussian width of bonds in 1D-RISM calculation

        REAL(DP) :: rism1d_dielectric = -1.0_DP
          ! dielectric constant for DRISM

        REAL(DP) :: rism1d_molesize = 2.0_DP
          ! size of solvent molecule for DRISM (in bohr)

        INTEGER :: rism1d_nproc = 128
          ! number of processes to calculate 1D-RISM

        INTEGER :: rism1d_nproc_switch = 16
          ! number of processes to calculate 1D-RISM

        REAL(DP) :: rism3d_conv_level = -1.0_DP
          ! convergence level of 3D-RISM

        LOGICAL :: rism3d_planar_average = .FALSE.
          ! calculate planar average of solvents after 3D-RISM calculation, or not

        INTEGER :: laue_nfit = 4
          ! number of fitting points in Laue-RISM calculation

        REAL(DP) :: laue_expand_right = -1.0_DP
          ! expanding length on right-hand side in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_expand_left = -1.0_DP
          ! expanding length on left-hand side in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_starting_right = 0.0_DP
          ! starting position on right-hand side in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_starting_left = 0.0_DP
          ! starting position on left-hand side in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_buffer_right = -1.0_DP
          ! buffering length on right-hand side in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_buffer_right_solu = -1.0_DP
          ! additional buffering length on right-hand side
          ! of solute-ward in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_buffer_right_solv = -1.0_DP
          ! additional buffering length on right-hand side
          ! of solvent-ward in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_buffer_left = -1.0_DP
          ! buffering length on left-hand side in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_buffer_left_solu = -1.0_DP
          ! additional buffering length on left-hand side
          ! of solute-ward in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_buffer_left_solv = -1.0_DP
          ! additional buffering length on left-hand side
          ! of solvent-ward in Laue-RISM calculation (in bohr)

        LOGICAL :: laue_both_hands = .FALSE.
          ! use both-hands method in Laue-RISM calculation, or not

        CHARACTER(len=80) :: laue_reference = 'none'
          ! laue_reference = 'none'* | 'average' | 'right' | 'left'
          ! reference of electrostatic potential in Laue-RISM calculation
          ! (used to evaluate Fermi energy and to calculate FCP)
          ! 'none'     explicit reference is not defined
          ! 'average'  average of right-hand side and left-hand side
          ! 'right'    right-hand side
          ! 'left'     left-hand side
        CHARACTER(len=80) :: laue_reference_allowed(4)
        DATA laue_reference_allowed / 'none', 'average', 'right', 'left' /

        CHARACTER(len=80) :: laue_wall = 'auto'
          ! laue_wall = 'none' | 'auto'* | 'manual'
          ! define repulsive wall in Laue-RISM calculation
          ! 'none'    wall is not defined
          ! 'auto'    edge position of wall is defined automatically
          ! 'manual'  edge position of wall is defined manually

        CHARACTER(len=80) :: laue_wall_allowed(3)
        DATA laue_wall_allowed / 'none', 'manual', 'auto' /

        REAL(DP) :: laue_wall_z = 0.0_DP
          ! edge position of repulsive wall in Laue-RISM calculation (in bohr)

        REAL(DP) :: laue_wall_rho = 0.01_DP
          ! density of repulsive wall in Laue-RISM calculation (in 1/bohr^3)

        REAL(DP) :: laue_wall_epsilon = 0.1_DP
          ! Lennard-Jones parameters `epsilon' for repulsive wall
          ! in Laue-RISM calculation (in kcal/mol)

        REAL(DP) :: laue_wall_sigma = 4.0_DP
          ! Lennard-Jones parameters `sigma' for repulsive wall
          ! in Laue-RISM calculation (in angstrom)

        LOGICAL :: laue_wall_lj6 = .FALSE.
          ! use attractive term of Lennard-Jones: -(1/r)^6, or not

        NAMELIST / rism / nsolv, closure, tempv, ecutsolv, solute_lj, &
                          solute_epsilon, solute_sigma, rmax_lj, rmax1d, &
                          starting1d, starting3d, smear1d, smear3d, &
                          rism1d_maxstep, rism3d_maxstep, rism1d_conv_thr, rism3d_conv_thr, &
                          mdiis1d_size, mdiis3d_size, mdiis1d_step, mdiis3d_step, &
                          rism1d_bond_width, rism1d_dielectric, rism1d_molesize, &
                          rism1d_nproc, rism1d_nproc_switch, &
                          rism3d_conv_level, rism3d_planar_average, &
                          laue_nfit, laue_expand_right, laue_expand_left, &
                          laue_starting_right, laue_starting_left, &
                          laue_buffer_right, laue_buffer_right_solu, laue_buffer_right_solv, &
                          laue_buffer_left, laue_buffer_left_solu, laue_buffer_left_solv, &
                          laue_both_hands, laue_reference, laue_wall, laue_wall_z, laue_wall_rho, &
                          laue_wall_epsilon, laue_wall_sigma, laue_wall_lj6
!  END manual
! ----------------------------------------------------------------------


! ----------------------------------------------------------------
! BEGIN manual
!
!=----------------------------------------------------------------------------=!
!  CARDS parameters
!=----------------------------------------------------------------------------=!
!
!  Note: See file read_cards.f90 for card syntax and usage
!
!    ATOMIC_SPECIES
!
        CHARACTER(len=6)  :: atom_label(nsx) = 'XX'
        !! label of the atomic species being read
        CHARACTER(len=80) :: atom_pfile(nsx) = 'YY'
        !! pseudopotential file name
        REAL(DP)          :: atom_mass(nsx)  = 0.0_DP
        !! atomic mass of the i-th atomic species in atomic mass units:  
        !! \(1\text{a.m.u.}=1822.9\text{a.u.}=1.6605\cdot 10^{-27}kg\)
        LOGICAL   :: taspc = .false.
        LOGICAL   :: tkpoints = .false.
        LOGICAL   :: tforces = .false.
        LOGICAL   :: tocc = .false.
        LOGICAL   :: tcell = .false.
        LOGICAL   :: tionvel = .false.
        LOGICAL   :: tconstr = .false.
        LOGICAL   :: tksout = .false.
        LOGICAL   :: ttemplate = .false.
        LOGICAL   :: twannier = .false.
        LOGICAL   :: tsolvents = .false.
        LOGICAL   :: ttotcharge = .false.

!
!    ATOMIC_POSITIONS
!
        REAL(DP), ALLOCATABLE :: rd_pos(:,:)
        !! unsorted positions from input
        INTEGER,  ALLOCATABLE :: sp_pos(:)
        INTEGER,  ALLOCATABLE :: rd_if_pos(:,:)
        INTEGER,  ALLOCATABLE :: na_inp(:)
        LOGICAL  :: tapos = .false.
        LOGICAL  :: lsg   = .false.
        CHARACTER(len=80) :: atomic_positions = 'crystal'
        !! atomic\_positions = 'bohr' | 'angstrom' | 'crystal' | 'alat'  
        !! Select the units for the atomic positions being read from stdin

!
!    ION_VELOCITIES
!
        REAL(DP), ALLOCATABLE :: rd_vel(:,:)
        !! unsorted velocities from input
        INTEGER,  ALLOCATABLE :: sp_vel(:)
        LOGICAL  :: tavel          = .false.
!
!    ATOMIC_FORCES
!
        REAL(DP), ALLOCATABLE :: rd_for(:,:)
        !! external forces applied to single atoms

!
!    KPOINTS
!
! ...   k-points inputs
        LOGICAL :: tk_inp = .false.
        REAL(DP), ALLOCATABLE :: xk(:,:), wk(:)
        CHARACTER(len=50), ALLOCATABLE :: labelk(:)
        INTEGER :: nkstot = 0, nk1 = 0, nk2 = 0, nk3 = 0, k1 = 0, k2 = 0, k3 = 0
        CHARACTER(len=80) :: k_points = 'gamma'
        !! select the k points mesh. Available options:  
        !! 'automatic':  k points mesh is generated automatically
        !!               with Monkhorst-Pack algorithm;  
        !! 'crystal':    k points mesh is given in stdin in scaled units;  
        !! 'tpiba':      k points mesh is given in stdin in units of (2PI/alat);  
        !! 'gamma':      only gamma point is used (default in CPMD simulation);  
        !! 'crystal\_b', 'tpiba\_b': postfix '\_b' means that a band input is given.  
        !! The weights is a integer number that gives the number of points between
        !! the present point and the next. The weight of the last point is not used.
!
!    OCCUPATIONS
!
        REAL(DP), ALLOCATABLE :: f_inp(:,:)
        LOGICAL   :: tf_inp = .false.

!
!    CELL_PARAMETERS
!
       REAL(DP) :: rd_ht(3,3) = 0.0_DP
       CHARACTER(len=80) :: cell_units = 'none'
       LOGICAL   :: trd_ht = .false.

!
!    REFERENCE_CELL_PARAMETERS
!
       REAL(DP) :: rd_ref_ht(3,3) = 0.0_DP
       CHARACTER(len=80) :: ref_cell_units = 'alat'
       LOGICAL   :: ref_cell = .false.

!
!    CONSTRAINTS
!
      INTEGER :: nc_fields = 4
      !! max number of fields that is allowed to
      !! define a constraint

      INTEGER  :: nconstr_inp    = 0
      REAL(DP) :: constr_tol_inp = 1.E-6_DP
      !
      CHARACTER(len=20), ALLOCATABLE :: constr_type_inp(:)
      REAL(DP),          ALLOCATABLE :: constr_inp(:,:)
      REAL(DP),          ALLOCATABLE :: constr_target_inp(:)
      LOGICAL,           ALLOCATABLE :: constr_target_set(:)

!
!    KOHN_SHAM
!
      INTEGER, ALLOCATABLE :: iprnks( :, : )
      INTEGER :: nprnks( nspinx ) = 0
      !! logical mask used to specify which kohn sham orbital should be
      !! written to files 'KS.'

!
!   PLOT_WANNIER
!

      INTEGER, PARAMETER :: nwf_max = 1000
      !
      INTEGER :: wannier_index( nwf_max )

!
!   WANNIER_NEW
!
      TYPE (wannier_data) :: wan_data(nwanx,2)

!
!    SOLVENTS
!
      CHARACTER(len=10) :: solv_label(nsolx) = 'XX'    
      !! label of the solvents
      CHARACTER(len=80) :: solv_mfile(nsolx) = 'YY'    
      !! molecular file name
      REAL(DP)          :: solv_dens1(nsolx) = 0.0_DP  
      !! solvent's density (for the right-hand side)
      REAL(DP)          :: solv_dens2(nsolx) = 0.0_DP  
      !! solvent's density (for the left-hand side)
      CHARACTER(len=80) :: solvents_unit = '1/cell'
      !! solvents_unit = '1/cell' | 'mol/L' | 'g/cm^3'
      !! select the units for the solvent's densities being read from stdin
!   HUBBARD
!
      LOGICAL  :: tahub = .false.

!  END manual
! ----------------------------------------------------------------------

      LOGICAL :: xmloutput = .false.
      !! if TRUE PW produce an xml output

CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE reset_input_checks()
  !-----------------------------------------------------------------------------
  !! This routine sets to FALSE flags used to check whether some variables
  !! have been read. If called before reading, allows to read a different
  !! input file without triggering bogus error messages - useful for NEB.
  !
  IMPLICIT NONE
  !
  tapos = .false.
  tkpoints = .false.
  taspc = .false.
  twannier = .false.
  tconstr = .false.
  tforces = .false.
  tocc = .false.
  tksout = .false.
  tionvel = .false.
  tcell = .false.
  tsolvents = .false.
  ttotcharge = .false.
  !
  END SUBROUTINE reset_input_checks
  !
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE allocate_input_ions( ntyp, nat )
  !-----------------------------------------------------------------------------
    !
    INTEGER, INTENT(in) :: ntyp, nat
    !
    IF ( allocated( rd_pos ) ) DEALLOCATE( rd_pos )
    IF ( allocated( sp_pos ) ) DEALLOCATE( sp_pos )
    IF ( allocated( rd_if_pos ) ) DEALLOCATE( rd_if_pos )
    IF ( allocated( na_inp ) ) DEALLOCATE( na_inp )
    IF ( allocated( rd_vel ) ) DEALLOCATE( rd_vel )
    IF ( allocated( sp_vel ) ) DEALLOCATE( sp_vel )
    IF ( allocated( rd_for ) ) DEALLOCATE( rd_for )
    !
    ALLOCATE( rd_pos( 3, nat ) )
    ALLOCATE( sp_pos( nat)   )
    ALLOCATE( rd_if_pos( 3, nat ) )
    ALLOCATE( na_inp( ntyp)  )
    ALLOCATE( rd_vel( 3, nat ) )
    ALLOCATE( sp_vel( nat)   )
    ALLOCATE( rd_for( 3, nat ) )
    !
    rd_pos = 0.0_DP
    sp_pos = 0
    rd_if_pos = 1
    na_inp = 0
    rd_vel = 0.0_DP
    sp_vel = 0
    rd_for = 0.0_DP
    !
    RETURN
    !
  END SUBROUTINE allocate_input_ions

  !-----------------------------------------------------------------------------
  SUBROUTINE allocate_input_constr()
  !-----------------------------------------------------------------------------
    !
    IF ( allocated( constr_type_inp ) )   DEALLOCATE( constr_type_inp )
    IF ( allocated( constr_inp ) )        DEALLOCATE( constr_inp )
    IF ( allocated( constr_target_inp ) ) DEALLOCATE( constr_target_inp )
    IF ( allocated( constr_target_set ) ) DEALLOCATE( constr_target_set )
    !
    ALLOCATE( constr_type_inp(   nconstr_inp ) )
    ALLOCATE( constr_target_inp( nconstr_inp ) )
    ALLOCATE( constr_target_set( nconstr_inp ) )
    !
    ALLOCATE( constr_inp( nc_fields, nconstr_inp ) )
    !
    constr_type_inp   = ' '
    constr_inp        = 0.0_DP
    constr_target_inp = 0.0_DP
    constr_target_set = .false.
    !
    RETURN
    !
  END SUBROUTINE allocate_input_constr

  !-----------------------------------------------------------------------------
  SUBROUTINE allocate_input_iprnks( nksx, nspin )
  !-----------------------------------------------------------------------------
    !
    INTEGER, INTENT(in) :: nksx, nspin
    !
    IF( allocated( iprnks ) ) DEALLOCATE( iprnks )
    !
    ALLOCATE( iprnks( max( 1, nksx), nspin ) )
    !
    iprnks = 0
    !
    RETURN
    !
  END SUBROUTINE allocate_input_iprnks

  !-----------------------------------------------------------------------------
  SUBROUTINE deallocate_input_parameters()
  !-----------------------------------------------------------------------------
    !
    IF ( allocated( xk ) ) DEALLOCATE( xk )
    IF ( allocated( wk ) ) DEALLOCATE( wk )
    IF ( allocated( labelk ) ) DEALLOCATE( labelk )
    IF ( allocated( rd_pos ) ) DEALLOCATE( rd_pos )
    IF ( allocated( sp_pos ) ) DEALLOCATE( sp_pos )
    IF ( allocated( rd_if_pos ) ) DEALLOCATE( rd_if_pos )
    IF ( allocated( na_inp ) ) DEALLOCATE( na_inp )
    IF ( allocated( rd_vel ) ) DEALLOCATE( rd_vel )
    IF ( allocated( sp_vel ) ) DEALLOCATE( sp_vel )
    IF ( allocated( rd_for ) ) DEALLOCATE( rd_for )
    !
    IF ( allocated( f_inp ) ) DEALLOCATE( f_inp )
    !
    IF ( allocated( constr_type_inp ) )   DEALLOCATE( constr_type_inp )
    IF ( allocated( constr_inp ) )        DEALLOCATE( constr_inp )
    IF ( allocated( constr_target_inp ) ) DEALLOCATE( constr_target_inp )
    IF ( allocated( constr_target_set ) ) DEALLOCATE( constr_target_set )
    !
    IF ( allocated( iprnks ) )       DEALLOCATE( iprnks )
    !
    RETURN
    !
  END SUBROUTINE deallocate_input_parameters
  !
!=----------------------------------------------------------------------------=!
!
END MODULE input_parameters
!
!=----------------------------------------------------------------------------=!
