!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------------=!
!
MODULE input_parameters
!
!=----------------------------------------------------------------------------=!
!
!  this module contains
!  1) the definitions of all input parameters
!     (both those read from namelists and those read from cards)
!  2) the definitions of all namelists
!  3) routines that allocate data needed in input
!  Note that all values are initialized, but the default values should be
!  set in the appropriate routines contained in module "read_namelists"
!  Documentation of input variabe is in INPUT_PW
!  Originally written by Carlo Cavazzoni for FPMD
!
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nsx, npkx, nspinx, lqmax, nhclm, max_nconstr
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
!     &PHONON
!       phonon_parameter_1,
!       phonon_parameter_Lastone
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
!       ( when appropriate default values are marked with a "*" )
!

        CHARACTER(LEN=80) :: title = ' '
          ! title = 'title for this simulation'

        CHARACTER(LEN=80) :: calculation = 'none'
          ! Specify the type of the simulation
          ! 'scf'      = electron minimization/selfconsistency
          ! 'nscf'     = nonselfconsistent calculation of electron states
          ! 'bands'    = as above, band structure calculation only
          ! 'phonon'   = as above, plus data needed for a phonon calculation
          ! 'relax'    = ionic minimization
          ! 'cp'       = Car-Parrinello MD
          ! 'md'       = Car-Parrinello MD
          ! 'vc-relax' = ionic + cell minimization
          ! 'vc-cp'    = variable-cell Car-Parrinello (-Rahman) dynamics
          ! 'vc-md'    = variable-cell Car-Parrinello (-Rahman) dynamics
          ! 'neb'      = NEB    Method search of the Minimum Energy Path (MEP)
          ! 'smd'      = String Method search of the Minimum Energy Path (MEP)
          ! 'cp-wf'    = Car-Parrinello with wannier functions
          ! 'fpmd'     = Compatibility with the old FPMD code
          ! 'fpmd-neb' = NEB using fpmd as scf engine
          ! 'metadyn'  = meta-dynamics (Laio-Parrinello dynamics)

        CHARACTER(LEN=80) :: calculation_allowed(16)
        DATA calculation_allowed / 'scf', 'nscf', 'relax', 'md', 'cp', &
          'vc-relax', 'vc-md', 'vc-cp', 'phonon', 'bands', 'neb', 'smd', &
          'cp-wf', 'fpmd', 'metadyn', 'fpmd-neb' /
          ! Allowed value for calculation parameters


        CHARACTER(LEN=80) :: verbosity = 'default'
          ! verbosity = 'high' | 'default'* | 'low' | 'minimal' | 'default+projwfc'
          ! define the verbosity of the code output
        CHARACTER(LEN=80) :: verbosity_allowed(5)
        DATA verbosity_allowed / 'high' , 'default' , 'low' , 'minimal' , 'default+projwfc' /

        CHARACTER(LEN=80) :: restart_mode = 'restart'
          ! restart_mode = 'from_scratch' | 'restart'* | 'reset_counters'
          ! specify how to start/restart the simulation
          !   'from_scratch'    start a new simulation from scratch,
          !                     with random wave functions
          !                     for path only : a simulation is started
          !                     from scratch and the initial guess for the
          !                     path is the linear interpolation between
          !                     the initial and the final images
          !   'restart'         continue a previous simulation,
          !                     and performs  "nstep" new steps,
          !   'reset_counters'  continue a previous simulation,
          !                     performs  "nstep" new steps, resetting
          !                     the counter and averages
        CHARACTER(LEN=80) :: restart_mode_allowed(3)
        DATA restart_mode_allowed / 'from_scratch', 'restart', 'reset_counters' /

        INTEGER :: nstep = 10
          ! number of simulation steps, see "restart_mode"

        INTEGER :: iprint = 10
          ! number of steps between successive writings of relevant physical
          ! quantities to standard output and to files "fort.3?" or "prefix.???"
          ! depending on "prefix" parameter.
          ! In PW iprint is compared also with convergence iterations, and
          ! physical quantities are written every iprint iterations

        INTEGER :: isave = 100
          ! number of steps between successive savings of
          ! information needed to restart the run (see "ndr", "ndw")
          ! relevant only for CP and FPMD

        LOGICAL :: tstress = .TRUE.
          ! This flag controls the printing of the stress, its value is overwritten
          ! and set to ".TRUE." when the cell is moving
          !  .TRUE.  write the stress tensor to standard output every "iprint" steps
          !  .FALSE. do not write the stress tensor stdout

        LOGICAL :: tprnfor = .TRUE.
          ! This flag controls the printing of the interatomic forces,
          ! its value is overwritten and set to ".TRUE." the ions are moving
          !  .TRUE.  write the atomic forces to standard output every "iprint" steps
          !  .FALSE. do not write atomic forces to stdout

        REAL(DP) :: dt = 1.0_DP
          ! time step of the CP molecular dynamics simulation,
          ! in atomic units ( 1 a.u. of time = 2.4189 * 10^-17 s ),
          ! for non CP calculations, this represents the time advancing parameter.
          ! Note: typical values for CP simulations are between 1 and 10 a.u.
          ! In PW dt is used for Born-Oppenheimer molecular dynamics, and
          ! its value is usually larger than for CP dynamics, since it is related
          ! only to the mass of ions.

        INTEGER :: ndr = 50
          ! Fortran unit from which the code read the restart file
          ! at the beginning of the simulation, its value should be greather than 50
          ! and it is opened in the running directory.

        INTEGER :: ndw = 50
          ! Fortran unit to which the code writes the restart file
          ! at the end of the simulation, its value should be greather than 50
          ! and it is opened in the running directory.

        CHARACTER(LEN=256) :: outdir = './'
          ! specify the directory where the code opens output and restart
          ! files. When possible put this directory in the fastest available
          ! filesystem ( not NFS! )

        CHARACTER(LEN=256) :: prefix = 'prefix'
          ! specify the prefix for the output file, if not specified the
          ! files are opened as standard fortran units.
          ! The prefix does _NOT_ apply to restart file

        CHARACTER(LEN=256) :: pseudo_dir = './'
          ! specify the directory containing the pseudopotentials

        REAL(DP) :: refg = 0.05_DP
          ! Accurancy of the interpolation table, interval between
          ! table values in Rydberg

        CHARACTER(LEN=256) :: wfcdir = 'undefined'
          ! scratch directory that is hopefully local to the node
          ! to store large, usually temporary files.

        REAL(DP) :: max_seconds = 1.0E+7_DP
          ! smoothly terminate program after the specified number of seconds
          ! this parameter is typically used to prevent an hard kill from
          ! the queuing system.

        REAL(DP) :: ekin_conv_thr = 1.0E-5_DP
          ! convergence criterion, minimizing the electrons this criterion is met
          ! when "ekin < ekin_conv_thr"
          ! convergence is achieved when all criteria are met

        REAL(DP) :: etot_conv_thr = 1.0E-4_DP
          ! convergence criterion, minimizing the ions this criterion is met
          ! when "etot(n+1)-etot(n) < etot_conv_thr", where "n" is the step index,
          ! and "etot" the DFT energy
          ! convergence is achieved when all criteria are met

        REAL(DP) :: forc_conv_thr = 1.0E-3_DP
          ! convergence criterion, minimizing the ions this criterion is met
          ! when "MAXVAL(fion) < forc_conv_thr", where fion are the atomic forces
          ! convergence is achieved when all criteria are met

        CHARACTER(LEN=80) :: disk_io = 'default'
          ! disk_io = 'high', 'default', 'low', 'minimal'
          ! Specify the amount of I/O activities

        LOGICAL :: tefield  = .FALSE.
          ! if .TRUE. a finite electric field is added to the local potential
          ! only used in PW

          LOGICAL :: tefield2  = .FALSE.
          ! if .TRUE. a second finite electric field is added to the local potential
          ! only used in PW


        LOGICAL :: dipfield = .FALSE.
          ! if .TRUE. the dipole field is subtracted
          ! only used in PW

        LOGICAL :: lberry = .FALSE.
          ! if .TRUE., calculate polarization

        INTEGER :: gdir = 0
          ! G-vector for polarization calculation ( related to lberry )
          ! only used in PW

        INTEGER :: nppstr = 0
          ! number of k-points (parallel vector) ( related to lberry )
          ! only used in PW
        LOGICAL :: lelfield = .FALSE.
          !if true a static homogeneous electric field is present


        INTEGER  :: nberrycyc = 1
          !number of covergence cycles on electric field


        LOGICAL :: wf_collect = .FALSE.
          ! This flag controls the way wavefunctions are stored to disk:
          !  .TRUE.  collect all wavefunctions and store them in a single
          !          restart file ( outdir/prefix.save )
          !  .FALSE. do not collect wavefunctions, leave them in temporary
          !          local file

        INTEGER :: printwfc=1
          ! if <0 do nothing, if==0 print rho and fort.47, if == nband print band

        LOGICAL :: saverho = .TRUE.
          ! This flag controls the saving of charge density in CP codes:
          !  .TRUE.  save charge density to restart dir
          !  .FALSE. do not save charge density

        LOGICAL :: tabps = .FALSE. ! for ab-initio pressure and/or surface
                                   ! calculations

        LOGICAL :: lkpoint_dir = .TRUE. ! opens a directory for each k point

        NAMELIST / control / title, calculation, verbosity, restart_mode, &
          nstep, iprint, isave, tstress, tprnfor, dt, ndr, ndw, outdir,   &
          prefix, wfcdir, max_seconds, ekin_conv_thr, etot_conv_thr,      &
          forc_conv_thr, pseudo_dir, disk_io, tefield, dipfield, lberry,  &
          gdir, nppstr, wf_collect, printwfc, lelfield, nberrycyc, refg,  &
          tefield2, saverho, tabps, lkpoint_dir

!
!=----------------------------------------------------------------------------=!
!  SYSTEM Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        INTEGER :: ibrav = 14
          ! index the Bravais lattice:
          !    ibrav     structure                     point groups
          !     ---    --------------               -------------------
          !      1      cubic P (sc)              m3m, 432, m3, <4>3m, 23
          !      2      cubic F (fcc)             m3m, 432, m3, <4>3m, 23
          !      3      cubic I (bcc)             m3m, 432, m3, <4>3m, 23
          !      4       Hexagonal P        6, 6mm, 6/m, <6>, 622, 6/mmm, <6>2m
          !      4       Trigonal P                3, 3m, <3>, 32, <3>m
          !      5       Trigonal R                3, 3m, <3>, 32, <3>m
          !      6    Tetragonal P (st)     4, 4mm, 4/m, <4>, 422, 4/mmm, <4>2m
          !      7    Tetragonal I (bct)    4, 4mm, 4/m, <4>, 422, 4/mmm, <4>2m
          !      8     Orthorhombic P                  2mm, 222, mmm
          !     12      Monoclinic P                     2, 2/m, m
          !     14      Triclinic P                        1, <1>
          !
          ! Note: in variable cell CP molecular dynamics, usually one do not want
          !       to put constraints on the cell symmetries, therefore an
          !       ibrav = 14 is used

        REAL(DP) :: celldm(6) = 0.0_DP
          ! dimensions of the cell
          !   celldm(1) = a
          !   celldm(2) = b/a
          !   celldm(3) = c/a
          !   celldm(4) = cos(bc)
          !   celldm(5) = cos(ac)
          !   celldm(6) = cos(ab)
          ! only the celldm's that are meaningful for the Bravais lattice
          ! considered need be specified (others are ignored):
          !   ibrav = 1,2,3 : celldm(1)
          !   ibrav = 4,6,7 : celldm(1,3)
          !   ibrav =  5    : celldm(1,4)
          !   ibrav =  8    : celldm(1,2,3)
          !   ibrav = 12    : celldm(1,2,3,4)
          !   ibrav = 14    : celldm(1,2,3,4,5,6)

        INTEGER :: nat = 0
          ! total number of atoms

        INTEGER :: ntyp = 0
          ! number of atomic species

        INTEGER :: nbnd = 0
          ! number of electronic states, this parameter is MANDATORY in FPMD

        REAL(DP):: nelec = 0.0_DP
          ! number of electrons, this parameter is MANDATORY in FPMD
          ! may be fractionary in PW, but not in CP and FPMD !

        REAL(DP):: tot_charge = 0.0_DP
          ! total system charge

        INTEGER :: multiplicity = 0
          ! spin multiplicity (2s+1), 1 for singlet, 2 for doublet etc.
          ! when multiplicity = 0, it is unspecified

        INTEGER :: tot_magnetization = -1
          ! majority - minority spin.
          ! A value < 0 ==> unspecified
!
! A comment about variables nelup, neldw, multiplicity and tot_magnetization:
! All these variables contain the same information and must be kept harmonized.
! Variables nelup and neldw will be removed in future versions of the code.
! Variables multiplicity and tot_magnetization, though redundent will probably
! coexist since multiplicity is the more natural way (?)for defining the spin
! configuratio in the quantum-chemistry community while tot_magnetization is
! more natural (?) when dealing with extended systems.
!
        REAL(DP) :: ecutwfc = 0.0_DP
          ! energy cutoff for wave functions in k-space ( in Rydbergs )
          ! this parameter is MANDATORY in FPMD

        REAL(DP) :: ecutrho = 0.0_DP
          ! energy cutoff for charge density in k-space ( in Rydbergs )
          ! by default its value is "4 * ecutwfc"

        INTEGER :: nr1 = 0
        INTEGER :: nr2 = 0
        INTEGER :: nr3 = 0
          ! dimensions of the real space grid for charge and potentials
          ! presently NOT used in FPMD-N

        INTEGER :: nr1s = 0
        INTEGER :: nr2s = 0
        INTEGER :: nr3s = 0
          ! dimensions of the real space grid for wavefunctions
          ! presently NOT used in FPMD-N

        INTEGER :: nr1b = 0
        INTEGER :: nr2b = 0
        INTEGER :: nr3b = 0
          ! dimensions of the "box" grid for Ultrasoft pseudopotentials
          ! presently NOT used in FPMD-N

        CHARACTER(LEN=80) :: occupations = 'fixed'
          ! occupations = 'smearing' | 'tetrahedra' | 'fixed'* | 'from_input'
          ! select the electronic states filling mode
          ! 'smearing'    smearing function is used around Fermi Level
          !               (see "ngauss" and "dgauss")
          !               NOT used in FPMD-N
          ! 'tetrahedra'  tetrahedron method
          !               NOT used in FPMD-N
          ! 'fixed'       fixed occupations automatically calculated
          ! 'from_input'  fixed occupations given in the input
          !               ( see card 'OCCUPATIONS' )

        CHARACTER(LEN=80) :: smearing = 'gaussian'
          ! smearing = 'gaussian', 'methfessel-paxton', 'marzari-vanderbilt', 'fermi-dirac'

        REAL(DP) :: degauss = 0.0_DP
          ! parameter for the smearing functions
          ! NOT used in FPMD-N

        INTEGER :: nspin = 1
          ! number of spinors
          ! "nspin = 1" for LDA simulations
          ! "nspin = 2" for LSD simulations
          ! "nspin = 4" for NON COLLINEAR simulations

        REAL(DP) :: nelup = 0.0_DP, neldw = 0.0_DP
          ! meaningful only if "nspin = 2",
          ! "nelup" is the number of electrons with spin up
          ! "neldw" is the number of electrons with spin down
          ! Remember the following relation hold "nelec = nelup + neldw"

        LOGICAL :: nosym = .TRUE., noinv = .FALSE.
          ! (do not) use symmetry, q => -q symmetry in k-point generation

        REAL(DP) :: ecfixed = 0.0_DP, qcutz = 0.0_DP, q2sigma = 0.0_DP
          ! parameters for constant cut-off simulations
          ! "ecfixed" is the value (in Rydbergs) of the constant-cutoff
          ! "qcutz" and "q2sigma" are the height and the width (in Rydbergs)
          !   of the energy step for reciprocal vector whose square modulus
          !   is greater than  "ecfixed"

        CHARACTER(LEN=80) :: input_dft = 'none'
          ! Variable used to overwrite dft definition contained in
          ! pseudopotential files .
          ! Default value is 'none' meaning that DFT is take from pseudos.
          ! Allowed values: any legal DFT value as defined in pseudopotentials.
          !

        CHARACTER(LEN=80) :: xc_type = 'none'
          ! xc_type = 'BLYP' | 'BP' | 'PBE' | 'PZ' | 'PW' | 'LDA'
          ! select the exchange and correlation functionals
          ! 'BLYP'  use Becke-Lee-Yang-Parr GCC-XC Functional
          ! 'BP'    use Becke-Perdew GCC-XC Functionals
          ! 'PBE'   use Perdew-Burke-Ernzerhof GCC-XC Functionals
          ! 'PZ'    use Slater X, and Perdew-Zunger C Functionals
          ! 'PW'    use Slater X, and Perdew-Wang C Functionals
          ! 'LDA'   use LDA xc functional: the xc potential is
          !         calculated through an interpolation table

        REAL(DP) :: starting_magnetization( nsx ) = 0.0_DP
          ! ONLY PWSCF

        LOGICAL :: lda_plus_u = .FALSE.
          ! ONLY PWSCF

        REAL(DP) :: starting_ns_eigenvalue(lqmax,nspinx,nsx) = -1.0_DP
          ! ONLY PWSCF

        REAL(DP) :: hubbard_u(nsx) = 0.0_DP
          ! ONLY PWSCF

        REAL(DP) :: hubbard_alpha(nsx) = 0.0_DP
          ! ONLY PWSCF

        CHARACTER(LEN=80) :: U_projection_type = 'atomic'
          ! ONLY PWSCF
        LOGICAL :: la2F = .FALSE.
          !
#if defined (EXX)
        LOGICAL :: x_gamma_extrapolation = .true.
          ! ONLY PWSCF
        INTEGER :: nqx1 = 1, nqx2 = 1, nqx3=1
          ! ONLY PWSCF
!        REAL(DP) :: yukawa = 0.0_DP
          ! ONLY PWSCF
#endif

        REAL(DP) :: a = 0.0_DP

        REAL(DP) :: c = 0.0_DP

        REAL(DP) :: b = 0.0_DP

        REAL(DP) :: cosab = 0.0_DP

        REAL(DP) :: cosac = 0.0_DP

        REAL(DP) :: cosbc = 0.0_DP

        INTEGER   :: edir = 0

        REAL(DP) :: emaxpos = 0.0_DP

        REAL(DP) :: eopreg = 0.0_DP

        REAL(DP) :: eamp = 0.0_DP

        LOGICAL :: noncolin = .FALSE.

        LOGICAL :: lspinorb = .FALSE.

        REAL(DP) :: lambda = 1.0_DP

        REAL(DP) :: fixed_magnetization(3) = 0.0_DP

        REAL(DP) :: angle1(nsx) = 0.0_DP

        REAL(DP) :: angle2(nsx) = 0.0_DP

        INTEGER :: report = 1

        CHARACTER(LEN=80) :: constrained_magnetization = 'none'
!
! Used to perform constrained calculations in magnetic systems
! Currently available choices:
! `none` : no constraint
! `total`: total magmetization is constrained
!          If nspin=4 (noncolin=.True.) constraint is imposed by
!          adding a penalty functional to the total energy:
!  - LAMBDA * SUM_{i} ( magnetization(i) - fixed_magnetization(i) )**2
!          where the sum over i runs over the three components of
!          the magnetization. Lambda is a real number (see below).
!          If nspin=2 constraint is imposed by defining two Fermi
!          energies for spin up and down.
!          Only fixed_magnetization(3) can be defined in this case.
! `atomic`: atomic magmetization are constrained to the defined
!          starting magnetization adding a penalty
!  - LAMBDA * SUM_{i,itype} ( magnetic_moment(i,itype) - mcons(i,itype) )**2
!          where i runs over the components (1-3) and itype over
!          the types (1-ntype).
!          mcons(:,:) array is defined from starting_magnetization,
!          angle1 and angle2 variables. lambda is a real number
! `atomic direction`: not all the components of the atomic
!          magnetic moment are constrained but only the cosinus
!          of angle1, and the penalty functional is:
!  - LAMBDA * SUM_{itype} ( mag_mom(3,itype)/mag_mom_tot - cos(angle1(ityp) )**2
!
        REAL(DP) :: B_field(3) = 0.0_DP
!
! A fixed magnetic field defined by the vector B_field is added
! to the exchange and correlation magnetic field.
! The three components of the magnetic field are given in Ry.
! Only B_field(3) can be used if nspin=2.
!
        CHARACTER(LEN=80) :: sic = 'none'
          ! sic = 'none' | 'sic_mac'
          ! 'none'          no SIC
          ! 'sic_mac'       SIC correction: Exc[rhoup,rhodown]
          !                         - sic_alpha * ( Exc[rhoup, rhodwn] - Exc[ rhodwn, rhodwn ] )
          !                         + Uhartree[rhoup+rhodown] - sic_epsilon * Uhartree[rhoup-rhodown]

        REAL(DP) :: sic_epsilon = 0.0_DP
        REAL(DP) :: sic_alpha   = 0.0_DP

        LOGICAL   :: force_pairing = .FALSE.
          !  FORCEPAIRING
          !      paires electrons forcedly in two spin channels calculation
          !      works for stuctures with at most one unpaired eletron
          ! just a comment: nelu is the number of el with spin up
          ! while neld== number of el with spin down
          ! define the unpaired el with spin up

        LOGICAL :: assume_isolated = .FALSE.

        ! use spline interpolation for pseudopotential
        LOGICAL :: spline_ps = .false.

! DCC
        ! add electrostatic embedding part (details in the EE namelist)
        LOGICAL :: do_ee  = .false.
!
        NAMELIST / system / ibrav, celldm, a, b, c, cosab, cosac, cosbc, nat, &
             ntyp, nbnd, nelec, ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s,  &
             nr3s, nr1b, nr2b, nr3b, nosym, noinv, starting_magnetization,    &
             occupations, degauss, nelup, neldw, nspin, ecfixed,              &
             qcutz, q2sigma, xc_type, lda_plus_U, Hubbard_U, Hubbard_alpha,   &
             edir, emaxpos, eopreg, eamp, smearing, starting_ns_eigenvalue,   &
             U_projection_type, input_dft, la2F, assume_isolated,             &
#if defined (EXX)
             x_gamma_extrapolation, nqx1, nqx2, nqx3,                         &
#endif
             noncolin, lspinorb, lambda, angle1, angle2, report,              &
             constrained_magnetization, B_field, fixed_magnetization,         &
             sic, sic_epsilon, force_pairing, sic_alpha,                      &
             tot_charge, multiplicity, tot_magnetization,                     &
             spline_ps,                                                       &
! DCC
             do_ee
!
!=----------------------------------------------------------------------------=!
!  EE Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
! type of electrostatic embedding used
        CHARACTER( LEN = 256 ) :: which_compensation = 'none'
! kinetic energy cutoff for the coarse (MultiGrid) grid
        REAL(DP) :: ecutcoarse = 100.0d0
! amount of "new" correction introduced when mixing 
        REAL(DP) :: mixing_charge_compensation = 1.0
! error tolerance for the multigrid solver
        REAL(DP) :: errtol = 1.d-22
! how early in scf itarations should the corrective pot start being calculated
        REAL(DP) :: comp_thr = 1.d-2
! nlev number of grid levels in the multigrid solver
        INTEGER :: nlev = 2
! itmax maximum number of iterations in the multigrid solver
        INTEGER :: itmax = 1000
! whichbc 0 if aperiodic
        INTEGER :: whichbc(3) = 0
! sets after how many scf cycles the corrective potential should be calculated
        INTEGER :: n_charge_compensation = 5
!
        INTEGER :: ncompx = 1
        INTEGER :: ncompy = 1
        INTEGER :: ncompz = 1
          ! ONLY PWSCF
! 
        INTEGER :: mr1 = 0
        INTEGER :: mr2 = 0
        INTEGER :: mr3 = 0

        REAL(DP) :: cellmin( 3 ) = 0.D0
          ! ONLY PWSCF

        REAL(DP) :: cellmax( 3 ) = 1.D0

        NAMELIST / ee / which_compensation,comp_thr,    &
             ncompx,n_charge_compensation,              &
             ncompy, ncompz,mixing_charge_compensation, &
             mr1, mr2, mr3, ecutcoarse,                 &
             errtol, nlev, itmax, whichbc,              &
             cellmin, cellmax                                                

!=----------------------------------------------------------------------------=!
!  ELECTRONS Namelist Input Parameters
!=----------------------------------------------------------------------------=!

        REAL(DP) :: emass = 0.0_DP
          ! effective electron mass in the CP Lagrangian,
          ! in atomic units ( 1 a.u. of mass = 1/1822.9 a.m.u. = 9.10939 * 10^-31 kg )
          ! Typical values in CP simulation are between 100. and 1000.

        REAL(DP) :: emass_cutoff = 0.0_DP
          ! mass cut-off (in Rydbergs) for the Fourier acceleration
          ! effective mass is rescaled for "G" vector components with kinetic
          ! energy above "emass_cutoff"
          ! Use a value grether than "ecutwfc" to disable Fourier acceleration.

        CHARACTER(LEN=80) :: orthogonalization = 'ortho'
          ! orthogonalization = 'Gram-Schmidt' | 'ortho'*
          ! selects the orthonormalization method for electronic wave functions
          !  'Gram-Schmidt'  use Gram-Schmidt algorithm
          !  'ortho'         use iterative algorithm

        REAL(DP) :: ortho_eps = 1.E-8_DP
          ! meaningful only if orthogonalization = 'ortho'
          ! tolerance for iterative orthonormalization,
          ! a value of 1.d-8 is usually sufficent

        INTEGER   :: ortho_max = 20
          ! meaningful only if orthogonalization = 'ortho'
          ! maximum number of iterations for orthonormalization
          ! usually between 15 and 30.

        INTEGER   :: ortho_para = 0
          ! meaningful only if orthogonalization = 'ortho' and parallel build
          ! Suggested number of processors to be used for distributing
          ! lambda matrixes and for parallel diagonalization

        INTEGER :: electron_maxstep = 1000
          ! maximum number of steps in electronic minimization
          ! This parameter apply only when using 'cg' electronic or
          ! ionic dynamics

        CHARACTER(LEN=80) :: electron_dynamics = 'none'
          ! electron_dynamics = 'default' | 'sd' | 'cg' | 'damp' | 'none' | 'verlet'
          ! set how electrons shold be moved
          ! 'none'    electronic degrees of fredom (d.o.f.) are kept fixed
          ! 'sd'      steepest descent algorithm is used to minimize electronic d.o.f.
          ! 'cg'      conjugate gradient algorithm is used to minimize electronic d.o.f.
          ! 'damp'    damped dynamics is used to propagate electronic d.o.f.
          ! 'verlet'  standard Verlet algorithm is used to propagate electronic d.o.f.
          ! 'default' the value depends on variable "calculation"

        CHARACTER(LEN=80) :: electron_dynamics_allowed(6)
        DATA electron_dynamics_allowed &
          / 'default', 'sd', 'cg', 'damp', 'verlet', 'none' /

        REAL(DP) :: electron_damping = 0.0_DP
          ! meaningful only if " electron_dynamics = 'damp' "
          ! damping frequency times delta t, optimal values could be
          ! calculated with the formula
          !        sqrt(0.5*log((E1-E2)/(E2-E3)))
          ! where E1 E2 E3 are successive values of the DFT total energy
          ! in a steepest descent simulations

        CHARACTER(LEN=80) :: electron_velocities = 'default'
          ! electron_velocities = 'zero' | 'default'*
          ! 'zero'    restart setting electronic velocities to zero
          ! 'default' restart using electronic velocities of the previous run

        CHARACTER(LEN=80) :: electron_temperature = 'not_controlled'
          ! electron_temperature = 'nose' | 'not_controlled'* | 'rescaling'
          ! 'nose'           control electronic temperature using Nose thermostat
          !                  see parameter "fnosee" and "ekincw"
          ! 'rescaling'      control electronic temperature via velocities rescaling
          ! 'not_controlled' electronic temperature is not controlled

        REAL(DP) :: ekincw = 0.0_DP
          ! meaningful only with "electron_temperature /= 'not_controlled' "
          ! value of the average kinetic energy (in atomic units) forced
          ! by the temperature control

        REAL(DP) :: fnosee = 0.0_DP
          ! meaningful only with "electron_temperature = 'nose' "
          ! oscillation frequency of the nose thermostat (in terahertz)

        CHARACTER(LEN=80) :: startingwfc = 'random'
          ! startingwfc = 'atomic' | 'atomic+random' | 'random' | 'file'
          ! define how the code should initialize the wave function
          ! 'atomic'   start from superposition of atomic wave functions
          ! 'atomic+random' as above, plus randomization
          ! 'random'   start from random wave functions
          ! 'file'     read wavefunctions from file

        REAL(DP) :: ampre = 0.0_DP
          ! meaningful only if "startingwfc = 'random'", amplitude of the
          ! randomization ( allowed values: 0.0 - 1.0 )

        REAL(DP) :: grease = 0.0_DP
          ! a number <= 1, very close to 1: the damping in electronic
          ! damped dynamics is multiplied at each time step by "grease"
          ! (avoids overdamping close to convergence: Obsolete ?)
          ! grease = 1 : normal damped dynamics
          ! NOT used in FPMD

        INTEGER :: empty_states_nbnd = 0
          ! number of empty states to be calculated every iprint steps
          ! default value is zero

        INTEGER :: empty_states_maxstep = 100
          ! meaningful only with "empty_states_nbnd > 0 "
          ! maximum number of iteration in the empty states calculation
          ! default is 100

        REAL(DP) :: empty_states_ethr = 1.E-4_DP
          ! meaningful only with "empty_states_nbnd > 0 "
          ! wave function gradient threshold, for convergence of empty states
          ! default value is ekin_conv_thr

        INTEGER :: diis_size = 0
          ! meaningful only with " electron_dynamics = 'diis' "
          ! size of the matrix used for the inversion in the iterative subspace
          ! default is 4, allowed value 1-5

        INTEGER :: diis_nreset = 0
          ! meaningful only with " electron_dynamics = 'diis' "
          ! number of steepest descendent step after a reset of the diis
          ! iteration, default value is 3

        REAL(DP) :: diis_hcut = 0.0_DP
          ! meaningful only with " electron_dynamics = 'diis' "
          ! energy cutoff (a.u.), above which an approximate diagonal
          ! hamiltonian is used in finding the direction to the minimum
          ! default is "1.0"

        REAL(DP) :: diis_wthr = 1.E-4_DP
          ! meaningful only with " electron_dynamics = 'diis' "
          ! convergence threshold for wave function
          ! this criterion is satisfied when the maximum change
          ! in the wave functions component between two diis steps
          ! is less than this threshold
          ! default value is ekin_conv_thr

        REAL(DP) :: diis_delt = 1.0_DP
          ! meaningful only with " electron_dynamics = 'diis' "
          ! electronic time step used in the steepest descendent step
          ! default is "dt"

        INTEGER :: diis_maxstep = 100
          ! meaningful only with " electron_dynamics = 'diis' "
          ! maximum number of iteration in the diis minimization
          ! default is electron_maxstep

        LOGICAL :: diis_rot = .FALSE.
          ! meaningful only with " electron_dynamics = 'diis' "
          ! if "diis_rot = .TRUE." enable diis with charge mixing and rotations
          ! default is "diis_rot = .FALSE."

        REAL(DP) :: diis_fthr = 1.E-3_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! convergence threshold for ionic force
          ! this criterion is satisfied when the maximum change
          ! in the atomic force between two diis steps
          ! is less than this threshold
          ! default value is "0.0"

        REAL(DP) :: diis_temp = 0.0_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! electronic temperature, significant only if ???

        REAL(DP) :: diis_achmix  = 0.0_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! "A" parameter in the charge mixing formula
          ! chmix = A * G^2 / (G^2 + G0^2) , G represents reciprocal lattice vectors

        REAL(DP) :: diis_g0chmix  = 0.0_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! "G0^2" parameter in the charge mixing formula

        INTEGER :: diis_nchmix = 0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! dimension of the charge mixing

        REAL(DP) :: diis_g1chmix = 0.0_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! "G1^2" parameter in the charge mixing formula
          ! metric = (G^2 + G1^2) / G^2 , G represents reciprocal lattice vectors

        INTEGER :: diis_nrot(3) = 0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! start upgrading the charge density every "diis_nrot(1)" steps,
          ! then every "diis_nrot(2)", and at the end every "diis_nrot(3)",
          ! depending on "diis_rothr"

        REAL(DP) :: diis_rothr(3) = 1.E-4_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! threshold on the charge difference between two diis step
          ! when max charge difference is less than "diis_rothr(1)", switch
          ! between the "diis_nrot(1)" upgrade frequency to "diis_nrot(2)",
          ! then when the max charge difference is less than "diis_rothr(2)",
          ! switch between "diis_nrot(2)" and "diis_nrot(3)", upgrade frequency,
          ! finally when the max charge difference is less than "diis_nrot(3)"
          ! convergence is achieved

        REAL(DP) :: diis_ethr = 1.E-4_DP
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! convergence threshold for energy
          ! this criterion is satisfied when the change
          ! in the energy between two diis steps
          ! is less than this threshold
          ! default value is etot_conv_thr

        LOGICAL :: diis_chguess = .FALSE.
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! if "diis_chguess = .TRUE." enable charge density guess
          ! between two diis step, defaut value is "diis_chguess = .FALSE."

        CHARACTER(LEN=80) :: mixing_mode = 'default'
          ! mixing = ????
          ! define how to mix wave functions
          ! NOT used in FPMD

        REAL(DP) :: mixing_beta = 0.0_DP
          ! parameter for wave function mixing
          ! NOT used in FPMD

        INTEGER :: mixing_ndim = 0
          ! dimension of wave function mixing
          ! NOT used in FPMD

        CHARACTER(LEN=80) :: diagonalization = 'cg'
          ! diagonalization = 'cg' | 'david' | 'david-serial'
          ! algorithm used by PWscf for iterative diagonalization

        REAL(DP) :: diago_thr_init = 0.0_DP
          ! convergence threshold for the firts iterative diagonalization.
          ! NOT used in FPMD

        INTEGER :: diago_cg_maxiter = 100
          ! NOT used in FPMD

        INTEGER :: diago_david_ndim = 10
          ! NOT used in FPMD

        INTEGER :: diago_diis_ndim = 10
          ! NOT used in FPMD

        LOGICAL :: diago_full_acc = .FALSE.

        REAL(DP) :: conv_thr = 1.E-6_DP
          ! convergence threshold in electronic ONLY minimizations
          ! NOT used in FPMD

        INTEGER :: mixing_fixed_ns  = 0
          ! PWSCF only
          ! NOT used in FPMD

        CHARACTER(LEN=80) :: startingpot = 'potfile'
          ! specify the file containing the DFT potential of the system
          ! NOT used in FPMD

        INTEGER :: n_inner = 2
          ! number of inner loop per CG iteration.

        INTEGER :: niter_cold_restart = 1
          !frequency of full cold smearing inner cycle (in iterations)

        REAL(DP) :: lambda_cold
         !step for not complete cold smearing inner cycle


        LOGICAL :: tgrand = .FALSE.
          ! whether to do grand-canonical calculations.

        REAL(DP) :: fermi_energy = 0.0_DP
          ! chemical potential of the grand-canonical ensemble.

        CHARACTER(LEN=80) :: rotation_dynamics = "line-minimization"
          ! evolution the rotational degrees of freedom.

        CHARACTER(LEN=80) :: occupation_dynamics = "line-minimization"
          ! evolution of the occupational degrees of freedom.

        REAL(DP) :: rotmass = 0
          ! mass for the rotational degrees of freedom.

        REAL(DP) :: occmass = 0
          ! mass for the occupational degrees of freedom.

        REAL(DP) :: occupation_damping = 0
          ! damping for the rotational degrees of freedom.

        REAL(DP) :: rotation_damping = 0
          ! damping for the occupational degrees of freedom.

        LOGICAL :: tcg = .true.
          ! if true perform in cpv conjugate gradient minimization of electron energy

        INTEGER :: maxiter = 40
          ! man number of conjugate gradient iterations

        REAL(DP)  :: etresh =1.0E-7_DP
          ! treshhold on energy

        REAL(DP) :: passop =0.3_DP
          ! small step for parabolic interpolation

        INTEGER :: niter_cg_restart
          !frequency of restart for the conjugate gradient algorithm in iterations

        INTEGER  :: epol = 3
          ! electric field direction

        REAL(DP) :: efield =0.0_DP
          ! electric field intensity in atomic units

        REAL(DP) :: efield_cart(3)!electric field vector in carthesian system of reference

       INTEGER  :: epol2 = 3
          ! electric field direction

        REAL(DP) :: efield2 =0.0_DP
          ! electric field intensity in atomic units

        LOGICAL :: tqr = .FALSE.
          ! US contributions are added in real space

        LOGICAL :: occupation_constraints = .FALSE.
          ! If true perform CP dynamics with constrained occupations
          ! to be used together with penalty functional ...

        NAMELIST / electrons / emass, emass_cutoff, orthogonalization, &
          electron_maxstep, ortho_eps, ortho_max, electron_dynamics,   &
          electron_damping, electron_velocities, electron_temperature, &
          ekincw, fnosee, ampre, grease, empty_states_nbnd,            &
          empty_states_maxstep,                                        &
          empty_states_ethr, diis_size, diis_nreset, diis_hcut,        &
          diis_wthr, diis_delt, diis_maxstep, diis_rot, diis_fthr,     &
          diis_temp, diis_achmix, diis_g0chmix, diis_g1chmix,          &
          diis_nchmix, diis_nrot, diis_rothr, diis_ethr, diis_chguess, &
          mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns,      &
          tqr, diago_cg_maxiter, diago_david_ndim, diagonalization,    &
          startingpot, startingwfc , conv_thr, diago_diis_ndim,        &
          diago_thr_init, n_inner, fermi_energy, rotmass, occmass,     &
          rotation_damping, occupation_damping, rotation_dynamics,     &
          occupation_dynamics, tcg, maxiter, etresh, passop, epol,     &
          efield, epol2, efield2, diago_full_acc,                      &
          occupation_constraints, ortho_para, niter_cg_restart,        &
          niter_cold_restart, lambda_cold, efield_cart

!
!=----------------------------------------------------------------------------=!
!  IONS Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
        CHARACTER(LEN=80) :: phase_space = 'full'
          ! phase_space = 'full' | 'coarse-grained'
          ! 'full'             the full phase-space is used for the ionic
          !                    dynamics
          ! 'coarse-grained'   a coarse-grained phase-space, defined by a set
          !                    of constraints, is used for the ionic dynamics

        CHARACTER(LEN=80) :: phase_space_allowed(2)
        DATA phase_space_allowed / 'full', 'coarse-grained' /

        CHARACTER(LEN=80) :: ion_dynamics = 'none'
          ! ion_dynamics = 'sd' | 'cg' | 'damp' | 'verlet' | 'bfgs' | 'none'*
          ! set how ions shold be moved
          ! 'none'     ions are kept fixed
          ! 'bfgs'     a BFGS algorithm is used to minimize ionic configuration
          ! 'sd'       steepest descent algorithm is used to minimize ionic configuration
          ! 'cg'       conjugate gradient algorithm is used to minimize ionic configuration
          ! 'damp'     damped dynamics is used to propagate ions
          ! 'verlet'   standard Verlet algorithm is used to propagate ions

        CHARACTER(LEN=80) :: ion_dynamics_allowed(8)
        DATA ion_dynamics_allowed / 'none', 'sd', 'cg', 'langevin', &
                                    'damp', 'verlet', 'bfgs', 'beeman' /

        REAL(DP) :: ion_radius(nsx) = 0.5_DP
          ! pseudo-atomic radius of the i-th atomic species
          ! (for Ewald summation), values between 0.5 and 2.0 are usually used.

        REAL(DP) :: ion_damping = 0.2_DP
          ! meaningful only if " ion_dynamics = 'damp' "
          ! damping frequency times delta t, optimal values could be
          ! calculated with the formula
          !        sqrt(0.5*log((E1-E2)/(E2-E3)))
          ! where E1 E2 E3 are successive values of the DFT total energy
          ! in a ionic steepest descent simulation

        CHARACTER(LEN=80) :: ion_positions = 'default'
          ! ion_positions = 'default'* | 'from_input'
          ! 'default'    restart the simulation with atomic positions read
          !              from the restart file
          ! 'from_input' restart the simulation with atomic positions read
          !              from standard input ( see the card 'ATOMIC_POSITIONS' )

        CHARACTER(LEN=80) :: ion_velocities = 'default'
          ! ion_velocities = 'zero' | 'default'* | 'random' | 'from_input'
          ! 'default'    restart the simulation with atomic velocities read
          !              from the restart file
          ! 'random'     start the simulation with random atomic velocities
          ! 'from_input' restart the simulation with atomic velocities read
          !              from standard input (see the card 'ATOMIC_VELOCITIES' )
          ! 'zero'       restart the simulation with atomic velocities set to zero

        CHARACTER(LEN=80) :: ion_temperature = 'not_controlled'
          ! ion_temperature = 'nose' | 'not_controlled'* | 'rescaling' |
          !    'berendsen' | 'andersen' | 'rescale-v' | 'rescale-T' | 'reduce-T'
          !
          ! 'nose'           control ionic temperature using Nose thermostat
          !                  see parameters "fnosep" and "tempw"
          ! 'rescaling'      control ionic temperature via velocity rescaling
          !                  see parameters "tempw" and "tolp"
          ! 'rescale-v'      control ionic temperature via velocity rescaling
          !                  see parameters "tempw" and "nraise"
          ! 'rescale-T'      control ionic temperature via velocity rescaling
          !                  see parameter "delta_t"
          ! 'reduce-T'       reduce ionic temperature
          !                  see parameters "nraise", delta_t"
          ! 'berendsen'      control ionic temperature using "soft" velocity
          !                  rescaling - see parameters "tempw" and "nraise"
          ! 'andersen'       control ionic temperature using Andersen thermostat
          !                  see parameters "tempw" and "nraise"
          ! 'not_controlled' ionic temperature is not controlled

        REAL(DP) :: tempw = 300.0_DP
          ! meaningful only with "ion_temperature /= 'not_controlled' "
          ! value of the ionic temperature (in Kelvin) forced
          ! by the temperature control

        REAL(DP) :: fnosep( nhclm )  = 50.0_DP
          ! meaningful only with "ion_temperature = 'nose' "
          ! oscillation frequency of the nose thermostat (in terahertz)
          ! nhclm is a max length for the chain

        INTEGER   ::  nhpcl = 0
          ! non-zero only with "ion_temperature = 'nose' "
          ! this defines the length of the Nose-Hoover chain

        INTEGER   :: nhptyp = 0
        ! this parameter set the nose hoover thermostat to more than one

        INTEGER   ::  nhgrp(nsx)=0
          ! this is the array to assign thermostats to atomic types
          ! allows to use various thermostat setups

        INTEGER   ::  ndega = 0
          ! this is the parameter to control active degrees of freedom
          ! used for temperature control and the Nose-Hoover chains

        REAL(DP) :: tolp = 50.0_DP
          ! meaningful only with "ion_temperature = 'rescaling' "
          ! tolerance (in Kelvin) of the rescaling. When ionic temperature
          ! differs from "tempw" more than "tolp" apply rescaling.

        REAL(DP)  ::  fnhscl(nsx)=-1.0_DP
        ! this is to scale the target energy, in case there are constraints
        ! the dimension is the same as nhgrp, meaning that atomic type
        ! i with a group nhgrp(i) is scaled by fnhscl(i)
        
        LOGICAL   :: tranp(nsx) = .FALSE.
          ! tranp(i) control the randomization of the i-th atomic specie
          ! .TRUE.   randomize ionic positions ( see "amprp" )
          ! .FALSE.  do nothing

        REAL(DP) :: amprp(nsx) = 0.0_DP
          ! amprp(i) meaningful only if "tranp(i) = .TRUE.", amplitude of the
          ! randomization ( allowed values: 0.0 - 1.0 ) for the i-th atomic specie.
          ! Add to the positions a random displacements vector ( in bohr radius )
          ! defined as:  amprp( i ) * ( X, Y, Z )
          ! where X, Y, Z are pseudo random number in the interval [ -0.5 , 0.5 ]

        REAL(DP) :: greasp = 0.0_DP
          ! same as "grease", for ionic damped dynamics
          ! NOT used in FPMD

        INTEGER   :: ion_nstepe = 1
          ! number of electronic steps for each ionic step

        INTEGER   :: ion_maxstep = 1000
          ! maximum number of step in ionic minimization

        REAL(DP) :: upscale = 1.0_DP
          ! This variable is NOT used in FPMD

        CHARACTER(LEN=80) :: pot_extrapolation = 'default', &
                             wfc_extrapolation = 'default'
          !  These variables are used only by PWSCF

        LOGICAL :: refold_pos
        LOGICAL :: remove_rigid_rot = .FALSE.

        !
        ! ... delta_T, nraise, tolp are used to change temperature in PWscf
        !

        REAL(DP) :: delta_t = 1.0_DP

        INTEGER :: nraise = 1

        !
        ! ... variables added for "path" calculations
        !

        !
        ! ... these are two auxiliary variables used in read_cards to
        ! ... distinguish among neb and smd done in the full phase-space
        ! ... or in the coarse-grained phase-space
        !
        LOGICAL :: full_phs_path_flag = .FALSE.
        LOGICAL :: cg_phs_path_flag   = .FALSE.
        !
        INTEGER :: input_images = 0
        !
        INTEGER :: num_of_images = 0
        !
        CHARACTER(LEN=80) :: CI_scheme = 'no-CI'
          ! CI_scheme = 'no-CI' | 'auto' | 'manual'
          ! set the Climbing Image scheme
          ! 'no-CI'       Climbing Image is not used
          ! 'auto'        Standard Climbing Image
          ! 'manual'      the image is selected by hand
        !
        CHARACTER(LEN=80) :: CI_scheme_allowed(3)
        DATA CI_scheme_allowed / 'no-CI', 'auto', 'manual' /
        !
        LOGICAL :: first_last_opt = .FALSE.
        LOGICAL :: use_masses     = .FALSE.
        LOGICAL :: use_freezing   = .FALSE.
        LOGICAL :: fixed_tan      = .FALSE.
        !
        CHARACTER(LEN=80) :: opt_scheme = 'quick-min'
          ! minimization_scheme = 'quick-min' | 'damped-dyn' |
          !                       'mol-dyn'   | 'sd'
          ! set the minimization algorithm
          ! 'quick-min'   projected molecular dynamics
          ! 'sd'          steepest descent
          ! 'broyden'     broyden acceleration
          ! 'langevin'    langevin dynamics
        !
        CHARACTER(LEN=80) :: opt_scheme_allowed(4)
        DATA opt_scheme_allowed / 'quick-min', 'broyden', 'sd', 'langevin' /
        !
        REAL (DP)  :: temp_req = 0.0_DP
          ! meaningful only when minimization_scheme = 'sim-annealing'
        REAL (DP)  :: ds = 1.0_DP
        !
        REAL (DP)  :: k_max = 0.1_DP, k_min = 0.1_DP
        !
        REAL (DP)  :: path_thr = 0.05_DP

        !
        ! ... variables added for new BFGS algorithm
        !

        INTEGER ::  bfgs_ndim = 1

        REAL(DP)  :: trust_radius_max = 0.8_DP
        REAL(DP)  :: trust_radius_min = 1.E-3_DP
        REAL(DP)  :: trust_radius_ini = 0.5_DP

        REAL(DP)  :: w_1 = 0.5E-1_DP
        REAL(DP)  :: w_2 = 0.5_DP

        REAL(DP)  :: sic_rloc = 0.0_DP

        !
        ! ... variable for meta-dynamics
        !

        INTEGER  :: fe_nstep = 100
        INTEGER  :: sw_nstep = 10
        INTEGER  :: eq_nstep = 0
        REAL(DP) :: g_amplitude = 0.005_DP
        !
        REAL(DP) :: fe_step( max_nconstr ) = 0.4_DP
        !
        NAMELIST / ions / phase_space, ion_dynamics, ion_radius, ion_damping,  &
                          ion_positions, ion_velocities, ion_temperature,      &
                          tempw, fnosep, nhgrp, fnhscl, nhpcl, nhptyp, ndega, tranp,   &
                          amprp, greasp, tolp, ion_nstepe, ion_maxstep,        &
                          refold_pos, upscale, delta_t, pot_extrapolation,     &
                          wfc_extrapolation, nraise, remove_rigid_rot,         &
                          num_of_images, CI_scheme, opt_scheme, use_masses,    &
                          first_last_opt, ds, k_max, k_min, temp_req,          &
                          path_thr, fixed_tan, use_freezing,                   &
                          trust_radius_max, trust_radius_min,                  &
                          trust_radius_ini, w_1, w_2, bfgs_ndim, sic_rloc,     &
                          fe_step, fe_nstep, sw_nstep, eq_nstep, g_amplitude

!=----------------------------------------------------------------------------=!
!  CELL Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!

        CHARACTER(LEN=80) :: cell_parameters = 'default'
          ! cell_parameters = 'default'* | 'from_input'
          ! 'default'    restart the simulation with cell parameters read
          !              from the restart file or "celldm" if
          !              "restart = 'from_scratch'"
          ! 'from_input' restart the simulation with cell parameters
          !              from standard input ( see the card 'CELL_PARAMETERS' )

        CHARACTER(LEN=80) :: cell_dynamics  = 'none'
          ! cell_dynamics = 'sd' | 'pr' | 'none'*
          ! set how cell shold be moved
          ! 'none'   cell is kept fixed
          ! 'sd'     steepest descent algorithm is used to minimize the cell
          ! 'pr'     standard Verlet algorithm is used to propagate the cell

        CHARACTER(LEN=80) :: cell_dynamics_allowed(7)
        DATA cell_dynamics_allowed / 'sd', 'pr', 'none', 'w', 'damp-pr', 'damp-w', 'bfgs'  /

        CHARACTER(LEN=80) :: cell_velocities = 'default'
          ! cell_velocities = 'zero' | 'default'*
          ! 'zero'    restart setting cell velocitiy to zero
          ! 'default' restart using cell velocity of the previous run

        REAL(DP) :: press = 0.0_DP
          ! external pressure (in GPa, remember 1 kbar = 10^8 Pa)

        REAL(DP) :: wmass = 0.0_DP
          ! effective cell mass in the Parrinello-Rahman Lagrangian (in atomic units)
          ! of the order of magnitude of the total atomic mass
          ! (sum of the mass of the atoms) within the simulation cell.
          ! if you do not specify this parameters, the code will compute
          ! its value based on some physical consideration

        CHARACTER(LEN=80) :: cell_temperature  = 'not_controlled'
          ! cell_temperature = 'nose' | 'not_controlled'* | 'rescaling'
          ! 'nose'           control cell temperature using Nose thermostat
          !                  see parameters "fnoseh" and "temph"
          ! 'rescaling'      control cell temperature via velocities rescaling
          ! 'not_controlled' cell temperature is not controlled
          ! NOT used in FPMD

        REAL(DP) :: temph = 0.0_DP
          ! meaningful only with "cell_temperature /= 'not_controlled' "
          ! value of the cell temperature (in Kelvin) forced
          ! by the temperature control

        REAL(DP) :: fnoseh = 1.0_DP
          ! meaningful only with "cell_temperature = 'nose' "
          ! oscillation frequency of the nose thermostat (in terahertz)

        REAL(DP) :: greash = 0.0_DP
          ! same as "grease", for cell damped dynamics

        CHARACTER(LEN=80) :: cell_dofree = 'all'
          ! cell_dofree = 'all'* | 'volume' | 'x' | 'y' | 'z' | 'xy' | 'xz' | 'yz' | 'xyz'
          ! select which of the cell parameters should be moved
          ! 'all'    all axis and angles are propagated (default)
          ! 'volume' the cell is simply rescaled, without changing the shape
          ! 'x'      only the "x" axis is moved
          ! 'y'      only the "y" axis is moved
          ! 'z'      only the "z" axis is moved
          ! 'xy'     only the "x" and "y" axis are moved, angles are unchanged
          ! 'xz'     only the "x" and "z" axis are moved, angles are unchanged
          ! 'yz'     only the "y" and "z" axis are moved, angles are unchanged
          ! 'xyz'    "x", "y" and "z" axis are moved, angles are unchanged

        REAL(DP) :: cell_factor = 0.0_DP
          ! NOT used in FPMD

        INTEGER   :: cell_nstepe = 1
          ! number of electronic steps for each cell step

        REAL(DP) :: cell_damping = 0.0_DP
          ! meaningful only if " cell_dynamics = 'damp' "
          ! damping frequency times delta t, optimal values could be
          ! calculated with the formula
          !        sqrt(0.5*log((E1-E2)/(E2-E3)))
          ! where E1 E2 E3 are successive values of the DFT total energy
          ! in a ionic steepest descent simulation

        REAL(DP) :: press_conv_thr = 0.5_DP

        NAMELIST / cell / cell_parameters, cell_dynamics, cell_velocities, &
                          press, wmass, cell_temperature, temph, fnoseh,   &
                          cell_dofree, greash, cell_factor, cell_nstepe,   &
                          cell_damping, press_conv_thr

!
!=----------------------------------------------------------------------------=!!
! PRESS_AI Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
!
      LOGICAL  :: abivol = .FALSE.
      LOGICAL  :: abisur = .FALSE.
      LOGICAL  :: pvar   = .FALSE.
      LOGICAL  :: fill_vac=.FALSE.
      LOGICAL  :: scale_at=.FALSE.
      LOGICAL  :: t_gauss =.FALSE.
      LOGICAL  :: jellium= .FALSE.
      LOGICAL  :: cntr(nsx)=.FALSE.
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

!
!=----------------------------------------------------------------------------=!
!  PHONON Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!

        INTEGER :: modenum = 0

        REAL(DP) :: xqq(3) = 0.0_DP
          ! coordinates of q point for phonon calculation

        NAMELIST / phonon / modenum, xqq

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
          LOGICAL :: wfsd
          !
          REAL(DP) :: wfdt
          REAL(DP) :: maxwfdt
          REAL(DP) :: wf_q
          REAL(DP) :: wf_friction
          !
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
          !
          NAMELIST / wannier / wf_efield, wf_switch, sw_len, efx0, efy0, efz0, &
                               efx1, efy1, efz1, wfsd, wfdt, maxwfdt, wf_q,    &
                               wf_friction, nit, nsd, nsteps, tolw, adapt,     &
                               calwf, nwf, wffort, writev

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
        CHARACTER(LEN=3)  :: atom_label(nsx) = 'XX'   ! label of the atomic species being read
        CHARACTER(LEN=80) :: atom_pfile(nsx) = 'YY'   ! pseudopotential file name
        REAL(DP)          :: atom_mass(nsx)  = 0.0_DP ! atomic mass of the i-th atomic species
          ! in atomic mass units: 1 a.m.u. = 1822.9 a.u. = 1.6605 * 10^-27 kg
        LOGICAL   :: taspc = .FALSE.

!
!    ATOMIC_POSITIONS
!
        REAL(DP), ALLOCATABLE :: rd_pos(:,:)  ! unsorted positions from input
        INTEGER,  ALLOCATABLE :: sp_pos(:)
        INTEGER,  ALLOCATABLE :: if_pos(:,:)
        INTEGER,  ALLOCATABLE :: id_loc(:)
        INTEGER,  ALLOCATABLE :: na_inp(:)
        LOGICAL  :: tapos = .FALSE.
        CHARACTER(LEN=80) :: atomic_positions = 'crystal'
          ! atomic_positions = 'bohr' | 'angstrong' | 'crystal' | 'alat'
          ! select the units for the atomic positions being read from stdin

        !
        ! ... variable added for NEB  ( C.S. 17/10/2003 )
        !
        REAL(DP), ALLOCATABLE :: pos(:,:)
        !
        ! ... workaround for IBM xlf bug, compiler can't manage large
        !     array initialization

!
!    ION_VELOCITIES
!
        REAL(DP), ALLOCATABLE :: rd_vel(:,:)   ! unsorted velocities from input
        INTEGER,  ALLOCATABLE :: sp_vel(:)
        LOGICAL  :: tavel          = .FALSE.

!
!    KPOINTS
!
! ...   k-points inputs
        LOGICAL :: tk_inp = .FALSE.
        REAL(DP) :: xk(3,npkx) = 0.0_DP, wk(npkx) = 0.0_DP
        INTEGER :: nkstot = 0, nk1 = 0, nk2 = 0, nk3 = 0, k1 = 0, k2 = 0, k3 = 0
        CHARACTER(LEN=80) :: k_points = 'gamma'
          ! k_points = 'automatic' | 'crystal' | 'tpiba' | 'gamma'*
          ! select the k points mesh
          ! 'automatic'  k points mesh is generated automatically
          !              with Monkhorst-Pack algorithm
          ! 'crystal'    k points mesh is given in stdin in scaled units
          ! 'tpiba'      k points mesh is given in stdin in units of ( 2 PI / alat )
          ! 'gamma'      only gamma point is used ( default in CPMD simulation )
!
!    NEWNFI
!
        LOGICAL :: tnewnfi_card = .FALSE.
        INTEGER :: newnfi_card = 0

!
!    2DPROCMESH
!
        LOGICAL :: t2dpegrid_inp = .FALSE.

!
!    OCCUPATIONS
!
        REAL(DP), ALLOCATABLE :: f_inp(:,:)
        LOGICAL   :: tf_inp = .FALSE.

!
!    VHMEAN
!
! ...   card planar mean of the Hartree potential
        LOGICAL :: tvhmean_inp = .FALSE.
        INTEGER :: vhnr_inp = 0, vhiunit_inp = 0
        REAL(DP)  :: vhrmin_inp = 0.0_dP, vhrmax_inp = 0.0_DP
        CHARACTER :: vhasse_inp = 'X'

!
!    DIPOLE
!
        LOGICAL :: tdipole_card = .FALSE.

!
!    ESR
!
       INTEGER :: iesr_inp = 1

!
!    NEIGHBOURS
!
       LOGICAL :: tneighbo = .FALSE.
       REAL(DP) :: neighbo_radius = 0.0_DP

!
!    CELL_PARAMETERS
!
       REAL(DP) :: rd_ht(3,3) = 0.0_DP
       CHARACTER(len=80) :: cell_symmetry = 'none'
       CHARACTER(LEN=80) :: cell_units = 'alat'
       LOGICAL   :: trd_ht = .FALSE.

!
!    TURBO
!
      LOGICAL :: tturbo_inp = .FALSE.
      INTEGER :: nturbo_inp = 0

!
!    CONSTRAINTS
!
      INTEGER :: nc_fields = 4   ! max number of fields that is allowed to
                                 ! define a constraint

      INTEGER  :: nconstr_inp    = 0
      REAL(DP) :: constr_tol_inp = 1.E-6_DP
      !
      CHARACTER(LEN=20), ALLOCATABLE :: constr_type_inp(:)
      REAL(DP),          ALLOCATABLE :: constr_inp(:,:)
      REAL(DP),          ALLOCATABLE :: constr_target(:)
      LOGICAL,           ALLOCATABLE :: constr_target_set(:)

!
!    COLLECTIVE_VARS
!
      INTEGER  :: ncolvar_inp    = 0
      REAL(DP) :: colvar_tol_inp = 1.E-6_DP
      !
      CHARACTER(LEN=20), ALLOCATABLE :: colvar_type_inp(:)
      REAL(DP),          ALLOCATABLE :: colvar_inp(:,:)
      REAL(DP),          ALLOCATABLE :: colvar_target(:)
      LOGICAL,           ALLOCATABLE :: colvar_target_set(:)

!
!    KOHN_SHAM
!
      INTEGER, ALLOCATABLE :: iprnks( :, : )
      INTEGER :: nprnks( nspinx ) = 0
        ! logical mask used to specify which kohn sham orbital should be
        ! written to files 'KS.'
      INTEGER, ALLOCATABLE :: iprnks_empty( :, : )
      INTEGER :: nprnks_empty( nspinx ) = 0
        ! logical mask used to specify which empty kohn sham orbital should be
        ! written to files 'KS_EMP.'

!
!    CHI2
!
      LOGICAL :: tchi2_inp = .FALSE.

!
!   CLIMBING_IMAGES
!
      !
      ! ... variable added for NEB  ( C.S. 20/11/2003 )
      !
      LOGICAL, ALLOCATABLE :: climbing( : )

!
!   PLOT_WANNIER
!

      INTEGER, PARAMETER :: nwf_max = 1000
      !
      INTEGER :: wannier_index( nwf_max )

!  END manual
! ----------------------------------------------------------------------

CONTAINS

  SUBROUTINE allocate_input_ions( ntyp, nat )
    !
    INTEGER, INTENT(IN) :: ntyp, nat
    !
    IF ( ALLOCATED( rd_pos ) ) DEALLOCATE( rd_pos )
    IF ( ALLOCATED( sp_pos ) ) DEALLOCATE( sp_pos )
    IF ( ALLOCATED( if_pos ) ) DEALLOCATE( if_pos )
    IF ( ALLOCATED( id_loc ) ) DEALLOCATE( id_loc )
    IF ( ALLOCATED( na_inp ) ) DEALLOCATE( na_inp )
    IF ( ALLOCATED( rd_vel ) ) DEALLOCATE( rd_vel )
    IF ( ALLOCATED( sp_vel ) ) DEALLOCATE( sp_vel )
    !
    ALLOCATE( rd_pos( 3, nat ) )
    ALLOCATE( sp_pos( nat)   )
    ALLOCATE( if_pos( 3, nat ) )
    ALLOCATE( id_loc( nat)   )
    ALLOCATE( na_inp( ntyp)  )
    ALLOCATE( rd_vel( 3, nat ) )
    ALLOCATE( sp_vel( nat)   )
    !
    rd_pos = 0.0_DP
    sp_pos = 0
    if_pos = 1
    id_loc = 0
    na_inp = 0
    rd_vel = 0.0_DP
    sp_vel = 0
    !
    RETURN
    !
  END SUBROUTINE allocate_input_ions

  SUBROUTINE allocate_input_constr()
    !
    IF ( ALLOCATED( constr_type_inp ) )   DEALLOCATE( constr_type_inp )
    IF ( ALLOCATED( constr_inp ) )        DEALLOCATE( constr_inp )
    IF ( ALLOCATED( constr_target ) )     DEALLOCATE( constr_target )
    IF ( ALLOCATED( constr_target_set ) ) DEALLOCATE( constr_target_set )
    !
    ALLOCATE( constr_type_inp(   nconstr_inp ) )
    ALLOCATE( constr_target(     nconstr_inp ) )
    ALLOCATE( constr_target_set( nconstr_inp ) )
    !
    ALLOCATE( constr_inp( nc_fields, nconstr_inp ) )
    !
    constr_type_inp   = ' '
    constr_inp        = 0.0_DP
    constr_target     = 0.0_DP
    constr_target_set = .FALSE.
    !
    RETURN
    !
  END SUBROUTINE allocate_input_constr

  SUBROUTINE allocate_input_colvar()
    !
    IF ( ALLOCATED( colvar_type_inp ) )   DEALLOCATE( colvar_type_inp )
    IF ( ALLOCATED( colvar_inp ) )        DEALLOCATE( colvar_inp )
    IF ( ALLOCATED( colvar_target ) )     DEALLOCATE( colvar_target )
    IF ( ALLOCATED( colvar_target_set ) ) DEALLOCATE( colvar_target_set )
    !
    ALLOCATE( colvar_type_inp(   ncolvar_inp ) )
    ALLOCATE( colvar_target(     ncolvar_inp ) )
    ALLOCATE( colvar_target_set( ncolvar_inp ) )
    !
    ALLOCATE( colvar_inp( nc_fields, ncolvar_inp ) )
    !
    colvar_type_inp   = ' '
    colvar_inp        = 0.0_DP
    colvar_target     = 0.0_DP
    colvar_target_set = .FALSE.
    !
    RETURN
    !
  END SUBROUTINE allocate_input_colvar
  !
  SUBROUTINE allocate_input_iprnks( nksx, nspin )
    !
    INTEGER, INTENT(IN) :: nksx, nspin
    !
    IF( ALLOCATED( iprnks ) ) DEALLOCATE( iprnks )
    !
    ALLOCATE( iprnks( MAX( 1, nksx), nspin ) )
    !
    iprnks = 0
    !
    RETURN
    !
  END SUBROUTINE allocate_input_iprnks

  SUBROUTINE allocate_input_iprnks_empty( nksx, nspin )
    !
    INTEGER, INTENT(IN) :: nksx, nspin
    !
    IF( ALLOCATED( iprnks_empty ) ) DEALLOCATE( iprnks_empty )
    !
    ALLOCATE( iprnks_empty( MAX( 1, nksx), nspin ) )
    !
    iprnks_empty = 0
    !
    RETURN
    !
  END SUBROUTINE allocate_input_iprnks_empty

  SUBROUTINE deallocate_input_parameters()
    !
    IF ( ALLOCATED( rd_pos ) ) DEALLOCATE( rd_pos )
    IF ( ALLOCATED( sp_pos ) ) DEALLOCATE( sp_pos )
    IF ( ALLOCATED( if_pos ) ) DEALLOCATE( if_pos )
    IF ( ALLOCATED( id_loc ) ) DEALLOCATE( id_loc )
    IF ( ALLOCATED( na_inp ) ) DEALLOCATE( na_inp )
    IF ( ALLOCATED( rd_vel ) ) DEALLOCATE( rd_vel )
    IF ( ALLOCATED( sp_vel ) ) DEALLOCATE( sp_vel )
    !
    IF ( ALLOCATED( pos )    )   DEALLOCATE( pos )
    IF ( ALLOCATED( climbing ) ) DEALLOCATE( climbing )
    !
    IF ( ALLOCATED( constr_type_inp ) )   DEALLOCATE( constr_type_inp )
    IF ( ALLOCATED( constr_inp ) )        DEALLOCATE( constr_inp )
    IF ( ALLOCATED( constr_target ) )     DEALLOCATE( constr_target )
    IF ( ALLOCATED( constr_target_set ) ) DEALLOCATE( constr_target_set )
    !
    IF ( ALLOCATED( colvar_type_inp ) )   DEALLOCATE( colvar_type_inp )
    IF ( ALLOCATED( colvar_inp ) )        DEALLOCATE( colvar_inp )
    IF ( ALLOCATED( colvar_target ) )     DEALLOCATE( colvar_target )
    IF ( ALLOCATED( colvar_target_set ) ) DEALLOCATE( colvar_target_set )
    !
    IF ( ALLOCATED( iprnks ) )       DEALLOCATE( iprnks )
    IF ( ALLOCATED( iprnks_empty ) ) DEALLOCATE( iprnks_empty )
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
