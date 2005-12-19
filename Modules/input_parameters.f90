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
!  this module contains the definitions of all input parameters for CP, FPMD & 
!  PWSCF codes.
!  Originally written by Carlo Cavazzoni for FPMD
!
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nsx, natx, npkx, nbndxx, nspinx, &
                         lqmax, nhclm, max_num_of_images, max_nconstr
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
          ! calculation = 'scf', 'relax', 'md', 'cp'*, 'vc-relax', 'vc-md', 
          !               'vc-cp', 'neb', 'smd', 'cp-wf', 'fpmd', 'metadyn'
          ! Specify the type of the simulation
          ! 'scf'      = electron minimization
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
          ! 'metadyn'  = meta-dynamics (Laio-Parrinello dynamics)


        CHARACTER(LEN=80) :: calculation_allowed(15)
        DATA calculation_allowed / 'scf', 'nscf', 'relax', 'md', 'cp', &
          'vc-relax', 'vc-md', 'vc-cp', 'phonon', 'raman', 'neb', 'smd', &
          'cp-wf', 'fpmd', 'metadyn' /
          ! Allowed value for calculation parameters


        CHARACTER(LEN=80) :: verbosity = 'default'
          ! verbosity = 'high' | 'default'* | 'low' | 'minimal'
          ! define the verbosity of the code output

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

        REAL(DP) :: dt = 1.0d0
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
          ! specify the directory where the code opens output files
          ! _NOT_ for restart file

        CHARACTER(LEN=256) :: prefix = 'prefix'
          ! specify the prefix for the output file, if not specified the
          ! files are opened as standard fortran units.
          ! The prefix does _NOT_ apply to restart file

        CHARACTER(LEN=256) :: pseudo_dir = './'
          ! specify the directory containing the pseudopotentials

        REAL(DP) :: refg = 0.05
          ! Accurancy of the interpolation table, interval between
          ! table values in Rydberg

        CHARACTER(LEN=256) :: scradir = './'
          ! specify the scratch directory to store large, usually temporary, 
          ! files. When possible put this directory in the fastest available
          ! filesystem ( not NFS! )


        CHARACTER(LEN=256) :: wfcdir = 'undefined'
          ! scratch directory that is hopefully local to the node
          ! to store large, usually temporary files. Note that the 'scradir'
          ! above is quite confusing and does not do what is claimed

        REAL(DP) :: max_seconds = 1.0d+7
          ! smoothly terminate program after the specified number of seconds
          ! this parameter is typically used to prevent an hard kill from
          ! the queuing system.

        REAL(DP) :: ekin_conv_thr = 1.0d-5
          ! convergence criterion, minimizing the electrons this criterion is met 
          ! when "ekin < ekin_conv_thr"
          ! convergence is achieved when all criteria are met

        REAL(DP) :: etot_conv_thr = 1.0d-4
          ! convergence criterion, minimizing the ions this criterion is met
          ! when "etot(n+1)-etot(n) < etot_conv_thr", where "n" is the step index,
          ! and "etot" the DFT energy
          ! convergence is achieved when all criteria are met

        REAL(DP) :: forc_conv_thr = 1.0d-3
          ! convergence criterion, minimizing the ions this criterion is met
          ! when "MAXVAL(fion) < forc_conv_thr", where fion are the atomic forces
          ! convergence is achieved when all criteria are met

        CHARACTER(LEN=80) :: disk_io = 'default' 
          ! disk_io = 'high', 'default', 'low', 'minimal'
          ! Specify the amount of I/O activities ( not used in FPMD )

        LOGICAL :: tefield  = .FALSE. 
          ! if .TRUE. a finite electric field is added to the local potential
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
          ! This flag is effective only with PW code, and controls the way 
          ! in which wave functions are stored to disk, 
          !  .TRUE.  collect all wave functions and store them in restart file 
          !          ( .save )
          !  .FALSE. do not collect wave function they are left in temporary
          !          local file

        LOGICAL :: tqr = .FALSE.


        INTEGER :: printwfc=1
          ! if <0 do nothing, if==0 print rho and fort.47, if == nband print band

        NAMELIST / control / title, calculation, verbosity, restart_mode, &
          nstep, iprint, isave, tstress, tprnfor, dt, ndr, ndw, outdir, prefix, &
          wfcdir, max_seconds, ekin_conv_thr, etot_conv_thr, forc_conv_thr, &
          pseudo_dir, disk_io, tefield, dipfield, lberry, gdir, nppstr, &
          wf_collect, printwfc, scradir,lelfield, nberrycyc, refg, tqr


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

        REAL(DP) :: celldm(6) = 0.0d0
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

        REAL(DP):: nelec = 0.d0
          ! number of electrons, this parameter is MANDATORY in FPMD
          ! may be fractionary in PW, but not in CP and FPMD !

        REAL(DP):: tot_charge = 0.0d0
	  ! total system charge 

        INTEGER :: multiplicity = 0
          ! spin multiplicity (2s+1), 1 for singlet, 2 for doublet etc.
          ! when multiplicity = 0, it is unspecified

        REAL(DP) :: ecutwfc = 0.0d0
          ! energy cutoff for wave functions in k-space ( in Rydbergs )
          ! this parameter is MANDATORY in FPMD

        REAL(DP) :: ecutrho = 0.0d0
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

        REAL(DP) :: degauss = 0.0d0
          ! parameter for the smearing functions
          ! NOT used in FPMD-N

        INTEGER :: nspin = 1
          ! number of spinors
          ! "nspin = 1" for LDA simulations
          ! "nspin = 2" for LSD simulations
          ! "nspin = 4" for NON COLLINEAR simulations

        REAL(DP) :: nelup = 0.d0, neldw = 0.d0
          ! meaningful only if "nspin = 2", 
          ! "nelup" is the number of electrons with spin up
          ! "neldw" is the number of electrons with spin down
          ! Remember the following relation hold "nelec = nelup + neldw"

        LOGICAL :: nosym = .TRUE.
          ! do not use symmetry
          ! NOT used in FPMD-N

        REAL(DP) :: ecfixed = 0.0d0, qcutz = 0.0d0, q2sigma = 0.0d0
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

        REAL(DP) :: starting_magnetization( nsx ) = 0.0d0
          ! ONLY PWSCF

        LOGICAL :: lda_plus_u = .FALSE.
          ! ONLY PWSCF

        REAL(DP) :: starting_ns_eigenvalue(lqmax,nspinx,nsx)=-1.0
          ! ONLY PWSCF

        REAL(DP) :: hubbard_u(nsx) = 0.0d0
          ! ONLY PWSCF

        REAL(DP) :: hubbard_alpha(nsx) = 0.0d0
          ! ONLY PWSCF

        CHARACTER(LEN=80) :: U_projection_type = 'atomic'
          ! ONLY PWSCF
#if defined (EXX)
        INTEGER :: nqx1 = 1, nqx2 = 1, nqx3=1
          ! ONLY PWSCF
!        REAL(DP) :: yukawa = 0.d0
          ! ONLY PWSCF
#endif

        REAL(DP) :: a = 0.0d0

        REAL(DP) :: c = 0.0d0

        REAL(DP) :: b = 0.0d0

        REAL(DP) :: cosab = 0.0d0

        REAL(DP) :: cosac = 0.0d0

        REAL(DP) :: cosbc = 0.0d0

        INTEGER   :: edir = 0
 
        REAL(DP) :: emaxpos = 0.0d0

        REAL(DP) :: eopreg = 0.0d0

        REAL(DP) :: eamp = 0.0d0 

        LOGICAL :: noncolin = .FALSE.

        LOGICAL :: lspinorb = .FALSE.

        REAL(DP) :: lambda = 1.0D0

        REAL(DP) :: fixed_magnetization(3) = 0.0D0 

        REAL(DP) :: angle1(nsx) = 0.0D0

        REAL(DP) :: angle2(nsx) = 0.0D0

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
        REAL(DP) :: B_field(3) = 0.0D0 
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

        REAL(DP) :: sic_epsilon = 0.0d0
        REAL(DP) :: sic_alpha   = 0.0d0

        LOGICAL   :: force_pairing = .FALSE.
          !  FORCEPAIRING 
          !      paires electrons forcedly in two spin channels calculation
          !      works for stuctures with at most one unpaired eletron
          ! just a comment: nelu is the number of el with spin up
          ! while neld== number of el with spin down
          ! define the unpaired el with spin up



        NAMELIST / system / ibrav, celldm, a, b, c, cosab, cosac, cosbc, nat,&
             ntyp, nbnd, nelec, ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s, &
             nr3s, nr1b, nr2b, nr3b, nosym, starting_magnetization, &
             occupations, degauss, nelup, neldw, nspin, ecfixed, &
             qcutz, q2sigma, xc_type, lda_plus_U, Hubbard_U, Hubbard_alpha, &
             edir, emaxpos, eopreg, eamp, smearing, starting_ns_eigenvalue, &
             U_projection_type,  &
             input_dft, &
#if defined (EXX)
             nqx1, nqx2, nqx3, &
#endif
             noncolin, lspinorb, lambda, angle1, angle2, report, &
             constrained_magnetization, B_field, fixed_magnetization, &
             sic, sic_epsilon, force_pairing, sic_alpha, &
	     tot_charge, multiplicity


!=----------------------------------------------------------------------------=!
!  ELECTRONS Namelist Input Parameters
!=----------------------------------------------------------------------------=!

        REAL(DP) :: emass = 0.0d0
          ! effective electron mass in the CP Lagrangian, 
          ! in atomic units ( 1 a.u. of mass = 1/1822.9 a.m.u. = 9.10939 * 10^-31 kg )
          ! Typical values in CP simulation are between 100. and 1000.

        REAL(DP) :: emass_cutoff = 0.0d0
          ! mass cut-off (in Rydbergs) for the Fourier acceleration
          ! effective mass is rescaled for "G" vector components with kinetic 
          ! energy above "emass_cutoff" 
          ! Use a value grether than "ecutwfc" to disable Fourier acceleration.

        CHARACTER(LEN=80) :: orthogonalization = 'ortho'  
          ! orthogonalization = 'Gram-Schmidt' | 'ortho'*
          ! selects the orthonormalization method for electronic wave functions
          !  'Gram-Schmidt'  use Gram-Schmidt algorithm
          !  'ortho'         use iterative algorithm 

        REAL(DP) :: ortho_eps = 1.d-8
          ! meaningful only if orthogonalization = 'ortho'
          ! tolerance for iterative orthonormalization,
          ! a value of 1.d-8 is usually sufficent

        INTEGER   :: ortho_max = 20
          ! meaningful only if orthogonalization = 'ortho'
          ! maximum number of iterations for orthonormalization
          ! usually between 15 and 30.

        INTEGER :: electron_maxstep = 1000
          ! maximum number of steps in electronic minimization
          ! This parameter apply only when using 'cg' electronic or
          ! ionic dynamics

        CHARACTER(LEN=80) :: electron_dynamics = 'none'
          ! electron_dynamics = 'default' | 'sd' | 'cg' | 'damp' | 'md' | 'none' | 'diis' 
          ! set how electrons shold be moved
          ! 'none'    electronic degrees of fredom (d.o.f.) are kept fixed 
          ! 'sd'      steepest descent algorithm is used to minimize electronic d.o.f. 
          ! 'cg'      conjugate gradient algorithm is used to minimize electronic d.o.f. 
          ! 'diis'    DIIS algorithm is used to minimize electronic d.o.f. 
          ! 'damp'    damped dynamics is used to propagate electronic d.o.f. 
          ! 'verlet'  standard Verlet algorithm is used to propagate electronic d.o.f. 
          ! 'default' the value depends on variable "calculation"

        CHARACTER(LEN=80) :: electron_dynamics_allowed(7)
        DATA electron_dynamics_allowed &
          / 'default', 'sd', 'cg', 'damp', 'verlet', 'none', 'diis' /

        REAL(DP) :: electron_damping = 0.0d0
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

        REAL(DP) :: ekincw = 0.0d0
          ! meaningful only with "electron_temperature /= 'not_controlled' "
          ! value of the average kinetic energy (in atomic units) forced
          ! by the temperature control

        REAL(DP) :: fnosee = 0.0d0
          ! meaningful only with "electron_temperature = 'nose' "
          ! oscillation frequency of the nose thermostat (in terahertz) 

        CHARACTER(LEN=80) :: startingwfc = 'random'
          ! startingwfc = 'atomic' | 'random'* | 'none'
          ! define how the code should initialize the wave function
          ! 'atomic'   start from superposition of atomic wave functions
          ! 'random'   start from random wave functions

        REAL(DP) :: ampre = 0.0d0
          ! meaningful only if "startingwfc = 'random'", amplitude of the
          ! randomization ( allowed values: 0.0 - 1.0 )

        REAL(DP) :: grease = 0.0d0
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

        REAL(DP) :: empty_states_delt = 1.0d0
          ! meaningful only with "empty_states_nbnd > 0 "
          ! fictitious time step to be used in the empty states iteration
          ! default value is "dt"

        REAL(DP) :: empty_states_emass = 500.0d0
          ! meaningful only with "empty_states_nbnd > 0 "
          ! fictitious electronic mass to be used in the empty states iteration
          ! default value is "emass"

        REAL(DP) :: empty_states_ethr = 1.d-4
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

        REAL(DP) :: diis_hcut = 0.d0
          ! meaningful only with " electron_dynamics = 'diis' "
          ! energy cutoff (a.u.), above which an approximate diagonal
          ! hamiltonian is used in finding the direction to the minimum
          ! default is "1.0"

        REAL(DP) :: diis_wthr = 1.d-4
          ! meaningful only with " electron_dynamics = 'diis' "
          ! convergence threshold for wave function 
          ! this criterion is satisfied when the maximum change
          ! in the wave functions component between two diis steps
          ! is less than this threshold
          ! default value is ekin_conv_thr

        REAL(DP) :: diis_delt = 1.0d0
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

        REAL(DP) :: diis_fthr = 1.d-3
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! convergence threshold for ionic force 
          ! this criterion is satisfied when the maximum change
          ! in the atomic force between two diis steps
          ! is less than this threshold
          ! default value is "0.0"

        REAL(DP) :: diis_temp = 0.0d0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! electronic temperature, significant only if ???

        REAL(DP) :: diis_achmix  = 0.0d0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! "A" parameter in the charge mixing formula
          ! chmix = A * G^2 / (G^2 + G0^2) , G represents reciprocal lattice vectors

        REAL(DP) :: diis_g0chmix  = 0.0d0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! "G0^2" parameter in the charge mixing formula

        INTEGER :: diis_nchmix = 0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! dimension of the charge mixing

        REAL(DP) :: diis_g1chmix = 0.0d0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! "G1^2" parameter in the charge mixing formula
          ! metric = (G^2 + G1^2) / G^2 , G represents reciprocal lattice vectors

        INTEGER :: diis_nrot(3) = 0
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! start upgrading the charge density every "diis_nrot(1)" steps,
          ! then every "diis_nrot(2)", and at the end every "diis_nrot(3)",
          ! depending on "diis_rothr"

        REAL(DP) :: diis_rothr(3) = 1.d-4
          ! meaningful only with "electron_dynamics='diis' " and "diis_rot=.TRUE."
          ! threshold on the charge difference between two diis step
          ! when max charge difference is less than "diis_rothr(1)", switch
          ! between the "diis_nrot(1)" upgrade frequency to "diis_nrot(2)",
          ! then when the max charge difference is less than "diis_rothr(2)",
          ! switch between "diis_nrot(2)" and "diis_nrot(3)", upgrade frequency,
          ! finally when the max charge difference is less than "diis_nrot(3)"
          ! convergence is achieved

        REAL(DP) :: diis_ethr = 1.d-4
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

        REAL(DP) :: mixing_beta = 0.0d0
          ! parameter for wave function mixing
          ! NOT used in FPMD

        INTEGER :: mixing_ndim = 0
          ! dimension of wave function mixing
          ! NOT used in FPMD

        CHARACTER(LEN=80) :: diagonalization = 'cg'
          ! diagonalization = 'cg' | 'david' | 'david_overlap' | 'diis'
          ! NOTA: 'david' e 'david_overlap' per eliminare la variabile "loverlap"
          ! NOT used in FPMD

        REAL(DP) :: diago_thr_init = 0.D0
          ! convergence threshold for the firts iterative diagonalization.
          ! NOT used in FPMD

        INTEGER :: diago_cg_maxiter = 100
          ! NOT used in FPMD

        INTEGER :: diago_david_ndim = 10
          ! NOT used in FPMD

        INTEGER :: diago_diis_ndim = 10
          ! NOT used in FPMD

        REAL(DP) :: conv_thr = 1.d-6
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

        LOGICAL :: tgrand = .FALSE.
          ! whether to do grand-canonical calculations.

        REAL(DP) :: fermi_energy = 0.0d0
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
          ! falg per gradiente coniugato

        INTEGER :: maxiter = 40
          ! numero massime iterazioni 

        REAL(DP)  :: etresh =0.1d-6 
          ! treshhold su energia 

        REAL(DP) :: passop =0.3 
          ! passetto per ricerca minimo

        INTEGER  :: epol = 3
          ! direzione campo elettrico

        REAL(DP) :: efield =0.d0 
          ! electric field intensity in atomic units


        NAMELIST / electrons / emass, emass_cutoff, orthogonalization, &
          electron_maxstep, ortho_eps, ortho_max, electron_dynamics,   &
          electron_damping, electron_velocities, electron_temperature, &
          ekincw, fnosee, ampre, grease, empty_states_nbnd,            &
          empty_states_maxstep, empty_states_delt, empty_states_emass, &
          empty_states_ethr, diis_size, diis_nreset, diis_hcut,        &
          diis_wthr, diis_delt, diis_maxstep, diis_rot, diis_fthr,     &
          diis_temp, diis_achmix, diis_g0chmix, diis_g1chmix,          &
          diis_nchmix, diis_nrot, diis_rothr, diis_ethr, diis_chguess, &
          mixing_mode, mixing_beta, mixing_ndim, mixing_fixed_ns,      &
          diago_cg_maxiter, diago_david_ndim, diagonalization,         &
          startingpot, startingwfc , conv_thr, diago_diis_ndim,        &
          diago_thr_init, n_inner, fermi_energy, rotmass, occmass,     &
          rotation_damping, occupation_damping, rotation_dynamics,     &
          occupation_dynamics, tcg, maxiter, etresh, passop, epol,     &
          efield


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
          ! ion_dynamics = 'sd' | 'cg' | 'damp' | 'verlet' | 'bfgs' | 
          !                'old-bfgs' | 'none'*  
          ! set how ions shold be moved
          ! 'none'     ions are kept fixed 
          ! 'bfgs'     a new BFGS algorithm is used to minimize ionic configuration
          ! 'old-bfgs' old BFGS algorithm is used to minimize ionic configuration
          ! 'sd'       steepest descent algorithm is used to minimize ionic configuration
          ! 'cg'       conjugate gradient algorithm is used to minimize ionic configuration
          ! 'damp'     damped dynamics is used to propagate ions
          ! 'verlet'   standard Verlet algorithm is used to propagate ions

        CHARACTER(LEN=80) :: ion_dynamics_allowed(10)
        DATA ion_dynamics_allowed / 'sd', 'cg', 'damp', 'verlet', 'none',   &
                                    'bfgs', 'old-bfgs', 'constrained-damp', &
                                    'constrained-verlet', 'beeman' /

        REAL(DP) :: ion_radius(nsx) = 0.5d0
          ! pseudo-atomic radius of the i-th atomic species
          ! (for Ewald summation), values between 0.5 and 2.0 are usually used.

        REAL(DP) :: ion_damping = 0.2d0
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
          ! ion_temperature = 'nose' | 'not_controlled'* | 'rescaling'
          ! 'nose'           control ionic temperature using Nose thermostat
          !                  see parameters "fnosep" and "tempw"
          ! 'rescaling'      control ionic temperature via velocities rescaling 
          !                  see parameter "tolp"
          ! 'not_controlled' ionic temperature is not controlled

        REAL(DP) :: tempw = 300.0d0
          ! meaningful only with "ion_temperature /= 'not_controlled' "
          ! value of the ionic temperature (in Kelvin) forced
          ! by the temperature control

        REAL(DP) :: fnosep( nhclm )  = 50.0d0
          ! meaningful only with "ion_temperature = 'nose' "
          ! oscillation frequency of the nose thermostat (in terahertz) 
          ! nhclm is a max length for the chain

        INTEGER   ::  nhpcl = 0
          ! non-zero only with "ion_temperature = 'nose' "
          ! this defines the length of the Nose-Hoover chain

        INTEGER   :: nhptyp = 0
        ! this parameter set the nose hoover thermostat to more than one

        INTEGER   ::  ndega = 0
          ! this is the parameter to control active degrees of freedom
          ! used for temperature control and the Nose-Hoover chains

        REAL(DP) :: tolp = 50.0d0
          ! meaningful only with "ion_temperature = 'rescaling' "
          ! tolerance (in Kelvin) of the rescaling. When ionic temperature
          ! differs from "tempw" more than "tolp" apply rescaling.

        LOGICAL   :: tranp(nsx) = .FALSE.
          ! tranp(i) control the randomization of the i-th atomic specie
          ! .TRUE.   randomize ionic positions ( see "amprp" )
          ! .FALSE.  do nothing

        REAL(DP) :: amprp(nsx) = 0.0d0
          ! amprp(i) meaningful only if "tranp(i) = .TRUE.", amplitude of the
          ! randomization ( allowed values: 0.0 - 1.0 ) for the i-th atomic specie.
          ! Add to the positions a random displacements vector ( in bohr radius )
          ! defined as:  amprp( i ) * ( X, Y, Z )
          ! where X, Y, Z are pseudo random number in the interval [ -0.5 , 0.5 ]

        REAL(DP) :: greasp = 0.0d0
          ! same as "grease", for ionic damped dynamics
          ! NOT used in FPMD

        INTEGER   :: ion_nstepe = 1
          ! number of electronic steps for each ionic step

        INTEGER   :: ion_maxstep = 1000
          ! maximum number of step in ionic minimization

        REAL(DP) :: upscale = 0.0d0
          ! This variable is NOT used in FPMD

        CHARACTER(LEN=80) :: pot_extrapolation = 'default', &
                             wfc_extrapolation = 'default'
          !  These variables are used only by PWSCF
          ! NOT used in FPMD
        
        !
        ! ... delta_T, nraise are used to change the temperature in PWscf
        !
        
        REAL(DP) :: delta_t = 1.D0
        
        INTEGER :: nraise
        
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
        
        INTEGER :: num_of_images = 0
        
        INTEGER :: init_num_of_images = 0

        CHARACTER(LEN=80) :: CI_scheme = 'no-CI' 
          ! CI_scheme = 'no-CI' | 'highest-TS' | 'all-SP' | 'manual'
          ! set the Climbing Image scheme
          ! 'no-CI'       Climbing Image is not used
          ! 'highest-TS'  Standard Climbing Image
          ! 'manual'      the image is selected by hand
        
        CHARACTER(LEN=80) :: CI_scheme_allowed(4)
        DATA CI_scheme_allowed / 'no-CI', 'highest-TS', 'all-SP', 'manual' /
        
        LOGICAL :: first_last_opt = .FALSE.
        
        LOGICAL :: use_masses     = .FALSE.
        
        LOGICAL :: use_fourier    = .FALSE.
        
        LOGICAL :: use_multistep  = .FALSE.
        
        LOGICAL :: use_freezing   = .FALSE.
        
        LOGICAL :: fixed_tan      = .FALSE.
        
        LOGICAL :: free_energy    = .FALSE.
                
        LOGICAL :: write_save     = .FALSE.
        
        CHARACTER(LEN=80) :: opt_scheme = 'quick-min' 
          ! minimization_scheme = 'quick-min' | 'damped-dyn' | 
          !                       'mol-dyn'   | 'sd'
          ! set the minimization algorithm
          ! 'quick-min'   projected molecular dynamics
          ! 'damped-dyn'  damped molecular dynamics
          ! 'mol-dyn'     constant temperature molecular dynamics
          ! 'sd'          steepest descent
          ! 'broyden'     broyden acceleration
          ! 'langevin'    langevin dynamics

        CHARACTER(LEN=80) :: opt_scheme_allowed(6)
        DATA opt_scheme_allowed / 'quick-min', 'damped-dyn', 'broyden', &
                                  'mol-dyn', 'sd', 'langevin' /  

        REAL (DP)  :: damp = 1.D0
          ! meaningful only when minimization_scheme = 'damped-verlet'
        
        REAL (DP)  :: temp_req = 0.D0
          ! meaningful only when minimization_scheme = 'sim-annealing'

        REAL (DP)  :: ds = 1.5D0

        REAL (DP)  :: k_max = 0.1D0, k_min = 0.1D0

        REAL (DP)  :: path_thr = 0.05D0

        !
        ! ... variables added for new BFGS algorithm
        !
        
        INTEGER ::  bfgs_ndim = 1
                            
        REAL(DP)  :: trust_radius_max = 0.8D0
        REAL(DP)  :: trust_radius_min = 1.D-3
        REAL(DP)  :: trust_radius_ini = 0.5D0
                               
        REAL(DP)  :: w_1 = 0.5D-1
        REAL(DP)  :: w_2 = 0.5D0 

        REAL(DP)  :: sic_rloc = 0.0d0

        !
        ! ...  variable added for SMD  ( Y.K. 04/15/2004 )
        !
  
        LOGICAL :: smd_polm = .FALSE.
          ! polynomial interpolation for the initial path.
        
        INTEGER :: smd_kwnp = 2
          ! # of known points for polm. interp.
  
        LOGICAL :: smd_linr = .FALSE.
          ! linear interpolation for the initial path. 
  
        LOGICAL :: smd_stcd = .FALSE.
          ! for molecular reactions, make the 6 degrees of
          ! the same to prevent rotation and translation
          ! at the initial state.  
  
        INTEGER :: smd_stcd1, smd_stcd2, smd_stcd3  
          ! for STCD = .true.
  
        INTEGER :: smd_codf = 50, smd_forf = 50, smd_smwf = 1
          ! Frequency of writing coord, force, replica
          ! files
  
        INTEGER :: smd_lmfreq = 1
          ! Frequencey of computing the Lag. mul.
  
        REAL (DP) :: smd_tol = 0.0001
          ! Tolerance in constraints in units of
          ! [alpha(k)-alpha(k-1)] - 1/SM_P      
          !
        
        INTEGER :: smd_maxlm = 10
          ! Max number of iterations for calculating
          ! Lag. mul.
  
        LOGICAL :: smd_smcp = .TRUE.
        LOGICAL :: smd_smopt = .FALSE.
        LOGICAL :: smd_smlm = .FALSE.
        LOGICAL :: smd_splc = .FALSE.
         ! smcp  for regular elec min.
         ! smopt for the minimization of initial & final
         !       state
         ! smlm  for SMD.
         ! splc  for future sp locater using SMD.
  
        REAL (DP) :: smd_ene_ini = 1.D0
        REAL (DP) :: smd_ene_fin = 1.D0
        !
        REAL (DP) :: smd_spal = 1.D0
        
        !
        ! ... variable for meta-dynamics
        !
        
        INTEGER  :: fe_nstep    = 50
        INTEGER  :: shake_nstep = 5
        REAL(DP) :: g_amplitude = 0.01D0
        !
        REAL(DP) :: fe_step( max_nconstr ) = 0.4D0
        !
        NAMELIST / ions / phase_space, ion_dynamics, ion_radius, ion_damping,  &
                          ion_positions, ion_velocities, ion_temperature,      &
                          tempw, fnosep, nhpcl, nhptyp, ndega, tranp, amprp,   &
                          greasp, tolp, ion_nstepe, ion_maxstep, upscale,      &
                          delta_t, pot_extrapolation, wfc_extrapolation,       &
                          nraise, num_of_images, CI_scheme, opt_scheme,        &
                          use_masses, first_last_opt, use_multistep,           &
                          write_save, damp, temp_req, ds, k_max, k_min,        &
                          path_thr, fixed_tan, init_num_of_images,             &
                          free_energy, use_freezing, use_fourier,              &
                          trust_radius_max, trust_radius_min,                  &
                          trust_radius_ini, w_1, w_2, bfgs_ndim, sic_rloc,     &
                          smd_polm, smd_kwnp, smd_linr, smd_stcd, smd_stcd1,   &
                          smd_stcd2, smd_stcd3, smd_codf, smd_forf, smd_smwf,  &
                          smd_lmfreq, smd_tol, smd_maxlm, smd_smcp, smd_smopt, &
                          smd_smlm, smd_ene_ini, smd_ene_fin,                  &
                          fe_step, fe_nstep, shake_nstep, g_amplitude

!=----------------------------------------------------------------------------=!
!  CELL Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!

        CHARACTER(LEN=80) :: cell_parameters = 'default' 
          ! cell_parameters = 'default'* | 'from_input' 
          ! 'default'    restart the simulation with cell parameters read
          !              from the restart file or "celldm" if "restart = 'from_scratch'"
          ! 'from_input' restart the simulation with cell parameters
          !              from standard input ( see the card 'CELL_PARAMETERS' )

        CHARACTER(LEN=80) :: cell_dynamics  = 'none'
          ! cell_dynamics = 'sd' | 'pr' | 'none'*  
          ! set how cell shold be moved
          ! 'none'   cell is kept fixed 
          ! 'sd'     steepest descent algorithm is used to minimize the cell
          ! 'pr'     standard Verlet algorithm is used to propagate the cell

        CHARACTER(LEN=80) :: cell_dynamics_allowed(6)
        DATA cell_dynamics_allowed / 'sd', 'pr', 'none', 'w', 'damp-pr', 'damp-w'  /

        CHARACTER(LEN=80) :: cell_velocities = 'default'
          ! cell_velocities = 'zero' | 'default'*
          ! 'zero'    restart setting cell velocitiy to zero
          ! 'default' restart using cell velocity of the previous run

        REAL(DP) :: press = 0.0d0
          ! external pressure (in GPa, remember 1 kbar = 10^8 Pa)

        REAL(DP) :: wmass = 0.0d0
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

        REAL(DP) :: temph = 0.0d0
          ! meaningful only with "cell_temperature /= 'not_controlled' "
          ! value of the cell temperature (in Kelvin) forced
          ! by the temperature control

        REAL(DP) :: fnoseh = 1.0d0
          ! meaningful only with "cell_temperature = 'nose' "
          ! oscillation frequency of the nose thermostat (in terahertz) 

        REAL(DP) :: greash = 0.0d0
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

        REAL(DP) :: cell_factor = 0.0d0
          ! NOT used in FPMD

        INTEGER   :: cell_nstepe = 1
          ! number of electronic steps for each cell step

        REAL(DP) :: cell_damping = 0.0d0
          ! meaningful only if " cell_dynamics = 'damp' "
          ! damping frequency times delta t, optimal values could be
          ! calculated with the formula
          !        sqrt(0.5*log((E1-E2)/(E2-E3))) 
          ! where E1 E2 E3 are successive values of the DFT total energy 
          ! in a ionic steepest descent simulation

        NAMELIST / cell / cell_parameters, cell_dynamics, cell_velocities, &
                          press, wmass, cell_temperature, temph, fnoseh,   &
                          cell_dofree, greash, cell_factor, cell_nstepe,   &
                          cell_damping
!
!=----------------------------------------------------------------------------=!
!  PHONON Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!

        INTEGER :: modenum = 0
        
        INTEGER :: nq1, nq2, nq3
        ! number of q points in each direction
        
        REAL(DP) :: xqq(3) = 0.0d0
          ! coordinates of q point for phonon calculation

        REAL(DP) :: tr2_ph

        NAMELIST / phonon / modenum, xqq, nq1, nq2, nq3, tr2_ph

!=----------------------------------------------------------------------------=!
!  RAMAN Namelist Input Parameters
!=----------------------------------------------------------------------------=!
 
         REAL(DP) :: b_length
           ! length of the b-vector
         LOGICAL :: lcart
           ! cartesian directions
 
         NAMELIST / raman / b_length, lcart

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
        CHARACTER(LEN=3)  :: atom_label(nsx) = 'XX'  ! label of the atomic species being read
        CHARACTER(LEN=80) :: atom_pfile(nsx) = 'YY'  ! pseudopotential file name
        REAL(DP)         :: atom_mass(nsx)  = 0.0d0 ! atomic mass of the i-th atomic species
          ! in atomic mass units: 1 a.m.u. = 1822.9 a.u. = 1.6605 * 10^-27 kg
        LOGICAL   :: taspc = .FALSE.

!
!    ATOMIC_POSITIONS
!
        REAL(DP) :: rd_pos(3,natx) = 0.D0 ! unsorted position from input
        INTEGER  :: sp_pos(natx)   = 0
        INTEGER  :: if_pos(3,natx) = 1
        INTEGER  :: id_loc(natx)   = 0
        INTEGER  :: na_inp(nsx)    = 0    ! number of atom for each specie
        LOGICAL  :: tapos = .FALSE.
        CHARACTER(LEN=80) :: atomic_positions = 'crystal'
          ! atomic_positions = 'bohr' | 'angstrong' | 'crystal' | 'alat'
          ! select the units for the atomic positions being read from stdin

        !
        ! ... variable added for NEB  ( C.S. 17/10/2003 )
        !
        REAL(DP) :: pos( 3*natx, max_num_of_images ) = 0.D0

!
!    ION_VELOCITIES
!
        REAL(DP) :: rd_vel(3,natx) = 0.0d0   ! unsorted velocities from input
        INTEGER  :: sp_vel(natx)   = 0 
        LOGICAL  :: tavel          = .FALSE.

!
!    KPOINTS
!
! ...   k-points inputs
        LOGICAL :: tk_inp = .FALSE.
        REAL(DP) :: xk(3,npkx) = 0.0d0, wk(npkx) = 0.0d0
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
        REAL(DP) :: f_inp(nbndxx, nspinx) = 0.0d0
        LOGICAL   :: tf_inp = .FALSE.

!
!    VHMEAN
!
! ...   card planar mean of the Hartree potential
        LOGICAL :: tvhmean_inp = .FALSE.
        INTEGER :: vhnr_inp = 0, vhiunit_inp = 0
        REAL(DP)  :: vhrmin_inp = 0.0d0, vhrmax_inp = 0.0d0
        CHARACTER :: vhasse_inp = 'X'

!
!    OPTICAL
!
        LOGICAL :: toptical_card = .FALSE.
        REAL(DP) :: woptical = 0.0d0, boptical = 0.0d0
        INTEGER   :: noptical = 0

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
       REAL(DP) :: neighbo_radius = 0.0d0

!
!    CELL_PARAMETERS
!
       REAL(DP) :: rd_ht(3,3) = 0.0d0
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
      INTEGER   :: nconstr_inp               = 0
      REAL (DP) :: constr_tol_inp            = 0.D0
      INTEGER   :: constr_type_inp(natx)     = 0
      REAL (DP) :: constr_inp(4,natx)        = 0
      REAL (DP) :: constr_target(natx)       = 0.D0
      LOGICAL   :: constr_target_set(natx)   = .FALSE.

!
!    KOHN_SHAM
!
      LOGICAL :: tprnks( nbndxx, nspinx ) = .FALSE.
        ! logical mask used to specify which kohn sham orbital should be
        ! written to files 'KS.'
      LOGICAL :: tprnks_empty( nbndxx, nspinx ) = .FALSE.
        ! logical mask used to specify which empty kohn sham orbital should be
        ! written to files 'KS_EMP.'
      CHARACTER(LEN=256) :: ks_path = './'
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
      LOGICAL :: climbing( max_num_of_images ) = .FALSE.

!
!   PLOT_WANNIER
!

      INTEGER, PARAMETER :: nwf_max = 1000
      !
      INTEGER :: wannier_index(nwf_max)

!  END manual
! ----------------------------------------------------------------------
  !
!=----------------------------------------------------------------------------=!
!
END MODULE input_parameters
!
!=----------------------------------------------------------------------------=!
