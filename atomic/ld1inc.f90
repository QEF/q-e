!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module ld1inc
  use kinds, only : dp
  use ld1_parameters
  use atomic_paw, only : paw_t
  integer, parameter :: lmx=3, lmx2=2*lmx
  !
  !    variables for the all-electron calculation
  !
  integer  ::      &
       nn(nwfx),   &   ! the main quantum number
       ll(nwfx),   &   ! the orbital angular momentum
       nwf,        &   ! the number of wavefunctions
       isw(nwfx),  &   ! spin of the wfc. if(.not.lsd) all 1 (default)
       nspin           ! 1 (default) or 2 (if lsd=true)

  character(len=2) ::&
       el(nwfx)        !  the label of the states

  real(DP) ::   &
       jj(nwfx),     & ! the total angular momentum
       oc(nwfx),     & ! the occupations of the all-electron atom
       zed,          & ! the ionic charge 
       enne,         & ! the number of electrons
       sl3(0:lmx2,0:lmx2,0:lmx2)

  real(DP)::          &
       enl(nwfx),          & ! the energies of the all-electron atom
       psi(ndm,2,nwfx),    & ! the all-electron (dirac) wavefunctions
                             ! psi(:,1,n) = major component for state n 
                             ! psi(:,2,n) = minor component for state n
       rho(ndm,2),         & ! the all-electron density
                             ! rho(:,1) = spin-up, rho(:,2) = spin-down
       zeta(ndm)             ! the all-electron magnetization 

  logical :: &
       core_state(nwfx)   ! if true the state is in the core
  !
  !    the parameters of the logarithmic mesh
  !
  integer :: &
       mesh     ! the number of mesh points

  real(DP) :: &
       r(ndm),     & ! the radial mesh
       r2(ndm),    & ! the square of the radial mesh 
       rab(ndm),   & ! d r(x) / d x where x is the linear grid
       sqr(ndm),   & ! the square root of the radial mesh
       xmin,       & ! the minimum x
       rmax,       & ! the maximum radial point
       zmesh,      & ! the ionic charge used for the mesh
       dx            ! the deltax of the linear mesh
  !
  !    the variables for computing logarithmic derivatives
  !
  integer :: &
       nld,  &  ! computes the log der of the last nld wavefunctions
       npte     ! number of energy points

  real(DP) :: &
       rlderiv,    & ! the radius of logarithmic derivatives
       eminld,     & ! the minimum energy
       emaxld,     & ! the maximum energy
       deld          ! the deltae of energy
  !
  !   the variables which define the pseudopotential
  !
  integer ::       &
       nns(nwfsx), & ! the main quantum number of pseudopotential
       lls(nwfsx), & ! the angular momentum of pseudopotential
       isws(nwfsx),& ! the spin of each pseudo-wavefunctions (not used)
       ikk(nwfsx), & ! the maximum ik of each wavefunctions
       nwfs,       & ! the number of pseudo wavefunctions 
       nbeta,      & ! the number of projectors
       nsloc,      & ! the wavefunction which correspond to the loc pot
       lloc,       & ! the l component considered as local
       pseudotype, &  ! the type of pseudopotential
       nstoae(nwfsx)  ! for each pseudo the all-electron

  character(len=2) :: &
       els(nwfsx)       !  the label of the states

  real(DP) ::       &
       enls(nwfsx),      & ! the energies of the pseudo atom
       jjs(nwfsx),       & ! the j of each wavefunction (only rel=2)
       ocs(nwfsx),       & ! the occupations of the pseudo atom
       rcut(nwfsx),      & ! the cut-off radius for pseudowavefunctions
       rcutus(nwfsx),    & ! the cut-off radius for us-pseudowavefunctions
       rcloc,            & ! cut-off for local potential
       zval,             & ! the ionic pseudo charge
       phis(ndm,nwfsx),  & ! the pseudo wavefunctions
       psipsus(ndm,nwfx),& ! the all-electron wavefunctions for us pseudo
       rhos(ndm,2),      & ! the pseudo density
       zetas(ndm),       & ! the pseudo magnetization
       vnl(ndm,0:3,2),   & ! the pseudopotential in semilocal form
       betas(ndm,nwfsx), & ! the projector functions
       chis(ndm,nwfsx),  & ! auxiliary functions
       rho0,             & ! value of the charge at the origin
       bmat(nwfsx,nwfsx), &! the pseudo coefficients (unscreened D)
       ddd(nwfsx,nwfsx,2),&! the screened D
       qq(nwfsx,nwfsx),   &! the integrals of the qvan
       qvan(ndm,nwfsx,nwfsx) ! the augmentation functions

  logical :: &
       tm,            &!  if true use Troullier-Martins for norm-conserving PP
       new(nwfsx)      !  if true the fit is on arbitrary energy
  !
  !    the variable for multiconfigurations
  !
  integer ::                 &
       nconf,                & ! number of configuration
       nstoaec(nwfsx,ncmax1),& ! correspondence all-electron test
       nwftsc(ncmax1),       & ! number of wavefunctions for each config.
       nntsc(nwfsx,ncmax1),lltsc(nwfsx,ncmax1),& ! the quantum numbers of
                                ! each configuration
       iswtsc(nwfsx,ncmax1)    ! the spin index

  character(len=2) ::  &
       eltsc(nwfsx,ncmax1)     !  the labels for each configuration

  real(DP) ::              &
       rcuttsc(nwfsx,ncmax1),   & ! the cut-off radius of each configuration
       rcutustsc(nwfsx,ncmax1), & ! cut-off radius for us
       jjtsc(nwfsx,ncmax1),     & ! the j of a configuration 
       octsc(nwfsx,ncmax1)        ! the occupations of each configuration
  !
  ! for tests
  !
  integer ::        &
       nnts(nwfsx),  &   ! the main quantum number of pseudopotential
       llts(nwfsx),  &   ! the angular momentum of pseudopotential
       iswts(nwfsx), &   ! spin of the wfc. if(.not.lsd) all 1 (default)
       nwfts             ! the number of pseudo wavefunctions

  real(DP) ::        &
       enlts(nwfsx),       & ! the energies for the test configuration
       phits(ndm,nwfsx),   & ! the pseudo wavefunctions
       rcutts(nwfsx),      & ! cut-off radius for test wavefunction
       rcutusts(nwfsx),    & ! us cut-off radii for test wavefunct.
       jjts(nwfsx),        & ! jj of the test function (rel=2)
       octs(nwfsx)           ! the occupation numbers

  character(len=2) ::  &
       elts(nwfsx)           ! the label of the states
  !
  !    The control of the run
  !
  integer ::      &
       iter,      &  ! iteration conter
       lsd,       &  ! if true lsd calculation
       isic,      &  ! if true uses self-interaction correction
       latt,      &  ! if true Latter's correction is applied
       iswitch,   &  ! control the type of run
       rel           ! 0 nonrelativistic calculation
  ! 1 scalar relativistic calculation
  ! 2 calculation with the full dirac equation

  real(DP) :: &
       beta,       &   ! the mixing parameter
       tr2,        &   ! the required precision of the scf
       eps0            ! the reached precision of the scf
  !
  !    parameters for the old type pseudopotential
  !
  integer ::   &
       lmin,   &  ! the minimum angular momentum
       lmax,   &  ! the maximum angular momentum
       nlc,    &  ! number of core functions
       nnl        ! number of angular momentum functions

  real(DP) ::     &
       cc(2),          & ! the coeffients of the core part
       alpc(2),        & ! the alpha parameters of the core
       alc(6,0:3),     & ! the coefficients of the pseudopotential
       alps(3,0:3)       ! the alpha parameters
  !
  !   the energy parameters
  !
  real(DP) :: &
       etot,       &    ! total energy
       etot0,      &    ! saved value of the total energy
       ekin,       &    ! kinetic energy
       encl,       &    ! nuclear Coulomb energy
       ehrt,       &    ! Hartree energy
       ecxc,       &    ! exchange-correlation energy
       ecc,        &    ! core-only contribution to the energy
       evxt,       &    ! external field energy 
       epseu,      &    ! pseudopotential energy
       dhrsic,     &    ! Hartree sic energy
       dxcsic,     &    ! exchange sic energy
       etots,      &    ! total pseudopotential energy
       etots0           ! saved value of the total pseudopotential energy
  !
  !  variable for nlcc
  !
  real(DP) :: &
       rcore,      &  ! the points where core charge is smooth
       rhoc(ndm)      ! the core charge

  logical :: &
       nlcc    ! if true nlcc pseudopotential
  !
  !  the potential for the scf
  !
  real(DP) ::   &
       vpot(ndm,2),  & ! the all-electron scf potential
       vxt(ndm),     & ! the external potential
       vh(ndm),      & ! the hartree potential
       vpstot(ndm,2),& ! the total local pseudopotential
       vpsloc(ndm)  ,& ! the local pseudopotential
       vx(ndm,2)    ,& ! the OEP-X potential (when needed)
       enzero(2)
  !
  !  variables needed for PAW dataset generation and test
  !
  logical :: &
       lpaw        ! if true generate or test a PAW dataset
  type(paw_t) :: &
       pawsetup    ! the PAW dataset
  real(DP) ::       &
       psipaw(ndm,nwfsx),& ! the all-electron wavefunctions for any beta
       aeccharge(ndm),   & ! true, not smoothened, AE core charge for PAW
       psccharge(ndm)      ! smoothened core charge for PAW
  !
  !  file names
  !
  character(len=75)  :: title  ! the title of the run
  character(len=240) :: prefix ! prefix for file names
  character(len=256) ::      & ! 
       file_pseudo,          & ! input file containing the pseudopotential
       file_pseudopw           ! output file where the pseudopot is written
  character(len=256) ::      & ! output filenames read from input, containing:
       file_chi,             & ! chi functions
       file_beta,            & ! beta functions
       file_qvan,            & ! qvan functions
       file_screen,          & ! screening potential
       file_core,            & ! core charge
       file_recon              ! information for paw reconstruction
  ! the following filenames are determined by "prefix", not read from input
  character(len=256) ::      & ! output files, conatining:
       file_wavefunctions,   & ! all-electron results for orbitals
       file_wavefunctionsps, & ! pseudopotential results for orbitals
       file_logder,          & ! all-electron logarithmic derivatives
       file_logderps,        & ! pseudopotential logarithmic derivatives
       file_tests              ! results of pseudopotential tests
  !
  ! vdw calculation
  !
  logical :: vdw        ! optional variable
  real(DP) :: um,     & ! maximum frequency
              du,     & ! step of frequency
              tr_s    ! threshold for scf solution of modified Sternheimer equation
  !
  ! test on ghosts and convergences with spherical Bessel functions
  !
  real(DP) :: ecutmin, & ! min kinetic energy cutoff for j_l(qr)
              ecutmax, & ! max energy cutoff
              decut,   & ! step: ecut = ecutmin, ecutmin+decut, ... , ecutmax
              rm         ! radius of the box
  !
end module ld1inc
