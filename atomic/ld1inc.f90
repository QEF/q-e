module ld1inc
  use kinds, only : dp
  use ld1_parameters
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

  real(kind=dp) ::   &
       jj(nwfx),     & ! the total angular momentum
       oc(nwfx),     & ! the occupations of the all-electron atom
       oc_old(nwfx), & ! saves the occupations of the all-electron 
       zed,          & ! the ionic charge 
       enne            ! the number of electrons

  real(kind=dp)::          &
       enl(nwfx),          & ! the energies of the all-electron atom
       psi(ndm,nwfx),      & ! the all-electron wavefunctions
       psi_dir(ndm,2,nwfx),& ! all-electron dirac wavefunctions
       rho(ndm,2),         & ! the all-electron density
       zeta(ndm)             ! the all-electron magnetization 

  logical :: &
       core_state(nwfx)   ! if true the state is in the core
  !
  !    the parameters of the logaritmic mesh
  !
  integer :: &
       mesh     ! the number of mesh points

  real(kind=dp) :: &
       r(ndm),     & ! the radial mesh
       r2(ndm),    & ! the square of the radial mesh 
       sqr(ndm),   & ! the square root of the radial mesh
       xmin,       & ! the minimum x
       rmax,       & ! the maximum radial point
       zmesh,      & ! the ionic charge used for the mesh
       dx            ! the deltax of the linear mesh
  !
  !    the variables for computing logaritmic derivatives
  !
  integer :: &
       nld,  &  ! computes the log der of the last nld wavefunctions
       npte     ! number of energy points

  real(kind=dp) :: &
       rlderiv,    & ! the radius of logaritmic derivatives
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
       nstoaes(nwfsx), &  ! for each pseudo the all-electron (save variable)
       nstoae(nwfsx)  ! for each pseudo the all-electron

  character(len=2) :: &
       els(nwfsx)       !  the label of the states

  real(kind=dp) ::       &
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
       vnl(ndm,0:3),     & ! the nonlocal pseudopotential
       vnlo(ndm,0:3,2),  & ! the nonlocal pseudopotential (spin-orbit)
       betas(ndm,nwfsx), & ! the projector functions
       chis(ndm,nwfsx),  & ! auxiliary functions
       rho0,             & ! charge value at the origin
       bmat(nwfsx,nwfsx), &! the pseudo coefficients (unscreened D)
       ddd(nwfsx,nwfsx,2),&! the screened D
       qq(nwfsx,nwfsx),   &! the integrals of the qvan
       qvan(ndm,nwfsx,nwfsx) ! the augmentation functions

  logical :: &
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

  real(kind=dp) ::              &
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

  real(kind=dp) ::        &
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

  real(kind=dp) :: &
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

  real(kind=dp) ::     &
       cc(2),          & ! the coeffients of the core part
       alpc(2),        & ! the alpha parameters of the core
       alc(6,0:3),     & ! the coefficients of the pseudopotential
       alps(3,0:3)       ! the alpha parameters
  !
  !   the energy parameters
  !
  real(kind=dp) :: &
       etot,       &    ! total energy
       etot0,      &    ! save the total energy
       ekin,       &    ! kinetic energy
       encl,       &    ! nuclear Coulomb energy
       ehrt,       &    ! Hartree energy
       ecxc,       &    ! exchange energy
       evxt,       &    ! external field energy 
       epseu,      &    ! pseudopotential energy
       dhrsic,     &    ! Hartree sic energy
       dxcsic,     &    ! exchange sic energy
       etots,      &    ! total pseudopotential energy
       etots0           ! save the total pseudopotential energy
  !
  !  variable for nlcc
  !
  real(kind=dp) :: &
       rcore,      &  ! the points where core charge is smooth
       rhoc(ndm)    ! the core charge

  logical :: &
       nlcc    ! if true nlcc pseudopotential
  !
  !  the potential for the scf
  !
  real(kind=dp) ::   &
       vpot(ndm,2),  & ! the all-electron scf potential
       vxt(ndm),     & ! the external potential
       vh(ndm),      & ! the hartree potential
       vpstot(ndm,2),& ! the total local pseudopotential
       vpsloc(ndm)     ! the local pseudopotential
  !
  !   the files and relative names
  !
  character ::                  &
       title*75,                & ! the title of the run
       file_pseudo*80,          & ! the file of input pseudopotential
       file_wavefunctions*80,   & ! the file with output wavefunctions
       file_wavefunctionsps*80, & ! the file with output pseudowavefunctions
       file_pseudopw*80,        & ! the file with output pseudopotential
       file_logderae*80,        & ! the file with all-electron log der 
       file_logderps*80,        & ! the file with pseudo log der 
       file_screen*80,          & ! the file with the screening potential
       file_core*80,            & ! the file with the core charge
       file_tests*80,           & ! the file with the tests
       file_recon*80              ! the file for paw reconstruction

end module ld1inc
