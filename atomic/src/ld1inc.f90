!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ld1inc
  USE kinds, ONLY : dp
  USE ld1_parameters
  USE radial_grids, ONLY: radial_grid_type, ndmx
  USE paw_type, ONLY : paw_t
  IMPLICIT NONE
  SAVE
  PRIVATE :: nwfx, nwfsx, ncmax1
  INTEGER, PARAMETER :: lmx=3, lmx2=2*lmx
  !
  !    variables for the all-electron calculation
  !
  INTEGER  ::      &
       nn(nwfx),   &   ! the main quantum number
       ll(nwfx),   &   ! the orbital angular momentum
       nwf,        &   ! the number of wavefunctions
       isw(nwfx),  &   ! spin of the wfc. if(.not.lsd) all 1 (default)
       nspin           ! 1 (default) or 2 (if lsd=true)

  CHARACTER(len=2) ::&
       el(nwfx)        !  the label of the states

  real(DP) ::   &
       jj(nwfx),     & ! the total angular momentum
       oc(nwfx),     & ! the occupations of the all-electron atom
       zed,          & ! the ionic charge
       enne,         & ! the number of electrons
       sl3(0:lmx2,0:lmx2,0:lmx2)

  real(DP)::          &
       enl(nwfx),          & ! the energies of the all-electron atom
       psi(ndmx,2,nwfx),    & ! the all-electron (dirac) wavefunctions
                             ! psi(:,1,n) = major component for state n
                             ! psi(:,2,n) = minor component for state n
       rho(ndmx,2),         & ! the all-electron density
                             ! rho(:,1) = spin-up, rho(:,2) = spin-down
       zeta(ndmx),           & ! the all-electron magnetization
       ! relativistic perturbative terms
       evel(nwfx),       & ! p^4 ("velocity") correction
       edar(nwfx),       & ! Darwin term
       eso(nwfx)           ! spin-orbit splitting

  LOGICAL :: &
       core_state(nwfx)   ! if true the state is in the core
  !
  !    the parameters of the logarithmic mesh
  !
  TYPE(radial_grid_type) :: grid
  !
  !    the variables for computing logarithmic derivatives
  !
  INTEGER :: &
       nld,  &  ! computes the log der of the last nld wavefunctions
       npte     ! number of energy points

  real(DP) :: &
       rlderiv,    & ! the radius of logarithmic derivatives
       eminld,     & ! the minimum energy
       emaxld,     & ! the maximum energy
       deld,       & ! the deltae of energy
       rpwe          ! the radius of the partial wave expansion
  !
  !   the variables which define the pseudopotential
  !
  INTEGER ::       &
       nns(nwfsx), & ! the main quantum number of pseudopotential
       lls(nwfsx), & ! the angular momentum of pseudopotential
       isws(nwfsx),& ! the spin of each pseudo-wavefunctions (not used)
       ikk(nwfsx), & ! the maximum ik of each beta functions
       ik(nwfsx),  & ! the ik that correspond to rcut
       ikus(nwfsx), & ! the ik that corresponds to rcutus
       nwfs,       & ! the number of pseudo wavefunctions
       nbeta,      & ! the number of projectors
       nsloc,      & ! the wavefunction which correspond to the loc pot
       lloc,       & ! the l component considered as local
       pseudotype, &  ! the type of pseudopotential
       nstoae(nwfsx)  ! for each pseudo the all-electron

  CHARACTER(len=2) :: &
       els(nwfsx)       !  the label of the states

  real(DP) ::       &
       enls(nwfsx),      & ! the energies of the pseudo atom
       jjs(nwfsx),       & ! the j of each wavefunction (only rel=2)
       ocs(nwfsx),       & ! the occupations of the pseudo atom
       rcut(nwfsx),      & ! the cut-off radius for pseudowavefunctions
       rcutus(nwfsx),    & ! the cut-off radius for us-pseudowavefunctions
       rcloc,            & ! cut-off for local potential
       ecutrho,          & ! suggested cut-off for the change
       ecutwfc,          & ! suggested cut-off for the wavefunctions
       zval,             & ! the ionic pseudo charge
       phis(ndmx,nwfsx),  & ! the pseudo wavefunctions
       psipsus(ndmx,nwfx),& ! the all-electron wavefunctions for us pseudo
       rhos(ndmx,2),      & ! the pseudo density
       zetas(ndmx),       & ! the pseudo magnetization
       vnl(ndmx,0:3,2),   & ! the pseudopotential in semilocal form
       betas(ndmx,nwfsx), & ! the projector functions
       chis(ndmx,nwfsx),  & ! auxiliary functions
       rho0,             & ! value of the charge at the origin
       bmat(nwfsx,nwfsx), &! the pseudo coefficients (unscreened D)
       ddd(nwfsx,nwfsx,2),&! the screened D
       qq(nwfsx,nwfsx),   &! the integrals of the qvan
       qvan(ndmx,nwfsx,nwfsx), & ! the augmentation functions
       qvanl(ndmx,nwfsx,nwfsx,0:lmx2) ! the augmentation functions, l dependent

  LOGICAL :: &
       tm,            &!  if true use Troullier-Martins for norm-conserving PP
       new(nwfsx)      !  if true the fit is on arbitrary energy
  !
  !    the variable for multiconfigurations
  !
  INTEGER ::                 &
       nconf,                & ! number of configuration
       nstoaec(nwfsx,ncmax1),& ! correspondence all-electron test
       lsdts(ncmax1),        & ! for each configuration the lsd
       nwftsc(ncmax1),       & ! number of wavefunctions for each config.
       nntsc(nwfsx,ncmax1),lltsc(nwfsx,ncmax1),& ! the quantum numbers of
                                ! each configuration
       iswtsc(nwfsx,ncmax1)    ! the spin index

  CHARACTER(len=2) ::  &
       eltsc(nwfsx,ncmax1)     !  the labels for each configuration

  real(DP) ::              &
       rcuttsc(nwfsx,ncmax1),   & ! the cut-off radius of each configuration
       rcutustsc(nwfsx,ncmax1), & ! cut-off radius for us
       jjtsc(nwfsx,ncmax1),     & ! the j of a configuration
       octsc(nwfsx,ncmax1),     & ! the occupations of each configuration
       enltsc(nwfsx,ncmax1)       ! the energies of each configuration
  !
  ! for tests
  !
  INTEGER ::        &
       nnts(nwfsx),  &   ! the main quantum number of pseudopotential
       llts(nwfsx),  &   ! the angular momentum of pseudopotential
       iswts(nwfsx), &   ! spin of the wfc. if(.not.lsd) all 1 (default)
       nstoaets(nwfsx), & ! for each test wavefunction the all-electron
       nwfts             ! the number of pseudo wavefunctions

  real(DP) ::        &
       enlts(nwfsx),       & ! the energies for the test configuration
       phits(ndmx,nwfsx),   & ! the pseudo wavefunctions
       rcutts(nwfsx),      & ! cut-off radius for test wavefunction
       rcutusts(nwfsx),    & ! us cut-off radii for test wavefunct.
       jjts(nwfsx),        & ! jj of the test function (rel=2)
       octs(nwfsx)           ! the occupation numbers

  CHARACTER(len=2) ::  &
       elts(nwfsx)           ! the label of the states
  !
  ! for LDA-1/2
  !
  real(DP) :: rcutv ! CUT for LDA-1/2
  !
  !    The control of the run
  !
  INTEGER ::      &
       iter,      &  ! iteration conter
       lsd,       &  ! if true lsd calculation
       isic,      &  ! if true uses self-interaction correction
       latt,      &  ! if true Latter's correction is applied
       iswitch,   &  ! control the type of run
       max_out_wfc, & ! maximum number of wavefunctions written on the
                     !   wavefunctions file.
       rel           ! 0 nonrelativistic calculation
  ! 1 scalar relativistic calculation
  ! 2 calculation with the full dirac equation
  LOGICAL ::      &
       lsmall,    &     ! if true writes the small component on file
       relpert,   &     ! compute relativistic perturbative corrections
       frozen_core, &   ! if true the all-electron calculation is frozen core
       write_coulomb, & ! if true write a fake UPF pseudopotential file named
                        ! X.coul (X=atomic symbol) - for usage in special cases
                        ! when the bare coulomb potential is required
       noscf            ! if true a hydrogenic atom is solved. The charge
                        ! density is not computed.



  CHARACTER(len=4) :: &
       verbosity     ! if 'high' writes more information on output


  real(DP) :: &
       beta,       &   ! the mixing parameter
       tr2,        &   ! the required precision of the scf
       eps0            ! the reached precision of the scf
  !
  !    parameters for the old type pseudopotential
  !
  INTEGER ::   &
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
       ekinc,      &    ! core kinetic energy
       ekinc0,     &    ! core kinetic energy
       ekinv,      &    ! valence kinetic energy
       enclv, enclc,  & ! nuclear Coulomb energy of valence and core
       ehrtvv,     &    ! valence-valence Hartree energy
       ehrtcv,     &    ! core-valence Hartree energy
       ehrtcc,     &    ! core-core Hartree energy
       ae_fc_energy, &  ! frozen core energy calculated with all-electron char
       dhrsic,     &    ! Hartree sic energy
       dxcsic,     &    ! exchange sic energy
       etots,      &    ! total pseudopotential energy
       etots0           ! saved value of the total pseudopotential energy
  !
  !  variable for nlcc
  !
  real(DP) :: &
       rcore,      &  ! the points where core charge is smooth
       rhoc(ndmx)      ! the core charge

  LOGICAL :: &
       new_core_ps, & ! if true pseudize the core charge with bessel functions
       nlcc    ! if true nlcc pseudopotential
  !
  !  the potential for the scf
  !
  real(DP) ::   &
       v0(ndmx),      & ! the coulomb potential
       vpot(ndmx,2),  & ! the all-electron scf potential
       vxt(ndmx),     & ! the external potential
       vh(ndmx),      & ! the hartree potential
       vxc(ndmx,2),   & ! the exchange and correlation potential
       exc(ndmx),     & ! the exchange and correlation energy
       excgga(ndmx),  & ! the GGA exchange and correlation energy
       vxcts(ndmx,2), & ! the pseudo exchange and correlation potential
       excts(ndmx),   & ! the pseudo exchange and correlation energy
       excggats(ndmx),& ! the GGA exchange and correlation energy
       vpstot(ndmx,2),& ! the total local pseudopotential
       vpsloc(ndmx)  ,& ! the local pseudopotential
       vx(ndmx,2)    ,& ! the OEP-X potential (when needed)
       enzero(2)
  real(DP) ::  &
       tau(ndmx,2),   & ! kinetic energy density for metaGGA
       vtau(ndmx)       ! potential for metaGGA
  real(DP), ALLOCATABLE ::  &
       vsic(:,:), vsicnew(:), vhn1(:), egc(:) ! potentials for SIC
  !
  LOGICAL :: lsave_wfc  ! if true, wfcs (AE and PS) are saved to the UFP file
  !
  !  variables needed for PAW dataset generation and test
  !
  LOGICAL :: &
       lpaw,      &! if true generate or test a PAW dataset
       lnc2paw, &  ! if true the PAW dataset is generate from the NC one
       rmatch_augfun_nc, &  ! if .true. the norm conserving radii are
                            ! used to pseudize the q functions
       use_paw_as_gipaw ! if true, PAW data will be used for GIPAW
  TYPE(paw_t) :: &
       pawsetup    ! the PAW dataset
  real(DP) ::       &
       rmatch_augfun,     & ! define the matching radius for paw aug.fun.
       psipaw(ndmx,nwfsx),& ! the all-electron wavefunctions for any beta
       psipaw_rel(ndmx,nwfsx),& ! the all-electron wfc small component
       aeccharge(ndmx),   & ! true, not smoothened, AE core charge for PAW
       psccharge(ndmx),   & ! smoothened core charge for PAW
       paw_energy(5,3)

   CHARACTER(len=20) ::&
       which_augfun     ! choose shape of paw fun. (GAUSS, BESSEL..)
  !
  ! conversion factor
  !
  real(DP) :: &
             rytoev_fact    ! Conversion from Ry and eV. A value
                            ! different from default can be used
                            ! to reproduce results of old papers.
  real(DP) :: &
             cau_fact       ! speed of light in atomic units.
  !
  !  Auxiliary quantities for verbose output
  !
  real(DP) ::       &
       aevcharge(ndmx,2)     ! the all-electron valence charge

  !
  !  file names
  !
  CHARACTER(len=75)  :: title  ! the title of the run
  CHARACTER(len=75)  :: author ! the author of the pseudopotential
  CHARACTER(len=240) :: prefix ! prefix for file names
  CHARACTER(len=256) ::      & !
       file_pseudo,          & ! input file containing the pseudopotential
       file_pseudopw           ! output file where the pseudopot is written
  CHARACTER(len=256) ::      & ! output filenames read from input, containing:
       file_charge,          & ! all-electron total charge only
       file_chi,             & ! chi functions
       file_beta,            & ! beta functions
       file_qvan,            & ! qvan functions
       file_screen,          & ! screening potential
       file_core,            & ! core charge
       file_recon              ! information for paw reconstruction
  ! the following filenames are determined by "prefix", not read from input
  CHARACTER(len=256) ::      & ! output files, containing:
       file_wfcaegen,        & ! all-electron wavefunctions for generation
       file_wfcncgen,        & ! norm-conserving wavefunctions for generation
       file_wfcusgen,        & ! ultra-soft wavefunctions for generation
       file_potscf,          & ! scf potential at each iteration
       file_wavefunctions,   & ! all-electron results for orbitals
       file_wavefunctionsps, & ! pseudopotential results for orbitals
       file_logder,          & ! all-electron logarithmic derivatives
       file_logderps,        & ! pseudopotential logarithmic derivatives
       file_pawexp,          & ! quality index of partial wave expansion
       file_tests              ! results of pseudopotential tests
  !
  ! vdw calculation
  !
  LOGICAL :: vdw        ! optional variable
  !
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
  ! (GI)PAW reconstruction
  !
  LOGICAL :: lgipaw_reconstruction
  REAL ( dp ) :: wfc_ae_recon(ndmx,nwfx)
  REAL ( dp ) :: wfc_ps_recon(ndmx,nwfsx)
  REAL ( dp ) :: wfc_us_recon(ndmx,nwfsx)
  !
END MODULE ld1inc
