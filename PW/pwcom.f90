!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
!    Structure of commons
!
!    statically  allocated variables are dimensioned in commons
!    dynamically allocated variables are dimensioned in type declaration
!
!
module brilz
  use parameters
  !
  !    The variables needed to describe the lattice
  !
  real(kind=DP) ::                                                  &
       celldm(6),      &! dimensions of the unit cell
       at(3,3),        &! direct lattice vectors
       bg(3,3),        &! reciprocal lattice vectors
       alat,           &! cell parameter
       omega,          &! volume of the unit cell
       tpiba,          &! 2 times pi / alat
       tpiba2           ! the square of tpiba
  integer       ::                                                  &
       ibrav           ! index of the bravais lattice
  character(len=9)::                                                &
       symm_type       ! 'cubic' or 'hexagonal' when ibrav=0
end module brilz
!
module basis
  use parameters
  !
  !    The variables needed to describe the atoms in the unit cell
  !
  !
  integer     ::       &
       nat,            &! number of atoms in the unit cell
       ntyp,           &! number of different types of atoms
       natomwfc         ! number of starting wavefunctions
  integer, allocatable :: &
       ityp(:)               ! the type of each atom
  real(kind=DP), allocatable ::  &
       tau(:,:)         ! the positions of each atom
  logical         ::                                                &
       ltaucry          ! if true atom. cryst. coord. accepted in input
  character(len=6)::   &
       startingwfc ,   &! 'random' or 'atomic' or 'file'
       startingpot ,   &! 'atomic' or 'file'
       startingconfig   ! 'input' or 'file'

end module basis

!
module dynam
  use parameters
  !
  !    Variables needed for the dynamics
  !
  real(kind=8) ::                                                    &
       amass(ntypx),   &! mass of atoms
       dt,             &! time step
       temperature,    &! starting temperature
       delta_T          ! rate of thermalization
  !
  integer      ::                                                    &
       nraise           ! the frequency of temperature raising
  !
end module dynam
!
module gvect
  use parameters
  !
  !    The variables describing the reciprocal lattice vectors
  !
  integer                ::                                          &
       ngm,            &! number of g vectors
       ngm_g,          &! global number of g vectors (sum over all processors)
       ngm_l,          &! the local number of g vectors (only present processor)
       gstart,         &! first nonzero g vector
       nr1,            &! fft dimension along x
       nr2,            &! fft dimension along y
       nr3,            &! fft dimension along z
       nrx1,           &! maximum fft dimension along x
       nrx2,           &! maximum fft dimension along y
       nrx3,           &! maximum fft dimension along z
       nrxx,           &! maximum total fft
       ngl              ! number of |g| shells

  integer , allocatable :: &
       ig_l2g(:)        !"l2g" means local to global, this array convert a local
                        ! G-vector index into the global index, in other words
                        ! the index of the G-v. in the overall array of G-vectors

  integer , allocatable, target :: &
       nl(:),          &! correspondence fft <-> array of G vectors
       igtongl(:)       ! correspondence shells of G <-> G
  real(kind=DP) , allocatable, target  ::                                    &
       g(:,:),         &! coordinates of G vectors
       gg(:)            ! modulus G^2 of G vectors
                        ! G vectors are in order of increasing |G|
  real(kind=DP) ::    ecutwfc ! energy cut-off

  real(kind=DP), pointer  ::                                    &
       gl(:)            ! the modulus of g in each shell
  real (kind=DP) :: &
       gcutm,          &! cut-off for G vectors
       dual,           &! link between G of wavefunctions and charge
       ecfixed,        &!
       qcutz,          &! For the modified Ekin functional
       q2sigma          !
  complex(kind=DP), allocatable :: &
       eigts1(:,:),    &!
       eigts2(:,:),    &! the phases e^{-iG*tau_s}
       eigts3(:,:)      !
  integer, allocatable  ::         &
       ig1(:),       &!
       ig2(:),       &! the indices of G components
       ig3(:)         !
end module gvect
!
module gsmooth
  use parameters
  !
  !     the variables for the smooth mesh of the wavefunction. It can
  !     be different from the large mesh if dual > 4
  !
  !
  integer :: &
       ngms,       &! the number of smooth G vectors
       ngms_g,     &! the global number of smooth G vectors (sum over all processors)
       ngms_l,     &! the local number of smooth G vectors (only present processor)
       nr1s,       &!
       nr2s,       &! the dimension of the smooth grid
       nr3s,       &!                                              &
       nrx1s,      &! maximum dimension of the smooth grid
       nrx2s,      &! maximum dimension of the smooth grid
       nrx3s,      &! maximum dimension of the smooth grid
       nrxxs        ! the total dimension of the smooth grid
  integer, pointer :: &
       &       nls(:)   ! the correspondence  G <-> smooth mesh
  !
  logical                                                          &
       doublegrid  ! true if we use a double grid
  !
  real(kind=DP)                                                           &
       gcutms      ! the cut-off of the smooth mesh
end module gsmooth
!
module klist
  use parameters
  !
  !    The variables for the k-points
  !
  real(kind=DP) ::                                                  &
       xk(3,npk),      &! coordinates of k points
       wk(npk),        &! weight of k points
       xqq(3),         &! coordinates of q point (used with iswitch=-2)
       degauss,        &! smearing parameter
       nelec            ! number of electrons
  integer       ::                                                  &
       ngk(npk),       &! number of plane waves for each k point
       nks,            &! number of k points in this pool
       nkstot,         &! total number of k points
       ngauss          ! type of smearing technique
  logical       ::                                                  &
       lgauss,         &! if true: use gaussian broadening
       lxkcry           ! if true:k-pnts in cryst. basis accepted in input
end module klist
!
module lsda_mod
  use parameters
  !
  !    The variables needed for the lsda calculation
  !
  logical :: lsda
  real(kind=DP)  ::     &
       starting_magnetization(ntypx) ! the magnetization used to start with
  integer   ::                                           &
       nspin,           &! number of spin polarization: 2 if lsda, 1 other
       current_spin,    &! spin of the current kpoint
       isk(npk)          ! for each k-point: 1=spin up, 2=spin down
end module lsda_mod
!
module ktetra
  use parameters
  !
  !    The variables for the tetrahedron method
  !
  integer                ::                                         &
       nk1, nk2, nk3, &! the special-point grid
       k1, k2, k3,    &! the offset from the origin
       ntetra          ! number of tetrahedra
  integer, allocatable :: &
       tetra(:,:)      ! index of k-points in a given tetrahedron
  !shape (4,ntetra)
  logical              :: &
       ltetra          ! if true: use tetrahedron method
end module ktetra
!
module symme
  use parameters
  !
  !    The variables needed to describe the symmetry properties
  !
  integer   :: &
       s(3,3,48),              &! simmetry matrices
       ftau(3,48),             &! fractional translations
       nsym                     ! number of symmetries
  integer, allocatable :: &
       irt(:,:)                 ! symmetric atom for each atom and sym.op.

  logical  ::    invsym         ! if true the system has inversion symmetry
end module symme
!
module  atom
  use parameters
  !
  !    The variables needed to describe the atoms and related quantities
  !
  real(kind=DP) ::                                                        &
       zmesh(npsx),          &! the atomic charge for mesh generation
       xmin(npsx),           &! initial linear mesh point
       dx(npsx),             &! linear interval for logaritmic mesh
       r(0:ndm,npsx),        &! radial logaritmic mesh
       rab(0:ndm,npsx),      &! derivative of the radial mesh
       vnl(0:ndm,0:lmaxx,npsx),    &! non local radial potential (KB type)
       chi(0:ndm,nchix,npsx),      &! radial atomic orbitals
       oc(nchix,npsx),             &! atomic level occupation
       rho_at(0:ndm,npsx),         &! radial atomic charge density
       rho_atc(0:ndm,npsx)          ! radial core charge density
  integer     :: &
       mesh(npsx),                 &! number of mesh points
       msh(npsx),                  &! the point at rcut
       nchi(npsx),                 &! number of atomic orbitals
       lchi(nchix,npsx)             ! angular momentum of atomic orbitals
  logical     :: &
       numeric(npsx)         ! if true the potential is in numeric form
  !
end module atom
!
module pseud
  use parameters
  !
  !    The variables needed to compute the BHS pseudopotentials
  !
  real(kind=DP) ::     &
       cc(2,npsx),             &! the coefficients of the erf functions
       alpc(2,npsx),           &! the alpha of the erf functions
       zp(npsx),               &! the charge of the pseudopotential
       aps(6,0:3,npsx),        &! the a_l coefficient
       alps(3,0:3,npsx),       &! the b_l coefficient
       zv(ntypx)                ! the valence charge of the atom
  integer        ::    &
       nlc(npsx),              &! number of erf functions
       nnl(npsx),              &! number of the gaussian functions
       lmax(npsx),             &! maximum angular momentum of the pseudopot
       lloc(npsx)               ! angular momentum of the part taken as local
  logical        ::     &
       bhstype(npsx)            ! if true the parameter are from the BHS table
  !
end module pseud
!
module nl_c_c
  use parameters
  !
  !    The variable needed for the Non Linear Core Correction
  !
  real(kind=DP)     :: &
       a_nlcc(npsx),         &! the a_c coefficient of the gaussian
       b_nlcc(npsx),         &! the b_c coefficient of the gaussian
       alpha_nlcc(npsx)       ! the alpha coefficient of the gaussian
  logical           :: &
       nlcc(npsx)             ! if true the atom has nlcc
  !
end module nl_c_c

module vlocal
  use parameters
  !
  !   The variables needed for the local potential in reciprocal space
  !
  complex(kind=DP) , allocatable :: &
       strf(:,:)              ! the structure factor
  real(kind=DP), allocatable     :: &
       vloc(:,:),            &! the local potential for each atom type
       vnew(:,:)              ! V_out - V_in, needed in scf force correction
end module vlocal
!
module wvfct
  use parameters
  !
  !    The variables needed to compute the band structure
  !
  integer                  ::  &
       npwx,             &! maximum number of PW for wavefunctions
       nbndx,            &! maximum number of bands
       nbnd,             &! number of bands
       npw                ! the number of plane waves
  integer, allocatable, target :: &
       igk(:),           &! correspondence k+G <-> G
       igk_l2g(:,:)       ! correspondence local index k+G <-> global G index
                          ! see also ig_l2g
  real(kind=DP) , allocatable :: &
       et(:,:),          &! eigenvalues of the hamiltonian
       wg(:,:),          &! the weight of each k point and band
       g2kin(:)           ! kinetic energy
  complex(kind=DP) , allocatable, target :: &
       evc(:,:)           ! wavefunctions in the PW basis
  logical :: &
       gamma_only         ! if .true. only half G vectors are used
end module wvfct
!
module ener
  use parameters
  !
  !    The variables needed to compute the energies
  !
  real(kind=DP)           ::  &
       etot,           &! the total energy of the solid
       eband,          &! the band energy
       deband,         &! correction for variational energy
       ehart,          &! the hartree energy
       etxc,           &! the exchange and correlation energy
       vtxc,           &! another exchange-correlation energy
       etxcc,          &! the nlcc exchange and correlation
       ewld,           &! the ewald energy
       demet,          &! correction for metals
       ef               ! the fermi energy
end module ener
!
module force_mod
  use parameters
  !
  !    The variables for the first derivative of the energy
  !
  real(kind=DP) , allocatable :: &
       force(:,:)       ! the force on each atom
  real(kind=DP)           :: &
       sigma(3,3)       ! the stress acting on the system
  logical                 :: &
       lforce,         &! if true compute the forces
       lstres           ! if true compute the stress
end module force_mod
!
module scf
  use parameters
  !
  !    The variables needed to define the self-consistent cycle
  !
  real(kind=DP), allocatable :: &
       rho(:,:),  &! the charge density in real space
       rho_save(:,:), &! another charge density in real space
       vr(:,:),   &! the Hartree + xc potential in real space
       vltot(:),  &! the local potential in real space
       vrs(:,:),  &! the total pot. in real space (smooth grig)
       rho_core(:) ! the core charge in real space
end module scf
module workspace
  use parameters
  !
  !   additional memory needed in h_psi
  !
  complex(kind=DP) , allocatable, target :: &
       psic(:)     ! additional memory for FFT
end module workspace
!
module varie
  use parameters
  !
  !   Several variables controlling the run
  !
  real(kind=DP)  :: &
       mixing_beta,      &! the mixing parameter
       tr2,              &! the convergence threshold for potential
       upscale,          &! maximum reduction of convergence threshold
       time_max,         &! maximum allowed cpu time in sec
       ethr,             &! the convergence threshold for eigenvalues
       alpha0,           &! the mixing parameters for the extrapolation
       beta0              ! of the starting potential
  integer                    :: &
       ngm0,             &! used in mix_rho
       niter,            &! the maximum number of iteration
       nmix,             &! the number of iteration kept in the history
       imix,             &! the type of mixing (0=plain,1=TF,2=localTF)
       iprint,           &! the interval between full writing of results
       iverbosity,       &! type of printing ( 0 few, 1 all )
       david,            &! used on Davidson diagonalization
       nstep,            &! number of minimization steps
       istep,            &! current minimization step
       isolve,           &! Davidson or CG diagonalization
       iswitch,          &! general switch for the calculation type
       modenum,          &! used with iswitch=-4
       max_cg_iter,      &! maximum number of iterations in a CG di
       diis_buff,        &! dimension of the buffer in diis
       diis_start_dav,   &! number of initial davidson iterations
       order             ! type of potential updating ( see update_pot )
  !
  logical                  :: &
       loverlap,         &!if .t. the diagonalization with overlap
       lscf,             &! if true the calculation is selfconsiste
       conv_elec,        &! if true electron convergence has been reached
       conv_ions,        &! if true    ionic convergence has been reached
       nosym,            &! if true no symmetry is used
       newpseudo(npsx),  &! if true done with the new pseudopotentials
       noinv,            &! if true eliminates inversion symmetry
       diis_wfc_keep,    &! if true keeps old wfc for starting
       restart,          &! if true start for a given configuration
       reduce_io         ! if true reduce the I/O to the strict minimum
  !
end module varie
!
module relax
  use parameters
  !
  !   The variables used to control ionic relaxations
  !
  integer     :: &
       fixatom           ! last "fixatom" are kept fixed
  real(kind=DP)   :: &
       epse,             &! threshold on total energy
       epsf,             &! threshold on forces
       dtau_ref,         &! estimation of dtau
       starting_diag_threshold, &! self-explanatory
       starting_scf_threshold   ! as above
  !
end module relax
!
module cellmd
  use parameters
  !
  !   The variables used to control cell relaxation
  !
  real(kind=DP)    :: &
       press, cmass,     &! target pressure and cell mass,
       ttol,             &! tollerance for temperature rescaling
       at_old(3,3),      &! the lattice vectors at the previous ste
       omega_old,        &! the cell volume at the previous step
       cell_factor        ! maximum expected (linear) cell contraction
  ! during relaxation/MD
  integer              :: &
       nzero,            &! iteration # of last thermalization
       ntimes,           &! number of thermalization steps to be performed
       ntcheck            ! # of steps between thermalizations
  !
  logical :: lmovecell    ! used in cell relaxation
  !
  character(len=2) ::  &
       calc               ! main switch for variable cell shape MD
  ! see readin, vcsmd and/or INPUT files
  !
end module cellmd
!
module units
  use parameters
  !
  !    The units where various variables are saved
  !
  !
  integer               :: &
       iunpun,           &! unit for saving the final results
       iunwfc,           &! unit with wavefunctions
       iunat,            &! unit for saving orthogonal atomic wfcs
       iunocc,           &! unit for saving the atomic n_{ij}
       iunoldwfc,        &! unit with old wavefunctions (molecular dynamics)
       iunigk,           &! unit for saving indices
       iunres,           &! unit for the restart of the run
       nwordwfc,         &! lenght of record in wavefunction file
       nwordatwfc        ! lenght of record in atomic wfc file
  !
end module units
!
module char
  use parameters
  !
  !    The names of the atoms, of the solid and of the symmetries
  !
  character(len=75)   ::  title       ! title of the run
  character(len=20)   ::  crystal     ! type of the solid
  character(len=3 )   ::  atm(ntypx)  ! name of the type of the atoms
  character(len=2 )   ::  psd(npsx)   ! name of the pseudopotential
  character(len=45)   ::  sname(48)   ! name of the symmetries
  !
end module char
!
module filnam
  use parameters
  !
  !    The name of the files
  !
  character(len=14)   :: &
       filpun,           &! name of the punch file
       input_drho,       &! name of the file with the input drho
       output_drho        ! name of the file with the output drho
  !
end module filnam
!
module us
  use parameters
  !
  ! These parameters are needed with the US pseudopotentials
  !
  integer, parameter  :: &
       nbrx = 6,           &! maximum number of beta functions
       lqmax= 2*lmaxx+1,   &! maximum number of angular momenta of Q
       nlx  = (lmaxx+1)**2,&! maximum number of combined angular momentum
       mx   = 2*lqmax-1,   &! maximum magnetic angular momentum of Q
       nqfm = 8             ! maximum number of coefficients in Q smoothing
  real(kind=DP), parameter:: &
       dq = 0.01d0          ! space between points in the pseudopotential tab.
  !
  real(kind=DP)                :: &
       dion(nbrx,nbrx,npsx),       &! D_{mu,nu} parameters (in the atomic case)
       betar(0:ndm,nbrx,npsx),     &! radial beta_{mu} functions
       qqq(nbrx,nbrx,npsx),        &! q_{mu,nu} parameters (in the atomic case)
       qfunc(0:ndm,nbrx,nbrx,npsx),&! Q_{mu,nu}(|r|) function for |r|> r_L
       qfcoef(nqfm,lqmax,nbrx,nbrx,npsx),&! coefficients for Q in region |r|<r_L
       rinner(lqmax,npsx)! values of r_L
  integer                  :: &
       nh(npsx),             &! number of beta functions per atomic type
       nbeta(npsx),          &! number of beta functions
       kkbeta(npsx),         &! point where the beta are zero
       nqf(npsx),            &! number of coefficients for Q
       nqlc(npsx),           &! number of angular momenta in Q
       ifqopt(npsx),         &! level of q optimization
       lll(nbrx,npsx),       &! angular momentum of the beta function
       iver(3,npsx)           ! version of the atomic code
  !
  integer     ::    &
       nhm,              &! max number of different beta functions per atom
       nkb,              &! total number of beta functions, with struct.fact.
       nqxq,             &! size of interpolation table
       lpx(nlx,nlx),     &! for each limi,ljmj gives the maximum LM
       lpl(nlx,nlx,mx),  &! for each limi,ljmj gives the list of LM
       lmaxkb,           &! max angular momentum
       lqx,              &! max angular momentum + 1 for Q functions
       nqx                ! number of interpolation points
  integer, allocatable ::&
       indv(:,:),        &! correspondence of betas atomic <-> soli
       nhtol(:,:),       &! correspondence n <-> angular momentum
       nhtom(:,:)         ! correspondence n <-> magnetic angular m

  complex(kind=DP) , allocatable, target :: &
       vkb(:,:),     &! all beta functions in reciprocal space
       qgm(:)         ! complete fourier transform of Q
  real(kind=DP) , allocatable          :: &
       qq(:,:,:),    &! the q functions in the solid
       dvan(:,:,:),  &! the D functions of the solid
       deeq(:,:,:,:),&! the integral of V_eff and Q_{nm}
       becsum(:,:,:),&! the sum of bec functions
       qrad(:,:,:,:),&! radial FT of Q functi
       tab(:,:,:)! interpolation table for PPs
  real(kind=DP)   :: &
       ap(lqmax*lqmax,nlx,nlx) ! Clebsch-Gordan coefficients for spher.harm.
  !
  logical                :: &
       tvanp(npsx),       &! if true the atom is of Vanderbilt type
       okvan               ! if true at least one pseudo is Vanderbilt
  !
end module us
!
module ldaU
  use parameters
  !
  !    The quantities needed in lda+U calculations
  !
  complex(kind=DP), allocatable  :: &
       swfcatom(:,:)! orthogonalized atomic wfcs
  real(kind=DP), allocatable     :: &
       ns(:,:,:,:),   &! the occupation matrix used in h_psi
       nsnew(:,:,:,:)  ! the occupation matrix computed by at
  real(kind=DP) :: &       
       d1(3,3,48),          &! matrices for rotating spherical     
       d2(5,5,48),          &! harmonics                           
       d3(7,7,48),          &! 
       eth,                 &! the (corrected) Hubbard contribution
       Hubbard_U(ntypx),    &! the Hubbard U
       Hubbard_alpha(ntypx)  ! the Hubbard alpha (used to calculate U)
  integer           :: &                                                
       niter_with_fixed_ns, &! no. of iterations with fixed ns
       Hubbard_l(ntypx),    &! the agular momentum of Hubbard states
       Hubbard_lmax          ! maximum agular momentum of Hubbard states
  logical           :: &                                                
       lda_plus_u,          &! .true. if lda+u calculation is performed
       conv_ns              ! .true. if ns are converged
end module ldaU

module pwcom
  use constants, only: e2, degspin, rytoev, amconv, uakbar, pi, tpi, fpi
  use brilz
  use basis
  use dynam
  use gvect
  use gsmooth
  use klist
  use lsda_mod
  use ktetra
  use symme
  use atom
  use pseud
  use nl_c_c
  use vlocal
  use wvfct
  use ener
  use force_mod
  use scf
  use workspace
  use varie
  use relax
  use cellmd
  use units
  use char
  use filnam
  use us
  use ldaU
end module pwcom
!
