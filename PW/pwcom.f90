!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE basis
  !
  ! ... The variables needed to describe the atoms in the unit cell
  !
  USE ions_base,  ONLY : &
       zv,               &! the valence charge of the atom
       nat,              &! number of atoms in the unit cell
       ntyp => nsp,      &! number of different types of atoms
       tau,              &! the positions of each atom
       atm,              &! name of the type of the atoms
       ityp              ! the type of each atom
  !
  SAVE
  !
  INTEGER :: &
       natomwfc            ! number of starting wavefunctions
  CHARACTER(LEN=30) :: &   ! 'alat', 'crystal', 'angstrom', 'bohr'
       atomic_positions    ! specifies how input coordinates are given
  CHARACTER(LEN=6)  :: &
       startingwfc ,      &! 'random' or 'atomic' or 'file'
       startingpot ,      &! 'atomic' or 'file'
       startingconfig      ! 'input' or 'file'
  !
END MODULE basis
!
!
MODULE dynam
  !
  ! ... Variables needed for the dynamics
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : amass
  !
  SAVE
  !
  REAL(KIND=DP) :: &
       dt,            &! time step
       temperature,   &! starting temperature
       delta_T         ! rate of thermalization
  INTEGER :: &
       nraise          ! the frequency of temperature raising
  !
END MODULE dynam
!
!
MODULE gvect
  !
  ! ...The variables describing the reciprocal lattice vectors
  !
  USE kinds,              ONLY : DP
  USE reciprocal_vectors, ONLY : ig_l2g, sortedig_l2g
  !
  SAVE
  !
  INTEGER :: &
       ngm,           &! number of g vectors
       ngm_g,         &! global number of g vectors (sum over all processors)
       ngm_l,         &! the local number of g vectors (only present processor)
       gstart,        &! first nonzero g vector
       nr1,           &! fft dimension along x
       nr2,           &! fft dimension along y
       nr3,           &! fft dimension along z
       nrx1,          &! maximum fft dimension along x
       nrx2,          &! maximum fft dimension along y
       nrx3,          &! maximum fft dimension along z
       nrxx,          &! maximum total fft
       ngl             ! number of |g| shells
  !INTEGER, ALLOCATABLE :: &
  !     ig_l2g(:)      !"l2g" means local to global, this array convert a local
  !                    ! G-vector index into the global index, in other words
  !                    ! the index of the G-v. in the overall array of G-vectors
  INTEGER, ALLOCATABLE, TARGET :: &
       nl(:),         &! correspondence fft <-> array of G vectors
       nlm(:),        &! same for gamma point calculation
       igtongl(:)      ! correspondence shells of G <-> G
  REAL(KIND=DP), ALLOCATABLE, TARGET :: &
       g(:,:),        &! coordinates of G vectors
       gg(:)           ! modulus G^2 of G vectors
                       ! G vectors are in order of increasing |G|
  REAL(KIND=DP) :: &
       ecutwfc         ! energy cut-off
  REAL(KIND=DP), POINTER :: &
       gl(:)           ! the modulus of g in each shell
  REAL (KIND=DP) :: &
       gcutm,         &! cut-off for G vectors
       dual,          &! link between G of wavefunctions and charge
       ecfixed,       &!
       qcutz,         &! For the modified Ekin functional
       q2sigma         !
  complex(KIND=DP), ALLOCATABLE :: &
       eigts1(:,:),   &!
       eigts2(:,:),   &! the phases e^{-iG*tau_s}
       eigts3(:,:)     !
  INTEGER, ALLOCATABLE  :: &
       ig1(:),        &!
       ig2(:),        &! the indices of G components
       ig3(:)          !
  !
END MODULE gvect
!
!
MODULE gsmooth
  !
  ! ... the variables for the smooth mesh of the wavefunction. It can
  ! ... be different from the large mesh if dual > 4
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  INTEGER :: &
       ngms,        &! the number of smooth G vectors
       ngms_g,      &! the global number of smooth G vectors 
                     !  (sum over all processors)
       ngms_l,      &! the local number of smooth G vectors 
                     !  (only present processor)
       nr1s,        &!
       nr2s,        &! the dimension of the smooth grid
       nr3s,        &!                                              
       nrx1s,       &! maximum dimension of the smooth grid
       nrx2s,       &! maximum dimension of the smooth grid
       nrx3s,       &! maximum dimension of the smooth grid
       nrxxs         ! the total dimension of the smooth grid
  INTEGER, POINTER :: &
       nls(:),      &! the correspondence  G <-> smooth mesh
       nlsm(:)       ! the same for gamma point calculation
  LOGICAL :: &
       doublegrid    ! .TRUE. if we use a double grid
  REAL(KIND=DP) :: &
       gcutms        ! the cut-off of the smooth mesh
  !
END MODULE gsmooth
!
!
MODULE klist
  !
  ! ... The variables for the k-points
  !  
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npk
  !
  SAVE
  !
  REAL(KIND=DP) :: &
       xk(3,npk),      &! coordinates of k points
       wk(npk),        &! weight of k points
       xqq(3),         &! coordinates of q point (used with iswitch=-2)
       degauss,        &! smearing parameter
       nelec            ! number of electrons
  INTEGER :: &
       ngk(npk),       &! number of plane waves for each k point
       nks,            &! number of k points in this pool
       nkstot,         &! total number of k points
       ngauss           ! type of smearing technique
  LOGICAL :: &
       lgauss,         &! if .TRUE.: use gaussian broadening
       lxkcry           ! if .TRUE.:k-pnts in cryst. basis accepted in input
  !
END MODULE klist
!
!
MODULE lsda_mod
  !
  ! ... The variables needed for the lsda calculation
  !  
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx, npk
  !
  SAVE
  !
  LOGICAL :: &
       lsda
  REAL(KIND=DP) :: &
       starting_magnetization(ntypx) ! the magnetization used to start with
  INTEGER :: &
       nspin,           &! number of spin polarization: 2 if lsda, 1 other
       current_spin,    &! spin of the current kpoint
       isk(npk)          ! for each k-point: 1=spin up, 2=spin down
  !     
END MODULE lsda_mod
!
!
MODULE ktetra
  !
  ! ... The variables for the tetrahedron method
  !  
  SAVE
  !
  INTEGER :: &
       nk1, nk2, nk3,   &! the special-point grid
       k1, k2, k3,      &! the offset from the origin
       ntetra            ! number of tetrahedra
  INTEGER, ALLOCATABLE :: &
       tetra(:,:)        ! index of k-points in a given tetrahedron
                         ! shape (4,ntetra)
  LOGICAL :: &
       ltetra            ! if .TRUE.: use tetrahedron method
  !
END MODULE ktetra
!
!
MODULE symme
  !
  ! ... The variables needed to describe the symmetry properties
  !  
  SAVE
  !
  INTEGER :: &
       s(3,3,48),              &! simmetry matrices
       ftau(3,48),             &! fractional translations
       nsym                     ! number of symmetries
  INTEGER, ALLOCATABLE :: &
       irt(:,:)                 ! symmetric atom for each atom and sym.op.
  LOGICAL :: &
       invsym                   ! if .TRUE. the system has inversion symmetry
  !
END MODULE symme
!
!
MODULE atom
  !
  ! ... The variables needed to describe the atoms and related quantities
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npsx, ndmx, nchix
  !
  SAVE
  !
  REAL(KIND=DP) :: &
       zmesh(npsx),              &! the atomic charge for mesh generation
       xmin(npsx),               &! initial linear mesh point
       dx(npsx),                 &! linear interval for logaritmic mesh
       r(ndmx,npsx),             &! radial logaritmic mesh
       rab(ndmx,npsx),           &! derivative of the radial mesh
       chi(ndmx,nchix,npsx),     &! radial atomic orbitals
       oc(nchix,npsx),           &! atomic level occupation
       rho_at(ndmx,npsx),        &! radial atomic charge density
       rho_atc(ndmx,npsx)        ! radial core charge density
  INTEGER :: &
       mesh(npsx),               &! number of mesh points
       msh(npsx),                &! the point at rcut
       nchi(npsx),               &! number of atomic orbitals
       lchi(nchix,npsx)           ! angular momentum of atomic orbitals
  LOGICAL :: &
       numeric(npsx),            &! if .TRUE. the potential is in numeric form
       nlcc(npsx)                 ! if .TRUE. the atom has nlcc
  !
END MODULE atom
!
!
MODULE pseud
  !
  ! ... The variables describing pseudopotentials in analytical form
  !  
  USE kinds,      ONLY : DP
  USE parameters, ONLY : npsx
  !
  SAVE
  !
  REAL(KIND=DP) :: &
       cc(2,npsx),            &! the coefficients of the erf functions
       alpc(2,npsx),          &! the alpha of the erf functions
       zp(npsx),              &! the charge of the pseudopotential
       aps(6,0:3,npsx),       &! the a_l coefficient
       alps(3,0:3,npsx)        ! the b_l coefficient
  REAL(KIND=DP) :: &
       a_nlcc(npsx),         &! nonlinear core correction coefficients:
       b_nlcc(npsx),         &! rho_c(r) = (a_c + b_c*r^2) exp(-alpha_c*r^2)
       alpha_nlcc(npsx)       ! 
  INTEGER :: &
       nlc(npsx),             &! number of erf functions
       nnl(npsx),             &! number of the gaussian functions
       lmax(npsx),            &! maximum angular momentum of the pseudopot
       lloc(npsx)              ! angular momentum of the part taken as local
  !
END MODULE pseud
!
!
MODULE vlocal
  !
  ! ... The variables needed for the local potential in reciprocal space
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       strf(:,:)              ! the structure factor
  REAL(KIND=DP), ALLOCATABLE :: &
       vloc(:,:),            &! the local potential for each atom type
       vnew(:,:)              ! V_out - V_in, needed in scf force correction
  !
END MODULE vlocal
!
!
MODULE wvfct
  !
  ! ... The variables needed to compute the band structure
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  INTEGER ::  &
       npwx,             &! maximum number of PW for wavefunctions
       nbndx,            &! max number of bands use in iterative diag
       nbnd,             &! number of bands
       npw                ! the number of plane waves
  INTEGER, ALLOCATABLE, TARGET :: &
       igk(:),           &! correspondence k+G <-> G
       igk_l2g(:,:)       ! correspondence local index k+G <-> global G index
                          ! see also ig_l2g
  REAL(KIND=DP), ALLOCATABLE :: &
       et(:,:),          &! eigenvalues of the hamiltonian
       wg(:,:),          &! the weight of each k point and band
       g2kin(:)           ! kinetic energy
  LOGICAL :: &
       gamma_only         ! if .TRUE. only half G vectors are used
  !
  INTEGER :: ibm_baco2
  !
END MODULE wvfct
!
!
MODULE ener
  !
  ! ... The variables needed to compute the energies
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(KIND=DP) :: &
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
  !
END MODULE ener
!
!
MODULE force_mod
  !
  ! ... The variables for the first derivative of the energy
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(KIND=DP), ALLOCATABLE :: &
       force(:,:)       ! the force on each atom
  REAL(KIND=DP) :: &
       sigma(3,3)       ! the stress acting on the system
  LOGICAL :: &
       lforce,         &! if .TRUE. compute the forces
       lstres           ! if .TRUE. compute the stress
  !
END MODULE force_mod
!
!
MODULE scf
  !
  ! ... The variables needed to define the self-consistent cycle
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(KIND=DP), ALLOCATABLE :: &
       rho(:,:),       &! the charge density in real space
       rho_save(:,:),  &! another charge density in real space
       vr(:,:),        &! the Hartree + xc potential in real space
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grig)
       rho_core(:)      ! the core charge in real space
  !
END MODULE scf
!
!
MODULE relax
  !
  ! ... The variables used to control ionic relaxations
  !  
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : &
       if_pos,                  &! if 0 that coordinate will be kept fixed
       fixatom                   ! last "fixatom" are kept fixed
  !
  SAVE
  !
  LOGICAL :: &                   ! if .TRUE. start the structural optimization
       restart_bfgs              ! from the results of a previous run
  REAL(KIND=DP) :: &
       epse,                    &! threshold on total energy
       epsf,                    &! threshold on forces
       dtau_ref,                &! estimation of dtau
       starting_diag_threshold, &! self-explanatory
       starting_scf_threshold    ! as above
  !
END MODULE relax
!
!
MODULE cellmd
  !
  ! ... The variables used to control cell relaxation
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(KIND=DP) :: &
       press, cmass,     &! target pressure and cell mass,
       ttol,             &! tollerance for temperature rescaling
       at_old(3,3),      &! the lattice vectors at the previous ste
       omega_old,        &! the cell volume at the previous step
       cell_factor        ! maximum expected (linear) cell contraction
                          ! during relaxation/MD
  INTEGER :: &
       nzero,            &! iteration # of last thermalization
       ntimes,           &! number of thermalization steps to be performed
       ntcheck            ! # of steps between thermalizations
  LOGICAL :: lmovecell    ! used in cell relaxation
  !
  CHARACTER(LEN=2) :: &
       calc               ! main switch for variable cell shape MD
                          ! see readin, vcsmd and/or INPUT files
  !
END MODULE cellmd
!
!
!
MODULE char
  !
  ! ... The names of the atoms, of the solid and of the symmetries
  !  
  USE parameters, ONLY : npsx
  !
  SAVE
  !
  CHARACTER(LEN=75) ::  title       ! title of the run
  CHARACTER(LEN=20) ::  crystal     ! type of the solid
  CHARACTER(LEN=2 ) ::  psd(npsx)   ! name of the pseudopotential
  CHARACTER(LEN=45) ::  sname(48)   ! name of the symmetries
  !
END MODULE char
!
!
!
!
MODULE us
  !
  ! ... These parameters are needed with the US pseudopotentials
  !  
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lqmax, nbrx, npsx, nqfx, ndmx
  !
  SAVE
  !
  REAL(KIND=DP), PARAMETER:: &
       dq = 0.01D0           ! space between points in the pseudopotential tab.
  REAL(KIND=DP) :: &
       dion(nbrx,nbrx,npsx),              &! D_{mu,nu} parameters (in the 
                                           !  atomic case)
       betar(ndmx,nbrx,npsx),            &! radial beta_{mu} functions
       qqq(nbrx,nbrx,npsx),               &! q_{mu,nu} parameters (in the 
                                           !  atomic case)
       qfunc(ndmx,nbrx,nbrx,npsx),       &! Q_{mu,nu}(|r|) function for 
                                           !  |r|> r_L
       qfcoef(nqfx,lqmax,nbrx,nbrx,npsx), &! coefficients for Q in region 
                                           !  |r|<r_L
       vloc_at(ndmx,npsx),                &! local potential
       rinner(lqmax,npsx)                  ! values of r_L
  INTEGER :: &
       nh(npsx),             &! number of beta functions per atomic type
       nbeta(npsx),          &! number of beta functions
       kkbeta(npsx),         &! point where the beta are zero
       nqf(npsx),            &! number of coefficients for Q
       nqlc(npsx),           &! number of angular momenta in Q
       ifqopt(npsx),         &! level of q optimization
       lll(nbrx,npsx),       &! angular momentum of the beta function
       iver(3,npsx)           ! version of the atomic code
  INTEGER :: &
       nhm,              &! max number of different beta functions per atom
       nkb,              &! total number of beta functions, with struct.fact.
       nqxq,             &! size of interpolation table
       lmaxkb,           &! max angular momentum
       lqx,              &! max angular momentum + 1 for Q functions
       nqx                ! number of interpolation points
  INTEGER, ALLOCATABLE ::&
       indv(:,:),        &! correspondence of betas atomic <-> soli
       nhtol(:,:),       &! correspondence n <-> angular momentum
       nhtom(:,:)         ! correspondence n <-> magnetic angular m
  COMPLEX(KIND=DP), ALLOCATABLE, TARGET :: &
       vkb(:,:),              &! all beta functions in reciprocal space
       qgm(:)                  ! complete fourier transform of Q
  REAL(KIND=DP), ALLOCATABLE :: &
       qq(:,:,:),             &! the q functions in the solid
       dvan(:,:,:),           &! the D functions of the solid
       deeq(:,:,:,:),         &! the integral of V_eff and Q_{nm}
       becsum(:,:,:),         &! the sum of bec functions
       qrad(:,:,:,:),         &! radial FT of Q functions
       tab(:,:,:),            &! interpolation table for PPs
       tab_at(:,:,:)           ! interpolation table for atomic wfc
  LOGICAL :: &
       tvanp(npsx),           &! if .TRUE. the atom is of Vanderbilt type
       newpseudo(npsx),       &! if .TRUE. RRKJ3 US pseudopotentials
       okvan                   ! if .TRUE. at least one pseudo is Vanderbilt
  !
END MODULE us
!
!
MODULE ldaU
  !
  ! ... The quantities needed in lda+U calculations
  !  
  USE kinds,      ONLY : DP
  USE parameters, ONLY : lqmax, nspinx, ntypx
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       swfcatom(:,:)          ! orthogonalized atomic wfcs
  REAL(KIND=DP), ALLOCATABLE :: &
       ns(:,:,:,:),          &! the occupation matrix used in h_psi
       nsnew(:,:,:,:)         ! the occupation matrix computed by at
  REAL(KIND=DP) :: &       
       d1(3,3,48),           &! matrices for rotating spherical     
       d2(5,5,48),           &! harmonics                           
       d3(7,7,48),           &! 
       eth,                  &! the (corrected) Hubbard contribution
       Hubbard_U(ntypx),     &! the Hubbard U
       Hubbard_alpha(ntypx), &! the Hubbard alpha (used to calculate U)
       starting_ns(lqmax,nspinx,ntypx) !
  INTEGER :: &                                                
       niter_with_fixed_ns,  &! no. of iterations with fixed ns
       Hubbard_l(ntypx),     &! the agular momentum of Hubbard states
       Hubbard_lmax           ! maximum agular momentum of Hubbard states
  LOGICAL :: &                                                
       lda_plus_u,           &! .TRUE. if lda+u calculation is performed
       conv_ns                ! .TRUE. if ns are converged
  CHARACTER(LEN=30) :: &      ! 'atomic', 'ortho-atomic', 'file'
       U_projection           ! specifies how input coordinates are given
  !
END MODULE ldaU
!
!
MODULE extfield
  !
  ! ... The quantities needed in calculations with external field
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: &
       tefield,      &! if .TRUE. a finite electric field is added to the
                      ! local potential
       dipfield       ! if .TRUE. the dipole field is subtracted
  INTEGER :: &
       edir           ! direction of the field
  REAL(KIND=DP) :: &
      emaxpos,       &! position of the maximum of the field (0<emaxpos<1)
      eopreg,        &! amplitude of the inverse region (0<eopreg<1)
      eamp,          &! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
      etotefield      ! energy correction due to the field
  REAL(KIND=DP), ALLOCATABLE :: &
      forcefield(:,:)
  !
END MODULE extfield
!
!
MODULE sticks
  !
  ! ... data structure containing all information
  ! ... about fft data distribution for a given 
  ! ... potential grid, and its wave functions sub-grid.
  !  
  USE fft_types, ONLY : fft_dlay_descriptor
  !
  SAVE
  !
  TYPE ( fft_dlay_descriptor ) :: dfftp   ! dense grid
  TYPE ( fft_dlay_descriptor ) :: dffts   ! smooth grid
  !
END MODULE
!  
!  
MODULE bp
  !
  ! ... The variables needed for the Berry phase polarization calculation
  !
  SAVE
  !
  LOGICAL :: &
       lberry        ! if .TRUE., calculate polarization 
  INTEGER :: &
       gdir,        &! G-vector for polarization calculation
       nppstr        ! number of k-points (parallel vector)
  !
END MODULE bp
!
!
MODULE fixed_occ
  !
  ! ... The quantities needed in calculations with fixed occupations
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nbndxx, nspinx
  !
  SAVE
  !
  REAL(KIND=DP) :: &
       f_inp(nbndxx,nspinx)   ! the occupations for each spin
  LOGICAL :: &
       tfixed_occ             ! if .TRUE. the occupations are fixed.
  !
END MODULE fixed_occ
!
!
MODULE pwcom
  !
  USE constants, ONLY : e2, rytoev, amconv, uakbar, pi, tpi, fpi
  USE cell_base, ONLY : celldm, at, bg, alat, omega, tpiba, tpiba2, &
                        ibrav, symm_type
  USE basis
  USE dynam
  USE gvect
  USE gsmooth
  USE klist
  USE lsda_mod
  USE ktetra
  USE symme
  USE atom
  USE pseud
  USE vlocal
  USE wvfct
  USE ener
  USE force_mod
  USE scf
  USE relax
  USE cellmd
  USE char
  USE us
  USE ldaU
  USE extfield
  USE sticks
  USE bp
  USE fixed_occ
  !
END MODULE pwcom
