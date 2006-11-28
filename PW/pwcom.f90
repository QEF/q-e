!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
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
  REAL(DP), ALLOCATABLE, TARGET :: &
       g(:,:),        &! coordinates of G vectors
       gg(:)           ! modulus G^2 of G vectors
                       ! G vectors are in order of increasing |G|
  REAL(DP) :: &
       ecutwfc         ! energy cut-off
  REAL(DP), POINTER :: &
       gl(:)           ! the modulus of g in each shell
  REAL (DP) :: &
       gcutm,         &! cut-off for G vectors
       dual,          &! link between G of wavefunctions and charge
       ecfixed,       &!
       qcutz,         &! For the modified Ekin functional
       q2sigma         !
  complex(DP), ALLOCATABLE :: &
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
  REAL(DP) :: &
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
  REAL(DP) :: &
       xk(3,npk),         &! coordinates of k points
       wk(npk),           &! weight of k points
       xqq(3),            &! coordinates of q point (used with iswitch=-2)
       degauss,           &! smearing parameter
       nelec,             &! number of electrons
       nelup,             &! number of spin-up electrons (if two_fermi_energies=t)
       neldw,             &! number of spin-dw electrons (if two_fermi_energies=t)
       tot_charge
  INTEGER, ALLOCATABLE :: &
       ngk(:)              ! number of plane waves for each k point
  INTEGER :: &
       nks,               &! number of k points in this pool
       nkstot,            &! total number of k points
       ngauss,            &! type of smearing technique
       tot_magnetization, &! nelup-neldw >= 0 (negative value means unspecified)
       multiplicity        ! spin multiplicity
  LOGICAL :: &
       lgauss,         &! if .TRUE.: use gaussian broadening
       lxkcry,         &! if .TRUE.:k-pnts in cryst. basis accepted in input
       two_fermi_energies ! if .TRUE.: nelup and neldw set ef_up and ef_dw 
                          ! separately
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
  REAL(DP) :: &
       magtot,                       &! total magnetization
       absmag,                       &! total absolute magnetization  
       starting_magnetization(ntypx)  ! the magnetization used to start with
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
  INTEGER :: &
       t_rev(48) = 0            ! time reversal flag, for noncolinear magnetisation
  INTEGER, ALLOCATABLE :: &
       irt(:,:)                 ! symmetric atom for each atom and sym.op.
  LOGICAL :: &
       invsym                   ! if .TRUE. the system has inversion symmetry
  !
END MODULE symme

MODULE rap_point_group
   !
   USE kinds,      ONLY : DP
   !
   INTEGER :: &
          code_group,  &   ! The code of the point group
          nclass,  &       ! The number of classes of the point group
          nelem(12),   &   ! The elements of each class
          elem(8,12),  &   ! Which elements in the smat list for each class
          which_irr(12)    ! For each class gives its position in the 
                           ! character table.
   !
   COMPLEX(DP) :: char_mat(12,12)       ! the character tables

   CHARACTER(LEN=15) :: name_rap(12)  ! the name of the representation
   CHARACTER(LEN=11) :: gname           ! the name of the group
   CHARACTER(LEN=5) :: name_class(12)  ! the name of the class
   !
END MODULE rap_point_group

MODULE rap_point_group_so
   !
   USE kinds,      ONLY : DP
   !
   INTEGER :: &
          nrap,    &       ! The number of classes of the point group
          nelem_so(24),   &! The elements of each class
          elem_so(12,24),  &! Which elements in the smat list for each class
          has_e(12,24),  & ! if -1 the smat is multiplied by -E
          which_irr_so(24) ! For each class gives its position in the 
                           ! character table.
   !
   COMPLEX(DP) :: char_mat_so(12,24),  &   ! the character tables
                  d_spin(2,2,48)           ! the rotation in spin space

   CHARACTER(LEN=15) :: name_rap_so(12)  ! the name of the representation
   CHARACTER(LEN=5) :: name_class_so(24), &  ! the name of the class
                       name_class_so1(24)  ! the name of the class
   !
END MODULE rap_point_group_so
!
MODULE rap_point_group_is
   !
   USE kinds,      ONLY : DP
   !
   INTEGER :: &
          ftau_is(3,48), & ! The fractional transl. of the invariant subgroup
          nsym_is,       & ! The number of operations of the invariant subgroup
          code_group_is    ! The code of the point invariant subgroup

   REAL(DP) :: &
          sr_is(3,3,48)    ! The matrices of the invariant subgroup

   COMPLEX(DP) :: &
                d_spin_is(2,2,48)      ! the rotation in spin space

   CHARACTER(LEN=45) :: sname_is(48)   ! name of the symmetries
   CHARACTER(LEN=11) :: gname_is       ! the name of the invariant group
   !
END MODULE rap_point_group_is
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
  REAL(DP) :: &
       cc(2,npsx),            &! the coefficients of the erf functions
       alpc(2,npsx),          &! the alpha of the erf functions
       zp(npsx),              &! the charge of the pseudopotential
       aps(6,0:3,npsx),       &! the a_l coefficient
       alps(3,0:3,npsx)        ! the b_l coefficient
  REAL(DP) :: &
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
  COMPLEX(DP), ALLOCATABLE :: &
       strf(:,:)              ! the structure factor
  REAL(DP), ALLOCATABLE :: &
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
       npw,              &! the number of plane waves
       current_k          ! the index of k-point under consideration
  INTEGER, ALLOCATABLE, TARGET :: &
       igk(:)             ! index of G corresponding to a given index of k+G
  REAL(DP), ALLOCATABLE :: &
       et(:,:),          &! eigenvalues of the hamiltonian
       wg(:,:),          &! the weight of each k point and band
       g2kin(:)           ! kinetic energy
  INTEGER, ALLOCATABLE :: &
       btype(:,:)         ! one if the corresponding state has to be
                          ! converged to full accuracy, zero otherwise
  LOGICAL :: &
       gamma_only         ! if .TRUE. only half G vectors are used
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
  REAL(DP) :: &
       etot,           &! the total Kohn-Sham energy of the solid
       hwf_energy,     &! this is the Harris-Weinert-Foulkes energy
       eband,          &! the band energy
       deband,         &! scf correction to have variational energy
       ehart,          &! the hartree energy
       etxc,           &! the exchange and correlation energy
       vtxc,           &! another exchange-correlation energy
       etxcc,          &! the nlcc exchange and correlation
       ewld,           &! the ewald energy
       demet,          &! variational correction ("-TS") for metals
       ef, ef_up, ef_dw ! the fermi energy (up and dw if two_fermi_energies=.T.)
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
  REAL(DP), ALLOCATABLE :: &
       force(:,:)       ! the force on each atom
  REAL(DP) :: &
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
  REAL(DP), ALLOCATABLE :: &
       rho(:,:),       &! the charge density in real space
       vr(:,:),        &! the Hartree + xc potential in real space
       vltot(:),       &! the local potential in real space
       vrs(:,:),       &! the total pot. in real space (smooth grig)
       rho_core(:)      ! the core charge in real space
  COMPLEX(DP), ALLOCATABLE :: &
       rhog(:,:),      &! the charge density in reciprocal space
       rhog_core(:)     ! the core charge in reciprocal space
  !
END MODULE scf
!
!
MODULE relax
  !
  ! ... The variables used to control ionic relaxations
  !  
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: &                   ! if .TRUE. start the structural optimization
       restart_bfgs              ! from the results of a previous run
  REAL(DP) :: &
       epse,                    &! threshold on total energy
       epsf,                    &! threshold on forces
       epsp,                    &! threshold on pressure
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
  REAL(DP) :: &
       press, cmass,     &! target pressure and cell mass,
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
  ! ... The names of the system and of the symmetries
  !  
  USE printout_base, ONLY: title    ! title of the run
  !
  SAVE
  !
  CHARACTER(LEN=20) ::  crystal     ! type of the solid
  CHARACTER(LEN=45) ::  sname(48)   ! name of the symmetries
  !
END MODULE char
!
!
MODULE us
  !
  ! ... These parameters are needed with the US pseudopotentials
  !  
  USE kinds,      ONLY : DP
  !
  SAVE
  !
  INTEGER :: &
       nqxq,             &! size of interpolation table
       nqx                ! number of interpolation points
  REAL(DP), PARAMETER:: &
       dq = 0.01D0           ! space between points in the pseudopotential tab.
  REAL(DP), ALLOCATABLE :: &
       qrad(:,:,:,:),         &! radial FT of Q functions
       tab(:,:,:),            &! interpolation table for PPs
       tab_at(:,:,:)           ! interpolation table for atomic wfc
#ifdef USE_SPLINES
  REAL(DP), ALLOCATABLE :: &
       tab_d2y(:,:,:)            ! for cubic splines
#endif
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
  COMPLEX(DP), ALLOCATABLE :: &
       swfcatom(:,:)           ! orthogonalized atomic wfcs
  REAL(DP), ALLOCATABLE :: &
       ns(:,:,:,:),          &! the occupation matrix used in h_psi
       nsnew(:,:,:,:)         ! the occupation matrix computed by at
  REAL(DP) :: &       
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
  REAL(DP) :: &
      emaxpos,       &! position of the maximum of the field (0<emaxpos<1)
      eopreg,        &! amplitude of the inverse region (0<eopreg<1)
      eamp,          &! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
      etotefield      ! energy correction due to the field
  REAL(DP), ALLOCATABLE :: &
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
END MODULE sticks
!  
!  
!
MODULE fixed_occ
  !
  ! ... The quantities needed in calculations with fixed occupations
  !
  USE kinds,      ONLY : DP
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: &
       f_inp(:,:)             ! the occupations for each spin
  LOGICAL :: &
       tfixed_occ             ! if .TRUE. the occupations are fixed.
  !
END MODULE fixed_occ

MODULE spin_orb
  
  USE kinds, ONLY: DP
  USE parameters, ONLY : lmaxx, npsx
  
  SAVE

  LOGICAL :: &
      lspinorb,  &       ! if .TRUE. this is a spin-orbit calculation
      domag,     &       ! if .TRUE. magnetization is computed
      so(npsx)           ! for each pseudo tells if it is spin-orbit


  COMPLEX (DP) :: rot_ylm(2*lmaxx+1,2*lmaxx+1)  ! transform real
                         ! spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:) ! function needed to
                         ! account for spinors.
END MODULE spin_orb
!
!
MODULE bp
  USE kinds, ONLY: DP
  !
  ! ... The variables needed for the Berry phase polarization calculation
  !
  SAVE
  !
  LOGICAL :: &
       lberry, &        ! if .TRUE., calculate polarization
       lelfield         ! if .TRUE. static homogeneous electric field
  INTEGER :: &
       gdir,        &! G-vector for polarization calculation
       nppstr,      &! number of k-points (parallel vector)
       nberrycyc     !numer of cycles for cobergence in electric field without changing the selfconsistent charge
  REAL(DP) :: efield ! electric field intensity in a.u.
  COMPLEX(DP), ALLOCATABLE , TARGET :: evcel(:,:) ! wave function for calculating the electric field operator

!
END MODULE bp
!
MODULE pwcom
  !
  USE constants, ONLY : e2, rytoev, amconv, pi, tpi, fpi
  USE cell_base, ONLY : celldm, at, bg, alat, omega, tpiba, tpiba2, &
                        ibrav, symm_type
  ! TEMP
  USE uspp
  ! TEMP
  USE basis
  USE gvect
  USE gsmooth
  USE klist
  USE lsda_mod
  USE ktetra
  USE symme
  USE rap_point_group
  USE rap_point_group_so
  USE rap_point_group_is
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
  USE spin_orb
  USE bp
  !
END MODULE pwcom
