!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE pseudo_types

  !  this module contains the definitions of several TYPE structures,
  !  together with their allocation/deallocation routines

  USE upf_kinds,    ONLY: DP

  IMPLICIT NONE
  SAVE
  !
  ! Additional data to make a PAW setup out of an US pseudo,
  ! they are all stored on a radial grid:
  TYPE paw_in_upf
     REAL(DP),ALLOCATABLE :: ae_rho_atc(:) ! AE core charge (pseudo ccharge
     ! is already included in upf)
     REAL(DP),ALLOCATABLE :: pfunc(:,:,:),&! Psi_i(r)*Psi_j(r)
          pfunc_rel(:,:,:), & ! Psi_i(r)*Psi_j(r) small component
          ptfunc(:,:,:), & ! as above, but for pseudo
          aewfc_rel(:,:) ! as above, but for pseudo
     REAL(DP),ALLOCATABLE :: ae_vloc(:)    ! AE local potential (pseudo vloc
     ! is already included in upf)
     REAL(DP),ALLOCATABLE :: oc(:)         ! starting occupation used to init becsum
     ! they differ from US ones because they
     ! are indexed on BETA functions, non on WFC
     REAL(DP),ALLOCATABLE :: augmom(:,:,:) ! multipole AE-pseudo (i,j,l=0:2*lmax)
     REAL(DP)         :: raug          ! augfunction max radius
     INTEGER          :: iraug         ! index on rgrid closer to, and >, raug
     INTEGER          :: lmax_aug      ! max angmom of augmentation functions, it is ==
     ! to 2* max{l of pseudized wavefunctions}
     ! note that nqlc of upf also includes the angmom of
     ! empty virtual channel used to generate local potential
     REAL(DP)         :: core_energy   ! constant to add in order to get all-electron energy
     CHARACTER(len=12):: augshape      ! shape of augmentation charge
  END TYPE paw_in_upf


  TYPE pseudo_upf
     CHARACTER(LEN=80):: generated=' '! generator software
     CHARACTER(LEN=80):: author=' '   ! pseudopotential's author
     CHARACTER(LEN=80):: date=' '     ! generation date
     CHARACTER(LEN=80):: comment=' '  ! author's comment
     CHARACTER(LEN=2) :: psd=' '      ! Element label
     CHARACTER(LEN=4) :: typ=' '      ! Pseudo type:    NC, SL, US, PAW, 1/r
                                      ! NB: many files have USPP instead of US
     CHARACTER(len=6) :: rel=' '      ! relativistic: {no|scalar|full}
     LOGICAL :: tvanp              ! .true. if Ultrasoft
     LOGICAL :: tcoulombp          ! .true. if Coulomb 1/r potential
     LOGICAL :: nlcc               ! Non linear core corrections
     LOGICAL :: is_gth             ! .true. if Goedecker-Teter-Hutter
     LOGICAL :: is_multiproj       ! .true. if multiple projectors per l
     ! (for NC PP only; US-PP and PAW are assumed to be multi-projector)
     CHARACTER(LEN=25) :: dft      ! Exch-Corr type
     REAL(DP) :: zp                ! z valence
     REAL(DP) :: etotps            ! total energy
     REAL(DP) :: ecutwfc           ! suggested cut-off for wfc
     REAL(DP) :: ecutrho           ! suggested cut-off for rho
     !
     CHARACTER(len=11) :: nv       ! UPF file three-digit version i.e. 2.0.0
     INTEGER :: lmax               ! maximum l component in beta
     INTEGER :: lmax_rho           ! max l component in charge (should be 2*lmax)
     REAL(DP), ALLOCATABLE :: vnl(:,:,:) ! vnl(i,l,s) = V(r_i)_{ls}
     ! (semilocal form) only for single-channel NC PP
     ! Wavefunctions and projectors
     INTEGER :: nwfc               ! number of atomic wavefunctions
     INTEGER :: nbeta              ! number of projectors
     INTEGER,  ALLOCATABLE :: kbeta(:) ! kbeta(nbeta) see below
     INTEGER :: kkbeta             ! kkbeta=max(kbeta(:))
     !  kbeta<=mesh is the number of grid points for each beta function
     !              beta(r,nb) = 0 for r > r(kbeta(nb))
     ! kkbeta<=mesh is the largest of such number so that for all beta
     !              beta(r,nb) = 0 for r > r(kkbeta)
     !
     INTEGER,  ALLOCATABLE :: lll(:)     ! lll(nbeta) l of each projector
     REAL(DP), ALLOCATABLE :: beta(:,:)  ! beta(mesh,nbeta) projectors
     !
     CHARACTER(LEN=2), ALLOCATABLE :: els(:)  ! els(nwfc) label of wfc
     CHARACTER(LEN=2), ALLOCATABLE :: els_beta(:)  ! els(nbeta) label of beta
     INTEGER, ALLOCATABLE  :: nchi(:)    ! lchi(nwfc) value of pseudo-n for wavefcts
     INTEGER, ALLOCATABLE  :: lchi(:)    ! lchi(nwfc) value of l for wavefcts
     REAL(DP), ALLOCATABLE :: oc(:)      ! oc(nwfc) occupancies for wavefcts
     REAL(DP), ALLOCATABLE :: epseu(:)   ! pseudo one-particle energy (nwfc)
     REAL(DP), ALLOCATABLE :: rcut_chi(:)! rcut_chi(nwfc) cutoff inner radius
     REAL(DP), ALLOCATABLE :: rcutus_chi(:)! rcutus_chi(nwfc) ultrasoft outer radius
     ! Chi and rho_at are only used for initial density and initial wfcs:
     REAL(DP), ALLOCATABLE :: chi(:,:)   ! chi(mesh,nwfc) atomic wavefcts
     REAL(DP), ALLOCATABLE :: rho_at(:)  ! rho_at(mesh) atomic charge
     ! Minimal radial grid:
     INTEGER :: mesh               ! number of points in the radial mesh
     REAL(DP) :: xmin              ! the minimum x of the linear mesh
     REAL(DP) :: rmax              ! the maximum radius of the mesh
     REAL(DP) :: zmesh             ! the nuclear charge used for mesh
     REAL(DP) :: dx                ! the deltax of the linear mesh
     REAL(DP), ALLOCATABLE :: r(:)     ! r(mesh)  radial grid
     REAL(DP), ALLOCATABLE :: rab(:)   ! rab(mesh) dr(x)/dx (x=linear grid)
     ! Pseudized core charge
     REAL(DP), ALLOCATABLE :: rho_atc(:) ! rho_atc(mesh) atomic core charge
     ! Local potential
     INTEGER :: lloc                 ! L of channel used to generate local potential
     ! (if < 0 it was generated by smoothing AE potential)
     REAL(DP) :: rcloc               ! vloc = v_ae for r > rcloc
     REAL(DP), ALLOCATABLE :: vloc(:)    ! vloc(mesh) local atomic potential
     !
     REAL(DP), ALLOCATABLE :: dion(:,:)  ! dion(nbeta,nbeta) atomic D_{mu,nu}
     ! Augmentation
     LOGICAL :: q_with_l              ! if .true. qfunc is pseudized in
     ! different ways for different l
     INTEGER :: nqf                  ! number of Q coefficients
     INTEGER :: nqlc                 ! number of angular momenta in Q
     REAL(DP):: qqq_eps              ! qfunc is null if its norm is .lt. qqq_eps
     REAL(DP), ALLOCATABLE :: rinner(:)  ! rinner(0:2*lmax) r_L
     REAL(DP), ALLOCATABLE :: qqq(:,:)   ! qqq(nbeta,nbeta) q_{mu,nu}
     ! Augmentation without L dependecy
     REAL(DP), ALLOCATABLE :: qfunc(:,:) ! qfunc(mesh,nbeta*(nbeta+1)/2)
     ! Q_{mu,nu}(|r|) function for |r|> r_L
     ! Augmentation depending on L (optional, compulsory for PAW)
     REAL(DP), ALLOCATABLE :: qfuncl(:,:,:)!  qfuncl(mesh,nbeta*(nbeta+1)/2,l)
     ! Q_{mu,nu}(|r|) function for |r|> r_L
     ! Analitycal coeffs cor small r expansion of qfunc (Vanderbilt's code)
     REAL(DP), ALLOCATABLE :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
     ! coefficients for Q for |r|<r_L
     ! All electron and pseudo wavefunction, pswfc differ from chi as they are
     ! one for each beta, not just some choosen for initial conditions
     LOGICAL           :: has_wfc    ! if true, UPF contain AE and PS wfc for each beta
     REAL(DP), ALLOCATABLE :: aewfc(:,:) ! wfc(mesh,nbeta) all-electron wfc
     REAL(DP), ALLOCATABLE :: pswfc(:,:) ! wfc(mesh,nbeta) pseudo wfc

     LOGICAL :: has_so             ! if .true. includes spin-orbit
     INTEGER, ALLOCATABLE :: nn(:)     ! nn(nwfc) quantum number of wfc
     REAL(DP), ALLOCATABLE :: rcut(:)  ! cut-off radius(nbeta)
     REAL(DP), ALLOCATABLE :: rcutus(:)! ultrasoft cut-off radius (nbeta)
     REAL(DP), ALLOCATABLE :: jchi(:)  ! jchi(nwfc) j=l+1/2 or l-1/2 of wfc
     REAL(DP), ALLOCATABLE :: jjj(:)   ! jjj(nbeta) j=l+1/2 or l-1/2 of beta

     ! PAW:
     INTEGER :: paw_data_format      ! The version of the format
     LOGICAL  :: tpawp               ! true if atom is PAW, PAW data must be present
     TYPE(paw_in_upf) :: paw         ! additional data for PAW (see above)
     ! GIPAW:
     LOGICAL  :: has_gipaw           ! Whether GIPAW data is included
     LOGICAL  :: paw_as_gipaw        !EMINE
     INTEGER  :: gipaw_data_format   ! The version of the format
     INTEGER  :: gipaw_ncore_orbitals
     REAL(DP), ALLOCATABLE :: gipaw_core_orbital_n(:)
     REAL(DP), ALLOCATABLE :: gipaw_core_orbital_l(:)
     CHARACTER(LEN=2), ALLOCATABLE :: gipaw_core_orbital_el(:)
     REAL(DP), ALLOCATABLE :: gipaw_core_orbital(:,:)
     REAL(DP), ALLOCATABLE :: gipaw_vlocal_ae(:)
     REAL(DP), ALLOCATABLE :: gipaw_vlocal_ps(:)
     INTEGER :: gipaw_wfs_nchannels
     CHARACTER(LEN=2), ALLOCATABLE :: gipaw_wfs_el(:)
     INTEGER, ALLOCATABLE :: gipaw_wfs_ll(:)
     REAL(DP), ALLOCATABLE :: gipaw_wfs_ae(:,:)
     REAL(DP), ALLOCATABLE :: gipaw_wfs_rcut(:)
     REAL(DP), ALLOCATABLE :: gipaw_wfs_rcutus(:)
     REAL(DP), ALLOCATABLE :: gipaw_wfs_ps(:,:)
     !
     !  MD5 checksum ... used to verify integrity of the information contained
     !  in the pseudopotential file w.r.t previous run
     !
     CHARACTER(len=32) :: MD5_cksum = 'NOT SET'

  END TYPE pseudo_upf

  TYPE pseudo_config
     INTEGER :: nwfs
     CHARACTER(len=32)        :: pseud
     CHARACTER(len=2),ALLOCATABLE :: els(:)        ! label
     INTEGER,ALLOCATABLE          :: nns(:)        ! n
     INTEGER,ALLOCATABLE          :: lls(:)        ! l
     REAL(DP),ALLOCATABLE         :: ocs(:)        ! occupation
     REAL(DP),ALLOCATABLE         :: rcut(:)       ! NC cutoff radius
     REAL(DP),ALLOCATABLE         :: rcutus(:)     ! US cutoff radius
     REAL(DP),ALLOCATABLE         :: enls(:)       ! energy
  END TYPE pseudo_config

CONTAINS

  SUBROUTINE deallocate_paw_in_upf( paw )
    TYPE( paw_in_upf ), INTENT(INOUT) :: paw
    IF( ALLOCATED( paw%ae_rho_atc ) ) DEALLOCATE ( paw%ae_rho_atc )
    IF( ALLOCATED( paw%aewfc_rel ) )  DEALLOCATE (paw%aewfc_rel )
    IF( ALLOCATED( paw%pfunc ) )      DEALLOCATE ( paw%pfunc )
    IF( ALLOCATED( paw%pfunc_rel ) )  DEALLOCATE ( paw%pfunc_rel )
    IF( ALLOCATED( paw%ptfunc ) )     DEALLOCATE ( paw%ptfunc )
    IF( ALLOCATED( paw%ae_vloc )  )   DEALLOCATE ( paw%ae_vloc )
    IF( ALLOCATED( paw%augmom ) )     DEALLOCATE ( paw%augmom )
    IF( ALLOCATED( paw%oc ) )         DEALLOCATE ( paw%oc )
    paw%raug = 0.0
    paw%core_energy = 0.0
    paw%iraug = 0
    paw%lmax_aug = 0
    paw%augshape = ' '
  END SUBROUTINE deallocate_paw_in_upf
  !
   SUBROUTINE deallocate_pseudo_config(conf)
      TYPE(pseudo_config),INTENT(INOUT) :: conf
      IF ( ALLOCATED(conf%els)   ) DEALLOCATE(conf%els)
      IF ( ALLOCATED(conf%nns)   ) DEALLOCATE(conf%nns)
      IF ( ALLOCATED(conf%lls)   ) DEALLOCATE(conf%lls)
      IF ( ALLOCATED(conf%ocs)   ) DEALLOCATE(conf%ocs)
      IF ( ALLOCATED(conf%rcut)  ) DEALLOCATE(conf%rcut)
      IF ( ALLOCATED(conf%rcutus)) DEALLOCATE(conf%rcutus)
      IF ( ALLOCATED(conf%enls)  ) DEALLOCATE(conf%enls)
   END SUBROUTINE deallocate_pseudo_config



  SUBROUTINE deallocate_pseudo_upf( upf )
    TYPE( pseudo_upf ), INTENT(INOUT) :: upf
    CALL deallocate_paw_in_upf( upf%paw )
    IF( ALLOCATED( upf%els ) )     DEALLOCATE( upf%els )
    IF( ALLOCATED( upf%lchi ) )    DEALLOCATE( upf%lchi )
    IF( ALLOCATED( upf%nchi ) )    DEALLOCATE( upf%nchi )
    IF( ALLOCATED( upf%jchi ) )    DEALLOCATE( upf%jchi )
    IF( ALLOCATED( upf%oc ) )      DEALLOCATE( upf%oc )
    !
    IF( ALLOCATED( upf%r ) ) DEALLOCATE( upf%r )
    IF( ALLOCATED( upf%rab ) ) DEALLOCATE( upf%rab )
    !
    IF( ALLOCATED( upf%nn ) )      DEALLOCATE( upf%nn )
    IF( ALLOCATED( upf%els_beta ) )DEALLOCATE( upf%els_beta )
    IF( ALLOCATED( upf%rcut_chi ) )  DEALLOCATE( upf%rcut_chi )
    IF( ALLOCATED( upf%rcutus_chi ) )DEALLOCATE( upf%rcutus_chi )
    IF( ALLOCATED( upf%rcut ) )    DEALLOCATE( upf%rcut )
    IF( ALLOCATED( upf%rcutus ) )  DEALLOCATE( upf%rcutus )
    IF( ALLOCATED( upf%epseu ) )   DEALLOCATE( upf%epseu )
    IF( ALLOCATED( upf%rho_atc ) ) DEALLOCATE( upf%rho_atc )
    IF( ALLOCATED( upf%vloc ) )    DEALLOCATE( upf%vloc )
    IF( ALLOCATED( upf%lll ) )     DEALLOCATE( upf%lll )
    IF( ALLOCATED( upf%jjj ) )     DEALLOCATE( upf%jjj )
    IF( ALLOCATED( upf%kbeta ) )   DEALLOCATE( upf%kbeta )
    IF( ALLOCATED( upf%beta ) )    DEALLOCATE( upf%beta )
    IF( ALLOCATED( upf%vnl ) )     DEALLOCATE( upf%vnl )
    IF( ALLOCATED( upf%aewfc ) )   DEALLOCATE( upf%aewfc )
    IF( ALLOCATED( upf%pswfc ) )   DEALLOCATE( upf%pswfc )
    IF( ALLOCATED( upf%dion ) )    DEALLOCATE( upf%dion )
    IF( ALLOCATED( upf%rinner ) )  DEALLOCATE( upf%rinner )
    IF( ALLOCATED( upf%qqq ) )     DEALLOCATE( upf%qqq )
    IF( ALLOCATED( upf%qfunc ) )   DEALLOCATE( upf%qfunc )
    IF( ALLOCATED( upf%qfuncl ) )  DEALLOCATE( upf%qfuncl )
    IF( ALLOCATED( upf%qfcoef ) )  DEALLOCATE( upf%qfcoef )
    IF( ALLOCATED( upf%chi ) )     DEALLOCATE( upf%chi )
    IF( ALLOCATED( upf%rho_at ) )  DEALLOCATE( upf%rho_at )
    IF ( ALLOCATED ( upf%gipaw_core_orbital_n ) ) &
         DEALLOCATE ( upf%gipaw_core_orbital_n )
    IF ( ALLOCATED ( upf%gipaw_core_orbital_l ) ) &
         DEALLOCATE ( upf%gipaw_core_orbital_l )
    IF ( ALLOCATED ( upf%gipaw_core_orbital_el ) ) &
         DEALLOCATE ( upf%gipaw_core_orbital_el )
    IF ( ALLOCATED ( upf%gipaw_core_orbital ) ) &
         DEALLOCATE ( upf%gipaw_core_orbital )
    IF ( ALLOCATED ( upf%gipaw_vlocal_ae ) ) &
         DEALLOCATE ( upf%gipaw_vlocal_ae )
    IF ( ALLOCATED ( upf%gipaw_vlocal_ps ) ) &
         DEALLOCATE ( upf%gipaw_vlocal_ps )
    IF ( ALLOCATED ( upf%gipaw_wfs_el ) ) &
         DEALLOCATE ( upf%gipaw_wfs_el )
    IF ( ALLOCATED ( upf%gipaw_wfs_ll ) ) &
         DEALLOCATE ( upf%gipaw_wfs_ll )
    IF ( ALLOCATED ( upf%gipaw_wfs_ae ) ) &
         DEALLOCATE ( upf%gipaw_wfs_ae )
    IF ( ALLOCATED ( upf%gipaw_wfs_rcut ) ) &
         DEALLOCATE ( upf%gipaw_wfs_rcut )
    IF ( ALLOCATED ( upf%gipaw_wfs_rcutus ) ) &
         DEALLOCATE ( upf%gipaw_wfs_rcutus )
    IF ( ALLOCATED ( upf%gipaw_wfs_ps ) ) &
         DEALLOCATE ( upf%gipaw_wfs_ps )
  !
    upf%tvanp = .false.
    upf%tcoulombp = .false.
    upf%nlcc = .false.
    upf%dft = ' '
    upf%zp = 0.0
    upf%etotps = 0.0
    upf%ecutwfc = 0.0
    upf%ecutrho = 0.0
    upf%nv = ' '
    upf%lmax = 0
    upf%lmax_rho = 0
    upf%nwfc = 0
    upf%nbeta = 0
    upf%kkbeta = 0
    upf%mesh = 0
    upf%xmin = 0.0
    upf%rmax = 0.0
    upf%zmesh = 0.0
    upf%dx = 0.0
    upf%lloc = 0
    upf%rcloc = 0.0
    upf%q_with_l = .false.
    upf%nqf = 0
    upf%nqlc = 0
    upf%qqq_eps = 0.0
    upf%has_wfc = .false.
    upf%paw_data_format = 0
    upf%tpawp = .false.
    upf%has_gipaw = .false.
    upf%paw_as_gipaw = .false.
    upf%gipaw_data_format = 0
    upf%gipaw_ncore_orbitals = 0
    upf%gipaw_wfs_nchannels = 0

  END SUBROUTINE deallocate_pseudo_upf

END MODULE pseudo_types
