!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE rism
  !--------------------------------------------------------------------------
  !
  ! ... this module keeps data for 1D- or 3D-RISM calculations.
  ! ... also Laue-RISM is incluted.
  !
  USE fft_types,   ONLY : fft_type_descriptor, fft_type_deallocate
  USE gvec_3drism, ONLY : gvec_type, deallocate_gvec_3drism
  USE kinds,       ONLY : DP
  USE lauefft,     ONLY : lauefft_type, allocate_lauefft, deallocate_lauefft, &
                        & set_lauefft_offset, set_lauefft_offset0, set_lauefft_barrier
  USE mp_rism,     ONLY : mp_rism_site, mp_rism_task, mp_set_index_rism_site, &
                        & mp_set_index_rism_task, mp_start_rism_task_and_site, &
                        & mp_start_rism_task_on_site, mp_end_rism
  USE radfft,      ONLY : radfft_type, allocate_radfft, init_mpi_radfft, &
                        & deallocate_radfft
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  ! ..... type of `rism_type'
  INTEGER, PARAMETER :: ITYPE_NULL     = 0
  INTEGER, PARAMETER :: ITYPE_1DRISM   = 1
  INTEGER, PARAMETER :: ITYPE_3DRISM   = 2
  INTEGER, PARAMETER :: ITYPE_LAUERISM = 3
  !
  ! ..... type of Closure eqation
  INTEGER, PARAMETER :: CLOSURE_NULL = 0
  INTEGER, PARAMETER :: CLOSURE_HNC  = 1
  INTEGER, PARAMETER :: CLOSURE_KH   = 2
  !
  ! ..... type of chemical potential
  INTEGER, PARAMETER :: CHEMPOT_NULL = 0
  INTEGER, PARAMETER :: CHEMPOT_HNC  = 1
  INTEGER, PARAMETER :: CHEMPOT_KH   = 2
  INTEGER, PARAMETER :: CHEMPOT_GF   = 3
  !
  ! ... define data of 1D- or 3D-RISM calculations
  ! ..... units: length      -> [bohr]
  ! .....        energy      -> [Ry]
  ! .....        temperature -> [kelvin]
  ! .....
  ! ..... for Laue-RISM,
  ! .....        csgz, vlgz, hgz, hsgz, hlgz, rhog, vpot
  ! .....        are not in G-space but in Laue-representation.
  ! .....
  TYPE rism_type
    !
    ! ... control variables
    LOGICAL              :: avail   = .FALSE.       ! this rism_type is available ?
    INTEGER              :: itype   = ITYPE_NULL    ! data for 1D-RISM or 3D-RISM or Laue-RISM ?
    INTEGER              :: closure = CLOSURE_NULL  ! type of Closure equation
    REAL(DP)             :: temp    = 300.0_DP      ! temperature of solvent system
    REAL(DP)             :: tau     = 1.0_DP        ! coulomb smearing radius
    !
    ! ... dimensions of data
    INTEGER              :: nsite        ! number of solvent sites
    INTEGER              :: nr           ! number of meshes in R-space
    INTEGER              :: nrzs         ! number of meshes along short Z-stick in R-space (for Laue-RISM)
    INTEGER              :: nrzl         ! number of meshes along long Z-stick in R-space (for Laue-RISM)
    INTEGER              :: ng           ! number of meshes in G-space
    INTEGER              :: ngs          ! number of shells in G-space (in Gxy-space for Laue-RISM)
    INTEGER              :: ngxy         ! number of meshes on XY-plane in G-space (for Laue-RISM)
    !
    ! ... data to calculate RISM
    REAL(DP),    POINTER :: csr  (:,:)   ! short-range direct correlations in R-space
    REAL(DP),    POINTER :: csg  (:,:)   ! short-range direct correlations in G-space
    COMPLEX(DP), POINTER :: csgz (:,:)   ! short-range direct correlations in G-space or Laue-rep. (complex)
    REAL(DP),    POINTER :: csg0 (:,:)   ! short-range direct correlations when Gxy = 0 (Laue-RISM)
    REAL(DP),    POINTER :: csdr (:,:)   ! short-range + dipole of direct correlations in R-space (Laue-RISM)
    REAL(DP),    POINTER :: csdg0(:,:)   ! short-range + dipole of direct correlations when Gxy = 0 (Laue-RISM)
    REAL(DP),    POINTER :: cda  (:)     ! dipole-parts of direct correlations in Laue-rep. (amplitude)
    REAL(DP),    POINTER :: cdz  (:)     ! dipole-parts of direct correlations in Laue-rep. (step function)
    REAL(DP),    POINTER :: uljr (:,:)   ! Lennard-Jones potential functions in R-space
    REAL(DP),    POINTER :: uwr  (:,:)   ! repulsive-wall potential functions in R-space (Laue-RISM)
    REAL(DP),    POINTER :: usr  (:,:)   ! short-range potential functions in R-space
    REAL(DP),    POINTER :: usg0 (:,:)   ! short-range potential functions when Gxy = 0 (Laue-RISM)
    REAL(DP),    POINTER :: vsr  (:)     ! short-range coulomb potential in R-space
    REAL(DP),    POINTER :: ulr  (:,:)   ! long-range potential functions in R-space
    REAL(DP),    POINTER :: vlr  (:)     ! long-range coulomb potential in R-space
    REAL(DP),    POINTER :: ulg  (:,:)   ! long-range potential functions in G-space
    COMPLEX(DP), POINTER :: ulgz (:,:)   ! long-range potential functions in G-space (complex)
    COMPLEX(DP), POINTER :: vlgz (:)     ! long-range coulomb potential in G-space or Laue-rep. (complex)
    COMPLEX(DP), POINTER :: vright(:)    ! long-range coulomb potential coefficient for z > zright (Laue-RISM)
    COMPLEX(DP), POINTER :: vleft (:)    ! long-range coulomb potential coefficient for z < zleft  (Laue-RISM)
    LOGICAL,     POINTER :: do_vright(:) ! to consider vright or not (Laue-RISM)
    LOGICAL,     POINTER :: do_vleft (:) ! to consider vleft  or not (Laue-RISM)
    REAL(DP),    POINTER :: hr   (:,:)   ! total correlations in R-space
    REAL(DP),    POINTER :: hg   (:,:)   ! total correlations in G-space
    COMPLEX(DP), POINTER :: hgz  (:,:)   ! total correlations in G-space or Laue-rep. (complex)
    REAL(DP),    POINTER :: hg0  (:,:)   ! total correlations when Gxy = 0 (Laue-RISM)
    COMPLEX(DP), POINTER :: hsgz (:,:)   ! short-range total correlations in Laue-rep. (w/ dipole parts) (complex)
    COMPLEX(DP), POINTER :: hlgz (:,:)   ! long-range total correlations in Laue-rep. (complex)
    REAL(DP),    POINTER :: hdz  (:,:,:) ! dipole-parts of total correlations in Laue-rep. (integrated intermediate)
    REAL(DP),    POINTER :: gr   (:,:)   ! distribution functions in R-space
    REAL(DP),    POINTER :: gg0  (:,:)   ! distribution functions when Gxy = 0 (Laue-RISM)
    REAL(DP),    POINTER :: wg   (:,:)   ! intra-molecular correlations in G-space
    REAL(DP),    POINTER :: zg   (:,:)   ! intra-molecular dielectric bridge in G-space (for DRISM)
    REAL(DP),    POINTER :: xgs  (:,:,:) ! inter-site susceptibility in G-shell or Laue-rep.
    REAL(DP),    POINTER :: xgs0 (:,:,:) ! integrated xgs       (Laue-RISM).
    REAL(DP),    POINTER :: xgs1 (:,:,:) ! integrated (z * xgs) (Laue-RISM).
    REAL(DP),    POINTER :: ygs  (:,:,:) ! ygs is left-hand version of xgs (Laue-RISM).
    REAL(DP),    POINTER :: ygs0 (:,:,:) ! ygs0 is left-hand version of xgs0 (Laue-RISM).
    REAL(DP),    POINTER :: ygs1 (:,:,:) ! ygs1 is left-hand version of xgs1 (Laue-RISM).
    !
    ! ... results from RISM
    REAL(DP),    POINTER :: nsol(:)      ! solvent's number for each site
    REAL(DP),    POINTER :: qsol(:)      ! solvent's charge for each site
    REAL(DP)             :: qtot         ! solvent's total charge
    REAL(DP),    POINTER :: usol(:)      ! solvation's chemical potential for each site
    REAL(DP),    POINTER :: usol_GF(:)   ! solvation's chemical potential for each site (by G.F.)
    REAL(DP)             :: esol         ! solvation's energy
    REAL(DP)             :: vsol         ! potential shifting (for Laue-RISM)
    COMPLEX(DP), POINTER :: rhog(:)      ! charge density of solvent in G-space or Laue-rep.
    COMPLEX(DP), POINTER :: vpot(:)      ! coulomb potential of solvent in G-space or Laue-rep.
    INTEGER              :: pbc_nfit     ! number of fitting points for rhog_pbc and vpot_pbc. (for Laue-RISM)
    COMPLEX(DP), POINTER :: rhog_pbc(:)  ! effective charge density of solvent in G-space (for Laue-RISM)
    COMPLEX(DP), POINTER :: vpot_pbc(:)  ! effective coulomb potential of solvent in G-space (for Laue-RISM)
    !
    ! ... for MPI
    INTEGER              :: super_comm   ! parent group communicator
    INTEGER              :: super_root   ! root rank of parent group
    LOGICAL              :: is_intra     ! this process is in intra group or not
    INTEGER              :: intra_comm   ! intra group communicator
    TYPE(mp_rism_site)   :: mp_site      ! MPI-data for site parallel
    TYPE(mp_rism_task)   :: mp_task      ! MPI-data for task parallel
    !
    ! ... for FFT
    TYPE(radfft_type)         :: rfft    ! radial FFT for 1D-RISM
    TYPE(fft_type_descriptor) :: dfft    ! 3D-FFT(custom) for 3D-RISM
    TYPE(gvec_type)           :: gvec    ! G-vectors(custom) for 3D-RISM
    TYPE(lauefft_type)        :: lfft    ! 1D- and 2D-FFT for Laue-RISM
    !
  END TYPE rism_type
  !
  ! ... public components
  PUBLIC :: ITYPE_NULL
  PUBLIC :: ITYPE_1DRISM
  PUBLIC :: ITYPE_3DRISM
  PUBLIC :: ITYPE_LAUERISM
  PUBLIC :: CLOSURE_NULL
  PUBLIC :: CLOSURE_HNC
  PUBLIC :: CLOSURE_KH
  PUBLIC :: CHEMPOT_NULL
  PUBLIC :: CHEMPOT_HNC
  PUBLIC :: CHEMPOT_KH
  PUBLIC :: CHEMPOT_GF
  PUBLIC :: rism_type
  PUBLIC :: allocate_1drism
  PUBLIC :: allocate_3drism
  PUBLIC :: allocate_lauerism
  PUBLIC :: refresh_suscept_3drism
  PUBLIC :: refresh_suscept_lauerism
  PUBLIC :: deallocate_rism
  PUBLIC :: clean_rism_data
  PUBLIC :: get_chempot_type
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_1drism(rismt, nv, ngrid, rmax, mpi_radfft, &
                           & super_comm, super_root, is_intra, intra_comm)
    !--------------------------------------------------------------------------
    !
    ! ... initialize rism_type for 1D-RISM
    ! ...
    ! ... Variables:
    ! ...   nv:         number of solvent's sites
    ! ...   ngrid:      number of grids
    ! ...   rmax:       maximum radius of R-space (in bohr)
    ! ...   mpi_radfft: use MPI for radial FFT
    ! ...   super_comm: MPI's parent communicator
    ! ...   super_root: root rank of MPI's parent
    ! ...   is_intra:   this process is in MPI's intra group
    ! ...   intra_comm: MPI's intra communicator
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    INTEGER,         INTENT(IN)    :: nv
    INTEGER,         INTENT(IN)    :: ngrid
    REAL(DP),        INTENT(IN)    :: rmax
    LOGICAL,         INTENT(IN)    :: mpi_radfft
    INTEGER,         INTENT(IN)    :: super_comm
    INTEGER,         INTENT(IN)    :: super_root
    LOGICAL,         INTENT(IN)    :: is_intra
    INTEGER,         INTENT(IN)    :: intra_comm
    !
    INTEGER :: nvv
    INTEGER :: nsite
    INTEGER :: mgrid
    !
    nvv   = nv * (nv + 1) / 2
    nsite = nvv
    !
    ! ... check parameters
    IF (nsite <= 0) THEN
      CALL errore(' allocate_1drism ', ' too small nsite ', 1)
    END IF
    !
    IF (ngrid <= 0) THEN
      CALL errore(' allocate_1drism ', ' too small ngrid ', 1)
    END IF
    !
    IF (rmax <= 0.0_DP) THEN
      CALL errore(' allocate_1drism ', ' too small rmax ', 1)
    END IF
    !
    ! ... initialize MPI
    rismt%super_comm = super_comm
    rismt%super_root = super_root
    rismt%is_intra   = is_intra
    rismt%intra_comm = intra_comm
    CALL mp_start_rism_task_and_site(rismt%mp_site, rismt%mp_task, intra_comm)
    CALL mp_set_index_rism_site(rismt%mp_site, nsite)
    CALL mp_set_index_rism_task(rismt%mp_task, ngrid)
    !
    ! ... initialize radial FFT
    CALL allocate_radfft(rismt%rfft, ngrid, rmax)
    !
    IF (mpi_radfft) THEN
      CALL init_mpi_radfft(rismt%rfft, rismt%mp_task%itask_comm, &
                         & rismt%mp_task%ivec_start, rismt%mp_task%ivec_end)
    END IF
    !
    ! ... initialize rismt
    mgrid = rismt%mp_task%ivec_end - rismt%mp_task%ivec_start + 1
    CALL allocate_rism(rismt, ITYPE_1DRISM, nsite, nsite, mgrid, 0, 0, mgrid, mgrid, 0, .FALSE.)
    !
  END SUBROUTINE allocate_1drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_3drism(rismt, nv, ecutv, itask_comm, intra_comm)
    !--------------------------------------------------------------------------
    !
    ! ... initialize rism_type for 3D-RISM
    ! ...
    ! ... Variables:
    ! ...   nv:         number of solvent's sites
    ! ...   ecutv:      energy cutoff on G-space (in Ry)
    ! ...   itask_comm: MPI's communicator for task parallel
    ! ...   intra_comm: MPI's communicator, which is global in 3D-RISM calculation
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    INTEGER,         INTENT(IN)    :: nv
    REAL(DP),        INTENT(IN)    :: ecutv
    INTEGER,         INTENT(IN)    :: itask_comm
    INTEGER,         INTENT(IN)    :: intra_comm
    !
    INTEGER :: nsite
    INTEGER :: msite
    INTEGER :: nr
    INTEGER :: ng
    INTEGER :: ngs
    !
    nsite = nv
    !
    ! ... check parameters (solvent)
    IF (nsite <= 0) THEN
      CALL errore(' allocate_3drism ', ' too small nsite ', 1)
    END IF
    !
    ! ... initialize MPI
    rismt%super_comm = intra_comm
    rismt%super_root = 0
    rismt%is_intra   = .TRUE.
    rismt%intra_comm = intra_comm
    CALL mp_start_rism_task_on_site(rismt%mp_site, rismt%mp_task, itask_comm, intra_comm)
    CALL mp_set_index_rism_site(rismt%mp_site, nsite)
    !
    ! ... initialize 3D-FFT(custom)
    CALL allocate_fft_3drism(rismt%dfft, rismt%gvec, ecutv, .FALSE., rismt%mp_task)
    !
    nr  = rismt%dfft%nnr
    ng  = rismt%gvec%ngm
    ngs = rismt%gvec%ngl
    !
    ! ... check parameters (3D-FFT)
    IF (nr <= 0) THEN
      CALL errore(' allocate_3drism ', ' too small nr ', 1)
    END IF
    !
    IF (ng <= 0) THEN
      CALL errore(' allocate_3drism ', ' too small ng ', 1)
    END IF
    !
    IF (ngs <= 0) THEN
      CALL errore(' allocate_3drism ', ' too small ngs ', 1)
    END IF
    !
    ! ... initialize rismt
    msite = rismt%mp_site%isite_end - rismt%mp_site%isite_start + 1
    CALL allocate_rism(rismt, ITYPE_3DRISM, msite, nsite, nr, 0, 0, ng, ngs, 0, .FALSE.)
    !
  END SUBROUTINE allocate_3drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_lauerism(rismt, nv, ecutv, nfit, &
                             & dzright, dzleft, zright, zleft, &
                             & zrightu, zrightv, zleftu, zleftv, zrightb, zleftb, &
                             & lboth, itask_comm, intra_comm)
    !--------------------------------------------------------------------------
    !
    ! ... initialize rism_type for Laue-RISM
    ! ...
    ! ... Variables:
    ! ...   nv:         number of solvent's sites
    ! ...   ecutv:      energy cutoff on G-space (in Ry)
    ! ...   nfit:       #fitting points for effective solvation potential
    ! ...   dzright:    |zright - z0| (in alat units)
    ! ...   dzleft:     |zleft  + z0| (in alat units)
    ! ...   zright:     offset of solvent region on right-hand side (in alat units)
    ! ...   zleft:      offset of solvent region on left-hand side (in alat units)
    ! ...   zrightu:    offset of solvent region on right-hand side for Gxy=0 (in alat units)
    ! ...   zrightv:    -> solute-side and solvent-side, respectively.
    ! ...   zleftu:     offset of solvent region on left-hand side for Gxy=0 (in alat units)
    ! ...   zleftv:     -> solute-side and solvent-side, respectively.
    ! ...   zrightb:    offset of solvent barrier on right-hand side (in alat units)
    ! ...   zleftb:     offset of solvent barrier on left-hand side (in alat units)
    ! ...   lboth:      apply both-hands calculation, or not
    ! ...   itask_comm: MPI's communicator for task parallel
    ! ...   intra_comm: MPI's communicator, which is global in Laue-RISM calculation
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    INTEGER,         INTENT(IN)    :: nv
    REAL(DP),        INTENT(IN)    :: ecutv
    INTEGER,         INTENT(IN)    :: nfit
    REAL(DP),        INTENT(IN)    :: dzright
    REAL(DP),        INTENT(IN)    :: dzleft
    REAL(DP),        INTENT(IN)    :: zright
    REAL(DP),        INTENT(IN)    :: zleft
    REAL(DP),        INTENT(IN)    :: zrightu
    REAL(DP),        INTENT(IN)    :: zrightv
    REAL(DP),        INTENT(IN)    :: zleftu
    REAL(DP),        INTENT(IN)    :: zleftv
    REAL(DP),        INTENT(IN)    :: zrightb
    REAL(DP),        INTENT(IN)    :: zleftb
    LOGICAL,         INTENT(IN)    :: lboth
    INTEGER,         INTENT(IN)    :: itask_comm
    INTEGER,         INTENT(IN)    :: intra_comm
    !
    INTEGER :: nsite
    INTEGER :: msite
    INTEGER :: nr
    INTEGER :: nrzs
    INTEGER :: nrzl
    INTEGER :: ng
    INTEGER :: ngs
    INTEGER :: ngxy
    !
    nsite = nv
    !
    ! ... check parameters (solvent)
    IF (nsite <= 0) THEN
      CALL errore(' allocate_lauedrism ', ' too small nsite ', 1)
    END IF
    !
    ! ... initialize MPI
    rismt%super_comm = intra_comm
    rismt%super_root = 0
    rismt%is_intra   = .TRUE.
    rismt%intra_comm = intra_comm
    CALL mp_start_rism_task_on_site(rismt%mp_site, rismt%mp_task, itask_comm, intra_comm)
    CALL mp_set_index_rism_site(rismt%mp_site, nsite)
    !
    ! ... initialize 3D-FFT(custom)
    CALL allocate_fft_3drism(rismt%dfft, rismt%gvec, ecutv, .TRUE., rismt%mp_task)
    !
    ! ... initialize Laue-FFT
    CALL allocate_lauefft(rismt%lfft, rismt%dfft, dzright, dzleft, &
    & rismt%gvec%ngm, rismt%gvec%ig1, rismt%gvec%ig2, rismt%gvec%ig3, &
    & rismt%gvec%g, rismt%gvec%gcutm, rismt%mp_task%itask_comm)
    !
    CALL set_lauefft_offset( rismt%lfft, zright,  zleft)
    CALL set_lauefft_offset0(rismt%lfft, zrightu, zrightv, zleftu, zleftv)
    CALL set_lauefft_barrier(rismt%lfft, zrightb, zleftb)
    !
    ! ... set #fitting points
    rismt%pbc_nfit = nfit
    !
    nr   = rismt%dfft%nnr
    nrzs = rismt%dfft%nr3
    nrzl = rismt%lfft%nrz
    ng   = rismt%gvec%ngm
    ngs  = rismt%lfft%nglxy
    ngxy = rismt%lfft%ngxy
    !
    ! ... check parameters (3D-FFT and Laue-FFT)
    IF (nr <= 0) THEN
      CALL errore(' allocate_lauerism ', ' too small nr ', 1)
    END IF
    !
    IF (nrzs <= 0) THEN
      CALL errore(' allocate_lauerism ', ' too small nrzs ', 1)
    END IF
    !
    IF (nrzl <= 0) THEN
      CALL errore(' allocate_lauerism ', ' too small nrzl ', 1)
    END IF
    !
    IF (ng <= 0) THEN
      CALL errore(' allocate_lauerism ', ' too small ng ', 1)
    END IF
    !
    IF (ngs <= 0) THEN
      CALL errore(' allocate_lauerism ', ' too small ngs ', 1)
    END IF
    !
    IF (ngxy <= 0) THEN
      CALL errore(' allocate_lauerism ', ' too small ngxy ', 1)
    END IF
    !
    IF (rismt%pbc_nfit < 0) THEN
      CALL errore(' allocate_lauerism ', ' negative pbc_nfit ', 1)
    END IF
    !
    ! ... initialize rismt
    msite = rismt%mp_site%isite_end - rismt%mp_site%isite_start + 1
    CALL allocate_rism(rismt, ITYPE_LAUERISM, msite, nsite, nr, nrzs, nrzl, ng, ngs, ngxy, lboth)
    !
  END SUBROUTINE allocate_lauerism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_rism(rismt, itype, nsite, nsite_t, nr, nrzs, nrzl, ng, ngs, ngxy, lboth)
    !--------------------------------------------------------------------------
    !
    ! ... initialize rism_type
    ! ...
    ! ... Variables:
    ! ...   itype:   data type (ITYPE_1DRISM, ITYPE_3DRISM or ITYPE_LAUERISM)
    ! ...   nsite:   number of solvent's sites (i.e. degrees of freedom for correlation functions)
    ! ...   nsite_t: sum of nsite over all processies
    ! ...   nr:      dimension of R-space grid
    ! ...   nrzs:    dimension of R-space grid along short Z-stick
    ! ...   nrzl:    dimension of R-space grid along long Z-stick
    ! ...   ng:      dimension of G-space grid
    ! ...   ngs:     number of shells in G-space (in Gxy-space for Laue-RISM)
    ! ...   ngxy:    dimension of G-space grid on XY-plane
    ! ...   lboth:   apply both-hands calculation, or not (for Laue-RISM)
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    INTEGER,         INTENT(IN)    :: itype
    INTEGER,         INTENT(IN)    :: nsite
    INTEGER,         INTENT(IN)    :: nsite_t
    INTEGER,         INTENT(IN)    :: nr
    INTEGER,         INTENT(IN)    :: nrzs
    INTEGER,         INTENT(IN)    :: nrzl
    INTEGER,         INTENT(IN)    :: ng
    INTEGER,         INTENT(IN)    :: ngs
    INTEGER,         INTENT(IN)    :: ngxy
    LOGICAL,         INTENT(IN)    :: lboth
    !
    ! ... set variables
    rismt%avail = .FALSE.
    rismt%itype = itype
    rismt%nsite = nsite
    rismt%nr    = nr
    rismt%nrzs  = nrzs
    rismt%nrzl  = nrzl
    rismt%ng    = ng
    rismt%ngs   = ngs
    rismt%ngxy  = ngxy
    rismt%esol  = 0.0_DP
    rismt%vsol  = 0.0_DP
    !
    ! ... allocate arrays
    ! ..... R-space
    IF (itype == ITYPE_1DRISM) THEN
      IF ((nr * nsite) > 0) THEN
        ALLOCATE(rismt%csr(nr, nsite))
        ALLOCATE(rismt%usr(nr, nsite))
        ALLOCATE(rismt%ulr(nr, nsite))
        ALLOCATE(rismt%hr( nr, nsite))
        ALLOCATE(rismt%gr( nr, nsite))
      END IF
      !
    ELSE IF (itype == ITYPE_3DRISM) THEN
      IF ((nr * nsite) > 0) THEN
        ALLOCATE(rismt%csr( nr, nsite))
        ALLOCATE(rismt%usr( nr, nsite))
        ALLOCATE(rismt%ulr( nr, nsite))
        ALLOCATE(rismt%hr(  nr, nsite))
        ALLOCATE(rismt%gr(  nr, nsite))
        ALLOCATE(rismt%uljr(nr, nsite))
      END IF
      IF (nr > 0) THEN
        ALLOCATE(rismt%vsr(nr))
        ALLOCATE(rismt%vlr(nr))
      END IF
      !
    ELSE IF (itype == ITYPE_LAUERISM) THEN
      IF ((nr * nsite) > 0) THEN
        ALLOCATE(rismt%csr( nr, nsite))
        ALLOCATE(rismt%csdr(nr, nsite))
        ALLOCATE(rismt%usr( nr, nsite))
        ALLOCATE(rismt%hr(  nr, nsite))
        ALLOCATE(rismt%gr(  nr, nsite))
        ALLOCATE(rismt%uljr(nr, nsite))
        ALLOCATE(rismt%uwr( nr, nsite))
      END IF
      IF (nr > 0) THEN
        ALLOCATE(rismt%vsr(nr))
        ALLOCATE(rismt%vlr(nr))
      END IF
    END IF
    !
    ! ..... G-space (or Laue-rep.)
    IF (itype == ITYPE_1DRISM) THEN
      IF ((ng * nsite) > 0) THEN
        ALLOCATE(rismt%csg(ng, nsite))
        ALLOCATE(rismt%ulg(ng, nsite))
        ALLOCATE(rismt%hg( ng, nsite))
        ALLOCATE(rismt%wg( ng, nsite))
        ALLOCATE(rismt%zg( ng, nsite))
      END IF
      !
    ELSE IF (itype == ITYPE_3DRISM) THEN
      IF ((ng * nsite) > 0) THEN
        ALLOCATE(rismt%csgz(ng, nsite))
        ALLOCATE(rismt%ulgz(ng, nsite))
        ALLOCATE(rismt%hgz( ng, nsite))
      END IF
      IF ((ng) > 0) THEN
        ALLOCATE(rismt%vlgz(ng))
      END IF
      !
    ELSE IF (itype == ITYPE_LAUERISM) THEN
      IF ((nrzs * ngxy * nsite) > 0) THEN
        ALLOCATE(rismt%csgz(nrzs * ngxy, nsite))
        ALLOCATE(rismt%hgz( nrzs * ngxy, nsite))
      END IF
      IF ((nrzl * ngxy) > 0) THEN
        ALLOCATE(rismt%vlgz(nrzl * ngxy))
      END IF
      IF (ngxy > 0) THEN
        ALLOCATE(rismt%vright(ngxy))
        ALLOCATE(rismt%vleft( ngxy))
        ALLOCATE(rismt%do_vright(ngxy))
        ALLOCATE(rismt%do_vleft( ngxy))
      END IF
      IF ((nrzl * ngxy * nsite) > 0) THEN
        ALLOCATE(rismt%hsgz(nrzl * ngxy, nsite))
        ALLOCATE(rismt%hlgz(nrzl * ngxy, nsite))
      END IF
      IF (nsite > 0) THEN
        ALLOCATE(rismt%cda(nsite))
      END IF
      IF (nrzl > 0) THEN
        ALLOCATE(rismt%cdz(nrzl))
      END IF
      IF ((nrzl * nsite * nsite_t) > 0) THEN
        ALLOCATE(rismt%hdz(nrzl, nsite, nsite_t))
      END IF
      IF ((nrzl * nsite) > 0) THEN
        ALLOCATE(rismt%usg0( nrzl, nsite))
        ALLOCATE(rismt%csg0( nrzl, nsite))
        ALLOCATE(rismt%csdg0(nrzl, nsite))
        ALLOCATE(rismt%hg0(  nrzl, nsite))
        ALLOCATE(rismt%gg0(  nrzl, nsite))
      END IF
    END IF
    !
    ! ..... susceptibility, x12(g)
    CALL allocate_suscept(rismt, itype, nsite, nsite_t, nrzl, ngs, lboth)
    !
    ! ..... charges
    IF (itype == ITYPE_3DRISM .OR. itype == ITYPE_LAUERISM) THEN
      IF (nsite > 0) THEN
        ALLOCATE(rismt%nsol(nsite))
        ALLOCATE(rismt%qsol(nsite))
      END IF
    END IF
    !
    ! ..... chemical potentials
    IF (nsite > 0) THEN
      ALLOCATE(rismt%usol(nsite))
      ALLOCATE(rismt%usol_GF(nsite))
    END IF
    !
    ! ..... fields of solvent
    IF (itype == ITYPE_3DRISM) THEN
      IF (ng > 0) THEN
        ALLOCATE(rismt%rhog(ng))
        ALLOCATE(rismt%vpot(ng))
      END IF
      !
    ELSE IF (itype == ITYPE_LAUERISM) THEN
      IF ((nrzl * ngxy) > 0) THEN
        ALLOCATE(rismt%rhog(nrzl * ngxy))
        ALLOCATE(rismt%vpot(nrzl * ngxy))
      END IF
      IF (ng > 0) THEN
        ALLOCATE(rismt%rhog_pbc(ng))
        ALLOCATE(rismt%vpot_pbc(ng))
      END IF
    END IF
    !
  END SUBROUTINE allocate_rism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_suscept(rismt, itype, nsite, nsite_t, nrzl, ngs, lboth)
    !--------------------------------------------------------------------------
    !
    ! ... initialize rism_type
    ! ...
    ! ... Variables:
    ! ...   itype:   data type (ITYPE_1DRISM, ITYPE_3DRISM or ITYPE_LAUERISM)
    ! ...   nsite:   number of solvent's sites (i.e. degrees of freedom for correlation functions)
    ! ...   nsite_t: sum of nsite over all processies
    ! ...   nrzl:    dimension of R-space grid along long Z-stick
    ! ...   ngs:     number of shells in G-space (in Gxy-space for Laue-RISM)
    ! ...   lboth:   apply both-hands calculation, or not (for Laue-RISM)
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    INTEGER,         INTENT(IN)    :: itype
    INTEGER,         INTENT(IN)    :: nsite
    INTEGER,         INTENT(IN)    :: nsite_t
    INTEGER,         INTENT(IN)    :: nrzl
    INTEGER,         INTENT(IN)    :: ngs
    LOGICAL,         INTENT(IN)    :: lboth
    !
    ! ... deallocate memory, if needed
    IF (ASSOCIATED(rismt%xgs))  DEALLOCATE(rismt%xgs)
    IF (ASSOCIATED(rismt%xgs0)) DEALLOCATE(rismt%xgs0)
    IF (ASSOCIATED(rismt%xgs1)) DEALLOCATE(rismt%xgs1)
    IF (ASSOCIATED(rismt%ygs))  DEALLOCATE(rismt%ygs)
    IF (ASSOCIATED(rismt%ygs0)) DEALLOCATE(rismt%ygs0)
    IF (ASSOCIATED(rismt%ygs1)) DEALLOCATE(rismt%ygs1)
    !
    ! ... susceptibility, x12(g)
    IF (itype == ITYPE_3DRISM) THEN
      IF ((ngs * nsite * nsite_t) > 0) THEN
        ALLOCATE(rismt%xgs(ngs, nsite, nsite_t))
      END IF
      !
    ELSE IF (itype == ITYPE_LAUERISM) THEN
      IF ((nrzl * ngs * nsite * nsite_t) > 0) THEN
        ALLOCATE(rismt%xgs(nrzl * ngs, nsite, nsite_t))
        IF (lboth) THEN
          ALLOCATE(rismt%ygs(nrzl * ngs, nsite, nsite_t))
        END IF
      END IF
      IF ((nrzl * nsite * nsite_t) > 0) THEN
        ALLOCATE(rismt%xgs0(nrzl, nsite, nsite_t))
        ALLOCATE(rismt%xgs1(nrzl, nsite, nsite_t))
        IF (lboth) THEN
          ALLOCATE(rismt%ygs0(nrzl, nsite, nsite_t))
          ALLOCATE(rismt%ygs1(nrzl, nsite, nsite_t))
        END IF
      END IF
    END IF
    !
  END SUBROUTINE allocate_suscept
  !
  !--------------------------------------------------------------------------
  SUBROUTINE refresh_suscept_3drism(rismt)
    !--------------------------------------------------------------------------
    !
    ! ... refresh inter-site susceptibility for 3D-RISM
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    !
    INTEGER :: nsite
    INTEGER :: msite
    INTEGER :: ngs
    !
    ! ... get dimensions
    nsite = rismt%mp_site%nsite
    msite = rismt%nsite
    ngs   = rismt%gvec%ngl
    !
    ! ... check parameters
    IF (nsite <= 0) THEN
      CALL errore(' refresh_suscept_3drism ', ' too small nsite ', 1)
    END IF
    !
    IF (msite < 0) THEN
      CALL errore(' refresh_suscept_3drism ', ' msite is negative ', 1)
    END IF
    !
    IF (ngs <= 0) THEN
      CALL errore(' refresh_suscept_3drism ', ' too small ngs ', 1)
    END IF
    !
    ! ... refresh susceptibility
    rismt%ngs = ngs
    CALL allocate_suscept(rismt, ITYPE_3DRISM, msite, nsite, 0, ngs, .FALSE.)
    !
  END SUBROUTINE refresh_suscept_3drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE refresh_suscept_lauerism(rismt, lboth)
    !--------------------------------------------------------------------------
    !
    ! ... refresh inter-site susceptibility for Laue-RISM
    ! ...
    ! ... Variables:
    ! ...   lboth: apply both-hands calculation, or not
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    LOGICAL,         INTENT(IN)    :: lboth
    !
    INTEGER :: nsite
    INTEGER :: msite
    INTEGER :: nrzl
    INTEGER :: ngs
    !
    ! ... get dimensions
    nsite = rismt%mp_site%nsite
    msite = rismt%nsite
    nrzl  = rismt%nrzl
    ngs   = rismt%lfft%nglxy
    !
    ! ... check parameters
    IF (nsite <= 0) THEN
      CALL errore(' refresh_suscept_lauerism ', ' too small nsite ', 1)
    END IF
    !
    IF (msite < 0) THEN
      CALL errore(' refresh_suscept_lauerism ', ' msite is negative ', 1)
    END IF
    !
    IF (nrzl <= 0) THEN
      CALL errore(' refresh_suscept_lauerism ', ' too small nrzl ', 1)
    END IF
    !
    IF (ngs <= 0) THEN
      CALL errore(' refresh_suscept_lauerism ', ' too small ngs ', 1)
    END IF
    !
    ! ... refresh susceptibility
    rismt%ngs = ngs
    CALL allocate_suscept(rismt, ITYPE_LAUERISM, msite, nsite, nrzl, ngs, lboth)
    !
  END SUBROUTINE refresh_suscept_lauerism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_rism(rismt, lall)
    !--------------------------------------------------------------------------
    !
    ! ... finalize rism_type
    !
    ! ... if lall=.TRUE., deallocate all data.
    ! ... if lall=.FALSE., deallocate data, which depend on FFT box.
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    LOGICAL,         INTENT(IN)    :: lall
    !
    ! ... finalize MPI
    IF (lall) THEN
      CALL mp_end_rism(rismt%mp_site, rismt%mp_task)
    END IF
    !
    IF (rismt%itype == ITYPE_1DRISM) THEN
      !
      ! ... finalize radial FFT
      CALL deallocate_radfft(rismt%rfft)
      !
    ELSE IF (rismt%itype == ITYPE_3DRISM) THEN
      !
      ! ... finalize 3D-FFT(custom)
      CALL fft_type_deallocate(rismt%dfft)
      !
      ! ... finalize G-vectors(custom)
      CALL deallocate_gvec_3drism(rismt%gvec)
      !
    ELSE IF (rismt%itype == ITYPE_LAUERISM) THEN
      !
      ! ... finalize 3D-FFT(custom)
      CALL fft_type_deallocate(rismt%dfft)
      !
      ! ... finalize G-vectors(custom)
      CALL deallocate_gvec_3drism(rismt%gvec)
      !
      ! ... finalize Laue-FFT
      CALL deallocate_lauefft(rismt%lfft)
      !
    END IF
    !
    ! ... set variables to be 0
    IF (lall) THEN
      rismt%avail    = .FALSE.
      rismt%itype    = ITYPE_NULL
      rismt%closure  = CLOSURE_NULL
      rismt%temp     = 0.0_DP
      rismt%tau      = 0.0_DP
      rismt%nsite    = 0
      rismt%qtot     = 0.0_DP
      rismt%esol     = 0.0_DP
      rismt%vsol     = 0.0_DP
      rismt%pbc_nfit = 0
    END IF
    rismt%nr   = 0
    rismt%nrzs = 0
    rismt%nrzl = 0
    rismt%ng   = 0
    rismt%ngs  = 0
    rismt%ngxy = 0
    !
    ! ... deallocate arrays
    IF (lall) THEN
      IF (ASSOCIATED(rismt%nsol))    DEALLOCATE(rismt%nsol)
      IF (ASSOCIATED(rismt%qsol))    DEALLOCATE(rismt%qsol)
      IF (ASSOCIATED(rismt%usol))    DEALLOCATE(rismt%usol)
      IF (ASSOCIATED(rismt%usol_GF)) DEALLOCATE(rismt%usol_GF)
    END IF
    IF (ASSOCIATED(rismt%csr))       DEALLOCATE(rismt%csr)
    IF (ASSOCIATED(rismt%csg))       DEALLOCATE(rismt%csg)
    IF (ASSOCIATED(rismt%csgz))      DEALLOCATE(rismt%csgz)
    IF (ASSOCIATED(rismt%csg0))      DEALLOCATE(rismt%csg0)
    IF (ASSOCIATED(rismt%csdr))      DEALLOCATE(rismt%csdr)
    IF (ASSOCIATED(rismt%csdg0))     DEALLOCATE(rismt%csdg0)
    IF (ASSOCIATED(rismt%cda))       DEALLOCATE(rismt%cda)
    IF (ASSOCIATED(rismt%cdz))       DEALLOCATE(rismt%cdz)
    IF (ASSOCIATED(rismt%uljr))      DEALLOCATE(rismt%uljr)
    IF (ASSOCIATED(rismt%uwr))       DEALLOCATE(rismt%uwr)
    IF (ASSOCIATED(rismt%usr))       DEALLOCATE(rismt%usr)
    IF (ASSOCIATED(rismt%usg0))      DEALLOCATE(rismt%usg0)
    IF (ASSOCIATED(rismt%vsr))       DEALLOCATE(rismt%vsr)
    IF (ASSOCIATED(rismt%ulr))       DEALLOCATE(rismt%ulr)
    IF (ASSOCIATED(rismt%vlr))       DEALLOCATE(rismt%vlr)
    IF (ASSOCIATED(rismt%ulg))       DEALLOCATE(rismt%ulg)
    IF (ASSOCIATED(rismt%ulgz))      DEALLOCATE(rismt%ulgz)
    IF (ASSOCIATED(rismt%vlgz))      DEALLOCATE(rismt%vlgz)
    IF (ASSOCIATED(rismt%vright))    DEALLOCATE(rismt%vright)
    IF (ASSOCIATED(rismt%vleft))     DEALLOCATE(rismt%vleft)
    IF (ASSOCIATED(rismt%do_vright)) DEALLOCATE(rismt%do_vright)
    IF (ASSOCIATED(rismt%do_vleft))  DEALLOCATE(rismt%do_vleft)
    IF (ASSOCIATED(rismt%hr ))       DEALLOCATE(rismt%hr)
    IF (ASSOCIATED(rismt%hg ))       DEALLOCATE(rismt%hg)
    IF (ASSOCIATED(rismt%hgz ))      DEALLOCATE(rismt%hgz)
    IF (ASSOCIATED(rismt%hg0))       DEALLOCATE(rismt%hg0)
    IF (ASSOCIATED(rismt%hsgz))      DEALLOCATE(rismt%hsgz)
    IF (ASSOCIATED(rismt%hlgz))      DEALLOCATE(rismt%hlgz)
    IF (ASSOCIATED(rismt%hdz))       DEALLOCATE(rismt%hdz)
    IF (ASSOCIATED(rismt%gr ))       DEALLOCATE(rismt%gr)
    IF (ASSOCIATED(rismt%gg0))       DEALLOCATE(rismt%gg0)
    IF (ASSOCIATED(rismt%wg ))       DEALLOCATE(rismt%wg)
    IF (ASSOCIATED(rismt%zg ))       DEALLOCATE(rismt%zg)
    IF (ASSOCIATED(rismt%xgs))       DEALLOCATE(rismt%xgs)
    IF (ASSOCIATED(rismt%xgs0))      DEALLOCATE(rismt%xgs0)
    IF (ASSOCIATED(rismt%xgs1))      DEALLOCATE(rismt%xgs1)
    IF (ASSOCIATED(rismt%ygs))       DEALLOCATE(rismt%ygs)
    IF (ASSOCIATED(rismt%ygs0))      DEALLOCATE(rismt%ygs0)
    IF (ASSOCIATED(rismt%ygs1))      DEALLOCATE(rismt%ygs1)
    IF (ASSOCIATED(rismt%rhog))      DEALLOCATE(rismt%rhog)
    IF (ASSOCIATED(rismt%vpot))      DEALLOCATE(rismt%vpot)
    IF (ASSOCIATED(rismt%rhog_pbc))  DEALLOCATE(rismt%rhog_pbc)
    IF (ASSOCIATED(rismt%vpot_pbc))  DEALLOCATE(rismt%vpot_pbc)
    !
  END SUBROUTINE deallocate_rism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE clean_rism_data(rismt)
    !--------------------------------------------------------------------------
    !
    ! ... set zero to correlation functions,
    ! ... which are variable through iterations of RISM.
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    !
    REAL(DP),    PARAMETER :: R_ZERO = 0.0_DP
    COMPLEX(DP), PARAMETER :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    !
    rismt%qtot = R_ZERO
    rismt%esol = R_ZERO
    rismt%vsol = R_ZERO
    IF (ASSOCIATED(rismt%csr))      rismt%csr      = R_ZERO
    IF (ASSOCIATED(rismt%csg))      rismt%csg      = R_ZERO
    IF (ASSOCIATED(rismt%csgz))     rismt%csgz     = C_ZERO
    IF (ASSOCIATED(rismt%csg0))     rismt%csg0     = R_ZERO
    IF (ASSOCIATED(rismt%csdr))     rismt%csdr     = R_ZERO
    IF (ASSOCIATED(rismt%csdg0))    rismt%csdg0    = R_ZERO
    IF (ASSOCIATED(rismt%cda))      rismt%cda      = R_ZERO
    IF (ASSOCIATED(rismt%hr ))      rismt%hr       = R_ZERO
    IF (ASSOCIATED(rismt%hg ))      rismt%hg       = R_ZERO
    IF (ASSOCIATED(rismt%hgz ))     rismt%hgz      = C_ZERO
    IF (ASSOCIATED(rismt%hg0))      rismt%hg0      = R_ZERO
    IF (ASSOCIATED(rismt%hsgz))     rismt%hsgz     = C_ZERO
    IF (ASSOCIATED(rismt%hlgz))     rismt%hlgz     = C_ZERO
    IF (ASSOCIATED(rismt%gr ))      rismt%gr       = R_ZERO
    IF (ASSOCIATED(rismt%gg0))      rismt%gg0      = R_ZERO
    IF (ASSOCIATED(rismt%nsol))     rismt%nsol     = R_ZERO
    IF (ASSOCIATED(rismt%qsol))     rismt%qsol     = R_ZERO
    IF (ASSOCIATED(rismt%usol))     rismt%usol     = R_ZERO
    IF (ASSOCIATED(rismt%usol_GF))  rismt%usol_GF  = R_ZERO
    IF (ASSOCIATED(rismt%rhog))     rismt%rhog     = C_ZERO
    IF (ASSOCIATED(rismt%vpot))     rismt%vpot     = C_ZERO
    IF (ASSOCIATED(rismt%rhog_pbc)) rismt%rhog_pbc = C_ZERO
    IF (ASSOCIATED(rismt%vpot_pbc)) rismt%vpot_pbc = C_ZERO
    !
  END SUBROUTINE clean_rism_data
  !
  !--------------------------------------------------------------------------
  FUNCTION get_chempot_type(rismt) result(ichempot)
    !--------------------------------------------------------------------------
    !
    ! ... set zero to correlation functions
    !
    IMPLICIT NONE
    !
    INTEGER                        :: ichempot
    TYPE(rism_type), INTENT(INOUT) :: rismt
    !
    SELECT CASE (rismt%closure)
    CASE (CLOSURE_HNC)
      ichempot = CHEMPOT_HNC
      !
    CASE (CLOSURE_KH)
      ichempot = CHEMPOT_KH
      !
    CASE DEFAULT
      ichempot = CHEMPOT_NULL
      !
    END SELECT
    !
  END FUNCTION get_chempot_type
  !
END MODULE rism

