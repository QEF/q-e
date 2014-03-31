!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE paw_gipaw

  USE kinds,        ONLY: DP
  USE parameters,   ONLY: npsx
  use radial_grids, ONLY: ndmx
  !!use gipaw_module, ONLY: nbrx
  !
  ! ... These parameters are needed for the paw variables
  !
  SAVE
  !
!  REAL(DP) :: &
!       paw_betar(ndmx,nbrx,npsx)  ! radial beta_{mu} functions
!  INTEGER :: &
!       paw_nh(npsx),             &! number of beta functions per atomic type
!       paw_nbeta(npsx),          &! number of beta functions
!       paw_kkbeta(npsx),         &! point where the beta are zero
!       paw_lll(nbrx,npsx)         ! angular momentum of the beta function
!  INTEGER :: &
!       paw_nhm,              &! max number of different beta functions per atom
!       paw_nkb,              &! total number of beta functions, with st.fact.
!       paw_nqxq,             &! size of interpolation table
!       paw_lmaxkb,           &! max angular momentum
!       paw_lmaxq,              &! max angular momentum + 1 for Q functions
!       paw_nqx                ! number of interpolation points
!  INTEGER, ALLOCATABLE ::&
!       paw_indv(:,:),        &! correspondence of betas atomic <-> soli
!       paw_nhtol(:,:),       &! correspondence n <-> angular momentum
!       paw_nhtom(:,:),       &! correspondence n <-> magnetic angular m
!       paw_nl(:,:),          &! number of projectors for each l
!       paw_iltonh(:,:,:)        ! corresp l, num <--> n for each type
!  complex(DP), ALLOCATABLE, TARGET :: &
!       paw_vkb(:,:),         &   ! all beta functions in reciprocal space
!       paw_becp(:,:)             !  products of wavefunctions and proj
!  REAL(DP), ALLOCATABLE :: &
!       paw_tab(:,:,:)              ! interpolation table for PPs
!  !<ceres>
!  REAL(DP), ALLOCATABLE :: &
!       paw_tab_d2y(:,:,:)          ! for cubic splines
!  !</ceres>

  type wfc_label
     integer  :: na , &   ! Atom number
          nt ,        &   ! Type
          n  ,        &   ! Chi index
          l  ,        &   ! l
          m  ,        &   ! m
          nrc,        &   ! index of core radius on mesh
          nrs             ! index of inner core radius (where "f_step" starts)
     real(DP) :: rc  ! paw core radius
  end type wfc_label

  type at_wfc
     type(wfc_label)          :: label
     integer                  :: kkpsi
     real(DP)  , pointer :: psi(:)
  end type at_wfc

!  TYPE ( at_wfc ), POINTER :: aephi(:,:)
!  TYPE ( at_wfc ), POINTER :: psphi(:,:)
!
!  LOGICAL, ALLOCATABLE :: vloc_present ( : )
!  REAL(dp), ALLOCATABLE :: gipaw_ae_vloc ( :, : )
!  REAL(dp), ALLOCATABLE :: gipaw_ps_vloc ( :, : )
!
!  LOGICAL, ALLOCATABLE :: gipaw_data_in_upf_file ( : )
!
!  INTEGER, ALLOCATABLE :: gipaw_ncore_orbital ( : )
!  REAL(dp), ALLOCATABLE :: gipaw_core_orbital ( :, :, : )

  integer, parameter :: nbrx = 20

  INTEGER :: &
       paw_nkb,            & ! total number of beta functions, with st.fact.
       paw_nqxq,           & ! size of interpolation table
       paw_lmaxkb,         & ! max angular momentum
       paw_lmaxq,          & ! max angular momentum + 1 for Q functions
       paw_nqx               ! number of interpolation points

  complex(DP), ALLOCATABLE :: &
       paw_vkb(:,:),       & ! all beta functions in reciprocal space
       paw_becp(:,:)         !  products of wavefunctions and proj

  TYPE paw_recon_type
     REAL(DP) :: &
          paw_betar(ndmx,nbrx)  ! radial beta_{mu} functions
     INTEGER :: &
          paw_nh,             & ! number of beta functions per atomic type
          paw_nbeta,          & ! number of beta functions
          paw_kkbeta,         & ! point where the beta are zero
          paw_lll(nbrx)         ! angular momentum of the beta function
     INTEGER, POINTER ::&
          paw_indv(:),        & ! correspondence of betas atomic <-> soli
          paw_nhtol(:),       & ! correspondence n <-> angular momentum
          paw_nhtom(:),       & ! correspondence n <-> magnetic angular m
          paw_nl(:),          & ! number of projectors for each l
          paw_iltonh(:,:)       ! corresp l, num <--> n for each type
     REAL(DP), POINTER :: &
          paw_tab(:,:)          ! interpolation table for PPs

     REAL(DP), POINTER :: &
          paw_tab_d2y(:,:)      ! for cubic splines

     TYPE ( at_wfc ), POINTER :: aephi(:)
     TYPE ( at_wfc ), POINTER :: psphi(:)

     LOGICAL :: vloc_present
     REAL(dp), POINTER :: gipaw_ae_vloc(:)
     REAL(dp), POINTER :: gipaw_ps_vloc(:)

     LOGICAL :: gipaw_data_in_upf_file

     INTEGER :: gipaw_ncore_orbital
     REAL(dp), POINTER :: gipaw_core_orbital(:,:)
     INTEGER, POINTER :: gipaw_core_orbital_l(:)

  END TYPE paw_recon_type

  TYPE ( paw_recon_type ), ALLOCATABLE :: paw_recon(:)

CONTAINS

  SUBROUTINE paw_wfc_init(phi)
    !
    ! Initialize default values for labe end kkpsi
    !
    implicit none
    ! Argument
    TYPE ( at_wfc ), INTENT ( INOUT ) :: phi(:)

    phi(:)%label%na  = 0
    phi(:)%label%nt  = 0
    phi(:)%label%n   = 0
    phi(:)%label%l   = -99
    phi(:)%label%m   = -99
    phi(:)%label%nrc = 0
    phi(:)%kkpsi     = 0

  END SUBROUTINE paw_wfc_init

  !****************************************************************************

  subroutine read_recon ( filerec_sp, jtyp, paw_recon_sp )

    !
    ! Read all-electron and pseudo atomic wavefunctions
    !  needed for PAW reconstruction
    !

    use read_upf_v1_module,  only : scan_begin, scan_end
    use atom,               only : rgrid
    USE io_global,          ONLY : stdout

    IMPLICIT NONE

    ! Arguments
    CHARACTER ( LEN = 256 ), INTENT ( IN ) :: filerec_sp
    INTEGER, INTENT ( IN ) :: jtyp
    TYPE ( paw_recon_type ), INTENT ( INOUT ) :: paw_recon_sp

    ! Local
    INTEGER :: j, i, kkphi

    ! If the data has already been read from a UPF file
    IF ( paw_recon_sp%gipaw_data_in_upf_file ) RETURN

    OPEN ( 14, FILE = filerec_sp )
    CALL scan_begin ( 14, 'PAW', .true. )
    READ(14,*) paw_recon_sp%paw_nbeta
    CALL scan_end ( 14, 'PAW' )
    CLOSE ( 14 )

    ALLOCATE ( paw_recon_sp%psphi(paw_recon_sp%paw_nbeta) )
    ALLOCATE ( paw_recon_sp%aephi(paw_recon_sp%paw_nbeta) )

    CALL paw_wfc_init ( paw_recon_sp%psphi )
    CALL paw_wfc_init ( paw_recon_sp%aephi )

    OPEN ( 14, FILE = filerec_sp )

    WRITE (stdout,*) "N_AEwfc atom",jtyp,":", paw_recon_sp%paw_nbeta

    recphi_loop: DO i = 1, paw_recon_sp%paw_nbeta
       ALLOCATE ( paw_recon_sp%aephi(i)%psi(rgrid(jtyp)%mesh) )
       paw_recon_sp%aephi(i)%label%nt=jtyp
       paw_recon_sp%aephi(i)%label%n=i
       CALL scan_begin(14,'REC',.false.)
       CALL scan_begin(14,'kkbeta',.false.)
       read(14,*)  kkphi
       CALL scan_end(14,'kkbeta')
       paw_recon_sp%aephi(i)%kkpsi=kkphi
       CALL scan_begin(14,'L',.false.)
       read(14,*) paw_recon_sp%aephi(i)%label%l
       CALL scan_end(14,'L')
       CALL scan_begin(14,'REC_AE',.false.)
       read(14,*) ( paw_recon_sp%aephi(i)%psi(j),j=1,kkphi)
       CALL scan_end(14,'REC_AE')
       ALLOCATE ( paw_recon_sp%psphi(i)%psi(rgrid(jtyp)%mesh) )
       paw_recon_sp%psphi(i)%label%nt = jtyp
       paw_recon_sp%psphi(i)%label%n = i
       paw_recon_sp%psphi(i)%label%l = paw_recon_sp%aephi(i)%label%l
       paw_recon_sp%psphi(i)%kkpsi = kkphi
       CALL scan_begin(14,'REC_PS',.false.)
       READ(14,*) ( paw_recon_sp%psphi(i)%psi(j),j=1,kkphi)
       CALL scan_end(14,'REC_PS')
       CALL scan_end(14,'REC')
    END DO recphi_loop
    CLOSE(14)

  END SUBROUTINE read_recon

  !****************************************************************************


  subroutine set_paw_upf ( is, upf)
  !
  ! interface between the UPF pseudo type and the internal representation
  ! of the PAW-related variables
  !<apsi>
  USE pseudo_types
  USE ions_base, ONLY :  nsp
  USE atom, ONLY: rgrid
  implicit none
  !
  INTEGER, INTENT(IN) :: is
  TYPE (pseudo_upf) :: upf
  !
  INTEGER :: nb
  !
  IF ( upf%has_gipaw ) THEN
     IF ( .NOT. ALLOCATED ( paw_recon ) ) THEN !CG
        ALLOCATE ( paw_recon(nsp) )
        paw_recon(:)%gipaw_data_in_upf_file = .FALSE.
     END IF

     paw_recon(is)%paw_nbeta = upf%gipaw_wfs_nchannels
     paw_recon(is)%vloc_present = .TRUE.
     paw_recon(is)%gipaw_data_in_upf_file = .TRUE.

     paw_recon(is)%gipaw_ncore_orbital = upf%gipaw_ncore_orbitals
     ALLOCATE ( paw_recon(is)%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
     ALLOCATE ( paw_recon(is)%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
     paw_recon(is)%gipaw_core_orbital(:upf%mesh,:upf%gipaw_ncore_orbitals) &
          = upf%gipaw_core_orbital(:upf%mesh,:upf%gipaw_ncore_orbitals)
     paw_recon(is)%gipaw_core_orbital_l(:upf%gipaw_ncore_orbitals) &
          = upf%gipaw_core_orbital_l(:upf%gipaw_ncore_orbitals)

     ALLOCATE ( paw_recon(is)%gipaw_ae_vloc(upf%mesh) )
     ALLOCATE ( paw_recon(is)%gipaw_ps_vloc(upf%mesh) )
     paw_recon(is)%gipaw_ae_vloc(:upf%mesh) = upf%gipaw_vlocal_ae(:upf%mesh)
     paw_recon(is)%gipaw_ps_vloc(:upf%mesh) = upf%gipaw_vlocal_ps(:upf%mesh)

     ALLOCATE ( paw_recon(is)%aephi(upf%gipaw_wfs_nchannels) )
     ALLOCATE ( paw_recon(is)%psphi(upf%gipaw_wfs_nchannels) )

     DO nb = 1, upf%gipaw_wfs_nchannels
        ALLOCATE ( paw_recon(is)%aephi(nb)%psi(rgrid(is)%mesh) )
        paw_recon(is)%aephi(nb)%label%nt = is
        paw_recon(is)%aephi(nb)%label%n = nb
        paw_recon(is)%aephi(nb)%label%l = upf%gipaw_wfs_ll(nb)
        !paw_recon(is)%aephi(nb)%label%m =
        paw_recon(is)%aephi(nb)%label%nrc = upf%mesh
        paw_recon(is)%aephi(nb)%kkpsi = upf%mesh
        paw_recon(is)%aephi(nb)%label%rc = upf%gipaw_wfs_rcut(nb)
        paw_recon(is)%aephi(nb)%psi(:upf%mesh) = upf%gipaw_wfs_ae(:upf%mesh,nb)

        ALLOCATE ( paw_recon(is)%psphi(nb)%psi(rgrid(is)%mesh) )
        paw_recon(is)%psphi(nb)%label%nt = is
        paw_recon(is)%psphi(nb)%label%n = nb
        paw_recon(is)%psphi(nb)%label%l = upf%gipaw_wfs_ll(nb)
        !paw_recon(is)%psphi(nb)%label%m =
        paw_recon(is)%psphi(nb)%label%nrc = upf%mesh
        paw_recon(is)%psphi(nb)%kkpsi = upf%mesh
        paw_recon(is)%psphi(nb)%label%rc = upf%gipaw_wfs_rcutus(nb)
        paw_recon(is)%psphi(nb)%psi(:upf%mesh) = upf%gipaw_wfs_ps(:upf%mesh,nb)
     END DO
  END IF
  !</apsi>
  END SUBROUTINE set_paw_upf

END MODULE paw_gipaw
