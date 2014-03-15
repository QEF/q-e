!
! Copyright (C) 2008-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
      MODULE read_upf_v2_module
!=----------------------------------------------------------------------------=!
!  this module handles the reading of pseudopotential data

! ...   declare modules
        USE kinds,        ONLY: DP
        USE pseudo_types, ONLY: pseudo_upf
        USE radial_grids, ONLY: radial_grid_type
        USE parser,       ONLY : version_compare
        USE iotk_module
        !
        PRIVATE
        PUBLIC :: read_upf_v2
 CONTAINS

!------------------------------------------------+
SUBROUTINE read_upf_v2(u, upf, grid, ierr)             !
   !---------------------------------------------+
   ! Read pseudopotential in UPF format version 2, uses iotk
   !
   USE pseudo_types, ONLY: nullify_pseudo_upf, deallocate_pseudo_upf
   USE radial_grids, ONLY: radial_grid_type, nullify_radial_grid
   IMPLICIT NONE
   INTEGER,INTENT(IN)             :: u         ! i/o unit
   TYPE(pseudo_upf),INTENT(INOUT) :: upf       ! the pseudo data
   TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
   !
   INTEGER,OPTIONAL,INTENT(OUT):: ierr      ! /= 0 if something went wrong
   CHARACTER(len=iotk_namlenx) :: root
   CHARACTER(len=iotk_attlenx) :: attr
   INTEGER :: ierr_
   LOGICAL :: found
   LOGICAL,EXTERNAL :: matches
   CHARACTER(len=6),PARAMETER :: max_version = '2.0.1'
   !
   ! Prepare the type .  Should be done where upf is instantiated
   ! CALL deallocate_pseudo_upf(upf)
   ! CALL nullify_pseudo_upf(upf)
   !
   ! IF(present(grid)) call nullify_radial_grid(grid)
   ! nullify(upf%grid)
   !
   ! Initialize the file
   CALL iotk_open_read(u, attr=attr, root=root, ierr=ierr_)
   !
   IF((abs(ierr_)>0) .or. .not. matches('UPF',root) ) THEN
       !
       CALL iotk_close_read(u,ierr=ierr_)
       IF(.not. present(ierr)) &
         CALL errore('read_upf_v2','Cannot open UPF file.',1)
       ierr = 1
       RETURN
   ENDIF
   CALL iotk_scan_attr(attr, 'version', upf%nv)
   IF (version_compare(upf%nv, max_version) == 'newer') &
       CALL errore('read_upf_v2',&
                   'Unknown UPF format version: '//TRIM(upf%nv),1)
   !
   ! Skip human-readable header
   CALL iotk_scan_begin(u,'PP_INFO',found=found)
   if(found) CALL iotk_scan_end(u,'PP_INFO')
   !
   ! Read machine-readable header
   CALL read_upf_header(u, upf)
   IF(upf%tpawp .and. .not. present(grid)) &
      CALL errore('read_upf_v2', 'PAW requires a radial_grid_type.', 1)
   !
   ! CHECK for bug in version 2.0.0 of UPF file
   IF ( version_compare(upf%nv, '2.0.1') == 'older' .and. upf%tvanp .and.  &
        .not. upf%tpawp ) CALL errore('read_upf_v2',&
                   'Ultrasoft pseudopotentials in UPF format v.2.0.0 are &
                  & affected by a bug compromising their quality. Please &
                  & regenerate pseudopotential file for '//TRIM(upf%psd), 1)

   ! Read radial grid mesh
   CALL read_upf_mesh(u, upf, grid)
   ! Read non-linear core correction charge
   ALLOCATE( upf%rho_atc(upf%mesh) )
   IF(upf%nlcc) THEN
      CALL iotk_scan_dat(u, 'PP_NLCC',  upf%rho_atc)
   ELSE
      ! A null core charge simplifies several functions, mostly in PAW
      upf%rho_atc(1:upf%mesh) = 0._dp
   ENDIF

   ! Read local potential
   IF(.not. upf%tcoulombp) THEN
      ALLOCATE( upf%vloc(upf%mesh) )
      CALL iotk_scan_dat(u, 'PP_LOCAL', upf%vloc)
   ENDIF
   ! Read nonlocal components: projectors, augmentation, hamiltonian elements
   CALL read_upf_nonlocal(u, upf)
   ! Read initial pseudo wavefunctions
   ! (usually only wfcs with occupancy > 0)
   CALL read_upf_pswfc(u, upf)
   ! Read all-electron and pseudo wavefunctions
   CALL read_upf_full_wfc(u, upf)
   ! Read valence atomic density (used for initial density)
   ALLOCATE( upf%rho_at(upf%mesh) )
   CALL iotk_scan_dat(u, 'PP_RHOATOM', upf%rho_at)
   ! Read additional info for full-relativistic calculation
   CALL read_upf_spin_orb(u, upf)
   ! Read additional data for PAW (All-electron charge, wavefunctions, vloc..)
   CALL read_upf_paw(u, upf)
   ! Read data for gipaw reconstruction
   CALL read_upf_gipaw(u, upf)
   !
   ! Close the file (not the unit!)
   CALL iotk_close_read(u)
   !
   IF( present(ierr) ) ierr=0

   RETURN

   CONTAINS
   !
   SUBROUTINE read_upf_header(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)             :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER                     :: ierr ! /= 0 if something went wrong
      CHARACTER(len=iotk_attlenx) :: attr
      CHARACTER(len=256) :: dft_buffer  ! needed to allow the string defining the
                                        ! DFT flavor to be longer than upf%dft 
                                        ! (currntly 25) without getting iotk upset. 
                                        ! An error message is issued if trimmed 
                                        ! dft_buffer exceeds upf%dft size.
      INTEGER :: len_buffer
      !
      INTEGER :: nw
      !
      ! Read HEADER section with some initialization data
      CALL iotk_scan_empty(u, 'PP_HEADER', attr=attr)
         CALL iotk_scan_attr(attr, 'generated',      upf%generated, default=' ')
         CALL iotk_scan_attr(attr, 'author',         upf%author,    default='anonymous')
         CALL iotk_scan_attr(attr, 'date',           upf%date,      default=' ')
         CALL iotk_scan_attr(attr, 'comment',        upf%comment,   default=' ')
         !
         CALL iotk_scan_attr(attr, 'element',        upf%psd)
         CALL iotk_scan_attr(attr, 'pseudo_type',    upf%typ)
         CALL iotk_scan_attr(attr, 'relativistic',   upf%rel)
         !
         CALL iotk_scan_attr(attr, 'is_ultrasoft',   upf%tvanp)
         CALL iotk_scan_attr(attr, 'is_paw',         upf%tpawp)
         CALL iotk_scan_attr(attr, 'is_coulomb',     upf%tcoulombp, default=.false.)
         !
         CALL iotk_scan_attr(attr, 'has_so',         upf%has_so,    default=.false.)
         CALL iotk_scan_attr(attr, 'has_wfc',        upf%has_wfc,   default=upf%tpawp)
         CALL iotk_scan_attr(attr, 'has_gipaw',      upf%has_gipaw, default=.false.)
         !EMINE
         CALL iotk_scan_attr(attr, 'paw_as_gipaw',      upf%paw_as_gipaw, default=.false.)
         !
         CALL iotk_scan_attr(attr, 'core_correction',upf%nlcc)
!         CALL iotk_scan_attr(attr, 'functional',     upf%dft)
         CALL iotk_scan_attr(attr, 'functional',     dft_buffer)
         len_buffer=len_trim(dft_buffer)
         if (len_buffer > len(upf%dft)) &
            call errore('read_upf_v2','String defining DFT is too long',len_buffer)
         upf%dft=TRIM(dft_buffer)

         CALL iotk_scan_attr(attr, 'z_valence',      upf%zp)
         CALL iotk_scan_attr(attr, 'total_psenergy', upf%etotps,    default=0._dp)
         CALL iotk_scan_attr(attr, 'wfc_cutoff',     upf%ecutwfc,   default=0._dp)
         CALL iotk_scan_attr(attr, 'rho_cutoff',     upf%ecutrho,   default=0._dp)
         CALL iotk_scan_attr(attr, 'l_max',          upf%lmax,      default=0)
         CALL iotk_scan_attr(attr, 'l_max_rho',      upf%lmax_rho,  default=2*upf%lmax)
         CALL iotk_scan_attr(attr, 'l_local',        upf%lloc,      default=0)
         CALL iotk_scan_attr(attr, 'mesh_size',      upf%mesh)
         CALL iotk_scan_attr(attr, 'number_of_wfc',  upf%nwfc)
         CALL iotk_scan_attr(attr, 'number_of_proj', upf%nbeta)
      !
      !CALL iotk_scan_end(u, 'PP_HEADER')
      !CALL debug_pseudo_upf(upf)
      !
      RETURN
   END SUBROUTINE read_upf_header
   !
   SUBROUTINE read_upf_mesh(u, upf, grid)
      USE radial_grids, ONLY: allocate_radial_grid
      IMPLICIT NONE
      INTEGER,INTENT(IN)             :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
      !
      INTEGER                     :: ierr ! /= 0 if something went wrong
      CHARACTER(len=iotk_attlenx) :: attr
      LOGICAL :: found
      !
      CALL iotk_scan_begin(u, 'PP_MESH', attr=attr)

      CALL iotk_scan_attr(attr, 'dx',   upf%dx,    default=0._dp)
      CALL iotk_scan_attr(attr, 'mesh', upf%mesh,  default=upf%mesh)
      CALL iotk_scan_attr(attr, 'xmin', upf%xmin,  default=0._dp)
      CALL iotk_scan_attr(attr, 'rmax', upf%rmax,  default=0._dp)
      CALL iotk_scan_attr(attr, 'zmesh',upf%zmesh, default=0._dp)
      IF (present(grid)) THEN
         CALL allocate_radial_grid(grid, upf%mesh)
         !
         grid%dx    = upf%dx
         grid%mesh  = upf%mesh
         grid%xmin  = upf%xmin
         grid%rmax  = upf%rmax
         grid%zmesh = upf%zmesh
         !
         upf%grid => grid
         upf%r    => upf%grid%r
         upf%rab  => upf%grid%rab
      ELSE
         ALLOCATE( upf%r( upf%mesh ), upf%rab( upf%mesh ) )
      ENDIF
      !
      CALL iotk_scan_dat(u, 'PP_R',   upf%r(1:upf%mesh))
      CALL iotk_scan_dat(u, 'PP_RAB', upf%rab(1:upf%mesh))
      !
      IF (present(grid)) THEN
         ! Reconstruct additional grids
         upf%grid%r2 =  upf%r**2
         upf%grid%sqr = sqrt(upf%r)
         upf%grid%rm1 = upf%r**(-1)
         upf%grid%rm2 = upf%r**(-2)
         upf%grid%rm3 = upf%r**(-3)
      ENDIF

      CALL iotk_scan_end(u, 'PP_MESH')
      !
      RETURN
   END SUBROUTINE read_upf_mesh
   !
   SUBROUTINE read_upf_nonlocal(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)             :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nb,mb,ln,lm,l,nmb,ierr=0
      !INTEGER :: nb_=-1,mb_=-1,l_=-1,nmb_=-1
      REAL(DP):: zeros(upf%mesh)
      LOGICAL :: isnull, found
      zeros=0._dp
      !
      ! modified by AF
      !IF (upf%tcoulombp) RETURN
      IF (upf%tcoulombp) upf%nbeta = 0
      !
      ! Allocate space for non-local part
      IF ( upf%nbeta == 0) then
         upf%nqf = 0
         upf%nqlc= 0
         upf%qqq_eps= -1._dp
         upf%kkbeta = 0  
         ALLOCATE( upf%kbeta(1),         &
                   upf%lll(1),           &
                   upf%beta(upf%mesh,1), &
                   upf%dion(1,1),        &
                   upf%rinner(1),        &
                   upf%qqq(1,1),         &
                   upf%qfunc(upf%mesh,1),&
                   upf%qfcoef(1,1,1,1),  &
                   upf%rcut(1),          &
                   upf%rcutus(1),        &
                   upf%els_beta(1) )
         ! <AF>
         !CALL iotk_scan_end(u, 'PP_NONLOCAL')
         RETURN
      END IF
      !
      ! <AF>
      CALL iotk_scan_begin(u, 'PP_NONLOCAL')
      !
      ALLOCATE( upf%kbeta(upf%nbeta),          &
                upf%lll(upf%nbeta),            &
                upf%beta(upf%mesh, upf%nbeta), &
                upf%dion(upf%nbeta, upf%nbeta),&
                upf%rcut(upf%nbeta),           &
                upf%rcutus(upf%nbeta),         &
                upf%els_beta(upf%nbeta) )

      !
      ! Read the projectors:
      DO nb = 1,upf%nbeta
         CALL iotk_scan_dat(u, 'PP_BETA'//iotk_index( nb ), &
                             upf%beta(:,nb), attr=attr)
            CALL iotk_scan_attr(attr, 'label',                  upf%els_beta(nb), default='Xn')
            CALL iotk_scan_attr(attr, 'angular_momentum',       upf%lll(nb))
            CALL iotk_scan_attr(attr, 'cutoff_radius_index',    upf%kbeta(nb),    default=upf%mesh)
            CALL iotk_scan_attr(attr, 'cutoff_radius',          upf%rcut(nb),     default=0._dp)
            CALL iotk_scan_attr(attr, 'ultrasoft_cutoff_radius', upf%rcutus(nb),   default=0._dp)
!
!    Old version of UPF PPs v.2 contained an error in the tag. 
!    To be able to read the old PPs we need the following
!
            IF ( upf%rcutus(nb)==0._DP) &
            CALL iotk_scan_attr(attr,'norm_conserving_radius',upf%rcutus(nb), &
                                default=0._dp)
      ENDDO
      !
      ! Read the hamiltonian terms D_ij
      CALL iotk_scan_dat(u, 'PP_DIJ', upf%dion, attr=attr)
      !   CALL iotk_scan_attr(attr, 'non_zero_elements', upf%nd)
      !
      ! Read the augmentation charge section
      augmentation : &
      IF(upf%tvanp .or. upf%tpawp) THEN
      !
      CALL iotk_scan_begin(u, 'PP_AUGMENTATION', attr=attr)
         CALL iotk_scan_attr(attr, 'q_with_l', upf%q_with_l)
         CALL iotk_scan_attr(attr, 'nqf',      upf%nqf)
         CALL iotk_scan_attr(attr, 'nqlc',     upf%nqlc, default=2*upf%lmax+1)
         IF (upf%tpawp) THEN
            CALL iotk_scan_attr(attr,'shape',          upf%paw%augshape, default='UNKNOWN')
            CALL iotk_scan_attr(attr,'cutoff_r',       upf%paw%raug,     default=0._dp)
            CALL iotk_scan_attr(attr,'cutoff_r_index', upf%paw%iraug,    default=upf%mesh)
            CALL iotk_scan_attr(attr,'l_max_aug',      upf%paw%lmax_aug, default=upf%lmax_rho)
         ENDIF
         ! a negative number means that all qfunc are stored
         CALL iotk_scan_attr(attr,'augmentation_epsilon',upf%qqq_eps, default=-1._dp)
      !
      ALLOCATE( upf%rinner( upf%nqlc ) )
      ALLOCATE( upf%qqq   ( upf%nbeta, upf%nbeta ) )
      IF ( upf%q_with_l ) THEN
        ALLOCATE( upf%qfuncl ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2, 0:2*upf%lmax ) )
        upf%qfuncl=0._dp
      ELSE
        ALLOCATE( upf%qfunc (upf%mesh, upf%nbeta*(upf%nbeta+1)/2) )
      ENDIF
      !
      ! Read the integrals of the Q functions
      CALL iotk_scan_dat(u, 'PP_Q',upf%qqq )
      !
      ! read charge multipoles (only if PAW)
      IF( upf%tpawp ) THEN
         ALLOCATE(upf%paw%augmom(upf%nbeta,upf%nbeta, 0:2*upf%lmax))
         CALL iotk_scan_dat(u, 'PP_MULTIPOLES', upf%paw%augmom)
      ENDIF
      !
      ! Read polinomial coefficients for Q_ij expansion at small radius
      IF(upf%nqf <= 0) THEN
         upf%rinner(:) = 0._dp
         ALLOCATE( upf%qfcoef(1,1,1,1) )
         upf%qfcoef = 0._dp
      ELSE
         ALLOCATE( upf%qfcoef( MAX( upf%nqf,1 ), upf%nqlc, upf%nbeta, upf%nbeta ) )
         CALL iotk_scan_dat(u, 'PP_QFCOEF',upf%qfcoef, attr=attr)
         CALL iotk_scan_dat(u, 'PP_RINNER',upf%rinner, attr=attr)
      ENDIF
      !
      ! Read augmentation charge Q_ij
      ultrasoft_or_paw : &
      IF( upf%tvanp) THEN
         DO nb = 1,upf%nbeta
         ln = upf%lll(nb)
         DO mb = nb,upf%nbeta
         lm = upf%lll(mb)
         nmb = mb * (mb-1) /2 + nb
            q_with_l : &
            IF( upf%q_with_l ) THEN
               DO l = abs(ln-lm),ln+lm,2 ! only even terms
                  CALL iotk_scan_dat(u, 'PP_QIJL'//iotk_index((/nb,mb,l/)),&
                                    upf%qfuncl(:,nmb,l),default=zeros,attr=attr)
                  IF( upf%tpawp) upf%qfuncl(upf%paw%iraug+1:,nmb,l) = 0._dp
               ENDDO
            ELSE q_with_l 
               CALL iotk_scan_dat(u, 'PP_QIJ'//iotk_index((/nb,mb/)),&
                                 upf%qfunc(:,nmb),attr=attr,default=zeros)
            ENDIF q_with_l
         ENDDO
         ENDDO
      !
      ENDIF ultrasoft_or_paw
      !
      CALL iotk_scan_end(u, 'PP_AUGMENTATION')
      !
      ENDIF augmentation
      !
      ! Maximum radius of beta projector: outer radius to integrate
      upf%kkbeta = MAXVAL(upf%kbeta(1:upf%nbeta))
      ! For PAW augmentation charge may extend a bit further:
      IF(upf%tpawp) upf%kkbeta = MAX(upf%kkbeta, upf%paw%iraug)
      !
      CALL iotk_scan_end(u, 'PP_NONLOCAL')
      !
      RETURN
   END SUBROUTINE read_upf_nonlocal
   !
   SUBROUTINE read_upf_pswfc(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nw
      !
      CALL iotk_scan_begin(u, 'PP_PSWFC')
      !
      ALLOCATE( upf%chi(upf%mesh,upf%nwfc) )
      ALLOCATE( upf%els(upf%nwfc), &
                upf%oc(upf%nwfc), &
                upf%lchi(upf%nwfc), &
                upf%nchi(upf%nwfc), &
                upf%rcut_chi(upf%nwfc), &
                upf%rcutus_chi(upf%nwfc), &
                upf%epseu(upf%nwfc) &
              )
      !
      DO nw = 1,upf%nwfc
         CALL iotk_scan_dat(u, 'PP_CHI'//iotk_index(nw), &
                              upf%chi(:,nw), attr=attr)
         CALL iotk_scan_attr(attr, 'label',         upf%els(nw),     default='Xn')
         CALL iotk_scan_attr(attr, 'l',             upf%lchi(nw))
         CALL iotk_scan_attr(attr, 'occupation',    upf%oc(nw))
         CALL iotk_scan_attr(attr, 'n',             upf%nchi(nw),    default=upf%lchi(nw)-1)
         CALL iotk_scan_attr(attr, 'pseudo_energy', upf%epseu(nw),   default=0._dp)
         CALL iotk_scan_attr(attr, 'cutoff_radius', upf%rcut_chi(nw),default=0._dp)
         CALL iotk_scan_attr(attr, 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw),default=0._dp)
      ENDDO
      !
      CALL iotk_scan_end(u, 'PP_PSWFC')
      !
      RETURN
   END SUBROUTINE read_upf_pswfc

   SUBROUTINE read_upf_full_wfc(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      CHARACTER(len=iotk_attlenx) :: attr
      LOGICAL :: exst
      !
      INTEGER :: nb
      !
      IF(.not. upf%has_wfc) RETURN
      !
      CALL iotk_scan_begin(u, 'PP_FULL_WFC')
      !
      ALLOCATE( upf%aewfc(upf%mesh, upf%nbeta) )
      DO nb = 1,upf%nbeta
         CALL iotk_scan_dat(u, 'PP_AEWFC'//iotk_index(nb), &
                              upf%aewfc(:,nb), attr=attr)
      ENDDO

      IF (upf%has_so .and. upf%tpawp) THEN
         ALLOCATE( upf%paw%aewfc_rel(upf%mesh, upf%nbeta) )
nb_loop: DO nb = 1,upf%nbeta
            CALL iotk_scan_dat(u, 'PP_AEWFC_REL'//iotk_index(nb), &
                              upf%paw%aewfc_rel(:,nb), attr=attr, found=exst)
            IF (.not.exst) THEN
               upf%paw%aewfc_rel=0.0_DP
               EXIT nb_loop
            ENDIF
         ENDDO nb_loop
      ENDIF

      ALLOCATE( upf%pswfc(upf%mesh, upf%nbeta) )
      DO nb = 1,upf%nbeta
         CALL iotk_scan_dat(u, 'PP_PSWFC'//iotk_index(nb), &
                              upf%pswfc(:,nb), attr=attr)
      ENDDO
      CALL iotk_scan_end(u, 'PP_FULL_WFC')
      !
   END SUBROUTINE read_upf_full_wfc

   !
   SUBROUTINE read_upf_spin_orb(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nw, nb
      !
      IF (.not. upf%has_so) RETURN
      !
      CALL iotk_scan_begin(u, 'PP_SPIN_ORB')
      !
      ALLOCATE (upf%nn(upf%nwfc))
      ALLOCATE (upf%jchi(upf%nwfc))
      !
      DO nw = 1,upf%nwfc
         CALL iotk_scan_empty(u, 'PP_RELWFC'//iotk_index(nw),&
                               attr=attr)
            !CALL iotk_scan_attr(attr, 'els',   upf%els(nw))  ! already read
            CALL iotk_scan_attr(attr, 'nn',    upf%nn(nw))
            !CALL iotk_scan_attr(attr, 'lchi',  upf%lchi(nw)) ! already read
            CALL iotk_scan_attr(attr, 'jchi',  upf%jchi(nw))
            !CALL iotk_scan_attr(attr, 'oc',    upf%oc(nw))   ! already read
      ENDDO
      !
      ALLOCATE(upf%jjj(upf%nbeta))
      !
      DO nb = 1,upf%nbeta
         CALL iotk_scan_empty(u, 'PP_RELBETA'//iotk_index(nb),&
                               attr=attr)
            CALL iotk_scan_attr(attr, 'lll',   upf%lll(nb))
            CALL iotk_scan_attr(attr, 'jjj',   upf%jjj(nb))
      ENDDO
      !
      CALL iotk_scan_end(u, 'PP_SPIN_ORB')
      !
      RETURN
   END SUBROUTINE read_upf_spin_orb
   !
   SUBROUTINE read_upf_paw(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      !
      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nb,nb1

      IF (.not. upf%tpawp ) RETURN

      CALL iotk_scan_begin(u, 'PP_PAW', attr=attr)
      CALL iotk_scan_attr(attr, 'paw_data_format', upf%paw_data_format)
      IF(upf%paw_data_format /= 2) &
         CALL errore('read_upf_v2::paw',&
                     'Unknown format of PAW data.',1)
      CALL iotk_scan_attr(attr, 'core_energy', upf%paw%core_energy, default=0._dp)
      !
      ! Full occupation (not only > 0 ones)
      ALLOCATE( upf%paw%oc(upf%nbeta) )
      CALL iotk_scan_dat(u, 'PP_OCCUPATIONS',upf%paw%oc)
      !
      ! All-electron core charge
      ALLOCATE( upf%paw%ae_rho_atc(upf%mesh) )
      CALL iotk_scan_dat(u, 'PP_AE_NLCC', upf%paw%ae_rho_atc)
      !
      ! All-electron local potential
      ALLOCATE( upf%paw%ae_vloc(upf%mesh) )
      CALL iotk_scan_dat(u, 'PP_AE_VLOC', upf%paw%ae_vloc)
      !
      ALLOCATE(upf%paw%pfunc(upf%mesh, upf%nbeta,upf%nbeta) )
      upf%paw%pfunc(:,:,:) = 0._dp
      IF (upf%has_so) THEN
         ALLOCATE(upf%paw%pfunc_rel(upf%mesh, upf%nbeta,upf%nbeta) )
         upf%paw%pfunc_rel(:,:,:) = 0._dp
      ENDIF
      DO nb=1,upf%nbeta
         DO nb1=1,nb
            upf%paw%pfunc (1:upf%mesh, nb, nb1) = &
                  upf%aewfc(1:upf%mesh, nb) * upf%aewfc(1:upf%mesh, nb1)
            IF (upf%has_so) THEN
               upf%paw%pfunc_rel (1:upf%paw%iraug, nb, nb1) =  &
                        upf%paw%aewfc_rel(1:upf%paw%iraug, nb) *   &
                        upf%paw%aewfc_rel(1:upf%paw%iraug, nb1)
!
!    The small component is added to pfunc. pfunc_rel is useful only
!    to add a small magnetic contribution
!
               upf%paw%pfunc (1:upf%paw%iraug, nb, nb1) = &
                        upf%paw%pfunc (1:upf%paw%iraug, nb, nb1) + &
                        upf%paw%pfunc_rel (1:upf%paw%iraug, nb, nb1)  
            ENDIF
            upf%paw%pfunc(upf%paw%iraug+1:,nb,nb1) = 0._dp
            !
            upf%paw%pfunc (1:upf%mesh, nb1, nb) = upf%paw%pfunc (1:upf%mesh, nb, nb1)
            IF (upf%has_so) upf%paw%pfunc_rel (1:upf%mesh, nb1, nb) =  &
                                upf%paw%pfunc_rel (1:upf%mesh, nb, nb1)   
         ENDDO
      ENDDO
      !
      ! Pseudo wavefunctions (not only the ones for oc > 0)
      ! All-electron wavefunctions
      ALLOCATE(upf%paw%ptfunc(upf%mesh, upf%nbeta,upf%nbeta) )
      upf%paw%ptfunc(:,:,:) = 0._dp
      DO nb=1,upf%nbeta
         DO nb1=1,upf%nbeta
            upf%paw%ptfunc (1:upf%mesh, nb, nb1) = &
                  upf%pswfc(1:upf%mesh, nb) * upf%pswfc(1:upf%mesh, nb1)
            upf%paw%ptfunc(upf%paw%iraug+1:,nb,nb1) = 0._dp
            !
            upf%paw%ptfunc (1:upf%mesh, nb1, nb) = upf%paw%ptfunc (1:upf%mesh, nb, nb1)
         ENDDO
      ENDDO
      !
      ! Finalize
      CALL iotk_scan_end(u, 'PP_PAW')

      RETURN
   END SUBROUTINE read_upf_paw
!
   SUBROUTINE read_upf_gipaw(u, upf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      CHARACTER(len=iotk_attlenx) :: attr
      !
      INTEGER :: nb
      IF (.not. upf%has_gipaw ) RETURN
      
      CALL iotk_scan_begin(u, 'PP_GIPAW', attr=attr)
         CALL iotk_scan_attr(attr, 'gipaw_data_format', upf%gipaw_data_format)
      IF(upf%gipaw_data_format /= 2) &
         CALL infomsg('read_upf_v2::gipaw','Unknown format version')
      !
      CALL iotk_scan_begin(u, 'PP_GIPAW_CORE_ORBITALS', attr=attr)
         CALL iotk_scan_attr(attr, 'number_of_core_orbitals', upf%gipaw_ncore_orbitals)
      ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
      DO nb = 1,upf%gipaw_ncore_orbitals
         CALL iotk_scan_dat(u, 'PP_GIPAW_CORE_ORBITAL'//iotk_index(nb), &
                              upf%gipaw_core_orbital(:,nb), attr=attr)
            CALL iotk_scan_attr(attr, 'label', upf%gipaw_core_orbital_el(nb))
            CALL iotk_scan_attr(attr, 'n',     upf%gipaw_core_orbital_n(nb))
            CALL iotk_scan_attr(attr, 'l',     upf%gipaw_core_orbital_l(nb))
      ENDDO
      CALL iotk_scan_end(u, 'PP_GIPAW_CORE_ORBITALS')

      !
      ! Read valence all-electron and pseudo orbitals and their labels
      !
      IF (upf%paw_as_gipaw) THEN
         !READ PAW DATA INSTEAD OF GIPAW 
         upf%gipaw_wfs_nchannels = upf%nbeta
         ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
         DO nb = 1,upf%gipaw_wfs_nchannels
            upf%gipaw_wfs_el(nb) = upf%els_beta(nb)
            upf%gipaw_wfs_ll(nb) = upf%lll(nb)
            upf%gipaw_wfs_ae(:,nb) = upf%aewfc(:,nb)
         ENDDO
         DO nb = 1,upf%gipaw_wfs_nchannels
            upf%gipaw_wfs_ps(:,nb) = upf%pswfc(:,nb) 
         ENDDO
         ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
         ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
         upf%gipaw_vlocal_ae(:)= upf%paw%ae_vloc(:)  
         upf%gipaw_vlocal_ps(:)= upf%vloc(:)
         DO nb = 1,upf%gipaw_wfs_nchannels
            upf%gipaw_wfs_rcut(nb)=upf%rcut(nb)
            upf%gipaw_wfs_rcutus(nb)=upf%rcutus(nb)
         ENDDO
      ELSEIF (upf%tcoulombp) THEN
         upf%gipaw_wfs_nchannels = 1
         ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
         DO nb = 1,upf%gipaw_wfs_nchannels
            upf%gipaw_wfs_el(nb) = "1S"
            upf%gipaw_wfs_ll(nb) = 0
            upf%gipaw_wfs_ae(:,nb) = 0.0d0
            upf%gipaw_wfs_ps(:,nb) = 0.0d0
         ENDDO
         ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
         ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
         upf%gipaw_vlocal_ae(:)=  0.0d0
         upf%gipaw_vlocal_ps(:)=  0.0d0
         DO nb = 1,upf%gipaw_wfs_nchannels
            upf%gipaw_wfs_rcut(nb)=1.0d0
            upf%gipaw_wfs_rcutus(nb)=1.0d0
         ENDDO
      ELSE
         CALL iotk_scan_begin(u, 'PP_GIPAW_ORBITALS', attr=attr)
         CALL iotk_scan_attr(attr, 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels)
         ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
      
         DO nb = 1,upf%gipaw_wfs_nchannels
            CALL iotk_scan_begin(u, 'PP_GIPAW_ORBITAL'//iotk_index(nb), attr=attr)
            CALL iotk_scan_attr(attr, 'label', upf%gipaw_wfs_el(nb))
            CALL iotk_scan_attr(attr, 'l',     upf%gipaw_wfs_ll(nb))
            CALL iotk_scan_attr(attr, 'cutoff_radius',           upf%gipaw_wfs_rcut(nb))
            CALL iotk_scan_attr(attr, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb),&
                                                                 default=upf%gipaw_wfs_rcut(nb))
        ! read all-electron orbital
            CALL iotk_scan_dat(u, 'PP_GIPAW_WFS_AE', upf%gipaw_wfs_ae(:,nb))
        ! read pseudo orbital
            CALL iotk_scan_dat(u, 'PP_GIPAW_WFS_PS', upf%gipaw_wfs_ps(:,nb))
         !
            CALL iotk_scan_end(u, 'PP_GIPAW_ORBITAL'//iotk_index(nb))
         ENDDO
         CALL iotk_scan_end(u, 'PP_GIPAW_ORBITALS')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read all-electron and pseudo local potentials
         ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
         ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
         CALL iotk_scan_begin(u, 'PP_GIPAW_VLOCAL')
         CALL iotk_scan_dat(u, 'PP_GIPAW_VLOCAL_AE',upf%gipaw_vlocal_ae(:))
         CALL iotk_scan_dat(u, 'PP_GIPAW_VLOCAL_PS',upf%gipaw_vlocal_ps(:))
         CALL iotk_scan_end(u, 'PP_GIPAW_VLOCAL')

      ENDIF
      CALL iotk_scan_end(u, 'PP_GIPAW')
      RETURN
   END SUBROUTINE read_upf_gipaw
!
END SUBROUTINE read_upf_v2
!
END MODULE read_upf_v2_module
