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
        USE FoX_DOM
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
   TYPE(Node),POINTER,INTENT(IN)  :: u         ! pointer to root DOM Node 
   TYPE(pseudo_upf),INTENT(INOUT) :: upf       ! the pseudo data
   TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
   !
   INTEGER,OPTIONAL,INTENT(OUT):: ierr      ! /= 0 if something went wrong
   INTEGER :: ierr_
   TYPE(DOMException)   :: ex 
   TYPE(Node), POINTER  :: auxNode
   LOGICAL :: found
   LOGICAL,EXTERNAL :: matches
   CHARACTER(len = 256)  :: root
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
   root = getTagname(u, EX = ex)
   ierr_ = getExceptionCode(ex)  
   !
   IF((abs(ierr_)>0)  ) THEN
       !
       IF(.not. present(ierr)) &
         CALL errore('read_upf_v2','Cannot open UPF file.',1)
       ierr = 1
       RETURN
   ENDIF
   IF ( .not. matches('UPF',root) ) THEN
      IF (PRESENT (ierr) ) THEN 
         CALL infomsg( 'read_upf_v2', 'tagname is '//TRIM(root)//' instead of UPF' )
         ierr = 2 
         RETURN
      ELSE 
         CALL errore('read_upf_v2', 'tagname is '//TRIM(root)//' instead of UPF',2)
      END IF
   END IF
   CALL extractDataAttribute(u, 'version', upf%nv)
   IF (version_compare(upf%nv, max_version) == 'newer') &
       CALL errore('read_upf_v2',&
                   'Unknown UPF format version: '//TRIM(upf%nv),1)
   !
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
      auxNode => item(getElementsByTagname(u, 'PP_NLCC'), 0)
      CALL extractDataContent(auxNode, upf%rho_atc)
   ELSE
      ! A null core charge simplifies several functions, mostly in PAW
      upf%rho_atc(1:upf%mesh) = 0._dp
   ENDIF
   ! Read local potential
   IF(.not. upf%tcoulombp) THEN
      ALLOCATE( upf%vloc(upf%mesh) )
      auxNode => item( getElementsByTagname( u, 'PP_LOCAL'), 0)
      CALL extractDataContent(auxNode, upf%vloc)
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
   auxNode => item(getElementsByTagname(u, 'PP_RHOATOM'), 0) 
   CALL extractDataContent(auxNode, upf%rho_at)

   ! Read additional info for full-relativistic calculation
   CALL read_upf_spin_orb(u, upf)

   ! Read additional data for PAW (All-electron charge, wavefunctions, vloc..)
   CALL read_upf_paw(u, upf)

   ! Read data for gipaw reconstruction
   CALL read_upf_gipaw(u, upf)

   !
   ! Close the file (not the unit!)
   CALL destroy(u)
   !
   IF( present(ierr) ) ierr=0

   RETURN

   CONTAINS
   !
   SUBROUTINE read_upf_header(u, upf)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN)  :: u    ! parent node pointer
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER                     :: ierr, ios  ! /= 0 if something went wrong
      CHARACTER(len=256) :: dft_buffer  ! needed to allow the string defining the
                                        ! DFT flavor to be longer than upf%dft 
                                        ! (currntly 25) without getting iotk upset. 
                                        ! An error message is issued if trimmed 
                                        ! dft_buffer exceeds upf%dft size.
      INTEGER :: len_buffer
      !
      INTEGER :: nw
      TYPE(Node), POINTER  :: hdrNode
      CHARACTER(LEN=256)   :: attr
      TYPE(DOMException)   :: ex 
      !
      ! Read HEADER section with some initialization data
      hdrNode  => item( getElementsByTagname(u, 'PP_HEADER'), 0 )  
      IF ( hasAttribute( hdrNode, 'generated') ) THEN 
         CALL extractDataAttribute(hdrNode, 'generated', upf%generated) 
      ELSE 
        upf%generated = ' '
      END IF 
      IF ( hasAttribute( hdrNode, 'author') ) THEN 
         CALL extractDataAttribute(hdrNode, 'author', upf%author) 
      ELSE 
        upf%author = 'anonymous'
      END IF 
      IF ( hasAttribute( hdrNode, 'date') ) THEN 
         CALL extractDataAttribute(hdrNode, 'date', upf%date) 
      ELSE 
        upf%date = ' '
      END IF 
      IF ( hasAttribute( hdrNode, 'comment') ) THEN 
         CALL extractDataAttribute(hdrNode, 'comment', upf%comment) 
      ELSE 
        upf%comment = ' '
      END IF 
      !
      CALL extractDataAttribute(hdrNode, 'element', upf%psd)
      CALL extractDataAttribute(hdrNode, 'pseudo_type', upf%typ)
      
      CALL extractDataAttribute(hdrNode, 'relativistic', upf%rel)
      
      !
      CALL extractDataAttribute(hdrNode, 'is_ultrasoft', upf%tvanp, iostat = ios )
      IF ( ios /= 0 ) THEN 
          CALL extractDataAttribute(hdrNode, 'is_ultrasoft', attr) 
          upf%tvanp = ( INDEX (attr, 'T') > 0 )  
      END IF
      CALL extractDataAttribute(hdrNode, 'is_paw', upf%tpawp, iostat = ios)
      IF ( ios /= 0 ) THEN
          CALL extractDataAttribute(hdrNode, 'is_paw', attr)
          upf%tpawp = ( INDEX (attr, 'T') > 0 ) 
      END IF

      IF ( hasAttribute ( hdrNode, 'is_coulomb')) THEN
         CALL extractDataAttribute(hdrNode, 'is_coulomb', upf%tcoulombp, iostat = ios)
         IF ( ios /= 0 ) THEN 
            CALL extractDataAttribute ( hdrNode, 'is_coulomb', attr) 
            upf%tcoulombp = ( INDEX ( attr, 'T') > 0 ) 
         END IF
      ELSE 
         upf%tcoulombp = .FALSE.
      END IF
      !
      IF ( hasAttribute (hdrNode, 'has_so') ) THEN 
         CALL extractDataAttribute(hdrNode, 'has_so',         upf%has_so , IOSTAT = ios )
         IF ( ios /=0) THEN 
            CALL extractDataAttribute(hdrNode, 'has_so',  attr) 
            upf%has_so = ( INDEX ( attr, 'T') > 0 )
         END IF
      ELSE 
         upf%has_so = .false.
      END IF
      IF ( hasAttribute( hdrNode, 'has_wfc') ) THEN 
         CALL extractDataAttribute(hdrNode, 'has_wfc',        upf%has_wfc, IOSTAT = ios)
         IF ( ios /= 0 ) THEN 
            CALL extractDataAttribute (hdrNode, 'has_wfc', attr) 
            upf%has_wfc = ( INDEX(attr, 'T' ) > 0 ) 
         END IF 
      ELSE 
         upf%has_wfc = upf%tpawp
      END IF 
      IF ( hasAttribute ( hdrNode, 'has_gipaw' )) THEN  
         CALL extractDataAttribute(hdrNode, 'has_gipaw',      upf%has_gipaw, IOSTAT = ios )
         IF ( ios /= 0 ) THEN
           CALL extractDataAttribute(hdrNode, 'has_gipaw', attr ) 
           upf%has_gipaw = ( INDEX ( attr, 'T') > 0 ) 
         END IF 
      ELSE 
        upf%has_gipaw = .false.
      END IF 
      !EMINE
      IF ( hasAttribute ( hdrNode, 'paw_as_gipaw') ) THEN 
         CALL extractDataAttribute(hdrNode, 'paw_as_gipaw',      upf%paw_as_gipaw, IOSTAT = ios )
         IF ( ios /= 0 ) THEN 
            CALL extractDataAttribute(hdrNode, 'paw_as_gipaw', attr ) 
            upf%paw_as_gipaw = ( INDEX(attr, 'T') > 0 ) 
         END IF 
      ELSE 
        upf%paw_as_gipaw  = .false.
      END IF
      !
      CALL extractDataAttribute(hdrNode, 'core_correction',upf%nlcc, IOSTAT = ios)
      IF ( ios /= 0 ) THEN 
         CALL extractDataAttribute(hdrNode, 'core_correction', attr ) 
         upf%nlcc = ( INDEX( attr, 'T') > 0 ) 
      END IF
!        
      CALL extractDataAttribute(hdrNode, 'functional',  dft_buffer)
         len_buffer=len_trim(dft_buffer)
         if (len_buffer > len(upf%dft)) &
            call errore('read_upf_v2','String defining DFT is too long',len_buffer)
         upf%dft=TRIM(dft_buffer)

         CALL extractDataAttribute (hdrNode, 'z_valence',      upf%zp)
         IF ( hasAttribute (hdrNode,  'total_psenergy') ) THEN 
            CALL extractDataAttribute (hdrNode, 'total_psenergy', upf%etotps) 
         ELSE 
           upf%etotps = 0._dp
         END IF 
         IF ( hasAttribute (hdrNode, 'wfc_cutoff') ) THEN 
            CALL extractDataAttribute (hdrNode, 'wfc_cutoff',     upf%ecutwfc )
         ELSE 
           upf%ecutwfc = 0._dp 
         END IF
         IF  ( hasAttribute (hdrNode, 'rho_cutoff') ) THEN 
            CALL extractDataAttribute (hdrNode, 'rho_cutoff',     upf%ecutrho) 
         ELSE 
           upf%ecutrho = 0._dp 
         END IF
         IF ( hasAttribute ( hdrNode, 'l_max' ) ) THEN 
            CALL extractDataAttribute (hdrNode, 'l_max',          upf%lmax) 
         ELSE      
           upf%lmax =0 
         END IF
         IF ( hasAttribute ( hdrNode, 'l_max_rho') ) THEN 
            CALL extractDataAttribute (hdrNode, 'l_max_rho',      upf%lmax_rho) 
         ELSE 
           upf%lmax_rho = 2*upf%lmax
         END IF
         IF ( hasAttribute ( hdrNode, 'l_local') ) THEN 
            CALL extractDataAttribute (hdrNode, 'l_local',        upf%lloc ) 
         ELSE 
            upf%lloc = 0
         END IF
         CALL extractDataAttribute (hdrNode, 'mesh_size',      upf%mesh)
         CALL extractDataAttribute (hdrNode, 'number_of_wfc',  upf%nwfc)
         CALL extractDataAttribute (hdrNode, 'number_of_proj', upf%nbeta)
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
      TYPE (Node),POINTER,INTENT(IN)       :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT)         :: upf  ! the pseudo data
      TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
      !
      INTEGER                     :: ierr ! /= 0 if something went wrong
      TYPE (Node),POINTER         :: mshNode, locNode
      !
      LOGICAL :: found
      !
      mshNode => item( getElementsByTagname(u, 'PP_MESH'),0 )
      IF ( hasAttribute(mshNode, 'dx')) THEN 
         CALL extractDataAttribute(mshNode, 'dx',   upf%dx ) 
      ELSE 
        upf%dx  = 0._dp
      END IF
      IF ( hasAttribute (mshNode, 'mesh')) &
             CALL extractDataAttribute(mshNode, 'mesh', upf%mesh )
      IF ( hasAttribute ( mshNode, 'mesh') ) THEN 
         CALL extractDataAttribute(mshNode, 'xmin', upf%xmin ) 
      ELSE 
          upf%xmin = 0._dp
      END IF
      IF ( hasAttribute ( mshNode, 'rmax') ) THEN
          CALL extractDataAttribute(mshNode, 'rmax', upf%rmax )
      ELSE
          upf%rmax = 0._dp 
      END IF
      IF ( hasAttribute ( mshNode, 'zmesh') ) THEN
          CALL extractDataAttribute(mshNode, 'zmesh',upf%zmesh ) 
      ELSE 
          upf%zmesh = 0._dp 
      END IF
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
      locNode => item( getElementsByTagname( mshNode, 'PP_R'), 0 ) 
      CALL extractDataContent(locNode, upf%r(1:upf%mesh))  
      !
      locNode => item(getElementsByTagname( mshNode, 'PP_RAB'), 0)
      CALL extractDataContent(locNode, upf%rab(1:upf%mesh)) 
      !
      IF (present(grid)) THEN
         ! Reconstruct additional grids
         upf%grid%r2 =  upf%r**2
         upf%grid%sqr = sqrt(upf%r)
         upf%grid%rm1 = upf%r**(-1)
         upf%grid%rm2 = upf%r**(-2)
         upf%grid%rm3 = upf%r**(-3)
      ENDIF
      !
      RETURN
   END SUBROUTINE read_upf_mesh
   !
   SUBROUTINE read_upf_nonlocal(u, upf)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN)         :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      !
      !
      TYPE (Node),POINTER     :: nlcNode, locNode,locNode1, locNode2, locNode3
      TYPE (nodeList),POINTER :: tmpList 
      INTEGER :: nb,mb,ln,lm,l,nmb,ierr=0
      !INTEGER :: nb_=-1,mb_=-1,l_=-1,nmb_=-1
      REAL(DP):: zeros(upf%mesh)
      REAL(DP), ALLOCATABLE :: tmp_dbuffer(:)
      LOGICAL :: isnull, found
      CHARACTER(LEN = 256 )   :: attr 
      INTEGER :: ios 
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
         RETURN
      END IF
      !
      ! <AF>
      
      nlcNode => item ( getElementsByTagname( u, 'PP_NONLOCAL'), 0 )  
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
      locNode2 => getFirstChild(nlcNode) 
      nb = 0 
      DO 
         IF (.NOT. ASSOCIATED( locNode2) )  EXIT
         locNode => locNode2
         locNode2 => getNextSibling(locNode) 
         IF (getNodeType( locNode) .NE. ELEMENT_NODE ) CYCLE
         IF ( INDEX(getTagName(locNode), 'PP_BETA') .LE. 0 ) CYCLE
         nb = nb + 1
         CALL extractDataContent(locNode, upf%beta(:, nb)) 
         IF ( hasAttribute( locNode, 'label') ) THEN 
            CALL extractDataAttribute(locNode, 'label',                  upf%els_beta(nb))
         ELSE 
            upf%els_beta(nb) ='Xn'
         END IF
         CALL extractDataAttribute(locNode, 'angular_momentum',       upf%lll(nb))
         IF ( hasAttribute( locNode,'cutoff_radius_index' )) THEN 
            CALL extractDataAttribute(locNode, 'cutoff_radius_index',    upf%kbeta(nb)) 
         ELSE    
           upf%kbeta = upf%mesh
         END IF 
         IF ( hasAttribute( locNode,'cutoff_radius' )) THEN 
            CALL extractDataAttribute(locNode, 'cutoff_radius',          upf%rcut(nb) ) 
         ELSE  
            upf%rcut(nb) = 0._dp
         END IF
         IF ( hasAttribute( locNode,'ultrasoft_cutoff_radius' )) THEN
            CALL extractDataAttribute(locNode, 'ultrasoft_cutoff_radius', upf%rcutus(nb)) 
         ELSE    
           upf%rcutus(nb)  = 0._dp
         END IF
!
!    Old version of UPF PPs v.2 contained an error in the tag. 
!    To be able to read the old PPs we need the following
!
         IF ( upf%rcutus(nb)==0._DP) THEN 
            IF ( hasAttribute( locNode,'norm_conserving_radius' )) THEN
               CALL extractDataAttribute(locNode,'norm_conserving_radius',upf%rcutus(nb))  
            ELSE
               upf%rcutus(nb)  = 0._dp
            END IF
         END IF  
      ENDDO 
      !
      ! Read the hamiltonian terms D_ij
      locNode => item( getElementsByTagname(nlcNode, 'PP_DIJ'),0)    
      CALL extractDataContent(locNode, upf%dion)
      !   CALL iotk_scan_attr(attr, 'non_zero_elements', upf%nd)
      !
      ! Read the augmentation charge section
      augmentation : &
      IF(upf%tvanp .or. upf%tpawp) THEN
      !
      locNode => item(getElementsByTagname(nlcNode, 'PP_AUGMENTATION'),0) 
         CALL extractDataAttribute(locNode, 'q_with_l', upf%q_with_l, IOSTAT = ios )
         IF ( ios /= 0) THEN 
            CALL extractDataAttribute(locNode, 'q_with_l', attr )
            upf%q_with_l = ( INDEX ( attr, 'T') > 0)
         END IF
         CALL extractDataAttribute(locNode, 'nqf',      upf%nqf)
         IF (hasAttribute(locNode, 'nqlc') ) THEN 
            CALL extractDataAttribute(locNode, 'nqlc',     upf%nqlc) 
         ELSE  
            upf%nqlc =2*upf%lmax+1
         END IF
         IF (upf%tpawp) THEN
            IF (hasAttribute(locNode, 'shape') ) THEN 
               CALL extractDataAttribute(locNode,'shape',          upf%paw%augshape) 
            ELSE 
               upf%paw%augshape ='UNKNOWN'
            END IF
            IF (hasAttribute(locNode, 'cutoff_r') ) THEN
               CALL extractDataAttribute(locNode,'cutoff_r',       upf%paw%raug ) 
            ELSE 
               upf%paw%raug = 0._dp
            END IF
            IF (hasAttribute(locNode, 'cutoff_r_index') ) THEN
               CALL extractDataAttribute(locNode,'cutoff_r_index', upf%paw%iraug)  
            ELSE 
              upf%paw%iraug =upf%mesh 
            END IF
            IF (hasAttribute(locNode, 'l_max_aug') ) THEN
               CALL extractDataAttribute(locNode,'l_max_aug',      upf%paw%lmax_aug) 
            ELSE 
              upf%paw%lmax_aug   =upf%lmax_rho
            END IF
         ENDIF
         ! a negative number means that all qfunc are stored
         IF (hasAttribute(locNode,   'augmentation_epsilon'  ) ) THEN
            CALL extractDataAttribute(locNode,'augmentation_epsilon',upf%qqq_eps) 
         ELSE 
            upf%qqq_eps = -1._dp
         END IF
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
      locNode2 => item( getElementsByTagname( locNode, 'PP_Q'), 0) 
      CALL extractDataContent(locNode2, upf%qqq )
      !
      ! read charge multipoles (only if PAW)
      IF( upf%tpawp ) THEN   
         ALLOCATE(upf%paw%augmom(upf%nbeta,upf%nbeta, 0:2*upf%lmax))
         ALLOCATE( tmp_dbuffer(upf%nbeta*upf%nbeta*(2*upf%lmax+1)) )
         locNode2 => item( getElementsByTagname(locNode,'PP_MULTIPOLES'), 0)
         CALL extractDataContent(locNode2, tmp_dbuffer)
         upf%paw%augmom=reshape(tmp_dbuffer, [upf%nbeta,upf%nbeta,2*upf%lmax+1])
         DEALLOCATE (tmp_dbuffer)
      ENDIF
      !
      ! Read polinomial coefficients for Q_ij expansion at small radius
      IF(upf%nqf <= 0) THEN
         upf%rinner(:) = 0._dp
         ALLOCATE( upf%qfcoef(1,1,1,1) )
         upf%qfcoef = 0._dp
      ELSE
         ALLOCATE( upf%qfcoef( MAX( upf%nqf,1 ), upf%nqlc, upf%nbeta, upf%nbeta ) )
         ALLOCATE(tmp_dbuffer(MAX( upf%nqf,1 )*upf%nqlc*upf%nbeta*upf%nbeta))
         locNode2=> item(getElementsByTagname(locNode, 'PP_QFCOEFF'),0) 
         CALL extractDataContent(locNode2, tmp_dbuffer)
         upf%qfcoef = reshape(tmp_dbuffer,[size(upf%qfcoef,1),size(upf%qfcoef,2),&
                                           size(upf%qfcoef,3),size(upf%qfcoef,4)])
         DEALLOCATE(tmp_dbuffer)
         locNode2 => item(getElementsByTagname(locNode, 'PP_RINNER'),0)
         CALL extractDataContent(locNode2, upf%rinner)
      ENDIF
      !
      ! Read augmentation charge Q_ij
      ultrasoft_or_paw : &
      IF( upf%tvanp) THEN
         locNode3 => getFirstChild(locNode)
         IF (upf%q_with_l) THEN 
            upf%qfuncl = 0._dp
         ELSE 
            upf%qfunc = 0._dp
         END IF
         search_for_qij: DO 
           IF ( .NOT. ASSOCIATED(locNode3) ) EXIT search_for_qij
           locNode2 => locNode3
           locNode3 => getNextSibling(locNode2)
           IF (getNodeType(locNode2) .NE. ELEMENT_NODE) CYCLE search_for_qij
           !
           IF ( INDEX( getTagName(locNode2), 'PP_QIJ') .LE. 0) CYCLE search_for_qij
           CALL extractDataAttribute(locNode2, 'composite_index', nmb)
           IF (upf%q_with_l) THEN 
              CALL extractDataAttribute(locNode2, 'angular_momentum', l)
              CALL extractDataContent( locNode2, upf%qfuncl(:, nmb,l))
              IF (upf%tpawp) upf%qfuncl(upf%paw%iraug+1:,nmb,l) = 0._DP
           ELSE
              CALL extractDataContent ( locNode2, upf%qfunc(:,nmb))
           END IF
         END DO search_for_qij     
      !
      ENDIF ultrasoft_or_paw
      !
      !
      ENDIF augmentation
      !
      ! Maximum radius of beta projector: outer radius to integrate
      upf%kkbeta = MAXVAL(upf%kbeta(1:upf%nbeta))
      ! For PAW augmentation charge may extend a bit further:
      IF(upf%tpawp) upf%kkbeta = MAX(upf%kkbeta, upf%paw%iraug)
      !
      !
      RETURN
   END SUBROUTINE read_upf_nonlocal
   !
   SUBROUTINE read_upf_pswfc(u, upf)
      IMPLICIT NONE
      TYPE(Node), POINTER, INTENT(IN) :: u    ! pointer to root node
      TYPE(pseudo_upf),INTENT(INOUT)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      !
      INTEGER :: nw
      TYPE(Node),POINTER                :: pswfcNode, locNode, locNode2
      !
      pswfcNode  => item(getElementsByTagname(u, 'PP_PSWFC'), 0) 
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
      locNode2 => getFirstChild ( pswfcNode ) 
      nw = 0 
      DO  
         IF (.NOT. ASSOCIATED( locNode2) ) EXIT
         locNode => locNode2 
         locNode2 => getNextSibling(locNode) 
         IF (getNodeType(locNode) .NE. ELEMENT_NODE ) CYCLE 
         IF ( INDEX ( getTagName(locNode),'PP_CHI') .LE. 0 ) CYCLE
         nw = nw + 1 
         IF ( nw .GT. upf%nwfc ) THEN 
            CALL infomsg('pseudo '//trim(upf%psd), "too many chi found in pswfc" ) 
            EXIT
         END IF
         IF ( hasAttribute (locNode, 'label')) THEN 
            CALL extractDataAttribute(locNode, 'label', upf%els(nw) )
         ELSE 
            upf%els(nw) = 'Xn'
         END IF 
         CALL extractDataAttribute(locNode, 'l', upf%lchi(nw)) 
         CALL extractDataAttribute(locNode, 'occupation',    upf%oc(nw))
         IF ( hasAttribute(locNode, 'n')) THEN 
            CALL extractDataAttribute(locNode, 'n',             upf%nchi(nw)) 
         ELSE 
             upf%nchi = upf%lchi(nw)-1 
         END IF
         IF ( hasAttribute(locNode, 'pseudo_energy') ) THEN 
            CALL extractDataAttribute(locNode, 'pseudo_energy', upf%epseu(nw) ) 
         ELSE 
            upf%epseu(nw) = 0._dp
         END IF
         IF ( hasAttribute( locNode,'cutoff_radius') ) THEN 
            CALL extractDataAttribute(locNode, 'cutoff_radius', upf%rcut_chi(nw) ) 
         ELSE 
            upf%rcut_chi(nw) =0._dp 
         END IF
         IF ( hasAttribute(locNode, 'ultrasoft_cutoff_radius') ) THEN 
            CALL extractDataAttribute(locNode, 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw)) 
         ELSE 
           upf%rcutus_chi(nw) =0._dp 
         END IF 
         CALL extractDataContent(locNode, upf%chi(:,nw) )
      ENDDO
      !
      RETURN
   END SUBROUTINE read_upf_pswfc

   SUBROUTINE read_upf_full_wfc(u, upf)
      IMPLICIT NONE
      TYPE(Node),POINTER, INTENT(IN) :: u    ! parent node
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      LOGICAL :: exst
      !
      INTEGER :: nbae, nbae_rel, nbps
      TYPE(Node),POINTER       :: fllwfNode, locNode, locNode2
      !
      IF(.not. upf%has_wfc) RETURN
      !
      fllwfNode => item(getElementsByTagname (u, 'PP_FULL_WFC'), 0)
      !
      ALLOCATE( upf%aewfc(upf%mesh, upf%nbeta) )
      ALLOCATE( upf%pswfc(upf%mesh, upf%nbeta) )
      IF (upf%has_so .and. upf%tpawp) THEN 
        ALLOCATE( upf%paw%aewfc_rel(upf%mesh, upf%nbeta) )
        upf%paw%aewfc_rel = 0._dp
      END IF 
      locNode2 => getFirstChild(fllwfNode) 
      nbae=0 
      nbae_rel = 0 
      nbps = 0 
      DO 
         IF ( .NOT. ASSOCIATED ( locNode2) ) EXIT 
         locNode => locNode2
         locNode2=> getNextSibling(locNode) 
         IF (getNodeType(locNode)  .NE. ELEMENT_NODE) CYCLE 
         IF (INDEX(getTagName(locNode), 'PP_AEWFC_REL') .GT. 0 ) THEN 
            nbae_rel = nbae_rel+1
            IF (nbae_rel .GT. upf%nbeta ) THEN 
               CYCLE
            ELSE
              CALL extractDataContent(locNode, upf%paw%aewfc_rel(:,nbae_rel)) 
            END IF 
         ELSE IF (INDEX(getTagName(locNode),'PP_AEWFC') .GT. 0) THEN 
            nbae  = nbae +1 
            IF (nbae .GT. upf%nbeta) THEN 
               CYCLE
            ELSE 
               CALL extractDataContent(locNode, upf%aewfc(:,nbae)) 
            END IF 
         ELSE IF (INDEX(getTagName(locNode), 'PP_PSWFC') .GT. 0) THEN 
            nbps = nbps + 1 
            IF ( nbps .LE. upf%nbeta ) THEN 
               CALL extractDataContent(locNode, upf%pswfc(:, nbps) )
            END IF 
         END IF 
      ENDDO
   END SUBROUTINE read_upf_full_wfc

   !
   SUBROUTINE read_upf_spin_orb(u, upf)
      IMPLICIT NONE
      TYPE(Node), POINTER, INTENT(IN) :: u    ! parent node pointer
      TYPE(pseudo_upf),INTENT(INOUT)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      !
      !
      INTEGER :: nw, nb
      TYPE(Node), POINTER  :: soNode, locNode, locNode2
      !
      IF (.not. upf%has_so) RETURN
      !
      soNode => item(getElementsByTagName(u, 'PP_SPIN_ORB'), 0) 
      !
      ALLOCATE (upf%nn(upf%nwfc))
      ALLOCATE (upf%jchi(upf%nwfc))
      !
      ALLOCATE(upf%jjj(upf%nbeta))
      !
      locNode2=> getFirstChild(soNode)
      nw = 0 
      nb = 0
      DO 
        IF (.NOT. ASSOCIATED(locNode2)) EXIT
           locNode => locNode2
           locNode2 => getNextSibling(locNode) 
           IF ( getNodeType(locNode) .NE. ELEMENT_NODE ) CYCLE
           select_tag: IF ( INDEX(getTagName(locNode),'PP_RELWFC') .GT. 0) THEN 
              nw = nw + 1 
              IF (nw .LE. upf%nwfc ) THEN 
                 CALL extractDataAttribute (locNode, 'nn', upf%nn(nw))   
                 CALL extractDataAttribute (locNode, 'jchi',  upf%jchi(nw))
                 !extraxtDataAttribute(attr, 'els',   upf%els(nw))  ! already read 
                 !extraxtDataAttribute(attr, 'lchi',  upf%lchi(nw)) ! already read
                 !extraxtDataAttribute(attr, 'oc',    upf%oc(nw))   ! already read
              END IF
           ELSE IF (INDEX(getTagName(locNode),'PP_RELBETA') .GT. 0) THEN
              nb = nb + 1 
              IF (nb .LE. upf%nbeta ) THEN 
                 CALL extractDataAttribute(locNode, 'lll',   upf%lll(nb))
                 CALL extractDataAttribute(locNode, 'jjj',   upf%jjj(nb))
              END IF 
           END IF select_tag
      ENDDO
      !
      RETURN
   END SUBROUTINE read_upf_spin_orb
   !
   SUBROUTINE read_upf_paw(u, upf)
      IMPLICIT NONE
      TYPE(Node), POINTER, INTENT(IN) :: u    ! 
      TYPE(pseudo_upf),INTENT(INOUT)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      TYPE(Node), POINTER      :: pawNode, locNode
      INTEGER :: nb,nb1

      IF (.not. upf%tpawp ) RETURN

      pawNode => item(getElementsByTagname(u, 'PP_PAW'), 0)
      CALL extractDataAttribute(pawNode, 'paw_data_format', upf%paw_data_format)
      IF(upf%paw_data_format /= 2) &
         CALL errore('read_upf_v2::paw',&
                     'Unknown format of PAW data.',1)
      IF (hasAttribute(pawNode, 'core_energy')) THEN
         CALL extractDataAttribute(pawNode, 'core_energy', upf%paw%core_energy) 
      ELSE 
        upf%paw%core_energy = 0._dp 
      END IF
      !
      ! Full occupation (not only > 0 ones)
      ALLOCATE( upf%paw%oc(upf%nbeta) )
      locNode => item(getElementsByTagname(pawNode, 'PP_OCCUPATIONS'), 0)
      CALL  extractDataContent(locNode, upf%paw%oc)
      
      !
      ! All-electron core charge
      ALLOCATE( upf%paw%ae_rho_atc(upf%mesh) )
      locNode => item(getElementsByTagname(pawNode, 'PP_AE_NLCC'), 0)
      CALL extractDataContent(locNode, upf%paw%ae_rho_atc)
      
      !
      ! All-electron local potential
      ALLOCATE( upf%paw%ae_vloc(upf%mesh) )
      locNode => item(getElementsByTagname( pawNode, 'PP_AE_VLOC'), 0)
      CALL extractDataContent(locNode, upf%paw%ae_vloc)

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
      RETURN
   END SUBROUTINE read_upf_paw
!
   SUBROUTINE read_upf_gipaw(u, upf)
      IMPLICIT NONE
      TYPE(Node), POINTER, INTENT(IN)  :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT)   :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      INTEGER :: nb
      TYPE(Node), POINTER   :: gpawNode, locNode, locNode2, locNode3, locNode4
      TYPE(NodeList), POINTER  :: locList
      IF (.not. upf%has_gipaw ) RETURN
      !
      gpawNode => item(getElementsByTagname(u, 'PP_GIPAW'), 0)
         CALL extractDataAttribute(gpawNode, 'gipaw_data_format', upf%gipaw_data_format)
      IF(upf%gipaw_data_format /= 2) &
         CALL infomsg('read_upf_v2::gipaw','Unknown format version')
      !
      locNode => item(getElementsByTagname(gpawNode, 'PP_GIPAW_CORE_ORBITALS'), 0) 
         CALL extractDataAttribute(locNode, 'number_of_core_orbitals', upf%gipaw_ncore_orbitals)
      ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
      locNode2  => getFirstChild( locNode )  
      nb = 0 
      DO 
        IF ( .NOT. ASSOCIATED(locNode2) ) EXIT
        IF ( getNodeType( locNode2 )  == ELEMENT_NODE ) THEN 
           IF ( INDEX(getTagName(locNode2), 'PP_GIPAW_CORE_ORBITAL' ) > 0) THEN 
              nb = nb + 1 
              CALL extractDataContent(locNode2, upf%gipaw_core_orbital(:,nb) )
              CALL extractDataAttribute(locNode2, 'label', upf%gipaw_core_orbital_el(nb))
              CALL extractDataAttribute(locNode2, 'n',     upf%gipaw_core_orbital_n(nb))
              CALL extractDataAttribute(locNode2, 'l',     upf%gipaw_core_orbital_l(nb))
           END IF 
        END IF 
        locNode3 => locNode2 
        locNode2 => getNextSibling(locNode3) 
      ENDDO
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
         locNode => item(getElementsByTagname(gpawNode, 'PP_GIPAW_ORBITALS'), 0)
         CALL extractDataAttribute(locNode, 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels)

         ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )

         locNode2 => getFirstChild(locNode) 
         nb = 0  
         DO 
            IF (.NOT. ASSOCIATED( locNode2 ))  EXIT
            !
            IF ( getNodeType( locNode2 ) == ELEMENT_NODE)  THEN 
               IF ( index( getTagName( locNode2 ), "PP_GIPAW_ORBITAL." )> 0) THEN  
                  nb = nb + 1 
                  CALL extractDataAttribute(locNode2, 'label', upf%gipaw_wfs_el(nb))
                  CALL extractDataAttribute(locNode2, 'l',     upf%gipaw_wfs_ll(nb))
                  CALL extractDataAttribute(locNode2, 'cutoff_radius',           upf%gipaw_wfs_rcut(nb))
                  IF (hasAttribute(locNode, 'ultrasoft_cutoff_radius') ) THEN
                     CALL extractDataAttribute(locNode, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb) )
                  ELSE
                     upf%gipaw_wfs_rcutus(nb) = upf%gipaw_wfs_rcut(nb) 
                  END IF
                  ! read all-electron orbital
                  locNode4 => item( getElementsByTagname(locNode2, 'PP_GIPAW_WFS_AE'), 0)
                  CALL extractDataContent(locNode4, upf%gipaw_wfs_ae(:,nb))
                  ! read pseudo orbital
                  locNode4 => item( getElementsByTagname( locNode2, 'PP_GIPAW_WFS_PS'), 0)
                  CALL extractDataContent(locNode4, upf%gipaw_wfs_ps(:,nb))
               END IF 
           END IF 
           locNode3 => locNode2 
           locNode2 => getNextSibling(locNode3) 
         !
         ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read all-electron and pseudo local potentials
         ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
         ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
         locNode => item(getElementsByTagname( gpawNode, 'PP_GIPAW_VLOCAL'), 0)
           !
           locNode2 => item(getElementsBytagname( locNode, 'PP_GIPAW_VLOCAL_AE'), 0)
              CALL extractDataContent (locNode2, upf%gipaw_vlocal_ae(:)) 
           !
           locNode2 => item(getElementsByTagname( locNode, 'PP_GIPAW_VLOCAL_PS'), 0)
              CALL extractDataContent(locNode2, upf%gipaw_vlocal_ps(:))
      ENDIF
      RETURN
   END SUBROUTINE read_upf_gipaw
!
END SUBROUTINE read_upf_v2
!
END MODULE read_upf_v2_module
