!
! Copyright (C) 2008-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
      MODULE read_upf_schema_module
!=----------------------------------------------------------------------------=!
!  this module handles the reading of pseudopotential data

! ...   declare modules
        USE kinds,        ONLY: DP
        USE pseudo_types, ONLY: pseudo_upf
        USE radial_grids, ONLY: radial_grid_type
        USE parser,       ONLY : version_compare
        USE FoX_Dom
        !
        PRIVATE
        PUBLIC :: read_upf_schema

INTERFACE searchData
    MODULE PROCEDURE searchStringData, searchBooleanData, searchRealData, searchIntegerData
END INTERFACE searchData
 CONTAINS

!------------------------------------------------+
SUBROUTINE read_upf_schema(pseudo, upf, grid, ierr )             !
   !---------------------------------------------+
   ! Read pseudopotential in UPF schema, uses FoX libs
   !
   USE pseudo_types, ONLY: nullify_pseudo_upf, deallocate_pseudo_upf
   USE radial_grids, ONLY: radial_grid_type, nullify_radial_grid
   USE FoX_dom  
   IMPLICIT NONE
   TYPE(Node), POINTER, INTENT(IN)   :: pseudo   ! pointer to root node
   TYPE(pseudo_upf),INTENT(INOUT)    :: upf       ! the pseudo data
   TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
   !
   INTEGER,OPTIONAL,INTENT(OUT):: ierr      ! /= 0 if something went wrong
   type(Node), pointer :: root, header, info, pNode
   type(NodeList), pointer :: pointList
   INTEGER :: ierr_
   LOGICAL :: found
   LOGICAL,EXTERNAL :: matches
   CHARACTER(len=6),PARAMETER :: max_version = '0.1'
   CHARACTER(LEN=80)          :: fileRef
   INTEGER                    :: nw 
   !
   ! Initialize the file
   !
   ! header and info elements, check version extract main info  
   header => item ( getElementsByTagname(pseudo, "pp_header" ),0)
   info   => item ( getElementsByTagname(pseudo, "pp_info" ), 0)  
   call searchData ( 'xsd_version', upf%nv, pseudo)   
   IF (version_compare(upf%nv, max_version) == 'newer') &
       CALL errore('read_upf_schema',&
                   'Unknown UPF format version: '//TRIM(upf%nv),1)
   CALL read_upf_header( header, info, upf)
   pointList => getElementsByTagname(info, "valence_orbital") 
   !ALLOCATE(upf%els(upf%nwfc), upf%nchi(upf%nwfc), upf%lchi(upf%nwfc), upf%oc(upf%nwfc), &
   !         upf%rcut_chi(upf%nwfc), upf%rcutus_chi(upf%nwfc), upf%epseu(upf%nwfc))
   !DO nw = 1, upf%nwfc
   !   pNode => item( pointList, nw -1 ) 
   !   CALL extractDataAttribute(pNode, "nl", upf%els(nw))
   !   CALL extractDataAttribute(pNode, "pn", upf%nchi(nw))
   !   CALL extractDataAttribute(pNode, "l", upf%lchi(nw)) 
   !   CALL searchData ('occupation', upf%oc(nw), pNode) 
   !   CALL searchData('Rcut', upf%rcut_chi(nw), pNode) 
   !   CALL searchData('RcutUS', upf%rcutus_chi(nw), pNode)
   !   CALL searchData('Epseu', upf%epseu(nw), pNode )
   !END DO    
   !! no need to read these data here these field will be filled later in read_pswfc 
   fileRef = TRIM(upf%psd)//'-'//TRIM(upf%typ)
   IF(upf%tpawp .and. .not. present(grid)) &
      CALL errore('read_upf_v2', 'PAW requires a radial_grid_type.', 1)
   !
   ! Read radial grid mesh
   pNode => item ( getElementsByTagname (pseudo, "pp_mesh"), 0 ) 
   CALL read_upf_mesh(pNode, upf, grid)
   ! Read non-linear core correction charge
   ALLOCATE( upf%rho_atc(upf%mesh) )
   IF(upf%nlcc) THEN
      pNode => item(getElementsByTagname( pseudo, 'pp_nlcc' ), 0) 
      CALL extractDataContent (pNode, upf%rho_atc)
   ELSE
      ! A null core charge simplifies several functions, mostly in PAW
      upf%rho_atc(1:upf%mesh) = 0._dp
   ENDIF

   ! Read local potential
   IF(.not. upf%tcoulombp) THEN
      ALLOCATE( upf%vloc(upf%mesh) )
      pNode => item ( getElementsByTagname ( pseudo, 'pp_local'),0)
      CALL extractDataContent  (pNode, upf%vloc)
   ENDIF
   ! Read nonlocal components: projectors, augmentation, hamiltonian elements
   pNode => item ( getElementsByTagname ( pseudo, 'pp_nonlocal'),0)
   CALL read_upf_nonlocal( pNode, upf)
   ! Read initial pseudo wavefunctions
   ! (usually only wfcs with occupancy > 0)
   pNode => item(getElementsByTagname (pseudo, 'pp_pswfc'),0)
   ! Close the file (not the unit!)
   IF ( .NOT. ASSOCIATED ( pNode))  CALL errore ( 'read_upf_schema', 'pp_pswfc not found in '//TRIM(fileRef), 5)
   CALL read_upf_pswfc(pNode, upf)
   ! Read all-electron and pseudo wavefunctions
   pNode => item (getElementsByTagname (pseudo, 'pp_full_wfc'),0) 
   IF (ASSOCIATED ( pNode ) )  CALL read_upf_full_wfc( pNode, upf)
   ! Read valence atomic density (used for initial density)
   ALLOCATE( upf%rho_at(upf%mesh) )
   pNode => item(getElementsByTagname(pseudo, 'pp_rhoatom'), 0) 
   IF (.NOT. ASSOCIATED( pNode))  CALL errore ('read_upf_schema', 'pp_rhoatom not found in '//TRIM(fileRef), 6 )
   CALL extractDataContent (pNode, upf%rho_at, IOSTAT = ierr_ ) 
   IF ( ierr_ /= 0 ) &
         CALL errore ( 'read_upf_schema', 'error reading rho_atom not in '//TRIM(fileRef), ierr_ )  
   ! Read additional info for full-relativistic calculation
   IF ( upf%has_so ) THEN
      pNode => item ( getElementsByTagname (pseudo, 'pp_spinorb'), 0 ) 
      IF ( .NOT. ASSOCIATED(pNode) ) &
         CALL errore ('read_upf_schema', 'pp_spinorb not found in '//TRIM(fileRef), 7 )   
      CALL read_upf_spin_orb(pNode, upf)
   END IF 
   ! Read additional data for PAW (All-electron charge, wavefunctions, vloc..)
   IF (upf%tpawp ) THEN 
      pNode => item ( getElementsByTagname (pseudo, 'pp_paw'),0)
      IF (.NOT. ASSOCIATED(pNode) ) CALL errore ('read_upf_schema', 'pp_paw not found in '//TRIM(fileRef),8)
      CALL read_upf_paw(pNode, upf)
   END IF
   ! Read data for gipaw reconstruction
   IF (upf%has_gipaw) THEN 
      pNode => item( getElementsByTagname(pseudo, "pp_gipaw"), 0) 
      IF (.NOT. ASSOCIATED ( pNode) ) CALL errore ('read_upf_schema', 'pp_gipaw not found in '//TRIM(fileRef),9)
      CALL read_upf_gipaw(pNode, upf)
   END IF
   !
   !
   IF( present(ierr) ) ierr=0

   RETURN

   END SUBROUTINE read_upf_schema
   !
   SUBROUTINE read_upf_header(header, info,  upf)
      IMPLICIT NONE
      TYPE(Node),POINTER             :: header,&    ! XML pp_header node
                                        info        ! XML pp_info node
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER                     :: ierr ! /= 0 if something went wrong
      CHARACTER(len=256) :: dft_buffer   ! needed to allow the string defining the
                                         ! DFT flavor to be longer than upf%dft 
                                         ! (currntly 25) without getting iotk upset. 
                                         ! An error message is issued if trimmed 
                                         ! dft_buffer exceeds upf%dft size.
                                         ! also use somewhere else as generic buffer ... 
      INTEGER :: len_buffer
      !
      INTEGER :: nw
      TYPE(Node), POINTER           :: tmpNode
      !
      ! Read HEADER section with some initialization data
      CALL searchData ( 'generated', upf%generated  , info ) 
      CALL searchData('creator',      upf%author,  info)
      tmpNode => item( getElementsByTagname( info, 'created'), 0)
      CALL extractDataAttribute( tmpNode, 'DATE', upf%date)  
      upf%comment =' ' 
         !
      CALL searchData('element',     upf%psd,     info)
      CALL searchData('type',        upf%typ,     info)
      CALL searchData('relativistic',upf%rel,     header)
         !
      CALL searchData ('is_ultrasoft',   upf%tvanp, header)
      CALL searchData('is_paw',          upf%tpawp, header)
      CALL searchData( 'is_coulomb',     upf%tcoulombp, header)
         !
      CALL searchData('has_so',         upf%has_so,   header)
      CALL searchData('has_wfc',        upf%has_wfc,  header)
      CALL searchData('has_gipaw',      upf%has_gipaw,header)
         !EMINE
      CALL searchData('paw_as_gipaw',   upf%paw_as_gipaw, header)
         !
      CALL searchData('core_correction',upf%nlcc, header )
      CALL searchData('functional',  dft_buffer, info)
         len_buffer=len_trim(dft_buffer)
         if (len_buffer > len(upf%dft)) &
            call errore('read_upf_schema','String defining DFT is too long',len_buffer)
         upf%dft=TRIM(dft_buffer)

      CALL searchData('z_valence',      upf%zp, header)
      CALL searchData('total_psenergy', upf%etotps,  header)
      CALL searchData('wfc_cutoff',     upf%ecutwfc,   header)
      CALL searchData('rho_cutoff',     upf%ecutrho,   header)
      CALL searchData('l_max',          upf%lmax,      header)
      CALL searchData('l_max_rho',      upf%lmax_rho,  header)
      CALL searchData('l_local',        upf%lloc,      header)
      CALL searchData('mesh_size',      upf%mesh,      header)
      CALL searchData('number_of_wfc',  upf%nwfc,      header) 
      CALL searchData('number_of_proj', upf%nbeta,     header)
      !
      RETURN
   END SUBROUTINE read_upf_header
   !
   SUBROUTINE read_upf_mesh(u, upf, grid)
      USE radial_grids, ONLY: allocate_radial_grid
      IMPLICIT NONE
      TYPE (Node),POINTER, INTENT(IN) :: u    ! pointer to XML node corresponding to the 
      TYPE(pseudo_upf),INTENT(INOUT)  :: upf  ! the pseudo data
      TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
      !
      TYPE ( Node )    , POINTER      :: locNode
      TYPE ( NodeList ), POINTER      :: nList 
      INTEGER                         :: ierr ! /= 0 if something went wrong
      LOGICAL :: found
      !

      CALL extractDataAttribute(u, 'dx',   upf%dx)
      CALL extractDataAttribute(u, 'mesh', upf%mesh)
      CALL extractDataAttribute(u, 'xmin', upf%xmin)
      CALL extractDataAttribute(u, 'rmax', upf%rmax)
      CALL extractDataAttribute(u, 'zmesh',upf%zmesh)
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
      locNode => item(getElementsByTagname ( u, 'pp_r' ),0)  
      CALL extractDataContent (locNode,   upf%r(1:upf%mesh))
      nList =>  getElementsByTagname ( u, 'pp_rab')
      IF ( getLength ( nList) .GT. 0 ) THEN
        locNode => item ( nList, 0 )  
        CALL extractDataContent(locNode, upf%rab(1:upf%mesh))
      END IF
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
   SUBROUTINE read_upf_nonlocal( u, upf)
      IMPLICIT NONE
      TYPE ( Node), POINTER, INTENT(IN)     :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT)        :: upf  ! the pseudo data
      !
      TYPE ( NodeList ), POINTER            :: beta_list, qij_list
      TYPE ( Node ),     POINTER            :: locNode, locNode2
      INTEGER :: nb,mb,ln,lm,l,nmb,ierr=0, iItem, index
      !INTEGER :: nb_=-1,mb_=-1,l_=-1,nmb_=-1
      INTEGER                               :: i,j,k
      CHARACTER(256)                        :: buf, temp
      LOGICAL :: isnull, found
      TYPE (DomException)                    :: ex_obj
      REAL(DP),ALLOCATABLE                   :: aux(:)
      !
      ! modified by AF
      !IF (upf%tcoulombp) RETURN
      beta_list => getElementsByTagname ( u, 'pp_beta') 
      upf%nbeta = getLength ( beta_list )  
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
         RETURN
      END IF
      !
      ! <AF>
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
       
      DO nb = 1, upf%nbeta
         locNode => item ( beta_list , nb-1 )
         CALL extractDataAttribute ( locNode,"index", index)
         CALL extractDataContent (locNode,  upf%beta(:,index) )
         IF ( hasAttribute( locNode, 'label') )&
            CALL extractDataAttribute(locNode, 'label',              upf%els_beta(index) )
         CALL extractDataAttribute(locNode, 'angular_momentum',       upf%lll(index))
         CALL extractDataAttribute(locNode, 'cutoff_radius_index',    upf%kbeta(index) )
         CALL extractDataAttribute(locNode, 'cutoff_radius',          upf%rcut(index) )
         IF (hasAttribute ( locNode, 'ultrasoft_cutoff_radius') ) THEN 
           CALL extractDataAttribute(locNode, 'ultrasoft_cutoff_radius', upf%rcutus(index) )
         ELSE IF (hasAttribute ( locNode, 'norm_conserving_radius') ) THEN
           CALL extractDataAttribute(locNode, 'ultrasoft_cutoff_radius', upf%rcutus(index) ) 
         END IF 
      ENDDO
      !
      ! Read the hamiltonian terms D_ij
      locNode => item ( getElementsByTagname(u, 'pp_dij'),0, EX = ex_obj ) 
      ierr = getExceptionCode (ex_obj) 
      IF ( ierr /= 0 ) CALL errore ('read_upf_schema', 'error reading pp_dij', ierr) 
      CALL extractDataContent ( locNode, upf%dion, IOSTAT = ierr )
      IF ( ierr /= 0 ) CALL errore ('read_upf_schema', 'format error in pp_dij element', ierr ) 
      !
      ! Read the augmentation charge section
      augmentation : IF(upf%tvanp .or. upf%tpawp) THEN
      !      
         locNode => item ( getElementsByTagname(u,'pp_augmentation'),0 )
         IF (.NOT. ASSOCIATED( locNode) ) CALL errore ('read_upf_schema', &
                                           'augmentation part not found in '// upf%typ //' pseudo '//upf%psd, 3)  
         CALL searchData('q_with_l', upf%q_with_l, locNode)
         CALL searchData('nqf', upf%nqf, locNode)
         IF (getLength( getElementsByTagname (locNode, 'nqlc')) .GT.0 ) THEN  
            CALL searchData('nqlc', upf%nqlc, locNode ) 
         ELSE 
            upf%nqlc = 2*upf%lmax+1 
         END IF 
         IF (upf%tpawp) THEN
            IF (getLength( getElementsByTagname(locNode, 'shape')) .GT. 0 ) THEN 
               CALL searchData('shape',   buf , locNode)
               temp=TRIM(buf)
               upf%paw%augshape = temp(1:12)
            ELSE 
               upf%paw%augshape = 'UNKNOWN'
            END IF
            IF (getLength( getElementsByTagname(locNode, 'cutoff_r')) .GT. 0 ) THEN
               CALL searchData('cutoff_r',       upf%paw%raug, locNode)
            ELSE 
               upf%paw%raug = 0
            END IF
            IF (getLength( getElementsByTagname(locNode, 'cutoff_r_index')) .GT. 0 ) THEN
               CALL searchData('cutoff_r_index', upf%paw%iraug,  locNode)
            ELSE 
               upf%paw%iraug = upf%mesh
            END IF 
            IF (getLength( getElementsByTagname(locNode, 'l_max_aug')) .GT. 0 ) THEN
               CALL searchData('l_max_aug',      upf%paw%lmax_aug, locNode) 
            ELSE 
               upf%paw%lmax_aug=upf%lmax_rho
            END IF 
         ENDIF
         ! a negative number means that all qfunc are stored
         IF ( getLength ( getElementsByTagname ( locNode, 'augmentation_epsilon') ) > 0 ) THEN 
            CALL searchData('augmentation_epsilon',upf%qqq_eps, locNode) 
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
         locNode2 => item(getElementsByTagname ( locNode, 'pp_q' ),0)
         IF (.NOT. ASSOCIATED(locNode2) ) CALL errore ('read_upf_schema', 'pp_q not found in pseudo '//upf%psd, 4)  
         CALL extractDataContent (locNode2 ,upf%qqq )
         !
         ! read charge multipoles (only if PAW)
         IF( upf%tpawp ) THEN
            ALLOCATE(upf%paw%augmom(upf%nbeta,upf%nbeta, 0:2*upf%lmax))
            locNode2 => item(getElementsByTagname( locNode, 'pp_multipoles'),0) 
         IF (.NOT. ASSOCIATED(locNode2) ) CALL errore ( 'read_upf_schema', &
                                               'multipoles not found in pseudo '//upf%psd, 5)
         ALLOCATE( aux(size(upf%paw%augmom,1)*size(upf%paw%augmom,2)*size(upf%paw%augmom,3)) )
         CALL extractDataContent(locNode2, aux, IOSTAT = ierr)
         IF ( ierr /= 0 ) &
            CALL errore ('read_upf_schema', 'format error reading multipoles for pseudo '//upf%psd, ierr)
         DO k =0,2*upf%nbeta
            DO j =1,upf%nbeta
               DO i = 1,2*upf%nbeta
                  upf%paw%augmom(i,j,k)=aux(i+(j-1)*upf%nbeta+k*(upf%nbeta*upf%nbeta))
               END DO
            END DO
         END DO
         DEALLOCATE(aux)
      ENDIF
      !
      ! Read polinomial coefficients for Q_ij expansion at small radius
      IF(upf%nqf <= 0) THEN
         upf%rinner(:) = 0._dp
         ALLOCATE( upf%qfcoef(1,1,1,1) )
         upf%qfcoef = 0._dp
      ELSE
         CALL errore ('read_upf_schema', 'found positive nqf for pseudo '//upf%psd, 6) 
      ENDIF
      !
      ! Read augmentation charge Q_ij
      ultrasoft_or_paw : IF( upf%tvanp) THEN
         q_with_l : IF( upf%q_with_l ) THEN
            upf%qfuncl = 0._dp
            qij_list => getElementsByTagname ( locNode, 'pp_qijl')
            DO iItem = 1, getLength(qij_list) 
               locNode2 => item(qij_list,iItem - 1) 
               CALL extractDataAttribute (locNode2,'composite_index', nmb) 
               CALL extractDataAttribute(locNode2,'angular_momentum', l ) 
               CALL extractDataContent  (locNode2, upf%qfuncl(:,nmb,l) )
               IF( upf%tpawp) upf%qfuncl(upf%paw%iraug+1:,nmb,l) = 0._dp
            ENDDO
         ELSE q_with_l
            upf%qfunc=0._dp
            qij_list =>  getElementsByTagname ( locNode, 'pp_qij' )
            DO iItem  = 1, getLength ( qij_list) 
               locNode2 => item( qij_list, iItem - 1)
               CALL extractDataAttribute( locNode2, 'composite_index',nmb)    
               CALL extractDataContent(locNode2, upf%qfunc(:,nmb))
            END DO
         ENDIF q_with_l
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
      TYPE( Node),POINTER, INTENT(IN) :: u    ! xml node with pseudo wfc data
      TYPE(pseudo_upf),INTENT(INOUT)  :: upf  ! the pseudo data
      INTEGER                         :: ierr ! /= 0 if something went wrong
      !
      INTEGER                         :: nw, nw_, nwfc_
      TYPE( NodeList),POINTER         :: locList 
      TYPE ( Node),POINTER            :: locNode
      !
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
      locList => getElementsByTagname ( u, "pp_chi") 
      nwfc_ = getLength(locList) 
      IF (nwfc_ /= upf%nwfc) CALL errore ( 'read_upf_schema', &
                                 'npwfc in header inconsistent with chi funcs found in '//upf%psd,5) 
      DO nw_ = 1,upf%nwfc
         locNode => item( locList, nw_ -1) 
         CALL extractDataAttribute ( locNode, "index", nw) 
         CALL extractDataContent (  locNode, upf%chi(:,nw), IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi function', ierr )
         CALL extractDataAttribute( locNode, 'label',  upf%els(nw) , IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi label ', ierr ) 
         CALL extractDataAttribute( locNode, 'l',  upf%lchi(nw) , IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi angular momentum l', ierr )
         CALL extractDataAttribute( locNode, 'occupation',  upf%oc(nw) , IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi occupation', ierr )
         CALL extractDataAttribute( locNode, 'n',  upf%nchi(nw) , IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi n number', ierr )
         CALL extractDataAttribute( locNode, 'pseudo_energy',  upf%epseu(nw) , IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi pseudo energy', ierr )
         CALL extractDataAttribute( locNode, 'cutoff_radius',  upf%rcut_chi(nw) , IOSTAT = ierr )
         CALL errore ('upf_read_schema', 'error reading chi rcut', ierr )
         IF (hasAttribute( locNode, 'ultrasoft_cutoff_radius')) THEN 
            CALL extractDataAttribute( locNode, 'ultrasoft_cutoff_radius',  upf%rcutus_chi(nw) , IOSTAT = ierr )
            CALL errore ('upf_read_schema', 'error reading chi US_rcut', ierr )
         END IF
      END DO
      !
      !
      RETURN
   END SUBROUTINE read_upf_pswfc

   SUBROUTINE read_upf_full_wfc(u, upf)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN)  :: u    ! XML Node containing full ae wave functions.
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      LOGICAL :: exst(15) = .FALSE. 
      !
      INTEGER :: nb,index, l 
      TYPE(NodeList), POINTER           :: locList
      TYPE(Node) ,    POINTER           :: locNode
      !
      !
      !
      ALLOCATE( upf%aewfc(upf%mesh, upf%nbeta) )
      locList => getElementsByTagname(u, "pp_aewfc")
      IF ( getLength(locList) /= upf%nbeta) CALL errore ('read_upf_schema', 'number of AE wfc not equal to nbeta',1)   
      DO nb = 1,upf%nbeta
         locNode => item( locList, nb -1)
         CALL extractDataAttribute(locNode, "index", index) 
         CALL extractDataContent ( locNode, upf%aewfc(:,index))  
      ENDDO

      IF (upf%has_so .and. upf%tpawp) THEN
         ALLOCATE( upf%paw%aewfc_rel(upf%mesh, upf%nbeta) )
         locList => getElementsByTagname(u,"pp_aewfc_rel")
        
         list_loop: DO nb = 1,getLength( locList) 
            locNode => item( locList, nb -1 )
            CALL extractDataAttribute ( locNode, "index", index)
            CALL extractDataContent ( locNode, upf%paw%aewfc_rel(:,index) )
            exst(index) = .TRUE. 
         END DO list_loop
         beta_loop: DO nb =1, upf%nbeta
            IF (.NOT. exst(nb) ) upf%paw%aewfc_rel=0.0_DP
         END DO beta_loop 
      END IF 

      ALLOCATE( upf%pswfc(upf%mesh, upf%nbeta) )
      locList => getElementsByTagname(u, 'pp_pswfc') 
      IF ( getLength(locList) /= upf%nbeta) CALL errore ('read_upf_schema', 'number of PS wfc not equal to nbeta',1)  
      DO nb = 1,upf%nbeta
         locNode => item( locList, nb -1 ) 
         CALL extractDataAttribute( locNode, "index", index)
         CALL extractDataContent( locNode, upf%pswfc(:,index) )
      ENDDO
      !
   END SUBROUTINE read_upf_full_wfc

   !
   SUBROUTINE read_upf_spin_orb(u, upf)
      IMPLICIT NONE
      TYPE (Node),POINTER,INTENT(IN) :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      TYPE (NodeList),POINTER    :: locList
      TYPE( Node ), POINTER      :: locNode
      INTEGER :: nw, nb, index
      !
      IF (.not. upf%has_so) RETURN
      !
      ALLOCATE (upf%nn(upf%nwfc))
      ALLOCATE (upf%jchi(upf%nwfc))
      !
      locList => getElementsByTagname(u, 'pp_relwfc') 
      IF ( getLength( locList) /= upf%nwfc ) &
         CALL errore ('read_upf_schema', 'in pp_spinorb section relwfc labels are less than nwfc '//upf%psd,-1)  
      DO nw = 1,upf%nwfc
         locNode => item( locList, nw -1)
         CALL extractDataAttribute( locNode, "index", index)
         CALL extractDataAttribute( locNode, 'nn', upf%nn(index))               
         CALL extractDataAttribute( locNode, 'jchi', upf%jchi(index)) 
      ENDDO
      !
      ALLOCATE(upf%jjj(upf%nbeta))
      !
      locList => getElementsByTagname(u, 'pp_relbeta')
      IF ( getLength ( locList) /= upf%nbeta )&
         CALL errore ('read_upf_schema', 'in pp_spinorb section relbeta labels are less than nbeta '//upf%psd,-1)
      DO nb = 1,upf%nbeta
         locNode => item(locList, nb -1 ) 
         CALL extractDataAttribute ( locNode, "index", index ) 
         CALL extractDataAttribute ( locNode, "lll", upf%lll(index))
         CALL extractDataAttribute ( locNode, "jjj", upf%jjj(index)) 
      ENDDO
      !
      RETURN
   END SUBROUTINE read_upf_spin_orb
   !
   SUBROUTINE read_upf_paw(u, upf)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN)  :: u    ! XML node with paw info
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      INTEGER :: nb,nb1
      TYPE(NodeList ),POINTER        :: locList
      TYPE(Node),POINTER             :: locNode
      ! 
      IF (.not. upf%tpawp ) RETURN
      CALL extractDataAttribute (u, "paw_data_format", upf%paw_data_format)
      IF(upf%paw_data_format /= 2) &
         CALL errore('read_upf_schema::paw',&
                     'Unknown format of PAW data.',1)
      IF (hasAttribute(u,'core_energy')) THEN 
         CALL extractDataAttribute( u, 'core_energy', upf%paw%core_energy)
      ELSE 
        upf%paw%core_energy = 0._DP 
      END IF
      ! Full occupation (not only > 0 ones)
      ALLOCATE( upf%paw%oc(upf%nbeta) )
      locNode = item ( getElementsByTagname(u, 'pp_occupations'),0) 
      IF (.NOT. ASSOCIATED(locNode)) CALL errore ('read_schema_upf', 'pp_occupations not found '//upf%psd, -1)
      CALL extractDataContent ( locNode, upf%paw%oc, IOSTAT = ierr ) 
      IF (ierr /= 0 ) CALL errore ('read_schema_upf', 'error reading pp_occupations '//upf%psd, ierr)  
      !
      ! All-electron core charge
      ALLOCATE( upf%paw%ae_rho_atc(upf%mesh) )
      locNode = item ( getElementsByTagname(u, 'pp_ae_nlcc'), 0) 
      IF (.NOT. ASSOCIATED(locNode)) CALL errore ('read_schema_upf', 'pp_ae_nlcc not found '//upf%psd, -1)
      CALL extractDataContent(locNode, upf%paw%ae_rho_atc, IOSTAT = ierr)
      IF (ierr /= 0 ) CALL errore ('read_schema_upf', 'error reading pp_ae_nlcc '//upf%psd, ierr)
      !
      ! All-electron local potential
      ALLOCATE( upf%paw%ae_vloc(upf%mesh) )
      locNode = item (getElementsByTagname(u, 'pp_ae_vloc'), 0)
      IF (.NOT. ASSOCIATED(locNode)) CALL errore ('read_schema_upf', 'pp_ae_vloc not found '//upf%psd, -1)
      CALL extractDataContent(locNode, upf%paw%ae_vloc, IOSTAT = ierr)
      IF (ierr /= 0 ) CALL errore ('read_schema_upf', 'error reading pp_ae_nlcc '//upf%psd, ierr)
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
      RETURN
   END SUBROUTINE read_upf_paw
   !
   SUBROUTINE read_upf_gipaw(u, upf)
      IMPLICIT NONE
      TYPE(Node),POINTER,INTENT(IN)  :: u    ! XML node with gipaw data
      TYPE(pseudo_upf),INTENT(INOUT) :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      TYPE (NodeList), POINTER       :: locList
      TYPE (Node), POINTER           :: locNode, locNode1, locNode2
      INTEGER :: nb, index
      IF (.not. upf%has_gipaw ) RETURN
      
      CALL extractDataAttribute ( u, 'gipaw_data_format', upf%gipaw_data_format )
      IF(upf%gipaw_data_format /= 2) &
         CALL infomsg('read_upf_schema::gipaw','Unknown format version')
      !
      CALL searchData ('number_of_core_orbitals', upf%gipaw_ncore_orbitals, u)
      CALL searchData ('number_of_valence_orbitals', upf%gipaw_wfs_nchannels, u) 
      ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
      ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
      locList => getElementsByTagname ( u, 'pp_gipaw_core_orbital' )
      IF ( getLength (locList) /= upf%gipaw_ncore_orbitals ) &
         CALL errore( 'read_upf_schema::gipaw', 'wrong number_of_core_orbitals', 1)  
      DO nb = 1, upf%gipaw_ncore_orbitals
            locNode => item (locList,nb -1)    
            CALL extractDataAttribute ( locNode, 'index', index)
            CALL extractDataContent (locNode, upf%gipaw_core_orbital(:,index))
            CALL extractDataAttribute (locNode, 'label', upf%gipaw_core_orbital_el(index))
            CALL extractDataAttribute (locNode, 'n', upf%gipaw_core_orbital_n(index))
            CALL extractDataAttribute (locNode, 'l', upf%gipaw_core_orbital_l(index))
      END DO

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
         locList => getElementsByTagname(u, 'pp_gipaw_orbital') 
         IF ( getLength ( locList ) /=  upf%gipaw_wfs_nchannels)   &
            CALL errore( 'read_upf_schema::gipaw', 'wrong number_of_valence_orbitals',1) 
         ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
         ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
         DO nb = 1,upf%gipaw_wfs_nchannels
            locNode => item(locList, nb -1) 
            CALL extractDataAttribute(locNode, 'index', index)
            CALL extractDataAttribute(locNode, 'label', upf%gipaw_wfs_el(index) )
            CALL extractDataAttribute(locNode, 'l', upf%gipaw_wfs_ll(index) )
            CALL extractDataAttribute(locNode, 'cutoff_radius', upf%gipaw_wfs_rcut(index))
            IF ( hasAttribute ( locNode, 'ultrasoft_cutoff_radius') ) THEN 
               CALL extractDataAttribute(locNode, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(index)) 
            ELSE
                upf%gipaw_wfs_rcutus(index)  = upf%gipaw_wfs_rcut(index)
            END IF
            ! read all-electron orbital
            locNode2 => item( getElementsByTagname( locNode, 'pp_gipaw_wfs_ae'),0) 
               IF (.NOT. ASSOCIATED( locNode2) ) CALL errore ('read_upf_schema::gipaw', 'pp_gipaw_wfs_ae not found',-1)
               CALL extractDataContent(locNode2, upf%gipaw_wfs_ae(:,index))
            !
            ! read pseudo orbital
            locNode2 => item( getElementsByTagname( locNode, 'pp_gipaw_wfs_ps'),0)
              IF (.NOT. ASSOCIATED( locNode2) ) CALL errore ('read_upf_schema::gipaw', 'pp_gipaw_wfs_ps not found',-1)
              CALL extractDataContent(locNode2, upf%gipaw_wfs_ps(:,index))
            !
         ENDDO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read all-electron and pseudo local potentials
         ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
         ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
         locNode => item ( getElementsByTagname(u, 'pp_gipaw_vlocal'),0) 
            IF ( .NOT. ASSOCIATED (locNode) ) CALL errore( 'read_upf_schema::gipaw', 'pp_gipaw_vlocal not found', -1)
            locNode2 => item( getElementsByTagname(locNode, 'pp_gipaw_vlocal_ae'),0)
               IF ( .NOT. ASSOCIATED( locNode2) ) &
                  CALL errore( 'read_upf_schema::gipaw', 'pp_gipaw_vlocal_ae not found', -1)
               CALL extractDataContent(locNode2, upf%gipaw_vlocal_ae(:))
           !
           locNode2 => item( getElementsByTagname(locNode, 'pp_gipaw_vlocal_ps'),0)
               IF ( .NOT. ASSOCIATED( locNode2) ) &
                  CALL errore( 'read_upf_schema::gipaw', 'pp_gipaw_vlocal_ps not found', -1)
               CALL extractDataContent(locNode2, upf%gipaw_vlocal_ps(:))
           !
      ENDIF
      RETURN
   END SUBROUTINE read_upf_gipaw
   !
   !
   subroutine searchStringData (tagname_, outstr, node_point, error)
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(OUT)          :: outstr
   CHARACTER(LEN=*), INTENT(IN)           :: tagname_
   TYPE ( NODE), POINTER, INTENT(IN)      :: node_point
   INTEGER, OPTIONAL,INTENT(OUT)          :: error
   !
   INTEGER                                :: ierr, io_err, lenstr 
   TYPE(Node), POINTER                    :: point
   TYPE (DomException)                    :: exception_obj
   point => item(getElementsByTagname(node_point, trim(tagname_)),0)
   outstr = getTextContent(point, ex = exception_obj)
   ierr =   getExceptionCode(exception_obj)
   IF ( PRESENT (error) ) THEN 
      error = ierr 
   ELSE IF (ierr /= 0 ) THEN 
       CALL errore ('read_upf_schema', 'error getting  '// tagname_, ierr)
   END IF
   END SUBROUTINE searchStringData
   
   SUBROUTINE searchBooleanData(tagname_, out_bool, node_point, error) 
   IMPLICIT NONE
   LOGICAL, INTENT(OUT)                   :: out_bool 
   CHARACTER(LEN=*), INTENT(IN)           :: tagname_
   TYPE ( Node), POINTER, INTENT(IN)      :: node_point
   INTEGER, OPTIONAL, INTENT(OUT)         :: error 
   INTEGER                                :: ierr, io_err  
   TYPE(Node),POINTER                     :: point
   TYPE (DomException)                    :: exception_obj
   !
   point => item(getElementsByTagname(node_point, trim(tagname_) ),0)
   call extractDataContent(point, out_bool, EX = exception_obj, IOSTAT = io_err)
   ierr = getExceptionCode(exception_obj)
   IF ( PRESENT(error) ) THEN 
      IF ( ierr /=0 ) error = ierr
   ELSE IF ( ierr /= 0 ) THEN
       CALL errore ('read_upf_schema','error reading '//tagname_, ierr) 
   END IF
   IF ( io_err /= 0) CALL errore ('read_upf_schema','content error reading '//tagname_, ierr)
   END SUBROUTINE searchBooleanData
   
   
   SUBROUTINE searchRealData(tagname_, out_real, node_point, error )
   IMPLICIT NONE
   real(8), intent(out)                   :: out_real
   character(len=*), intent(in)           :: tagname_
   type ( Node), pointer, intent(in)      :: node_point
   INTEGER, OPTIONAL, INTENT(OUT)         :: error
   !
   INTEGER                                :: ierr, io_err
   TYPE(Node),POINTER                     :: point
   TYPE (DomException)                    :: exception_obj
   !
   point => item(getElementsByTagname(node_point, trim(tagname_) ),0)
   call extractDataContent(point, out_real, ex = exception_obj, IOSTAT = io_err)
   ierr = getExceptionCode(exception_obj)
   IF ( PRESENT ( error) ) THEN 
      error = ierr 
   ELSE IF ( ierr /= 0 ) THEN 
       CALL errore ( 'read_upf_schema', 'error reading '// tagname_, ierr ) 
   END IF 
   IF (io_err /= 0 ) CALL errore ( 'read_upf_schema', 'format error in '// tagname_, io_err) 
   END SUBROUTINE searchRealData
   
   SUBROUTINE searchIntegerData(tagname_, out_int, node_point, error)
   IMPLICIT NONE
   integer, intent(out)        :: out_int
   character(len=*), intent(in)    :: tagname_
   type ( Node), pointer, intent(in)      :: node_point
   INTEGER, OPTIONAL, INTENT(OUT)         :: error 
   !            
   INTEGER                                :: ierr, io_err                           
   TYPE(Node),POINTER                     :: point
   TYPE (DomException)                    :: exception_obj
   
   
   point => item(getElementsByTagname(node_point, trim(tagname_) ),0)
   call extractDataContent(point, out_int, EX = exception_obj, IOSTAT = io_err )
   ierr = getExceptionCode(exception_obj)
   IF (PRESENT ( error )) THEN 
       error = ierr 
   ELSE IF ( ierr /=0 ) THEN 
       CALL errore ('read_upf_schema', 'error reading '// tagname_, ierr ) 
   END IF
   IF ( io_err /= 0 ) CALL errore ( 'read_upf_schema', 'error reading '// tagname_, ierr )
   END SUBROUTINE searchIntegerData
   !
END MODULE read_upf_schema_module
