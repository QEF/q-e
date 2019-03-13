!--------------------------------------------------------------------
!
! 
! Program written by Yang Jiao, Nov 2018, GPL, No warranties.
!    partly ported from QE
!
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
!-----------------------------------------------------------------------
MODULE vasp_xml
!----------------------------------------------------------------------
  ! 
  ! ... this module contains subroutine to read data produced by VASP 
  ! ...    tested on VASP.5.3.3
  !
  !
  ! vasp_read_*                read variables from xml file
  ! vasp_readschema_*          read variables into internal varables
  !
USE kinds, ONLY : DP

IMPLICIT NONE

INTEGER          :: iunvasp, iunvaspchg 

PRIVATE
PUBLIC        :: readxmlfile_vasp
!
  !---------------------------------------------------------

! types for vasp input data
!
  TYPE :: vasp_kpoints_type
    !
    CHARACTER(len=100)             :: tagname
    INTEGER                        :: nk
    REAL(DP), ALLOCATABLE          :: xk(:,:)
    REAL(DP), ALLOCATABLE          :: wk(:)
    !
  END TYPE vasp_kpoints_type
  !
  TYPE :: vasp_parameters_type
    !
    CHARACTER(len=100) :: tagname
    CHARACTER(len=2)               :: gga
    INTEGER                        :: nbands
    INTEGER                        :: ispin
    INTEGER                        :: ngx
    INTEGER                        :: ngy
    INTEGER                        :: ngz
    INTEGER                        :: ngxf
    INTEGER                        :: ngyf
    INTEGER                        :: ngzf
    REAL(DP)                       :: enmax
    REAL(DP)                       :: aldax
    REAL(DP)                       :: aggax
    REAL(DP)                       :: aldac
    REAL(DP)                       :: aggac
    REAL(DP)                       :: zab_vdw
    REAL(DP)                       :: param1
    REAL(DP)                       :: param2
    REAL(DP)                       :: param3
    LOGICAL                        :: lnoncollinear
    LOGICAL                        :: lmetagga
    LOGICAL                        :: luse_vdw
    !
  END TYPE vasp_parameters_type
  !
  TYPE :: vasp_atominfo_type
    !
    CHARACTER(len=100)             :: tagname
    INTEGER                        :: nat          ! number of atoms
    INTEGER                        :: nsp          ! number of species
    INTEGER,ALLOCATABLE            :: ityp(:)      ! number of atoms per type
    REAL(DP),ALLOCATABLE           :: zv(:)
    CHARACTER(LEN=3),ALLOCATABLE   :: atm(:)
    CHARACTER(LEN=10)              :: pseudo
    !
  END TYPE vasp_atominfo_type
  !
  TYPE :: vasp_structure_type
    !
    CHARACTER(len=100)             :: tagname
    INTEGER                        :: nat
    REAL(DP)                       :: at(3,3)
    REAL(DP)                       :: volume
    REAL(DP)                       :: bg(3,3)
    REAL(DP),ALLOCATABLE           :: tau(:,:)
    !
  END TYPE vasp_structure_type
  !


CONTAINS
!----------------------------------------------------------------------
SUBROUTINE readxmlfile_vasp(iexch,icorr,igcx,igcc,inlc,ierr)
  !----------------------------------------------------------------------
  USE ions_base,            ONLY : nat, nsp, ityp, tau, atm
  USE cell_base,            ONLY : tpiba2, alat,omega, at, bg, ibrav
  USE funct,                ONLY : set_dft_from_indices
  USE klist,                ONLY : nkstot, nks, xk, wk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, nbndx, et, wg
  USE symm_base,            ONLY : irt, d1, d2, d3, checkallsym, nsym
  USE extfield,             ONLY : forcefield, tefield, gate, forcegate
  USE cellmd,               ONLY : cell_factor, lmovecell
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE fft_types,            ONLY : fft_type_allocate
  USE recvec_subs,          ONLY : ggen, ggens
  USE gvect,                ONLY : gg, ngm, g, gcutm, mill, ngm_g, ig_l2g, &
                                   eigts1, eigts2, eigts3, gstart, gshells
  USE fft_base,             ONLY : dfftp, dffts
  USE gvecs,                ONLY : ngms, gcutms
  USE spin_orb,             ONLY : lspinorb, domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  USE wavefunctions,        ONLY : psic
  USE vlocal,               ONLY : strf
  USE io_files,             ONLY : tmp_dir, prefix, iunpun, nwordwfc, iunwfc
  USE io_global,            ONLY : stdout
  USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda, nspin_mag, nspin_gga
  USE pw_restart_new,       ONLY :  pw_readschema_file, init_vars_from_schema
  USE qes_types_module,     ONLY :  output_type, parallel_info_type, general_info_type, input_type
  USE qes_libs_module,      ONLY :  qes_reset 
  USE io_rho_xml,           ONLY : read_scf
  USE fft_rho,              ONLY : rho_g2r
  USE uspp,                 ONLY : becsum
  USE uspp_param,           ONLY : upf
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_init,             ONLY : paw_init_onecenter, allocate_paw_internals
  USE control_flags,        ONLY : gamma_only
  USE funct,                ONLY : get_inlc, get_dft_name
  USE kernel_table,         ONLY : initialize_kernel_table
  USE esm,                  ONLY : do_comp_esm, esm_init
  USE mp_bands,             ONLY : intra_bgrp_comm, nyfft
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_fact
  USE vasp_read_chgcar,     ONLY : vaspread_rho
  !
  IMPLICIT NONE

  INTEGER,           INTENT(out)  :: iexch,icorr,igcx,igcc,inlc
  INTEGER,           INTENT(out)  :: ierr
  INTEGER  :: i, is, ik, ibnd, nb, nt, ios, isym
  INTEGER  :: ngx_, ngy_, ngz_
  INTEGER  :: ngxf_, ngyf_, ngzf_
  INTEGER  :: imeta = 0
  REAL(DP) :: rdum(1,1), ehart, etxc, vtxc, etotefield, charge
  REAL(DP) :: sr(3,3,48)
  CHARACTER(LEN=20) dft_name
  TYPE(vasp_kpoints_type)               :: vasp_kpoints_obj
  TYPE(vasp_parameters_type)            :: vasp_parameters_obj
  TYPE(vasp_atominfo_type)              :: vasp_atominfo_obj
  TYPE(vasp_structure_type)             :: vasp_structure_obj
  !
  !
  CALL vasp_readschema_file ( ierr, vasp_kpoints_obj, vasp_parameters_obj,    &
                                    vasp_atominfo_obj, vasp_structure_obj)
  IF ( ierr /= 0 ) CALL errore ( 'read_schema', 'unable to read xml file', ierr )
  !------------------------------------------------------------------
  ngx_=vasp_parameters_obj%ngx
  ngy_=vasp_parameters_obj%ngy
  ngz_=vasp_parameters_obj%ngz
  ngxf_=vasp_parameters_obj%ngxf
  ngyf_=vasp_parameters_obj%ngyf
  ngzf_=vasp_parameters_obj%ngzf
  !------------------------------------------------------------------
  CALL vasp_init_xc(vasp_parameters_obj,vasp_atominfo_obj,iexch,icorr,igcx,igcc,inlc,ierr)
  !------------------------------------------------------------------
  !
  CALL vasp_init_vars_from_schema( 'dim', ierr , vasp_kpoints_obj, vasp_parameters_obj, &
                                    vasp_atominfo_obj, vasp_structure_obj)
  CALL errore( 'read_xml_file ', 'problem reading file ' // TRIM( tmp_dir ) //'vasprun.xml', ierr )
  !
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau( 3, nat ) )

  CALL vasp_init_vars_from_schema( 'atom', ierr , vasp_kpoints_obj, vasp_parameters_obj, &
                                    vasp_atominfo_obj, vasp_structure_obj)
  CALL errore( 'read_xml_file ', 'problem reading file ' // TRIM( tmp_dir ) //'vasprun.xml', ierr )
  !
  !
  CALL vasp_init_vars_from_schema( 'kpoint', ierr , vasp_kpoints_obj, vasp_parameters_obj, &
                                    vasp_atominfo_obj, vasp_structure_obj)
  CALL errore( 'read_xml_file ', 'problem reading file ' // TRIM( tmp_dir ) //'vasprun.xml', ierr )
  !
  CALL set_dft_from_indices(iexch, icorr, igcx, igcc, inlc)
  WRITE( stdout, '(5X,"Exchange-correlation      = ", &
        &  " (",I2,3I3,2I2,")")') iexch,icorr,igcx,igcc,inlc,imeta
  !
  ! ... read the vdw kernel table if needed
  !
  IF (inlc > 0 ) THEN
     call initialize_kernel_table(inlc) 
  END IF
  !
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  CALL set_dimensions()  ! input: ecutwfc dual output: gcutm, gcutms
  CALL fft_type_setdim( dffts, ngxf_, ngyf_, ngzf_)
  CALL fft_type_setdim( dfftp, ngxf_, ngyf_, ngzf_)
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  gamma_only=.FALSE.   ! gamma_only tricks not used 
  CALL data_structure ( gamma_only ) 
  CALL allocate_fft()
  CALL ggen ( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
              g, gg, mill, ig_l2g, gstart )
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )

  CALL gshells ( .False. ) 

  CALL allocate_wfc()
     
  CALL vaspread_rho(ierr)
  !
  ! ignore core charge
  rhog_core(:) = 0.0_DP
  rho_core(:)  = 0.0_DP
  !
  RETURN

END SUBROUTINE readxmlfile_vasp
!
SUBROUTINE vasp_readschema_file (ierr, vasp_kpoints, vasp_parameters, vasp_atominfo, vasp_structure)
  !----------------------------------------------------------------
  USE io_files,             ONLY : tmp_dir
  USE FoX_dom,              ONLY : parseFile, item, getElementsByTagname, destroy, nodeList, Node
  USE FoX_dom,              ONLY : getLength
  USE FoX_dom,              ONLY : hasAttribute, getAttributes, extractDataAttribute
  IMPLICIT NONE
  !
  INTEGER                                               :: ierr, io_err
  TYPE(vasp_kpoints_type), OPTIONAL,       INTENT(OUT)  :: vasp_kpoints
  TYPE(vasp_parameters_type), OPTIONAL,    INTENT(OUT)  :: vasp_parameters
  TYPE(vasp_atominfo_type), OPTIONAL,      INTENT(OUT)  :: vasp_atominfo
  TYPE(vasp_structure_type), OPTIONAL,     INTENT(OUT)  :: vasp_structure
  !
  TYPE(Node), POINTER       :: root, nodePointer
  TYPE(nodeList), POINTER   :: listPointer
  LOGICAL                   :: found
  CHARACTER(LEN=80)         :: errmsg = ' '
  CHARACTER(LEN=320)        :: filename
  INTEGER,EXTERNAL          :: find_free_unit
  INTEGER                   :: listPointer_size, iitem
  CHARACTER(LEN=100)        :: attr_name
  !
  !
  ierr = 0
  io_err = 0
  !
  iunvasp = find_free_unit()
  IF (iunvasp .LT. 0 ) &
     CALL errore ("vasp_readschema_file", "could not fine a free unit to open vasprun.xml", 1)
  filename = TRIM( tmp_dir ) // "vasprun.xml"
  INQUIRE ( file=filename, exist=found )
  IF (.NOT. found ) ierr = ierr + 1
  IF ( ierr /= 0 ) THEN
     errmsg='xml data file not found'
     GOTO 100
  END IF
  !
  root => parseFile(filename)
  !
  IF ( PRESENT (vasp_kpoints ) ) THEN
     nodePointer => item ( getElementsByTagname(root, "kpoints"),0)
     CALL vasp_read_kpoints( nodePointer, vasp_kpoints, ierr)
     IF ( ierr /= 0 ) THEN
        errmsg='error kpoints of xml data file'
        GOTO 100
     END IF
  END IF
  !
  IF ( PRESENT (vasp_parameters ) ) THEN
     nodePointer => item ( getElementsByTagname(root, "parameters"),0)
     CALL vasp_read_parameters( nodePointer, vasp_parameters, ierr)
     IF ( ierr /= 0 ) THEN
        errmsg='error parameters of xml data file'
        GOTO 100
     END IF
  END IF
  !
  IF ( PRESENT ( vasp_atominfo ) ) THEN
     nodePointer => item ( getElementsByTagname(root, "atominfo"),0)
     CALL vasp_read_atominfo( nodePointer, vasp_atominfo, ierr)
     IF ( ierr /= 0 ) THEN
        errmsg='error atominfo of xml data file'
        GOTO 100
     END IF
  END IF
  !
  IF ( PRESENT ( vasp_structure ) ) THEN
     listPointer => getElementsByTagname(root, "structure")
     listPointer_size = getLength(listPointer)
     DO iitem = 0, listPointer_size-1
        nodePointer => item ( listPointer, iitem )
        IF (hasAttribute(nodePointer, "name")) THEN
           CALL extractDataAttribute(nodePointer, "name", attr_name)
           IF(attr_name == "finalpos") THEN
              CALL vasp_read_structure( nodePointer, vasp_structure, ierr)
           END IF
        END IF
     END DO
     IF ( ierr /= 0 ) THEN
        errmsg='error structure of xml data file'
        GOTO 100
     END IF
  END IF
  !
  CALL destroy(root)
  !
 100  CALL errore('vasp_readschemafile',TRIM(errmsg),ierr)
  !
END SUBROUTINE vasp_readschema_file
!
  !---------------------------------------------------------
  SUBROUTINE vasp_read_atominfo(xml_node, obj, ierr)
    !
    USE FoX_dom,              ONLY : item, getElementsByTagname, nodeList, Node
    USE FoX_dom,              ONLY : getTagName, getLength, extractDataContent
    USE FoX_dom,              ONLY : hasAttribute, getAttributes, extractDataAttribute
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER          :: xml_node
    TYPE(vasp_atominfo_type), INTENT(OUT)    :: obj
    INTEGER, OPTIONAL, INTENT(OUT)           :: ierr
    !
    TYPE(Node), POINTER        :: tmp_node, sub_tmp_node, ssub_tmp_node, sssub_tmp_node
    TYPE(NodeList),  POINTER   :: tmp_node_list, sub_tmp_node_list, ssub_tmp_node_list, sssub_tmp_node_list
    INTEGER :: tmp_node_list_size, sub_tmp_node_list_size, ssub_tmp_node_list_size, sssub_tmp_node_list_size
    INTEGER :: index, iostat_
    INTEGER :: iitem, i_sitem, i_ssitem
    CHARACTER(LEN=100) :: attr_name, pseudofile_ 
    CHARACTER(LEN=10)  :: pseudo_
    !
    obj%tagname = getTagName(xml_node)
    !
    tmp_node_list => getElementsByTagname(xml_node, "atoms")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1 ) THEN
       IF (PRESENT(ierr) ) THEN
          CALL infomsg("vasp_read_atominfo","atoms: wrong number of occurrences")
          ierr = ierr + 1
       ELSE
          CALL errore("vasp_read_atominfo","atoms: wrong number of occurrences",10)
       END IF
    END IF
    ! 
    tmp_node => item(tmp_node_list,0)
    IF (ASSOCIATED(tmp_node)) THEN
       CALL extractDataContent(tmp_node, obj%nat )
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "types")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1 ) THEN
       IF (PRESENT(ierr) ) THEN
          CALL infomsg("vasp_read_atominfo","types: wrong number of occurrences")
          ierr = ierr + 1
       ELSE
          CALL errore("vasp_read_atominfo","types: wrong number of occurrences",10)
       END IF
    END IF
    !
    tmp_node => item(tmp_node_list,0)
    IF (ASSOCIATED(tmp_node)) THEN
       CALL extractDataContent(tmp_node, obj%nsp )
    END IF
    !
    ALLOCATE(obj%atm(1:obj%nsp), obj%ityp(1:obj%nat), obj%zv(1:obj%nsp))
    tmp_node_list => getElementsByTagname(xml_node, "array")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    DO iitem = 0, tmp_node_list_size-1
       tmp_node => item(tmp_node_list,iitem)
       IF (hasAttribute(tmp_node, "name")) THEN
          CALL extractDataAttribute(tmp_node, "name", attr_name)
          IF (attr_name == "atoms") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"set")
             sub_tmp_node => item(sub_tmp_node_list,0)
             ssub_tmp_node_list => getElementsByTagname(sub_tmp_node,"rc")
             ssub_tmp_node_list_size = getLength(ssub_tmp_node_list)
             IF (ssub_tmp_node_list_size /= obj%nat) THEN
                CALL errore("vasp_read_atominfo","atoms: wrong number of occurrences",10)
             END IF
             DO i_sitem = 0, ssub_tmp_node_list_size-1
                ssub_tmp_node => item(ssub_tmp_node_list, i_sitem)
                sssub_tmp_node_list => getElementsByTagname(ssub_tmp_node,"c")
                sssub_tmp_node_list_size = getLength(sssub_tmp_node_list)
             !   sssub_tmp_node => item(sssub_tmp_node_list,0)
             !   CALL extractDataContent(sssub_tmp_node,obj%atm(i_sitem+1))
                sssub_tmp_node => item(sssub_tmp_node_list,1)
                CALL extractDataContent(sssub_tmp_node,obj%ityp(i_sitem+1))
             END DO
          ELSE IF (attr_name == "atomtypes") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"set")
             sub_tmp_node => item(sub_tmp_node_list,0)
             ssub_tmp_node_list => getElementsByTagname(sub_tmp_node,"rc")
             ssub_tmp_node_list_size = getLength(ssub_tmp_node_list)
             IF (ssub_tmp_node_list_size /= obj%nsp) THEN
                CALL errore("vasp_read_atominfo","atomtypes: wrong number of occurrences",10)
             END IF
             DO i_sitem = 0, ssub_tmp_node_list_size-1
                ssub_tmp_node => item(ssub_tmp_node_list, i_sitem)
                sssub_tmp_node_list => getElementsByTagname(ssub_tmp_node,"c")
                sssub_tmp_node_list_size = getLength(sssub_tmp_node_list)
                IF (sssub_tmp_node_list_size /= 5) THEN
                   CALL errore("vasp_read_atominfo","atomtypes set: wrong number of occurrences",10)
                END IF
                sssub_tmp_node => item(sssub_tmp_node_list,1)
                CALL extractDataContent(sssub_tmp_node,obj%atm(i_sitem+1))
                sssub_tmp_node => item(sssub_tmp_node_list,3)
                CALL extractDataContent(sssub_tmp_node,obj%zv(i_sitem+1))
                sssub_tmp_node => item(sssub_tmp_node_list,4)
                CALL extractDataContent(sssub_tmp_node,pseudofile_)
                READ(pseudofile_,*) pseudo_
                IF (i_sitem == 0) THEN
                   obj%pseudo = pseudo_
                ELSE IF (pseudo_ /= obj%pseudo) THEN
                   CALL errore("vasp_read_atominfo", "atomtypes: inconsistent pseudo files",10)
                END IF
             END DO
          END IF
       END IF
    END DO
    !
    RETURN         
    ! 
  END SUBROUTINE vasp_read_atominfo
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  SUBROUTINE vasp_read_kpoints(xml_node, obj, ierr)
    !
    USE FoX_dom,              ONLY : item, getElementsByTagname, nodeList, Node
    USE FoX_dom,              ONLY : getTagName, getLength, extractDataContent
    USE FoX_dom,              ONLY : hasAttribute, getAttributes, extractDataAttribute
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER          :: xml_node
    TYPE(vasp_kpoints_type), INTENT(OUT)     :: obj
    INTEGER, OPTIONAL, INTENT(OUT)           :: ierr
    !
    TYPE(Node), POINTER        :: tmp_node, sub_tmp_node
    TYPE(NodeList),  POINTER   :: tmp_node_list, sub_tmp_node_list
    INTEGER :: tmp_node_list_size, sub_tmp_node_list_size
    INTEGER :: index, iostat_
    INTEGER :: iitem, i_sitem
    CHARACTER(LEN=100) :: attr_name
    !
    obj%tagname = getTagName(xml_node)
    !
    tmp_node_list => getElementsByTagname(xml_node, "varray")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    DO iitem = 0, tmp_node_list_size-1
       tmp_node => item(tmp_node_list,iitem)
       IF (hasAttribute(tmp_node, "name")) THEN
          CALL extractDataAttribute(tmp_node, "name", attr_name)
          IF (attr_name == "kpointlist") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"v")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             obj%nk=sub_tmp_node_list_size
             ALLOCATE(obj%xk(1:3,1:obj%nk))
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                CALL extractDataContent(sub_tmp_node, obj%xk(1:3,i_sitem+1))
             END DO
          END IF
          IF (attr_name == "weights") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"v")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             IF(obj%nk/=sub_tmp_node_list_size) THEN
                CALL errore("vasp_read_kpoints","weights: wrong number of occurrences",sub_tmp_node_list_size )
             END IF
             ALLOCATE(obj%wk(1:obj%nk))
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                CALL extractDataContent(sub_tmp_node, obj%wk(i_sitem+1))
             END DO
          END IF
       END IF
    END DO
    !
    RETURN
    !
  END SUBROUTINE vasp_read_kpoints
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  SUBROUTINE vasp_read_parameters(xml_node, obj, ierr)
    !
    USE FoX_dom,              ONLY : item, getElementsByTagname, nodeList, Node
    USE FoX_dom,              ONLY : getTagName, getLength, extractDataContent
    USE FoX_dom,              ONLY : hasAttribute, getAttributes, extractDataAttribute
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER          :: xml_node
    TYPE(vasp_parameters_type), INTENT(OUT)  :: obj
    INTEGER, OPTIONAL, INTENT(OUT)           :: ierr
    !
    TYPE(Node), POINTER        :: tmp_node, sub_tmp_node, ssub_tmp_node, sssub_tmp_node
    TYPE(NodeList),  POINTER   :: tmp_node_list, sub_tmp_node_list, ssub_tmp_node_list, sssub_tmp_node_list
    INTEGER :: tmp_node_list_size, sub_tmp_node_list_size, ssub_tmp_node_list_size, sssub_tmp_node_list_size
    INTEGER :: index, iostat_
    INTEGER :: iitem, i_sitem, i_ssitem
    CHARACTER(LEN=100) :: attr_name
    CHARACTER(LEN=3)   :: logical_
    !
    obj%tagname = getTagName(xml_node)
    !
    tmp_node_list => getElementsByTagname(xml_node, "separator")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    DO iitem = 0, tmp_node_list_size-1
       tmp_node => item(tmp_node_list,iitem)
       IF (hasAttribute(tmp_node, "name")) THEN
          CALL extractDataAttribute(tmp_node, "name", attr_name)
          IF (attr_name == "electronic") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"i")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                IF(hasAttribute(sub_tmp_node, "name")) THEN
                   CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
                   IF (attr_name == "ENMAX") THEN
                      CALL extractDataContent(sub_tmp_node,obj%enmax)
                   ELSE IF (attr_name == "NBANDS") THEN
                      CALL extractDataContent(sub_tmp_node,obj%nbands)
                   END IF
                END IF
             END DO
             sub_tmp_node_list => getElementsByTagname(tmp_node,"separator")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                IF(hasAttribute(sub_tmp_node, "name")) THEN
                   CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
                   IF (attr_name == "electronic spin") THEN
                      ssub_tmp_node_list => getElementsByTagname(sub_tmp_node,"i")
                      ssub_tmp_node_list_size = getLength(ssub_tmp_node_list)
                      DO i_ssitem = 0, ssub_tmp_node_list_size-1
                         ssub_tmp_node => item(ssub_tmp_node_list,i_ssitem)
                         IF(hasAttribute(ssub_tmp_node, "name")) THEN
                            CALL extractDataAttribute(ssub_tmp_node,"name",attr_name)
                            IF (attr_name=="ISPIN") THEN
                               CALL extractDataContent(ssub_tmp_node,obj%ispin)
                            ELSE IF (attr_name=="LNONCOLLINEAR") THEN
                               CALL extractDataContent(ssub_tmp_node,logical_)
                               IF(TRIM(logical_)=="F") THEN
                                  obj%lnoncollinear = .FALSE.
                               ELSE IF(TRIM(logical_)=="T") THEN
                                  obj%lnoncollinear = .TRUE.
                               ELSE
                                  CALL errore("vasp_read_parameters","LNONCOLLINEAR: wrong occurrence",10)
                               END IF
                            END IF
                         END IF
                      END DO
                   ELSE IF (attr_name == "electronic exchange-correlation") THEN
                      ssub_tmp_node_list => getElementsByTagname(sub_tmp_node,"i")
                      ssub_tmp_node_list_size = getLength(ssub_tmp_node_list)
                      DO i_ssitem = 0, ssub_tmp_node_list_size-1
                         ssub_tmp_node => item(ssub_tmp_node_list,i_ssitem)
                         IF(hasAttribute(ssub_tmp_node, "name")) THEN
                            CALL extractDataAttribute(ssub_tmp_node,"name",attr_name)
                            IF (attr_name=="LMETAGGA") THEN
                               CALL extractDataContent(ssub_tmp_node,logical_)
                               IF(TRIM(logical_)=="F") THEN
                                  obj%lmetagga = .FALSE.
                               ELSE IF(TRIM(logical_)=="T") THEN
                                  obj%lmetagga = .TRUE.
                                  CALL errore("vasp_read_parameters","LMETAGGA: metaGGA not implemented",10)  
                               ELSE
                                  CALL errore("vasp_read_parameters","LMETAGGA: wrong occurrence",10)
                               END IF
                            END IF 
                         END IF
                      END DO
                   END IF
                END IF
             END DO
          ELSE IF(attr_name == "grids") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"i")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                IF(hasAttribute(sub_tmp_node, "name")) THEN
                   CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
                   IF (attr_name == "NGX") THEN
                      CALL extractDataContent(sub_tmp_node,obj%ngx)
                   ELSE IF (attr_name == "NGY") THEN
                      CALL extractDataContent(sub_tmp_node,obj%ngy)
                   ELSE IF (attr_name == "NGZ") THEN
                      CALL extractDataContent(sub_tmp_node,obj%ngz)
                   ELSE IF (attr_name == "NGXF") THEN
                      CALL extractDataContent(sub_tmp_node,obj%ngxf)
                   ELSE IF (attr_name == "NGYF") THEN
                      CALL extractDataContent(sub_tmp_node,obj%ngyf)
                   ELSE IF (attr_name == "NGZF") THEN
                      CALL extractDataContent(sub_tmp_node,obj%ngzf)
                   END IF
                END IF
             END DO
          ELSE IF(attr_name == "electronic exchange-correlation") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"i")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                IF(hasAttribute(sub_tmp_node, "name")) THEN
                   CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
                   IF (attr_name == "GGA") THEN
                      CALL extractDataContent(sub_tmp_node,obj%gga)
                   ELSE IF (attr_name == "ALDAX") THEN
                      CALL extractDataContent(sub_tmp_node,obj%aldax)
                   ELSE IF (attr_name == "AGGAX") THEN
                      CALL extractDataContent(sub_tmp_node,obj%aggax)
                   ELSE IF (attr_name == "ALDAC") THEN
                      CALL extractDataContent(sub_tmp_node,obj%aldac)
                   ELSE IF (attr_name == "AGGAC") THEN
                      CALL extractDataContent(sub_tmp_node,obj%aggac)
                   END IF
                END IF
             END DO
          ELSE IF(attr_name == "vdW DFT") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node,"i")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                IF(hasAttribute(sub_tmp_node, "name")) THEN
                   CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
                   IF (attr_name == "LUSE_VDW") THEN
                      CALL extractDataContent(sub_tmp_node,logical_)
                      IF(TRIM(logical_)=="F") THEN
                         obj%luse_vdw = .FALSE.
                      ELSE IF(TRIM(logical_)=="T") THEN
                         obj%luse_vdw = .TRUE.
                      ELSE
                         CALL errore("vasp_read_parameters","LUSE_VDW: wrong occurrence",10)
                      END IF
                   ELSE IF (attr_name == "Zab_VDW") THEN
                      CALL extractDataContent(sub_tmp_node,obj%zab_vdw)
                   ELSE IF (attr_name == "PARAM1") THEN 
                      CALL extractDataContent(sub_tmp_node,obj%param1) 
                   ELSE IF (attr_name == "PARAM2") THEN 
                      CALL extractDataContent(sub_tmp_node,obj%param2) 
                   ELSE IF (attr_name == "PARAM3") THEN 
                      CALL extractDataContent(sub_tmp_node,obj%param3) 
                   END IF 
                END IF
             END DO
          END IF
       END IF
    END DO
    !
    RETURN
    !
  END SUBROUTINE vasp_read_parameters
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  SUBROUTINE vasp_read_structure(xml_node, obj, ierr)
    !
    USE FoX_dom,              ONLY : item, getElementsByTagname, nodeList, Node
    USE FoX_dom,              ONLY : getTagName, getLength, extractDataContent
    USE FoX_dom,              ONLY : hasAttribute, getAttributes, extractDataAttribute
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER          :: xml_node
    TYPE(vasp_structure_type), INTENT(OUT)   :: obj
    INTEGER, OPTIONAL, INTENT(OUT)           :: ierr
    !
    TYPE(Node), POINTER        :: tmp_node, sub_tmp_node, ssub_tmp_node, sssub_tmp_node
    TYPE(NodeList),  POINTER   :: tmp_node_list, sub_tmp_node_list, ssub_tmp_node_list, sssub_tmp_node_list
    INTEGER :: tmp_node_list_size, sub_tmp_node_list_size, ssub_tmp_node_list_size, sssub_tmp_node_list_size
    INTEGER :: index, iostat_
    INTEGER :: iitem, i_sitem, i_ssitem
    CHARACTER(LEN=100) :: attr_name
    !
    obj%tagname = getTagName(xml_node)
    !
    tmp_node_list => getElementsByTagname(xml_node, "crystal")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
       IF (PRESENT(ierr)) THEN
          CALL infomsg("vasp_read_structure","crystal: wrong number of occurrences")
          ierr = ierr + 1
       ELSE
          CALL errore("vasp_read_structure","crystal: wrong number of occurrences",10)
       END IF
    END IF
    !
    tmp_node => item(tmp_node_list,0)
    sub_tmp_node_list => getElementsByTagname(tmp_node, "varray")
    sub_tmp_node_list_size = getLength(sub_tmp_node_list)
    DO i_sitem = 0, sub_tmp_node_list_size-1
       sub_tmp_node => item(sub_tmp_node_list,i_sitem)
       IF (hasAttribute(sub_tmp_node, "name")) THEN
          CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
          IF (attr_name == "basis") THEN
             ssub_tmp_node_list => getElementsByTagname(sub_tmp_node, "v")
             ssub_tmp_node_list_size = getLength(ssub_tmp_node_list)
             IF (ssub_tmp_node_list_size /= 3) THEN
                CALL errore("vasp_read_structure","basis: wrong number of occurrences",10)
             END IF
             DO i_ssitem = 0, ssub_tmp_node_list_size-1
                ssub_tmp_node => item(ssub_tmp_node_list,i_ssitem)
                CALL extractDataContent(ssub_tmp_node,obj%at(1:3,i_ssitem+1))
             END DO
          ELSE IF (attr_name == "rec_basis") THEN
             ssub_tmp_node_list => getElementsByTagname(sub_tmp_node, "v")
             ssub_tmp_node_list_size = getLength(ssub_tmp_node_list)
             IF (ssub_tmp_node_list_size /= 3) THEN
                CALL errore("vasp_read_structure","rec_basis: wrong number of occurrences",10)
             END IF
             DO i_ssitem = 0, ssub_tmp_node_list_size-1
                ssub_tmp_node => item(ssub_tmp_node_list,i_ssitem)
                CALL extractDataContent(ssub_tmp_node,obj%bg(1:3,i_ssitem+1))
             END DO
          END IF
       END IF
    END DO
    sub_tmp_node_list => getElementsByTagname(tmp_node, "i")
    sub_tmp_node_list_size = getLength(sub_tmp_node_list)
    DO i_sitem = 0, sub_tmp_node_list_size-1
       sub_tmp_node => item(sub_tmp_node_list,i_sitem)
       IF (hasAttribute(sub_tmp_node,"name")) THEN
          CALL extractDataAttribute(sub_tmp_node,"name",attr_name)
          IF(attr_name=="volume") THEN
             CALL extractDataContent(sub_tmp_node,obj%volume)
          END IF
       END IF
    END DO
    ! 
    tmp_node_list => getElementsByTagname(xml_node, "varray")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    DO iitem = 0, tmp_node_list_size-1
       tmp_node => item(tmp_node_list,iitem)
       IF(hasAttribute(tmp_node,"name")) THEN
          CALL extractDataAttribute(tmp_node,"name",attr_name)
          IF(attr_name=="positions") THEN
             sub_tmp_node_list => getElementsByTagname(tmp_node, "v")
             sub_tmp_node_list_size = getLength(sub_tmp_node_list)
             obj%nat=sub_tmp_node_list_size
             ALLOCATE(obj%tau(1:3,obj%nat))
             DO i_sitem = 0, sub_tmp_node_list_size-1
                sub_tmp_node => item(sub_tmp_node_list,i_sitem)
                CALL extractDataContent(sub_tmp_node,obj%tau(1:3,i_sitem+1))
             END DO
          END IF
       END IF
    END DO
    
    
  END SUBROUTINE vasp_read_structure
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
SUBROUTINE vasp_init_xc(vasp_parameters,vasp_atominfo,iexch,icorr,igcx,igcc,inlc,ierr)
  !---------------------------------------------------------
  USE constants,            ONLY : eps4
  USE vdW_DF,               ONLY : vdw_type
  IMPLICIT NONE
  !
  TYPE(vasp_parameters_type), INTENT(IN)      :: vasp_parameters
  TYPE(vasp_atominfo_type), INTENT(IN)        :: vasp_atominfo
  INTEGER,INTENT(OUT)   :: iexch, icorr, igcx, igcc, inlc, ierr
  !
  ierr = 0
  iexch = -1; icorr = -1; igcx = -1; igcc = -1
  !
  ! set xc from POTCAR
  !
  IF (vasp_atominfo%pseudo=='PAW') THEN
     iexch = 1; icorr = 1; igcx = 0; igcc = 0
  ELSE IF (vasp_atominfo%pseudo=='PAW_GGA') THEN
     iexch = 1; icorr = 4; igcx = 2; igcc = 2
  ELSE IF (vasp_atominfo%pseudo=='PAW_PBE') THEN
     iexch = 1; icorr = 4; igcx = 3; igcc = 4
  ELSE IF (vasp_atominfo%pseudo=='US') THEN
     iexch = 1; icorr = 4; igcx = 2; igcc = 2
  ELSE 
     CALL errore ("vasp_init_xc", "unknown potential file", 1)
  END IF
  !
  ! overwrite xc if GGA is set
  !
  IF(vasp_parameters%gga=='CA') THEN
     iexch = 1; icorr = 1; igcx = 0; igcc = 0
  ELSE IF(vasp_parameters%gga=='91') THEN
     iexch = 1; icorr = 4; igcx = 2; igcc = 2
  ELSE IF(vasp_parameters%gga=='PE') THEN
     iexch = 1; icorr = 4; igcx = 3; igcc = 4
  ELSE IF(vasp_parameters%gga=='CX') THEN
     iexch = 1; icorr = 4; igcx = 27
  ELSE IF(vasp_parameters%gga=='RE') THEN
     iexch = 1; icorr = 4; igcx = 4; igcc = 4
  ELSE IF (vasp_parameters%gga/='--') THEN
     CALL errore ("vasp_init_xc", "GGA type not implemented", 1)
  ENDIF
  !
  !
  !
  IF(ABS(vasp_parameters%aldax) .LT. eps4) THEN
     iexch = 0
  ELSE IF(ABS(vasp_parameters%aldax-1._DP) .GT. eps4) THEN
     CALL errore ("vasp_init_xc", "hybrid calculations not implemented", 1)
  ENDIF
  !
  IF(ABS(vasp_parameters%aldac) .LT. eps4) THEN
     icorr = 0
  ELSE IF(ABS(vasp_parameters%aldac-1._DP) .GT. eps4) THEN
     CALL errore ("vasp_init_xc", "hybrid calculations not implemented", 1)
  ENDIF
  !
  IF(ABS(vasp_parameters%aggax) .LT. eps4) THEN
     igcx = 0
  ELSE IF(ABS(vasp_parameters%aggax-1._DP) .GT. eps4) THEN
     CALL errore ("vasp_init_xc", "hybrid calculations not implemented", 1)
  ENDIF
  !
  IF(ABS(vasp_parameters%aggac) .LT. eps4) THEN
     igcc = 0
  ELSE IF(ABS(vasp_parameters%aggac-1._DP) .GT. eps4) THEN
     CALL errore ("vasp_init_xc", "hybrid calculations not implemented", 1)
  ENDIF
  !
  !  set vdW_DF
  !
  IF(vasp_parameters%luse_vdw) THEN
     IF(ABS(vasp_parameters%zab_vdw-(-0.8491))<eps4) THEN
        vdw_type=1
        inlc = 1
     ELSEIF(ABS(vasp_parameters%zab_vdw-(-1.887))<eps4) THEN
        vdw_type=2
        inlc = 2
     ELSE
        CALL errore ('vasp_init_xc', 'Zab_vdW not implemented', vasp_parameters%zab_vdw)
     END IF
  ELSE
     inlc = 0
  END IF
  !
  RETURN
  ! 
END SUBROUTINE vasp_init_xc
!---------------------------------------------------------
SUBROUTINE vasp_init_vars_from_schema( what, ierr, vasp_kpoints, vasp_parameters, vasp_atominfo, vasp_structure)
  !---------------------------------------------------------
  USE scf,                 ONLY: rho
  USE lsda_mod,            ONLY: nspin
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN)                :: what
  TYPE(vasp_kpoints_type), INTENT(IN)         :: vasp_kpoints
  TYPE(vasp_parameters_type), INTENT(IN)      :: vasp_parameters
  TYPE(vasp_atominfo_type), INTENT(IN)        :: vasp_atominfo
  TYPE(vasp_structure_type), INTENT(IN)       :: vasp_structure
  INTEGER,INTENT (OUT)                        :: ierr
  !
  LOGICAL                :: ldim
  LOGICAL                :: latom
  LOGICAL                :: lkpoint
  LOGICAL                :: found

  !
  ierr = 0
  !
  ldim             = .FALSE.
  latom            = .FALSE.
  lkpoint          = .FALSE.
  found            = .FALSE.
  !
  SELECT CASE ( what )
  CASE( 'dim' )
     !
     ldim = .TRUE.
     !
  CASE( 'atom' )
     !
     latom = .TRUE.
     !
  CASE( 'kpoint' )
     !
     lkpoint = .TRUE.
     !
  END SELECT
  !
  !
  IF (ldim) THEN
     CALL vasp_readschema_dim( vasp_kpoints, vasp_parameters, vasp_atominfo, vasp_structure )
  END IF
  !
  IF (latom) THEN
     CALL vasp_readschema_atom( vasp_atominfo, vasp_structure )
  END IF
  !
  IF (lkpoint) THEN
     CALL vasp_readschema_kpoint( vasp_kpoints )
  END IF
  !
  RETURN
  !
END SUBROUTINE vasp_init_vars_from_schema
!---------------------------------------------------------
  SUBROUTINE vasp_readschema_dim( vasp_kpoints, vasp_parameters, vasp_atominfo, vasp_structure ) 
  !---------------------------------------------------------
    USE constants,        ONLY : RYTOEV
    USE constants,        ONLY : e2
    USE cell_base,        ONLY : at, bg, alat, omega, cell_base_init
    USE ions_base,        ONLY : nat, nsp
    USE symm_base,        ONLY : nsym
    USE gvect,            ONLY : ngm_g, ecutrho
    USE fft_base,         ONLY : dfftp
    USE gvecs,            ONLY : ngms_g, dual
    USE fft_base,         ONLY : dffts
    USE lsda_mod,         ONLY : nspin, lsda
    USE noncollin_module, ONLY : noncolin
    USE klist,            ONLY : nkstot, nelec
    USE wvfct,            ONLY : nbnd, npwx
    USE gvecw,            ONLY : ecutwfc
    USE control_flags,    ONLY : gamma_only
    USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                 nproc_image_file, ntask_groups_file, &
                                 nproc_bgrp_file, nproc_ortho_file
    ! 
    IMPLICIT NONE
    !
    TYPE (vasp_kpoints_type), INTENT(IN)     :: vasp_kpoints
    TYPE (vasp_parameters_type), INTENT(IN)  :: vasp_parameters
    TYPE (vasp_atominfo_type), INTENT(IN)    :: vasp_atominfo
    TYPE (vasp_structure_type), INTENT(IN)   :: vasp_structure
    !
    INTEGER                                  :: npwx_
    REAL(DP)                                 :: celldm_(6)
    !
    !---------------------------------------------------------
    !           ENERGY CUTOFF
    !---------------------------------------------------------
    ecutwfc=vasp_parameters%enmax/RYTOEV 
    dual = 4._DP
    !---------------------------------------------------------
    !           SPIN
    !---------------------------------------------------------
    nspin = vasp_parameters%ispin
    IF(nspin == 1) THEN
       lsda = .FALSE.  
    ELSE IF(nspin == 2) THEN 
       lsda = .TRUE.
    END IF
    noncolin = vasp_parameters%lnoncollinear
    IF (noncolin) THEN
       CALL errore ("vasp_readschema_dim", "noncollinear calculations not implemented", 1)
    END IF
    !---------------------------------------------------------
    !           BANDS
    !---------------------------------------------------------
    nbnd = vasp_parameters%nbands
    !---------------------------------------------------------
    !           ATOMS AND SPECIES
    !---------------------------------------------------------
    IF(vasp_atominfo%nat == vasp_structure%nat) THEN
       nat = vasp_atominfo%nat
    ELSE
       CALL errore ("vasp_readschema_dim", "wrong atom coordinate length", 1)
    END IF
    nsp = vasp_atominfo%nsp
    !---------------------------------------------------------
    !           CRYSTAL
    !---------------------------------------------------------
    ! set alat, omega, tpiba2 from at 
    at=vasp_structure%at
    celldm_=0._DP
    CALL cell_base_init( 0, celldm_, 0._DP, 0._DP, 0._DP, 0._DP, 0._DP, &
                         0._DP, .true., TRANSPOSE(at), 'angstrom')
    
    !---------------------------------------------------------
    !           KPOINTS
    !---------------------------------------------------------
    nkstot = vasp_kpoints%nk
    !
    !
    IF (lsda) THEN
       nkstot = nkstot*2
       nbnd = nbnd / 2
    END IF

    
  END SUBROUTINE vasp_readschema_dim
  !---------------------------------------------------------
  !---------------------------------------------------------
  SUBROUTINE vasp_readschema_atom( vasp_atominfo, vasp_structure) 
  !---------------------------------------------------------
    USE constants,        ONLY : e2, ANGSTROM_AU 
    USE cell_base,        ONLY : at, alat, omega
    USE ions_base,        ONLY : nat, nsp, ityp, tau, zv, atm
    USE symm_base,        ONLY : nsym
    USE gvect,            ONLY : ngm_g, ecutrho
    USE fft_base,         ONLY : dfftp
    USE gvecs,            ONLY : ngms_g, dual
    USE fft_base,         ONLY : dffts
    USE lsda_mod,         ONLY : lsda
    USE noncollin_module, ONLY : noncolin
    USE klist,            ONLY : nkstot, nelec
    USE wvfct,            ONLY : nbnd, npwx
    USE gvecw,            ONLY : ecutwfc
    USE control_flags,    ONLY : gamma_only
    USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                 nproc_image_file, ntask_groups_file, &
                                 nproc_bgrp_file, nproc_ortho_file
    ! 
    IMPLICIT NONE
    !
    TYPE (vasp_atominfo_type), INTENT(IN)    :: vasp_atominfo
    TYPE (vasp_structure_type), INTENT(IN)   :: vasp_structure
    !
    INTEGER                                  :: npwx_
    INTEGER                                  :: ia, i
    !
    !---------------------------------------------------------
    !           ATOMS AND SPECIES
    !---------------------------------------------------------
    atm(1:nsp) = vasp_atominfo%atm(1:nsp)
    ityp(1:nat) = vasp_atominfo%ityp(1:nat)
    tau(1:3,1:nat) = vasp_structure%tau(1:3,1:nat)
    zv(1:nsp) = vasp_atominfo%zv(1:nsp)
    DO ia = 1, nat
    !
       DO i = 1, 3
       !
       tau(i,ia) = at(i,1) * vasp_structure%tau(1,ia) + &
                   at(i,2) * vasp_structure%tau(2,ia) + &
                   at(i,3) * vasp_structure%tau(3,ia)
       !
       END DO
    ! 
    END DO
    ! 
  END SUBROUTINE vasp_readschema_atom
  !---------------------------------------------------------
  !---------------------------------------------------------
  SUBROUTINE vasp_readschema_kpoint( vasp_kpoints ) 
  !---------------------------------------------------------
    USE constants,        ONLY : e2
    USE cell_base,        ONLY : at, alat, omega
    USE ions_base,        ONLY : nat, nsp, ityp, tau, atm
    USE symm_base,        ONLY : nsym
    USE gvect,            ONLY : ngm_g, ecutrho
    USE fft_base,         ONLY : dfftp
    USE gvecs,            ONLY : ngms_g, dual
    USE fft_base,         ONLY : dffts
    USE lsda_mod,         ONLY : lsda, isk
    USE noncollin_module, ONLY : noncolin
    USE klist,            ONLY : nkstot, nks, xk, wk
    USE wvfct,            ONLY : nbnd, npwx
    USE gvecw,            ONLY : ecutwfc
    USE control_flags,    ONLY : gamma_only
    USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                 nproc_image_file, ntask_groups_file, &
                                 nproc_bgrp_file, nproc_ortho_file
    ! 
    IMPLICIT NONE
    !
    TYPE (vasp_kpoints_type), INTENT(IN)    :: vasp_kpoints
    !
    INTEGER                 :: ik, num_k_points
    !
    !---------------------------------------------------------
    !           KPOINTS
    !---------------------------------------------------------
    num_k_points = vasp_kpoints%nk
    xk(1:3,1:num_k_points) = vasp_kpoints%xk(1:3,1:num_k_points)
    wk(1:num_k_points) = vasp_kpoints%wk(1:num_k_points)
    DO ik = 1, num_k_points
       !
       isk = 1
       IF ( lsda ) THEN
          !
          xk(:,ik+num_k_points) = xk(:,ik)
          !
          wk(ik+num_k_points) = wk(ik)
          !
          isk(ik+num_k_points) = 2
          !
       END IF
       !
    END DO
    !
  END SUBROUTINE vasp_readschema_kpoint
  !---------------------------------------------------------
      !------------------------------------------------------------------------
      SUBROUTINE set_dimensions()
        !------------------------------------------------------------------------
        !
        USE constants, ONLY : pi, eps8
        USE cell_base, ONLY : alat, tpiba, tpiba2
        USE gvect,     ONLY : ecutrho, gcutm
        USE gvecs,     ONLY : gcutms, dual, doublegrid
        USE gvecw,     ONLY : gcutw, ecutwfc
        !
        !
        ! ... Set the units in real and reciprocal space
        !
        ! tpiba  = 2.D0 * pi / alat
        ! tpiba2 = tpiba**2
        !
        ! ... Compute the cut-off of the G vectors
        !
        gcutw =        ecutwfc / tpiba2
        gcutm = dual * ecutwfc / tpiba2
        ecutrho=dual * ecutwfc
        !
        doublegrid = ( dual > 4.D0 + eps8 )
        IF ( doublegrid ) THEN
           gcutms = 4.D0 * ecutwfc / tpiba2
        ELSE
           gcutms = gcutm
        END IF
        !
      END SUBROUTINE set_dimensions
      !
      SUBROUTINE fft_type_setdim( desc, nr1, nr2, nr3 )
         USE fft_types,          ONLY: fft_type_descriptor
         USE fft_support,        ONLY: good_fft_order, good_fft_dimension
         TYPE (fft_type_descriptor) :: desc
         INTEGER, INTENT(IN) :: nr1, nr2, nr3
         IF (desc%nr1 /= 0 .OR. desc%nr1 /= 0 .OR. desc%nr1 /= 0 ) &
            CALL fftx_error__(' fft_type_setdim ', ' fft dimensions already set ', 1)
         desc%nr1 = nr1
         desc%nr2 = nr2
         desc%nr3 = nr3
         desc%nr1 = good_fft_order( desc%nr1 )
         desc%nr2 = good_fft_order( desc%nr2 )
         desc%nr3 = good_fft_order( desc%nr3 )
         desc%nr1x  = good_fft_dimension( desc%nr1 )
         desc%nr2x  = desc%nr2
         desc%nr3x  = good_fft_dimension( desc%nr3 )
      END SUBROUTINE fft_type_setdim
      !
END MODULE vasp_xml
