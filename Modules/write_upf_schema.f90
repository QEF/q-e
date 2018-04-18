!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
MODULE write_upf_schema_module
  !----------------------------------------------------------------------------=!
  !  this module handles the writing of pseudopotential data
  ! ...   declare modules
  USE kinds,        ONLY: DP
  USE pseudo_types, ONLY: pseudo_upf, pseudo_config, deallocate_pseudo_config
  USE radial_grids, ONLY: radial_grid_type
  USE global_version, ONLY: version_number, svn_revision 
  USE Fox_wxml
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: write_upf_schema 

CONTAINS

  !-------------------------------+
  SUBROUTINE write_upf_schema(xf, upf, conf, u_input)
    !----------------------------+
    ! Write pseudopotential in UPF format version 2, uses iotk
    !
    IMPLICIT NONE
    TYPE(xmlf_t),INTENT(INOUT)  :: xf   ! xmlfile for writing
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    ! optional: configuration used to generate the pseudopotential
    TYPE(pseudo_config),OPTIONAL,INTENT(IN) :: conf
    ! optional: unit pointing to input file containing generation data
    INTEGER, OPTIONAL, INTENT(IN):: u_input
    INTEGER                      :: irow
    CHARACTER(LEN=*),PARAMETER   :: QE_PP_URI="http://www.quantum-espresso.org/ns/qes/qe_pp-1.0",&
                                    QE_PP_LOC=QE_PP_URI//"  "//&
                                              "http://www.quantum-espresso.org/ns/qes/qe_pp-1.0.xsd",&
                                    XSD_VERSION = "QE_PP-1.0"
    !
    ! Initialize the file
    CALL xml_DeclareNameSpace(XF = xf, PREFIX = "xsi", nsUri="http://www.w3.org/2001/XMLSchema-instance") 
    CALL xml_DeclareNameSpace(XF = xf, PREFIX = "qe_pp", nsUri=QE_PP_URI)
    CALL xml_NewElement(XF=xf, name="qe_pp:pseudo")
    CALL xml_addAttribute(XF=xf, name="xsi:schemalocation", VALUE = QE_PP_LOC)
    !
    CALL xml_NewElement(xf,"xsd_version")
    CALL xml_AddCharacters(xf, XSD_VERSION)
    CALL xml_endElement(xf, "xsd_version")
    !
    !
    CALL write_info(xf, upf, conf, u_input)
    ! Write machine-readable header
    CALL write_header(xf, upf)
    ! Write radial grid mesh
    CALL write_mesh(xf, upf)
    ! Write non-linear core correction charge
    IF(upf%nlcc) THEN 
      CALL xml_newElement(xf, 'pp_nlcc')
         CALL xml_addAttribute(xf, 'size', upf%mesh)
         DO irow = 1, upf%mesh, 4
            CALL xml_addNewLine(xf)
            CALL xml_addCharacters(xf, upf%rho_atc(irow:min(irow-1+4,upf%mesh)) , fmt = 's16')
         END DO
         CALL xml_addNewLine(xf)
      CALL xml_endElement(xf, 'pp_nlcc') 
    END IF
    ! Write local potential
    IF(.not. upf%tcoulombp) THEN
       CALL xml_newElement(xf, 'pp_local')
          CALL xml_addAttribute(xf,'size', upf%mesh)
          DO irow = 1, upf%mesh, 4
             CALL xml_addNewLine(xf)
             CALL xml_addCharacters(xf, upf%vloc(irow:min(irow-1+4,upf%mesh)) , fmt = 's16')
          END DO
          CALL xml_addNewLine(xf)
       CALL xml_endElement(xf, 'pp_local')    
    ENDIF
    ! Write potentials in semilocal form (if existing)
    IF ( upf%typ == "SL" ) CALL write_semilocal(xf, upf)
    ! Write nonlocal components: projectors, augmentation, hamiltonian elements
    CALL write_nonlocal(xf, upf)
    ! Write initial pseudo wavefunctions
    ! (usually only wfcs with occupancy > 0)
    CALL write_pswfc(xf, upf)
    ! If included, write all-electron and pseudo wavefunctions
    CALL write_full_wfc(xf, upf)
    ! Write valence atomic density (used for initial density)
    CALL xml_NewElement(xf, 'pp_rhoatom') 
       CALL xml_addAttribute(xf,'size', upf%mesh) 
       DO irow = 1, upf%mesh, 4
          CALL xml_addNewLine(xf)
          CALL xml_addCharacters(xf, upf%rho_at(irow:min(irow-1+4,upf%mesh)) , fmt = 's16')
       END DO
       CALL xml_addNewLine(xf)
    CALL xml_endElement(xf, 'pp_rhoatom')
    ! Write additional data for PAW (All-electron charge, wavefunctions, vloc..)
    CALL write_paw(xf, upf)
    ! Write additional data for GIPAW reconstruction
    CALL write_gipaw(xf, upf)
    !
    ! Close the file (not the unit!)
    CALL xml_endElement(xf,'qe_pp:pseudo')
    CALL xml_Close(xf)

  CONTAINS
    !
    SUBROUTINE write_info(u, upf, conf, u_input)
      ! Write human-readable header
      ! The header is written directly, not via iotk
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)  :: u    ! write to xml file u
      TYPE(pseudo_upf),INTENT(IN) :: upf  ! the pseudo data
      ! optional: configuration used to generate the pseudopotential
      TYPE(pseudo_config),OPTIONAL,INTENT(IN) :: conf
      INTEGER, OPTIONAL, INTENT(IN) :: u_input ! read input data from u_input
      !
      INTEGER :: nb ! aux counter
      INTEGER :: ierr ! /= 0 if something went wrong
      CHARACTER(LEN=4096) :: char_buf
      CHARACTER(len=256) :: line
      LOGICAL :: opnd
      !
      CALL xml_NewElement( u, "pp_info") 
        CALL xml_NewElement(u,"generated")
          CALL xml_addCharacters(u, TRIM(upf%generated), parsed = .FALSE. )
        CALL xml_EndElement(u,"generated")
        ! 
        CALL xml_NewElement(u,"creator")
          CALL xml_addAttribute(u,name="NAME",value="QE Atomic Code")
          CALL xml_addAttribute(u,name= "VERSION", value = version_number // ':'//svn_revision) 
          CALL xml_addCharacters(u,TRIM(upf%author))
        CALL xml_EndElement(u, 'creator')
        !
        CALL xml_NewElement(u,"created")
          CALL xml_addAttribute(u, name="DATE", VALUE=TRIM(upf%date))
        CALL xml_endElement(u,'created')
        IF ( PRESENT(u_input) ) THEN
         !
         ! copy content of input file used in pseudopotential generation
         !
         INQUIRE (unit=u_input, opened=opnd)
         IF (opnd) THEN
            CALL xml_NewElement(u,"input")
            CALL xml_addAttribute(u, name = "program", value = "ld1.x") 
            
            REWIND (unit=u_input)
            char_buf=""
            read_write_loop: DO
               READ (u_input, '(A)',end=20,err=25) line
               char_buf = TRIM(char_buf) // new_line('a') // trim(line)
               CYCLE read_write_loop
25             CALL infomsg('write_upf_schema::write_inputfile', 'problem writing input data')
20             EXIT read_write_loop
            END DO read_write_loop
            char_buf = TRIM(char_buf)// new_line('a')
            CALL xml_addCharacters(u, TRIM(char_buf), parsed = .FALSE.)
            CALL xml_endElement(u, 'input')
         ELSE
            CALL infomsg('write_upf_v2::write_inputfile', 'input file not open')
         END IF
         !
      END IF
 
        
        CALL xml_NewElement(u,'type')
          CALL xml_addCharacters(u,TRIM(upf%typ))
        CALL xml_EndElement(u,'type')
        CALL xml_NewElement(u,'relativistic_effects')
        IF (TRIM(upf%rel)=='full') THEN
           CALL xml_addCharacters(u,'full')
        ELSE IF (TRIM(upf%rel)=='scalar') THEN
           CALL xml_addCharacters(u,'scalar')
        ELSE
          CALL xml_addCharacters(u,'none')
        ENDIF
        CALL xml_EndElement(u, 'relativistic_effects')
        CALL xml_NewElement(u,'element')
          CALL xml_addCharacters(u,TRIM(upf%psd))
        CALL xml_endElement(u,'element')
        CALL xml_NewElement(u,'functional')
          CALL xml_addCharacters(u,TRIM(upf%dft))
        CALL xml_endElement(u,'functional')
        CALL xml_newElement(u,'suggested_basis')
          CALL xml_addAttribute(u,name='ecutwfc',value = upf%ecutwfc)
          IF (upf%tpawp .OR. upf%tvanp ) CALL xml_addAttribute(u,name='ecutrho',value=upf%ecutrho)
        CALL xml_endElement(u,'suggested_basis')
        DO nb =1, upf%nwfc
           IF( upf%oc(nb) >= 0._dp) THEN 
             CALL xml_newElement(u, name = "valence_orbital")
                CALL xml_addAttribute(u, name='nl', value = upf%els(nb))
                CALL xml_addAttribute(u, name = 'pn', value = upf%nchi(nb) )
                CALL xml_addAttribute(u, name = 'l', value = upf%lchi(nb) )
                !
                CALL xml_newElement (u, name = "occupation")
                   CALL xml_addCharacters(u, upf%oc(nb) )
                CALL xml_endElement(u, "occupation") 
                CALL xml_newElement(u, "Rcut") 
                   CALL xml_addCharacters(u, upf%rcut_chi(nb))
                CALL xml_endElement(u, "Rcut")
                IF (upf%rcutus_chi(nb) > 0.d0) THEN
                   CALL xml_newElement(u, "RcutUS")
                      CALL xml_addCharacters(u, upf%rcutus_chi(nb))
                   CALL xml_endElement(u, "RcutUS")
                END IF
                CALL xml_newElement(u, "Epseu")
                   CALL xml_addCharacters(u, upf%epseu(nb))
                CALL xml_endElement(u,"Epseu")
             CALL xml_endElement(u,"valence_orbital")
          END IF
       END DO
             
      
      IF( present(conf) ) THEN
         CALL xml_newElement(u, "generation_configuration")
         char_buf = ""
         DO nb = 1,conf%nwfs
            WRITE(line, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
                 conf%els(nb), conf%nns(nb), &
                 conf%lls(nb), conf%ocs(nb), conf%rcut(nb), &
                 conf%rcutus(nb), conf%enls(nb)
            char_buf = TRIM(char_buf) // new_line('a') // TRIM(line)
         ENDDO
         WRITE(line,'(4x,2a)') 'Pseudization used: ',TRIM(conf%pseud)
         char_buf = TRIM(char_buf) // new_line('a') // TRIM(line) // new_line('a') 
         CALL xml_addCharacters(u, TRIM(char_buf), parsed = .FALSE.)
         CALL xml_endElement(u,"generation_configuration")
      ENDIF

      IF(TRIM(upf%comment) /= ' ') THEN
         CALL xml_addComment (u, TRIM(upf%comment), WS_SIGNIFICANT=.TRUE. )
      END IF
      !
      !
      CALL xml_endElement(u,"pp_info")
      !
      RETURN
100   CALL errore('write_upf_schema::write_info', 'Writing pseudo file', 1)
      !
    END SUBROUTINE write_info
    !
    !
    SUBROUTINE write_header(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t), INTENT(INOUT)  :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      !
      INTEGER :: nw
      !
      ! Write HEADER section with some initialization data
      !
      CALL xml_newElement(u, 'pp_header')
         CALL xml_newElement(u,'element')
            print *, (upf%psd)
            CALL xml_addCharacters(u, TRIM(upf%psd))         
         CALL xml_endElement(u,'element')                    
         !                                                   
         CALL xml_newElement(u,'z_valence')                  
            CALL xml_addCharacters(u, upf%zp)                
         CALL xml_endElement(u,'z_valence')                  
         !                                                   
         CALL xml_newElement(u, 'type')                      
            CALL xml_addCharacters(u, TRIM(upf%typ))         
         CALL xml_endElement(u,'type')                       
         CALL xml_newElement(u,'functional')                 
            CALL xml_addCharacters(u, TRIM(upf%dft))         
         CALL xml_endElement(u,'functional')                 
         !                                                   
         CALL xml_newElement(u, 'relativistic')              
            CALL xml_addCharacters(u, TRIM(upf%rel))         
         CALL xml_endElement(u,'relativistic')              
         !                                                   
         CALL xml_newElement(u,'is_ultrasoft')               
            CALL xml_addCharacters(u, upf%tvanp)             
         CALL xml_endElement(u,'is_ultrasoft')               
         !                                                   
         CALL xml_newElement(u,'is_paw')                     
            CALL xml_addCharacters(u, upf%tpawp)              
         CALL xml_endElement(u,'is_paw')                     
         !                                                   
         CALL xml_newElement(u,'is_coulomb')                 
            CALL xml_addCharacters(u, upf%tcoulombp)         
         CALL xml_endElement(u,'is_coulomb')                 
         !                                                   
         CALL xml_newElement(u,'has_so')                     
            CALL xml_addCharacters(u, upf%has_so)            
         CALL xml_endElement(u,'has_so')                     
         !                                                   
         CALL xml_newElement(u,'has_wfc')                    
            CALL xml_addCharacters(u, upf%has_wfc)           
         CALL xml_endElement(u,'has_wfc')                    
         !                                                   
         CALL xml_newElement(u,'has_gipaw')                  
            CALL xml_addCharacters(u, upf%has_gipaw)         
         CALL xml_endElement(u,'has_gipaw')                  
         !                                                   
         !Emine                                              
         CALL xml_newElement(u,'paw_as_gipaw')               
            CALL xml_addCharacters(u, upf%paw_as_gipaw)      
         CALL xml_endElement(u,'paw_as_gipaw')               
         !                                                   
         CALL xml_newElement(u,'core_correction')            
            CALL xml_addCharacters(u, upf%nlcc)              
         CALL xml_endElement(u,'core_correction')            
         !                                                   
         CALL xml_newElement(u,'total_psenergy')             
            CALL xml_addCharacters(u, upf%etotps)            
         CALL xml_endElement(u,'total_psenergy')             
         !                                                   
         CALL xml_newElement(u,'wfc_cutoff')                 
            CALL xml_addCharacters(u, upf%ecutwfc)           
         CALL xml_endElement(u,'wfc_cutoff')                 
         !                                                   
         CALL xml_newElement(u,'rho_cutoff')                 
            CALL xml_addCharacters(u, upf%ecutrho)           
         CALL xml_endElement(u,'rho_cutoff')                 
         !                                                   
         CALL xml_newElement(u,'l_max')                      
            CALL xml_addCharacters(u, upf%lmax)              
         CALL xml_endElement(u,'l_max')                      
         !                                                   
         CALL xml_newElement(u,'l_max_rho')                  
            CALL xml_addCharacters(u, upf%lmax_rho)          
         CALL xml_endElement(u,'l_max_rho')                  
         !                                                   
         CALL xml_newElement(u,'l_local')                    
            CALL xml_addCharacters(u, upf%lloc)              
         CALL xml_endElement(u,'l_local')                    
         !                                                   
         CALL xml_newElement(u,'mesh_size')                  
            CALL xml_addCharacters(u, upf%mesh)              
         CALL xml_endElement(u,'mesh_size')                  
         !                                                   
         CALL xml_newElement(u,'number_of_wfc')              
            CALL xml_addCharacters(u, upf%nwfc)              
         CALL xml_endElement(u,'number_of_wfc')              
         !                                                   
         CALL xml_newElement(u,'number_of_proj')             
            CALL xml_addCharacters(u, upf%nbeta)             
         CALL xml_endElement(u,'number_of_proj')             
      !
      CALL xml_endElement(u, 'pp_header') 
      RETURN
    END SUBROUTINE write_header
    !
    SUBROUTINE write_mesh(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t), INTENT(INOUT)  :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong

      !
      !
      CALL xml_NewElement(u, 'pp_mesh')
      IF ( upf%dx .GT. 0.d0)   CALL xml_addAttribute(u, name='dx', value = upf%dx , fmt = 's16')
      IF ( upf%mesh .GT. 0 )   CALL xml_addAttribute(u, name='mesh', value = upf%mesh)
      IF (upf%dx  .GT. 0.d0) CALL xml_addAttribute(u, name='xmin', value = upf%xmin,  fmt = 's16')
      IF (upf%rmax  .GT. 0.d0) CALL xml_addAttribute(u, name='rmax', value = upf%rmax,  fmt = 's16')
      IF (upf%zmesh .GT.0.d0 )CALL xml_addAttribute(u, name='zmesh',value = upf%zmesh, fmt = 's16')

      CALL xml_NewElement(u, 'pp_r' ) 
         DO irow =1, upf%mesh, 4
            CALL xml_addNewLine(u)
            CALL xml_addCharacters(u, upf%r(irow:min(irow-1+4,upf%mesh) ) , fmt='s16')
         END DO
         CALL xml_addNewLine(xf)
      CALL xml_endElement(u,'pp_r')
      CALL xml_NewElement(u, 'pp_rab')
         DO irow = 1, upf%mesh, 4 
            CALL xml_addNewLine(u) 
            CALL xml_addCharacters(u, upf%rab(irow:min(irow-1+4, upf%mesh)), fmt = 's16')
         END DO
         CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_rab')
      !

      CALL xml_endElement(u, 'pp_mesh') 
      !
      RETURN
    END SUBROUTINE write_mesh
    !
    SUBROUTINE write_semilocal(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)   :: u    ! i/o unit
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      INTEGER :: nb, l, ind
      !
      CALL xml_newElement(u, 'pp_semilocal')
      !
      ! Write V_l(r)
      DO nb = 1,upf%nbeta
         l = upf%lll(nb)
         ind = 1
         CALL xml_newElement(u, 'vnl')
         CALL xml_addAttribute(u, name = 'l', value = l)
         IF ( upf%has_so ) THEN
            CALL xml_addAttribute(u, 'j', upf%jjj(nb))
            IF ( l > 0 .AND. ABS (upf%jjj(nb)-l-0.5_dp) < 0.001_dp) ind = 2
         ENDIF
         DO irow = 1, upf%mesh, 4
            CALL xml_addNewLine(u)
            CALL xml_addCharacters(u, upf%vnl(irow:min(irow-1+4, upf%mesh),l, ind), fmt = 's16')
         END DO
         CALL xml_addNewLine(xf)
         CALL xml_endElement(u, 'vnl')
      END DO
      !
      CALL xml_endElement(u, 'pp_semilocal')
      !
    END SUBROUTINE write_semilocal
    !
    SUBROUTINE write_nonlocal(u, upf)
      IMPLICIT NONE
      TYPE(xmlf_t),INTENT(INOUT)      :: u    ! xml file descriptor
      TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
      INTEGER :: ierr ! /= 0 if something went wrong
      !
      INTEGER :: nb,mb,ln,lm,l,nmb
      LOGICAL :: isnull
      REAL(DP),ALLOCATABLE :: tmp_dbuffer(:)
      !
      IF (upf%tcoulombp) RETURN
      !
      CALL xml_NewElement(u, 'pp_nonlocal') 
      !
      ! Write the projectors:
      DO nb = 1,upf%nbeta
         CALL xml_NewElement(u, 'pp_beta')
            CALL xml_addAttribute(u, 'index',                  nb )
            CALL xml_addAttribute(u, 'label',                  upf%els_beta(nb))
            CALL xml_addAttribute(u, 'angular_momentum',       upf%lll(nb))
            IF (upf%has_so) &
                          CALL xml_addAttribute(u, 'tot_ang_mom', upf%jjj(nb))
            CALL xml_addAttribute(u, 'cutoff_radius_index',    upf%kbeta(nb))
            CALL xml_addAttribute(u, 'cutoff_radius',          upf%rcut(nb), fmt = 's16')
            CALL xml_addAttribute(u, 'ultrasoft_cutoff_radius',upf%rcutus(nb), fmt = 's16')
            DO irow = 1, upf%mesh, 4
               CALL xml_addNewLine(u) 
               CALL xml_addCharacters(u, upf%beta(irow:min(irow-1+4,upf%mesh), nb), fmt = 's16')
            END DO
            CALL xml_addNewLine(xf)
         CALL xml_endElement(u, 'pp_beta') 
      ENDDO
      !
      ! Write the hamiltonian terms D_ij
      
      CALL xml_newElement(u, 'pp_dij') 
         CALL xml_addAttribute(u,'columns', upf%nbeta)
         CALL xml_addAttribute(u, 'rows', upf%nbeta )
         nb = upf%nbeta*upf%nbeta 
         ALLOCATE(tmp_dbuffer(nb))
         tmp_dbuffer = reshape(upf%dion,[nb])
         DO irow = 1, nb, 4
            CALL xml_addNewLine(u)
            CALL xml_addCharacters(u, tmp_dbuffer(irow:min(irow-1+4,nb)), fmt ='s16')
         END DO
         CALL xml_addNewLine(xf)
         DEALLOCATE(tmp_dbuffer)
      CALL xml_endElement(u, 'pp_dij')
      !
      ! Write the augmentation charge section
      augmentation : &
           IF(upf%tvanp .or. upf%tpawp) THEN
           CALL xml_newElement(u, 'pp_augmentation')
           CALL xml_newElement(u, 'nqf')
             CALL xml_addCharacters(u, upf%nqf)
           CALL xml_endElement(u, 'nqf')
           CALL xml_newElement(u, 'q_with_l') 
              CALL xml_addCharacters(u, upf%q_with_l)
           CALL xml_endElement(u, 'q_with_l')
           CALL xml_newElement(u, 'nqlc')
              CALL xml_addCharacters(u, upf%nqlc) 
           CALL xml_endElement(u, 'nqlc')            
           IF (upf%tpawp) THEN
              CALL xml_newElement(u, 'shape')
                 CALL xml_addCharacters(u, upf%paw%augshape) 
              CALL xml_endElement(u, 'shape')
              CALL xml_newElement(u, 'cutoff_r')
                 CALL xml_addCharacters(u, upf%paw%raug) 
              CALL xml_endElement(u, 'cutoff_r')
              CALL xml_newElement(u, 'cutoff_r_index')
                 CALL xml_addCharacters(u, upf%paw%iraug) 
              CALL xml_endElement(u, 'cutoff_r_index')
              CALL xml_newElement(u, 'l_max_aug')
                 CALL xml_addCharacters(u, upf%paw%lmax_aug) 
              CALL xml_endElement(u, 'l_max_aug')
              CALL xml_newElement(u, 'augmentation_epsilon')
                 CALL xml_addCharacters(u, upf%qqq_eps)
              CALL xml_endElement(u, 'augmentation_epsilon')
      ENDIF
      !
      ! Write the integrals of the Q functions
      CALL xml_newElement(u, 'pp_q')
          nb = upf%nbeta*upf%nbeta
          CALL xml_addAttribute(u, 'size',nb)
          ALLOCATE(tmp_dbuffer (nb))
          tmp_dbuffer = reshape(upf%qqq, [nb])
          DO irow =1, nb, 4
             CALL xml_addNewLine(u)
             CALL xml_addCharacters(u, tmp_dbuffer(irow:min(irow-1+4, nb)), fmt='s16')
          END DO
          CALL xml_addNewLine(xf)
          DEALLOCATE(tmp_dbuffer)
      CALL xml_endElement(u, 'pp_q')    
      !
      ! Write charge multipoles (only if PAW)
      IF ( upf%tpawp ) THEN
         CALL xml_addComment(u, 'augmentation charge multipoles ( only for PAW) ' //&
                                'multipole array dims = (nbeta,nbeta,2*lmax+1)')
         CALL xml_newElement(u, 'pp_multipoles')
              CALL xml_addAttribute(u, 'nbeta', upf%nbeta)
              CALL xml_addAttribute(u, 'lmax', upf%lmax)
              nb = upf%nbeta*upf%nbeta*(2*upf%lmax+1)
              ALLOCATE (tmp_dbuffer(nb))
              tmp_dbuffer = reshape(upf%paw%augmom,[nb]) 
              DO irow = 1, nb, 4
                 CALL xml_addNewLine(u)
                 CALL xml_addCharacters(u, tmp_dbuffer(irow: min(irow-1+4, nb) ), fmt='s16')
              END DO
              CALL xml_addNewLine(u)
              DEALLOCATE(tmp_dbuffer)
         CALL xml_endElement(u, 'pp_multipoles') 
      ENDIF
      !
      ! Write polinomial coefficients for Q_ij expansion at small radius
      !IF ( upf%nqf > 0) THEN
      !   CALL iotk_write_comment(u, ' polinomial expansion of Q_ij at small radius ')
      !   CALL iotk_write_dat(u, 'PP_QFCOEF',upf%qfcoef, attr=attr, columns=4)
      !   CALL iotk_write_dat(u, 'PP_RINNER',upf%rinner, attr=attr, columns=4)
      !ENDIF
      !
      ! Write augmentation charge Q_ij
      loop_on_nb: DO nb = 1,upf%nbeta
         ln = upf%lll(nb)
         loop_on_mb: DO mb = nb,upf%nbeta
            lm = upf%lll(mb)
            nmb = mb * (mb-1) /2 + nb
            IF( upf%q_with_l ) THEN
               loop_on_l: DO l = abs(ln-lm),ln+lm,2 ! only even terms
                  isnull = .FALSE. 
                  IF( upf%tpawp ) isnull = (abs(upf%paw%augmom(nb,mb,l)) < upf%qqq_eps)
                  IF( isnull) CYCLE loop_on_l
                  CALL xml_NewElement(u, 'pp_qijl')
                      CALL xml_addAttribute(u, 'first_index',  nb )
                      CALL xml_addAttribute(u, 'second_index', mb)
                      CALL xml_addAttribute(u, 'composite_index', nmb)
                      CALL xml_addAttribute(u, 'angular_momentum', l)
                      CALL xml_addAttribute(u, 'size', upf%mesh ) 
                      DO irow = 1, upf%mesh, 4 
                         CALL xml_addNewLine(u) 
                         CALL xml_addCharacters(u,upf%qfuncl(irow:min(irow-1+4,upf%mesh),nmb,l), fmt = 's16')
                      END DO
                      CALL xml_addNewLine(xf)
                  CALL xml_EndElement(u, 'pp_qijl')  
               ENDDO loop_on_l
            ELSE
               isnull = .FALSE. 
               IF  ( upf%tpawp ) isnull = ( abs(upf%qqq(nb,mb)) < upf%qqq_eps )
               IF (isnull) CYCLE loop_on_mb
               CALL xml_NewElement(u, 'pp_qij')   
                   CALL xml_addAttribute(u, 'size', upf%mesh)
                   CALL xml_addAttribute(u, 'first_index',  nb )
                   CALL xml_addAttribute(u, 'second_index', mb)
                   CALL xml_addAttribute(u, 'composite_index', nmb)
                   DO irow = 1, upf%mesh, 4
                     CALL xml_addNewLine(u)
                     CALL xml_addCharacters(u, upf%qfunc(irow:min(irow-1+4, upf%mesh), nmb), fmt = 's16')
                   END DO
                   CALL xml_addNewLine(xf)
               CALL xml_endElement(u, 'pp_qij') 
               !
            ENDIF
         ENDDO loop_on_mb
      ENDDO  loop_on_nb 
      !
      CALL xml_endElement(u, 'pp_augmentation') 
      !
   ENDIF augmentation
   !
   CALL xml_endElement(u, 'pp_nonlocal')
   !
   RETURN
 END SUBROUTINE write_nonlocal
 !
 SUBROUTINE write_pswfc(u, upf)
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)      :: u    ! i/o unit descriptor
   TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
   INTEGER :: ierr ! /= 0 if something went wrong

   !
   INTEGER :: nw
   !
   CALL xml_newElement(u, 'pp_pswfc')
   !
   DO nw = 1,upf%nwfc
      CALL xml_newElement(u, 'pp_chi') 
          CALL xml_addAttribute(u, 'size', upf%mesh)
          CALL xml_addAttribute(u, 'index',         nw )
          CALL xml_addAttribute(u, 'label',         upf%els(nw))
          CALL xml_addAttribute(u, 'l',             upf%lchi(nw))
          IF ( upf%has_so) THEN 
             CALL xml_addAttribute(u, 'nn', upf%nn(nw)) 
             CALL xml_addAttribute (u, 'jchi', upf%jchi(nw))
          END IF
          CALL xml_addAttribute(u, 'occupation',    upf%oc(nw))
          CALL xml_addAttribute(u, 'n',             upf%nchi(nw))
          CALL xml_addAttribute(u, 'pseudo_energy', upf%epseu(nw))
          CALL xml_addAttribute(u, 'cutoff_radius', upf%rcut_chi(nw))
          CALL xml_addAttribute(u, 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw))
          DO irow =1, upf%mesh, 4
             CALL xml_addNewLine(u)
             CALL xml_addCharacters(u, upf%chi(irow:min(irow-1+4,upf%mesh),nw) , fmt = 's16')
          END DO 
          CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_chi') 
   ENDDO
   !
   CALL xml_endElement(u, 'pp_pswfc')
   !
   RETURN
 END SUBROUTINE write_pswfc
 !
! SUBROUTINE write_spin_orb(u, upf)
!   IMPLICIT NONE
!   INTEGER,INTENT(IN)           :: u    ! i/o unit
!   TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
!   INTEGER :: ierr ! /= 0 if something went wrong
!
!   CHARACTER(len=iotk_attlenx) :: attr
!   !
!   INTEGER :: nw, nb
!   !
!   IF (.not. upf%has_so) RETURN
!   !
!   CALL iotk_write_begin(u, 'PP_SPIN_ORB')
!   !
!   DO nw = 1,upf%nwfc
!      CALL iotk_write_attr(attr, 'index', nw, first=.true.)
!      CALL iotk_write_attr(attr, 'els',   upf%els(nw))
!      CALL iotk_write_attr(attr, 'nn',    upf%nn(nw))
!      CALL iotk_write_attr(attr, 'lchi',  upf%lchi(nw))
!      CALL iotk_write_attr(attr, 'jchi',  upf%jchi(nw))
!      CALL iotk_write_attr(attr, 'oc',    upf%oc(nw))
!      CALL iotk_write_empty(u, 'PP_RELWFC'//iotk_index(nw),&
!           attr=attr)
!   ENDDO
!   !
!   DO nb = 1,upf%nbeta
!      CALL iotk_write_attr(attr, 'index', nb, first=.true.)
!      CALL iotk_write_attr(attr, 'lll',   upf%lll(nb))
!      CALL iotk_write_attr(attr, 'jjj',   upf%jjj(nb))
!      CALL iotk_write_empty(u, 'PP_RELBETA'//iotk_index(nb),&
!           attr=attr)
!   ENDDO
!   !
!   CALL iotk_write_end(u, 'PP_SPIN_ORB')
!   !
!   RETURN
! END SUBROUTINE write_spin_orb
 !

 SUBROUTINE write_full_wfc(u, upf)
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)      :: u    ! i/o unit descriptor
   TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
   INTEGER :: ierr ! /= 0 if something went wrong
   !
   INTEGER :: nb

   IF(.not. upf%has_wfc) RETURN
   CALL xml_NewElement(u, 'pp_full_wfc')
   CALL xml_addAttribute(u, 'number_of_wfc', upf%nbeta )
   ! All-electron wavefunctions corresponding to beta functions
   DO nb = 1,upf%nbeta
      CALL xml_NewElement(u, 'pp_aewfc')
          CALL xml_AddAttribute(u, 'index',      nb )
          CALL xml_AddAttribute(u, 'label',      upf%els_beta(nb))
          CALL xml_AddAttribute(u, 'l',          upf%lll(nb))
          DO irow = 1, upf%mesh, 4
             CALL xml_addNewLine(u)
             CALL xml_addCharacters(u, upf%aewfc(irow:min(irow-1+4, upf%mesh), nb), fmt = 's16')
          END DO
          CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_aewfc')
   ENDDO
   IF (upf%has_so.and.upf%tpawp) THEN
      DO nb = 1,upf%nbeta
         CALL xml_NewElement(u, 'pp_aewfc_rel')
             CALL xml_AddAttribute(u, 'size', upf%mesh) 
             CALL xml_addAttribute(u, 'index',      nb )
             CALL xml_addAttribute(u, 'label',      upf%els_beta(nb))
             CALL xml_addAttribute(u, 'l',          upf%lll(nb))
             CALL xml_addAttribute(u, 'j',          upf%jjj(nb))
             DO irow = 1, upf%mesh, 4
                CALL xml_addNewLine(u)
                CALL xml_addCharacters(u, upf%paw%aewfc_rel(irow:min(irow-1+4, upf%mesh), nb) , fmt = 's16')
             END DO
             CALL xml_addNewLine(xf)
         CALL xml_endElement(u, 'pp_aewfc_rel') 
      ENDDO
   ENDIF
   ! Pseudo wavefunctions 
   DO nb = 1,upf%nbeta
      CALL xml_newElement(u, 'pp_pswfc') 
          CALL xml_addAttribute(u, 'size', upf%mesh) 
          CALL xml_addAttribute(u, 'index',      nb )
          CALL xml_addAttribute(u, 'label',      upf%els_beta(nb))
          CALL xml_addAttribute(u, 'l',          upf%lll(nb))
          DO irow = 1, upf%mesh, 4
             CALL xml_addNewLine(u)
             CALL xml_AddCharacters(u, upf%pswfc(irow:min(irow-1+4, upf%mesh), nb) , fmt = 's16')
          END DO
          CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_pswfc')
   ENDDO
   ! Finalize
   CALL xml_endELement(u, 'pp_full_wfc')
 END SUBROUTINE write_full_wfc

 SUBROUTINE write_paw(u, upf)
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)      :: u    ! i/o unit descriptor
   TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
   INTEGER :: ierr ! /= 0 if something went wrong
   !
   INTEGER :: nb

   IF (.not. upf%tpawp ) RETURN
   CALL xml_newElement(u, 'pp_paw')
      CALL xml_addAttribute(u, 'paw_data_format', upf%paw_data_format )
      CALL xml_addAttribute(u, 'core_energy',     upf%paw%core_energy)
      ! Full occupation (not only > 0 ones)
      CALL xml_newElement(u, 'pp_occupations') 
         CALL xml_addAttribute(u, 'size', upf%nbeta)
         DO irow = 1, upf%nbeta, 4
            CALL xml_addNewLine(u)
            CALL xml_addCharacters(u, upf%paw%oc(irow:min(irow-1+4,upf%nbeta)), fmt = 's16')
         END DO 
         CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_occupations')  
      ! All-electron core charge
      CALL xml_newElement(u, 'pp_ae_nlcc')
          CALL xml_AddAttribute(u, 'size' , upf%mesh)
          DO irow = 1, upf%mesh, 4
             CALL xml_addNewLine(u)
             CALL xml_addCharacters(u, upf%paw%ae_rho_atc(irow:min(irow-1+4, upf%mesh)), fmt ='s16')
          END DO
          CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_ae_nlcc')
      ! All-electron local potential
      CALL xml_newElement(u, 'pp_ae_vloc')
         CALL xml_addAttribute(u, 'size', upf%mesh) 
         DO irow = 1, upf%mesh, 4
             CALL xml_addNewLine(u)
             CALL xml_addCharacters(u, upf%paw%ae_vloc(irow:min(irow-1+4, upf%mesh)), fmt ='s16')
         END DO
         CALL xml_addNewLine(xf)
      CALL xml_endElement(u, 'pp_ae_vloc') 
   !
   CALL xml_endElement(u, 'pp_paw')
   !
   RETURN
 END SUBROUTINE write_paw
 !
 SUBROUTINE write_gipaw(u, upf)
   IMPLICIT NONE
   TYPE(xmlf_t),INTENT(INOUT)      :: u    ! i/o unit descriptor
   TYPE(pseudo_upf),INTENT(IN)  :: upf  ! the pseudo data
   INTEGER :: ierr ! /= 0 if something went wrong
   !
   INTEGER :: nb
   IF (.not. upf%has_gipaw) RETURN

   CALL xml_newElement(u, 'pp_gipaw') 
   CALL xml_addAttribute(u, 'gipaw_data_format', upf%gipaw_data_format ) 
   DO nb = 1,upf%gipaw_ncore_orbitals
      CALL xml_NewElement(u, 'pp_gipaw_core_orbital') 
          CALL xml_addAttribute(u, 'size', upf%mesh)
          CALL xml_addAttribute(u, 'index', nb )
          CALL xml_addAttribute(u, 'label', upf%gipaw_core_orbital_el(nb))
          CALL xml_addAttribute(u, 'n',     upf%gipaw_core_orbital_n(nb))
          CALL xml_addAttribute(u, 'l',     upf%gipaw_core_orbital_l(nb))
          DO irow = 1, upf%mesh, 4
             CALL xml_addNewLine(u)
             CALL xml_addCharacters(u, upf%gipaw_core_orbital(irow:min(irow-1+4,upf%mesh), nb), fmt = 's16')
          END DO
          CALL xml_addNewLine(xf)
      CALL xml_EndElement(u, 'pp_gipaw_core_orbital') 
   ENDDO
   !
   ! Only write core orbitals in the PAW as GIPAW case
   IF (upf%paw_as_gipaw) THEN
      CALL xml_EndElement(u, 'pp_gipaw')
      RETURN
   ENDIF
   !
   ! Write valence all-electron and pseudo orbitals
   !
   DO nb = 1,upf%gipaw_wfs_nchannels
      CALL xml_NewElement(u, 'pp_gipaw_orbital')
          CALL xml_addAttribute(u, 'index', nb )
          CALL xml_addAttribute(u, 'label', upf%gipaw_wfs_el(nb))
          CALL xml_addAttribute(u, 'l',     upf%gipaw_wfs_ll(nb))
          CALL xml_addAttribute(u, 'cutoff_radius',           upf%gipaw_wfs_rcut(nb), fmt = 's16')
          CALL xml_addAttribute(u, 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb), fmt = 's16')
      !
          CALL xml_NewElement(u, 'pp_gipaw_wfs_ae') 
              CALL xml_addAttribute(u, 'size', upf%mesh)
              DO irow = 1, upf%mesh, 4 
                 CALL xml_addNewLine(u)
                 CALL xml_addCharacters(u, upf%gipaw_wfs_ae(irow:min(irow-1+4, upf%mesh), nb), fmt = 's16')
              END DO
              CALL xml_addNewLine(xf)
          CALL xml_endElement (u, 'pp_gipaw_wfs_ae')
          CALL xml_NewElement(u, 'pp_gipaw_wfs_ps') 
              CALL xml_addAttribute(u, 'size', upf%mesh)
              DO irow = 1, upf%mesh, 4 
                 CALL xml_addNewLine(u)
                 CALL xml_addCharacters(u, upf%gipaw_wfs_ps(irow:min(irow-1+4, upf%mesh), nb), fmt = 's16')
              END DO
              CALL xml_addNewLine(xf)
          CALL xml_endElement (u, 'pp_gipaw_wfs_ps')
      !
      CALL xml_endElement(u, 'pp_gipaw_orbital' )
   ENDDO
   !
   ! Write all-electron and pseudo local potentials
   CALL xml_NewElement(u, 'pp_gipaw_vlocal')
       CALL xml_NewElement(u, 'pp_gipaw_vlocal_ae')
            CALL xml_addAttribute(u, 'size', upf%mesh)
            DO irow = 1, upf%mesh, 4
               CALL xml_addNewLine(u)
               CALL xml_addCharacters(u, upf%gipaw_vlocal_ae(irow:min(irow-1+4, upf%mesh)), fmt = 's16')
            END DO
            CALL xml_addNewLine(xf)
       CALL xml_endElement(u, 'pp_gipaw_vlocal_ae')
       CALL xml_NewElement(u, 'pp_gipaw_vlocal_ps')
            CALL xml_addAttribute(u, 'size', upf%mesh)
            DO irow = 1, upf%mesh, 4
               CALL xml_addNewLine(u)
               CALL xml_addCharacters(u, upf%gipaw_vlocal_ps(irow:min(irow-1+4, upf%mesh)), fmt = 's16')
            END DO
            CALL xml_addNewLine(xf)
       CALL xml_endElement(u, 'pp_gipaw_vlocal_ps')
   CALL xml_endElement(u, 'pp_gipaw_vlocal')
   !
   CALL xml_endElement(u, 'pp_gipaw')

   RETURN
 END SUBROUTINE write_gipaw
 !
 ! Remove '<' and '>' from string, replacing them with '/', necessary
 ! or iotk will complain while read-skipping PP_INFO section.
END SUBROUTINE write_upf_schema

END MODULE write_upf_schema_module
