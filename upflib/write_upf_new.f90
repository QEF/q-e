!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
MODULE  write_upf_new
  !-----------------------------------------------------
  !! this module contains the simplified code for writing
  !! pseudopotential files in either UPF v.2 or xml
  ! 
  USE xmltools
  USE upf_kinds, ONLY: dp
  USE pseudo_types, ONLY: pseudo_upf, pseudo_config
  !
  LOGICAL :: v2
  !! true if UPF v.2 version, false if new UPF with xml schema
  INTEGER :: iun
  !! unit for writing data
  PRIVATE
  PUBLIC :: write_upf
  !
CONTAINS
  
  SUBROUTINE write_upf ( filename, upf, schema, conf, u_input)
    !
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)   :: filename
    !! name of the output file
    TYPE(pseudo_upf),INTENT(IN)   :: upf
    !! pseudo_upf structure containing all the pseudo data
    CHARACTER(LEN=*),INTENT(IN), OPTIONAL :: schema
    !!  optional, character flag which selects what schema will be used on writing 
    TYPE(pseudo_config), OPTIONAL, INTENT(IN)   :: conf 
    !!  optional, pseudo_conf data structure containing the atomic configuration used
    !!  to generate the pseudo
    INTEGER,OPTIONAL                        :: u_input
    !!  optional: unit of stdin for  the generation program, used to write the 
    !!  generation input in the upf file
    !
    CHARACTER(LEN=5) :: schema_='qe_pp'
    CHARACTER(LEN=*), PARAMETER   :: QE_PP_URI = &
         "http://www.quantum-espresso.org/ns/qes/qe_pp-1.0", &
         XSI = "http://www.w3.org/2001/XMLSchema-instance", &
         XSD_VERSION = "QE_PP-1.0"
    !
    IF ( PRESENT(schema) ) schema_ = schema
    SELECT CASE (TRIM(schema_))
    CASE ('qe_pp', 'QE_PP')
       v2 = .false.
    CASE ('V2', 'v2' ,'upf', 'UPF') 
       v2 = .true.
    END SELECT
    !
    iun = xml_open_file ( filename )
    IF ( iun == -1 ) CALL upf_error('write_upf', 'cannot open file',1)
    !
    ! The starting line, header and info sections differ a lot
    ! between UPF v.2 and UPF with schema so we call different routines
    !
    IF ( v2 ) THEN
       !
       CALL add_attr ('version', upf%nv )
       CALL xmlw_opentag ( 'UPF')
       !
       ! pp_info
       !
       CALL write_pp_info_v2 ( upf, conf, u_input )
       !
       ! pp_header
       !
       CALL write_pp_header_v2 ( upf )
       !
    ELSE
       !
       call add_attr( 'version','1.0')
       call add_attr( 'encoding','UTF-8')
       CALL xmlw_writetag ( 'xml', '?' )
       call add_attr( 'xsi:schemalocation', QE_PP_URI//' '//QE_PP_URI//'.xsd')
       call add_attr( 'xmlns:xsi',XSI)
       call add_attr( 'xmlns:qe_pp', QE_PP_URI)
       CALL xmlw_opentag ( 'qe_pp:pseudo' )
       CALL xmlw_writetag ( 'xsd_version', XSD_VERSION )
       !
       ! pp_info
       !
       CALL write_pp_info_schema ( upf, conf, u_input )
       !
       ! pp_header
       !
       CALL write_pp_header_schema ( upf )
       !
    END IF
    !
    ! From here on the format of v2 and schema do not differ much:
    ! the most frequent difference is capitalization of tags
    ! (see function capitalize_if_v2)
    !
    CALL write_pp_mesh ( upf )
    !
    IF( upf%nlcc ) THEN
       CALL add_attr( 'size', upf%mesh )
       CALL xmlw_writetag(capitalize_if_v2('pp_nlcc'), upf%rho_atc(1:upf%mesh))
    END IF
    IF( .NOT. upf%tcoulombp ) THEN
       CALL add_attr( 'size', upf%mesh )
       CALL xmlw_writetag( capitalize_if_v2('pp_local'), upf%vloc(1:upf%mesh))
    END IF
    !
    CALL write_pp_semilocal ( upf )
    !
    CALL write_pp_nonlocal ( upf )
    !
    CALL write_pp_pswfc ( upf )
    !
    CALL write_pp_full_wfc ( upf )
    !
    CALL add_attr( 'size', upf%mesh )
    CALL xmlw_writetag( capitalize_if_v2('pp_rhoatom'), upf%rho_at(1:upf%mesh))
    !
    CALL write_pp_metagga ( upf )
    !
    CALL write_pp_spinorb ( upf )
    !
    CALL write_pp_paw ( upf )
    !
    CALL write_pp_gipaw ( upf )
    !
    ! close initial tag, qe_pp:pseudo or UPF
    !
    CALL xmlw_closetag ( )
    !
    CALL xml_closefile ( )
    !
    RETURN
    !
  END SUBROUTINE write_upf
  !
  FUNCTION capitalize_if_v2 ( strin ) RESULT ( strout )
    !
    ! returns a capitalized string for UPF v.2, the same string otherwise
    ! (UPF v.2 uses capitalized tags, UPF with schema use lowercase)
    !
    USE upf_utils, ONLY: capital
    IMPLICIT NONE
    CHARACTER(LEN=*) :: strin
    !
    INTEGER :: n
    CHARACTER(LEN=:), ALLOCATABLE :: strout
    !
    IF ( v2 ) THEN
       strout = ''
       DO n = 1,LEN_TRIM(strin)
          strout = strout // capital(strin(n:n))
       END DO
    ELSE
       strout = TRIM(strin)
    END IF
    !
  END FUNCTION capitalize_if_v2
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_info_schema ( upf, conf, u_input )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    ! optional: configuration used to generate the pseudopotential
    TYPE(pseudo_config), OPTIONAL, INTENT(IN) :: conf
    ! optional: unit pointing to input file containing generation data
    INTEGER, OPTIONAL, INTENT(IN):: u_input
    !
#include "qe_version.h"
    INTEGER :: nw, nb
    !
    CALL xmlw_opentag ( 'pp_info' )
    CALL xmlw_writetag ( 'generated', xml_protect(upf%generated) )
    CALL add_attr( 'NAME', 'QE Atomic Code' )
    CALL add_attr( 'VERSION', version_number )
    CALL xmlw_writetag ( 'creator', xml_protect(upf%author) )
    CALL add_attr( 'DATE', upf%date )
    CALL xmlw_writetag ( 'created', '' )
    !
    IF ( PRESENT(u_input) ) CALL copy_input_data ( u_input )
    !
    CALL xmlw_writetag ( 'type', upf%typ )
    IF (TRIM(upf%rel)=='full' .OR. TRIM(upf%rel)=='scalar' ) THEN
       CALL xmlw_writetag ( 'relativistic_effects', upf%rel )
    ELSE
       CALL xmlw_writetag ( 'relativistic_effects', 'none' )
    ENDIF
    CALL xmlw_writetag ( 'element', upf%psd )
    CALL xmlw_writetag ( 'functional', upf%dft )
    CALL add_attr( 'ecutwfc', upf%ecutwfc )
    IF (upf%tpawp .OR. upf%tvanp ) THEN
       CALL add_attr( 'ecutrho', upf%ecutrho )
       CALL xmlw_writetag ( 'suggested_basis', '' )
    ELSE
       CALL xmlw_writetag ( 'suggested_basis', '' )
    END IF
    DO nw =1, upf%nwfc
       IF( upf%oc(nw) >= 0.0_dp) THEN 
          CALL add_attr( 'nl', upf%els(nw) )
          CALL add_attr( 'pn', upf%nchi(nw) )
          CALL add_attr( 'l', upf%lchi(nw) )
          CALL xmlw_opentag ( "valence_orbital" )
          CALL xmlw_writetag ( "occupation", upf%oc(nw) )
          CALL xmlw_writetag ( "Rcut", upf%rcut_chi(nw) )
          IF (upf%rcutus_chi(nw) > 0.0_dp) &
               CALL xmlw_writetag (  "RcutUS", upf%rcutus_chi(nw) )
          CALL xmlw_writetag ( "Epseu", upf%epseu(nw) )
          CALL xmlw_closetag ( )
       END IF
    END DO
    IF( present(conf) ) THEN
       CALL xmlw_opentag ( "generation_configuration" )
       DO nb = 1,conf%nwfs
          WRITE(iun, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
               conf%els(nb), conf%nns(nb), &
               conf%lls(nb), conf%ocs(nb), conf%rcut(nb), &
               conf%rcutus(nb), conf%enls(nb)
       ENDDO
       WRITE(iun,'(4x,2a)') 'Pseudization used: ',TRIM(conf%pseud)
       CALL xmlw_closetag ( )
    ENDIF
    IF( TRIM(upf%comment) /= ' ') &
       WRITE(iun,'("<!--",a,"-->")') TRIM(upf%comment)
    CALL xmlw_closetag ( ) ! end pp_info
    !
  END SUBROUTINE write_pp_info_schema
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_header_schema ( upf )
    !--------------------------------------------------------
    !
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    CALL xmlw_opentag ( 'pp_header')
    CALL xmlw_writetag( 'element', upf%psd )
    CALL xmlw_writetag( 'z_valence', upf%zp )
    CALL xmlw_writetag( 'type', upf%typ )
    CALL xmlw_writetag( 'functional', upf%dft )
    CALL xmlw_writetag( 'relativistic', upf%rel )
    CALL xmlw_writetag( 'is_ultrasoft', upf%tvanp )
    CALL xmlw_writetag( 'is_paw', upf%tpawp )
    CALL xmlw_writetag( 'is_coulomb', upf%tcoulombp )
    CALL xmlw_writetag( 'has_so', upf%has_so )
    CALL xmlw_writetag( 'has_wfc', upf%has_wfc )
    CALL xmlw_writetag( 'has_gipaw', upf%has_gipaw )
    CALL xmlw_writetag( 'paw_as_gipaw', upf%paw_as_gipaw)
    CALL xmlw_writetag( 'core_correction', upf%nlcc)
    CALL xmlw_writetag( 'with_metagga_info', upf%with_metagga_info)
    CALL xmlw_writetag( 'total_psenergy', upf%etotps )
    CALL xmlw_writetag( 'wfc_cutoff', upf%ecutwfc )
    CALL xmlw_writetag( 'rho_cutoff', upf%ecutrho )
    CALL xmlw_writetag( 'l_max', upf%lmax )
    CALL xmlw_writetag( 'l_max_rho', upf%lmax_rho )
    CALL xmlw_writetag( 'l_local', upf%lloc )
    CALL xmlw_writetag( 'mesh_size', upf%mesh )
    CALL xmlw_writetag( 'number_of_wfc', upf%nwfc )
    CALL xmlw_writetag( 'number_of_proj', upf%nbeta )
    CALL xmlw_closetag( ) ! end pp_header
    !
  END SUBROUTINE write_pp_header_schema
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_info_v2 ( upf, conf, u_input )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    ! optional: configuration used to generate the pseudopotential
    TYPE(pseudo_config), OPTIONAL, INTENT(IN) :: conf
    ! optional: unit pointing to input file containing generation data
    INTEGER, OPTIONAL, INTENT(IN):: u_input
    !
    INTEGER :: nw, nb
    CALL xmlw_opentag ( 'PP_INFO' )
    WRITE(iun,'(4x,a)') TRIM(upf%generated)
    WRITE(iun,'(4x,a)') 'Author: '//TRIM(upf%author)
    WRITE(iun,'(4x,a)') 'Generation date: '//TRIM(upf%date) 
    WRITE(iun,'(4x,a)') 'Pseudopotential type: '//TRIM(upf%typ) 
    WRITE(iun,'(4x,a)') 'Element: '//TRIM(upf%psd) 
    WRITE(iun,'(4x,a)') 'Functional: '//TRIM(upf%dft) 
    !
    ! Cutoff Information
    WRITE(iun,'(4x,a,f5.0,a)') &
         'Suggested minimum cutoff for wavefunctions:',upf%ecutwfc,' Ry'
    WRITE(iun, '(4x,a,f5.0,a)') &
         'Suggested minimum cutoff for charge density:',upf%ecutrho,' Ry'
    ! Write relativistic information
    IF (TRIM(upf%rel)=='full') THEN
       WRITE(iun, '(4x,a)') &
            "The Pseudo was generated with a Fully-Relativistic Calculation"
    ELSE IF (TRIM(upf%rel)=='scalar') THEN
       WRITE(iun, '(4x,a)') &
            "The Pseudo was generated with a Scalar-Relativistic Calculation"
    ELSE
       WRITE(iun, '(4x,a)') &
            "The Pseudo was generated with a Non-Relativistic Calculation"
    ENDIF
    !
    ! Write local potential information
    IF (upf%lloc >= 0 ) THEN
       WRITE(iun, '(4x,a,i3,f9.4)') &
            "L component and cutoff radius for Local Potential:", upf%lloc, upf%rcloc
    ELSE IF (upf%lloc == -1 ) THEN
       WRITE(iun, '(4x,a,f9.4)') &
            "Local Potential by smoothing AE potential with Bessel fncs, cutoff radius:", upf%rcloc
    ELSE IF (upf%lloc == -2 ) THEN
       WRITE(iun, '(4x,a,f9.4)') &
            "Local Potential according to Troullier-Martins recipe, cutoff radius:", upf%rcloc
    ELSE
       WRITE(iun, '(4x,a,i3,f9.4)') &
            "Local Potential: unknown format, L component and cutoff radius:",upf%lloc, upf%rcloc
    ENDIF
    !
    IF (upf%has_so) WRITE(iun, '(4x,a,i3,f9.4)') &
            "Pseudopotential contains additional information for spin-orbit calculations."
    IF (upf%has_gipaw) WRITE(iun, '(4x,a,i3,f9.4)') &
         "Pseudopotential contains additional information for GIPAW reconstruction."
    !
    ! Write valence orbitals information
    WRITE(iun, '(4x,a)') 'Valence configuration: '
    WRITE(iun, '(4x,a2,2a3,a6,2a11,1a13)') &
         "nl"," pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    DO nb = 1, upf%nwfc
       IF(upf%oc(nb) >= 0._dp) THEN
          WRITE(iun, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
               upf%els(nb), upf%nchi(nb), &
               upf%lchi(nb), upf%oc(nb), upf%rcut_chi(nb), &
               upf%rcutus_chi(nb), upf%epseu(nb)
       ENDIF
    END DO
    IF( present(conf) ) THEN
       WRITE(iun, '(4x,a)') 'Generation configuration:'
       DO nb = 1,conf%nwfs
          WRITE(iun, '(4x,a2,2i3,f6.2,2f11.3,1f13.6)') &
               conf%els(nb), conf%nns(nb), &
               conf%lls(nb), conf%ocs(nb), conf%rcut(nb), &
               conf%rcutus(nb), conf%enls(nb)
       ENDDO
       WRITE(iun,'(4x,2a)') 'Pseudization used: ',TRIM(conf%pseud)
    ELSE
       WRITE(iun, '(4x,a)') 'Generation configuration: not available.'
    ENDIF

    IF(TRIM(upf%comment) /= ' ') WRITE(iun, '(4x,"Comment:",2x,a)') &
         xml_protect(TRIM(upf%comment))
    !
    IF ( PRESENT(u_input) ) CALL copy_input_data ( u_input )
    !
    ! end PP_INFO
    !
    CALL xmlw_closetag ( )
    WRITE(iun, '("    <!-- END OF HUMAN READABLE SECTION -->")')
    !
  END SUBROUTINE write_pp_info_v2
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_header_v2 ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    call add_attr("generated", xml_protect(upf%generated) )
    call add_attr("author", xml_protect(upf%author) )
    call add_attr("date", upf%date )
    call add_attr("comment", xml_protect(upf%comment) )
    call add_attr("element", upf%psd )
    call add_attr("pseudo_type", upf%typ )
    call add_attr("relativistic", upf%rel )
    call add_attr("is_ultrasoft", upf%tvanp )
    call add_attr("is_paw", upf%tpawp )
    call add_attr("is_coulomb", upf%tcoulombp )
    call add_attr("has_so", upf%has_so )
    call add_attr("has_wfc", upf%has_wfc )
    call add_attr("has_gipaw", upf%has_gipaw )
    call add_attr("paw_as_gipaw", upf%paw_as_gipaw )
    call add_attr("core_correction", upf%nlcc )
    call add_attr("with_metagga_info", upf%with_metagga_info)
    call add_attr("functional", upf%dft )
    call add_attr("z_valence", upf%zp )
    call add_attr("total_psenergy", upf%etotps )
    call add_attr("wfc_cutoff", upf%ecutwfc )
    call add_attr("rho_cutoff", upf%ecutrho )
    call add_attr("l_max", upf%lmax )
    call add_attr("l_max_rho", upf%lmax_rho )
    call add_attr("l_local", upf%lloc )
    call add_attr("mesh_size", upf%mesh )
    call add_attr("number_of_wfc", upf%nwfc )
    call add_attr("number_of_proj", upf%nbeta )
    !
    CALL xmlw_writetag ( "PP_HEADER", '')
    !
  END SUBROUTINE write_pp_header_v2
  !
  !--------------------------------------------------------
  SUBROUTINE copy_input_data ( u_input )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: u_input
    CHARACTER(len=256) :: line
    LOGICAL :: opnd
    !
    ! copy content of input file used in pseudopotential generation
    !
    INQUIRE (unit=u_input, opened=opnd)
    IF (opnd) THEN
       IF ( v2 ) THEN
          CALL xmlw_opentag ( 'PP_INPUTFILE' )
       ELSE
          CALL add_attr ('program', 'ld1.x' )
          CALL xmlw_opentag ( 'input' ) 
       END IF
       REWIND (unit=u_input)
       read_write_loop: DO
          READ (u_input, '(A)',end=20,err=25) line
          WRITE (iun, '(A)') xml_protect(line)
          CYCLE read_write_loop
25        WRITE(*,'(5X,"write_upf::copy_input_data warning: problem writing input data")')
20        EXIT read_write_loop
       END DO read_write_loop
       CALL xmlw_closetag ( ) ! 'input'
    ELSE
       WRITE(*,'(5X,"write_upf::copy_input_data warning: input file not open")')
    END IF
    !
  END SUBROUTINE copy_input_data
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_mesh ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    IF ( upf%dx > 0.d0) THEN
       CALL add_attr( 'mesh', upf%mesh )
       CALL add_attr( 'dx', upf%dx )
       CALL add_attr( 'xmin', upf%xmin )
       CALL add_attr( 'rmax', upf%rmax )
       CALL add_attr( 'zmesh', upf%zmesh )
       CALL xmlw_opentag( capitalize_if_v2('pp_mesh') )
    ELSE
       CALL xmlw_opentag( capitalize_if_v2('pp_mesh') )
    END IF
    CALL xmlw_writetag( capitalize_if_v2('pp_r'), upf%r(1:upf%mesh) )
    CALL xmlw_writetag( capitalize_if_v2('pp_rab'), upf%rab(1:upf%mesh) )
    !
    CALL xmlw_closetag( ) ! end pp_mesh
    !
  END SUBROUTINE write_pp_mesh
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_semilocal ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    INTEGER :: nb, ind, l
    CHARACTER(LEN=8) :: tag
    !
    IF ( upf%typ == "SL" ) THEN
       CALL xmlw_opentag( capitalize_if_v2('pp_semilocal') )
       !
       DO nb = 1,upf%nbeta
          l = upf%lll(nb)
          ind = 1
          IF ( upf%has_so ) THEN
             IF ( l > 0 .AND. ABS(upf%jjj(nb)-l-0.5_dp) < 0.001_dp ) ind = 2
          END IF
          IF ( v2 ) THEN
             tag = 'PP_VNL.'//i2c(l)
          ELSE
             tag = 'vnl'
          END IF
          CALL add_attr( 'l', l )
          IF ( upf%has_so ) THEN
             CALL add_attr( 'j', upf%jjj(nb) )
             CALL xmlw_writetag( tag, upf%vnl(1:upf%mesh,l,ind) )
          ELSE
             CALL xmlw_writetag( tag, upf%vnl(1:upf%mesh,l,ind) )
          END IF
       END DO
       CALL xmlw_closetag( ) ! end pp_semilocal
    END IF
    !
  END SUBROUTINE write_pp_semilocal
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_nonlocal ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    LOGICAL :: isnull
    INTEGER :: nb, ind, l, ln, lm, mb, nmb
    CHARACTER(LEN=15) :: tag
    !
    CALL xmlw_opentag( capitalize_if_v2('pp_nonlocal') )
    !
    DO nb = 1,upf%nbeta
       !
       IF ( v2 ) THEN
          tag = 'PP_BETA.'//i2c(nb)
       ELSE
          tag = 'pp_beta'
       END IF
       call add_attr( 'index', nb )
       call add_attr( 'label', upf%els_beta(nb) )
       call add_attr( 'angular_momentum', upf%lll(nb) )
       call add_attr( 'cutoff_radius_index', upf%kbeta(nb) )
       call add_attr( 'cutoff_radius', upf%rcut(nb) )
       call add_attr( 'ultrasoft_cutoff_radius', upf%rcutus(nb) )
       IF ( .NOT. v2 .AND. upf%has_so ) THEN
          call add_attr( 'tot_ang_mom', upf%jjj(nb) )
          CALL xmlw_writetag( tag, upf%beta(1:upf%mesh,nb) )
       ELSE
          CALL xmlw_writetag( tag, upf%beta(1:upf%mesh,nb) )
       END IF
    END DO
    !
    ! pp_dij (D_lm matrix)
    !
    call add_attr( 'columns',  upf%nbeta )
    call add_attr( 'rows', upf%nbeta )
    CALL xmlw_opentag( capitalize_if_v2 ('pp_dij') )
    DO nb = 1,upf%nbeta
       WRITE(iun,*) upf%dion(1:upf%nbeta,nb)
    END DO
    CALL xmlw_closetag( ) 
    !
    ! pp_augmentation
    !
    IF (upf%tvanp .or. upf%tpawp) THEN
       IF ( v2 ) THEN
          call add_attr('q_with_l', upf%q_with_l )
          call add_attr('nqf', upf%nqf )
          call add_attr('nqlc', upf%nqlc )
          IF (upf%tpawp) THEN
             CALL add_attr( 'shape', upf%paw%augshape )
             CALL add_attr( 'cutoff_r', upf%paw%raug )
             CALL add_attr( 'cutoff_r_index', upf%paw%iraug )
             CALL add_attr( 'augmentation_epsilon', upf%qqq_eps )
             CALL add_attr( 'l_max_aug', upf%paw%lmax_aug )
          ENDIF
       END IF
       CALL xmlw_opentag( capitalize_if_v2('pp_augmentation') )
       !
       IF ( .NOT. v2 ) THEN
          CALL xmlw_writetag( 'q_with_l', upf%q_with_l )
          CALL xmlw_writetag( 'nqf', upf%nqf )
          CALL xmlw_writetag( 'nqlc', upf%nqlc )
          IF (upf%tpawp) THEN
             CALL xmlw_writetag( 'shape', upf%paw%augshape )
             CALL xmlw_writetag( 'cutoff_r', upf%paw%raug )
             CALL xmlw_writetag( 'cutoff_r_index', upf%paw%iraug )
             CALL xmlw_writetag( 'augmentation_epsilon', upf%qqq_eps )
             CALL xmlw_writetag( 'l_max_aug', upf%paw%lmax_aug )
          ENDIF
       END IF
       !
       nb = upf%nbeta*upf%nbeta
       call add_attr( 'size', nb )
       CALL xmlw_opentag( capitalize_if_v2('pp_q') )
       DO nb = 1,upf%nbeta
          WRITE(iun,*) upf%qqq(1:upf%nbeta,nb)
       END DO
       CALL xmlw_closetag( )
       !
       IF ( upf%tpawp ) THEN
          WRITE(iun,"('<!--augmentation charge multipoles ( only for PAW) ',/,&
               &      'multipole array dims = (nbeta,nbeta,2*lmax+1)-->')")
          call add_attr( 'nbeta', upf%nbeta )
          call add_attr( 'lmax', upf%lmax )
          CALL xmlw_opentag( capitalize_if_v2('pp_multipoles') )
          DO l = 0,2*upf%lmax
             DO nb = 1,upf%nbeta
                WRITE(iun,*) upf%paw%augmom(1:upf%nbeta,nb,l)
             END DO
          END DO
          CALL xmlw_closetag ()
       ENDIF
       !
    END IF
    !
    ! Write polinomial coefficients for Q_ij expansion at small radius
    IF ( v2 .AND. upf%nqf > 0) THEN
       WRITE(iun,"('<!--polinomial expansion of Q_ij at small radius-->')")
       CALL xmlw_opentag('PP_QFCOEF')
       WRITE(iun,*) upf%qfcoef
       CALL xmlw_closetag ()
       CALL xmlw_opentag('PP_RINNER')
       WRITE(iun,*) upf%rinner
       CALL xmlw_closetag ()
    ENDIF
    !
    IF ( upf%tpawp .or. upf%tvanp ) THEN
       !
       ! Write augmentation charge Q_ij
       !
       loop_on_nb: DO nb = 1,upf%nbeta
          ln = upf%lll(nb)
          loop_on_mb: DO mb = nb,upf%nbeta
             lm = upf%lll(mb)
             nmb = mb * (mb-1) /2 + nb
             IF( upf%q_with_l ) THEN
                loop_on_l: DO l = abs(ln-lm),ln+lm,2 ! only even terms
                   isnull = .FALSE. 
                   IF( upf%tpawp ) isnull = (abs(upf%paw%augmom(nb,mb,l)) < upf%qqq_eps)
                   IF(isnull) CYCLE loop_on_l
                   IF ( v2 ) THEN
                      tag = 'PP_QIJL.'//i2c(nb)//'.'//i2c(mb)//'.'//i2c(l)
                   ELSE
                      tag = 'pp_qijl'
                   END IF
                   call add_attr( 'first_index', nb )
                   call add_attr( 'second_index', mb )
                   call add_attr( 'composite_index', nmb )
                   call add_attr( 'angular_momentum', l )
                   call add_attr( 'size', upf%mesh )
                   CALL xmlw_writetag( tag, upf%qfuncl(1:upf%mesh,nmb,l) )
                ENDDO loop_on_l
             ELSE
                isnull = .FALSE. 
                IF  ( upf%tpawp ) isnull = ( abs(upf%qqq(nb,mb)) < upf%qqq_eps )
                IF (isnull) CYCLE loop_on_mb
                IF ( v2 ) THEN
                   tag = 'PP_QIJ.'//i2c(nb)//'.'//i2c(mb)
                ELSE
                   tag = 'pp_qij'
                END IF
                call add_attr( 'size', upf%mesh )
                call add_attr( 'first_index', nb )
                call add_attr( 'second_index', mb )
                call add_attr( 'composite_index', nmb )
                CALL xmlw_writetag( tag, upf%qfunc(1:upf%mesh,nmb) )
                !
             ENDIF
          ENDDO loop_on_mb
       ENDDO  loop_on_nb
       !
       CALL xmlw_closetag( ) ! end pp_augmentation
       !
    END IF
    CALL xmlw_closetag( ) ! end pp_nonlocal
    !
  END SUBROUTINE write_pp_nonlocal
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_pswfc ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    INTEGER :: nw, ind, l
    CHARACTER(LEN=8) :: tag
    !
    CALL xmlw_opentag( capitalize_if_v2('pp_pswfc') )
    DO nw =1, upf%nwfc
       call add_attr( 'size', upf%mesh )
       call add_attr( 'index', nw )
       call add_attr( 'label', upf%els(nw) )
       call add_attr( 'l', upf%lchi(nw) )
       IF ( upf%has_so) THEN
          call add_attr( 'nn', upf%nchi(nw) )
          call add_attr( 'jchi', upf%jchi(nw) )
       END IF
       call add_attr( 'occupation', upf%oc(nw) )
       IF ( upf%nchi(nw) > upf%lchi(nw) ) call add_attr( 'n', upf%nchi(nw) )
       IF ( upf%epseu(nw) > 0.0_dp ) &
            call add_attr( 'pseudo_energy',upf%epseu(nw) )
       IF ( upf%rcut_chi(nw) > 0.0_dp ) &
            call add_attr( 'cutoff_radius',upf%rcut_chi(nw) )
       IF ( upf%rcutus_chi(nw) > 0.0_dp ) &
            call add_attr( 'ultrasoft_cutoff_radius',upf%rcutus_chi(nw) )
       !
       IF ( v2 ) THEN
          tag = 'PP_CHI.'//i2c(nw)
       ELSE
          tag = 'pp_chi'
       END IF
       CALL xmlw_writetag( tag, upf%chi(1:upf%mesh,nw) )
    END DO
    CALL xmlw_closetag( ) ! end pp_pswfc
    !
  END SUBROUTINE write_pp_pswfc
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_full_wfc ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    INTEGER :: nb
    CHARACTER(LEN=15) :: tag
    !
    IF ( upf%has_wfc ) THEN
       !
       call add_attr( 'number_of_wfc', upf%nbeta )
       CALL xmlw_opentag( capitalize_if_v2('pp_full_wfc') )
       !
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             tag = 'PP_AEWFC.'//i2c(nb)
          ELSE
             tag = 'pp_aewfc'
          END IF
          call add_attr( 'index', nb )
          call add_attr( 'label', upf%els_beta(nb) )
          call add_attr( 'l', upf%lll(nb) )
          CALL xmlw_writetag( tag, upf%aewfc(1:upf%mesh,nb) )
       END DO
       !
       IF ( upf%has_so .AND. upf%tpawp ) THEN
          DO nb = 1, upf%nbeta
             IF ( v2 ) THEN
                tag = 'PP_AEWFC_REL.'//i2c(nb)
             ELSE
                tag = 'pp_aewfc_rel'
             END IF
             call add_attr( 'index', nb )
             call add_attr( 'label', upf%els_beta(nb) )
             call add_attr( 'l', upf%lll(nb) )
             CALL xmlw_writetag(tag, upf%paw%aewfc_rel(1:upf%mesh,nb) )
          END DO
       END IF
       !
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             tag = 'PP_PSWFC.'//i2c(nb)
          ELSE
             tag = 'pp_pswfc'
          END IF
          call add_attr( 'size', upf%mesh )
          call add_attr( 'index', nb )
          call add_attr( 'label', upf%els_beta(nb) )
          call add_attr( 'l', upf%lll(nb) )
          CALL xmlw_writetag( tag, upf%pswfc(1:upf%mesh,nb) )
       END DO
       !
       CALL xmlw_closetag( )
       !
    END IF
    !
  END SUBROUTINE write_pp_full_wfc
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_metagga ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    if ( .NOT. upf%with_metagga_info ) RETURN
    !
    CALL xmlw_writetag( capitalize_if_v2('pp_taumod'), upf%tau_core(:) )
    CALL xmlw_writetag( capitalize_if_v2('pp_tauatom'), upf%tau_atom(:) )
    !
  END SUBROUTINE write_pp_metagga
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_spinorb ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    INTEGER :: nw, nb
    !
    IF ( .NOT. v2 .OR. .NOT. upf%has_so ) RETURN
    !
    CALL xmlw_opentag( 'PP_SPIN_ORB' )
    DO nw = 1,upf%nwfc
       CALL add_attr( 'index' , nw )
       CALL add_attr( 'els',   upf%els(nw) )
       CALL add_attr( 'nn',    upf%nchi(nw) )
       CALL add_attr( 'lchi',  upf%lchi(nw) )
       CALL add_attr( 'jchi',  upf%jchi(nw) )
       CALL add_attr( 'oc',    upf%oc(nw) )
       CALL xmlw_writetag( 'PP_RELWFC.'//i2c(nw), '' )
    ENDDO
    !
    DO nb = 1,upf%nbeta
       CALL add_attr( 'index' , nb )
       CALL add_attr( 'lll',  upf%lll(nb) )
       CALL add_attr( 'jjj',  upf%jjj(nb) )
       CALL xmlw_writetag( 'PP_RELBETA.'//i2c(nb), '' )
    ENDDO
    CALL xmlw_closetag () ! end pp_spin_orb
    !
  END SUBROUTINE write_pp_spinorb
  !
  !--------------------------------------------------------
  SUBROUTINE write_pp_paw ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    IF ( upf%tpawp ) THEN
       call add_attr( 'paw_data_format', upf%paw_data_format )
       call add_attr( 'core_energy', upf%paw%core_energy )
       CALL xmlw_opentag( capitalize_if_v2('pp_paw') )
       ! Full occupation (not only > 0 ones)
       call add_attr( 'size', upf%nbeta )
       CALL xmlw_writetag( capitalize_if_v2('pp_occupations'), &
            upf%paw%oc(1:upf%nbeta) )
       ! All-electron core charge
       call add_attr( 'size', upf%mesh )
       CALL xmlw_writetag( capitalize_if_v2('pp_ae_nlcc'), &
            upf%paw%ae_rho_atc(1:upf%mesh) )
       ! All-electron local potential
       call add_attr( 'size', upf%mesh )
       CALL xmlw_writetag( capitalize_if_v2('pp_ae_vloc'), &
            upf%paw%ae_vloc(1:upf%mesh) )
       CALL xmlw_closetag () ! end pp_paw
    END IF
    !
  END SUBROUTINE write_pp_paw
  !--------------------------------------------------------
  SUBROUTINE write_pp_gipaw ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    INTEGER :: nb
    CHARACTER(LEN=24) :: tag
    !
    IF (upf%has_gipaw) THEN
       call add_attr( 'gipaw_data_format', upf%gipaw_data_format )
       CALL xmlw_opentag( capitalize_if_v2('pp_gipaw') )
       IF ( v2 ) THEN
          call add_attr( 'number_of_core_orbitals', upf%gipaw_ncore_orbitals )
          CALL xmlw_opentag( 'PP_GIPAW_CORE_ORBITALS' )
       ELSE 
          CALL xmlw_writetag('number_of_core_orbitals', upf%gipaw_ncore_orbitals )
          IF ( .NOT. upf%paw_as_gipaw ) &
             CALL xmlw_writetag('number_of_valence_orbitals', upf%gipaw_wfs_nchannels) 
       END IF
       DO nb = 1,upf%gipaw_ncore_orbitals
          IF ( v2 ) THEN
             tag = "PP_GIPAW_CORE_ORBITAL."//i2c(nb)
          ELSE
             tag = 'pp_gipaw_core_orbital'
          END IF
          call add_attr( 'size', upf%mesh )
          call add_attr( 'index', nb )
          call add_attr( 'label', upf%gipaw_core_orbital_el(nb) )
          call add_attr( 'n', upf%gipaw_core_orbital_n(nb) )
          call add_attr( 'l', upf%gipaw_core_orbital_l(nb) )
          CALL xmlw_writetag( tag, upf%gipaw_core_orbital(1:upf%mesh,nb) )
       END DO
       IF ( v2 ) CALL xmlw_closetag ( )
       !
       ! Only core orbitals are written in the PAW as GIPAW case
       !
       IF ( .NOT. upf%paw_as_gipaw) THEN
          !
          ! Write valence all-electron and pseudo orbitals
          !
          IF ( v2 ) THEN
             call add_attr( 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels )
             CALL xmlw_opentag( 'PP_GIPAW_ORBITALS' )
          END IF
          DO nb = 1,upf%gipaw_wfs_nchannels
             IF ( v2 ) THEN
                tag = "PP_GIPAW_ORBITAL."//i2c(nb)
             ELSE
                tag = 'pp_gipaw_orbital'
             END IF
             call add_attr( 'index', nb )
             call add_attr( 'label', upf%gipaw_wfs_el(nb) )
             call add_attr( 'l', upf%gipaw_wfs_ll(nb) )
             call add_attr( 'cutoff_radius', upf%gipaw_wfs_rcut(nb) )
             call add_attr( 'ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb) )
             CALL xmlw_opentag( tag)
             !
             call add_attr( 'size', upf%mesh )
             CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_wfs_ae'), &
                  upf%gipaw_wfs_ae(1:upf%mesh,nb) )
             call add_attr( 'size', upf%mesh )
             CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_wfs_ps'),&
                  upf%gipaw_wfs_ps(1:upf%mesh,nb) )
             CALL xmlw_closetag ()
          END DO
          IF ( v2 ) CALL xmlw_closetag( )
          !
          ! Write all-electron and pseudo local potentials
          CALL xmlw_opentag( capitalize_if_v2('pp_gipaw_vlocal') )
          call add_attr( 'size', upf%mesh )
          CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_vlocal_ae'), &
               upf%gipaw_vlocal_ae(1:upf%mesh) )
          call add_attr( 'size', upf%mesh )
          CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_vlocal_ps'), &
               upf%gipaw_vlocal_ps(1:upf%mesh) )
          CALL xmlw_closetag ()
       END IF
       CALL xmlw_closetag () ! end pp_gipaw
    END IF
    !
  END SUBROUTINE write_pp_gipaw
  
END MODULE write_upf_new
  
