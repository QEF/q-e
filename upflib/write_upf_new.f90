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
    iun = xml_openfile ( filename )
    IF ( iun == -1 ) CALL upf_error('write_upf', 'cannot open file',1)
    !
    ! The starting line, header and info sections differ a lot
    ! between UPF v.2 and UPF with schema so we call different routines
    !
    IF ( v2 ) THEN
       !
       CALL xmlw_opentag ( 'UPF', 'version', upf%nv )
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
       CALL xmlw_writetag ( 'xml', '?', 'version,encoding','1.0,UTF-8')
       CALL xmlw_opentag ( 'qe_pp:pseudo', &
            'xsi:schemalocation,xmlns:xsi,xmlns:qe_pp', &
            QE_PP_URI//' '//QE_PP_URI//'.xsd,'//XSI//','//QE_PP_URI )
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
    IF(upf%nlcc) CALL xmlw_writetag( capitalize_if_v2('pp_nlcc'), &
         upf%rho_atc(1:upf%mesh), 'size', i2c(upf%mesh) )
    IF( .NOT. upf%tcoulombp) &
         CALL xmlw_writetag( capitalize_if_v2('pp_local'), &
         upf%vloc(1:upf%mesh), 'size', i2c(upf%mesh) )
    !
    CALL write_pp_semilocal ( upf )
    !
    CALL write_pp_nonlocal ( upf )
    !
    CALL write_pp_pswfc ( upf )
    !
    CALL write_pp_full_wfc ( upf )
    !
    CALL xmlw_writetag( capitalize_if_v2('pp_rhoatom'), &
         upf%rho_at(1:upf%mesh), 'size', i2c(upf%mesh) )
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
#include "version.h"
    INTEGER :: nw, nb
    !
    CALL xmlw_opentag ( 'pp_info' )
    CALL xmlw_writetag ( 'generated', xml_protect(upf%generated) )
    CALL xmlw_writetag ( 'creator', upf%author, &
         'NAME,VERSION', 'QE Atomic Code,'//version_number )
    CALL xmlw_writetag ( 'created', '', 'DATE', upf%date )
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
    IF (upf%tpawp .OR. upf%tvanp ) THEN
       CALL xmlw_writetag ( 'suggested_basis', '', &
            'ecutwfc,ecutrho', &
            r2c(upf%ecutwfc)//','//r2c(upf%ecutrho) )
    ELSE
       CALL xmlw_writetag ( 'suggested_basis', '', &
            'ecutwfc', r2c(upf%ecutwfc) )
    END IF
    DO nw =1, upf%nwfc
       IF( upf%oc(nw) >= 0.0_dp) THEN 
          CALL xmlw_opentag ( "valence_orbital", &
               'nl,pn,l',  &
               upf%els(nw)//','//i2c(upf%nchi(nw))//','//i2c(upf%lchi(nw)) )
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
    CHARACTER(LEN=:), ALLOCATABLE :: attr_list, attr_vals
    !
    attr_list="generated,author,date,comment,element,pseudo_type,relativistic," // &
         & "is_ultrasoft,is_paw,is_coulomb,has_so,has_wfc,has_gipaw," // &
         & "paw_as_gipaw,core_correction,functional,z_valence," // &
         & "total_psenergy,wfc_cutoff,rho_cutoff,l_max,l_max_rho,l_local," // &
         & "mesh_size,number_of_wfc,number_of_proj"

    attr_vals=xml_protect(upf%generated) //','// TRIM(upf%author) //','// &
         &TRIM(upf%date) //','// xml_protect(upf%comment) //','// &
         &TRIM(upf%psd) //','// TRIM(upf%typ) //','// TRIM(upf%rel) //','// &
         &l2c(upf%tvanp)//','// l2c(upf%tpawp)//','// l2c(upf%tcoulombp)//','//&
         &l2c(upf%has_so)//','// l2c(upf%has_wfc)//','//l2c(upf%has_gipaw) &
         &//','//l2c(upf%paw_as_gipaw) //','//l2c(upf%nlcc)//','// &
         &TRIM(upf%dft)//','//r2c(upf%zp)//','//r2c(upf%etotps)//','//&
         &r2c(upf%ecutwfc)//','//r2c(upf%ecutrho)//','//i2c(upf%lmax)//&
         & i2c(upf%lmax_rho)//i2c(upf%lloc)//i2c(upf%mesh)//i2c(upf%nwfc)&
         &//','//i2c(upf%nbeta)
    !
    CALL xmlw_writetag ( "PP_HEADER", '', attr_list, attr_vals )
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
          CALL xmlw_opentag ( 'input', 'program', 'ld1.x' )
       END IF
       REWIND (unit=u_input)
       read_write_loop: DO
          READ (u_input, '(A)',end=20,err=25) line
          WRITE (iun, '(A)') xml_protect(line)
          CYCLE read_write_loop
25        CALL upf_error('write_upf::write_inputfile', 'problem writing input data',-1)
20        EXIT read_write_loop
       END DO read_write_loop
    ELSE
       CALL upf_error('write_upf::write_inputfile', 'input file not open',-1)
    END IF
    CALL xmlw_closetag ( ) ! 'input'
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
       CALL xmlw_opentag( capitalize_if_v2('pp_mesh'), &
            'mesh,dx,xmin,rmax,zmesh', &
            & i2c(upf%mesh)//','//r2c(upf%dx)//','//r2c(upf%xmin)//','&
            & //r2c(upf%rmax)//','//r2c(upf%zmesh) )
    ELSE
       CALL xmlw_opentag( capitalize_if_v2('pp_mesh'), 'mesh', i2c(upf%mesh) )
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
             tag = 'PP_VNL.'//i2c(ind)
          ELSE
             tag = 'vnl'
          END IF
          IF ( upf%has_so ) THEN
             CALL xmlw_writetag( tag, upf%vnl(1:upf%mesh,l,ind), &
                  'l,j', i2c(l)//','//r2c(upf%jjj(nb)) )
          ELSE
             CALL xmlw_writetag( tag, upf%vnl(1:upf%mesh,l,ind), &
                  'l', i2c(l) )
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
       IF ( .NOT. v2 .AND. upf%has_so ) THEN
          CALL xmlw_writetag( tag, upf%beta(1:upf%mesh,nb), &
     &      'index,label,angular_momentum,tot_ang_mom,cutoff_radius_index,' &
     &      //'cutoff_radius,ultrasoft_cutoff_radius',&
     &      i2c(nb)//','//trim(upf%els_beta(nb))//',' &
     &               //i2c(upf%lll(nb))//','//r2c(upf%jjj(nb))//','&
     &               //i2c(upf%kbeta(nb))//','//r2c(upf%rcut(nb))//','&
     &               //r2c(upf%rcutus(nb)) )
       ELSE
          CALL xmlw_writetag( tag, upf%beta(1:upf%mesh,nb), &
     &      'index,label,angular_momentum,cutoff_radius_index,' &
     &      //'cutoff_radius,ultrasoft_cutoff_radius',&
     &      i2c(nb)//','//trim(upf%els_beta(nb))//',' &
     &               //i2c(upf%lll(nb))//','//i2c(upf%kbeta(nb))//',' &
     &               //r2c(upf%rcut(nb))//','//r2c(upf%rcutus(nb)) )
       END IF
    END DO
    !
    ! pp_dij (D_lm matrix)
    !
    CALL xmlw_opentag( capitalize_if_v2 ('pp_dij'), 'columns,rows', &
         i2c(upf%nbeta)//','//i2c(upf%nbeta) )
    WRITE(iun,*) upf%dion(1:upf%nbeta,1:upf%nbeta)
    CALL xmlw_closetag( ) 
    !
    ! pp_augmentation
    !
    IF (upf%tvanp .or. upf%tpawp) THEN
       CALL xmlw_opentag( capitalize_if_v2('pp_augmentation') )
       !
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
       !
       nb = upf%nbeta*upf%nbeta
       CALL xmlw_opentag( capitalize_if_v2('pp_q'), 'size', i2c(nb) )
       WRITE(iun,*) upf%qqq(1:upf%nbeta,1:upf%nbeta)
       CALL xmlw_closetag( )
       !
       IF ( upf%tpawp ) THEN
          WRITE(iun,"('<!--augmentation charge multipoles ( only for PAW) ',/,&
               &      'multipole array dims = (nbeta,nbeta,2*lmax+1)-->')")
          CALL xmlw_opentag( capitalize_if_v2('pp_multipoles'), 'nbeta,lmax', &
               i2c(upf%nbeta)//','//i2c(upf%lmax) )
          WRITE(iun,*) upf%paw%augmom(1:upf%nbeta,1:upf%nbeta,0:2*upf%lmax)
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
                   CALL xmlw_writetag( tag, upf%qfuncl(1:upf%mesh,nmb,l), &
                        'first_index,second_index,composite_index,angular_momentum,size', &
                        i2c(nb)//','//i2c(mb)//','//i2c(nmb)//','//i2c(l)//','//i2c(upf%mesh) )
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
                CALL xmlw_writetag( tag, upf%qfunc(1:upf%mesh,nmb), &
                     'size,first_index,second_index,composite_index', &
                     i2c(upf%mesh)//','//i2c(nb)//','//i2c(mb)//','//i2c(nmb) )
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
    CHARACTER(LEN=:), ALLOCATABLE :: attr_list, attr_vals
    !
    CALL xmlw_opentag( capitalize_if_v2('pp_pswfc') )
    DO nw =1, upf%nwfc
       attr_list = 'size,index,label,l'
       attr_vals = i2c(upf%mesh)//','//i2c(nw)//','//trim(upf%els(nw))// &
            & ','//i2c(upf%lchi(nw))
       IF ( upf%has_so) THEN
          attr_list = attr_list//',nn,jchi' 
          attr_vals = attr_vals//','//i2c(upf%nn(nw))//',' //r2c(upf%jchi(nw))
       END IF
       attr_list = attr_list//',occupation'
       attr_vals = attr_vals//','//r2c(upf%oc(nw))
       IF ( upf%nchi(nw) > upf%lchi(nw) ) THEN
          attr_list = attr_list//',n'
          attr_vals = attr_vals//','//i2c(upf%nchi(nw))
       END IF
       IF ( upf%epseu(nw) > 0.0_dp ) THEN
          attr_list = attr_list//',pseudo_energy'
          attr_vals = attr_vals//','//r2c(upf%epseu(nw))
       END IF
       IF ( upf%rcut_chi(nw) > 0.0_dp) THEN
          attr_list = attr_list//',cutoff_radius'
          attr_vals = attr_vals//','//r2c(upf%rcut_chi(nw))
       END IF
       IF ( upf%rcut_chi(nw) > 0.0_dp)THEN
          attr_list = attr_list//',ultrasoft_cutoff_radius' 
          attr_vals = attr_vals//','//r2c(upf%rcutus_chi(nw))
       END IF
       IF ( v2 ) THEN
          tag = 'PP_CHI.'//i2c(nw)
       ELSE
          tag = 'pp_chi'
       END IF
       CALL xmlw_writetag( tag, upf%chi(1:upf%mesh,nw), &
            attr_list, attr_vals )
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
       CALL xmlw_opentag( capitalize_if_v2('pp_full_wfc'), &
            'number_of_wfc', i2c(upf%nbeta) )
       !
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             tag = 'PP_AEWFC.'//i2c(nb)
          ELSE
             tag = 'pp_aewfc'
          END IF
          CALL xmlw_writetag( tag, upf%aewfc(1:upf%mesh,nb), 'index,label,l',&
               i2c(nb)//','//trim(upf%els_beta(nb))//','//i2c(upf%lll(nb)) )
       END DO
       !
       IF ( upf%has_so .AND. upf%tpawp ) THEN
          DO nb = 1, upf%nbeta
             IF ( v2 ) THEN
                tag = 'PP_AEWFC_rel.'//i2c(nb)
             ELSE
                tag = 'pp_aewfc_rel'
             END IF
             CALL xmlw_writetag(tag, upf%aewfc(1:upf%mesh,nb), 'index,label,l',&
                  i2c(nb)//','//trim(upf%els_beta(nb))//','//i2c(upf%lll(nb)) )
          END DO
       END IF
       !
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             tag = 'PP_PSWFC.'//i2c(nb)
          ELSE
             tag = 'pp_pswfc'
          END IF
          CALL xmlw_writetag( tag, upf%pswfc(1:upf%mesh,nb), &
              & 'size,index,label,l', i2c(upf%mesh)//','// &
              & i2c(nb)//','//trim(upf%els_beta(nb))//','//i2c(upf%lll(nb)) )
       END DO
       !
       CALL xmlw_closetag( )
       !
    END IF
    !
  END SUBROUTINE write_pp_full_wfc
  !--------------------------------------------------------
  SUBROUTINE write_pp_paw ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(IN) :: upf ! the pseudo data
    !
    IF ( upf%tpawp ) THEN
       CALL xmlw_opentag( capitalize_if_v2('pp_paw'), &
            'paw_data_format,core_energy', &
            i2c(upf%paw_data_format)//','//r2c(upf%paw%core_energy) ) 
       ! Full occupation (not only > 0 ones)
       CALL xmlw_writetag( capitalize_if_v2('pp_occupations'), &
            upf%paw%oc(1:upf%nbeta), 'size', i2c(upf%nbeta) )
       ! All-electron core charge
       CALL xmlw_writetag( capitalize_if_v2('pp_ae_nlcc'), &
            upf%paw%ae_rho_atc(1:upf%mesh), 'size', i2c(upf%mesh) )
       ! All-electron local potential
       CALL xmlw_writetag( capitalize_if_v2('pp_ae_vloc'), &
            upf%paw%ae_vloc(1:upf%mesh), 'size', i2c(upf%mesh) )
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
       CALL xmlw_opentag( capitalize_if_v2('pp_gipaw'), &
            'gipaw_data_format', i2c(upf%gipaw_data_format) ) 
       IF ( v2 ) CALL xmlw_opentag( 'PP_GIPAW_CORE_ORBITALS', &
            'number_of_core_orbitals', i2c(upf%gipaw_ncore_orbitals) )
       DO nb = 1,upf%gipaw_ncore_orbitals
          IF ( v2 ) THEN
             tag = "PP_GIPAW_CORE_ORBITAL."//i2c(nb)
          ELSE
             tag = 'pp_gipaw_core_orbital'
          END IF
          CALL xmlw_writetag( tag, &
               & upf%gipaw_core_orbital(1:upf%mesh,nb), &
               & 'size,index,label,n,l', &
               & i2c(upf%mesh) // ',' // i2c(nb) // ',' // &
               & trim(upf%gipaw_core_orbital_el(nb)) // ',' // & 
               & r2c(upf%gipaw_core_orbital_n(nb)) // ',' // & 
               & r2c(upf%gipaw_core_orbital_l(nb)) )
       END DO
       IF ( v2 ) CALL xmlw_closetag ( )
       ! Only write core orbitals in the PAW as GIPAW case
       IF ( .NOT. upf%paw_as_gipaw) THEN
          !
          ! Write valence all-electron and pseudo orbitals
          !
          IF ( v2 ) CALL xmlw_opentag( 'PP_GIPAW_ORBITALS', &
               'number_of_valence_orbitals', i2c(upf%gipaw_wfs_nchannels) )
          DO nb = 1,upf%gipaw_wfs_nchannels
             IF ( v2 ) THEN
                tag = "PP_GIPAW_ORBITAL."//i2c(nb)
             ELSE
                tag = 'pp_gipaw_orbital'
             END IF
             CALL xmlw_opentag( tag, &
               & 'index,label,l,cutoff_radius,ultrasoft_cutoff_radius', &
               & i2c(nb) // ',' // &
               & trim(upf%gipaw_wfs_el(nb)) // ',' // & 
               & i2c(upf%gipaw_wfs_ll(nb)) // ',' // &
               & r2c(upf%gipaw_wfs_rcut(nb)) // ',' // &
               & r2c(upf%gipaw_wfs_rcutus(nb)) )
             !
             CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_wfs_ae'), &
                  upf%gipaw_wfs_ae(1:upf%mesh,nb), 'size', i2c(upf%mesh) )
             CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_wfs_ps'),&
                  upf%gipaw_wfs_ps(1:upf%mesh,nb), 'size', i2c(upf%mesh) )
             CALL xmlw_closetag ()
          END DO
          IF ( v2 ) CALL xmlw_closetag( )
          !
          ! Write all-electron and pseudo local potentials
          CALL xmlw_opentag( capitalize_if_v2('pp_gipaw_vlocal') )
          CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_vlocal_ae'), &
               upf%gipaw_vlocal_ae(1:upf%mesh), 'size', i2c(upf%mesh) )
          CALL xmlw_writetag( capitalize_if_v2('pp_gipaw_vlocal_ps'), &
               upf%gipaw_vlocal_ps(1:upf%mesh), 'size', i2c(upf%mesh) )
          CALL xmlw_closetag ()
       END IF
       CALL xmlw_closetag () ! end pp_gipaw
    END IF
    !
  END SUBROUTINE write_pp_gipaw
  
END MODULE write_upf_new
  
