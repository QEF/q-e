!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
MODULE  read_upf_new__
  !-----------------------------------------------------
  !! this module contains the simplified code for reading
  !! pseudopotential files in either UPF v.2 or xml
  !
  USE xmltools
  USE upf_kinds, ONLY: dp
  USE pseudo_types, ONLY: pseudo_upf, pseudo_config
  !
  LOGICAL :: v2
  !! true if UPF v.2 version, false if new UPF with xml schema
  INTEGER :: iun
  !! unit for reading data
  !
  PUBLIC
  !
CONTAINS
  !
  !------------------------------------------------+
  SUBROUTINE read_upf_new_ (filename, upf, ierr)         !
    !---------------------------------------------+
    !! Reads pseudopotential in UPF format (either v.2 or upf_schema).
    !! Derived-type variable *upf* store in output the data read from file. 
    !! File *filename* is opened and closed inside the routine
    !
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: filename  
    !! i/o filename
    TYPE(pseudo_upf),INTENT(OUT) :: upf
    !! the derived type storing the pseudo data
    INTEGER, INTENT(OUT) :: ierr
    !! ierr=0  : xml schema, ierr=-2: UPF v.2
    !! ierr=-81: error reading PP file
    !
    iun = xml_openfile ( filename )
    IF ( iun == -1 ) CALL upf_error('read_upf', 'cannot open file',1)
    call xmlr_opentag ( 'qe_pp:pseudo', IERR = ierr )
    if ( ierr == 0 ) then
       v2 =.false.
    else if ( ierr == -1 ) then
       rewind (iun) 
       call xmlr_opentag ( 'UPF', IERR = ierr )
       if ( ierr == 0 ) then
          v2 =.true.
          ierr = -2
          CALL get_attr ( 'version', upf%nv )
       end if
    end if
    if ( ierr /= 0 .and. ierr /= -2 ) then
       ierr = -81
       return
    end if
    !
    ! The header sections differ a lot between UPF v.2 and UPF with schema
    !
    IF ( v2 ) THEN
       CALL read_pp_header_v2 ( upf )
    ELSE
       CALL read_pp_header_schema ( upf )
    END IF
    !
    ! From here on the format of v2 and schema do not differ much:
    ! the most frequent difference is capitalization of tags
    ! (see function capitalize_if_v2)
    !
    CALL read_pp_mesh ( upf )
    !
    IF(upf%nlcc) then
       allocate ( upf%rho_atc(1:upf%mesh) )
       CALL xmlr_readtag( capitalize_if_v2('pp_nlcc'), &
            upf%rho_atc(1:upf%mesh) )
    end if
    IF( .NOT. upf%tcoulombp) then
       allocate ( upf%vloc(1:upf%mesh) )
       CALL xmlr_readtag( capitalize_if_v2('pp_local'), &
            upf%vloc(1:upf%mesh) )
    end if
    !
    CALL read_pp_semilocal ( upf )
    !
    !CALL read_pp_nonlocal ( upf )
    !
    !CALL read_pp_pswfc ( upf )
    !
    !CALL read_pp_full_wfc ( upf )
    !
    allocate( upf%rho_at(1:upf%mesh) )
    CALL xmlr_readtag( capitalize_if_v2('pp_rhoatom'), &
         upf%rho_at(1:upf%mesh) )
    !
    !CALL read_pp_paw ( upf )
    !
    !CALL read_pp_gipaw ( upf )
    !
    ! close initial tag, qe_pp:pseudo or UPF
    !
    CALL xmlr_closetag ( )
    !
    CALL xml_closefile ( )
    !
  END SUBROUTINE read_upf_new_
  !--------------------------------------------------------
  SUBROUTINE read_pp_header_schema ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf), INTENT(INOUT) :: upf ! the pseudo data
    !
    CALL xmlr_opentag( capitalize_if_v2('pp_header') )
    !
    CALL xmlr_readtag( 'element', upf%psd )
    CALL xmlr_readtag( 'z_valence', upf%zp )
    CALL xmlr_readtag( 'type', upf%typ )
    CALL xmlr_readtag( 'functional', upf%dft )
    CALL xmlr_readtag( 'relativistic', upf%rel )
    CALL xmlr_readtag( 'is_ultrasoft', upf%tvanp )
    CALL xmlr_readtag( 'is_paw', upf%tpawp )
    CALL xmlr_readtag( 'is_coulomb', upf%tcoulombp )
    CALL xmlr_readtag( 'has_so', upf%has_so )
    CALL xmlr_readtag( 'has_wfc', upf%has_wfc )
    CALL xmlr_readtag( 'has_gipaw', upf%has_gipaw )
    CALL xmlr_readtag( 'paw_as_gipaw', upf%paw_as_gipaw)
    CALL xmlr_readtag( 'core_correction', upf%nlcc)
    CALL xmlr_readtag( 'total_psenergy', upf%etotps )
    CALL xmlr_readtag( 'wfc_cutoff', upf%ecutwfc )
    CALL xmlr_readtag( 'rho_cutoff', upf%ecutrho )
    CALL xmlr_readtag( 'l_max', upf%lmax )
    CALL xmlr_readtag( 'l_max_rho', upf%lmax_rho )
    CALL xmlr_readtag( 'l_local', upf%lloc )
    CALL xmlr_readtag( 'mesh_size', upf%mesh )
    CALL xmlr_readtag( 'number_of_wfc', upf%nwfc )
    CALL xmlr_readtag( 'number_of_proj', upf%nbeta )
    !
    CALL xmlr_closetag( )
    !
  END SUBROUTINE read_pp_header_schema
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_header_v2 ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf), INTENT(INOUT) :: upf ! the pseudo data
    !
    CHARACTER(LEN=1) :: dummy
    !
    CALL xmlr_readtag ( capitalize_if_v2('pp_header'), dummy )
    CALL get_attr ('generated', upf%generated)
    CALL get_attr ('author', upf%author)
    CALL get_attr ('date', upf%date)
    CALL get_attr ('comment', upf%comment)
    CALL get_attr ('element', upf%psd)
    CALL get_attr ('pseudo_type', upf%typ)
    CALL get_attr ('relativistic', upf%rel)
    CALL get_attr ('is_ultrasoft', upf%tvanp)
    CALL get_attr ('is_paw', upf%tpawp)
    CALL get_attr ('is_coulomb', upf%tcoulombp)
    CALL get_attr ('has_so', upf%has_so)
    CALL get_attr ('has_wfc', upf%has_wfc)
    CALL get_attr ('has_gipaw', upf%has_gipaw)
    CALL get_attr ('paw_as_gipaw', upf%paw_as_gipaw)
    CALL get_attr ('core_correction', upf%nlcc)
    CALL get_attr ('functional', upf%dft)
    CALL get_attr ('z_valence', upf%zp)
    CALL get_attr ('total_psenergy', upf%etotps)
    CALL get_attr ('wfc_cutoff', upf%ecutwfc)
    CALL get_attr ('rho_cutoff', upf%ecutrho)
    CALL get_attr ('l_max', upf%lmax)
    CALL get_attr ('l_max_rho', upf%lmax_rho)
    CALL get_attr ('l_local', upf%lloc)
    CALL get_attr ('mesh_size', upf%mesh)
    CALL get_attr ('number_of_wfc', upf%nwfc)
    CALL get_attr ('number_of_proj', upf%nbeta )
    !
  END SUBROUTINE read_pp_header_v2
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_mesh ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    integer :: mesh
    !
    CALL xmlr_opentag( capitalize_if_v2('pp_mesh') )
    CALL get_attr ( 'mesh', mesh )
    if ( mesh /= upf%mesh ) call upf_error('read_pp_mesh','mismatch in mesh',mesh)
    CALL get_attr ( 'dx'  , upf%dx   )
    CALL get_attr ( 'xmin', upf%xmin )
    CALL get_attr ( 'rmax', upf%rmax )
    CALL get_attr ( 'zmesh', upf%zmesh )
    allocate ( upf%r(1:upf%mesh) )
    CALL xmlr_readtag( capitalize_if_v2('pp_r'), upf%r(1:upf%mesh) )
    allocate ( upf%rab(1:upf%mesh) )
    CALL xmlr_readtag( capitalize_if_v2('pp_rab'), upf%rab(1:upf%mesh) )
    !
    CALL xmlr_closetag( ) ! end pp_mesh
    !
  END SUBROUTINE read_pp_mesh
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_semilocal ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    !
    INTEGER :: nb, ind, l, j, ierr
    CHARACTER(LEN=8) :: tag
    real(dp), allocatable :: vnl(:)
    !
    IF ( upf%typ == "SL" ) THEN
       !
       IF ( upf%has_so ) then
          ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,2))
       else
          ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,1))
       end if
       allocate ( vnl(1:upf%mesh) )
       CALL xmlr_opentag( capitalize_if_v2('pp_semilocal') )       
       !
       IF ( v2 ) THEN
          tag = 'PP_VNL.1'
       ELSE
          tag = 'vnl'
       END IF
       DO nb = 1,upf%nbeta
          CALL xmlr_readtag( tag, vnl, ierr )
          if ( ierr /= 0 ) then
             if ( v2 ) then
                go to 10
             else
                call upf_error('read_pp_semilocal','error reading SL PPs',1)
             end if
          end if
          CALL get_attr ( 'l', l)
          ind = 1
          IF ( upf%has_so ) then
             CALL get_attr ( 'j', j)
             IF ( l > 0 .AND. ABS(j-l-0.5_dp) < 0.001_dp ) ind = 2
             if ( v2 .and. ind == 2 ) &
                  call upf_error('read_pp_semilocal','inconsistency in SL',1)
          END IF
          upf%vnl(:,l,ind) = vnl(:)
       END DO
       !
       CALL xmlr_closetag( ) ! end pp_semilocal
       !
10     IF ( v2 .and. upf%has_so ) then
          rewind ( iun )
          CALL xmlr_opentag( capitalize_if_v2('pp_semilocal') )
          ind = 2
          tag = 'PP_VNL.2'
          DO nb = 1,upf%nbeta
             CALL xmlr_readtag( tag, vnl, ierr )
             if ( ierr /= 0 ) exit
             CALL get_attr ( 'l', l)
             CALL get_attr ( 'j', j)
             IF ( .not. (l > 0 .AND. ABS(j-l-0.5_dp) < 0.001_dp) ) ind = 1
             if ( v2 .and. ind == 1 ) &
                  call upf_error('read_pp_semilocal','inconsistency in SL',2)
             upf%vnl(:,l,ind) = vnl(:)
          END DO
          CALL xmlr_closetag( ) ! end pp_semilocal
       END IF
       deallocate ( vnl )
    END IF
    !
  END SUBROUTINE read_pp_semilocal
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_pswfc ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    !
    INTEGER :: nw, ind, l
    CHARACTER(LEN=8) :: tag
    !
    allocate ( upf%chi(1:upf%mesh,upf%nwfc) ) 
    CALL xmlr_opentag( capitalize_if_v2('pp_pswfc') )
    DO nw=1,upf%nwfc
       IF ( v2 ) THEN
          tag = 'PP_CHI.'//i2c(nw)
       ELSE
          tag = 'pp_chi'
       END IF
       CALL xmlr_readtag( tag, upf%chi(1:upf%mesh,nw) )
       call get_attr('index', ind)
       if ( ind /= nw ) &
            call upf_error('read_pp_pswfc','mismatch reading PSWFC', nw)
       call get_attr( 'label', upf%els(nw) )
       call get_attr( 'l', upf%lchi(nw) )
       IF ( upf%has_so) THEN
          call get_attr( 'nn', upf%nn(nw) )
          call get_attr( 'jchi', upf%jchi(nw) )
       END IF
       call get_attr( 'occupation', upf%oc(nw) )
       call get_attr( 'n', upf%nchi(nw) )
       call get_attr( 'pseudo_energy', upf%epseu(nw) )
       call get_attr( 'cutoff_radius', upf%rcut_chi(nw) )
       call get_attr( 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw) )
    END DO
    CALL xmlw_closetag( ) ! end pp_pswfc
    !
  END SUBROUTINE read_pp_pswfc
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
END MODULE read_upf_new__
