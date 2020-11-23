!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
MODULE  read_upf_new_module
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
  SUBROUTINE read_upf_new (filename, upf, ierr)         !
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
       call xml_closefile( )
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
    ! compatibility
    upf%is_gth = .false.
    upf%is_multiproj = .true.
    !
    ! From here on the format of v2 and schema do not differ much:
    ! the most frequent difference is capitalization of tags
    ! (see function capitalize_if_v2)
    !
    CALL read_pp_mesh ( upf )
    !
    allocate ( upf%rho_atc(upf%mesh) )
    IF(upf%nlcc) then
       CALL xmlr_readtag( capitalize_if_v2('pp_nlcc'), &
            upf%rho_atc(:) )
    else
       upf%rho_atc(:) = 0.0_dp
    end if
    IF( .NOT. upf%tcoulombp) then
       allocate ( upf%vloc(upf%mesh) )
       CALL xmlr_readtag( capitalize_if_v2('pp_local'), &
            upf%vloc(:), ierr )
       !
       ! existing PP files may have pp_nlcc first, pp_local later,
       ! but also the other way round - check that everything was right
       !
       if ( ierr /= 0 ) then
          ierr = -81
          return
       end if
    end if
    !
    CALL read_pp_semilocal ( upf )
    !
    CALL read_pp_nonlocal ( upf )
    !
    CALL read_pp_pswfc ( upf )
    !
    CALL read_pp_full_wfc ( upf )
    !
    allocate( upf%rho_at(1:upf%mesh) )
    CALL xmlr_readtag( capitalize_if_v2('pp_rhoatom'), &
         upf%rho_at(1:upf%mesh) )
    !
    CALL read_pp_spinorb ( upf )
    !
    CALL read_pp_paw ( upf )
    !
    CALL read_pp_gipaw ( upf )
    !
    ! close initial tag, qe_pp:pseudo or UPF
    !
    CALL xmlr_closetag ( )
    !
    CALL xml_closefile ( )
    !
  END SUBROUTINE read_upf_new
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
    if ( mesh == 0 ) THEN
       call upf_error('read_pp_mesh',&
         'mesh size missing, using the one in header',-1)
    else if ( mesh /= upf%mesh ) THEN
       call upf_error('read_pp_mesh',&
         'mismatch in mesh size, discarding the one in header',-1)
       upf%mesh = mesh
    end if
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
       tag = 'vnl'
       DO nb = 1,upf%nbeta
          IF ( v2 ) THEN
             ! NOTA BENE: v2 format follows available PP files, written 
             ! using original write_upf_v2; not FoX-based write_upf_v2
             IF ( nb - 1 == upf%lloc ) CYCLE
             tag = 'PP_VNL.'//i2c(nb-1)
          END IF
          CALL xmlr_readtag( tag, vnl, ierr )
          if ( ierr /= 0 ) &
               call upf_error('read_pp_semilocal','error reading SL PPs',1)
          CALL get_attr ( 'l', l)
          ind = 1
          IF ( upf%has_so ) then
             CALL get_attr ( 'j', j)
             IF ( l > 0 .AND. ABS(j-l-0.5_dp) < 0.001_dp ) ind = 2
             ! FIXME: what about spin-orbit case for v.2 upf?
             if ( v2 ) &
                  call upf_error('read_pp_semilocal','check spin-orbit',1)
          END IF
          upf%vnl(:,l,ind) = vnl(:)
       END DO
       deallocate ( vnl )
       !
       CALL xmlr_closetag( ) ! end pp_semilocal
       !
    END IF
    !
  END SUBROUTINE read_pp_semilocal
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_nonlocal ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    !
    LOGICAL :: isnull
    INTEGER :: nb, ind, l, l_, ln, lm, mb, nmb
    CHARACTER(LEN=15) :: tag
    REAL(dp), ALLOCATABLE :: aux(:)
    !
    nb = upf%nbeta
    IF ( nb == 0 ) nb = 1
    ALLOCATE (upf%beta(upf%mesh,nb) )
    ALLOCATE (upf%els_beta(nb), &
              upf%lll(nb),      &
              upf%kbeta(nb),    &
              upf%rcut(nb),     &
              upf%rcutus(nb),   &
              upf%dion(nb,nb),  &
              upf%qqq(nb,nb)    )
    !
    IF (upf%has_so) ALLOCATE( upf%jjj(upf%nbeta)) 
    !
    IF ( upf%nbeta == 0 ) THEN
       upf%nqf = 0
       upf%nqlc= 0
       upf%kkbeta = 0
       upf%qqq_eps=-1.0_dp
       RETURN
    END IF
    !
    CALL xmlr_opentag( capitalize_if_v2('pp_nonlocal') )
    !
    DO nb = 1,upf%nbeta
       !
       IF ( v2 ) THEN
          tag = 'PP_BETA.'//i2c(nb)
       ELSE
          tag = 'pp_beta'
       END IF
       CALL xmlr_readtag( tag, upf%beta(1:upf%mesh,nb) )
       CALL get_attr('index', mb)
       ! not-so-strict test: index is absent or incorrect in some UPF v.2 files
       IF ( .NOT. v2 .AND. nb /= mb ) &
            CALL upf_error('read_pp_nonlocal','mismatch',nb)
       CALL get_attr('label', upf%els_beta(nb))
       CALL get_attr('angular_momentum', upf%lll(nb))
       IF ( .NOT. v2 .AND. upf%has_so ) &
            CALL get_attr('tot_ang_mom', upf%jjj(nb))
       CALL get_attr('cutoff_radius_index', upf%kbeta(nb))
       CALL get_attr('cutoff_radius', upf%rcut(nb))
       CALL get_attr('ultrasoft_cutoff_radius', upf%rcutus(nb))
       !
    END DO
    !
    ! pp_dij (D_lm matrix)
    !
    CALL xmlr_opentag( capitalize_if_v2 ('pp_dij') )
    READ(iun,*) upf%dion(1:upf%nbeta,1:upf%nbeta)
    CALL xmlr_closetag( ) 
    !
    ! pp_augmentation
    !
    IF (upf%tvanp .or. upf%tpawp) THEN
       CALL xmlr_opentag( capitalize_if_v2('pp_augmentation') )
       !
       IF ( v2 ) THEN
          CALL get_attr ( 'q_with_l', upf%q_with_l )
          CALL get_attr ( 'nqf', upf%nqf )
          CALL get_attr ( 'nqlc', upf%nqlc )
          IF (upf%tpawp) THEN
             CALL get_attr ( 'shape', upf%paw%augshape )
             CALL get_attr ( 'cutoff_r', upf%paw%raug )
             CALL get_attr ( 'cutoff_r_index', upf%paw%iraug )
             CALL get_attr ( 'augmentation_epsilon', upf%qqq_eps )
             CALL get_attr ( 'l_max_aug', upf%paw%lmax_aug )
          ENDIF
       ELSE
          CALL xmlr_readtag( 'q_with_l', upf%q_with_l )
          CALL xmlr_readtag( 'nqf', upf%nqf )
          CALL xmlr_readtag( 'nqlc', upf%nqlc )
          IF (upf%tpawp) THEN
             CALL xmlr_readtag( 'shape', upf%paw%augshape )
             CALL xmlr_readtag( 'cutoff_r', upf%paw%raug )
             CALL xmlr_readtag( 'cutoff_r_index', upf%paw%iraug )
             CALL xmlr_readtag( 'augmentation_epsilon', upf%qqq_eps )
             CALL xmlr_readtag( 'l_max_aug', upf%paw%lmax_aug )
          ENDIF
       ENDIF
       !
       CALL xmlr_opentag( capitalize_if_v2('pp_q') )
       READ(iun,*) upf%qqq(1:upf%nbeta,1:upf%nbeta)
       CALL xmlr_closetag( )
       !
       IF ( upf%tpawp ) THEN
          CALL xmlr_opentag( capitalize_if_v2('pp_multipoles') )
          ALLOCATE ( upf%paw%augmom(1:upf%nbeta,1:upf%nbeta,0:2*upf%lmax) )
          READ(iun,*) upf%paw%augmom(1:upf%nbeta,1:upf%nbeta,0:2*upf%lmax)
          CALL xmlr_closetag ()
       ENDIF
       !
       ! read polinomial coefficients for Q_ij expansion at small radius
       !
       IF ( upf%nqlc == 0 ) upf%nqlc = 2*upf%lmax+1
       ALLOCATE( upf%rinner( upf%nqlc ) )
       IF ( v2 .AND. upf%nqf > 0) THEN
          ALLOCATE ( upf%qfcoef(upf%nqf, upf%nqlc, upf%nbeta, upf%nbeta) )
          CALL xmlr_opentag('PP_QFCOEF')
          READ(iun,*) upf%qfcoef
          CALL xmlr_closetag ()
          CALL xmlr_readtag('PP_RINNER',upf%rinner)
       ELSE IF ( upf%nqf == 0 ) THEN
          ALLOCATE( upf%qfcoef(1,1,1,1) )
          upf%qfcoef =0.0_dp
       ENDIF
       !
       ! Read augmentation charge Q_ij
       !
       IF( upf%q_with_l ) THEN
          ALLOCATE( upf%qfuncl(upf%mesh,upf%nbeta*(upf%nbeta+1)/2,0:2*upf%lmax) )
          upf%qfuncl(:,:,:) = 0.0_dp
          ! NOTE: it would be wiser to dimension qfuncl as (:,:,0:upf%lmax)
          ! and store the q_l(r) with index l=L/2 (see loop_on_l below)
          ! This would save some storage and avoid "holes" in the array
          ! that may be a source of trouble if not initialized to zero 
       ELSE
          ALLOCATE ( upf%qfunc(upf%mesh,upf%nbeta*(upf%nbeta+1)/2) )
          upf%qfunc (:,:) = 0.0_dp
       END IF
       ALLOCATE ( aux(upf%mesh) )
       loop_on_nb: DO nb = 1,upf%nbeta
          ln = upf%lll(nb)
          loop_on_mb: DO mb = nb,upf%nbeta
             lm = upf%lll(mb)
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
                   CALL xmlr_readtag( tag, aux )
                   CALL get_attr ('composite_index', nmb)
                   IF ( nmb /= mb*(mb-1)/2 + nb ) &
                        CALL upf_error ('read_pp_nonlocal','mismatch',1)
                   CALL get_attr ('angular_momentum', l_)
                   IF ( l /= l_ ) CALL upf_error ('read_pp_nonlocal','mismatch',2)                 
                   upf%qfuncl(:,nmb,l) = aux(:)
                   IF (upf%tpawp) upf%qfuncl(upf%paw%iraug+1:,nmb,l) = 0._DP
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
                CALL xmlr_readtag( tag, aux )
                CALL get_attr ('composite_index', nmb)
                IF ( nmb /= mb*(mb-1)/2 + nb ) &
                     CALL upf_error ('read_pp_nonlocal','mismatch',3)
                upf%qfunc(:,nmb) = aux(:)
                !
             ENDIF
          ENDDO loop_on_mb
       ENDDO  loop_on_nb
       !
       DEALLOCATE (aux)
       CALL xmlr_closetag( ) ! end pp_augmentation
       !
    END IF
    CALL xmlr_closetag( ) ! end pp_nonlocal
    !
    ! Maximum radius of beta projector: outer radius to integrate
    upf%kkbeta = MAXVAL(upf%kbeta(1:upf%nbeta))
    ! For PAW, augmentation charge may extend a bit further:
    IF(upf%tpawp) upf%kkbeta = MAX(upf%kkbeta, upf%paw%iraug)
    !
  END SUBROUTINE read_pp_nonlocal
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
    allocate ( upf%els(upf%nwfc), &
                upf%oc(upf%nwfc), &
                upf%lchi(upf%nwfc), &
                upf%nchi(upf%nwfc), &
                upf%rcut_chi(upf%nwfc), &
                upf%rcutus_chi(upf%nwfc), &
                upf%epseu(upf%nwfc) )
    IF ( upf%has_so ) THEN
       allocate ( upf%nn(upf%nwfc) )
       allocate ( upf%jchi(upf%nwfc) )
    END IF
    !
    CALL xmlr_opentag( capitalize_if_v2('pp_pswfc') )
    DO nw=1,upf%nwfc
       IF ( v2 ) THEN
          tag = 'PP_CHI.'//i2c(nw)
       ELSE
          tag = 'pp_chi'
       END IF
       CALL xmlr_readtag( tag, upf%chi(1:upf%mesh,nw) )
       call get_attr('index', ind)
       ! not-so-strict test: index is absent or incorrect in some UPF v.2 files
       if ( .NOT. v2 .AND. ind /= nw ) &
            call upf_error('read_pp_pswfc','mismatch reading PSWFC', nw)
       call get_attr( 'label', upf%els(nw) )
       call get_attr( 'l', upf%lchi(nw) )
       IF ( .not. v2 .and. upf%has_so ) THEN
          call get_attr( 'nn', upf%nn(nw) )
          call get_attr( 'jchi', upf%jchi(nw) )
       END IF
       call get_attr( 'occupation', upf%oc(nw) )
       call get_attr( 'n', upf%nchi(nw) )
       call get_attr( 'pseudo_energy', upf%epseu(nw) )
       call get_attr( 'cutoff_radius', upf%rcut_chi(nw) )
       call get_attr( 'ultrasoft_cutoff_radius', upf%rcutus_chi(nw) )
    END DO
    CALL xmlr_closetag( ) ! end pp_pswfc
    !
  END SUBROUTINE read_pp_pswfc
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_full_wfc ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    !
    INTEGER :: nb, mb
    CHARACTER(LEN=15) :: tag
    !
    IF ( upf%has_wfc ) THEN
       !
       ALLOCATE (upf%aewfc(1:upf%mesh,upf%nbeta) )
       CALL xmlr_opentag( capitalize_if_v2('pp_full_wfc') )
       !
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             tag = 'PP_AEWFC.'//i2c(nb)
          ELSE
             tag = 'pp_aewfc'
          END IF
          CALL xmlr_readtag( tag, upf%aewfc(1:upf%mesh,nb) )
          CALL get_attr ('index',mb)
          ! not-so-strict test (and two more below):
          ! index may be absent or incorrect in some UPF v.2 files
          IF ( .NOT. v2 .AND. nb /= mb ) CALL upf_error('read_pp_full_wfc','mismatch',1)
       END DO
       !
       IF ( upf%has_so .AND. upf%tpawp ) THEN
          ALLOCATE (upf%paw%aewfc_rel(1:upf%mesh,upf%nbeta) )
          DO nb = 1, upf%nbeta
             IF ( v2 ) THEN
                tag = 'PP_AEWFC_rel.'//i2c(nb)
             ELSE
                tag = 'pp_aewfc_rel'
             END IF
             CALL xmlr_readtag(tag, upf%paw%aewfc_rel(1:upf%mesh,nb) )
             CALL get_attr ('index',mb)
             IF ( .NOT. v2 .AND. nb /= mb ) CALL upf_error('read_pp_full_wfc','mismatch',2)
          END DO
       END IF
       !
       ALLOCATE (upf%pswfc(1:upf%mesh,upf%nbeta) )
       DO nb = 1, upf%nbeta
          IF ( v2 ) THEN
             tag = 'PP_PSWFC.'//i2c(nb)
          ELSE
             tag = 'pp_pswfc'
          END IF
          CALL xmlr_readtag(tag, upf%pswfc(1:upf%mesh,nb) )
          CALL get_attr ('index',mb)
          IF ( .NOT. v2 .AND. nb /= mb ) CALL upf_error('read_pp_full_wfc','mismatch',3)
       END DO
       !
       CALL xmlr_closetag( )
       !
    END IF
    !
  END SUBROUTINE read_pp_full_wfc
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_spinorb ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER :: nw, nb, ierr
    CHARACTER(LEN=1) :: dummy
    !
    IF ( .NOT. v2 .OR. .NOT. upf%has_so ) RETURN
    !
    CALL xmlr_opentag( 'PP_SPIN_ORB' )
    DO nw = 1,upf%nwfc
       CALL xmlr_readtag( 'PP_RELWFC.'//i2c(nw), dummy )
       CALL get_attr( 'index' , nb )
       ! not-so-strict test: index absent or incorrect in some UPF v.2 files
       IF ( .NOT. v2 .AND. nb /= nw ) CALL upf_error('read_pp_spinorb','mismatch',1)
       CALL get_attr( 'els',   upf%els(nw) )
       CALL get_attr( 'nn',    upf%nn(nw) )
       CALL get_attr( 'lchi',  upf%lchi(nw) )
       CALL get_attr( 'jchi',  upf%jchi(nw) )
       CALL get_attr( 'oc',    upf%oc(nw) )
    ENDDO
    !
    DO nb = 1,upf%nbeta
       CALL xmlr_readtag( 'PP_RELBETA.'//i2c(nb), dummy, ierr )
       !
       ! existing PP files may have pp_relbeta first, pp_relwfc later,
       ! but also the other way round - check that everything was right
       !
       if ( ierr /= 0 ) then
          ierr = -81
          return
       end if
       CALL get_attr( 'index' , nw )
       IF ( .NOT.v2 .AND. nb /= nw ) CALL upf_error('read_pp_spinorb','mismatch',2)
       CALL get_attr( 'lll',  upf%lll(nb) )
       CALL get_attr( 'jjj',  upf%jjj(nb) )
    ENDDO
    CALL xmlr_closetag () ! end pp_spin_orb
    !
  END SUBROUTINE read_pp_spinorb
  !
  !--------------------------------------------------------
  SUBROUTINE read_pp_paw ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    INTEGER :: nb, mb
    !
    IF ( .NOT. upf%tpawp ) RETURN
    !
    CALL xmlr_opentag( capitalize_if_v2('pp_paw') )
    CALL get_attr ('paw_data_format', upf%paw_data_format)
    CALL get_attr ('core_energy', upf%paw%core_energy) 
    ! Full occupation (not only > 0 ones)
    ALLOCATE (upf%paw%oc(upf%nbeta) )
    ALLOCATE (upf%paw%ae_rho_atc(upf%mesh) )
    ALLOCATE (upf%paw%ae_vloc(upf%mesh) )
    CALL xmlr_readtag( capitalize_if_v2('pp_occupations'), &
         upf%paw%oc(1:upf%nbeta) )
    ! All-electron core charge
    CALL xmlr_readtag( capitalize_if_v2('pp_ae_nlcc'), &
         upf%paw%ae_rho_atc(1:upf%mesh) )
    ! All-electron local potential
    CALL xmlr_readtag( capitalize_if_v2('pp_ae_vloc'), &
         upf%paw%ae_vloc(1:upf%mesh) )
    CALL xmlr_closetag () ! end pp_paw
    !
    ALLOCATE(upf%paw%pfunc(upf%mesh, upf%nbeta,upf%nbeta) )
    upf%paw%pfunc(:,:,:) = 0._dp
    IF (upf%has_so) THEN
       ALLOCATE(upf%paw%pfunc_rel(upf%mesh, upf%nbeta,upf%nbeta) )
       upf%paw%pfunc_rel(:,:,:) = 0._dp
    ENDIF
    DO nb=1,upf%nbeta
       DO mb=1,nb
          upf%paw%pfunc (1:upf%mesh, nb, mb) = &
               upf%aewfc(1:upf%mesh, nb) * upf%aewfc(1:upf%mesh, mb)
          IF (upf%has_so) THEN
             upf%paw%pfunc_rel (1:upf%paw%iraug, nb, mb) =  &
                  upf%paw%aewfc_rel(1:upf%paw%iraug, nb) *   &
                  upf%paw%aewfc_rel(1:upf%paw%iraug, mb)
!
!    The small component is added to pfunc. pfunc_rel is useful only
!    to add a small magnetic contribution
!
             upf%paw%pfunc (1:upf%paw%iraug, nb, mb) = &
                        upf%paw%pfunc (1:upf%paw%iraug, nb, mb) + &
                        upf%paw%pfunc_rel (1:upf%paw%iraug, nb, mb)
          ENDIF
          upf%paw%pfunc(upf%paw%iraug+1:,nb,mb) = 0._dp
          !
          upf%paw%pfunc (1:upf%mesh, mb, nb) = upf%paw%pfunc (1:upf%mesh, nb, mb)
          IF (upf%has_so) upf%paw%pfunc_rel (1:upf%mesh, mb, nb) =  &
               upf%paw%pfunc_rel (1:upf%mesh, nb, mb)
       ENDDO
    ENDDO
    !
    ! Pseudo wavefunctions (not only the ones for oc > 0)
    ! All-electron wavefunctions
    ALLOCATE(upf%paw%ptfunc(upf%mesh, upf%nbeta,upf%nbeta) )
    upf%paw%ptfunc(:,:,:) = 0._dp
    DO nb=1,upf%nbeta
       DO mb=1,upf%nbeta
          upf%paw%ptfunc (1:upf%mesh, nb, mb) = &
               upf%pswfc(1:upf%mesh, nb) * upf%pswfc(1:upf%mesh, mb)
          upf%paw%ptfunc(upf%paw%iraug+1:,nb,mb) = 0._dp
          !
          upf%paw%ptfunc (1:upf%mesh, mb, nb) = upf%paw%ptfunc (1:upf%mesh, nb, mb)
       ENDDO
    ENDDO
    !
  END SUBROUTINE read_pp_paw
  !--------------------------------------------------------
  SUBROUTINE read_pp_gipaw ( upf )
    !--------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(pseudo_upf),INTENT(INOUT) :: upf ! the pseudo data
    !
    INTEGER :: nb, mb
    CHARACTER(LEN=24) :: tag
    !
    IF (.NOT. upf%has_gipaw) RETURN
    !
    CALL xmlr_opentag( capitalize_if_v2('pp_gipaw') )
    CALL get_attr ('gipaw_data_format', upf%gipaw_data_format ) 
    IF ( v2 ) THEN
       CALL xmlr_opentag( 'PP_GIPAW_CORE_ORBITALS')
       CALL get_attr ('number_of_core_orbitals', upf%gipaw_ncore_orbitals)
    ELSE
       CALL xmlr_readtag ('number_of_core_orbitals', upf%gipaw_ncore_orbitals) 
       IF ( .NOT. upf%paw_as_gipaw) & 
          CALL xmlr_readtag( 'number_of_valence_orbitals', upf%gipaw_wfs_nchannels)  
    END IF
    ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
    ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
    ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
    ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
    DO nb = 1,upf%gipaw_ncore_orbitals
       IF ( v2 ) THEN
          tag = "PP_GIPAW_CORE_ORBITAL."//i2c(nb)
       ELSE
          tag = 'pp_gipaw_core_orbital'
       END IF
       CALL xmlr_readtag( tag, upf%gipaw_core_orbital(1:upf%mesh,nb) )
       CALL get_attr ('index', mb)
       IF ( nb /= mb ) CALL upf_error('read_pp_gipaw','mismatch',1)
       CALL get_attr ('label', upf%gipaw_core_orbital_el(nb) )
       CALL get_attr ('n', upf%gipaw_core_orbital_n(nb) )
       CALL get_attr ('l', upf%gipaw_core_orbital_l(nb) )
    END DO
    IF ( v2 ) CALL xmlr_closetag ( )
    !
    IF ( upf%paw_as_gipaw) THEN
       !
       !    PAW as GIPAW case: all-electron and pseudo-orbitals not read here
       !
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
       !
    ELSE
       !
       ! Read valence all-electron and pseudo orbitals
       !
       IF ( v2 ) THEN
          CALL xmlr_opentag( 'PP_GIPAW_ORBITALS' )
          CALL get_attr( 'number_of_valence_orbitals', &
               upf%gipaw_wfs_nchannels )
       END IF
       ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
       ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
       ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
       DO nb = 1,upf%gipaw_wfs_nchannels
          IF ( v2 ) THEN
             tag = "PP_GIPAW_ORBITAL."//i2c(nb)
          ELSE
             tag = 'pp_gipaw_orbital'
          END IF
          CALL xmlr_opentag( tag )
          CALL get_attr ('index', mb)
          IF ( nb /= mb ) CALL upf_error('read_pp_gipaw','mismatch',2)
          CALL get_attr ('label', upf%gipaw_wfs_el(nb) )
          CALL get_attr ('l',     upf%gipaw_wfs_ll(nb) )
          CALL get_attr ('cutoff_radius', upf%gipaw_wfs_rcut(nb) )
          CALL get_attr ('ultrasoft_cutoff_radius', upf%gipaw_wfs_rcutus(nb) )
          CALL xmlr_readtag( capitalize_if_v2('pp_gipaw_wfs_ae'), &
               upf%gipaw_wfs_ae(1:upf%mesh,nb) )
          CALL xmlr_readtag( capitalize_if_v2('pp_gipaw_wfs_ps'),&
               upf%gipaw_wfs_ps(1:upf%mesh,nb) )
          CALL xmlr_closetag ()
       END DO
       IF ( v2 ) CALL xmlr_closetag( )
       !
       ! Read all-electron and pseudo local potentials
       !
       CALL xmlr_opentag( capitalize_if_v2('pp_gipaw_vlocal') )
       CALL xmlr_readtag( capitalize_if_v2('pp_gipaw_vlocal_ae'), &
            upf%gipaw_vlocal_ae(1:upf%mesh) )
       CALL xmlr_readtag( capitalize_if_v2('pp_gipaw_vlocal_ps'), &
            upf%gipaw_vlocal_ps(1:upf%mesh) )
       CALL xmlr_closetag ()
    END IF
    CALL xmlr_closetag () ! end pp_gipaw
    !
  END SUBROUTINE read_pp_gipaw
  !
END MODULE read_upf_new_module
