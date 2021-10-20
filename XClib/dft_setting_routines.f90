!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE dft_setting_routines
  !--------------------------------------------------------------------------
  !! Routines to set and/or recover DFT names, parameters and flags.
  !
  USE xclib_utils_and_para,  ONLY: stdout
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
#endif  
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: xclib_set_dft_from_name, xclib_set_dft_IDs,           &
            xclib_set_auxiliary_flags, xclib_set_threshold,       &
            xclib_set_exx_fraction, xclib_set_finite_size_volume, &
            set_screening_parameter, set_gau_parameter
  PUBLIC :: xclib_get_name, xclib_get_ID,             & 
            xclib_get_dft_short, xclib_get_dft_long,  &
            xclib_get_exx_fraction, xclib_get_finite_size_cell_volume, &
            get_screening_parameter, get_gau_parameter
  PUBLIC :: xclib_dft_is, xclib_dft_is_libxc, xclib_init_libxc,  &
            start_exx, stop_exx, dft_has_finite_size_correction, &
            exx_is_active, igcc_is_lyp, xclib_reset_dft, dft_force_hybrid, &
            xclib_finalize_libxc
  PUBLIC :: set_libxc_ext_param, get_libxc_ext_param
  PUBLIC :: capital
  !
CONTAINS
  !
  !-------------------------------------------------------------------------
  SUBROUTINE xclib_set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !! Translates a string containing the exchange-correlation name
    !! into internal indices iexch, icorr, igcx, igcc, inlc, imeta.
    !
    USE xclib_utils_and_para,ONLY: nowarning
    USE dft_setting_params,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                                   discard_input_dft, is_libxc, dft, scan_exx, notset
    USE qe_dft_list,         ONLY: nxc, ncc, ngcx, ngcc, nmeta, get_IDs_from_shortname, &
                                   dft_LDAx_name, dft_LDAc_name, dft_GGAx_name,         &
                                   dft_GGAc_name, dft_MGGA_name
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !! DFT full name
    !
    ! ... local variables
    !
    INTEGER :: leng, l, i
    CHARACTER(len=150):: dftout
    LOGICAL :: dft_defined
    LOGICAL :: check_libxc
    CHARACTER(len=1) :: lxc
    INTEGER :: ID_vec(6)
    INTEGER :: save_iexch, save_icorr, save_igcx, save_igcc, save_meta, &
               save_metac
    !
    ! Exit if set to discard further input dft
    !
    IF ( discard_input_dft ) RETURN
    !
    is_libxc(:) = .FALSE.
    !
    ! save current status of XC indices
    !
    dft_defined = .FALSE.
    !
    save_iexch = iexch
    save_icorr = icorr
    save_igcx  = igcx
    save_igcc  = igcc
    save_meta  = imeta
    save_metac = imetac
    !
    ! convert to uppercase
    !
    leng = LEN_TRIM(dft_)
    dftout = ' '
    !
    DO l = 1, leng
       dftout(l:l) = capital( dft_(l:l) )
    ENDDO
    !
    ! ----------------------------------------------
    ! NOW WE CHECK ALL THE SHORT NAMES
    ! Note: comparison is done via exact matching
    ! ----------------------------------------------
    !
    CALL get_IDs_from_shortname( dftout, ID_vec(1:6) )
    !
    IF ( ALL(ID_vec(1:6)/=notset) ) THEN
       !
       iexch = ID_vec(1) ;  icorr = ID_vec(2)
       igcx  = ID_vec(3) ;  igcc  = ID_vec(4)
       imeta = ID_vec(5) ;  imetac= ID_vec(6)
       dft_defined = .TRUE.
       !
    ELSE
       !
       ! ----------------------------------------------------------------------
       ! CHECK LIBXC FUNCTIONALS BY INDEX NOTATION, IF PRESENT
       ! ----------------------------------------------------------------------
       !
       IF (dftout(1:3) .EQ. 'XC-') THEN
          !
          ! ... short notation with libxc DFTs: 'XC-000i-000i-000i-000i-000i-000i'
          !
          READ( dftout(4:6), * ) iexch
          READ( dftout(7:7), '(a)' ) lxc
          IF (lxc == 'L') is_libxc(1) = .TRUE.
          READ( dftout(9:11), * ) icorr
          READ( dftout(12:12), '(a)' ) lxc
          IF (lxc == 'L') is_libxc(2) = .TRUE.
          READ( dftout(14:16), * ) igcx
          READ( dftout(17:17), '(a)' ) lxc
          IF (lxc == 'L') is_libxc(3) = .TRUE.
          READ( dftout(19:21), * ) igcc
          READ( dftout(22:22), '(a)' ) lxc
          IF (lxc == 'L') is_libxc(4) = .TRUE.
          READ( dftout(24:26), * ) imeta
          READ( dftout(27:27), '(a)' ) lxc
          IF (lxc == 'L') is_libxc(5) = .TRUE.
          READ( dftout(29:31), * ) imetac
          READ( dftout(32:32), '(a)' ) lxc
          IF (lxc == 'L') is_libxc(6) = .TRUE.
          !
          dft_defined = .TRUE.
          !
#if !defined(__LIBXC)
          IF (ANY(is_libxc(:))) THEN
            CALL xclib_error( 'set_dft_from_name', 'libxc functionals needed, but &
                                            &libxc is not active', 1 )
          ENDIF
#endif
          !
       ENDIF
       !
    ENDIF
    !
    !----------------------------------------------------------------
    ! If the DFT was not yet defined, check every part of the string
    !----------------------------------------------------------------
    !
    IF (.NOT. dft_defined) THEN
       !
       iexch = matching( dftout, nxc,   dft_LDAx_name )
       icorr = matching( dftout, ncc,   dft_LDAc_name )
       igcx  = matching( dftout, ngcx,  dft_GGAx_name )
       igcc  = matching( dftout, ngcc,  dft_GGAc_name )
       imeta = matching( dftout, nmeta, dft_MGGA_name )
       imetac = 0
       !
    ENDIF
    !
#if defined(__LIBXC)
    IF (.NOT. dft_defined) CALL matching_libxc( dftout )
#endif
    !
    ! Back compatibility - TO BE REMOVED
    !
    IF (igcx == 14) igcx = 3 ! PBE -> PBX
    IF (igcc ==  9) igcc = 4 ! PBE -> PBC
    !
    IF (igcx == 6 .AND. .NOT.nowarning ) CALL xclib_infomsg( 'set_dft_from_name', 'OPTX &
                                                             &untested! please test' )
    !
    ! ... Check for conflicts with MGGA functionals of QE
    !
    IF (imeta/=0 .AND. (.NOT.is_libxc(5)) .AND. (iexch+icorr+igcx+igcc)>0 ) THEN
      WRITE(stdout,'(/5X,"WARNING: the MGGA functional with ID ",I4," has been ",&
                    &/5X,"read together with other dft terms, but it should ",   &
                    &/5X,"stand alone in order to work properly. The other ",    &
                    &/5x,"terms will be ignored.")' ) imeta
      iexch=0 ; igcx=0
      icorr=0 ; igcc=0
    ENDIF
    !
    ! ... A workaround to keep the q-e shortname notation for SCAN and TB09
    !     functionals valid.
#if defined(__LIBXC)
    IF (imeta==5 .OR. imeta==6) THEN
      imeta = 263
      IF (imeta==6) THEN  ! SCAN0
        scan_exx = .TRUE.
        imeta = 264
      ENDIF
      imetac = 267
      is_libxc(5:6) = .TRUE.
    ELSEIF (imeta==3) THEN
      imeta  = 208
      imetac = 231
      is_libxc(5:6) = .TRUE.
    ENDIF
#else
    IF (imeta==3 .OR. imeta==5 .OR. imeta==6) &
      CALL xclib_error( 'set_dft_from_name', 'libxc needed for this functional', 2 )
#endif
    !
    ! Fill variables and exit
    !
    dft = dftout
    !
    dft_defined = .TRUE.
    !
#if defined(__LIBXC)
    dft_defined = xclib_set_dft_IDs( iexch, icorr, igcx, igcc, imeta, imetac )
#else
    dft_defined = xclib_set_dft_IDs( iexch, icorr, igcx, igcc, imeta, 0 )
#endif
    !
    !dft_longname = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
    !     &//gradc (igcc) //'-'// nonlocc(inlc)
    !
    ! check dft has not been previously set differently
    !
    IF (save_iexch /= notset .AND. save_iexch /= iexch) THEN
       WRITE(stdout,*) iexch, save_iexch
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for iexch', 2 )
    ENDIF
    IF (save_icorr /= notset .AND. save_icorr /= icorr) THEN
       WRITE(stdout,*) icorr, save_icorr
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for icorr', 3 )
    ENDIF
    IF (save_igcx /= notset  .AND. save_igcx /= igcx)   THEN
       WRITE(stdout,*) igcx, save_igcx
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for igcx',  4 )
    ENDIF
    IF (save_igcc /= notset  .AND. save_igcc /= igcc)   THEN
       WRITE(stdout,*) igcc, save_igcc
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for igcc',  5 )
    ENDIF
    IF (save_meta /= notset  .AND. save_meta /= imeta)  THEN
       WRITE(stdout,*) imeta, save_meta
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for imeta', 6 )
    ENDIF
    IF (save_metac /= notset  .AND. save_metac /= imetac)  THEN
       WRITE(stdout,*) imetac, save_metac
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for imetac', 7 )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE xclib_set_dft_from_name
  !
  !----------------------------------------------------------------------------------
  LOGICAL FUNCTION xclib_set_dft_IDs( iexch_, icorr_, igcx_, igcc_, imeta_, imetac_ )
    !--------------------------------------------------------------------------------
    !! Set XC functional IDs. It can be easily extended to include libxc functionals
    !! by adding the 'is_libxc_' array as argument. 
    !
    USE dft_setting_params,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iexch_, icorr_
    INTEGER, INTENT(IN) :: igcx_, igcc_
    INTEGER, INTENT(IN) :: imeta_, imetac_
    !LOGICAL, OPTIONAL   :: is_libxc_(6)
    !
    iexch = iexch_
    icorr = icorr_
    igcx  = igcx_
    igcc  = igcc_
    imeta = imeta_
    imetac = imetac_
    !
    !IF ( PRESENT(is_libxc_) ) is_libxc = is_libxc_
    !
    xclib_set_dft_IDs = .TRUE.
    !
    RETURN
    !
  END FUNCTION xclib_set_dft_IDs
  !
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION matching( dft_, n, name )
    !-----------------------------------------------------------------
    !! Looks for matches between the names of each single term of the 
    !! xc-functional and the input dft string.
    !
    USE dft_setting_params,  ONLY: notset
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN):: n
    CHARACTER(LEN=*), INTENT(IN):: name(0:n)
    CHARACTER(LEN=*), INTENT(IN):: dft_
    !
    INTEGER :: i
    !
    matching = notset
    !
    DO i = n, 0, -1
       IF ( matches(name(i), TRIM(dft_)) ) THEN
          !
          IF (matching==notset .OR. name(i)=='REVX') THEN
             !WRITE(*, '("matches",i2,2X,A,2X,A)') i, name(i), TRIM(dft)
             matching = i
          ELSE
#if defined(__LIBXC)
             IF (name(i)=='B88' .OR. name(i)=='CX0') CYCLE
#else
             IF (name(i)/='B88' .AND. name(i)/='CX0') THEN
                WRITE(stdout, '(2(2X,i2,2X,A))') i, TRIM(name(i)), &
                                     matching, TRIM(name(matching))
                CALL xclib_error( 'set_dft', 'two conflicting matching values', 1 )
             ENDIF
#endif
          ENDIF
       ENDIF
    ENDDO
    !
    IF (matching == notset) matching = 0
    !
  END FUNCTION matching
  !
  !
#if defined(__LIBXC)
  !--------------------------------------------------------------------------------
  SUBROUTINE matching_libxc( dft_ )
    !------------------------------------------------------------------------------
    !! It spans the libxc functionals and looks for matches with the input dft
    !! string. Then stores the corresponding indices.  
    !! It also makes some compatibility checks.
    !
    USE dft_setting_params,   ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                                    is_libxc, exx_fraction, xc_kind_error
    USE xclib_utils_and_para, ONLY: nowarning
    ! 
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !
    CHARACTER(LEN=256) :: name
    INTEGER :: i, l, prev_len(6), fkind, fkind_v(3), family
    INTEGER, PARAMETER :: ID_MAX_LIBXC=999
    TYPE(xc_f03_func_t) :: xc_func
    TYPE(xc_f03_func_info_t) :: xc_info
#if (XC_MAJOR_VERSION>5)
    !workaround to keep compatibility with libxc develop version
    INTEGER, PARAMETER :: XC_FAMILY_HYB_GGA  = -10 
    INTEGER, PARAMETER :: XC_FAMILY_HYB_MGGA = -11 
#endif
    !
    prev_len(:) = 1
    !
    DO i = 1, ID_MAX_LIBXC
       !
       name = xc_f03_functional_get_name( i )
       !
       DO l = 1, LEN_TRIM(name)
          name(l:l) = capital( name(l:l) )
       ENDDO
       !
       IF ( TRIM(name) == '' ) CYCLE
       !
       IF ( matches(TRIM(name), TRIM(dft_)) ) THEN
          !
          !WRITE(*, '("matches libxc",i2,2X,A,2X,A)') i, TRIM(name), TRIM(dft)
          !
          fkind=-100 ; family=-100
          CALL xc_f03_func_init( xc_func, i, 1 )
          xc_info = xc_f03_func_get_info( xc_func )
          fkind = xc_f03_func_info_get_kind( xc_info )
          family = xc_f03_func_info_get_family( xc_info )
          CALL xc_f03_func_end( xc_func )
          !   
          SELECT CASE( family )
          CASE( XC_FAMILY_LDA )
             IF (fkind==XC_EXCHANGE) THEN
                IF ( LEN(TRIM(name)) > prev_len(1) ) iexch = i
                is_libxc(1) = .TRUE.
                prev_len(1) = LEN(TRIM(name))
             ELSEIF (fkind==XC_CORRELATION .OR. fkind==XC_EXCHANGE_CORRELATION) THEN
                IF ( LEN(TRIM(name)) > prev_len(2) ) icorr = i
                is_libxc(2) = .TRUE.
                prev_len(2) = LEN(TRIM(name))
                IF (fkind==XC_EXCHANGE_CORRELATION) THEN
                  iexch = 0
                  is_libxc(1) = .FALSE.
                ENDIF
             ELSE
                xc_kind_error = .TRUE.
             ENDIF
             fkind_v(1) = fkind
             !
          CASE( XC_FAMILY_GGA, XC_FAMILY_HYB_GGA )
             IF (fkind==XC_EXCHANGE) THEN
                IF ( LEN(TRIM(name)) > prev_len(3) ) igcx = i
                is_libxc(3) = .TRUE.
                prev_len(3) = LEN(TRIM(name))
             ELSEIF (fkind==XC_CORRELATION .OR. fkind==XC_EXCHANGE_CORRELATION) THEN
                IF ( LEN(TRIM(name)) > prev_len(4) ) igcc = i
                is_libxc(4) = .TRUE.
                prev_len(4) = LEN(TRIM(name))
             ELSE
                xc_kind_error = .TRUE.
             ENDIF
             fkind_v(2) = fkind
             !
          CASE( XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA )
             IF (fkind==XC_EXCHANGE) THEN
                IF ( LEN(TRIM(name)) > prev_len(5) ) imeta = i
                is_libxc(5) = .TRUE.
                prev_len(5) = LEN(TRIM(name))
             ELSEIF (fkind==XC_CORRELATION .OR. fkind==XC_EXCHANGE_CORRELATION) THEN
                IF ( LEN(TRIM(name)) > prev_len(6) ) imetac = i
                is_libxc(6) = .TRUE.
                prev_len(6) = LEN(TRIM(name))
             ELSE
                xc_kind_error = .TRUE.
             ENDIF
             fkind_v(3) = fkind
             !
          END SELECT
          !
       ENDIF
       !
    ENDDO
    !
    ! ... overlaps check (between qe and libxc names)
    !
    IF (ANY(.NOT.is_libxc(:)).AND.ANY(is_libxc(:))) CALL check_overlaps_qe_libxc(dft_)
    !
    ! ... Compatibility checks
    !
    IF ( xc_kind_error .AND. .NOT.nowarning ) &
       CALL xclib_error( 'matching_libxc', 'a Libxc functional of a kind not &
                         &usable in QE has been found', 1 )
    !
    ! LDA:
    IF (iexch/=0 .AND. is_libxc(2).AND. fkind_v(1)==XC_EXCHANGE_CORRELATION)  &
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an EXCHANGE+CORRELATION &
                           &functional has been found together with an exchange&
                           & one (LDA)' )
    ! GGA:
    IF (igcx/=0 .AND. is_libxc(4).AND. fkind_v(2)==XC_EXCHANGE_CORRELATION)   &
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an EXCHANGE+CORRELATION &
                           &functional has been found together with an exchange&
                           & one (GGA)' )
    !
    IF ( (is_libxc(3).AND.iexch/=0) .OR. (is_libxc(4).AND. icorr/=0) ) &
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an LDA functional has bee&
                           &n found, but libxc GGA functionals already include t&
                           &he LDA part' )
    ! mGGA:
    ! (imeta defines both exchange and correlation term for q-e mGGA functionals)
    IF (imeta/=0 .AND. (.NOT. is_libxc(5)) .AND. imetac/=0)   &
       CALL xclib_error( 'matching_libxc', 'Two conflicting metaGGA functionals &
                         &have been found', 2 )
    !
    IF (imeta/=0 .AND. is_libxc(6).AND. fkind_v(3)==XC_EXCHANGE_CORRELATION)  &   
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an EXCHANGE+CORRELATION f&
                           &unctional has been found together with an exchange o&
                           &ne (mGGA)' )
    !   
  END SUBROUTINE matching_libxc
  !
  !--------------------------------------------------------------------------
  SUBROUTINE check_overlaps_qe_libxc( dft_ )
    !------------------------------------------------------------------------
    !! It fixes eventual overlap issues between qe and libxc names when qe and
    !! libxc functionals are used together.
    !
    USE dft_setting_params, ONLY: iexch, icorr, igcx, igcc, imeta, imetac, is_libxc
    USE qe_dft_list,        ONLY: nxc, ncc, ngcx, ngcc, nmeta, dft_LDAx_name,  &
                                  dft_LDAc_name, dft_GGAx_name, dft_GGAc_name, &
                                  dft_MGGA_name
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !! DFT full name
    !
    CHARACTER(LEN=4) :: qe_name
    CHARACTER(LEN=256) :: lxc_name
    INTEGER :: i, ii, l, ch, qedft, nlxc
    INTEGER :: ID_vec(6)
    !
    ID_vec(1)=iexch ; ID_vec(2)=icorr
    ID_vec(3)=igcx  ; ID_vec(4)=igcc
    ID_vec(5)=imeta ; ID_vec(6)=imetac
    !
    DO ch = 1, 5
       IF (.NOT.is_libxc(ch)) THEN
          !
          SELECT CASE( ch )
          CASE( 1 )
             qe_name = dft_LDAx_name(iexch)
          CASE( 2 )
             qe_name = dft_LDAc_name(icorr)
          CASE( 3 )
             qe_name = dft_GGAx_name(igcx)
          CASE( 4 )
             qe_name = dft_GGAc_name(igcc)
          CASE( 5 )
             qe_name = dft_MGGA_name(imeta)
          END SELECT
          !
          qedft = 0
          i = 0
          DO WHILE ( i < LEN_TRIM(dft_) )
            i = i + 1
            IF ( matches( TRIM(qe_name), TRIM(dft_(i:i+1)) ) ) THEN
               qedft = qedft + 1
               i = i + 1
            ELSEIF (matches( TRIM(qe_name), TRIM(dft_(i:i+2)) ) ) THEN
               qedft = qedft + 1
               i = i + 2
            ELSEIF (matches( TRIM(qe_name), TRIM(dft_(i:i+3)) ) ) THEN
               qedft = qedft + 1
               i = i + 3
            ENDIF
          ENDDO
          !
          nlxc = 0
          DO i = 1, 6
            IF (is_libxc(i)) THEN
              lxc_name = xc_f03_functional_get_name( ID_vec(i) )
              DO l = 1, LEN_TRIM(lxc_name)
                 lxc_name(l:l) = capital( lxc_name(l:l) )
              ENDDO
              nlxc = 0
              ii = 0
              DO WHILE ( ii < LEN_TRIM(lxc_name) )
                ii = ii + 1
                IF ( matches( TRIM(qe_name), TRIM(lxc_name(ii:ii+1)) ) ) THEN
                  nlxc = nlxc + 1
                  ii = ii + 1
                ELSEIF (matches( TRIM(qe_name), TRIM(lxc_name(ii:ii+2)) ) ) THEN
                  nlxc = nlxc + 1
                  ii = ii + 2
                ELSEIF (matches( TRIM(qe_name), TRIM(lxc_name(ii:ii+3)) ) ) THEN
                  nlxc = nlxc + 1
                  ii = ii + 3
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          !
          IF (qedft == nlxc) ID_vec(ch) = 0  
          !
       ENDIF
    ENDDO
    !
    iexch = ID_vec(1) ;  icorr  = ID_vec(2)
    igcx  = ID_vec(3) ;  igcc   = ID_vec(4)
    imeta = ID_vec(5) ;  imetac = ID_vec(6)
    !
  END SUBROUTINE
#endif
  !
  !
  !---------------------------------------------------------------------------
  SUBROUTINE xclib_set_auxiliary_flags( isnonlocc )
    !-------------------------------------------------------------------------
    !! Set logical flags describing the complexity of the xc functional
    !! define the fraction of exact exchange used by hybrid fuctionals.
    !
    USE kind_l,              ONLY: DP
    USE dft_setting_params,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                                   islda, isgradient, ismeta, exx_fraction, &
                                   screening_parameter, gau_parameter,      & 
                                   ishybrid, has_finite_size_correction, is_libxc
    !  
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: isnonlocc
    !! The non-local part, for now, is not included in xc_lib, but this variable
    !! is needed to establish 'isgradient'.
    LOGICAL :: is_libxc13, is_libxc12
    !
    ismeta    = (imeta+imetac > 0)
    isgradient= (igcx > 0) .OR.  (igcc > 0)  .OR. ismeta .OR. isnonlocc
    islda     = (iexch> 0) .AND. (icorr > 0) .AND. .NOT. isgradient
    is_libxc13 = is_libxc(1) .OR. is_libxc(3)
    !
    ! PBE0/DF0
    IF ( iexch==6 .AND. .NOT.is_libxc(1) ) exx_fraction = 0.25_DP
    IF ( igcx==8  .AND. .NOT.is_libxc(3) ) exx_fraction = 0.25_DP
    ! CX0P
    IF ( iexch==6 .AND. igcx==31 .AND. .NOT.is_libxc13 ) exx_fraction = 0.20_DP
    ! B86BPBEX
    IF ( iexch==6 .AND. igcx==41 .AND. .NOT.is_libxc13 ) exx_fraction = 0.25_DP
    ! BHANDHLYP
    IF ( iexch==6 .AND. igcx==42 .AND. .NOT.is_libxc13 ) exx_fraction = 0.50_DP
    ! HSE
    IF ( igcx ==12 .AND. .NOT.is_libxc(3) ) THEN
       exx_fraction = 0.25_DP
       screening_parameter = 0.106_DP
    ENDIF
    ! gau-pbe
    IF ( igcx ==20 .AND. .NOT.is_libxc(3) ) THEN
       exx_fraction = 0.24_DP
       gau_parameter = 0.150_DP
    ENDIF
    ! HF or OEP
    IF ( iexch==4 .AND. .NOT.is_libxc(1)) exx_fraction = 1.0_DP
    IF ( iexch==5 .AND. .NOT.is_libxc(1)) exx_fraction = 1.0_DP
    ! B3LYP or B3LYP-VWN-1-RPA
    IF ( iexch == 7 .AND. .NOT.is_libxc(3) ) exx_fraction = 0.2_DP
    ! X3LYP
    IF ( iexch == 9 .AND. .NOT.is_libxc(3) ) exx_fraction = 0.218_DP
    !
    ishybrid = ( exx_fraction /= 0.0_DP )
    !
    has_finite_size_correction = ( (iexch==8 .AND. .NOT.is_libxc(1)) .OR. &
                                   (icorr==10.AND. .NOT.is_libxc(2)) )
    !
    RETURN
    !
  END SUBROUTINE xclib_set_auxiliary_flags
  !
  !======================= EXX-hybrid ======================================
  ! 
!   !-----------------------------------------------------------------------
!   SUBROUTINE enforce_dft_exxrpa( )
!     !---------------------------------------------------------------------
!     !
!     USE dft_setting_params
!     !
!     IMPLICIT NONE
!     !
!     !character(len=*), intent(in) :: dft_
!     !logical, intent(in), optional :: nomsg
!     !
!     iexch = 0; icorr = 0; igcx = 0; igcc = 0
!     exx_fraction = 1.0_DP
!     ishybrid = ( exx_fraction /= 0.0_DP )
!     !
!     WRITE(*,'(/,5x,a)') "XC functional enforced to be EXXRPA"
!     CALL write_dft_name
!     WRITE(*,'(5x,a)') "!!! Any further DFT definition will be discarded"
!     WRITE(*,'(5x,a/)') "!!! Please, verify this is what you really want !"
!     !
!     RETURN
!     !
!   END SUBROUTINE enforce_dft_exxrpa
!   !
!   !
!   !-----------------------------------------------------------------------
!   SUBROUTINE init_dft_exxrpa( )
!     !-----------------------------------------------------------------------
!     !
!     USE dft_setting_params
!     !
!     IMPLICIT NONE
!     !
!     exx_fraction = 1.0_DP
!     ishybrid = ( exx_fraction /= 0.0_DP )
!     !
!     WRITE(*,'(/,5x,a)') "Only exx_fraction is set to 1.d0"
!     WRITE(*,'(5x,a)') "XC functional still not changed"
!     !
!     CALL write_dft_name
!     !
!     RETURN
!     !
!   END SUBROUTINE init_dft_exxrpa
  !
  SUBROUTINE start_exx
    !! Activate exact exchange (exx_started=TRUE)
    USE dft_setting_params, ONLY: ishybrid, exx_started
    IMPLICIT NONE
    IF (.NOT. ishybrid) &
       CALL xclib_error( 'start_exx', 'dft is not hybrid, wrong call', 1 )
    exx_started = .TRUE.
  END SUBROUTINE start_exx
  !-----------------------------------------------------------------------
  SUBROUTINE stop_exx
    !! Deactivate exact exchange (exx_started=FALSE)
    USE dft_setting_params, ONLY: ishybrid, exx_started
    IMPLICIT NONE
    IF (.NOT. ishybrid) &
       CALL xclib_error( 'stop_exx', 'dft is not hybrid, wrong call', 1 )
    exx_started = .FALSE.
  END SUBROUTINE stop_exx
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_set_exx_fraction( exx_fraction_ )
    !! Impose input parameter as exact exchange fraction value
    USE kind_l,             ONLY: DP
    USE dft_setting_params, ONLY: exx_fraction
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: exx_fraction_
    !! Imposed value of exact exchange fraction
    exx_fraction = exx_fraction_
    WRITE( stdout,'(5x,a,f6.2)') 'EXX fraction changed: ', exx_fraction
    RETURN
  END SUBROUTINE xclib_set_exx_fraction
  !-----------------------------------------------------------------------
  SUBROUTINE dft_force_hybrid( request )
    !! Impose hybrid condition.
    USE dft_setting_params, ONLY: ishybrid
    IMPLICIT NONE
    LOGICAL, OPTIONAL, INTENT(INOUT) :: request
    !! Impose input request as hybrid condition and return output request
    !! as previous hybrid condition.
    LOGICAL :: aux
    IF (PRESENT(request)) THEN
      aux = ishybrid
      ishybrid = request
      request = aux
    ELSE
      ishybrid= .TRUE.
    ENDIF
  END SUBROUTINE dft_force_hybrid
  !-----------------------------------------------------------------------
  FUNCTION exx_is_active()
     !! TRUE if exact exchange is active.
     USE dft_setting_params, ONLY: exx_started
     IMPLICIT NONE
     LOGICAL :: exx_is_active
     exx_is_active = exx_started
  END FUNCTION exx_is_active
  !-----------------------------------------------------------------------
  FUNCTION xclib_get_exx_fraction()
     !! Recover exact exchange fraction.
     USE kind_l,             ONLY: DP
     USE dft_setting_params, ONLY: exx_fraction
     IMPLICIT NONE
     REAL(DP) :: xclib_get_exx_fraction
     xclib_get_exx_fraction = exx_fraction
     RETURN
  END FUNCTION xclib_get_exx_fraction
  !-----------------------------------------------------------------------
  !
  !
  !============ PBE gau-screening ========================================
  !
  !-----------------------------------------------------------------------
  SUBROUTINE set_screening_parameter( scrparm_ )
    !! Impose input parameter as screening parameter (for pbexsr)
    USE kind_l,             ONLY: DP
    USE dft_setting_params, ONLY: screening_parameter
    IMPLICIT NONE
    REAL(DP):: scrparm_
    !! Value to impose as screening parameter
    screening_parameter = scrparm_
    WRITE(stdout,'(5x,a,f12.7)') 'EXX Screening parameter changed: ', &
         & screening_parameter
  END SUBROUTINE set_screening_parameter
  !-----------------------------------------------------------------------
  FUNCTION get_screening_parameter()
     !! Recover screening parameter (for pbexsr)
     USE kind_l,             ONLY: DP
     USE dft_setting_params, ONLY: screening_parameter
     IMPLICIT NONE
     REAL(DP):: get_screening_parameter
     get_screening_parameter = screening_parameter
     RETURN
  END FUNCTION get_screening_parameter
  !-----------------------------------------------------------------------
  SUBROUTINE set_gau_parameter( gauparm_ )
    !! Impose input parameter as gau parameter (for gau-pbe)
    USE kind_l,             ONLY: DP
    USE dft_setting_params, ONLY: gau_parameter
    IMPLICIT NONE
    REAL(DP):: gauparm_
    !! Value to impose as gau parameter
    gau_parameter = gauparm_
    WRITE(stdout,'(5x,a,f12.7)') 'EXX Gau parameter changed: ', &
         & gau_parameter
  END SUBROUTINE set_gau_parameter
  !-----------------------------------------------------------------------
  FUNCTION get_gau_parameter()
    !! Recover gau parameter (for gau-pbe)
    USE kind_l,             ONLY: DP
    USE dft_setting_params, ONLY: gau_parameter
    IMPLICIT NONE
    REAL(DP):: get_gau_parameter
    get_gau_parameter = gau_parameter
    RETURN
  END FUNCTION get_gau_parameter
  !-----------------------------------------------------------------------
  !
  !
  !============ DFT NAME & ID SETTING AND RECOVERY =======================
  !
  !-----------------------------------------------------------------------
  FUNCTION xclib_get_ID( family, kindf )
     !--------------------------------------------------------------------
     !! Get functionals index of \(\text{family}\) and \(\text{kind}\).
     !
     USE dft_setting_params, ONLY: iexch, icorr, igcx, igcc, imeta, imetac
     !
     IMPLICIT NONE
     !
     INTEGER :: xclib_get_ID
     CHARACTER(len=*), INTENT(IN) :: family
     !! LDA, GGA or MGGA
     CHARACTER(len=*), INTENT(IN) :: kindf
     !! EXCH or CORR
     !
     CHARACTER(len=4) :: cfamily, ckindf
     INTEGER :: i, ln
     !
     ln = LEN_TRIM(family)
     !
     DO i = 1, ln
       cfamily(i:i) = capital(family(i:i))
     ENDDO
     DO i = 1, 4
       ckindf(i:i) = capital(kindf(i:i))
     ENDDO
     !
     SELECT CASE( cfamily(1:ln) )
     CASE( 'LDA' )
       IF (ckindf=='EXCH') xclib_get_ID = iexch
       IF (ckindf=='CORR') xclib_get_ID = icorr
     CASE( 'GGA' )
       IF (ckindf=='EXCH') xclib_get_ID = igcx
       IF (ckindf=='CORR') xclib_get_ID = igcc
     CASE( 'MGGA' )
       IF (ckindf=='EXCH') xclib_get_ID = imeta
       IF (ckindf=='CORR') xclib_get_ID = imetac
     CASE DEFAULT
       CALL xclib_error( 'xclib_get_id', 'input not recognized', 1 )
     END SELECT
     !
     RETURN
     !
  END FUNCTION xclib_get_ID
  !
  !-------------------------------------------------------------------       
  SUBROUTINE xclib_get_name( family, kindf, name )
     !----------------------------------------------------------------
     !! Gets QE name for 'family'-'kind' term of the XC functional.
     !
     USE dft_setting_params, ONLY: iexch, icorr, igcx, igcc, imeta, imetac
     USE qe_dft_list,        ONLY: dft_LDAx_name, dft_LDAc_name, dft_GGAx_name, &
                                   dft_GGAc_name, dft_MGGA_name
     !
     IMPLICIT NONE
     !
     CHARACTER(len=4) :: name
     CHARACTER(len=*), INTENT(IN) :: family
     !! LDA, GGA or MGGA
     CHARACTER(len=*), INTENT(IN) :: kindf
     !! EXCH or CORR
     !
     CHARACTER(len=4) :: cfamily, ckindf
     INTEGER :: i, ln
     !
     ln = LEN_TRIM(family)
     !
     DO i = 1, ln
       cfamily(i:i) = capital(family(i:i))
     ENDDO
     DO i = 1, 4
       ckindf(i:i) = capital(kindf(i:i))
     ENDDO
     !
     SELECT CASE( cfamily(1:ln) )
     CASE( 'LDA' )
       IF (ckindf=='EXCH') name = dft_LDAx_name(iexch)
       IF (ckindf=='CORR') name = dft_LDAc_name(icorr)
     CASE( 'GGA' )
       IF (ckindf=='EXCH') name = dft_GGAx_name(igcx)
       IF (ckindf=='CORR') name = dft_GGAc_name(igcc)
     CASE( 'MGGA' )
       IF (ckindf=='EXCH') name = dft_MGGA_name(imeta) 
     CASE DEFAULT
       CALL xclib_error( 'get_name', 'input not recognized', 1 )
     END SELECT
     !
     RETURN
     !
  END SUBROUTINE xclib_get_name
  !
  !--------------------------------------------------------------------
  FUNCTION xclib_dft_is_libxc( family, kindf )
     !-----------------------------------------------------------------
     !! Establish if the XC term family-kind is Libxc or not.
     !
     USE dft_setting_params,  ONLY: is_libxc
     !
     IMPLICIT NONE
     !
     LOGICAL :: xclib_dft_is_libxc
     CHARACTER(len=*), INTENT(IN) :: family
     !! LDA, GGA or MGGA
     CHARACTER(len=*), INTENT(IN), OPTIONAL :: kindf
     !! EXCH or CORR
     !
     CHARACTER(len=4) :: cfamily='', ckindf
     INTEGER :: i, ln
     !
     xclib_dft_is_libxc = .FALSE.
     !
     ln = LEN_TRIM(family)
     !
     DO i = 1, ln
       cfamily(i:i) = capital(family(i:i))
     ENDDO
     IF ( PRESENT(kindf) ) THEN
       DO i = 1, 4
         ckindf(i:i) = capital(kindf(i:i))
       ENDDO
       !
       SELECT CASE( cfamily(1:ln) )
       CASE( 'LDA' )
         IF (ckindf=='EXCH') xclib_dft_is_libxc = is_libxc(1)
         IF (ckindf=='CORR') xclib_dft_is_libxc = is_libxc(2)
       CASE( 'GGA' )
         IF (ckindf=='EXCH') xclib_dft_is_libxc = is_libxc(3)
         IF (ckindf=='CORR') xclib_dft_is_libxc = is_libxc(4)
       CASE( 'MGGA' )
         IF (ckindf=='EXCH') xclib_dft_is_libxc = is_libxc(5)
         IF (ckindf=='CORR') xclib_dft_is_libxc = is_libxc(6)
       CASE DEFAULT
         CALL xclib_error( 'xclib_dft_is_libxc', 'input not recognized', 1 )
       END SELECT
     ELSE
       IF (family=='ANY'.AND.ANY(is_libxc(:))) xclib_dft_is_libxc=.TRUE.
     ENDIF
     !
     RETURN
     !
  END FUNCTION
  !
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_reset_dft()
    !---------------------------------------------------------------------
    !! Unset DFT indexes and main parameters.
    USE dft_setting_params
    IMPLICIT NONE
    dft = 'not set'
    iexch  = notset ; icorr  = notset
    igcx   = notset ; igcc   = notset
    imeta  = notset ; imetac = notset
    exx_fraction = 0.d0
    is_libxc(:) = .FALSE.
    exx_started = .FALSE.
    exx_fraction = 0.0_DP
    finite_size_cell_volume = -1._DP
    rho_threshold_lda = 1.E-10_DP
    rho_threshold_gga = 1.E-6_DP   ; grho_threshold_gga = 1.E-10_DP
    rho_threshold_mgga = 1.E-12_DP ; grho2_threshold_mgga = 1.E-24_DP
    tau_threshold_mgga = 1.0E-12_DP
    islda = .FALSE. ; isgradient  = .FALSE.
    has_finite_size_correction = .FALSE.
    finite_size_cell_volume_set = .FALSE.
    ismeta = .FALSE.
    ishybrid = .FALSE.
    scan_exx = .FALSE.
    beeftype = -1 ; beefvdw = 0
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  FUNCTION get_dft_name()
     !---------------------------------------------------------------------
     !! Get full DFT name
     USE dft_setting_params, ONLY: dft
     IMPLICIT NONE
     CHARACTER(LEN=32) :: get_dft_name
     get_dft_name = dft
     RETURN
  END FUNCTION get_dft_name
  !
  !-----------------------------------------------------------------------
  FUNCTION xclib_dft_is( what )
     !---------------------------------------------------------------------
     !! Find if DFT has gradient correction, meta or hybrid.
     !
     USE dft_setting_params,  ONLY: isgradient, ismeta, ishybrid
     !
     IMPLICIT NONE
     !
     LOGICAL :: xclib_dft_is
     CHARACTER(len=*) :: what
     !! gradient, meta or hybrid
     !
     CHARACTER(len=15) :: cwhat
     INTEGER :: i, ln
     !
     ln = LEN_TRIM(what)
     !
     DO i = 1, ln
       cwhat(i:i) = capital(what(i:i))
     ENDDO
     !
     SELECT CASE( cwhat(1:ln) )
     CASE( 'GRADIENT' )
       xclib_dft_is = isgradient
     CASE( 'META' )
       xclib_dft_is = ismeta
     CASE( 'HYBRID' )
       xclib_dft_is = ishybrid
     CASE DEFAULT
       CALL xclib_error( 'xclib_dft_is', 'wrong input', 1 )
     END SELECT
     !
     RETURN
     !
  END FUNCTION xclib_dft_is
  !
  !-----------------------------------------------------------------------
  FUNCTION igcc_is_lyp()
     !! Find if correlation GGA is Lee-Yang-Parr.
     USE dft_setting_params,  ONLY: igcc
     IMPLICIT NONE
     LOGICAL :: igcc_is_lyp
     igcc_is_lyp = (igcc==3 .OR. igcc==7 .OR. igcc==13)
     RETURN
  END FUNCTION igcc_is_lyp
  !-----------------------------------------------------------------------
  FUNCTION dft_has_finite_size_correction()
     !! TRUE if finite size correction present
     USE dft_setting_params, ONLY: has_finite_size_correction
     IMPLICIT NONE
     LOGICAL :: dft_has_finite_size_correction
     dft_has_finite_size_correction = has_finite_size_correction
     RETURN
  END FUNCTION dft_has_finite_size_correction
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_set_finite_size_volume( volume )
     !! Set value for finite size cell volume.
     USE dft_setting_params, ONLY: has_finite_size_correction, &
                                   finite_size_cell_volume,    &
                                   finite_size_cell_volume_set
     IMPLICIT NONE
     REAL, INTENT(IN) :: volume
     !! finite size cell volume
     IF (.NOT. has_finite_size_correction) &
         CALL xclib_error( 'set_finite_size_volume', &
                      'dft w/o finite_size_correction, wrong call', 1 )
     IF (volume <= 0.d0) &
         CALL xclib_error( 'set_finite_size_volume', &
                      'volume is not positive, check omega and/or nk1,nk2,nk3', 1 )
     finite_size_cell_volume = volume
     finite_size_cell_volume_set = .TRUE.
  END SUBROUTINE xclib_set_finite_size_volume
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_get_finite_size_cell_volume( is_present, volume )
     !! Recover value for finite size cell volume.
     USE kind_l,             ONLY: DP
     USE dft_setting_params, ONLY: finite_size_cell_volume, finite_size_cell_volume_set
     IMPLICIT NONE
     LOGICAL, INTENT(OUT) :: is_present
     !! TRUE if finite size correction present
     REAL(DP), INTENT(OUT) :: volume
     !! finite size cell volume
     is_present = finite_size_cell_volume_set
     volume = -1.d0
     IF (is_present) volume = finite_size_cell_volume
  END SUBROUTINE xclib_get_finite_size_cell_volume
  !
  !--------------------------------------------------------------------------
  SUBROUTINE xclib_init_libxc( xclib_nspin, domag )
    !------------------------------------------------------------------------
    !! Initialize Libxc functionals, if present.
    USE dft_setting_params,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                                   is_libxc, libxc_initialized
#if defined(__LIBXC)
    USE xclib_utils_and_para,ONLY: nowarning
    USE dft_setting_params,  ONLY: n_ext_params, xc_func, xc_info, par_list, &
                                   libxc_flags, n_ext_params, exx_fraction,  &
                                   ishybrid
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: xclib_nspin
    LOGICAL, INTENT(IN) :: domag
    !! 1: unpolarized case; 2: polarized
    INTEGER :: iid, ip, p0, pn, ips, nspin0, iflag, family
    INTEGER :: id_vec(6), flags_tot    
    !
#if defined(__LIBXC)
    !call xclib_init_libxc_print_version_info
    !
    nspin0 = xclib_nspin
    IF ( xclib_nspin==4 ) THEN
      nspin0 = 1
      IF ( domag ) nspin0 = 2
    ENDIF  
    !
    id_vec(1)=iexch ; id_vec(2)=icorr
    id_vec(3)=igcx  ; id_vec(4)=igcc
    id_vec(5)=imeta ; id_vec(6)=imetac
    !
    DO iid = 1, 6
      IF (libxc_initialized(iid)) THEN
        CALL xc_f03_func_end( xc_func(iid) )
        libxc_initialized(iid) = .FALSE.
      ENDIF
      IF (is_libxc(iid)) THEN
        CALL xc_f03_func_init( xc_func(iid), id_vec(iid), nspin0 )
        xc_info(iid) = xc_f03_func_get_info( xc_func(iid) )
        family = xc_f03_func_info_get_family( xc_info(iid) )
        !
        flags_tot = xc_f03_func_info_get_flags( xc_info(iid) )
        DO iflag = 15, 0, -1
          libxc_flags(iid,iflag) = 0
          IF ( flags_tot-2**iflag < 0 ) CYCLE
          libxc_flags(iid,iflag) = 1
          flags_tot = flags_tot-2**iflag
        ENDDO
        !
        IF ( family==XC_FAMILY_HYB_GGA .OR. family==XC_FAMILY_HYB_MGGA ) THEN
           exx_fraction = xc_f03_hyb_exx_coef( xc_func(iid) )
           ishybrid = ( exx_fraction /= 0.d0 )
        ENDIF
        !
        n_ext_params(iid) = xc_f03_func_info_get_n_ext_params( xc_info(iid) )
#if (XC_MAJOR_VERSION<=5)
        p0 = 0 ;  pn = n_ext_params(iid)-1 ;  ips = 1
#else
        p0 = 1 ;  pn = n_ext_params(iid)   ;  ips = 0
#endif
        DO ip = p0, pn
          par_list(iid,ip+ips) = xc_f03_func_info_get_ext_params_default_value( &
                                                           xc_info(iid), ip )
        ENDDO
        libxc_initialized(iid) = .TRUE.
        !
        IF ( .NOT. nowarning ) THEN
          IF ( n_ext_params(iid) /= 0 ) &
            WRITE(stdout,'(/5X,"WARNING: libxc functional with ID ",I4," depends",&
                          &/5X," on external parameters: check the user_guide of",&
                          &/5X," QE if you need to modify them or to check their",&
                          &/5x," default values.")' ) id_vec(iid)
          IF ( libxc_flags(iid,0) == 0 ) &
            WRITE(stdout,'(/5X,"WARNING: libxc functional with ID ",I4," does not ",&
                          &/5X,"provide Exc.")' ) id_vec(iid)
          IF ( libxc_flags(iid,1) == 0 ) &
            WRITE(stdout,'(/5X,"WARNING: libxc functional with ID ",I4," does not ",&
                          &/5X,"provide Vxc.")' ) id_vec(iid)
          IF ( libxc_flags(iid,2) == 0 ) &
            WRITE(stdout,'(/5X,"WARNING: libxc functional with ID ",I4," does not ", &
                          &/5X,"provide Vxc derivative.")' ) id_vec(iid)
          IF ( libxc_flags(iid,15) == 1 ) &
            WRITE(stdout,'(/5X,"WARNING: libxc functional with ID ",I4," depends on", &
                          &/5X," the laplacian of the density, which is currently set",&
                          &/5X," to zero.")' ) id_vec(iid)
        ENDIF
      ENDIF  
    ENDDO
    !
#endif
    RETURN
  END SUBROUTINE xclib_init_libxc
  !
  !--------------------------------------------------------------------------
#if defined(__LIBXC)
  SUBROUTINE xclib_init_libxc_print_version_info()
    !------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER :: major, minor, micro
    CALL xc_f03_version(major, minor, micro)
    WRITE(stdout, '(3X,"Using LIBXC version       = ",3I4)') major, minor, micro
    RETURN
  END SUBROUTINE xclib_init_libxc_print_version_info
#endif
  !
  !--------------------------------------------------------------------------
  SUBROUTINE xclib_finalize_libxc()
    !------------------------------------------------------------------------
    !! Finalize Libxc functionals, if present.
    USE dft_setting_params,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                                   is_libxc, libxc_initialized, notset
#if defined(__LIBXC)
    USE dft_setting_params,  ONLY: xc_func, xc_kind_error, libxc_flags, &
                                   n_ext_params, par_list
#endif
    IMPLICIT NONE
    INTEGER :: iid
    INTEGER :: id_vec(6)
    !
#if defined(__LIBXC)
    id_vec(1)=iexch ; id_vec(2)=icorr
    id_vec(3)=igcx  ; id_vec(4)=igcc
    id_vec(5)=imeta ; id_vec(6)=imetac
    !
    DO iid = 1, 6
      IF (is_libxc(iid)) THEN
        CALL xc_f03_func_end( xc_func(iid) )
        libxc_initialized(iid) = .FALSE.
        is_libxc(iid) = .FALSE.
        libxc_flags(:,:) = notset
        n_ext_params(:) = 0
        par_list(:,:) = 0.d0
      ENDIF  
    ENDDO
    xc_kind_error = .FALSE.
#endif
    RETURN
  END SUBROUTINE xclib_finalize_libxc
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_libxc_ext_param( sid, i_param, param )
    !------------------------------------------------------------------------
    !! Routine to set external parameters of some Libxc functionals.  
    !! In order to get a list and description of all the available parameters
    !! for a given Libxc functional you can use the \(\texttt{xclib_test} 
    !! program with input \(\text{test}=\text{'dft-info'}\).
    USE kind_l,       ONLY: DP
#if defined(__LIBXC)
    USE dft_setting_params,  ONLY: xc_func, par_list
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sid
    !! Index for family-kind.  
    !! 1:lda-exch, 2:lda-corr, 3:gga-exch, ...
    INTEGER, INTENT(IN) :: i_param
    !! Index of the chosen Libxc parameter
    REAL(DP), INTENT(IN) :: param
    !! Input value of the parameter
#if defined(__LIBXC)    
    par_list(sid,i_param) = param
    CALL xc_f03_func_set_ext_params( xc_func(sid), par_list(sid,:) )
#else
    CALL xclib_infomsg( 'set_libxc_ext_param', 'WARNING: an external parameter&
                         &was enforced into Libxc, but Libxc is not active' )
#endif
    RETURN
  END SUBROUTINE
  !
  !----------------------------------------------------------------------------
  FUNCTION get_libxc_ext_param( sid, i_param )
    !--------------------------------------------------------------------------
    !! Get the value of the i-th external parameter of Libxc functional with
    !! \(ID = \text{func_id}\)
    USE kind_l,       ONLY: DP
#if defined(__LIBXC)
    USE dft_setting_params,  ONLY: par_list
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: sid
    !! ID of the libxc functional
    INTEGER, INTENT(IN) :: i_param
    !! Index of the chosen parameter
    REAL(DP) :: get_libxc_ext_param
    !! Value of the parameter
#if defined(__LIBXC)
    get_libxc_ext_param = par_list(sid,i_param)
#else
    CALL xclib_infomsg( 'get_libxc_ext_param', 'WARNING: an external parameter&
                         &was sought in Libxc, but Libxc is not linked' )
    get_libxc_ext_param = 0.d0
#endif
    RETURN
  END FUNCTION
  !
  !-------------------------------------------------------------------------
  FUNCTION xclib_get_dft_short()
    !---------------------------------------------------------------------
    !! Get DFT name in short notation.
    !
    USE dft_setting_params, ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                                  is_libxc, scan_exx, notset
    USE qe_dft_list,        ONLY: dft_LDAc_name, get_shortname_from_IDs
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=32) :: xclib_get_dft_short
    CHARACTER(LEN=32) :: shortname
    INTEGER :: ID_vec(6)
    !
    shortname = 'no shortname'
    !
    ID_vec(1) = iexch  ;  ID_vec(2) = icorr
    ID_vec(3) = igcx   ;  ID_vec(4) = igcc
    ID_vec(5) = imeta  ;  ID_vec(6) = imetac
    !
    CALL get_shortname_from_IDs( ID_vec, shortname )
    !
    IF ( shortname/='no shortname' .AND. iexch==1 .AND. igcx==0 .AND. igcc==0) THEN
       shortname = TRIM(dft_LDAc_name(icorr))
    ENDIF
    !
    !
    IF ( ANY(is_libxc(5:6)) ) THEN
       IF (imeta==263 .AND. imetac==267) THEN
          shortname = 'SCAN'
          IF (scan_exx) shortname = 'SCAN0'
       ELSEIF (imeta == 208 .AND. imetac==231) THEN
          shortname = 'TB09'
       ENDIF  
    ENDIF
    !
    IF ( TRIM(shortname)=='no shortname' ) THEN
       shortname = 'XC-000I-000I-000I-000I-000I-000I'
       WRITE( shortname(4:6),   '(i3.3)' ) iexch
       IF ( is_libxc(1) ) WRITE( shortname(7:7),   '(a)' ) 'L'
       WRITE( shortname(9:11),  '(i3.3)' ) icorr
       IF ( is_libxc(2) ) WRITE( shortname(12:12), '(a)' ) 'L'
       WRITE( shortname(14:16), '(i3.3)' ) igcx
       IF ( is_libxc(3) ) WRITE( shortname(17:17), '(a)' ) 'L'
       WRITE( shortname(19:21), '(i3.3)' ) igcc
       IF ( is_libxc(4) ) WRITE( shortname(22:22), '(a)' ) 'L'
       WRITE( shortname(24:26), '(i3.3)' ) imeta
       IF ( is_libxc(5) ) WRITE( shortname(27:27), '(a)' ) 'L'
       WRITE( shortname(29:31), '(i3.3)' ) imetac
       IF ( is_libxc(6) ) WRITE( shortname(32:32), '(a)' ) 'L'
    ENDIF
    !
    xclib_get_dft_short = shortname
    !
  END FUNCTION xclib_get_dft_short
  !
  !
  !---------------------------------------------------------------------
  FUNCTION xclib_get_dft_long()
    !---------------------------------------------------------------------
    !! Get DFT name in long notation.
    !
    USE dft_setting_params, ONLY: iexch, icorr, igcx, igcc, imeta
    USE qe_dft_list,        ONLY: dft_LDAx_name, dft_LDAc_name, dft_GGAx_name, &
                                  dft_GGAc_name, dft_MGGA_name
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=25) :: xclib_get_dft_long
    CHARACTER(LEN=25) :: longname
    !
    WRITE(longname,'(4a5)') dft_LDAx_name(iexch), dft_LDAc_name(icorr), &
                            dft_GGAx_name(igcx),  dft_GGAc_name(igcc)
    !
    IF ( imeta > 0 )  longname = longname(1:20)//TRIM(dft_MGGA_name(imeta))
    !
    xclib_get_dft_long = longname
    !
  END FUNCTION xclib_get_dft_long
  !
  !---------------------------------------------------------------------------
  SUBROUTINE xclib_set_threshold( family, rho_threshold_, grho_threshold_, tau_threshold_ )
   !--------------------------------------------------------------------------
   !! Set input threshold for \(\text{family}\)-term of XC functional.
   !
   USE kind_l,             ONLY: DP
   USE dft_setting_params, ONLY: rho_threshold_lda,  rho_threshold_gga,    &
                                 rho_threshold_mgga, grho2_threshold_mgga, &
                                 grho_threshold_gga, tau_threshold_mgga
   !
   IMPLICIT NONE
   !
   CHARACTER(len=*), INTENT(IN) :: family
   !! LDA, GGA or MGGA
   REAL(DP), INTENT(IN) :: rho_threshold_
   !! Density threshold
   REAL(DP), INTENT(IN), OPTIONAL :: grho_threshold_
   !! Threshold for density gradient
   REAL(DP), INTENT(IN), OPTIONAL :: tau_threshold_
   !! Threshold for density laplacian
   !
   CHARACTER(len=4) :: cfamily
   INTEGER :: i, ln
   !
   ln = LEN_TRIM(family)
   !
   DO i = 1, ln
     cfamily(i:i) = capital(family(i:i))
   ENDDO
   !
   SELECT CASE( cfamily(1:ln) )
   CASE( 'LDA' )
     rho_threshold_lda = rho_threshold_
   CASE( 'GGA' )
     rho_threshold_gga = rho_threshold_
     IF ( PRESENT(grho_threshold_) ) grho_threshold_gga = grho_threshold_
   CASE( 'MGGA' )
     rho_threshold_mgga = rho_threshold_
     IF ( PRESENT(grho_threshold_) ) grho2_threshold_mgga = grho_threshold_
     IF ( PRESENT(tau_threshold_)  ) tau_threshold_mgga   = tau_threshold_
   END SELECT
   !
   RETURN
   !
  END SUBROUTINE xclib_set_threshold
  !
  !-----------------------------------------------------------------------
  FUNCTION matches( string1, string2 )  
   !-----------------------------------------------------------------------
   !! TRUE if string1 is contained in string2, .FALSE. otherwise
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: string1, string2
   LOGICAL :: matches
   INTEGER :: len1, len2, l  
   !
   len1 = LEN_TRIM( string1 )  
   len2 = LEN_TRIM( string2 )  
   !
   DO l = 1, ( len2 - len1 + 1 )  
     !   
     IF ( string1(1:len1) == string2(l:(l+len1-1)) ) THEN  
        !
        matches = .TRUE.  
        !
        RETURN  
        !
     END IF
     !
   END DO
   !
   matches = .FALSE.
   ! 
   RETURN
   !
  END FUNCTION matches
  !
  !-----------------------------------------------------------------------
  FUNCTION capital( in_char )  
   !-----------------------------------------------------------------------
   !! Converts character to capital if lowercase.  
   !! Copies character to output in all other cases.
   !
   IMPLICIT NONE  
   !
   CHARACTER(LEN=1), INTENT(IN) :: in_char
   CHARACTER(LEN=1)             :: capital
   CHARACTER(LEN=26), PARAMETER :: lower = 'abcdefghijklmnopqrstuvwxyz', &
                                   upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   INTEGER                      :: i
   !
   DO i=1, 26
     IF ( in_char == lower(i:i) ) THEN
       capital = upper(i:i)
       RETURN
     ENDIF
   ENDDO
   !
   capital = in_char
   !
   RETURN 
   !
  END FUNCTION capital
  !
  !
END MODULE dft_setting_routines
