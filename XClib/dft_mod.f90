!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE dft_mod
  !--------------------------------------------------------------------------
  !! Routines to set and/or recover DFT names, parameters and flags.
  !
#if defined(__LIBXC)
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
#if defined(__LIBXC)
  PUBLIC :: get_libxc_flags_exc
#endif
  !
CONTAINS
  !
  !-------------------------------------------------------------------------
  SUBROUTINE xclib_set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !! Translates a string containing the exchange-correlation name
    !! into internal indices iexch, icorr, igcx, igcc, inlc, imeta.
    !
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                            discard_input_dft, is_libxc, dft, exc,   &
                            corr, gradx, gradc, meta, nxc, ncc, ngcx,&
                            ngcc, nmeta, scan_exx, notset
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
#if defined(__LIBXC)
    INTEGER :: ii, id_vec(6), n_ext_params
    INTEGER :: flag_v(16), exp2, ftot, ftotx
    TYPE(xc_f03_func_t) :: xc_func03
    TYPE(xc_f03_func_info_t) :: xc_info03
#endif
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
    SELECT CASE( TRIM(dftout) )
    ! special cases : PZ  (LDA is equivalent to PZ)
    CASE( 'PZ', 'LDA' )
       dft_defined = xclib_set_dft_IDs(1,1,0,0,0,0)
    ! speciale cases : PW ( LDA with PW correlation )
    CASE( 'PW' )
       dft_defined = xclib_set_dft_IDs(1,4,0,0,0,0)
    ! special cases : VWN-RPA
    CASE( 'VWN-RPA' )
       dft_defined = xclib_set_dft_IDs(1,11,0,0,0,0)
    ! special cases : OEP no GC part (nor LDA...) and no correlation by default
    CASE( 'OEP' )
       dft_defined = xclib_set_dft_IDs(4,0,0,0,0,0)
    !
    CASE( 'KLI' )
       dft_defined = xclib_set_dft_IDs(10,0,0,0,0,0)
    ! special cases : HF no GC part (nor LDA...) and no correlation by default
    CASE( 'HF' )
       dft_defined = xclib_set_dft_IDs(5,0,0,0,0,0)
    ! special case : PBE
    CASE( 'PBE' )
       dft_defined = xclib_set_dft_IDs(1,4,3,4,0,0)
    ! special case : B88
    CASE( 'B88' )
       dft_defined = xclib_set_dft_IDs(1,1,1,0,0,0)
    ! special case : BP = B88 + P86
    CASE( 'BP' )
       dft_defined = xclib_set_dft_IDs(1,1,1,1,0,0)
    ! special case : PW91 = GGX + GGC
    CASE( 'PW91' )
       dft_defined = xclib_set_dft_IDs(1,4,2,2,0,0)
    ! special case : revPBE
    CASE( 'REVPBE' )
       dft_defined = xclib_set_dft_IDs(1,4,4,4,0,0)
    ! special case : PBEsol
    CASE( 'PBESOL' )
       dft_defined = xclib_set_dft_IDs(1,4,10,8,0,0)
    ! special cases : BLYP (note, BLYP=>B88)
    CASE( 'BLYP' )
       dft_defined = xclib_set_dft_IDs(1,3,1,3,0,0)
    ! Special case optB88
    CASE( 'OPTBK88' )
       dft_defined = xclib_set_dft_IDs(1,4,23,1,0,0)
    ! Special case optB86b
    CASE( 'OPTB86B' )
       dft_defined = xclib_set_dft_IDs(1,4,24,1,0,0)
    ! special case : PBC  = PW + PBC
    CASE( 'PBC' )
       dft_defined = xclib_set_dft_IDs(1,4,0,4,0,0)
    ! special case : HCTH
    CASE( 'HCTH' )
       dft_defined = xclib_set_dft_IDs(0,0,5,5,0,0)
       CALL xclib_error( 'set_dft_from_name', 'HCTH yields suspicious results', 1 )
    ! special case : OLYP = OPTX + LYP
    CASE( 'OLYP' )
       dft_defined = xclib_set_dft_IDs(0,3,6,3,0,0)
       CALL xclib_error( 'set_dft_from_name', 'OLYP yields suspicious results', 1 )
    ! special case : Wu-Cohen
    CASE( 'WC' )
       dft_defined = xclib_set_dft_IDs(1,4,11,4,0,0)
    ! special case : PW86PBE
    CASE( 'PW86PBE' )
       dft_defined = xclib_set_dft_IDs(1,4,21,4,0,0)
    ! special case : B86BPBE
    CASE( 'B86BPBE' )
       dft_defined = xclib_set_dft_IDs(1,4,22,4,0,0)
    ! special case : PBEQ2D
    CASE( 'PBEQ2D', 'Q2D' )
       dft_defined = xclib_set_dft_IDs(1,4,19,12,0,0)
    ! special case : SOGGA = SOX + PBEc
    CASE( 'SOGGA' )
       dft_defined = xclib_set_dft_IDs(1,4,17,4,0,0)
    ! special case : Engel-Vosko
    CASE( 'EV93' )
       dft_defined = xclib_set_dft_IDs(1,4,25,0,0,0)
    ! special case : RPBE
    CASE( 'RPBE' )
       CALL xclib_error( 'set_dft_from_name', &
                    'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)', 1 )
    ! special case : PBE0
    CASE( 'PBE0' )
       dft_defined = xclib_set_dft_IDs(6,4,8,4,0,0)
    ! special case : B86BPBEX
    CASE( 'B86BPBEX' )
       dft_defined = xclib_set_dft_IDs(6,4,41,4,0,0)
    !
    CASE( 'BHAHLYP', 'BHANDHLYP' )
       dft_defined = xclib_set_dft_IDs(6,4,42,3,0,0)
    ! special case : HSE
    CASE( 'HSE' )
       dft_defined = xclib_set_dft_IDs(1,4,12,4,0,0)
    ! special case : GAUPBE
    CASE( 'GAUP', 'GAUPBE' )
       dft_defined = xclib_set_dft_IDs(1,4,20,4,0,0)
    ! special case : B3LYP
    CASE( 'B3LYP' )
       dft_defined = xclib_set_dft_IDs(7,12,9,7,0,0)
    ! special case : B3LYP-VWN-1-RPA hybrid
    CASE( 'B3LYP-V1R' )
       dft_defined = xclib_set_dft_IDs(7,13,9,7,0,0)
    ! special case : X3LYP hybrid
    CASE( 'X3LYP' )
       dft_defined = xclib_set_dft_IDs(9,14,28,13,0,0)
    ! special case : TPSS meta-GGA Exc
    CASE( 'TPSS' )
       dft_defined = xclib_set_dft_IDs(1,4,7,6,1,0)
    ! special case : TPSS meta-GGA - mgga term only
    CASE( 'TPSS-only' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,1,0)
    ! special case : M06L Meta GGA
    CASE( 'M06L' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,2,0)
    ! special case : TB09 meta-GGA Exc
    CASE( 'TB09' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,3,0)
    ! special case : SCAN Meta GGA
    CASE( 'SCAN' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,5,0)
    ! special case : SCAN0
    CASE( 'SCAN0' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,6,0)
    ! special case : PZ/LDA + null meta-GGA
    CASE( 'PZ+META', 'LDA+META' )
       dft_defined = xclib_set_dft_IDs(1,1,0,0,4,0)
    ! special case : PBE + null meta-GGA
    CASE( 'PBE+META' )
       dft_defined = xclib_set_dft_IDs(1,4,3,4,4,0)
    !
    CASE DEFAULT
       !
       ! ----------------------------------------------------------------------
       ! CHECK LIBXC FUNCTIONALS BY INDEX NOTATION, IF PRESENT (USED in PHonon)
       ! ----------------------------------------------------------------------
       !
       IF (dftout(1:3) .EQ. 'XC-') THEN
#if defined(__LIBXC)
          is_libxc = .FALSE.
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
#else
          CALL xclib_error( 'set_dft_from_name', 'libxc functionals needed, but &
                                            &libxc is not active', 1 )
#endif
       ENDIF
       !
    END SELECT
    !
    !
    ! ... A temporary fix to keep the q-e input notation for SCAN-functionals
    !     valid.
#if defined(__LIBXC)
    IF (imeta==5 .OR. imeta==6) THEN
       IF (imeta==6) scan_exx = .TRUE.
       imeta  = 263 
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
    !----------------------------------------------------------------
    ! If the DFT was not yet defined, check every part of the string
    !----------------------------------------------------------------
    !
    IF (.NOT. dft_defined) THEN
       !
       is_libxc(:) = .FALSE.
       !
       iexch = matching( dftout, nxc,   exc   )
       icorr = matching( dftout, ncc,   corr  )
       igcx  = matching( dftout, ngcx,  gradx )
       igcc  = matching( dftout, ngcc,  gradc )
       imeta = matching( dftout, nmeta, meta  )
       imetac = 0
       !
    ENDIF
    !
#if defined(__LIBXC)
    IF (.NOT. dft_defined) CALL matching_libxc( dftout )
    !
    !------------------------------------------------------------------
    ! Checks whether external parameters are required by the libxc
    ! functionals (if present)
    !------------------------------------------------------------------
    !
    id_vec(1) = iexch  ;  id_vec(2) = icorr
    id_vec(3) = igcx   ;  id_vec(4) = igcc
    id_vec(5) = imeta  ;  id_vec(6) = imetac
    !
    n_ext_params = 0
    DO ii = 1, 6
      IF (is_libxc(ii)) THEN
        CALL xc_f03_func_init( xc_func03, id_vec(ii), 1 )
        xc_info03 = xc_f03_func_get_info(xc_func03)
        n_ext_params = xc_f03_func_info_get_n_ext_params(xc_info03)
        ftot = xc_f03_func_info_get_flags(xc_info03)
        flag_v(1:16) = 0
        exp2 = 16
        DO WHILE (ftot > 0)
          exp2 = exp2 - 1
          ftotx = ftot - 2**exp2
          IF (ftotx >= 0) THEN
            flag_v(exp2+1) = 1
            ftot = ftotx
          ENDIF
        ENDDO
        !
        IF ( n_ext_params /= 0 ) THEN       
          WRITE( *, '(/5X,"WARNING: libxc functional with ID ",I4," depends",&
                     &/5X," on external parameters: the correct operation in",&
                     &/5X," QE is not guaranteed with default values.")' ) id_vec(ii)
        ENDIF
        IF ( flag_v(1) == 0 ) THEN
          WRITE( *, '(/5X,"WARNING: libxc functional with ID ",I4," does not ",&
                     &/5X,"provide Exc: its correct operation in QE is not ",&
                     &/5X,"guaranteed.")' ) id_vec(ii)
        ENDIF
        IF ( flag_v(2) == 0 ) THEN
          WRITE( *, '(/5X,"WARNING: libxc functional with ID ",I4," does not ",&
                     &/5X,"provide Vxc: its correct operation in QE is not ",&
                     &/5X,"guaranteed.")' ) id_vec(ii)
        ENDIF
        IF (dftout(1:3) .EQ. 'XC-' .AND. flag_v(3) == 0 ) THEN
          WRITE( *, '(/5X,"WARNING: libxc functional with ID ",I4," does not ",&
                     &/5X,"provide Vxc derivative: its correct operation in QE is",&
                     &/5X," not guaranteed when derivative is needed.")' ) id_vec(ii)
        ENDIF
        CALL xc_f03_func_end( xc_func03 )
      ENDIF
    ENDDO
    !
#endif
    !
    ! Back compatibility - TO BE REMOVED
    !
    IF (igcx == 14) igcx = 3 ! PBE -> PBX
    IF (igcc ==  9) igcc = 4 ! PBE -> PBC
    !
    IF (igcx == 6) CALL xclib_infomsg( 'set_dft_from_name', 'OPTX untested! please test' )
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
       WRITE(*,*) iexch, save_iexch
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for iexch', 1 )
    ENDIF
    IF (save_icorr /= notset .AND. save_icorr /= icorr) THEN
       WRITE(*,*) icorr, save_icorr
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for icorr', 1 )
    ENDIF
    IF (save_igcx /= notset  .AND. save_igcx /= igcx)   THEN
       WRITE(*,*) igcx, save_igcx
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for igcx',  1 )
    ENDIF
    IF (save_igcc /= notset  .AND. save_igcc /= igcc)   THEN
       WRITE (*,*) igcc, save_igcc
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for igcc',  1 )
    ENDIF
    IF (save_meta /= notset  .AND. save_meta /= imeta)  THEN
       WRITE (*,*) imeta, save_meta
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for imeta', 1 )
    ENDIF
    IF (save_metac /= notset  .AND. save_metac /= imetac)  THEN
       WRITE (*,*) imetac, save_metac
       CALL xclib_error( 'set_dft_from_name', ' conflicting values for imetac', 1 )
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
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac
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
    USE dft_par_mod,  ONLY: notset
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
#if defined(__LIBXC)
          IF ( matching == notset ) matching = i
#else
          IF ( matching == notset ) THEN
             !WRITE(*, '("matches",i2,2X,A,2X,A)') i, name(i), TRIM(dft)
             matching = i
          ELSE
             WRITE(*, '(2(2X,i2,2X,A))') i, TRIM(name(i)), &
                                  matching, TRIM(name(matching))
             CALL xclib_error( 'set_dft', 'two conflicting matching values', 1 )
          ENDIF
#endif
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
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, is_libxc, &
                            exx_fraction
    ! 
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !
    CHARACTER(LEN=256) :: name
    INTEGER :: i, l, prev_len(6), fkind, fkind_v(3), family
    INTEGER, PARAMETER :: ID_MAX_LIBXC=600
    TYPE(xc_f03_func_t) :: xc_func
    TYPE(xc_f03_func_info_t) :: xc_info
#if (XC_MAJOR_VERSION > 5)
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
          IF ( matches('HYB_', TRIM(name)) ) THEN
            exx_fraction = xc_f03_hyb_exx_coef( xc_func )
          ENDIF
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
    ! LDA:
    IF (iexch/=0 .AND. fkind_v(1)==XC_EXCHANGE_CORRELATION)  &
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an EXCHANGE+CORRELATION &
                           &functional has been found together with an exchange&
                           & one (LDA)' )
    ! GGA:
    IF (igcx/=0 .AND. fkind_v(2)==XC_EXCHANGE_CORRELATION)   &
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an EXCHANGE+CORRELATION &
                           &functional has been found together with an exchange&
                           & one (GGA)' )
    !
    IF ( (is_libxc(3).AND.iexch/=0) .OR. (is_libxc(4).AND. icorr/=0) )    &
       CALL xclib_infomsg( 'matching_libxc', 'WARNING: an LDA functional has bee&
                           &n found, but libxc GGA functionals already include t&
                           &he LDA part' )
    ! mGGA:
    ! (imeta defines both exchange and correlation term for q-e mGGA functionals)
    IF (imeta/=0 .AND. (.NOT. is_libxc(5)) .AND. imetac/=0)   &
       CALL xclib_error( 'matching_libxc', 'Two conflicting metaGGA functionals &
                         &have been found', 1 )
    !
    IF (imeta/=0 .AND. fkind_v(3)==XC_EXCHANGE_CORRELATION)  &   
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
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, is_libxc, &
                            exc, corr, gradx, gradc, meta, nxc, ncc, ngcx,     &
                            ngcc, nmeta
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !! DFT full name
    !
    CHARACTER(LEN=4) :: qe_name
    CHARACTER(LEN=256) :: lxc_name
    INTEGER :: i, l, ch, qedft, nlxc
    INTEGER :: id_vec(6)
    !
    id_vec(1)=iexch ; id_vec(2)=icorr
    id_vec(3)=igcx  ; id_vec(4)=igcc
    id_vec(5)=imeta ; id_vec(6)=imetac
    !
    DO ch = 1, 5
       IF (.NOT.is_libxc(ch)) THEN
          !
          SELECT CASE( ch )
          CASE( 1 )
             qe_name = exc(iexch)
          CASE( 2 )
             qe_name = corr(icorr)
          CASE( 3 )
             qe_name = gradx(igcx)
          CASE( 4 )
             qe_name = gradc(igcc)
          CASE( 5 )
             qe_name = meta(imeta)
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
              lxc_name = xc_f03_functional_get_name( id_vec(i) )
              DO l = 1, LEN_TRIM(lxc_name)
                 lxc_name(l:l) = capital( lxc_name(l:l) )
              ENDDO
              IF (matches( TRIM(qe_name), TRIM(lxc_name))) nlxc = nlxc + 1
            ENDIF
          ENDDO
          !
          IF (qedft == nlxc) id_vec(ch) = 0  
          !
       ENDIF
    ENDDO
    !
    iexch = id_vec(1) ;  icorr  = id_vec(2)
    igcx  = id_vec(3) ;  igcc   = id_vec(4)
    imeta = id_vec(5) ;  imetac = id_vec(6)
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
    USE kind_l,       ONLY: DP
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                            islda, isgradient, ismeta, exx_fraction, &
                            screening_parameter, gau_parameter,      & 
                            ishybrid, has_finite_size_correction
    !  
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: isnonlocc
    !! The non-local part, for now, is not included in xc_lib, but this variable
    !! is needed to establish 'isgradient'.
    !
    ismeta    = (imeta+imetac > 0)
    isgradient= (igcx > 0) .OR.  (igcc > 0)  .OR. ismeta .OR. isnonlocc
    islda     = (iexch> 0) .AND. (icorr > 0) .AND. .NOT. isgradient
    ! PBE0/DF0
    IF ( iexch==6 .OR.  igcx == 8 ) exx_fraction = 0.25_DP
    ! CX0P
    IF ( iexch==6 .AND. igcx ==31 ) exx_fraction = 0.20_DP
    ! B86BPBEX
    IF ( iexch==6 .AND. igcx ==41 ) exx_fraction = 0.25_DP
    ! BHANDHLYP
    IF ( iexch==6 .AND. igcx ==42 ) exx_fraction = 0.50_DP
    ! HSE
    IF ( igcx ==12 ) THEN
       exx_fraction = 0.25_DP
       screening_parameter = 0.106_DP
    ENDIF
    ! gau-pbe
    IF ( igcx ==20 ) THEN
       exx_fraction = 0.24_DP
       gau_parameter = 0.150_DP
    ENDIF
    ! HF or OEP
    IF ( iexch==4 .OR. iexch==5 ) exx_fraction = 1.0_DP
    ! B3LYP or B3LYP-VWN-1-RPA
    IF ( iexch == 7 ) exx_fraction = 0.2_DP
    ! X3LYP
    IF ( iexch == 9 ) exx_fraction = 0.218_DP
    !
    ishybrid = ( exx_fraction /= 0.0_DP )
    !
    has_finite_size_correction = ( iexch==8 .OR. icorr==10)
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
!     USE dft_par_mod
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
!     USE dft_par_mod
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
    USE dft_par_mod, ONLY: ishybrid, exx_started
    IMPLICIT NONE
    IF (.NOT. ishybrid) &
       CALL xclib_error( 'start_exx', 'dft is not hybrid, wrong call', 1 )
    exx_started = .TRUE.
  END SUBROUTINE start_exx
  !-----------------------------------------------------------------------
  SUBROUTINE stop_exx
    !! Deactivate exact exchange (exx_started=FALSE)
    USE dft_par_mod, ONLY: ishybrid, exx_started
    IMPLICIT NONE
    IF (.NOT. ishybrid) &
       CALL xclib_error( 'stop_exx', 'dft is not hybrid, wrong call', 1 )
    exx_started = .FALSE.
  END SUBROUTINE stop_exx
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_set_exx_fraction( exx_fraction_ )
    !! Impose input parameter as exact exchange fraction value
    USE kind_l,      ONLY: DP
    USE dft_par_mod, ONLY: exx_fraction
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: exx_fraction_
    !! Imposed value of exact exchange fraction
    exx_fraction = exx_fraction_
    WRITE( *,'(5x,a,f6.2)') 'EXX fraction changed: ', exx_fraction
    RETURN
  END SUBROUTINE xclib_set_exx_fraction
  !-----------------------------------------------------------------------
  SUBROUTINE dft_force_hybrid( request )
    !! Impose hybrid condition.
    USE dft_par_mod, ONLY: ishybrid
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
     USE dft_par_mod, ONLY: exx_started
     IMPLICIT NONE
     LOGICAL :: exx_is_active
     exx_is_active = exx_started
  END FUNCTION exx_is_active
  !-----------------------------------------------------------------------
  FUNCTION xclib_get_exx_fraction()
     !! Recover exact exchange fraction.
     USE kind_l,      ONLY: DP
     USE dft_par_mod, ONLY: exx_fraction
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
    USE kind_l,      ONLY: DP
    USE dft_par_mod, ONLY: screening_parameter
    IMPLICIT NONE
    REAL(DP):: scrparm_
    !! Value to impose as screening parameter
    screening_parameter = scrparm_
    WRITE( *,'(5x,a,f12.7)') 'EXX Screening parameter changed: ', &
         & screening_parameter
  END SUBROUTINE set_screening_parameter
  !-----------------------------------------------------------------------
  FUNCTION get_screening_parameter()
     !! Recover screening parameter (for pbexsr)
     USE kind_l,      ONLY: DP
     USE dft_par_mod, ONLY: screening_parameter
     IMPLICIT NONE
     REAL(DP):: get_screening_parameter
     get_screening_parameter = screening_parameter
     RETURN
  END FUNCTION get_screening_parameter
  !-----------------------------------------------------------------------
  SUBROUTINE set_gau_parameter( gauparm_ )
    !! Impose input parameter as gau parameter (for gau-pbe)
    USE kind_l,      ONLY: DP
    USE dft_par_mod, ONLY: gau_parameter
    IMPLICIT NONE
    REAL(DP):: gauparm_
    !! Value to impose as gau parameter
    gau_parameter = gauparm_
    WRITE( *,'(5x,a,f12.7)') 'EXX Gau parameter changed: ', &
         & gau_parameter
  END SUBROUTINE set_gau_parameter
  !-----------------------------------------------------------------------
  FUNCTION get_gau_parameter()
    !! Recover gau parameter (for gau-pbe)
    USE kind_l,      ONLY: DP
    USE dft_par_mod, ONLY: gau_parameter
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
  FUNCTION xclib_get_id( family, kindf )
     !--------------------------------------------------------------------
     !! Get functionals index of \(\text{family}\) and \(\text{kind}\).
     !
     USE dft_par_mod, ONLY: iexch, icorr, igcx, igcc, imeta, imetac
     !
     IMPLICIT NONE
     !
     INTEGER :: xclib_get_id
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
       IF (ckindf=='EXCH') xclib_get_id = iexch
       IF (ckindf=='CORR') xclib_get_id = icorr
     CASE( 'GGA' )
       IF (ckindf=='EXCH') xclib_get_id = igcx
       IF (ckindf=='CORR') xclib_get_id = igcc
     CASE( 'MGGA' )
       IF (ckindf=='EXCH') xclib_get_id = imeta
       IF (ckindf=='CORR') xclib_get_id = imetac
     CASE DEFAULT
       CALL xclib_error( 'xclib_get_id', 'input not recognized', 1 )
     END SELECT
     !
     RETURN
     !
  END FUNCTION xclib_get_id
  !
  !-------------------------------------------------------------------       
  SUBROUTINE xclib_get_name( family, kindf, name )
     !----------------------------------------------------------------
     !! Gets QE name for 'family'-'kind' term of the XC functional.
     !
     USE dft_par_mod, ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                            exc, corr, gradx, gradc, meta
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
       IF (ckindf=='EXCH') name = exc(iexch)
       IF (ckindf=='CORR') name = corr(icorr)
     CASE( 'GGA' )
       IF (ckindf=='EXCH') name = gradx(igcx)
       IF (ckindf=='CORR') name = gradc(igcc)
     CASE( 'MGGA' )
       IF (ckindf=='EXCH') name = meta(imeta) 
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
     USE dft_par_mod,  ONLY: is_libxc
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
       IF (TRIM(cfamily)=='ANY'.AND.ANY(is_libxc(:))) xclib_dft_is_libxc=.TRUE.
     ENDIF
     !
     RETURN
     !
  END FUNCTION
  !
  ! ... previous format:
  !
!   FUNCTION get_iexch()
!      INTEGER get_iexch
!      get_iexch = iexch
!      RETURN
!   END FUNCTION get_iexch
!   !-----------------------------------------------------------------------
!   FUNCTION get_icorr()
!      INTEGER get_icorr
!      get_icorr = icorr
!      RETURN
!   END FUNCTION get_icorr
!   !-----------------------------------------------------------------------
!   FUNCTION get_igcx()
!      INTEGER get_igcx
!      get_igcx = igcx
!      RETURN
!   END FUNCTION get_igcx
!   !-----------------------------------------------------------------------
!   FUNCTION get_igcc()
!      INTEGER get_igcc
!      get_igcc = igcc
!      RETURN
!   END FUNCTION get_igcc
!   !-----------------------------------------------------------------------
!   FUNCTION get_meta()
!      INTEGER get_meta
!      get_meta = imeta
!      RETURN
!   END FUNCTION get_meta
!   !
!   FUNCTION get_metac()
!     INTEGER get_metac
!     get_metac = imetac
!     RETURN
!   END FUNCTION get_metac
  !
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_reset_dft()
    !---------------------------------------------------------------------
    !! Unset DFT indexes.
    USE dft_par_mod, ONLY: iexch, icorr, igcx, igcc, imeta, imetac, notset
    IMPLICIT NONE
    iexch  = notset ; icorr  = notset
    igcx   = notset ; igcc   = notset
    imeta  = notset ; imetac = notset
  END SUBROUTINE
  !
  !------------------------------------------------------------------------
  FUNCTION get_dft_name()
     !---------------------------------------------------------------------
     !! Get full DFT name
     USE dft_par_mod, ONLY: dft
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
     USE dft_par_mod,  ONLY: isgradient, ismeta, ishybrid
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
  ! ... previous format:
  !
!   FUNCTION dft_is_gradient()
!      LOGICAL :: dft_is_gradient
!      dft_is_gradient = isgradient
!      RETURN
!   END FUNCTION dft_is_gradient
!   !-----------------------------------------------------------------------
!   FUNCTION dft_is_meta()
!      LOGICAL :: dft_is_meta
!      dft_is_meta = ismeta
!      RETURN
!   END FUNCTION dft_is_meta
!   !-----------------------------------------------------------------------
!   FUNCTION dft_is_hybrid()
!      LOGICAL :: dft_is_hybrid
!      dft_is_hybrid = ishybrid
!      RETURN
!   END FUNCTION dft_is_hybrid
  !
  !-----------------------------------------------------------------------
  FUNCTION igcc_is_lyp()
     !! Find if correlation GGA is Lee-Yang-Parr.
     USE dft_par_mod,  ONLY: igcc
     IMPLICIT NONE
     LOGICAL :: igcc_is_lyp
     igcc_is_lyp = (igcc==3 .OR. igcc==7 .OR. igcc==13)
     RETURN
  END FUNCTION igcc_is_lyp
  !-----------------------------------------------------------------------
  FUNCTION dft_has_finite_size_correction()
     !! TRUE if finite size correction present
     USE dft_par_mod, ONLY: has_finite_size_correction
     IMPLICIT NONE
     LOGICAL :: dft_has_finite_size_correction
     dft_has_finite_size_correction = has_finite_size_correction
     RETURN
  END FUNCTION dft_has_finite_size_correction
  !-----------------------------------------------------------------------
  SUBROUTINE xclib_set_finite_size_volume( volume )
     !! Set value for finite size cell volume.
     USE dft_par_mod, ONLY: has_finite_size_correction, finite_size_cell_volume,&
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
     USE kind_l,      ONLY: DP
     USE dft_par_mod, ONLY: finite_size_cell_volume, finite_size_cell_volume_set
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
  SUBROUTINE xclib_init_libxc( xclib_nspin )
    !------------------------------------------------------------------------
    !! Initialize Libxc functionals, if present.
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                            is_libxc, libxc_initialized
#if defined(__LIBXC)
    USE dft_par_mod,  ONLY: n_ext_params, xc_func, xc_info, par_list
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: xclib_nspin
    !! 1: unpolarized case; 2: polarized
    INTEGER :: iid, ip
    INTEGER :: id_vec(6)
    !
#if defined(__LIBXC)
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
        CALL xc_f03_func_init( xc_func(iid), id_vec(iid), xclib_nspin )
        xc_info(iid) = xc_f03_func_get_info( xc_func(iid) )
        n_ext_params(iid) = xc_f03_func_info_get_n_ext_params( xc_info(iid) )
        DO ip = 1, n_ext_params(iid)
          par_list(iid,ip) = xc_f03_func_info_get_ext_params_default_value( &
                                                           xc_info(iid), ip )
        ENDDO
        libxc_initialized(iid) = .TRUE.
      ENDIF
    ENDDO
#endif
    RETURN
  END SUBROUTINE xclib_init_libxc
  !
  !--------------------------------------------------------------------------
  SUBROUTINE xclib_finalize_libxc()
    !------------------------------------------------------------------------
    !! Finalize Libxc functionals, if present.
    USE dft_par_mod,  ONLY: iexch, icorr, igcx, igcc, imeta, imetac, &
                            is_libxc
#if defined(__LIBXC)
    USE dft_par_mod,  ONLY: xc_func
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
      IF (is_libxc(iid)) CALL xc_f03_func_end( xc_func(iid) )
    ENDDO
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
    USE dft_par_mod,  ONLY: xc_func, par_list
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
    USE dft_par_mod,  ONLY: par_list
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
                         &was sought in Libxc, but Libxc is not active' )
#endif
    RETURN
  END FUNCTION
  !
#if defined(__LIBXC)
  !------------------------------------------------------------------------
  SUBROUTINE get_libxc_flags_exc( xc_info, eflag )
     !--------------------------------------------------------------------
     !! Checks whether Exc is present or not in the output of a libxc 
     !! functional (e.g. TB09 and a few others)
     IMPLICIT NONE
     TYPE(xc_f03_func_info_t) :: xc_info
     INTEGER :: ii, flags_tot
     INTEGER, INTENT(OUT) :: eflag 
     flags_tot = xc_f03_func_info_get_flags(xc_info)
     eflag = 0
     DO ii = 15, 0, -1
       IF ( flags_tot-2**ii<0 ) CYCLE
       flags_tot = flags_tot-2**ii
       IF ( ii==0 ) eflag = 1
     ENDDO
     RETURN
  END SUBROUTINE
#endif
  !
  !-------------------------------------------------------------------------
  FUNCTION xclib_get_dft_short()
    !---------------------------------------------------------------------
    !! Get DFT name in short notation.
    !
    USE dft_par_mod, ONLY: iexch, icorr, igcx, igcc, imeta, imetac, corr, &
                           is_libxc, scan_exx, notset
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=32) :: xclib_get_dft_short
    CHARACTER(LEN=32) :: shortname
    !
    shortname = 'no shortname'
    !
    IF ( iexch==1 .AND. igcx==0 .AND. igcc==0) THEN
       shortname = TRIM(corr(icorr))
    ELSEIF (iexch==4 .AND. icorr==0  .AND. igcx==0  .AND. igcc== 0) THEN
       shortname = 'OEP'
    ELSEIF (iexch==1 .AND. icorr==11 .AND. igcx==0  .AND. igcc== 0) THEN
       shortname = 'VWN-RPA'
    ELSEIF (iexch==1 .AND. icorr==3  .AND. igcx==1  .AND. igcc== 3) THEN
       shortname = 'BLYP'
    ELSEIF (iexch==1 .AND. icorr==1  .AND. igcx==1  .AND. igcc== 0) THEN
       shortname = 'B88'
    ELSEIF (iexch==1 .AND. icorr==1  .AND. igcx==1  .AND. igcc== 1) THEN
       shortname = 'BP'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==2  .AND. igcc== 2) THEN
       shortname = 'PW91'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==3  .AND. igcc== 4) THEN
       shortname = 'PBE'
    ELSEIF (iexch==6 .AND. icorr==4  .AND. igcx==8  .AND. igcc== 4) THEN
       shortname = 'PBE0'
    ELSEIF (iexch==6 .AND. icorr==4  .AND. igcx==41 .AND. igcc== 4) THEN
       shortname = 'B86BPBEX'
    ELSEIF (iexch==6 .AND. icorr==4  .AND. igcx==42 .AND. igcc== 3) THEN
       shortname = 'BHANDHLYP'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==4  .AND. igcc== 4) THEN
       shortname = 'revPBE'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==10 .AND. igcc== 8) THEN
       shortname = 'PBESOL'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==19 .AND. igcc==12) THEN
       shortname = 'Q2D'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==12 .AND. igcc== 4) THEN
       shortname = 'HSE'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==20 .AND. igcc== 4) THEN
       shortname = 'GAUPBE'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==21 .AND. igcc== 4) THEN
       shortname = 'PW86PBE'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==22 .AND. igcc== 4) THEN
       shortname = 'B86BPBE'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==11 .AND. igcc== 4) THEN
       shortname = 'WC'
    ELSEIF (iexch==7 .AND. icorr==12 .AND. igcx==9  .AND. igcc== 7) THEN
       shortname = 'B3LYP'
    ELSEIF (iexch==7 .AND. icorr==13 .AND. igcx==9  .AND. igcc== 7) THEN
       shortname = 'B3LYP-V1R'
    ELSEIF (iexch==9 .AND. icorr==14 .AND. igcx==28 .AND. igcc==13) THEN
       shortname = 'X3LYP'
    ELSEIF (iexch==0 .AND. icorr==3  .AND. igcx==6  .AND. igcc== 3) THEN
       shortname = 'OLYP'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==17 .AND. igcc== 4) THEN
       shortname = 'SOGGA'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==23 .AND. igcc== 1) THEN
       shortname = 'OPTBK88'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==24 .AND. igcc== 1) THEN
       shortname = 'OPTB86B'
    ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==25 .AND. igcc== 0) THEN
       shortname = 'EV93'
    ELSEIF (iexch==5 .AND. icorr==0  .AND. igcx==0  .AND. igcc== 0) THEN
       shortname = 'HF'
    ENDIF
    !
    IF (imeta==1) THEN
       shortname = 'TPSS'
    ELSEIF (imeta == 2) THEN
       shortname = 'M06L'
    ELSEIF (imeta == 4) THEN
       IF ( iexch == 1 .AND. icorr == 1) THEN
          shortname = 'PZ+META'
       ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==3 .AND. igcc==4) THEN
          shortname = 'PBE+META'
       ENDIF
    ENDIF
    !
#if defined(__LIBXC)
    IF ( ANY(is_libxc(:)) ) THEN
       IF (imeta==263 .AND. imetac==267) THEN
          shortname = 'SCAN'
          IF (scan_exx) shortname = 'SCAN0'
       ELSEIF (imeta == 208 .AND. imetac==231) THEN
          shortname = 'TB09'
       ELSE
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
    ENDIF
#endif
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
    USE dft_par_mod, ONLY: iexch, icorr, igcx, igcc, imeta, exc, corr, gradx,&
                           gradc, meta
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=25) :: xclib_get_dft_long
    CHARACTER(LEN=25) :: longname
    !
    WRITE(longname,'(4a5)') exc(iexch), corr(icorr), gradx(igcx), gradc(igcc)
    !
    IF ( imeta > 0 )  longname = longname(1:20)//TRIM(meta(imeta))
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
   USE kind_l,      ONLY: DP
   USE dft_par_mod, ONLY: rho_threshold_lda, rho_threshold_gga, grho2_threshold_mgga, &
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
     rho_threshold_gga = rho_threshold_
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
END MODULE dft_mod
