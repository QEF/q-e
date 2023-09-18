!
! Copyright (C) 2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!============================================================================
!============================================================================
PROGRAM xc_infos
  !==========================================================================
  !! Provides infos on the input DFTs (both QE and Libxc).  
  !! Currently does not cover vdW functionals.
  !
  ! --- To be run on a single processor ---
  !
  USE kind_l,               ONLY: DP
  USE xc_lib,               ONLY: xclib_set_dft_from_name, xclib_get_ID, &
                                  xclib_dft_is_libxc, xclib_init_libxc,  &
                                  xclib_finalize_libxc, xclib_set_auxiliary_flags
  USE qe_dft_list
  USE qe_dft_refs
  USE dft_setting_params,   ONLY: ishybrid, exx_fraction, screening_parameter, &
                                  gau_parameter
  USE xclib_utils_and_para, ONLY: stdout, nowarning
#if defined(__LIBXC)
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, xc_info, xc_kind_error, n_ext_params, &
                                  par_list, libxc_flags
#endif
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=120) :: lxc_kind, lxc_family
  CHARACTER(LEN=150) :: dft_r
  CHARACTER(LEN=100) :: dft_w
  CHARACTER(LEN=10)  :: dft_n
  INTEGER :: n_ext, id(6), idfull, fkind
  INTEGER :: i, ii
  !
  !-------- Input var -----------------------
  CHARACTER(LEN=80) :: dft
  !
  !---------- DFT infos -------------------------
  INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac, idx
  LOGICAL :: is_libxc(6)
  CHARACTER(LEN=80) :: name_check
  !
  dft = 'none'
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Insert DFT name:  "
  READ(*,'(A)') dft
  !
  !==========================================================================
  ! PRINT DFT INFOS
  !==========================================================================
  !
  nowarning = .TRUE.
  !
  CALL xclib_set_dft_from_name( dft )
  !
#if defined(__LIBXC)
  IF ( xc_kind_error ) WRITE(stdout,*) 'WARNING: This functional includes Libxc &
                                        &terms that are currently not usable in &
                                        &QE (kinetic) and they will be ignored.'
#endif
  !
  iexch = xclib_get_ID('LDA','EXCH')
  is_libxc(1) = xclib_dft_is_libxc('LDA','EXCH')
  icorr = xclib_get_ID('LDA','CORR')
  is_libxc(2) = xclib_dft_is_libxc('LDA','CORR')
  igcx = xclib_get_ID('GGA','EXCH')
  is_libxc(3) = xclib_dft_is_libxc('GGA','EXCH')
  igcc = xclib_get_ID('GGA','CORR')
  is_libxc(4) = xclib_dft_is_libxc('GGA','CORR')
  imeta = xclib_get_ID('MGGA','EXCH')
  is_libxc(5) = xclib_dft_is_libxc('MGGA','EXCH')
  imetac = xclib_get_ID('MGGA','CORR')
  is_libxc(6) = xclib_dft_is_libxc('MGGA','CORR')
  !
  id(1) = iexch ; id(2) = icorr
  id(3) = igcx  ; id(4) = igcc
  id(5) = imeta ; id(6) = imetac
  !
  name_check = 'noshortname'
  CALL get_shortname_from_IDs( id, name_check, idfull )
  !
  WRITE(stdout,*) " "
  WRITE(stdout,*) "=================================== "//CHAR(10)//" "
  !
  IF (TRIM(name_check)/='noshortname') THEN
    WRITE(stdout,*) dft_full_descr(idfull)
    WRITE(stdout,*) CHAR(10)
  ENDIF
  !
  WRITE(stdout,*) "The selected XC functional is a composition of the &
                  &following terms:"
  WRITE(stdout,*) CHAR(10)//"LDA"
  WRITE(stdout,121) iexch, TRIM(xc_library(is_libxc(1),iexch)), &
                    icorr, TRIM(xc_library(is_libxc(2),icorr))
  WRITE(stdout,*) CHAR(10)//"GGA"
  WRITE(stdout,121) igcx,  TRIM(xc_library(is_libxc(3),igcx)),  &
                    igcc,  TRIM(xc_library(is_libxc(4),igcc))
  WRITE(stdout,*) CHAR(10)//"MGGA"
  WRITE(stdout,121) imeta, TRIM(xc_library(is_libxc(5),imeta)), &
                    imetac,TRIM(xc_library(is_libxc(6),imetac))
  WRITE(stdout,*) " "
  WRITE(stdout,*) "============== "
  !
  CALL xclib_set_auxiliary_flags( .FALSE. )
  !
#if defined(__LIBXC)
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( 1, .FALSE. )
#endif
  !
  DO i = 1, 6
    idx = id(i)
    !
    IF (.NOT.is_libxc(i) .AND. idx/=0) THEN
      !
      SELECT CASE( i )
      CASE( 1 )
        WRITE(lxc_kind, '(a)') 'EXCHANGE'
        WRITE(lxc_family,'(a)') "LDA"
        dft_n = dft_LDAx_name(idx)
        dft_r = dft_LDAx(idx)%ref
        dft_w = dft_LDAx(idx)%wrn
      CASE( 2 )
        WRITE(lxc_kind, '(a)') 'CORRELATION'
        WRITE(lxc_family,'(a)') "LDA"
        dft_n = dft_LDAc_name(idx)
        dft_r = dft_LDAc(idx)%ref
        dft_w = dft_LDAc(idx)%wrn
      CASE( 3 )
        WRITE(lxc_kind, '(a)') 'EXCHANGE'
        IF (ishybrid) WRITE(lxc_family,'(a)') "Hybrid GGA"
        IF (.NOT. ishybrid) WRITE(lxc_family,'(a)') "GGA"
        dft_n = dft_GGAx_name(idx)
        dft_r = dft_GGAx(idx)%ref
        dft_w = dft_GGAx(idx)%wrn
      CASE( 4 )
        WRITE(lxc_kind, '(a)') 'CORRELATION'
        IF (ishybrid) WRITE(lxc_family,'(a)') "Hybrid GGA"
        IF (.NOT. ishybrid) WRITE(lxc_family,'(a)') "GGA"
        dft_n = dft_GGAc_name(idx)
        dft_r = dft_GGAc(idx)%ref
        dft_w = dft_GGAc(idx)%wrn
      CASE( 5 )
        WRITE(lxc_kind, '(a)') 'EXCHANGE+CORRELATION'
        IF (ishybrid) WRITE(lxc_family,'(a)') "Hybrid MGGA"
        IF (.NOT. ishybrid) WRITE(lxc_family,'(a)') "MGGA"
        dft_n = dft_MGGA_name(idx)
        dft_r = dft_MGGA(idx)%ref
        dft_w = dft_MGGA(idx)%wrn
      END SELECT
      !
      WRITE(stdout,*) CHAR(10)
      WRITE(*,'(i1,". Functional with ID:", i3 )') i, idx
      IF (TRIM(dft_n)/='xxxx') THEN
        WRITE(stdout, '(" - Name:   ",a)') TRIM(dft_n)
      ELSE
        WRITE(stdout, '(" - Name:   xxxx [",a,"]")') TRIM(dft)
      ENDIF
      WRITE(stdout, '(" - Family: ",a)') TRIM(lxc_family)
      WRITE(stdout, '(" - Kind:   ",a)') TRIM(lxc_kind)
      !
      WRITE(stdout, '(" - Warnings:")')
      WRITE(stdout,'(a,2a)') '    ', TRIM(dft_w)
      !
      n_ext = 0
      IF ( (ishybrid .AND. (i==1 .OR. i==3)) .OR. &
           (i==3 .AND. idx==12) .OR. (i==3 .AND. idx==20) ) n_ext = 1
      WRITE(stdout, '(" - External parameters:")')
      IF ( n_ext/=0 ) THEN
        IF ( ishybrid .AND. (i==1 .OR. i==3)) WRITE(stdout,*) '   exx_fraction (default)= ', exx_fraction
        IF ( i==3 .AND. idx==12 ) WRITE(stdout,*) '   screening_parameter (default)= ',&
                                                  screening_parameter
        IF ( i==3 .AND. idx==20 ) WRITE(stdout,*) '   gau_parameter (default)= ',      &
                                                  gau_parameter
      ELSE
        WRITE(stdout, '("    none")')
      ENDIF
      WRITE(stdout, '(" - Reference(s):")')
      WRITE(stdout,'(a,i1,2a)') '    [',1,'] ', TRIM(dft_r)
      !
#if defined(__LIBXC)
      !
    ELSEIF (is_libxc(i)) THEN
      !
      fkind = xc_f03_func_info_get_kind(xc_info(i))
      SELECT CASE( fkind )
      CASE( XC_EXCHANGE )
        WRITE(lxc_kind, '(a)') 'EXCHANGE'
      CASE( XC_CORRELATION )
        WRITE(lxc_kind, '(a)') 'CORRELATION'
      CASE( XC_EXCHANGE_CORRELATION )
        WRITE(lxc_kind, '(a)') 'EXCHANGE+CORRELATION'
      CASE( XC_KINETIC )
        WRITE(lxc_kind, '(a)') 'KINETIC ENERGY FUNCTIONAL - currently NOT usable&
                               & in QE.'
      CASE DEFAULT
        WRITE(lxc_kind, '(a)') 'UNKNOWN'
      END SELECT
      !
      SELECT CASE( xc_f03_func_info_get_family(xc_info(i)) )
      CASE( XC_FAMILY_LDA )
        WRITE(lxc_family,'(a)') "LDA"
      CASE( XC_FAMILY_HYB_LDA )
        WRITE(lxc_family,'(a)') "Hybrid LDA" 
      CASE( XC_FAMILY_GGA )
        WRITE(lxc_family,'(a)') "GGA"
      CASE( XC_FAMILY_HYB_GGA )
        WRITE(lxc_family,'(a)') "Hybrid GGA"
      CASE( XC_FAMILY_MGGA )
        WRITE(lxc_family,'(a)') "MGGA"
      CASE( XC_FAMILY_HYB_MGGA )
        WRITE(lxc_family,'(a)') "Hybrid MGGA"
      CASE DEFAULT
        WRITE(lxc_family,'(a)') "UNKNOWN"
      END SELECT
      !
      WRITE(stdout,*) CHAR(10)
      WRITE(*,'(i1,". Functional with ID: ", i3 )') i, idx
      WRITE(stdout, '(" - Name:   ",a)') TRIM(xc_f03_func_info_get_name(xc_info(i)))
      WRITE(stdout, '(" - Family: ",a)') TRIM(lxc_family)
      WRITE(stdout, '(" - Kind:   ",a)') TRIM(lxc_kind)
      !
      IF (lxc_family(1:6)=="Hybrid" .AND. (MOD(i,2)==1 .OR. ((MOD(i,2)==0) &
                                    .AND.fkind==XC_EXCHANGE_CORRELATION))) THEN
         WRITE(stdout,*) '- Default exx fraction: ', exx_fraction
         IF (screening_parameter/=0.d0) &
             WRITE(stdout,*) '- Default screening parameter: ', screening_parameter
      ENDIF
      IF ( n_ext_params(i)/=0 ) THEN
        WRITE(stdout, '(" - External parameters: ",i3)') n_ext_params(i)
        DO ii = 0, n_ext_params(i)-1
          WRITE(stdout, '("  ",i3,") ",a)') ii,&
            TRIM(xc_f03_func_info_get_ext_params_description(xc_info(i), ii))  
          WRITE(stdout,*) '      Default value: ', par_list(i,ii+1)
        ENDDO
      ELSE
        WRITE(stdout, '(" - External parameters: NONE")')
      ENDIF
      !
      WRITE(stdout, '(" - Special warnings: ")')
      IF ( libxc_flags(i,0)  == 0 ) &
        WRITE(stdout,'(4X,"[w00] libxc functional with ID ",I4," does not ",&
                      &/4X,"provide Exc.")' ) idx
      IF ( libxc_flags(i,1)  == 0 ) &
        WRITE(stdout,'(4X,"[w01] libxc functional with ID ",I4," does not ",&
                      &/4X,"provide Vxc.")' ) idx
      IF ( libxc_flags(i,2)  == 0 ) &
        WRITE(stdout,'(4X,"[w02] libxc functional with ID ",I4," does not ", &
                      &/4X,"provide Vxc derivative.")' ) idx
      IF ( libxc_flags(i,8) == 1 ) &
            WRITE(stdout,'(/5X,"WARNING: libxc functional with ID ",I4," is CAM, ", &
                          &/5X,"long range exx is not available yet (if needed).")' ) idx
      IF ( libxc_flags(i,14) == 1 ) &
        WRITE(stdout,'(4X,"[w02] libxc functional with ID ",I4," is still ", &
                      &/4X,"in development.")' ) idx
      IF ( libxc_flags(i,15) == 1 ) &
        WRITE(stdout,'(4X,"[w15] libxc functional with ID ",I4," depends on", &
                      &/4X," the laplacian of the density, which is currently set",&
                      &/4X," to zero.")' ) idx
      IF ( ALL(libxc_flags(i,0:2)==1) .AND. libxc_flags(i,14)==0 .AND. libxc_flags(i,15)==0 ) &
        WRITE(stdout, '(4X,"NONE")')
      !
      WRITE(stdout, '(" - Reference(s):")') 
      ii = 0  
      DO WHILE( ii >= 0 )  
        WRITE(stdout,'(a,i1,2a)') '    [',ii+1,'] ',TRIM(xc_f03_func_reference_get_ref( &  
                                  xc_f03_func_info_get_references(xc_info(i), ii)))  
      ENDDO
#endif
      !
    ENDIF
  ENDDO  
  !  
#if defined(__LIBXC)
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_finalize_libxc()  
#endif
  !
  WRITE(stdout,*) CHAR(10)//" "
  !
  121 FORMAT( 'Exchange ID: ', i3, ', Library: ', a, ' ;  Correlation ID: ', i3, ', Library: ',a )
  !
  STOP
  !
 CONTAINS
  !
  CHARACTER(11) FUNCTION xc_library( islibxc, idxc )
    !
    LOGICAL, INTENT(IN) :: islibxc
    INTEGER, INTENT(IN) :: idxc
    !
    xc_library = ''
    IF (idxc /= 0) THEN
      IF ( islibxc ) THEN
        xc_library = 'Libxc'
      ELSE
        xc_library = 'QE_internal'
      ENDIF
    ELSE
      xc_library = 'none'
    ENDIF
    !
    RETURN
    !
  END FUNCTION
  !
END PROGRAM xc_infos
