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
  !
  ! --- To be run on a single processor ---
  !
  USE kind_l,               ONLY: DP
  USE xc_lib,               ONLY: xclib_set_dft_from_name, xclib_get_ID, &
                                  xclib_dft_is_libxc, xclib_init_libxc,  &
                                  xclib_finalize_libxc
  USE qe_dft_list
  USE xclib_utils_and_para, ONLY: stdout
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, xc_info
#endif
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=120) :: lxc_kind, lxc_family
  CHARACTER(LEN=150) :: dft_r
  CHARACTER(LEN=10) :: dft_n
  INTEGER :: n_ext, id(6)
  INTEGER :: i, ii
#if defined(__LIBXC)
#if (XC_MAJOR_VERSION>5)
  !workaround to keep compatibility with libxc develop version
  INTEGER, PARAMETER :: XC_FAMILY_HYB_GGA  = -10
  INTEGER, PARAMETER :: XC_FAMILY_HYB_MGGA = -11 
#endif
#endif
  !
  !-------- Input var -----------------------
  CHARACTER(LEN=80) :: dft
  !
  !---------- DFT infos -------------------------
  INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac, idx
  LOGICAL :: is_libxc(6)
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
  CALL xclib_set_dft_from_name( dft )
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
  WRITE(stdout,*) " "  
  WRITE(stdout,*) "=================================== "//CHAR(10)//" "  
  WRITE(stdout,*) "XC functional IDs:"  
  WRITE(stdout,*) CHAR(10)//"LDA IDs"  
  WRITE(stdout,121) iexch, is_libxc(1), icorr, is_libxc(2)  
  WRITE(stdout,*) CHAR(10)//"GGA IDs"  
  WRITE(stdout,121) igcx, is_libxc(3), igcc, is_libxc(4)  
  WRITE(stdout,*) CHAR(10)//"MGGA IDs"  
  WRITE(stdout,121) imeta, is_libxc(5), imetac, is_libxc(6) 
  WRITE(stdout,*) " "  
  WRITE(stdout,*) "============== "//CHAR(10)//" " 
  ! 
  !  
#if defined(__LIBXC)
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( 1 )  
#endif
  !  
  !WRITE(stdout,*) CHAR(10)//"LIBXC functional infos:"  
  !  
  id(1) = iexch ; id(2) = icorr
  id(3) = igcx  ; id(4) = igcc
  id(5) = imeta ; id(6) = imetac
  !  
  
  DO i = 1, 6  
    idx = id(i)
    IF (.NOT.is_libxc(i) .AND. idx/=0) THEN
      WRITE(stdout,*) CHAR(10)//"Functional with ID:", idx
      !
      SELECT CASE( i )
      CASE( 1 ) 
        WRITE(lxc_kind, '(a)') 'Exchange functional'
        WRITE(lxc_family,'(a)') "LDA"
        dft_n = dft_LDAx_name(idx)
        dft_r = dft_LDAx_ref(idx)
      CASE( 2 )
        WRITE(lxc_kind, '(a)') 'Correlation functional'
        WRITE(lxc_family,'(a)') "LDA"
        dft_n = dft_LDAc_name(idx)
        dft_r = dft_LDAc_ref(idx)
      CASE( 3 )
        WRITE(lxc_kind, '(a)') 'Exchange functional'
        WRITE(lxc_family,'(a)') "GGA"
        dft_n = dft_GGAx_name(idx)
        dft_r = dft_GGAx_ref(idx)
      CASE( 4 )
        WRITE(lxc_kind, '(a)') 'Correlation functional'
        WRITE(lxc_family,'(a)') "GGA"
        dft_n = dft_GGAc_name(idx)
        dft_r = dft_GGAc_ref(idx)
      CASE( 5 )
        WRITE(lxc_kind, '(a)') 'Exchange+Correlation functional'
        WRITE(lxc_family,'(a)') "MGGA"
        dft_n = dft_MGGA_name(idx)
        dft_r = dft_MGGA_ref(idx)
      !CASE( 6 )
      !  WRITE(lxc_kind, '(a)') 'Correlation functional'
      !  WRITE(lxc_family,'(a)') "MGGA"
      !  dft_n = dft_MGGA_name(idx)
      !  dft_r = dft_MGGA_ref(idx)
      END SELECT
      !  
      WRITE(*,'("The functional ''", a, "'' is a ", a, ", it belongs to &  
             &the ''", a, "'' family and is defined in the reference(s): &  
             &")') TRIM(dft_n), TRIM(lxc_kind)&  
             ,TRIM(lxc_family)  
      
      !DO WHILE( ii >= 0 )  
        WRITE(*,'(a,i1,2a)') '[',1,'] ',TRIM(dft_r)  
      !ENDDO
      !
#if defined(__LIBXC)
      !
    ELSEIF (is_libxc(i)) THEN  
      !
      WRITE(stdout,*) CHAR(10)//"Functional with ID:", id(i)  
      !  
      SELECT CASE( xc_f03_func_info_get_kind(xc_info(i)) )  
      CASE( XC_EXCHANGE )  
        WRITE(lxc_kind, '(a)') 'Exchange functional'  
      CASE( XC_CORRELATION )  
        WRITE(lxc_kind, '(a)') 'Correlation functional'  
      CASE( XC_EXCHANGE_CORRELATION )  
        WRITE(lxc_kind, '(a)') 'Exchange+Correlation functional'  
      CASE( XC_KINETIC )  
        WRITE(lxc_kind, '(a)') 'Kinetic energy functional - not implemented&  
                               &in QE.'  
      CASE DEFAULT  
        WRITE(lxc_kind, '(a)') 'Unknown kind'  
      END SELECT  
      !  
      SELECT CASE( xc_f03_func_info_get_family(xc_info(i)) )  
      CASE( XC_FAMILY_LDA )  
        WRITE(lxc_family,'(a)') "LDA"  
      CASE( XC_FAMILY_GGA )  
        WRITE(lxc_family,'(a)') "GGA"  
      CASE( XC_FAMILY_HYB_GGA )  
        WRITE(lxc_family,'(a)') "Hybrid GGA"  
      CASE( XC_FAMILY_MGGA )  
        WRITE(lxc_family,'(a)') "MGGA"  
      CASE( XC_FAMILY_HYB_MGGA )  
        WRITE(lxc_family,'(a)') "Hybrid MGGA"  
      CASE DEFAULT  
        WRITE(lxc_family,'(a)') "unknown"  
      END SELECT  
      !  
      WRITE(*,'("The functional ''", a, "'' is a ", a, ", it belongs to &  
             &the ''", a, "'' family and is defined in the reference(s): &  
             &")') TRIM(xc_f03_func_info_get_name(xc_info(i))), TRIM(lxc_kind)&  
             ,TRIM(lxc_family)  
      ii = 0  
      DO WHILE( ii >= 0 )  
       WRITE(*,'(a,i1,2a)') '[',ii+1,'] ',TRIM(xc_f03_func_reference_get_ref( &  
                                 xc_f03_func_info_get_references(xc_info(i), ii)))  
      ENDDO  
      !  
      WRITE(stdout,*)  
      n_ext = xc_f03_func_info_get_n_ext_params( xc_info(i) )  
      WRITE(stdout,*) 'Number of external parameters: ', n_ext  
      !  
      IF ( n_ext/=0 ) THEN  
        DO ii = 0, n_ext-1  
          WRITE(stdout,*) &  
            TRIM(xc_f03_func_info_get_ext_params_description(xc_info(i), ii))  
          WRITE(stdout,*) 'Default value: ', &  
                 xc_f03_func_info_get_ext_params_default_value(xc_info(i), ii)  
        ENDDO  
      ENDIF
#endif
      !
    ENDIF
  ENDDO  
  !  
#if defined(__LIBXC)
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_finalize_libxc()  
#endif
  !
  121 FORMAT('Exch: ',I3,' is libxc: ',L1,';  Corr: ',I3,' is libxc: ',L1 )
  !
  STOP
  !
END PROGRAM xc_infos
