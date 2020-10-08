
!=---------------------------------------------------------------------------=!
MODULE xc_interfaces
  !
  IMPLICIT NONE
  PRIVATE
  !
  ! LDA
  PUBLIC :: XC, DMXC
  ! GGA
  PUBLIC :: XC_GCX, DGCXC
  ! MGGA
  PUBLIC :: XC_METAGCX
  ! 
  PUBLIC :: XCLIB_SET_THRESHOLD
  !
  PUBLIC :: XCLIB_SET_DFT_FROM_NAME, XCLIB_SET_DFT_IDs, XCLIB_GET_ID, &
            XCLIB_SET_AUXILIARY_FLAGS, XCLIB_DFT_IS, XCLIB_IS_LIBXC, &
            START_EXX, STOP_EXX, EXX_IS_ACTIVE, XCLIB_SET_EXX_FRACTION, &
            XCLIB_GET_EXX_FRACTION, XCLIB_RESET_DFT, &
            XCLIB_SET_FINITE_SIZE_VOLUME, XCLIB_GET_FINITE_SIZE_CELL_VOLUME, &
            XCLIB_GET_DFT_SHORT, XCLIB_GET_DFT_LONG, &
            get_screening_parameter, get_gau_parameter, &
            set_screening_parameter, set_gau_parameter, dft_force_hybrid, &
            igcc_is_lyp, xclib_get_name, dft_has_finite_size_correction
  
  !
  !
  !
  INTERFACE XCLIB_SET_DFT_FROM_NAME
    !
    SUBROUTINE set_dft_from_name( dft_ )
      !
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: dft_
      !
    END SUBROUTINE
    !
  END INTERFACE
  !
  !
  INTERFACE XCLIB_SET_DFT_IDs
     !
     LOGICAL FUNCTION set_dft_IDs( iexch_, icorr_, igcx_, igcc_, imeta_, imetac_) !, is_libxc_ )
       !
       USE dft_par_mod
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iexch_, icorr_
       INTEGER, INTENT(IN) :: igcx_, igcc_
       INTEGER, INTENT(IN) :: imeta_, imetac_
       !LOGICAL, OPTIONAL, INTENT(IN) :: is_libxc_(6)
       !
     END FUNCTION
     !
  END INTERFACE
  !
  
  INTERFACE XCLIB_GET_ID
    !
    FUNCTION get_id( family, kindf )
      !
      IMPLICIT NONE
      INTEGER :: get_id
      CHARACTER(len=*), INTENT(IN) :: family, kindf
      !
    END FUNCTION
    !
  END INTERFACE
  
  
  INTERFACE XCLIB_GET_NAME
    !
    FUNCTION get_name( family, kindf )
      !
      IMPLICIT NONE
      CHARACTER(len=*) :: get_name
      CHARACTER(len=*), INTENT(IN) :: family, kindf
      !
    END FUNCTION
    !
  END INTERFACE
  
  
  INTERFACE XCLIB_IS_LIBXC
    !
    FUNCTION is_lxc( family, kindf )
      !
      IMPLICIT NONE
      LOGICAL :: is_lxc
      CHARACTER(len=*), INTENT(IN) :: family, kindf
      !
    END FUNCTION
    !
  END INTERFACE
  
  
  !
  INTERFACE XCLIB_SET_AUXILIARY_FLAGS
     !
     SUBROUTINE set_auxiliary_flags
       !
       IMPLICIT NONE
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE XCLIB_GET_DFT_SHORT
     !
     FUNCTION get_dft_short()
       !
       USE dft_par_mod
       IMPLICIT NONE
       CHARACTER(LEN=32) :: get_dft_short
       !
     END FUNCTION
     !
  END INTERFACE
  !
  !
  INTERFACE XCLIB_GET_DFT_LONG
     !
     FUNCTION get_dft_long()
       !
       USE dft_par_mod
       IMPLICIT NONE
       CHARACTER(LEN=25) :: get_dft_long
       !
     END FUNCTION
     !
  END INTERFACE
  !
  
  INTERFACE XCLIB_RESET_DFT
     !
     SUBROUTINE reset_dft()
       !
       USE dft_par_mod
       IMPLICIT NONE
       !
     END SUBROUTINE
     !
  END INTERFACE
  
  !
  INTERFACE get_screening_parameter
    FUNCTION get_screening_parameter()
      USE kind_l, ONLY: DP
      IMPLICIT NONE
      REAL(DP):: get_screening_parameter
    END FUNCTION
  END INTERFACE

  INTERFACE get_gau_parameter
    FUNCTION get_gau_parameter()
      USE kind_l, ONLY: DP
      IMPLICIT NONE
      REAL(DP):: get_gau_parameter
    END FUNCTION
  END INTERFACE
  
  INTERFACE set_screening_parameter
    SUBROUTINE set_screening_parameter(scrparm_)
      USE kind_l, ONLY: DP
      IMPLICIT NONE
      REAL(DP), INTENT(IN):: scrparm_
    END SUBROUTINE
  END INTERFACE

  INTERFACE set_gau_parameter
    FUNCTION set_gau_parameter()
      USE kind_l, ONLY: DP
      IMPLICIT NONE
      REAL(DP):: set_gau_parameter
    END FUNCTION
  END INTERFACE
  
  
  INTERFACE igcc_is_lyp
    FUNCTION igcc_is_lyp()
      USE kind_l, ONLY: DP
      IMPLICIT NONE
      LOGICAL :: igcc_is_lyp
    END FUNCTION
  END INTERFACE
  
  
!==================================
  INTERFACE START_EXX
     !
     SUBROUTINE start_exx()
       !
       IMPLICIT NONE
       !
     END SUBROUTINE
     !
  END INTERFACE
  
  !
  INTERFACE STOP_EXX
     !
     SUBROUTINE stop_exx()
       !
       IMPLICIT NONE
       !
     END SUBROUTINE
     !
  END INTERFACE 
!===============================================
  !
  INTERFACE EXX_IS_ACTIVE
     !
     LOGICAL FUNCTION exx_is_active()
       !
       IMPLICIT NONE
       !
     END FUNCTION
     !
  END INTERFACE 
  !
  INTERFACE XCLIB_GET_EXX_FRACTION
     !
     FUNCTION get_exx_fraction()
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP) :: get_exx_fraction
       !
     END FUNCTION
     !
  END INTERFACE
  !
  INTERFACE XCLIB_SET_EXX_FRACTION
     !
     SUBROUTINE set_exx_fraction( exx_fraction_ )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: exx_fraction_
       !
     END SUBROUTINE
     !
  END INTERFACE  
  !
  !
  INTERFACE XCLIB_SET_FINITE_SIZE_VOLUME
     !
     SUBROUTINE set_finite_size_volume( volume )
       IMPLICIT NONE
       REAL, INTENT(IN) :: volume
     END SUBROUTINE set_finite_size_volume
     !
  END INTERFACE
  !
  INTERFACE XCLIB_GET_FINITE_SIZE_CELL_VOLUME
     !
     SUBROUTINE get_finite_size_cell_volume( is_present, volume )
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       LOGICAL, INTENT(OUT) :: is_present
       REAL(DP), INTENT(OUT) :: volume
     END SUBROUTINE get_finite_size_cell_volume
     !
  END INTERFACE
  
  
  INTERFACE dft_force_hybrid
     !
     SUBROUTINE dft_force_hybrid( request )
       IMPLICIT NONE
       LOGICAL,OPTIONAL,INTENT(INOUT) :: request
     END SUBROUTINE
     !
  END INTERFACE
  
  INTERFACE dft_has_finite_size_correction
     !
     LOGICAL FUNCTION dft_has_finite_size_correction()
       IMPLICIT NONE
     END FUNCTION
     !
  END INTERFACE
  
  !
  
  INTERFACE XCLIB_DFT_IS
     !
     FUNCTION dft_is( what )
       !
       IMPLICIT NONE
       LOGICAL :: dft_is
       CHARACTER(len=*) :: what
       !
     END FUNCTION
     !
  END INTERFACE
  !
  !
  INTERFACE XCLIB_SET_THRESHOLD
     !
     SUBROUTINE set_threshold_l( fkind, rho_threshold_, grho_threshold_, tau_threshold_ )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       CHARACTER(len=*), INTENT(IN) :: fkind
       REAL(DP), INTENT(IN) :: rho_threshold_
       REAL(DP), INTENT(IN), OPTIONAL :: grho_threshold_
       REAL(DP), INTENT(IN), OPTIONAL :: tau_threshold_
       !
     END SUBROUTINE
     !
  END INTERFACE
  
  !--------------------------------------------------------
  !
  INTERFACE XC
     !
     SUBROUTINE xc_l( length, sr_d, sv_d, rho_in, ex_out, ec_out, vx_out, vc_out )
       !
       USE dft_par_mod
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, sr_d, sv_d
       REAL(DP), INTENT(IN) :: rho_in(length,sr_d)
       REAL(DP), INTENT(OUT) :: ex_out(length), ec_out(length)
       REAL(DP), INTENT(OUT) :: vx_out(length,sv_d), vc_out(length,sv_d)
       !
     END SUBROUTINE xc_l
     !
  END INTERFACE
  !
  INTERFACE XC_GCX
     !
     SUBROUTINE xc_gcx_l( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
       !
       USE kind_l,        ONLY: DP  
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, ns
       REAL(DP), INTENT(IN) :: rho(:,:), grho(:,:,:)
       REAL(DP), INTENT(OUT) :: ex(:), ec(:)
       REAL(DP), INTENT(OUT) :: v1x(:,:), v2x(:,:)
       REAL(DP), INTENT(OUT) :: v1c(:,:), v2c(:,:)
       REAL(DP), OPTIONAL, INTENT(OUT) :: v2c_ud(:)
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE XC_METAGCX
     !
     SUBROUTINE xc_metagcx_l( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
       !
       USE kind_l,        ONLY: DP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: length, ns, np
       REAL(DP), INTENT(IN) :: rho(length,ns), grho(3,length,ns), tau(length,ns)
       REAL(DP), INTENT(OUT) :: ex(length), ec(length)
       REAL(DP), INTENT(OUT) :: v1x(length,ns), v2x(length,ns), v3x(length,ns)
       REAL(DP), INTENT(OUT) :: v1c(length,ns), v2c(np,length,ns), v3c(length,ns)
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE DMXC
     !
     SUBROUTINE dmxc_l( length, sr_d, rho_in, dmuxc )
       !
       USE kind_l,            ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, sr_d
       REAL(DP), INTENT(IN) :: rho_in(length,sr_d)
       REAL(DP), INTENT(OUT) :: dmuxc(length,sr_d,sr_d)
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE DGCXC
     !
     SUBROUTINE dgcxc_l( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
       !
       USE kind_l,           ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, sp
       REAL(DP), INTENT(IN) :: r_in(length,sp), g_in(length,3,sp)
       REAL(DP), INTENT(OUT) :: dvxc_rr(length,sp,sp), dvxc_sr(length,sp,sp), &
                                dvxc_ss(length,sp,sp)
       !                    
     END SUBROUTINE
     !
  END INTERFACE
  !
END MODULE xc_interfaces
!=---------------------------------------------------------------------------=!
