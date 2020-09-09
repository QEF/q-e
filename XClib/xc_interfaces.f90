!
! ... LDAlib
!
!=---------------------------------------------------------------------------=!
MODULE xc_interfaces
  !
  IMPLICIT NONE
  PRIVATE
  !
  ! LDA
  PUBLIC :: XC, DMXC
  PUBLIC :: SLATER, SLATER_SPIN, PW, PW_SPIN, LYP, &
            LSD_LYP
  ! GGA
  PUBLIC :: XC_GCX, GCXC, GCX_SPIN, GCC_SPIN, DGCXC
  PUBLIC :: LSD_GLYP
  ! MGGA
  PUBLIC :: XC_METAGCX, TAU_XC, TAU_XC_SPIN
  PUBLIC :: TPSSCXC
  ! 
  PUBLIC :: XCLIB_GET_IDs, XCLIB_GET_EXX, XCLIB_GET_FINITE_SIZE_CELL_VOL, &
            XCLIB_SET_THRESHOLD, XCLIB_GET_GAU_SCR_PARAM
  !
  !
  INTERFACE XCLIB_GET_IDs
     !
     SUBROUTINE get_IDs( iexch_, icorr_, igcx_, igcc_, imeta_, imetac_, is_libxc_ )
       !
       USE dft_par_mod
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iexch_, icorr_
       INTEGER, INTENT(IN) :: igcx_, igcc_
       INTEGER, INTENT(IN) :: imeta_, imetac_
       LOGICAL, OPTIONAL, INTENT(IN) :: is_libxc_(6)
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE XCLIB_GET_EXX
     !
     SUBROUTINE get_exx_started_l( exx_started_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       LOGICAL, INTENT(IN) :: exx_started_
       !
     END SUBROUTINE
     !
     SUBROUTINE get_exx_fraction_l( exx_fraction_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: exx_fraction_
       !
     END SUBROUTINE
     !
  END INTERFACE  
  !
  INTERFACE XCLIB_GET_FINITE_SIZE_CELL_VOL
     !
     SUBROUTINE get_finite_size_cell_l( finite_size_cell_volume_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: finite_size_cell_volume_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE XCLIB_GET_GAU_SCR_PARAM
     !
     SUBROUTINE get_gau_scr_par_l( gau_scr_par_ )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: gau_scr_par_
       !
     END SUBROUTINE
     !
  END INTERFACE
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
  INTERFACE GCXC
     !
     SUBROUTINE gcxc_l( length, rho_in, grho_in, sx_out, sc_out, v1x_out, &
                                          v2x_out, v1c_out, v2c_out )
       USE kind_l, ONLY: DP
       USE dft_par_mod
       USE exch_gga_l
       USE corr_gga_l
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length 
       REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in, grho_in
       REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_out, sc_out, v1x_out, &
                                                   v2x_out, v1c_out, v2c_out
     END SUBROUTINE gcxc_l
     !
  END INTERFACE
  !
  !
  INTERFACE GCX_SPIN
     !
     SUBROUTINE gcx_spin_l( length, rho_in, grho2_in, sx_tot, v1x_out, v2x_out )
       !
       USE kind_l, ONLY: DP
       USE dft_par_mod
       USE exch_gga_l
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: length
       REAL(DP), INTENT(IN),  DIMENSION(length,2) :: rho_in, grho2_in
       REAL(DP), INTENT(OUT), DIMENSION(length) :: sx_tot
       REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1x_out, v2x_out
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE GCC_SPIN
     !
     SUBROUTINE gcc_spin_l( length, rho_in, zeta_io, grho_in, sc_out, v1c_out, v2c_out )
       !
       USE kind_l, ONLY: DP
       USE dft_par_mod
       USE corr_gga_l
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: length
       REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
       REAL(DP), INTENT(INOUT), DIMENSION(length) :: zeta_io
       REAL(DP), INTENT(IN), DIMENSION(length) :: grho_in
       REAL(DP), INTENT(OUT), DIMENSION(length) :: sc_out
       REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c_out
       REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_out
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
  !
  INTERFACE TAU_XC
     !
     SUBROUTINE tau_xc_l( length, rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
       !
       USE kind_l
       USE dft_par_mod
       USE metagga_l
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: length
       REAL(DP), DIMENSION(length) :: rho, grho2, tau, &
                                      ex, ec, v1x, v2x, v3x, v1c, v2c, v3c
       !
     END SUBROUTINE tau_xc_l
     !
  END INTERFACE
  !
  INTERFACE TAU_XC_SPIN
     !
     SUBROUTINE tau_xc_spin_l( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
       !
       USE kind_l
       USE dft_par_mod
       USE metagga_l
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       REAL(DP), INTENT(IN) :: rho(length,2), tau(length,2)
       REAL(DP), INTENT(IN) :: grho(3,length,2)
       REAL(DP), INTENT(OUT) :: ex(length), ec(length), v1x(length,2), v2x(length,2), &
                                v3x(length,2), v1c(length,2), v3c(length,2)
       REAL(DP), INTENT(OUT) :: v2c(3,length,2)
       !
     END SUBROUTINE tau_xc_spin_l
     !
  END INTERFACE
  !
  !
  !---PROVISIONAL .. for cases when functional routines are called outside xc-drivers---
  !
  INTERFACE SLATER
     !
     SUBROUTINE slater_ext( rs, ex, vx )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ex
       REAL(DP), INTENT(OUT) :: vx
       !
     END SUBROUTINE slater_ext
     !
  END INTERFACE
  !
  INTERFACE SLATER_SPIN
     !
     SUBROUTINE slater_spin_ext( rho, zeta, ex, vx_up, vx_dw )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rho, zeta
       REAL(DP), INTENT(OUT) :: ex, vx_up, vx_dw
       !
     END SUBROUTINE slater_spin_ext
     !
  END INTERFACE
  !
  !
  INTERFACE PW
     !
     SUBROUTINE pw_ext( rs, iflag, ec, vc )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ec, vc
       INTEGER,  INTENT(IN)  :: iflag
       !
     END SUBROUTINE pw_ext
     !
  END INTERFACE
  !
  INTERFACE PW_SPIN
     !
     SUBROUTINE pw_spin_ext( rs, zeta, ec, vc_up, vc_dw )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs, zeta
       REAL(DP), INTENT(OUT) :: ec, vc_up, vc_dw
       !
     END SUBROUTINE pw_spin_ext
     !
  END INTERFACE
  !
  !
  INTERFACE LYP
     !
     SUBROUTINE lyp_ext( rs, ec, vc )
       !
       USE kind_l,      ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rs
       REAL(DP), INTENT(OUT) :: ec, vc
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE LSD_LYP
     !
     SUBROUTINE lsd_lyp_ext( rho, zeta, elyp, vlyp_up, vlyp_dw )
       !
       USE kind_l,       ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rho, zeta
       REAL(DP), INTENT(OUT) :: elyp, vlyp_up, vlyp_dw
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE LSD_GLYP
     !
     SUBROUTINE lsd_glyp_ext( rho_in_up, rho_in_dw, grho_up, grho_dw, grho_ud, sc, v1c_up, v1c_dw, v2c_up, v2c_dw, v2c_ud )                     !<GPU:DEVICE>
       !
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rho_in_up, rho_in_dw
       REAL(DP), INTENT(IN) :: grho_up, grho_dw
       REAL(DP), INTENT(IN) :: grho_ud
       REAL(DP), INTENT(OUT) :: sc
       REAL(DP), INTENT(OUT) :: v1c_up, v1c_dw
       REAL(DP), INTENT(OUT) :: v2c_up, v2c_dw
       REAL(DP), INTENT(OUT) :: v2c_ud
       !
     END SUBROUTINE
     !
  END INTERFACE   
  !
  !
  INTERFACE TPSSCXC
     !
     SUBROUTINE tpsscxc_ext( rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )
       !
       USE kind_l,      ONLY : DP
       USE metagga_l  , ONLY : metax, metac
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rho, grho, tau
       REAL(DP), INTENT(OUT) :: sx, sc
       REAL(DP), INTENT(OUT) :: v1x, v2x, v3x
       REAL(DP), INTENT(OUT) :: v1c, v2c, v3c
     END SUBROUTINE
     !
  END INTERFACE   
  !
END MODULE xc_interfaces
!=---------------------------------------------------------------------------=!
