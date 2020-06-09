!
! ... LDAlib
!
!=---------------------------------------------------------------------------=!
MODULE xc_interfaces
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: XC_LDA, XC_LSDA, DMXC_LDA, DMXC_LSDA, DMXC_NC
  PUBLIC :: SLATER, SLATER_SPIN, PZ, PZ_POLARIZED, PW, PW_SPIN, LYP, &
            LSD_LYP
  PUBLIC :: GCXC, GCX_SPIN, GCC_SPIN, GCC_SPIN_MORE, DGCXC_UNPOL
  PUBLIC :: PBEC, PBEC_SPIN, LSD_GLYP
  PUBLIC :: GET_XC_INDEXES, GET_LDAXC_PARAM, GET_LDA_THRESHOLD, &
            GET_GGAXC_PARAM, GET_GGA_THRESHOLD
  !
  !
  INTERFACE GET_XC_INDEXES
     !
     SUBROUTINE get_xclib_IDs( iexch_, icorr_, igcx_, igcc_ )
       !
       USE dft_par_mod
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iexch_, icorr_
       INTEGER, INTENT(IN) :: igcx_, igcc_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE GET_LDAXC_PARAM
     !
     SUBROUTINE get_ldaxcparlib( finite_size_cell_volume_, exx_started_, exx_fraction_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       REAL(DP), OPTIONAL, INTENT(IN) :: finite_size_cell_volume_
       LOGICAL , OPTIONAL, INTENT(IN) :: exx_started_
       REAL(DP), OPTIONAL, INTENT(IN) :: exx_fraction_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE GET_GGAXC_PARAM
     !
     SUBROUTINE get_ggaxcparlib( gau_scr_par_, exx_started_, exx_fraction_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: gau_scr_par_
       LOGICAL , OPTIONAL, INTENT(IN) :: exx_started_
       REAL(DP), OPTIONAL, INTENT(IN) :: exx_fraction_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE GET_LDA_THRESHOLD
     !
     SUBROUTINE get_lda_threshold( rho_threshold_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rho_threshold_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE GET_GGA_THRESHOLD
     !
     SUBROUTINE get_GGA_threshold( rho_threshold_, grho_threshold_ )
       !
       USE kind_l,  ONLY: DP
       USE dft_par_mod
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rho_threshold_
       REAL(DP), INTENT(IN) :: grho_threshold_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !--------------------------------------------------------
  !
  INTERFACE XC_LDA
     !
     SUBROUTINE xc_lda_l( length, rho_in, ex_out, ec_out, vx_out, vc_out )
       !
       USE dft_par_mod
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN)  :: length
       REAL(DP), INTENT(IN)  :: rho_in(length)
       REAL(DP), INTENT(OUT) :: ec_out(length), vc_out(length), &
                                ex_out(length), vx_out(length)
       !
     END SUBROUTINE xc_lda_l
     !
  END INTERFACE
  !
  !
  INTERFACE XC_LSDA
     !
     SUBROUTINE xc_lsda_l( length, rho_in, zeta_in, ex_out, ec_out, vx_out, vc_out )
       !
       USE dft_par_mod
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN)  :: length
       REAL(DP), INTENT(IN)  :: rho_in(length), zeta_in(length)
       REAL(DP), INTENT(OUT) :: ex_out(length), ec_out(length), &
                                vx_out(length,2), vc_out(length,2)
       !
     END SUBROUTINE xc_lsda_l
     !
  END INTERFACE
  !
  !
  INTERFACE DMXC_LDA
     !
     SUBROUTINE dmxc_lda_l( length, rho_in, dmuxc )
       !
       USE dft_par_mod
       USE exch_lda_l,   ONLY: slater_l
       USE kind_l,       ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
       REAL(DP), INTENT(OUT), DIMENSION(length) :: dmuxc
       !
     END SUBROUTINE dmxc_lda_l
     !
  END INTERFACE
  !
  INTERFACE DMXC_LSDA
     !
     SUBROUTINE dmxc_lsda_l( length, rho_in, dmuxc )
       !
       USE dft_par_mod
       USE exch_lda_l,   ONLY: slater_l
       USE corr_lda_l,   ONLY: pz_l, pz_polarized_l
       USE kind_l,       ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho_in
       REAL(DP), INTENT(OUT), DIMENSION(length,2,2) :: dmuxc
       !
     END SUBROUTINE dmxc_lsda_l
     !
  END INTERFACE
  !
  INTERFACE DMXC_NC
     !
     SUBROUTINE dmxc_nc_l( length, rho_in, m, dmuxc )
       !
       USE dft_par_mod
       USE kind_l,       ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
       REAL(DP), INTENT(IN), DIMENSION(length,3) :: m
       REAL(DP), INTENT(OUT), DIMENSION(length,4,4) :: dmuxc
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
  INTERFACE GCC_SPIN_MORE
     !
     SUBROUTINE gcc_spin_more_l( length, rho_in, grho_in, grho_ud_in, &
                                                 sc, v1c, v2c, v2c_ud )
       !
       USE kind_l, ONLY: DP
       USE dft_par_mod
       USE corr_gga_l
       IMPLICIT NONE 
       INTEGER, INTENT(IN) :: length
       REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho_in, grho_in
       REAL(DP), INTENT(IN), DIMENSION(length) :: grho_ud_in
       REAL(DP), INTENT(OUT), DIMENSION(length) :: sc
       REAL(DP), INTENT(OUT), DIMENSION(length,2) :: v1c, v2c
       REAL(DP), INTENT(OUT), DIMENSION(length) :: v2c_ud
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE DGCXC_UNPOL
     !
     SUBROUTINE dgcxc_unpol_l( length, r_in, s2_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
       USE kind_l,         ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       REAL(DP), INTENT(IN),  DIMENSION(length) :: r_in, s2_in
       REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrx, vsrx, vssx
       REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrc, vsrc, vssc
       !
     END SUBROUTINE
     !
  END INTERFACE
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
  INTERFACE PZ
     !
     SUBROUTINE pz_ext( rs, iflag, ec, vc )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ec, vc
       INTEGER,  INTENT(IN)  :: iflag
       !
     END SUBROUTINE pz_ext
     !
  END INTERFACE
  !
  INTERFACE PZ_POLARIZED
     !
     SUBROUTINE pz_polarized_ext( rs, ec, vc )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ec, vc
       !
     END SUBROUTINE pz_polarized_ext
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
  INTERFACE PBEC
     !
     SUBROUTINE pbec_ext( rho, grho, iflag, sc, v1c, v2c )
       !
       USE kind_l,    ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: iflag
       REAL(DP), INTENT(IN) :: rho, grho
       REAL(DP), INTENT(OUT) :: sc, v1c, v2c
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  !
  INTERFACE PBEC_SPIN
     !
     SUBROUTINE pbec_spin_ext( rho, zeta, grho, iflag, sc, v1c_up, v1c_dw, v2c )                    !<GPU:DEVICE>
       !
       USE corr_lda_l, ONLY: pw_spin_l
       USE kind_l,     ONLY : DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: iflag        !<GPU:VALUE>
       REAL(DP), INTENT(IN) :: rho, zeta, grho
       REAL(DP), INTENT(OUT) :: sc, v1c_up, v1c_dw, v2c
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
END MODULE xc_interfaces
!=---------------------------------------------------------------------------=!
