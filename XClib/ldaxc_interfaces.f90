!
! ... LDAlib
!
!=---------------------------------------------------------------------------=!
MODULE ldaxc_interfaces
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: XC_LDA, XC_LSDA, SLATER, SLATER_SPIN, PZ, PZ_POLARIZED, PW, PW_SPIN, LYP, LSD_LYP
  PUBLIC :: GET_LDAXC_INDEXES, GET_LDAXC_PARAM, GET_LDA_THRESHOLD
  !
  !
  INTERFACE GET_LDAXC_INDEXES
     !
     SUBROUTINE get_ldaxclib( iexch_, icorr_ )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: iexch_, icorr_
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
       IMPLICIT NONE
       REAL(DP), OPTIONAL, INTENT(IN) :: finite_size_cell_volume_
       LOGICAL , OPTIONAL, INTENT(IN) :: exx_started_
       REAL(DP), OPTIONAL, INTENT(IN) :: exx_fraction_
       !
     END SUBROUTINE
     !
  END INTERFACE
  !
  INTERFACE GET_LDA_THRESHOLD
     !
     SUBROUTINE get_threshold( rho_threshold_ )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: rho_threshold_
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
END MODULE ldaxc_interfaces
!=---------------------------------------------------------------------------=!
