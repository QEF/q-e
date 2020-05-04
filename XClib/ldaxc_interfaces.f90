!
! ... LDAlib
!
!=---------------------------------------------------------------------------=!
MODULE ldaxc_interfaces
  !
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: XC_LDA, SLATER, PZ, PW
  PUBLIC :: GET_LDAXC_INDEXES, GET_LDAXC_PARAM, GET_LDA_THRESHOLD
  !
  !
  INTERFACE GET_LDAXC_INDEXES
     !
     SUBROUTINE get_ldaxclib( iexch_, icorr_ )
       !
       !USE dft_par_mod
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
       !USE dft_par_mod
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
       !USE dft_par_mod
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
  !--------------------------------------------------------
  !
  INTERFACE SLATER
     !
     SUBROUTINE slater_l( rs, ex, vx )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ex
       REAL(DP), INTENT(OUT) :: vx
       !
     END SUBROUTINE slater_l
     !
  END INTERFACE
  !
  !
  INTERFACE PZ
     !
     SUBROUTINE pz_l( rs, iflag, ec, vc )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ec, vc
       INTEGER,  INTENT(IN)  :: iflag
       !
     END SUBROUTINE pz_l
     !
  END INTERFACE
  !
  !
  INTERFACE PW
     !
     SUBROUTINE pw_l( rs, iflag, ec, vc )
       !
       USE kind_l,  ONLY: DP
       IMPLICIT NONE
       REAL(DP), INTENT(IN)  :: rs
       REAL(DP), INTENT(OUT) :: ec, vc
       INTEGER,  INTENT(IN)  :: iflag
       !
     END SUBROUTINE pw_l
     !
  END INTERFACE
  !
  !
END MODULE ldaxc_interfaces
!=---------------------------------------------------------------------------=!
