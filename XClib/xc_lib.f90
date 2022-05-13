!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------
MODULE xc_lib
  !----------------------------------------------------
  !! Interface module for \(\texttt{xc_lib}\) library.
  !
  USE dft_setting_routines
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  !
  PUBLIC :: xc, dmxc          !LDA
  !
  PUBLIC :: xc_gcx, dgcxc     !GGA
  !
  PUBLIC :: xc_metagcx        !MGGA
  !
  PUBLIC :: xclib_set_dft_from_name,           &
            xclib_set_dft_IDs,                 &
            xclib_set_auxiliary_flags,         &
            xclib_set_threshold,               &
            xclib_set_exx_fraction,            &
            xclib_set_finite_size_volume,      &
            set_screening_parameter,           &
            set_gau_parameter
  !
  PUBLIC :: xclib_get_name,                    &
            xclib_get_ID,                      & 
            xclib_get_dft_short,               &
            xclib_get_dft_long,                &
            xclib_get_exx_fraction,            &
            xclib_get_finite_size_cell_volume, &
            get_screening_parameter,           &
            get_gau_parameter        
  !
  PUBLIC :: xclib_dft_is,                      &
            xclib_dft_is_libxc,                &
            xclib_init_libxc,                  &
            xclib_finalize_libxc,              &
            set_libxc_ext_param,               &
            get_libxc_ext_param,               &
            start_exx, stop_exx,               &
            dft_has_finite_size_correction,    &
            exx_is_active,                     &
            igcc_is_lyp,                       &
            xclib_reset_dft,                   &
            dft_force_hybrid
  !
  !
  INTERFACE xc
     SUBROUTINE xc( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out, gpu_args_ )
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, srd, svd
       REAL(DP), INTENT(IN) :: rho_in(length,srd)
       REAL(DP), INTENT(OUT) :: ex_out(length), ec_out(length)
       REAL(DP), INTENT(OUT) :: vx_out(length,svd), vc_out(length,svd)
       LOGICAL,  OPTIONAL, INTENT(IN) :: gpu_args_
     END SUBROUTINE
  END INTERFACE
  !
  !
  INTERFACE xc_gcx
     SUBROUTINE xc_gcx( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud, &
                        gpu_args_ )
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, ns
       REAL(DP), INTENT(IN) :: rho(length,ns), grho(3,length,ns)
       REAL(DP), INTENT(OUT) :: ex(length), ec(length)
       REAL(DP), INTENT(OUT) :: v1x(length,ns), v2x(length,ns)
       REAL(DP), INTENT(OUT) :: v1c(length,ns), v2c(length,ns)
       REAL(DP), OPTIONAL, INTENT(OUT) :: v2c_ud(length)
       LOGICAL,  OPTIONAL, INTENT(IN)  :: gpu_args_
     END SUBROUTINE
  END INTERFACE
  !
  !
  INTERFACE xc_metagcx
     SUBROUTINE xc_metagcx( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, &
                            v1c, v2c, v3c, gpu_args_ )
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, ns, np
       REAL(DP), INTENT(IN) :: rho(length,ns), grho(3,length,ns), tau(length,ns)
       REAL(DP), INTENT(OUT) :: ex(length), ec(length)
       REAL(DP), INTENT(OUT) :: v1x(length,ns), v2x(length,ns), v3x(length,ns)
       REAL(DP), INTENT(OUT) :: v1c(length,ns), v2c(np,length,ns), v3c(length,ns)
       LOGICAL,  OPTIONAL, INTENT(IN) :: gpu_args_
     END SUBROUTINE
  END INTERFACE
  !
  !
  INTERFACE dmxc
     SUBROUTINE dmxc( length, srd, rho_in, dmuxc, gpu_args_ )
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       INTEGER,  INTENT(IN) :: srd
       REAL(DP), INTENT(IN) :: rho_in(length,srd)
       REAL(DP), INTENT(OUT) :: dmuxc(length,srd,srd)
       LOGICAL,  OPTIONAL, INTENT(IN) :: gpu_args_
     END SUBROUTINE
  END INTERFACE
  !
  !
  INTERFACE dgcxc
     SUBROUTINE dgcxc( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss, gpu_args_ )
       USE kind_l, ONLY: DP
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length
       INTEGER,  INTENT(IN) :: sp
       REAL(DP), INTENT(IN) :: r_in(length,sp)
       REAL(DP), INTENT(IN) :: g_in(length,3,sp)
       REAL(DP), INTENT(OUT) :: dvxc_rr(length,sp,sp), dvxc_sr(length,sp,sp), &
                                dvxc_ss(length,sp,sp)
       LOGICAL, OPTIONAL, INTENT(IN) :: gpu_args_
     END SUBROUTINE
  END INTERFACE
  !
END MODULE xc_lib
!--------------------------------------------------
