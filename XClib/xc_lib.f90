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
  USE dft_mod
  !
  IMPLICIT NONE
  !
  SAVE 
  ! 
  PRIVATE
  !
  !
  !PUBLIC :: xc, dmxc         !LDA
  !
  PUBLIC :: xc_gcx !, dgcxc   !GGA
  !
  !PUBLIC :: xc_metagcx       !MGGA
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
  INTERFACE xc_gcx
     SUBROUTINE xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
       USE kind_l,        ONLY: DP  
       IMPLICIT NONE
       INTEGER,  INTENT(IN) :: length, ns
       REAL(DP), INTENT(IN) :: rho(:,:), grho(:,:,:)
       REAL(DP), INTENT(OUT) :: ex(:), ec(:)
       REAL(DP), INTENT(OUT) :: v1x(:,:), v2x(:,:)
       REAL(DP), INTENT(OUT) :: v1c(:,:), v2c(:,:)
       REAL(DP), OPTIONAL, INTENT(OUT) :: v2c_ud(:)
     END SUBROUTINE
  END INTERFACE
  !
  !
END MODULE xc_lib
!--------------------------------------------------
