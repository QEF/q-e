!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
!  This subroutine allocates all of the fft stuff for the custom defined grid
! 
SUBROUTINE allocate_fft_custom(fc)
  
  USE kinds,              ONLY : DP
  USE gvect,              ONLY : g, mill
  USE cell_base,          ONLY : at, bg, tpiba2
  USE control_flags,      ONLY : gamma_only
  USE fft_custom,         ONLY : fft_cus, set_custom_grid, ggent
  USE grid_subroutines,   ONLY : realspace_grid_init_custom
  IMPLICIT NONE
  
  TYPE (fft_cus) :: fc
  
  INTEGER :: ng,n1t,n2t,n3t
  
  IF(fc%initalized) RETURN
  !
  fc%gcutmt = fc%dual_t*fc%ecutt / tpiba2
  !
  CALL realspace_grid_init_custom(fc%dfftt, at, bg, fc%gcutmt)
  !
  CALL data_structure_custom(fc, .TRUE.)
  !
  fc%initalized = .true.
  !
  CALL ggent(fc)
  
  RETURN
END SUBROUTINE allocate_fft_custom
