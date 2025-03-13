!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE convert_kpts_names()
  !---------------------------------------------------------------------
  !
  !! This routine redefines the names of some of the variables
  !! read from the card K_POINTS in the kcw_ham program, to avoid
  !! conflicts between different variables, coming from klist and 
  !! input_parameters, having the same name.
  !
  USE input_parameters,     ONLY : xk, wk, nk1, nk2, nk3, nkstot, k_points
  USE control_kcw,          ONLY : xk_bands, wk_bands, nks_bands 
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k
  INTEGER :: kn
  !
  !
  IF ( k_points == 'automatic' ) THEN
    !
    nks_bands = nk1 * nk2 * nk3
    !
    ALLOCATE( xk_bands(3,nks_bands) )
    ALLOCATE( wk_bands(nks_bands) )
    !
    kn = 0
    !
    DO i = 1, nk1
      DO j = 1, nk2
        DO k = 1, nk3
          !
          kn = kn + 1
          !
          xk_bands(:,kn) = DBLE(i-1) / nk1 * (/1,0,0/) + &
                           DBLE(j-1) / nk2 * (/0,1,0/) + &
                           DBLE(k-1) / nk3 * (/0,0,1/)
          !
        ENDDO
      ENDDO
    ENDDO
    !
    wk_bands(:) = 1.0
    !
  ELSE
    !
    nks_bands = nkstot
    ALLOCATE( xk_bands(3,nks_bands) )
    ALLOCATE( wk_bands(nks_bands) )
    xk_bands = xk
    wk_bands = wk
    !
  ENDIF
  !
  !
END SUBROUTINE convert_kpts_names
