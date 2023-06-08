!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE deriv_drhoc( nt, ngl, gl, tpiba2, drhocg )
  !--------------------------------------------------------------------------
  !! Calculates the Fourier transform of \(d\text{Rho}_c/dG\).
  !
  USE kinds
  USE uspp_data, ONLY : tab_rhc, dq
  !
  IMPLICIT NONE
  !
  INTEGER :: nt
  !! input: atomic type
  INTEGER :: ngl
  !! input: the number of g shell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! input: the number of G shells
  REAL(DP), INTENT(IN) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP), INTENT(OUT) :: drhocg(ngl)
  !! Fourier transform of d Rho_c/dG
  !
  ! ... local variables
  !
    REAL(DP) :: gx, px, ux, vx, wx
  ! the modulus of g for a given shell
  ! variables used for interpolation
  INTEGER :: igl, i0, i1, i2, i3
  ! counters
  !
  !$acc data present_or_copyin(gl) present_or_copyout(drhocg)
  !$acc parallel loop
  DO igl = 1, ngl
     gx = SQRT( gl(igl) * tpiba2 )
     px = gx / dq - int (gx/dq)
     ux = 1.d0 - px
     vx = 2.d0 - px
     wx = 3.d0 - px
     i0 = INT(gx/dq) + 1
     i1 = i0 + 1
     i2 = i0 + 2
     i3 = i0 + 3
     drhocg (igl) = (- tab_rhc(i0, nt) * (ux*vx + vx*wx + ux*wx) / 6.0_dp &
                     + tab_rhc(i1, nt) * (wx*vx - px*wx - px*vx) / 2.0_dp &
                     - tab_rhc(i2, nt) * (wx*ux - px*wx - px*ux) / 2.0_dp &
                     + tab_rhc(i3, nt) * (ux*vx - px*ux - px*vx) / 6.0_dp ) / dq
  ENDDO
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE deriv_drhoc

