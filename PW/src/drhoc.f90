!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE drhoc( nt, ngl, gl, tpiba2, rhocg )
  !-----------------------------------------------------------------------
  !! Calculates the radial Fourier transform of the core charge.
  !
  USE kinds,     ONLY : dp
  USE uspp_data, ONLY : tab_rhc, dq
  !
  IMPLICIT NONE
  !
  INTEGER :: nt
  !! input: atomic type
  INTEGER :: ngl
  !! input: the number of g shell
  REAL(DP) :: gl(ngl)
  !! input: the number of G shells
  REAL(DP) :: tpiba2
  !! input: 2 times pi / alat
  REAL(DP) :: rhocg(ngl)
  !! output: the Fourier transform of the core charge
  !
  ! ... local variables
  !
  REAL(DP) :: gx, px, ux, vx, wx
  ! the modulus of g for a given shell
  ! variables used for interpolation
  INTEGER :: igl, i0, i1, i2, i3
  ! counters
  !
  !$acc data present_or_copyin(gl) present_or_copyout(rhocg) present(tab_rhc)
  !$acc parallel loop
  DO igl = 1, ngl
     gx = SQRT(gl(igl) * tpiba2)
     px = gx / dq - int (gx/dq)
     ux = 1.d0 - px
     vx = 2.d0 - px
     wx = 3.d0 - px
     i0 = INT(gx/dq) + 1
     i1 = i0 + 1
     i2 = i0 + 2
     i3 = i0 + 3
     rhocg (igl) = tab_rhc(i0, nt) * ux * vx * wx / 6.d0 + &
                   tab_rhc(i1, nt) * px * vx * wx / 2.d0 - &
                   tab_rhc(i2, nt) * px * ux * wx / 2.d0 + &
                   tab_rhc(i3, nt) * px * ux * vx / 6.d0

  ENDDO
  !$acc end data
  !
END SUBROUTINE drhoc
