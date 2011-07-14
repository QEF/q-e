!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE ggen1d (ngm1d, g1d, gg1d, ig1dto3d, nl1d, igtongl1d)
  !
  !    This subroutine is the one dimensional equivalent of the ggen
  !    routine. It is used to produce the G vectors in the z directions
  !    to compute and output planar averages of the charge density
  !    of single levels. It should work in parallel machines.
  !
  USE kinds, ONLY: DP
  USE cell_base, ONLY : at
  USE gvect, ONLY: ngm, g, gg, igtongl
  USE fft_base, ONLY: dfftp
  IMPLICIT NONE
  INTEGER :: ngm1d, ig1dto3d (dfftp%nr3), igtongl1d (dfftp%nr3), nl1d (dfftp%nr3)
  ! output: the number of 1D G vectors on this processor
  ! output: correspondence 1D with 3D G vectors
  ! output: the correspondence with the shells
  ! output: correspondence 1D FFT mesh G with array

  real(DP) :: g1d (3, dfftp%nr3), gg1d (dfftp%nr3)
  ! output: ngm1d 3D vectors with the 1D G of this
  ! output: ngm1d scalars with the modulus of 1D G
  !
  !   local variables
  !

  INTEGER :: ig, ig1d
  ! counter on 3D vectors
  ! counter on 1D vectors

  real(DP), PARAMETER :: eps = 1.d-12


  g1d(:,:) = 0.d0
  gg1d(:) = 0.d0

  ig1d = 0
  DO ig = 1, ngm
     IF ( (abs(g(1,ig))<eps) .and. (abs(g(2,ig)) <eps) ) THEN
        !
        !    a vector of the 1D grid has been found
        !
        ig1d = ig1d+1
        IF (ig1d>dfftp%nr3) CALL errore ('ggen1d', 'too many G', 1)
        g1d (3, ig1d) = g (3, ig)
        gg1d (ig1d) = gg (ig)
        ig1dto3d (ig1d) = ig
        nl1d (ig1d) = nint (g (3, ig) * at (3, 3) ) + 1
        IF (nl1d (ig1d) <1) nl1d (ig1d) = nl1d (ig1d) + dfftp%nr3
     ENDIF
  ENDDO

  ngm1d = ig1d
  DO ig1d = 1, ngm1d
     igtongl1d (ig1d) = igtongl (ig1dto3d (ig1d) )
  ENDDO
  RETURN
END SUBROUTINE ggen1d
