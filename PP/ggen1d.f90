!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine ggen1d (ngm1d, g1d, gg1d, ig1dto3d, nl1d, igtongl1d)
  !
  !    This subroutine is the one dimensional equivalent of the ggen
  !    routine. It is used to produce the G vectors in the z directions
  !    to compute and output planar averages of the charge density
  !    of single levels. It should work in parallel machines.
  !

  use pwcom

  integer :: ngm1d, ig1dto3d (nr3), igtongl1d (nr3), nl1d (nr3)
  ! output: the number of 1D G vectors on this proce
  ! output: correspondence 1D with 3D G vectors
  ! output: the correspondence with the shells
  ! output: correspondence 1D FFT mesh G with array

  real(kind=DP) :: g1d (3, nr3), gg1d (nr3)
  ! output: ngm1d 3D vectors with the 1D G of this
  ! output: ngm1d scalars with the modulus of 1D G
  !
  !   local variables
  !

  integer :: ig, ig1d
  ! counter on 3D vectors
  ! counter on 1D vectors

  real(kind=DP) :: eps
  ! a small number
  parameter (eps = 1.d-12)


  call setv (nr3 * 3, 0.d0, g1d, 1)
  call setv (nr3, 0.d0, gg1d, 1)

  ig1d = 0
  do ig = 1, ngm
     if ( (abs (g (1, ig) ) .lt.eps) .and. (abs (g (2, ig) ) .lt.eps) ) &
          then
        !
        !    a vector of the 1D grid has been found
        !
        ig1d = ig1d+1
        if (ig1d.gt.nr3) call error ('ggen1d', 'too many G', 1)
        g1d (3, ig1d) = g (3, ig)
        gg1d (ig1d) = gg (ig)
        ig1dto3d (ig1d) = ig
        nl1d (ig1d) = nint (g (3, ig) * at (3, 3) ) + 1
        if (nl1d (ig1d) .lt.1) nl1d (ig1d) = nl1d (ig1d) + nr3
     endif
  enddo

  ngm1d = ig1d
  do ig1d = 1, ngm1d
     igtongl1d (ig1d) = igtongl (ig1dto3d (ig1d) )
  enddo
  return
end subroutine ggen1d
