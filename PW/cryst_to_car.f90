!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cryst_to_cart (nvec, vec, trmat, iflag)
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to cartesian coordinates ( iflag=
  !     and viceversa ( iflag=-1 ).
  !     Output cartesian coordinates are stored in the input ('vec') array
  !
  !
  use parameters
  implicit none
  !
  !     first the dummy variables
  !
  integer :: nvec, iflag
  ! input: number of vectors (atom. pos. or k-p
  !        to be transf. from cryst. to cart. a
  ! input: gives the sense of the transformatio
  real(kind=DP) :: vec (3, nvec), trmat (3, 3)
  ! input/output: cryst./cart. coord. of the ve
  !               (atom. pos. or k-points)
  ! input: transformation matrix
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space la
  !                  for atoms   or
  !          = bg ,  basis of the rec.-space la
  !                  for k-points
  ! if iflag=-1: the opposite
  !
  !    here the local variables
  !

  integer :: nv, kpol
  ! counter on vectors
  ! counter on polarizations
  real(kind=DP) :: vau (3)
  ! auxil. vector (containing the temp. transf.
  !
  !     Compute the cartesian coordinates of each vectors
  !     (atomic positions or k-points components)
  !
  do nv = 1, nvec
     if (iflag.eq.1) then
        do kpol = 1, 3
           vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
                * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
        enddo
     else
        do kpol = 1, 3
           vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
                * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
        enddo
     endif
     do kpol = 1, 3
        vec (kpol, nv) = vau (kpol)
     enddo
  enddo
  !
  return
end subroutine cryst_to_cart
