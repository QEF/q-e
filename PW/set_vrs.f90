!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)
  !--------------------------------------------------------------------
  ! set the total local potential vrs on the smooth mesh to be used in h_p
  ! adding the (spin dependent) scf (H+xc) part and the sum of all the loc
  ! pseudopotential contributions.
  !
  USE kinds
  implicit none

  integer :: nspin, nrxx, ir, is
  ! input: number of spin components: 1 if lda, 2 if lsd
  ! input: the fft grid dimension
  ! counter for fft grid
  ! counter for spin polarizations

  real(kind=DP) :: vrs (nrxx, nspin), vltot (nrxx), vr (nrxx, nspin)
  ! output: total local potential on the smooth grid
  !         vrs=vltot+vr
  ! input: the total local pseudopotential
  ! input: the scf(H+xc) part of the local po
  ! input: true if a doublegrid is used

  logical :: doublegrid
  do is = 1, nspin
     !
     ! define the total local potential (external + scf) for each spin ...
     !
     do ir = 1, nrxx
        vrs (ir, is) = vltot (ir) + vr (ir, is)
     enddo
     !
     ! ... and interpolate it on the smooth mesh if necessary
     !
     if (doublegrid) call interpolate (vrs (1, is), vrs (1, is), - 1)
  enddo
  return

end subroutine set_vrs
