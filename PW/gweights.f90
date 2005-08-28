!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine gweights (nks, wk, nbnd, nelec, degauss, ngauss, &
     et, ef, demet, wg, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  USE kinds
  implicit none
  !
  integer, intent(in) :: nks, nbnd, ngauss
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss
  real(DP), intent(out) :: wg (nbnd, nks), ef, demet
  integer, intent(in) :: is, isk(nks)
  !
  integer :: kpoint, ibnd
  real(DP) , external :: wgauss, w1gauss

  ! Calculate the Fermi energy ef


  call efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, ef, is, isk)
  demet = 0.d0
  do kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        ! Calculate the gaussian weights
        wg (ibnd, kpoint) = wk (kpoint) * &
                            wgauss ( (ef-et(ibnd,kpoint)) / degauss, ngauss)
        ! The correct form of the band energy is  \int e n(e) de   for e<ef
        ! demet is the correction to add to the sum of eigenvalues
        demet = demet + wk (kpoint) * &
                 degauss * w1gauss ( (ef-et(ibnd,kpoint)) / degauss, ngauss)
     enddo

  enddo
  return
end subroutine gweights
