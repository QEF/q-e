!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------


subroutine gweights (nks, wk, nbndx, nbnd, nelec, degauss, ngauss, &
     et, ef, demet, wg)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  use parameters
  implicit none  
  !
  integer :: nks, nbndx, nbnd, ngauss  
  real(kind=DP) :: wk (nks), et (nbndx, nks), nelec, degauss  

  real(kind=DP) :: wg (nbnd, nks), ef, demet  
  real(kind=DP) :: wgauss, w1gauss  
  integer :: kpoint, ibnd  



  external efermig, wgauss, w1gauss  
  ! Calculate the Fermi energy ef

  call efermig (et, nbndx, nbnd, nks, nelec, wk, degauss, ngauss, &
       ef)
  demet = 0.d0  
  do kpoint = 1, nks  
     do ibnd = 1, nbnd  
        ! Calculate the gaussian weights


        wg (ibnd, kpoint) = wk (kpoint) * wgauss ( (ef - et (ibnd, kpoint) &
             ) / degauss, ngauss)
        ! The correct form of the band energy is  \int e n(e) de   for e<ef
        ! demet is the correction to add to the sum of eigenvalues
        demet = demet + wk (kpoint) * degauss * w1gauss ( (ef - et (ibnd, &
             kpoint) ) / degauss, ngauss)
     enddo

  enddo
  return  
end subroutine gweights
