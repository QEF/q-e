!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function sumkg (et, nbndx, nbnd, nks, wk, degauss, ngauss, e)
!-----------------------------------------------------------------------
!
!     This function computes the number of states under a given energy e
!
!
use parameters
implicit none
!
!    First the dummy variables.
!
integer :: nks, nbndx, nbnd, ngauss
                         ! input: the total number of K points
                         ! input: the maximum number of bands
                         ! input: the number of bands
                         ! input: the type of smearing
real(kind=DP) :: wk (nks), et (nbndx, nks), degauss, sumkg, e
                         ! input: the weight of the k points
                          !input: the energy eigenvalues
                         ! input: the smearing function
                         ! output: the sum of the eigenvalues
                         ! input: the energy to check
!
!   the local variables
!
real(kind=DP) :: wgauss, sum1
                         ! function which compute the smearing
                         ! auxiliary function
integer :: ik, ibnd
                         ! counter on k points
                         ! counter on the band energy
!
sumkg = 0.d0
do ik = 1, nks
sum1 = 0.d0
do ibnd = 1, nbnd
sum1 = sum1 + wgauss ( (e-et (ibnd, ik) ) / degauss, ngauss)
enddo
sumkg = sumkg + wk (ik) * sum1
enddo
#ifdef PARA
call poolreduce (1, sumkg)
#endif
return
end function sumkg

