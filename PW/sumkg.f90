!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function sumkg (et, nbnd, nks, wk, degauss, ngauss, e, is, isk)
  !-----------------------------------------------------------------------
  !
  !     This function computes the number of states under a given energy e
  !
  !
  USE kinds
  implicit none
  ! Output variable
  real(kind=DP) :: sumkg
  ! Input variables
  integer, intent(in) :: nks, nbnd, ngauss
  ! input: the total number of K points
  ! input: the number of bands
  ! input: the type of smearing
  real(kind=DP), intent(in) :: wk (nks), et (nbnd, nks), degauss, e
  ! input: the weight of the k points
  ! input: the energy eigenvalues
  ! input: gaussian broadening
  ! input: the energy to check
  integer, intent(in) :: is, isk(nks)
  !
  ! local variables
  !
  real(kind=DP), external :: wgauss
  ! function which compute the smearing
  real(kind=DP) ::sum1
  integer :: ik, ibnd
  ! counter on k points
  ! counter on the band energy
  !
  sumkg = 0.d0
  do ik = 1, nks
     sum1 = 0.d0
     if (is /= 0) then
        if (isk(ik).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        sum1 = sum1 + wgauss ( (e-et (ibnd, ik) ) / degauss, ngauss)
     enddo
     sumkg = sumkg + wk (ik) * sum1
  enddo
#ifdef __PARA
  call poolreduce (1, sumkg)
#endif
  return
end function sumkg

