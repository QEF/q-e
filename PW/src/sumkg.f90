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
  USE mp_pools, ONLY : inter_pool_comm
  USE mp,       ONLY : mp_sum
  implicit none
  ! Output variable
  real(DP) :: sumkg
  ! Input variables
  integer, intent(in) :: nks, nbnd, ngauss
  ! input: the total number of K points
  ! input: the number of bands
  ! input: the type of smearing
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), degauss, e
  ! input: the weight of the k points
  ! input: the energy eigenvalues
  ! input: gaussian broadening
  ! input: the energy to check
  integer, intent(in) :: is, isk(nks)
  !
  ! local variables
  !
  real(DP), external :: wgauss
  ! function which compute the smearing
  real(DP) ::sum1
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
#if defined(__MPI)
  call mp_sum ( sumkg, inter_pool_comm )
#endif
  return
end function sumkg

