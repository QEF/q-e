!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine iweights (nks, wk, nbnd, nelec, wg)
  !--------------------------------------------------------------------
  !     calculates weights for semiconductors and insulators
  !     (bands are either empty or filled)
  use parameters
  implicit none
  !
  integer :: nks, nbnd
  real(kind=DP) :: wk (nks), nelec
  real(kind=DP) :: wg (nbnd, nks)
  real(kind=DP), parameter :: degspin = 2.d0
  integer :: kpoint, ibnd

  do kpoint = 1, nks
     do ibnd = 1, nbnd
        if (ibnd.le.nint (nelec) / degspin) then
           wg (ibnd, kpoint) = wk (kpoint)
        else
           wg (ibnd, kpoint) = 0.d0
        endif
     enddo

  enddo
  return
end subroutine iweights
