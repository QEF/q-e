!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine iweights (nks, wk, nbnd, nelec, et, Ef, wg, is, isk)
  !--------------------------------------------------------------------
  !     calculates weights for semiconductors and insulators
  !     (bands are either empty or filled)
  !     On output, Ef is the highest occupied Kohn-Sham level
  USE kinds
  USE noncollin_module, ONLY: noncolin
  USE mp,       ONLY : mp_max
  USE mp_pools, ONLY : inter_pool_comm
  implicit none
  !
  integer, intent(in)  :: nks, nbnd, is, isk(nks)
  real(DP), intent(in) :: wk (nks), et(nbnd, nks), nelec
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  real(DP), intent(out) :: Ef
  integer :: kpoint, ibnd

  CALL iweights_only (nks, wk, is, isk, nbnd, nelec, wg )
  !
  Ef = - 1.0d+20
  do kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint) .ne.  is ) cycle
     end if
     do ibnd = 1, nbnd
        if (wg (ibnd, kpoint) > 0.d0 ) Ef = MAX (Ef, et (ibnd, kpoint) )
     enddo
  enddo
  !
  ! find max across pools
  !
  CALL mp_max( ef, inter_pool_comm )

  return
end subroutine iweights
!
!--------------------------------------------------------------------
subroutine iweights_only (nks, wk, is, isk, nbnd, nelec, wg )
  !--------------------------------------------------------------------
  !     calculates weights for semiconductors and insulators
  !     (bands are either empty or filled)

  USE kinds
  USE noncollin_module, ONLY: noncolin
  implicit none
  !
  integer, intent(in) :: nks, nbnd, is, isk(nks)
  real(DP), intent(in) :: wk (nks), nelec
  real(DP), intent(out) :: wg (nbnd, nks)
  real(DP) :: degspin 
  integer :: kpoint, ibnd

  degspin=2.d0
  if (noncolin) degspin=1.d0
  if (is /= 0)  degspin=1.d0
  do kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint) .ne.  is ) cycle
     end if
     do ibnd = 1, nbnd
        if (ibnd <= nint (nelec) / degspin) then
           wg (ibnd, kpoint) = wk (kpoint)
        else
           wg (ibnd, kpoint) = 0.d0
        endif
     enddo
  enddo

  return
end subroutine iweights_only
