!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
FUNCTION efermig (et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk)
  !--------------------------------------------------------------------
  !
  !     Finds the Fermi energy - Gaussian Broadening
  !     (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  USE mp,        ONLY : mp_max, mp_min
  USE mp_pools,  ONLY : inter_pool_comm
  implicit none
  !  I/O variables
  integer, intent(in) :: nks, nbnd, Ngauss, is, isk(nks)
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), Degauss, nelec
  real(DP) :: efermig
  !
  real(DP), parameter :: eps= 1.0d-10
  integer, parameter :: maxiter = 300
  ! internal variables
  real(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  real(DP), external::  sumkg
  integer :: i, kpoint
  !
  !      find bounds for the Fermi energy. Very safe choice!
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  do kpoint = 2, nks
     Elw = min (Elw, et (1, kpoint) )
     Eup = max (Eup, et (nbnd, kpoint) )
  enddo
  Eup = Eup + 2 * Degauss
  Elw = Elw - 2 * Degauss
  !
  ! find min and max across pools
  !
  call mp_max( eup, inter_pool_comm )
  call mp_min( elw, inter_pool_comm )
  !
  !      Bisection method
  !
  sumkup = sumkg (et, nbnd, nks, wk, Degauss, Ngauss, Eup, is, isk)
  sumklw = sumkg (et, nbnd, nks, wk, Degauss, Ngauss, Elw, is, isk)
  if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
       call errore ('efermig', 'internal error, cannot bracket Ef', 1)
  do i = 1, maxiter
     Ef = (Eup + Elw) / 2.d0
     sumkmid = sumkg (et, nbnd, nks, wk, Degauss, Ngauss, Ef, is, isk)
     if (abs (sumkmid-nelec) < eps) then
        efermig = Ef
        return
     elseif ( (sumkmid-nelec) < -eps) then
        Elw = Ef
     else
        Eup = Ef
     endif
  enddo
  if (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
       Ef * rytoev, sumkmid
  !
  efermig = Ef
  return
end FUNCTION efermig

