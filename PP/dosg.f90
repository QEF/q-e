!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine dos_g (et, nspin, nbnd, nks, wk, Degauss, ngauss, E, dosg)
  !--------------------------------------------------------------------
  !
  use parameters, only : DP
  implicit none
  integer :: nspin, nks, nbnd, ngauss

  real(kind=DP) :: wk (nks), et (nbnd, nks), Degauss, E, dosg (2)
  real(kind=DP) :: w0gauss
  integer :: n, ns, nk0, nk, ik
  external w0gauss
  !
  if (nspin.eq.1) then
     nk = nks
  else
     nk = nks / 2
  endif
  !
  do ns = 1, nspin
     if (ns.eq.1) then
        nk0 = 1
     else
        nk0 = nks / 2 + 1
     endif
     dosg (ns) = 0.0
     do ik = nk0, nk0 + nk-1
        do n = 1, nbnd
           dosg (ns) = dosg (ns) + wk (ik) * w0gauss ( (E-et (n, ik) ) &
                / Degauss, ngauss)
        enddo
     enddo
     !
     dosg (ns) = dosg (ns) / Degauss
     !
  enddo
  !
  return
end subroutine dos_g
