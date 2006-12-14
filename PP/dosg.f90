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
  USE kinds, only : DP
  implicit none
  integer :: nspin, nks, nbnd, ngauss

  real(DP) :: wk (nks), et (nbnd, nks), Degauss, E, dosg (2)
  real(DP) :: w0gauss
  integer :: n, ns, nk0, nk, ik
  integer :: nspin0
  external w0gauss
  !
  if (nspin == 1 .or. nspin == 4) then
     nk = nks
  else
     nk = nks / 2
  endif
  nspin0=nspin
  if (nspin==4) nspin0=1
  !
  do ns = 1, nspin0
     if (ns.eq.1) then
        nk0 = 1
     else
        nk0 = nks / 2 + 1
     endif
     dosg (ns) = 0.0d0
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
