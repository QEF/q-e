!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine efermit (et, nbnd, nks, nelec, nspin, ntetra, tetra, ef, is, isk)
  !--------------------------------------------------------------------
  !
  !     Finds the Fermi energy - tetrahedron method (Bloechl)
  !
  USE io_global, ONLY : stdout
  USE kinds
  implicit none
  integer, intent(in)  :: nks, nbnd, nspin, ntetra, tetra (4, ntetra)
  ! input: the number of k points
  ! input: the number of bands
  ! input: the number of spin components
  ! input: the number of tetrahedra
  ! input: the vertices of a tetrahedron
  real(DP), intent(in) :: et (nbnd, nks), nelec
  ! input: the eigenvalues
  ! input: the number of electrons
  real(DP), intent(out) :: ef
  ! output: the fermi energy
  integer, intent(in) :: is, isk(nks)
  !
  !     two parameters
  !
  integer :: maxiter
  ! the maximum number of iterations in

  real(DP) :: rydtoev, eps
  ! the transformation Ry to eV
  ! a small quantity
  parameter (maxiter = 300, rydtoev = 13.6058d0, eps = 1.0d-10)
  !
  !     here the local variables
  !
  integer :: nlw, ik, iter
  ! the minimum energy band
  ! counter on k points
  ! counter on iterations

  real(DP) :: elw, eup, sumkup, sumklw, sumkt, sumkmid
  ! the lower limit of the fermi ener
  ! the upper limit of the fermi ener
  ! the number of states with eup
  ! the number of states with elw
  ! the sum over the states below E_f
  ! the number of states with ef

  real(DP) :: efbetter, better
  external sumkt
  !
  !      find bounds for the Fermi energy.
  !
  nlw = max (1, nint (nelec / 2.0 - 5.0) )
  elw = et (nlw, 1)
  eup = et (nbnd, 1)
  do ik = 2, nks
     elw = min (elw, et (nlw, ik) )
     eup = max (eup, et (nbnd, ik) )
  enddo
  !
  !      Bisection method
  !
  sumkup = sumkt (et, nbnd, nks, nspin, ntetra, tetra, eup, is, isk)
  sumklw = sumkt (et, nbnd, nks, nspin, ntetra, tetra, elw, is, isk)
  if ( (sumkup - nelec) .lt. - eps.or. (sumklw - nelec) .gt.eps) &
       call errore ('efermit', 'unexpected error', 1)
  better = 1.0e+10
  do iter = 1, maxiter
     ef = (eup + elw) / 2.0
     sumkmid = sumkt (et, nbnd, nks, nspin, ntetra, tetra, ef, is, isk)
     if (abs (sumkmid-nelec) .lt.better) then
        better = abs (sumkmid-nelec)
        efbetter = ef
     endif
     ! converged
     if (abs (sumkmid-nelec) .lt. eps) then
        goto 100
     elseif ( (sumkmid-nelec) .lt. -eps) then
        elw = ef
     else
        eup = ef
     endif

  enddo
  !     unconverged exit:
  !     the best available ef is used . Needed in some difficult cases
  ef = efbetter
  sumkmid = sumkt (et, nbnd, nks, nspin, ntetra, tetra, ef, is, isk )

  if (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, 9010) ef * rydtoev, sumkmid
  !     converged exit:
100 continue
  !     Check if Fermi level is above any of the highest eigenvalues
  do ik = 1, nks
     if (is /= 0) then
        if (isk(ik) .ne. is ) cycle
     end if
     if (ef.gt.et (nbnd, ik) + 1.d-4) WRITE( stdout, 9020) ef * rydtoev, ik, &
          et (nbnd, ik) * rydtoev
  enddo

  return
9010 format (/5x,'Warning: too many iterations in bisection'/ &
       &          5x,'ef = ',f10.6,' sumk = ',f10.6,' electrons')

9020 format (/5x,'Warning: ef =',f10.6, &
       &     ' is above the highest band at k-point',i4,/5x,9x, &
       &     'e  = ',f10.6)
end subroutine efermit
