!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine efermit (et, nbndx, nbnd, nks, nelec, nspin, ntetra, &
     tetra, ef)
  !--------------------------------------------------------------------
  !
  !     Finds the Fermi energy - tetrahedron method (Bloechl)
  !
  use parameters
  implicit none
  integer :: nks, nbndx, nbnd, nspin, ntetra, tetra (4, ntetra)
  ! input: the number of k points
  ! input: the maximum number of bands
  ! input: the number of bands
  ! input: the number of spin components
  ! input: the number of tetrahedra
  ! input: the vertices of a tetrahedron
  real(kind=DP) :: et (nbndx, nks), nelec, ef
  ! input: the eigenvalues
  ! input: the number of electrons
  ! output: the fermi energy
  !
  !     two parameters
  !
  integer :: maxiter
  ! the maximum number of iterations in

  real(kind=DP) :: rydtoev, eps
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

  real(kind=DP) :: elw, eup, sumkup, sumklw, sumkt, sumkmid
  ! the lower limit of the fermi ener
  ! the upper limit of the fermi ener
  ! the number of states with eup
  ! the number of states with elw
  ! the sum over the states below E_f
  ! the number of states with ef

  real(kind=DP) :: efbetter, better
  external sumkt, error
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
  sumkup = sumkt (et, nbndx, nbnd, nks, nspin, ntetra, tetra, eup)
  sumklw = sumkt (et, nbndx, nbnd, nks, nspin, ntetra, tetra, elw)
  if ( (sumkup - nelec) .lt. - eps.or. (sumklw - nelec) .gt.eps) &
       call error ('efermit', 'unexpected error', 1)
  better = 1.0e+10
  do iter = 1, maxiter
     ef = (eup + elw) / 2.0
     sumkmid = sumkt (et, nbndx, nbnd, nks, nspin, ntetra, tetra, ef)
     if (abs (sumkmid-nelec) .lt.better) then
        better = abs (sumkmid-nelec)
        efbetter = ef
     endif
     ! converged
     if (abs (sumkmid-nelec) .lt.eps) then
        goto 100
     elseif ( (sumkmid-nelec) .lt. - eps) then
        elw = ef
     else
        eup = ef
     endif


  enddo
  !     unconverged exit:
  !     the best available ef is used . Needed in some difficult cases
  ef = efbetter
  sumkmid = sumkt (et, nbndx, nbnd, nks, nspin, ntetra, tetra, ef)


  write (6, 9010) ef * rydtoev, sumkmid
  !     converged exit:
100 continue
  !     Check if Fermi level is above any of the highest eigenvalues
  do ik = 1, nks
     if (ef.gt.et (nbnd, ik) + 1.d-4) write (6, 9020) ef * rydtoev, ik, &
          et (nbnd, ik) * rydtoev

  enddo

  return
9010 format (/5x,'Warning: too many iterations in bisection'/ &
       &          5x,'ef = ',f10.6,' sumk = ',f10.6,' electrons')

9020 format (/5x,'Warning: ef =',f10.6, &
       &     ' is above the highest band at k-point',i4,/5x,9x, &
       &     'e  = ',f10.6)
end subroutine efermit
