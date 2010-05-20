!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
FUNCTION efermit (et, nbnd, nks, nelec, nspin, ntetra, tetra, is, isk)
  !--------------------------------------------------------------------
  !
  !     Finds the Fermi energy - tetrahedron method
  !     (see P. E. Bloechl et al, PRB49, 16223 (1994))
  !
  USE io_global, ONLY : stdout
  USE kinds, ONLY: DP
  USE constants, ONLY: rytoev
  implicit none
  integer, intent(in)  :: nks, nbnd, nspin, ntetra, tetra (4, ntetra)
  ! nks   : the number of k points
  ! nbnd  : the number of bands
  ! nspin : the number of spin components
  ! ntetra: the number of tetrahedra
  ! tetra : the vertices of a tetrahedron
  real(DP), intent(in) :: et (nbnd, nks), nelec
  ! input: the eigenvalues
  ! input: the number of electrons
  real(DP):: efermit
  ! output: the fermi energy
  integer, intent(in) :: is, isk(nks)
  !
  !     two parameters
  !
  integer, parameter :: maxiter = 300
  ! the maximum number of iterations in bisection

  real(DP), parameter :: eps= 1.0d-10
  ! a small quantity
  !
  !     here the local variables
  !
  integer :: nlw, ik, iter
  ! the minimum energy band
  ! counter on k points
  ! counter on iterations

  real(DP) :: ef, elw, eup, sumkup, sumklw, sumkmid
  ! elw, eup: lower and upper bounds for fermi energy (ef)
  ! sumklw, sumkup: number of states for ef=elw, ef=eup resp.
  ! sumkmid:        number of states for ef=(elw+eup)/2
  real(DP), external :: sumkt

  real(DP) :: efbetter, better
  !
  !      find bounds for the Fermi energy.
  !
  nlw = max (1, nint (nelec / 2.0d0 - 5.0d0) )
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
  better = 1.0d+10
  if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps)  then
     !
     ! this is a serious error and the code should stop here
     ! we don't stop because this may occasionally happen in nonscf
     ! calculations where it may be completely irrelevant
     !
     call infomsg ('efermit', 'internal error, cannot bracket Ef')
     efermit = better
     return
  end if
  do iter = 1, maxiter
     ef = (eup + elw) / 2.d0
     sumkmid = sumkt (et, nbnd, nks, nspin, ntetra, tetra, ef, is, isk)
     if (abs (sumkmid-nelec) < better) then
        better = abs (sumkmid-nelec)
        efbetter = ef
     endif
     ! converged
     if (abs (sumkmid-nelec) < eps) then
        goto 100
     elseif ( (sumkmid-nelec) < -eps) then
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
  WRITE( stdout, 9010) ef * rytoev, sumkmid
  !     converged exit:
100 continue
  !     Check if Fermi level is above any of the highest eigenvalues
  do ik = 1, nks
     if (is /= 0) then
        if (isk(ik) /= is ) cycle
     end if
     if (ef > et (nbnd, ik) + 1.d-4) &
          WRITE( stdout, 9020) ef * rytoev, ik, et (nbnd, ik) * rytoev
  enddo

  efermit = ef
  return
9010 format (/5x,'Warning: too many iterations in bisection'/ &
       &          5x,'ef = ',f10.6,' sumk = ',f10.6,' electrons')

9020 format (/5x,'Warning: ef =',f10.6, &
       &     ' is above the highest band at k-point',i4,/5x,9x, &
       &     'e  = ',f10.6)
end FUNCTION efermit

