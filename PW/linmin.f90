!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine linmin (xold, eold, deold, xtot, etot, detot, xnew, &
     minimum_ok)
  !-----------------------------------------------------------------------
  !
  ! as an approximation to the line minimization find the minimum of a
  ! cubic interpolation of the available data at the end points:
  !
  !    e(x) = a + b x + c xx + d xxx ; de/dx = b + 2c x + 3d xx
  !
  ! If a reliable minimum cannot be found make a reasonable guess for
  ! an additional step and set minimum_ok=.false.
  !
  ! on input: eold  = e (xold)
  !           etot  = e (xtot)
  !           deold = de/dx(xold)
  !           detot = de/dx(xtot)
  !           dtau  = tau(xtot)-tau(xold)
  ! deold must be negative (i.e. the searched direction must be down-hill)
  !
  ! several cases apply:
  ! (detot > 0) case: a nice minimum should exist with xold < xnew < xtot
  ! (detot < 0) cases:
  !   ...     (d > 0) case: a nice minimum should exist with  xnew > xtot
  !   ...     (d < 0) cases:
  !   ...       ...   dbc2 > 1: No local minimum exists
  !                      c < 0: No local minimum exists for xnew > xold
  !                             in both cases better take an additional st
  !                             Next step size is changed according to
  !                             force ratio if not too large
  !   ...       ...  (dbc2 < 1) AND (c > 0): a minimum exists for positive
  !                             if (xnew>=xtot): let's hope it is a good o
  !                             if (xnew< xtot): something strange is happ
  !                                        let's go there and check
  !
  ! In any case, if the new estimated position is too far, better take an
  ! additional reduced step and see what happens
  !
  USE io_global,  ONLY : stdout
  USE kinds
  implicit none
  !  Input
  real(DP) :: xold, eold, deold, xtot, etot, detot
  !  Output
  logical :: minimum_ok
  real(DP) :: xnew
  !  Local
  ! maximum allowed step size
  real(DP) :: xmax
  parameter (xmax = 10.d0)


  real(DP) :: b, c, c2, d, dbc2, dx, x, enew
  WRITE( stdout, 100) eold, etot, deold, detot
  if (deold.gt.0) call errore ('linmin', 'search direction is up-hill &
       &', 1)

100 format (5x,'searching for next position (pslinmin)...', &
       &          /5x,'Eold  = ',f15.8, &
       &          /5x,'Etot  = ',f15.8, &
       &          /5x,'DEold = ',f15.8 &
       &          /5x,'DEtot = ',f15.8)
  minimum_ok = .true.
  dx = xtot - xold
  !
  ! the variable x used in the following is x = (xnew-xold)/(xtot-xold)
  !
  b = deold * dx
  c = 3.d0 * etot - detot * dx - 3.d0 * eold-2.d0 * deold * dx
  d = - 2.d0 * etot + detot * dx + 2.d0 * eold+deold * dx
  dbc2 = 3.d0 * d * b / c**2
  c2 = (detot - deold) / 2.d0 * dx
  if (dbc2.gt.1.d0) then
     !
     ! NB: since deold.lt.0, dbc2.gt.1.d0 may occour only if detot.lt.0 !!
     !
     WRITE( stdout, '(5x,"linmin: no 3rd order solution")')
     x = - sign (999.d0, d)
     !
     ! 2nd order solution: x = -b/(2.d0*c2)
     !
  elseif (abs (dbc2) .lt.1.d-2) then
     WRITE( stdout, 110)
110  format    (5x,'linmin: 2nd order interpolation', &
          &              ' plus 3rd order corrections')
     x = - b / (2.d0 * c) * (1.d0 + dbc2 / 2.d0)
  else
     !
     ! x is the solution corresponding to the minimum
     !
     x = c * ( - 1.d0 + sign (1.d0, c) * sqrt (1.d0 - dbc2) ) &
          / 3.d0 / d
     WRITE( stdout, '(5x,"linmin: 3rd order interpolation")')


  endif
  !c      WRITE( stdout,'(5x,"b, c, d, dbc2 =",4f12.6)') b, c, d, dbc2
  !c      WRITE( stdout,'(5x,"x ",f12.6)') x
  if (detot.gt.0.d0) then
     !
     ! (detot > 0) case: a nice minimum should exist with 0 < x < 1
     !
     if (x.lt.0.d0.or.x.gt.1.d0) call errore ('linmin', 'unexpected error', 1)
  else
     !
     ! (detot < 0) cases:
     !
     if (d.gt.0.d0) then
        !
        !   ...     (d > 0) case: a nice minimum should exist with x>1
        !
        if (x.lt.1.d0) call errore ('linmin', 'unexpected error', 2)
     else
        !
        !   ...     (d < 0) cases:
        !
        if (dbc2.gt.1.d0.or.c.lt.0.d0) then
           !
           !   ...       ...   dbc2 > 1: No local minimum exists
           !                      c < 0: No local minimum exists for x > 0
           !                             in both cases better take an additional st
           !                             Next step size (x-1) is changed according
           !                             force ratio if not too large
           !
           minimum_ok = .false.
           x = 1.d0 + min (2.d0, abs (detot / deold) )
           WRITE( stdout, '(5x,"linmin: no reliable minimum found")')
        else
           !
           !   ...       ...  (dbc2 < 1) AND (c > 0): a minimum exists for positive
           !                             if (x>=1): let's hope it is a good one.
           !                             if (x< 1): something strange is happening:
           !                                        let's go there and check
           !
           if (x.lt.0) call errore ('linmin', 'unexpected error', 3)
           if (x.lt.1) then
              minimum_ok = .false.
              WRITE( stdout, '(5x,"linmin: new pos. on the wrong side")')
           endif
        endif
     endif
  endif
  !
  ! In any case, if the new estimated position is too far, better take an
  ! additional reduced step and see what happens
  !
  if (x.gt.xmax) then
     x = xmax
     minimum_ok = .false.

  endif
  if (minimum_ok) then
     enew = eold+b * x + c * x**2 + d * x**3
     WRITE( stdout, '(/5x,"Enext = ",f15.8," Xnext=",f12.6)') enew, &
          x
  else
     WRITE( stdout, '(5x,"linmin: take another downhill step")')

  endif

  xnew = xold+x * dx
  return
end subroutine linmin
