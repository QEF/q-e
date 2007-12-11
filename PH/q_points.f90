!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------
SUBROUTINE q_points ( )
!----------========------------------------------

  USE kinds, only : dp
  USE io_global,  ONLY :  stdout, ionode
  USE disp,  ONLY : nqmax, nq1, nq2, nq3, x_q, nqs
  USE disp,  ONLY : iq1, iq2, iq3
  USE output, ONLY : fildyn
  USE symme, ONLY : nsym, s, time_reversal, t_rev
  USE cell_base, ONLY : bg

  implicit none
  
  integer :: i, iq, ierr, iudyn = 26
  logical :: exist_gamma, single_q
  logical, external :: is_equivalent
  real(DP), allocatable, dimension(:) :: wq  

  !
  !  calculate the Monkhorst-Pack grid
  !

  if( nq1 <= 0 .or. nq2 <= 0 .or. nq3 <= 0 ) &
       call errore('q_points','nq1 or nq2 or nq3 <= 0',1)

  allocate (wq(nqmax))
  allocate (x_q(3,nqmax))
  call kpoint_grid( nsym, time_reversal, s, t_rev, bg, nqmax, &
                         0,0,0, nq1,nq2,nq3, nqs, x_q, wq )
  deallocate (wq)
  !
  !  if a single q-point of the grid requested
  !
  IF ( iq1 < 0 .or. iq2 < 0 .or. iq3 < 0 ) &
     CALL errore('q_points','iq1 or iq2 or iq3 < 0',1)
  IF ( iq1 > nq1 .or. iq2 > nq2 .or. iq3 > iq3 ) &
     CALL errore('q_points','iq1 or iq2 or iq3 > nq1 or nq2 or nq3',1)
  single_q = iq1 > 0 .AND. iq2 > 0 .AND. iq3 > 0
  IF ( single_q ) THEN
     DO iq = 1, nqs
        IF ( is_equivalent ( iq1, iq2, iq3, nq1, nq2, nq3, x_q(1,iq), bg, &
                             time_reversal, nsym, s, t_rev ) ) THEN
           x_q(:,1) = x_q(1,iq)
           nqs = 1
           GO TO 10
        END IF
     END DO
     CALL errore('q_points','could not find required q-point',1)
10   CONTINUE
  END IF
  !
  ! Check if the Gamma point is one of the points and put
  ! it in the first position (it should already be the first)
  !
  exist_gamma = .false.
  do iq = 1, nqs
     if ( abs(x_q(1,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(2,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(3,iq)) .lt. 1.0e-10_dp ) then
        exist_gamma = .true.
        if (iq .ne. 1) then
           do i = 1, 3
              x_q(i,iq) = x_q(i,1)
              x_q(i,1) = 0.0_dp 
           end do
        end if
     end if
  end do
  !
  ! Write the q points in the output
  !
  write(stdout, '(//5x,"Calculation of the dynamical matrices for (", & 
       &3(i2,","),") uniform grid of q-points")') nq1, nq2, nq3
  if ( single_q ) write(stdout, '(5x, "with only (", 3(i2,","), &
                                    & ") point requested")') iq1, iq2, iq3
  write(stdout, '(5x,"(",i4,"q-points):")') nqs
  write(stdout, '(5x,"  N       xq(1)       xq(2)       xq(3) " )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 3f12.5)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
  end do
  !
  IF ( .NOT. single_q .AND. .NOT. exist_gamma) &
     CALL errore('q_points','Gamma is not a q point',1)
  !
  ! ... write the information on the grid of q-points to file
  !
  IF (ionode) THEN
     OPEN (unit=iudyn, file=TRIM(fildyn)//'0', status='unknown', iostat=ierr)
     IF ( ierr > 0 ) CALL errore ('phonon','cannot open file ' &
          & // TRIM(fildyn) // '0', ierr)
     WRITE (iudyn, '(3i4)' ) nq1, nq2, nq3
     WRITE (iudyn, '( i4)' ) nqs
     DO  iq = 1, nqs
        WRITE (iudyn, '(3e24.15)') x_q(1,iq), x_q(2,iq), x_q(3,iq)
     END DO
     CLOSE (unit=iudyn)
  END IF
  return
end subroutine q_points
!
!-----------------------------------------------------------------------
LOGICAL FUNCTION is_equivalent ( ik1, ik2, ik3, nk1, nk2, nk3, xk, bg, &
                                 time_reversal, nsym, s, t_rev )
!-----------------------------------------------------------------------
  !
  ! ... Check if q-point defined by ik1, ik2, ik3 in the uniform grid
  ! ... is equivalent to xk (cartesian) - used for single-q calculation
  ! 
  USE kinds, ONLY: DP
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in) :: ik1,ik2,ik3, nk1,nk2,nk3, nsym, t_rev(48), s(3,3,48)
  LOGICAL,  INTENT(in) :: time_reversal 
  REAL(DP), INTENT(in) :: bg(3,3)
  REAL(DP), INTENT(in) :: xk(3)
  !
  REAL(DP), PARAMETER :: eps=1.0d-5
  REAL(DP) :: xk_(3), xkr(3)
  INTEGER :: ns
  !
  xk_(1) = DBLE(ik1-1)/nk1
  xk_(2) = DBLE(ik2-1)/nk2
  xk_(3) = DBLE(ik3-1)/nk3
  xk_(:) = xk_(:)-NINT(xk_(:))
  !
  DO ns=1,nsym
     !
     xkr(:) = s(:,1,ns) * xk_(1) &
            + s(:,2,ns) * xk_(2) &
            + s(:,3,ns) * xk_(3)
     xkr(:) = xkr(:) - NINT( xkr(:) )
     IF (t_rev(ns) == 1) xkr(:) = -xkr(:)
     ! 
     CALL cryst_to_cart (1, xkr, bg, 1)
     !
     is_equivalent = ABS( xkr(1)-xk(1) ) < eps .AND. &
                     ABS( xkr(2)-xk(2) ) < eps .AND. &
                     ABS( xkr(3)-xk(3) ) < eps
     IF ( is_equivalent)  RETURN 
     !
  END DO

  is_equivalent = .false.
  RETURN 

END FUNCTION is_equivalent
