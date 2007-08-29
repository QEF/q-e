!
! Copyright (C) 2001-2004 PWSCF group
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
  USE symme, ONLY : nsym, s
  USE cell_base, ONLY : bg

  implicit none
  
  integer :: i, iq, ierr, iudyn = 26

  real(DP), allocatable, dimension(:) :: wq  

  logical :: exist_gamma

  !
  !  calculates the Monkhorst-Pack grid
  !

  if( nq1 .le. 0 .or. nq2 .le. 0 .or. nq3 .le. 0 ) &
       call errore('q_points','nq1 or nq2 or nq3 .le. 0',1)
  if( iq1 .lt. 0 .or. iq2 .lt. 0 .or. iq3 .lt. 0 ) &
       call errore('q_points','iq1 or iq2 or iq3 .le. 0',1)
  if( iq1 .gt. nq1 .or. iq2 .gt. nq2 .or. iq3 .gt. iq3 ) &
       call errore('q_points','iq1 or iq2 or iq3 .le. 0',1)

  allocate (wq(nqmax))
  allocate (x_q(3,nqmax))
  call kpoint_grid( nsym, s, bg, nqmax, 2*iq1, 2*iq2, 2*iq3, &
                         nq1, nq2, nq3, nqs, x_q, wq )
  deallocate (wq)
  !
  ! Check if the Gamma point is one of the points and put
  ! it in the first position (it should already be the first)
  !
  exist_gamma = .false.
  do iq = 1, nqs
     if (abs(x_q(1,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(2,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(3,iq)) .lt. 1.0e-10_dp) then
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
  if (.not.(iq1.eq.0.and.iq2.eq.0.and.iq3.eq.0)) write(stdout, '(5x, &
       "with only (", 3(i2,","),") point requested")') iq1, iq2, iq3
  write(stdout, '(5x,"(",i4,"q-points):")') nqs
  write(stdout, '(5x,"  N       xq(1)       xq(2)       xq(3) " )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 3f12.5)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
  end do
  !
  if ((iq1.eq.0.and.iq2.eq.0.and.iq3.eq.0).and.&
       (.not. exist_gamma)) call errore('q_points','Gamma is not a q point',1)
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
