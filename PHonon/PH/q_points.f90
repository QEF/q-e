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
  USE io_global,  ONLY :  stdout, ionode, ionode_id
  USE disp,  ONLY : nq1, nq2, nq3, x_q, nqs, lgamma_iq
  USE output, ONLY : fildyn
  USE symm_base, ONLY : nsym, s, time_reversal, t_rev, invs
  USE cell_base, ONLY : at, bg
  USE control_ph, ONLY : search_sym
  USE mp_images,  ONLY : intra_image_comm
  USE mp,         ONLY : mp_bcast
  USE elph_tetra_mod, ONLY : lshift_q

  implicit none

  integer :: i, iq, ierr, iudyn = 26
  logical :: exist_gamma, check, skip_equivalence=.FALSE.
  logical, external :: check_q_points_sym
  real(DP), allocatable :: xq(:,:), wq(:)

  INTEGER :: nqmax
  !
  !  calculate the Monkhorst-Pack grid
  !

  if( nq1 <= 0 .or. nq2 <= 0 .or. nq3 <= 0 ) &
       call errore('q_points','nq1 or nq2 or nq3 <= 0',1)

  nqmax= nq1 * nq2 * nq3

  allocate (wq(nqmax))
  allocate (xq(3,nqmax))
  if(lshift_q) then
     call kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, bg, nqmax,&
     &                  1,1,1, nq1,nq2,nq3, nqs, xq, wq )
  else
     call kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, bg, nqmax,&
     &                  0,0,0, nq1,nq2,nq3, nqs, xq, wq )
  end if
  allocate(x_q(3,nqs))
  allocate(lgamma_iq(nqs))
  x_q(:,:)=xq(:,1:nqs)
  deallocate (xq)
  deallocate (wq)
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
  lgamma_iq=.FALSE.
  IF(.NOT. lshift_q) lgamma_iq(1)=.TRUE.
  !
  ! Write the q points in the output
  !
  write(stdout, '(//5x,"Dynamical matrices for (", 2(i2,","),i2,") &
           & uniform grid of q-points")') nq1, nq2, nq3
  IF (lshift_q) write(stdout,'(a)') "     With a half shift" 
  write(stdout, '(5x,"(",i4,"q-points):")') nqs
  write(stdout, '(5x,"  N         xq(1)         xq(2)         xq(3) " )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 3f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
  end do
  !
  IF ( (.NOT. exist_gamma) .AND. (.NOT. lshift_q) ) &
     CALL errore('q_points','Gamma is not a q point',1)
!
!  Check that the q point grid is compatible with the symmetry.
!  If this test is not passed, q2r will stop in any case.
!
  IF (lshift_q) THEN
     WRITE(stdout,'(a)') "     Because shifted q grid is used, q2r will not work !"
  ELSE
     IF (search_sym) THEN
        check=check_q_points_sym(nqs, x_q, at, bg, nsym, s, invs, nq1, nq2, nq3)
        IF (.NOT.check) THEN
           WRITE(stdout, '(/,5x,"This q-mesh breaks symmetry!")')
           WRITE(stdout, '(5x,"Try to choose different nq1, nq2, nq3")')
           WRITE(stdout, '(5x,"You can also continue by setting &
                       &search_sym=.false.")')
           WRITE(stdout, '(5x,"but be careful because q2r will not work")')
           CALL errore('q_points', 'q-mesh breaks symmetry', 1)
        ENDIF
     ENDIF
  END IF
  !
  ! ... write the information on the grid of q-points to file
  !
  IF (ionode) &
     OPEN (unit=iudyn, file=TRIM(fildyn)//'0', status='unknown', iostat=ierr)
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF ( ierr > 0 ) CALL errore ('q_points','cannot open file ' &
                   & // TRIM(fildyn) // '0', ierr)
  IF (ionode) THEN
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
