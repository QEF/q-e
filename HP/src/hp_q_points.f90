!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------
SUBROUTINE hp_q_points ( )
!------------------------------------------------
  !
  ! Generate the q points
  !
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE symm_base,    ONLY : nsym, s, time_reversal, t_rev, invs
  USE cell_base,    ONLY : at, bg
  USE ldaU_hp,      ONLY : skip_equivalence_q, nq1, nq2, nq3, &
                           x_q, nqs, search_sym

  implicit none

  integer :: i, iq, ierr
  logical :: exist_gamma, check
  logical, external :: check_q_points_sym
  real(DP), allocatable :: xq(:,:), wq(:)
  INTEGER :: nqmax
  !
  if ( nq1 <= 0 .or. nq2 <= 0 .or. nq3 <= 0 ) &
       call errore('hp_q_points','nq1 or nq2 or nq3 <= 0',1)
  !
  nqmax = nq1 * nq2 * nq3
  !
  allocate (wq(nqmax))
  allocate (xq(3,nqmax))
  !
  CALL kpoint_grid( nsym, time_reversal, skip_equivalence_q, s, t_rev, bg, nqmax,&
                         0,0,0, nq1,nq2,nq3, nqs, xq, wq )
  !
  allocate(x_q(3,nqs))
  x_q(:,:) = xq(:,1:nqs)
  deallocate (xq)
  !
  ! Check if the Gamma point is one of the points and put
  ! it in the first position (it should already be the first)
  !
  exist_gamma = .false.
  !
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
  write(stdout, '(//5x,"The grid of q-points (", 2(i2,","),i2,")",2x, &
                          & "(",i3," q-points ) :")') nq1, nq2, nq3, nqs
  write(stdout, '(5x,"  N       xq(1)         xq(2)         xq(3)       wq" )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 4f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq), wq(iq)
  enddo
  !
  IF (.NOT. exist_gamma) CALL errore('hp_q_points','Gamma is not a q point', 1)
  !
  !  Check that the q point grid is compatible with the symmetry.
  !
  IF (search_sym) THEN
     !
     check = check_q_points_sym(nqs, x_q, at, bg, nsym, s, invs, nq1, nq2, nq3)
     !
     IF (.NOT.check) THEN
        WRITE(stdout, '(/,5x,"This q-mesh breaks symmetry!")')
        WRITE(stdout, '(/,5x,"Try to disable the symmetry (nosym=.true. and noinv=.true. in PWscf)!")')
        WRITE(stdout, '(5x,"Or try to choose different nq1, nq2, nq3")')
        CALL errore('hp_q_points', 'q-mesh breaks symmetry', 1)
     ENDIF
     !
  ENDIF
  !
  deallocate(wq)
  !
  return
  !
END SUBROUTINE hp_q_points
!
