!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------
SUBROUTINE q_points_wannier ( )
!----------========------------------------------

  USE kinds,          ONLY : dp
  USE io_global,      ONLY : stdout, ionode, ionode_id
  USE io_files,       ONLY : prefix
  USE mp,             ONLY : mp_bcast
  USE mp_world,       ONLY : world_comm
  USE disp,           ONLY : nq1, nq2, nq3, x_q, nqs, lgamma_iq
  USE output,         ONLY : fildyn
  USE el_phon,        ONLY : wan_index_dyn
  USE dfile_autoname, ONLY : dfile_get_qlist
  USE dfile_star,     ONLY : dvscf_star
  USE io_files,       ONLY : prefix
  USE control_ph,     ONLY : last_q
  USE mp_images,      ONLY : intra_image_comm

  implicit none

  integer :: i, iq, ierr, iudyn = 26
  logical :: exist_gamma
  logical :: exst
  real(DP), allocatable :: wq(:)

  !
  !  calculate the Monkhorst-Pack grid
  !

  if( nq1 <= 0 .or. nq2 <= 0 .or. nq3 <= 0 ) &
       call errore('q_points_wannier','nq1 or nq2 or nq3 <= 0',1)

  nqs=nq1*nq2*nq3
 
  if(last_q.lt.nqs.and.last_q.gt.0) nqs=last_q
  allocate (lgamma_iq(nqs))
  allocate (x_q(3,nqs))
  allocate(wan_index_dyn(nqs))

! here read q_points
  CALL dfile_get_qlist(x_q, nqs, dvscf_star%ext, TRIM(dvscf_star%dir)//prefix, wan_index_dyn )

  call mp_bcast(x_q,ionode_id, world_comm)
  call mp_bcast(wan_index_dyn, ionode_id, world_comm)
  
  !
  ! Check if the Gamma point is one of the points and put
  ! 
  exist_gamma = .false.
  do iq = 1, nqs
     if ( abs(x_q(1,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(2,iq)) .lt. 1.0e-10_dp .and. &
          abs(x_q(3,iq)) .lt. 1.0e-10_dp ) then
        exist_gamma = .true.
        if (iq .ne. 1) then
           call errore('q_points_wannier','first q in dirfile must be Gamma',1)
        end if
     end if
  end do
  lgamma_iq=.FALSE.
  lgamma_iq(1)=.TRUE.
  !
  ! Write the q points in the output
  !
  write(stdout, '(//5x,"Dynamical matrices for (", 3(i2,","),") &
           & uniform grid of q-points")') nq1, nq2, nq3
  write(stdout, '(5x,"(",i4,"q-points):")') nqs
  write(stdout, '(5x,"  N         xq(1)         xq(2)         xq(3) " )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 3f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
  end do
  !
  IF ( .NOT. exist_gamma) &
     CALL errore('q_points_wannier','Gamma is not a q point',1)

  !
  ! ... write the information on the grid of q-points to file
  !
  IF (ionode) &
     OPEN (unit=iudyn, file=TRIM(fildyn)//'0_qstar', status='unknown', iostat=ierr)
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF ( ierr > 0 ) CALL errore ('q_point_wannier','cannot open file ' &
          & // TRIM(fildyn) // '0_qstar', ierr)
  IF (ionode) THEN
     WRITE (iudyn, '(3i4)' ) nq1, nq2, nq3
     WRITE (iudyn, '( i4)' ) nqs
     DO  iq = 1, nqs
        WRITE (iudyn, '(3e24.15)') x_q(1,iq), x_q(2,iq), x_q(3,iq)
     END DO
     CLOSE (unit=iudyn)
  END IF
  return
end subroutine q_points_wannier
!
