!------------------------------------------------
subroutine q_points
!----------========------------------------------

  use kinds, only : dp
  USE io_global,  ONLY :  stdout
  use disp, only : nqmax, nq1, nq2, nq3, x_q, nqs
  USE symme, only : nsym, s
  USE cell_base, ONLY : bg

  implicit none
  
  integer :: i, iq

  real(kind = dp), allocatable, dimension(:) :: wq  

  logical :: exist_gamma

  !
  !  calculates the Monkhorst-Pack grid
  !

  if( nq1 .le. 0 .or. nq2 .le. 0 .or. nq3 .le. 0 ) &
       call errore('q_points','nq1 or nq2 or nq3 .le. 0',1)

  allocate (wq(nqmax))
  allocate (x_q(3,nqmax))
  call kpoint_grid( nsym, s, bg, nqmax, 0, 0, 0, &
                         nq1, nq2, nq3, nqs, x_q, wq )
  deallocate (wq)
  !
  ! Check if the Gamma point is one of the points and put
  ! it in the first position
  ! (It should already be the first)
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
  write(stdout, '(//5x,"Calculation of the dynamical matrices for the following points:")')
  write(stdout, '(5x,"Number of q points:",i4)') nqs
  write(stdout, '(5x," Nr:      xq(1)       xq(2)       xq(3) " )')
  do iq = 1, nqs
     write(stdout, '(5x,i3, 3f12.5)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq)
  end do

  if(.not. exist_gamma) call errore('q_points','Gamma is not a q point',1)

  return
end subroutine q_points
