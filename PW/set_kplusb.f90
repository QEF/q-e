!-----------------------------------------------------------------------
subroutine set_kplusb (ibrav, xk, wk, b_length, nks, npk, lcart)
  !-----------------------------------------------------------------------
  !    This subroutine sets the k and k+b points (with zero weights) 
  !    used in the preparatory run for a raman calculation
  !
  !    on input: xk and wk contain k-points and corresponding weights
  !
  !    on output: the number of points is enlarged, the first are the 
  !               original ones
  !  
  use kinds, only : DP
  implicit none
  
  integer :: npk, nks, ibrav
  
  real (kind = dp), intent(in) :: b_length
  ! input: length of the b_vector

  real(kind=DP), intent(inout):: xk (3, npk), wk (npk)
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  
  logical, intent(in) :: lcart
  ! input: if .true. construction of additional k-points in cart. directions
  !        else in directions of the nearest neighbor k-points

  integer :: i, ik
  real(kind = dp), allocatable, dimension(:,:) :: b 

  if (lcart) then
     if (7 * nks.gt.npk) call errore ('set_kplusb', 'too many k points', &
          & nks)
     allocate( b(3,6))
     b = 0.d0
     do i = 1, 3
        b(i,i) = 1.0_dp
        b(i,i+3) = -1.0_dp
     end do
     b = b * b_length

     ! Add k+b to the k-points 

     do ik = nks, 1, -1
        xk(:,( 7*(ik-1) )+1) = xk(:,ik)
        wk( ( 7*(ik-1) )+1 ) = wk(ik)
     end do
     do ik = 1, nks
        do i = 1, 6
           xk(:,7*(ik-1)+i+1) = xk(:,7*(ik-1) +1) + b(:,i)
           wk( 7 *(ik-1)+i+1 ) = 0.0_dp
        end do
     end do
     nks = 7 * nks
     deallocate(b)
  else
     select case (ibrav)
     case(1)
        if (7 * nks.gt.npk) call errore ('set_kplusb', 'too many k points', &
             & nks)
        allocate( b(3,6))
        b = 0.d0
        do i = 1, 3
           b(i,i) = 1.0_dp
           b(i,i+3) = -1.0_dp
        end do
        b = b * b_length
        ! Add k+b to the k-points 
        do ik = nks, 1, -1
           xk(:,( 7*(ik-1) )+1) = xk(:,ik)
           wk( ( 7*(ik-1) )+1 ) = wk(ik)
        end do
        do ik = 1, nks
           do i = 1, 6
              xk(:,7*(ik-1)+i+1) = xk(:,7*(ik-1) +1) + b(:,i)
              wk( 7 *(ik-1)+i+1 ) = 0.0_dp
           end do
        end do
        nks = 7 * nks
        deallocate(b)
    case (2)
        if ( 9 * nks.gt.npk) call errore ('set_kplusb', 'too many k points', &
             & nks)
        ! set the vector b for this case
        allocate(b(3,8))
        b = 1.0_dp
        b(:,2) = - b(:,2)
        b(1,7) = - b(1,7)
        b(:,8) = - b(:,7)
        b(1,5) = - b(1,5)
        b(2,5) = - b(2,5)
        b(3,6) = - b(3,6)
        b(2,3) = - b(2,3)
        b(:,4) = - b(:,3)
        b = b * b_length 
        ! Add k+b to the k-points

        do ik = nks, 1, -1
           xk(:,( 9*(ik-1) )+1) = xk(:,ik)
           wk( ( 9*(ik-1) )+1 ) = wk(ik)
        end do
        do ik = 1, nks
           do i = 1, 8
              xk(:,9*(ik-1)+i+1) = xk(:,( 9*(ik-1) )+1 ) + b(:,i)
              wk( 9 *(ik-1)+i+1 ) = 0.0_dp
           end do
        end do
        nks = 9 * nks
        deallocate(b)
     case(3)
        if (13 * nks.gt.npk) call errore ('set_kplusb', 'too many k points', &
             & nks)
        allocate(b(3,12))
        b = 0.0_dp
        b(1:2,1:4) = 1.0_dp
        b(1,2) = -b(1,1)
        b(2,3) = -b(2,1)
        b(1,4) = -b(1,1)
        b(2,4) = -b(2,1)
        b(1,5:8) = 1.0_dp
        b(3,5:8) = 1.0_dp
        b(1,6) = -b(1,5)
        b(3,7) = -b(3,5)
        b(1,8) = -b(1,5)
        b(3,8) = -b(3,5)
        b(2:3,9:12) = 1.0_dp
        b(2,10) = -b(2,9)
        b(3,11) = -b(3,9)
        b(2,12) = -b(2,9)
        b(3,12) = -b(3,9)
        b = b * b_length
        do ik = nks, 1, -1
           xk(:,( 13*(ik-1) )+1) = xk(:,ik)
           wk( ( 13*(ik-1) )+1 ) = wk(ik)
        end do
        do ik = nks, 1, -1
           do i = 1, 12
              xk(:,13*(ik-1)+i+1) = xk(:,( 13*(ik-1) )+1 ) + b(:,i)
              wk( 13 *(ik-1)+i+1 ) = 0.0_dp
           end do
        end do
        nks = 13 * nks
        deallocate(b)
     case default
        call errore ('set_kplusb', 'IBRAV not implemented', ibrav)
     end select
  end if

  return
end subroutine set_kplusb
