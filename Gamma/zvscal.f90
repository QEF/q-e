!
subroutine zvscal(n,lda,m,v,zin,zout)
  implicit none
  integer :: n, lda, m
  real(kind=8) :: v(n), zin(2,lda,m), zout(2,lda,m)
  integer :: i,j
  !
  do j = 1,m
     do i = 1,n
        zout(1,i,j) = zin(1,i,j)*v(i)
        zout(2,i,j) = zin(2,i,j)*v(i)
     end do
  end do
  !
  return
end subroutine zvscal
