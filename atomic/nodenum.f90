!
!----------------------------------------------------------------------------
subroutine nodeno(snlo,jj1,jj2,nodes,idim1)
!----------------------------------------------------------------------------
!
!   routine counts the number of nodes of the wavefunction snlo
!   between the points jj1 and jj2
!
implicit none
integer,parameter :: dp = kind(1.d0)
integer :: jj1,jj2,nodes,idim1
!
!   wavefunction array
!
real(kind=dp) :: snlo(idim1)

integer :: i
!
nodes = 0
!
do i = jj1+1,jj2
   if ( snlo(i-1) * snlo(i) .lt. 0.0d0 ) nodes = nodes + 1
enddo
return
end
