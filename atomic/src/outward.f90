!---------------------------------------------------------------
subroutine outward(y,f,g,mesh,imatch,ncross)
!---------------------------------------------------------------
   ! outward integration. numerov method.
   !
   use kinds, only: DP
   use radial_grids, only: ndmx
   implicit none
   !
   ! I/O variables
   !
   integer,intent(in) :: mesh 
   integer, intent(in) :: imatch
   integer, intent(out) ::  ncross  ! num of axis crosses
   real(DP),intent(in) :: f(ndmx), g(ndmx)
   real (DP),intent(out) :: y(ndmx)

   !
   ! local variables
   !
   integer :: n ! lopp variable
   real (DP) :: ymx ! max value of the function over the grid 

   if (ndmx.lt.mesh) stop ' outward : ndmx .lt. mesh !!!!'
   !
   ncross=0
   ymx=0.d0
   do n=2,imatch-1
      y(n+1)=((12.d0-10.d0*f(n))*y(n)-f(n-1)*y(n-1)+g(n))/f(n+1)
      if ( y(n) .ne. sign(y(n),y(n+1)) ) ncross = ncross + 1
      ymx=max(ymx,abs(y(n)))
   end do
   if(ymx.ge.1.0d10) write (*,*) ' ******** ymx.ge.1.0e10 ********'
!
   return
end subroutine outward

