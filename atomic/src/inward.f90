!
!---------------------------------------------------------------
subroutine inward(y,f,g,mesh,imatch)
!---------------------------------------------------------------
   ! inward integration, charlotte froese can j phys 41,1895(1963)
   !
   use kinds,  only:DP
   use radial_grids, only:ndmx
   implicit none
   !
   ! I/O variables
   !
   integer :: mesh, imatch
   real(DP) :: y(ndmx), f(ndmx), g(ndmx)
   !
   ! local variables
   !
   integer :: imp1, imp2, nm1, n
   real(DP) :: di, ff, expn
   real(DP) :: el(ndmx),c(ndmx)

   if (ndmx.lt.mesh) stop ' inward : ndmx .lt. mesh !!!!'
! set up el, and c vectors
   imp1=imatch+1
   el(imp1)=10.d0*f(imp1)-12.d0
   c(imp1)=g(imp1)-f(imatch)*y(imatch)
   imp2=imatch+2
   nm1=mesh-1
   do n=imp2,nm1
      di=10.d0*f(n)-12.d0
      el(n)=di-f(n)*f(n-1)/el(n-1)
      c(n)=g(n)-c(n-1)*f(n-1)/el(n-1)
   end do
! start with y=aexp(-ax)-b, with b="g"/f
   ff=12.d0*abs(1.d0-f(nm1))
   expn=exp(-sqrt(ff))
   y(mesh)=(c(nm1)*expn+g(nm1)*el(nm1)*(expn-1.d0)/ff)/(el(nm1)+f(mesh)*expn)
! inward integration
   do n=nm1,imp1,-1
      y(n)=(c(n)-f(n+1)*y(n+1))/el(n)
   end do
!
   return
end subroutine inward
