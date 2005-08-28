!-----------------------------------------------------------------------
subroutine set_sl3(sl3,lmax)
!-----------------------------------------------------------------------
   ! compute sl3(l1,l2,l3)=\int pl(x,l1) pl(x,l2) pl(x,l3) dx
   !
   use kinds, only: DP
   implicit none
   !
   ! I/O variables
   !
   integer :: lmax
   real (DP) :: sl3(0:2*lmax,0:2*lmax,0:2*lmax)
   !
   ! local variables
   !
   integer :: lmax2, l1, l2, l3, l1i, l1f

   lmax2 = 2*lmax
   do l1=0,lmax2
      do l2=0,lmax2
         do l3=0,lmax2
            sl3(l1,l2,l3) = 0.0d0
         end do
      end do
   end do

   do l3=0,lmax2
      sl3(l3,0,l3) = 2.0d0/(2.0d0*l3+1.0d0)
      sl3(0,l3,l3) = 2.0d0/(2.0d0*l3+1.0d0)
      do l2=1,lmax2
         l1i = max(abs(l2-l3),1)
         l1f = lmax2-abs(lmax2-l2-l3)
         if (l1f.eq.lmax2) then
            sl3(lmax2,l2,l3) = (2.0d0*l2-1.0d0)/(2.0d0*l1f+1.0d0) *  &
                               l1f*sl3(l1f-1,l2-1,l3)/l2
            l1f=l1f-1
         end if
         do l1=l1i,l1f
            sl3(l1,l2,l3) = (2.0d0*l2-1.0d0)/(2.0d0*l1+1.0d0) * &
               ((l1+1.0d0)*sl3(l1+1,l2-1,l3)+l1*sl3(l1-1,l2-1,l3))/l2
            if (l2.gt.1) sl3(l1,l2,l3) = sl3(l1,l2,l3) -(l2-1.0d0)/l2 * &
                                        sl3(l1,l2-2,l3)
         end do
      end do
   end do
   return
end subroutine set_sl3
