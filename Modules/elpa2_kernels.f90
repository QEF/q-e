! --------------------------------------------------------------------------------------------------
!
! This file contains the compute intensive kernels for the Householder transformations.
! It should be compiled with the highest possible optimization level.
!
! On Intel use -O3 -xSSE4.2 (or the SSE level fitting to your CPU)
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
! --------------------------------------------------------------------------------------------------

subroutine double_hh_trafo(q, hh, nb, nq, ldq, ldh)

   implicit none

   integer, intent(in) :: nb, nq, ldq, ldh
   real*8, intent(inout) :: q(ldq,*)
   real*8, intent(in) :: hh(ldh,*)

   real*8 s
   integer i

   ! Safety only:

   if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

   ! Calculate dot product of the two Householder vectors

   s = hh(2,2)*1
   do i=3,nb
      s = s+hh(i,2)*hh(i-1,1)
   enddo

   ! Do the Householder transformations

   ! Always a multiple of 4 Q-rows is transformed, even if nq is smaller

   do i=1,nq-8,12
      call hh_trafo_kernel_12(q(i,1),hh, nb, ldq, ldh, s)
   enddo

   ! i > nq-8 now, i.e. at most 8 rows remain

   if(nq-i+1 > 4) then
      call hh_trafo_kernel_8(q(i,1),hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 0) then
      call hh_trafo_kernel_4(q(i,1),hh, nb, ldq, ldh, s)
   endif

end

! --------------------------------------------------------------------------------------------------
! The following kernels perform the Householder transformation on Q for 12/8/4 rows.
! Please note that Q is declared complex*16 here.
! This is a hint for compilers that packed arithmetic can be used for Q
! (relevant for Intel SSE and BlueGene double hummer CPUs).
! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_12(q, hh, nb, ldq, ldh, s)

   implicit none

   integer, intent(in) :: nb, ldq, ldh
   complex*16, intent(inout) :: q(ldq/2,*)
   real*8, intent(in) :: hh(ldh,*), s

   complex*16 x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6
   real*8 h1, h2, tau1, tau2
   integer i


   x1 = q(1,2)
   x2 = q(2,2)
   x3 = q(3,2)
   x4 = q(4,2)
   x5 = q(5,2)
   x6 = q(6,2)

   y1 = q(1,1) + q(1,2)*hh(2,2)
   y2 = q(2,1) + q(2,2)*hh(2,2)
   y3 = q(3,1) + q(3,2)*hh(2,2)
   y4 = q(4,1) + q(4,2)*hh(2,2)
   y5 = q(5,1) + q(5,2)*hh(2,2)
   y6 = q(6,1) + q(6,2)*hh(2,2)

!DEC$ VECTOR ALIGNED
   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      x1 = x1 + q(1,i)*h1
      y1 = y1 + q(1,i)*h2
      x2 = x2 + q(2,i)*h1
      y2 = y2 + q(2,i)*h2
      x3 = x3 + q(3,i)*h1
      y3 = y3 + q(3,i)*h2
      x4 = x4 + q(4,i)*h1
      y4 = y4 + q(4,i)*h2
      x5 = x5 + q(5,i)*h1
      y5 = y5 + q(5,i)*h2
      x6 = x6 + q(6,i)*h1
      y6 = y6 + q(6,i)*h2
   enddo

   x1 = x1 + q(1,nb+1)*hh(nb,1)
   x2 = x2 + q(2,nb+1)*hh(nb,1)
   x3 = x3 + q(3,nb+1)*hh(nb,1)
   x4 = x4 + q(4,nb+1)*hh(nb,1)
   x5 = x5 + q(5,nb+1)*hh(nb,1)
   x6 = x6 + q(6,nb+1)*hh(nb,1)

   tau1 = hh(1,1)
   tau2 = hh(1,2)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   h1 = -tau2
   h2 = -tau2*s
   y1 = y1*h1 + x1*h2
   y2 = y2*h1 + x2*h2
   y3 = y3*h1 + x3*h2
   y4 = y4*h1 + x4*h2
   y5 = y5*h1 + x5*h2
   y6 = y6*h1 + x6*h2

   q(1,1) = q(1,1) + y1
   q(2,1) = q(2,1) + y2
   q(3,1) = q(3,1) + y3
   q(4,1) = q(4,1) + y4
   q(5,1) = q(5,1) + y5
   q(6,1) = q(6,1) + y6
   q(1,2) = q(1,2) + x1 + y1*hh(2,2)
   q(2,2) = q(2,2) + x2 + y2*hh(2,2)
   q(3,2) = q(3,2) + x3 + y3*hh(2,2)
   q(4,2) = q(4,2) + x4 + y4*hh(2,2)
   q(5,2) = q(5,2) + x5 + y5*hh(2,2)
   q(6,2) = q(6,2) + x6 + y6*hh(2,2)

!DEC$ VECTOR ALIGNED
   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      q(1,i) = q(1,i) + x1*h1 + y1*h2
      q(2,i) = q(2,i) + x2*h1 + y2*h2
      q(3,i) = q(3,i) + x3*h1 + y3*h2
      q(4,i) = q(4,i) + x4*h1 + y4*h2
      q(5,i) = q(5,i) + x5*h1 + y5*h2
      q(6,i) = q(6,i) + x6*h1 + y6*h2
   enddo

   q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
   q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)
   q(3,nb+1) = q(3,nb+1) + x3*hh(nb,1)
   q(4,nb+1) = q(4,nb+1) + x4*hh(nb,1)
   q(5,nb+1) = q(5,nb+1) + x5*hh(nb,1)
   q(6,nb+1) = q(6,nb+1) + x6*hh(nb,1)

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_8(q, hh, nb, ldq, ldh, s)

   implicit none

   integer, intent(in) :: nb, ldq, ldh
   complex*16, intent(inout) :: q(ldq/2,*)
   real*8, intent(in) :: hh(ldh,*), s

   complex*16 x1, x2, x3, x4, y1, y2, y3, y4
   real*8 h1, h2, tau1, tau2
   integer i


   x1 = q(1,2)
   x2 = q(2,2)
   x3 = q(3,2)
   x4 = q(4,2)

   y1 = q(1,1) + q(1,2)*hh(2,2)
   y2 = q(2,1) + q(2,2)*hh(2,2)
   y3 = q(3,1) + q(3,2)*hh(2,2)
   y4 = q(4,1) + q(4,2)*hh(2,2)

!DEC$ VECTOR ALIGNED
   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      x1 = x1 + q(1,i)*h1
      y1 = y1 + q(1,i)*h2
      x2 = x2 + q(2,i)*h1
      y2 = y2 + q(2,i)*h2
      x3 = x3 + q(3,i)*h1
      y3 = y3 + q(3,i)*h2
      x4 = x4 + q(4,i)*h1
      y4 = y4 + q(4,i)*h2
   enddo

   x1 = x1 + q(1,nb+1)*hh(nb,1)
   x2 = x2 + q(2,nb+1)*hh(nb,1)
   x3 = x3 + q(3,nb+1)*hh(nb,1)
   x4 = x4 + q(4,nb+1)*hh(nb,1)

   tau1 = hh(1,1)
   tau2 = hh(1,2)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   h1 = -tau2
   h2 = -tau2*s
   y1 = y1*h1 + x1*h2
   y2 = y2*h1 + x2*h2
   y3 = y3*h1 + x3*h2
   y4 = y4*h1 + x4*h2

   q(1,1) = q(1,1) + y1
   q(2,1) = q(2,1) + y2
   q(3,1) = q(3,1) + y3
   q(4,1) = q(4,1) + y4
   q(1,2) = q(1,2) + x1 + y1*hh(2,2)
   q(2,2) = q(2,2) + x2 + y2*hh(2,2)
   q(3,2) = q(3,2) + x3 + y3*hh(2,2)
   q(4,2) = q(4,2) + x4 + y4*hh(2,2)

!DEC$ VECTOR ALIGNED
   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      q(1,i) = q(1,i) + x1*h1 + y1*h2
      q(2,i) = q(2,i) + x2*h1 + y2*h2
      q(3,i) = q(3,i) + x3*h1 + y3*h2
      q(4,i) = q(4,i) + x4*h1 + y4*h2
   enddo

   q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
   q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)
   q(3,nb+1) = q(3,nb+1) + x3*hh(nb,1)
   q(4,nb+1) = q(4,nb+1) + x4*hh(nb,1)

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_4(q, hh, nb, ldq, ldh, s)

   implicit none

   integer, intent(in) :: nb, ldq, ldh
   complex*16, intent(inout) :: q(ldq/2,*)
   real*8, intent(in) :: hh(ldh,*), s

   complex*16 x1, x2, y1, y2
   real*8 h1, h2, tau1, tau2
   integer i


   x1 = q(1,2)
   x2 = q(2,2)

   y1 = q(1,1) + q(1,2)*hh(2,2)
   y2 = q(2,1) + q(2,2)*hh(2,2)

!DEC$ VECTOR ALIGNED
   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      x1 = x1 + q(1,i)*h1
      y1 = y1 + q(1,i)*h2
      x2 = x2 + q(2,i)*h1
      y2 = y2 + q(2,i)*h2
   enddo

   x1 = x1 + q(1,nb+1)*hh(nb,1)
   x2 = x2 + q(2,nb+1)*hh(nb,1)

   tau1 = hh(1,1)
   tau2 = hh(1,2)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   h1 = -tau2
   h2 = -tau2*s
   y1 = y1*h1 + x1*h2
   y2 = y2*h1 + x2*h2

   q(1,1) = q(1,1) + y1
   q(2,1) = q(2,1) + y2
   q(1,2) = q(1,2) + x1 + y1*hh(2,2)
   q(2,2) = q(2,2) + x2 + y2*hh(2,2)

!DEC$ VECTOR ALIGNED
   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      q(1,i) = q(1,i) + x1*h1 + y1*h2
      q(2,i) = q(2,i) + x2*h1 + y2*h2
   enddo

   q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
   q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)

end

! --------------------------------------------------------------------------------------------------

subroutine single_hh_trafo_complex(q, hh, nb, nq, ldq)

   implicit none

   integer, intent(in) :: nb, nq, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   integer i

   ! Safety only:

   if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

   ! Do the Householder transformations

   ! Always a multiple of 4 Q-rows is transformed, even if nq is smaller

   do i=1,nq-8,12
      call hh_trafo_complex_kernel_12(q(i,1),hh, nb, ldq)
   enddo

   ! i > nq-8 now, i.e. at most 8 rows remain

   if(nq-i+1 > 4) then
      call hh_trafo_complex_kernel_8(q(i,1),hh, nb, ldq)
   else if(nq-i+1 > 0) then
      call hh_trafo_complex_kernel_4(q(i,1),hh, nb, ldq)
   endif

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_complex_kernel_12(q, hh, nb, ldq)

   implicit none

   integer, intent(in) :: nb, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   complex*16 x1, x2, x3, x4, x5, x6, x7, x8, x9, xa, xb, xc
   complex*16 h1, tau1
   integer i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)
   x5 = q(5,1)
   x6 = q(6,1)
   x7 = q(7,1)
   x8 = q(8,1)
   x9 = q(9,1)
   xa = q(10,1)
   xb = q(11,1)
   xc = q(12,1)

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = conjg(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
      x5 = x5 + q(5,i)*h1
      x6 = x6 + q(6,i)*h1
      x7 = x7 + q(7,i)*h1
      x8 = x8 + q(8,i)*h1
      x9 = x9 + q(9,i)*h1
      xa = xa + q(10,i)*h1
      xb = xb + q(11,i)*h1
      xc = xc + q(12,i)*h1
   enddo

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   x7 = x7*h1
   x8 = x8*h1
   x9 = x9*h1
   xa = xa*h1
   xb = xb*h1
   xc = xc*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4
   q(5,1) = q(5,1) + x5
   q(6,1) = q(6,1) + x6
   q(7,1) = q(7,1) + x7
   q(8,1) = q(8,1) + x8
   q(9,1) = q(9,1) + x9
   q(10,1) = q(10,1) + xa
   q(11,1) = q(11,1) + xb
   q(12,1) = q(12,1) + xc

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
      q(5,i) = q(5,i) + x5*h1
      q(6,i) = q(6,i) + x6*h1
      q(7,i) = q(7,i) + x7*h1
      q(8,i) = q(8,i) + x8*h1
      q(9,i) = q(9,i) + x9*h1
      q(10,i) = q(10,i) + xa*h1
      q(11,i) = q(11,i) + xb*h1
      q(12,i) = q(12,i) + xc*h1
   enddo

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_complex_kernel_8(q, hh, nb, ldq)

   implicit none

   integer, intent(in) :: nb, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   complex*16 x1, x2, x3, x4, x5, x6, x7, x8
   complex*16 h1, tau1
   integer i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)
   x5 = q(5,1)
   x6 = q(6,1)
   x7 = q(7,1)
   x8 = q(8,1)

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = conjg(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
      x5 = x5 + q(5,i)*h1
      x6 = x6 + q(6,i)*h1
      x7 = x7 + q(7,i)*h1
      x8 = x8 + q(8,i)*h1
   enddo

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1
   x5 = x5*h1
   x6 = x6*h1
   x7 = x7*h1
   x8 = x8*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4
   q(5,1) = q(5,1) + x5
   q(6,1) = q(6,1) + x6
   q(7,1) = q(7,1) + x7
   q(8,1) = q(8,1) + x8

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
      q(5,i) = q(5,i) + x5*h1
      q(6,i) = q(6,i) + x6*h1
      q(7,i) = q(7,i) + x7*h1
      q(8,i) = q(8,i) + x8*h1
   enddo

end

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_complex_kernel_4(q, hh, nb, ldq)

   implicit none

   integer, intent(in) :: nb, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   complex*16 x1, x2, x3, x4
   complex*16 h1, tau1
   integer i


   x1 = q(1,1)
   x2 = q(2,1)
   x3 = q(3,1)
   x4 = q(4,1)

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = conjg(hh(i))
      x1 = x1 + q(1,i)*h1
      x2 = x2 + q(2,i)*h1
      x3 = x3 + q(3,i)*h1
      x4 = x4 + q(4,i)*h1
   enddo

   tau1 = hh(1)

   h1 = -tau1
   x1 = x1*h1
   x2 = x2*h1
   x3 = x3*h1
   x4 = x4*h1

   q(1,1) = q(1,1) + x1
   q(2,1) = q(2,1) + x2
   q(3,1) = q(3,1) + x3
   q(4,1) = q(4,1) + x4

!DEC$ VECTOR ALIGNED
   do i=2,nb
      h1 = hh(i)
      q(1,i) = q(1,i) + x1*h1
      q(2,i) = q(2,i) + x2*h1
      q(3,i) = q(3,i) + x3*h1
      q(4,i) = q(4,i) + x4*h1
   enddo

end

! --------------------------------------------------------------------------------------------------
