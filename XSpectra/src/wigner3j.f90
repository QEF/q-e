
! Copyright (C) 2002-2004 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: wigner3j
! !INTERFACE:
real(8) function wigner3j(j1,j2,j3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Wigner $3j$-symbol. There are many equivalent definitions for
!   the $3j$-symbols, the following provides high accuracy for $j\le 50$
!   \begin{align}
!    &\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}=(-1)^
!    {j1+j2+m3}\nonumber\\
!    &\times\sqrt{\frac{(j_1+m_1)!(j_2+m_2)!(j_3+m_3)!(j_3-m_3)!(j_1-m_1)!
!    (j_2-m_2)!}{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!(1+j_1+j_2+j_3)!}}
!    \times\sum_{\max(0,j_2-j_3-m_1,j_1-j_3+m_2)}^
!    {\min(j_1+j_2-j_3,j_1-m_1,j_2+m_2)}\nonumber\\
!    &(-1)^k\frac{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!}
!    {(j_3-j_1-m_2+k)!(j_3-j_2+m_1+k)!(j_1+j_2-j_3-k)!k!(j_1-m_1-k)!
!    (j_2+m_2-k)}.\nonumber
!   \end{align}
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: j1
integer, intent(in) :: j2
integer, intent(in) :: j3
integer, intent(in) :: m1
integer, intent(in) :: m2
integer, intent(in) :: m3
! local variables
integer k,k1,k2,l1,l2,l3,n1,n2
real(8) sgn,sum,t1
! external functions
real(8) factnm,factr
external factnm,factr
! check input variables
if ((j1.lt.0).or.(j2.lt.0).or.(j3.lt.0).or.(abs(m1).gt.j1).or.(abs(m2).gt.j2) &
 .or.(abs(m3).gt.j3)) then
  write(*,*)
  write(*,'("Error(wigner3j): non-physical arguments :")')
  write(*,'("j1 = ",I8," j2 = ",I8," j3 = ",I8)') j1,j2,j3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((j1.eq.0).and.(j2.eq.0).and.(j3.eq.0)) then
  wigner3j=1.d0
  return
end if
if ((j1.gt.50).or.(j2.gt.50).or.(j3.gt.50)) then
  write(*,*)
  write(*,'("Error(wigner3j): angular momenta out of range : ",3I8)') j1,j2,j3
  write(*,*)
  stop
end if
l1=j2-j1+j3
l2=j1-j2+j3
l3=j1+j2-j3
if ((m1+m2+m3.ne.0).or.(l1.lt.0).or.(l2.lt.0).or.(l3.lt.0)) then
  wigner3j=0.d0
  return
end if
n1=j1-m1
n2=j2+m2
k1=max(0,j2-j3-m1,j1-j3+m2)
k2=min(l3,n1,n2)
sgn=dble((-1)**(k1+j1+j2+m3))
sum=0.d0
do k=k1,k2
  t1=sgn*factr(l1,l1-n2+k)*factr(l2,l2-n1+k)*factr(l3,l3-k)
  sum=sum+t1/(factnm(k,1)*factnm(n1-k,1)*factnm(n2-k,1))
  sgn=-sgn
end do
t1=factr(j1+m1,l1)*factr(j2+m2,l2)*factr(j3+m3,l3)
t1=t1*factr(j3-m3,1+j1+j2+j3)*factnm(j1-m1,1)*factnm(j2-m2,1)
wigner3j=sum*sqrt(t1)
return
end function
!EOC
