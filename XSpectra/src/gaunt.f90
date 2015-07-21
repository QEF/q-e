
! Copyright (C) 2002-2004 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gaunt
! !INTERFACE:
real(8) function gaunt(l1,l2,l3,m1,m2,m3)
! !INPUT/OUTPUT PARAMETERS:
!   l1, l2, l3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Gaunt coefficient given by
!   $$  \langle Y^{l_1}_{m_1}|Y^{l_2}_{m_2}|Y^{l_3}_{m_3} \rangle 
!    = (-1)^{m_1}\left[\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi} \right]
!    ^{\frac{1}{2}}
!    \begin{pmatrix} l_1 & l_2 & l_3 \\  0   & 0   & 0   \end{pmatrix}
!    \begin{pmatrix} l_1 & l_2 & l_3 \\ -m_1 & m_2 & m_3 \end{pmatrix}. $$
!   Suitable for $l_i$ less than 50.
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: l1
integer, intent(in) :: l2
integer, intent(in) :: l3
integer, intent(in) :: m1
integer, intent(in) :: m2
integer, intent(in) :: m3
! local variables
integer j,j1,j2,j3,jh
real(8) t1
! real constant 1/sqrt(4*pi)
real(8), parameter :: c1=0.28209479177387814347d0
! external functions
real(8) wigner3j,factr,factnm
external wigner3j,factr,factnm
if ((l1.lt.0).or.(l2.lt.0).or.(l3.lt.0).or.(abs(m1).gt.l1).or.(abs(m2).gt.l2) &
 .or.(abs(m3).gt.l3)) then
  write(*,*)
  write(*,'("Error(gaunt): non-physical arguments :")')
  write(*,'("l1 = ",I8," l2 = ",I8," l3 = ",I8)') l1,l2,l3
  write(*,'("m1 = ",I8," m2 = ",I8," m3 = ",I8)') m1,m2,m3
  write(*,*)
  stop
end if
if ((l1.gt.50).or.(l2.gt.50).or.(l3.gt.50)) then
  write(*,*)
  write(*,'("Error(gaunt): angular momenta out of range : ",3I8)') l1,l2,l3
  write(*,*)
  stop
end if
if (m1-m2-m3.ne.0) then
  gaunt=0.d0
  return
end if
j1=l2-l1+l3
j2=l1-l2+l3
j3=l1+l2-l3
if ((j1.lt.0).or.(j2.lt.0).or.(j3.lt.0)) then
  gaunt=0.d0
  return
end if
j=l1+l2+l3
if (mod(j,2).ne.0) then
  gaunt=0.d0
  return
end if
jh=j/2
t1=dble((-1)**(m1+jh))
t1=t1*sqrt(dble((2*l1+1)*(2*l2+1)*(2*l3+1))*factr(j1,j+1)*factnm(j2,1) &
 *factnm(j3,1))
t1=t1*factr(jh,jh-l1)/(factnm(jh-l2,1)*factnm(jh-l3,1))
gaunt=t1*c1*wigner3j(l1,l2,l3,-m1,m2,m3)
return
end function
!EOC
