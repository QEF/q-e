!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine series(f,r,r2,b)
  !---------------------------------------------------------------
  !
  use kinds, only : DP
  implicit none
  real(kind=dp):: dr21,dr31,dr32,dr41,dr42,dr43,df21,df32,df43, &
       ddf42,ddf31 
  real(kind=dp):: f(4),r(4),r2(4),b(0:3)
  dr21=r(2)-r(1)
  dr31=r(3)-r(1)
  dr32=r(3)-r(2)
  dr41=r(4)-r(1)
  dr42=r(4)-r(2)
  dr43=r(4)-r(3)
  df21=(f(2)-f(1))/dr21
  df32=(f(3)-f(2))/dr32
  df43=(f(4)-f(3))/dr43
  ddf42=(df43-df32)/dr42
  ddf31=(df32-df21)/dr31
  b(3)=(ddf42-ddf31)/dr41
  b(2)=ddf31-b(3)*(r(1)+r(2)+r(3))
  b(1)=df21-b(2)*(r(2)+r(1))-b(3)*(r2(1)+r2(2)+r(1)*r(2))
  b(0)=f(1)-r(1)*(b(1)+r(1)*(b(2)+r(1)*b(3)))
  return
end subroutine series
