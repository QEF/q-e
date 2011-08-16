!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function bessj(n,x)
!
! It computes the Bessel functions J(n,x)
!
  USE kinds, only : DP

  implicit none
  integer, parameter :: iacc=40
  integer :: n, j, jsum, m
  real(DP), parameter :: bigno=1.d10, bigni=1.d-10
  real(DP) :: x, bessj, bessj0, bessj1, bj, bjm, bjp, &
                   sum, tox, ans

  if (n.lt.2) then
     if (n.eq.0) bessj=bessj0(x)
     if (n.eq.1) bessj=bessj1(x)
  else
     if (x.eq.0.d0) then
       ans=0.d0
     elseif (abs(x).gt.1.d0*n) then
        tox=2.d0/abs(x)
        bjm=bessj0(abs(x))
        bj=bessj1(abs(x))
        do j=1, n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
        enddo
        ans=bj
     else
        tox=2.d0/abs(x)
        m=2*((n+nint(sqrt(1.d0*(iacc*n))))/2)
        ans=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.d0
        do j=m, 1, -1
           bjm=j*tox*bj-bjp
           bjp=bj
           bj=bjm
           if (abs(bj).gt.bigno) then
              bj=bj*bigni
              bjp=bjp*bigni
              ans=ans*bigni
              sum=sum*bigni
           endif
           if (jsum.ne.0) sum=sum+bj
           jsum=1-jsum
           if (j.eq.n) ans=bjp
        enddo
        sum=2.d0*sum-bj
        ans=ans/sum
     endif
     if (0.d0.gt.x.and.mod(n,2).ne.n) ans=-ans
     bessj=ans
  endif

  return
end function bessj
!---------------------------------------

function bessj0(x)
  USE kinds, only : DP
  IMPLICIT NONE
  real(DP) :: x, ax, xx, z, y, ans, ans1, ans2, bessj0

  if (abs(x).lt.8.d0) then
     y=x**2
     ans1=57568490574.d0+y*(-13362590354.d0+y*(651619640.7d0+   &
      y*(-11214424.18d0+y*(77392.33017d0+y*(-184.9052456d0)))))
     ans2=57568490411.d0+y*(1029532985.d0+y*(9494680.718d0+     &
      y*(59272.64853d0+y*(267.8532712d0+y*1.d0))))
     bessj0=ans1/ans2
  else
     ax=abs(x)
     z=8.d0/ax
     y=z**2
     xx=ax-0.785398164d0
     ans1=1.d0+y*(-0.1098628627d-2+y*(0.2734510407d-4+         &
             y*(-0.2073370639d-5+y*0.2093887211d-6)))
     ans2=-0.1562499995d-1+y*(0.1430488765d-3+                 &
             y*(-0.6911147651d-5+y*(0.7621095161d-6-           &
             y*0.934945152d-7)))
     ans=sqrt(0.636619772d0/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
     bessj0=ans
  endif

  return
end function bessj0
!-----------------------------------

function bessj1(x)
  USE kinds, only : DP
  real(DP) :: x, ax, xx, y, z, ans, ans1, ans2, bessj1

  if (abs(x).le.8.d0) then
     y=x**2
     ans1=x*(72362614232.d0+y*(-7895059235.d0+                &
             y*(242396853.1d0+                                &
             y*(-2972611.439d0+y*(15704.48260d0+              &
             y*(-30.16036606d0))))))
     ans2=144725228442.d0+y*(2300535178.d0+y*(18583304.74d0+  &
             y*(99447.43394d0+y*(376.9991397d0+y*1.d0))))
     bessj1=ans1/ans2
  else
     ax=abs(x)
     z=8.d0/ax
     y=z**2
     xx=ax-2.356194491d0
     ans1=1.d0+y*(0.183105d-2+y*(-0.3516396496d-4+           &
            y*(0.2457520174d-5+y*(-0.240337019d-6))))
     ans2=0.04687499995d0+y*(-0.2002690873d-3+               &
            y*(0.8449199096d-5+                              &
            y*(-0.88228987d-6+y*0.105787412d-6)))
     ans=sqrt(0.636619772d0/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
     if (x.le.0.d0) ans=-ans
     bessj1=ans
  endif

  return
end function bessj1

