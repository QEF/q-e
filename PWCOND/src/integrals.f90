!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Optimized Aug. 2004 (ADC)
!
!
function int1d(fun, zk, dz, dz1, nz1, tpiba, sign)
!
! This function computes the integral of beta function with the
! exponential
!
  USE kinds, only : DP
  implicit none
  integer ::  &
      ik,     &        ! counter on slab points
      nz1,    &        ! input: the number of integration points
      sign             ! input: the sign of the exponential
  real(DP), parameter :: eps=1.d-8
  real(DP) :: tpi, dz, dz1, tpiba
  complex(DP), parameter :: cim = (0.d0,1.d0)
  complex(DP) ::  &
      zk,              & ! the exponential k
      fun(nz1),        & ! the beta function on the slab points
      fact,fact0,      & ! auxiliary
      arg,             & ! auxiliary
      int1d              ! output: the value of the integral

  tpi = 8.d0*atan(1.d0)

  int1d = (0.d0,0.d0)
  arg = sign*tpi*cim*zk*dz1
  fact0=exp(arg)
  fact=fact0
  do ik=1, nz1
    int1d = int1d+CONJG(fun(ik))*fact
    fact=fact*fact0
  enddo
  if (abs(DBLE(zk))+abs(AIMAG(zk)).gt.eps) then
    int1d =-sign*cim*int1d*(1.d0-exp(-arg))/(zk*tpiba)
    if (sign.lt.0) int1d=int1d*exp(tpi*cim*zk*dz)
  else
    int1d = int1d*dz1/tpiba*tpi
  endif

  return
end function int1d
!-----------------------------------
!
function int2d(fun1, fun2, int1, int2, fact1, fact2, zk, dz1, tpiba, nz1 )
!
! This function computes the 2D integrals of beta functions with
! exponential
!
  USE kinds, only : DP
  USE constants, ONLY : tpi
  implicit none
  integer ::    &
     nz1,       &  ! number of points for the slab integration
     ik       ! counters on the slab points
  real(DP), parameter :: eps=1.d-8
  real(DP) :: dz1, tpiba
  complex(DP), parameter :: cim=(0.d0,1.d0), one=(1.d0,0.d0)
  complex(DP) ::       &
     fun1(nz1), fun2(nz1),  &  ! the two arrays to be integrated
     int1(nz1), int2(nz1),  &  ! auxiliary arrays for integration
     fact1(nz1), fact2(nz1),&
     s1, s2, s3, ff,        &  ! auxiliary for integration
     fact,fact0,            &  ! auxiliary
     f1, f2, zk,            &  ! the complex k of the exponent
     int2d                     ! output: the result of the integration

  s1=(0.d0,0.d0)
  s2=(0.d0,0.d0)
  s3=(0.d0,0.d0)
!
! integral for i > = j
!
  fact=fact1(1)
  fact0=fact2(1)
  do ik=1, nz1
    ff=CONJG(fun1(ik))
    s1=s1+int1(ik)*ff*fact1(ik)
    s2=s2+int2(ik)*ff*fact2(ik)
    s3=s3+fun2(ik)*ff
  enddo
!
! complete integral
!
  f1=cim*zk*dz1*tpi
  f2=one/(zk*tpiba)**2
  if (abs(f1).gt.eps) then
     int2d=((1.d0-fact+f1)*s3*2.d0+(2.d0-fact-fact0)*(s1+s2))*f2
  else
     int2d=(s1+s2+s3)*(dz1*tpi/tpiba)**2
  endif

  return
end function int2d

subroutine setint(fun,int1,int2,fact1,fact2,nz1)

  USE kinds, only : DP
  implicit none
  integer ::    &
     nz1,       &  ! number of points for the slab integration
     ik            ! counters on the slab points
  complex(DP) ::       &
     fun(nz1),              &  ! the arrays to be integrated
     int1(nz1), int2(nz1),  &  ! auxiliary arrays for integration
     fact1(nz1), fact2(nz1)    !
!
  int1(1)=(0.d0, 0.d0)
  int2(nz1)=(0.d0, 0.d0)
  do ik=2, nz1
     int1(ik)=int1(ik-1)+fun(ik-1)*fact2(ik-1)
  enddo
  do ik=nz1-1,1,-1
     int2(ik)=int2(ik+1)+fun(ik+1)*fact1(ik+1)
  enddo

  return
end subroutine setint
