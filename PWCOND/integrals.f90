!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function int1d(fun, zk, dz, dz1, nz1, tpiba, sign)
!
! This function computes the integral of beta function with the
! exponential
!
  use parameters, only : DP
  implicit none
  integer ::  &
      ik,     &        ! counter on slab points
      nz1,    &        ! input: the number of integration points
      sign             ! input: the sign of the exponential
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: tpi, dz, dz1, tpiba 
  complex(kind=DP), parameter :: cim = (0.d0,1.d0)
  complex(kind=DP) ::  &
      zk,              & ! the exponential k
      fun(nz1),        & ! the beta function on the slab points
      arg,             & ! auxiliary
      int1d              ! output: the value of the integral

  tpi = 8.d0*atan(1.d0)       

  int1d = (0.d0,0.d0)
  arg = sign*tpi*cim*zk*dz1
  do ik=1, nz1
    int1d = int1d+conjg(fun(ik))*exp(arg*ik)    
  enddo
  if (abs(real(zk))+abs(DIMAG(zk)).gt.eps) then
    int1d =-sign*cim*int1d*(1.d0-exp(-arg))/(zk*tpiba)
    if (sign.lt.0) int1d=int1d*exp(tpi*cim*zk*dz)
  else
    int1d = int1d*dz1/tpiba*tpi
  endif

  return
end function int1d
!-----------------------------------

function int2d(fun1, fun2, zk, dz1, nz1, tpiba)
!
! This function computes the 2D integrals of beta functions with
! exponential
!
  use parameters, only : DP 
  implicit none
  integer ::    &
     nz1,       &  ! number of points for the slab integration 
     ik, ik1       ! counters on the slab points
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: dz1, tpiba, tpi
  complex(kind=DP), parameter :: cim=(0.d0,1.d0)
  complex(kind=DP) ::       &
     fun1(nz1), fun2(nz1),  &  ! the two arrays to be integrated
     arg, s1, s2, s3, s,    &  ! auxyliary for integration
     zk,                    &  ! the complex k of the exponent
     int2d                     ! output: the result of the integration

  int2d=(0.d0, 0.d0)
  tpi=8.d0*atan(1.d0)
  s1=0.d0
  s2=0.d0
  s3=0.d0 
  arg=tpi*cim*zk*dz1
!
! integral for i > = j
!
  do ik=1, nz1
    s=0.d0
    do ik1=1, ik-1
       s=s+fun2(ik1)*exp(arg*(ik-ik1))
    enddo
    s1=s1+s*conjg(fun1(ik))
    s3=s3+conjg(fun1(ik))*fun2(ik)
  enddo
!
! integral for i < j
!
  do ik=1, nz1
    s=0.d0
    do ik1=ik+1, nz1
      s=s+fun2(ik1)*exp(arg*(ik1-ik))
    enddo
    s2=s2+s*conjg(fun1(ik))
  enddo 
!
! complete integral
!
  if (abs(DREAL(zk))+abs(DIMAG(zk)).gt.eps) then
     int2d=((1.d0-exp(arg)+cim*zk*dz1*tpi)*s3*2.d0+      &
           (2.d0-exp(arg)-exp(-arg))*(s1+s2))/((zk*tpiba)**2)
  else
     int2d=(s1+s2+s3)*(dz1*tpi/tpiba)**2
  endif

  return
end function int2d
