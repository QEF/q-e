!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine starting_potential &
     (ndm,mesh,zval,zed,nwf,oc,nn,ll,r,enl, &
     v0,vxt,vpot,enne,nspin)
  !---------------------------------------------------------------
  !
  ! starting potential: generalized thomas-fermi atomic potential
  ! it is assumed that the effective charge seen by the reference
  ! electron cannot be smaller than 1 (far from the core)
  !
  use kinds, only : DP
  use ld1inc, only : frozen_core, noscf
  implicit none
  integer :: nwf, nn(nwf), ll(nwf), ndm, mesh, n, i, nspin
  real(DP) :: r(ndm), vpot(ndm,2), v0(ndm), vxt(ndm), enl(nwf), oc(nwf), &
       zed, zval, zz, zen, enne, t,x, vext, oce
  real(DP), parameter :: e2 = 2.0_dp
  external vext
  !
  enne = 0.0_dp
  zz = max(zed,zval)
  do  n=1,nwf
     oce=max(0.0_dp, oc(n))
     enne = enne + oce
     zen= 0.0_dp
     do  i=1,nwf
        oce=max(0.0_dp, oc(i))
        if(nn(i).lt.nn(n)) zen=zen+oce
        if(nn(i).eq.nn(n).and.ll(i).le.ll(n)) zen=zen+oce
     end do
     zen = max(zz-zen+1.0_dp,1.0_dp)
     if (ABS(enl(n))<1.d-7.or..not.frozen_core) enl(n) =-(zen/nn(n))**2
  end do
  !
  do  i=1,mesh
     vxt(i)=vext(r(i))
     x =r(i)*enne**(1.0_dp/3.0_dp)/0.885_dp
     t= zz/(1.0_DP+sqrt(x)*(0.02747_dp-x*(0.1486_dp-0.007298_dp*x)) &
          + x*(1.243_dp+x*(0.2302_dp+0.006944_dp*x)))
     t = max(1.0_dp,t)
     v0(i)= -e2 * zed / r(i)
     if (noscf) then
        vpot(i,1) = v0(i)+vxt(i)
     else
        vpot(i,1) = -e2*t/r(i) + vxt(i)
     endif
  enddo
  !
  if (nspin.eq.2) then
     do i=1,mesh
        vpot(i,2)=vpot(i,1)
     enddo
  endif
  !
  return
end subroutine starting_potential
