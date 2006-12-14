!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine esic
  !---------------------------------------------------------------
  !
  use ld1inc
  implicit none
  ! output
  ! local
  integer:: n, i
  real(DP) :: int_0_inf_dr,vxup,vxdw,vcup,vcdw,ex,ec,deksic  
  real(DP) :: work1(ndm),v(ndm),vsic(ndm)
  real(DP) :: egc(ndm)
  external int_0_inf_dr
  !
  deksic = 0.0_DP
  dhrsic = 0.0_DP
  dxcsic = 0.0_DP
  do n=1,nwf
     call sic_correction(n,v,vsic,work1)
     if (rel.eq.2) then
        do i=1,mesh
           v(i)=v(i)*(psi(i,1,n)**2+psi(i,2,n)**2)
           vsic(i)=vsic(i)*(psi(i,1,n)**2+psi(i,2,n)**2)
        end do
     else
        do i=1,mesh
           v(i)=v(i)*psi(i,1,n)**2
           vsic(i)=vsic(i)*psi(i,1,n)**2
        enddo
     endif
     deksic = deksic +  &
          &      oc(n)*int_0_inf_dr(vsic,r,r2,dx,mesh,2*(ll(n)+1))
     dhrsic = dhrsic - 0.5_DP*oc(n)*int_0_inf_dr(v    ,r,r2,dx,mesh,2)
     dxcsic = dxcsic - oc(n)*int_0_inf_dr(work1,r,r2,dx,mesh,2)
  enddo
  !
  ekin=ekin+deksic
  ehrt=ehrt
  ecxc=ecxc+dxcsic+dhrsic
  etot=etot+dhrsic+dxcsic+deksic
  !
  return
end subroutine esic
