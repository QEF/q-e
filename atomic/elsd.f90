!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine elsd (mesh,zed,r,r2,dx,rho,zeta,vxt,vh,nlcc,   &
     nwf,enl,ll,lsd,nspin,oc,ndm, vnl, &
     etot,ekin,encl,epseu,ehrt,ecxc,evxt)    
  !---------------------------------------------------------------
  !
  !   atomic total energy in the local-spin-density scheme
  !   atomic pseudopotentials with nonlinear core correction are allowed
  !   gradient correction allowed (A. Dal Corso fecit AD 1993)
  !
  use kinds, only : DP
  use constants, only: fpi
  use funct, only: get_iexch, dft_is_gradient
  use ld1inc, only: vx
  implicit none
  integer:: ndm,mesh,nwf,i,n,ll(nwf),lam,lmax,lsd,nspin,is
  logical:: nlcc, gga, oep
  real(DP):: zed, int_0_inf_dr, rh(2),rhc,vxc,exc,vxcp(2), &
       etot,encl,epseu,ekin,ehrt,ecxc,evxt
  real(DP):: enl(nwf),oc(nwf), rhotot, exc_t, &
       r(ndm),r2(ndm),dx, rho(ndm,2),zeta(ndm), &
       vxt(ndm),vnl(ndm,0:3),vh(ndm)
  real(DP),allocatable :: f1(:), f2(:), f3(:), f4(:)
  real(DP),allocatable :: vgc(:,:), egc(:), rhoc(:)
  integer:: mgcx,mgcc,ierr

  gga=dft_is_gradient() 
  oep=get_iexch().eq.4
  allocate(vgc(ndm,2),stat=ierr)
  allocate(rhoc(ndm),stat=ierr)
  allocate(egc(ndm),stat=ierr)
  allocate(f1(mesh),stat=ierr)
  allocate(f2(mesh),stat=ierr)
  allocate(f3(mesh),stat=ierr)
  allocate(f4(mesh),stat=ierr)

  rhoc=0.0_DP
  call vxcgc(ndm,mesh,nspin,r,r2,rho,rhoc,vgc,egc)

  rhc=0.0_DP
  do i=1,mesh
     rhotot=rho(i,1)
     if (lsd.eq.1) rhotot=rhotot+rho(i,2) 
     f1(i)=-2.0_DP*zed/r(i) * rhotot
     f4(i)= vxt(i)       * rhotot
     vh(i)= vh (i)       * rhotot
     do is=1, nspin
        rh(is) = rho(i,is)/r2(i)/fpi
     enddo
     call vxc_t(rh,rhc,lsd,vxcp)
     if (gga) then
        f2(i) =-(vxcp(1)+vgc(i,1))*rho(i,1)-f1(i)-vh(i)-f4(i)
        f3(i) = exc_t(rh,rhc,lsd)*rhotot+egc(i)*r2(i)*fpi
        if (lsd.eq.1) f2(i)=f2(i)-(vxcp(2)+vgc(i,2))*rho(i,2)
     else
        f2(i) =-vxcp(1)*rho(i,1)-f1(i)-vh(i)-f4(i)
        f3(i) = exc_t(rh,rhc,lsd) * rhotot
        if (lsd.eq.1) f2(i) =f2(i)-vxcp(2)*rho(i,2)
     endif
     if (oep) then
        do is = 1, nspin
           f2(i) = f2(i) - vx(i,is)*rho(i,is)
        end do
     end if
  enddo

  encl=    int_0_inf_dr(f1,r,r2,dx,mesh,1)
  ehrt=0.5_DP*int_0_inf_dr(vh,r,r2,dx,mesh,2)
  ecxc=    int_0_inf_dr(f3,r,r2,dx,mesh,2)
  evxt=    int_0_inf_dr(f4,r,r2,dx,mesh,2)
  !
  epseu=0.0_DP
  ekin = int_0_inf_dr(f2,r,r2,dx,mesh,1)

  do n=1,nwf
     ekin=ekin+oc(n)*enl(n)
  enddo

  if (oep) call add_exchange (ecxc)

  etot= ekin + encl + ehrt + ecxc + evxt

  deallocate(f4)
  deallocate(f3)
  deallocate(f2)
  deallocate(f1)
  deallocate(egc)
  deallocate(vgc)
  deallocate(rhoc)

  return
end subroutine elsd
