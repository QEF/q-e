!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine new_potential &
     (ndm,mesh,r,r2,sqr,dx,zed,vxt,lsd, &
     nlcc,latt,enne,rhoc,rho,vh,vnew)
  !---------------------------------------------------------------
  !   set up the selfconsistent atomic potential
  !
  use constants, only: fpi, e2
  use kinds, only : DP
  use funct, only : get_iexch, get_igcx, get_igcc
  use ld1inc, only : nwf, vx
  implicit none
  logical :: nlcc, gga, oep
  integer :: ndm,mesh,lsd,latt,i,is,nu, nspin, ierr
  real(DP):: rho(ndm,2),r(ndm),r2(ndm),vxcp(2), &
       sqr(ndm),vnew(ndm,2),vxt(ndm),vh(ndm), rhoc(ndm)
  real(DP):: zed,enne,dx,rh(2),rhc
  real(DP),allocatable:: vgc(:,:), egc(:), rhotot(:)
!  real(DP),allocatable:: vx(:,:)
  real(DP),allocatable:: dchi0(:,:)

  gga=get_igcx().ne.0.or.get_igcc().ne.0
  oep=get_iexch().eq.4
  nspin=1
  if (lsd.eq.1) nspin=2
  !
  !   compute hartree potential with the total charge
  !
  allocate(rhotot(ndm),stat=ierr)
  do i=1,ndm
     rhotot(i)=rho(i,1)
  enddo
  if (lsd.eq.1) then
     do i=1,ndm
        rhotot(i)=rhotot(i)+rho(i,2)
     enddo
  endif
  call hartree(0,2,mesh,r,r2,sqr,dx,rhotot,vh)
  !
  ! add exchange and correlation potential: LDA or LSDA only
  !
  rhc=0.0_DP
  do i=1,mesh
     vh(i) = e2*vh(i)
     do is=1,nspin
        rh(is) = rho(i,is)/r2(i)/fpi
     enddo
     if (nlcc) rhc = rhoc(i)/r2(i)/fpi
     call vxc_t(rh,rhc,lsd,vxcp)
     do is=1,nspin
        vnew(i,is)= - zed*e2/r(i)+vxt(i)+vh(i)+vxcp(is)
     enddo
  end do

  deallocate(rhotot)
  !
  ! add exchange and correlation potential: GGA only
  !
  if (gga) then
     allocate(vgc(ndm,2),stat=ierr)
     allocate(egc(ndm),stat=ierr)
     call errore('new_potential','allocating vgc and egc',ierr)

     call vxcgc(ndm,mesh,nspin,r,r2,rho,rhoc,vgc,egc)
     do is=1,nspin
        do i=1,mesh
           vnew(i,is)=vnew(i,is)+vgc(i,is)
        enddo
     enddo
     deallocate(egc)
     deallocate(vgc)
  end if

  !
  ! add OEP exchange 
  !
  if (oep) then
!     write (*,*) ndm, nwf
     allocate(dchi0(ndm,nwf))

     do nu=1,nwf
        call dvex(nu,dchi0(1,nu))
     end do

     call dfx_new(dchi0, vx)
     do is=1,nspin
        do i=1,mesh
! ADD OEP VX
           vnew(i,is)= vnew(i,is)  + vx(i,is)
        end do
     end do
     deallocate(dchi0)
  end if 
  !
  ! latter correction
  !
  if (latt.ne.0) then
     do is=1,nspin
        do i=1,mesh
           vnew(i,is)= min(vnew(i,is),-e2*(zed-enne+1.0_DP)/r(i))
        enddo
     enddo
  end if

  return
end subroutine new_potential
