!
! Copyright (C) 2004-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine new_potential &
     (ndm,mesh,grid,zed,vxt,lsd,nlcc,latt,enne,rhoc,rho,vh,vnew,iflag)
  !---------------------------------------------------------------
  !   set up the selfconsistent atomic potential
  !
  use constants, only: fpi, e2
  use radial_grids, only: radial_grid_type, hartree
  use kinds, only : DP
  use funct, only : get_iexch, dft_is_meta, dft_is_gradient
  use ld1inc, only : nwf, vx, vxc, exc, excgga, tau, vtau
  implicit none
  type(radial_grid_type),intent(in):: grid
  integer, intent(in) :: iflag
  logical :: nlcc, gga, oep, meta
  integer :: ndm,mesh,lsd,latt,i,is,nu, nspin, ierr
  real(DP):: rho(ndm,2),vxcp(2),vnew(ndm,2),vxt(ndm),vh(ndm), rhoc(ndm)
  real(DP):: zed,enne,rh(2),rhc, excp
  real(DP),allocatable:: vgc(:,:), egc(:), rhotot(:)
!  real(DP),allocatable:: vx(:,:)
  real(DP),allocatable:: dchi0(:,:)

  if (mesh.ne.grid%mesh) &
       call errore('new_potential','mesh dimension is not as expected',1)
  gga=dft_is_gradient()
  meta=dft_is_meta()
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
  call hartree(0,2,mesh,grid,rhotot,vh)
  deallocate(rhotot)
  !
  ! add exchange and correlation potential: LDA or LSDA only
  !
  rhc=0.0_DP
  do i=1,mesh
     vh(i) = e2*vh(i)
     do is=1,nspin
        rh(is) = rho(i,is)/grid%r2(i)/fpi
     enddo
     if (nlcc) rhc = rhoc(i)/grid%r2(i)/fpi
     if (meta) then
        !
        ! Workaround: the meta-GGA XC functional already contains the LDA part
        !
        vxcp(:)=0.0_dp
        exc(i) =0.0_dp
     else
        call vxc_t(lsd,rh,rhc,excp,vxcp)
        exc(i)=excp
     endif
     do is=1,nspin
        vxc(i,is)=vxcp(is)
        vnew(i,is)= - zed*e2/grid%r(i)+vxt(i)+vh(i)+vxcp(is)
     enddo
  end do
  !
  ! add exchange and correlation potential: GGA only
  !
  if (gga) then
     allocate(vgc(ndm,2),stat=ierr)
     allocate(egc(ndm),stat=ierr)
     call errore('new_potential','allocating vgc and egc',ierr)

     call vxcgc (ndm, mesh, nspin, grid%r, grid%r2, rho, rhoc, &
          vgc, egc, tau, vtau, iflag)
     do is=1,nspin
        do i=1,mesh
           vxc(i,is)=vxc(i,is)+vgc(i,is)
           vnew(i,is)=vnew(i,is)+vgc(i,is)
           excgga(i) =egc(i)*fpi*grid%r2(i)
        enddo
     enddo
     deallocate(egc)
     deallocate(vgc)
  else
     excgga=0.0_DP
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
           vnew(i,is)= min(vnew(i,is),-e2*(zed-enne+1.0_DP)/grid%r(i))
        enddo
     enddo
  end if

  return
end subroutine new_potential
