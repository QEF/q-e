!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine v_of_rho_at (rho,rhoc,vh,vxc,exc,excgga,vnew,nlcc,iflag)
  !---------------------------------------------------------------
  !   set up the selfconsistent atomic potential
  !
  use kinds, only : DP
  use constants, only: fpi, e2
  use radial_grids, only: ndmx, hartree
  use funct, only : get_iexch, dft_is_gradient
  use ld1inc, only : nwf, grid, vx, vxt, lsd, zed, enne, latt, nspin
  implicit none
  integer, intent(in) :: iflag
  real(DP), intent(in):: rho(ndmx,2), rhoc(ndmx) ! valence and core charges
  real(DP), intent(out) :: vh(ndmx), vxc(ndmx,2), exc(ndmx), excgga(ndmx), &
                           vnew(ndmx,2)
  logical, intent(in) :: nlcc
  REAL(dp) :: & ! compatibility with metaGGA - not yet used
       tau(ndmx) = 0.0_dp, vtau(ndmx) = 0.0_dp
  ! Hartree potential, exchange and correlation potential and energy
  ! gga exchange and correlation energy, 
  ! vnew is in output potential vnew = vh+vxc+vxc
  !
  logical :: gga, oep
  integer :: i,is,nu,ierr
  real(DP):: vxcp(2),rh(2),rhc, excp
  real(DP),allocatable:: vgc(:,:), egc(:), rhotot(:)
  real(DP),allocatable:: dchi0(:,:)

  gga=dft_is_gradient()
  oep=get_iexch().eq.4
  !
  !   compute hartree potential with the total charge
  !
  allocate(rhotot(ndmx),stat=ierr)
  rhotot=0.0_DP
  do i=1,grid%mesh
     rhotot(i)=rho(i,1)
     if (nspin==2) rhotot(i)=rhotot(i)+rho(i,2)
  enddo
  call hartree(0,2,grid%mesh,grid,rhotot,vh)
  vh(1:grid%mesh) = e2*vh(1:grid%mesh)
  deallocate(rhotot)
  !
  ! calculate exchange and correlation potential: LDA or LSDA only
  !
  rhc=0.0_DP
  exc=0.0_DP
  vxc=0.0_DP
  do i=1,grid%mesh
     do is=1,nspin
        rh(is) = rho(i,is)/grid%r2(i)/fpi
     end do
     if (nlcc) rhc = rhoc(i)/grid%r2(i)/fpi
     call vxc_t(lsd,rh,rhc,excp,vxcp)
     do is=1,nspin
        vxc(i,is)=vxcp(is)
     end do
     exc(i)=excp
  end do
  !
  ! if gga add gga exchange and correlation potential
  !
  excgga=0.0_DP
  if (gga) then
     allocate(vgc(ndmx,2),stat=ierr)
     call errore('new_potential','allocating vgc',ierr)
     allocate(egc(ndmx),stat=ierr)
     call errore('new_potential','allocating egc',ierr)

     call vxcgc (ndmx, grid%mesh, nspin, grid%r, grid%r2, rho, rhoc, &
          vgc, egc, tau, vtau, iflag)
     do is=1,nspin
        do i=1,grid%mesh
           vxc(i,is)=vxc(i,is)+vgc(i,is)
           excgga(i) =egc(i)*fpi*grid%r2(i)
        enddo
     enddo
     deallocate(egc)
     deallocate(vgc)
  end if
!
!  Calculate the self consistent potential
!
  do is=1,nspin
     do i=1,grid%mesh
        vnew(i,is)=vh(i)+vxc(i,is)
     enddo
  enddo
  !
  ! add OEP exchange 
  !
  if (oep) then
!     write (*,*) ndmx, nwf
     allocate(dchi0(ndmx,nwf))

     do nu=1,nwf
        call dvex(nu,dchi0(1,nu))
     end do

     call dfx_new(dchi0, vx)
     do is=1,nspin
        do i=1,grid%mesh
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
        do i=1,grid%mesh
           vnew(i,is)= min(vnew(i,is),-e2*(1.0_DP-enne)/grid%r(i))
        enddo
     enddo
  end if

  return
end subroutine v_of_rho_at
