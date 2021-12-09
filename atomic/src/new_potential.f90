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
  !   from the density and the KS wavefunctions
  !
  use constants, only: fpi, e2
  use radial_grids, only: radial_grid_type, hartree
  use kinds, only : DP
  use xc_lib, only: xclib_get_id, xclib_dft_is, xclib_dft_is_libxc
  use ld1inc, only : nwf, vx, vxc, exc, excgga, tau, vtau
  use kli, only : compute_kli_potential
  implicit none
  type(radial_grid_type),intent(in):: grid
  integer, intent(in) :: iflag
  logical :: nlcc, gga, oep, meta, kli_
  integer :: ndm,mesh,lsd,latt,i,is,nu, nspin, ierr
  real(DP):: rho(ndm,2),vxcp(2),vnew(ndm,2),vxt(ndm),vh(ndm), rhoc(ndm)
  real(DP):: zed,enne,rh(2),rhc, excp
  real(DP),allocatable:: vgc(:,:), egc(:), rhotot(:)
!  real(DP),allocatable:: vx(:,:)


  real(DP),allocatable:: dchi0(:,:) ! 

  if (mesh.ne.grid%mesh) &
       call errore('new_potential','mesh dimension is not as expected',1)
  gga  = xclib_dft_is('gradient')
  meta = xclib_dft_is('meta')

  oep = xclib_get_id('LDA','EXCH')== 4 .and. (.not.xclib_dft_is_libxc('LDA','EXCH'))
  kli_= xclib_get_id('LDA','EXCH')==10 .and. (.not.xclib_dft_is_libxc('LDA','EXCH'))

  nspin = 1
  if (lsd.eq.1) nspin = 2
  !
  !   compute hartree potential with the total charge
  !
  allocate(rhotot(ndm),stat=ierr)

  ! get total density
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
  rhc = 0.0_DP
  do i= 1, mesh
     vh(i) = e2 * vh(i)
     do is = 1, nspin
        rh(is) = rho(i,is)/grid%r2(i)/fpi
     enddo
     if (nlcc) rhc = rhoc(i)/grid%r2(i)/fpi
     if (meta) then
        !
        ! Workaround: the meta-GGA XC functional already contains the LDA part
        !
        vxcp(:)=0.0_dp
        exc(i) =0.0_dp
        print *, "meta gga"
     else
        vxcp = 0
        call vxc_t(lsd,rh,rhc,excp,vxcp)
        exc(i) = excp
     endif
     
     do is =1, nspin
        vxc(i, is) = vxcp(is)
        vnew(i,is)= - zed*e2/grid%r(i)+vxt(i) + vh(i) + vxcp(is)
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
           vxc(i,is) = vxc(i,is) + vgc(i,is)
           vnew(i,is) = vnew(i,is) + vgc(i,is)
           excgga(i) = egc(i)*fpi*grid%r2(i)
        enddo
     enddo
     deallocate(egc)
     deallocate(vgc)
  else
     excgga = 0.0_DP
  end if
  !
  ! add OEP exchange 
  !
  if (oep) then
!     write (*,*) ndm, nwf
     allocate(dchi0(ndm,nwf))
    
     do nu=1,nwf ! num wave functions
        call dvex(nu,dchi0(1,nu))
     end do

     
     call dfx_new(dchi0, vx)
     ! vx contains the oep term
     ! ADD OEP VX
     vnew(:,1:nspin)= vnew(:,1:nspin)  + vx(:,1:nspin)
    
     deallocate(dchi0)
  end if 

  if(kli_) then
      
      call compute_kli_potential(grid%mesh,vx)
      vnew(:, 1:nspin) = vnew(:, 1:nspin ) + vx(:,1:nspin)
  endif


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
