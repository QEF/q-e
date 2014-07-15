!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine sic_correction(n,vhn1,vhn2,egc) 
  !---------------------------------------------------------------
  !   set up the orbital-dependent selfconsistent potential generated
  !   by the n-th wavefunction - for self-interaction correction
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use constants, only: e2, fpi
  use ld1inc, only : nspin, lsd, rel, nlcc, rhoc, grid, psi
  use funct, only: dft_is_gradient
  use radial_grids, only: hartree
  implicit none
  integer :: n
  real(DP):: vhn1(ndmx),vhn2(ndmx), egc(ndmx)
  REAL(dp) :: & ! compatibility with metaGGA - not yet used
       tau(ndmx) = 0.0_dp, vtau(ndmx) = 0.0_dp
  !
  integer :: i
  real(DP):: rh(2), rhc, vxcp(2), excp
  real(DP):: vgc(ndmx,2),  egc0(ndmx), rhotot(ndmx,2)
  logical :: gga

  vhn1=0.0_dp
  vhn2=0.0_dp
  gga=dft_is_gradient()
  nspin=1
  if (lsd.eq.1) nspin=2
  !
  !   compute hartree potential with the charge of orbital n
  !
  rhotot=0.0_dp
  if (rel.eq.2) then
     do i=1,grid%mesh
        rhotot(i,1)=psi(i,1,n)**2+psi(i,2,n)**2
     enddo
  else
     do i=1,grid%mesh
        rhotot(i,1)=psi(i,1,n)**2
     enddo
  endif
  !call hartree(0,2*(ll(n)+1),grid%mesh,grid,rhotot,vhn1)
  call hartree(0,2,grid%mesh,grid,rhotot,vhn1)
  !
  !    add exchange and correlation potential: LDA or LSDA terms
  !
  rhc=0.0_dp
  rh=0.0_dp
  do i=1,grid%mesh
     vhn1(i) = e2*vhn1(i)
     rh(1) = rhotot(i,1)/grid%r2(i)/fpi
     if (nlcc) rhc = rhoc(i)/grid%r2(i)/fpi
     call vxc_t(lsd,rh,rhc,excp,vxcp)
     vhn2(i)= vhn1(i)+vxcp(1)
     egc(i)= excp*rhotot(i,1)
  end do

  if (.not.gga) return
  !
  !   add gradient-correction terms to exchange-correlation potential
  !
  egc0=egc
  call vxcgc ( ndmx, grid%mesh, nspin, grid%r, grid%r2, rhotot, rhoc, &
       vgc, egc, tau, vtau, 1)
  do i=1,grid%mesh
     vhn2(i)=vhn2(i)+vgc(i,1)
     egc(i)=egc(i)*grid%r2(i)*fpi+egc0(i)
  enddo
  return
end subroutine sic_correction
