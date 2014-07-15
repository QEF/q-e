!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine elsdps_paw( )
  !---------------------------------------------------------------
  !
  !   total paw energy in the local-spin-density scheme
  !
  use kinds, only: DP
  use constants, only: fpi
  use radial_grids, only : ndmx
  use ld1_parameters, only : nwfsx
  use ld1inc, only : nlcc, grid, nspin, rhoc, lsd, &
       encl, ehrt, ecxc, evxt, ekin, ecc, epseu,  &
       nwfts, enlts, octs, paw_energy
  use funct, only : dft_is_gradient
  implicit none
  real(DP) :: &
       excc, vxcc(2), &   ! exch-corr energy from core charge
       int_0_inf_dr,  &   ! the integral function
       rh0(2),        &   ! the charge in a given point
       rhc,           &   ! core charge in a given point
       edcts              ! auxiliary energy

  real(DP),allocatable :: &
       vgc(:,:),   &   ! the gga potential
       egc(:),     &   ! the gga energy
       rho_aux(:,:), & ! auxiliary space
       exccc(:)        ! the exchange and correlation energy of the core
  REAL(dp) :: & ! compatibility with metaGGA - not yet used
       tau(ndmx) = 0.0_dp, vtau(ndmx) = 0.0_dp

  integer :: &
       i,ns,ierr

  !
  !  If there is NLCC we calculate here also the exchange and correlation
  !  energy of the pseudo core charge.
  !  This quantity is printed but not added to the total energy
  !
  ecc=0.0_DP
  if (nlcc) then
     allocate(exccc(ndmx), stat=ierr)
     exccc=0.0_DP
     rh0(1)=0.0_DP
     rh0(2)=0.0_DP
     do i=1,grid%mesh
        rhc= rhoc(i)/grid%r2(i)/fpi
        call vxc_t(lsd,rh0,rhc,excc,vxcc)
        exccc(i) = excc*rhoc(i) 
     enddo
     if (dft_is_gradient()) then
        allocate(rho_aux(ndmx,2), stat=ierr)
        allocate(vgc(ndmx,2),stat=ierr)
        allocate(egc(ndmx),stat=ierr)
        vgc=0.0_DP
        egc=0.0_DP
        rho_aux=0.0_DP
        call vxcgc ( ndmx, grid%mesh, nspin, grid%r, grid%r2, rho_aux, &
             rhoc, vgc, egc, tau, vtau, 1)
        do i=1,grid%mesh
           exccc(i) = exccc(i) + egc(i)*fpi*grid%r2(i)
        enddo
        deallocate(egc)
        deallocate(vgc)
        deallocate(rho_aux)
     endif
     ecc=  int_0_inf_dr(exccc,grid,grid%mesh,2)
     deallocate(exccc)
  endif
  !
  !  Add the three contributions for each energy
  !
  encl= paw_energy(5,1)+paw_energy(5,2)-paw_energy(5,3)
  ehrt= paw_energy(2,1)+paw_energy(2,2)-paw_energy(2,3)
  ecxc= paw_energy(3,1)+paw_energy(3,2)-paw_energy(3,3)
  edcts=  paw_energy(4,1)+paw_energy(4,2)-paw_energy(4,3)
  !
  !  The nonlocal pseudopotential energy is not computed.
  !
  epseu=0.0_DP
  !
  !  Now compute the kinetic energy
  !
  ekin = -encl-edcts
  do ns=1,nwfts
     if (octs(ns) > 0.0_DP) then
        ekin=ekin+octs(ns)*enlts(ns)
     endif
  end do
  return
end subroutine elsdps_paw
