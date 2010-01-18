!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine elsd_highv (nc)
  !---------------------------------------------------------------
  !
  !   additional information on the atomic total energy. The valence
  !   and core contributions to the different terms are calculated
  !
  use kinds, only : DP
  use constants, only : e2
  use radial_grids, only: ndmx, hartree
  use ld1inc, only: grid, aeccharge, aevcharge, nwf, nspin, enl, oc, v0, &
                    vxcts, excts, excggats, nlcc, enclc, enclv, ehrtvv, &
                    ehrtcv, ehrtcc, ekinc, ekinv, ecxc, ae_fc_energy, &
                    core_state, ekinc0, etot, frozen_core, iswitch
  implicit none
  integer, intent(in) :: nc
  real(DP),allocatable :: f2vv(:), f2cv(:), f2vc(:), f2cc(:), f1c(:), f1v(:),  &
                          f5c(:), f5v(:), vhval(:), vhcore(:), vnew(:,:)

  real(DP) :: int_0_inf_dr, rhotot, fact
  integer:: i,n,ierr
!
!  when iswitch=1 this routine cannot work because the code does not
!  know which are the core and valence states.
!
  if (iswitch==1) return

  allocate(f1c(grid%mesh),stat=ierr)
  allocate(f1v(grid%mesh),stat=ierr)

  allocate(f2vv(grid%mesh),stat=ierr)
  allocate(f2cv(grid%mesh),stat=ierr)
  allocate(f2vc(grid%mesh),stat=ierr)
  allocate(f2cc(grid%mesh),stat=ierr)

  allocate(vnew(ndmx,2),stat=ierr)
  allocate(vhval(ndmx),stat=ierr)
  allocate(vhcore(ndmx),stat=ierr)

  allocate(f5c(grid%mesh),stat=ierr)
  allocate(f5v(grid%mesh),stat=ierr)

  call set_rc_rv()
  call v_of_rho_at (aevcharge,aeccharge,vhval,vxcts,excts,excggats, &
                   vnew,.true.,1)
  call hartree(0,2,grid%mesh,grid,aeccharge,vhcore)
  vhcore=e2*vhcore

  fact=1.0_DP/nspin
  do i=1,grid%mesh
     rhotot=aevcharge(i,1)
     if (nspin==2) rhotot=rhotot+aevcharge(i,2) 
!
!   The integral for the energy due to the interaction with the nuclei
!
     f1c(i)= v0(i) * aeccharge(i)
     f1v(i)= v0(i) * rhotot
!
!   The integrals for the Hartree energy
!
     f2vv(i)= vhval (i) * rhotot
     f2cv(i)= vhcore(i) * rhotot
     f2vc(i)= vhval(i) *  aeccharge(i)
     f2cc(i)= vhcore(i) * aeccharge(i)  
!
!   The integral to be subtracted to the sum of the eigenvalues to
!   get the kinetic energy
!
     f5c(i) =-vxcts(i,1)*aeccharge(i)*fact-f1c(i)-f2vc(i)-f2cc(i)
     if (nspin==2) f5c(i) =f5c(i)-vxcts(i,2)*aeccharge(i)*fact
     f5v(i) =-vxcts(i,1)*aevcharge(i,1)-f1v(i)-f2cv(i)-f2vv(i)
     if (nspin==2) f5v(i) =f5v(i)-vxcts(i,2)*aevcharge(i,2)
  enddo
!
!  Now compute the integrals
!
  enclc=      int_0_inf_dr(f1c,grid,grid%mesh,1)
  enclv=      int_0_inf_dr(f1v,grid,grid%mesh,1)
  ehrtvv=0.5_DP*int_0_inf_dr(f2vv,grid,grid%mesh,2)
  ehrtcc=0.5_DP*int_0_inf_dr(f2cc,grid,grid%mesh,2)
  ehrtcv=     int_0_inf_dr(f2cv,grid,grid%mesh,2)

  ekinc=      int_0_inf_dr(f5c,grid,grid%mesh,1)
  ekinv=      int_0_inf_dr(f5v,grid,grid%mesh,1)
  do n=1,nwf
     if (oc(n)>0.0_DP.and.core_state(n)) ekinc=ekinc+oc(n)*enl(n)
     if (oc(n)>0.0_DP.and..not.core_state(n)) ekinv=ekinv+oc(n)*enl(n)
  enddo
  if (nc==1) ekinc0=ekinc
  if (frozen_core.and.nc>1) then
     etot=etot-ekinc+ekinc0
     ekinc=ekinc0
  endif

  ae_fc_energy=ekinv+ehrtvv+ehrtcv+ecxc+enclv

  deallocate(f5c)
  deallocate(f5v)
  deallocate(f2vc)
  deallocate(f2cv)
  deallocate(f2cc)
  deallocate(f2vv)
  deallocate(f1c)
  deallocate(f1v)
  deallocate(vnew)
  deallocate(vhval)
  deallocate(vhcore)

  return
end subroutine elsd_highv
