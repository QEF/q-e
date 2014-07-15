!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine elsdps 
  !---------------------------------------------------------------
  !
  !   atomic total energy in the local-spin-density scheme
  !   atomic pseudopotentials with nonlinear core correction are allowed
  !   gradient correction allowed (A. Dal Corso fecit AD 1993)
  !
  use kinds, only: DP
  use constants, only: fpi
  use radial_grids, only : ndmx
  use ld1_parameters, only : nwfsx
  use ld1inc, only : nlcc, grid, nspin, rhoc, rhos, lsd, vpsloc, vxt, vh, &
       encl, ehrt, ecxc, evxt, ekin, ecc, epseu, vnl, &
       etots, pseudotype, phits, ikk, nbeta, betas, bmat, &
       nwfts, rel, jjts, llts, octs, enlts, jjs, lls, &
       vxc, exc, excgga
  use funct, only : dft_is_gradient
  implicit none
  real(DP) :: &
       excc, vxcc(2), &   ! exch-corr energy from core charge
       int_0_inf_dr,  &   ! the integral function
       rh0(2),        &   ! the charge in a given point
       rhc,           &   ! core charge in a given point
       rho_tot,       &   ! the total charge in one point
       work(nwfsx)        ! auxiliary space (similar to becp)

  real(DP),allocatable :: &
       f1(:),      &   ! auxiliary
       f2(:),      &   ! auxiliary
       f3(:),      &   ! auxiliary
       f4(:),      &   ! auxiliary
       f5(:),      &   ! auxiliary
       vgc(:,:),   &   ! the gga potential
       egc(:),     &   ! the gga energy
       rho_aux(:,:), & ! auxiliary space
       exccc(:)        ! the exchange and correlation energy of the core
  REAL(dp) :: & ! compatibility with metaGGA - not yet used
       tau(ndmx) = 0.0_dp, vtau(ndmx) = 0.0_dp

  integer :: &
       n,i,ns,nst,lam,n1,n2,ikl,ierr,ind

  allocate(f1(grid%mesh), stat=ierr)
  allocate(f2(grid%mesh), stat=ierr)
  allocate(f3(grid%mesh), stat=ierr)
  allocate(f4(grid%mesh), stat=ierr)
  allocate(f5(grid%mesh), stat=ierr)
  allocate(exccc(ndmx), stat=ierr)
  !
  !  If there is NLCC we calculate here also the exchange and correlation
  !  energy of the pseudo core charge.
  !  This quantity is printed but not added to the total energy
  !
  exccc=0.0_DP
  ecc=0.0_DP
  if (nlcc) then
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
  endif
  !
  !  Now prepare the integrals
  !
  do i=1,grid%mesh
     rho_tot=rhos(i,1)
     if (lsd.eq.1) rho_tot=rho_tot+rhos(i,2)
     !
     !    The integral for the interaction with the local potential
     !
     f1(i)= vpsloc(i) * rho_tot
     !
     !    The integral for the Hartree energy
     !
     f2(i)= vh(i) * rho_tot
     !
     !    The integral for the exchange and correlation energy
     !
     f3(i)= exc(i) * (rho_tot+rhoc(i)) + excgga(i)
     !
     !    The integral for the interaction with the external potential
     !
     f4(i)= vxt(i)*rho_tot
     !
     !    The integral to add to the sum of the eigenvalues to have the
     !    kinetic energy.
     !
     f5(i) =-vxc(i,1)*rhos(i,1)-f1(i)-f2(i)-f4(i)
     if (nspin==2) f5(i)=f5(i)-vxc(i,2)*rhos(i,2)
  enddo
  !
  !  And now compute the integrals
  !
  encl=       int_0_inf_dr(f1,grid,grid%mesh,1)
  ehrt=0.5_DP*int_0_inf_dr(f2,grid,grid%mesh,2)
  ecxc=       int_0_inf_dr(f3,grid,grid%mesh,2)
  evxt=       int_0_inf_dr(f4,grid,grid%mesh,2)
  !
  !  Now compute the nonlocal pseudopotential energy. There are two cases:
  !  The potential in semilocal form or in fully separable form
  !
  epseu=0.0_DP
  if (pseudotype == 1) then
     !
     !   Semilocal form
     !
     do ns=1,nwfts
        if (octs(ns)>0.0_DP) then
           if ( rel < 2 .or. llts(ns) == 0 .or. &
                abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
              ind=1
           else if ( rel == 2 .and. llts(ns) > 0 .and. &
                abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
              ind=2
           endif
           f1=0.0_DP
           lam=llts(ns)
           do n=1, grid%mesh
              f1(n) = f1(n) + phits(n,ns)**2 * octs(ns)
           enddo
           do n=1,grid%mesh
              f1(n) = f1(n) * vnl(n,lam,ind)
           end do
           if (ikk(ns) > 0) &
                epseu = epseu + int_0_inf_dr(f1,grid,ikk(ns),2*(lam+1))
        endif
     enddo
  else
     !
     !  Fully separable form
     !
     do ns=1,nwfts
        if (octs(ns).gt.0.0_DP) then
           do n1=1,nbeta
              if ( llts(ns).eq.lls(n1).and.   &
                   abs(jjts(ns)-jjs(n1)).lt.1.e-7_DP) then
                 nst=(llts(ns)+1)*2
                 ikl=ikk(n1)
                 do n=1,ikl
                    f1(n)=betas(n,n1)*phits(n,ns)
                 enddo
                 work(n1)=int_0_inf_dr(f1,grid,ikl,nst)
              else
                 work(n1)=0.0_DP
              endif
           enddo
           do n1=1,nbeta
              do n2=1,nbeta
                 epseu=epseu  &
                      + bmat(n1,n2)*work(n1)*work(n2)*octs(ns)
              enddo
           enddo
        endif
     enddo
  endif
  !
  !  Now compute the kinetic energy
  !
  ekin = int_0_inf_dr(f5,grid,grid%mesh,2) - epseu
  do ns=1,nwfts
     if (octs(ns).gt.0.0_DP) then
        ekin=ekin+octs(ns)*enlts(ns)
     endif
  end do
  !
  !  And the total energy
  !
  etots= ekin + encl + epseu + ehrt + ecxc + evxt

  deallocate(f5)
  deallocate(f4)
  deallocate(f3)
  deallocate(f2)
  deallocate(f1)
  deallocate(exccc)

  return
end subroutine elsdps
