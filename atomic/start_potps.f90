!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine start_potps
  !---------------------------------------------------------------
  !
  !     This routine computes an initial estimate of the screening
  !     potential
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use io_global, only : stdout
  use ld1inc, only : grid, nspin, lsd, nlcc, latt, enne, rhos, rhoc, &
                     nwfts, llts, jjts, octs, iswts, phits, &
                     vxt, vh, vpstot, vpsloc
  implicit none

  integer :: &
       is,   &   ! counter on spin
       n         ! counter on mesh points

  real(DP) ::    &
       vnew(ndmx,2)    ! the potential
  !
  !    compute an initial estimate of the potential
  !
  call chargeps(rhos,phits,nwfts,llts,jjts,octs,iswts)
  call new_potential(ndmx,grid%mesh,grid,0.0_dp,vxt,lsd,nlcc, &
       latt,enne,rhoc,rhos,vh,vnew,1)

  do is=1,nspin
     do n=1,grid%mesh
        vpstot(n,is)=vpsloc(n)+vnew(n,is)
        !      if (is.eq.1) &
        !           write(stdout,'(3f25.16)') grid%r(n), rhos(n,1),vpstot(n,1)
        !      if (is.eq.2.and.nspin.eq.2)  &
        !   write(stdout,'(3f25.16)') 2.0_dp*rhos(n,1),vpstot(n,1),vpstot(n,2)
     enddo
  enddo
  !
  !    screening the D coefficients
  !
  call newd_at ( )
  !
  return
end subroutine start_potps


subroutine guess_initial_wfc()
!
!  This subroutine guess some initial wavefunction for the test configuration
!
use kinds, only : DP
use radial_grids, only : ndmx
use ld1inc, only: nwfts, octs, llts, jjts, nstoaets, rcutts, rcutusts, grid, &
                  lpaw, psi, lpaw, tm, enlts, elts, phits, pseudotype, iswitch
implicit none

integer ::    &
       ns,n,     & ! counters 
       lam,      & ! angular momentum
       nwf0,     & ! all-electron state
       ik,ikus     ! points on the mesh

real(DP) :: psi_in(ndmx), xc(8)

do ns=1,nwfts
   if (octs(ns).gt.0.0_dp) then
      lam=llts(ns)
      nwf0=nstoaets(ns)
      !
      !        compute the ik closer to r_cut
      !
      ik=0
      ikus=0
!      write(6,*) ns, lam, rcutts(ns), rcutusts(ns)
      do n=1,grid%mesh
         if (grid%r(n).lt.rcutts(ns)) ik=n
         if (grid%r(n).lt.rcutusts(ns)) ikus=n
      enddo
      if (mod(ik,2).eq.0) ik=ik+1
      if (mod(ikus,2).eq.0) ikus=ikus+1
      if (ikus.gt.grid%mesh) &
           call errore('starting potential','ik is wrong ',1)
      !
      !    compute the phi functions
      !
      if (lpaw) then
          phits(:,ns)=psi(:,1,nwf0)
      else
         if (tm.or.(pseudotype<3.and.iswitch==2)) then
            call compute_phi_tm(lam,ik,psi(1,1,nwf0),phits(1,ns),0,xc,   &
                                                     enlts(ns),elts(ns))
         else
            call compute_phi(lam,ik,psi(1,1,nwf0),phits(1,ns),xc,0,octs(ns), &
                                 enlts(ns),'  ')
         endif
      endif
      if (pseudotype.eq.3) then
         !
         !   US only on the components where ikus <> ik
         !
         psi_in(:)=phits(:,ns)
         if (ikus.ne.ik.or.lpaw) call compute_phius(lam,ikus,psi_in, &
                                                     phits(1,ns),xc,0,'  ')
      endif
      call normalize(phits(1,ns),llts(ns),jjts(ns), ns)
   else
      phits(:,ns)=0.0_dp
   endif
enddo

return
end subroutine guess_initial_wfc
