!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine allocate_cond
!
!  This subroutine allocates the variables for PWCOND  
!
#include "f_defs.h"
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau
  use pwcom
  USE uspp_param, ONLY : nbrx, nbeta, lll, betar, tvanp
  use atom, only: mesh, r
  use cond 
  implicit none
  integer :: mmax, nt, k, iw, ib, ir, na, naux 
  real(kind=DP), parameter :: eps=1.d-8 
  real(kind=DP) :: & 
     bmax,         & ! maximum of beta function   
     ledge,        & ! left edge of the orbital   
     redge,        & ! right edge of the orbital  
     dz1, dz2,     & ! auxiliary  
     epsbeta       
  real(kind=DP), allocatable :: dwid(:)

!
! FFT parameters for local potential
!
  nrx=nr1s
  nry=nr2s
  nrz=nr3s+dnslab
  if(nrz/2*2.eq.nrz) nrz=nrz+1

  allocate ( z(nrz+1) )
  allocate ( rsph(nbrx, ntypx) )
  allocate ( dwid(6) )              

!
! slab parameters
!
  zl=at(3,3)
  sarea=abs(at(1,1)*at(2,2)-at(2,1)*at(1,2))*alat**2
  if (abs(ecut2d).le.eps) ecut2d=ecutwfc                   

  if (ikind.gt.2.or.ikind.le.-1) &
       call errore ('init_cond','ikind is not in allowed range',1)   
  if (ikind.eq.0) then
     bdr1=bdl2
     bdr2=bdl2
  endif     
  if (ikind.eq.1) then
    bdr1 = zl
    bdr2 = zl
  endif

  dz1=zl/nrz
  dwid(1)=0.d0
  dwid(2)=bdl1
  dwid(3)=bdl2
  dwid(4)=bdr1
  dwid(5)=bdr2
  dwid(6)=zl
  nrz=0
  do iw=2, 6
     naux=(dwid(iw)-dwid(1)+dz1*0.5d0)/dz1-nrz
     if (naux.gt.0) then
       dz2=(dwid(iw)-dwid(iw-1))/naux
       do k=1, naux
         z(nrz+k)=dwid(iw-1)+(k-1)*dz2
       enddo
       nrz=nrz+naux
     endif
  enddo
  z(nrz+1)=zl   

!
! to determine radii of nonlocal spheres
!
  epsbeta=1.d-4
  mmax = 0
  do nt=1, ntyp
    do ib=1, nbeta(nt)
      mmax = max(mmax, lll(ib, nt))
      bmax=0.d0
      do ir=2, mesh(nt)
         bmax=max(bmax, betar(ir,ib,nt)/r(ir,nt)) 
      enddo  
      ir=mesh(nt)
      do while (abs(betar(ir,ib,nt)/r(ir,nt)).le.epsbeta*bmax)
        ir=ir-1
      enddo
      rsph(ib,nt)=r(ir,nt)/alat
    enddo
  enddo

  if (mmax.gt.3) call errore ('allocate','for l>3 -orbitals  &
        &                  the method is not implemented',1)
!
! We set up the radii of US spheres to be the same (to avoid 
! the problem with the spheres crossing or not the boundaries)
!
  do nt=1, ntyp
   if (tvanp(nt)) then
      bmax=0.d0
      do ib=1, nbeta(nt)
         bmax=max(bmax, rsph(ib,nt))
      enddo
      do ib=1, nbeta(nt)
         rsph(ib,nt)=bmax
      enddo  
   endif
  enddo              
!
! Move all atoms into the unit cell 
!
  do na=1, nat
    tau(3,na) = tau(3,na)/zl
    tau(3,na) = tau(3,na) - int(tau(3,na))
    if(tau(3,na).le.eps) tau(3,na)=tau(3,na)+1.d0 
    tau(3,na) = tau(3,na)*zl
  enddo

!
! Compute the number of orbitals inside the cell and 
! crossing the boundary
!
  noinsl=0
  nocrosl=0
  noinss=0
  noinsr=0
  nocrosr=0
  norbs=0
!
! for left tip
!    
  do na=1, nat
     nt=ityp(na)
     do ib=1, nbeta(nt)
        ledge=tau(3,na)-rsph(ib,nt)
        redge=tau(3,na)+rsph(ib,nt)
        if (ledge.le.bdl1.and.redge.gt.bdl2) &
            call errore ('init_cond','The unit cell in left',1)
        if (tau(3,na).gt.bdl1+eps.and.tau(3,na).le.bdl2+eps) then
          if (ledge.gt.bdl1.and.redge.le.bdl2) then
             noinsl=noinsl+2*lll(ib,nt)+1
          else
             nocrosl=nocrosl+2*lll(ib,nt)+1
          endif
        endif   
     enddo
  enddo

  norb=2*nocrosl+noinsl
  norbf=norb
!
! If scattering problem
!
  if (ikind.ne.0) then
!
!       Scatt. region
!
    do na=1, nat
       nt=ityp(na)
       do ib=1, nbeta(nt)
          ledge=tau(3,na)-rsph(ib,nt)
          redge=tau(3,na)+rsph(ib,nt)
          if (ledge.le.bdl2.and.redge.gt.bdr1) &
             call errore ('init_cond','The unit cell in scatt.',1) 
          if (ledge.gt.bdl2.and.redge.le.bdr1) & 
             noinss=noinss+2*lll(ib,nt)+1
       enddo
    enddo
    norb=norb+noinss
!
!       Right tip
!
    if (ikind.eq.1) then
!
!       If the tips are equal
!   
      nocrosr=nocrosl
      noinsr=noinsl
      norb=norb+nocrosr 
!
!       If the tips are different
! 
    else
      do na=1, nat
         nt=ityp(na)
         do ib=1, nbeta(nt)
            ledge=tau(3,na)-rsph(ib,nt)
            redge=tau(3,na)+rsph(ib,nt)
            if (ledge.le.bdr1.and.redge.gt.bdr2) &
               call errore ('init_cond','The unit cell in right',1)
            if (ledge.gt.bdr1.and.redge.le.bdr2) &
               noinsr=noinsr+2*lll(ib,nt)+1
            if (ledge.le.bdr1.and.redge.gt.bdr1) &
               nocrosr=nocrosr+2*lll(ib,nt)+1
         enddo
      enddo
      norb=norb+2*nocrosr+noinsr     
      norbf=max(norbf, 2*nocrosr+noinsr)
    endif

    norbs=nocrosl+noinss+nocrosr
    norbf=max(norbf, norbs)
  endif 

  allocate( vppot(nrz, nrx * nry) )
  allocate( itnew(norb) )
  allocate( nbnew(norb) )
  allocate( natih(norb, 2) )
  allocate( ls(norb) )
  allocate( mnew(norb) )
  allocate( cross(norb, nrz) )
  allocate( taunew(3, norb) )
  allocate( zpseu(norb, norb, nspin) )
  
  allocate( nkofz(nrz) ) 

  deallocate(dwid) 

  return 
end subroutine allocate_cond
