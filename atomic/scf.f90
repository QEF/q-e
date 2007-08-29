!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine scf(ic)
  !---------------------------------------------------------------
  !
  !   this routine performs the atomic self-consistent procedure
  !   self-interaction-correction allowed
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use constants, only: e2
  use ld1inc, only : grid, zed, psi, isic, vpot, vh, vxt, rho, iter, &
                     lsd, rel, latt, enne, vh, beta, nspin, tr2, eps0, &
                     nwf, nn, ll, jj, enl, oc, isw, core_state, frozen_core
  implicit none

  integer, intent(in) :: ic

  logical:: conv
  integer:: nerr, nstop, n, i, is, id, nin, mch
  real(DP) ::  vnew(ndmx,2), rhoc1(ndmx), ze2
  real(DP), allocatable ::  vsic(:,:), vsicnew(:), vhn1(:), egc(:)
  integer, parameter :: maxter=200
  real(DP), parameter :: thresh=1.0e-10_dp
  !
  ! 
  ze2 = - zed * e2
  rhoc1=0.0_dp
  id=3
  if (.not.frozen_core.or.ic==1) psi=0.0_dp
  !
  if (isic /= 0) then
     allocate(vsic(ndmx,nwf), vsicnew(ndmx), vhn1(ndmx), egc(ndmx))
     vsic=0.0_dp
     ! id=1
  endif
  do iter=1,maxter
     nerr=0
     vnew=vpot
     do n=1,nwf
        if (oc(n) >= 0.0_dp) then
           if (ic==1.or..not.frozen_core.or..not.core_state(n)) then
              is=isw(n)
              if (isic /= 0 .and. iter > 1) vnew(:,is)=vpot(:,is)-vsic(:,n)
              if (rel == 0) then
                 call ascheq (nn(n),ll(n),enl(n),grid%mesh,grid,vnew(1,is),&
                      ze2,thresh,psi(1,1,n),nstop)
              elseif (rel == 1) then
                 call lschps (1,zed,grid,nin,mch,nn(n),ll(n),enl(n),&
                             psi(1,1,n),vnew(1,is))
                 nstop=0
              elseif (rel == 2) then
                 call dirsol (ndmx,grid%mesh,nn(n),ll(n),jj(n),iter,enl(n), &
                      thresh,grid,psi(1,1,n),vnew(1,is))
                 nstop=0
              else
                 call errore('scf','relativistic not programmed',1)
              endif
              !      write(6,*) el(n),enl(n)
              ! if (nstop /= 0) write(6,'(4i6)') iter,nn(n),ll(n),nstop
              nerr=nerr+nstop
           endif
        else
           enl(n)=0.0_dp
           psi(:,:,n)=0.0_dp
        endif
     enddo
     !
     ! calculate charge density (spherical approximation)
     !
     rho=0.0_dp
     do n=1,nwf
        do i=1,grid%mesh
           rho(i,isw(n))=rho(i,isw(n))+oc(n)*(psi(i,1,n)**2+psi(i,2,n)**2)
        enddo
     enddo
     !
     ! calculate new potential
     !
     call new_potential(ndmx,grid%mesh,grid,zed,vxt,&
          lsd,.false.,latt,enne,rhoc1,rho,vh,vnew,0)
     !
     ! calculate SIC correction potential (if present)
     !
     if (isic /= 0) then
        do n=1,nwf
           if (oc(n) >= 0.0_dp) then
              is=isw(n)
              call sic_correction(n,vhn1,vsicnew,egc)
              !
              ! use simple mixing for SIC correction
              !
              vsic(:,n) = (1.0_dp-beta)*vsic(:,n)+beta*vsicnew(:)
           end if
        enddo
     endif
     !
     ! mix old and new potential
     !
     call vpack(grid%mesh,ndmx,nspin,vnew,vpot,1)
     call dmixp(grid%mesh*nspin,vnew,vpot,beta,tr2,iter,id,eps0,conv)
     call vpack(grid%mesh,ndmx,nspin,vnew,vpot,-1)
!        write(6,*) iter, eps0
     !
     if (conv) then
        if (nerr /= 0) call infomsg ('scf','errors in KS equations')
        goto 45
     endif
  enddo
  call infomsg('scf','warning: convergence not achieved')
45 if (isic /= 0) then
     deallocate(egc, vhn1, vsicnew, vsic)
  endif

  return
end subroutine scf
