!---------------------------------------------------------------
subroutine scf
  !---------------------------------------------------------------
  !
  !   this routine performs the atomic self-consistent procedure
  !   self-interaction-correction allowed
  !
  use constants, only: e2
  use ld1inc
  implicit none

  logical:: conv
  integer:: nerr, nstop, n, i, is, id, nin, mch
  real(kind=dp) ::  vnew(ndm,2), rhoc1(ndm), vhn1(ndm), egc(ndm), &
       rab(ndm), ze2
  real(kind=dp), allocatable ::  vsic(:,:), vsicnew(:)
  integer, parameter :: maxter=200
  real(kind=dp), parameter :: thresh=1.0e-10_dp
  !
  ! 
  ze2 = - zed * e2
  rhoc1=0.0_dp
  id=3
  psi_dir=0.0_dp
  rab=dx*r
  !
  if (isic /= 0) then
     allocate(vsic(ndm,nwf))
     allocate(vsicnew(ndm))
     vsic=0.0_dp
     id=1
  endif
  do iter=1,maxter
     nerr=0
     vnew = vpot
     do n=1,nwf
        if (oc(n) >= 0.0_dp) then
           is=isw(n)
           if (isic /= 0 .and. iter > 1) vnew(:,is)=vsic(:,n)
           if (rel == 0) then
              call ascheq (nn(n),ll(n),enl(n),mesh,dx,r,r2, &
                   sqr,vnew(1,is),ze2,thresh,psi(1,n),nstop)
           elseif (rel == 1) then
              call lschps (1,zed,exp(dx),dx,mesh,nin,mch, &
                   nn(n),ll(n),enl(n),psi(1,n),r,vnew(1,is))
              nstop=0
           elseif (rel == 2) then
              call dirsol (ndm,mesh,nn(n),ll(n),jj(n),iter,enl(n), &
                   thresh,dx,psi_dir(1,1,n),r,rab,vnew(1,is))
              psi(:,n)=psi_dir(:,2,n)
              nstop=0
           else
              call errore('scf','relativistic not programmed',1)
           endif
           !      write(6,*) el(n),enl(n)
           if (nstop /= 0) write(6,'(4i6)') iter,nn(n),ll(n),nstop
           nerr=nerr+nstop
        else
           enl(n)=0.0_dp
           psi(:,n)=0.0_dp
        endif
     enddo

     if (rel==2) then
        call charge(ndm,mesh,nwf,2,oc,psi_dir,rho,isw)
     else
        call charge(ndm,mesh,nwf,1,oc,psi,rho,isw)
     endif

     call new_potential(ndm,mesh,r,r2,sqr,dx,zed,vxt,    &
          lsd,.false.,latt,enne,rhoc1,rho,vh,vnew)

     if (isic /= 0) then
        do n=1,nwf
           if (oc(n) >= 0.0_dp) then
              is=isw(n)
              call sic_correction(n,vhn1,vsicnew,egc)
              vsicnew(:) = vnew(:,is)-vsicnew(:)
              call dmixp(mesh,vsicnew,vsic(1,n),beta,tr2,iter,id,eps0,conv)
           end if
        enddo
     endif
     call vpack(mesh,ndm,nspin,vnew,vpot,1)
     call dmixp(mesh*nspin,vnew,vpot,beta,tr2,iter,id,eps0,conv)
     call vpack(mesh,ndm,nspin,vnew,vpot,-1)
     !   write(6,*) iter, eps0

     if (conv) then
        if (nerr /= 0) call errore('scf','errors in KS equations',-1)
        goto 45
     endif
  enddo
  call errore('scf','warning: convergence not achieved',-1)
45 if (isic /= 0) then
     deallocate(vsicnew)
     deallocate(vsic)
  endif

  return
end subroutine scf
