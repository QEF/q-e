!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine scatter_forw(zin, zfin, orbin, orbfin)
!
! This subroutine computes local Phi_n and partial nonlocal Phi_alpha
! solutions of the Schrodinger equation in the region zin<z<zfin
! on which the general solution may be constructed:
!        \sum a_n Phi_n + \sum a_alpha Phi_alpha
! It computes also the integrals (intw1, intw2) of  Phi_n and 
! Phi_alpha over beta-functions inside the unit cell. 
!
#include "machine.h"
  use pwcom
  use para
  use cond
implicit none

  integer ::      &
        orbin,    & ! starting orbital         in zin<z<zfin 
        orbfin,   & ! final orbital            in zin<z<zfin
        norbnow,  & ! total number of orbitals in zin<z<zfin
        kin,      & ! first slab               in zin<z<zfin 
        kfin,     & ! last slab                in zin<z<zfin 
        startk,   & ! first slab for a given CPU 
        lastk,    & ! last slab  for a given CPU 
        k, kz, n, lam, ig, lam1, mdim, itt, nbb, iorb, iorb1,   &
        iorbs, iorb1s
  integer :: info
  integer, allocatable :: inslab(:)
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: DDOT, zin, zfin, dz, tr, tr1, dz1 
  complex(kind=DP), parameter :: cim=(0.d0,1.d0) 
  complex(kind=DP) :: int1d, int2d, c, d, e, f, s1, s2, s3, s4, arg,&
                      f1p, ZDOTC
  complex(kind=DP), allocatable ::   &
     psigper(:,:), & ! psigper(g,lam)=newbg(g,lam1) psiper(lam1,lam)
     w0(:,:,:),  &   ! w0(z,g,m) are 2D Fourier components (see four.f)
     w(:,:,:),   &   ! w(z,lam,m)=psigper(g,lam)^* \exp{-igr^m_perp} 
                     !                            w0(z,g,m) 
     ci(:,:,:),  &   ! ci(m,lam,k)=\int_{z(k)}^{z(k+1)} dz 
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z-z(k))}
     di(:,:,:),  &   ! di(m,lam,k)=\int_{z(k)}^{z(k+1)} dz
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z(k+1)-z)}   
     bf(:,:), an(:,:), bn(:,:),        & 
     app(:,:), bpp(:,:), al(:,:),      & 
     bl(:,:), af(:,:),                 &
     bnlf(:,:), anln(:,:), bnln(:,:),  &
     anlp(:,:), bnlp(:,:), anll(:,:),  &
     ff(:,:), fl(:,:)                
                               
  norbnow = orbfin - orbin + 1
  orbin = orbin - 1
!
! Find first and last slab (kin, kfin)
!
  do k=1, nrz
    if (z(k).le.zin+eps)  kin=k
    if (z(k).le.zfin-eps) kfin=k
  enddo
  dz=z(kin+1)-z(kin) 
  dz1=dz/nz1
!
! Divide the slabs among CPU
!
  call divide(kfin-kin+1,startk,lastk)
  startk=kin+startk-1
  lastk=kin+lastk-1

!------------------------
! Start of 2D Fourier components calculations and depending
! variables
!
  allocate( psigper( ngper, n2d ) )
  allocate( w( nz1, n2d, norbnow ) )
  allocate( w0( nz1, ngper, 5 ) )
  allocate( ci( norbnow, n2d, nrzp ) )
  allocate( di( norbnow, n2d, nrzp ) )     
  allocate( inslab( norbnow ) )
  intw1=(0.d0,0.d0)
  intw2=(0.d0,0.d0)
!
! some orbitals relations
!                    
  do iorb=1, norbnow
    inslab(iorb)=0
  enddo    
  do iorb=1, norbnow
    iorbs=orbin+iorb 
    if(inslab(iorb).eq.0.and.mnew(iorbs).eq.1) then
      itt=itnew(iorbs)
      nbb=nbnew(iorbs)
      do iorb1=iorb, norbnow
       iorb1s=orbin+iorb1
       if(mnew(iorb1s).eq.1.and.nbnew(iorb1s).eq.nbb) then
         tr=abs(taunew(3,iorbs)-taunew(3,iorb1s))
         if(itnew(iorb1s).eq.itt.and.tr.le.eps) then
           inslab(iorb1)=iorb
         endif
       endif
      enddo
    endif
  enddo 
!
! The loop over slabs to compute ci, di, and initial intw2
!
  do k=startk, lastk

    c=(1.d0,0.d0)
    d=(0.d0,0.d0)  
    call ZGEMM('n', 'n', ngper, n2d, n2d, c, newbg, ngper,  &
               psiper(1,1,nkofz(k)), n2d, d, psigper, ngper)
    w=(0.d0,0.d0)
    do iorb=1, norbnow
      iorbs=orbin+iorb
      if(cross(iorbs,k).eq.1.and.inslab(iorb).eq.iorb) then  
       mdim=2*ls(iorbs)+1
       call four(iorbs, w0, k, dz)
       do iorb1=1, norbnow
        iorb1s=orbin+iorb1
        if(inslab(iorb1).eq.iorb) then
         do ig=1, ngper
          tr=-tpi*DDOT(2,gper(1,ig),1,taunew(1,iorb1s),1)
          c=cmplx(cos(tr),sin(tr))
          do lam=1, n2d
           d=CONJG(psigper(ig,lam))*c
           do n=1, mdim
            do kz=1, nz1
             w(kz,lam,iorb1-1+n)=w(kz,lam,iorb1-1+n)+d*w0(kz,ig,n)
            enddo
           enddo
          enddo
         enddo                
        endif
       enddo
      endif 
    enddo
    do iorb=1, norbnow
     iorbs = orbin + iorb
     if (cross(iorbs,k).eq.1) then   
        do lam=1, n2d
           ci(iorb,lam,nkofz(k))=int1d(w(1,lam,iorb),   &
                        zk(lam,nkofz(k)),dz,dz1,nz1,tpiba,1)
           di(iorb,lam,nkofz(k))=int1d(w(1,lam,iorb),   &
                        zk(lam,nkofz(k)),dz,dz1,nz1,tpiba,-1)
        enddo                                                        
        do iorb1=1, norbnow
           iorb1s = orbin + iorb1
           if (cross(iorb1s,k).eq.1) then
             do lam=1, n2d
                intw2(iorb,iorb1)=intw2(iorb,iorb1)-       &
                   cim*int2d(w(1,lam,iorb),w(1,lam,iorb1), & 
                         zk(lam,nkofz(k)),dz1,nz1,tpiba)/  &
                       (2.d0*zk(lam,nkofz(k))*tpiba) 
             enddo
           endif 
        enddo
     endif
    enddo      
  enddo 

  deallocate(psigper)
  deallocate(w)
  deallocate(w0)
  deallocate(inslab)

!-----------------------------------

!
! Some allocation for iterative process
!
  allocate( bf( n2d, n2d ) )
  allocate( an( n2d, n2d ) )
  allocate( bn( n2d, n2d ) )
  allocate( app( n2d, n2d ) )
  allocate( bpp( n2d ,n2d ) )
  allocate( al( n2d, n2d ) )
  allocate( bl( n2d, n2d ) )
  allocate( af( n2d, n2d ) )
  allocate( bnlf( n2d, norbnow ) )
  allocate( anln( n2d, norbnow ) )
  allocate( bnln( n2d, norbnow ) )
  allocate( anlp( n2d, norbnow ) )
  allocate( bnlp( n2d, norbnow ) )
  allocate( anll( n2d, norbnow ) )
  allocate( ff( n2d, norbnow ) )
  allocate( fl( n2d, norbnow ) )

!
!    We set up the starting values
!                                       
  app=(0.d0,0.d0)
  an=(0.d0,0.d0)
  bpp=(0.d0,0.d0)
  bn=(0.d0,0.d0)
  bf=(0.d0,0.d0)
  anlp=(0.d0,0.d0)
  anln=(0.d0,0.d0)
  bnlp=(0.d0,0.d0)
  bnln=(0.d0,0.d0)
  bnlf=(0.d0,0.d0)
  ff=(0.d0,0.d0)
  fl=(0.d0,0.d0)

  do lam=1, n2d
    bf(lam,lam)=(1.d0,0.d0)
    bpp(lam,lam)=(1.d0,0.d0)
  enddo                                               
!
! To compute intw1, ff, fl for the first slab
!
  kz=nkofz(startk)
  do iorb=1, norbnow
    iorbs = orbin + iorb
    if (cross(iorbs, startk).eq.1) then
       do lam=1, n2d
          intw1(iorb,lam)=di(iorb, lam, kz)
          arg=tpi*zk(lam, kz)*dz
          if (abs(DIMAG(zk(lam, kz))).lt.eps) then
           ff(lam,iorb)=-cim*EXP(cim*arg)*CONJG(di(iorb,lam,kz))/ &
                             (2.d0*zk(lam,kz)*tpiba) 
           fl(lam,iorb)=-cim*EXP(cim*arg)*CONJG(ci(iorb,lam,kz))/ &
                             (2.d0*zk(lam,kz)*tpiba)
          else
           ff(lam,iorb)=-cim*CONJG(ci(iorb,lam,kz))/              &
                             (2.d0*zk(lam,kz)*tpiba)
           fl(lam,iorb)=-cim*CONJG(di(iorb,lam,kz))/              &
                             (2.d0*zk(lam,kz)*tpiba)
          endif
       enddo
    endif
  enddo

!------------------------------------
! The main loop over slabs 
!
  do k=startk+1, lastk
    kz=nkofz(k) 
    do lam=1, n2d
       do lam1=1,n2d
          c=ZDOTC(n2d,psiper(1,lam,kz),1,psiper(1,lam1,kz-1),1)
          s1=(zk(lam,kz)+zk(lam1,kz-1))/zk(lam,kz)*c
          s2=(zk(lam,kz)-zk(lam1,kz-1))/zk(lam,kz)*c
          c=EXP(tpi*cim*zk(lam1,kz-1)*dz)
          s3=s1*c
          s4=s2*c
          do n=1, n2d
             an(lam,n)= an(lam,n)+s3*app(lam1,n)
             bn(lam,n)= bn(lam,n)+s4*app(lam1,n)   
          enddo
          an(lam,lam1)= an(lam,lam1)+s2
          bn(lam,lam1)= bn(lam,lam1)+s1
          do iorb=1, norbnow
            iorbs = orbin + iorb
            tr=taunew(3,iorbs)-rsph(nbnew(iorbs),itnew(iorbs))
            if (z(k)+dz.gt.tr) then
               anln(lam,iorb)=anln(lam,iorb)+s3*anlp(lam1,iorb)+ &
                       s1*fl(lam1,iorb)      
               bnln(lam,iorb)=bnln(lam,iorb)+s4*anlp(lam1,iorb)+ &
                       s2*fl(lam1,iorb)      
            endif
          enddo
       enddo 
       do n=1, n2d
          an(lam,n)=an(lam,n)*0.5d0
          bn(lam,n)=bn(lam,n)*0.5d0*EXP(-tpi*cim*zk(lam,kz)*dz)
       enddo
       do iorb=1, norbnow
          iorbs = orbin + iorb
          tr=taunew(3,iorbs)-rsph(nbnew(iorbs),itnew(iorbs))
          if (z(k)+dz.gt.tr) then 
             anln(lam,iorb)=anln(lam,iorb)*0.5d0
             bnln(lam,iorb)=bnln(lam,iorb)*0.5d0                &
                                    *EXP(-tpi*cim*zk(lam,kz)*dz)
          endif
       enddo
    enddo
    call DCOPY(2*n2d*n2d, an, 1, app, 1)
    call DCOPY(2*n2d*n2d, bn, 1, bpp, 1)
    an=(0.d0,0.d0)
    bn=(0.d0,0.d0)
    call DCOPY(2*norbnow*n2d, anln, 1, anlp, 1)
    call DCOPY(2*norbnow*n2d, bnln, 1, bnlp, 1)
    anln=(0.d0,0.d0)
    bnln=(0.d0,0.d0)
    fl=(0.d0,0.d0)
!
    do iorb=1, norbnow
       iorbs = orbin + iorb 
       if(cross(iorbs, k).eq.1) then
          do lam=1, n2d
             arg=tpi*zk(lam,kz)*dz
             do n=1, n2d   
              intw1(iorb,n)=intw1(iorb,n)+app(lam,n)*ci(iorb,lam,kz)+ &
                                    bpp(lam,n)*di(iorb,lam,kz)       
             enddo
             if (abs(DIMAG(zk(lam,kz))).lt.eps) then
              f1p=-cim*EXP(cim*arg)*CONJG(di(iorb,lam,kz))/        &
                          (2.d0*zk(lam,kz)*tpiba) 
              fl(lam,iorb)=-cim*EXP(cim*arg)*CONJG(ci(iorb,lam,kz))/ &
                          (2.d0*zk(lam,kz)*tpiba) 
             else
              f1p=-cim*CONJG(ci(iorb,lam,kz))/(2.d0*zk(lam,kz)*tpiba)
              fl(lam,iorb)=-cim*CONJG(di(iorb,lam,kz))/              &
                                      (2.d0*zk(lam,kz)*tpiba) 
             endif
             bnlp(lam,iorb)=bnlp(lam,iorb)-f1p*EXP(-cim*arg)
          enddo
       endif   
    enddo
    do iorb=1, norbnow
     iorbs = orbin + iorb
     if(cross(iorbs, k).eq.1) then 
          do iorb1=1, norbnow
             iorb1s = orbin + iorb1
             tr=taunew(3,iorb1s)-rsph(nbnew(iorb1s),itnew(iorb1s))
             if (z(k)+dz.gt.tr) then 
               c=(0.d0, 0.d0)
               do lam=1, n2d
                  c=c+anlp(lam,iorb1)*ci(iorb,lam,kz)+              &
                      bnlp(lam,iorb1)*di(iorb,lam,kz)   
               enddo
               intw2(iorb,iorb1)=intw2(iorb,iorb1)+c
             endif
          enddo
       endif
    enddo
!
! Rotation of linear solutions
!
    call rotatef(app, bpp, bf, anlp, bnlp, bnlf, intw1, intw2,     &
                 n2d, norbf, norbnow)
  enddo
!---------------------------------------------

  call DCOPY(2*n2d*n2d, app, 1, al, 1)
!
! To compute the 2nd half of linear solutions
!
  call scatter_back(app, bpp, an, bn, af, ci, di, startk, lastk, &
                    norbnow, orbin, dz)

  call DCOPY(2*n2d*n2d, bpp, 1, bl, 1)
  call DCOPY(2*n2d*norbnow, anlp, 1, anll, 1)
  call DSCAL(2*norbf*2*n2d, sarea, intw1, 1)
  call DSCAL(2*norbf*norbf, sarea, intw2, 1)

!
! To construct functions and derivetives on the boundaries
!
! local solutions 
  fun0=(0.d0,0.d0)
  fun1=(0.d0,0.d0)
  fund0=(0.d0,0.d0)
  fund1=(0.d0,0.d0)
  k=nkofz(startk)
  kz=nkofz(lastk)
  do n=1, n2d
    do lam=1, n2d
      s1=bf(lam,n)*EXP(cim*tpi*zk(lam,k)*dz)
      s2=al(lam,n)*EXP(cim*tpi*zk(lam,kz)*dz)
      if (lam.eq.n) s2=s2+(1.d0,0.d0)
      s3=-cim*zk(lam,k)*s1*tpiba
      s4=cim*zk(lam,kz)*s2*tpiba
      if (lam.eq.n) s4=s4-2.d0*cim*zk(lam,kz)*tpiba
      c=bl(lam,n)*EXP(cim*tpi*zk(lam,k)*dz)
      if (lam.eq.n) c=c+(1.d0,0.d0)
      d=af(lam,n)*EXP(cim*tpi*zk(lam,kz)*dz)
      e=-cim*zk(lam,k)*c*tpiba
      if (lam.eq.n) e=e+2.d0*cim*zk(lam,k)*tpiba
      f=cim*zk(lam,kz)*d*tpiba
      do ig=1, n2d
        fun0(ig,n)=fun0(ig,n)+psiper(ig,lam,k)*s1
        fun1(ig,n)=fun1(ig,n)+psiper(ig,lam,kz)*s2
        fund0(ig,n)=fund0(ig,n)+psiper(ig,lam,k)*s3
        fund1(ig,n)=fund1(ig,n)+psiper(ig,lam,kz)*s4
        fun0(ig,n+n2d)=fun0(ig,n+n2d)+psiper(ig,lam,k)*c
        fun1(ig,n+n2d)=fun1(ig,n+n2d)+psiper(ig,lam,kz)*d
        fund0(ig,n+n2d)=fund0(ig,n+n2d)+psiper(ig,lam,k)*e
        fund1(ig,n+n2d)=fund1(ig,n+n2d)+psiper(ig,lam,kz)*f
       enddo
     enddo
  enddo                                                               
! nonlocal solutions
  funl0=(0.d0,0.d0)
  funl1=(0.d0,0.d0)
  fundl0=(0.d0,0.d0)
  fundl1=(0.d0,0.d0)
  do iorb=1, norbnow
    do lam=1, n2d
      s1=ff(lam,iorb)+bnlf(lam,iorb)*                   &
                            exp(cim*tpi*zk(lam,k)*dz)
      s2=fl(lam,iorb)+anll(lam,iorb)*                   &
                            exp(cim*tpi*zk(lam,kz)*dz)
      s3=-cim*zk(lam,k)*s1*tpiba
      s4=cim*zk(lam,kz)*s2*tpiba
      do ig=1, n2d
       funl0(ig,iorb)=funl0(ig,iorb)+psiper(ig,lam,k)*s1
       funl1(ig,iorb)=funl1(ig,iorb)+psiper(ig,lam,kz)*s2
       fundl0(ig,iorb)=fundl0(ig,iorb)+psiper(ig,lam,k)*s3
       fundl1(ig,iorb)=fundl1(ig,iorb)+psiper(ig,lam,kz)*s4
      enddo
    enddo
  enddo               

  deallocate(ci)
  deallocate(di)   
  deallocate(bf)
  deallocate(an)
  deallocate(bn)
  deallocate(app)
  deallocate(bpp)
  deallocate(al)
  deallocate(bl)
  deallocate(af)
  deallocate(bnlf)
  deallocate(anln)
  deallocate(bnln)
  deallocate(anlp)
  deallocate(bnlp)
  deallocate(anll)
  deallocate(ff)
  deallocate(fl)

!
! To construct the functions in the whole rigion zin<z<zfin in the
! case of multiparallel running
!
#ifdef __PARA
  call rotproc(fun0, fund0, fun1, fund1, funl0, fundl0, funl1, &
               fundl1, intw1, intw2, n2d, norbf, norbnow)                    
#endif       

  return
end subroutine scatter_forw
