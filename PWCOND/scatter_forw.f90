!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Optimized Aug. 2004 (ADC)
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
!
#include "f_defs.h"
!
subroutine scatter_forw(nrz, nrzp, z, psiper, zk, norb, tblm, cros, &
                        taunew, r, rab, betar)
!
! This subroutine computes local Phi_n and partial nonlocal Phi_alpha
! solutions of the Schrodinger equation in the region zin<z<zfin
! on which the general solution may be constructed:
!        \sum a_n Phi_n + \sum a_alpha Phi_alpha
! It computes also the integrals (intw1, intw2) of  Phi_n and 
! Phi_alpha over beta-functions inside the unit cell. 
!
  USE constants, ONLY : tpi
  USE parameters, only : ndmx, nbrx, npsx
  USE cell_base, ONLY : tpiba
  USE noncollin_module, ONLY : npol
  USE cond
  !
  IMPLICIT NONE

  INTEGER ::      &
        nrz, nrzp, norb,& 
        cros(norb,nrz), &
        tblm(4,norb),   &
        k, kz, n, lam, ig, lam1, mdim, itt, nbb, iorb, iorb1,   &
        iorba, iorb1a, is, kp, nok, k1, nt, nb, kin, kfin
  INTEGER :: info
  INTEGER, ALLOCATABLE :: inslab(:)
  real(kind=DP) :: z(nrz+1), r(1:ndmx,npsx), rab(1:ndmx,npsx),        &
                   betar(1:ndmx,nbrx,npsx), taunew(4,norb) 
  REAL(kind=DP), PARAMETER :: eps=1.d-8
  REAL(kind=DP) :: DDOT, dz, tr, tr1, dz1 
  COMPLEX(kind=DP), PARAMETER :: cim=(0.d0,1.d0), one=(1.d0, 0.d0), &
                                 zero=(0.d0,0.d0) 
  COMPLEX(kind=DP) :: int1d, int2d, c, d, e, f, s1, s2, s3, s4, arg,&
                      f1p, ZDOTC, fact, factm, psiper(n2d,n2d,nrzp), &
                      zk(n2d,nrzp)
  COMPLEX(kind=DP), ALLOCATABLE ::   &
     psigper(:,:), & ! psigper(g,lam)=newbg(g,lam1) psiper(lam1,lam)
     w0(:,:,:),  &   ! w0(z,g,m) are 2D Fourier components (see four.f)
     w(:,:,:),   &   ! w(z,lam,m)=psigper(g,lam)^* \exp{-igr^m_perp} 
                     !                            w0(z,g,m) 
     ci(:,:,:),  &   ! ci(m,lam,k)=\int_{z(k)}^{z(k+1)} dz 
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z-z(k))}
     di(:,:,:),  &   ! di(m,lam,k)=\int_{z(k)}^{z(k+1)} dz
                     !          w(z,lam,m)^*\exp{izk(lam,k)(z(k+1)-z)}   
     cix(:,:,:), &   !
     dix(:,:,:), &   !
     bf(:,:), an(:,:), bn(:,:),        & 
     app(:,:), bpp(:,:), al(:,:),      & 
     bl(:,:), af(:,:),                 &
     bnlf(:,:), anln(:,:), bnln(:,:),  &
     anlp(:,:), bnlp(:,:), anll(:,:),  &
     ff(:,:), fl(:,:), ezk(:,:), emzk(:,:), zk2(:,:), s1m(:,:), s2m(:,:), &
     s3m(:,:), s4m(:,:), s5m(:,:), s6m(:,:), s7m(:,:), s8m(:,:),  &
     ezk1(:,:), emzk1(:,:)
                               
  CALL start_clock('scatter_forw')

!
! Find first and last slab (kin, kfin)
!
  dz = z(2) - z(1)
  dz1 = dz/nz1
!
! Divide the slabs among CPU
!
  call divide(nrz,kin,kfin)


!------------------------
! Start of 2D Fourier components calculations and depending
! variables
!
  ALLOCATE( psigper( ngper*npol, n2d ) )
  ALLOCATE( w( nz1, n2d, norb*npol ) )
  ALLOCATE( w0( nz1, ngper, 5 ) )
  ALLOCATE( cix( nz1, n2d, norb*npol ) )
  ALLOCATE( dix( nz1, n2d, norb*npol ) )
  ALLOCATE( ci( norb*npol, n2d, nrzp ) )
  ALLOCATE( di( norb*npol, n2d, nrzp ) )     
  ALLOCATE( inslab( norb ) )
  ALLOCATE( ezk( n2d, nrzp ) )
  ALLOCATE( emzk( n2d, nrzp ) )
  ALLOCATE( ezk1( nz1, n2d ) )
  ALLOCATE( emzk1( nz1, n2d ) )
  ALLOCATE( zk2( n2d, nrzp ) )
  intw1=(0.d0,0.d0)
  intw2=(0.d0,0.d0)

  do k=kin,kfin
     kp = k-kin+1
     do lam=1,n2d
        arg=cim*tpi*zk(lam, kp)*dz
        ezk(lam,kp)=exp(arg)
        emzk(lam,kp)=exp(-arg)
        zk2(lam,kp)=cim/(2.d0*zk(lam,kp)*tpiba)
     enddo
  enddo

!
! some orbitals relations
!                    
  do iorb=1, norb
    inslab(iorb) = 0
  enddo
  do iorb = 1, norb
    if(inslab(iorb).eq.0.and.tblm(4,iorb).eq.1) then
      itt = tblm(1,iorb)
      nbb = tblm(2,iorb)
      do iorb1 = iorb, norb
       if(tblm(4,iorb1).eq.1.and.tblm(2,iorb1).eq.nbb) then
         tr = abs(taunew(3,iorb)-taunew(3,iorb1))
         if(tblm(1,iorb1).eq.itt.and.tr.le.eps) then
           inslab(iorb1) = iorb
         endif
       endif
      enddo
    endif
  enddo

  CALL start_clock('integrals')
!
! The loop over slabs to compute ci, di, and initial intw2
!
  do k = kin, kfin
    kp = k-kin+1

!    write(6,*) 'integrals k=', k

    DO lam=1,n2d
       arg=cim*zk(lam,kp)*dz1*tpi
       fact=EXP(arg)
       factm=EXP(-arg)
       ezk1(1,lam)=fact
       emzk1(1,lam)=factm
       DO k1=2,nz1
          ezk1(k1,lam)=ezk1(k1-1,lam)*fact
          emzk1(k1,lam)=emzk1(k1-1,lam)*factm
       ENDDO
    ENDDO

    if(ewind.le.100.d0) then
        CALL ZGEMM('n', 'n', ngper*npol, n2d, n2d, one, newbg,   &
                   ngper*npol, psiper(1,1,kp), n2d, zero,        &
                   psigper, ngper*npol)
    else
        psigper(:,:) = psiper(:,:,kp)
    endif

    w=(0.d0,0.d0)
    DO iorb=1, norb
      IF(cros(iorb,k).EQ.1.AND.inslab(iorb).EQ.iorb) THEN  
       mdim = 2*tblm(3,iorb)+1
       nt = tblm(1,iorb)
       nb = tblm(2,iorb)
       call four(w0, z(k), dz, tblm(1,iorb), taunew(1,iorb),  &
                 r(1,nt), rab(1,nt), betar(1,nb,nt))
       DO iorb1=1, norb
        IF(inslab(iorb1).EQ.iorb) THEN
         DO ig=1, ngper
          tr=-tpi*DDOT(2,gper(1,ig),1,taunew(1,iorb1),1)
          c=CMPLX(COS(tr),SIN(tr))
          DO lam=1, n2d
           DO is=1,npol
            d=CONJG(psigper(ngper*(is-1)+ig,lam))*c
            DO n=1, mdim
             DO kz=1, nz1
              w(kz,lam,npol*(iorb1-2+n)+is)= &
                 w(kz,lam,npol*(iorb1-2+n)+is)+d*w0(kz,ig,n)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO                
        ENDIF
       ENDDO
      ENDIF 
    ENDDO
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       IF (cros(iorba,k).EQ.1) THEN
          DO lam=1, n2d
             CALL setint(w(1,lam,iorb),cix(1,lam,iorb),dix(1,lam,iorb),  &
                                ezk1(1,lam), emzk1(1,lam), nz1)
          ENDDO
       ENDIF
    ENDDO
    DO iorb=1, norb*npol
     iorba=iorb
     IF (npol.EQ.2) iorba=(iorb+1)/2
     IF (cros(iorba,k).EQ.1) THEN   
        DO lam=1, n2d
           ci(iorb,lam,kp)=int1d(w(1,lam,iorb),   &
                        zk(lam,kp),dz,dz1,nz1,tpiba,1)
           di(iorb,lam,kp)=int1d(w(1,lam,iorb),   &
                        zk(lam,kp),dz,dz1,nz1,tpiba,-1)
        ENDDO         
        DO iorb1=1, norb*npol
           iorb1a=iorb1
           IF (npol.EQ.2) iorb1a=(iorb1+1)/2
           IF (cros(iorb1a,k).EQ.1) THEN
             DO lam=1, n2d
                intw2(iorb,iorb1)=intw2(iorb,iorb1)-                     &
                   int2d(w(1,lam,iorb),w(1,lam,iorb1),cix(1,lam,iorb1),  &
                         dix(1,lam,iorb1),ezk1(1,lam),emzk1(1,lam),      &
                         zk(lam,kp),dz1,tpiba,nz1)*zk2(lam,kp)
             ENDDO
           ENDIF 
        ENDDO
     ENDIF
    ENDDO      
  ENDDO 

  CALL stop_clock('integrals')


  DEALLOCATE(psigper)
  DEALLOCATE(cix)
  DEALLOCATE(dix)
  DEALLOCATE(w)
  DEALLOCATE(w0)
  DEALLOCATE(inslab)

!-----------------------------------

!
! Some allocation for iterative process
!
  ALLOCATE( bf( n2d, n2d ) )
  ALLOCATE( an( n2d, n2d ) )
  ALLOCATE( bn( n2d, n2d ) )
  ALLOCATE( app( n2d, n2d ) )
  ALLOCATE( bpp( n2d, n2d ) )
  ALLOCATE( al( n2d, n2d ) )
  ALLOCATE( bl( n2d, n2d ) )
  ALLOCATE( af( n2d, n2d ) )
  ALLOCATE( s1m( n2d, n2d ) )
  ALLOCATE( s2m( n2d, n2d ) )
  ALLOCATE( s3m( n2d, n2d ) )
  ALLOCATE( s4m( n2d, n2d ) )
  ALLOCATE( bnlf( n2d, norb*npol ) )
  ALLOCATE( anln( n2d, norb*npol ) )
  ALLOCATE( bnln( n2d, norb*npol ) )
  ALLOCATE( anlp( n2d, norb*npol ) )
  ALLOCATE( bnlp( n2d, norb*npol ) )
  ALLOCATE( anll( n2d, norb*npol ) )
  ALLOCATE( ff( n2d, norb*npol ) )
  ALLOCATE( fl( n2d, norb*npol ) )


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

  DO lam=1, n2d
    bf(lam,lam)=(1.d0,0.d0)
    bpp(lam,lam)=(1.d0,0.d0)
  ENDDO                                               
!
! To compute intw1, ff, fl for the first slab
!
  DO iorb=1, norb*npol
    iorba=iorb
    IF (npol.EQ.2) iorba=(iorb+1)/2
    IF (cros(iorba,kin).EQ.1) THEN
       DO lam=1, n2d
          intw1(iorb,lam)=di(iorb, lam, 1)
          arg=ezk(lam, 1)
          IF (ABS(AIMAG(zk(lam, 1))).LT.eps) THEN
           ff(lam,iorb)=-arg*CONJG(di(iorb,lam,1))*zk2(lam,1) 
           fl(lam,iorb)=-arg*CONJG(ci(iorb,lam,1))*zk2(lam,1)
          ELSE
           ff(lam,iorb)=-CONJG(ci(iorb,lam,1))*zk2(lam,1)
           fl(lam,iorb)=-CONJG(di(iorb,lam,1))*zk2(lam,1)
          ENDIF
       ENDDO
    ENDIF
  ENDDO

!------------------------------------
! The main loop over slabs 
!
  DO k=kin+1, kfin
    kp = k-kin+1
    nok = 0
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       tr=taunew(3,iorba)-taunew(4,iorba)
       IF (z(k)+dz.GT.tr) nok=iorb
    ENDDO
    DO lam=1, n2d
       DO lam1=1,n2d
          c=ZDOTC(n2d,psiper(1,lam,kp),1,psiper(1,lam1,kp-1),1)
          s1m(lam,lam1)=(zk(lam,kp)+zk(lam1,kp-1))/zk(lam,kp)*c
          s2m(lam,lam1)=(zk(lam,kp)-zk(lam1,kp-1))/zk(lam,kp)*c
          c=ezk(lam1,kp-1)
          s3m(lam,lam1)=s1m(lam,lam1)*c
          s4m(lam,lam1)=s2m(lam,lam1)*c
       ENDDO
    ENDDO

    CALL ZGEMM('n','n',n2d,n2d,n2d,one,s3m,n2d,app,n2d,one,an,n2d)
    CALL ZGEMM('n','n',n2d,n2d,n2d,one,s4m,n2d,app,n2d,one,bn,n2d)
    an= an+s2m
    bn= bn+s1m
    CALL ZGEMM('n','n',n2d,nok,n2d,one,s3m,n2d,anlp,n2d,one,anln,n2d)
    CALL ZGEMM('n','n',n2d,nok,n2d,one,s1m,n2d,fl,n2d,one,anln,n2d)
    CALL ZGEMM('n','n',n2d,nok,n2d,one,s4m,n2d,anlp,n2d,one,bnln,n2d)
    CALL ZGEMM('n','n',n2d,nok,n2d,one,s2m,n2d,fl,n2d,one,bnln,n2d)
    an=an*0.5d0
    bn=bn*0.5d0
    anln=anln*0.5d0
    bnln=bnln*0.5d0
    DO lam=1, n2d
       DO n=1, n2d
          bn(lam,n)=bn(lam,n)*emzk(lam,kp)
       ENDDO
       DO iorb=1, nok
          bnln(lam,iorb)=bnln(lam,iorb)*emzk(lam,kp)
       ENDDO
    ENDDO
    CALL DCOPY(2*n2d*n2d, an, 1, app, 1)
    CALL DCOPY(2*n2d*n2d, bn, 1, bpp, 1)
    an=(0.d0,0.d0)
    bn=(0.d0,0.d0)
    CALL DCOPY(2*norb*npol*n2d, anln, 1, anlp, 1)
    CALL DCOPY(2*norb*npol*n2d, bnln, 1, bnlp, 1)
    anln=(0.d0,0.d0)
    bnln=(0.d0,0.d0)
    fl=(0.d0,0.d0)

!
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       IF(cros(iorba,k).EQ.1) THEN
          DO lam=1, n2d
             arg=ezk(lam,kp)
             DO n=1, n2d   
              intw1(iorb,n)=intw1(iorb,n)+app(lam,n)*ci(iorb,lam,kp)+ &
                                    bpp(lam,n)*di(iorb,lam,kp)       
             ENDDO
             IF (ABS(AIMAG(zk(lam,kp))).LT.eps) THEN
              f1p=-arg*CONJG(di(iorb,lam,kp))*zk2(lam,kp)
              fl(lam,iorb)=-arg*CONJG(ci(iorb,lam,kp))*zk2(lam,kp)
             ELSE
              f1p=-CONJG(ci(iorb,lam,kp))*zk2(lam,kp)
              fl(lam,iorb)=-CONJG(di(iorb,lam,kp))*zk2(lam,kp)
             ENDIF
             bnlp(lam,iorb)=bnlp(lam,iorb)-f1p*emzk(lam,kp)
          ENDDO
       ENDIF   
    ENDDO
    DO iorb=1, norb*npol
       iorba=iorb
       IF (npol.EQ.2) iorba=(iorb+1)/2
       IF(cros(iorba,k).EQ.1) THEN 
          DO iorb1=1, norb*npol
             iorb1a=iorb1
             IF (npol.EQ.2) iorb1a=(iorb1+1)/2
             tr=taunew(3,iorb1a)-taunew(4,iorb1a)
             IF (z(k)+dz.GT.tr) THEN 
               c=(0.d0, 0.d0)
               DO lam=1, n2d
                  c=c+anlp(lam,iorb1)*ci(iorb,lam,kp)+              &
                      bnlp(lam,iorb1)*di(iorb,lam,kp)   
               ENDDO
               intw2(iorb,iorb1)=intw2(iorb,iorb1)+c
             ENDIF
          ENDDO
       ENDIF
    ENDDO
!
! Rotation of linear solutions
!

    CALL rotatef(app, bpp, bf, anlp, bnlp, bnlf, intw1, intw2,     &
                 n2d, norbf, norb, npol)
!    write(6,*) 'done k', k

  ENDDO

!---------------------------------------------

  CALL DCOPY(2*n2d*n2d, app, 1, al, 1)
!
! To compute the 2nd half of linear solutions
!
  CALL scatter_back(psiper, zk, app, bpp, an, bn, af, ci, di, &
                    ezk, emzk, s1m, s2m, s3m, s4m, kin, kfin, &
                    nrz, nrzp, norb, cros, dz)

  CALL DCOPY(2*n2d*n2d, bpp, 1, bl, 1)
  CALL DCOPY(2*n2d*norb*npol, anlp, 1, anll, 1)
  CALL DSCAL(2*norbf*npol*2*n2d, sarea, intw1, 1)
  CALL DSCAL(2*norbf*npol*norbf*npol, sarea, intw2, 1)

!
! To construct functions and derivetives on the boundaries
!
! local solutions 
  ALLOCATE( s5m( n2d, n2d ) )
  ALLOCATE( s6m( n2d, n2d ) )
  ALLOCATE( s7m( n2d, n2d ) )
  ALLOCATE( s8m( n2d, n2d ) )
!  fun0=(0.d0,0.d0)
!  fun1=(0.d0,0.d0)
!  fund0=(0.d0,0.d0)
!  fund1=(0.d0,0.d0)
  k = 1
  kp = nrzp
  do n=1, n2d
    do lam=1, n2d
      s1m(lam,n)=bf(lam,n)*ezk(lam,k)
      s2m(lam,n)=al(lam,n)*ezk(lam,kp)
      if (lam.eq.n) s2m(lam,n)=s2m(lam,n)+(1.d0,0.d0)
      s3m(lam,n)=-cim*zk(lam,k)*s1m(lam,n)*tpiba
      s4m(lam,n)=cim*zk(lam,kp)*s2m(lam,n)*tpiba
      if (lam.eq.n) s4m(lam,n)=s4m(lam,n)-2.d0*cim*zk(lam,kp)*tpiba
      s5m(lam,n)=bl(lam,n)*ezk(lam,k)
      if (lam.eq.n) s5m(lam,n)=s5m(lam,n)+(1.d0,0.d0)
      s6m(lam,n)=af(lam,n)*ezk(lam,kp)
      s7m(lam,n)=-cim*zk(lam,k)*s5m(lam,n)*tpiba
      if (lam.eq.n) s7m(lam,n)=s7m(lam,n)+2.d0*cim*zk(lam,k)*tpiba
      s8m(lam,n)=cim*zk(lam,kp)*s6m(lam,n)*tpiba
    enddo
  enddo
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,k),n2d,s1m,n2d,zero,fun0,n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,kp),n2d,s2m,n2d,zero,fun1,n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,k),n2d,s3m,n2d,zero,fund0,n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,kp),n2d,s4m,n2d,zero,fund1,n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,k),n2d,s5m,n2d,zero,  &
                                     fun0(1,n2d+1),n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,kp),n2d,s6m,n2d,zero, &
                                     fun1(1,n2d+1),n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,k),n2d,s7m,n2d,zero,  &
                                     fund0(1,n2d+1),n2d)
  call ZGEMM('n','n',n2d,n2d,n2d,one,psiper(1,1,kp),n2d,s8m,n2d,zero, &
                                     fund1(1,n2d+1),n2d)

! nonlocal solutions
  funl0=(0.d0,0.d0)
  funl1=(0.d0,0.d0)
  fundl0=(0.d0,0.d0)
  fundl1=(0.d0,0.d0)
  do iorb = 1, norb*npol
    do lam = 1, n2d
      s1=ff(lam,iorb)+bnlf(lam,iorb)*ezk(lam,k)
      s2=fl(lam,iorb)+anll(lam,iorb)*ezk(lam,kp)
      s3=-cim*zk(lam,k)*s1*tpiba
      s4=cim*zk(lam,kp)*s2*tpiba
      do ig=1, n2d
       funl0(ig,iorb)=funl0(ig,iorb)+psiper(ig,lam,k)*s1
       funl1(ig,iorb)=funl1(ig,iorb)+psiper(ig,lam,kp)*s2
       fundl0(ig,iorb)=fundl0(ig,iorb)+psiper(ig,lam,k)*s3
       fundl1(ig,iorb)=fundl1(ig,iorb)+psiper(ig,lam,kp)*s4
      enddo
    enddo
  enddo

!
! To construct the functions in the whole rigion zin<z<zfin in the
! case of multiparallel running
!
#ifdef __PARA
  CALL rotproc(fun0, fund0, fun1, fund1, funl0, fundl0, funl1, &
               fundl1, intw1, intw2, n2d, norbf, norb)                    
#endif       
  CALL stop_clock('scatter_forw')

  DEALLOCATE(ci)
  DEALLOCATE(di)
  DEALLOCATE(bf)
  DEALLOCATE(an)
  DEALLOCATE(bn)
  DEALLOCATE(app)
  DEALLOCATE(bpp)
  DEALLOCATE(al)
  DEALLOCATE(bl)
  DEALLOCATE(af)
  DEALLOCATE(bnlf)
  DEALLOCATE(anln)
  DEALLOCATE(bnln)
  DEALLOCATE(anlp)
  DEALLOCATE(bnlp)
  DEALLOCATE(anll)
  DEALLOCATE(ff)
  DEALLOCATE(fl)
  DEALLOCATE(s1m)
  DEALLOCATE(s2m)
  DEALLOCATE(s3m)
  DEALLOCATE(s4m)
  DEALLOCATE(s5m)
  DEALLOCATE(s6m)
  DEALLOCATE(s7m)
  DEALLOCATE(s8m)
  DEALLOCATE(ezk)
  DEALLOCATE(emzk)
  DEALLOCATE(ezk1)
  DEALLOCATE(emzk1)
  DEALLOCATE(zk2)

  RETURN
END SUBROUTINE scatter_forw
