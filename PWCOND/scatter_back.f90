!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine scatter_back(app, bpp, an, bn, af, ci, di,       &
                        startk, lastk, norbnow, orbin, dz)
!
! This subroutine computes the second half of the local functions
! and their integrals.
!
#include "machine.h"
  use pwcom
  use cond
  implicit none  
  
  integer :: kz, n, lam, lam1, k, iorb, startk, lastk, norbnow, &
             orbin
  real(kind=DP) :: dz
  complex(kind=DP), parameter :: cim=(0.d0, 1.d0) 
  complex(kind=DP) :: c, s1, s2, s3, s4, ZDOTC, CONJG,          & 
                      app(n2d,n2d), an(n2d,n2d), bn(n2d,n2d),   &
                      bpp(n2d,n2d), af(n2d,n2d),                &
                      ci(norbnow,n2d,nrzp), di(norbnow,n2d,nrzp)

  app=(0.d0,0.d0)
  an=(0.d0,0.d0)
  bpp=(0.d0,0.d0)
  bn=(0.d0,0.d0)
  af=(0.d0,0.d0)
!
! The last slab calculations
!
  do lam=1, n2d
     af(lam,lam)=(1.d0,0.d0)
     app(lam,lam)=(1.d0,0.d0)
  enddo
  kz=nkofz(lastk)
  do iorb=1, norbnow
    if (cross(orbin+iorb,lastk).eq.1) then
      do n=1, n2d
        intw1(iorb,n2d+n)=ci(iorb,n,kz)
      enddo
    endif
  enddo

!
! The loop over slabs
!
  do k=lastk-1, startk, -1
     kz=nkofz(k)
     do lam=1, n2d
       do lam1=1,n2d
          c=ZDOTC(n2d,psiper(1,lam,kz),1,psiper(1,lam1,kz+1),1)  
          s1=(zk(lam,kz)+zk(lam1,kz+1))/zk(lam,kz)*c
          s2=(zk(lam,kz)-zk(lam1,kz+1))/zk(lam,kz)*c
          c=EXP(tpi*cim*zk(lam1,kz+1)*dz) 
          s3=s1*c
          s4=s2*c
          do n=1, n2d
            an(lam,n)= an(lam,n)+s4*bpp(lam1,n)
            bn(lam,n)= bn(lam,n)+s3*bpp(lam1,n)
          enddo
          an(lam,lam1)= an(lam,lam1)+s1
          bn(lam,lam1)= bn(lam,lam1)+s2                         
       enddo
       do n=1, n2d
          an(lam,n)=an(lam,n)*0.5d0*EXP(-tpi*cim*zk(lam,kz)*dz)
          bn(lam,n)=bn(lam,n)*0.5d0
       enddo
     enddo
     do iorb=1, norbnow
       if (cross(orbin+iorb,k).eq.1) then
         do lam=1, n2d
           do n=1, n2d   
             intw1(iorb,n2d+n)=intw1(iorb,n2d+n)+   &
                        an(lam,n)*ci(iorb,lam,kz)+    &
                        bn(lam,n)*di(iorb,lam,kz)       
           enddo
         enddo
       endif
     enddo
     call DCOPY(2*n2d*n2d, an, 1, app, 1)
     call DCOPY(2*n2d*n2d, bn, 1, bpp, 1)
     an=(0.d0,0.d0)
     bn=(0.d0,0.d0)
     call rotateb(app, bpp, af, intw1, n2d, norbf, norbnow)
  enddo 

  return
end subroutine scatter_back
