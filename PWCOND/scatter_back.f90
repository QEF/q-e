!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generalized to spinor wavefunctions and spin-orbit Oct. 2004 (ADC).
!
!
subroutine scatter_back(psiper, zk, app, bpp, an, bn, af, ci, di, &
                        ezk, emzk, s1m, s2m, s3m, s4m, kin, kfin, &
                        nrz, nrzp, norb, cros, dz)
!
! This subroutine computes the second half of the local functions
! and their integrals.
!
#include "f_defs.h"
  use noncollin_module, ONLY : npol
  use cond
  implicit none  
  
  integer :: kp, n, lam, lam1, k, iorb, kin, kfin, norb,        &
             nrz, nrzp, iorba, cros(norb,nrz)
  real(kind=DP) :: dz
  complex(kind=DP), parameter :: cim=(0.d0, 1.d0), one=(1.d0,0.d0)
  complex(kind=DP) :: c, s1, s2, s3, s4, ZDOTC, CONJG,          & 
                      psiper(n2d,n2d,nrzp), zk(n2d,nrzp),       &
                      app(n2d,n2d), an(n2d,n2d), bn(n2d,n2d),   &
                      bpp(n2d,n2d), af(n2d,n2d),                &
                      ci(norb*npol,n2d,nrzp), di(norb*npol,n2d,nrzp), &
                      ezk(n2d,nrzp), emzk(n2d,nrzp), &
                      s1m(n2d,n2d),s2m(n2d,n2d),s3m(n2d,n2d),s4m(n2d,n2d)

  call start_clock('scatter_back')
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
  do iorb=1, norb*npol
    iorba=iorb
    if (npol.eq.2) iorba=(iorb+1)/2
    if (cros(iorba,kfin).eq.1) then
      do n=1, n2d
        intw1(iorb,n2d+n)=ci(iorb,n,nrzp)
      enddo
    endif
  enddo

!
! The loop over slabs
!
  do k=kfin-1, kin, -1
     kp = k-kin+1
     do lam=1, n2d
       do lam1=1,n2d
          c=ZDOTC(n2d,psiper(1,lam,kp),1,psiper(1,lam1,kp+1),1)  
          s1m(lam,lam1)=(zk(lam,kp)+zk(lam1,kp+1))/zk(lam,kp)*c
          s2m(lam,lam1)=(zk(lam,kp)-zk(lam1,kp+1))/zk(lam,kp)*c
          c=ezk(lam1,kp+1)
          s3m(lam,lam1)=s1m(lam,lam1)*c
          s4m(lam,lam1)=s2m(lam,lam1)*c
       enddo
     enddo
     call ZGEMM('n','n',n2d,n2d,n2d,one,s4m,n2d,bpp,n2d,one,an,n2d)
     call ZGEMM('n','n',n2d,n2d,n2d,one,s3m,n2d,bpp,n2d,one,bn,n2d)
     an=an+s1m
     bn=bn+s2m                         
     an=an*0.5d0
     bn=bn*0.5d0
     do lam=1,n2d
       do n=1, n2d
          an(lam,n)=an(lam,n)*emzk(lam,kp)
       enddo
     enddo
     do iorb=1, norb*npol
       iorba=iorb
       if (npol.eq.2) iorba=(iorb+1)/2
       if (cros(iorba,k).eq.1) then
         do lam=1, n2d
           do n=1, n2d   
             intw1(iorb,n2d+n)=intw1(iorb,n2d+n)+   &
                        an(lam,n)*ci(iorb,lam,kp)+    &
                        bn(lam,n)*di(iorb,lam,kp)       
           enddo
         enddo
       endif
     enddo
     call DCOPY(2*n2d*n2d, an, 1, app, 1)
     call DCOPY(2*n2d*n2d, bn, 1, bpp, 1)
     an=(0.d0,0.d0)
     bn=(0.d0,0.d0)
     call rotateb(app, bpp, af, intw1, n2d, norbf, norb, npol)

!     write(6,*) k, kin

  enddo 

  call stop_clock('scatter_back')

  return
end subroutine scatter_back
