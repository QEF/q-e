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
subroutine rotatef(app, bpp, bf, anlp, bnlp, bnlf, intw1, intw2,    &
                   n2d, norbf, norbnow, npol)
!
! This subroutine makes a linear combination of the solutions 
! in such a way that bpp at this slab becomes a delta symbol.
! It works for forward iterative process. 
!
#include "f_defs.h"
  USE kinds, only : DP
  implicit none
  integer :: norbf, n2d, norbnow, lam, n, n1, iorb, iorb1, npol, info
  integer, allocatable :: ipiv(:) 
  complex(DP) ::                                   &
          app(n2d, n2d),       &      ! coeff. of local functions 
          bpp(n2d, n2d),       &      !          -- 
          bf(n2d, n2d),        &      !          -- 
          anlp(n2d, norbnow*npol),  &      ! coeff. of nonloc. functions 
          bnlp(n2d, norbnow*npol),  &      !          -- 
          bnlf(n2d, norbnow*npol),  &      !          --  
          intw1(norbf*npol, 2*n2d), &      ! integral of loc. fun. 
          intw2(norbf*npol, norbf*npol)    ! integral of nonloc. fun.   
  complex(DP), allocatable :: h(:,:), aux(:,:)
  complex(DP), parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)

  call start_clock('rotatef')
  allocate( h( n2d, n2d ) )
  allocate( ipiv( n2d ) )

!
! To find the needed matrix h of the linear transformation
!
  h=(0.d0,0.d0)
  do lam=1, n2d
     h(lam,lam)=(1.d0,0.d0)
  enddo
  call ZGESV(n2d,n2d,bpp,n2d,ipiv,h,n2d,info)

!
! To rotate app, bf, bpp
!
  allocate( aux( n2d, n2d ) )
  call ZGEMM('n','n',n2d,n2d,n2d,one,app,n2d,h,n2d,zero,aux,n2d)
  app=aux
  call ZGEMM('n','n',n2d,n2d,n2d,one,bf,n2d,h,n2d,zero,aux,n2d)
  bf=aux
  bpp=(0.d0,0.d0)
  do lam=1, n2d
     bpp(lam,lam)=(1.d0,0.d0)
  enddo   
  deallocate(aux)

!
! To recalculate intw1 with new functions
!
  if (norbnow==0) goto 100

  allocate( aux( norbf*npol, n2d ) )
  call ZGEMM('n','n',norbnow*npol,n2d,n2d,one,intw1,norbf*npol,h,n2d,zero,&
                 aux,norbf*npol)
  do iorb=1,norbnow*npol
     do n=1,n2d
        intw1(iorb,n)=aux(iorb,n)
     enddo
  enddo
  deallocate(aux)

!
! To reobtain nonlinear functions and the integrals 
! intw2 on them.
!
  call ZGEMM('n','n',n2d,norbnow*npol,n2d,-one,app,n2d,bnlp,n2d,one,&
                 anlp,n2d)
  call ZGEMM('n','n',n2d,norbnow*npol,n2d,-one,bf,n2d,bnlp,n2d,one,&
                 bnlf,n2d)
  call ZGEMM('n','n',norbnow*npol,norbnow*npol,n2d,-one,intw1,norbf*npol, &
                     bnlp,n2d,one,intw2,norbf*npol)
  bnlp=(0.d0,0.d0)

100 continue

  deallocate(h)
  deallocate(ipiv)
     
  call stop_clock('rotatef')
  return
end subroutine rotatef
!------------------------------------------

subroutine rotateb (app, bpp, af, intw1, n2d, norbf, norbnow, npol)
!
! This subroutine makes a linear combination of the solutions
! in such a way that app at this slab becomes a delta symbol.
! It works for backward iterative process.
!
#include "f_defs.h"
  USE kinds, only : DP
  implicit none 

  integer :: norbf, n2d, norbnow, lam, n, n1, iorb, npol, info
  integer, allocatable :: ipiv(:) 
  complex(DP) :: app(n2d,n2d), af(n2d,n2d), bpp(n2d,n2d), &
                      intw1(norbf*npol,2*n2d) 
  complex(DP), allocatable :: h(:,:), aux(:,:), aux1(:,:)
  complex(DP), parameter :: one=(1.d0,0.d0), zero=(0.d0,0.d0)

  call start_clock('rotateb')
  allocate( h( n2d, n2d ) )
  allocate( ipiv( n2d ) )  

!
! To find the needed matrix h of the linear transformation
!                                                           
  h=(0.d0,0.d0)
  do lam=1, n2d
     h(lam,lam)=(1.d0,0.d0)
  enddo
  call ZGESV(n2d, n2d, app, n2d, ipiv, h, n2d, info)
!
! To rotate app, bpp, af
!                           
  allocate( aux( n2d, n2d ) )
  call ZGEMM('n','n',n2d,n2d,n2d,one,bpp,n2d,h,n2d,zero,aux,n2d)
  bpp=aux
  call ZGEMM('n','n',n2d,n2d,n2d,one,af,n2d,h,n2d,zero,aux,n2d)
  af=aux

  app=(0.d0,0.d0)
  do lam=1, n2d
     app(lam,lam)=(1.d0,0.d0)
  enddo   
  deallocate(aux)
!
! To recalculate intw1 with new functions
!                                         
  if (norbnow==0) goto 100

  allocate( aux( norbf*npol, n2d ) )
  allocate( aux1( norbf*npol, n2d ) )
  do iorb=1,norbnow*npol
     do n=1,n2d
        aux1(iorb,n)= intw1(iorb,n2d+n)
     enddo
  enddo
  call ZGEMM('n','n',norbnow*npol,n2d,n2d,one,aux1,norbf*npol, &
                     h,n2d,zero,aux,norbf*npol)
  do iorb=1,norbnow*npol
     do n=1,n2d
        intw1(iorb,n2d+n)= aux(iorb,n)
     enddo
  enddo
  deallocate(aux)
  deallocate(aux1)

100 continue

  deallocate(h)
  deallocate(ipiv)

  call stop_clock('rotateb')

  return
end subroutine rotateb
