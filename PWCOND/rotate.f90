!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine rotatef(app, bpp, bf, anlp, bnlp, bnlf, intw1, intw2,    &
                   n2d, norbf, norbnow)
!
! This subroutine makes a linear combination of the solutions 
! in such a way that bpp at this slab becomes a delta symbol.
! It works for forward iterative process. 
!
#include "f_defs.h"
  USE kinds, only : DP
  implicit none
  integer :: norbf, n2d, norbnow, lam, n, n1, iorb, iorb1, info
  integer, allocatable :: ipiv(:) 
  complex(kind=DP) ::                                   &
          app(n2d, n2d),       &      ! coeff. of local functions 
          bpp(n2d, n2d),       &      !          -- 
          bf(n2d, n2d),        &      !          -- 
          anlp(n2d, norbnow),  &      ! coeff. of nonloc. functions 
          bnlp(n2d, norbnow),  &      !          -- 
          bnlf(n2d, norbnow),  &      !          --  
          intw1(norbf, 2*n2d), &      ! integral of loc. fun. 
          intw2(norbf, norbf)         ! integral of nonloc. fun.   
  complex(kind=DP), allocatable :: x(:), y(:), h(:,:)

  allocate( x( n2d ) )
  allocate( y( n2d ) )
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
  bpp=(0.d0,0.d0)
  do lam=1, n2d
     x=(0.d0,0.d0)
     y=(0.d0,0.d0)
     do n=1, n2d
        do n1=1, n2d
           x(n)=x(n)+h(n1,n)*app(lam,n1)
           y(n)=y(n)+h(n1,n)*bf(lam,n1)
        enddo
     enddo
     do n=1, n2d
        app(lam,n)=x(n)
        bf(lam,n)=y(n)
     enddo
     bpp(lam,lam)=(1.d0,0.d0)
  enddo   

!
! To recalculate intw1 with new functions
!
  do iorb=1, norbnow
     x=(0.d0,0.d0)
     do n=1, n2d
        do n1=1, n2d
           x(n)=x(n)+h(n1,n)*intw1(iorb,n1)
        enddo
     enddo
     do n=1, n2d
        intw1(iorb,n)=x(n)
     enddo
  enddo  

!
! To reobtain nonlinear functions and the integrals 
! intw2 on them.
!
  do iorb=1, norbnow
    do lam=1, n2d
      do n=1, n2d
       anlp(lam,iorb)=anlp(lam,iorb)-bnlp(n,iorb)*app(lam,n)
       bnlf(lam,iorb)=bnlf(lam,iorb)-bnlp(n,iorb)*bf(lam,n) 
      enddo
    enddo
    do iorb1=1, norbnow
       do n=1, n2d
         intw2(iorb,iorb1)=intw2(iorb,iorb1)-           &
                     bnlp(n,iorb1)*intw1(iorb,n)
       enddo
    enddo
  enddo 
  bnlp=(0.d0,0.d0)

  deallocate(x)
  deallocate(y)
  deallocate(h)
  deallocate(ipiv)
     
  return
end subroutine rotatef
!------------------------------------------

subroutine rotateb (app, bpp, af, intw1, n2d, norbf, norbnow)
!
! This subroutine makes a linear combination of the solutions
! in such a way that app at this slab becomes a delta symbol.
! It works for backward iterative process.
!
#include "f_defs.h"
  USE kinds, only : DP
  implicit none 

  integer :: norbf, n2d, norbnow, lam, n, n1, iorb, info
  integer, allocatable :: ipiv(:) 
  complex(kind=DP) :: app(n2d,n2d), af(n2d,n2d), bpp(n2d,n2d), &
                      intw1(norbf,2*n2d) 
  complex(kind=DP), allocatable :: x(:), y(:), h(:,:)

  allocate( x( n2d ) )
  allocate( y( n2d ) )
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
  app=(0.d0,0.d0)
  do lam=1, n2d
     x=(0.d0,0.d0)
     y=(0.d0,0.d0)
     do n=1, n2d
        do n1=1, n2d
           x(n)=x(n)+h(n1,n)*bpp(lam,n1)
           y(n)=y(n)+h(n1,n)*af(lam,n1)
        enddo
     enddo
     do n=1, n2d
        bpp(lam,n)=x(n)
        af(lam,n)=y(n)
     enddo
     app(lam,lam)=(1.d0,0.d0)
  enddo   

!
! To recalculate intw1 with new functions
!                                         
  do iorb=1, norbnow
     x=(0.d0,0.d0)
     do n=1, n2d
        do n1=1, n2d
           x(n)=x(n)+h(n1,n)*intw1(iorb,n2d+n1)
        enddo
     enddo
     do n=1, n2d
        intw1(iorb,n2d+n)=x(n)
     enddo
  enddo

  deallocate(x)
  deallocate(y)
  deallocate(h)
  deallocate(ipiv)

  return
end subroutine rotateb
