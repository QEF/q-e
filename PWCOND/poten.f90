!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine poten
!
! This subroutine computes the 2D Fourier components of the
! local potential in each slab.
!
#include "machine.h"
  use pwcom
  use para
  use cond

  implicit none
  integer :: i, j, k, n, p, il, ik, kstart, klast, &
             ix, jx, kx, ir, ir1, is, ixy, info
  integer, allocatable :: ipiv(:) 
  real(kind=DP) :: arg, bet, alph   
  real(kind=DP), parameter :: eps=1.d-8
  complex(kind=DP), parameter :: cim=(0.d0,1.d0)
  real(kind=DP), allocatable :: gz(:), auxr(:,:)
  logical :: lg
  complex(kind=DP) :: xfact, auxy1, auxy2 
  complex(kind=DP), allocatable :: aux(:), amat(:,:) 

  allocate( ipiv( nrz ) )
  allocate( gz( nrz ) )
  allocate( aux( nrx1*nrx2*nrx3 ) )
  allocate( auxr( nrxx, nspin ) )
  allocate( amat( nrz, nrz ) )

!
!  Compute the Gz vectors in the z direction
!
  do k=1,nrz
     il=k-1
     if (il.gt.nrz/2) il=il-nrz
       gz(k)=il*bg(3,3)
  enddo

!
!     To form local potential on the real space mesh
!
  do is = 1, nspin
    do ir = 1, nrxx
      auxr (ir, is) = vltot (ir) + vr (ir, is)
    enddo
  enddo

!
! To collect the potential from different CPUs
!
  kstart=1
  klast=nr3 
#ifdef __PARA
  aux=(0.d0,0.d0)
  do i=1, me-1
    kstart=kstart+npp(i)
  enddo
  klast=kstart+npp(me)-1 
#endif
  do i=1,nr1*nr2
    do k=1,nr3
      lg=.true.
      ir=i+(k-1)*nr2*nr1
      ir1=ir       
#ifdef __PARA
      if(k.ge.kstart.and.k.le.klast) then
        lg=.true.
        ir1=i+(k-kstart)*nr2*nr1 
      else
        lg=.false.
      endif
#endif                       
      if (lg) then
        aux(ir)=auxr(ir1, iofspin)
      endif 
    enddo
  enddo
#ifdef __PARA
  call reduce (2*nrx1*nrx2*nrx3,aux)
#endif    

!
! To find FFT of the local potential
!
#ifdef __PARA
!
!  This FFT is needed to make a non-parallel FFT in the parallel case
!
  call cft3sp(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
#else
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
#endif

  do i=1,nrx
    ix=i
    if(ix.gt.nrx/2+1) ix=nr1-(nrx-i)
    do j=1,nry
      jx=j
      if(jx.gt.nry/2+1) jx=nr2-(nry-j) 
      ixy=i+(j-1)*nrx
      do k=1,nrz
        if(abs(gz(k)).gt.gz(nr3s/2)) then
          vppot(k, ixy)=0.d0
        else
          kx=k
          if(kx.gt.nrz/2+1) kx=nr3-(nrz-k)
          ir=ix+(jx-1)*nr1+(kx-1)*nr2*nr1
          vppot(k, ixy)=aux(ir)
        endif
      enddo
    enddo
  enddo

!
! set up the matrix for the linear system
!
  do n=1,nrz
    do p=1,nrz
      arg=gz(n)*z(p)*tpi
      bet=gz(n)*(z(p+1)-z(p))*tpi             
      if (abs(gz(n)).gt.eps) then
        xfact=cim*(CMPLX(cos(bet),-sin(bet))-(1.d0,0.d0))  &
                                    /zl/gz(n)/tpi
      else
        xfact=(z(p+1)-z(p))/zl
      endif
      amat(n,p)=CMPLX(cos(arg),-sin(arg))*xfact
    enddo
  enddo

!
! solve the linear system
!
  call ZGESV(nrz, nrx*nry, amat, nrz, ipiv, vppot, nrz, info)
  call errore ('poten','info different from zero',abs(info))


  deallocate(ipiv) 
  deallocate(gz) 
  deallocate(aux) 
  deallocate(auxr) 
  deallocate(amat) 

  return
end subroutine poten
