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
#include "f_defs.h"
  use pwcom
  use cond
#ifdef __PARA
  use para
#endif
  implicit none

  integer ::                                                & 
             i, j, ij, ijx, k, n, p, il, ik, kstart, klast, &
             ix, jx, kx, ir, ir1, is, ixy, info
  integer :: ionode_id
  integer, allocatable :: ipiv(:) 

  real(kind=DP), parameter :: eps = 1.d-8
  real(kind=DP) :: zlen, dz, arg, bet   
  real(kind=DP), allocatable :: gz(:), auxr(:)

  complex(kind=DP), parameter :: cim = (0.d0,1.d0)
  complex(kind=DP) :: caux
  complex(kind=DP), allocatable :: aux(:), amat(:,:)

  logical :: lg

  allocate( ipiv( nrz ) )
  allocate( gz( nrz ) )
  allocate( aux( nrx1*nrx2*nrx3 ) )
  allocate( auxr( nrxx ) )
  allocate( amat( nrz, nrz ) )

!
!  Compute the Gz vectors in the z direction
!
  do k = 1, nrz
     il = k-1
     if (il.gt.nrz/2) il = il-nrz
     gz(k) = il*bg(3,3)
  enddo

!
!     To form local potential on the real space mesh
!
  auxr(:) = vltot(:) + vr(:,iofspin) 

!
! To collect the potential from different CPUs
!
  aux(:) = (0.d0,0.d0)
#ifdef __PARA
  kstart = 1
  do i=1, me-1
    kstart = kstart+npp(i)
  enddo
  klast = kstart+npp(me)-1 
#endif
  do i = 1, nrx1*nrx2
    do k = 1, nr3
      lg = .true.
      ir = i+(k-1)*nrx2*nrx1
      ir1 = ir       
#ifdef __PARA
      if(k.ge.kstart.and.k.le.klast) then
        lg = .true.
        ir1 = i+(k-kstart)*nrx2*nrx1 
      else
        lg = .false.
      endif
#endif                       
      if (lg) then
        aux(ir) = auxr(ir1)
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

  vppot(:,:) = 0.d0
  do i = 1, nrx
    if(i.gt.nrx/2+1) then
        ix = nr1-(nrx-i) 
    else
        ix = i
    endif
    do j = 1, nry
      if(j.gt.nry/2+1) then
         jx = nr2-(nry-j)
      else
         jx = j
      endif 
      ij = i+(j-1)*nrx
      ijx = ix+(jx-1)*nrx1 

      do k = 1, nrz
        il = k-1
        if (il.gt.nrz/2) il = il-nrz
        if(il.le.nr3/2.and.il.ge.-(nr3-1)/2) then

         if(k.gt.nrz/2+1) then 
            kx = nr3-(nrz-k)  
         else
            kx = k
         endif 
         vppot(k, ij) = aux(ijx+(kx-1)*nrx1*nrx2)

        endif
      enddo
    enddo
  enddo

!
! set up the matrix for the linear system
!


  zlen = z(nrz+1)-z(1)
  dz = z(2)-z(1)

  do n = 1, nrz
    bet = gz(n)*dz*tpi
    if (abs(gz(n)).gt.eps) then
      caux = cim*(CMPLX(cos(bet),-sin(bet))-(1.d0,0.d0))  &
                                  /(zlen*gz(n)*tpi)
    else
      caux = dz/zlen
    endif
    do p = 1, nrz
      arg = gz(n)*z(p)*tpi
      amat(n,p) = CMPLX(cos(arg),-sin(arg))*caux
    enddo
  enddo

!
! solve the linear system
!
  call ZGESV(nrz, nrx*nry, amat, nrz, ipiv, vppot, nrz, info)
  call errore ('poten','info different from zero',abs(info))

 ! do p = 1, nrz
 !   write(6,'(i5,2f12.6)') p, real(vppot(p,2)), imag(vppot(p,2))
 ! enddo

  deallocate(ipiv) 
  deallocate(gz) 
  deallocate(aux) 
  deallocate(auxr) 
  deallocate(amat) 

  return
end subroutine poten
