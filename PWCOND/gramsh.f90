!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine gramsh (n, nvec, nstart, nfinish,         & 
                    psibase, psiprob, ndim, epsproj)
!
! This routine orthogonalizes a set of vectors with respect 
! to the basis set and supplies the latter. It uses the 
! Gram-Schmidt method. 
!
#include "f_defs.h"
  USE kinds, only : DP
  implicit none
  integer ::         & 
     n,              &  ! input: physical dimension
     nvec,           &  ! input: number of vectors  
     nstart,         &  ! input: first vector to orthogonalize
     nfinish,        &  ! input: last vector to orthogonalize
     ndim,           &  ! inp/out: dimension of psibase old/new 
     ivec,           &  ! counter on vectors
     ic,             &  ! coordinates
     ivecp              ! counter on vectors    
  real(kind=DP) ::   &
     epsproj,        &  ! accuracy
     norm,           &  ! the norm of a vector
     DDOT               ! to compute the dot product of two vectors  
  real(kind=DP), parameter :: eps=1.d-8
  complex(kind=DP) :: &
     psibase(n,n),    & ! i/o:basis vector set
     psiprob(n,nvec), & ! i/o:vectors to be orthog. and added to psibas
     ZDOTC              ! to compute scalar products
  complex(kind=DP), allocatable ::  &
     ps(:)                        ! the scalar products  

  allocate( ps( n ) )

  if (ndim.eq.n) return

  do ivec = nstart, nfinish
   if(ndim.lt.n) then
!
! To find orthogonal to psibase projection of psiprob
!
     do ivecp=1, ndim
       ps(ivecp)=ZDOTC(n,psibase(1,ivecp),1,psiprob(1,ivec),1)
     enddo
     do ivecp=1,ndim
       call ZAXPY (n,-ps(ivecp),psibase(1,ivecp),1,psiprob(1,ivec),1)
     enddo
! and its norm
     norm=DDOT(2*n,psiprob(1,ivec),1,psiprob(1,ivec),1)
!
! adding (or not) psiprob to psibase
!
     if (norm.le.-eps) then
       print*,'norma = ',norm,ivec
       call errore ('gramsh',' negative norm in S ',1) 
     endif        
     if (abs(norm).gt.epsproj) then
       ndim=ndim+1 
       if (ndim.eq.n) then
          psibase=(0.d0,0.d0)
          do ivecp=1, n
             psibase(ivecp,ivecp)=(1.d0,0.d0)
          enddo                      
       else
         norm = 1.d0/sqrt(norm)
         call DSCAL( 2*n,norm,psiprob(1,ivec),1 )
         call DCOPY(2*n,psiprob(1,ivec),1,psibase(1,ndim),1) 
       endif
     endif
   endif
  enddo

  deallocate(ps)

  return
end subroutine gramsh
