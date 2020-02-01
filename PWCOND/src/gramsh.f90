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
  USE kinds, only : DP
  implicit none
  integer ::         &
     n,              &  ! input: physical dimension
     nvec,           &  ! input: number of vectors
     nstart,         &  ! input: first vector to orthogonalize
     nfinish,        &  ! input: last vector to orthogonalize
     ndim,           &  ! inp/out: dimension of psibase old/new
     ivec,           &  ! counter on vectors
     ivecp              ! counter on vectors
  real(DP) ::   &
     epsproj,        &  ! accuracy
     norm,           &  ! the norm of a vector
     ddot               ! to compute the dot product of two vectors
  real(DP), parameter :: eps=1.d-8
  complex(DP) :: &
     psibase(n,n),    & ! i/o:basis vector set
     psiprob(n,nvec), & ! i/o:vectors to be orthog. and added to psibas
     zdotc              ! to compute scalar products
  complex(DP), allocatable ::  &
     ps(:)                        ! the scalar products

  allocate( ps( n ) )

  if (ndim.eq.n) return

  do ivec = nstart, nfinish
   if(ndim.lt.n) then
!
! To find orthogonal to psibase projection of psiprob
!
     do ivecp=1, ndim
       ps(ivecp)=zdotc(n,psibase(1,ivecp),1,psiprob(1,ivec),1)
     enddo
     do ivecp=1,ndim
       call zaxpy (n,-ps(ivecp),psibase(1,ivecp),1,psiprob(1,ivec),1)
     enddo
! and its norm
     norm=ddot(2*n,psiprob(1,ivec),1,psiprob(1,ivec),1)
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
         call dscal( 2*n,norm,psiprob(1,ivec),1 )
         call dcopy(2*n,psiprob(1,ivec),1,psibase(1,ndim),1)
       endif
     endif
   endif
  enddo

  deallocate(ps)

  return
end subroutine gramsh
