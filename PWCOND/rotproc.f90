!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine rotproc (fun0, fund0, fun1, fund1, funl0, fundl0, funl1,  & 
                   fundl1, intw1, intw2, n2d, norbf, norbnow)
!
! This subroutine implements a matching procedure to construct 
! local and nonlocal functions on the whole region from those computed
! by different CPU.
! It works well with 1, 2, 4, 8, 16... CPU.
! The matching scheme with 8 CPU looks like:

! |   1   |   2   |   3   |   4   |   5   |   6   |   7   |   8   |
!         +               +               +               +
! |               |               |               |               | 
!                 +                               +                
! |                               |                               | 
!                                 +                                 
! |                                                               |
!
! So in this case there are 3 matching steps. 
!
#include "f_defs.h"
  USE kinds, only : DP
#ifdef __PARA
  USE mp_global, ONLY: nproc
  use para
  USE parallel_include
  
  implicit none 


  integer :: k, ig, n, lam, lam1, iorb, iorb1, norbf, norbnow, n2d,  &
             ibound, numb, ninsl, ib, icolor, ikey, new_comm, info
  integer, allocatable :: ipiv(:) 
  complex(kind=DP) :: fun0(n2d, 2*n2d),    & ! phi_n(0) 
                      fund0(n2d, 2*n2d),   & ! phi'_n(0)  
                      fun1(n2d, 2*n2d),    & ! phi_n(d) 
                      fund1(n2d, 2*n2d),   & ! phi'_n(d) 
                      funl0(n2d, norbf),   & ! phi_alpha(0) 
                      fundl0(n2d, norbf),  & ! phi'_alpha(0) 
                      funl1(n2d, norbf),   & ! phi_alpha(d) 
                      fundl1(n2d, norbf),  & ! phi'_alpha(d) 
                      intw1(norbf, 2*n2d), & ! integrals on phi_n 
                      intw2(norbf, norbf)    ! integrals on phi_alpha

  complex(kind=DP), allocatable :: x(:), y(:), amat(:,:), vec(:,:),  &
                               amat_aux(:,:), vec_aux(:,:)

  if(nproc.eq.1) return

  allocate( x( n2d ) )
  allocate( y( n2d ) )    
  allocate( amat( 2*n2d, 2*n2d ) )
  allocate( amat_aux( 2*n2d, 2*n2d ) ) 
  allocate( vec( 2*n2d, 2*n2d+norbnow ) )
  allocate( vec_aux( 2*n2d, 2*n2d+norbnow ) )
  allocate( ipiv( 2*n2d ) )

  numb=0
  ibound=nproc/2
  ninsl=1
  ib=2*me-me/2*2    

  do while(ibound.gt.0) 

!
!   To find the matching coefficients for a group of CPU 
!
    icolor=(ib+ninsl-1)/(2*ninsl)
    ikey=(me+ninsl)-ib               
    call mpi_barrier (MPI_COMM_WORLD, info)      
    call mpi_comm_split(MPI_COMM_WORLD, icolor, ikey, new_comm, info)
    amat=(0.d0,0.d0)
    vec=(0.d0,0.d0)

    if(me.eq.ib) then
      do lam=1, n2d
        do lam1=1, n2d
          amat(lam, n2d+lam1)=-fun0(lam,n2d+lam1)
          amat(n2d+lam, n2d+lam1)=-fund0(lam,n2d+lam1)
          vec(lam, lam1)=fun0(lam, lam1)
          vec(n2d+lam, lam1)=fund0(lam, lam1)
        enddo
        do iorb=1, norbnow
          vec(lam, 2*n2d+iorb)=funl0(lam, iorb)
          vec(n2d+lam, 2*n2d+iorb)=fundl0(lam, iorb)
        enddo
      enddo              
      numb=numb+1
    endif
    if(me.eq.ib-1) then
      do lam=1, n2d
        do lam1=1, n2d
          amat(lam, lam1)=fun1(lam,lam1)
          amat(n2d+lam, lam1)=fund1(lam,lam1)
          vec(lam, n2d+lam1)=-fun1(lam, n2d+lam1)
          vec(n2d+lam, n2d+lam1)=-fund1(lam, n2d+lam1)
        enddo
        do iorb=1, norbnow
          vec(lam, 2*n2d+iorb)=-funl1(lam, iorb)
          vec(n2d+lam, 2*n2d+iorb)=-fundl1(lam, iorb)
        enddo
      enddo                      
      numb=numb+1
    endif
    call mpi_allreduce(amat, amat_aux, 2*2*n2d*2*n2d, MPI_REAL8,     &
                       MPI_SUM, new_comm, info)
    call mpi_allreduce(vec, vec_aux, 2*2*n2d*(2*n2d+norbnow),        &
                       MPI_REAL8, MPI_SUM, new_comm, info)
    call DCOPY(2*2*n2d*2*n2d, amat_aux, 1, amat, 1)  
    call DCOPY(2*2*n2d*(2*n2d+norbnow), vec_aux, 1, vec, 1) 
    call ZGESV(2*n2d, 2*n2d+norbnow, amat, 2*n2d, ipiv,              &
               vec, 2*n2d, info)          
!
!   recalculate the functions for CPU which is left to matching 
!   boundary
!
    if(numb.le.1.and.me/2*2.eq.me) then
     do ig=1, n2d
      do n=1, n2d
       do lam=1, n2d
        fun1(ig, n)=fun1(ig, n)+  &
                           vec(n2d+lam, n)*fun1(ig, n2d+lam)
        fund1(ig, n)=fund1(ig, n)+  &
                           vec(n2d+lam, n)*fund1(ig, n2d+lam)
       enddo
      enddo
      do iorb=1, norbnow
       do lam=1, n2d
        funl1(ig, iorb)=funl1(ig, iorb)+  &
                vec(n2d+lam, 2*n2d+iorb)*fun1(ig, n2d+lam)
        fundl1(ig, iorb)=fundl1(ig, iorb)+  &
                vec(n2d+lam, 2*n2d+iorb)*fund1(ig, n2d+lam)
       enddo
      enddo
     enddo
     do ig=1, n2d
       x=(0.d0,0.d0)
       y=(0.d0,0.d0)
       do n=1, n2d
         do lam=1, n2d
           x(n)=x(n)+vec(n2d+lam, n2d+n)*fun1(ig, n2d+lam)
           y(n)=y(n)+vec(n2d+lam, n2d+n)*fund1(ig, n2d+lam)
         enddo
       enddo
       do n=1, n2d
         fun1(ig, n2d+n)=x(n)
         fund1(ig, n2d+n)=y(n)
       enddo
     enddo                                   
    endif
!
!   recalculate the functions for CPU which is right to matching
!   boundary  
!
    if(numb.le.1.and.me/2*2.ne.me) then
     do ig=1, n2d
      do n=1, n2d
       do lam=1, n2d
        fun0(ig, n2d+n)=fun0(ig, n2d+n)+  &
                           vec(lam, n2d+n)*fun0(ig, lam)
        fund0(ig, n2d+n)=fund0(ig, n2d+n)+  &
                           vec(lam, n2d+n)*fund0(ig, lam)
       enddo
      enddo
      do iorb=1, norbnow
       do lam=1, n2d
        funl0(ig, iorb)=funl0(ig, iorb)+  &
                     vec(lam, 2*n2d+iorb)*fun0(ig, lam)
        fundl0(ig, iorb)=fundl0(ig, iorb)+  &
                     vec(lam, 2*n2d+iorb)*fund0(ig, lam)
       enddo
      enddo
     enddo
     do ig=1, n2d
       x=(0.d0,0.d0)
       y=(0.d0,0.d0)
       do n=1, n2d
         do lam=1, n2d
           x(n)=x(n)+vec(lam, n)*fun0(ig, lam)
           y(n)=y(n)+vec(lam, n)*fund0(ig, lam)
         enddo
       enddo
       do n=1, n2d
         fun0(ig, n)=x(n)
         fund0(ig, n)=y(n)
       enddo
     enddo                                             
    endif  
!
!  to recalculate the integrals for a given group of CPU
!
    if(me.ge.ib) then
     do iorb=1, norbnow
      do iorb1=1, norbnow
       do lam=1, n2d
        intw2(iorb,iorb1)=intw2(iorb,iorb1)+           &
           vec(n2d+lam, 2*n2d+iorb1)*intw1(iorb, n2d+lam)
       enddo
      enddo
     enddo
     do iorb=1, norbnow
       x=(0.d0,0.d0)
       do n=1, n2d
         do lam=1, n2d
           x(n)=x(n)+vec(n2d+lam, n2d+n)*intw1(iorb, n2d+lam)
           intw1(iorb, n)=intw1(iorb, n)+   &
               vec(n2d+lam, n)*intw1(iorb, n2d+lam)
         enddo
       enddo
       do n=1, n2d
         intw1(iorb, n2d+n)=x(n)
       enddo
     enddo                            
    else
     do iorb=1, norbnow
      do iorb1=1, norbnow
       do lam=1, n2d
        intw2(iorb, iorb1)=intw2(iorb, iorb1)+           &
             vec(lam, 2*n2d+iorb1)*intw1(iorb, lam)
       enddo
      enddo
     enddo
     do iorb=1, norbnow
       x=(0.d0,0.d0)
       do n=1, n2d
         do lam=1, n2d
           x(n)=x(n)+vec(lam, n)*intw1(iorb, lam)
           intw1(iorb, n2d+n)=intw1(iorb, n2d+n)+   &
               vec(lam, n2d+n)*intw1(iorb, lam)
         enddo
       enddo
       do n=1, n2d
         intw1(iorb, n)=x(n)
       enddo
     enddo                                       
    endif 

!
!  to next matching step
!
    n=(ib+ninsl-1)/(2*ninsl)
    if(n/2*2.eq.n) then
      ib=ib-ninsl
    else
      ib=ib+ninsl
    endif
    ninsl=ninsl*2
    ibound=ibound/2              

    call mpi_comm_free(new_comm, info) 
  enddo

!
! Broadcast of the functions and the integrals to all CPU
!
  call mpi_bcast(fun0, 2*n2d*2*n2d, MPI_REAL8,  0,          &
                 MPI_COMM_WORLD, info)
  call mpi_bcast(fund0, 2*n2d*2*n2d, MPI_REAL8, 0,          &
                 MPI_COMM_WORLD, info)   
  call mpi_bcast(funl0, 2*n2d*norbf, MPI_REAL8,  0,         &
                 MPI_COMM_WORLD, info)   
  call mpi_bcast(fundl0, 2*n2d*norbf, MPI_REAL8, 0,         &
                 MPI_COMM_WORLD, info)   

  call mpi_bcast(fun1, 2*n2d*2*n2d, MPI_REAL8,  nproc-1,    &
                 MPI_COMM_WORLD, info)
  call mpi_bcast(fund1, 2*n2d*2*n2d, MPI_REAL8, nproc-1,    &
                 MPI_COMM_WORLD, info)
  call mpi_bcast(funl1, 2*n2d*norbf, MPI_REAL8,  nproc-1,   &
                 MPI_COMM_WORLD, info)
  call mpi_bcast(fundl1, 2*n2d*norbf, MPI_REAL8, nproc-1,   &
                 MPI_COMM_WORLD, info)                         

!
! Gathering of the integrals
!
  call reduce(2*norbf*2*n2d, intw1)
  call reduce(2*norbf*norbf, intw2)

  deallocate(x)
  deallocate(y)
  deallocate(amat)
  deallocate(amat_aux)
  deallocate(vec)
  deallocate(vec_aux)
  deallocate(ipiv)
     
#endif
  return
end subroutine rotproc


