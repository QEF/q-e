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
SUBROUTINE rotproc (fun0, fund0, fun1, fund1, funl0, fundl0, funl1,  &
                   fundl1, intw1, intw2, n2d, norbf, norb, nrzp)
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

  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE parallel_include
  USE mp_world,         ONLY : nproc
  USE mp_pools,         ONLY : me_pool, intra_pool_comm
  USE mp,               ONLY : mp_sum
  use cond,             ONLY : lorb, funz0



  IMPLICIT NONE


  INTEGER :: ig, n, lam, lam1, iorb, iorb1, norbf, norb, n2d,  &
             ibound, numb, ninsl, ib, icolor, ikey, new_comm, nrzp, info
  INTEGER, ALLOCATABLE :: ipiv(:)
  COMPLEX(DP), PARAMETER :: one=(1.d0, 0.d0), zero=(0.d0,0.d0)
  COMPLEX(DP) :: fun0(n2d, 2*n2d),    & ! phi_n(0)
                      fund0(n2d, 2*n2d),   & ! phi'_n(0)
                      fun1(n2d, 2*n2d),    & ! phi_n(d)
                      fund1(n2d, 2*n2d),   & ! phi'_n(d)
                      funl0(n2d, npol*norbf), & ! phi_alpha(0)
                      fundl0(n2d, npol*norbf),& ! phi'_alpha(0)
                      funl1(n2d, npol*norbf), & ! phi_alpha(d)
                      fundl1(n2d, npol*norbf),  & ! phi'_alpha(d)
                      intw1(norbf*npol, 2*n2d), & ! integrals on phi_n
                      intw2(norbf*npol, norbf*npol) ! integrals on phi_alpha

  COMPLEX(DP), ALLOCATABLE :: x(:), y(:), amat(:,:), vec(:,:),  &
                               amat_aux(:,:), vec_aux(:,:)

#if defined(__MPI)

  IF(nproc.EQ.1) RETURN

  ALLOCATE( x( n2d ) )
  ALLOCATE( y( n2d ) )
  ALLOCATE( amat( 2*n2d, 2*n2d ) )
  ALLOCATE( amat_aux( 2*n2d, 2*n2d ) )
  ALLOCATE( vec( 2*n2d, 2*n2d+npol*norb ) )
  ALLOCATE( vec_aux( 2*n2d, 2*n2d+npol*norb ) )
  ALLOCATE( ipiv( 2*n2d ) )

  numb=0
  ibound=nproc/2
  ninsl=1
  ib=2*(me_pool+1)-(me_pool+1)/2*2

  DO WHILE(ibound.GT.0)

!
!   To find the matching coefficients for a group of CPU
!
    icolor=(ib+ninsl-1)/(2*ninsl)
    ikey=((me_pool+1)+ninsl)-ib
    CALL mpi_barrier (MPI_COMM_WORLD, info)
    CALL mpi_comm_split(MPI_COMM_WORLD, icolor, ikey, new_comm, info)
    amat=(0.d0,0.d0)
    vec=(0.d0,0.d0)

    IF((me_pool+1).EQ.ib) THEN
      DO lam=1, n2d
        DO lam1=1, n2d
          amat(lam, n2d+lam1)=-fun0(lam,n2d+lam1)
          amat(n2d+lam, n2d+lam1)=-fund0(lam,n2d+lam1)
          vec(lam, lam1)=fun0(lam, lam1)
          vec(n2d+lam, lam1)=fund0(lam, lam1)
        ENDDO
        DO iorb=1, npol*norb
          vec(lam, 2*n2d+iorb)=funl0(lam, iorb)
          vec(n2d+lam, 2*n2d+iorb)=fundl0(lam, iorb)
        ENDDO
      ENDDO
      numb=numb+1
    ENDIF
    IF((me_pool+1).EQ.ib-1) THEN
      DO lam=1, n2d
        DO lam1=1, n2d
          amat(lam, lam1)=fun1(lam,lam1)
          amat(n2d+lam, lam1)=fund1(lam,lam1)
          vec(lam, n2d+lam1)=-fun1(lam, n2d+lam1)
          vec(n2d+lam, n2d+lam1)=-fund1(lam, n2d+lam1)
        ENDDO
        DO iorb=1, npol*norb
          vec(lam, 2*n2d+iorb)=-funl1(lam, iorb)
          vec(n2d+lam, 2*n2d+iorb)=-fundl1(lam, iorb)
        ENDDO
      ENDDO
      numb=numb+1
    ENDIF
    CALL mpi_allreduce(amat, amat_aux, 2*2*n2d*2*n2d, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, new_comm, info)
    CALL mpi_allreduce(vec, vec_aux, 2*2*n2d*(2*n2d+npol*norb),        &
                       MPI_DOUBLE_PRECISION, MPI_SUM, new_comm, info)
    CALL dcopy(2*2*n2d*2*n2d, amat_aux, 1, amat, 1)
    CALL dcopy(2*2*n2d*(2*n2d+npol*norb), vec_aux, 1, vec, 1)
    CALL ZGESV(2*n2d, 2*n2d+npol*norb, amat, 2*n2d, ipiv,              &
               vec, 2*n2d, info)
!
!   recalculate the functions for CPU which is left to matching
!   boundary
!
    IF(numb.LE.1.AND.(me_pool+1)/2*2.EQ.(me_pool+1)) THEN
     DO ig=1, n2d
      DO n=1, n2d
       DO lam=1, n2d
        fun1(ig, n)=fun1(ig, n)+  &
                           vec(n2d+lam, n)*fun1(ig, n2d+lam)
        fund1(ig, n)=fund1(ig, n)+  &
                           vec(n2d+lam, n)*fund1(ig, n2d+lam)
       ENDDO
      ENDDO
      DO iorb=1, npol*norb
       DO lam=1, n2d
        funl1(ig, iorb)=funl1(ig, iorb)+  &
                vec(n2d+lam, 2*n2d+iorb)*fun1(ig, n2d+lam)
        fundl1(ig, iorb)=fundl1(ig, iorb)+  &
                vec(n2d+lam, 2*n2d+iorb)*fund1(ig, n2d+lam)
       ENDDO
      ENDDO
     ENDDO
     DO ig=1, n2d
       x=(0.d0,0.d0)
       y=(0.d0,0.d0)
       DO n=1, n2d
         DO lam=1, n2d
           x(n)=x(n)+vec(n2d+lam, n2d+n)*fun1(ig, n2d+lam)
           y(n)=y(n)+vec(n2d+lam, n2d+n)*fund1(ig, n2d+lam)
         ENDDO
       ENDDO
       DO n=1, n2d
         fun1(ig, n2d+n)=x(n)
         fund1(ig, n2d+n)=y(n)
       ENDDO
     ENDDO
    ENDIF
!
!   recalculate the functions for CPU which is right to matching
!   boundary
!
    IF(numb.LE.1.AND.(me_pool+1)/2*2.NE.(me_pool+1)) THEN
     DO ig=1, n2d
      DO n=1, n2d
       DO lam=1, n2d
        fun0(ig, n2d+n)=fun0(ig, n2d+n)+  &
                           vec(lam, n2d+n)*fun0(ig, lam)
        fund0(ig, n2d+n)=fund0(ig, n2d+n)+  &
                           vec(lam, n2d+n)*fund0(ig, lam)
       ENDDO
      ENDDO
      DO iorb=1, npol*norb
       DO lam=1, n2d
        funl0(ig, iorb)=funl0(ig, iorb)+  &
                     vec(lam, 2*n2d+iorb)*fun0(ig, lam)
        fundl0(ig, iorb)=fundl0(ig, iorb)+  &
                     vec(lam, 2*n2d+iorb)*fund0(ig, lam)
       ENDDO
      ENDDO
     ENDDO
     DO ig=1, n2d
       x=(0.d0,0.d0)
       y=(0.d0,0.d0)
       DO n=1, n2d
         DO lam=1, n2d
           x(n)=x(n)+vec(lam, n)*fun0(ig, lam)
           y(n)=y(n)+vec(lam, n)*fund0(ig, lam)
         ENDDO
       ENDDO
       DO n=1, n2d
         fun0(ig, n)=x(n)
         fund0(ig, n)=y(n)
       ENDDO
     ENDDO
    ENDIF
!
!  to recalculate the integrals for a given group of CPU
!
    IF((me_pool+1).GE.ib) THEN
     DO iorb=1, npol*norb
      DO iorb1=1, npol*norb
       DO lam=1, n2d
        intw2(iorb,iorb1)=intw2(iorb,iorb1)+           &
           vec(n2d+lam, 2*n2d+iorb1)*intw1(iorb, n2d+lam)
       ENDDO
      ENDDO
     ENDDO
     DO iorb=1, npol*norb
       x=(0.d0,0.d0)
       DO n=1, n2d
         DO lam=1, n2d
           x(n)=x(n)+vec(n2d+lam, n2d+n)*intw1(iorb, n2d+lam)
           intw1(iorb, n)=intw1(iorb, n)+   &
               vec(n2d+lam, n)*intw1(iorb, n2d+lam)
         ENDDO
       ENDDO
       DO n=1, n2d
         intw1(iorb, n2d+n)=x(n)
       ENDDO
     ENDDO

     IF (lorb) THEN
        DO n = 1, nrzp
          CALL zgemm('n','n',n2d,n2d,n2d,one,funz0(1,n2d+1,n),n2d,&
                     vec(n2d+1,1),2*n2d,one,funz0(1,1,n),n2d)
          CALL zgemm('n','n',n2d,npol*norb,n2d,one,funz0(1,n2d+1,n),n2d,&
                     vec(n2d+1,2*n2d+1),2*n2d,one,funz0(1,2*n2d+1,n),n2d)
          CALL zgemm('n','n',n2d,n2d,n2d,one,funz0(1,n2d+1,n),n2d,&
                      vec(n2d+1,n2d+1),2*n2d,zero,vec_aux(1,1),2*n2d)
          do ig = 1, n2d
            do lam = 1, n2d
              funz0(ig,n2d+lam,n) = vec_aux(ig,lam)
            enddo
          enddo
        END DO
     ENDIF

    ELSE
     DO iorb=1, npol*norb
      DO iorb1=1, npol*norb
       DO lam=1, n2d
        intw2(iorb, iorb1)=intw2(iorb, iorb1)+           &
             vec(lam, 2*n2d+iorb1)*intw1(iorb, lam)
       ENDDO
      ENDDO
     ENDDO
     DO iorb=1, npol*norb
       x=(0.d0,0.d0)
       DO n=1, n2d
         DO lam=1, n2d
           x(n)=x(n)+vec(lam, n)*intw1(iorb, lam)
           intw1(iorb, n2d+n)=intw1(iorb, n2d+n)+   &
               vec(lam, n2d+n)*intw1(iorb, lam)
         ENDDO
       ENDDO
       DO n=1, n2d
         intw1(iorb, n)=x(n)
       ENDDO
     ENDDO

     IF (lorb) THEN
        DO n = 1, nrzp
          CALL zgemm('n','n',n2d,n2d,n2d,one,funz0(1,1,n),n2d,&
                     vec(1,n2d+1),2*n2d,one,funz0(1,n2d+1,n),n2d)
          CALL zgemm('n','n',n2d,npol*norb,n2d,one,funz0(1,1,n),n2d,&
                     vec(1,2*n2d+1),2*n2d,one,funz0(1,2*n2d+1,n),n2d)
          CALL zgemm('n','n',n2d,n2d,n2d,one,funz0(1,1,n),n2d,&
                      vec(1,1),2*n2d,zero,vec_aux(1,1),2*n2d)
          do ig = 1, n2d
            do lam = 1, n2d
              funz0(ig,lam,n) = vec_aux(ig,lam)
            enddo
          enddo
        END DO
     ENDIF

    ENDIF

!
!  to next matching step
!
    n=(ib+ninsl-1)/(2*ninsl)
    IF(n/2*2.EQ.n) THEN
      ib=ib-ninsl
    ELSE
      ib=ib+ninsl
    ENDIF
    ninsl=ninsl*2
    ibound=ibound/2

    CALL mpi_comm_free(new_comm, info)
  ENDDO

!
! Broadcast of the functions and the integrals to all CPU
!
  CALL mpi_bcast(fun0, 2*n2d*2*n2d, MPI_DOUBLE_PRECISION,  0,          &
                 MPI_COMM_WORLD, info)
  CALL mpi_bcast(fund0, 2*n2d*2*n2d, MPI_DOUBLE_PRECISION, 0,          &
                 MPI_COMM_WORLD, info)
  CALL mpi_bcast(funl0, 2*n2d*npol*norbf, MPI_DOUBLE_PRECISION,  0,    &
                 MPI_COMM_WORLD, info)
  CALL mpi_bcast(fundl0, 2*n2d*npol*norbf, MPI_DOUBLE_PRECISION, 0,    &
                 MPI_COMM_WORLD, info)

  CALL mpi_bcast(fun1, 2*n2d*2*n2d, MPI_DOUBLE_PRECISION,  nproc-1,    &
                 MPI_COMM_WORLD, info)
  CALL mpi_bcast(fund1, 2*n2d*2*n2d, MPI_DOUBLE_PRECISION, nproc-1,    &
                 MPI_COMM_WORLD, info)
  CALL mpi_bcast(funl1, 2*n2d*npol*norbf, MPI_DOUBLE_PRECISION,  nproc-1,   &
                 MPI_COMM_WORLD, info)
  CALL mpi_bcast(fundl1, 2*n2d*npol*norbf, MPI_DOUBLE_PRECISION, nproc-1,   &
                 MPI_COMM_WORLD, info)

!
! Gathering of the integrals
!
  CALL mp_sum( intw1, intra_pool_comm )
  CALL mp_sum( intw2, intra_pool_comm )

  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(amat)
  DEALLOCATE(amat_aux)
  DEALLOCATE(vec)
  DEALLOCATE(vec_aux)
  DEALLOCATE(ipiv)

#endif
  RETURN
END SUBROUTINE rotproc
