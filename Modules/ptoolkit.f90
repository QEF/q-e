!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!==----------------------------------------------==!
    MODULE parallel_toolkit
!==----------------------------------------------==!

    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE parallel_include

    IMPLICIT NONE
    SAVE
    PRIVATE

    PUBLIC :: pdspev_drv, dspev_drv, &
              diagonalize, pzhpev_drv, zhpev_drv, cdiagonalize
    PUBLIC :: rep_matmul_drv

!==----------------------------------------------==!
    CONTAINS
!==----------------------------------------------==!

    SUBROUTINE ptredv( a, lda, d, e, v, ldv, nrl, n, nproc, me, comm )

!
!     Parallel version of the famous HOUSEHOLDER tridiagonalization
!     Algorithm for simmetric matrix.
!
!     AUTHOR : Carlo Cavazzoni - SISSA 1997
!              comments and suggestions to : carlo.cavazzoni@cineca.it
!
! REFERENCES :
!
! NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
! W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
! CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.  
!
! PARALLEL NUMERICAL ALGORITHMS,
! T.L. FREEMAN AND C.PHILLIPS,
! PRENTICE HALL INTERNATIONAL (1992). 
!
!
!
!     INPUTS :
!
!     A(NRL,N) Local part of the global matrix A(N,N) to be reduced,
!              only the upper triangle is needed.
!              The rows of the matrix are distributed amoung processors
!              with blocking factor 1. 
!              Example for NPROC = 4 :
!              ROW | PE
!              1   | 0
!              2   | 1
!              3   | 2
!              4   | 3
!              5   | 0
!              6   | 1
!              ..  | ..
!
!     LDA      LEADING DIMENSION OF MATRIX A.
!
!     LDV      LEADING DIMENSION OF MATRIX V.
!
!     NRL      NUMBER OF ROWS BELONGING TO THE LOCAL PROCESSOR.
!
!     N        DIMENSION OF THE GLOBAL MATRIX.
!
!     NPROC    NUMBER OF PROCESSORS.
!
!     ME       INDEX OF THE LOCAL PROCESSOR (Starting from 0).
!
!
!     OUTPUTS :
!
!     V(NRL,N) Ortogonal transormation that tridiagonalize A,
!              this matrix is distributed omong processor
!              in the same way as A.
!
!     D(N)     Diagonal elements of the tridiagonal matrix
!              this vector is equal on all processors.
!
!     E(N)     Subdiagonal elements of the tridiagonal matrix
!              this vector is equal on all processors. 
!
!

      USE kinds
      IMPLICIT NONE

      INTEGER, PARAMETER :: PTRED_WORK_SIZE = 2333

#if defined __PARA
#  if defined __SHMEM
      include "mpp/shmem.fh"
      INTEGER PSYNC_BC(SHMEM_BCAST_SYNC_SIZE)
      INTEGER PSYNC_B(SHMEM_BARRIER_SYNC_SIZE)
      INTEGER PSYNC_STA(SHMEM_REDUCE_SYNC_SIZE)
      REAL(DP)  pWrk(MAX(PTRED_WORK_SIZE,SHMEM_REDUCE_MIN_WRKDATA_SIZE))
      DATA PSYNC_BC /SHMEM_BCAST_SYNC_SIZE*SHMEM_SYNC_VALUE/
      DATA PSYNC_B /SHMEM_BARRIER_SYNC_SIZE*SHMEM_SYNC_VALUE/
      DATA PSYNC_STA /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/
      SAVE PSYNC_B, PSYNC_STA, PSYNC_BC,pWrk
#  endif 
#endif

      INTEGER :: N, NRL, LDA, LDV
      INTEGER :: NPROC, ME, comm
      REAL(DP) :: A(LDA,N), D(N), E(N), V(LDV,N)
!
      REAL(DP) :: DDOT
!
      REAL(DP) :: g, scale, sigma, kappa, f, h, tmp
      REAL(DP) :: u(PTRED_WORK_SIZE)
      REAL(DP) :: p(PTRED_WORK_SIZE)
      REAL(DP) :: vtmp(PTRED_WORK_SIZE)

      REAL(DP) :: tu, tp, one_over_h
      REAL(DP) :: one_over_scale
      REAL(DP) :: ul(PTRED_WORK_SIZE)
      REAL(DP) :: pl(PTRED_WORK_SIZE)
#if (defined __SHMEM && defined __ALTIX) || (defined __SHMEM && defined __ORIGIN)
      REAL(DP), SAVE :: d_tmp(PTRED_WORK_SIZE)
#endif
      integer :: l, i, j, k, t, tl, ierr
      integer :: kl, jl, ks, lloc
      integer :: is(PTRED_WORK_SIZE)
      integer :: ri(PTRED_WORK_SIZE)
     
      save :: g,h,scale,sigma,kappa,f,u,p,tmp,vtmp
      save :: tu,tp,one_over_h,one_over_scale,ul,pl
      save :: l,i,j,k,t,tl,ierr,kl,jl,is,ri,ks,lloc

      
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........

!      ttot = 0
!      mem1 = 0

      IF( N == 0 ) THEN
        RETURN
      ELSE IF( (N > PTRED_WORK_SIZE) .OR. (N < 0) ) THEN
        WRITE( stdout,*) ' *** ERROR IN PTREDV'
        WRITE( stdout,*) ' N OUT OF RANGE : ',N
        STOP
      END IF
      
      DO I = N, 1, -1
        IS(I)  = (I-1)/NPROC
        RI(I)  = MOD((I-1),NPROC)
        IF(ME .le. RI(I) ) then
          IS(I) = IS(I) + 1
        END IF
!        WRITE( stdout,100) 'ME,I,RI,IS',ME,I,RI(I),IS(I)
100   FORMAT(A,2X,5I4)
      END DO

      DO I = N, 2, -1 

         L     = I - 1         ! primo elemeto
         H     = 0.0d0

         IF ( L .GT. 1 ) THEN

           SCALE = 0.0D0
           DO K = 1, is(l)
              SCALE = SCALE + DABS( A(K,I) )
           END DO

#if defined __PARA
#  if defined __SHMEM
           call shmem_barrier_all
           CALL SHMEM_REAL8_SUM_TO_ALL(SCALE, SCALE, 1, 0, 0, nproc, pWrk, pSync_sta)
#  elif defined __MPI
           CALL MPI_ALLREDUCE(SCALE, TMP, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
           SCALE = TMP
#  endif
#endif

           IF ( SCALE .EQ. 0.0D0 )  THEN
             IF (RI(L).EQ.ME) THEN
               E(I) = A(is(L),I) 
             END IF
           ELSE 

!  ......  CALCOLO DI SIGMA E DI H

             ONE_OVER_SCALE = 1.0d0/SCALE
             SIGMA = 0.0D0
             DO k = 1,is(L)
               A(k,I) = A(k,I) * ONE_OVER_SCALE
               SIGMA  = SIGMA + A(k,I)**2 
             END DO

             IF(RI(L) .eq. ME ) THEN
               F = A(is(L),I)
             END IF

#if defined __PARA
#  if defined __SHMEM
             call shmem_barrier_all
             CALL SHMEM_REAL8_SUM_TO_ALL(SIGMA, SIGMA, 1, 0, 0, nproc, pWrk, pSync_sta)
             call shmem_barrier_all
#    if defined __ALTIX || defined __ORIGIN
             CALL SHMEM_BROADCAST8(F,F,1,RI(L),0,0,nproc,pSync_bc)
#    else
             CALL SHMEM_BROADCAST(F,F,1,RI(L),0,0,nproc,pSync_bc)
#    endif
#  elif defined __MPI
             CALL MPI_ALLREDUCE(SIGMA, TMP, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
             SIGMA = TMP
             CALL MPI_BCAST(F, 1, MPI_DOUBLE_PRECISION, RI(L), comm, IERR)
#  endif
#endif

             G          = -SIGN(SQRT(SIGMA),F)
             H          = SIGMA - F*G
             ONE_OVER_H = 1.0d0/H
             E(I)       = SCALE*G

!  ......  COSTRUZIONE DEL VETTORE U

#if defined __T3E
!DIR$ UNROLL 8
#endif
             DO k = 1,L          
               vtmp(k) = 0.0d0 
             END DO

             k = ME + 1

#if defined __T3E
!DIR$ UNROLL 4
#endif
             DO kl = 1,is(l)          
               vtmp(k)   = A(kl,I)
               UL(kl) = A(kl,I)
               k      = k + NPROC
             END DO
             IF(RI(L) .eq. ME ) THEN
               vtmp(L)      = F - G 
               UL(is(l))  = F - G
               A(is(l),I) = F - G 
             END IF

#if defined __PARA
#  if defined __SHMEM
             call shmem_barrier_all
             CALL SHMEM_REAL8_SUM_TO_ALL(VTMP, U, L, 0, 0, nproc, pWrk, pSync_sta)
#  elif defined __MPI
             CALL MPI_ALLREDUCE(VTMP,U,L,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
             u(1:l) = vtmp(1:l)
#endif

!  ......  COSTRUZIONE DEL VETTORE P

             KAPPA = 0.0d0

             DO J = 1,L

               vtmp(j) = 0.0d0

#if defined __T3E
!DIR$ UNROLL 8
#endif
               DO KL = 1, IS(J)
                 vtmp(J) = vtmp(J) + A(KL,J) * UL(KL)
               END DO

               IF(L.GT.J .AND. ME.EQ.RI(J)) then

#if defined __T3E
!DIR$ UNROLL 8
#endif
                 DO K = J+1,L
                   vtmp(J) = vtmp(J) + A(IS(J),K) * U(K)
                 END DO
               END IF

               vtmp(J) = vtmp(J) * ONE_OVER_H
               KAPPA = KAPPA + vtmp(J) * U(J) * 0.5d0 * ONE_OVER_H
             END DO 

#if defined __PARA
#  if defined __SHMEM
             call shmem_barrier_all
             CALL SHMEM_REAL8_SUM_TO_ALL(KAPPA, KAPPA, 1, 0, 0, nproc, pWrk, pSync_sta)
             call shmem_barrier_all
             CALL SHMEM_REAL8_SUM_TO_ALL(vtmp, P, L, 0, 0, nproc, pWrk, pSync_sta)
             call shmem_barrier_all
#  elif defined __MPI
             CALL MPI_ALLREDUCE(KAPPA,TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM, comm,IERR)
             KAPPA = TMP
             CALL MPI_ALLREDUCE(vtmp,p,L,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
             p(1:l) = vtmp(1:l)
#endif

             CALL DAXPY(l, -kappa, u, 1, p, 1)

             k = me + 1
#if defined __T3E
!DIR$ UNROLL 8
#endif
             DO kl = 1,is(l)          
               PL(kl) = P(k)
               k      = k + NPROC
             END DO

             DO J = 1,L
               tu= U(J)
               tp= P(J)
#if defined __T3E
!DIR$ UNROLL 12
#endif
               DO KL = 1,is(l)
                 A(KL,j) = A(KL,j) - UL(KL) * tp 
                 A(KL,j) = A(KL,j) - PL(KL) * tu 
               END DO
             END DO
           END IF

         ELSE 

           IF(RI(L).EQ.ME) THEN
             G = A(is(l),I)
           END IF

#if defined __PARA
#  if defined __SHMEM
           call shmem_barrier_all
#    if defined __ALTIX || defined __ORIGIN
           CALL SHMEM_BROADCAST8(G,G,1,RI(L),0,0,nproc,pSync_bc)
#    else
           CALL SHMEM_BROADCAST(G,G,1,RI(L),0,0,nproc,pSync_bc)
#    endif
#  elif defined __MPI
           CALL MPI_BCAST(G,1,MPI_DOUBLE_PRECISION,RI(L),comm,IERR)
#  endif
#endif
           E(I) = G

         END IF

         D(I) = H

#if defined __PARA
#  if defined __SHMEM
         call shmem_barrier_all
#  endif
#endif

      END DO
      
      E(1) = 0.0d0
      D(1) = 0.0d0


!      t1 = mclock()

      DO J = 1,N
#if defined __T3E
!DIR$ UNROLL 12
#endif
        DO I = 1,NRL
          V(I,J) = 0.0d0
        END DO
        IF(RI(J).EQ.ME) THEN
          V(IS(J),J) = 1.0d0
        END IF
      END DO

      DO I = 2,N
        L = I - 1
        LLOC = IS(L)
        IF(D(I).NE.0.0d0) THEN
          ONE_OVER_H = 1.0d0/D(I)
          DO J = 1, L
            P(J) = DDOT(LLOC,V(1,J),1,A(1,I),1)
          END DO

#if defined __PARA
#  if defined __SHMEM
          call shmem_barrier_all
          CALL SHMEM_REAL8_SUM_TO_ALL(P, VTMP, L, 0, 0, nproc, pWrk, pSync_sta)
#  elif defined __MPI
          CALL MPI_ALLREDUCE(P,VTMP,L,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
          vtmp(1:l) = p(1:l)
#endif

          DO J = 1, L
            call DAXPY(LLOC,-VTMP(J)*ONE_OVER_H,A(1,I),1,V(1,J),1)
          END DO

        END IF

#if defined __PARA
#  if defined __SHMEM
         call shmem_barrier_all
#  endif
#endif

      END DO 


      DO I = 1,N
        U(I) = 0.0d0
        IF(RI(I).eq.ME) then
          U(I) = A(IS(I),I) 
        END IF
      END DO

#if defined __PARA
#  if defined __SHMEM
#if defined __ALTIX || defined __ORIGIN
      d_tmp(1:n) = d(1:n)
      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(U, d_tmp, N, 0, 0, nproc, pWrk,       &
     &                            pSync_sta)
#else
      call shmem_barrier_all
      CALL SHMEM_REAL8_SUM_TO_ALL(U, D, N, 0, 0, nproc, pWrk, pSync_sta)
#endif
      call shmem_barrier_all
#  elif defined __MPI
      CALL MPI_ALLREDUCE(U,D,N,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
      D(1:N) = U(1:N)
#endif

    RETURN
    END SUBROUTINE ptredv

!==----------------------------------------------==!


    SUBROUTINE ptqliv( d, e, n, z, ldz, nrl )

!
! Modified QL algorithm for CRAY T3E PARALLEL MACHINE
! calculate the eigenvectors and eigenvalues of a matrix reduced to 
! tridiagonal form by PTREDV. 
!
! AUTHOR : Carlo Cavazzoni - SISSA 1997
!          comments and suggestions to : cava@sissa.it
!
! REFERENCES :
!
! NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
! W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
! CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.  
!
! PARALLEL NUMERICAL ALGORITHMS,
! T.L. FREEMAN AND C.PHILLIPS,
! PRENTICE HALL INTERNATIONAL (1992). 
!
! NOTE : the algorithm tha finds the eigenvalues is not parallelized
!        ( it scales as O(N^2) ), I prefere to parallelize only the
!        updating of the eigenvectors because its the most costly
!        part of the algorithm ( it scales as O(N^3) ).
!        For large matrix in practice all the time is spent in the updating
!        that in this routine scales linearly with the number of processors,
!        in fact there is no communication at all.
!
!  
!     INPUTS :
!
!     D(N)     Diagonal elements of the tridiagonal matrix
!              this vector is equal on all processors.
!
!     E(N)     Subdiagonal elements of the tridiagonal matrix
!              this vector is equal on all processors.
!
!     N        DIMENSION OF THE GLOBAL MATRIX.
!
!     NRL      NUMBER OF ROWS OF Z BELONGING TO THE LOCAL PROCESSOR.
!
!     LDZ      LEADING DIMENSION OF MATRIX Z. 
!
!     Z(LDZ,N) Ortogonal transormation that tridiagonalize the original
!              matrix A.
!              The rows of the matrix are distributed amoung processors
!              with blocking factor 1.
!              Example for NPROC = 4 :
!              ROW | PE
!              1   | 0
!              2   | 1
!              3   | 2
!              4   | 3
!              5   | 0
!              6   | 1
!              ..  | ..
!
!
!
!     OUTPUTS :
!
!     Z(LDZ,N) EIGENVECTORS OF THE ORIGINAL MATRIX.
!              THE Jth COLUMN of Z contains the eigenvectors associated
!              with the jth eigenvalue.
!              The eigenvectors are scattered among processors (4PE examp. )
!              eigenvector | PE
!               elements   |
!                 V(1)     | 0
!                 V(2)     | 1
!                 V(3)     | 2
!                 V(4)     | 3
!                 V(5)     | 0
!                 V(6)     | 1
!                 ....       ..
! 
!     D(N)     Eigenvalues of the original matrix, 
!              this vector is equal on all processors.
!
!
!
!
      USE kinds
      IMPLICIT NONE

      INTEGER n,nrl,ldz
      REAL(DP) d(n),e(n)
      REAL(DP) z(ldz,n)

      INTEGER FV_WORK_SIZE
      INTEGER CSV_WORK_SIZE
      PARAMETER ( FV_WORK_SIZE = 2000 )
      PARAMETER ( CSV_WORK_SIZE = 2000 )

      INTEGER i,iter,mk,k,l,m
      REAL(DP) b,dd,f,g,p,r,c,s
      REAL(DP) cv(CSV_WORK_SIZE)
      REAL(DP) sv(CSV_WORK_SIZE)
      REAL(DP) fv1(FV_WORK_SIZE)
      REAL(DP) fv2(FV_WORK_SIZE)


      save cv,sv,fv1,fv2,b,dd,f,g,p,r,c,s
      save i,iter,mk,k,l,m


      if(n.gt.CSV_WORK_SIZE) then
        print *,' ptqli CSV_WORK_SIZE too small '
        stop
      end if
      if(nrl.gt.FV_WORK_SIZE) then
        print *,' ptqli FV_WORK_SIZE too small '
        stop
      end if

      do l = 2,n
        e(l-1) = e(l)
      end do

      do 15 l=1,n
        iter=0
1       do m=l,n-1
          dd = abs(d(m))+abs(d(m+1))
          if ( abs(e(m))+dd .eq. dd ) goto 2
        end do
        m=n

2       if(m.ne.l)then
          if(iter.eq.200) then
            call errore(' tqli ',' too many iterations ', iter)
          end if
          iter=iter+1
          g=(d(l+1)-d(l))/(2.0d0*e(l))
          r=pythag(g,1.0d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0d0
          c=1.0d0
          p=0.0d0
          do i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.0d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.0d0
              goto 1
            endif
            c=g/r
            g=d(i+1)-p
            s=f/r
            r=(d(i)-g)*s+2.0d0*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b

!            do k=1,nrl
!              f=z(k,i+1)
!              z(k,i+1)=s*z(k,i)+c*f
!              z(k,i)=c*z(k,i)-s*f
!            end do

            cv(i) = c
            sv(i) = s

          end do

          do i=m-1,l,-1
#if defined __T3E
!DIR$ UNROLL 12
#endif
            do k=1,nrl
              fv2(k)  =z(k,i+1)
            end do
#if defined __T3E
!DIR$ UNROLL 12
#endif
            do k=1,nrl
              fv1(k)  =z(k,i)
            end do
#if defined __T3E
!DIR$ UNROLL 12
#endif
            do k=1,nrl
              z(k,i+1)  =sv(i)*fv1(k) + cv(i)*fv2(k)
              z(k,i)    =cv(i)*fv1(k) - sv(i)*fv2(k)
            end do
          end do

          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0d0
          goto 1

        endif
15    continue

      return
      END SUBROUTINE ptqliv

!==----------------------------------------------==!


      SUBROUTINE peigsrtv(d,v,ldv,n,nrl)

!
!     This routine sort eigenvalues and eigenvectors 
!     generated by PTREDV and PTQLIV.  
!
!     AUTHOR : Carlo Cavazzoni - SISSA 1997
!              comments and suggestions to : cava@sissa.it
!

      USE kinds
      IMPLICIT NONE
      INTEGER n,ldv,nrl
      REAL(DP) d(n),v(ldv,n)

      INTEGER i,j,k
      REAL(DP) p
      save i,j,k
      save p
!      REAL(DP) fv(nrl)

      do 13 i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
          if(d(j).le.p)then
            k=j
            p=d(j)
          endif
        end do
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
!
!         Exchange local elements of eigenvectors.
!
#if defined __T3E
!DIR$ UNROLL 12
#endif
          do j=1,nrl
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
          END DO

        endif
13    continue
      return
      END SUBROUTINE peigsrtv

      FUNCTION pythag(a,b)
      USE kinds
      IMPLICIT NONE
      REAL(DP) a,b,pythag
      REAL(DP) absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          pythag=0.0d0
        else
          pythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END FUNCTION pythag
   !
   !-------------------------------------------------------------------------
   SUBROUTINE diagonalize( iopt, a, d, ev, n, nproc, mpime, comm_in )
     !-------------------------------------------------------------------------
     !
     ! This subroutine calculate eigenvalues and optionally eigenvectors
     ! of a real symmetric matrix "a" distributed or replicated across
     ! processors
     !
     ! iopt     input, choose the distributions of "a"
     !          iopt = 0 matrix "a" and "ev" are distributed across procs
     !          iopt = 1 matrix "a" and "ev" are replicated across procs
     ! a(:,:)   input, matrix to be diagonalized, overwritten on output 
     ! d(:)     output, array of the eigenvalues
     ! ev(:,:)  output, matrix of the eigenvectors
     ! n        input, dimension of matrix a
     ! nproc    input, number of processors
     ! mpime    input, index of this processors (starting from 0!)
     !
     ! When iopt = 0 the matrix "a" shoud have the following layout:
     !
     !     A(NRL,N) Local part of the global matrix A(N,N) to be reduced,
     !              The rows of the matrix are cyclically distributed among 
     !              processors with blocking factor 1.
     !              Example for NPROC = 4 :
     !              ROW | PE
     !              1   | 0
     !              2   | 1
     !              3   | 2
     !              4   | 3 
     !              5   | 0
     !              6   | 1
     !              ..  | ..
     !
     IMPLICIT NONE
     !   
     INTEGER,           INTENT(IN)    :: iopt, n, nproc, mpime
     REAL(DP),          INTENT(OUT)   :: d(n)
     REAL(DP),          INTENT(INOUT) :: a(n,n)
     REAL(DP),          INTENT(OUT)   :: ev(n,n)
     INTEGER, OPTIONAL, INTENT(IN)    :: comm_in
     !
     REAL(DP), ALLOCATABLE :: aloc(:,:)
     REAL(DP), ALLOCATABLE :: evloc(:,:)
     REAL(DP), ALLOCATABLE :: e(:)
     INTEGER               :: nrl, i, j, jl, ierr, comm
     !
     !
     IF ( PRESENT( comm_in ) ) THEN
        !
        comm = comm_in
        !
     ELSE
        !
        comm = MPI_COMM_WORLD
        !
     END IF
     !
     nrl = n / nproc
     !
     IF ( mpime < MOD( n, nproc ) ) nrl = nrl + 1
     !
     ALLOCATE( evloc( nrl, n ) )
     ALLOCATE( e( n ) )
     !
     IF ( iopt == 1 ) THEN
        !
        ALLOCATE( aloc( nrl, n ) )
        !
        DO i = 1, n
           DO jl = 1, nrl
              aloc( jl, i ) = a( (jl-1)*nproc + mpime + 1, i )
           END DO
        END DO
        !
        CALL ptredv( aloc, nrl, d, e, evloc, nrl, nrl, n, nproc, mpime, comm )
        !
        DEALLOCATE( aloc )
        !
      ELSE
         !
         CALL ptredv( a, nrl, d, e, evloc, nrl, nrl, n, nproc, mpime, comm )
         !
      END IF
      !
      CALL ptqliv( d, e, n, evloc, nrl, nrl )
      CALL peigsrtv( d, evloc, nrl, n, nrl )
      !
      IF ( iopt == 1 ) THEN
         !
         DO i = 1,n
            DO j = 1,n
               a(j,i) = 0.D0
            END DO
            DO jl = 1,nrl
                a((jl-1)*nproc + mpime + 1,i) = evloc(jl,i)
            END DO
         END DO
         !
#if defined (__PARA)
#  if defined (__MPI)
         CALL MPI_ALLREDUCE( a, ev, n*n, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
#  endif
#else
         ev = a
#endif
         !
      ELSE
         !
         DO i = 1, n
            DO jl = 1, nrl
              ev(jl,i) = evloc(jl,i)
            END DO
         END DO
         !
      END IF
      !
      DEALLOCATE( evloc )
      DEALLOCATE( e )
      !
      RETURN
      !
   END SUBROUTINE diagonalize

!==----------------------------------------------==!

   SUBROUTINE pdspev_drv( jobz, ap, lda, w, z, ldz, &
                          nrl, n, nproc, mpime, comm_in )
     USE kinds
     IMPLICIT NONE
     CHARACTER :: JOBZ
     INTEGER   :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, OPTIONAL, INTENT(IN) :: comm_in
     REAL(DP) :: ap( lda, * ), w( * ), z( ldz, * )
     REAL(DP) :: sd( n )
     INTEGER   :: comm
     !
     IF ( PRESENT( comm_in ) ) THEN
        comm = comm_in
     ELSE
        comm = MPI_COMM_WORLD
     END IF
     CALL ptredv(ap, lda, w, sd, z, ldz, nrl, n, nproc, mpime, comm)
     CALL ptqliv(w, sd, n, z, ldz, nrl)
     CALL peigsrtv(w, z, ldz, n, nrl)
     RETURN
   END SUBROUTINE pdspev_drv
 

!==----------------------------------------------==!

      SUBROUTINE dspev_drv( JOBZ, UPLO, N, AP, W, Z, LDZ )
        USE kinds
        IMPLICIT NONE
        CHARACTER ::       JOBZ, UPLO
        INTEGER   ::       IOPT, INFO, LDZ, N
        REAL(DP) ::  AP( * ), W( * ), Z( LDZ, * )
        REAL(DP), ALLOCATABLE :: WORK(:)

        ALLOCATE( work( 3*n ) )

#if defined __AIX

          IOPT = 0
          IF((JOBZ .EQ. 'V') .OR. (JOBZ .EQ. 'v') ) iopt = iopt + 1
          IF((UPLO .EQ. 'U') .OR. (UPLO .EQ. 'u') ) iopt = iopt + 20
          CALL DSPEV(IOPT, ap, w, z, ldz, n, work, 3*n)

#else 

          CALL DSPEV(jobz, uplo, n, ap(1), w(1), z(1,1), ldz, work, INFO)
          IF( info .NE. 0 ) THEN
            CALL errore( ' dspev_drv ', ' diagonalization failed ',info )
          END IF

#endif

        DEALLOCATE( work )
 
        RETURN
      END SUBROUTINE dspev_drv

   !-------------------------------------------------------------------------
   SUBROUTINE cdiagonalize( iopt, a, d, ev, n, nproc, mpime, comm_in )
     !-------------------------------------------------------------------------
     !
     ! ... this routine calls the appropriate Lapack routine for diagonalizing
     ! ... a complex Hermitian matrix
     !
     IMPLICIT NONE
     !
     INTEGER,           INTENT(IN)    :: iopt, n, nproc, mpime
     REAL(DP),          INTENT(OUT)   :: d(n)
     COMPLEX(DP),       INTENT(INOUT) :: a(n,n)
     COMPLEX(DP),       INTENT(OUT)   :: ev(n,n)
     INTEGER, OPTIONAL, INTENT(IN)    :: comm_in
     !
     INTEGER :: i, j, k, nrl, comm, ierr
     INTEGER :: info = 0
     !
     COMPLEX(DP), ALLOCATABLE :: aloc(:)
     COMPLEX(DP), ALLOCATABLE :: ap(:,:)
     COMPLEX(DP), ALLOCATABLE :: vp(:,:)
     REAL(DP),    ALLOCATABLE :: rwork(:)
     COMPLEX(DP), ALLOCATABLE :: cwork(:)
     !
     !
     IF ( PRESENT( comm_in ) ) THEN
        !
        comm = comm_in
        !
     ELSE
        !
        comm = MPI_COMM_WORLD
        !
     END IF
     !
     IF ( ( nproc == 1 ) .OR. ( n < nproc ) ) THEN
        !
        ! ... here use the standard scalar lapack drivers
        !
        ALLOCATE( aloc( n*(n+1)/2 ) )
        ALLOCATE( rwork( 4*n ) )
        ALLOCATE( cwork( 4*n ) )
        !
        ! ... copy the lower-diagonal part of the matrix according to the
        ! ... Lapack packed storage scheme for Hermitian matrices
        !
        k = 0
        DO j = 1, n
           DO i = j, n
              k = k + 1
              aloc(k) = a(i,j)
           END DO
        END DO
        !
        ! ... call the Lapack routine
        !
#if defined (__AIX)
        !
        iopt = 1
        !
        CALL ZHPEV( iopt, aloc, d, ev, n, n, cwork, 4*n )
        !
#else
        CALL ZHPEV( 'V', 'L', n, aloc, d, ev, n, cwork, rwork, info )
#endif
        !
        DEALLOCATE( cwork )
        DEALLOCATE( rwork )
        DEALLOCATE( aloc )
        !                                                                       
        IF ( info /=  0 ) &
           CALL errore( 'cdiagonalize', 'info', ABS( info ) )
        !
     ELSE
        !
        nrl = n / nproc
        !
        IF ( mpime < MOD( n, nproc ) )  nrl = nrl + 1
        !
        ALLOCATE( vp( nrl, n ) )
        ALLOCATE( rwork( n ) )
        ALLOCATE( cwork( n ) )
        !
        IF ( iopt == 1 ) THEN
           !
           ALLOCATE( ap( nrl, n ) )
           !
           DO j = 1, n
              DO i = 1, nrl
                 ap(i,j) = a(mpime+(i-1)*nproc+1,j)
              END DO
           END DO
           !
           CALL pzhptrd( n, nrl, ap, nrl, d, rwork, cwork, nproc, mpime, comm )
           CALL pzupgtr( n, nrl, ap, nrl, cwork, vp, nrl, nproc, mpime, comm)
           !
           DEALLOCATE( ap )
           !
        ELSE
           !
           CALL pzhptrd( n, nrl, a, nrl, d, rwork, cwork, nproc, mpime, comm )
           CALL pzupgtr( n, nrl, a, nrl, cwork, vp, nrl, nproc, mpime, comm )
           !
        END IF
        !
        CALL pzsteqr( 'V', n, nrl, d, rwork, vp, nrl, nproc, mpime, comm )
        !
        IF ( iopt == 1 ) THEN
           !
           DO j = 1,n
              DO i = 1,n
                 a(i,j) = 0.D0
              END DO
              DO i = 1, nrl
                 a((i-1)*nproc+mpime+1,j) = vp(i,j)
              END DO
           END DO
           !
#if defined (__PARA)
#  if defined (__MPI)
           CALL MPI_ALLREDUCE( a, ev, 2*n*n, &
                               MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
#  endif
#else
           ev(:,:) = a(:,:)
#endif
        ELSE
           !
           DO j = 1,n
              DO i = 1, nrl
                 ev(i,j) = vp(i,j)
              END DO
           END DO
           !
        END IF
        !
        DEALLOCATE( vp )
        DEALLOCATE( rwork )
        DEALLOCATE( cwork )
        !
     END IF
     !
     RETURN
     !
   END SUBROUTINE cdiagonalize
   !
   !-------------------------------------------------------------------------
   SUBROUTINE pzhptrd( n, nrl, ap, lda, d, e, tau, nproc, me, comm )
     !-------------------------------------------------------------------------
      !
      !  Parallel MPI version of the LAPACK routine ZHPTRD
      !
      !     Carlo Cavazzoni (carlo.cavazzoni@cineca.it) -- CINECA
      !     Dicember 12, 1999
      !
      !  REFERENCES :
      !
      !     NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
      !     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
      !     CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
      !
      !     PARALLEL NUMERICAL ALGORITHMS,
      !     T.L. FREEMAN AND C.PHILLIPS,
      !     PRENTICE HALL INTERNATIONAL (1992).
      !
      !     LAPACK routine (version 2.0) --
      !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
      !     Courant Institute, Argonne National Lab, and Rice University
      !

      IMPLICIT NONE

!     .. __SCALAR Arguments ..
      INTEGER            LDA, N, NRL, NPROC, ME, comm
!     ..
!     .. Array Arguments ..
      REAL(DP)             D( * ), E( * )
      COMPLEX(DP)         AP(LDA, * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  PZHPTRD reduces a complex Hermitian distributed matrix AP  to
!  real symmetric tridiagonal form T by a unitary similarity
!  transformation: Q**H * A * Q = T.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the mglobal atrix AP.  N >= 0.
!
!  NRL     (input) INTEGER
!          The number of local rows of the matrix AP. NRL >= 0.
!
!  AP      (input/output) COMPLEX(DP) array, dimension (LDA,N)
!          On entry, the Hermitian matrix AP.
!          The rows of the matrix are distributed among processors
!          with blocking factor 1.
!              Example for NPROC = 4 :
!              ROW | PE
!              1   | 0
!              2   | 1
!              3   | 2
!              4   | 3
!              5   | 0
!              6   | 1
!              ..  | ..

!          On exit, the diagonal and first subdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements below the first
!          subdiagonal, with the array TAU, represent the unitary
!          matrix Q as a product of elementary reflectors; 
!
!  LDA     (input) INTEGER
!          Leading dimension of the local matrix AP, LDA > NRL
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = AP(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i+1,i) 
!
!  TAU     (output) COMPLEX(DP) array, dimension (N-1)
!          The __SCALAR factors of the elementary reflectors (see Further
!          Details).
!
!  NPROC   (input) INTEGER
!          Number of processors
!
!  ME      (input) INTEGER
!          Index of the local processor  ( 0, 1, 2, ..., NPROC-1 )

!
!  Further Details
!  ===============
!
!  the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex __SCALAR, and v is a complex vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..

      COMPLEX(DP)  ONE, ZERO, HALF
      PARAMETER   ( ONE = ( 1.0D+0, 0.0D+0 ),ZERO = ( 0.0D+0, 0.0D+0 ),  &
     &             HALF = ( 0.5D+0, 0.0D+0 ) )
      REAL(DP)      RONE, RZERO
      PARAMETER   ( RONE = 1.0D+0, RZERO = 0.0D+0 ) 


      INTEGER QI
      INTEGER IL(N+1)
      INTEGER OW(N+1)  
      COMPLEX(DP) CTMP
      COMPLEX(DP) CTMPV(N+1)
      COMPLEX(DP) TAUL(N+1)
      COMPLEX(DP) APKI(N+1)
      REAL(DP)     TMP
      REAL(DP)     TMPV(N+1)

!     ..
!     .. Local __SCALARs ..
      INTEGER            J, I, I1, K, I2, NI1, JL
      INTEGER            KL, J1
      COMPLEX(DP)         ALPHA, TAUI
      INTEGER            KNT, IERR
      REAL(DP)             ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZAXPY
      EXTERNAL           ZDSCAL, ZSCAL                                          
!     ..
!     .. External Functions ..
      COMPLEX(DP)         ZDOTC
      EXTERNAL           ZDOTC
      REAL(DP)             DLAMCH, DLAPY3, DZNRM2
      COMPLEX(DP)         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DABS, DBLE, AIMAG, SIGN
!     cmplx removed because preprocessed
!
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF(N.LE.0) THEN
        RETURN
      END IF


      DO I = 1,N+1
        QI     = (I-1)/NPROC
        OW(I)  = MOD((I-1),NPROC) 
        IF(ME .le. OW(I) ) then
          IL(I) = QI + 1
        ELSE
          IL(I) = QI
        END IF
      END DO                                                                    
!
!        Reduce the lower triangle of A. 
!
         IF (OW(1).EQ.ME) THEN
           AP( IL(1), 1 ) = DBLE( AP( IL(1), 1 ) )
         END IF                                                             

         DO I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            IF (OW(I+1).EQ.ME) THEN
              ALPHA = AP( IL(I+1), I ) 
            END IF                                                             

#if defined __PARA
#  if defined __MPI
            CALL MPI_BCAST(ALPHA,1,MPI_DOUBLE_COMPLEX,OW(I+1),comm,IERR)
#  endif
#endif

            IF( (N-I).LE.0 ) THEN
              TAUI = RZERO
            ELSE
              IF(OW(I+2).EQ.ME) THEN
                I2 = IL(I+2)
              ELSE
                I2 = IL(I+2) + 1          ! I+2
              ENDIF
              NI1 = NRL - I2 + 1          ! N-I-1

              IF((N-I-1).GT.0) THEN
                XNORM = DZNRM2( NI1, AP( I2, I ), 1 )
#if defined __PARA
                XNORM = XNORM ** 2 
#  if defined __MPI
                CALL MPI_ALLREDUCE(XNORM,TMP,1,MPI_DOUBLE_PRECISION, &
     &                             MPI_SUM, comm,IERR)
#  endif
                XNORM = SQRT(TMP)
#endif
              ELSE
                XNORM = 0.0D0
              ENDIF

              ALPHR = DBLE( ALPHA )
              ALPHI = AIMAG( ALPHA )
              IF( XNORM.EQ.RZERO .AND. ALPHI.EQ.RZERO ) THEN
                TAUI = RZERO
              ELSE  
                BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
                SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
                RSAFMN = RONE / SAFMIN   
                IF( DABS( BETA ).LT.SAFMIN ) THEN
                  KNT = 0
   10             CONTINUE
                  KNT = KNT + 1

                  IF(NI1.GT.0) THEN
                    CALL ZDSCAL( NI1, RSAFMN, AP( I2, I ), 1 )
                  ENDIF

                  BETA = BETA*RSAFMN
                  ALPHI = ALPHI*RSAFMN
                  ALPHR = ALPHR*RSAFMN
                  IF( DABS( BETA ).LT.SAFMIN ) GO TO 10 

                  IF((N-I-1).GT.0) THEN
                    XNORM = DZNRM2( NI1, AP( I2, I ), 1 )
#if defined __PARA
                    XNORM = XNORM ** 2 
#  if defined __MPI
                    CALL MPI_ALLREDUCE(XNORM,TMP,1,MPI_DOUBLE_PRECISION, &
     &                                 MPI_SUM, comm,IERR)
#  endif
                    XNORM = SQRT(TMP)
#endif
                  ELSE
                    XNORM = 0.0D0
                  ENDIF

                  ALPHA = CMPLX( ALPHR, ALPHI )
                  BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
                  TAUI = CMPLX( (BETA-ALPHR)/BETA, -ALPHI/BETA )
                  ALPHA = ZLADIV( ONE, ALPHA-BETA )

                  IF(NI1.GT.0) THEN
                    CALL ZSCAL( NI1, ALPHA, AP( I2, I ), 1 )
                  ENDIF

                  ALPHA = BETA
                  DO J = 1, KNT
                    ALPHA = ALPHA*SAFMIN
                  END DO   

                ELSE

                  TAUI = CMPLX( (BETA-ALPHR)/BETA, -ALPHI/BETA )
                  ALPHA = ZLADIV( ONE, ALPHA-BETA )

                  IF(NI1.GT.0) THEN
                    CALL ZSCAL( NI1, ALPHA, AP( I2, I ), 1 )
                  ENDIF

                  ALPHA = BETA
                END IF
              END IF    
            ENDIF
!
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               ! ... AP( I+1, I ) = ONE
               IF (OW(I+1).EQ.ME) THEN
                 AP( IL(I+1), I ) = ONE
               END IF
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!

               ! ... broadcast A(K,I)
               IF(OW(I+1).EQ.ME) THEN
                 I1 = IL(I+1)
               ELSE
                 I1 = IL(I+1) + 1          ! I+2
               ENDIF

#if defined __PARA
               DO J = I+1, N
                 CTMPV(J) = ZERO
               END DO
               DO JL = I1, NRL
                 J = ME + (JL-1)*NPROC + 1
                 CTMPV(J) = AP(JL,I)
               END DO
#  if defined __MPI
               CALL MPI_ALLREDUCE(CTMPV(I+1),APKI(I+1),N-I, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, comm,IERR)
#  endif
#else
               DO J = I+1,N
                 APKI(J) = AP(J,I)
               ENDDO
#endif
               DO J = I+1, N+1 
                 TAU(J-1) = ZERO
               END DO
               DO JL = I1, NRL 
                 J = ME + (JL-1)*NPROC + 1
                 TAU(J-1) = ZERO
                 DO K = I+1, J
                   TAU(J-1) = TAU(J-1) + TAUI * AP(JL,K) * APKI(K)
                 END DO
               END DO
               DO J = I+1, N 
                 IF(OW(J+1).EQ.ME) THEN
                   J1 = IL(J+1)
                 ELSE
                   J1 = IL(J+1) + 1          ! I+2
                 ENDIF
                 DO KL = J1, NRL
                   K = ME + (KL-1)*NPROC + 1
                   TAU(J-1) = TAU(J-1) + TAUI * CONJG(AP(KL,J)) * APKI(K)
                 END DO
               END DO


#if defined __PARA
               ! ... parallel sum TAU
#  if defined __MPI
               CALL MPI_ALLREDUCE(TAU(I),CTMPV(I),N-I+1, &
     &              MPI_DOUBLE_COMPLEX,MPI_SUM,comm,IERR)
#  endif
               DO J = I, N
                 TAU(J) = CTMPV(J)
               END DO
#endif
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ! ... ALPHA = -HALF*TAUI*ZDOTC(N-I,TAU(I),1,AP(I+1,I),1)

               JL = 1
               DO J = I, N
                 IF(OW(J+1).EQ.ME) THEN
                   TAUL(JL) = TAU(J)
                   JL = JL + 1
                 END IF
               END DO
               IF(OW(I+1).EQ.ME) THEN
                 I1 = IL(I+1)
               ELSE
                 I1 = IL(I+1) + 1          ! I+1
               ENDIF
               NI1 = NRL - I1 + 1          ! N-I
               ALPHA = -HALF*TAUI*ZDOTC(NI1,TAUL(1),1,AP(I1,I),1)

#if defined __PARA
#  if defined __MPI
               CALL MPI_ALLREDUCE(ALPHA,CTMP,1, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, comm,IERR)
#  endif
               ALPHA = CTMP
#endif


#if defined __PARA
               CALL ZAXPY(NI1,ALPHA,AP(I1,I),1,TAUL(1),1)
               JL = 1
               DO J = I, N
                 CTMPV(J) = ZERO
                 IF(OW(J+1).EQ.ME) THEN
                   CTMPV(J) = TAUL(JL)
                   JL = JL + 1
                 END IF
               END DO
#  if defined __MPI
               CALL MPI_ALLREDUCE(CTMPV(I),TAU(I),N-I+1, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, comm,IERR)
#  endif
#else
               CALL ZAXPY(N-I,ALPHA,AP(I+1,I),1,TAU(I),1)
#endif

!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               ! ... broadcast A(K,I)
               IF(OW(I+1).EQ.ME) THEN
                 I1 = IL(I+1)
               ELSE
                 I1 = IL(I+1) + 1          ! I+2
               ENDIF

#if defined __PARA
               DO J = I+1, N
                 CTMPV(J) = ZERO
               END DO
               DO JL = I1, NRL
                 J = ME + (JL-1)*NPROC + 1
                 CTMPV(J) = AP(JL,I)
               END DO
#  if defined __MPI
               CALL MPI_ALLREDUCE(CTMPV(I+1),APKI(I+1),N-I, &
     &         MPI_DOUBLE_COMPLEX, MPI_SUM, comm,IERR)
#  endif
#else
               DO J = I+1, N
                 APKI(J) = AP(J,I)
               END DO
#endif

               DO K = I+1,N
                 DO JL = I1,NRL
                   J = ME + (JL-1)*NPROC + 1
                   AP(JL,K) = AP(JL,K) - ONE * AP(JL,I) * CONJG(TAU(K-1)) - &
     &                CONJG(ONE) * TAU(J-1) * CONJG(APKI(K))
                 END DO
               END DO
!
            END IF
            IF(OW(I+1).EQ.ME) THEN
              AP(IL(I+1),I) = E( I )
            END IF
            IF(OW(I).EQ.ME) THEN
              D( I ) = DBLE(AP( IL(I),I ))
            END IF
#if defined __PARA
#  if defined __MPI
            CALL MPI_BCAST(D(I),1,MPI_DOUBLE_PRECISION,OW(I), comm,IERR)
#  endif
#endif
            TAU( I ) = TAUI
         END DO
         IF(OW(I).EQ.ME) THEN
            D( N ) = DBLE(AP( IL(I),I ))
         END IF
#if defined __PARA
#  if defined __MPI
         CALL MPI_BCAST(D(N),1,MPI_DOUBLE_PRECISION,OW(I), comm,IERR)
#  endif
#endif
!
      RETURN

!
!     End of ZHPTRD
!
      END SUBROUTINE pzhptrd

!==----------------------------------------------==!

   SUBROUTINE pzupgtr( n, nrl, ap, lda, tau, q, ldq, nproc, me, comm)

!
!  Parallel MPI version of the LAPACK routine ZUPGTR
!
!     Carlo Cavazzoni (carlo.cavazzoni@cineca.it) -- CINECA
!     Dicember 12, 1999
!
!  REFERENCES :
!
!     NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
!     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
!     CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
!
!     PARALLEL NUMERICAL ALGORITHMS,
!     T.L. FREEMAN AND C.PHILLIPS,
!     PRENTICE HALL INTERNATIONAL (1992).
!
!     LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University


      USE kinds
      IMPLICIT NONE

!
!     .. __SCALAR Arguments ..

      INTEGER            INFO, LDQ, N, LDA, NRL, NPROC, ME, comm
!     ..
!     .. Array Arguments ..
      COMPLEX(DP)         AP(LDA, * ), Q( LDQ, * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  PZUPGTR generates a complex unitary matrix Q which is defined as the
!  product of n-1 elementary reflectors H(i) of order n, as returned by
!  PZHPTRD :
!
!  Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the mglobal atrix AP.  N >= 0.
!
!  NRL     (input) INTEGER
!          The number of local rows of the matrix AP. NRL >= 0.
!
!  AP      (input) COMPLEX(DP) array, dimension (LDA,N)
!          The vectors which define the elementary reflectors, as
!          returned by PZHPTRD.
!          The rows of the matrix are distributed among processors
!          with blocking factor 1.
!              Example for NPROC = 4 :
!              ROW | PE
!              1   | 0
!              2   | 1
!              3   | 2
!              4   | 3
!              5   | 0
!              6   | 1
!              ..  | ..
!
!  LDA     (input) INTEGER
!          Leading dimension of the local matrix AP, LDA > NRL
!
!  TAU     (input) COMPLEX(DP) array, dimension (N-1)
!          TAU(i) must contain the __SCALAR factor of the elementary
!          reflector H(i), as returned by PZHPTRD.
!
!  Q       (output) COMPLEX(DP) array, dimension (LDQ,N)
!          The N-by-N unitary matrix Q.
!          The rows of the matrix are distributed among processors
!          in the same way of the matrix AP
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= max(1,NRL).
!
!  NPROC   (input) INTEGER
!          Number of processors
!
!  ME      (input) INTEGER
!          Index of the local processor  ( 0, 1, 2, ..., NPROC-1 )
!
!  =====================================================================
!
!     .. Parameters ..

      COMPLEX(DP)         ONE, ZERO
      PARAMETER          ( ONE = (1.0D+0,0.0D+0), ZERO = (0.0D+0,0.0D+0) ) 

      INTEGER QI
      INTEGER IL(N+1)
      INTEGER OW(N+1)
      COMPLEX(DP) CTMP
      COMPLEX(DP) CTMPV(N+1)
      COMPLEX(DP) WORK(N+1)

!     ..
!     .. Local __SCALARs ..
      INTEGER            I, IINFO, J, K, JL, KL, J1, I1, I2, NI1, L, IERR
!     ..

!     .. Executable Statements ..
!
!     Test the input arguments
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
        RETURN
      END IF

!
      DO I = 1,N+1
        QI     = (I-1)/NPROC
        OW(I)  = MOD((I-1),NPROC)
        IF(ME .le. OW(I) ) then
          IL(I) = QI + 1
        ELSE
          IL(I) = QI
        END IF
      END DO
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
      IF(OW(1).EQ.ME) THEN
        Q( IL(1), 1 ) = ONE
        DO KL = 2, NRL
          Q( KL, 1 ) = ZERO
        END DO
        DO J = 2, N
          Q( IL(1), J ) = ZERO
        END DO
      ELSE
        DO KL = 1, NRL
          Q( KL, 1 ) = ZERO
        END DO
      ENDIF

      DO J = 2, N
        IF(OW(J+1).EQ.ME) THEN
          J1 = IL(J+1)
        ELSE
          J1 = IL(J+1) + 1 
        ENDIF
        DO KL = J1, NRL
          Q( KL, J ) = AP( KL, J-1 )
        END DO
      END DO

      IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
        DO I = N-1, 1, -1
!
!         Apply H(i) to A(i:m,i:n) from the left
!
          IF( I.LT.(N-1) ) THEN

            IF(OW(I+1).EQ.ME) THEN
              Q( IL(I+1), I+1 ) = ONE
            END IF
!
!           Form  H * C
!
            IF( TAU(I).NE.ZERO ) THEN
!
!             w := C' * v
!
              IF(OW(I+1).EQ.ME) THEN
                I1 = IL(I+1)
              ELSE
                I1 = IL(I+1) + 1 
              ENDIF
              DO J = 1, N-1-I
                CTMP = ZERO
                DO KL = I1, NRL
                  CTMP = CTMP + CONJG( Q( KL, J+I+1 ) ) * Q( KL,I+1 )  
                END DO
                WORK(J) = CTMP
              END DO 

#if defined __PARA
#  if defined __MPI
              CALL MPI_ALLREDUCE(WORK,CTMPV,N-1-I,MPI_DOUBLE_COMPLEX, &
     &             MPI_SUM, comm,IERR)
#  endif
              DO J = 1,N-1-I
                WORK(J) = CTMPV(J)
              END DO
#endif
!
!             C := C - v * w'
!
              DO J = 1, N-1-I
                CTMP = -TAU(I) * CONJG( WORK( J ) ) 
                DO KL = I1, NRL
                  Q( KL, J+I+1 ) = Q( KL, J+I+1 ) + CTMP * Q(KL, I+1)    
                END DO
              END DO 
            END IF
          END IF

          IF( I.LT.(N-1) ) THEN
            IF(OW(I+2).EQ.ME) THEN
              I2 = IL(I+2)              ! I+2
            ELSE
              I2 = IL(I+2) + 1          ! local ind. of the first element > I+2 
            ENDIF
            NI1 = NRL - I2 + 1          ! N-I-1
            CALL ZSCAL( NI1, -TAU( I ), Q( I2, I+1 ), 1 )
          END IF

          IF(OW(I+1).EQ.ME) THEN
            Q( IL(I+1), I+1 ) = ONE - TAU( I )
          END IF
!
!             Set A(1:i-1,i) to zero
!
          DO L = 1, I - 1
            IF(OW(L+1).EQ.ME) THEN
              Q( IL(L+1), I+1 ) = ZERO
            END IF
          END DO
        END DO   
      END IF


      RETURN

!
!     End of ZUPGTR
!
      END SUBROUTINE pzupgtr

!==----------------------------------------------==!

      SUBROUTINE pzsteqr( compz, n, nrl, d, e, z, ldz, nproc, me, comm )
!
!  Parallel MPI version of the LAPACK routine ZHPTRD
!
!     Carlo Cavazzoni (carlo.cavazzoni@cineca.it) -- CINECA
!     Dicember 12, 1999
!
!  REFERENCES :
!
!     NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
!     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
!     CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
!
!     PARALLEL NUMERICAL ALGORITHMS,
!     T.L. FREEMAN AND C.PHILLIPS,
!     PRENTICE HALL INTERNATIONAL (1992).
!
!     LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!

      USE kinds
      IMPLICIT NONE

!     .. __SCALAR Arguments ..
      CHARACTER          COMPZ
      INTEGER            LDZ, N, NRL, NPROC, ME, comm
!     ..
!     .. Array Arguments ..
      REAL(DP)   D( * ), E( * )
      COMPLEX(DP)         Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  PZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band complex Hermitian matrix can also
!  be found if PZHPTRD has been used to reduce this
!  matrix to tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  Hermitian matrix.  On entry, Z must contain the
!                  unitary matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the mglobal atrix AP.  N >= 0.
!
!  NRL     (input) INTEGER
!          The number of local rows of the matrix AP. NRL >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) COMPLEX(DP) array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the unitary
!          matrix used in the reduction to tridiagonal form.
!          On exit if COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original Hermitian matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!          The rows of the matrix are distributed among processors
!          with blocking factor 1, i.e. for NPROC = 4 :
!              ROW Index | Processor index owning the row 
!                    1   |    0
!                    2   |    1
!                    3   |    2
!                    4   |    3
!                    5   |    0
!                    6   |    1
!                    ..  |    ..
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,NRL).
!
!  NPROC   (input) INTEGER
!          Number of processors
!
!  ME      (input) INTEGER
!          Index of the local processor  ( 0, 1, 2, ..., NPROC-1 )
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(DP)   ZERO, ONE, TWO, THREE, CTEMP, STEMP
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
     &                   THREE = 3.0D0 )
      COMPLEX(DP)         CZERO, CONE,ZTEMP
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..

      INTEGER QI, KL, INFO
      INTEGER IL(N+1)
      INTEGER OW(N+1)
      REAL(DP)  WORK(2*N)

!     .. Local __SCALARs ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
     &                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,   &
     &                   NM1, NMAXIT
      REAL(DP)   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,   &
     &                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL(DP)   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA, &
     &                   ZLASET, ZLASR, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( (LDZ.LT.1) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX(1,NRL) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF(N.LE.0) THEN
        RETURN
      END IF
!
      DO I = 1,N+1
        QI     = (I-1)/NPROC
        OW(I)  = MOD((I-1),NPROC)
        IF(ME .le. OW(I) ) then
          IL(I) = QI + 1
        ELSE
          IL(I) = QI
        END IF
      END DO

      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 .AND. OW(1).EQ.ME ) Z( IL(1), 1 ) = CONE
         RETURN
      END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF( ICOMPZ.EQ.2 ) THEN
        CALL ZLASET( 'Full', NRL, N, CZERO, CZERO, Z, LDZ )
        DO J = 1, N
          IF(OW(J).EQ.ME) THEN
            Z( IL(J), J ) = CONE
          END IF
        END DO
      END IF
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF( L1.GT.N )   GO TO 160
      IF( L1.GT.1 )   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = DABS( E( M ) )
            IF( TST.EQ.ZERO )        GO TO 30
            IF( TST.LE.( SQRT(DABS(D(M)))*SQRT(DABS(D(M+1))))*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )  GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO )
      END IF
!
!     Choose between QL and QR iteration
!
      IF( DABS( D( LEND ) ).LT.DABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = DABS( E( M ) )**2
               IF( TST.LE.( EPS2*DABS(D(M)) )*DABS(D(M+1))+ SAFMIN )GO TO 60
   50       CONTINUE
         END IF
!
         M = LEND
!
   60    CONTINUE
         IF( M.LT.LEND )  E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )     GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CTEMP = WORK( L )
               STEMP = WORK( N-1+L )
               IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                  DO KL = 1, NRL
                     ZTEMP = Z( KL, 1+L )
                     Z( KL, 1+L ) = CTEMP*ZTEMP - STEMP*Z( KL, L )
                     Z( KL, L )   = STEMP*ZTEMP + CTEMP*Z( KL, L )
                  END DO
               END IF
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )     GO TO 40
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT )   GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )  E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
!
   70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
           DO J = M - L + 1 - 1, 1, -1
             CTEMP =  WORK( L + J -1)
             STEMP =  WORK( N-1+L + J-1)
             IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
               DO KL = 1, NRL
                 ZTEMP = Z( KL, J+1+L-1 )
                 Z( KL, J+1+L-1 ) = CTEMP*ZTEMP - STEMP*Z( KL, J+L-1 )
                 Z( KL, J+L-1 ) = STEMP*ZTEMP + CTEMP*Z( KL, J+L-1 )
               END DO
             END IF
           END DO                                                         
         END IF
!
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
!
!        Eigenvalue found.
!
   80    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L.LE.LEND )   GO TO 40
         GO TO 140
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = DABS( E( M-1 ) )**2
               IF( TST.LE.(EPS2*DABS(D(M)))*DABS(D(M-1))+ SAFMIN )GO TO 110
  100       CONTINUE
         END IF
!
         M = LEND
!
  110    CONTINUE
         IF( M.GT.LEND )   E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )      GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CTEMP = WORK( M )
               STEMP = WORK( N-1+M )
               IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                  DO KL = 1, NRL
                     ZTEMP = Z( KL, L)
                     Z( KL, L )   = CTEMP*ZTEMP - STEMP*Z( KL, L-1 )
                     Z( KL, L-1 ) = STEMP*ZTEMP + CTEMP*Z( KL, L-1 )
                  END DO
               END IF
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )    GO TO 90
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT )  GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )     E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
!
  120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            DO J = 1, L - M
               CTEMP = WORK( M+J-1 )
               STEMP = WORK( N-1+M+J-1 )
               IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                  DO KL = 1, NRL
                     ZTEMP = Z( KL, J+M )
                     Z( KL, J+M )   = CTEMP*ZTEMP - STEMP*Z(KL, J+M-1)
                     Z( KL, J+M-1 ) = STEMP*ZTEMP + CTEMP*Z(KL, J+M-1)
                  END DO
               END IF
            END DO                                                         
         END IF
!
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L.GE.LEND )   GO TO 90
         GO TO 140
!
      END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
     &                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
     &                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
     &                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
     &                N, INFO )
      END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO )  INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( NRL, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
!
!     End of ZSTEQR
!
      END SUBROUTINE pzsteqr

!==----------------------------------------------==!

   SUBROUTINE zhpev_drv( JOBZ, UPLO, N, AP, W, Z, LDZ )

        USE kinds
        IMPLICIT NONE

        CHARACTER ::       JOBZ, UPLO
        INTEGER   ::       IOPT, INFO, LDZ, N
        COMPLEX(DP) ::  AP( * ), Z( LDZ, * )
        REAL(DP) ::  W( * )
        REAL(DP), ALLOCATABLE :: RWORK(:)
        COMPLEX(DP), ALLOCATABLE :: ZWORK(:)

#if defined __AIX

        IOPT = 0
        IF((JOBZ .EQ. 'V') .OR. (JOBZ .EQ. 'v') ) iopt = iopt + 1
        IF((UPLO .EQ. 'U') .OR. (UPLO .EQ. 'u') ) iopt = iopt + 20
        ALLOCATE( rwork( 4*n ) )
        CALL ZHPEV(IOPT, ap, w, z, ldz, n, rwork, 4*n)
        DEALLOCATE( rwork )

#else

        ALLOCATE( rwork( MAX(1, 3*n-2) ), zwork( MAX(1, 2*n-1)) )
        CALL ZHPEV(jobz, uplo, n, ap, w, z, ldz, zwork, rwork, INFO)
        DEALLOCATE( rwork, zwork )
        IF( info .NE. 0 ) THEN
          CALL errore( ' dspev_drv ', ' diagonalization failed ',info )
        END IF

#endif

        RETURN
   END SUBROUTINE zhpev_drv

!==----------------------------------------------==!

   SUBROUTINE pzhpev_drv( jobz, ap, lda, w, z, ldz, &
                          nrl, n, nproc, mpime, comm_in )
     USE kinds
     IMPLICIT NONE
     CHARACTER :: JOBZ
     INTEGER   :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, OPTIONAL, INTENT(IN) :: comm_in
     COMPLEX(DP) :: ap( lda, * ), z( ldz, * )
     REAL(DP) :: w( * )
     REAL(DP) :: rwork( n )
     COMPLEX(DP) :: cwork( n )
     INTEGER   :: comm
     !
     IF ( PRESENT( comm_in ) ) THEN
        comm = comm_in
     ELSE
        comm = MPI_COMM_WORLD
     END IF
     CALL pzhptrd( n, nrl, ap, lda, w, rwork, cwork, nproc, mpime, comm)
     IF( jobz == 'V' .OR. jobz == 'v' ) THEN
        CALL pzupgtr( n, nrl, ap, lda, cwork, z, ldz, nproc, mpime, comm)
     END IF
     CALL pzsteqr( jobz, n, nrl, w, rwork, z, ldz, nproc, mpime, comm)
     RETURN
   END SUBROUTINE pzhpev_drv



!==----------------------------------------------==!
!
!  My parallel blas
!
!==----------------------------------------------==!


SUBROUTINE mattr_drv( m, k, a, lda, b, ldb, nb, dims, coor, comm )
  !
  !  Compute B as the transpose of matrix A 
  !  A and B are distributed on a 2D cartesian processor
  !  grid in a block cyclic way (as in scalapack),
  !  using a block size of NB
  !
  !     B :=  A'
  !
  !  A is a K by M matrix
  !  B is an M by K matrix
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, k
  INTEGER, INTENT(IN) :: lda, ldb
  REAL*8              :: a(lda,*), b(ldb,*)
  INTEGER, INTENT(IN) :: nb, dims(2), coor(2), comm
  !
#if defined __MPI

  include 'mpif.h'

  !
  integer ierr
  integer ndims, rowid, colid
  integer coosrc(2), coodst(2), ipsrc, ipdst, mpime
  integer ihsnd, ihrcv
  logical periods(2)
  !
  integer :: iu
  integer :: i, j, nk, nm
  integer :: ii, jj
  integer :: isrc, jsrc
  integer :: idst, jdst
  integer :: itag
  integer :: nmb, nkb 
  integer :: istatus( MPI_STATUS_SIZE )
  real*8, allocatable :: abuf(:,:)
  !
  integer :: numroc
  integer :: indxg2l
  external :: numroc, indxg2l

  !
  CALL GRID2D_RANK( dims(1), dims(2), coor(1), coor(2), mpime )
  !
  iu = 200 + mpime
  !
  !  Compute the global number of blocks for matrix dimension
  !
  nmb = ( m + nb - 1 ) / nb  
  nkb = ( k + nb - 1 ) / nb
  !
  ! WRITE(iu,*) 'TR',coor(1),coor(2),' m,k,nmb,nkb,nb = ',m,k,nmb,nkb,nb
  ! WRITE(iu,*) 'TR',coor(1),coor(2),' ---------------- '
  !
  ALLOCATE( abuf( nb, nb ) )
  !
  DO i = 1, nmb
     DO j = 1, nkb
        !
        itag = j + nkb * (i-1)
        !
        coosrc(1) = MOD( (j-1), dims(1) )
        coosrc(2) = MOD( (i-1), dims(2) ) 
        !
        coodst(1) = MOD( (i-1), dims(1) )
        coodst(2) = MOD( (j-1), dims(2) ) 
        !
        CALL GRID2D_RANK( dims(1), dims(2), coosrc(1), coosrc(2), ipsrc )
        CALL GRID2D_RANK( dims(1), dims(2), coodst(1), coodst(2), ipdst )
        !
        jsrc = INDXG2L( 1 + (j-1)*nb, nb, coor(1), 0, dims(1) )
        isrc = INDXG2L( 1 + (i-1)*nb, nb, coor(2), 0, dims(2) )
        !
        jdst = INDXG2L( 1 + (j-1)*nb, nb, coor(2), 0, dims(2) )
        idst = INDXG2L( 1 + (i-1)*nb, nb, coor(1), 0, dims(1) )
        !
        nk = MIN( nb, k - (j-1)*nb ) !  number of element in the block
        nm = MIN( nb, m - (i-1)*nb ) !  number of element in the block
        !
        ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' cccsrc,ipsrc   = ', coosrc(1), coosrc(2), ipsrc
        ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' cccdst,ipdst   = ', coodst(1), coodst(2), ipdst
        ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' i,j,nm,nk      = ', i,j,nm,nk 
        ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' isrc,jsrc      = ', isrc,jsrc
        ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' idst,jdst      = ', idst,jdst
        !
        IF( ipsrc == ipdst ) THEN
          IF( ipsrc == mpime ) THEN
            DO ii = 1, nm
              DO jj = 1, nk
                b( idst + ii - 1, jdst + jj - 1 ) = a( jsrc + jj - 1, isrc + ii - 1 )
              END DO
            END DO
          END IF
        ELSE
          IF( ipsrc == mpime ) THEN
            DO ii = 1, nm
              DO jj = 1, nk
                abuf( ii, jj ) = a( jsrc + jj - 1, isrc + ii - 1 )
                !
                ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' ii,jj,abuf_snd = ', ii,jj,abuf( ii, jj )
                !
              END DO
            END DO
            CALL MPI_ISEND( abuf, nb*nb, MPI_DOUBLE_PRECISION, ipdst, itag, comm, ihsnd, ierr )
            CALL mpi_wait(ihsnd, istatus, ierr)
          ELSE IF( ipdst == mpime ) THEN
            CALL MPI_IRECV( abuf, nb*nb, MPI_DOUBLE_PRECISION, ipsrc, itag, comm, ihrcv, ierr )
            CALL mpi_wait(ihrcv, istatus, ierr)
            DO jj = 1, nk
              DO ii = 1, nm
                !
                ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' ii,jj,abuf_rcv = ', ii,jj,abuf( ii, jj )
                !
                b( idst + ii - 1, jdst + jj - 1 ) = abuf( ii, jj )
              END DO
            END DO
          END IF
        END IF 
        !
        ! WRITE(iu,*) 'TR',coor(1),coor(2),itag,' ---------------- '
        !
     END DO
  END DO

#else

  INTEGER :: i, j

  DO j = 1, k
     DO i = 1, m
        B( i, j ) = A( j, i )
     END DO
  END DO

#endif

  RETURN

END SUBROUTINE mattr_drv


! ---------------------------------------------------------------------------------

SUBROUTINE matsplit_drv( m, k, ar, ldar, a, lda, nb, dims, coor, comm )
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, k
  INTEGER, INTENT(IN) :: ldar
  REAL*8              :: ar(ldar,*)  !  matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL*8              :: a(lda,*)
  INTEGER, INTENT(IN) :: nb, coor(2), dims(2), comm
  !
  INTEGER :: i, j, nra, nca, ii, jj
  !
  INTEGER  :: numroc, INDXL2G
  EXTERNAL :: numroc, INDXL2G

  nra = NUMROC( m, nb, coor(1), 0, dims(1) )  !  total number of local row for matrix A, C
  nca = NUMROC( k, nb, coor(2), 0, dims(2) )  !  total number of local columns of A

  do j = 1, nca
     jj = INDXL2G( j, NB, coor(2), 0, dims(2) )
     do i = 1, nra
        ii = INDXL2G( i, NB, coor(1), 0, dims(1) )
        a( i, j ) = ar( ii, jj )
     end do
  end do

  RETURN

END SUBROUTINE matsplit_drv



! ---------------------------------------------------------------------------------

SUBROUTINE matmerge_drv( m, k, a, lda, ar, ldar, nb, dims, coor, comm )
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, k
  INTEGER, INTENT(IN) :: ldar
  REAL*8              :: ar(ldar,*)  !  matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL*8              :: a(lda,*)
  INTEGER, INTENT(IN) :: nb, coor(2), dims(2), comm
  !
  INTEGER :: i, j, ii, jj, ierr

#if defined __MPI

  include 'mpif.h'
  !

  INTEGER :: jsrc, isrc, ipsrc, coosrc(2)
  INTEGER :: nmb, nkb, nk, nm, mpime

  REAL*8, ALLOCATABLE :: buf(:,:)
  !
  INTEGER  :: INDXG2L
  EXTERNAL :: INDXG2L
  !
  CALL GRID2D_RANK( dims(1), dims(2), coor(1), coor(2), mpime )

  nmb = ( m + nb - 1 ) / nb  
  nkb = ( k + nb - 1 ) / nb

  ALLOCATE( buf( nb, nb ) )

  DO j = 1, nkb
     DO i = 1, nmb
        !
        coosrc(1) = MOD( (i-1), dims(1) )
        coosrc(2) = MOD( (j-1), dims(2) )
        !
        CALL GRID2D_RANK( dims(1), dims(2), coosrc(1), coosrc(2), ipsrc )
        !
        isrc = INDXG2L( 1 + (i-1)*nb, nb, coor(1), 0, dims(1) )
        jsrc = INDXG2L( 1 + (j-1)*nb, nb, coor(2), 0, dims(2) )
        !
        nm = MIN( nb, m - (i-1)*nb ) !  number of element in the block
        nk = MIN( nb, k - (j-1)*nb ) !  number of element in the block

        IF( ipsrc == mpime ) THEN
           DO jj = 1, nk
              DO ii = 1, nm
                 buf( ii, jj ) = a( isrc + ii - 1, jsrc + jj - 1 )
              END DO
           END DO
        ENDIF
        !
        CALL MPI_BCAST( buf, nb*nb, MPI_DOUBLE_PRECISION, ipsrc, comm, ierr )
        !
        do jj = 1, nk
           do ii = 1, nm
              ar( ii + (i-1)*nb, jj + (j-1)*nb ) = buf( ii, jj )
           end do
        end do
        !
     END DO
  END DO
  !
  DEALLOCATE( buf )

#else

  DO j = 1, k
     DO i = 1, m
        ar( i, j ) = a( i, j )
     END DO
  END DO

#endif

  RETURN
END SUBROUTINE matmerge_drv



! ---------------------------------------------------------------------------------


SUBROUTINE matscal_drv( m, n, beta, c, ldc, nb, dims, coor, comm )
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: m, n
  REAL*8,  INTENT(IN) :: beta
  INTEGER, INTENT(IN) :: ldc
  REAL*8              :: c(ldc,*)
  INTEGER, INTENT(IN) :: nb, coor(2), dims(2), comm
  !
  INTEGER :: i, j, nr, nc, ierr
  !
  INTEGER  :: numroc
  EXTERNAL :: numroc

  nr  = NUMROC( m, nb, coor(1), 0, dims(1) )  ! local row of C
  nc  = NUMROC( n, nb, coor(2), 0, dims(2) )  ! local colum of C 

  IF( beta == 0.0d0 ) THEN
    do j = 1, nc
      do i = 1, nr
        c(i,j) = 0.0d0
      end do
    end do 
  ELSE
    do j = 1, nc
      do i = 1, nr
        c(i,j) = beta * c(i,j)
      end do
    end do 
  END IF

  RETURN

END SUBROUTINE


! ---------------------------------------------------------------------------------


SUBROUTINE matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, nb, dims, coor, comm )
  !
  implicit none
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER, INTENT(IN) :: m, n, k
  REAL*8, INTENT(IN) :: alpha, beta
  INTEGER, INTENT(IN) :: lda, ldb, ldc
  REAL*8 :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: nb, dims(2), coor(2), comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI

  include 'mpif.h'
  !
  integer ierr
  integer ndims, rowid, colid
  integer comm_row, comm_col
  !
  integer :: ib, jb, kb, ibl, kbl, jbl
  integer :: i, j, kk, ni, nj, nk, nm, il, jl
  integer :: nnb, nmb, nkb 
  integer :: nr, nra, nca, nc, nrb, ncb, ii, jj
  integer :: nrt, ncat, nct, nrbt
  real*8, allocatable :: abuf(:,:), bbuf(:,:)
  real*8, allocatable :: at(:,:)
  real*8, allocatable :: bt(:,:)
  !
  integer :: numroc
  integer :: indxg2l
  external :: numroc, indxg2l
  !
  IF( dims(1) * dims(2) == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     call dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

     RETURN

  END IF
  !

  CALL MPI_COMM_SPLIT( COMM, coor(2), coor(1), COMM_COL, IERR )
  CALL MPI_COMM_RANK( COMM_COL, rowid, IERR )
  !
  CALL MPI_COMM_SPLIT( COMM, coor(1), coor(2), COMM_ROW, IERR )
  CALL MPI_COMM_RANK( COMM_ROW, colid, IERR )
  !
  !  Compute the global number of blocks for matrix dimension
  !
  nmb = ( m + nb - 1 ) / nb  
  !
  nnb = ( n + nb - 1 ) / nb
  !
  nkb = ( k + nb - 1 ) / nb
  !
  !  Compute the total number of local row for matrix A, C
  !
  nr  = NUMROC( m, nb, coor(1), 0, dims(1) )  ! local row of C
  !
  nra = NUMROC( m, nb, coor(1), 0, dims(1) )  ! local row of OP( A )
  nca = NUMROC( k, nb, coor(2), 0, dims(2) )  ! local columns of OP( A )
  !
  nrb = NUMROC( k, nb, coor(1), 0, dims(1) )  ! local row of OP( B )
  ncb = NUMROC( n, nb, coor(2), 0, dims(2) )  ! local colum of OP( B )
  !
  nc  = NUMROC( n, nb, coor(2), 0, dims(2) )  ! local colum of C 
  !
  IF( transa == 'T' .OR. transa == 't' ) THEN
    !
    ALLOCATE( at( nra, nca ) )
    !
    CALL mattr_drv( m, k, a, lda, at, nra, nb, dims, coor, comm )
    !
  END IF
  !
  IF( transb == 'T' .OR. transb == 't' ) THEN
    !
    ALLOCATE( bt( nrb, ncb ) )
    !
    CALL mattr_drv( k, n, b, ldb, bt, nrb, nb, dims, coor, comm )
    !
  END IF
  !
  !  Scale matrix C
  !
  CALL matscal_drv( m, n, beta, c, ldc, nb, dims, coor, comm )
  !
  !  loop over the rows/columns blocks of matrix OP(A)/OP(B)
  !
  do kb = 1, nkb
    !
    kk  = ( kb - 1 ) * nb + 1  !  first element of the block (global index)
    nk = MIN( nb, k - kk + 1 ) !  number of element in the block

    colid = MOD( (kb-1), dims(2) )  ! processor owning the block
    rowid = MOD( (kb-1), dims(1) )

    allocate( abuf( nr, nk ) )

    if( colid == coor(2) ) then
      nrt = 0
      ibl = 0
      kbl = INDXG2L( 1 + (kb-1)*nb, nb, coor(2), 0, dims(2) )
      do ib = 1 + coor(1), nmb, dims(1)
        i = ( ib - 1 ) * nb + 1
        ni = MIN( nb, m - i + 1 )
        IF( transa == 'T' .OR. transa == 't' ) THEN
          do jj = 1, nk
            do ii = 1, ni
              abuf( ii + nrt, jj ) = at( ii + ibl*nb, jj + kbl - 1 )
            end do
          end do
        ELSE
          do jj = 1, nk
            do ii = 1, ni
              abuf( ii + nrt, jj ) = a( ii + ibl*nb, jj + kbl - 1 )
            end do
          end do
        END IF
        nrt = nrt + ni
        ibl = ibl + 1
      end do
    end if
    CALL MPI_BCAST( abuf(1,1), nr*nk, MPI_DOUBLE_PRECISION, colid, COMM_ROW, IERR )

    allocate( bbuf( nk, nc ) )

    if( rowid == coor(1) ) then
      nct = 0 
      jbl = 0
      kbl = INDXG2L( 1 + (kb-1)*nb, nb, coor(1), 0, dims(1) )
      do jb = 1 + coor(2), nnb, dims(2)
        j = ( jb - 1 ) * nb + 1
        nj = MIN( nb, n - j + 1 )
        IF( transb == 'T' .OR. transb == 't' ) THEN
          do jj = 1, nj
            do ii = 1, nk
              bbuf( ii, jj + nct ) = bt( ii + kbl - 1, jj + jbl*nb )
            end do
          end do
        ELSE
          do jj = 1, nj
            do ii = 1, nk
              bbuf( ii, jj + nct ) = b( ii + kbl - 1, jj + jbl*nb )
            end do
          end do
        END IF
        nct = nct + nj
        jbl = jbl + 1
      end do
    end if

    CALL MPI_BCAST( bbuf(1,1), nk*nc, MPI_DOUBLE_PRECISION, rowid, COMM_COL, IERR )

    ii = 1
    do ib = 1 + coor(1), nmb, dims(1)
      i = ( ib - 1 ) * nb + 1
      il = INDXG2L( i, nb, coor(1), 0, dims(1) )
      ni = MIN( nb, m - i + 1 )
      jj = 1
      do jb = 1 + coor(2), nnb, dims(2)
        j = ( jb - 1 ) * nb + 1
        jl = INDXG2L( j, nb, coor(2), 0, dims(2) )
        nj = MIN( nb, n - j + 1 )
        call dgemm( 'n', 'n', ni, nj, nk, alpha, abuf( ii, 1 ), nra, bbuf( 1, jj ), nk, 1.0d0, c( il, jl ), ldc )
        jj = jj + nj
      end do
      ii = ii + ni
    end do

    deallocate( abuf )
    deallocate( bbuf )

  end do

  IF( ALLOCATED( at ) ) DEALLOCATE( at )
  IF( ALLOCATED( bt ) ) DEALLOCATE( bt )


#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     call dgemm( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

#endif


  RETURN

END SUBROUTINE


!==----------------------------------------------==!
!
! Copyright (C) 2005 Carlo Cavazzoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


SUBROUTINE rep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !  Parallel matrix multiplication with replicated matrix
  !
  implicit none
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER, INTENT(IN) :: m, n, k
  REAL*8, INTENT(IN) :: alpha, beta
  INTEGER, INTENT(IN) :: lda, ldb, ldc
  REAL*8 :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI

  include 'mpif.h'
  !

  INTEGER :: ME, I, II, J, JJ, IP, SOUR, DEST, INFO, IERR, ioff, ldx
  INTEGER :: NB, IB_S, NB_SOUR, IB_SOUR, IBUF
  INTEGER :: nproc, mpime, q, r

  REAL*8, ALLOCATABLE :: auxa( : )
  REAL*8, ALLOCATABLE :: auxc( : )

  !
  ! ... BODY
  !

  CALL MPI_COMM_SIZE(comm, NPROC, IERR)
  CALL MPI_COMM_RANK(comm, MPIME, IERR)

  IF ( NPROC == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     CALL DGEMM(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

     RETURN

  END IF

  ME = MPIME + 1
  Q = INT( m / NPROC )
  R = MOD( m , NPROC )

  ! ... Find out the number of elements in the local block
  !     along "M" first dimension os matrix A

  NB = Q
  IF( ME <= R ) NB = NB + 1

  ! ... Find out the global index of the local first row

  IF( ME <= R ) THEN
     ib_s = (Q+1)*(ME-1) + 1
  ELSE
     ib_s = Q*(ME-1) + R + 1
  END IF

  ldx = m / nproc + 1

  ALLOCATE( auxa( MAX( n, m ) * ldx ) )
  ALLOCATE( auxc( MAX( n, m ) * ldx ) )

  IF( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  ELSE
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( J, I + ioff )
        END DO
        ibuf = ibuf + ldx
     END DO
     !ioff = ib_s - 1
     !call mytranspose( A( 1, ioff + 1 ), lda, auxa(1), ldx, m, nb)
  END IF

  IF( beta /= 0.0d0 ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, n
        DO I = 1, NB
           auxc( ibuf + I ) = C( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  END IF

  CALL DGEMM( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

  ! ... Here processors exchange blocks

  DO IP = 0, NPROC-1

     ! ...    Find out the number of elements in the block of processor SOUR

     NB_SOUR = q
     IF( (IP+1) .LE. r ) NB_SOUR = NB_SOUR+1

     ! ...    Find out the global index of the first row owned by SOUR

     IF( (IP+1) .LE. r ) THEN
        ib_sour = (Q+1)*IP + 1
     ELSE
        ib_sour = Q*IP + R + 1
     END IF

     IF( mpime == ip ) auxa(1:n*ldx) = auxc(1:n*ldx)

     CALL MPI_BCAST( auxa(1), ldx*n, mpi_double_precision, ip, comm, IERR)

     IBUF = 0
     ioff = IB_SOUR - 1
     DO J = 1, N
        DO I = 1, NB_SOUR
           C( I + ioff, J ) = AUXA( IBUF + I )
        END DO
        IBUF = IBUF + ldx
     END DO

  END DO

  DEALLOCATE( auxa, auxc )


#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     CALL DGEMM(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

#endif



  RETURN

END SUBROUTINE rep_matmul_drv

!==----------------------------------------------==!
    END MODULE parallel_toolkit
!==----------------------------------------------==!
