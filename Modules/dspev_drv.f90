!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!

MODULE dspev_module


    IMPLICIT NONE
    SAVE
    PRIVATE

    PUBLIC :: pdspev_drv, dspev_drv, diagonalize

CONTAINS


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
      USE kinds,     ONLY : DP
      USE io_global, ONLY : stdout
      USE parallel_include

      IMPLICIT NONE

      INTEGER :: N, NRL, LDA, LDV
      INTEGER :: NPROC, ME, comm
      REAL(DP) :: A(LDA,N), D(N), E(N), V(LDV,N)
!
      REAL(DP) :: DDOT
!
      REAL(DP) :: g, scalef, sigma, kappa, f, h, tmp
      REAL(DP), ALLOCATABLE :: u(:)
      REAL(DP), ALLOCATABLE :: p(:)
      REAL(DP), ALLOCATABLE :: vtmp(:)

      REAL(DP) :: tu, tp, one_over_h
      REAL(DP) :: one_over_scale
      REAL(DP) :: redin(3), redout(3)
      REAL(DP), ALLOCATABLE :: ul(:)
      REAL(DP), ALLOCATABLE :: pl(:)
      integer :: l, i, j, k, t, tl, ierr
      integer :: kl, jl, ks, lloc
      integer, ALLOCATABLE :: is(:)
      integer, ALLOCATABLE :: ri(:)
     
      
      !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........

      IF( N == 0 ) THEN
        RETURN
      END IF

      ALLOCATE( u( n ), p( n+1 ), vtmp( n+1 ), ul( n ), pl( n ), is( n ), ri( n ) )

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
         H     = 0.0_DP

         IF ( L > 1 ) THEN

           SCALEF = 0.0_DP
           DO K = 1, is(l)
              SCALEF = SCALEF + DABS( A(K,I) )
           END DO

#if defined __PARA
#  if defined __MPI
           redin(1) = scalef
           CALL MPI_ALLREDUCE(redin, redout, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
           SCALEF = redout(1)
#  endif
#endif

           IF ( SCALEF .EQ. 0.0_DP )  THEN
             IF (RI(L).EQ.ME) THEN
               E(I) = A(is(L),I) 
             END IF
           ELSE 

             !  ......  CALCOLO DI SIGMA E DI H

             ONE_OVER_SCALE = 1.0_DP/SCALEF
             SIGMA = 0.0_DP
             DO k = 1,is(L)
               A(k,I) = A(k,I) * ONE_OVER_SCALE
               SIGMA  = SIGMA + A(k,I)**2 
             END DO

             IF( ri(l) .eq. me ) THEN
               F = A( is(l), i )
             ELSE
               F = 0.0_DP
             END IF

#if defined __PARA
#  if defined __MPI
             redin(1) = sigma
             redin(2) = f
             CALL MPI_ALLREDUCE( redin, redout, 2, MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR)
             SIGMA = redout(1)
             f     = redout(2)
#  endif
#endif

             G          = -SIGN(SQRT(SIGMA),F)
             H          = SIGMA - F*G
             ONE_OVER_H = 1.0_DP/H
             E(I)       = SCALEF*G

!  ......  COSTRUZIONE DEL VETTORE U

             DO k = 1,L          
               vtmp(k) = 0.0_DP 
             END DO

             k = ME + 1

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
#  if defined __MPI
             CALL MPI_ALLREDUCE(VTMP,U,L,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
             u(1:l) = vtmp(1:l)
#endif

!  ......  COSTRUZIONE DEL VETTORE P

             KAPPA = 0.0_DP

             DO J = 1,L

               vtmp(j) = 0.0_DP

               DO KL = 1, IS(J)
                 vtmp(J) = vtmp(J) + A(KL,J) * UL(KL)
               END DO

               IF(L.GT.J .AND. ME.EQ.RI(J)) then

                 DO K = J+1,L
                   vtmp(J) = vtmp(J) + A(IS(J),K) * U(K)
                 END DO
               END IF

               vtmp(J) = vtmp(J) * ONE_OVER_H
               KAPPA = KAPPA + vtmp(J) * U(J) * 0.5_DP * ONE_OVER_H
             END DO 

#if defined __PARA
#  if defined __MPI
             vtmp( l + 1 ) = kappa
             CALL MPI_ALLREDUCE( vtmp, p, l+1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, IERR )
             kappa = p( l + 1 )
#  endif
#else
             p(1:l) = vtmp(1:l)
#endif

             CALL DAXPY( l, -kappa, u, 1, p, 1 )
             CALL DGER( is(l), l, -1.0_DP, ul, 1, p, 1, a, lda )
             CALL DGER( is(l), l, -1.0_DP, p( me + 1 ), nproc, u, 1, a, lda )

             ! k = me + 1
             ! DO kl = 1,is(l)          
             !   PL(kl) = P(k)
             !   k      = k + NPROC
             ! END DO


             ! DO J = 1,L
             !   tu= U(J)
             !   tp= P(J)
             !   DO KL = 1,is(l)
             !     A(KL,j) = A(KL,j) - UL(KL) * tp 
             !     A(KL,j) = A(KL,j) - PL(KL) * tu 
             !   END DO
             ! END DO

           END IF

         ELSE 

           IF(RI(L).EQ.ME) THEN
             G = A(is(l),I)
           END IF

#if defined __PARA
#  if defined __MPI
           CALL MPI_BCAST(G,1,MPI_DOUBLE_PRECISION,RI(L),comm,IERR)
#  endif
#endif
           E(I) = G

         END IF

         D(I) = H

      END DO
      
      E(1) = 0.0_DP
      D(1) = 0.0_DP

      DO J = 1,N
        DO I = 1,NRL
          V(I,J) = 0.0_DP
        END DO
        IF(RI(J).EQ.ME) THEN
          V(IS(J),J) = 1.0_DP
        END IF
      END DO

      DO I = 2,N
        L = I - 1
        LLOC = IS(L)
        IF(D(I).NE.0.0_DP) THEN
          ONE_OVER_H = 1.0_DP/D(I)
          DO J = 1, L
            P(J) = DDOT(LLOC,V(1,J),1,A(1,I),1)
          END DO

#if defined __PARA
#  if defined __MPI
          CALL MPI_ALLREDUCE(P,VTMP,L,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
          vtmp(1:l) = p(1:l)
#endif

          DO J = 1, L
            call DAXPY(LLOC,-VTMP(J)*ONE_OVER_H,A(1,I),1,V(1,J),1)
          END DO

        END IF

      END DO 


      DO I = 1,N
        U(I) = 0.0_DP
        IF(RI(I).eq.ME) then
          U(I) = A(IS(I),I) 
        END IF
      END DO

#if defined __PARA
#  if defined __MPI
      CALL MPI_ALLREDUCE(U,D,N,MPI_DOUBLE_PRECISION,MPI_SUM,comm,IERR)
#  endif
#else
      D(1:N) = U(1:N)
#endif

      DEALLOCATE( u, p, vtmp, ul, pl, is, ri )

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
      USE kinds,     ONLY : DP
      USE io_global, ONLY : stdout

      IMPLICIT NONE

      INTEGER  :: n, nrl, ldz
      REAL(DP) :: d(n), e(n)
      REAL(DP) :: z(ldz,n)

      INTEGER  :: i, iter, mk, k, l, m
      REAL(DP) :: b, dd, f, g, p, r, c, s
      REAL(DP), ALLOCATABLE :: cv(:)
      REAL(DP), ALLOCATABLE :: sv(:)
      REAL(DP), ALLOCATABLE :: fv1(:)
      REAL(DP), ALLOCATABLE :: fv2(:)

      ALLOCATE( cv( n ) )
      ALLOCATE( sv( n ) )
      ALLOCATE( fv1( nrl ) )
      ALLOCATE( fv2( nrl ) )

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
          g=(d(l+1)-d(l))/(2.0_DP*e(l))
          r=pythag(g,1.0_DP)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0_DP
          c=1.0_DP
          p=0.0_DP
          do i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.0_DP)then
              d(i+1)=d(i+1)-p
              e(m)=0.0_DP
              goto 1
            endif
            c=g/r
            g=d(i+1)-p
            s=f/r
            r=(d(i)-g)*s+2.0_DP*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b

            cv(i) = c
            sv(i) = s

          end do

          do i=m-1,l,-1
            do k=1,nrl
              fv2(k)  =z(k,i+1)
            end do
            do k=1,nrl
              fv1(k)  =z(k,i)
            end do
            do k=1,nrl
              z(k,i+1)  =sv(i)*fv1(k) + cv(i)*fv2(k)
              z(k,i)    =cv(i)*fv1(k) - sv(i)*fv2(k)
            end do
          end do

          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0_DP
          goto 1

        endif
15    continue

      DEALLOCATE( cv )
      DEALLOCATE( sv )
      DEALLOCATE( fv1 )
      DEALLOCATE( fv2 )

      return
      END SUBROUTINE ptqliv

!==----------------------------------------------==!


      SUBROUTINE peigsrtv(d,v,ldv,n,nrl)

      USE kinds,     ONLY : DP
!
!     This routine sort eigenvalues and eigenvectors 
!     generated by PTREDV and PTQLIV.  
!
!     AUTHOR : Carlo Cavazzoni - SISSA 1997
!              comments and suggestions to : cava@sissa.it
!

      IMPLICIT NONE
      INTEGER n,ldv,nrl
      REAL(DP) d(n),v(ldv,n)

      INTEGER i,j,k
      REAL(DP) p
      save i,j,k
      save p

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
          do j=1,nrl
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
          END DO

        endif
13    continue
      return
      END SUBROUTINE peigsrtv

   !
   !-------------------------------------------------------------------------
   FUNCTION pythag(a,b)
      USE kinds,     ONLY : DP
      IMPLICIT NONE
      REAL(DP) :: a, b, pythag
      REAL(DP) :: absa, absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.0_DP+(absb/absa)**2)
      else
        if(absb.eq.0.0_DP)then
          pythag=0.0_DP
        else
          pythag=absb*sqrt(1.0_DP+(absa/absb)**2)
        endif
      endif
      return
   END FUNCTION pythag
   !
   !
   !
   !-------------------------------------------------------------------------
   SUBROUTINE diagonalize( iopt, a, lda, d, ev, ldv, n, nproc, mpime, comm_in )
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
     USE kinds,     ONLY : DP
     USE io_global, ONLY : stdout
     USE parallel_include
     !
     IMPLICIT NONE
     !   
     INTEGER,           INTENT(IN)    :: iopt, n, nproc, mpime, lda, ldv
     REAL(DP),          INTENT(OUT)   :: d(n)
     REAL(DP),          INTENT(IN)    :: a(lda,*)
     REAL(DP),          INTENT(OUT)   :: ev(ldv,*)
     INTEGER, OPTIONAL, INTENT(IN)    :: comm_in
     !
     REAL(DP), ALLOCATABLE :: aloc(:,:)
     REAL(DP), ALLOCATABLE :: evloc(:,:)
     REAL(DP), ALLOCATABLE :: e(:)
     INTEGER               :: nrl, i, j, jl, ierr, comm
     INTEGER               :: ip, nrl_ip
     !
     !
#if defined (__PARA) && defined (__MPI)
     IF ( PRESENT( comm_in ) ) THEN
        !
        comm = comm_in
        !
     ELSE
        !
        comm = MPI_COMM_WORLD
        !
     END IF
#endif
     !
     nrl = n / nproc
     !
     IF ( mpime < MOD( n, nproc ) ) nrl = nrl + 1
     !
     ALLOCATE( evloc( nrl, n ) )
     ALLOCATE( e( n ) )
     ALLOCATE( aloc( nrl, n ) )
     !
     IF ( iopt == 1 ) THEN
        !
        DO i = 1, n
           DO jl = 1, nrl
              aloc( jl, i ) = a( (jl-1)*nproc + mpime + 1, i )
           END DO
        END DO
        !
        CALL ptredv( aloc, nrl, d, e, evloc, nrl, nrl, n, nproc, mpime, comm )
        !
     ELSE
        !
        DO i = 1, n
           DO jl = 1, nrl
              aloc( jl, i ) = a( jl, i )
           END DO
        END DO
        !
        CALL ptredv( aloc, nrl, d, e, evloc, nrl, nrl, n, nproc, mpime, comm )
        !
     END IF
     !
     DEALLOCATE( aloc )
     !
     CALL ptqliv( d, e, n, evloc, nrl, nrl )
     CALL peigsrtv( d, evloc, nrl, n, nrl )
     !
     IF ( iopt == 1 ) THEN
        !

#if defined (__PARA)
#  if defined (__MPI)

           DO ip = 0, nproc - 1

              nrl_ip = n / nproc
              !
              IF ( ip < MOD( n, nproc ) )  nrl_ip = nrl_ip + 1

              ALLOCATE( aloc( nrl_ip, n ) )

              IF( mpime == ip ) THEN
                 CALL MPI_BCAST( evloc, nrl_ip*n, MPI_DOUBLE_PRECISION, ip, comm, ierr )
              ELSE
                 CALL MPI_BCAST( aloc, nrl_ip*n, MPI_DOUBLE_PRECISION, ip, comm, ierr )
              END IF

              IF( mpime == ip ) THEN
                 DO j = 1,n
                    DO i = 1, nrl_ip
                       ev((i-1)*nproc+ip+1,j) = evloc(i,j)
                    END DO
                 END DO
              ELSE
                 DO j = 1,n
                    DO i = 1, nrl_ip
                       ev((i-1)*nproc+ip+1,j) = aloc(i,j)
                    END DO
                 END DO
              END IF
              !
              DEALLOCATE( aloc )

           END DO

#  endif
#else

           ev( 1:n, 1:n ) = evloc( 1:n, 1:n )

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
     USE kinds,     ONLY : DP
     USE io_global, ONLY : stdout
     USE parallel_include
     IMPLICIT NONE
     CHARACTER :: JOBZ
     INTEGER   :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, OPTIONAL, INTENT(IN) :: comm_in
     REAL(DP) :: ap( lda, * ), w( * ), z( ldz, * )
     REAL(DP) :: sd( n )
     INTEGER   :: comm
     !
#if defined (__PARA) && defined (__MPI)
     IF ( PRESENT( comm_in ) ) THEN
        comm = comm_in
     ELSE
        comm = MPI_COMM_WORLD
     END IF
#endif
     CALL ptredv(ap, lda, w, sd, z, ldz, nrl, n, nproc, mpime, comm)
     CALL ptqliv(w, sd, n, z, ldz, nrl)
     CALL peigsrtv(w, z, ldz, n, nrl)
     RETURN
   END SUBROUTINE pdspev_drv
 
!==----------------------------------------------==!

      SUBROUTINE dspev_drv( JOBZ, UPLO, N, AP, W, Z, LDZ )
        USE kinds,     ONLY : DP
        USE io_global, ONLY : stdout
        USE parallel_include
        IMPLICIT NONE
        CHARACTER ::       JOBZ, UPLO
        INTEGER   ::       IOPT, INFO, LDZ, N
        REAL(DP) ::  AP( * ), W( * ), Z( LDZ, * )
        REAL(DP), ALLOCATABLE :: WORK(:)

        ALLOCATE( work( 3*n ) )

#if defined __ESSL
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


END MODULE dspev_module
