!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

MODULE dspev_module


    IMPLICIT NONE
    SAVE
    PRIVATE

    PUBLIC :: pdspev_drv, dspev_drv

#if defined __SCALAPACK
    PUBLIC :: pdsyevd_drv
#endif


CONTAINS


    SUBROUTINE ptredv( tv, a, lda, d, e, v, ldv, nrl, n, nproc, me, comm )

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
!     TV       if it is true compute eigrnvectors "v"
!
!     A(NRL,N) Local part of the global matrix A(N,N) to be reduced,
!              only the upper triangle is needed.
!              The rows of the matrix are distributed among processors
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
!     V(NRL,N) Orthogonal transformation that tridiagonalize A,
!              this matrix is distributed among processor
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

      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: tv
      INTEGER, intent(in) :: N, NRL, LDA, LDV
      INTEGER, intent(in) :: NPROC, ME, comm
      REAL(DP) :: A(LDA,N), D(N), E(N), V(LDV,N)
!
      REAL(DP), external ::ddot
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

      ALLOCATE( u( n+2 ), p( n+1 ), vtmp( n+2 ), ul( n ), pl( n ), is( n ), ri( n ) )

      DO I = N, 1, -1
        IS(I)  = (I-1)/NPROC
        RI(I)  = MOD((I-1),NPROC)  !  owner of I-th row
        IF(ME .le. RI(I) ) then
          IS(I) = IS(I) + 1
        END IF
      END DO

      DO I = N, 2, -1 

         L     = I - 1         ! first element
         H     = 0.0_DP

         IF ( L > 1 ) THEN

           SCALEF = 0.0_DP
           DO K = 1, is(l)
              SCALEF = SCALEF + DABS( A(K,I) )
           END DO

#if defined __MPI
           CALL reduce_base_real( 1, scalef, comm, -1 )
#endif

           IF ( SCALEF .EQ. 0.0_DP )  THEN
             !
             IF (RI(L).EQ.ME) THEN
               E(I) = A(is(L),I) 
             END IF
             !
           ELSE 

             !  ......  CALCULATION OF SIGMA AND H

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

             !  CONSTRUCTION OF VECTOR U

             vtmp( 1:l ) = 0.0_DP 

             k = ME + 1
             DO kl = 1,is(l)          
               vtmp(k)   = A(kl,I)
               k         = k + NPROC
             END DO

             DO kl = 1,is(l)          
               UL(kl)    = A(kl,I)
             END DO

#if defined __MPI
             vtmp( l + 1 ) = sigma
             vtmp( l + 2 ) = f
             CALL reduce_base_real_to( L + 2, vtmp, u, comm, -1 )
             sigma = u( l + 1 )
             f     = u( l + 2 )
#else
             u(1:l) = vtmp(1:l)
#endif

             G          = -SIGN(SQRT(SIGMA),F)
             H          = SIGMA - F*G
             ONE_OVER_H = 1.0_DP/H
             E(I)       = SCALEF*G

             U(L)       = F - G 

             IF( RI(L) ==  ME ) THEN
               UL(is(l))  = F - G
               A(is(l),I) = F - G 
             END IF

             !  CONSTRUCTION OF VECTOR P

             DO J = 1,L

               vtmp(j) = 0.0_DP

               DO KL = 1, IS(J)
                 vtmp(J) = vtmp(J) + A(KL,J) * UL(KL)
               END DO

               IF( L > J .AND. ME == RI(J) ) then
                 DO K = J+1,L
                   vtmp(J) = vtmp(J) + A(IS(J),K) * U(K)
                 END DO
               END IF

               vtmp(J) = vtmp(J) * ONE_OVER_H

             END DO 

             KAPPA = 0.5_DP * ONE_OVER_H * ddot( l, vtmp, 1, u, 1 )

#if defined __MPI
             vtmp( l + 1 ) = kappa
             CALL reduce_base_real_to( L + 1, vtmp, p, comm, -1 )
             kappa = p( l + 1 )
#else
             p(1:l) = vtmp(1:l)
#endif

             CALL daxpy( l, -kappa, u, 1, p, 1 )
             CALL DGER( is(l), l, -1.0_DP, ul, 1, p, 1, a, lda )
             CALL DGER( is(l), l, -1.0_DP, p( me + 1 ), nproc, u, 1, a, lda )

           END IF

         ELSE 

           IF(RI(L).EQ.ME) THEN
             G = A(is(l),I)
           END IF

#if defined __MPI
           CALL bcast_real( g, 1, ri( L ), comm )
#endif
           E(I) = G

         END IF

         D(I) = H

      END DO
      
      E(1) = 0.0_DP
      D(1) = 0.0_DP

      IF( tv ) THEN
        DO J = 1,N
          V(1:nrl,J) = 0.0_DP
          IF(RI(J).EQ.ME) THEN
            V(IS(J),J) = 1.0_DP
          END IF
        END DO

        DO I = 2,N
          L = I - 1
          LLOC = IS(L)
          !
          IF( D(I) .NE. 0.0_DP ) THEN
            !
            ONE_OVER_H = 1.0_DP/D(I)
            !
            IF( lloc > 0 ) THEN
               CALL DGEMV( 't', lloc, l, 1.0d0, v(1,1), ldv, a(1,i), 1, 0.0d0, p(1), 1 )
            ELSE
               P(1:l) = 0.0d0
            END IF
           

#if defined __MPI
            CALL reduce_base_real_to( L, p, vtmp, comm, -1 )
#else
            vtmp(1:l) = p(1:l)
#endif

            IF( lloc > 0 ) THEN
               CALL DGER( lloc, l, -ONE_OVER_H, a(1,i), 1, vtmp, 1, v, ldv )
            END IF

          END IF

        END DO 

      END IF


      DO I = 1,N
        U(I) = 0.0_DP
        IF(RI(I).eq.ME) then
          U(I) = A(IS(I),I) 
        END IF
      END DO

#if defined __MPI
      CALL reduce_base_real_to( n, u, d, comm, -1 )
#else
      D(1:N) = U(1:N)
#endif

      DEALLOCATE( u, p, vtmp, ul, pl, is, ri )

    RETURN
    END SUBROUTINE ptredv

!==----------------------------------------------==!

    SUBROUTINE ptqliv( tv, d, e, n, z, ldz, nrl, mpime, comm )

!
! Modified QL algorithm for CRAY T3E PARALLEL MACHINE
! calculate the eigenvectors and eigenvalues of a matrix reduced to 
! tridiagonal form by PTREDV. 
!
! AUTHOR : Carlo Cavazzoni - SISSA 1997
!          comments and suggestions to : carlo.cavazzoni@cineca.it
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
! NOTE : the algorithm that finds the eigenvalues is not parallelized
!        ( it scales as O(N^2) ), I preferred to parallelize only the
!        updating of the eigenvectors because it is the most costly
!        part of the algorithm ( it scales as O(N^3) ).
!        For large matrix in practice all the time is spent in the updating
!        that in this routine scales linearly with the number of processors,
!        in fact there is no communication at all.
!
!  
!     INPUTS :
!
!     TV       if it is true compute eigrnvectors "z"
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
!     Z(LDZ,N) Orthogonal transformation that tridiagonalizes the original
!              matrix A.
!              The rows of the matrix are distributed among processors
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

      IMPLICIT NONE

      LOGICAL, INTENT(IN)  :: tv
      INTEGER, INTENT(IN)  :: n, nrl, ldz, mpime, comm
      REAL(DP) :: d(n), e(n)
      REAL(DP) :: z(ldz,n)

      INTEGER  :: i, iter, mk, k, l, m, ierr
      REAL(DP) :: b, dd, f, g, p, r, c, s
      REAL(DP), ALLOCATABLE :: cv(:,:)
      REAL(DP), ALLOCATABLE :: fv1(:)
      REAL(DP), ALLOCATABLE :: fv2(:)

      ALLOCATE( cv( 2,n ) )
      ALLOCATE( fv1( nrl ) )
      ALLOCATE( fv2( nrl ) )

      do l = 2,n
        e(l-1) = e(l)
      end do
      do l=1,n
        iter=0
1       do m=l,n-1
          dd = abs(d(m))+abs(d(m+1))
          if ( abs(e(m))+dd .eq. dd ) goto 2
        end do
        m=n

2       if ( m /= l ) then
          if ( iter == 200 ) then
             call errore(' tqli ',' too many iterations ', iter)
          end if
          iter=iter+1
          !
          ! iteration is performed on one processor and results broadcast 
          ! to all others to prevent potential problems if all processors
          ! do not behave in exactly the same way (even with the same data!)
          !
          if ( mpime == 0 ) then
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
                if ( r == 0.0_DP) then
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
                !
                cv(1,i-l+1) = c
                cv(2,i-l+1) = s
                !cv(1,i) = c
                !cv(2,i) = s
             end do
             !
             d(l)=d(l)-p
             e(l)=g
             e(m)=0.0_DP
           end if
#if defined __MPI
           CALL bcast_real( cv, 2*(m-l), 0, comm )
           CALL bcast_real( d(l), m-l+1, 0, comm )
           CALL bcast_real( e(l), m-l+1, 0, comm )
#endif

          if( tv ) then
            do i=m-1,l,-1
              do k=1,nrl
                fv2(k)  =z(k,i+1)
              end do
              do k=1,nrl
                fv1(k)  =z(k,i)
              end do
              c = cv(1,i-l+1)
              s = cv(2,i-l+1)
              do k=1,nrl
                z(k,i+1)  =s*fv1(k) + c*fv2(k)
                z(k,i)    =c*fv1(k) - s*fv2(k)
              end do
            end do
          end if

          goto 1

        endif
      end do

      DEALLOCATE( cv )
      DEALLOCATE( fv1 )
      DEALLOCATE( fv2 )

      RETURN
      END SUBROUTINE ptqliv

!==----------------------------------------------==!


      SUBROUTINE peigsrtv(tv,d,v,ldv,n,nrl)

      USE kinds,     ONLY : DP
!
!     This routine sorts eigenvalues and eigenvectors 
!     generated by PTREDV and PTQLIV.  
!
!     AUTHOR : Carlo Cavazzoni - SISSA 1997
!              comments and suggestions to : carlo.cavazzoni@cineca.it
!

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: tv
      INTEGER, INTENT (IN) :: n,ldv,nrl
      REAL(DP), INTENT(INOUT) :: d(n),v(ldv,n)

      INTEGER :: i,j,k
      REAL(DP):: p

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
          if( tv ) then
            do j=1,nrl
              p=v(j,i)
              v(j,i)=v(j,k)
              v(j,k)=p
            END DO
          end if

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
!==----------------------------------------------==!

   SUBROUTINE pdspev_drv( jobz, ap, lda, w, z, ldz, &
                          nrl, n, nproc, mpime, comm )
     USE kinds,     ONLY : DP
     IMPLICIT NONE
     CHARACTER, INTENT(IN) :: JOBZ
     INTEGER, INTENT(IN) :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, INTENT(IN) :: comm
     REAL(DP) :: ap( lda, * ), w( * ), z( ldz, * )
     REAL(DP), ALLOCATABLE :: sd( : )
     LOGICAL :: tv
     !
     IF( n < 1 ) RETURN
     !
     tv = .false.
     IF( jobz == 'V' .OR. jobz == 'v' ) tv = .true.

     ALLOCATE ( sd ( n ) )
     CALL ptredv( tv, ap, lda, w, sd, z, ldz, nrl, n, nproc, mpime, comm)
     CALL ptqliv( tv, w, sd, n, z, ldz, nrl, mpime, comm)
     DEALLOCATE ( sd )
     CALL peigsrtv( tv, w, z, ldz, n, nrl)

     RETURN
   END SUBROUTINE pdspev_drv
 
!==----------------------------------------------==!

      SUBROUTINE dspev_drv( JOBZ, UPLO, N, AP, W, Z, LDZ )
        USE kinds,     ONLY : DP
        IMPLICIT NONE
        CHARACTER ::       JOBZ, UPLO
        INTEGER   ::       IOPT, INFO, LDZ, N
        REAL(DP) ::  AP( * ), W( * ), Z( LDZ, * )
        REAL(DP), ALLOCATABLE :: WORK(:)

        IF( n < 1 ) RETURN

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


#if defined __SCALAPACK

  SUBROUTINE pdsyevd_drv( tv, n, nb, s, lds, w, ortho_cntx )
     USE kinds,     ONLY : DP
     USE mp_bands,  ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, root_bgrp
     USE mp_diag,   ONLY: ortho_comm
     USE mp,        ONLY: mp_comm_free
#if defined(__ELPA)
     USE elpa1
#endif
     IMPLICIT NONE

     LOGICAL, INTENT(IN)  :: tv  
       ! if tv is true compute eigenvalues and eigenvectors (not used)
     INTEGER, INTENT(IN)  :: nb, n, ortho_cntx 
       ! nb = block size, n = matrix size, ortho_cntx = BLACS context
     INTEGER, INTENT(IN)  :: lds
       ! lds = leading dim of s
     REAL(DP) :: s(:,:), w(:)    
       ! input:  s = matrix to be diagonalized
       ! output: s = eigenvectors, w = eigenvalues

     INTEGER     :: desch( 10 )
     REAL(DP)    :: rtmp( 4 )
     INTEGER     :: itmp( 4 )
     REAL(DP), ALLOCATABLE :: work(:)
     REAL(DP), ALLOCATABLE :: vv(:,:)
     INTEGER,  ALLOCATABLE :: iwork(:)
     INTEGER     :: LWORK, LIWORK, info
     CHARACTER   :: jobv
     INTEGER     :: i
#if defined(__ELPA)
     INTEGER     :: nprow,npcol,my_prow, my_pcol,mpi_comm_rows, mpi_comm_cols
#endif 

     IF( SIZE( s, 1 ) /= lds ) &
        CALL errore( ' pdsyevd_drv ', ' wrong matrix leading dimension ', 1 )
     !
     IF( tv ) THEN
        ALLOCATE( vv( SIZE( s, 1 ), SIZE( s, 2 ) ) )
        jobv = 'V'
     ELSE
        CALL errore('pdsyevd_drv','PDSYEVD does not compute eigenvalue only',1)
     END IF

     CALL descinit( desch, n, n, nb, nb, 0, 0, ortho_cntx, SIZE( s, 1 ) , info )

     IF( info /= 0 ) CALL errore( ' pdsyevd_drv ', ' desckinit ', ABS( info ) )

     lwork = -1
     liwork = 1
     itmp = 0
     rtmp = 0.0_DP

#if defined(__ELPA)
     CALL BLACS_Gridinfo(ortho_cntx,nprow, npcol, my_prow,my_pcol)
     CALL get_elpa_row_col_comms(ortho_comm, my_prow, my_pcol,mpi_comm_rows, mpi_comm_cols)
     CALL solve_evp_real(n,  n,   s, lds,    w,  vv, lds     ,nb  ,mpi_comm_rows, mpi_comm_cols)
     
     IF( tv )  s = vv
     IF( ALLOCATED( vv ) ) DEALLOCATE( vv )

     CALL mp_comm_free ( mpi_comm_rows )
     CALL mp_comm_free ( mpi_comm_cols )
#else
     CALL PDSYEVD( jobv, 'L', n, s, 1, 1, desch, w, vv, 1, 1, desch, rtmp, lwork, itmp, liwork, info )

     IF( info /= 0 ) CALL errore( ' pdsyevd_drv ', ' PDSYEVD ', ABS( info ) )

     lwork  = MAX( 131072, 2*INT( rtmp(1) ) + 1 )
     liwork = MAX( 8*n , itmp(1) + 1 )

     ALLOCATE( work( lwork ) )
     ALLOCATE( iwork( liwork ) )

     CALL PDSYEVD( jobv, 'L', n, s, 1, 1, desch, w, vv, 1, 1, desch, work, lwork, iwork, liwork, info )

     IF( info /= 0 ) CALL errore( ' pdsyevd_drv ', ' PDSYEVD ', ABS( info ) )

     IF( tv ) s = vv

     IF( ALLOCATED( vv ) ) DEALLOCATE( vv )
     DEALLOCATE( work )
     DEALLOCATE( iwork )
#endif 

     RETURN
  END SUBROUTINE pdsyevd_drv

#endif


END MODULE dspev_module
