! Slightly modified version of LINPACK routines zgefa and zgedi

      SUBROUTINE ZGEFA(A,LDA,N,IPVT,INFO)
      USE kinds, ONLY : DP
      INTEGER LDA,N,IPVT(*),INFO
      COMPLEX(DP) A(LDA,*)
!
!     ZGEFA FACTORS A COMPLEX(DP) MATRIX BY GAUSSIAN ELIMINATION.
!
!     ZGEFA IS USUALLY CALLED BY ZGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR ZGECO) = (1 + 9/N)*(TIME FOR ZGEFA) .
!
!     ON ENTRY
!
!        A       COMPLEX(DP)(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT ZGESL OR ZGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN ZGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS ZAXPY,ZSCAL,IZAMAX
!     FORTRAN DABS
!
!     INTERNAL VARIABLES
!
      COMPLEX(DP) T
      INTEGER IZAMAX,J,K,KP1,L,NM1
!
      COMPLEX(DP) ZDUM
      REAL(DP) CABS1
      REAL(DP) REAL,AIMAG
      COMPLEX(DP) ZDUMR,ZDUMI
      REAL(ZDUMR) = ZDUMR
      AIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CABS1(ZDUM) = DABS(REAL(ZDUM)) + DABS(AIMAG(ZDUM))
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
!
!        FIND L = PIVOT INDEX
!
         L = IZAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
         IF (CABS1(A(L,K)) .EQ. 0.0D0) GO TO 40
!
!           INTERCHANGE IF NECESSARY
!
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
!
!           COMPUTE MULTIPLIERS
!
            T = -(1.0D0,0.0D0)/A(K,K)
            CALL ZSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL ZAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (CABS1(A(N,N)) .EQ. 0.0D0) INFO = N
      RETURN
      END SUBROUTINE ZGEFA

      SUBROUTINE ZGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      USE kinds, ONLY : DP
      INTEGER LDA,N,IPVT(*),JOB
      COMPLEX(DP) A(LDA,*),DET(2),WORK(*)
!
!     ZGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
!     USING THE FACTORS COMPUTED BY ZGECO OR ZGEFA.
!
!     ON ENTRY
!
!        A       COMPLEX(DP)(LDA, N)
!                THE OUTPUT FROM ZGECO OR ZGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM ZGECO OR ZGEFA.
!
!        WORK    COMPLEX(DP)(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     COMPLEX(DP)(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF ZGECO HAS SET RCOND .GT. 0.0 OR ZGEFA HAS SET
!        INFO .EQ. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS ZAXPY,ZSCAL,ZSWAP
!     FORTRAN DABS,CMPLX,MOD
!
!     INTERNAL VARIABLES
!
      COMPLEX(DP) T
      REAL(DP) TEN
      INTEGER I,J,K,KB,KP1,L,NM1
!
      COMPLEX(DP) ZDUM
      REAL(DP) CABS1
      REAL(DP) REAL,AIMAG
      COMPLEX(DP) ZDUMR,ZDUMI
      REAL(ZDUMR) = ZDUMR
      AIMAG(ZDUMI) = (0.0D0,-1.0D0)*ZDUMI
      CABS1(ZDUM) = DABS(REAL(ZDUM)) + DABS(AIMAG(ZDUM))
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = (1.0D0,0.0D0)
         DET(2) = (0.0D0,0.0D0)
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
!        ...EXIT
            IF (CABS1(DET(1)) .EQ. 0.0D0) GO TO 60
   10       IF (CABS1(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = CMPLX(TEN,0.0D0,KIND=dp)*DET(1)
               DET(2) = DET(2) - (1.0D0,0.0D0)
            GO TO 10
   20       CONTINUE
   30       IF (CABS1(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/CMPLX(TEN,0.0D0,KIND=dp)
               DET(2) = DET(2) + (1.0D0,0.0D0)
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(U)
!
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = (1.0D0,0.0D0)/A(K,K)
            T = -A(K,K)
            CALL ZSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0D0,0.0D0)
               CALL ZAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM INVERSE(U)*INVERSE(L)
!
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = (0.0D0,0.0D0)
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL ZAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL ZSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END SUBROUTINE ZGEDI
