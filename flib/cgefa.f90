#include"../include/fpmd.h"

      SUBROUTINE FPMD_ZGEFA(A,LDA,N,IPVT,INFO)
      USE kinds
      INTEGER LDA,N,IPVT(1),INFO
      COMPLEX(dbl) A(LDA,1)
!
!     ZGEFA FACTORS A COMPLEX(dbl) MATRIX BY GAUSSIAN ELIMINATION.
!
!     ZGEFA IS USUALLY CALLED BY ZGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR ZGECO) = (1 + 9/N)*(TIME FOR ZGEFA) .
!
!     ON ENTRY
!
!        A       COMPLEX(dbl)(LDA, N)
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
      COMPLEX(dbl) T
      INTEGER IZAMAX,J,K,KP1,L,NM1
!
      COMPLEX(dbl) ZDUM
      REAL(dbl) CABS1
      REAL(dbl) REAL,AIMAG
      COMPLEX(dbl) ZDUMR,ZDUMI
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
         L = ICAMAX(N-K+1,A(K,K),1) + K - 1
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
            CALL FPMD_ZSCAL(N-K,T,A(K+1,K),1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL FPMD_ZAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
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
      END
