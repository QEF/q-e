      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. SCALAR ARGUMENTS ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
*
*     C := ALPHA*OP( A )*OP( B ) + BETA*C,
*
*  WHERE  OP( X ) IS ONE OF
*
*     OP( X ) = X   OR   OP( X ) = X',
*
*  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
*  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  TRANSA - CHARACTER*1.
*           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
*           THE MATRIX MULTIPLICATION AS FOLLOWS:
*
*              TRANSA = 'N' OR 'N',  OP( A ) = A.
*
*              TRANSA = 'T' OR 'T',  OP( A ) = A'.
*
*              TRANSA = 'C' OR 'C',  OP( A ) = A'.
*
*           UNCHANGED ON EXIT.
*
*  TRANSB - CHARACTER*1.
*           ON ENTRY, TRANSB SPECIFIES THE FORM OF OP( B ) TO BE USED IN
*           THE MATRIX MULTIPLICATION AS FOLLOWS:
*
*              TRANSB = 'N' OR 'N',  OP( B ) = B.
*
*              TRANSB = 'T' OR 'T',  OP( B ) = B'.
*
*              TRANSB = 'C' OR 'C',  OP( B ) = B'.
*
*           UNCHANGED ON EXIT.
*
*  M      - INTEGER.
*           ON ENTRY,  M  SPECIFIES  THE NUMBER  OF ROWS  OF THE  MATRIX
*           OP( A )  AND OF THE  MATRIX  C.  M  MUST  BE AT LEAST  ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY,  N  SPECIFIES THE NUMBER  OF COLUMNS OF THE MATRIX
*           OP( B ) AND THE NUMBER OF COLUMNS OF THE MATRIX C. N MUST BE
*           AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  K      - INTEGER.
*           ON ENTRY,  K  SPECIFIES  THE NUMBER OF COLUMNS OF THE MATRIX
*           OP( A ) AND THE NUMBER OF ROWS OF THE MATRIX OP( B ). K MUST
*           BE AT LEAST  ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, KA ), WHERE KA IS
*           K  WHEN  TRANSA = 'N' OR 'N',  AND IS  M  OTHERWISE.
*           BEFORE ENTRY WITH  TRANSA = 'N' OR 'N',  THE LEADING  M BY K
*           PART OF THE ARRAY  A  MUST CONTAIN THE MATRIX  A,  OTHERWISE
*           THE LEADING  K BY M  PART OF THE ARRAY  A  MUST CONTAIN  THE
*           MATRIX A.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. WHEN  TRANSA = 'N' OR 'N' THEN
*           LDA MUST BE AT LEAST  MAX( 1, M ), OTHERWISE  LDA MUST BE AT
*           LEAST  MAX( 1, K ).
*           UNCHANGED ON EXIT.
*
*  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, KB ), WHERE KB IS
*           N  WHEN  TRANSB = 'N' OR 'N',  AND IS  K  OTHERWISE.
*           BEFORE ENTRY WITH  TRANSB = 'N' OR 'N',  THE LEADING  K BY N
*           PART OF THE ARRAY  B  MUST CONTAIN THE MATRIX  B,  OTHERWISE
*           THE LEADING  N BY K  PART OF THE ARRAY  B  MUST CONTAIN  THE
*           MATRIX B.
*           UNCHANGED ON EXIT.
*
*  LDB    - INTEGER.
*           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. WHEN  TRANSB = 'N' OR 'N' THEN
*           LDB MUST BE AT LEAST  MAX( 1, K ), OTHERWISE  LDB MUST BE AT
*           LEAST  MAX( 1, N ).
*           UNCHANGED ON EXIT.
*
*  BETA   - DOUBLE PRECISION.
*           ON ENTRY,  BETA  SPECIFIES THE SCALAR  BETA.  WHEN  BETA  IS
*           SUPPLIED AS ZERO THEN C NEED NOT BE SET ON INPUT.
*           UNCHANGED ON EXIT.
*
*  C      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDC, N ).
*           BEFORE ENTRY, THE LEADING  M BY N  PART OF THE ARRAY  C MUST
*           CONTAIN THE MATRIX  C,  EXCEPT WHEN  BETA  IS ZERO, IN WHICH
*           CASE C NEED NOT BE SET ON ENTRY.
*           ON EXIT, THE ARRAY  C  IS OVERWRITTEN BY THE  M BY N  MATRIX
*           ( ALPHA*OP( A )*OP( B ) + BETA*C ).
*
*  LDC    - INTEGER.
*           ON ENTRY, LDC SPECIFIES THE FIRST DIMENSION OF C AS DECLARED
*           IN  THE  CALLING  (SUB)  PROGRAM.   LDC  MUST  BE  AT  LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 3 BLAS ROUTINE.
*
*  -- WRITTEN ON 8-FEBRUARY-1989.
*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
*     IAIN DUFF, AERE HARWELL.
*     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
*     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
*
*
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     .. LOCAL SCALARS ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     SET  NOTA  AND  NOTB  AS  TRUE IF  A  AND  B  RESPECTIVELY ARE NOT
*     TRANSPOSED AND SET  NROWA, NCOLA AND  NROWB  AS THE NUMBER OF ROWS
*     AND  COLUMNS OF  A  AND THE  NUMBER OF  ROWS  OF  B  RESPECTIVELY.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     AND IF  ALPHA.EQ.ZERO.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     START THE OPERATIONS.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           FORM  C := ALPHA*A*B + BETA*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           FORM  C := ALPHA*A'*B + BETA*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           FORM  C := ALPHA*A*B' + BETA*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           FORM  C := ALPHA*A'*B' + BETA*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     END OF DGEMM .
*
      END
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGEMV  PERFORMS ONE OF THE MATRIX-VECTOR OPERATIONS
*
*     Y := ALPHA*A*X + BETA*Y,   OR   Y := ALPHA*A'*X + BETA*Y,
*
*  WHERE ALPHA AND BETA ARE SCALARS, X AND Y ARE VECTORS AND A IS AN
*  M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  TRANS  - CHARACTER*1.
*           ON ENTRY, TRANS SPECIFIES THE OPERATION TO BE PERFORMED AS
*           FOLLOWS:
*
*              TRANS = 'N' OR 'N'   Y := ALPHA*A*X + BETA*Y.
*
*              TRANS = 'T' OR 'T'   Y := ALPHA*A'*X + BETA*Y.
*
*              TRANS = 'C' OR 'C'   Y := ALPHA*A'*X + BETA*Y.
*
*           UNCHANGED ON EXIT.
*
*  M      - INTEGER.
*           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
*           M MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
*           CONTAIN THE MATRIX OF COEFFICIENTS.
*           UNCHANGED ON EXIT.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCX ) ) WHEN TRANS = 'N' OR 'N'
*           AND AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCX ) ) OTHERWISE.
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE
*           VECTOR X.
*           UNCHANGED ON EXIT.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  BETA   - DOUBLE PRECISION.
*           ON ENTRY, BETA SPECIFIES THE SCALAR BETA. WHEN BETA IS
*           SUPPLIED AS ZERO THEN Y NEED NOT BE SET ON INPUT.
*           UNCHANGED ON EXIT.
*
*  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCY ) ) WHEN TRANS = 'N' OR 'N'
*           AND AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCY ) ) OTHERWISE.
*           BEFORE ENTRY WITH BETA NON-ZERO, THE INCREMENTED ARRAY Y
*           MUST CONTAIN THE VECTOR Y. ON EXIT, Y IS OVERWRITTEN BY THE
*           UPDATED VECTOR Y.
*
*  INCY   - INTEGER.
*           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           Y. INCY MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 2 BLAS ROUTINE.
*
*  -- WRITTEN ON 22-OCTOBER-1986.
*     JACK DONGARRA, ARGONNE NATIONAL LAB.
*     JEREMY DU CROZ, NAG CENTRAL OFFICE.
*     SVEN HAMMARLING, NAG CENTRAL OFFICE.
*     RICHARD HANSON, SANDIA NATIONAL LABS.
*
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     SET  LENX  AND  LENY, THE LENGTHS OF THE VECTORS X AND Y, AND SET
*     UP THE START POINTS IN  X  AND  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
*
*     FIRST FORM  Y := BETA*Y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        FORM  Y := ALPHA*A*X + Y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        FORM  Y := ALPHA*A'*X + Y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     END OF DGEMV .
*
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGER   PERFORMS THE RANK 1 OPERATION
*
*     A := ALPHA*X*Y' + A,
*
*  WHERE ALPHA IS A SCALAR, X IS AN M ELEMENT VECTOR, Y IS AN N ELEMENT
*  VECTOR AND A IS AN M BY N MATRIX.
*
*  PARAMETERS
*  ==========
*
*  M      - INTEGER.
*           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
*           M MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  N      - INTEGER.
*           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
*           N MUST BE AT LEAST ZERO.
*           UNCHANGED ON EXIT.
*
*  ALPHA  - DOUBLE PRECISION.
*           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
*           UNCHANGED ON EXIT.
*
*  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( M - 1 )*ABS( INCX ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE M
*           ELEMENT VECTOR X.
*           UNCHANGED ON EXIT.
*
*  INCX   - INTEGER.
*           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           X. INCX MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
*           ( 1 + ( N - 1 )*ABS( INCY ) ).
*           BEFORE ENTRY, THE INCREMENTED ARRAY Y MUST CONTAIN THE N
*           ELEMENT VECTOR Y.
*           UNCHANGED ON EXIT.
*
*  INCY   - INTEGER.
*           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
*           Y. INCY MUST NOT BE ZERO.
*           UNCHANGED ON EXIT.
*
*  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
*           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
*           CONTAIN THE MATRIX OF COEFFICIENTS. ON EXIT, A IS
*           OVERWRITTEN BY THE UPDATED MATRIX.
*
*  LDA    - INTEGER.
*           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
*           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
*           MAX( 1, M ).
*           UNCHANGED ON EXIT.
*
*
*  LEVEL 2 BLAS ROUTINE.
*
*  -- WRITTEN ON 22-OCTOBER-1986.
*     JACK DONGARRA, ARGONNE NATIONAL LAB.
*     JEREMY DU CROZ, NAG CENTRAL OFFICE.
*     SVEN HAMMARLING, NAG CENTRAL OFFICE.
*     RICHARD HANSON, SANDIA NATIONAL LABS.
*
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
*     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     END OF DGER  .
*
      END
