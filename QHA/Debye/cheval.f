      DOUBLE PRECISION FUNCTION CHEVAL(N,A,T)
C
C   This function evaluates a Chebyshev series, using the
C   Clenshaw method with Reinsch modification, as analysed
C   in the paper by Oliver.
C
C   INPUT PARAMETERS
C
C       N - INTEGER - The no. of terms in the sequence
C
C       A - DOUBLE PRECISION ARRAY, dimension 0 to N - The coefficients of
C           the Chebyshev series
C
C       T - DOUBLE PRECISION - The value at which the series is to be
C           evaluated
C
C
C   REFERENCES
C
C        "An error analysis of the modified Clenshaw method for
C         evaluating Chebyshev and Fourier series" J. Oliver,
C         J.I.M.A., vol. 20, 1977, pp379-391
C
C
C MACHINE-DEPENDENT CONSTANTS: NONE
C
C
C INTRINSIC FUNCTIONS USED;
C
C    ABS
C
C
C AUTHOR:  Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley ,
C          High St.,
C          PAISLEY,
C          SCOTLAND
C
C
C LATEST MODIFICATION:   21 December , 1992
C
C
      INTEGER I,N
      DOUBLE PRECISION A(0:N),D1,D2,HALF,T,TEST,TT,TWO,U0,U1,U2,ZERO
      INTRINSIC ABS
      DATA ZERO,HALF/ 0.0 D 0 , 0.5 D 0 /
      DATA TEST,TWO/ 0.6 D 0 , 2.0 D 0 /
      U1 = ZERO
C
C   If ABS ( T )  < 0.6 use the standard Clenshaw method
C
      IF ( ABS( T ) .LT. TEST ) THEN
         U0 = ZERO
         TT = T + T
         DO 100 I = N , 0 , -1
            U2 = U1
            U1 = U0
            U0 = TT * U1 + A( I ) - U2
 100     CONTINUE
         CHEVAL =  ( U0 - U2 ) / TWO
      ELSE
C
C   If ABS ( T )  > =  0.6 use the Reinsch modification
C
         D1 = ZERO
C
C   T > =  0.6 code
C
         IF ( T .GT. ZERO ) THEN
            TT =  ( T - HALF ) - HALF
            TT = TT + TT
            DO 200 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) + D2
               U1 = D1 + U2
 200        CONTINUE
            CHEVAL =  ( D1 + D2 ) / TWO
         ELSE
C
C   T < =  -0.6 code
C
            TT =  ( T + HALF ) + HALF
            TT = TT + TT
            DO 300 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) - D2
               U1 = D1 - U2
 300        CONTINUE
            CHEVAL =  ( D1 - D2 ) / TWO
         ENDIF
      ENDIF
      RETURN
      END

