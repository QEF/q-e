      DOUBLE PRECISION FUNCTION DEBYE3(XVALUE)
C
C
C   DEFINITION:
C
C      This program calculates the Debye function of order 3, defined as
C
C            DEBYE3(x) = 3*[Integral {0 to x} t^3/(exp(t)-1) dt] / (x^3)
C
C      The code uses Chebyshev series whose coefficients
C      are given to 20 decimal places.
C
C
C   ERROR RETURNS:
C
C      If XVALUE < 0.0 an error message is printed and the
C      function returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERMS - INTEGER - The no. of elements of the array ADEB3.
C                         The recommended value is such that
C                             ABS(ADEB3(NTERMS)) < EPS/100,
C                         subject to 1 <= NTERMS <= 18
C
C      XLOW - DOUBLE PRECISION - The value below which
C                    DEBYE3 = 1 - 3x/8 + x*x/20 to machine precision.
C                    The recommended value is
C                        SQRT(8*EPSNEG)
C
C      XUPPER - DOUBLE PRECISION - The value above which
C               DEBYE3 = (18*zeta(4)/x^3) - 3*exp(-x)(x^3+3x^2+6x+6)/x^3.
C                      The recommended value is
C                          -LOG(2*EPS)
C
C      XLIM1 - DOUBLE PRECISION - The value above which DEBYE3 = 18*zeta(4)/x^3
C                     The recommended value is
C                          -LOG(XMIN)
C
C      XLIM2 - DOUBLE PRECISION - The value above which DEBYE3 = 0.0 to machine
C                     precision. The recommended value is
C                          CUBE ROOT(19/XMIN)
C
C      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
C
C      The machine-dependent constants are computed internally by
C      using the D1MACH subroutine.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN, D1MACH
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      AINT , EXP , INT , LOG , SQRT
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley
C          High St.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:  macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:  23 January, 1996
C
      INTEGER I,NEXP,NTERMS
      DOUBLE PRECISION ADEB3(0:18),CHEVAL,DEBINF,EIGHT,EXPMX,FOUR,
     &     HALF,ONE,ONEHUN,PT375,RK,SEVP5,SIX,SUM,T,THREE,TWENTY,X,
     &     XK,XKI,XLIM1,XLIM2,XLOW,XUPPER,XVALUE,ZERO,D1MACH
C
C   OTHER MISCFUN SUBROUTINES USED:
C
       external CHEVAL, D1MACH
C
C
C   INTRINSIC FUNCTIONS USED:
C
       intrinsic  AINT , EXP , INT , LOG , SQRT, abs
C
c*****CHARACTER FNNAME*6,ERRMSG*17
c*****DATA FNNAME/'DEBYE3'/
c*****DATA ERRMSG/'ARGUMENT NEGATIVE'/
      DATA ZERO,PT375/0.0 D 0 , 0.375 D 0/
      DATA HALF,ONE/0.5 D 0 , 1.0 D 0/
      DATA THREE,FOUR,SIX/3.0 D 0 , 4.0 D 0 , 6.0 D 0/
      DATA SEVP5,EIGHT,TWENTY/7.5 D 0 , 8.0 D 0 , 20.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA DEBINF/0.51329 91127 34216 75946 D -1/
      DATA ADEB3/2.70773 70683 27440 94526  D    0,
     1           0.34006 81352 11091 75100  D    0,
     2          -0.12945 15018 44408 6863   D   -1,
     3           0.79637 55380 17381 64     D   -3,
     4          -0.54636 00095 90823 8      D   -4,
     5           0.39243 01959 88049        D   -5,
     6          -0.28940 32823 5386         D   -6,
     7           0.21731 76139 625          D   -7,
     8          -0.16542 09994 98           D   -8,
     9           0.12727 96189 2            D   -9,
     X          -0.98796 3459               D  -11,
     1           0.77250 740                D  -12,
     2          -0.60779 72                 D  -13,
     3           0.48075 9                  D  -14,
     4          -0.38204                    D  -15,
     5           0.3048                     D  -16,
     6          -0.244                      D  -17,
     7           0.20                       D  -18,
     8          -0.2                        D  -19/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
c********CALL ERRPRN(FNNAME,ERRMSG)
         DEBYE3 = ZERO
         RETURN
      ENDIF
C
C   Compute the machine-dependent constants.
C
      T = D1MACH(1)
      XLIM1 = - LOG( T )
      XK = ONE / THREE
      XKI = (ONE/DEBINF) ** XK
      RK = T ** XK
      XLIM2 = XKI / RK
      T = D1MACH(3)
      XLOW = SQRT ( T * EIGHT )
      XUPPER = - LOG( T + T )
      T = T / ONEHUN
      DO 10 NTERMS = 18 , 0 , -1
         IF ( ABS(ADEB3(NTERMS)) .GT. T ) GOTO 19
 10   CONTINUE
C
C   Code for x <= 4.0
C
 19   IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW ) THEN
            DEBYE3 = ( ( X - SEVP5 ) * X + TWENTY ) / TWENTY
         ELSE
            T = ( ( X * X / EIGHT ) - HALF ) - HALF
            DEBYE3 = CHEVAL ( NTERMS , ADEB3 , T ) - PT375 * X
         ENDIF
      ELSE
C
C   Code for x > 4.0
C
         IF ( X .GT. XLIM2 ) THEN
            DEBYE3 = ZERO
         ELSE
            DEBYE3 = ONE / ( DEBINF * X * X * X )
            IF ( X .LT. XLIM1 ) THEN
               EXPMX = EXP ( -X )
               IF ( X .GT. XUPPER ) THEN
                  SUM = (((X+THREE)*X+SIX)*X+SIX) / (X*X*X)
               ELSE
                  SUM = ZERO
                  RK = AINT ( XLIM1 / X )
                  NEXP = INT ( RK )
                  XK = RK * X
                  DO 100 I = NEXP,1,-1
                     XKI = ONE / XK
                     T =  (((SIX*XKI+SIX)*XKI+THREE)*XKI+ONE) / RK
                     SUM = SUM * EXPMX + T
                     RK = RK - ONE
                     XK = XK - X
 100              CONTINUE
               ENDIF
               DEBYE3 = DEBYE3 - THREE * SUM * EXPMX
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
