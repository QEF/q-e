      FUNCTION RRANF()
        USE kinds
        IMPLICIT NONE
        REAL(dbl) :: RRANF
        INTEGER :: KK, M, KONST
        DATA M/100001/, KONST/125/
        SAVE M, KONST
          M=M*KONST
          M=M-2796203 * (M/2796203)
          RRANF = REAL(M)/2796203.D0
        RETURN
      END FUNCTION
