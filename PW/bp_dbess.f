C
C-------------------------------------------------------------------------
      SUBROUTINE DBESS(XG,L,MMAX,R,DJL)
C-------------------------------------------------------------------------
C     CALCULATES DERIVATIVES OF SPHERICAL BESSEL FUNCTIONS  j_l(Gr)
C     WITH RESPECT TO h_alpha,beta (WITHOUT THE FACTOR GAGK(KK,IG)*HTM1)
C     I.E. -x * D(jl(x))/dx
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(EPS=1.E-8)
      REAL*8   DJL(MMAX),R(MMAX)
 
      IF(L.EQ.1) THEN                      !   S  PART
        IF(XG.LT.EPS) THEN
          DO IR=1,MMAX
            DJL(IR) = 0.D0
          END DO
        ELSE
          DJL(1) = 0.D0
          DO IR=2,MMAX
            XRG=R(IR)*XG
            DJL(IR) = SIN(XRG)/XRG-COS(XRG)
          END DO
        ENDIF
      ENDIF
 
      IF(L.EQ.2) THEN                      !   P  PART
        IF(XG.LT.EPS) THEN
          DO IR=1,MMAX
            DJL(IR) = 0.D0
          END DO
        ELSE
          DJL(1) = 0.D0
          DO IR=2,MMAX
            XRG=R(IR)*XG
            DJL(IR) = 2.D0*(SIN(XRG)/XRG-COS(XRG))/XRG - SIN(XRG)
          END DO
        ENDIF
      ENDIF
 
      IF(L.EQ.3) THEN                      !   D  PART
        IF(XG.LT.EPS) THEN
          DO IR=1,MMAX
            DJL(IR) = 0.D0
          END DO
        ELSE
          DJL(1) = 0.D0
          DO IR=2,MMAX
            XRG=R(IR)*XG
            DJL(IR) = ( SIN(XRG)*(9.D0/(XRG*XRG)-4.D0) -
     -                9.D0*COS(XRG)/XRG ) /XRG + COS(XRG)
          END DO
        ENDIF
      ENDIF

      IF(L.EQ.4) THEN                      !   F  PART
        IF(XG.LT.EPS) THEN
          DO IR=1,MMAX
            DJL(IR) = 0.D0
          END DO
        ELSE
          DJL(1) = 0.D0
          DO IR=2,MMAX
            XRG=R(IR)*XG
            XRG2=XRG*XRG
            DJL(IR)=SIN(XRG)*(60.D0/(XRG2*XRG2)-27.D0/XRG2+1.d0)
     $           -COS(XRG)*(60.D0/XRG2-7.D0)/XRG
          END DO
        ENDIF
      ENDIF

      IF(L.EQ.5) THEN                      !   G  PART
        IF(XG.LT.EPS) THEN
          DO IR=1,MMAX
            DJL(IR) = 0.D0
          END DO
        ELSE
          DJL(1) = 0.D0
          DO IR=2,MMAX
            XRG=R(IR)*XG
            XRG2=XRG*XRG
            DJL(IR)=SIN(XRG)*(525.D0/(XRG2*XRG2)-240.D0/XRG2+11.D0)/XRG
     $           -  COS(XRG)*(525.D0/(XRG2*XRG2)-65.D0/XRG2+1.D0)
          END DO
        ENDIF
      ENDIF
 
      IF(L.LE.0 .OR. L.GE.6) THEN
        CALL ERRORE('DBESS',' L NOT PROGRAMMED, L= ',L)
      END IF

      RETURN
      END
