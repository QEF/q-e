C-------------------------------------------------------------------------
      SUBROUTINE BESS(XG,L,MMAX,R,JL)
C-------------------------------------------------------------------------
C     CALCULATES SPHERICAL BESSEL FUNCTIONS  j_l(Gr)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(EPS=1.E-8)
      REAL*8   JL(MMAX),R(MMAX)
 
      IF(L.EQ.1) THEN                      !   S  PART
        IF(XG.LT.EPS) THEN
         DO 41 IR=1,MMAX
 41        JL(IR)=1.D0
        ELSE
          JL(1)=1.D0
          DO 42 IR=2,MMAX
            XRG=R(IR)*XG
            JL(IR)=SIN(XRG)/XRG
 42       CONTINUE
        ENDIF
      ENDIF
 
      IF(L.EQ.2) THEN                      !   P  PART
        IF(XG.LT.EPS) THEN
        DO 43 IR=1,MMAX
 43        JL(IR)=0.D0
        ELSE
          JL(1)=0.
          DO 44 IR=2,MMAX
            XRG=R(IR)*XG
            JL(IR)=(SIN(XRG)/XRG-COS(XRG))/XRG
 44       CONTINUE
        ENDIF
      ENDIF
 
      IF(L.EQ.3) THEN                      !   D  PART
        IF(XG.LT.EPS) THEN
        DO 45 IR=1,MMAX
 45       JL(IR)=0.D0
        ELSE
          JL(1)=0.D0
          DO 46 IR=2,MMAX
           XRG=R(IR)*XG
           JL(IR)=(SIN(XRG)*(3./(XRG*XRG)-1.)
     +             -3.*COS(XRG)/XRG) /XRG
 46       CONTINUE
        ENDIF
      ENDIF
 
      IF(L.EQ.4) THEN                      !   F  PART
        IF(XG.LT.EPS) THEN
        DO 47 IR=1,MMAX
 47        JL(IR)=0.D0
        ELSE
          JL(1)=0.D0
          DO 48 IR=2,MMAX
           XRG=R(IR)*XG
           XRG2=XRG*XRG
           JL(IR)=( SIN(XRG)*(15./(XRG2*XRG)-6./XRG)
     +             +COS(XRG)*(1.-15./XRG2) )/XRG
 48       CONTINUE
        ENDIF
      ENDIF
 
      IF(L.EQ.5) THEN                      !   G  PART
        IF(XG.LT.EPS) THEN
        DO 49 IR=1,MMAX
 49        JL(IR)=0.D0
        ELSE
          JL(1)=0.D0
          DO 50 IR=2,MMAX
           XRG=R(IR)*XG
           XRG2=XRG*XRG
           JL(IR)=( SIN(XRG)*(105./(XRG2*XRG2)-45./XRG2+1.)
     +             +COS(XRG)*(10./XRG-105./(XRG2*XRG)) )/XRG
 50       CONTINUE
        ENDIF
      ENDIF
 
      RETURN
      END
