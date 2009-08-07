        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [-1,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        N0=(N+1)/2
        DO 45 NR=1,N0
           Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10         Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
           F1=Z
           DO 20 K=2,N
              PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
              PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
              F0=F1
20            F1=PF
           IF (Z.EQ.0.0) GO TO 40
           FD=PF/P
           Q=0.0D0
           DO 35 I=1,NR
              WP=1.0D0
              DO 30 J=1,NR
                 IF (J.NE.I) WP=WP*(Z-X(J))
30            CONTINUE
35            Q=Q+WP
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40         X(NR)=Z
           X(N+1-NR)=-Z
           W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
45         W(N+1-NR)=W(NR)
        RETURN
        END

