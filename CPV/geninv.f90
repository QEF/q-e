!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sat Nov 13 11:04:11 MET 1999
!  ----------------------------------------------

#include "f_defs.h"

!  routines in this file:
!  SUBROUTINE geninv(a,ld,n,mrank,cond,u,v,work,toleig,info,iopt)
!  SUBROUTINE zgeninv(a,ld,n,mrank,cond,u,v,work,toleig,info,iopt)
!  ----------------------------------------------
!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE geninv( a, ld, n, mrank, cond, u, v, &
        work, toleig, info, iopt)

!  get a general inverse matrix
!
!  iopt = 0:   using calculated rank for pseudo inverse
!  iopt = 1 :  rank is assumed to be n-6 for pseudo inverse
!  iopt > 10:  scale matrix before decomposition
!  ----------------------------------------------
!  END manual
!  end of declarations
!  ----------------------------------------------

      USE kinds
      IMPLICIT NONE

      INTEGER :: ld, n, mrank, info, iopt
      REAL(DP) :: cond, toleig

      REAL(DP)  :: a(ld,n), u(ld,n), v(ld,n), work(4*n)
      REAL(DP)  :: zero=0.0d0
      REAL(DP)  :: one=1.0d0

      INTEGER :: n3, i, j, k, m

! ... scale matrix before inversion
      IF ( iopt >= 10 ) THEN
        n3 = 3 * n
        DO i = 1, n
          work( n3 + i ) = one
          IF ( ABS( a(i,i) ) >= 1.d-13 ) &
             work( n3 + i ) = one / SQRT( ABS( a(i,i) ) )
        END DO
        DO i = 1, n
          DO j = 1, n
            a(i,j) = a(i,j) * work( n3 + i ) * work( n3 + j )
          END DO
        END DO
      END IF

! ... get singular values
      CALL dsvdc( a, ld, n, n, work, work(n+1), u, ld, v, ld, work(2*n+1), 11, info)
      mrank = 0
      DO i = 1, n
        IF ( ABS( work(i) ) > toleig ) mrank = mrank + 1
      END DO

      m = mrank
      IF ( iopt == 1 .OR. iopt == 11 ) m = n - 6
      cond = work(1) / work(m)
      DO i = 1, m
        work(i) = one / work(i)
      END DO

      DO j = 1, n
        DO i = 1, n
          a(i,j) = zero
        END DO
        DO k = 1, m
          DO i = 1, n
            a(i,j) = a(i,j) + v(i,k) * work(k) * u(j,k)
          END DO
        END DO
      END DO

! ... rescale matrix after inversion
      IF ( iopt >= 10 ) THEN
        DO i = 1, n
          DO j = 1, n
            a(i,j) = a(i,j) * work( n3 + i ) * work( n3 + j )
          END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE geninv

!  ----------------------------------------------
!  ----------------------------------------------
!  BEGIN manual
      SUBROUTINE zgeninv(a,ld,n,mrank,cond,u,v,work,toleig,info,iopt)

!  get a general inverse matrix
!
!  iopt = 0:   using calculated rank for pseudo inverse
!  iopt = 1 :  rank is assumed to be n-6 for pseudo inverse
!  iopt > 10:  scale matrix before decomposition
!  ----------------------------------------------
!  END manual
!  end of declarations
!  ----------------------------------------------

      USE kinds
      IMPLICIT NONE

      INTEGER :: ld, n, mrank, info, iopt
      REAL(DP) :: cond, toleig

      COMPLEX(DP) a(ld,n),u(ld,n),v(ld,n),work(4*n)
      REAL(DP)  :: zero=0.0d0
      REAL(DP)  :: one=1.0d0

      INTEGER :: n3, i, j, k, m 

! ... scale matrix before inversion
      IF (iopt.GE.10) THEN
        n3=3*n
        DO i=1,n
          work(n3+i)=one
          IF (abs(a(i,i)).GE.1.d-13) &
             work(n3+i)=one/dsqrt(abs(a(i,i)))
        END DO
        DO i=1,n
          DO j=1,n
            a(i,j)=a(i,j)*work(n3+i)*work(n3+j)
          END DO
        END DO
      END IF

! ... get singular values
      CALL dsvdc(a,ld,n,n,work,work(n+1),u,ld,v,ld,work(2*n+1),11,info)
      mrank=0
      DO i=1,n
        IF (abs(work(i)).GT.toleig) mrank=mrank+1
      END DO

      m=mrank
      IF (iopt.EQ.1.OR.iopt.EQ.11) m=n-6
      cond=work(1)/work(m)
      DO i=1,m
        work(i)=one/work(i)
      END DO

      DO j=1,n
        DO i=1,n
          a(i,j)=zero
        END DO
        DO k=1,m
          DO i=1,n
            a(i,j)=a(i,j)+v(i,k)*work(k)*u(j,k)
          END DO
        END DO
      END DO

! ... rescale matrix after inversion
      IF (iopt.GE.10) THEN
        DO i=1,n
          DO j=1,n
            a(i,j)=a(i,j)*work(n3+i)*work(n3+j)
          END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE zgeninv

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)

!  (describe briefly what this routine does...)
!  ----------------------------------------------

!     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X
!     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
!     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
!     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
!     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
!
!     ON ENTRY
!
!         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N.
!                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
!                   DECOMPOSITION IS TO BE COMPUTED.  X IS
!                   DESTROYED BY DSVDC.
!
!         LDX       INTEGER.
!                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!
!         N         INTEGER.
!                   N IS THE NUMBER OF ROWS OF THE MATRIX X.
!
!         P         INTEGER.
!                   P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
!
!         LDU       INTEGER.
!                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
!                   (SEE BELOW).
!
!         LDV       INTEGER.
!                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
!                   (SEE BELOW).
!
!         WORK      DOUBLE PRECISION(N).
!                   WORK IS A SCRATCH ARRAY.
!
!         JOB       INTEGER.
!                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
!                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
!                   WITH THE FOLLOWING MEANING
!
!                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
!                                  VECTORS.
!                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
!                                  IN U.
!                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
!                                  VECTORS IN U.
!                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
!                                  VECTORS.
!                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
!                                  IN V.
!
!     ON RETURN
!
!         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P).
!                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
!                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
!                   ORDER OF MAGNITUDE.
!
!         E         DOUBLE PRECISION(P),
!                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
!                   DISCUSSION OF INFO FOR EXCEPTIONS.
!
!         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF
!                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2
!                                   THEN K.EQ.MIN(N,P).
!                   U CONTAINS THE MATRIX OF LEFT SINGULAR VECTORS.
!                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
!                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
!                   IN THE SUBROUTINE CALL.
!
!         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P.
!                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
!                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
!                   THEN V MAY BE IDENTIFIED WITH X IN THE
!                   SUBROUTINE CALL.
!
!         INFO      INTEGER.
!                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
!                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
!                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
!                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
!                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
!                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
!                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
!                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
!                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
!                   VALUES OF X AND B ARE THE SAME.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!              CORRECTION MADE TO SHIFT 2/84.
!     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!
!     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!
!     EXTERNAL DROT
!     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG
!     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT
!

      USE kinds
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      REAL(DP) X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)

!     INTERNAL VARIABLES
!
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT, &
     &        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      REAL(DP)  DDOT,T,R
      REAL(DP)  B,C,CS,EL,EMM1,F,G,DNRM2,SCALEF,SHIFT,SL,SM,SN, &
     &                 SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
!
!
!     SET THE MAXIMUM NUMBER OF ITERATIONS.
!
      MAXIT = 30
!
!     DETERMINE WHAT IS TO BE COMPUTED.
!
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
!
!     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
!     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
!
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
!
!           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
!           PLACE THE L-TH DIAGONAL IN S(L).
!
      s(l) = DNRM2(n-l+1,x(l,l),1)
            IF (S(L) .EQ. 0.0D0) GO TO 10
               IF (X(L,L) .NE. 0.0D0) S(L) = SIGN(S(L),X(L,L))
      call DSCAL(n-l+1,1.0d0/s(l),x(l,l),1)
               X(L,L) = 1.0D0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0D0) GO TO 30
!
!              APPLY THE TRANSFORMATION.
!
      t = - DDOT(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
      call DAXPY(n-l+1,t,x(l,l),1,x(l,j),1)
   30       CONTINUE
!
!           PLACE THE L-TH ROW OF X INTO  E FOR THE
!           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
!
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
!
!           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
!           MULTIPLICATION.
!
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
!
!           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
!           L-TH SUPER-DIAGONAL IN E(L).
!
      e(l) = DNRM2(p-l,e(lp1),1)
            IF (E(L) .EQ. 0.0D0) GO TO 80
               IF (E(LP1) .NE. 0.0D0) E(L) = SIGN(E(L),E(LP1))
      call DSCAL(p-l,1.0d0/e(l),e(lp1),1)
               E(LP1) = 1.0D0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
!
!              APPLY THE TRANSFORMATION.
!
               DO 90 I = LP1, N
                  WORK(I) = 0.0D0
   90          CONTINUE
               DO 100 J = LP1, P
      call DAXPY(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          CONTINUE
               DO 110 J = LP1, P
      call DAXPY(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
!
!              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
!              BACK MULTIPLICATION.
!
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
!
!     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
!
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
!
!     IF REQUIRED, GENERATE U.
!
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0D0
  180       CONTINUE
            U(J,J) = 1.0D0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
               t = - DDOT(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
               call DAXPY(n-l+1,t,u(l,l),1,u(l,j),1)
  210          CONTINUE
  220          CONTINUE
               call DSCAL(n-l+1,-1.0d0,u(l,l),1)
               U(L,L) = 1.0D0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0D0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0D0
  260          CONTINUE
               U(L,L) = 1.0D0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
!
!     IF IT IS REQUIRED, GENERATE V.
!
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO 310 J = LP1, P
               t = - DDOT(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
               call DAXPY(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0D0
  330       CONTINUE
            V(L,L) = 1.0D0
  340    CONTINUE
  350 CONTINUE
!
!     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
!
      MM = M
      ITER = 0
  360 CONTINUE
!
!        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
!
!     ...EXIT
         IF (M .EQ. 0) GO TO 620
!
!        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
!        FLAG AND RETURN.
!
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
!     ......EXIT
            GO TO 620
  370    CONTINUE
!
!        THIS SECTION OF THE PROGRAM INSPECTS FOR
!        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
!        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
!
!           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
!           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
!           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
!                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
!           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
!
         DO 390 LL = 1, M
            L = M - LL
!        ...EXIT
            IF (L .EQ. 0) GO TO 400
            TEST = DABS(S(L)) + DABS(S(L+1))
            ZTEST = TEST + DABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0D0
!        ......EXIT
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
!           ...EXIT
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + DABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + DABS(E(LS-1))
               ZTEST = TEST + DABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0D0
!           ......EXIT
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
!
!        PERFORM THE TASK INDICATED BY KASE.
!
         GO TO (490,520,540,570), KASE
!
!        DEFLATE NEGLIGIBLE S(M).
!
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
!
!        SPLIT AT NEGLIGIBLE S(L).
!
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO 530 K = L, M
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
!
!        PERFORM ONE QR STEP.
!
  540    CONTINUE
!
!           CALCULATE THE SHIFT.
!
            SCALEF = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)), &
     &                    DABS(S(L)),DABS(E(L)))
            SM = S(M)/SCALEF
            SMM1 = S(M-1)/SCALEF
            EMM1 = E(M-1)/SCALEF
            SL = S(L)/SCALEF
            EL = E(L)/SCALEF
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550
               SHIFT = DSQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) + SHIFT
            G = SL*EL
!
!           CHASE ZEROS.
!
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N) &
     &            CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
!
!        CONVERGENCE.
!
  570    CONTINUE
!
!           MAKE THE SINGULAR VALUE  POSITIVE.
!
            IF (S(L) .GE. 0.0D0) GO TO 580
               S(L) = -S(L)
               if(wantv)call DSCAL(p,-1.0d0,v(1,l),1)
  580       CONTINUE
!
!           ORDER THE SINGULAR VALUE.
!
  590       IF (L .EQ. MM) GO TO 600
!           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
      if(wantv.and.l.lt.p)call DSWAP(p,v(1,l),1,v(1,l+1),1)
      if(wantu.and.l.lt.n)call DSWAP(n,u(1,l),1,u(1,l+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END SUBROUTINE dsvdc
