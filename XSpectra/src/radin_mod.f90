!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
Module radin_mod

!  Use constants  !WARNING, here you should put other constants!!!
  use kinds, ONLY : DP
  Implicit None

Contains

  Function para_radin(y,x,n)

    Real(dp), Intent(in) :: y(n),x(n)
    Integer,  Intent(in) :: n
    Real(dp)             :: para_radin

    ! Local

    Integer :: ierr

    Real(dp) :: ans(1),up(1)
    Real(dp), Parameter :: zero = 0.d0
    Real(dp), Allocatable, Dimension(:) :: yp,ypp

    Allocate(yp(n),ypp(n))

    Call splift(x,y,yp,ypp,n,ierr,0,zero,zero,zero,zero)

    If(ierr.Ne.1) Stop 'error calling splift from para_radin'

    up(1) = x(n)

    Call spliq(x,y,yp,ypp,n,x(1),up,1,ans,ierr)

    para_radin=ans(1)

    If(ierr.Ne.1) Stop 'error calling spliq from para_radin'

    Deallocate(yp,ypp)

  End Function para_radin

  Function apply_ke_radial(l,y,x,n)

    Real(dp), Intent(in) :: y(:),x(:)
    Integer,  Intent(in) :: l,n
    Real(dp)             :: apply_ke_radial(n)

    ! Local

    Integer :: i,ierr

    Real(dp), Parameter :: zero = 0.d0
    Real(dp), Allocatable, Dimension(:) :: yp,ypp

    Allocate(yp(n),ypp(n))

    Call splift(x,y,yp,ypp,n,ierr,0,zero,zero,zero,zero)

    If(ierr.Ne.1) Stop 'error calling splift from radin'

    apply_ke_radial(1) = 0.d0

    Do i=2,n
       apply_ke_radial(i) = -ypp(i) + Real(l*(l+1))*y(i)/x(i)**2
    End Do

    Deallocate(yp,ypp)

  end Function apply_ke_radial

  Function apply_deriv_radial(y,x,n)

    ! Assumes y is passed as r x y, and return d y / d r

    Real(dp), Intent(in) :: y(:),x(:)
    Integer,  Intent(in) :: n
    Real(dp)             :: apply_deriv_radial(n)

    ! Local

    Integer :: i,ierr

    Real(dp), Parameter :: zero = 0.d0
    Real(dp), Allocatable, Dimension(:) :: yp,ypp

    Allocate(yp(n),ypp(n))

    Call splift(x,y,yp,ypp,n,ierr,0,zero,zero,zero,zero)

    If(ierr.Ne.1) Stop 'error calling splift from para_radin'

    Do i=2,n
       apply_deriv_radial(i) = (yp(i) - y(i)/x(i))/x(i)
    End Do

    ! Project back to the origin

    apply_deriv_radial(1) = apply_deriv_radial(2) - &
         (apply_deriv_radial(3)-apply_deriv_radial(2))*x(2)/(x(3)-x(2))

    Deallocate(yp,ypp)

  end Function apply_deriv_radial


  Subroutine SPLIFT (X,Y,YP,YPP,N,IERR,ISX,A1,B1,AN,BN)


    !
    !     SANDIA MATHEMATICAL PROGRAM LIBRARY
    !     APPLIED MATHEMATICS DIVISION 2613
    !     SANDIA LABORATORIES
    !     ALBUQUERQUE, NEW MEXICO  87185
    !     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !                    ISSUED BY SANDIA LABORATORIES
    !  *                   A PRIME CONTRACTOR TO THE
    !  *                UNITED STATES DEPARTMENT OF ENERGY
    !  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
    !  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
    !  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
    !  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
    !  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
    !  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
    !  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
    !  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
    !  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
    !  * OWNED RIGHTS.
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
    !  * PART IS SAND77-1441.
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !     WRITTEN BY RONDALL E. JONES
    !
    !     ABSTRACT
    !         SPLIFT FITS AN INTERPOLATING CUBIC SPLINE TO THE N DATA POINT
    !         GIVEN IN X AND Y AND RETURNS THE FIRST AND SECOND DERIVATIVES
    !         IN YP AND YPP.  THE RESULTING SPLINE (DEFINED BY X, Y, AND
    !         YPP) AND ITS FIRST AND SECOND DERIVATIVES MAY THEN BE
    !         EVALUATED USING SPLINT.  THE SPLINE MAY BE INTEGRATED USING
    !         SPLIQ.  FOR A SMOOTHING SPLINE FIT SEE SUBROUTINE SMOO.
    !
    !     DESCRIPTION OF ARGUMENTS
    !         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
    !         E.G.   X(N), Y(N), YP(N), YPP(N), W(3N)
    !
    !       --INPUT--
    !
    !         X    - ARRAY OF ABSCISSAS OF DATA (IN INCREASING ORDER)
    !         Y    - ARRAY OF ORDINATES OF DATA
    !         N    - THE NUMBER OF DATA POINTS.  THE ARRAYS X, Y, YP, AND
    !                YPP MUST BE DIMENSIONED AT LEAST N.  (N .GE. 4)
    !         ISX  - MUST BE ZERO ON THE INITIAL CALL TO SPLIFT.
    !                IF A SPLINE IS TO BE FITTED TO A SECOND SET OF DATA
    !                THAT HAS THE SAME SET OF ABSCISSAS AS A PREVIOUS SET,
    !                AND IF THE CONTENTS OF W HAVE NOT BEEN CHANGED SINCE
    !                THAT PREVIOUS FIT WAS COMPUTED, THEN ISX MAY BE
    !                SET TO ONE FOR FASTER EXECUTION.
    !         A1,B1,AN,BN - SPECIFY THE END CONDITIONS FOR THE SPLINE WHICH
    !                ARE EXPRESSED AS CONSTRAINTS ON THE SECOND DERIVATIVE
    !                OF THE SPLINE AT THE END POINTS (SEE YPP).
    !                THE END CONDITION CONSTRAINTS ARE
    !                        YPP(1) = A1*YPP(2) + B1
    !                AND
    !                        YPP(N) = AN*YPP(N-1) + BN
    !                WHERE
    !                        ABS(A1).LT. 1.0  AND  ABS(AN).LT. 1.0.
    !
    !                THE SMOOTHEST SPLINE (I.E., LEAST INTEGRAL OF SQUARE
    !                OF SECOND DERIVATIVE) IS OBTAINED BY A1=B1=AN=BN=0.
    !                IN THIS CASE THERE IS AN INFLECTION AT X(1) AND X(N).
    !                IF THE DATA IS TO BE EXTRAPOLATED (SAY, BY USING SPLIN
    !                TO EVALUATE THE SPLINE OUTSIDE THE RANGE X(1) TO X(N))
    !                THEN TAKING A1=AN=0.5 AND B1=BN=0 MAY YIELD BETTER
    !                RESULTS.  IN THIS CASE THERE IS AN INFLECTION
    !                AT X(1) - (X(2)-X(1)) AND AT X(N) + (X(N)-X(N-1)).
    !                IN THE MORE GENERAL CASE OF A1=AN=A  AND B1=BN=0,
    !                THERE IS AN INFLECTION AT X(1) - (X(2)-X(1))*A/(1.0-A)
    !                AND AT X(N) + (X(N)-X(N-1))*A/(1.0-A).
    !
    !                A SPLINE THAT HAS A GIVEN FIRST DERIVATIVE YP1 AT X(1)
    !                AND YPN AT Y(N) MAY BE DEFINED BY USING THE
    !                FOLLOWING CONDITIONS.
    !
    !                A1=-0.5
    !
    !                B1= 3.0*((Y(2)-Y(1))/(X(2)-X(1))-YP1)/(X(2)-X(1))
    !
    !                AN=-0.5
    !
    !                BN=-3.0*((Y(N)-Y(N-1))/(X(N)-X(N-1))-YPN)/(X(N)-X(N-1)
    !
    !       --OUTPUT--
    !
    !         YP   - ARRAY OF FIRST DERIVATIVES OF SPLINE (AT THE X(I))
    !         YPP  - ARRAY OF SECOND DERIVATIVES OF SPLINE (AT THE X(I))
    !         IERR - A STATUS CODE
    !              --NORMAL CODE
    !                 1 MEANS THAT THE REQUESTED SPLINE WAS COMPUTED.
    !              --ABNORMAL CODES
    !                 2 MEANS THAT N, THE NUMBER OF POINTS, WAS .LT. 4.
    !                 3 MEANS THE ABSCISSAS WERE NOT STRICTLY INCREASING.

    !
    !       --WORK--
    !
    !         W    - ARRAY OF WORKING STORAGE DIMENSIONED AT LEAST 3N.
    !

    ! INPUT

    Integer,  Intent(in) :: N,ISX
    Real(dp), Intent(in) :: X(N),Y(N),A1,B1,AN,BN

    ! OUTPUT

    Integer,  Intent(out) :: IERR
    Real(dp), Intent(out) :: YP(N),YPP(N)

    ! WORK

    Real(dp), Allocatable, Dimension(:,:) :: W

    ! LOCAL

    Real(dp), Parameter :: FOUR=4.D0

    Integer :: nm1,nm2,i,j
    Real(dp) :: dold,dnew


    If (N.Lt.4) Then
       IERR = 2
       Return
    Endif
    NM1  = N-1
    NM2  = N-2
    If (ISX.Gt.0) GO TO 40
    Do I=2,N
       If (X(I)-X(I-1) .Le. 0) Then
          IERR = 3
          Return
       Endif
    End Do

    ! Allocate the work space

    Allocate(W(N,3))

    !
    !     DEFINE THE TRIDIAGONAL MATRIX
    !
    W(1,3) = X(2)-X(1)
    Do I=2,NM1
       W(I,2) = W(I-1,3)
       W(I,3) = X(I+1)-X(I)
       W(I,1) = 2*(W(I,2)+W(I,3))
    End Do
    W(1,1) = FOUR
    W(1,3) =-4*A1
    W(N,1) = FOUR
    W(N,2) =-4*AN
    !
    !     L U DECOMPOSITION
    !
    Do  I=2,N
       W(I-1,3) = W(I-1,3)/W(I-1,1)
       W(I,1) = W(I,1) - W(I,2)*W(I-1,3)
    End Do
    !
    !     DEFINE *CONSTANT* VECTOR
    !
40  YPP(1) = 4*B1
    DOLD = (Y(2)-Y(1))/W(2,2)
    Do  I=2,NM2
       DNEW   = (Y(I+1) - Y(I))/W(I+1,2)
       YPP(I) = 6*(DNEW - DOLD)
       YP(I)  = DOLD
       DOLD = DNEW
    End Do
    DNEW = (Y(N)-Y(N-1))/(X(N)-X(N-1))
    YPP(NM1) = 6*(DNEW - DOLD)
    YPP(N) = 4*BN
    YP(NM1)= DOLD
    YP(N) = DNEW
    !
    !     FORWARD SUBSTITUTION
    !
    YPP(1) = YPP(1)/W(1,1)
    Do  I=2,N
       YPP(I) = (YPP(I) - W(I,2)*YPP(I-1))/W(I,1)
    End Do
    !
    !     BACKWARD SUBSTITUTION
    !
    Do  J=1,NM1
       I = N-J
       YPP(I) = YPP(I) - W(I,3)*YPP(I+1)
    End Do
    !
    !     COMPUTE FIRST DERIVATIVES
    !
    YP(1) = (Y(2)-Y(1))/(X(2)-X(1)) - (X(2)-X(1))*(2*YPP(1) + YPP(2))/6
    Do  I=2,NM1
       YP(I) = YP(I) + W(I,2)*(YPP(I-1) + 2*YPP(I))/6
    End Do
    YP(N) = YP(N) + (X(N)-X(NM1))*(YPP(NM1) + 2*YPP(N))/6

    IERR = 1

    Deallocate(W)

    Return
  End Subroutine SPLIFT

  Subroutine SPLIQ(X,Y,YP,YPP,N,XLO,XUP,NUP,ANS,IERR)

    ! INPUT

    Integer, Intent(in) :: N,NUP
    Real(dp), Intent(in) :: X(N),Y(N),YP(N),YPP(N),XLO,XUP(NUP)

    ! OUTPUT

    Integer, Intent(out) :: IERR
    Real(dp), Intent(out) :: ANS(NUP)

    !     
    !     SANDIA MATHEMATICAL PROGRAM LIBRARY
    !     APPLIED MATHEMATICS DIVISION 2613
    !     SANDIA LABORATORIES
    !     ALBUQUERQUE, NEW MEXICO  87185
    !     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !                    ISSUED BY SANDIA LABORATORIES
    !  *                   A PRIME CONTRACTOR TO THE
    !  *                UNITED STATES DEPARTMENT OF ENERGY
    !  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
    !  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
    !  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
    !  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
    !  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
    !  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
    !  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
    !  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
    !  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
    !  * OWNED RIGHTS.
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
    !  * PART IS SAND77-1441.
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    !     THIS ROUTINE WAS WRITTEN BY M. K. GORDON
    !
    !     ABSTRACT
    !
    !     SUBROUTINE SPLIQ INTEGRATES A CUBIC SPLINE (GENERATED BY
    !     SPLIFT, SMOO, ETC.) ON THE INTERVALS (XLO,XUP(I)), WHERE XUP
    !     IS A SEQUENCE OF UPPER LIMITS ON THE INTERVALS OF INTEGRATION.
    !     THE ONLY RESTRICTIONS ON XLO AND XUP(*) ARE
    !                XLO .LT. XUP(1),
    !                XUP(I) .LE. XUP(I+1)   FOR EACH I .
    !     ENDPOINTS BEYOND THE SPAN OF ABSCISSAS ARE ALLOWED.
    !     THE SPLINE OVER THE INTERVAL (X(I),X(I+1)) IS REGARDED
    !     AS A CUBIC POLYNOMIAL EXPANDED ABOUT X(I) AND IS INTEGRATED
    !     ANALYTICALLY.
    !
    !     DESCRIPTION OF ARGUMENTS
    !         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
    !         E.G.  X(N), Y(N), YP(N), YPP(N), XUP(NUP), ANS(NUP)
    !
    !      --INPUT--
    !
    !        X    - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
    !               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
    !        Y    - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
    !               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
    !        YP   - ARRAY OF FIRST DERIVATIVES OF THE SPLINE AT ABSCISSAS.
    !               USUALLY YP IS THE SAME AS YP IN SPLIFT OR R1 IN SMOO.
    !        YPP  - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
    !               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
    !        N    - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
    !        XLO  - LEFT ENDPOINT OF INTEGRATION INTERVALS.
    !        XUP  - RIGHT ENDPOINT OR ARRAY OF RIGHT ENDPOINTS OF
    !               INTEGRATION INTERVALS IN ASCENDING ORDER.
    !        NUP  - THE NUMBER OF RIGHT ENDPOINTS.  IF NUP IS GREATER THAN
    !               1, THEN XUP AND ANS MUST BE DIMENSIONED AT LEAST NUP.
    !
    !      --OUTPUT--
    !
    !        ANS -- ARRAY OF INTEGRAL VALUES, THAT IS,
    !               ANS(I) = INTEGRAL FROM XLO TO XUP(I)
    !        IERR -- ERROR STATUS
    !                = 1 INTEGRATION SUCCESSFUL
    !                = 2 IMPROPER INPUT - N.LT.4 OR NUP.LT.1
    !                = 3 IMPROPER INPUT - ABSCISSAS NOT IN
    !                        STRICTLY ASCENDING ORDER
    !                = 4 IMPROPER INPUT - RIGHT ENDPOINTS XUP NOT
    !                        IN ASCENDING ORDER
    !                = 5 IMPROPER INPUT - XLO.GT.XUP(1)
    !                = 6 INTEGRATION SUCCESSFUL BUT AT LEAST ONE ENDPOINT
    !                        NOT WITHIN SPAN OF ABSCISSAS
    !              ** NOTE.  ERRCHK PROCESSES DIAGNOSTICS FOR CODES 2,3,4,5

    ! *** Local

    Integer :: nm1,nm2,i,j,m
    Real(dp) :: hlo,hlo2,hi,hi2,hi3,hup,hup2,hup3,hup4,hsum,hdiff
    Real(dp) :: sum,sum0,sum1,sum2,sum3,psum0,psum1,psum2,psum3

    !
    !   CHECK FOR IMPROPER INPUT
    !
    IERR = 2
    If(N .Lt. 4  .Or.  NUP .Lt. 1) Then 
       Return
    Endif
    NM1 = N-1
    NM2 = N-2
    IERR = 3
    Do  I = 1,NM1
       If(X(I) .Ge. X(I+1)) Then
          Return
       Endif

    End Do
    If(NUP .Ne. 1) Then
       IERR = 4
       Do  I = 2,NUP
          If(XUP(I-1) .Gt. XUP(I)) Then
             Return
          Endif
       End Do
    Endif
    IERR = 5
    If(XLO .Gt. XUP(1)) Then
       Return
    Endif
    IERR = 1
    If(XLO .Lt. X(1)  .Or.  XUP(NUP) .Gt. X(N)) IERR = 6
    !     
    !     LOCATE XLO IN INTERVAL (X(I),X(I+1))
    !     
    Do  I = 1,NM2
       If(XLO .Lt. X(I+1)) GO TO 20
    End Do
    I = NM1
20  HLO = XLO-X(I)
    HLO2 = HLO*HLO
    HI = X(I+1)-X(I)
    HI2 = HI*HI
    Do  J = 1,NUP
       If(XUP(J) .Gt. X(I+1)  .And.  XLO .Lt. X(NM1)) GO TO 40
       !     
       !     COMPUTE SPECIAL CASES OF XUP IN INTERVAL WITH XLO
       !     
       HUP = XUP(J)-X(I)
       HSUM = HUP+HLO
       HDIFF = HUP-HLO
       HUP2 = HUP*HUP
       SUM = (YPP(I+1)-YPP(I))*HSUM*HDIFF*(HUP2+HLO2)/(24*HI)
       SUM = SUM + YPP(I)*HDIFF*(HUP2+HLO*HUP+HLO2)/6
       SUM = SUM + YP(I)*HDIFF*HSUM/2
       SUM = SUM + Y(I)*HDIFF
       ANS(J) = SUM
    End Do
    Return
    !     
    !     COMPUTE INTEGRAL BETWEEN XLO AND X(I+1) AS FOUR TERMS IN TAYLOR
    !     POLYNOMIAL AND ADVANCE I TO I+1
    !     
40  HDIFF = HI-HLO
    HSUM = HI+HLO
    SUM0 = Y(I)*HDIFF
    SUM1 = YP(I)*HDIFF*HSUM
    SUM2 = YPP(I)*HDIFF*(HI2+HI*HLO+HLO2)
    SUM3 = (YPP(I+1)-YPP(I))*HDIFF*HSUM*(HI2+HLO2)/HI
    I = I+1
    !     
    !     LOCATE EACH XUP(M) IN INTERVAL (X(I),X(I+1))
    !     
    Do  M = J,NUP
50     If(XUP(M) .Lt. X(I+1)  .Or.  I .Eq. NM1) GO TO 60
       !     
       !     AUGMENT INTEGRAL BETWEEN ABSCISSAS TO INCLUDE INTERVAL
       !     (X(I),X(I+1)) AND ADVANCE I TO I+1
       !     
       HI = X(I+1)-X(I)
       HI2 = HI*HI
       HI3 = HI2*HI
       SUM0 = SUM0 + Y(I)*HI
       SUM1 = SUM1 + YP(I)*HI2
       SUM2 = SUM2 + YPP(I)*HI3
       SUM3 = SUM3 + (YPP(I+1)-YPP(I))*HI3
       I = I+1
       GO TO 50
       !     
       !     INTEGRAL BETWEEN X(I) AND XUP(M) IS ZERO
       !     
60     If(XUP(M) .Ne. X(I)) Then
          !     
          !     COMPUTE INTEGRAL BETWEEN X(I) AND XUP(M) AND EVALUATE
          !     TAYLOR POLYNOMIAL IN REVERSE ORDER
          !     
          HUP = XUP(M)-X(I)
          HUP2 = HUP*HUP
          HUP3 = HUP2*HUP
          HUP4 = HUP3*HUP
          HI = X(I+1)-X(I)
          PSUM0 = Y(I)*HUP
          PSUM1 = YP(I)*HUP2
          PSUM2 = YPP(I)*HUP3
          PSUM3 = (YPP(I+1)-YPP(I))*HUP4/HI
          SUM = (SUM3+PSUM3)/24 + (SUM2+PSUM2)/6
          SUM = SUM + (SUM1+PSUM1)/2
          SUM = SUM + (SUM0+PSUM0)
       Else
          SUM = ((SUM3/24 + SUM2/6) + SUM1/2) + SUM0
       Endif
       ANS(M) = SUM
    End Do
    Return
  End Subroutine SPLIQ


  SUBROUTINE SPLINT (X,Y,YPP,N,XI,YI,YPI,YPPI,NI,KERR)
         
         ! INPUT
    
         Integer, Intent(in) :: N,NI
         Real(dp), Intent(in) :: X(N),Y(N),YPP(N),XI(NI)

         ! OUTPUT

         Integer, Intent(out) :: KERR
         Real(dp), Intent(out) :: YI(NI),YPI(NI),YPPI(NI)



!
 !C     SANDIA MATHEMATICAL PROGRAM LIBRARY
 !C     APPLIED MATHEMATICS DIVISION 2613
 !C     SANDIA LABORATORIES
 !C     ALBUQUERQUE, NEW MEXICO  87185
 !C     CONTROL DATA 6600/7600  VERSION 7.2  MAY 1978
 !C  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 !C                    ISSUED BY SANDIA LABORATORIES
 !C  *                   A PRIME CONTRACTOR TO THE
 !C  *                UNITED STATES DEPARTMENT OF ENERGY
 !C  * * * * * * * * * * * * * * * NOTICE  * * * * * * * * * * * * * * *
 !C  * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
 !C  * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
 !C  * UNITED STATES DEPARTMENT OF ENERGY NOR ANY OF THEIR EMPLOYEES,
 !C  * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR EMPLOYEES
 !C  * MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL
 !!  * LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS OR
!  * USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
!  * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
!  * OWNED RIGHTS.
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!  * THE PRIMARY DOCUMENT FOR THE LIBRARY OF WHICH THIS ROUTINE IS
!  * PART IS SAND77-1441.
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     WRITTEN BY RONDALL E. JONES
!
!     ABSTRACT
!
!         SPLINT EVALUATES A CUBIC SPLINE AND ITS FIRST AND SECOND
!         DERIVATIVES AT THE ABSCISSAS IN XI.  THE SPLINE (WHICH
!         IS DEFINED BY X, Y, AND YPP) MAY HAVE BEEN DETERMINED BY
!         SPLIFT OR SMOO OR ANY OTHER SPLINE FITTING ROUTINE THAT
!         PROVIDES SECOND DERIVATIVES.
!
!     DESCRIPTION OF ARGUMENTS
!         THE USER MUST DIMENSION ALL ARRAYS APPEARING IN THE CALL LIST
!         E.G.  X(N), Y(N), YPP(N), XI(NI), YI(NI), YPI(NI), YPPI(NI)
!
!       --INPUT--
!
!         X   - ARRAY OF ABSCISSAS (IN INCREASING ORDER) THAT DEFINE TH
!               SPLINE.  USUALLY X IS THE SAME AS X IN SPLIFT OR SMOO.
!         Y   - ARRAY OF ORDINATES THAT DEFINE THE SPLINE.  USUALLY Y I
!               THE SAME AS Y IN SPLIFT OR AS R IN SMOO.
!         YPP - ARRAY OF SECOND DERIVATIVES THAT DEFINE THE SPLINE.
!               USUALLY YPP IS THE SAME AS YPP IN SPLIFT OR R2 IN SMOO.
!         N   - THE NUMBER OF DATA POINTS THAT DEFINE THE SPLINE.
!               THE ARRAYS X, Y, AND YPP MUST BE DIMENSIONED AT LEAST N
!               N MUST BE GREATER THAN OR EQUAL TO 2.
!         XI  - THE ABSCISSA OR ARRAY OF ABSCISSAS (IN ARBITRARY ORDER)
!               AT WHICH THE SPLINE IS TO BE EVALUATED.
!               EACH XI(K) THAT LIES BETWEEN X(1) AND X(N) IS A CASE OF
!               INTERPOLATION.  EACH XI(K) THAT DOES NOT LIE BETWEEN
!               X(1) AND X(N) IS A CASE OF EXTRAPOLATION.  BOTH CASES
!               ARE ALLOWED.  SEE DESCRIPTION OF KERR.
!         NI  - THE NUMBER OF ABSCISSAS AT WHICH THE SPLINE IS TO BE
!               EVALUATED.  IF NI IS GREATER THAN 1, THEN XI, YI, YPI,
!               AND YPPI MUST BE ARRAYS DIMENSIONED AT LEAST NI.
!               NI MUST BE GREATER THAN OR EQUAL TO 1.
!
!       --OUTPUT--
!
!         YI  - ARRAY OF VALUES OF THE SPLINE (ORDINATES) AT XI.
!         YPI - ARRAY OF VALUES OF THE FIRST DERIVATIVE OF SPLINE AT XI

!         YPPI- ARRAY OF VALUES OF SECOND DERIVATIVES OF SPLINE AT XI.
!         KERR- A STATUS CODE
!             --NORMAL CODES
!                1 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
!                  IN XI USING ONLY INTERPOLATION.
!                2 MEANS THAT THE SPLINE WAS EVALUATED AT EACH ABSCISSA
!                  IN XI, BUT AT LEAST ONE EXTRAPOLATION WAS PERFORMED.
!             -- ABNORMAL CODE
!                3 MEANS THAT THE REQUESTED NUMBER OF EVALUATIONS, NI,
!                 WAS NOT POSITIVE.
!


   ! *** Local

    Integer :: nm1,i,K,IL,IR
     Real(dp) :: H,H2,XR,XR2,XR3,XL,XL2,XL3,XX
  
!
!    CHECK INPUT
!
    IF (NI) 1,1,2
1   CONTINUE
    !    1 CALL ERRCHK(67,67HIN SPLINT,  THE REQUESTED NUMBER OF INTERPOLATI
    !     1NS WAS NOT POSITIVE)
    KERR = 3
    RETURN
2   KERR = 1
    NM1= N-1

    !
    !    K IS INDEX ON VALUE OF XI BEING WORKED ON.  XX IS THAT VALUE.
    !    I IS CURRENT INDEX INTO X ARRAY.
    !
    K  = 1
    XX = XI(1)
    IF (XX.LT.X(1)) GO TO 90
    IF (XX.GT.X(N)) GO TO 80
    IL = 1
    IR = N
    !
    !   BISECTION SEARCH
    !
10  I  = (IL+IR)/2
    IF (I.EQ.IL) GO TO 100
    IF (XX-X(I)) 20,100,30
20  IR = I
    GO TO 10
30  IL = I
    GO TO 10
    !
    !     LINEAR FORWARD SEARCH
    !
50  IF (XX-X(I+1)) 100,100,60
60  IF (I.GE.NM1) GO TO 80
    I  = I+1
    GO TO 50
    !
    ! EXTRAPOLATION
    !
80  KERR = 2
    I  = NM1
    GO TO 100
90  KERR = 2
    I  = 1
    !
    !     INTERPOLATION
    !
100 H  = X(I+1) - X(I)
    H2 = H*H
    XR = (X(I+1)-XX)/H
    XR2= XR*XR
    XR3= XR*XR2
    XL = (XX-X(I))/H
    XL2= XL*XL
    XL3= XL*XL2
    YI(K) = Y(I)*XR + Y(I+1)*XL-H2*(YPP(I)*(XR-XR3) + YPP(I+1)*(XL-XL3))/6.0D0
    YPI(K) = (Y(I+1)-Y(I))/H+H*(YPP(I)*(1.0D0-3.0D0*XR2)-YPP(I+1)*&
         (1.0D0-3.0D0*XL2))/6.0D0
    YPPI(K) = YPP(I)*XR + YPP(I+1)*XL
    !
    !     NEXT POINT
    !
    IF (K.GE.NI) RETURN
    K = K+1
    XX = XI(K)
    IF (XX.LT.X(1)) GO TO 90
    IF (XX.GT.X(N)) GO TO 80
    IF (XX-XI(K-1)) 110,100,50
110 IL = 1
    IR = I+1
    GO TO 10
    !
  END subroutine SPLINT
     
End Module radin_mod
