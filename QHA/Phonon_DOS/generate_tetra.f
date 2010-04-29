!
! Copyright (C) 2006 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! GNU License
!
! Eyvaz Isaev

! Theoretical Physics Department,
! Moscow State Institute of Steel and Alloys
! (Technological University)
!
! Condensed Matter Theory Group,
! Uppsala University, Sweden
!
! Eyvaz.Isaev@fysik.uu.se, eyvaz_isaev@yahoo.com
!
! Adopted from a program published in a preprint (early 90th) of Lebedev Physical Institute (Moscow)
! Early versions of this program was used in ab initio psedopotentials code of E.I.Isaev
!

	subroutine generate_tetra(npt)
        include 'parameters.h'
        PARAMETER (KMAX=15000)
C
	dimension lstpnt(100000),iarr(4),pnt1(4,4)
	integer ttrin(4,27000)

        EQUIVALENCE(IARR(1),I),(IARR(2),J),(IARR(3),K),(IARR(4),L)
C
        COMMON /INTET/ NIN,TTRIN
        omg=0.D0
        NPNT=0
	nt0=ntet0
	n1=ndiv
        NTETMX=N1**3
	write(6,991)
991	format(22('*'),' generate_tetra ',21(('*')))
998   FORMAT(4I3)
      PRINT *,' NT0=',NT0,' NTETMX=',NTETMX
        GOTO 11
12      DO 58 KTET0=1,NT0
        LKTET=NTETMX*(KTET0-1)
      DO 141 I=1,4
      IPNT=TTR0(I,KTET0)
      DO 140 K=1,3
      PNT1(K,I)=PNT0(K,IPNT)
140     CONTINUE
      PNT1(4,I)=1.
141     CONTINUE
      PRINT 100,((PNT1(II,JJ),JJ=1,4),II=1,4)
100   FORMAT(1X,4F9.4)
       ot=dabs(det4(pnt1))/6.

	omg=omg+ot
        WRITE(6,99)  ot
99	format('  volume of tetrahedron =',f9.5)
        NUM=0
        DO 21 IJK=0,N1
        L=N1-IJK
        DO 21 JK=0,IJK
        I=IJK-JK
        DO 21 K=0,JK
        J=JK-K
        NUM=NUM+1
        NPNT=NPNT+1
C       WRITE(6,*)'I,J,K,L,PNT',I,J,K,L,NPNT
        IF(NPNT.GT.KMAX) then
	print*,'NPNT==',npnt,'  KMAX=',kmax
	STOP'MAXIMUM NUMBER OF K-POINTS EXCEEDED'
        endif
	DO 23 M1=1,3
        PNT(M1,NPNT)=0.
        DO 23 M2=1,4
23      PNT(M1,NPNT)=PNT(M1,NPNT)+PNT0(M1,TTR0(M2,KTET0))*IARR(M2)/N1
        IF(I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0)GOTO 22
        DO 1 IND=1,NPNT-1
        IF (NPNT.EQ.1)GOTO 1
        S=0.
        DO 2 KCHK =1,3
2       S=S+dABS(PNT(KCHK,IND)-PNT(KCHK,NPNT))
        IF(S-.0001)3,3,1
1       CONTINUE
        IND=0
3       IF(IND.EQ.0) GOTO 22
        NPNT=NPNT-1
        LSTPNT(NUM)=IND
        GOTO 21
22      LSTPNT(NUM)=NPNT
21      CONTINUE
20      DO 10 IT =1,NTETMX
        LIT=LKTET+IT
        DO 10 K=1,4
          I25=TTRIN(K,IT)
      TTR(K,LIT)=LSTPNT(I25)
!!10	print*,'ttr in trmsh7',lit,ttr(k,lit)
10	continue
        NTET(KTET0)=N1**3
58      CONTINUE
	npt=npnt
        omg48=1./omg*2.
	write(6,98) omg,omg48
c	write(10,98) omg,omg48
98	format('  total volume of BZ is =',f9.4,'  omg48=',f9.4)
	write(6,992) 
c	write(10,992)
992	format(18('*'),' end of generate_tetra ',18('*'))
       DO 777 JJ=1,NPNT
c       PRINT*,  (PNT(II,JJ),II=1,3)
C       PRINT*, ((PNT(II,JJ),II=1,3),JJ=1,NPNT)
777   CONTINUE
         RETURN
C       GET TTRIN
11      NIN=0
        DO 15 IJK=0,N1
        DO 15 JK=0,IJK
        I=IJK-JK
        DO 15 K=0,JK
        J=JK-K
        IF(IJK.EQ.N1)GOTO 15
        IF(IJK.EQ.N1-1)GOTO 16
        IF(IJK.EQ.N1-2)GOTO 17
        CALL TTRGEN(I,J,K,6)
17      DO 18 L=2,5
18      CALL TTRGEN(I,J,K,L)
16      CALL TTRGEN(I,J,K,1)
15      CONTINUE
        GOTO 12
C**************************************
C       GET LSTPNT
      END
        SUBROUTINE TTRGEN(I1,J1,K1,L)
        INTEGER TTRIN(4,27000),IND(3,4,6)
        COMMON/INTET/NIN,TTRIN
        DATA IND / 0,0,0, 1,0,0, 0,1,0, 0,0,1,
     1             1,0,0, 0,1,0, 0,0,1, 1,0,1,
     2             1,0,0, 0,1,0, 1,1,0, 1,0,1,
     3             0,0,1, 0,1,0, 0,1,1, 1,0,1,
     4             0,1,1, 0,1,0, 1,1,0, 1,0,1,
     5             1,1,0, 0,1,1, 1,0,1, 1,1,1/
        NUM(I,J,K)=I*(I+1)*(I+2)/6+J*(J+1)/2+K+1
        NIN=NIN+1
        DO 1 M=1,4
        K=IND(3,M,L)+K1
        JK=K+IND(2,M,L)+J1
        IJK=JK+IND(1,M,L)+I1
1       TTRIN(M,NIN)=NUM(IJK,JK,K)
        RETURN
        END 
