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

       SUBROUTINE Tetrahedra(E0,ns,IV)
       include 'parameters.h'

       DIMENSION IND(4),L(4),e0(4),et(4)
       LOGICAL GF
       GF(Xx)=EF.GT.Xx
       DF(Xx)=EF-Xx

	lv=ns
	if(iv.eq.0) lv=1
       DO 10 I=1,4
 10    IND(I)=1
       DO 1 I=1,3
       I1=I+1
       DO 1 J=I1 ,4
       IF(E0(I).GT.E0(J))GOTO 2
       IND(I)=IND(I)+1
       GOTO 1
 2     IND(J)=IND(J)+1
 1     CONTINUE
       DO 20 I=1,4
       JJ=IND(I)
       L(JJ)=I
 20   ET(JJ)=E0((I))
c	print*,et
       IF(GF(ET(4))) GOTO 4
       TS=0.
       DS=0.
       DO 24 K=1,LV
       ADS(K)=0.
 24    ATS(K)=0.
       RETURN
 4     IF(GF(ET(3))) GOTO 5
       EF4=DF(ET(4))
       E24=EF4/(ET(2)-ET(4))
       E14=EF4/(ET(1)-ET(4))
       E34=EF4/(ET(3)-ET(4))
       TMP=E24*E14*E34*OT
       DS=TMP*3.0/EF4
       TS=TMP
        IF(IV.EQ.1)THEN
       DO 25 K=1,LV
       A14=E14*( A0(L(1),K)- A0(L(4),K))
       A34=E34*(A0(L(3),K)-A0(L(4),K))
       A24=E24*( A0(L(2),K)- A0(L(4),K))
       ADS(K)=DS*(A0(L(4),K)+(A24+A14+A34)*.33333333)
 25    ATS(K)=TS*(A0(L(4),K)+(A24+A14+A34)*.25)
        ENDIF
       RETURN
 5     IF(GF(ET(2)))GOTO6
C SECTION IS A QUADRANGLE
        OM23=ET(2)-ET(3)
C DIVIDING INTO TWO TETRAHEDRA
C THE FIRST TETRAHEDRON
      DOM1=DF(ET(4))
      DOM11=ET(1)-ET(4)
      DOM12=ET(2)-ET(4)
      DOM13=DOM1
      F1=1./DOM11/DOM12*(ET(2)-EF)/OM23
C THE SECOND TETRAHEDRON
      DOM2=-DF(ET(1))
      DOM21=DOM2
      DOM22=ET(1)-ET(3)
      DOM23=ET(1)-ET(4)
      F2=1./DOM22/DOM23
      G2=(EF-ET(3))/OM23
C
        DS=(DOM1*F1+DOM2*F2*G2)*3.*OT
        TS=(DOM1**2*F1+(1.-DOM2**2*F2)*G2)*OT
        IF(IV.EQ.1)THEN
       DO 26 K=1,LV
        PP1=A0(L(4),K)
        PP2=A0(L(1),K)
        PMDL=A0(L(3),K)+(A0(L(2),K)-A0(L(3),K))*DF(ET(3))/(ET(2)-ET(3))
      PI11=(PP2-PP1)/DOM11
      PI12=(A0(L(2),K)-PP1)/DOM12
      PI13=(PMDL-PP1)/DOM13
        PIDS1=(PI11+PI12+PI13)*.333333333
        PITS1=(PI11+PI12+PI13)*.25
      PI21=-(PP2-PMDL)/DOM21
      PI22=-(PP2-A0(L(3),K))/DOM22
      PI23=-(PP2-PP1)/DOM23
        PIDS2=(PI21+PI22+PI23)*.333333333
        PITS2=(PI21+PI22+PI23)*.25
        ADS(K)=(DOM1*F1*(PP1+DOM1*PIDS1)
     1      +DOM2*F2*G2*(PP2+DOM2*PIDS2))*3.*OT
        ATS(K)=(DOM1**2*F1*(PP1+DOM1*PITS1)
     1      +((PMDL+PP1+PP2+A0(L(3),K))*.25
     2      -DOM2**2*F2*(PP2+DOM2*PITS2))*G2)*OT
 26    CONTINUE
        ENDIF
        RETURN
 6     IF(GF(ET(1)))GOTO7
       EF1=DF(ET(1))
       E14=EF1/(ET(1)-ET(4))
       E12=EF1/(ET(1)-ET(2))
       E13=EF1/(ET(1)-ET(3))
       TMP=E12*E13*E14*OT
C       !NB TMP<0,EF1<0
       DS=TMP*3./EF1
        TS=OT+TMP
        IF(IV.EQ.1)THEN
       DO 27 K=1,LV
       A14=E14*( A0(L(1),K)- A0(L(4),K))
       A13=E13*(A0(L(1),K)-A0(L(3),K))
       A12=E12*( A0(L(1),K)- A0(L(2),K))
       ADS(K)=DS*(A0(L(1),K)+(A12+A14+A13)*.33333333)
 27   ATS(K)=(A0(1,K)+A0(2,K)+A0(3,K)+A0(4,K))*.25*OT
     1      +TMP*(A0(L(1),K)+(A12+A14+A13)*.25)
        ENDIF
        RETURN
 7    TS=OT
       DS=0.
       DO 28 K=1,LV
       ADS(K)=0.
 28   ATS(K)=(A0(1,K)+A0(2,K)+A0(3,K)+A0(4,K))*.25*OT
       RETURN
       END 
