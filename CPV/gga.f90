!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      module gradient_corrections

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

! ... Gradient Correction & exchange and correlation

        logical :: tgc
        integer :: mfxcx,mfxcc,mgcx,mgcc
        REAL(dbl) :: salpha,bbeta,betapp

        ! PUBLIC :: gc_setup

      contains

        subroutine gc_setup(tgc_inp,type_inp)
          logical, intent(in)      :: tgc_inp
          character*20, intent(in) :: type_inp

          tgc = tgc_inp
 
          IF(tgc) THEN
            IF(type_inp.EQ.'BLYP') THEN
              mfxcx = 1
              mgcx  = 1
              mfxcc = 3
              mgcc  = 2
              bbeta = 0.0042d0
              betapp= 0.0d0
              salpha= 0.6666666666d0
            ELSE IF(type_inp.EQ.'BP') THEN
              mfxcx = 1
              mgcx  = 1
              mfxcc = 1
              mgcc  = 1
              bbeta = 0.0042d0
              betapp= 0.0d0
              salpha= 0.6666666666d0
            ELSE
              call errore(' module gradient corrections ', &
                         ' tgc is true but type is wrong ',0)
            ENDIF
          ENDIF
          return
        end subroutine gc_setup

!     ==================================================================
      SUBROUTINE XC(RHO,EX,EC,VX,VC)
!     ==--------------------------------------------------------------==
!     ==  LDA EXCHANGE AND CORRELATION FUNCTIONALS                    ==
!     ==                                                              ==
!     ==  EXCHANGE  :  SLATER alpha                                   ==
!     ==  CORRELATION : CEPERLEY & ALDER (PERDEW-ZUNGER PARAMETERS)   ==
!     ==                VOSKO, WILK & NUSSAIR                         ==
!     ==                LEE, YANG & PARR                              ==
!     ==                PERDEW & WANG                                 ==
!     ==                WIGNER                                        ==
!     ==                HEDIN & LUNDQVIST                             ==
!     ==                ORTIZ & BALLONE (PERDEW-ZUNGER FORMULA)       ==
!     ==                ORTIZ & BALLONE (PERDEW-WANG FORMULA)         ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)

! ... Gradient Correction & exchange and correlation
      REAL(dbl) :: SMALL, PI34, THIRD
      PARAMETER (SMALL=1.D-10)
      PARAMETER (PI34= 0.75D0 / 3.141592653589793D+00 ,THIRD=1.D0/3.D0)
!     ==--------------------------------------------------------------==
!..Exchange
      IF(MFXCX.EQ.1) THEN
        CALL SLATERX(RHO,EX,VX,SALPHA)
      ELSE
        EX=0.0D0
        VX=0.0D0
      ENDIF
      IF(RHO.LE.SMALL) THEN
        EC = 0.0D0
        VC = 0.0D0
        EX = 0.0D0
        VX = 0.0D0
      ELSE IF(MFXCC.EQ.1) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL PZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.2) THEN
        RS = (PI34/RHO)**THIRD
        CALL VWN(RS,EC,VC)
      ELSEIF(MFXCC.EQ.3) THEN
        CALL LYP(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.4) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL PW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.5) THEN
        CALL WIGNER(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.6) THEN
        CALL HEDIN(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.7) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL OBPZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.8) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL OBPW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.9) THEN
        RS=(PI34/RHO)**THIRD
        CALL PADE(RS,EC,VC)
      ELSE
        EC=0.0D0
        VC=0.0D0
      ENDIF
      RETURN
      END SUBROUTINE XC
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE GCXC(RHO,GRHO,SX,SC,V1X,V2X,V1C,V2C)
!     ==--------------------------------------------------------------==
!     ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
!     ==                                                              ==
!     ==  EXCHANGE  :  BECKE88                                        ==
!     ==               GGAX                                           ==
!     ==               PBEX                                           ==
!     ==  CORRELATION : PERDEW86                                      ==
!     ==                LEE, YANG & PARR                              ==
!     ==                GGAC                                          ==
!     ==                PBEC                                          ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: SMALL
      PARAMETER (SMALL=1.D-10)

! ... Gradient Correction & exchange and correlation

!     ==--------------------------------------------------------------==
!..Exchange
      IF(RHO.LE.SMALL) THEN
        SX  = 0.0D0
        V1X = 0.0D0
        V2X = 0.0D0
      ELSEIF(MGCX.EQ.1) THEN
        CALL BECKE88(bbeta,RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.2) THEN
        CALL GGAX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.3) THEN
        CALL PBEX(RHO,GRHO,SX,V1X,V2X)
      ELSE
        SX=0.0D0
        V1X=0.0D0
        V2X=0.0D0
      ENDIF
!..Correlation
      IF(RHO.LE.SMALL) THEN
        SC  = 0.0D0
        V1C = 0.0D0
        V2C = 0.0D0
      ELSEIF(MGCC.EQ.1) THEN
        CALL PERDEW86(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.2) THEN
        CALL GLYP(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.3) THEN
        CALL GGAC(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.4) THEN
        CALL PBEC(RHO,GRHO,SC,V1C,V2C)
      ELSE
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ENDIF
      RETURN
      END SUBROUTINE GCXC
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE SLATERX(RHO,EX,VX,ALPHA)
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: SMALL, F1, THIRD, F43
      PARAMETER (SMALL=1.D-10)
      PARAMETER (F1 = -1.10783814957303361D0)
      PARAMETER (THIRD=1.D0/3.D0,F43=4.D0/3.D0)
!     ==--------------------------------------------------------------==
      IF(RHO.LE.SMALL) THEN
        EX = 0.0D0
        VX = 0.0D0
      ELSE
        RS = RHO**THIRD
        EX = F1*ALPHA*RS
        VX = F43*F1*ALPHA*RS
      ENDIF
      RETURN
      END SUBROUTINE SLATERX
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE PZ(RS,EPZ,VPZ,IFLG)
!     ==--------------------------------------------------------------==
!     ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: A, B, C, D, GC, B1, B2
      PARAMETER (A=0.0311D0,B=-0.048D0,C=0.0020D0,D=-0.0116D0, &
                GC=-0.1423D0,B1=1.0529D0,B2=0.3334D0)
!     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
!..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+ (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
!..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7./6.*B1*RS1+4./3.*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
      RETURN
      END SUBROUTINE PZ
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE VWN(RS,EVWN,VVWN)
!     ==--------------------------------------------------------------==
!     ==  S.H VOSKO, L.WILK, AND M. NUSAIR,                           ==
!     ==                 CAN. J. PHYS. 58 1200  (1980)                ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: a, b, c, x0, two
      PARAMETER (A=0.0310907,B=3.72744,C=12.9352,X0=-0.10498)
      PARAMETER (TWO=2.0D0)
!     ==--------------------------------------------------------------==
      Q  = SQRT(4.D0*C-B*B)
      F1 = TWO*B/Q
      F2 = B*X0/(X0*X0+B*X0+C)
      F3 = TWO*(TWO*X0+B)/Q
      X  = SQRT(RS)
      FX = X*X+B*X+C
      QX = ATAN(Q/(TWO*X+B))
      EVWN=A*(LOG(RS/FX)+F1*QX-F2*(LOG((X-X0)**2/FX)+F3*QX))
      TXPB=TWO*X+B
      TTQQ=TXPB*TXPB+Q*Q
      VVWN=EVWN - X*A/6.*(TWO/X-TXPB/FX-4.*B/TTQQ-F2*(TWO/(X-X0) &
               -TXPB/FX-4.*(TWO*X0+B)/TTQQ))
      RETURN
      END SUBROUTINE VWN
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE LYP(RHO,ELYP,VLYP)
!     ==--------------------------------------------------------------==
!     ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
!     ==  THIS IS ONLY THE LDA PART                                   ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: a, b, c, d, cf
      PARAMETER (A=0.04918,B=0.132,C=0.2533,D=0.349)
      PARAMETER (CF=2.87123400018819108)
!     ==--------------------------------------------------------------==
      RS=RHO**(-1./3.)
      ECRS=B*CF*EXP(-C*RS)
      OX=1.D0/(1.D0+D*RS)
      ELYP=-A*OX*(1.D0+ECRS)
      VLYP=ELYP-RS/3.*A*OX*(D*OX+ECRS*(D*OX+C))
      RETURN
      END SUBROUTINE LYP
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE PW(RS,EPWC,VPWC,IFLG)
!     ==--------------------------------------------------------------==
!     ==  J.P. PERDEW AND YUE WANG PRB 45, 13244 (1992)               ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
      PARAMETER (A=0.031091,A1=0.21370,B1=7.5957,B2=3.5876,B3=1.6382, &
                 B4=0.49294,C0=A,C1=0.046644,C2=0.00664,C3=0.01043,   &
                 D0=0.4335,D1=1.4408)
!     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
!..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN- (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
!..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1+A1*RS)*OLOG
        VPWC=-2.*A*(1.D0+2./3.*A1*RS)*OLOG-2./3.*A*(1.+A1*RS)*DOM/(OM*(OM+1.))
      ELSEIF(IFLG.EQ.3) THEN
!..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5
      ENDIF
      RETURN
      END SUBROUTINE PW
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE WIGNER(RHO,EXC,FXC)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      RH=RHO
      X=RH**0.33333333333333333D0
      FXC=-X*((0.943656+8.8963*X)/(1.0+12.57*X)**2)
      EXC=-0.738*X*(0.959/(1.0+12.57*X))
      RETURN
      END SUBROUTINE WIGNER
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE HEDIN(RHO,ECP,FCP)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      SAVE RH
      IF(RH .EQ. 0.0D0) RETURN
      RH=RHO
      RSM1=0.62035049D0*RH**(0.3333333333333333D0)
      ALN=DLOG(1.0D0 + 21.0D0*RSM1)
      X=21.0D0/RSM1
      ECP = ALN+(X**3*ALN-X*X)+X/2-1.0D0/3.0D0
      ECP = -0.0225*ECP
      FCP = -0.0225*ALN
      RETURN
      END SUBROUTINE HEDIN
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE OBPZ(RS,EPZ,VPZ,IFLG)
!     ==--------------------------------------------------------------==
!     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
!     ==  PERDEW-ZUNGER FORMULA                                       ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER(I-N)
      REAL(dbl) :: a, b, c, d, gc, b1, b2
      PARAMETER (A=0.031091D0,B=-0.046644D0,C=0.00419D0,D=-0.00983D0, &
                GC=-0.103756D0,B1=0.56371D0,B2=0.27358D0)
!     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
!..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+  (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
!..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7./6.*B1*RS1+4./3.*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
      RETURN
      END SUBROUTINE OBPZ
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE OBPW(RS,EPWC,VPWC,IFLG)
!     ==--------------------------------------------------------------==
!     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
!     ==  PERDEW-WANG FORMULA                                         ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER(I-N)
      REAL(dbl) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
      PARAMETER (A=0.031091,A1=0.026481,B1=7.5957,B2=3.5876,B3=-0.46647, &
                 B4=0.13354,C0=A,C1=0.046644,C2=0.00664,C3=0.01043,      &
                 D0=0.4335,D1=1.4408)
!     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
!..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN- (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
!..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1+A1*RS)*OLOG
        VPWC=-2.*A*(1.D0+2./3.*A1*RS)*OLOG  &
             -2./3.*A*(1.+A1*RS)*DOM/(OM*(OM+1.))
      ELSEIF(IFLG.EQ.3) THEN
!..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5
      ENDIF
      RETURN
      END SUBROUTINE OBPW
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE PADE(RS,EC,VC)
!     ==--------------------------------------------------------------==
!     ==  PADE APPROXIMATION                                          ==
!     ==  S. GOEDECKER, M. TETER, J. HUTTER, PRB in press             ==
!     ==--------------------------------------------------------------==
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N) 
      REAL(dbl) :: a0, a1, a2, a3, b1, b2, b3, b4, O3
      PARAMETER (A0=0.4581652932831429D0,A1=2.217058676663745D0,  &
                 A2=0.7405551735357053D0,A3=0.01968227878617998D0)
      PARAMETER (B1=1.0000000000000000D0,B2=4.504130959426697D0,  &
                 B3=1.110667363742916D0,B4=0.02359291751427506D0)
      PARAMETER (O3=1.D0/3.D0)
!     ==--------------------------------------------------------------==
      TOP=A0+RS*(A1+RS*(A2+RS*A3))
      DTOP=A1+RS*(2.D0*A2+3.D0*A3*RS)
      BOT=RS*(B1+RS*(B2+RS*(B3+RS*B4)))
      DBOT=B1+RS*(2.D0*B2+RS*(3.D0*B3+RS*4.D0*B4))
      EC=-TOP/BOT
      VC=EC+RS*O3*(DTOP/BOT-TOP*DBOT/(BOT*BOT))
      RETURN
      END SUBROUTINE PADE
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE BECKE88(B1,RHO,GRHO,SX,V1X,V2X)
!     ==--------------------------------------------------------------==
! BECKE EXCHANGE: PRA 38, 3098 (1988)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: OB3
      PARAMETER(OB3=1.D0/3.D0)
!     ==--------------------------------------------------------------==
      TWO13 = 2.0D0**(1.D0/3.D0)
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      XS    = TWO13*A/BR4
      XS2   = XS*XS
      SA2B8 = SQRT(1.0D0+XS2)
      SHM1  = LOG(XS+SA2B8)
      DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
      DD2   = DD*DD
      EE    = 6.0D0*B1*XS2/SA2B8 - 1.D0
      SX    = TWO13*AA/BR4*(-B1/DD)
      V1X   = -(4./3.)/TWO13*XS2*B1*BR1*EE/DD2
      V2X   = TWO13*B1*(EE-DD)/(BR4*DD2)
      RETURN
      END SUBROUTINE BECKE88
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE GGAX(RHO,GRHO,SX,V1X,V2X)
!     ==--------------------------------------------------------------==
! J.P.PERDEW ET AL. PRB 46 6671 (1992)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: f1, f2, f3, f4, f5, pi
      PARAMETER(F1=0.19645,F2=7.7956,F3=0.2743,F4=0.1508,F5=0.004)
      PARAMETER(PI=3.141592653589793D+00)
!     ==--------------------------------------------------------------==
      FP1   = -3./(16.*PI)*(3.*PI*PI)**(-1./3.)
      FP2   = 0.5*(3.*PI*PI)**(-1./3.)
      AA    = GRHO
      A     = SQRT(AA)
      RR    = RHO**(-4./3.)
      S     = FP2*A*RR
      S2    = S*S
      S3    = S2*S
      S4    = S2*S2
      EXPS  = F4*EXP(-100.*S2)
      AS    = F3-EXPS-F5*S2
      SA2B8 = SQRT(1.0D0+F2*F2*S2)
      SHM1  = LOG(F2*S+SA2B8)
      BS    = 1.D0+F1*S*SHM1+F5*S4
      DAS   = 200.*S*EXPS-2.*S*F5
      DBS   = F1*(SHM1+F2*S/SA2B8)+4.*F5*S3
      DLS   = (DAS/AS-DBS/BS)
      SX    = FP1*AA*RR*AS/BS
      V1X   = -4./3.*SX/RHO*(1.+S*DLS)
      V2X   = FP1*RR*AS/BS*(2.+S*DLS)
      RETURN
      END SUBROUTINE GGAX
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE PBEX(RHO,GRHO,SX,V1X,V2X)
!     ==--------------------------------------------------------------==
! J.P.PERDEW ET AL. PRL XX XXXX (1996)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: us, ax, um, uk, ul
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0, &
                UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK)
!     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4./3.)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
      RETURN
      END SUBROUTINE PBEX
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE PERDEW86(RHO,GRHO,SC,V1C,V2C)
!     ==--------------------------------------------------------------==
! PERDEW CORRELATION: PRB 33, 8822 (1986)
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: p1, p2, p3, p4, pc1, pc2, pci, ob3, fpi
      PARAMETER(P1=0.023266D0,P2=7.389D-6,P3=8.723D0,P4=0.472D0)
      PARAMETER(PC1=0.001667D0,PC2=0.002568,PCI=PC1+PC2)
      PARAMETER(OB3=1.D0/3.D0, FPI=4.0D0*3.141592653589793D0)
!     ==--------------------------------------------------------------==
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      RS    = (3./(FPI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      CNA   = PC2+P1*RS+P2*RS2
      CNB   = 1.+P3*RS+P4*RS2+1.D4*P2*RS3
      CN    = PC1 + CNA/CNB
      DRS   = -OB3*(3./FPI)**OB3 / BR4
      DCNA  = (P1+2.*P2*RS)*DRS
      DCNB  = (P3+2.*P4*RS+3.D4*P2*RS2)*DRS
      DCN   = DCNA/CNB - CNA/(CNB*CNB)*DCNB
      PHI   = 0.192*PCI/CN*A*RHO**(-7./6.)
      EPHI  = EXP(-PHI)
      SC    = AA/BR4*CN*EPHI
      V1C   = SC*((1.+PHI)*DCN/CN -((4./3.)-(7./6.)*PHI)/RHO)
      V2C   = CN*EPHI/BR4*(2.-PHI)
      RETURN
      END SUBROUTINE PERDEW86
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE GLYP(RHO,GRHO,SC,V1C,V2C)
!     ==--------------------------------------------------------------==
! LEE, YANG PARR: GRADIENT CORRECTION PART
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: a, b, c, d
      PARAMETER(A=0.04918,B=0.132,C=0.2533,D=0.349)
!     ==--------------------------------------------------------------==
      AA    = GRHO
      R     = RHO**(-1.d0/3.d0)
      OM    = EXP(-C*R)/(1.d0+D*R)
      R5    = R**5
      XL    = 1.d0+(7.d0/3.d0)*(C*R + D*R/(1.d0+D*R))
      FF    = A*B*AA/24.d0
      SC    = FF*R5*OM*XL
      DR5   = 5.d0*R*R*R*R
      DOM   = -OM*(C+D+C*D*R)/(1.d0+D*R)
      DXL   = (7.d0/3.d0)*(C+D+2.d0*C*D*R+C*D*D*R*R)/(1.d0+D*R)**2
      V1C   = -FF*(R*R*R*R)/3.d0*( DR5*OM*XL + R5*DOM*XL + R5*OM*DXL)
      V2C   = A*B*R5*OM*XL/12.d0
      RETURN
      END SUBROUTINE GLYP
!     ==--------------------------------------------------------------==


!     ==================================================================
      SUBROUTINE GGAC(RHO,GRHO,SC,V1C,V2C)
!     ==--------------------------------------------------------------==
! PERDEW & WANG GGA CORRELATION PART      
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: al, pa, pb, pc, pd, cx, cxc0, cc0, ob3, pi
      PARAMETER(AL=0.09,PA=0.023266D0,PB=7.389D-6,pc=8.723D0,PD=0.472D0)
      PARAMETER(CX=-0.001667D0,CXC0=0.002568,CC0=-CX+CXC0)
      PARAMETER(OB3=1.D0/3.D0, PI=3.141592653589793D0)
!     ==--------------------------------------------------------------==
      XNU   = 16./PI*(3.*PI*PI)**OB3
      BE    = XNU*CC0
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3./(4.*PI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      XKF   = (9.*PI/4.)**OB3/RS
      XKS   = SQRT(4.*XKF/PI)
      T     = A/(2.*XKS*RHO)
      EXPE  = EXP(-2.*AL*EC/(BE*BE))
      AF    = 2.*AL/BE * (1./(EXPE-1.))
      BF    = EXPE*(VC-EC)
      Y     = AF*T*T
      XY    = (1.+Y)/(1.+Y+Y*Y)
      QY    = Y*Y*(2.+Y)/(1.+Y+Y*Y)**2
      S1    = 1.+2.*AL/BE*T*T*XY
      H0    = BE*BE/(2.*AL) * LOG(S1)
      DH0   = BE*T*T/S1*(-7./3.*XY-QY*(AF*BF/BE-7./3.))
      DDH0  = BE/(2.*XKS*XKS*RHO)*(XY-QY)/S1
      EE    = -100.*(XKS/XKF*T)**2
      CNA   = CXC0+PA*RS+PB*RS2
      DCNA  = -(PA*RS+2.*PB*RS2)/3.
      CNB   = 1.+pc*RS+PD*RS2+1.D4*PB*RS3
      DCNB  = -(pc*RS+2.*PD*RS2+3.D4*PB*RS3)/3.
      CN    = CNA/CNB - CX
      DCN   = DCNA/CNB - CNA*DCNB/(CNB*CNB)
      H1    = XNU*(CN-CC0-3./7.*CX)*T*T*EXP(EE)
      DH1   = -OB3*(H1*(7.+8.*EE)+XNU*T*T*EXP(EE)*DCN)
      DDH1  = 2.*H1*(1.+EE)*RHO/AA
      SC    = RHO*(H0+H1)
      V1C   = H0+H1+DH0+DH1
      V2C   = DDH0+DDH1
      RETURN
      END SUBROUTINE GGAC
!     ==--------------------------------------------------------------==



!     ==================================================================
      SUBROUTINE PBEC(RHO,GRHO,SC,V1C,V2C)
!     ==--------------------------------------------------------------==
! PBE Correlation functional
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER(I-N)
      REAL(dbl) :: be, ga, cx, cxc0, cc0, ob3, pi
      PARAMETER(BE=0.06672455060314922D0,GA=0.031090690869654895D0)
      PARAMETER(CX=-0.001667D0,CXC0=0.002568,CC0=-CX+CXC0)
      PARAMETER(OB3=1.D0/3.D0, PI=3.141592653589793D0)
!     ==--------------------------------------------------------------==
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3./(4.*PI*RHO))**OB3
      XKF   = (9.*PI/4.)**OB3/RS
      XKS   = SQRT(4.*XKF/PI)
      T     = A/(2.*XKS*RHO)
      EXPE  = EXP(-EC/GA)
      AF    = BE/GA * (1./(EXPE-1.))
      Y     = AF*T*T
      XY    = (1.+Y)/(1.+Y+Y*Y)
      S1    = 1.+BE/GA*T*T*XY
      H0    = GA * LOG(S1)
      DTDR  = -T*7./(6.*RHO)
      DADR  = AF*AF*EXPE/BE*(VC-EC)/RHO
      DSDA  = -BE/GA * AF * T**6 * (2.+Y) / (1.+Y+Y*Y)**2
      DSDT  = 2.*BE/GA * T * (1.+2.*Y) / (1.+Y+Y*Y)**2
      DSDR  = DSDA*DADR + DSDT*DTDR
      DHDT  = GA/S1*DSDT
      DHDR  = GA/S1*DSDR
      SC    = RHO*H0
      V1C   = H0+DHDR*RHO
      V2C   = RHO*DHDT*T/AA
      RETURN
      END  SUBROUTINE PBEC
!     ==--------------------------------------------------------------==


!     ==================================================================


      end module gradient_corrections
