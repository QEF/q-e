!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE local_density_approximation

        USE kinds

        IMPLICIT NONE
        SAVE

        INTEGER :: narray

        REAL(dbl) ::  ddro
        REAL(dbl) ::  rmxxc
        REAL(dbl), ALLOCATABLE :: eexc(:)
        REAL(dbl), ALLOCATABLE :: vvxc(:)

        PUBLIC :: lda_setup

      CONTAINS
  
        SUBROUTINE lda_setup(narray_inp, rmxxc_inp)
          INTEGER, INTENT(IN) :: narray_inp
          REAL(dbl), INTENT(IN) :: rmxxc_inp
            narray = narray_inp
            rmxxc  = rmxxc_inp
            CALL genxc()
          RETURN
        END SUBROUTINE lda_setup

        SUBROUTINE lda_print_info(iunit)
          INTEGER, INTENT(IN) :: iunit
            WRITE(iunit,10) RMXXC,NARRAY
 10 FORMAT(   3X,'Using Local Density Approximation',/ &
             ,3X,'functional type : SLATER, PERDEW & ZUNGER',/ &
             ,3X,'Using an interpolating table with : ',/ &
             ,3X,'  Density maximum value = ',E10.4,/ &
             ,3X,'  No. of table entries  = ',I6)
          RETURN
        END SUBROUTINE lda_print_info

        SUBROUTINE allocate_lda
          ALLOCATE(EEXC(0:NARRAY))
          ALLOCATE(VVXC(0:NARRAY))
          RETURN
        END SUBROUTINE allocate_lda

        SUBROUTINE deallocate_lda
          IF( ALLOCATED( EEXC ) ) DEALLOCATE(EEXC)
          IF( ALLOCATED( VVXC ) )  DEALLOCATE(VVXC)
          RETURN
        END SUBROUTINE deallocate_lda

        SUBROUTINE vofxc_lda(rhoe, vpot, sxc, vxc)
          USE io_global, ONLY: stdout
          REAL(dbl) :: RHOE(:,:,:)
          REAL(dbl) :: vpot(:,:,:)
          REAL(dbl) :: sxc, vxc
          REAL(dbl) :: ratio, dd, ddi, ee, vxc1, exc1, OBYDRO, roe
          INTEGER   :: i, j, k, i1
          LOGICAL   :: tmax, tmin

          sxc = 0.0d0
          vxc = 0.0d0
          tmax = .FALSE.
          tmin = .FALSE.

          OBYDRO = 1.D0/DDRO
          DO k = 1, SIZE(rhoe,3)
            DO j = 1, SIZE(rhoe,2)
              DO i = 1, SIZE(rhoe,1)
                ROE   = RHOE(i,j,k)
                RATIO = OBYDRO * ROE
                I1    = INT(RATIO)
                DDI   = REAL(I1)
                DD    = RATIO - DDI
                EE    = 1.D0 - DD
                IF(I1 .GE. NARRAY) THEN
                   I1 = NARRAY - 1
                   tmax = .TRUE.
                END IF
                IF(I1 .LT. 0) then
                   I1 = 0
                   tmin = .TRUE.
                END IF
                VXC1 = VVXC(I1)*EE + VVXC(I1+1)*DD
                EXC1 = EEXC(I1)*EE + EEXC(I1+1)*DD
                SXC = SXC + EXC1*ROE
                VXC = VXC + VXC1*ROE
                VPOT(i,j,k) = VXC1
              END DO
            END DO
          END DO

          IF( tmax ) THEN
            WRITE( stdout, fmt= &
            '(" ** vofxc_lda: WARNING charge values excede the table max value")')
          END IF
          IF( tmin ) THEN
            WRITE( stdout, fmt= &
            '(" ** vofxc_lda: WARNING negative charge values")')
          END IF

          RETURN
        END SUBROUTINE vofxc_lda

!
!---------------------------------------------------------------
      REAL(dbl) FUNCTION ECCA(RS,IFLG)
!---------------------------------------------------------------
!
      REAL(dbl), INTENT(IN) :: RS
      INTEGER :: IFLG
      REAL(dbl) :: A(2),B(2),C(2),D(2),G(2),B1(2),B2(2), RSL, RSQ
      DATA A/0.0622D0,0.0311D0/, B/-0.096D0,-0.0538D0/, C/0.0040D0,0.0014D0/,&
       D/-0.0232D0,-0.0096D0/, B1/1.0529D0,1.3981D0/, B2/0.3334D0,0.2611D0/,&
       G/-0.2846D0,-0.1686D0/
      IF(RS.LE.1.D0) THEN
        RSL=LOG(RS)
        ECCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
        RSQ=SQRT(RS)
        ECCA=G(IFLG)/(1.D0+B1(IFLG)*RSQ+B2(IFLG)*RS)
      END IF
      RETURN
      END FUNCTION ECCA
!
!---------------------------------------------------------------
      REAL(dbl) FUNCTION ECGL(RS,IFLG)
!---------------------------------------------------------------
!
      REAL(dbl), INTENT(IN) :: RS
      INTEGER :: IFLG
      REAL(dbl) :: C(2), R(2), X
      REAL(dbl), PARAMETER :: THIRD=1.D0/3.D0
      DATA C/0.0666D0,0.0406D0/, R/11.4D0,15.9D0/
      X=RS/R(IFLG)
      ECGL=-C(IFLG) * ( (1.D0+X**3)*LOG(1.D0+1.D0/X) - THIRD +X*(0.5D0-X) )
      RETURN
      END FUNCTION ECGL
!
!---------------------------------------------------------------
      REAL(dbl) FUNCTION ECHL(RS,IFLG)
!---------------------------------------------------------------
!
      REAL(dbl), INTENT(IN) :: RS
      INTEGER :: IFLG
      REAL(dbl) :: C(2), R(2), X
      REAL(dbl), PARAMETER :: THIRD=1.D0/3.D0
      DATA C/0.045D0,0.0225D0/, R/21.D0,52.917D0/
      X=RS/R(IFLG)
      ECHL=-C(IFLG) * ( (1.D0+X**3)*LOG(1.D0+1.D0/X) - THIRD +X*(0.5D0-X) )
      RETURN
      END FUNCTION ECHL
!
!--------------------------------------------------------------
      REAL(dbl) FUNCTION EXC(RHO,ZETA,IXC)
!---------------------------------------------------------------
!
!   INTERPOLATION FORMULA FOR THE EXCHANGE-CORRELATION ENERGY
!   DENSITY:
!      IXC =    0  EXCHANGE ONLY (KOHN & SHAM)
!           +/- 1  GUNNARSON & LUNDQVIST
!           +/- 2  HEDIN & LUNDQVIST
!           +/- 3  CEPERLEY & ALDER
!          >    0  EXCHANGE + CORRELATION
!          <    0  CORRELATION ONLY
!
 
      IMPLICIT REAL(dbl) (A-H,O-Z)
      REAL(dbl) :: rho, zeta
      INTEGER :: ixc
      INTEGER :: ic, ifl

      REAL(dbl), PARAMETER :: PI34 = 0.75D0 / 3.141592653589793D+00 
      REAL(dbl), PARAMETER :: ALPHA=0.5210617612D+00
      REAL(dbl), PARAMETER :: THIRD=1.D0/3.D0  
      REAL(dbl), PARAMETER :: TW43=2.519842099D0
      REAL(dbl), PARAMETER :: SMALL=1.0D-10
      IC=IABS(IXC)
      IF(IC.GT.3) STOP '     IXC NOT ALLOWED IN EXC'
      EXC=0.D0
      IF(RHO.LE.SMALL) RETURN
      RS=(PI34/RHO)**THIRD
      IFL=IC+1
      IF(IXC.GE.0) EXC=-2.D0*PI34/RS/ALPHA
      GO TO (40,10,20,30) IFL
  10  EXC=EXC+ECGL(RS,1)
      GO TO 40
  20  EXC=EXC+ECHL(RS,1)
      GO TO 40
  30  EXC=EXC+ECCA(RS,1)
  40  IF(ZETA.EQ.0.D0) RETURN
      EXCF=0.D0
      IF(IXC.GE.0) EXCF=-TW43*PI34/RS/ALPHA
      GO TO (41,11,21,31) IFL
  11  EXCF=EXCF+ECGL(RS,2)
      GO TO 41
  21  EXCF=EXCF+ECHL(RS,2)
      GO TO 41
  31  EXCF=EXCF+ECCA(RS,2)
  41  EXC=EXC+FZ(ZETA,RS)*(EXCF-EXC)
      RETURN
      END FUNCTION EXC
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FZ(ZETA,RS)
!---------------------------------------------------------------
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      INTEGER IFLAG
      SAVE IFLAG
      IFLAG = 0
      IF(IFLAG.EQ.0) THEN
        FZ=FZ0(ZETA)
      ELSE
        FZ=FZX(ZETA,RS)
      END IF
      RETURN
      END  FUNCTION FZ
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FZ0(ZETA)
!---------------------------------------------------------------
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      REAL(dbl), PARAMETER :: F43= 4.D0/3.D0
      REAL(dbl), PARAMETER :: DEN= 5.198420997897450D-01
      ZS=ABS(ZETA)
      IF( ZS-1.D0 .GT. 1.0D-10 ) STOP '     ZETA OUT OF RANGE IN FZ0'
      ZS=MIN(ZS,1.D0)
      IF(ZS.EQ.0.D0) THEN
        FZ0=0.D0
      ELSE IF(ZS.EQ.1.D0) THEN
        FZ0=1.D0
      ELSE
        FZ0=( (1.D0+ZS)**F43 + (1.D0-ZS)**F43 -2.D0 ) / DEN
      END IF
      RETURN
      END FUNCTION FZ0
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FZX(ZETA,RS)
!---------------------------------------------------------------
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)

      REAL(dbl) :: RSX(9),A0(9),A1(9),A2(9),A3(9),B0(9),B1(9),B2(9) &
         ,B3(9),D0(9),D1(9),D2(9),D3(9)

      REAL(dbl), PARAMETER :: FZ2=1.709920933D0
      REAL(dbl), PARAMETER :: RSXM=15.D0

      DATA RSX   / 5.0000D-01, 1.0000D+00, 2.0000D+00, 3.0000D+00, &
       4.0000D+00, 5.0000D+00, 6.0000D+00, 7.5000D+00, 1.0000D+01/ &
      ,A0        / 9.7980D+01, 7.7100D+01, 5.7860D+01, 4.7620D+01, &
       4.0920D+01, 3.6090D+01, 3.2390D+01, 2.8180D+01, 2.3290D+01/ &
      ,A1        /-5.3220D+01,-3.1618D+01,-1.2133D+01,-8.2905D+00, &
      -5.5253D+00,-4.1984D+00,-3.2712D+00,-2.4165D+00,-1.5615D+00/ &
      ,A2        / 2.5555D+01, 1.7649D+01, 1.8363D+00, 2.0062D+00, &
       7.5897D-01, 5.6795D-01, 3.5924D-01, 2.1057D-01, 1.3144D-01/ &
      ,A3        /-5.2708D+00,-5.2708D+00, 5.6645D-02,-4.1574D-01, &
      -6.3674D-02,-6.9568D-02,-3.3039D-02,-1.0550D-02,-1.0550D-02/ &
      ,B0        / 2.3870D-01, 1.9430D-01, 1.4660D-01, 1.1920D-01, &
       1.0080D-01, 8.7500D-02, 7.7400D-02, 6.6100D-02, 5.3300D-02/ &
      ,B1        /-1.0887D-01,-7.0850D-02,-3.3049D-02,-2.2253D-02, &
      -1.5338D-02,-1.1496D-02,-8.8781D-03,-6.4154D-03,-4.0142D-03/ &
      ,B2        / 4.4399D-02, 3.1650D-02, 6.1515D-03, 4.6443D-03, &
       2.2714D-03, 1.5702D-03, 1.0478D-03, 5.9401D-04, 3.6650D-04/ &
      ,B3        /-8.4994D-03,-8.4994D-03,-5.0241D-04,-7.9097D-04, &
      -2.3372D-04,-1.7416D-04,-1.0083D-04,-3.0335D-05,-3.0335D-05/  
      DATA D0    / 7.0977D+01, 5.3849D+01, 3.8797D+01, 3.1170D+01, &
       2.6345D+01, 2.2952D+01, 2.0408D+01, 1.7570D+01, 1.4348D+01/ &
      ,D1        /-4.4080D+01,-2.5573D+01,-9.0941D+00,-6.0874D+00, &
      -3.9121D+00,-2.9180D+00,-2.2268D+00,-1.6129D+00,-1.0129D+00/ &
      ,D2        / 2.1929D+01, 1.5084D+01, 1.3947D+00, 1.6120D+00, &
       5.6329D-01, 4.3084D-01, 2.6036D-01, 1.4891D-01, 9.1108D-02/ &
      ,D3        /-4.5632D+00,-4.5632D+00, 7.2443D-02,-3.4957D-01, &
      -4.4150D-02,-5.6827D-02,-2.4766D-02,-7.7072D-03,-7.7072D-03/  
      IF(RS.LT.RSX(1).OR.RS.GT.RSXM) THEN
      FZX=FZ0(ZETA)
      ELSE
      I=1
  10  I=I+1
      IF(RSX(I).LT.RS) GO TO 10
      I=I-1
      DRS=RS-RSX(I)
      A=A0(I)+DRS*(A1(I)+DRS*(A2(I)+DRS*A3(I)))
      B=B0(I)+DRS*(B1(I)+DRS*(B2(I)+DRS*B3(I)))
      D=D0(I)+DRS*(D1(I)+DRS*(D2(I)+DRS*D3(I)))
      DRPA=A*FZ0(ZETA)*(1.D0+B*ZETA**4)/FZ2
      FZX=DRPA/D
      END IF
      RETURN
      END FUNCTION FZX
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VCCA(RS,IFLG)
!---------------------------------------------------------------
!
!       CEPERLEY & ALDER'S CORRELATION POTENTIAL.
!       AFTER J.P. PERDEW & A. ZUNGER PRB 23, 5048 (1981)
!       IFLG = 1: PARAMAGNETIC  RESULTS
!            = 2: FERROMAGNETIC RESULTS
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      REAL(dbl) :: A(2),B(2),C(2),D(2),BT1(2),BT2(2)
      REAL(dbl) :: X76, X43, AP, BP, CP, DP, AF, BF, CF, DF
      REAL(dbl) :: BP1, CP1, DP1, BF1, CF1, DF1
      PARAMETER(X76=7.D0/6.D0, X43=4.D0/3.D0,  &
        AP=0.03110D0*2.D0, BP=-0.0480D0*2.D0,  &
        CP=0.0020D0*2.D0,  DP=-0.0116D0*2.D0,  &
        AF=0.01555D0*2.D0, BF=-0.0269D0*2.D0,  &
        CF=0.0007D0*2.D0,  DF=-0.0048D0*2.D0,  &
        BP1=BP-AP/3.D0, CP1=2.D0*CP/3.D0, DP1=(2.D0*DP-CP)/3.D0,      &
        BF1=BF-AF/3.D0, CF1=2.D0*CF/3.D0, DF1=(2.D0*DF-CF)/3.D0 )
      DATA A/AP ,AF /, B/BP1,BF1/, C/CP1,CF1/, D/DP1,DF1/,            &
           BT1/1.0529D0,1.3981D0/, BT2/0.3334D0,0.2611D0/
!
      IF(RS.LE.1.D0) THEN
        RSL=LOG(RS)
        VCCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
        RSQ=SQRT(RS)
        VCCA=ECCA(RS,IFLG)*(1.D0+X76*BT1(IFLG)*RSQ+X43*BT2(IFLG)*RS)/  &
                           (1.D0+    BT1(IFLG)*RSQ+    BT2(IFLG)*RS)
      END IF
      RETURN
      END FUNCTION VCCA
!
!***************************************************************
!*****                                                     *****
!*****   THIS SECTION CONTAINS SEVERAL PARAMETRIZATIONS    *****
!*****   FOR THE EXCHANGE-CORRELATION ENERGY AND POTENTIAL *****
!*****                                                     *****
!***************************************************************
!
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VCGL(RS,IFLG)
!---------------------------------------------------------------
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      REAL(dbl), PARAMETER :: CP=0.0666D0
      REAL(dbl), PARAMETER :: RP=11.4D0
      VCGL=-CP*LOG(1.D0+RP/RS)
      RETURN
      END FUNCTION VCGL
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VCHL(RS,IFLG)
!---------------------------------------------------------------
!
!   INTERPOLATION FORMULA FOR THE CORRELATION POTENTIAL AFTER
!   L. HEDIN & S. LUNDQVIST J.PHYS. C4, 2064 (1971).
!   IFLG=1 : POTENTIEL DE CORRELATION PARAMAGNETIQUE;
!   IFLG=2 : POTENTIEL DE CORRELATION FERROMAGNETIQUE (MISAWA).
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      REAL(dbl) :: R(2),C(2)
      DATA R/21.0D0,52.917D0/, C/0.045D0,0.0225D0/
      X=RS/R(IFLG)
      VCHL=-C(IFLG)*LOG(1.0D0+1.D0/X)
      RETURN
      END FUNCTION VCHL
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VX(RS,IFLG)
!---------------------------------------------------------------
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      REAL(dbl) :: CX(2)
      DATA CX/1.221774115D0, 1.539338926D0/
      VX=-CX(IFLG)/RS
      RETURN
      END FUNCTION VX
!
!---------------------------------------------------------------
      DOUBLE PRECISION FUNCTION VXC0(RHO,IXC)
!---------------------------------------------------------------
!
!   INTERPOLATION FORMULA FOR THE EXCHANGE-CORRELATION POTENTIAL
!   IN THE PARAMAGNETIC CASE:
!      IXC =    0  EXCHANGE ONLY (KOHN & SHAM)
!           +/- 1  GUNNARSON & LUNDQVIST
!           +/- 2  HEDIN & LUNDQVIST
!           +/- 3  CEPERLEY & ALDER
!          >    0  EXCHANGE + CORRELATION
!          <    0  CORRELATION ONLY
!
      IMPLICIT REAL(dbl) (A-H,O-Z) , INTEGER (I-N)
      REAL(dbl) :: PI34, ALPHA, THIRD, CX, SMALL
      PARAMETER(PI34= 0.75D0 / 3.141592653589793D+00    &
               ,ALPHA=0.5210617612D+00, THIRD=1.D0/3.D0  &
               ,CX=8.D0*PI34*THIRD/ALPHA, SMALL=1.0D-10)
!
      IC=IABS(IXC)
      IF(IC.GT.3) STOP '     IXC NOT ALLOWED IN VXC0'
      VXC0=0.D0
      IF(RHO.LE.SMALL) RETURN
      RS=(PI34/RHO)**THIRD
      IFL=IC+1
      IF(IXC.GE.0) VXC0=-CX/RS
      GO TO (40,10,20,30) IFL
  10  VXC0=VXC0+VCGL(RS,1)
      GO TO 40
  20  VXC0=VXC0+VCHL(RS,1)
      GO TO 40
  30  VXC0=VXC0+VCCA(RS,1)
  40  RETURN
      END  FUNCTION VXC0
!
!===============================================================
!
      SUBROUTINE XC(IRO,RO1,RO2,EX,EC,VX1,VX2,VC1,VC2,T58)
!
      IMPLICIT REAL(dbl) (A-H,O-Z), INTEGER (I-N)
      IF(IRO.NE.1) STOP '    WRONG IRO IN XC'
      EX=EXC(RO1,0.D0,0)/2.D0
      EC=EXC(RO1,0.D0,-3)/2.D0
      VX1=VXC0(RO1,0)/2.D0
      VX2=0.D0
      VC1=VXC0(RO1,-3)/2.D0
      VC2=0.D0
      RETURN
      END SUBROUTINE XC


!
!===============================================================
!
      SUBROUTINE GENXC
!
!........THIS PROGRAM GENERATES THE FILE TABXC. THE FILE TABXC
!........CONTAINS A TABULATION OF XC ENERGIES AND POTENTIALS ON
!........A GRID CONSISTING OF NARRAY ELEMENTS UNIFORMLY SPACED
!........BETWEEN 0 AND RMAX.
!........RECOMMENDED VALUES FOR SILICON: RMAX=0.2, NARRAY=10000
!........RECOMMENDED VALUES FOR GALLIUM: RMAX=0.2, NARRAY=10000
!........RECOMMENDED VALUES FOR CARBON : RMAX=2.0, NARRAY=40000
!........IBM VERSION (2.AUG.85).
!
      IMPLICIT NONE
      
      REAL(dbl)  RR,EX,EC,VX,VC,DUM
      INTEGER I

      call allocate_lda

!      WRITE( stdout,10) RMXXC,NARRAY
      DDRO=RMXXC/REAL(NARRAY-1)
!     WRITE( stdout,40)
      DO I=0,NARRAY-1
        RR=REAL(I)*DDRO
        CALL XC(1,RR,DUM,EX,EC,VX,DUM,VC,DUM,DUM)
        EEXC(I)=EX+EC
        VVXC(I)=VX+VC
!       IF(MOD(I,100).EQ.0) THEN
!         WRITE( stdout,30) RR,EX,EC,EEXC(I),VX,VC,VVXC(I)
!       END IF
      END DO
!
   10 FORMAT(//'   GENXC: RMAX, NARRAY = ',E10.4,3X,I5//)
   40 FORMAT(3X,'RHO',8X,'EX',8X,'EC',8X,'EEXC',8X,'VX',8X,'VC', 8X,'VVXC')
   30 FORMAT(1X,F9.6,6(1X,F9.6))
!
      RETURN
      END SUBROUTINE GENXC

      END MODULE local_density_approximation
