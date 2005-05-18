!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
     MODULE spherical_harmonics

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: set_dmqm, set_fmrm, set_pmtm

     CONTAINS

      SUBROUTINE mspharm( lm1x, ngw, gx, gg, ylm)
        INTEGER, INTENT(IN) :: lm1x, ngw
        REAL(dbl), INTENT(IN) :: gx(:,:), gg(:)
        REAL(dbl), INTENT(OUT) :: ylm(:,:)
        INTEGER :: l, m, iy
        iy  = 0
!        DO l = 0, lm1x
!          DO m = -l, l
!            iy = iy + 1
!            CALL spharm( ylm(:,iy), gx, gg, ngw, l, m )
!          END DO
!        END DO
!        DO l = 0, lm1x
!          DO m = 0, l
!            iy = iy + 1
!            CALL spharm( ylm(:,iy), gx, gg, ngw, l, m )
!            IF( m > 0 ) THEN
!              iy = iy + 1
!              CALL spharm( ylm(:,iy), gx, gg, ngw, l, -m )
!            END IF
!          END DO
!        END DO
        CALL ylmr2( (lm1x+1)**2, ngw, gx, gg, ylm )
        RETURN
      END SUBROUTINE mspharm


      SUBROUTINE spharm(S,G,GSQM,NG,L,M)

        USE constants, only: fpi

! ...   L = 0, 1, 2     angular momentum
! ...   M = -L, ..., L  magnetic quantum number
! ...   NG = number of plane wave
! ...   G(:,:)  = cartesian components of the reciprocal space vectors
! ...   GSQM(:) = square modulus of the reciprocal space vectors
! ...   S(:)    = spherical harmonic components

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NG, M, L
        REAL(dbl), INTENT(OUT) :: S(:)
        REAL(dbl), INTENT(IN) :: G(:,:), GSQM(:)
        REAL(dbl) :: x,y,z,r,r2
        INTEGER :: i, mm

        MM = M + L + 1
        SELECT CASE (L)
          CASE (0)
            S = 1.0d0
          CASE (1)
            DO i = 1, ng
              IF(GSQM(i).GE.1.0d-12) THEN
                SELECT CASE ( m )
                  CASE(  0 )  ! ok
                    S(i) =  G(3,i) / sqrt(GSQM(i))
                  CASE(  1 )  ! ok
                    S(i) = -G(1,i) / sqrt(GSQM(i))
                  CASE( -1 )  ! ok
                    S(i) = -G(2,i) / sqrt(GSQM(i))
                  CASE DEFAULT
                    CALL errore(' spharm ',' magnetic moment not implementent ',mm)
                END SELECT
              ELSE
                S(i) = 0.0d0
              END IF
            END DO
          CASE (2)
            DO i = 1, ng
              IF(GSQM(i).GE.1.0d-12) THEN
                x = G(1,i)
                y = G(2,i)
                z = G(3,i)
                r2 = GSQM(i)
                SELECT CASE (m)
                  CASE ( 0)
                    s(i) =  0.50d0 * (3.D0*Z*Z/r2-1.D0)
                  CASE ( 2)
                    s(i) =  0.50d0 * (X*X-Y*Y)/r2 * SQRT(3.D0)
                  CASE (-2)
                    s(i) =  X*Y/r2 * SQRT(3.D0)
                  CASE (-1)
                    s(i) = -Y*Z/r2 * SQRT(3.D0)
                  CASE ( 1)
                    s(i) = -Z*X/r2 * SQRT(3.D0)
                  CASE DEFAULT
                    CALL errore(' spharm ',' magnetic moment not implementent ',mm)
                END SELECT
              ELSE
                S(i) = 0.0d0
              END IF
            END DO
          CASE (3)
            DO i = 1, ng
              IF(GSQM(i).GE.1.0d-12) THEN
                x = G(1,i)
                y = G(2,i)
                z = G(3,i)
                r = SQRT(GSQM(i))
                SELECT CASE(M)
                  CASE( -3 )
                    S(i) = -0.5d0 * SQRT(5.0d0/2.0d0) * y * ( 3.0d0 * x**2 - y**2 ) / r**3
                  CASE( -2 )
                    S(i) = SQRT(15.0d0) * x * y * z / r**3
                  CASE( -1 )
                    S(i) = -0.5d0 * SQRT(3.0d0/2.0d0) * ( 5.0d0 * z**2 / r**2 - 1.0d0 ) * y / r
                  CASE(  0 )
                    S(i) = 0.5d0 * ( 5.0d0 * z**3 / r**3 - 3.0d0 * z / r )
                  CASE(  1 )
                    S(i) = -0.5d0 * SQRT(3.0d0/2.0d0) * ( 5.0d0 * z**2 / r**2 - 1.0d0 ) * x / r
                  CASE(  2 )
                    S(i) = 0.5d0 * SQRT(15.0d0) * z * ( x**2 - y**2 ) / r**3
                  CASE(  3 )
                    S(i) = -0.5d0 * SQRT(5.0d0/2.0d0) * x * ( x**2 - 3.0d0 * y**2 ) / r**3
                END SELECT
              ELSE
                S(i) = 0.0d0
              END IF
            END DO
          CASE DEFAULT
            CALL errore(' spharm ',' angular momuntum not implementent ',l)
        END SELECT

        s( 1:ng ) = s( 1:ng ) * sqrt ( dble(2*l+1) / fpi)

        RETURN
      END SUBROUTINE spharm

!=========================================================================
!== PM is the matrix used to construct the spherical harmonics with L=P ==
!==      Y_M(G) = PM(A,M) * G_A                         [A=1,3]         ==
!==---------------------------------------------------------------------==
!== PMTM is the product of PM * TM , TM being the "inverse" of TM, i.e. ==
!==      G_A = TM(A,M) * Y_M(G)                                         ==
!=========================================================================
!
      SUBROUTINE set_pmtm(pm,pmtm)

      IMPLICIT NONE
      REAL(dbl) P_M(3,3),T_M(3,3)
      REAL(dbl) PM(3,3),PMTM(6,3,3)
      INTEGER A,M1,M2,KK

      INTEGER, dimension(6), parameter :: ALPHA = (/ 1,2,3,2,3,3 /)
      INTEGER, dimension(6), parameter :: BETA  = (/ 1,1,1,2,2,3 /)

      P_M = 0.0d0
      T_M = 0.0d0

!======== PM ================
! M=1
      P_M(2,1)= -1.D0
! M=2
      P_M(3,2)=  1.D0
! M=3
      P_M(1,3)= -1.D0

      PM = P_M

!======== TM ================
! M=1
      T_M(2,1)= -1.0d0
! M=2
      T_M(3,2)=  1.0d0
! M=3
      T_M(1,3)= -1.0d0

!======= PMTM =============

      DO M1=1,3
        DO M2=1,3
          DO KK=1,6
            PMTM(KK,M1,M2)= P_M(BETA(KK),M1)*T_M(ALPHA(KK),M2)
          END DO
        END DO
      END DO


      RETURN
      END SUBROUTINE set_pmtm


!=========================================================================
!== DM is the matrix used to construct the spherical harmonics with L=D ==
!==      Y_M(G) = DM(A,B,M) * G_A * G_B      [A,B=1,3 => K=1,6]         ==
!==---------------------------------------------------------------------==
!== DMQM is the product of DM * QM , QM being the "inverse" of DM, i.e. ==
!==      G_A * G_B = QM(A,B,M) * Y_M(G) + DELTA_A,B / 3                 ==
!=========================================================================
!
      SUBROUTINE set_dmqm(dm,dmqm)

      IMPLICIT NONE
      REAL(dbl) D_M(3,3,5),Q_M(3,3,5)
      REAL(dbl) DM(6,5),DMQM(6,5,5)
      INTEGER A,M1,M2,KK

      INTEGER, dimension(6), parameter :: ALPHA = (/ 1,2,3,2,3,3 /)
      INTEGER, dimension(6), parameter :: BETA  = (/ 1,1,1,2,2,3 /)

      D_M = 0.0d0
      Q_M = 0.0d0

!======== DM ================
! M=1
      D_M(1,1,3)=-1.D0/2.D0
      D_M(2,2,3)=-1.D0/2.D0
      D_M(3,3,3)= 1.D0
! M=2
      D_M(1,1,5)= SQRT(3.D0)/2.D0
      D_M(2,2,5)=-SQRT(3.D0)/2.D0
! M=3
      D_M(1,2,1)=  SQRT(3.D0)/2.D0
      D_M(2,1,1)=  SQRT(3.D0)/2.D0
! M=4
      D_M(2,3,2)= -SQRT(3.D0)/2.D0
      D_M(3,2,2)= -SQRT(3.D0)/2.D0
! M=5
      D_M(1,3,4)= -SQRT(3.D0)/2.D0
      D_M(3,1,4)= -SQRT(3.D0)/2.D0

!======== QM ================
! M=1
      Q_M(1,1,3)=-1.D0/3.D0
      Q_M(2,2,3)=-1.D0/3.D0
      Q_M(3,3,3)= 2.D0/3.D0
! M=2
      Q_M(1,1,5)=  1.D0/SQRT(3.D0)
      Q_M(2,2,5)= -1.D0/SQRT(3.D0)
! M=3
      Q_M(1,2,1)= 1.D0/SQRT(3.D0)
      Q_M(2,1,1)= 1.D0/SQRT(3.D0)
! M=4
      Q_M(2,3,2)= -1.D0/SQRT(3.D0)
      Q_M(3,2,2)= -1.D0/SQRT(3.D0)
! M=5
      Q_M(1,3,4)= -1.D0/SQRT(3.D0)
      Q_M(3,1,4)= -1.D0/SQRT(3.D0)

!======= DMQM =============

      DO M1=1,5
        DO M2=1,5
          DO KK=1,6
            DMQM(KK,M1,M2)=0.D0
            DO A=1,3
              DMQM(KK,M1,M2)=DMQM(KK,M1,M2)+ &
                             D_M(BETA(KK),A,M1)*Q_M(A,ALPHA(KK),M2)
            END DO
          END DO
        END DO
      END DO

      DO M1=1,5
        DO KK=1,6
          DM(KK,M1)=D_M(ALPHA(KK),BETA(KK),M1)
        END DO
      END DO

      RETURN
      END SUBROUTINE set_dmqm




!=========================================================================
!== FM is the matrix used to construct the spherical harmonics with L=F ==
!==      Y_M(G) = FM(A,B,C,M) * G_A * G_B * G_C     [A,B,C=1,3]         ==
!==---------------------------------------------------------------------==
!== RM is the "inverse" of FM, i.e.                                     ==
!==      G_A * G_B * G_C = RM(A,B,C,M) * Y_M(G) +                       ==
!==                        DELTA_A,B * G_C / 5.0 +                      ==
!==                        DELTA_B,C * G_A / 5.0 +                      ==
!==                        DELTA_A,C * G_B / 5.0 +                      ==
!=========================================================================
!
      SUBROUTINE set_fmrm(fm,fmrm)

      IMPLICIT NONE
      REAL(dbl) :: F_M(3,3,3,7), R_M(3,3,3,7)
      REAL(dbl) :: FM(3,3,3,7), FMRM(6,7,7)
      INTEGER :: M1,M2,KK,s,k

      INTEGER, dimension(6), parameter :: ALPHA = (/ 1,2,3,2,3,3 /)
      INTEGER, dimension(6), parameter :: BETA  = (/ 1,1,1,2,2,3 /)

      F_M = 0.0d0
      R_M = 0.0d0

!======== FM ================
! M=-3
      F_M(1,1,2,1)= -1.D0/2.D0 * SQRT(5.0d0/2.0d0)
      F_M(1,2,1,1)= -1.D0/2.D0 * SQRT(5.0d0/2.0d0)
      F_M(2,1,1,1)= -1.D0/2.D0 * SQRT(5.0d0/2.0d0)
      F_M(2,2,2,1)= +1.D0/2.D0 * SQRT(5.0d0/2.0d0)
! M=-2
      F_M(1,2,3,2)= SQRT(15.0d0)/6.0d0
      F_M(1,3,2,2)= SQRT(15.0d0)/6.0d0
      F_M(2,1,3,2)= SQRT(15.0d0)/6.0d0
      F_M(2,3,1,2)= SQRT(15.0d0)/6.0d0
      F_M(3,1,2,2)= SQRT(15.0d0)/6.0d0
      F_M(3,2,1,2)= SQRT(15.0d0)/6.0d0
! M=-1
      F_M(1,1,2,3)= +SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(1,2,1,3)= +SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(2,1,1,3)= +SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(2,2,2,3)= +SQRT(3.0d0/2.0d0) / 2.0d0
      F_M(3,3,2,3)= -SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(3,2,3,3)= -SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(2,3,3,3)= -SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
! M=0
      F_M(1,1,3,4)= -1.0d0/2.0d0
      F_M(1,3,1,4)= -1.0d0/2.0d0
      F_M(3,1,1,4)= -1.0d0/2.0d0
      F_M(2,2,3,4)= -1.0d0/2.0d0
      F_M(2,3,2,4)= -1.0d0/2.0d0
      F_M(3,2,2,4)= -1.0d0/2.0d0      
      F_M(3,3,3,4)=  1.0d0
! M=1
      F_M(1,1,1,5)= +SQRT(3.0d0/2.0d0) / 2.0d0
      F_M(2,2,1,5)= +SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(2,1,2,5)= +SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(1,2,2,5)= +SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(3,3,1,5)= -SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(3,1,3,5)= -SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(1,3,3,5)= -SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
! M=2
      F_M(3,1,1,6)=  SQRT(15.0d0) / 6.0d0      
      F_M(1,3,1,6)=  SQRT(15.0d0) / 6.0d0
      F_M(1,1,3,6)=  SQRT(15.0d0) / 6.0d0
      F_M(3,2,2,6)= -SQRT(15.0d0) / 6.0d0
      F_M(2,3,2,6)= -SQRT(15.0d0) / 6.0d0
      F_M(2,2,3,6)= -SQRT(15.0d0) / 6.0d0
! M=3
      F_M(1,1,1,7)= -SQRT(5.0d0/2.0d0) / 2.0d0
      F_M(2,2,1,7)= +SQRT(5.0d0/2.0d0) / 2.0d0
      F_M(2,1,2,7)= +SQRT(5.0d0/2.0d0) / 2.0d0
      F_M(1,2,2,7)= +SQRT(5.0d0/2.0d0) / 2.0d0


      FM = F_M


!======== FM ================
! M=-3
      R_M(1,1,2,1)= -1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(1,2,1,1)= -1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,1,1,1)= -1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,2,2,1)= +1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
! M=-2
      R_M(1,2,3,2)= 1.0d0 / SQRT(15.0d0) 
      R_M(1,3,2,2)= 1.0d0 / SQRT(15.0d0)
      R_M(2,1,3,2)= 1.0d0 / SQRT(15.0d0)
      R_M(2,3,1,2)= 1.0d0 / SQRT(15.0d0)
      R_M(3,1,2,2)= 1.0d0 / SQRT(15.0d0)
      R_M(3,2,1,2)= 1.0d0 / SQRT(15.0d0)
! M=-1
      R_M(1,1,2,3)= +1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(1,2,1,3)= +1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,1,1,3)= +1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,2,2,3)= +3.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(3,3,2,3)= -2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(3,2,3,3)= -2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(2,3,3,3)= -2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
! M=0
      R_M(1,1,3,4)= -1.0d0/5.0d0 
      R_M(1,3,1,4)= -1.0d0/5.0d0 
      R_M(3,1,1,4)= -1.0d0/5.0d0 
      R_M(2,2,3,4)= -1.0d0/5.0d0 
      R_M(2,3,2,4)= -1.0d0/5.0d0 
      R_M(3,2,2,4)= -1.0d0/5.0d0      
      R_M(3,3,3,4)=  2.0d0/5.0d0
! M=1
      R_M(1,1,1,5)= +3.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,2,1,5)= +1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,1,2,5)= +1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(1,2,2,5)= +1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(3,3,1,5)= -2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(3,1,3,5)= -2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(1,3,3,5)= -2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
! M=2
      R_M(3,1,1,6)=  1.0d0 / SQRT(15.0d0)    
      R_M(1,3,1,6)=  1.0d0 / SQRT(15.0d0)
      R_M(1,1,3,6)=  1.0d0 / SQRT(15.0d0)
      R_M(3,2,2,6)= -1.0d0 / SQRT(15.0d0)
      R_M(2,3,2,6)= -1.0d0 / SQRT(15.0d0)
      R_M(2,2,3,6)= -1.0d0 / SQRT(15.0d0)
! M=3
      R_M(1,1,1,7)= -1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,2,1,7)= +1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,1,2,7)= +1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(1,2,2,7)= +1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0

      DO M1=1,7
        DO M2=1,7
          DO KK=1,6
            FMRM(KK,M1,M2)=0.D0
            DO s=1,3
              DO k=1,3
                FMRM(KK,M1,M2) = FMRM(KK,M1,M2) + &
                  F_M(BETA(KK),s,k,M1) * R_M(s,k,ALPHA(KK),M2)
              END DO
            END DO
          END DO
        END DO
      END DO


      RETURN
      END SUBROUTINE set_fmrm

     END MODULE
