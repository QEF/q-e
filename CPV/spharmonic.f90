!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

     module spherical_harmonics

        USE kinds
        implicit none
        save

        PRIVATE

        PUBLIC :: spharm, set_dmqm, gkl, set_fmrm

     contains

      SUBROUTINE spharm(S,G,GSQM,NG,L,M)
! ...   L = 0, 1, 2     angular momentum
! ...   M = -L, ..., L  magnetic quantum number
! ...   NG = number of plane wave
! ...   G(:,:)  = cartesian components of the reciprocal space vectors
! ...   GSQM(:) = square modulus of the reciprocal space vectors
! ...   S(:)    = spherical harmonic components
        IMPLICIT NONE
        REAL(dbl), INTENT(OUT) :: S(:)
        REAL(dbl), INTENT(IN) :: G(:,:), GSQM(:)
        REAL(dbl) :: x,y,z,r,r2
        INTEGER, INTENT(IN) :: NG, M, L
        INTEGER :: i, mm

        MM = M + L + 1
        SELECT CASE (L)
          CASE (0)
            S = 1.0d0
          CASE (1)
            DO i = 1, ng
              IF(GSQM(i).GE.1.0d-12) THEN
                S(i) = G(mm,i) / sqrt(GSQM(i))
              ELSE
                S(i) = 0.0d0
              END IF
            END DO
          CASE (2)
            DO i = 1, ng
              IF(GSQM(i).GE.1.0d-12) THEN
                ! S(i) = gkl(G(1,i),G(2,i),G(3,i),GSQM(i),mm)
                x = G(1,i)
                y = G(2,i)
                z = G(3,i)
                r2 = GSQM(i)
                SELECT CASE (mm)
                  CASE (1)
                    s(i) =  0.50d0 * (3.D0*Z*Z/r2-1.D0)
                  CASE (2)
                    s(i) =  0.50d0 * (X*X-Y*Y)/r2 * SQRT(3.D0)
                  CASE (3)
                    s(i) =  X*Y/r2 * SQRT(3.D0)
                  CASE (4)
                    s(i) = -Y*Z/r2 * SQRT(3.D0)
                  CASE (5)
                    s(i) = -Z*X/r2 * SQRT(3.D0)
                  CASE DEFAULT
                    CALL errore(' GKL ',' magnetic moment not implementent ',mm)
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
                    S(i) = 0.5d0 * SQRT(5.0d0/2.0d0) * y * ( 3.0d0 * x**2 - y**2 ) / r**3
                  CASE( -2 )
                    S(i) = SQRT(15.0d0) * x * y * z / r**3
                  CASE( -1 )
                    S(i) = 0.5d0 * SQRT(3.0d0/2.0d0) * ( 5.0d0 * z**2 / r**2 - 1.0d0 ) * y / r
                  CASE(  0 )
                    S(i) = 0.5d0 * ( 5.0d0 * z**3 / r**3 - 3.0d0 * z / r )
                  CASE(  1 )
                    S(i) = 0.5d0 * SQRT(3.0d0/2.0d0) * ( 5.0d0 * z**2 / r**2 - 1.0d0 ) * x / r
                  CASE(  2 )
                    S(i) = 0.5d0 * SQRT(15.0d0) * z * ( x**2 - y**2 ) / r**3
                  CASE(  3 )
                    S(i) = 0.5d0 * SQRT(5.0d0/2.0d0) * x * ( x**2 - 3.0d0 * y**2 ) / r**3
                END SELECT
              ELSE
                S(i) = 0.0d0
              END IF
            END DO
          CASE DEFAULT
            CALL errore(' spharm ',' angular momuntum not implementent ',l)
        END SELECT

        RETURN
      END SUBROUTINE spharm

!================================================================
!== Spherical harmonics (l=2) in cartesian coordinates
!================================================================
      REAL(dbl) FUNCTION GKL(X,Y,Z,SQM,M)
      IMPLICIT NONE
      REAL(dbl) :: X,Y,Z,SQM
      INTEGER :: M

      GKL = 0.0d0

      SELECT CASE (M)
        CASE (1)
          GKL=(3.D0*Z*Z/SQM-1.D0)/2.D0
        CASE (2)
          GKL=(X*X-Y*Y)/SQM*SQRT(3.D0)/2.D0
        CASE (3)
          GKL= X*Y/SQM*SQRT(3.D0)
        CASE (4)
          GKL=-Y*Z/SQM*SQRT(3.D0)
        CASE (5)
          GKL=-Z*X/SQM*SQRT(3.D0)
        CASE DEFAULT
          CALL errore(' GKL ',' magnetic moment not implementent ',m)
      END SELECT

      RETURN
      END FUNCTION GKL


!=========================================================================
!== DM is the matrix used to construct the spherical harmonics with L=D ==
!==      Y_M(G) = DM(A,B,M) * G_A * G_B      [A,B=1,3 => K=1,6]         ==
!==---------------------------------------------------------------------==
!== DMQM is the product of DM * QM , QM being the "inverse" of DM, i.e. ==
!==      G_A * G_B = QM(A,B,M) * Y_M(G) + DELTA_A,B / 3                 ==
!=========================================================================
!
      SUBROUTINE SET_DMQM(DM,DMQM)

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
      D_M(1,1,1)=-1.D0/2.D0
      D_M(2,2,1)=-1.D0/2.D0
      D_M(3,3,1)= 1.D0
! M=2
      D_M(1,1,2)= SQRT(3.D0)/2.D0
      D_M(2,2,2)=-SQRT(3.D0)/2.D0
! M=3
      D_M(1,2,3)=  SQRT(3.D0)/2.D0
      D_M(2,1,3)=  SQRT(3.D0)/2.D0
! M=4
      D_M(2,3,4)= -SQRT(3.D0)/2.D0
      D_M(3,2,4)= -SQRT(3.D0)/2.D0
! M=5
      D_M(1,3,5)= -SQRT(3.D0)/2.D0
      D_M(3,1,5)= -SQRT(3.D0)/2.D0

!======== QM ================
! M=1
      Q_M(1,1,1)=-1.D0/3.D0
      Q_M(2,2,1)=-1.D0/3.D0
      Q_M(3,3,1)= 2.D0/3.D0
! M=2
      Q_M(1,1,2)=  1.D0/SQRT(3.D0)
      Q_M(2,2,2)= -1.D0/SQRT(3.D0)
! M=3
      Q_M(1,2,3)= 1.D0/SQRT(3.D0)
      Q_M(2,1,3)= 1.D0/SQRT(3.D0)
! M=4
      Q_M(2,3,4)= -1.D0/SQRT(3.D0)
      Q_M(3,2,4)= -1.D0/SQRT(3.D0)
! M=5
      Q_M(1,3,5)= -1.D0/SQRT(3.D0)
      Q_M(3,1,5)= -1.D0/SQRT(3.D0)

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
      END SUBROUTINE SET_DMQM




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
      SUBROUTINE SET_FMRM(FM,FMRM)

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
      F_M(1,1,2,1)=  1.D0/2.D0 * SQRT(5.0d0/2.0d0)
      F_M(1,2,1,1)=  1.D0/2.D0 * SQRT(5.0d0/2.0d0)
      F_M(2,1,1,1)=  1.D0/2.D0 * SQRT(5.0d0/2.0d0)
      F_M(2,2,2,1)= -1.D0/2.D0 * SQRT(5.0d0/2.0d0)
! M=-2
      F_M(1,2,3,2)= SQRT(15.0d0)/6.0d0
      F_M(1,3,2,2)= SQRT(15.0d0)/6.0d0
      F_M(2,1,3,2)= SQRT(15.0d0)/6.0d0
      F_M(2,3,1,2)= SQRT(15.0d0)/6.0d0
      F_M(3,1,2,2)= SQRT(15.0d0)/6.0d0
      F_M(3,2,1,2)= SQRT(15.0d0)/6.0d0
! M=-1
      F_M(1,1,2,3)= -SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(1,2,1,3)= -SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(2,1,1,3)= -SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(2,2,2,3)= -SQRT(3.0d0/2.0d0) / 2.0d0
      F_M(3,3,2,3)=  SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(3,2,3,3)=  SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(2,3,3,3)=  SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
! M=0
      F_M(1,1,3,4)= -1.0d0/2.0d0
      F_M(1,3,1,4)= -1.0d0/2.0d0
      F_M(3,1,1,4)= -1.0d0/2.0d0
      F_M(2,2,3,4)= -1.0d0/2.0d0
      F_M(2,3,2,4)= -1.0d0/2.0d0
      F_M(3,2,2,4)= -1.0d0/2.0d0      
      F_M(3,3,3,4)=  1.0d0
! M=1
      F_M(1,1,1,5)= -SQRT(3.0d0/2.0d0) / 2.0d0
      F_M(2,2,1,5)= -SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(2,1,2,5)= -SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(1,2,2,5)= -SQRT(3.0d0/2.0d0) / 6.0d0
      F_M(3,3,1,5)=  SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(3,1,3,5)=  SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
      F_M(1,3,3,5)=  SQRT(3.0d0/2.0d0) * 2.0d0 / 3.0d0
! M=2
      F_M(3,1,1,6)=  SQRT(15.0d0) / 6.0d0      
      F_M(1,3,1,6)=  SQRT(15.0d0) / 6.0d0
      F_M(1,1,3,6)=  SQRT(15.0d0) / 6.0d0
      F_M(3,2,2,6)= -SQRT(15.0d0) / 6.0d0
      F_M(2,3,2,6)= -SQRT(15.0d0) / 6.0d0
      F_M(2,2,3,6)= -SQRT(15.0d0) / 6.0d0
! M=3
      F_M(1,1,1,7)=  SQRT(5.0d0/2.0d0) / 2.0d0
      F_M(2,2,1,7)= -SQRT(5.0d0/2.0d0) / 2.0d0
      F_M(2,1,2,7)= -SQRT(5.0d0/2.0d0) / 2.0d0
      F_M(1,2,2,7)= -SQRT(5.0d0/2.0d0) / 2.0d0


      FM = F_M


!======== FM ================
! M=-3
      R_M(1,1,2,1)=  1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(1,2,1,1)=  1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,1,1,1)=  1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,2,2,1)= -1.D0 / SQRT(5.0d0/2.0d0) / 2.0d0
! M=-2
      R_M(1,2,3,2)= 1.0d0 / SQRT(15.0d0) 
      R_M(1,3,2,2)= 1.0d0 / SQRT(15.0d0)
      R_M(2,1,3,2)= 1.0d0 / SQRT(15.0d0)
      R_M(2,3,1,2)= 1.0d0 / SQRT(15.0d0)
      R_M(3,1,2,2)= 1.0d0 / SQRT(15.0d0)
      R_M(3,2,1,2)= 1.0d0 / SQRT(15.0d0)
! M=-1
      R_M(1,1,2,3)= -1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(1,2,1,3)= -1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,1,1,3)= -1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,2,2,3)= -3.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(3,3,2,3)=  2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(3,2,3,3)=  2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(2,3,3,3)=  2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
! M=0
      R_M(1,1,3,4)= -1.0d0/5.0d0 
      R_M(1,3,1,4)= -1.0d0/5.0d0 
      R_M(3,1,1,4)= -1.0d0/5.0d0 
      R_M(2,2,3,4)= -1.0d0/5.0d0 
      R_M(2,3,2,4)= -1.0d0/5.0d0 
      R_M(3,2,2,4)= -1.0d0/5.0d0      
      R_M(3,3,3,4)=  2.0d0/5.0d0
! M=1
      R_M(1,1,1,5)= -3.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,2,1,5)= -1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(2,1,2,5)= -1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(1,2,2,5)= -1.0d0 / SQRT(3.0d0/2.0d0) / 10.0d0
      R_M(3,3,1,5)=  2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(3,1,3,5)=  2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
      R_M(1,3,3,5)=  2.0d0 / SQRT(3.0d0/2.0d0) / 5.0d0
! M=2
      R_M(3,1,1,6)=  1.0d0 / SQRT(15.0d0)    
      R_M(1,3,1,6)=  1.0d0 / SQRT(15.0d0)
      R_M(1,1,3,6)=  1.0d0 / SQRT(15.0d0)
      R_M(3,2,2,6)= -1.0d0 / SQRT(15.0d0)
      R_M(2,3,2,6)= -1.0d0 / SQRT(15.0d0)
      R_M(2,2,3,6)= -1.0d0 / SQRT(15.0d0)
! M=3
      R_M(1,1,1,7)=  1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,2,1,7)= -1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(2,1,2,7)= -1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0
      R_M(1,2,2,7)= -1.0d0 / SQRT(5.0d0/2.0d0) / 2.0d0

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
      END SUBROUTINE SET_FMRM

     end module
