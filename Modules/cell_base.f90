!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE cell_base
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl
!
      IMPLICIT NONE
      SAVE
!

! ...  periodicity box
! ...  In the matrix "a" every row is the vector of each side of 
! ...  the cell in the real space

        TYPE boxdimensions
          REAL(dbl) :: a(3,3)    ! direct lattice generators
          REAL(dbl) :: m1(3,3)   ! reciprocal lattice generators
          REAL(dbl) :: omega     ! cell volume = determinant of a
          REAL(dbl) :: g(3,3)    ! metric tensor
          REAL(dbl) :: pail(3,3) ! stress tensor
          REAL(dbl) :: hmat(3,3)
          REAL(dbl) :: h_inv(3,3)
          REAL(dbl) :: deth
          INTEGER :: perd(3)
        END TYPE boxdimensions

        INTERFACE cell_init
          MODULE PROCEDURE cell_init_ht, cell_init_a
        END INTERFACE

        INTERFACE pbcs
          MODULE PROCEDURE pbcs_components, pbcs_vectors
        END INTERFACE
!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!
!


        SUBROUTINE updatecell(box_tm2, box_tm1, box_t0, box_tp1)
          type (boxdimensions) :: box_tm2, box_tm1, box_t0, box_tp1
            box_tm2  = box_tm1
            box_tm1  = box_t0
            box_t0   = box_tp1
            box_t0%g = MATMUL( box_t0%a(:,:), TRANSPOSE( box_t0%a(:,:) ) )
            call gethinv(box_t0)
          RETURN
        END SUBROUTINE UPDATECELL

!------------------------------------------------------------------------------!

        SUBROUTINE dgcell(gcm1, gcdot, box_tm2, box_tm1, box_t0, delt)
          REAL(dbl) :: GCM1(3,3)
          REAL(dbl) :: GCDOT(3,3)
          REAL(dbl) :: delt
          type (boxdimensions), intent(in) :: box_tm2, box_tm1, box_t0
          REAL(dbl) :: DUM
            CALL matinv(box_t0%G,GCM1,DUM)
            GCDOT = (3.D0*box_t0%G - 4.D0*box_tm1%G + box_tm2%G)/ (2.0d0 * delt )
          RETURN
        END SUBROUTINE DGCELL

!------------------------------------------------------------------------------!
! ...     set box
! ...     box%m1(i,1) == b1(i)   COLUMN are B vectors
! ...     box%a(1,i)  == a1(i)   ROW are A vector
! ...     box%omega   == volume
! ...     box%g(i,j)  == metric tensor G
!------------------------------------------------------------------------------!

        SUBROUTINE cell_init_ht( box, ht )
          TYPE (boxdimensions) :: box
          REAL(dbl) :: ht(3,3)
            box%a = ht
            box%hmat = TRANSPOSE( ht )
            CALL gethinv(box)
            box%g = MATMUL(box%a(:,:),TRANSPOSE(box%a(:,:)))
            box%pail = 0.0d0
          RETURN
        END SUBROUTINE
          
!------------------------------------------------------------------------------!

        SUBROUTINE cell_init_a( box, a1, a2, a3 )
          TYPE (boxdimensions) :: box
          REAL(dbl) :: a1(3), a2(3), a3(3)
          INTEGER :: i
            DO i=1,3
              box%a(1,I) = A1(I)     ! this is HT: the row are the lattice vectors
              box%a(2,I) = A2(I)
              box%a(3,I) = A3(I)
              box%hmat(I,1) = A1(I)  ! this is H : the column are the lattice vectors
              box%hmat(I,2) = A2(I)
              box%hmat(I,3) = A3(I)
            END DO
            box%pail = 0.0d0
            CALL gethinv(box)
            box%g = MATMUL(box%a(:,:),TRANSPOSE(box%a(:,:)))
          RETURN
        END SUBROUTINE

!------------------------------------------------------------------------------!

        SUBROUTINE R_TO_S (R,S,box)
          REAL(dbl), intent(out) ::  S(3)
          REAL(dbl), intent(in) :: R(3)
          type (boxdimensions), intent(in) :: box
          integer i,j
          DO I=1,3
            S(I) = 0.D0
            DO J=1,3
              S(I) = S(I) + R(J)*box%m1(J,I)
            END DO
          END DO
          RETURN
        END SUBROUTINE R_TO_S

!------------------------------------------------------------------------------!

        SUBROUTINE S_TO_R (S,R,box)
          REAL(dbl), intent(in) ::  S(3)
          REAL(dbl), intent(out) :: R(3)
          type (boxdimensions), intent(in) :: box
          integer i,j
          DO I=1,3
            R(I) = 0.D0
            DO J=1,3
              R(I) = R(I) + S(J)*box%a(J,I)
            END DO
          END DO
          RETURN
        END SUBROUTINE S_TO_R

!------------------------------------------------------------------------------!
!  BEGIN manual

        SUBROUTINE RECIPS(alat, a1, a2, a3, b1, b2, b3, den)

!  this routine computes:
!   b1, b2, b3 the reciprocal lattice base vectors
!   in units of [2pi / alat]
!
!                 a2 x a3
!   b1 = alat --------------  [ 2pi / alat ] 
!             a1 ( a2 x a3 )
!
!  ----------------------------------------------
!  END manual

          REAL(dbl), intent(in)  :: alat
          REAL(dbl), intent(in)  :: a1(3), a2(3), a3(3)
          REAL(dbl), intent(out) :: b1(3), b2(3), b3(3)
          REAL(dbl), intent(out) :: den
          INTEGER I,J,K,IPERM,IR,L
          REAL(dbl)  S
            DEN=0.D0
            I=1; J=2; K=3; S=1.D0
            SIG: DO
              DO IPERM=1,3
                DEN=DEN+S*A1(I)*A2(J)*A3(K)
                L=I; I=J; J=K; K=L
              END DO
              I=2; J=1; K=3; S=-S
              IF(S.LT.0.D0) CYCLE SIG
              EXIT SIG 
            END DO SIG
            I=1; J=2; K=3
            DEN=ALAT/ABS(DEN)
            DO IR=1,3
              B1(IR)=DEN*(A2(J)*A3(K)-A2(K)*A3(J))
              B2(IR)=DEN*(A3(J)*A1(K)-A3(K)*A1(J))
              B3(IR)=DEN*(A1(J)*A2(K)-A1(K)*A2(J))
              L=I; I=J; J=K; K=L
            END DO
          RETURN
        END SUBROUTINE RECIPS

!------------------------------------------------------------------------------!

        SUBROUTINE LATGEN(IBRAV,CELLDM,A1,A2,A3,OMEGA)

!.........SETS UP THE CRYSTALLOGRAPHIC VECTORS A1,A2, AND A3.
!.........IBRAV AND CELLDM ARE DEFINED IN THE TABLE AT THE BEGINNING
!.........OF SUBROUTINE KSUM

          IMPLICIT NONE

          integer, intent(in) :: ibrav
          REAL(dbl), intent(in) :: celldm(6)
          REAL(dbl), intent(out) :: a1(3), a2(3), a3(3)
          REAL(dbl), intent(out) :: omega

          INTEGER IR,I,J,K,L,IPERM
          REAL(dbl) TERM,CBYA,SR3,TERM1,TERM2,SIN,S,SINGAM
          DATA SR3/1.732051D0/
          DO IR=1,3
            A1(IR)=0.D0
            A2(IR)=0.D0
            A3(IR)=0.D0
          END DO
          IF(IBRAV.LT.0 .OR. IBRAV.GT.14) GO TO 110
          GO TO (2,4,6,8,10,12,14,16,18,20,22,24,26,28),IBRAV
    2     A1(1)=CELLDM(1)
          A2(2)=CELLDM(1)
          A3(3)=CELLDM(1)
          GO TO 100
    4     TERM=CELLDM(1)/2.D0
          A1(1)=-TERM
          A1(3)=TERM
          A2(2)=TERM
          A2(3)=TERM
          A3(1)=-TERM
          A3(2)=TERM
          GO TO 100
    6     TERM=CELLDM(1)/2.D0
          DO IR=1,3
            A1(IR)=TERM
            A2(IR)=TERM
            A3(IR)=TERM
          END DO
          A2(1)=-TERM
          A3(1)=-TERM
          A3(2)=-TERM
          GO TO 100
    8     CBYA=CELLDM(3)
          A1(1)=CELLDM(1)
          A2(1)=-CELLDM(1)/2.D0
          A2(2)=CELLDM(1)*SR3/2.D0
          A3(3)=CELLDM(1)*CBYA
          GO TO 100
   10     TERM1=SQRT(1.D0+2.D0*CELLDM(4))
          TERM2=SQRT(1.D0-CELLDM(4))
          A1(2)=1.414214D0*CELLDM(1)*TERM2/SR3
          A1(3)=CELLDM(1)*TERM1/SR3
          A2(1)=CELLDM(1)*TERM2/1.414214D0
          A2(2)=-A2(1)/SR3
          A2(3)=A1(3)
          A3(1)=-A2(1)
          A3(2)=A2(2)
          A3(3)=A1(3)
          GO TO 100
   12     CBYA=CELLDM(3)
          A1(1)=CELLDM(1)
          A2(2)=CELLDM(1)
          A3(3)=CELLDM(1)*CBYA
          GO TO 100
   14     CBYA=CELLDM(3)
          A1(1)=CELLDM(1)/2.D0
          A1(2)=A1(1)
          A1(3)=CBYA*CELLDM(1)/2.D0
          A2(1)=A1(1)
          A2(2)=-A1(1)
          A2(3)=A1(3)
          A3(1)=-A1(1)
          A3(2)=-A1(1)
          A3(3)=A1(3)
          GO TO 100
   16     A1(1)=CELLDM(1)
          A2(2)=CELLDM(1)*CELLDM(2)
          A3(3)=CELLDM(1)*CELLDM(3)
          GO TO 100
   18     GO TO 110
   20     GO TO 110
   22     GO TO 110
   24     SIN=SQRT(1.D0-CELLDM(4)**2)
          A1(1)=CELLDM(1)
          A2(1)=CELLDM(1)*CELLDM(2)*CELLDM(4)
          A2(2)=CELLDM(1)*CELLDM(2)*SIN
          A3(3)=CELLDM(1)*CELLDM(3)
          GO TO 100
   26     GO TO 110
   28     SINGAM=SQRT(1.D0-CELLDM(6)**2)
          TERM=SQRT((1.D0+2.D0*CELLDM(4)*CELLDM(5)*CELLDM(6)  &
           -CELLDM(4)**2-CELLDM(5)**2-CELLDM(6)**2)/ (1.D0-CELLDM(6)**2))
          A1(1)=CELLDM(1)
          A2(1)=CELLDM(1)*CELLDM(2)*CELLDM(6)
          A2(2)=CELLDM(1)*CELLDM(2)*SINGAM
          A3(1)=CELLDM(1)*CELLDM(3)*CELLDM(5)
          A3(2)=CELLDM(1)*CELLDM(3)* (CELLDM(4)-CELLDM(5)*CELLDM(6))/SINGAM
          A3(3)=CELLDM(1)*CELLDM(3)*TERM
  100     OMEGA=0.D0
          S=1.D0; I=1; J=2; K=3
  101     DO 102 IPERM=1,3
          OMEGA=OMEGA+S*A1(I)*A2(J)*A3(K)
          L=I; I=J; J=K; K=L
  102     CONTINUE
          I=2; J=1; K=3; S=-S
          IF(S.LT.0.D0) GO TO 101
          OMEGA=ABS(OMEGA)
          RETURN
  110     WRITE(6,120) IBRAV
  120     FORMAT(' BRAVAIS LATTICE',I3,' NOT PROGRAMMED. STOPPING')
          STOP
        END SUBROUTINE LATGEN

!
!------------------------------------------------------------------------------!
!

      SUBROUTINE gethinv(box)
        IMPLICIT NONE
        TYPE (boxdimensions), INTENT (INOUT) :: box
        REAL (dbl), DIMENSION (3,3) :: hmat, hmati
        REAL (dbl) :: odet

        hmat = box%hmat
        box%deth = hmat(1,1)*(hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2)) + &
          hmat(1,2)*(hmat(2,3)*hmat(3,1)-hmat(2,1)*hmat(3,3)) + &
          hmat(1,3)*(hmat(2,1)*hmat(3,2)-hmat(2,2)*hmat(3,1))
        IF (box%deth<1.E-10) & 
          CALL error('gethinv', 'box determinant too small', 1)
        odet = 1._dbl/box%deth
        hmati(1,1) = (hmat(2,2)*hmat(3,3)-hmat(2,3)*hmat(3,2))*odet
        hmati(2,2) = (hmat(1,1)*hmat(3,3)-hmat(1,3)*hmat(3,1))*odet
        hmati(3,3) = (hmat(1,1)*hmat(2,2)-hmat(1,2)*hmat(2,1))*odet
        hmati(1,2) = (hmat(1,3)*hmat(3,2)-hmat(1,2)*hmat(3,3))*odet
        hmati(2,1) = (hmat(3,1)*hmat(2,3)-hmat(2,1)*hmat(3,3))*odet
        hmati(1,3) = (hmat(1,2)*hmat(2,3)-hmat(1,3)*hmat(2,2))*odet
        hmati(3,1) = (hmat(2,1)*hmat(3,2)-hmat(3,1)*hmat(2,2))*odet
        hmati(2,3) = (hmat(1,3)*hmat(2,1)-hmat(2,3)*hmat(1,1))*odet
        hmati(3,2) = (hmat(3,1)*hmat(1,2)-hmat(3,2)*hmat(1,1))*odet
        box%h_inv = hmati

        CALL matinv(box%a,box%m1,box%omega)

        IF(abs(box%omega-box%deth)/abs(box%omega+box%deth).gt.1.0d-12) THEN
          CALL error('gethinv', 'box determinants are different',2)
        END IF

      END SUBROUTINE gethinv
!
!------------------------------------------------------------------------------!
!
      FUNCTION pbc(rin,box,nl) RESULT (rout)
        IMPLICIT NONE
        TYPE (boxdimensions) :: box
        REAL (dbl) :: rin(3)
        REAL (dbl) :: rout(3), s(3)
        INTEGER, OPTIONAL :: nl(3)

        s = matmul(box%h_inv(:,:),rin)
        s = s - box%perd*nint(s)
        rout = matmul(box%hmat(:,:),s)
        IF (present(nl)) THEN
          s = float(nl)
          rout = rout + matmul(box%hmat(:,:),s)
        END IF
      END FUNCTION pbc
!
!------------------------------------------------------------------------------!
!
          SUBROUTINE get_cell_param(box,cell,ang)
          IMPLICIT NONE
          TYPE(boxdimensions), INTENT(in) :: box
          REAL(dbl), INTENT(out), DIMENSION(3) :: cell
          REAL(dbl), INTENT(out), DIMENSION(3), OPTIONAL :: ang
! This code gets the cell parameters given the h-matrix:
! a
          cell(1)=sqrt(box%hmat(1,1)*box%hmat(1,1)+box%hmat(2,1)*box%hmat(2,1) &
                      +box%hmat(3,1)*box%hmat(3,1))
! b
          cell(2)=sqrt(box%hmat(1,2)*box%hmat(1,2)+box%hmat(2,2)*box%hmat(2,2) &
                      +box%hmat(3,2)*box%hmat(3,2))
! c
          cell(3)=sqrt(box%hmat(1,3)*box%hmat(1,3)+box%hmat(2,3)*box%hmat(2,3) &
                      +box%hmat(3,3)*box%hmat(3,3))
          IF (PRESENT(ang)) THEN
! gamma
             ang(1)=acos((box%hmat(1,1)*box%hmat(1,2)+ &
                          box%hmat(2,1)*box%hmat(2,2) &
                      +box%hmat(3,1)*box%hmat(3,2))/(cell(1)*cell(2)))
! beta
             ang(2)=acos((box%hmat(1,1)*box%hmat(1,3)+ &
                          box%hmat(2,1)*box%hmat(2,3) &
                      +box%hmat(3,1)*box%hmat(3,3))/(cell(1)*cell(3)))
! alpha
           ang(3)=acos((box%hmat(1,2)*box%hmat(1,3)+ &
                        box%hmat(2,2)*box%hmat(2,3) &
                      +box%hmat(3,2)*box%hmat(3,3))/(cell(2)*cell(3)))
!           ang=ang*180.0_dbl/pi

          ENDIF
          END SUBROUTINE get_cell_param

!------------------------------------------------------------------------------!

      SUBROUTINE pbcs_components(x1, y1, z1, x2, y2, z2, m)
! ... This subroutine compute the periodic boundary conditions in the scaled
! ... variables system
        USE kinds
        INTEGER, INTENT(IN)  :: M
        REAL(dbl),  INTENT(IN)  :: X1,Y1,Z1
        REAL(dbl),  INTENT(OUT) :: X2,Y2,Z2
        REAL(dbl) MIC
        MIC = REAL(M)
        X2 = X1 - DNINT(X1/MIC)*MIC
        Y2 = Y1 - DNINT(Y1/MIC)*MIC
        Z2 = Z1 - DNINT(Z1/MIC)*MIC
        RETURN
      END SUBROUTINE


      SUBROUTINE pbcs_vectors(v, w, m)
! ... This subroutine compute the periodic boundary conditions in the scaled
! ... variables system
        USE kinds
        INTEGER, INTENT(IN)  :: m
        REAL(dbl),  INTENT(IN)  :: v(3)
        REAL(dbl),  INTENT(OUT) :: w(3)
        REAL(dbl) :: MIC
        MIC = REAL(M)
        w(1) = v(1) - DNINT(v(1)/MIC)*MIC
        w(2) = v(2) - DNINT(v(2)/MIC)*MIC
        w(3) = v(3) - DNINT(v(3)/MIC)*MIC
        RETURN
      END SUBROUTINE


!
!------------------------------------------------------------------------------!
   END MODULE cell_base
!------------------------------------------------------------------------------!
