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

        REAL(dbl) :: alat = 0.0d0   !  lattice parameter, often used to scale quantities
                                    !  or in combination to other parameters/constants
                                    !  to define new units

        !  celldm are che simulation cell parameters

        REAL(dbl) :: celldm(6) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)

        !  a1, a2 and a3 are the simulation cell base vector as calculated from celldm

        REAL(dbl) :: a1(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        REAL(dbl) :: a2(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        REAL(dbl) :: a3(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)

        REAL(dbl) :: ainv(3,3) = 0.0d0

        REAl(dbl) :: omega = 0.0d0  !  volume of the simulation cell

        REAL(dbl) :: tpiba  = 0.0d0   !  = 2 PI / alat
        REAL(dbl) :: tpiba2 = 0.0d0   !  = ( 2 PI / alat ) ** 2

        !  direct lattice vectors and reciprocal lattice vectors
        !  The folloving relations should alwais be kept valid
        !  at( :, 1 ) = a1( : ) / alat  ;  h( :, 1 ) = a1( : )
        !  at( :, 2 ) = a2( : ) / alat  ;  h( :, 2 ) = a2( : )
        !  at( :, 3 ) = a3( : ) / alat  ;  h( :, 3 ) = a3( : )
        !  ht = h^t ; ainv = h^(-1)
        !
        !  bg( :, 1 ) = b1( : )
        !  bg( :, 2 ) = b2( : )
        !  bg( :, 3 ) = b3( : )

        REAL(dbl) :: at(3,3) = RESHAPE( (/ 0.0d0 /), (/ 3, 3 /), (/ 0.0d0 /) )
        REAL(dbl) :: bg(3,3) = RESHAPE( (/ 0.0d0 /), (/ 3, 3 /), (/ 0.0d0 /) )

        INTEGER          :: ibrav      ! index of the bravais lattice
        CHARACTER(len=9) :: symm_type  ! 'cubic' or 'hexagonal' when ibrav=0

        REAL(dbl) :: h(3,3) = 0.0d0
        REAL(dbl) :: hold(3,3) = 0.0d0
        REAL(dbl) :: deth = 0.0d0

        REAL(dbl) :: wmass = 0.0d0

        LOGICAL :: tcell_base_init = .FALSE.

        INTERFACE cell_init
          MODULE PROCEDURE cell_init_ht, cell_init_a
        END INTERFACE

        INTERFACE pbcs
          MODULE PROCEDURE pbcs_components, pbcs_vectors
        END INTERFACE

        INTERFACE s_to_r
          MODULE PROCEDURE s_to_r1, s_to_r3
        END INTERFACE

        INTERFACE r_to_s
          MODULE PROCEDURE r_to_s1, r_to_s3
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
            CALL invmat(3,box_t0%G,GCM1,DUM)
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

        SUBROUTINE R_TO_S1 (R,S,box)
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
        END SUBROUTINE R_TO_S1

        SUBROUTINE R_TO_S3 ( R, S, na, nsp, hinv )
          REAL(dbl), intent(out) ::  S(:,:,:)
          INTEGER, intent(in) ::  na(:), nsp
          REAL(dbl), intent(in) :: R(:,:,:)
          REAL(dbl), intent(in) :: hinv(:,:)    ! hinv = TRANSPOSE( box%m1 )
          integer :: i, j, ia, is
          DO is = 1, nsp
            DO ia = 1, na(is)
              DO I=1,3
                S(I,ia,is) = 0.D0
                DO J=1,3
                  S(I,ia,is) = S(I,ia,is) + R(J,ia,is)*hinv(i,j)
                END DO
              END DO
            END DO
          END DO
          RETURN
        END SUBROUTINE R_TO_S3

!------------------------------------------------------------------------------!

        SUBROUTINE S_TO_R1 (S,R,box)
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
        END SUBROUTINE S_TO_R1

        SUBROUTINE S_TO_R3 ( S, R, na, nsp, h )
          REAL(dbl), intent(in) ::  S(:,:,:)
          INTEGER, intent(in) ::  na(:), nsp
          REAL(dbl), intent(out) :: R(:,:,:)
          REAL(dbl), intent(in) :: h(:,:)    ! h = TRANSPOSE( box%a )
          integer :: i, j, ia, is
          DO is = 1, nsp
            DO ia = 1, na(is)
              DO I = 1, 3
                R(I,ia,is) = 0.D0
                DO J = 1, 3
                  R(I,ia,is) = R(I,ia,is) + S(J,ia,is) * h(I,j)
                END DO
              END DO
            END DO
          END DO
          RETURN
        END SUBROUTINE S_TO_R3


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
          CALL errore('gethinv', 'box determinant too small', 1)
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

        CALL invmat(3,box%a,box%m1,box%omega)

        IF(abs(box%omega-box%deth)/abs(box%omega+box%deth).gt.1.0d-12) THEN
          CALL errore('gethinv', 'box determinants are different',2)
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
          s = dble(nl)
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

!------------------------------------------------------------------------------!

  SUBROUTINE cell_base_init( ibrav_ , celldm_ , trd_ht, cell_symmetry, rd_ht, &
               a, b, c, cosab, cosac, cosbc, alat_ )

    USE constants, ONLY: bohr_radius_angs

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ibrav_
    REAL(dbl), INTENT(IN) :: celldm_ (6)
    LOGICAL, INTENT(IN) :: trd_ht
    CHARACTER(LEN=*), INTENT(IN) :: cell_symmetry
    REAL(dbl), INTENT(IN) :: rd_ht (3,3)
    REAL(dbl), INTENT(IN) :: a, b, c, cosab, cosac, cosbc
    REAL(dbl), INTENT(OUT) :: alat_

    REAL(dbl) :: b1(3), b2(3), b3(3)
    !
    ! ... set up crystal lattice, and initialize cell_base module
    !

    celldm = celldm_
    ibrav  = ibrav_

    ! ... if celldm(1) /= 0  rd_ht should be in unit of alat

    IF ( trd_ht ) THEN
      symm_type = cell_symmetry
      at        = TRANSPOSE( rd_ht )
    END IF

    IF ( ibrav == 0 .AND. .NOT. trd_ht ) &
      CALL errore( ' cell_base_init ', ' ibrav=0: must read cell parameters', 1 )
    IF ( ibrav /= 0 .AND. trd_ht ) &
      CALL errore( ' cell_base_init ', ' redundant data for cell parameters', 2 )


    IF ( celldm(1) == 0.D0 .AND. a /= 0.D0 ) THEN
      IF ( ibrav == 0 ) ibrav = 14
      celldm(1) = a / bohr_radius_angs
      celldm(2) = b / a
      celldm(3) = c / a
      celldm(4) = cosab
      celldm(5) = cosac
      celldm(6) = cosbc
    ELSE IF ( celldm(1) /= 0.D0 .AND. a /= 0.D0 ) THEN
      CALL errore( ' cell_base_init ', ' do not specify both celldm and a,b,c!', 1 )
    END IF

    IF ( ibrav == 0 .AND. celldm(1) /= 0.D0 ) THEN
      !
      ! ... input at are in units of alat
      !
      alat = celldm(1)
    ELSE IF ( ibrav == 0 .AND. celldm(1) == 0.D0 ) THEN
      !
      ! ... input at are in atomic units: define alat
      !
      celldm(1) = SQRT( at(1,1)**2 + at(1,2)**2 + at(1,3)**2 )
      alat = celldm(1)
      !
      ! ... bring at to alat units
      !
      at(:,:) = at(:,:) / alat
    ELSE
      !
      ! ... generate at (atomic units)
      !
      CALL latgen( ibrav, celldm, at(1,1), at(1,2), at(1,3), omega )
      alat = celldm(1)
      !
      ! ... bring at to alat units
      !
      at(:,:) = at(:,:) / alat
    END IF

    !
    a1  =  at( :, 1 ) * alat
    a2  =  at( :, 2 ) * alat
    a3  =  at( :, 3 ) * alat

    CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )

    CALL recips( a1, a2, a3, b1, b2, b3 )
    ainv( 1, : ) = b1( : )
    ainv( 2, : ) = b2( : )
    ainv( 3, : ) = b3( : )

    bg( :, 1 ) = b1( : )
    bg( :, 2 ) = b2( : )
    bg( :, 3 ) = b3( : )

    tcell_base_init = .TRUE.
    alat_ = alat

    RETURN
  END SUBROUTINE
      


  SUBROUTINE cell_verlet( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, hnos )
    REAL(kind=8), INTENT(OUT) :: hnew(3,3)
    REAL(kind=8), INTENT(IN) :: h(3,3), hold(3,3), hnos(3,3), fcell(3,3)
    INTEGER,      INTENT(IN) :: iforceh(3,3)
    REAL(kind=8), INTENT(IN) :: frich, delt
    LOGICAL,      INTENT(IN) :: tnoseh

    REAL(kind=8) :: htmp(3,3)
    REAL(kind=8) :: verl1, verl2, verl3, dt2, ftmp
    INTEGER      :: i, j
  
    dt2 = delt * delt

    IF( tnoseh ) THEN
      ftmp = 0.0d0
      htmp = hnos
    ELSE
      ftmp = frich
      htmp = 0.0d0
    END IF

    verl1 = 2. / ( 1. + ftmp )
    verl2 = 1. - verl1
    verl3 = dt2 / ( 1. + ftmp )
  
    DO j=1,3
      DO i=1,3
        hnew(i,j) = h(i,j) + ( ( verl1 - 1.0d0 ) * h(i,j)               &
     &     + verl2 * hold(i,j) + &
             verl3 * ( fcell(i,j) - htmp(i,j) ) ) * iforceh(i,j)
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE

!
!------------------------------------------------------------------------------!
   END MODULE cell_base
!------------------------------------------------------------------------------!
