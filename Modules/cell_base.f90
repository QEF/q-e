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
          REAL(dbl) :: gvel(3,3) ! metric velocity
          REAL(dbl) :: pail(3,3) ! stress tensor
          REAL(dbl) :: hmat(3,3) ! cell parameters ( transpose of "a" )
          REAL(dbl) :: hvel(3,3) ! cell velocity
          REAL(dbl) :: hinv(3,3)
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

        REAL(dbl) :: h(3,3)    = 0.0d0 ! simulation cell at time t 
        REAL(dbl) :: hold(3,3) = 0.0d0 ! simulation cell at time t-delt
        REAL(dbl) :: deth      = 0.0d0 ! determinant of h ( cell volume )

        INTEGER   :: iforceh(3,3) = 1  ! if iforceh( i, j ) = 0 then h( i, j ) 
                                       ! is not allowed to move
        LOGICAL   :: thdiag = .FALSE.  ! True if only cell diagonal elements 
                                       ! should be updated

        REAL(dbl) :: wmass = 0.0d0     ! cell fictitious mass
        REAL(dbl) :: press = 0.0d0     ! external pressure 

        REAL(dbl) :: frich  = 0.0d0    ! firction parameter for cell damped dynamics
        REAL(dbl) :: greash = 1.0d0    ! greas parameter for damped dynamics

        LOGICAL :: tcell_base_init = .FALSE.

        INTERFACE cell_init
          MODULE PROCEDURE cell_init_ht, cell_init_a
        END INTERFACE

        INTERFACE pbcs
          MODULE PROCEDURE pbcs_components, pbcs_vectors
        END INTERFACE

        INTERFACE s_to_r
          MODULE PROCEDURE s_to_r1, s_to_r1b, s_to_r3
        END INTERFACE

        INTERFACE r_to_s
          MODULE PROCEDURE r_to_s1, r_to_s1b, r_to_s3
        END INTERFACE

!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!
!


        SUBROUTINE updatecell(box_tm1, box_t0, box_tp1)
          type (boxdimensions) :: box_tm1, box_t0, box_tp1
            box_tm1  = box_t0
            box_t0   = box_tp1
          RETURN
        END SUBROUTINE UPDATECELL

!------------------------------------------------------------------------------!

        SUBROUTINE dgcell( gcdot, box_tm1, box_t0, delt )
          REAL(dbl), INTENT(OUT) :: GCDOT(3,3)
          REAL(dbl), INTENT(IN) :: delt
          type (boxdimensions), intent(in) :: box_tm1, box_t0
            !
            GCDOT = 2.0d0 * ( box_t0%g - box_tm1%g ) / delt - box_t0%gvel
            !
            ! this is equivalent to:
            ! GCDOT = (3.D0*box_t0%G - 4.D0*box_tm1%G + box_tm2%G)/ (2.0d0 * delt )
            !
          RETURN 
        END SUBROUTINE dgcell

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
            CALL gethinv( box )
            box%g = MATMUL( box%a(:,:), box%hmat(:,:) )
            box%gvel = 0.0d0
            box%hvel = 0.0d0
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
            box%hvel = 0.0d0
            CALL gethinv(box)
            box%g    = MATMUL( box%a(:,:), box%hmat(:,:) )
            box%gvel = 0.0d0
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

!------------------------------------------------------------------------------!

        SUBROUTINE R_TO_S3 ( R, S, na, nsp, hinv )
          REAL(dbl), intent(out) ::  S(:,:)
          INTEGER, intent(in) ::  na(:), nsp
          REAL(dbl), intent(in) :: R(:,:)
          REAL(dbl), intent(in) :: hinv(:,:)    ! hinv = TRANSPOSE( box%m1 )
          integer :: i, j, ia, is, isa
          isa = 0
          DO is = 1, nsp
            DO ia = 1, na(is)
              isa = isa + 1
              DO I=1,3
                S(I,isa) = 0.D0
                DO J=1,3
                  S(I,isa) = S(I,isa) + R(J,isa)*hinv(i,j)
                END DO
              END DO
            END DO
          END DO
          RETURN
        END SUBROUTINE R_TO_S3

!------------------------------------------------------------------------------!

        SUBROUTINE R_TO_S1B ( R, S, hinv )
          REAL(dbl), intent(out) ::  S(:)
          REAL(dbl), intent(in) :: R(:)
          REAL(dbl), intent(in) :: hinv(:,:)    ! hinv = TRANSPOSE( box%m1 )
          integer :: i, j
          DO I=1,3
            S(I) = 0.D0
            DO J=1,3
              S(I) = S(I) + R(J)*hinv(i,j)
            END DO
          END DO
          RETURN
        END SUBROUTINE R_TO_S1B


!------------------------------------------------------------------------------!

        SUBROUTINE s_to_r1 (S,R,box)
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
        END SUBROUTINE s_to_r1

!------------------------------------------------------------------------------!

        SUBROUTINE s_to_r1b (S,R,h)
          REAL(dbl), intent(in) ::  S(3)
          REAL(dbl), intent(out) :: R(3)
          REAL(dbl), intent(in) :: h(:,:)    ! h = TRANSPOSE( box%a )
          integer i,j
          DO I=1,3
            R(I) = 0.D0
            DO J=1,3
              R(I) = R(I) + S(J)*h(I,j)
            END DO
          END DO
          RETURN
        END SUBROUTINE s_to_r1b

!------------------------------------------------------------------------------!

        SUBROUTINE s_to_r3 ( S, R, na, nsp, h )
          REAL(dbl), intent(in) ::  S(:,:)
          INTEGER, intent(in) ::  na(:), nsp
          REAL(dbl), intent(out) :: R(:,:)
          REAL(dbl), intent(in) :: h(:,:)    ! h = TRANSPOSE( box%a )
          integer :: i, j, ia, is, isa
          isa = 0
          DO is = 1, nsp
            DO ia = 1, na(is)
              isa = isa + 1
              DO I = 1, 3
                R(I,isa) = 0.D0
                DO J = 1, 3
                  R(I,isa) = R(I,isa) + S(J,isa) * h(I,j)
                END DO
              END DO
            END DO
          END DO
          RETURN
        END SUBROUTINE s_to_r3


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
        box%hinv = hmati

        CALL invmat( 3, box%a, box%m1, box%omega )

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

        s = matmul(box%hinv(:,:),rin)
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

!------------------------------------------------------------------------------!

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

  SUBROUTINE cell_base_init( ibrav_ , celldm_ , trd_ht, cell_symmetry, rd_ht,  &
               a_ , b_ , c_ , cosab, cosac, cosbc, wc_ , total_ions_mass , press_ ,  &
               frich_ , greash_ , cell_dofree )

    USE constants, ONLY: bohr_radius_angs, gpa_au, pi, uma_au 
    USE io_global, ONLY: stdout

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ibrav_
    REAL(dbl), INTENT(IN) :: celldm_ (6)
    LOGICAL, INTENT(IN) :: trd_ht
    CHARACTER(LEN=*), INTENT(IN) :: cell_symmetry
    REAL(dbl), INTENT(IN) :: rd_ht (3,3)
    REAL(dbl), INTENT(IN) :: a_ , b_ , c_ , cosab, cosac, cosbc
    CHARACTER(LEN=*), INTENT(IN) :: cell_dofree
    REAL(dbl),  INTENT(IN) :: wc_ , frich_ , greash_ , total_ions_mass
    REAL(dbl),  INTENT(IN) :: press_  ! external pressure from imput ( GPa )


    REAL(dbl) :: b1(3), b2(3), b3(3)
    REAL(dbl) :: a, b, c
    INTEGER   :: j

    !
    ! ... set up crystal lattice, and initialize cell_base module
    !
    celldm = celldm_
    a = a_
    b = b_
    c = c_
    ibrav  = ibrav_
    press  = press_ * gpa_au
    !  frich  = frich_   ! for the time being this is set elsewhere
    greash = greash_

    WRITE( stdout, 105 )
    WRITE( stdout, 110 ) press_
105 format(/,3X,'Simulation Cell Parameters (from input)')
110 format(  3X,'external pressure       = ',f15.2,' [GPa]')

    wmass  = wc_
    IF( wmass == 0.d0 ) THEN
      wmass = 3.d0 / (4.d0 * pi**2 ) * total_ions_mass
      wmass = wmass * UMA_AU
      WRITE( stdout,130) wmass
    ELSE
      WRITE( stdout,120) wmass
    END IF
120 format(3X,'wmass (read from input) = ',f15.2,' [AU]')
130 format(3X,'wmass (calculated)      = ',f15.2,' [AU]')


    ! ... if celldm(1) /= 0  rd_ht should be in unit of alat

    IF ( trd_ht ) THEN
      symm_type = cell_symmetry
      !
      !    The matrix "ht" in FPMD correspond to the transpose of matrix "at" in PW
      !
      at        = TRANSPOSE( rd_ht )
      WRITE( stdout, 210 )
      WRITE( stdout, 220 ) ( rd_ht( 1, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( rd_ht( 2, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( rd_ht( 3, j ), j = 1, 3 )
      IF ( ANY( celldm(1:6) /= 0 ) ) THEN
        WRITE( stdout, 230 )
        celldm(1:6) = 0.0d0
      END IF
      IF ( a /= 0 ) THEN
        WRITE( stdout, 240 )
        a = 0.0d0
        b = 0.0d0
        c = 0.0d0
      END IF
210   format(3X,'initial cell from CELL_PARAMETERS card')
220   format(3X,3F14.8)
230   format(3X,'celldm(1:6) are ignored')
240   format(3X,'a, b, c are ignored')
    END IF

    IF ( ibrav == 0 .AND. .NOT. trd_ht ) &
      CALL errore( ' cell_base_init ', ' ibrav=0: must read cell parameters', 1 )
    IF ( ibrav /= 0 .AND. trd_ht ) &
      CALL errore( ' cell_base_init ', ' redundant data for cell parameters', 2 )
    IF( wmass <= 0.0d0 ) &
      CALL errore(' cell_base_init ',' wmass out of range ',0)


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

    ! ...     The matrix "htm1" in FPMD correspond to the matrix "bg" in PW
    !           

    tcell_base_init = .TRUE.

    thdiag = .false.

    SELECT CASE ( TRIM( cell_dofree ) )

            CASE ( 'all', 'default' )
              iforceh = 1
            CASE ( 'volume' )
              CALL errore(' metric_setup ', &
                 ' cell_dofree = '//TRIM(cell_dofree)//' not yet implemented ', 1 )
            CASE ('x')
              iforceh      = 0
              iforceh(1,1) = 1
            CASE ('y')
              iforceh      = 0
              iforceh(2,2) = 1
            CASE ('z')
              iforceh      = 0
              iforceh(3,3) = 1
            CASE ('xy')
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(2,2) = 1
            CASE ('xz')
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(3,3) = 1
            CASE ('yz')
              iforceh      = 0
              iforceh(2,2) = 1
              iforceh(3,3) = 1
            CASE ('xyz')
              thdiag       = .true.
              iforceh      = 0
              iforceh(1,1) = 1
              iforceh(2,2) = 1
              iforceh(3,3) = 1
            CASE DEFAULT
              CALL errore(' metric_setup ',' unknown cell_dofree '//TRIM(cell_dofree), 1 )

    END SELECT

    WRITE( stdout, 300 ) ibrav
    WRITE( stdout, 305 ) alat
    WRITE( stdout, 310 ) a1
    WRITE( stdout, 320 ) a2
    WRITE( stdout, 330 ) a3
    WRITE( stdout, *   )
    WRITE( stdout, 350 ) b1
    WRITE( stdout, 360 ) b2
    WRITE( stdout, 370 ) b3
    WRITE( stdout, 340 ) omega
300 FORMAT( 3X, 'ibrav = ',I4)
305 FORMAT( 3X, 'alat  = ',F14.8)
310 FORMAT( 3X, 'a1    = ',3F14.8)
320 FORMAT( 3X, 'a2    = ',3F14.8)
330 FORMAT( 3X, 'a3    = ',3F14.8)
350 FORMAT( 3X, 'b1    = ',3F14.8)
360 FORMAT( 3X, 'b2    = ',3F14.8)
370 FORMAT( 3X, 'b3    = ',3F14.8)
340 FORMAT( 3X, 'omega = ',F14.8)


    RETURN
  END SUBROUTINE


      
!------------------------------------------------------------------------------!



  SUBROUTINE cell_base_reinit( ht )

    USE constants, ONLY: pi 
    USE io_global, ONLY: stdout
    USE control_flags, ONLY: iprsta

    IMPLICIT NONE
    REAL(dbl), INTENT(IN) :: ht (3,3)

    REAL(dbl) :: b1(3), b2(3), b3(3)
    INTEGER   :: j

    alat   =  sqrt( ht(1,1)*ht(1,1) + ht(1,2)*ht(1,2) + ht(1,3)*ht(1,3) )
    tpiba  = 2.d0 * pi / alat
    tpiba2 = tpiba * tpiba
    !
    !    The matrix "ht" in FPMD correspond to the transpose of matrix "at" in PW
    !
    at     = TRANSPOSE( ht )
    IF( iprsta > 3 ) THEN
      WRITE( stdout, 210 )
      WRITE( stdout, 220 ) ( ht( 1, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( ht( 2, j ), j = 1, 3 )
      WRITE( stdout, 220 ) ( ht( 3, j ), j = 1, 3 )
    END IF
210 format(3X,'Simulation cell parameters with the new cell:')
220 format(3X,3F14.8)

    !
    a1  =  at( :, 1 )
    a2  =  at( :, 2 )
    a3  =  at( :, 3 )

    at( :, : ) = at( :, : ) / alat

    CALL volume( alat, at(1,1), at(1,2), at(1,3), deth )
    omega = deth

    CALL recips( a1, a2, a3, b1, b2, b3 )
    ainv( 1, : ) = b1( : )
    ainv( 2, : ) = b2( : )
    ainv( 3, : ) = b3( : )

    bg( :, 1 ) = b1( : )
    bg( :, 2 ) = b2( : )
    bg( :, 3 ) = b3( : )

    ! ...     The matrix "htm1" in FPMD correspond to the matrix "bg" in PW

    IF( iprsta > 3 ) THEN
      WRITE( stdout, 305 ) alat
      WRITE( stdout, 310 ) a1
      WRITE( stdout, 320 ) a2
      WRITE( stdout, 330 ) a3
      WRITE( stdout, *   )
      WRITE( stdout, 350 ) b1
      WRITE( stdout, 360 ) b2
      WRITE( stdout, 370 ) b3
      WRITE( stdout, 340 ) omega
    END IF

300 FORMAT( 3X, 'ibrav = ',I4)
305 FORMAT( 3X, 'alat  = ',F14.8)
310 FORMAT( 3X, 'a1    = ',3F14.8)
320 FORMAT( 3X, 'a2    = ',3F14.8)
330 FORMAT( 3X, 'a3    = ',3F14.8)
350 FORMAT( 3X, 'b1    = ',3F14.8)
360 FORMAT( 3X, 'b2    = ',3F14.8)
370 FORMAT( 3X, 'b3    = ',3F14.8)
340 FORMAT( 3X, 'omega = ',F14.8)


    RETURN
  END SUBROUTINE



!------------------------------------------------------------------------------!

  SUBROUTINE cell_steepest( hnew, h, delt, iforceh, fcell )
    REAL(kind=8), INTENT(OUT) :: hnew(3,3)
    REAL(kind=8), INTENT(IN) :: h(3,3), fcell(3,3)
    INTEGER,      INTENT(IN) :: iforceh(3,3)
    REAL(kind=8), INTENT(IN) :: delt
    INTEGER      :: i, j
    REAL(kind=8) :: dt2
    dt2 = delt * delt
    DO j=1,3
      DO i=1,3
        hnew(i,j) = h(i,j) + dt2 * fcell(i,j) * iforceh(i,j)
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE


!------------------------------------------------------------------------------!

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

!------------------------------------------------------------------------------!

  subroutine cell_hmove( h, hold, delt, iforceh, fcell )
    real(kind=8), intent(out) :: h(3,3)
    real(kind=8), intent(in) :: hold(3,3), fcell(3,3)
    real(kind=8), intent(in) :: delt
    integer, intent(in) :: iforceh(3,3)
    real(kind=8) :: dt2by2, fac
    integer :: i, j
    dt2by2 = .5d0 * delt * delt
    fac = dt2by2
    do i=1,3
      do j=1,3
        h(i,j) = hold(i,j) + fac * iforceh(i,j) * fcell(i,j)
      end do
    end do
    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine cell_force( fcell, ainv, stress, omega, press, wmass )
    real(kind=8), intent(out) :: fcell(3,3)
    real(kind=8), intent(in) :: stress(3,3), ainv(3,3)
    real(kind=8), intent(in) :: omega, press, wmass
    integer      :: i, j
    do j=1,3
      do i=1,3
        fcell(i,j) = ainv(j,1)*stress(i,1) + ainv(j,2)*stress(i,2) + ainv(j,3)*stress(i,3)
      end do
    end do
    do j=1,3
      do i=1,3
        fcell(i,j) = fcell(i,j) - ainv(j,i) * press
      end do
    end do
    fcell = omega * fcell / wmass
    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine cell_move( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, vnhh, velh, tsdc )
    real(kind=8), intent(out) :: hnew(3,3)
    real(kind=8), intent(in) :: h(3,3), hold(3,3), fcell(3,3)
    real(kind=8), intent(in) :: vnhh(3,3), velh(3,3)
    integer,      intent(in) :: iforceh(3,3)
    real(kind=8), intent(in) :: frich, delt
    logical,      intent(in) :: tnoseh, tsdc

    real(kind=8) :: hnos(3,3)

    if( tnoseh ) then
      hnos = vnhh * velh
    else
      hnos = 0.0d0
    end if
!
    IF( tsdc ) THEN
      call cell_steepest( hnew, h, delt, iforceh, fcell )
    ELSE
      call cell_verlet( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, hnos )
    END IF

    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine cell_gamma( hgamma, ainv, h, velh )
    implicit none
    real(kind=8) :: hgamma(3,3)
    real(kind=8), intent(in) :: ainv(3,3), h(3,3), velh(3,3)
    integer :: i,j,k,l,m
         do i=1,3
            do j=1,3
               do k=1,3
                  do l=1,3
                     do m=1,3
                        hgamma(i,j)=hgamma(i,j)+ainv(i,l)*ainv(k,l)*    &
     &                       (velh(m,k)*h(m,j)+h(m,k)*velh(m,j))
                     enddo
                  enddo
               enddo
            enddo
         enddo
    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine cell_kinene( ekinh, temphh, velh )
    use constants, only: factem
    implicit none
    real(kind=8), intent(out) :: ekinh, temphh(3,3)
    real(kind=8), intent(in)  :: velh(3,3)
    integer :: i,j
    ekinh = 0.0d0
    do j=1,3
      do i=1,3
        ekinh=ekinh+0.5*wmass*velh(i,j)*velh(i,j)
        temphh(i,j)=factem*wmass*velh(i,j)*velh(i,j)
      end do
    end do
    return
  end subroutine

!------------------------------------------------------------------------------!

  function cell_alat( )
    real(dbl) :: cell_alat
    if( .NOT. tcell_base_init ) &
      call errore( ' cell_alat ', ' alat has not been set ', 1 )
    cell_alat = alat
    return 
  end function

!
!------------------------------------------------------------------------------!
   END MODULE cell_base
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
  MODULE cell_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl
!
      IMPLICIT NONE
      SAVE

      REAL(dbl) :: xnhh0(3,3) = 0.0d0
      REAL(dbl) :: xnhhm(3,3) = 0.0d0
      REAL(dbl) :: xnhhp(3,3) = 0.0d0
      REAL(dbl) :: vnhh(3,3)  = 0.0d0
      REAL(dbl) :: temph      = 0.0d0  !  Thermostat temperature (from input)
      REAL(dbl) :: fnoseh     = 0.0d0  !  Thermostat frequency (from input)
      REAL(dbl) :: qnh        = 0.0d0  !  Thermostat mass (computed)

CONTAINS

  subroutine cell_nose_init( temph_init, fnoseh_init )
     USE constants, ONLY: factem, pi, terahertz
     REAL(dbl), INTENT(IN) :: temph_init, fnoseh_init
     ! set thermostat parameter for cell
     qnh    = 0.0d0
     temph  = temph_init
     fnoseh = fnoseh_init
     if( fnoseh > 0.0d0 ) qnh = 2.d0 * ( 3 * 3 )*temph/factem/(fnoseh*(2.d0*pi)*terahertz)**2
    return
  end subroutine

  subroutine cell_nosezero( vnhh, xnhh0, xnhhm )
    real(dbl), intent(out) :: vnhh(3,3), xnhh0(3,3), xnhhm(3,3)
    xnhh0=0.0d0
    xnhhm=0.0d0
    vnhh =0.0d0
    return
  end subroutine

  subroutine cell_nosevel( vnhh, xnhh0, xnhhm, delt )
    implicit none
    real(kind=8), intent(inout) :: vnhh(3,3)
    real(kind=8), intent(in) :: xnhh0(3,3), xnhhm(3,3), delt
    vnhh(:,:)=2.*(xnhh0(:,:)-xnhhm(:,:))/delt-vnhh(:,:)
    return
  end subroutine

  subroutine cell_noseupd( xnhhp, xnhh0, xnhhm, delt, qnh, temphh, temph, vnhh )
    use constants, only: factem
    implicit none
    real(kind=8), intent(out) :: xnhhp(3,3), vnhh(3,3)
    real(kind=8), intent(in) :: xnhh0(3,3), xnhhm(3,3), delt, qnh, temphh(3,3), temph
    integer :: i, j
    do j=1,3
      do i=1,3
        xnhhp(i,j)=2.*xnhh0(i,j)-xnhhm(i,j)+ (delt**2/qnh)/factem*(temphh(i,j)-temph)
        vnhh(i,j) =(xnhhp(i,j)-xnhhm(i,j))/( 2.0d0 * delt )
      end do
    end do
    return
  end subroutine

  
  real(kind=8) function cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
    use constants, only: factem
    implicit none
    real(kind=8) :: qnh, vnhh( 3, 3 ), temph, xnhh0( 3, 3 )
    integer :: iforceh( 3, 3 )
    integer :: i, j
    real(kind=8) :: enij
    cell_nose_nrg = 0.0d0
    do i=1,3
      do j=1,3
        enij = 0.5*qnh*vnhh(i,j)*vnhh(i,j)+temph/factem*xnhh0(i,j)
        cell_nose_nrg = cell_nose_nrg + iforceh( i, j ) * enij
      enddo
    enddo
    return
  end function

  subroutine cell_nose_shiftvar( xnhhp, xnhh0, xnhhm )
    !  shift values of nose variables to start a new step
    implicit none
    real(kind=8), intent(out) :: xnhhm(3,3)
    real(kind=8), intent(inout) :: xnhh0(3,3)
    real(kind=8), intent(in) :: xnhhp(3,3)
      xnhhm = xnhh0
      xnhh0 = xnhhp
    return
  end subroutine


  SUBROUTINE cell_nose_info()

      use constants,     only: factem, terahertz, pi
      use time_step,     only: delt
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnoseh

      IMPLICIT NONE

      INTEGER   :: nsvar, i
      REAL(dbl) :: wnoseh

      IF( tnoseh ) THEN
        !
        IF( fnoseh <= 0.D0) &
          CALL errore(' cell_nose_info ', ' fnoseh less than zero ', 1)
        IF( delt <= 0.D0) &
          CALL errore(' cell_nose_info ', ' delt less than zero ', 1)

        wnoseh = fnoseh * ( 2.d0 * pi ) * terahertz
        nsvar  = ( 2.d0 * pi ) / ( wnoseh * delt )

        WRITE( stdout,563) temph, nsvar, fnoseh, qnh
      END IF

 563  format( //, &
            & 3X,'cell dynamics with nose` temperature control:', /, &
            & 3X,'Kinetic energy required   = ', f10.5, ' (Kelvin) ', /, &
            & 3X,'time steps per nose osc.  = ', i5, /, &
            & 3X,'nose` frequency           = ', f10.3, ' (THz) ', /, &
            & 3X,'nose` mass(es)            = ', 20(1X,f10.3),//)

    RETURN
  END SUBROUTINE cell_nose_info


!
!------------------------------------------------------------------------------!
   END MODULE cell_nose
!------------------------------------------------------------------------------!

