!
! Copyright (C) 2013 A. Dal Corso 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE bz_form
    !
    USE kinds, ONLY : DP
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    PUBLIC
    SAVE
    

TYPE bz
    INTEGER :: ind      ! number of the bz
    INTEGER :: nfaces   ! The number of faces
    INTEGER :: nvertices ! The number of vertices
    REAL(DP), ALLOCATABLE :: normal(:,:)  ! The G vector normal to each face
                                          ! in unit 2 pi / celldm(1)
    REAL(DP), ALLOCATABLE :: vertex_coord(:,:) ! coordinates of each vertex
                                           !(carthesian units 2 pi / celldm(1))
    INTEGER, ALLOCATABLE :: ivertex(:,:) ! for each vertex which faces define it
    INTEGER, ALLOCATABLE :: indsur(:,:) ! for each surface the vertex that
                                        ! define it
    INTEGER :: xaxis, yaxis, zaxis   ! the indices of the surfaces that
                                     ! intersect the x, y, and z axis
    REAL(DP) :: xi(3), yi(3), zi(3)  ! the actual coordinates of intersection
    INTEGER :: nlett                 ! number of letters for which the position
                                     ! in the BZ is known
    CHARACTER(LEN=3), ALLOCATABLE :: letter_list(:) ! list of each letter
    REAL(DP), ALLOCATABLE :: letter_coord(:,:) ! coordinates of each letter
    CHARACTER(LEN=20) :: letter_type  ! currently supported convention
                                      ! for the name of the points
    ! Setyawan-Curtarolo (SC) (default)  arXiv:1004.2974
    ! Bilbao    (BI)         Bilbao crystallographic server
    !                        see www.cryst.ehu.es/cryst/get_kvec.html
    ! Bradley-Cracknell (BC) The mathematical theory of symmetry in solids
                            ! is a subset of BI, so use the former

    INTEGER :: npx = 8
 
    INTEGER :: ibrav              ! the bravais lattice index
    REAL(DP) :: celldm(6)         ! the new celldm parameter. In the orthorombic
                                  ! case switches the axis so that a < b < c
    REAL(DP) :: omega             ! volume of the unit cell
    REAL(DP) :: at(3,3), bg(3,3)  ! the direct and reciprocal lattices
                                  ! used in this module. They are switched 
                                  ! in the orthorombic cases 
    LOGICAL :: switch_a_b = .FALSE., &  !  If true the axis are rotated
               switch_b_c = .FALSE., &  !  and the coordinates are switched
               rotate_a_b_c = .FALSE.

END TYPE

CONTAINS

SUBROUTINE allocate_bz(ibrav, ibz, bz_struc, celldm, at, bg )
IMPLICIT NONE
INTEGER, INTENT(IN) :: ibz, ibrav
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: celldm(6), at(3,3), bg(3,3)

bz_struc%ind=ibz
bz_struc%ibrav=ibrav
bz_struc%celldm=celldm
bz_struc%at=at
bz_struc%bg=bg

IF ( ibz ==1) THEN
!
!  simple cubic bz
!
   bz_struc%nfaces=6
   bz_struc%nvertices=8
   bz_struc%nlett=4
   IF (bz_struc%letter_type=='BI') bz_struc%nlett=5
ELSEIF (ibz==2) THEN
!
!  fcc bz
!
   bz_struc%nfaces=14
   bz_struc%nvertices=24
   bz_struc%nlett=6
   IF (bz_struc%letter_type=='BI') bz_struc%nlett=13
ELSEIF (ibz==3) THEN
!
!  bcc bz
!
   bz_struc%nfaces=12
   bz_struc%nvertices=14
   bz_struc%nlett=4
   IF (bz_struc%letter_type=='BI') bz_struc%nlett=5
ELSEIF (ibz==4) THEN
!
!  simple tetragonal bz
!
   bz_struc%nfaces=6
   bz_struc%nvertices=8
   bz_struc%nlett=6
ELSEIF (ibz==5) THEN
!
!  centered tetragonal (c<a) bz
!
   bz_struc%nfaces=12
   bz_struc%nvertices=18
   bz_struc%nlett=7
ELSEIF (ibz==6) THEN
!
!  centered tetragonal (c>a) bz
!
   bz_struc%nfaces=14
   bz_struc%nvertices=24
   bz_struc%nlett=9
   IF (bz_struc%letter_type=='BI') bz_struc%nlett=15
ELSEIF (ibz==7) THEN
!
!  simple orthorombic bz
!
   bz_struc%nfaces=6
   bz_struc%nvertices=8
   bz_struc%nlett=8
ELSEIF (ibz==8) THEN
!
!  face centered orthorombic (1/a^2 > 1/b^2 + 1/c^2) bz case 1
!
   bz_struc%nfaces=12
   bz_struc%nvertices=18
   bz_struc%nlett=9
   CALL check_orthorombic(bz_struc)
ELSEIF (ibz==9) THEN
!
!  face centered orthorombic (1/a^2 < 1/b^2 + 1/c^2) bz case 2
!
   bz_struc%nfaces=14
   bz_struc%nvertices=24
   bz_struc%nlett=11
   CALL check_orthorombic(bz_struc)
ELSEIF (ibz==10) THEN
!
!  face centered orthorombic (1/a^2 = 1/b^2 + 1/c^2) bz case 3
!
   bz_struc%nfaces=12
   bz_struc%nvertices=14
   bz_struc%nlett=8
   CALL check_orthorombic(bz_struc)
ELSEIF (ibz==11) THEN
!
!  body centered orthorombic bz
!
   bz_struc%nfaces=14
   bz_struc%nvertices=24
   bz_struc%nlett=13
   CALL check_orthorombic(bz_struc)
ELSEIF (ibz==12) THEN
!
!  one face centered orthorombic bz
!
   bz_struc%nfaces=8
   bz_struc%nvertices=12
   bz_struc%nlett=10
   CALL check_one_face_centerer_orthorombic(bz_struc)
ELSEIF (ibz==13) THEN
!
!  hexagonal
!
   bz_struc%nfaces=8
   bz_struc%nvertices=12
   bz_struc%nlett=6
ELSEIF (ibz==14) THEN
!
!  trigonal alpha < 90 bz
!
   bz_struc%nfaces=14
   bz_struc%nvertices=24
   bz_struc%nlett=12
ELSEIF (ibz==15) THEN
!
!  trigonal alpha > 90 bz
!
   bz_struc%nfaces=12
   bz_struc%nvertices=14
   bz_struc%nlett=8
ELSEIF (ibz==16) THEN
!
!  simple monoclinic lattice
!
   bz_struc%nfaces=8
   bz_struc%nvertices=12
   bz_struc%nlett=6
ELSE
   CALL errore('allocate_bz','Brillouin zone type not available',1)
ENDIF

ALLOCATE(bz_struc%normal(3,bz_struc%nfaces))
ALLOCATE(bz_struc%ivertex(3,bz_struc%nvertices))
ALLOCATE(bz_struc%vertex_coord(3,bz_struc%nvertices))
ALLOCATE(bz_struc%indsur(9,bz_struc%nfaces))
ALLOCATE(bz_struc%letter_list(bz_struc%nlett))
ALLOCATE(bz_struc%letter_coord(3,bz_struc%nlett))

RETURN
END SUBROUTINE allocate_bz

SUBROUTINE deallocate_bz(bz_struc)
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc

DEALLOCATE(bz_struc%normal)
DEALLOCATE(bz_struc%ivertex)
DEALLOCATE(bz_struc%vertex_coord)
DEALLOCATE(bz_struc%indsur)
DEALLOCATE(bz_struc%letter_list)
DEALLOCATE(bz_struc%letter_coord)

RETURN
END SUBROUTINE deallocate_bz

SUBROUTINE set_label_type(bz_struc, label)
IMPLICIT NONE
TYPE(bz) :: bz_struc
CHARACTER(LEN=*) :: label

IF (TRIM(label)/='SC'.AND.TRIM(label)/='BI'.AND.TRIM(label)/='BC') THEN
   WRITE(stdout,'("Label type not recognized, using default SC")')
   bz_struc%letter_type='SC'
   RETURN
ENDIF

bz_struc%letter_type=TRIM(label)

RETURN
END SUBROUTINE

SUBROUTINE init_bz(bz_struc)
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
INTEGER :: ibz

bz_struc%letter_list(1)='gG '
bz_struc%letter_coord(:,1)=0.0_DP
ibz=bz_struc%ind

IF ( ibz ==1) THEN
   CALL init_bz_1(bz_struc)
ELSEIF (ibz==2) THEN
   CALL init_bz_2(bz_struc)
ELSEIF (ibz==3) THEN
   CALL init_bz_3(bz_struc)
ELSEIF (ibz==4) THEN
   CALL init_bz_4(bz_struc)
ELSEIF (ibz==5) THEN
   CALL init_bz_5(bz_struc)
ELSEIF (ibz==6) THEN
   CALL init_bz_6(bz_struc)
ELSEIF (ibz==7) THEN
   CALL init_bz_7(bz_struc)
ELSEIF (ibz==8) THEN
   CALL init_bz_8(bz_struc)
ELSEIF (ibz==9) THEN
   CALL init_bz_9(bz_struc)
ELSEIF (ibz==10) THEN
   CALL init_bz_10(bz_struc)
ELSEIF (ibz==11) THEN
   CALL init_bz_11(bz_struc)
ELSEIF (ibz==12) THEN
   CALL init_bz_12(bz_struc)
ELSEIF (ibz==13) THEN
   CALL init_bz_13(bz_struc)
ELSEIF (ibz==14) THEN
   CALL init_bz_14(bz_struc)
ELSEIF (ibz==15) THEN
   CALL init_bz_15(bz_struc)
ELSEIF (ibz==16) THEN
   CALL init_bz_16(bz_struc)
ELSE
   CALL errore('init_bz','Brillouin zone type not available',1)
ENDIF

RETURN
END SUBROUTINE init_bz

SUBROUTINE init_bz_1(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  simple cubic bz
!
!
!  G vector normal to each surface
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1)
   bz_struc%normal(:,2)=bz_struc%bg(:,2)
   bz_struc%normal(:,3)=-bz_struc%bg(:,1)
   bz_struc%normal(:,4)=-bz_struc%bg(:,2)
   bz_struc%normal(:,5)=bz_struc%bg(:,3)
   bz_struc%normal(:,6)=-bz_struc%bg(:,3)
!
!  The number of vertice of each surface and its number
!
   bz_struc%indsur(:,1) = (/ 4, 1, 5, 8, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 1, 5, 6, 2 /)
   bz_struc%indsur(:,3) = (/ 4, 2, 6, 7, 3 /)
   bz_struc%indsur(:,4) = (/ 4, 3, 7, 8, 4 /)
   bz_struc%indsur(:,5) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,6) = (/ 4, 5, 6, 7, 8 /)

   CALL find_vertices(bz_struc)
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' M '
   bz_struc%letter_list(3)=' X '
   bz_struc%letter_list(4)=' R '


   bz_struc%letter_coord(:,2)=0.5_DP*( bz_struc%bg(:,1) + bz_struc%bg(:,2) )
   bz_struc%letter_coord(:,3)=0.5_DP*( bz_struc%bg(:,2) )
   bz_struc%letter_coord(:,4)=0.5_DP*( bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                     + bz_struc%bg(:,3) )
   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(5)=' X1'
      bz_struc%letter_coord(:,5)=0.5_DP*( bz_struc%bg(:,1) )
   ENDIF

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_1

SUBROUTINE init_bz_2(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  fcc bz
!
!  standard bg (-1,1,1) (1,-1,1) (1,1,-1)
!
!   bz_struc%normal(:,1)=bg(:,2) + bg(:,3)
!   bz_struc%normal(:,2)=bg(:,1) + bg(:,3)
!   bz_struc%normal(:,3)=-(bg(:,2)+bg(:,3))
!   bz_struc%normal(:,4)=-(bg(:,1) + bg(:,3))
!   bz_struc%normal(:,5)=bg(:,1) + bg(:,2)
!   bz_struc%normal(:,6)=-(bg(:,1)+bg(:,2))
!   bz_struc%normal(:,7)=bg(:,1) + bg(:,2) + bg(:,3)
!   bz_struc%normal(:,8)=bg(:,1) 
!   bz_struc%normal(:,9)=- bg(:,3)
!   bz_struc%normal(:,10)=bg(:,2) 
!   bz_struc%normal(:,11)=bg(:,3)
!   bz_struc%normal(:,12)=-bg(:,2)
!   bz_struc%normal(:,13)=-(bg(:,1) + bg(:,2) + bg(:,3))
!   bz_struc%normal(:,14)=-bg(:,1)
!
!  bg of pwscf (-1,-1,1) (1,1,1) (-1,1,-1)
!
   bz_struc%normal(:,1) =-(bz_struc%bg(:,1) + bz_struc%bg(:,3))
   bz_struc%normal(:,2) =  bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,3) =  bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,4) =-(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%normal(:,5) =  bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,6) =-(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%normal(:,7) =  bz_struc%bg(:,2)
   bz_struc%normal(:,8) = bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                            + bz_struc%bg(:,3)
   bz_struc%normal(:,9) =  bz_struc%bg(:,1)
   bz_struc%normal(:,10)=- bz_struc%bg(:,3) 
   bz_struc%normal(:,11)=- bz_struc%bg(:,1)
   bz_struc%normal(:,12)=  bz_struc%bg(:,3)
   bz_struc%normal(:,13)=- bz_struc%bg(:,2) 
   bz_struc%normal(:,14)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2) + &
                           bz_struc%bg(:,3))

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 5, 6, 7, 8 /)
   bz_struc%indsur(:,3) = (/ 4, 9, 10, 11, 12 /)
   bz_struc%indsur(:,4) = (/ 4, 13, 14, 15, 16 /)
   bz_struc%indsur(:,5) = (/ 4, 17, 18, 19, 20 /)
   bz_struc%indsur(:,6) = (/ 4, 21, 22, 23, 24 /)
   bz_struc%indsur(:,7) = (/ 6, 1, 2, 8, 5, 18, 17 /)
   bz_struc%indsur(:,8) = (/ 6, 5, 6, 12, 9, 19, 18 /)
   bz_struc%indsur(:,9) = (/ 6, 9, 10, 16, 13, 20, 19 /)
   bz_struc%indsur(:,10) = (/ 6, 1, 17, 20, 13, 14, 4 /)
   bz_struc%indsur(:,11) = (/ 6, 3, 2, 8, 7, 22, 21 /)
   bz_struc%indsur(:,12) = (/ 6, 22, 7, 6, 12, 11, 23 /)
   bz_struc%indsur(:,13) = (/ 6, 11, 10, 16, 15, 24, 23 /)
   bz_struc%indsur(:,14) = (/ 6, 3, 4, 14, 15, 24, 21 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 


   bz_struc%letter_list(2)=' K '
   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%vertex_coord(:,2) &
                                      +bz_struc%vertex_coord(:,8))
   bz_struc%letter_list(3)=' W '
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,8)
   bz_struc%letter_list(4)=' U '
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%vertex_coord(:,5)+ &
                                        bz_struc%vertex_coord(:,8))
!
!   conventional bg  
!
!   bz_struc%letter_coord(:,5) = 0.5_DP*(bg(:,2) + bg(:,3))
!   bz_struc%letter_coord(:,6) = 0.5_DP*(bg(:,1) + bg(:,2) + bg(:,3))

   bz_struc%letter_list(5)=' X '
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%letter_list(6)=' L '
   bz_struc%letter_coord(:,6) = 0.5_DP * bz_struc%bg(:,2) 

   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(7) =' X1'
      bz_struc%letter_coord(:,7) =-0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,3))
      bz_struc%letter_list(8) =' W1'
      bz_struc%letter_coord(:,8) =bz_struc%vertex_coord(:,2)
      bz_struc%letter_list(9) =' U1'
      bz_struc%letter_coord(:,9) =0.5_DP*(bz_struc%vertex_coord(:,1)+&
                                           bz_struc%vertex_coord(:,2))
      bz_struc%letter_list(10) =' X3'
      bz_struc%letter_coord(:,10) =0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2))
      bz_struc%letter_list(11) =' W3'
      bz_struc%letter_coord(:,11) =bz_struc%vertex_coord(:,17)
      bz_struc%letter_list(12) =' L1'
      bz_struc%letter_coord(:,12) =-0.5_DP*(bz_struc%bg(:,1))
      bz_struc%letter_list(13) =' M '
      bz_struc%letter_coord(:,13) =0.5_DP*(bz_struc%bg(:,2) - bz_struc%bg(:,1))
   ENDIF

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_2

SUBROUTINE init_bz_3(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc

!
!  bcc bz
!
!  standard vectors (0,1,1) (1,0,1) (1,1,0)
!
!   bz_struc%normal(:,1)=bg(:,3)
!   bz_struc%normal(:,2)=bg(:,1) - bg(:,2)
!   bz_struc%normal(:,3)=-bg(:,3)
!   bz_struc%normal(:,4)=-(bg(:,1) - bg(:,2))
!   bz_struc%normal(:,5)=bg(:,2)
!   bz_struc%normal(:,6)=bg(:,1)
!   bz_struc%normal(:,7)=bg(:,1) - bg(:,3)
!   bz_struc%normal(:,8)=bg(:,2) - bg(:,3)
!   bz_struc%normal(:,9)=-(bg(:,1) - bg(:,3))
!   bz_struc%normal(:,10)=-(bg(:,2) - bg(:,3))
!   bz_struc%normal(:,11)=-bg(:,2)
!   bz_struc%normal(:,12)=-bg(:,1)
!
!  QE vectors (1,0,1) (-1,1,0) (0,-1,1)
!
   bz_struc%normal(:,1) = bz_struc%bg(:,1) - bz_struc%bg(:,3) 
   bz_struc%normal(:,2) = bz_struc%bg(:,2) 
   bz_struc%normal(:,3) =-(bz_struc%bg(:,1) - bz_struc%bg(:,3))
   bz_struc%normal(:,4) =-bz_struc%bg(:,2)
   bz_struc%normal(:,5) = bz_struc%bg(:,1)
   bz_struc%normal(:,6) = bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,7) = bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,8) = bz_struc%bg(:,3) 
   bz_struc%normal(:,9) =-(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%normal(:,10)=-bz_struc%bg(:,3)
   bz_struc%normal(:,11)=-bz_struc%bg(:,1)
   bz_struc%normal(:,12)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2))

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 3, 10 /)
   bz_struc%indsur(:,2) = (/ 4, 3, 4, 5, 11 /)
   bz_struc%indsur(:,3) = (/ 4, 5, 6, 7, 12 /)
   bz_struc%indsur(:,4) = (/ 4, 1, 8, 7, 13 /)
   bz_struc%indsur(:,5) = (/ 4, 2, 9, 8, 1 /)
   bz_struc%indsur(:,6) = (/ 4, 2, 9, 4, 3 /)
   bz_struc%indsur(:,7) = (/ 4, 4, 5, 6, 9 /)
   bz_struc%indsur(:,8) = (/ 4, 6, 7, 8, 9 /)
   bz_struc%indsur(:,9) = (/ 4, 1, 10, 14, 13 /)
   bz_struc%indsur(:,10) = (/ 4, 3, 10, 14, 11 /)
   bz_struc%indsur(:,11) = (/ 4, 5, 12, 14, 11 /)
   bz_struc%indsur(:,12) = (/ 4, 12, 7, 13, 14 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' N '
   bz_struc%letter_list(3)=' P '
   bz_struc%letter_list(4)=' H '
!
!   standard bg
!
!   bz_struc%letter_coord(:,2) = 0.5_DP*bg(:,3)
   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%bg(:,1) - bz_struc%bg(:,3))
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,4) = bz_struc%vertex_coord(:,3)

   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(5)=' H1'
      bz_struc%letter_coord(:,5) = bz_struc%vertex_coord(:,1)
   ENDIF

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_3

SUBROUTINE init_bz_4(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc

!
!  simple tetragonal bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1)
   bz_struc%normal(:,2)=bz_struc%bg(:,2)
   bz_struc%normal(:,3)=-bz_struc%bg(:,1)
   bz_struc%normal(:,4)=-bz_struc%bg(:,2)
   bz_struc%normal(:,5)=bz_struc%bg(:,3)
   bz_struc%normal(:,6)=-bz_struc%bg(:,3)

   bz_struc%indsur(:,1) = (/ 4, 1, 5, 8, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 1, 5, 6, 2 /)
   bz_struc%indsur(:,3) = (/ 4, 2, 6, 7, 3 /)
   bz_struc%indsur(:,4) = (/ 4, 3, 7, 8, 4 /)
   bz_struc%indsur(:,5) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,6) = (/ 4, 5, 6, 7, 8 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' M '
   bz_struc%letter_list(3)=' X '
   bz_struc%letter_list(4)=' Z '
   bz_struc%letter_list(5)=' A '
   bz_struc%letter_list(6)=' R '

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%letter_coord(:,3) = 0.5_DP*bz_struc%bg(:,2)
   bz_struc%letter_coord(:,4) = 0.5_DP*bz_struc%bg(:,3)
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                     + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,6) = 0.5_DP*(bz_struc%bg(:,2) + bz_struc%bg(:,3))

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_4

SUBROUTINE init_bz_5(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc

!
!  centered tetragonal (c<a) bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,2) - bz_struc%bg(:,3)
   bz_struc%normal(:,2)=-bz_struc%bg(:,1)
   bz_struc%normal(:,3)=-(bz_struc%bg(:,2) - bz_struc%bg(:,3))
   bz_struc%normal(:,4)=bz_struc%bg(:,1)
   bz_struc%normal(:,5)=bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,6)=bz_struc%bg(:,2)
   bz_struc%normal(:,7)=bz_struc%bg(:,3)
   bz_struc%normal(:,8)=bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,9)=-bz_struc%bg(:,3)
   bz_struc%normal(:,10)=-(bz_struc%bg(:,1) + bz_struc%bg(:,3))
   bz_struc%normal(:,11)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%normal(:,12)=-bz_struc%bg(:,2)

   bz_struc%indsur(:,1) = (/ 6, 1, 2, 3, 12, 11, 10 /)
   bz_struc%indsur(:,2) = (/ 6, 3, 4, 5, 14, 13, 12 /)
   bz_struc%indsur(:,3) = (/ 6, 5, 6, 7, 16, 15, 14 /)
   bz_struc%indsur(:,4) = (/ 6, 7, 8, 1, 10, 17, 16 /)
   bz_struc%indsur(:,5) = (/ 4, 1, 2, 9, 8 /)
   bz_struc%indsur(:,6) = (/ 4, 2, 3, 4, 9 /)
   bz_struc%indsur(:,7) = (/ 4, 4, 5, 6, 9 /)
   bz_struc%indsur(:,8) = (/ 4, 6, 7, 8, 9 /)
   bz_struc%indsur(:,9) = (/ 4, 10, 11, 18, 17 /)
   bz_struc%indsur(:,10) = (/ 4, 11, 12, 13, 18 /)
   bz_struc%indsur(:,11) = (/ 4, 13, 14, 15, 18 /)
   bz_struc%indsur(:,12) = (/ 4, 15, 16, 17, 18 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' M '
   bz_struc%letter_list(3)=' X '
   bz_struc%letter_list(4)=' Z1'
   bz_struc%letter_list(5)=' N '
   bz_struc%letter_list(6)=' Z '
   bz_struc%letter_list(7)=' P '

   bz_struc%letter_coord(:,2) = 0.5_DP * ( bz_struc%vertex_coord(:,1) + &
                                           bz_struc%vertex_coord(:,10) )
   bz_struc%letter_coord(:,3) = 0.5_DP*( bz_struc%bg(:,2)- bz_struc%bg(:,3))
   bz_struc%letter_coord(:,4) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%vertex_coord(:,9) &
                                       +bz_struc%vertex_coord(:,1))
   bz_struc%letter_coord(:,6) = bz_struc%vertex_coord(:,9)
   bz_struc%letter_coord(:,7) = bz_struc%vertex_coord(:,2)

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_5

SUBROUTINE init_bz_6(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  centered tetragonal (c>a) bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,2) - bz_struc%bg(:,3)
   bz_struc%normal(:,2)=-bz_struc%bg(:,1)
   bz_struc%normal(:,3)=-(bz_struc%bg(:,2) - bz_struc%bg(:,3))
   bz_struc%normal(:,4)=bz_struc%bg(:,1)
   bz_struc%normal(:,5)=bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,6)=bz_struc%bg(:,2)
   bz_struc%normal(:,7)=bz_struc%bg(:,3)
   bz_struc%normal(:,8)=bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,9)=-bz_struc%bg(:,3)
   bz_struc%normal(:,10)=-(bz_struc%bg(:,1) + bz_struc%bg(:,3))
   bz_struc%normal(:,11)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%normal(:,12)=-bz_struc%bg(:,2)
   bz_struc%normal(:,13)=bz_struc%bg(:,1) + bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,14)=-( bz_struc%bg(:,1) + bz_struc%bg(:,2) + &
                            bz_struc%bg(:,3) )

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 5, 6, 7, 8 /)
   bz_struc%indsur(:,3) = (/ 4, 9, 10, 11, 12 /)
   bz_struc%indsur(:,4) = (/ 4, 13, 14, 15, 16 /)
   bz_struc%indsur(:,5) = (/ 6, 1, 2, 17, 20, 14, 15 /)
   bz_struc%indsur(:,6) = (/ 6, 2, 3, 5, 6, 18, 17 /)
   bz_struc%indsur(:,7) = (/ 6, 6, 7, 9, 10, 19, 18 /)
   bz_struc%indsur(:,8) = (/ 6, 10, 11, 13, 14, 20, 19 /)
   bz_struc%indsur(:,9) = (/ 6, 1, 4, 21, 24, 16, 15 /)
   bz_struc%indsur(:,10) = (/ 6, 3, 4, 21, 22, 8, 5 /)
   bz_struc%indsur(:,11) = (/ 6, 7, 8, 22, 23, 12, 9  /)
   bz_struc%indsur(:,12) = (/ 6, 11, 12, 23, 24, 16, 13 /)
   bz_struc%indsur(:,13) = (/ 4, 17, 18, 19, 20 /)
   bz_struc%indsur(:,14) = (/ 4, 21, 22, 23, 24 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)='gS '
   bz_struc%letter_list(3)=' N '
   bz_struc%letter_list(4)='gS1'
   bz_struc%letter_list(5)=' Z '
   bz_struc%letter_list(6)=' Y1'
   bz_struc%letter_list(7)=' P '
   bz_struc%letter_list(8)=' X '
   bz_struc%letter_list(9)=' Y '
   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(2)=' S0'
      bz_struc%letter_list(4)=' S '
      bz_struc%letter_list(5)=' M '
      bz_struc%letter_list(6)=' G '
      bz_struc%letter_list(9)=' R '
   ENDIF

   bz_struc%letter_coord(:,2) = 0.5_DP*( bz_struc%vertex_coord(:,1) &
                                       + bz_struc%vertex_coord(:,15))
   bz_struc%letter_coord(:,3) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%letter_coord(:,4) = 0.5_DP*( bz_struc%vertex_coord(:,20) &
                                       + bz_struc%vertex_coord(:,17))
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                      + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,6) = bz_struc%vertex_coord(:,17)
   bz_struc%letter_coord(:,7) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,8) = 0.5_DP*(bz_struc%bg(:,2)-bz_struc%bg(:,3))
   bz_struc%letter_coord(:,9) = bz_struc%vertex_coord(:,1)

   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(10)=' N0'
      bz_struc%letter_coord(:,10)= -0.5_DP*bz_struc%bg(:,3)
      bz_struc%letter_list(11)=' M2'
      bz_struc%letter_coord(:,11) =-0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                      + bz_struc%bg(:,3))
      bz_struc%letter_list(12)=' S2'
      bz_struc%letter_coord(:,12) = 0.5_DP*( bz_struc%vertex_coord(:,24) &
                                       + bz_struc%vertex_coord(:,21))
      bz_struc%letter_list(13)=' G0'
      bz_struc%letter_coord(:,13) = bz_struc%vertex_coord(:,21)
      bz_struc%letter_list(14)=' T '
      bz_struc%letter_coord(:,14)=0.0_DP
      bz_struc%letter_coord(3,14)=0.5_DP*(bz_struc%bg(3,1) + bz_struc%bg(3,2))
      bz_struc%letter_list(15)=' T4'
      bz_struc%letter_coord(:,15)=0.0_DP
      bz_struc%letter_coord(3,15)=-0.5_DP*(bz_struc%bg(3,1) + bz_struc%bg(3,2))
   ENDIF

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_6

SUBROUTINE init_bz_7(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  simple orthorombic bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1)
   bz_struc%normal(:,2)=bz_struc%bg(:,2)
   bz_struc%normal(:,3)=-bz_struc%bg(:,1)
   bz_struc%normal(:,4)=-bz_struc%bg(:,2)
   bz_struc%normal(:,5)=bz_struc%bg(:,3)
   bz_struc%normal(:,6)=-bz_struc%bg(:,3)

   bz_struc%indsur(:,1) = (/ 4, 1, 5, 8, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 1, 5, 6, 2 /)
   bz_struc%indsur(:,3) = (/ 4, 2, 6, 7, 3 /)
   bz_struc%indsur(:,4) = (/ 4, 3, 7, 8, 4 /)
   bz_struc%indsur(:,5) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,6) = (/ 4, 5, 6, 7, 8 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' S '
   bz_struc%letter_list(4)=' Y '
   bz_struc%letter_list(5)=' U '
   bz_struc%letter_list(6)=' R '
   bz_struc%letter_list(7)=' T '
   bz_struc%letter_list(8)=' Z '

   bz_struc%letter_coord(:,2) = 0.5_DP*bz_struc%bg(:,1)
   bz_struc%letter_coord(:,3) = 0.5_DP*(bz_struc%bg(:,1)+bz_struc%bg(:,2))
   bz_struc%letter_coord(:,4) = 0.5_DP*bz_struc%bg(:,2)
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,1)+bz_struc%bg(:,3))
   bz_struc%letter_coord(:,6) = 0.5_DP*(bz_struc%bg(:,1)+bz_struc%bg(:,3)+ &
                                        bz_struc%bg(:,2))
   bz_struc%letter_coord(:,7) = 0.5_DP*(bz_struc%bg(:,2)+bz_struc%bg(:,3))
   bz_struc%letter_coord(:,8) = 0.5_DP*bz_struc%bg(:,3)

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_7

SUBROUTINE init_bz_8(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  face centered orthorombic (1/a^2 > 1/b^2 + 1/c^2) bz
!
   bz_struc%normal(:,1) = bz_struc%bg(:,1) + bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,2) = bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,3) = bz_struc%bg(:,3)
   bz_struc%normal(:,4) = - bz_struc%bg(:,2)
   bz_struc%normal(:,5) = - ( bz_struc%bg(:,2) + bz_struc%bg(:,3) )
   bz_struc%normal(:,6) = bz_struc%bg(:,1)
   bz_struc%normal(:,7) = bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,8) = bz_struc%bg(:,2) 
   bz_struc%normal(:,9) = - bz_struc%bg(:,1)
   bz_struc%normal(:,10)= - ( bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                            + bz_struc%bg(:,3))
   bz_struc%normal(:,11)= - bz_struc%bg(:,3)
   bz_struc%normal(:,12)=-( bz_struc%bg(:,1) + bz_struc%bg(:,3) )

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,2) = (/ 6, 3, 4, 13, 15, 6, 5 /)
   bz_struc%indsur(:,3) = (/ 4, 5, 6, 7, 8 /)
   bz_struc%indsur(:,4) = (/ 4, 7, 8, 10, 9 /)
   bz_struc%indsur(:,5) = (/ 6, 9, 10, 11, 12, 18, 17 /)
   bz_struc%indsur(:,6) = (/ 4, 1, 2, 11, 12 /)
   bz_struc%indsur(:,7) = (/ 6, 2, 3, 5, 8, 10, 11 /)
   bz_struc%indsur(:,8) = (/ 4, 1, 4, 13, 14 /)
   bz_struc%indsur(:,9) = (/ 4, 6, 7, 16, 15 /)
   bz_struc%indsur(:,10) = (/ 4, 7, 9, 17, 16 /)
   bz_struc%indsur(:,11) = (/ 4, 14, 18, 12, 1 /)
   bz_struc%indsur(:,12) = (/ 6, 14, 13, 15, 16, 17, 18 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' A '
   bz_struc%letter_list(3)=' X '
   bz_struc%letter_list(4)=' Z '
   bz_struc%letter_list(5)=' L '
   bz_struc%letter_list(6)=' X1'
   bz_struc%letter_list(7)=' A1'
   bz_struc%letter_list(8)=' T '
   bz_struc%letter_list(9)=' Y '

   IF (bz_struc%switch_b_c) THEN
      bz_struc%letter_list(2)=' A1'
      bz_struc%letter_list(4)=' Y '
      bz_struc%letter_list(7)=' A '
      bz_struc%letter_list(9)=' Z '
   ENDIF
   IF ( bz_struc%switch_a_b ) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(2)=' T '
         bz_struc%letter_list(3)=' Z '
         bz_struc%letter_list(4)=' Y '
         bz_struc%letter_list(7)=' A '
         bz_struc%letter_list(8)=' A1'
         bz_struc%letter_list(9)=' X '
      ELSE
         bz_struc%letter_list(2)=' T '
         bz_struc%letter_list(3)=' Y '
         bz_struc%letter_list(8)=' A '
         bz_struc%letter_list(9)=' X '
      END IF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(2)=' A '
         bz_struc%letter_list(3)=' Z '
         bz_struc%letter_list(4)=' X '
         bz_struc%letter_list(7)=' T '
         bz_struc%letter_list(8)=' A1'
         bz_struc%letter_list(9)=' Y '
      ELSE
         bz_struc%letter_list(2)=' A1 '
         bz_struc%letter_list(3)=' Y '
         bz_struc%letter_list(4)=' X '
         bz_struc%letter_list(7)=' T '
         bz_struc%letter_list(8)=' A '
         bz_struc%letter_list(9)=' Z '
      END IF
   ENDIF

   bz_struc%letter_coord(:,2) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%vertex_coord(:,1) &
                                      + bz_struc%vertex_coord(:,3))
   bz_struc%letter_coord(:,6) = bz_struc%vertex_coord(:,3)
   bz_struc%letter_coord(:,7) = bz_struc%vertex_coord(:,4)
   bz_struc%letter_coord(:,8) = 0.5_DP*(bz_struc%vertex_coord(:,3) &
                                      + bz_struc%vertex_coord(:,5))
   bz_struc%letter_coord(:,9) = 0.5_DP*(bz_struc%vertex_coord(:,4) &
                                      + bz_struc%vertex_coord(:,6))

   CALL find_axis_coordinates(bz_struc)

   CALL adjust_orthorombic(bz_struc)

END SUBROUTINE init_bz_8

SUBROUTINE init_bz_9(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  face centered orthorombic (1/a^2 < 1/b^2 + 1/c^2) bz case 2
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,2)=bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,3)=-(bz_struc%bg(:,1)+bz_struc%bg(:,2))
   bz_struc%normal(:,4)=-(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%normal(:,5)=bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,6)=-(bz_struc%bg(:,1)+bz_struc%bg(:,3))
   bz_struc%normal(:,7)=bz_struc%bg(:,1) + bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,8)=bz_struc%bg(:,3) 
   bz_struc%normal(:,9)=- bz_struc%bg(:,2)
   bz_struc%normal(:,10)=bz_struc%bg(:,1) 
   bz_struc%normal(:,11)=bz_struc%bg(:,2)
   bz_struc%normal(:,12)=-bz_struc%bg(:,1)
   bz_struc%normal(:,13)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2) + &
                           bz_struc%bg(:,3))
   bz_struc%normal(:,14)=-bz_struc%bg(:,3)

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 5, 6, 7, 8 /)
   bz_struc%indsur(:,3) = (/ 4, 9, 10, 11, 12 /)
   bz_struc%indsur(:,4) = (/ 4, 13, 14, 15, 16 /)
   bz_struc%indsur(:,5) = (/ 4, 17, 18, 19, 20 /)
   bz_struc%indsur(:,6) = (/ 4, 21, 22, 23, 24 /)
   bz_struc%indsur(:,7) = (/ 6, 1, 2, 8, 5, 18, 17 /)
   bz_struc%indsur(:,8) = (/ 6, 5, 6, 12, 9, 19, 18 /)
   bz_struc%indsur(:,9) = (/ 6, 9, 10, 16, 13, 20, 19 /)
   bz_struc%indsur(:,10) = (/ 6, 1, 17, 20, 13, 14, 4 /)
   bz_struc%indsur(:,11) = (/ 6, 3, 2, 8, 7, 22, 21 /)
   bz_struc%indsur(:,12) = (/ 6, 22, 7, 6, 12, 11, 23 /)
   bz_struc%indsur(:,13) = (/ 6, 11, 10, 16, 15, 24, 23 /)
   bz_struc%indsur(:,14) = (/ 6, 3, 4, 14, 15, 24, 21 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' D '
   bz_struc%letter_list(4)=' C '
   bz_struc%letter_list(5)=' Y '
   bz_struc%letter_list(6)=' H1'
   bz_struc%letter_list(7)=' L '
   bz_struc%letter_list(8)=' H '
   bz_struc%letter_list(9)=' C1'
   bz_struc%letter_list(10)=' D1'
   bz_struc%letter_list(11)=' Z '

   IF (bz_struc%switch_b_c) THEN
      bz_struc%letter_list(3)=' H1'
      bz_struc%letter_list(4)=' C1'
      bz_struc%letter_list(5)=' Z '
      bz_struc%letter_list(6)=' D '
      bz_struc%letter_list(8)=' D1'
      bz_struc%letter_list(9)=' C '
      bz_struc%letter_list(10)=' H '
      bz_struc%letter_list(11)=' Y '
   ENDIF

   IF (bz_struc%switch_a_b) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(2)=' Z '
         bz_struc%letter_list(3)=' C1'
         bz_struc%letter_list(4)=' H1'
         bz_struc%letter_list(5)=' X '
         bz_struc%letter_list(6)=' D1'
         bz_struc%letter_list(8)=' D '
         bz_struc%letter_list(9)=' H '
         bz_struc%letter_list(10)=' C '
         bz_struc%letter_list(11)=' Y '
      ELSE
         bz_struc%letter_list(2)=' Y '
         bz_struc%letter_list(3)=' C '
         bz_struc%letter_list(4)=' D '
         bz_struc%letter_list(5)=' X '
         bz_struc%letter_list(6)=' H '
         bz_struc%letter_list(8)=' H1'
         bz_struc%letter_list(9)=' D1'
         bz_struc%letter_list(10)=' C1'
      ENDIF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(2)=' Z '
         bz_struc%letter_list(3)=' D1'
         bz_struc%letter_list(4)=' H '
         bz_struc%letter_list(5)=' Y '
         bz_struc%letter_list(6)=' C1'
         bz_struc%letter_list(8)=' C '
         bz_struc%letter_list(9)=' H1'
         bz_struc%letter_list(10)=' D '
         bz_struc%letter_list(11)=' X '
      ELSE
         bz_struc%letter_list(2)=' Y '
         bz_struc%letter_list(3)=' H '
         bz_struc%letter_list(4)=' D1'
         bz_struc%letter_list(6)=' C '
         bz_struc%letter_list(5)=' Z '
         bz_struc%letter_list(8)=' C1'
         bz_struc%letter_list(9)=' D '
         bz_struc%letter_list(10)=' H1'
         bz_struc%letter_list(11)=' X '
      ENDIF
   ENDIF

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,4) = bz_struc%vertex_coord(:,8)
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,6) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,7) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                              + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,8) = bz_struc%vertex_coord(:,5)
   bz_struc%letter_coord(:,9) = bz_struc%vertex_coord(:,17)
   bz_struc%letter_coord(:,10) = bz_struc%vertex_coord(:,18)
   bz_struc%letter_coord(:,11) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,3))


   CALL find_axis_coordinates(bz_struc)

   CALL adjust_orthorombic(bz_struc)

END SUBROUTINE init_bz_9

SUBROUTINE init_bz_10(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc

!
!  face centered orthorombic (1/a^2 = 1/b^2 + 1/c^2) bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1) + bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,2)=bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,3)=bz_struc%bg(:,3)
   bz_struc%normal(:,4)=-bz_struc%bg(:,2) 
   bz_struc%normal(:,5)=-(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%normal(:,6)=bz_struc%bg(:,1)
   bz_struc%normal(:,7)=bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,8)=bz_struc%bg(:,2) 
   bz_struc%normal(:,9)=- bz_struc%bg(:,1)
   bz_struc%normal(:,10)= -( bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                           + bz_struc%bg(:,3) )
   bz_struc%normal(:,11)=-bz_struc%bg(:,3)
   bz_struc%normal(:,12)=-(bz_struc%bg(:,1) + bz_struc%bg(:,3))

   bz_struc%indsur(:,1)  = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,2)  = (/ 4, 2, 3, 5, 11 /)
   bz_struc%indsur(:,3)  = (/ 4, 3, 5, 6, 7 /)
   bz_struc%indsur(:,4)  = (/ 4, 6, 7, 9, 8 /)
   bz_struc%indsur(:,5)  = (/ 4, 8, 9, 10, 12 /)
   bz_struc%indsur(:,6)  = (/ 4, 1, 4, 9, 10 /)
   bz_struc%indsur(:,7)  = (/ 4, 3, 4, 9, 7 /)
   bz_struc%indsur(:,8)  = (/ 4, 1, 2, 11, 13 /)
   bz_struc%indsur(:,9)  = (/ 4, 11, 5, 6, 14 /)
   bz_struc%indsur(:,10) = (/ 4, 6, 14, 12, 8 /)
   bz_struc%indsur(:,11) = (/ 4, 1, 10, 12, 13 /)
   bz_struc%indsur(:,12) = (/ 4, 11, 14, 12, 13 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' A '
   bz_struc%letter_list(4)=' L '
   bz_struc%letter_list(5)=' Z '
   bz_struc%letter_list(6)=' T '
   bz_struc%letter_list(7)=' A1'
   bz_struc%letter_list(8)=' Y '

   IF (bz_struc%switch_b_c) THEN
      bz_struc%letter_list(3)=' A1 '
      bz_struc%letter_list(5)=' Y '
      bz_struc%letter_list(7)=' A '
      bz_struc%letter_list(8)=' Z '
   ENDIF

   IF (bz_struc%switch_a_b) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(2)=' Z '
         bz_struc%letter_list(3)=' T '
         bz_struc%letter_list(5)=' Y '
         bz_struc%letter_list(6)=' A1'
         bz_struc%letter_list(7)=' A '
         bz_struc%letter_list(8)=' X '
      ELSE
         bz_struc%letter_list(3)=' T '
         bz_struc%letter_list(6)=' A '
         bz_struc%letter_list(8)=' X '
         bz_struc%letter_list(2)=' Y '
      ENDIF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(2)=' Z '
         bz_struc%letter_list(3)=' A '
         bz_struc%letter_list(5)=' X '
         bz_struc%letter_list(6)=' A1'
         bz_struc%letter_list(7)=' T'
         bz_struc%letter_list(8)=' Y '
      ELSE
         bz_struc%letter_list(2)=' Y '
         bz_struc%letter_list(3)=' A1'
         bz_struc%letter_list(5)=' X '
         bz_struc%letter_list(6)=' A '
         bz_struc%letter_list(7)=' T '
         bz_struc%letter_list(8)=' Z '
      ENDIF
   ENDIF


   bz_struc%letter_coord(:,2) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,4)
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                      + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,6) = bz_struc%vertex_coord(:,3)
   bz_struc%letter_coord(:,7) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,8) = 0.5_DP*(bz_struc%bg(:,2) + bz_struc%bg(:,3))

   CALL find_axis_coordinates(bz_struc)
   CALL adjust_orthorombic(bz_struc)

END SUBROUTINE init_bz_10

SUBROUTINE init_bz_11(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  body centered orthorombic bz
!
!  this is ok with the standard definition of vectors for bcc (see the 
!  corresponding bz)
!
!   bz_struc%normal(:,1)=bg(:,2) - bg(:,3)
!   bz_struc%normal(:,2)=-bg(:,1)
!   bz_struc%normal(:,3)=-(bg(:,2) - bg(:,3))
!   bz_struc%normal(:,4)=bg(:,1)
!   bz_struc%normal(:,5)=bg(:,1) + bg(:,2)
!   bz_struc%normal(:,6)=bg(:,2)
!   bz_struc%normal(:,7)=bg(:,3)
!   bz_struc%normal(:,8)=bg(:,1) + bg(:,3)
!   bz_struc%normal(:,9)=-bg(:,3)
!   bz_struc%normal(:,10)=-(bg(:,1) + bg(:,3))
!   bz_struc%normal(:,11)=-(bg(:,1) + bg(:,2))
!   bz_struc%normal(:,12)=-bg(:,2)
!   bz_struc%normal(:,13)=bg(:,1) + bg(:,2) + bg(:,3)
!   bz_struc%normal(:,14)=-( bg(:,1) + bg(:,2) + bg(:,3) )

   bz_struc%normal(:,1) = bz_struc%bg(:,1) - bz_struc%bg(:,3)
   bz_struc%normal(:,2) = bz_struc%bg(:,2)
   bz_struc%normal(:,3) =-(bz_struc%bg(:,1) - bz_struc%bg(:,3))
   bz_struc%normal(:,4) =-bz_struc%bg(:,2)
   bz_struc%normal(:,5) = bz_struc%bg(:,1)
   bz_struc%normal(:,6) = bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,7) = bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,8) = bz_struc%bg(:,3)
   bz_struc%normal(:,9) =-(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   bz_struc%normal(:,10)=-bz_struc%bg(:,3)
   bz_struc%normal(:,11)=-bz_struc%bg(:,1)
   bz_struc%normal(:,12)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%normal(:,13)= bz_struc%bg(:,1) + bz_struc%bg(:,2) + &
                          bz_struc%bg(:,3)
   bz_struc%normal(:,14)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2) + &
                           bz_struc%bg(:,3))

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 3, 4 /)
   bz_struc%indsur(:,2) = (/ 4, 5, 6, 7, 8 /)
   bz_struc%indsur(:,3) = (/ 4, 9, 10, 11, 12 /)
   bz_struc%indsur(:,4) = (/ 4, 13, 14, 15, 16 /)
   bz_struc%indsur(:,5) = (/ 6, 1, 2, 17, 20, 14, 15 /)
   bz_struc%indsur(:,6) = (/ 6, 2, 3, 5, 6, 18, 17 /)
   bz_struc%indsur(:,7) = (/ 6, 6, 7, 9, 10, 19, 18 /)
   bz_struc%indsur(:,8) = (/ 6, 10, 11, 13, 14, 20, 19 /)
   bz_struc%indsur(:,9) = (/ 6, 1, 4, 21, 24, 16, 15 /)
   bz_struc%indsur(:,10) = (/ 6, 3, 4, 21, 22, 8, 5 /)
   bz_struc%indsur(:,11) = (/ 6, 7, 8, 22, 23, 12, 9  /)
   bz_struc%indsur(:,12) = (/ 6, 11, 12, 23, 24, 16, 13 /)
   bz_struc%indsur(:,13) = (/ 4, 17, 18, 19, 20 /)
   bz_struc%indsur(:,14) = (/ 4, 21, 22, 23, 24 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' L '
   bz_struc%letter_list(4)=' T '
   bz_struc%letter_list(5)=' L1'
   bz_struc%letter_list(6)=' Y '
   bz_struc%letter_list(7)=' R '
   bz_struc%letter_list(8)=' W '
   bz_struc%letter_list(9)=' S '
   bz_struc%letter_list(10)=' X1'
   bz_struc%letter_list(11)=' L2'
   bz_struc%letter_list(12)=' Y1'
   bz_struc%letter_list(13)=' Z '

   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(2)='gS0'
      bz_struc%letter_list(4)=' R '
      bz_struc%letter_list(6)='gL0'
      bz_struc%letter_list(7)=' T '
      bz_struc%letter_list(10)=' F0'
      bz_struc%letter_list(12)=' G0'
      bz_struc%letter_list(13)=' G '
   ENDIF

   IF (bz_struc%switch_b_c) THEN
      bz_struc%letter_list(3)=' L '
      bz_struc%letter_list(4)=' R '
      bz_struc%letter_list(5)=' L2'
      bz_struc%letter_list(7)=' T '
      bz_struc%letter_list(8)=' W '
      bz_struc%letter_list(9)=' S '
      bz_struc%letter_list(10)=' X1'
      bz_struc%letter_list(11)=' L1'
      bz_struc%letter_list(12)=' Z1'
      bz_struc%letter_list(6)=' Z '
      bz_struc%letter_list(13)=' Y '
      IF (bz_struc%letter_type=='BI') THEN
         bz_struc%letter_list(6)='gL0'
         bz_struc%letter_list(13)=' G '
      ENDIF
   ENDIF

   IF (bz_struc%switch_a_b) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(3)=' L2 '
         bz_struc%letter_list(4)=' R '
         bz_struc%letter_list(5)=' L'
         bz_struc%letter_list(7)=' S '
         bz_struc%letter_list(9)=' T '
         bz_struc%letter_list(10)=' Z1'
         bz_struc%letter_list(12)=' X1'
         bz_struc%letter_list(2)=' Z '
         bz_struc%letter_list(6)=' X '
         bz_struc%letter_list(11)=' L1'
         bz_struc%letter_list(13)=' Y '
         IF (bz_struc%letter_type=='BI') THEN
            bz_struc%letter_list(2)='gL0'
            bz_struc%letter_list(6)='gS0'
            bz_struc%letter_list(13)=' G '
         ENDIF
      ELSE
         bz_struc%letter_list(2)=' Y '
         bz_struc%letter_list(3)=' L1'
         bz_struc%letter_list(5)=' L '
         bz_struc%letter_list(6)=' X '
         bz_struc%letter_list(13)=' Z '
         bz_struc%letter_list(7)=' S '
         bz_struc%letter_list(9)=' R '
         bz_struc%letter_list(10)=' Y1'
         bz_struc%letter_list(11)=' L2'
         bz_struc%letter_list(12)=' X1'
         IF (bz_struc%letter_type=='BI') THEN
            bz_struc%letter_list(2)=' G '
            bz_struc%letter_list(6)='gS0'
            bz_struc%letter_list(13)='gL0'
         ENDIF
      END IF
   ELSEIF (bz_struc%rotate_a_b_c) THEN
      IF (bz_struc%switch_b_c) THEN
         bz_struc%letter_list(3)=' L2'
         bz_struc%letter_list(4)=' S '
         bz_struc%letter_list(5)=' L1'
         bz_struc%letter_list(7)=' R '
         bz_struc%letter_list(8)=' W '
         bz_struc%letter_list(9)=' T '
         bz_struc%letter_list(10)=' Z1'
         bz_struc%letter_list(11)=' L '
         bz_struc%letter_list(12)=' Y1'
         bz_struc%letter_list(2)=' Z '
         bz_struc%letter_list(6)=' Y '
         bz_struc%letter_list(13)=' X '
         IF (bz_struc%letter_type=='BI') THEN
            bz_struc%letter_list(2)='gL0'
            bz_struc%letter_list(6)=' G '
            bz_struc%letter_list(13)='gS0'
         ENDIF
      ELSE
         bz_struc%letter_list(2)=' Y '
         bz_struc%letter_list(3)=' L1 '
         bz_struc%letter_list(4)=' S '
         bz_struc%letter_list(5)=' L2'
         bz_struc%letter_list(6)=' Z '
         bz_struc%letter_list(7)=' T '
         bz_struc%letter_list(9)=' R '
         bz_struc%letter_list(10)=' Y1'
         bz_struc%letter_list(11)=' L '
         bz_struc%letter_list(12)=' Z1'
         bz_struc%letter_list(13)=' X '
         IF (bz_struc%letter_type=='BI') THEN
            bz_struc%letter_list(2)=' G '
            bz_struc%letter_list(6)='gL0'
            bz_struc%letter_list(13)='gS0'
         ENDIF
      ENDIF
   END IF

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%vertex_coord(:,1)+ &
                                        bz_struc%vertex_coord(:,15)) 
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,1)
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%bg(:,1)-bz_struc%bg(:,3))
   bz_struc%letter_coord(:,5) = bz_struc%vertex_coord(:,3)
   bz_struc%letter_coord(:,6) = 0.5_DP*(bz_struc%vertex_coord(:,3)+ &
                                        bz_struc%vertex_coord(:,5)) 
   bz_struc%letter_coord(:,7) = 0.5_DP*bz_struc%bg(:,1)
   bz_struc%letter_coord(:,8) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,9) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2))
   bz_struc%letter_coord(:,10) = 0.5_DP*(bz_struc%vertex_coord(:,17)+ &
                                         bz_struc%vertex_coord(:,20)) 
   bz_struc%letter_coord(:,11) = bz_struc%vertex_coord(:,17)
   bz_struc%letter_coord(:,12) = 0.5_DP*(bz_struc%vertex_coord(:,17)+ &
                                         bz_struc%vertex_coord(:,18)) 
   bz_struc%letter_coord(:,13) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                       + bz_struc%bg(:,3))


   CALL find_axis_coordinates(bz_struc)
   CALL adjust_orthorombic(bz_struc)

END SUBROUTINE init_bz_11

SUBROUTINE init_bz_12(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc

!
!  one face centered orthorombic bz
!
   IF (bz_struc%ibrav==9) THEN
      bz_struc%normal(:,1) =  bz_struc%bg(:,1)
      bz_struc%normal(:,2) =  bz_struc%bg(:,1) + bz_struc%bg(:,2)
      bz_struc%normal(:,3) =  bz_struc%bg(:,2) 
      bz_struc%normal(:,4) = -bz_struc%bg(:,1)
      bz_struc%normal(:,5) = -( bz_struc%bg(:,1) + bz_struc%bg(:,2) )
      bz_struc%normal(:,6) = -bz_struc%bg(:,2)
   ELSE
      bz_struc%normal(:,1) =  bz_struc%bg(:,2)
      bz_struc%normal(:,2) =  bz_struc%bg(:,2) - bz_struc%bg(:,1)
      bz_struc%normal(:,3) = -bz_struc%bg(:,1) 
      bz_struc%normal(:,4) = -bz_struc%bg(:,2)
      bz_struc%normal(:,5) =  bz_struc%bg(:,1) - bz_struc%bg(:,2) 
      bz_struc%normal(:,6) =  bz_struc%bg(:,1)
   ENDIF
   bz_struc%normal(:,7) =  bz_struc%bg(:,3)
   bz_struc%normal(:,8) = -bz_struc%bg(:,3)

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 8, 7 /)
   bz_struc%indsur(:,2) = (/ 4, 2, 3, 9, 8 /)
   bz_struc%indsur(:,3) = (/ 4, 3, 4, 10, 9 /)
   bz_struc%indsur(:,4) = (/ 4, 4, 5, 11, 10 /)
   bz_struc%indsur(:,5) = (/ 4, 5, 6, 12, 11 /)
   bz_struc%indsur(:,6) = (/ 4, 6, 1, 7, 12 /)
   bz_struc%indsur(:,7) = (/ 6, 1, 2, 3, 4, 5, 6 /)
   bz_struc%indsur(:,8) = (/ 6, 7, 8, 9, 10, 11, 12 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' X '
   bz_struc%letter_list(3)=' S '
   bz_struc%letter_list(4)=' X1'
   bz_struc%letter_list(5)=' Y '
   bz_struc%letter_list(6)=' T '
   bz_struc%letter_list(7)=' A1'
   bz_struc%letter_list(8)=' R '
   bz_struc%letter_list(9)=' A '
   bz_struc%letter_list(10)=' Z '

   IF (bz_struc%switch_a_b) THEN
      bz_struc%letter_list(2)=' Y '
      bz_struc%letter_list(3)=' S '
      bz_struc%letter_list(4)=' Y1'
      bz_struc%letter_list(5)=' X '
      bz_struc%letter_list(6)=' A '
      bz_struc%letter_list(7)=' A1'
      bz_struc%letter_list(8)=' R '
      bz_struc%letter_list(9)=' T '
      bz_struc%letter_list(10)=' Z '
   ENDIF

   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(2)='gD0'
      bz_struc%letter_list(9)=' B0'
   ENDIF

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%vertex_coord(:,1)+ &
                                        bz_struc%vertex_coord(:,7))
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%vertex_coord(:,2)+ &
                                        bz_struc%vertex_coord(:,8))
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%vertex_coord(:,2)+ &
                                        bz_struc%vertex_coord(:,9))
   bz_struc%letter_coord(:,6) = 0.5_DP*(bz_struc%vertex_coord(:,2) + &
                               bz_struc%vertex_coord(:,9) + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,7) = bz_struc%vertex_coord(:,2) 
   bz_struc%letter_coord(:,9) = 0.5_DP*(bz_struc%vertex_coord(:,1) + &
                                        bz_struc%vertex_coord(:,7) + &
                                        bz_struc%bg(:,3))
   bz_struc%letter_coord(:,10) = 0.5_DP*bz_struc%bg(:,3)
   IF (bz_struc%ibrav==9) THEN
      bz_struc%letter_coord(:,3) = 0.5_DP*bz_struc%bg(:,1)
      bz_struc%letter_coord(:,8) = 0.5_DP*(bz_struc%bg(:,1) + bz_struc%bg(:,3))
   ELSE
      bz_struc%letter_coord(:,3) = 0.5_DP*bz_struc%bg(:,2)
      bz_struc%letter_coord(:,8) = 0.5_DP*(bz_struc%bg(:,2) + bz_struc%bg(:,3))
   ENDIF

   CALL find_axis_coordinates(bz_struc)

   CALL adjust_one_face_centered_orthorombic(bz_struc)

END SUBROUTINE init_bz_12

SUBROUTINE init_bz_13(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  hexagonal
!
   bz_struc%normal(:,1) = bz_struc%bg(:,1)
   bz_struc%normal(:,2) = bz_struc%bg(:,2)
   bz_struc%normal(:,3) = bz_struc%bg(:,2) - bz_struc%bg(:,1)
   bz_struc%normal(:,4) = -bz_struc%bg(:,1)
   bz_struc%normal(:,5) = -bz_struc%bg(:,2)
   bz_struc%normal(:,6) = -(bz_struc%bg(:,2) - bz_struc%bg(:,1))
   bz_struc%normal(:,7) = bz_struc%bg(:,3)
   bz_struc%normal(:,8) = -bz_struc%bg(:,3)

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 8, 7 /)
   bz_struc%indsur(:,2) = (/ 4, 2, 3, 9, 8 /)
   bz_struc%indsur(:,3) = (/ 4, 3, 4, 10, 9 /)
   bz_struc%indsur(:,4) = (/ 4, 4, 5, 11, 10 /)
   bz_struc%indsur(:,5) = (/ 4, 5, 6, 12, 11 /)
   bz_struc%indsur(:,6) = (/ 4, 6, 1, 7, 12 /)
   bz_struc%indsur(:,7) = (/ 6, 1, 2, 3, 4, 5, 6 /)
   bz_struc%indsur(:,8) = (/ 6, 7, 8, 9, 10, 11, 12 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' K '
   bz_struc%letter_list(3)=' M '
   bz_struc%letter_list(4)=' A '
   bz_struc%letter_list(5)=' L '
   bz_struc%letter_list(6)=' H '

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%vertex_coord(:,1) + &
                                        bz_struc%vertex_coord(:,7))
   bz_struc%letter_coord(:,3) = 0.5_DP*bz_struc%bg(:,1)
   bz_struc%letter_coord(:,4) = 0.5_DP*bz_struc%bg(:,3)
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%bg(:,1)+bz_struc%bg(:,3))
   bz_struc%letter_coord(:,6) = bz_struc%vertex_coord(:,1) 

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_13

SUBROUTINE init_bz_14(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
!
!  trigonal alpha < 90 bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1) + bz_struc%bg(:,2)
   bz_struc%normal(:,2)=bz_struc%bg(:,2)
   bz_struc%normal(:,3)=bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,4)=bz_struc%bg(:,3)
   bz_struc%normal(:,5)=bz_struc%bg(:,1) + bz_struc%bg(:,3)
   bz_struc%normal(:,6)=bz_struc%bg(:,1) 
   bz_struc%normal(:,7)=bz_struc%bg(:,1) + bz_struc%bg(:,2) + bz_struc%bg(:,3)
   bz_struc%normal(:,8)= -( bz_struc%bg(:,1) + bz_struc%bg(:,3) )
   bz_struc%normal(:,9)= -bz_struc%bg(:,1)
   bz_struc%normal(:,10)=-(bz_struc%bg(:,1) + bz_struc%bg(:,2) )
   bz_struc%normal(:,11)=-bz_struc%bg(:,2)
   bz_struc%normal(:,12)=-( bz_struc%bg(:,2) + bz_struc%bg(:,3) )
   bz_struc%normal(:,13)=-bz_struc%bg(:,3) 
   bz_struc%normal(:,14)=-( bz_struc%bg(:,1) + bz_struc%bg(:,2) + &
                            bz_struc%bg(:,3) )

   bz_struc%indsur(:,1) = (/ 4, 1, 6, 13, 14 /)
   bz_struc%indsur(:,2) = (/ 6, 1, 2, 15, 20, 19, 14 /)
   bz_struc%indsur(:,3) = (/ 4, 2, 3, 16, 15 /)
   bz_struc%indsur(:,4) = (/ 6, 3, 4, 17, 22, 21, 16 /)
   bz_struc%indsur(:,5) = (/ 4, 4, 5, 18, 17 /)
   bz_struc%indsur(:,6) = (/ 6, 5, 6, 13, 24, 23, 18 /)
   bz_struc%indsur(:,7) = (/ 6, 1, 2, 3, 4, 5, 6 /)
   bz_struc%indsur(:,8) = (/ 4, 7, 8, 20, 19/)
   bz_struc%indsur(:,9) = (/ 6, 8, 9, 21, 16, 15, 20 /)
   bz_struc%indsur(:,10) = (/ 4, 9, 10, 22, 21 /)
   bz_struc%indsur(:,11) = (/ 6, 17, 18, 23, 11, 10, 22 /)
   bz_struc%indsur(:,12) = (/ 4, 11, 12, 24, 23 /)
   bz_struc%indsur(:,13) = (/ 6, 12, 24, 13, 14, 19, 7  /)
   bz_struc%indsur(:,14) = (/ 6, 7, 8, 9, 10, 11, 12  /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' L1'
   bz_struc%letter_list(3)=' L '
   bz_struc%letter_list(4)=' X '
   bz_struc%letter_list(5)=' Q '
   bz_struc%letter_list(6)=' F '
   bz_struc%letter_list(7)=' P '
   bz_struc%letter_list(8)=' B '
   bz_struc%letter_list(9)=' P1'
   bz_struc%letter_list(10)=' Z '
   bz_struc%letter_list(11)=' B1'
   bz_struc%letter_list(12)=' P2'

   bz_struc%letter_coord(:,2) = -0.5_DP*bz_struc%bg(:,3)
   bz_struc%letter_coord(:,3) =  0.5_DP*bz_struc%bg(:,1)
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%vertex_coord(:,13)+ &
                                        bz_struc%vertex_coord(:,24))
   bz_struc%letter_coord(:,5) = 0.5_DP*(bz_struc%vertex_coord(:,13)+ &
                                        bz_struc%vertex_coord(:,6))
   bz_struc%letter_coord(:,6) = 0.5_DP*( bz_struc%bg(:,1) + bz_struc%bg(:,2) )
   bz_struc%letter_coord(:,7) = 0.5_DP*(bz_struc%vertex_coord(:,5)+ &
                                        bz_struc%vertex_coord(:,6))
   bz_struc%letter_coord(:,8) = bz_struc%vertex_coord(:,6)
   bz_struc%letter_coord(:,9) = 0.5_DP*(bz_struc%vertex_coord(:,6)+ &
                                        bz_struc%vertex_coord(:,1))
   bz_struc%letter_coord(:,10) = 0.5_DP*( bz_struc%bg(:,1) + bz_struc%bg(:,2) &
                                        + bz_struc%bg(:,3))
   bz_struc%letter_coord(:,11) = bz_struc%vertex_coord(:,13)
   bz_struc%letter_coord(:,12) = 0.5_DP*(bz_struc%vertex_coord(:,13)+ &
                                         bz_struc%vertex_coord(:,14))
   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_14

SUBROUTINE init_bz_15(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc

!
!  trigonal alpha > 90 bz
!
   bz_struc%normal(:,1)=bz_struc%bg(:,1) - bz_struc%bg(:,3)
   bz_struc%normal(:,2)=-(bz_struc%bg(:,3) - bz_struc%bg(:,2))
   bz_struc%normal(:,3)=-(bz_struc%bg(:,1) - bz_struc%bg(:,2))
   bz_struc%normal(:,4)=-(bz_struc%bg(:,1) - bz_struc%bg(:,3))
   bz_struc%normal(:,5)=bz_struc%bg(:,3) - bz_struc%bg(:,2)
   bz_struc%normal(:,6)=bz_struc%bg(:,1) - bz_struc%bg(:,2)
   bz_struc%normal(:,7)=bz_struc%bg(:,1) 
   bz_struc%normal(:,8)=bz_struc%bg(:,2) 
   bz_struc%normal(:,9)=bz_struc%bg(:,3)
   bz_struc%normal(:,10)=-bz_struc%bg(:,3) 
   bz_struc%normal(:,11)=-bz_struc%bg(:,1)
   bz_struc%normal(:,12)=-bz_struc%bg(:,2)

   bz_struc%indsur(:,1) = (/ 4, 1, 2, 9, 8 /)
   bz_struc%indsur(:,2) = (/ 4, 2, 3, 10, 9 /)
   bz_struc%indsur(:,3) = (/ 4, 3, 4, 11, 10 /)
   bz_struc%indsur(:,4) = (/ 4, 4, 5, 12, 11 /)
   bz_struc%indsur(:,5) = (/ 4, 5, 6, 13, 12 /)
   bz_struc%indsur(:,6) = (/ 4, 6, 1, 8, 13 /)
   bz_struc%indsur(:,7) = (/ 4, 1, 2, 7, 6 /)
   bz_struc%indsur(:,8) = (/ 4, 2, 3, 4, 7 /)
   bz_struc%indsur(:,9) = (/ 4, 4, 5, 6, 7 /)
   bz_struc%indsur(:,10) = (/ 4, 8, 9, 10, 14 /)
   bz_struc%indsur(:,11) = (/ 4, 10, 11, 12, 14 /)
   bz_struc%indsur(:,12) = (/ 4, 12, 13, 8, 14 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_list(2)=' P1'
   bz_struc%letter_list(3)=' Q1'
   bz_struc%letter_list(4)=' P '
   bz_struc%letter_list(5)=' Z '
   bz_struc%letter_list(6)=' F '
   bz_struc%letter_list(7)=' L '
   bz_struc%letter_list(8)=' Q '

   IF (bz_struc%letter_type=='BI') THEN
      bz_struc%letter_list(3)=' R0'
      bz_struc%letter_list(5)=' T '
      bz_struc%letter_list(6)=' FA'
      bz_struc%letter_list(8)=' P2'
   ENDIF

   bz_struc%letter_coord(:,2) = 0.5_DP*(bz_struc%vertex_coord(:,10)+ &
                                        bz_struc%vertex_coord(:,3))
   bz_struc%letter_coord(:,3) = bz_struc%vertex_coord(:,3)
   bz_struc%letter_coord(:,4) = 0.5_DP*(bz_struc%vertex_coord(:,2)+ &
                                        bz_struc%vertex_coord(:,9))
   bz_struc%letter_coord(:,5) = bz_struc%vertex_coord(:,2)
   bz_struc%letter_coord(:,6) = 0.5_DP * ( bz_struc%vertex_coord(:,2)+ &
                                           bz_struc%vertex_coord(:,10))
   bz_struc%letter_coord(:,7) = 0.5_DP * bz_struc%bg(:,2)
   bz_struc%letter_coord(:,8) = bz_struc%vertex_coord(:,7)

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_15

SUBROUTINE init_bz_16(bz_struc)
  IMPLICIT NONE
  TYPE(bz), INTENT(INOUT) :: bz_struc
  INTEGER :: n1(6), n2(6), i, idir, idir1
!
!  Simple monoclinic lattice
!
   IF (bz_struc%ibrav==12) THEN
      idir=2
      idir1=3
   ELSE
      idir=3
      idir1=2
   ENDIF
   CALL find_n1n2_monoclinic(n1, n2, idir, bz_struc)

   DO i=1,6 
      bz_struc%normal(:,i) = n1(i) * bz_struc%bg(:,1) + n2(i) *&
                             bz_struc%bg(:,idir)
   END DO
   bz_struc%normal(:,7) =   bz_struc%bg(:,idir1) 
   bz_struc%normal(:,8) = - bz_struc%bg(:,idir1) 

   bz_struc%indsur(:,1) = (/ 4, 1,  7,  8, 2 /)
   bz_struc%indsur(:,2) = (/ 4, 2,  8,  9, 3 /)
   bz_struc%indsur(:,3) = (/ 4, 3,  9, 10, 4 /)
   bz_struc%indsur(:,4) = (/ 4, 4, 10, 11, 5 /)
   bz_struc%indsur(:,5) = (/ 4, 5, 11, 12, 6 /)
   bz_struc%indsur(:,6) = (/ 4, 6, 12,  7, 1 /)
   bz_struc%indsur(:,7) = (/ 6, 1,  2,  3,  4,  5, 6 /)
   bz_struc%indsur(:,8) = (/ 6, 7,  8,  9, 10, 11, 12 /)

   CALL find_vertices(bz_struc) 
   CALL compute_vertices(bz_struc) 

   bz_struc%letter_coord(:,2) = 0.5_DP* (bz_struc%vertex_coord(:,1) + &
                                         bz_struc%vertex_coord(:,12) )
   bz_struc%letter_coord(:,3) = 0.5_DP* (bz_struc%vertex_coord(:,2) + &
                                         bz_struc%vertex_coord(:,9) )
   bz_struc%letter_coord(:,4) = 0.5_DP* bz_struc%bg(:,idir1) 
   bz_struc%letter_coord(:,5) = bz_struc%letter_coord(:,idir1) + &
                                         bz_struc%letter_coord(:,4)
   bz_struc%letter_coord(:,6) = bz_struc%letter_coord(:,idir) + &
                                         bz_struc%letter_coord(:,4)


   IF (bz_struc%ibrav==12) THEN
      bz_struc%letter_list(2)=' X '
      bz_struc%letter_list(3)=' Y '
      bz_struc%letter_list(4)=' Z '
      bz_struc%letter_list(5)=' D '
      bz_struc%letter_list(6)=' A '
   ELSE
      bz_struc%letter_list(2)=' X '
      bz_struc%letter_list(3)=' Z '
      bz_struc%letter_list(4)=' Y '
      bz_struc%letter_list(5)=' A '
      bz_struc%letter_list(6)=' D '
   ENDIF

   CALL find_axis_coordinates(bz_struc)

END SUBROUTINE init_bz_16

SUBROUTINE compute_vertices(bz_struc)
!
!  This routine finds the coordinates of the vertex of the BZ, given
!  the index of the three planes that define each vertex.
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: xk(3)
INTEGER :: i

DO i = 1, bz_struc%nvertices
   CALL find_intersection( bz_struc%ivertex(:,i), bz_struc%normal, &
                           bz_struc%nfaces, xk)
   bz_struc%vertex_coord(:,i)=xk(:)
ENDDO
RETURN
END SUBROUTINE compute_vertices

SUBROUTINE find_letter_coordinate(bz_struc, letter, xk_let)
!
!  This routine checks if among the labels of special points defined
!  for each BZ there is the label letter and in that case it 
!  returns the coordinates of the point with that label. It stops
!  if the letter is not recognized.
!
IMPLICIT NONE
TYPE(bz), INTENT(IN) :: bz_struc
REAL(DP), INTENT(OUT) :: xk_let(3)
CHARACTER(LEN=3), INTENT(IN) :: letter

INTEGER :: i

DO i=1, bz_struc%nlett
   IF ((letter(1:2) == bz_struc%letter_list(i)(2:3) .AND. &
         bz_struc%letter_list(i)(1:1)/='g') .OR. &
        (letter(1:3) == bz_struc%letter_list(i)(1:3) )) THEN
         xk_let(:) = bz_struc%letter_coord(:,i)
         RETURN
   ENDIF
ENDDO

CALL errore('find_letter_coordinate','Letter not recognized '//TRIM(letter),1)

RETURN
END SUBROUTINE find_letter_coordinate


SUBROUTINE find_intersection( ivertex, normal, nfaces, outputk) 
!
! This routine receives as input the number of the three faces that define
! a vertex of the BZ, the reciprocal vectors perpendicular to all the 
! faces and gives as output the intersection point.
!
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: nfaces, ivertex(3)
REAL(DP), INTENT(IN) :: normal(3,nfaces)
REAL(DP), INTENT(OUT) :: outputk(3)

REAL(DP) :: a(3,3)
REAL(DP) :: b(3)
INTEGER :: ipol, jpol, ipiv(3), info

DO ipol=1,3
   DO jpol=1,3
      a(ipol,jpol) = normal(jpol,ivertex(ipol))
   ENDDO
   b(ipol)=0.5_DP * ( normal(1,ivertex(ipol))**2  &
                   +  normal(2,ivertex(ipol))**2  &
                   +  normal(3,ivertex(ipol))**2 )
ENDDO
CALL DGESV(3,1,a,3,ipiv,b,3,info)
IF (ABS(info) > 0) CALL errore('find_intersection', &
                               'problem in linear system',ABS(info))
outputk=b

RETURN
END SUBROUTINE find_intersection

SUBROUTINE find_bz_type(ibrav, celldm, ibz)
!
!  This routine identifies the bz type that corresponds to the given
!  bravais lattice and structural parameters. In the orthorombic
!  case it exchanges b/a and c/a so that a < b < c and sets the 
!  three variables change_a_b, change_a_c, change_b_c. It sets
!  celldm_bz
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: ibrav
REAL(DP), INTENT(IN) :: celldm(6)
INTEGER, INTENT(OUT) :: ibz
REAL(DP) :: value, celldm_c(6)
LOGICAL :: ldummy

IF (ibrav==1) THEN
   ibz=1
ELSEIF (ibrav==2) THEN
   ibz=2
ELSEIF (ibrav==3) THEN
   ibz=3
ELSEIF (ibrav==4) THEN
   ibz = 13
ELSEIF (ibrav==5) THEN
   IF (celldm(4) > 0.0_DP) THEN
      ibz=14
   ELSE
      ibz=15
   ENDIF
ELSEIF (ibrav==6) THEN
   ibz=4
ELSEIF (ibrav==7) THEN
   IF (celldm(3) > 1.0_DP) THEN
      ibz=6
   ELSE
      ibz=5
   ENDIF
ELSEIF (ibrav==8) THEN
   ibz=7
ELSEIF ( ibrav==9 .OR. ibrav==-9 ) THEN
   ibz=12
ELSEIF (ibrav==10) THEN
   celldm_c=celldm
   IF ((celldm(2)< 1.0_DP) .OR. (celldm(3) < 1.0_DP) .OR. &
             (celldm(2) > celldm(3))) THEN
      CALL canonical_celldm(celldm, celldm_c, ldummy, ldummy, ldummy)
   ENDIF

   value=1.0_DP/celldm_c(2)**2 + 1.0_DP/celldm_c(3)**2
   IF (value < 1.0_DP-1.D-8) THEN
      ibz=8
   ELSEIF (value > 1.0_DP+1.D-8) THEN
      ibz=9
   ELSE
      ibz=10
   ENDIF
ELSEIF (ibrav==11) THEN
   ibz=11
ELSEIF (ibrav==12.OR.ibrav==-12) THEN
   ibz=16
ELSEIF (ibrav==13.OR.ibrav==14) THEN
   CALL errore('find_bz_type','This ibrav is not supported',1)
ELSE
   CALL errore('find_bz_type','Wrong ibrav',1)
ENDIF

RETURN
END SUBROUTINE find_bz_type

SUBROUTINE find_vertices(bz_struc) 
!
!  This routine uses the definition of the vertices of each face to 
!  identify, for each vertex, the three faces that define it.
!  At least three faces must contain the vertex
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc

INTEGER :: ivert, iface, i, iv

DO ivert = 1, bz_struc%nvertices
   iv=1
   DO iface=1, bz_struc%nfaces
      DO i=2, bz_struc%indsur(1,iface)+1
         IF (bz_struc%indsur(i,iface) == ivert) THEN
             bz_struc%ivertex(iv,ivert)= iface
             iv=iv+1
             IF (iv > 3) GOTO 100
             EXIT
         ENDIF
      ENDDO
   ENDDO
   CALL errore('find_vertices','face not found',ivert)
100 CONTINUE   
ENDDO

RETURN
END SUBROUTINE find_vertices

SUBROUTINE check_one_face_centerer_orthorombic(bz_struc)
!
!  This subroutine changes the order of a and b axis in one face centered
!  orthorhombic lattice so that a < b. The bg corresponding to
!  such a rotated orthorombic lattice are used to generate the BZ. 
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: buffer

IF (bz_struc%celldm(2) < 1.0_DP) THEN
   bz_struc%switch_a_b=.TRUE.
   buffer=bz_struc%celldm(1)
   bz_struc%celldm(1)=bz_struc%celldm(1)*bz_struc%celldm(2)
   bz_struc%celldm(3)=bz_struc%celldm(3)/bz_struc%celldm(2)
   bz_struc%celldm(2)=1.0_DP / bz_struc%celldm(2)
ENDIF

CALL direct_and_reciprocal_lattice(bz_struc)

RETURN
END SUBROUTINE check_one_face_centerer_orthorombic

SUBROUTINE check_orthorombic(bz_struc)
!
!  This subroutine changes the order of a, b, and c axis in an
!  orthorhombic lattice so that a < b < c. The bg corresponding to
!  such a rotated orthorombic are used to generate the BZ.
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
!
!   first exchange the axis if a, b, and c are not ordered in the
!   sequence a < b < c
!
CALL canonical_celldm(bz_struc%celldm, bz_struc%celldm, bz_struc%switch_a_b, &
                      bz_struc%switch_b_c, bz_struc%rotate_a_b_c)

CALL direct_and_reciprocal_lattice(bz_struc)

RETURN
END SUBROUTINE check_orthorombic

SUBROUTINE direct_and_reciprocal_lattice(bz_struc)
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: alat
!
!  generate direct lattice vectors
!
CALL latgen(bz_struc%ibrav,bz_struc%celldm,bz_struc%at(:,1), &
bz_struc%at(:,2),bz_struc%at(:,3),bz_struc%omega)
!
!  generate reciprocal lattice vectors
!
alat = bz_struc%celldm(1)
bz_struc%at(:,:) = bz_struc%at(:,:) / alat

CALL recips( bz_struc%at(:,1), bz_struc%at(:,2), bz_struc%at(:,3), &
             bz_struc%bg(:,1), bz_struc%bg(:,2), bz_struc%bg(:,3) )

RETURN
END SUBROUTINE direct_and_reciprocal_lattice

SUBROUTINE adjust_orthorombic(bz_struc)
!
!   This routine rotates the coordinates of the x,y,z points of bz_struc
!   so that they correspond to the original orthorombic lattice
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: buffer, buffer3(3)
INTEGER :: ivert, ilett

IF (bz_struc%switch_a_b) THEN
   DO ivert = 1, bz_struc%nvertices
      buffer = bz_struc%vertex_coord(1,ivert)
      bz_struc%vertex_coord(1,ivert)=bz_struc%vertex_coord(2,ivert)
      bz_struc%vertex_coord(2,ivert)=buffer
   ENDDO
   DO ilett = 1, bz_struc%nlett
      buffer = bz_struc%letter_coord(1,ilett)
      bz_struc%letter_coord(1,ilett)=bz_struc%letter_coord(2,ilett)
      bz_struc%letter_coord(2,ilett)=buffer
   ENDDO
   bz_struc%vertex_coord=bz_struc%vertex_coord*bz_struc%celldm(2)
   bz_struc%letter_coord=bz_struc%letter_coord*bz_struc%celldm(2)
   buffer3=bz_struc%xi
   bz_struc%xi=0.0_DP
   bz_struc%xi(1) = bz_struc%yi(2) * bz_struc%celldm(2)
   bz_struc%yi=0.0_DP
   bz_struc%yi(2) = buffer3(1) * bz_struc%celldm(2)
   bz_struc%zi = bz_struc%zi * bz_struc%celldm(2)
ELSEIF (bz_struc%rotate_a_b_c) THEN
   DO ivert = 1, bz_struc%nvertices
      buffer = bz_struc%vertex_coord(1,ivert)
      bz_struc%vertex_coord(1,ivert)=bz_struc%vertex_coord(3,ivert)
      bz_struc%vertex_coord(3,ivert)=bz_struc%vertex_coord(2,ivert)
      bz_struc%vertex_coord(2,ivert)=buffer
   ENDDO
   DO ilett = 1, bz_struc%nlett
      buffer = bz_struc%letter_coord(1,ilett)
      bz_struc%letter_coord(1,ilett)=bz_struc%letter_coord(3,ilett)
      bz_struc%letter_coord(3,ilett)=bz_struc%letter_coord(2,ilett)
      bz_struc%letter_coord(2,ilett)=buffer
   ENDDO
   bz_struc%vertex_coord=bz_struc%vertex_coord*bz_struc%celldm(2)
   bz_struc%letter_coord=bz_struc%letter_coord*bz_struc%celldm(2)
   buffer3 = bz_struc%xi 
   bz_struc%xi=0.0_DP
   bz_struc%xi(1) = bz_struc%zi(3) * bz_struc%celldm(2)
   bz_struc%zi=0.0_DP
   bz_struc%zi(3) = bz_struc%yi(2) * bz_struc%celldm(2)
   bz_struc%yi=0.0_DP
   bz_struc%yi(2) = buffer3(1) * bz_struc%celldm(2)
ENDIF

IF (bz_struc%switch_b_c) THEN
!
!  first adjust the coordinates of the vertices
!
   DO ivert = 1, bz_struc%nvertices
      buffer = bz_struc%vertex_coord(2,ivert)
      bz_struc%vertex_coord(2,ivert)=bz_struc%vertex_coord(3,ivert)
      bz_struc%vertex_coord(3,ivert)=buffer
   ENDDO
!
!  then adjust the coordinates of the letters
!
   DO ilett = 1, bz_struc%nlett
      buffer = bz_struc%letter_coord(2,ilett)
      bz_struc%letter_coord(2,ilett)=bz_struc%letter_coord(3,ilett)
      bz_struc%letter_coord(3,ilett)=buffer
   ENDDO
!
!  then adjust the crossing of the coordinate axis with the surfaces
!
  buffer3=bz_struc%yi
  bz_struc%yi=0.0_DP
  bz_struc%yi(2)=bz_struc%zi(3)
  bz_struc%zi=0.0_DP
  bz_struc%zi(3)=buffer3(2)
ENDIF

RETURN
END SUBROUTINE adjust_orthorombic

SUBROUTINE adjust_orthorombic_vect(bz_struc,vect)
!
!   This routine rotates the coordinates of a vector of bz_struc
!   so that they correspond to the original orthorombic lattice
!
IMPLICIT NONE
TYPE(bz), INTENT(IN) :: bz_struc
REAL(DP), INTENT(INOUT) :: vect(3)
REAL(DP) :: buffer

IF (bz_struc%switch_b_c) THEN
   buffer = vect(2)
   vect(2)=vect(3)
   vect(3)=buffer
ENDIF

IF (bz_struc%switch_a_b) THEN
   buffer = vect(1)
   vect(1)=vect(2)
   vect(2)=buffer
ELSEIF (bz_struc%rotate_a_b_c) THEN
   buffer = vect(3)
   vect(3)=vect(1)
   vect(1)=vect(2)
   vect(2)=buffer
ENDIF

RETURN
END SUBROUTINE adjust_orthorombic_vect

SUBROUTINE adjust_one_face_centered_orthorombic(bz_struc)
!
!   This routine rotates the coordinates of the x, y, z points of bz_struc
!   so that they correspond to the original one face centered orthorombic 
!   lattice
!
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: buffer, buffer3(3)
INTEGER :: ivert, ilett

IF (bz_struc%switch_a_b) THEN
   DO ivert = 1, bz_struc%nvertices
      buffer = bz_struc%vertex_coord(1,ivert)
      bz_struc%vertex_coord(1,ivert)=bz_struc%vertex_coord(2,ivert)
      bz_struc%vertex_coord(2,ivert)=buffer
   ENDDO
   DO ilett = 1, bz_struc%nlett
      buffer = bz_struc%letter_coord(1,ilett)
      bz_struc%letter_coord(1,ilett)=bz_struc%letter_coord(2,ilett)
      bz_struc%letter_coord(2,ilett)=buffer
   ENDDO
   bz_struc%vertex_coord=bz_struc%vertex_coord*bz_struc%celldm(2)
   bz_struc%letter_coord=bz_struc%letter_coord*bz_struc%celldm(2)
   buffer3=bz_struc%xi
   bz_struc%xi=0.0_DP
   bz_struc%xi(1) = bz_struc%yi(2) * bz_struc%celldm(2)
   bz_struc%yi=0.0_DP
   bz_struc%yi(2) = buffer3(1) * bz_struc%celldm(2)
   bz_struc%zi = bz_struc%zi * bz_struc%celldm(2)
ENDIF
RETURN
END SUBROUTINE adjust_one_face_centered_orthorombic

SUBROUTINE canonical_celldm(celldm, celldm_c, switch_a_b, switch_b_c, &
                                              rotate_a_b_c)
!
!  This subroutine changes the order of a, b, and c axis in an
!  orthorhombic lattice so that a < b < c. The bg corresponding to
!  such a rotated orthorombic are used to generate the BZ and then
!  the coordinates are rotated so as to recover the original order.
!
IMPLICIT NONE
REAL(DP), INTENT(IN) :: celldm(6)
REAL(DP), INTENT(OUT) :: celldm_c(6)
LOGICAL :: switch_a_b, switch_b_c, rotate_a_b_c
REAL(DP) :: buffer

celldm_c = celldm
switch_a_b = .FALSE.
switch_b_c = .FALSE.
rotate_a_b_c = .FALSE.
!
!  If b > c switch them
!
IF ( celldm_c(2) > celldm_c(3) ) THEN
   switch_b_c=.TRUE.
   buffer=celldm_c(2)
   celldm_c(2)=celldm_c(3)
   celldm_c(3)=buffer
ENDIF

IF ( celldm_c(2) < 1.0_DP .AND. celldm_c(3) > 1.0_DP ) THEN
!
!  In this case a > b but a < c, just switch a and b
!
   switch_a_b=.TRUE.
   celldm_c(1) = celldm_c(2) * celldm_c(1)
   celldm_c(3) = celldm_c(3) / celldm_c(2) 
   celldm_c(2) = 1.0_DP / celldm_c(2) 
ELSEIF ( celldm_c(2) < 1.0_DP ) THEN
!
! In this case  b < c < a, we rotate a, b, and c 
!
   rotate_a_b_c=.TRUE.
   buffer=celldm_c(3)
   celldm_c(1) = celldm_c(2) * celldm_c(1)
   celldm_c(3) = 1.0_DP / celldm_c(2) 
   celldm_c(2) = buffer / celldm_c(2) 
ENDIF     

RETURN
END SUBROUTINE canonical_celldm

SUBROUTINE find_axis_coordinates(bz_struc)
IMPLICIT NONE
TYPE(bz), INTENT(INOUT) :: bz_struc
REAL(DP) :: x0(3), vect(3), xi(3), xmin
INTEGER :: ifaces

x0 = 0.0_DP
vect=0.0_DP
vect(1)=1.0_DP
xmin=1.D20
DO ifaces=1, bz_struc%nfaces
   IF (ABS(bz_struc%normal(1,ifaces)) > 1.d-9) THEN
      CALL inter_plane_line(x0, vect, bz_struc%normal(:,ifaces), xi)
      IF (xi(1) > 0.0_DP .AND. xi(1) < xmin) THEN
         bz_struc%xi=xi
         bz_struc%xaxis=ifaces
         xmin=xi(1)
      ENDIF
   ENDIF
ENDDO

x0 = 0.0_DP
vect=0.0_DP
vect(2)=1.0_DP
xmin=1.D20
DO ifaces=1, bz_struc%nfaces
   IF (ABS(bz_struc%normal(2,ifaces)) > 1.d-9) THEN
      CALL inter_plane_line(x0, vect, bz_struc%normal(:,ifaces), xi)
      IF (xi(2) > 0.0_DP .AND. xi(2) < xmin) THEN
         bz_struc%yi=xi
         bz_struc%yaxis=ifaces
         xmin=xi(2)
      ENDIF
   ENDIF
ENDDO

x0 = 0.0_DP
vect=0.0_DP
vect(3)=1.0_DP
xmin=1.D20
DO ifaces=1, bz_struc%nfaces
   IF (ABS(bz_struc%normal(3,ifaces)) > 1.d-9) THEN
      CALL inter_plane_line(x0, vect, bz_struc%normal(:,ifaces), xi)
      IF (xi(3) > 0.0_DP .AND. xi(3) < xmin) THEN
         bz_struc%zi=xi
         bz_struc%zaxis=ifaces
         xmin=xi(3)
      ENDIF
   ENDIF
ENDDO

RETURN
END SUBROUTINE find_axis_coordinates

SUBROUTINE inter_plane_line(x0, vect, bplane, xk)
!
!  This routine finds the intersection between the line passing through
!  x0 and parallel to vect, and the plane passing through bplane/2 and
!  perpendicular to bplane
!
IMPLICIT NONE

REAL(DP), INTENT(IN) :: x0(3), vect(3), bplane(3)
REAL(DP), INTENT(OUT) :: xk(3)

REAL(DP) :: lambda
REAL(DP) :: x0b, vb, b2

vb=vect(1)*bplane(1)+vect(2)*bplane(2)+vect(3)*bplane(3)
IF (ABS(vb) < 1.d-9) &
   CALL errore('inter_plane_line','The line and the plane are parallel',1)

x0b=x0(1)*bplane(1)+x0(2)*bplane(2)+x0(3)*bplane(3)
b2=bplane(1)*bplane(1)+bplane(2)*bplane(2)+bplane(3)*bplane(3)

lambda=(b2-2.0_DP*x0b)/(2.0_DP * vb)

xk(:)=x0(:)+lambda*vect(:)

RETURN
END SUBROUTINE inter_plane_line

SUBROUTINE transform_label_coord(ibrav, celldm, xk, letter, label_list, &
                                 npk_label, nks, k_points, point_label_type )
!
!  This routine transforms the labels in the array letter into k points
!  coordinates that are put in the array xk in the position indicated
!  by label_list. If k_point='crystal' the coordinates are tranformed
!  in the basis of the crystal. point_label_type selects the type of
!  labels. npk_label is the size of the array letter and label_list,
!  while nks is the size of the array xk.
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: npk_label
INTEGER, INTENT(IN) :: nks
INTEGER, INTENT(IN) :: ibrav
INTEGER, INTENT(IN) :: label_list(npk_label)
REAL(DP), INTENT(IN) :: celldm(6)
REAL(DP), INTENT(INOUT) :: xk(3, nks)
CHARACTER(LEN=3), INTENT(IN) :: letter(npk_label)
CHARACTER(LEN=*), INTENT(IN) :: k_points, point_label_type
INTEGER :: bzt, i
REAL(DP) :: omega, at(3,3), bg(3,3), xk_buffer(3)
TYPE(bz) :: bz_struc
!
!    Find the brillouin zone type
!
   CALL find_bz_type(ibrav, celldm, bzt)
!
!    generate direct lattice vectors 
!
   CALL latgen(ibrav,celldm,at(:,1),at(:,2),at(:,3),omega)
!
!   we use at in units of celldm(1)
!
   at=at/celldm(1)
!
! generate reciprocal lattice vectors
!
   CALL recips( at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3) )
!
! load the information on the Brillouin zone
!
   CALL set_label_type(bz_struc, point_label_type)
   CALL allocate_bz(ibrav, bzt, bz_struc, celldm, at, bg )
   CALL init_bz(bz_struc)
!
! find for each label the corresponding coordinates and save them
! on the k point list
!
   DO i=1, npk_label
      CALL find_letter_coordinate(bz_struc, letter(i), xk_buffer)
!
!  The output of this routine is in cartesian coordinates. If the other
!  k points are in crystal coordinates we transform xk_buffer to the bg
!  base.
!
      IF (TRIM(k_points)=='crystal') &
         CALL cryst_to_cart( 1, xk_buffer, at, -1 ) 
      xk(:,label_list(i))=xk_buffer(:)
   ENDDO
  CALL deallocate_bz(bz_struc)
RETURN
END SUBROUTINE transform_label_coord

SUBROUTINE find_n1n2_monoclinic(n1, n2, idir, bz_struc)
!
!   This routine finds the six reciprocal lattices closest to the origin and
!   order them in order of increasing angle with the x-axis. These are the six
!   normals to the faces of the monoclinic Brillouin zone.
!
USE constants, ONLY : pi
IMPLICIT NONE
INTEGER, INTENT(IN) :: idir
INTEGER, INTENT(INOUT) :: n1(6), n2(6)
TYPE(bz), INTENT(INOUT) :: bz_struc
INTEGER :: in1, in2, ind(6), inaux(6), nfound, ivect, isub
REAL(DP) :: min_mod, max_mod, modul, vect(3), save_mod(6), save_angle(6), angle
INTEGER :: npx_
LOGICAL :: done
!
!  Search among (2*npx+1)**2 vectors. Not all cases are covered by this
!  search, but the largest part should be.
!
npx_=bz_struc%npx
min_mod=0.0_DP
nfound=0
DO in1=-npx_,npx_
   DO in2=-npx_,npx_
      IF ( in1 == 0 .AND. in2 == 0 ) CYCLE
      vect(1:3)= in1 * bz_struc%bg(1:3,1) + in2 * bz_struc%bg(1:3,idir)
      modul = SQRT( vect(1) ** 2 + vect(2) ** 2 + vect(3) ** 2 ) 
      angle = ACOS( vect(1) / modul )
      IF (vect(idir) < 0.0_DP ) angle= 2.0_DP * pi - angle
!
!    Check that if new vector is a multiple of one already found
!
      done=.FALSE.
      DO ivect=1, nfound
         done=done.OR.(ABS(angle-save_angle(ivect))< 1.D-7)
      ENDDO
      IF ( done ) THEN
!
!    In this case the new vector is a multiple of one vector already found.
!    Substitute if its modulus is lower
!
         DO ivect=1, nfound
            IF (ABS(angle-save_angle(ivect))< 1.D-7 .AND. &
                                          modul < save_mod(ivect)) THEN
               n1(ivect)=in1
               n2(ivect)=in2
               save_mod(ivect)=modul
               save_angle(ivect)=angle
            END IF
         ENDDO
!
!   ricompute the maximum modulus of the found vectors
!
         min_mod=0.0_DP
         DO ivect=1, nfound
            IF (save_mod(ivect) > min_mod) min_mod=save_mod(ivect)
         ENDDO
      ELSE
         IF ( nfound < 6 .OR. modul < min_mod ) THEN
            IF (nfound < 6 ) THEN
!
!    In this case we have not yet found six vectors. 
!
               nfound = nfound + 1
               n1(nfound)=in1
               n2(nfound)=in2
               save_mod(nfound) = modul
               save_angle(nfound) = angle
               IF (modul > min_mod) min_mod = modul
            ELSE
!
!    In this case we substitute the vector with maximum modulus with this one
!
               max_mod=0.0_DP
               isub=0
               DO ivect=1,6
                  IF (save_mod(ivect) > max_mod) THEN
                     isub = ivect
                     max_mod=save_mod(ivect)
                  ENDIF
               END DO
               IF (isub==0) CALL errore('find_n1n2_monoclinic','Problem with isub',1)
               n1(isub)=in1
               n2(isub)=in2
               save_mod(isub)=modul
               save_angle(isub)=angle
               min_mod=0.0_DP
               DO ivect=1, 6
                  IF (save_mod(ivect) > min_mod) min_mod=save_mod(ivect)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
ENDDO
IF (nfound /= 6) CALL errore('find_n1n2_monoclinic','Problem with nfound',1)
!
!  now order the six vectors, according to the angle they form with the x axis
!
!
!  If n1 or n2 is at the limit of the checked cell, tell the user to
!  double the parameter npx
!
DO ivect=1,6
   IF (n1(ivect)==npx_ .OR. n2(ivect)==npx_) &
      CALL errore('find_n1n2_monoclinic','Difficult monoclinic cell, &
                                                  &double npx',1)
ENDDO
!
!  now order the six vectors, according to the angle they form with the x axis
!
ind(1)=0
CALL hpsort(6, save_angle, ind)
!
!  order n1 and n2
!
inaux=n1
n1(:) =inaux(ind(:))
inaux=n2
n2(:) =inaux(ind(:))

RETURN
END SUBROUTINE find_n1n2_monoclinic

END MODULE bz_form
