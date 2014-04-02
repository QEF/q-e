!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE find_group(nrot,smat,gname,code_group)
!
!  Given a group of nrot rotation matrices smat (in cartesian coordinates)
!  this routine finds the name of the point group. It assumes but does not
!  check that: 
!  1) The nrot matrices smat are actually a group.
!  2) The group is one of the thirty-two point groups.
!

USE kinds, ONLY : DP

IMPLICIT NONE

INTEGER :: nrot, code_group
REAL(DP) :: smat(3,3,nrot)
CHARACTER (LEN=11) :: gname, group_name
INTEGER :: noperation(6), irot, ts, tipo_sym
!
! For each possible group operation the function tipo_sym gives a code
!   1 identity, 
!   2 inversion, 
!   3 proper rotation <> 180, 
!   4 proper rotation 180 degrees, 
!   5 mirror, 
!   6 improper rotation
! the variable noperation counts how many operations are present in the group.
!
noperation=0
DO irot=1,nrot
   ts=tipo_sym(smat(1,1,irot))
   noperation(ts)=noperation(ts)+1
END DO

IF (noperation(1).ne.1) call errore('find_group','the group has not identity',1)

code_group=0

IF (noperation(2)==0) THEN
!
!  There is not inversion
!
   IF (nrot==1) THEN
      code_group=1                                          ! C_1
   ELSEIF (nrot==2) THEN
      IF (noperation(4)==1) code_group=4                    ! C_2
      IF (noperation(5)==1) code_group=3                    ! C_s
   ELSEIF (nrot==3) THEN
      IF (noperation(3)==2) code_group=5                    ! C_3
   ELSEIF (nrot==4) THEN
      IF (noperation(6)>0)  code_group=26                   ! S_4
      IF (noperation(5)>0.and.code_group==0) code_group=12  ! C_2v
      IF (noperation(3)>0.and.code_group==0) code_group=6   ! C_4
      IF (noperation(4)>0.and.code_group==0) code_group=8   ! D_2
   ELSEIF (nrot==6) THEN
      IF (noperation(5)==3) code_group=13                   ! C_3v
      IF (noperation(5)==1) code_group=17                   ! C_3h
      IF (noperation(4)==3.and.code_group==0) code_group=9  ! D_3
      IF (noperation(3)>0.and.code_group==0) code_group=7   ! C_6
   ELSEIF (nrot==8) THEN
      IF (noperation(5)==4) code_group=14                   ! C_4v
      IF (noperation(5)==2) code_group=24                   ! D_2d
      IF (noperation(3)>0.and.code_group==0) code_group=10  ! D_4
   ELSEIF (nrot==12) THEN
      IF (noperation(5)==6) code_group=15                   ! C_6v
      IF (noperation(5)==4) code_group=21                   ! D_3h
      IF (noperation(4)>6.and.code_group==0) code_group=11  ! D_6
      IF (noperation(3)>0.and.code_group==0) code_group=28  ! T
   ELSEIF (nrot==24) THEN
      IF (noperation(5)>0) code_group=30                    ! T_d
      IF (noperation(5)==0) code_group=31                   ! O
   ELSE
      CALL errore('find_group','wrong number of elements',1)
   ENDIF

ELSEIF (noperation(2)==1) THEN
!
!  There is inversion
!
   IF (nrot==2) THEN
      code_group=2                                          ! C_i
   ELSEIF (nrot==4) THEN
      code_group=16                                         ! C_2h
   ELSEIF (nrot==6) THEN
      code_group=27                                         ! S_6
   ELSEIF (nrot==8) THEN
      IF (noperation(5)==3) code_group=20                   ! D_2h
      IF (noperation(5)==1) code_group=18                   ! C_4h
   ELSEIF (nrot==12) THEN
      IF (noperation(5)==3) code_group=25                   ! D_3d
      IF (noperation(5)==1) code_group=19                   ! C_6h
   ELSEIF (nrot==16) THEN
      IF (noperation(5)==5) code_group=22                   ! D_4h
   ELSEIF (nrot==24) THEN
      IF (noperation(5)>6) code_group=23                    ! D_6h
      IF (noperation(5)==3) code_group=29                   ! T_h
   ELSEIF (nrot==48) THEN
      code_group=32                                         ! O_h
   ELSE
      CALL errore('find_group','wrong number of elements',1)
   ENDIF
ELSE
   CALL errore('find_group','too many inversions',1)
ENDIF
IF (code_group==0) call errore('find_group','incompatible operations',1)

gname=group_name(code_group)

RETURN
END SUBROUTINE find_group

!--------------------------------------------------------------------------
FUNCTION group_name(code)
!--------------------------------------------------------------------------
! This function receives a code of the group and provides the name of the
! group. The order is the following:
!
!   1  "C_1 "     11 "D_6 "     21 "D_3h"     31 "O   " 
!   2  "C_i "     12 "C_2v"     22 "D_4h"     32 "O_h "  
!   3  "C_s "     13 "C_3v"     23 "D_6h" 
!   4  "C_2 "     14 "C_4v"     24 "D_2d" 
!   5  "C_3 "     15 "C_6v"     25 "D_3d" 
!   6  "C_4 "     16 "C_2h"     26 "S_4 " 
!   7  "C_6 "     17 "C_3h"     27 "S_6 " 
!   8  "D_2 "     18 "C_4h"     28 "T   " 
!   9  "D_3 "     19 "C_6h"     29 "T_h " 
!   10 "D_4 "     20 "D_2h"     30 "T_d "
!
IMPLICIT NONE

INTEGER :: code
CHARACTER(LEN=11) :: group_name

CHARACTER(LEN=11) :: gname(32)

data gname  / "C_1 (1)    ", "C_i (-1)   ", "C_s (m)    ", "C_2  (2)   ", &
              "C_3 (3)    ", "C_4 (4)    ", "C_6 (6)    ", "D_2  (222) ", &
              "D_3 (32)   ", "D_4 (422)  ", "D_6 (622)  ", "C_2v (mm2) ", &
              "C_3v (3m)  ", "C_4v (4mm) ", "C_6v (6mm) ", "C_2h (2/m) ", &
              "C_3h (-6)  ", "C_4h (4/m) ", "C_6h (6/m) ", "D_2h (mmm) ", &
              "D_3h (-62m)", "D_4h(4/mmm)", "D_6h(6/mmm)", "D_2d (-42m)", &
              "D_3d (-3m) ", "S_4 (-4)   ", "S_6 (-3)   ", "T    (23)  ", & 
              "T_h (m-3)  ", "T_d (-43m) ", "O   (432)  ", "O_h (m-3m) "  /

IF (code < 1 .OR. code > 32 ) CALL errore('group_name','code is out of range',1)

group_name=gname(code)

RETURN
END FUNCTION group_name


!--------------------------------------------------------------------------
FUNCTION tipo_sym(s)
!--------------------------------------------------------------------------
! This function receives a 3x3 orthogonal matrix which is a symmetry 
! operation of the point group of the crystal written in cartesian 
! coordinates and gives as output a code according to the following:
!
!  1   Identity
!  2   Inversion
!  3   Proper rotation of an angle <> 180 degrees
!  4   Proper rotation of 180 degrees
!  5   Mirror symmetry
!  6   Improper rotation
!
USE kinds, ONLY : DP

IMPLICIT NONE

REAL(DP), PARAMETER :: eps=1.d-7
REAL(DP) :: s(3,3), det, det1
INTEGER  :: tipo_sym
!
! Check for identity
!
IF ((ABS(s(1,1)-1.d0) < eps).AND. &
    (ABS(s(2,2)-1.d0) < eps).AND. &
    (ABS(s(3,3)-1.d0) < eps).AND. &
    (ABS(s(1,2)) < eps).AND.(ABS(s(2,1)) < eps).AND.(ABS(s(2,3)) < eps).AND. &
    (ABS(s(3,2)) < eps).AND.(ABS(s(1,3)) < eps).AND.(ABS(s(3,1)) < eps)) THEN
   tipo_sym=1
   RETURN
ENDIF
!
! Check for inversion
!
IF ((ABS(s(1,1)+1.d0) < eps).AND. &
    (ABS(s(2,2)+1.d0) < eps).AND. &
    (ABS(s(3,3)+1.d0) < eps).AND. &
    (ABS(s(1,2)) < eps).AND.(ABS(s(2,1)) < eps).AND.(ABS(s(2,3)) < eps).AND. &
    (ABS(s(3,2)) < eps).AND.(ABS(s(1,3)) < eps).AND.(ABS(s(3,1)) < eps)) THEN
   tipo_sym=2
   RETURN
ENDIF
!
! compute the determinant
!
det = s(1,1) * ( s(2,2) * s(3,3) - s(3,2) * s(2,3) )-   &
      s(1,2) * ( s(2,1) * s(3,3) - s(3,1) * s(2,3) )+   &
      s(1,3) * ( s(2,1) * s(3,2) - s(3,1) * s(2,2) ) 
!
! Determinant equal to 1: proper rotation
!
IF (abs(det-1.d0) < eps) THEN
!
!  check if an eigenvalue is equal to -1.d0 (180 rotation)
!
   det1=(s(1,1)+1.d0)*((s(2,2)+1.d0)*(s(3,3)+1.d0)-s(3,2)*s(2,3))-   &
         s(1,2)*       (s(2,1)*      (s(3,3)+1.d0)-s(3,1)*s(2,3))+   &
         s(1,3)*       (s(2,1)*s(3,2)             -s(3,1)*(s(2,2)+1.d0)) 

   IF (abs(det1) < eps) THEN
      tipo_sym=4     ! 180 proper rotation
   ELSE
      tipo_sym=3     ! proper rotation <> 180
   ENDIF
   RETURN
ENDIF
!
! Determinant equal to -1: mirror symmetry or improper rotation
!
IF (abs(det+1.d0) < eps) THEN
!
!  check if an eigenvalue is equal to 1.d0 (mirror symmetry)
!
   det1=(s(1,1)-1.d0)*((s(2,2)-1.d0)*(s(3,3)-1.d0)-s(3,2)*s(2,3))-   &
         s(1,2)*       (s(2,1)*      (s(3,3)-1.d0)-s(3,1)*s(2,3))+   &
         s(1,3)*       (s(2,1)*s(3,2)             -s(3,1)*(s(2,2)-1.d0)) 

   IF (abs(det1) < eps) THEN
      tipo_sym=5   ! mirror symmetry
   ELSE
      tipo_sym=6   ! improper rotation
   ENDIF
   RETURN
ELSE
   call errore('tipo_sym','symmetry not recognized',1)
ENDIF
END FUNCTION tipo_sym
!
!--------------------------------------------------------------------------
FUNCTION laue_class(code)
!--------------------------------------------------------------------------
! This function receives a code of the point group and provides the
! code of the point group that defines the Laue class (that is the point 
! group obtained by multipling by inversion).
! The order is the following:
!
!   1  "C_1 " -> 2   11 "D_6 " -> 23  21 "D_3h" -> 23  31 "O   " -> 32
!   2  "C_i " -> 2   12 "C_2v" -> 20  22 "D_4h" -> 22  32 "O_h " -> 32
!   3  "C_s " -> 16  13 "C_3v" -> 25  23 "D_6h" -> 23
!   4  "C_2 " -> 16  14 "C_4v" -> 22  24 "D_2d" -> 22
!   5  "C_3 " -> 27  15 "C_6v" -> 23  25 "D_3d" -> 25
!   6  "C_4 " -> 18  16 "C_2h" -> 16  26 "S_4 " -> 18
!   7  "C_6 " -> 19  17 "C_3h" -> 19  27 "S_6 " -> 27
!   8  "D_2 " -> 20  18 "C_4h" -> 18  28 "T   " -> 29
!   9  "D_3 " -> 25  19 "C_6h" -> 19  29 "T_h " -> 29
!   10 "D_4 " -> 22  20 "D_2h" -> 20  30 "T_d " -> 32
!
IMPLICIT NONE

INTEGER :: code
INTEGER :: laue_class

INTEGER :: laue(32)

DATA  laue  / 2,   2, 16, 16, 27, 18, 19, 20, 25, 22,  &
              23, 20, 25, 22, 23, 16, 19, 18, 19, 20,  &
              23, 22, 23, 22, 25, 18, 27, 29, 29, 32,  &
              32, 32 /

IF (code < 1 .OR. code > 32 ) CALL errore('laue_class','code is out of range',1)

laue_class=laue(code)

RETURN
END FUNCTION laue_class

