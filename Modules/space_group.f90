!
! Copyright (C) 2014 Federico Zadra
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE space_group
USE kinds, ONLY: DP
   IMPLICIT NONE

   SAVE
   PRIVATE

   PUBLIC sym_brav, find_equivalent_tau

   CONTAINS
   SUBROUTINE sym_brav(space_group_number,sym_n,ibrav)
   
   ! Sym_brav ->   
   ! input    spacegroup number
   ! output   sym_n = symmetries number
   !          ibrav = Bravais lattice number 

      INTEGER, INTENT(IN) :: space_group_number
      INTEGER, INTENT(OUT) :: sym_n,ibrav
      
      simmetria: SELECT CASE (space_group_number)
      !Triclinic 1-2
      CASE (1)
         sym_n=1
         ibrav=14
      CASE (2)
         sym_n=2
         ibrav=14
      !Monoclinic 3-15
      CASE (3) !P2
         sym_n=2
         ibrav=12
      CASE (4) !P2(1)
         sym_n=2
         ibrav=12
      CASE (5) !C2
         sym_n=2
         ibrav=13
      CASE (6) !PM
         sym_n=2
         ibrav=12
      CASE (7) !Pc
         sym_n=2
         ibrav=12
      CASE (8) !Cm
         sym_n=2
         ibrav=13
      CASE (9) !Cc
         sym_n=2
         ibrav=13
      CASE (10) !P2/m
         sym_n=4
         ibrav=12
      CASE (11) !P2(1)/m
         sym_n=4
         ibrav=12
      CASE (12) !C2/m
         sym_n=4
         ibrav=13
      CASE (13) !P2/c
         sym_n=4
         ibrav=12
      CASE (14) !P2(1)/c
         sym_n=4
         ibrav=12
      CASE (15) !C2/c
         sym_n=4
         ibrav=13
      !Orthorhombic
      CASE (16) !P222
         sym_n=4
         ibrav=8
      CASE (17) !P222(1)
         sym_n=4
         ibrav=8
      CASE (18) !P2(1)2(1)2
         sym_n=4
         ibrav=8
      CASE (19) !P2(1)2(1)2(1)
         sym_n=4
         ibrav=8
      CASE (20) !C222(1)
         sym_n=4
         ibrav=9
      CASE (21) !C222
         sym_n=4
         ibrav=9
      CASE (22) !F222
         sym_n=4
         ibrav=10
      CASE (23) !I222
         sym_n=4
         ibrav=11
      CASE (24) !I2(1)2(1)2(1)
         sym_n=4
         ibrav=11
      CASE (25) !Pmm2
         sym_n=4
         ibrav=8
      CASE (26) !Pmc2(1)
         sym_n=4
         ibrav=8
      CASE (27) !Pcc2
         sym_n=4
         ibrav=8
      CASE (28) !Pma2
         sym_n=4
         ibrav=8
      CASE (29) !Pca2(1)
         sym_n=4
         ibrav=8
      CASE (30) !Pnc2
         sym_n=4
         ibrav=8
      CASE (31) !Pmn2(1)
         sym_n=4
         ibrav=8
      CASE (32) !Pba2
         sym_n=4
         ibrav=8
      CASE (33) !Pna2(1)
         sym_n=4
         ibrav=8
      CASE (34) !Pnn2
         sym_n=4
         ibrav=8
      CASE (35) !Cmm2
         sym_n=4
         ibrav=9
      CASE (36) !Cmc2(1)
         sym_n=4
         ibrav=9
      CASE (37) !Ccc2
         sym_n=4
         ibrav=9
      CASE (38) !Amm2
         sym_n=4
         ibrav=91
      CASE (39) !Abm2
         sym_n=4
         ibrav=91
      CASE (40) !Ama2
         sym_n=4
         ibrav=91
      CASE (41) !Aba2
         sym_n=4
         ibrav=91
      CASE (42) !Fmm2
         sym_n=4
         ibrav=10
      CASE (43) !Fdd2
         sym_n=4
         ibrav=10
      CASE (44) !Imm2
         sym_n=4
         ibrav=11
      CASE (45) !Iba2
         sym_n=4
         ibrav=11
      CASE (46) !Ima2
         sym_n=4
         ibrav=11
      CASE (47) !Pmmm
         sym_n=8
         ibrav=8
      CASE (48) !Pnnn
         sym_n=8
         ibrav=8
      CASE (49) !Pccm
         sym_n=8
         ibrav=8
      CASE (50) !Pban
         sym_n=8
         ibrav=8
      CASE (51) !Pmma
         sym_n=8
         ibrav=8
      CASE (52) !Pnna
         sym_n=8
         ibrav=8
      CASE (53) !Pmna
         sym_n=8
         ibrav=8
      CASE (54) !Pcca
         sym_n=8
         ibrav=8
      CASE (55) !Pbam
         sym_n=8
         ibrav=8
      CASE (56) !Pccn
         sym_n=8
         ibrav=8
      CASE (57) !Pbcm
         sym_n=8
         ibrav=8
      CASE (58) !Pnnm
         sym_n=8
         ibrav=8
      CASE (59) !Pmmn
         sym_n=8
         ibrav=8
      CASE (60) !Pbcn
         sym_n=8
         ibrav=8
      CASE (61) !Pbca
         sym_n=8
         ibrav=8
      CASE (62) !Pnma
         sym_n=8
         ibrav=8
      CASE (63) !Cmcm
         sym_n=8
         ibrav=9
      CASE (64) !Cmca
         sym_n=8
         ibrav=9
      CASE (65) !Cmmm
         sym_n=8
         ibrav=9
      CASE (66) !Cccm
         sym_n=8
         ibrav=9
      CASE (67) !Cmma
         sym_n=8
         ibrav=9
      CASE (68) !Ccca
         sym_n=8
         ibrav=9
      CASE (69) !Fmmm
         sym_n=8
         ibrav=10
      CASE (70) !Fddd
         sym_n=8
         ibrav=10
      CASE (71) !Immm
         sym_n=8
         ibrav=11
      CASE (72) !Ibam
         sym_n=8
         ibrav=11
      CASE (73) !Ibca
         sym_n=8
         ibrav=11
      CASE (74) !Imma
         sym_n=8
         ibrav=11
      !Tetragonal
      CASE (75) !P4
         sym_n=4
         ibrav=6
      CASE (76) !P4(1)
         sym_n=4
         ibrav=6
      CASE (77) !P4(2)
         sym_n=4
         ibrav=6
      CASE (78) !P4(3)
         sym_n=4
         ibrav=6
      CASE (79) !I4
         sym_n=4
         ibrav=7
      CASE (80) !I4(1)
         sym_n=4
         ibrav=7
      CASE (81) !P-4
         sym_n=4
         ibrav=6
      CASE (82) !I-4
         sym_n=4
         ibrav=7
      CASE (83) !P4/m
         sym_n=8
         ibrav=6
      CASE (84) !P4(2)/m
         sym_n=8
         ibrav=6
      CASE (85) !P4/n
         sym_n=8
         ibrav=6
      CASE (86) !P4(2)/n
         sym_n=8
         ibrav=6
      CASE (87) !I4/m
         sym_n=8
         ibrav=7
      CASE (88) !I4(1)/a
         sym_n=8
         ibrav=7
      CASE (89) !P422
         sym_n=8
         ibrav=6
      CASE (90) !P42(1)2
         sym_n=8
         ibrav=6
      CASE (91) !P4(1)22
         sym_n=8
         ibrav=6
      CASE (92) !P4(1)2(1)2
         sym_n=8
         ibrav=6
      CASE (93) !P4(2)22
         sym_n=8
         ibrav=6
      CASE (94) !P4(2)2(1)2
         sym_n=8
         ibrav=6
      CASE (95) !P4(3)22
         sym_n=8
         ibrav=6
      CASE (96) !P4(3)2(1)2
         sym_n=8
         ibrav=6
      CASE (97) !I422
         sym_n=8
         ibrav=7
      CASE (98) !I4(1)22
         sym_n=8
         ibrav=7
      CASE (99) !P4mm
         sym_n=8
         ibrav=6
      CASE (100) !P4bm
         sym_n=8
         ibrav=6
      CASE (101) !P4(2)cm
         sym_n=8
         ibrav=6
      CASE (102) !P4(2)nm
         sym_n=8
         ibrav=6
      CASE (103) !P4cc
         sym_n=8
         ibrav=6
      CASE (104) !P4nc
         sym_n=8
         ibrav=6
      CASE (105) !P4(2)mc
         sym_n=8
         ibrav=6
      CASE (106) !P4(2)bc
         sym_n=8
         ibrav=6
      CASE (107) !I4mm
         sym_n=8
         ibrav=7
      CASE (108) !I4cm
         sym_n=8
         ibrav=7
      CASE (109) !I4(!)md
         sym_n=8
         ibrav=7
      CASE (110) !I4(1)cd
         sym_n=8
         ibrav=7
      CASE (111) !P-42m
         sym_n=8
         ibrav=6
      CASE (112) !P-42c
         sym_n=8
         ibrav=6
      CASE (113) !P-42(1)m
         sym_n=8
         ibrav=6
      CASE (114) !P-42(1)c
         sym_n=8
         ibrav=6
      CASE (115) !P-4m2
         sym_n=8
         ibrav=6
      CASE (116) !P-4c2
         sym_n=8
         ibrav=6
      CASE (117) !P-4b2
         sym_n=8
         ibrav=6
      CASE (118) !P-4n2
         sym_n=8
         ibrav=6
      CASE (119) !I-4m2
         sym_n=8
         ibrav=7
      CASE (120) !I-4c2
         sym_n=8
         ibrav=7
      CASE (121) !I-42m
         sym_n=8
         ibrav=7
      CASE (122) !I-42d
         sym_n=8
         ibrav=7
      CASE (123) !P4/mmm
         sym_n=16
         ibrav=6
      CASE (124) !P4/mcc
         sym_n=16
         ibrav=6
      CASE (125) !P4/nbm
         sym_n=16
         ibrav=6
      CASE (126) !P4/nnc
         sym_n=16
         ibrav=6
      CASE (127) !P4/mbm
         sym_n=16
         ibrav=6
      CASE (128) !P4/mnc
         sym_n=16
         ibrav=6
      CASE (129) !P4/nmm
         sym_n=16
         ibrav=6
      CASE (130) !P4/ncc
         sym_n=16
         ibrav=6
      CASE (131) !P4(2)/mmc
         sym_n=16
         ibrav=6
      CASE (132) !P4(2)/mcm
         sym_n=16
         ibrav=6
      CASE (133) !P4(2)nbc
         sym_n=16
         ibrav=6
      CASE (134) !P4(2)/nnm
         sym_n=16
         ibrav=6
      CASE (135) !P4(2)/mbc
         sym_n=16
         ibrav=6
      CASE (136) !P4(2)/mnm
         sym_n=16
         ibrav=6
      CASE (137) !P4(2)/nmc
         sym_n=16
         ibrav=6
      CASE (138) !P4(2)/ncm
         sym_n=16
         ibrav=6
      CASE (139) !I4/mmm
         sym_n=16
         ibrav=7
      CASE (140) !I4/mcm
         sym_n=16
         ibrav=7
      CASE (141) !I4(1)/amd
         sym_n=16
         ibrav=7
      CASE (142) !I4(1)/acd
         sym_n=16
         ibrav=7
      ! Trigonal
      CASE (143) !P3
         sym_n=3
         ibrav=4
      CASE (144)
         sym_n=3
         ibrav=4
      CASE (145)
         sym_n=3
         ibrav=4
      CASE (146) !R3
         sym_n=3
         ibrav=5
      CASE (147)
         sym_n=6
         ibrav=4
      CASE (148) !R-3
         sym_n=6
         ibrav=5
      CASE (149)
         sym_n=6
         ibrav=4
      CASE (150) 
         sym_n=6
         ibrav=4
      CASE (151)
         sym_n=6
         ibrav=4
      CASE (152)
         sym_n=6
         ibrav=4
      CASE (153)
         sym_n=6
         ibrav=4
      CASE (154)
         sym_n=6
         ibrav=4
      CASE (155) !R32
         sym_n=6
         ibrav=5
      CASE (156)
         sym_n=6
         ibrav=4
      CASE (157) 
         sym_n=6
         ibrav=4
      CASE (158)
         sym_n=6
         ibrav=4
      CASE (159)
         sym_n=6
         ibrav=4
      CASE (160) !R3m
         sym_n=6
         ibrav=5
      CASE (161) !R3c
         sym_n=6
         ibrav=5
      CASE (162)
         sym_n=12
         ibrav=4
      CASE (163)
         sym_n=12
         ibrav=4
      CASE (164)
         sym_n=12
         ibrav=4
      CASE (165)
         sym_n=12
         ibrav=4
      CASE (166) !R-3m
         sym_n=12
         ibrav=5
      CASE (167) !R-3c
         sym_n=12
         ibrav=5
      ! Exagonal
      CASE (168)
         sym_n=6
         ibrav=4
      CASE (169)
         sym_n=6
         ibrav=4
      CASE (170)
         sym_n=6
         ibrav=4
      CASE (171)
         sym_n=6
         ibrav=4
      CASE (172)
         sym_n=6
         ibrav=4
      CASE (173)
         sym_n=6
         ibrav=4
      CASE (174)
         sym_n=6
         ibrav=4
      CASE (175)
         sym_n=12
         ibrav=4
      CASE (176)
         sym_n=12
         ibrav=4
      CASE (177)
         sym_n=12
         ibrav=4
      CASE (178)
         sym_n=12
         ibrav=4
      CASE (179)
         sym_n=12
         ibrav=4
      CASE (180)
         sym_n=12
         ibrav=4
      CASE (181)
         sym_n=12
         ibrav=4
      CASE (182)
         sym_n=12
         ibrav=4
      CASE (183)
         sym_n=12
         ibrav=4
      CASE (184)
         sym_n=12
         ibrav=4
      CASE (185)
         sym_n=12
         ibrav=4
      CASE (186)
         sym_n=12
         ibrav=4
      CASE (187)
         sym_n=12
         ibrav=4
      CASE (188)
         sym_n=12
         ibrav=4
      CASE (189)
         sym_n=12
         ibrav=4
      CASE (190)
         sym_n=12
         ibrav=4
      CASE (191)
         sym_n=24
         ibrav=4
      CASE (192)
         sym_n=24
         ibrav=4
      CASE (193)
         sym_n=24
         ibrav=4
      CASE (194)
         sym_n=24
         ibrav=4
      !Cubic
      CASE (195)
         sym_n=12
         ibrav=1
      CASE (196)
         sym_n=12
         ibrav=2
      CASE (197)
         sym_n=12
         ibrav=3
      CASE (198)
         sym_n=12
         ibrav=1
      CASE (199)
         sym_n=12
         ibrav=3
      CASE (200)
         sym_n=24
         ibrav=1
      CASE (201)
         sym_n=24
         ibrav=1
      CASE (202)
         sym_n=24
         ibrav=2
      CASE (203)
         sym_n=24
         ibrav=2
      CASE (204)
         sym_n=24
         ibrav=3
      CASE (205)
         sym_n=24
         ibrav=1
      CASE (206)
         sym_n=24
         ibrav=3
      CASE (207)
         sym_n=24
         ibrav=1
      CASE (208)
         sym_n=24
         ibrav=1
      CASE (209)
         sym_n=24
         ibrav=2
      CASE (210)
         sym_n=24
         ibrav=2
      CASE (211)
         sym_n=24
         ibrav=3
      CASE (212)
         sym_n=24
         ibrav=1
      CASE (213)
         sym_n=24
         ibrav=1
      CASE (214)
         sym_n=24
         ibrav=3
      CASE (215)
         sym_n=24
         ibrav=1
      CASE (216)
         sym_n=24
         ibrav=2
      CASE (217)
         sym_n=24
         ibrav=3
      CASE (218)
         sym_n=24
         ibrav=1
      CASE (219)
         sym_n=24
         ibrav=2
      CASE (220)
         sym_n=24
         ibrav=3
      CASE (221)
         sym_n=48
         ibrav=1
      CASE (222)
         sym_n=48
         ibrav=1
      CASE (223)
         sym_n=48
         ibrav=1
      CASE (224)
         sym_n=48
         ibrav=1
      CASE (225)
         sym_n=48
         ibrav=2
      CASE (226)
         sym_n=48
         ibrav=2
      CASE (227)
         sym_n=48
         ibrav=2
      CASE (228)
         sym_n=48
         ibrav=2
      CASE (229)
         sym_n=48
         ibrav=3
      CASE (230)
         sym_n=48
         ibrav=3
      END SELECT simmetria
      RETURN
   END SUBROUTINE sym_brav

   SUBROUTINE find_equivalent_tau(space_group_number,inco,outco,i,unique)
      
   !sel_grup ->   input   space_group_number
   !         inco coordinate
   !         i element index
   !      output outco coordinates

      INTEGER, INTENT(IN) :: space_group_number,i
      REAL(DP),dimension(:,:), INTENT(IN) :: inco
      REAL(DP),dimension(:,:,:), INTENT(OUT) :: outco
      character(LEN=1), INTENT(IN) :: unique

      REAL(DP), PARAMETER :: unterz=(1.0_DP)/(3.0_DP)
      REAL(DP), PARAMETER :: duterz=(2.0_DP)/(3.0_DP)
      REAL(DP), PARAMETER :: unsest=(1.0_DP)/(6.0_DP)
      REAL(DP), PARAMETER :: cisest=(5.0_DP)/(6.0_DP)

      INTEGER :: k,j
      simmetria: SELECT CASE (space_group_number)
      !*****************************************
      !Triclinic 1-2
      CASE (1)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
      CASE (2)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         outco(k,2,i)=-inco(k,i)
         END DO
      !*****************************************
      !Monoclinic 3-15
      CASE (3)
         !x,y,z
         !-x,y,-z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO

         IF (unique=='2') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF

         IF (unique=='1') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF
      CASE (4)
         !x,y,z
         !-X,Y+1/2,-Z 
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO

         IF (unique=='2') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)+0.5_DP
         outco(3,2,i)=-inco(3,i)
         END IF

         IF (unique=='1') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         END IF

      CASE (5)
         !X,Y,Z identita
         !-X,Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2         
         IF (unique=='2') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF

         IF (unique=='1') THEN
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF
      CASE (6)
         !ID
         !x,-y,z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         
         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF

         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
      CASE (7)
         !ID
         !x,-y,1/2+z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO

         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         END IF

         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=0.5_DP+inco(2,i)
         outco(3,2,i)=-inco(3,i)         
         END IF
      CASE (8)
         !symmetry= X,Y,Z
         !symmetry= X,-Y,Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2

         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         END IF

         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         END IF
      CASE (9)
         !symmetry= X,Y,Z
         !symmetry= X,-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2

         IF (unique=='2') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         END IF
         
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)+0.5_DP
         outco(3,2,i)=-inco(3,i)
         END IF
      CASE (10)
         !symmetry= X,Y,Z
         !symmetry= X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= -X,-Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO

         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF

         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
      CASE (11)
         !symmetry= X,Y,Z
         !symmetry= -X,1/2+Y,-Z
         !symmetry= -X,-Y,-Z
         !symmetry= X,1/2-Y,Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=0.5_DP+inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=inco(3,i)
         END IF
         
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=inco(2,i)
         outco(3,4,i)=0.5_DP-inco(3,i)
         END IF
      CASE (12)
         !symmetry= X,Y,Z
         !symmetry= X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= -X,-Y,-Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
         
         IF (unique=='1') THEN
         outco(1,2,i)=inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
      CASE (13)
         !symmetry= X,Y,Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= -X,-Y,-Z
         !symmetry= X,-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO

         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=0.5_DP-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=0.5_DP+inco(3,i)
         END IF
         
         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=0.5_DP-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF
      CASE (14)
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,-Z
         !symmetry= -X,1/2+Y,1/2-Z
         !symmetry= X,1/2-Y,1/2+Z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO

         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=0.5_DP+inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=0.5_DP+inco(3,i)
         END IF

         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=0.5_DP-inco(2,i)
         outco(3,3,i)=0.5_DP+inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP+inco(2,i)
         outco(3,4,i)=0.5_DP-inco(3,i)
         END IF
      CASE (15)
         !symmetry= X,Y,Z
         !symmetry= -X,Y,1/2-Z 
         !symmetry= -X,-Y,-Z
         !symmetry= X,-Y,1/2+Z 
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         
         IF (unique=='2') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=inco(2,i)
         outco(3,2,i)=0.5_DP-inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=3
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=0.5_DP+inco(3,i)
         END IF

         IF (unique=='1') THEN
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=0.5_DP-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=3
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=0.5_DP+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         END IF

      !*****************************************
      !Orthorhombic 16-74
      CASE (16) !P222
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= X,-Y,-Z 
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
      CASE (17) !P222(1)
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,1/2+Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= X,-Y,-Z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
      CASE (18) !P2(1)2(1)2
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= 1/2-X,1/2+Y,-Z
         !symmetry= 1/2+X,1/2-Y,-Z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=0.5_DP-inco(1,i)
         outco(2,3,i)=0.5_DP+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=0.5_DP+inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=-inco(3,i)

      CASE (19) !P2(1)2(1)2(1)
         !symmetry= X,Y,Z
         !symmetry= 1/2-X,-Y,1/2+Z
         !symmetry= -X,1/2+Y,1/2-Z
         !symmetry= 1/2+X,1/2-Y,-Z 

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=0.5_DP-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=0.5_DP+inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=0.5_DP+inco(1,i)
         outco(2,4,i)=0.5_DP-inco(2,i)
         outco(3,4,i)=-inco(3,i)
      CASE (20) !C222(1)

         ! symmetry= X,Y,Z
         !symmetry= -X,-Y,1/2+Z
         !symmetry= -X,Y,1/2-Z
         !symmetry= X,-Y,-Z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=0.5_DP+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=0.5_DP-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
      
      CASE (21) !C222
         !symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= X,-Y,-Z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         

      CASE (22) !F222
         ! symmetry= X,Y,Z
         !symmetry= -X,-Y,Z
         !symmetry= -X,Y,-Z
         !symmetry= X,-Y,-Z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
      CASE (23) !I222
         !id
         !-x,-y,z
         !x,,y,-z
         !x,-y,-z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
      CASE (24) !I2(1)2(1)2(1)
         !id
         !-x+1/2,-y,z+1/2
         !-x,1/2+y,1/2-z
         !x+1/2,-y+1/2,-z
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         
      CASE (25) !Pmm2
         !id
         !-x,-y,z
         !+x,-y,+z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
      
      CASE (26) !Pmc2(1)
         !id
         !-x,-y,z+1/2
         !+x,-y,+z+1/2
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (27) !Pcc2
         !id
         !-x,-y,z
         !+x,-y,+z+1/2
         !-x,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (28) !Pma2
         !id
         !-x,-y,z
         !1/2+x,-y,z
         !1/2-x,y,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (29) !Pca2(1)
         !id
         !-x,-y,z+1/2
         !1/2+x,-y,z
         !1/2-x,y,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (30) !Pnc2
         !id
         !-x,-y,z
         !+x,1/2-y,z+1/2
         !-x,y+1/2,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (31) !Pmn2(1)
         !id
         !1/2-x,-y,z+1/2
         !1/2+x,-y,z+1/2
         !-x,y,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
      
      CASE (32) !Pba2
         !id
         !-x,-y,z
         !1/2+x,1/2-y,z
         !1/2-x,1/2+y,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
      
      CASE (33) !Pna2(1)
         !id
         !-x,-y,z+1/2
         !1/2+x,1/2-y,z
         !1/2-x,1/2+y,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (34) !Pnn2
         !id
         !-x,-y,z
         !1/2+x,1/2-y,1/2+z
         !1/2-x,1/2+y,1/2+z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (35) !Cmm2
         !id
         !-x,-y,z
         !+x,-y,z
         !-x,+y,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (36) !Cmc2(1)
         !id
         !-x,-y,z+1/2
         !+x,-y,z+1/2
         !-x,+y,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)
      
      CASE (37) !Ccc2
         !id
         !-x,-y,z
         !+x,-y,z+1/2
         !-x,+y,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (38) !Amm2
         !id
         !-x,-y,z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (39) !Abm2
         !id
         !-x,-y,z
         !x,-y+1/2,z
         !-x,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)

      CASE (40) !Ama2
         !id
         !-x,-y,z
         !x+1/2,-y,z
         !-x+1/2,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (41) !Aba2
         !id
         !-x,-y,z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)

      CASE (42) !Fmm2
         !id
         !-x,-y,z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (43) !Fdd2
         !id
         !-x,-y,z
         !x+1/4,-y+1/4,z+1/4
         !-x+1/4,y+1/4,z+1/4

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.25_DP
         outco(2,3,i)=-inco(2,i)+0.25_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=-inco(1,i)+0.25_DP
         outco(2,4,i)=+inco(2,i)+0.25_DP
         outco(3,4,i)=+inco(3,i)+0.25_DP

      CASE (44) !Imm2
         !id
         !-x,-y,z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (45) !Iba2
         !id
         !-x,-y,z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)

      CASE (46) !Ima2
         !id
         !-x,-y,z
         !x+1/2,-y,z
         !-x+1/2,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(1,i)+0.5_DP
         outco(2,3,i)=-inco(2,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)+0.5_DP
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=+inco(3,i)

      CASE (47) !Pmmm
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (48) !Pnnn
         
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x+1/2,-y+1/2,-z+1/2
         !x+1/2,y+1/2,-z+1/2
         !x+1/2,-y+/2,z+1/2
         !-x+1/2,y+1/2,z+1/2
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x+1/2,+y,-z+1/2
         !+x,-y+1/2,-z+1/2
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x+1/2,-y,z+1/2
         !-x,y+1/2,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF

      CASE (49) !Pccm
         !id
         !-x,-y,z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !x,-y,z+1/2
         !-x,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (50) !Pban
         
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x+1/2,+y,-z
         !+x,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x+1/2,-y,z
         !-x,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)
         END IF

      CASE (51) !Pmma
         !id
         !-x+1/2,-y,z
         !-x,+y,-z
         !+x+1/2,-y,-z
         !-x,-y,-z
         !x+1/2,y,-z
         !x,-y,z
         !-x+1/2,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (52) !Pnna
         !id
         !-x+1/2,-y,z
         !-x+1/2,+y+1/2,-z+1/2
         !+x,-y+1/2,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z
         !x+1/2,-y+1/2,z+1/2
         !-x,y+1/2,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (53) !Pmna
         !id
         !-x+1/2,-y,z+1/2
         !-x+1/2,+y,-z+1/2
         !+x,-y,-z
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x+1/2,-y,z+1/2
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (54) !Pcca
         !id
         !-x+1/2,-y,z
         !-x,+y,-z+1/2
         !+x+1/2,-y,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z
         !x,-y,z+1/2
         !-x+1/2,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (55) !Pbam
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x,y,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (56) !Pccn
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y+1/2,-z+1/2
         !+x+1/2,-y,-z+1/2
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x,-y+1/2,z+1/2
         !-x+1/2,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (57) !Pbcm
         !id
         !-x,-y,z+1/2
         !-x,+y+1/2,-z+1/2
         !+x,-y+1/2,-z
         !-x,-y,-z
         !x,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (58) !Pnnm
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !x+1/2,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (59) !Pmmn

         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y+1/2,-z
         !+x+1/2,-y,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !x,-y+1/2,z
         !-x+1/2,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)         
         END IF

      CASE (60) !Pbcn
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-x,+y,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z+1/2
         !x,-y,z+1/2
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (61) !Pbca
         !id
         !-x+1/2,-y,z+1/2
         !-x,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (62) !Pnma
         !id
         !-x+1/2,-y,z+1/2
         !-x,+y+1/2,-z
         !+x+1/2,-y+1/2,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z
         !-x+1/2,y+1/2,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP
      
      CASE (63) !Cmcm
         !id
         !-x,-y,z+1/2
         !-x,+y,-z+1/2
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z+1/2
         !x,-y,z+1/2
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (64) !Cmca
         !id
         !-x,-y+1/2,z+1/2
         !-x,+y,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z+1/2
         !x,-y,z+1/2
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (65) !Cmmm
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (66) !Cccm
         !id
         !-x,-y,z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !x,-y,z+1/2
         !-x,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (67) !Cmma
         !id
         !-x,-y+1/2,z
         !-x,+y,-z+1/2
         !+x,-y,-z
         !-x,-y,-z
         !x,y+1/2,-z
         !x,-y+1/2,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (68) !Ccca
         
         IF (unique=='1') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y,-z
         !+x+1/2,-y+1/2,-z
         !-x,-y+1/2,-z+1/2
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x+1/2,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-x,+y,-z
         !+x+1/2,-y,-z+1/2
         !-x,-y,-z
         !x+1/2,y,-z
         !x,-y+,z+1/2
         !-x+1/2,y,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP      
         END IF

      CASE (69) !Fmmm
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (70) !Fddd

         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x+1/4,-y+1/4,-z+1/4
         !x+1/4,y+1/4,-z+1/4
         !x+1/4,-y+1/4,z+1/4
         !-x+1/4,y+1/4,z+1/4

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.25_DP
         outco(2,5,i)=-inco(2,i)+0.25_DP
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.25_DP
         outco(2,6,i)=+inco(2,i)+0.25_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(1,i)+0.25_DP
         outco(2,7,i)=-inco(2,i)+0.25_DP
         outco(3,7,i)=+inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.25_DP
         outco(2,8,i)=+inco(2,i)+0.25_DP
         outco(3,8,i)=+inco(3,i)+0.25_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x+3/4,-y+3/4,z
         !-x+3/4,+y,-z+3/4
         !+x,-y+3/4,-z+3/4
         !-x,-y,-z
         !x+3/4,y+3/4,-z
         !x+3/4,-y,z+3/4
         !-x,y+3/4,z+3/4

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.75_DP
         outco(2,2,i)=-inco(2,i)+0.75_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.75_DP
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.75_DP
         outco(3,4,i)=-inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.75_DP
         outco(2,6,i)=+inco(2,i)+0.75_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.75_DP
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)+0.75_DP
         outco(3,8,i)=+inco(3,i)+0.75_DP
         END IF

      CASE (71) !Immm
         !id
         !-x,-y,z
         !-x,+y,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y,-z
         !x,-y,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      CASE (72) !Ibam
         !id
         !-x,-y,z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x,y,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)+0.5_DP
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (73) !Ibca
         !id
         !-x+1/2,-y,z+1/2
         !-x,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !x,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+0.5_DP
         outco(2,8,i)=+inco(2,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (74) !Imma
         !id
         !-x,-y+1/2,z
         !-x,+y+1/2,-z
         !+x,-y,-z
         !-x,-y,-z
         !x,y+1/2,-z
         !x,-y+1/2,z
         !-x,y,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=+inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(1,i)
         outco(2,7,i)=-inco(2,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)
         outco(2,8,i)=+inco(2,i)
         outco(3,8,i)=+inco(3,i)

      !*****************************************
      !Tetragonal 75-142

      CASE (75) !P4
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)

      CASE (76) !P4(1)
         !id
         !-x,-y,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP

      CASE (77) !P4(2)
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP

      CASE (78) !P4(3)
         !id
         !-x,-y,z+1/2
         !-y,x,z+3/4
         !y,-x,z+1/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.25_DP

      CASE (79) !I4
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)

      CASE (80) !I4(1)
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP

      CASE (81) !P-4
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,x,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)

      CASE (82) !I-4
         !id
         !-x,-y,z
         !+y,-x,-z
         !-y,x,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)

      CASE (83) !P4/m
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (84) !P(2)/m
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z+1/2
         !-y,x-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP

      CASE (85) !P4/n

         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z
         !y+1/2,-x+1/2,z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !y,-x,-z
         !-y,x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-y+1/2,x,z
         !y,-x+1/2,z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !y+1/2,-x,-z
         !-y,x+1/2,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         END IF

      CASE (86) !P4(2)/n
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z+1/2
         !y+1/2,-x+1/2,z+1/2
         !-x+1/2,-y+1/2,-z+1/2
         !x+1/2,y+1/2,-z+1/2
         !y,-x,-z
         !-y,x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-y,x+1/2,z+1/2
         !y+1/2,-x,z+1/2
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !y,-x+1/2,-z+1/2
         !-y+1/2,x,-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         END IF

      CASE (87) !I4/m
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=-inco(3,i)
      
      CASE (88) !I4(1)/a
         IF (unique=='1') THEN
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         !-x,-y+1/2,-z+1/4
         !x+1/2,y,-z+3/4
         !y,-x,-z
         !-y+1/2,x+1/2,-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y,z+1/2
         !-y+3/4,x+1/4,z+1/4
         !y+3/4,-x+3/4,z+3/4
         !-x,-y,-z
         !x+1/2,y,-z+1/2
         !y+1/4,-x+3/4,-z+3/4
         !-y+1/4,x+1/4,-z+1/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.75_DP
         outco(2,3,i)=+inco(1,i)+0.25_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.75_DP
         outco(2,4,i)=-inco(1,i)+0.75_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.25_DP
         outco(2,7,i)=-inco(1,i)+0.75_DP
         outco(3,7,i)=-inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.25_DP
         outco(2,8,i)=+inco(1,i)+0.25_DP
         outco(3,8,i)=-inco(3,i)+0.25_DP
         END IF

      CASE (89) !P422
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,+y,-z
         !x,-y,-z
         !y,x,-z
         !-y,-x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (90) !P42(1)2
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z
         !y+1/2,-x+1/2,z
         !-x+1/2,+y+1/2,-z
         !x+1/2,-y+1/2,-z
         !y,x,-z
         !-y,-x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (91) !P4(1)22
         !id
         !-x,-y,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+3/4
         !-y,-x-z+1/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.25_DP
      
      CASE (92) !P4(1)2(1)2
         !id
         !-x,-y,z+1/2
         !-y+1/2,x+1/2,z+1/4
         !y+1/2,-x+1/2,z+3/4
         !-x+1/2,+y+1/2,-z+1/4
         !x+1/2,-y+1/2,-z+3/4
         !y,x,-z
         !-y,-x-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP

      CASE (93) !P4(2)22
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !-x,+y,-z
         !x,-y,-z
         !y,x,-z+1/2
         !-y,-x-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
      
      CASE (94) !P4(2)2(1)2
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z+1/2
         !y+1/2,-x+1/2,z+1/2
         !-x+1/2,+y+1/2,-z+1/2
         !x+1/2,-y+1/2,-z+1/2
         !y,x,-z
         !-y,-x,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (95) !P4(3)22
         !id
         !-x,-y,z+1/2
         !-y,x,z+3/4
         !y,-x,z+1/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+1/4
         !-y,-x-z+3/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.25_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.75_DP

      CASE (96) !P4(3)2(1)2
         !id
         !-x,-y,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+3/4
         !-y,-x-z+1/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.25_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP

      CASE (97) !I422
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !-x,+y,-z
         !x,-y,-z
         !y,x,-z
         !-y,-x-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (98) !I4(1)22
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x,z+1/4
         !y,-x,z+3/4
         !-x,+y,-z
         !x,-y,-z+1/2
         !y,x,-z+3/4
         !-y,-x-z+1/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (99) !P4mm
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x,-y,+z
         !-x,+y,+z
         !-y,-x,+z
         !y,x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)

      CASE (100) !P4bm
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z
         !-y+1/2,-x+1/2,+z
         !y+1/2,x+1/2,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (101) !P4(2)cm
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !-y,-x,+z
         !y,x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)

      CASE (102) !P4(2)nm
         !id
         !-x,-y,z
         !-y+1/2,x+1/2,z+1/2
         !y+1/2,-x+1/2,z+1/2
         !+x+1/2,-y+1/2,+z+1/2
         !-x+1/2,+y+1/2,+z+1/2
         !-y,-x,+z
         !y,x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)

      CASE (103) !P4cc
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !-y,-x,+z+1/2
         !y,x,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (104) !P4nc
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x+1/2,-y+1/2,+z+1/2
         !-x+1/2,+y+1/2,+z+1/2
         !-y+1/2,-x+1/2,+z+1/2
         !y+1/2,x+1/2,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (105) !P4(2)mc
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !+x,-y,+z
         !-x,+y,+z
         !-y,-x,+z+1/2
         !y,x,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (106) !P4(2)bc
         !id
         !-x,-y,z
         !-y,x,z+1/2
         !y,-x,z+1/2
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z+1/2
         !-y+1/2,-x+1/2,+z+1/2
         !y+1/2,x+1/2,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (107) !I4mm
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x,-y,+z
         !-x,+y,+z
         !-y,-x,+z
         !y,x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)

      CASE (108) !I4cm
         !id
         !-x,-y,z
         !-y,x,z
         !y,-x,z
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z
         !-y+1/2,-x+1/2,+z
         !y+1/2,x+1/2,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (109) !I4(1)md
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         !+x,-y,+z
         !-x+1/2,+y+1/2,+z+1/2
         !-y,-x+1/2,+z+1/2
         !y+1/2,x,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.75_DP

      CASE (110) !I4(1)cd
         !id
         !-x+1/2,-y+1/2,z+1/2
         !-y,x+1/2,z+1/4
         !y+1/2,-x,z+3/4
         !+x,-y,+z+1/2
         !-x+1/2,+y+1/2,+z
         !-y,-x+1/2,+z+3/4
         !y+1/2,x,z+1/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.25_DP

      CASE (111) !P-42m
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x,+y,-z
         !+x,-y,-z
         !-y,-x,+z
         !y,x,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)

      CASE (112) !P-42c
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !-y,-x,+z+1/2
         !y,x,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (113) !P-42(1)m
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !-y+1/2,-x+1/2,+z
         !y+1/2,x+1/2,z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)

      CASE (114) !P-42(1)c
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x+1/2,+y+1/2,-z+1/2
         !+x+1/2,-y+1/2,-z+1/2
         !-y+1/2,-x+1/2,+z+1/2
         !y+1/2,x+1/2,z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)+0.5_DP
         outco(3,8,i)=+inco(3,i)+0.5_DP

      CASE (115) !P-4m2
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z
         !-x,+y,+z
         !+y,+x,-z
         !-y,-x,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (116) !P-4c2
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !+y,+x,-z+1/2
         !-y,-x,-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP

      CASE (117) !P-4b2
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x+1/2,-y+1/2,+z
         !-x+1/2,+y+1/2,+z
         !+y+1/2,+x+1/2,-z
         !-y+1/2,-x+1/2,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)

      CASE (118) !P-4n2
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x+1/2,-y+1/2,+z+1/2
         !-x+1/2,+y+1/2,+z+1/2
         !+y+1/2,+x+1/2,-z+1/2
         !-y+1/2,-x+1/2,-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)+0.5_DP
         outco(2,5,i)=-inco(2,i)+0.5_DP
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=+inco(2,i)+0.5_DP
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP

      CASE (119) !I-4m2
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z
         !-x,+y,+z
         !+y,+x,-z
         !-y,-x,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)

      CASE (120) !I-4c2
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z+1/2
         !-x,+y,+z+1/2
         !+y,+x,-z+1/2
         !-y,-x,-z+1/2
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP

      CASE (121) !I-42m
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !+x,-y,+z
         !-x,+y,+z
         !+y,+x,-z
         !-y,-x,-z
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)

      CASE (122) !I-42d
         !id
         !-x,-y,z
         !y,-x,-z
         !-y,+x,-z
         !-x+1/2,+y,-z+3/4
         !+x+1/2,-y,-z+3/4
         !-y+1/2,-x,+z+3/4
         !+y+1/2,+x,+z+3/4
         
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=-inco(2,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=+inco(2,i)+0.5_DP
         outco(2,8,i)=+inco(1,i)
         outco(3,8,i)=+inco(3,i)+0.75_DP

      CASE (123) !P4/mmm
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)

      CASE (124) !P4/mcc
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z+1/2
         !+x,-y,-z+1/2
         !+y,+x,-z+1/2
         !-y,-x,-z+1/2
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z+1/2
         !-x,y,z+1/2
         !-y,-x,z+1/2
         !y,x,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP

      CASE (125) !P4/nbm
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !y+1/2,-x+1/2,-z
         !-y+1/2,x+1/2,-z
         !x+1/2,-y+1/2,z
         !-x+1/2,y+1/2,z
         !-y+1/2,-x+1/2,z
         !y+1/2,x+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)+0.5_DP
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x+1/2,-y+1/2,z
         !-y+1/2,+x,+z
         !+y,-x+1/2,+z
         !-x+1/2,+y,-z
         !+x,-y+1/2,-z
         !+y,+x,-z
         !-y+1/2,-x+1/2,-z
         !-x,-y,-z
         !x+1/2,y+1/2,-z
         !y+1/2,-x,-z
         !-y,x+1/2,-z
         !x+1/2,-y,z
         !-x,y+1/2,z
         !-y,-x,z
         !y+1/2,x+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF

      CASE (126) !P4/nnc
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x+1/2,-y+1/2,-z+1/2
         !x+1/2,y+1/2,-z+1/2
         !y+1/2,-x+1/2,-z+1/2
         !-y+1/2,x+1/2,-z+1/2
         !x+1/2,-y+1/2,z+1/2
         !-x+1/2,y+1/2,z+1/2
         !-y+1/2,-x+1/2,z+1/2
         !y+1/2,x+1/2,z+1/2

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)+0.5_DP
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

      CASE (127) !P4/mbm
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)

      CASE (128) !P4/mnc
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP

      CASE (129)
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y+1/2,+x+1/2,+z
         !+y+1/2,-x+1/2,+z
         !-x+1/2,+y+1/2,-z
         !+x+1/2,-y+1/2,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x+1/2,-y+1/2,-z
         !x+1/2,y+1/2,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y+1/2,-x+1/2,z
         !y+1/2,x+1/2,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         END IF

      CASE (130) !P4/ncc
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF
      
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

      CASE (131) !P4(2)/mmc
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP

      CASE (132) !P4(2)mcm
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)

      CASE (133) !P4(2)/nbc
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

      CASE (134) !P4(2)/nnm
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF

      CASE (135) !P4(2)/mbc
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.5_DP
         outco(2,8,i)=-inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP

      CASE (136) !P4(2)mnm
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)+0.5_DP
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)

      CASE (137) !P4(2)/nmc
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP
         END IF

      CASE (138) !P4(2)/ncm
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)+0.5_DP
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         END IF
         
         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)+0.5_DP
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)+0.5_DP
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)+0.5_DP
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         END IF

      CASE (139) !I4/mmm
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)

      CASE (140) !I4/mcm
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=+inco(1,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=+inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)+0.5_DP

      CASE (141) !I4(1)amd
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.75_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.25_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.25_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.75_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.25_DP
         outco(2,3,i)=+inco(1,i)+0.75_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.25_DP
         outco(2,4,i)=-inco(1,i)+0.25_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)+0.25_DP
         outco(2,7,i)=+inco(1,i)+0.75_DP
         outco(3,7,i)=-inco(3,i)+0.25_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.25_DP
         outco(2,8,i)=-inco(1,i)+0.25_DP
         outco(3,8,i)=-inco(3,i)+0.75_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.75_DP
         outco(2,11,i)=-inco(1,i)+0.25_DP
         outco(3,11,i)=-inco(3,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=+inco(1,i)+0.75_DP
         outco(3,12,i)=-inco(3,i)+0.25_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)
         !S=15
         outco(1,15,i)=-inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=+inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         END IF

      CASE (142) !I4(1)/acd
         IF (unique=='1') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)
         outco(2,3,i)=+inco(1,i)+0.5_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.5_DP
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.25_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)+0.5_DP
         outco(3,6,i)=-inco(3,i)+0.75_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.5_DP
         outco(2,7,i)=+inco(1,i)+0.5_DP
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(2,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)+0.5_DP
         outco(3,9,i)=-inco(3,i)+0.25_DP
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.75_DP
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=-inco(1,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=+inco(1,i)+0.5_DP
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=+inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=+inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.75_DP
         END IF

         IF (unique=='2') THEN
         !id
         !-x,-y,z
         !-y,+x,+z
         !+y,-x,+z
         !-x,+y,-z
         !+x,-y,-z
         !+y,+x,-z
         !-y,-x,-z
         !-x,-y,-z
         !x,y,-z
         !y,-x,-z
         !-y,x,-z
         !x,-y,z
         !-x,y,z
         !-y,-x,z
         !y,x,z

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=+inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(2,i)+0.25_DP
         outco(2,3,i)=+inco(1,i)+0.75_DP
         outco(3,3,i)=+inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=+inco(2,i)+0.25_DP
         outco(2,4,i)=-inco(1,i)+0.25_DP
         outco(3,4,i)=+inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+0.5_DP
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)+0.25_DP
         outco(2,7,i)=+inco(1,i)+0.75_DP
         outco(3,7,i)=-inco(3,i)+0.75_DP
         !S=8
         outco(1,8,i)=-inco(2,i)+0.25_DP
         outco(2,8,i)=-inco(1,i)+0.25_DP
         outco(3,8,i)=-inco(3,i)+0.25_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)+0.5_DP
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.75_DP
         outco(2,11,i)=-inco(1,i)+0.25_DP
         outco(3,11,i)=-inco(3,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=+inco(1,i)+0.75_DP
         outco(3,12,i)=-inco(3,i)+0.25_DP
         !S=13
         outco(1,13,i)=+inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=+inco(3,i)
         !S=14
         outco(1,14,i)=-inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=+inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=-inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=+inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=+inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.75_DP
         END IF

      !*****************************************
      !Trigonal 143-167
      
      CASE (143) !P3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)

      CASE (144) !P3(1)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz

      CASE (145) !P3(2)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz

      CASE (146) !R3
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=+inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=inco(3,i)
         outco(3,3,i)=+inco(1,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         END IF

      CASE (147) !P-3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=-inco(3,i)
      CASE (148) !R-3
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=+inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(3,i)
         outco(2,5,i)=-inco(1,i)
         outco(3,5,i)=-inco(2,i)
         !S=6
         outco(1,6,i)=-inco(2,i)
         outco(2,6,i)=-inco(3,i)
         outco(3,6,i)=-inco(1,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=-inco(3,i)
         END IF
      CASE (149) !P312
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)

      CASE (150) !P321
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)

      CASE (151) !P3(1)12
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)+duterz
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+unterz
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)

      CASE (152) !P3(1)21
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+duterz
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+unterz

      CASE (153) !P3(2)12
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)+unterz
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+duterz
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)
      
      CASE (154) !P3(2)21
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)+inco(1,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+unterz
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+duterz

      CASE (155) !R32
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=inco(2,i)
         outco(2,3,i)=inco(3,i)
         outco(3,3,i)=inco(1,i)
         !S=4
         outco(1,4,i)=-inco(3,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(1,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=-inco(1,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(3,i)
         outco(3,6,i)=-inco(2,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         END IF

      CASE (156) !P3m1
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !s=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)
         outco(2,6,i)=inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)

      CASE (157) !P31m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=+inco(3,i)

      CASE (158) !P3c1
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !s=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)-inco(1,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)
         outco(2,6,i)=inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP

      CASE (159) !P31c
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP

      CASE (160) !R3m
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=inco(3,i)
         outco(2,4,i)=inco(2,i)
         outco(3,4,i)=inco(1,i)
         !S=5
         outco(1,5,i)=inco(2,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)
         outco(2,6,i)=inco(3,i)
         outco(3,6,i)=inco(2,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)
         END IF
         
      CASE (161) !R3c
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=inco(3,i)+0.5_DP
         outco(2,4,i)=inco(2,i)+0.5_DP
         outco(3,4,i)=inco(1,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(2,i)+0.5_DP
         outco(2,5,i)=inco(1,i)+0.5_DP
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)+0.5_DP
         outco(2,6,i)=inco(3,i)+0.5_DP
         outco(3,6,i)=inco(2,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         END IF

      CASE (162) !P-31m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(2,i)
         outco(2,10,i)=+inco(1,i)
         outco(3,10,i)=+inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)
      
      CASE (163) !P-31c
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(2,i)
         outco(2,4,i)=-inco(1,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(1,i)+inco(2,i)
         outco(2,5,i)=+inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)
         outco(2,6,i)=+inco(1,i)-inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(2,i)
         outco(2,10,i)=+inco(1,i)
         outco(3,10,i)=+inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP

      CASE (164) !P-3m1
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)

      CASE (165) !P-3c1
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP

      CASE (166) !R-3m
         IF (unique=='1') THEN
         !Rhombohedral
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=-inco(3,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(1,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=-inco(1,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(3,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=-inco(2,i)
         outco(2,9,i)=-inco(3,i)
         outco(3,9,i)=-inco(1,i)
         !S=10
         outco(1,10,i)=inco(3,i)
         outco(2,10,i)=inco(2,i)
         outco(3,10,i)=inco(1,i)
         !S=11
         outco(1,11,i)=+inco(2,i)
         outco(2,11,i)=+inco(1,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(3,i)
         outco(3,12,i)=+inco(2,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)
         END IF

      CASE (167) !R-3c
         IF (unique=='1') THEN
         !Rhombohedral
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=inco(3,i)
         outco(2,2,i)=inco(1,i)
         outco(3,2,i)=inco(2,i)
         !S=3
         outco(1,3,i)=+inco(2,i)
         outco(2,3,i)=+inco(3,i)
         outco(3,3,i)=+inco(1,i)
         !S=4
         outco(1,4,i)=-inco(3,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(1,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(2,i)+0.5_DP
         outco(2,5,i)=-inco(1,i)+0.5_DP
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+0.5_DP
         outco(2,6,i)=-inco(3,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=-inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=-inco(2,i)
         outco(2,9,i)=-inco(3,i)
         outco(3,9,i)=-inco(1,i)
         !S=10
         outco(1,10,i)=inco(3,i)+0.5_DP
         outco(2,10,i)=inco(2,i)+0.5_DP
         outco(3,10,i)=inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=+inco(2,i)+0.5_DP
         outco(2,11,i)=+inco(1,i)+0.5_DP
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)+0.5_DP
         outco(2,12,i)=+inco(3,i)+0.5_DP
         outco(3,12,i)=+inco(2,i)+0.5_DP
         END IF
   
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(2,i)
         outco(2,4,i)=+inco(1,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(1,i)-inco(2,i)
         outco(2,5,i)=-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)
         outco(2,6,i)=-inco(1,i)+inco(2,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=+inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
         END IF

      !*****************************************
      !Exagonal 168-194
      CASE (168) !P6
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)

      CASE (169) !P6(1)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+cisest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unsest

      CASE (170) !P6(5)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unsest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+cisest

      CASE (171) !P6(2)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+duterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unterz

      CASE (172) !P6(4)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+duterz

      CASE (173) !P6(3)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP

      CASE (174) !P-6
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)

      CASE (175) !P6/m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(2,i)
         outco(2,11,i)=+inco(1,i)-inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)+inco(2,i)
         outco(2,12,i)=-inco(1,i)
         outco(3,12,i)=-inco(3,i)
      
      CASE (176) !P6(3)/m
                  DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(1,i)
         outco(2,7,i)=-inco(2,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(2,i)
         outco(2,8,i)=-inco(1,i)+inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)-inco(2,i)
         outco(2,9,i)=inco(1,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=+inco(1,i)
         outco(2,10,i)=+inco(2,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(2,i)
         outco(2,11,i)=+inco(1,i)-inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)+inco(2,i)
         outco(2,12,i)=-inco(1,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP

      CASE (177) !P622
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)

      CASE (178) !P(1)22
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+cisest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unsest
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+unterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+duterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+cisest
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+unsest

      CASE (179) !P6(5)22
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unsest
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+cisest
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+duterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+unterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+unsest
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+cisest

      CASE (180) !P6(2)22
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+duterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+unterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+duterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+unterz
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+duterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+unterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+duterz
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+unterz

      CASE (181) !P6(4)22
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)+unterz
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)+duterz
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+unterz
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+duterz
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)+unterz
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+duterz
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+unterz
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+duterz

      CASE (182) !6(3)22
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=inco(2,i)
         outco(2,7,i)=inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP

      CASE (183) !P6mm
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)

      CASE (184) !P6cc
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP
      
      CASE (185) !P6(3)cm
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)

      CASE (186) !P(3)mc
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=+inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=+inco(3,i)+0.5_DP

      CASE (187) !P-6m2
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
      
      CASE (188) !P-6c2
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(2,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=+inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(1,i)+inco(2,i)
         outco(2,8,i)=inco(2,i)
         outco(3,8,i)=inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(1,i)
         outco(2,9,i)=inco(1,i)-inco(2,i)
         outco(3,9,i)=inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=+inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)

      CASE (189) !P-62m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=inco(3,i)
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=inco(3,i)

      CASE (190) !P-62c
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=+inco(1,i)
         outco(2,4,i)=+inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=-inco(2,i)
         outco(2,5,i)=+inco(1,i)-inco(2,i)
         outco(3,5,i)=-inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=-inco(1,i)+inco(2,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=inco(2,i)
         outco(2,10,i)=inco(1,i)
         outco(3,10,i)=inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(1,i)-inco(2,i)
         outco(2,11,i)=-inco(2,i)
         outco(3,11,i)=inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(1,i)
         outco(2,12,i)=-inco(1,i)+inco(2,i)
         outco(3,12,i)=inco(3,i)+0.5_DP

      CASE (191) !P6/mmm
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)

      CASE (192) !P6/mmc
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)+0.5_DP
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)+0.5_DP
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)+0.5_DP

      CASE (193) !P6(3)/mcm
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)+0.5_DP
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)+0.5_DP
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)+0.5_DP
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)+0.5_DP
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)

      CASE (194)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(2,i)
         outco(2,2,i)=inco(1,i)-inco(2,i)
         outco(3,2,i)=+inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+inco(2,i)
         outco(2,3,i)=-inco(1,i)
         outco(3,3,i)=+inco(3,i)
         !S=4
         outco(1,4,i)=-inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=+inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=+inco(2,i)
         outco(2,5,i)=-inco(1,i)+inco(2,i)
         outco(3,5,i)=+inco(3,i)+0.5_DP
         !S=6
         outco(1,6,i)=+inco(1,i)-inco(2,i)
         outco(2,6,i)=+inco(1,i)
         outco(3,6,i)=+inco(3,i)+0.5_DP
         !S=7
         outco(1,7,i)=+inco(2,i)
         outco(2,7,i)=+inco(1,i)
         outco(3,7,i)=-inco(3,i)
         !S=8
         outco(1,8,i)=+inco(1,i)-inco(2,i)
         outco(2,8,i)=-inco(2,i)
         outco(3,8,i)=-inco(3,i)
         !S=9
         outco(1,9,i)=-inco(1,i)
         outco(2,9,i)=-inco(1,i)+inco(2,i)
         outco(3,9,i)=-inco(3,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=-inco(1,i)
         outco(3,10,i)=-inco(3,i)+0.5_DP
         !S=11
         outco(1,11,i)=-inco(1,i)+inco(2,i)
         outco(2,11,i)=+inco(2,i)
         outco(3,11,i)=-inco(3,i)+0.5_DP
         !S=12
         outco(1,12,i)=inco(1,i)
         outco(2,12,i)=+inco(1,i)-inco(2,i)
         outco(3,12,i)=-inco(3,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=inco(2,i)
         outco(2,14,i)=-inco(1,i)+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=+inco(1,i)-inco(2,i)
         outco(2,15,i)=+inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=+inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(2,i)
         outco(2,17,i)=+inco(1,i)-inco(2,i)
         outco(3,17,i)=-inco(3,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+inco(2,i)
         outco(2,18,i)=-inco(1,i)
         outco(3,18,i)=-inco(3,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(2,i)
         outco(2,19,i)=-inco(1,i)
         outco(3,19,i)=+inco(3,i)
         !S=20
         outco(1,20,i)=-inco(1,i)+inco(2,i)
         outco(2,20,i)=+inco(2,i)
         outco(3,20,i)=+inco(3,i)
         !S=21
         outco(1,21,i)=inco(1,i)
         outco(2,21,i)=+inco(1,i)-inco(2,i)
         outco(3,21,i)=+inco(3,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=inco(1,i)
         outco(3,22,i)=inco(3,i)+0.5_DP
         !S=23
         outco(1,23,i)=inco(1,i)-inco(2,i)
         outco(2,23,i)=-inco(2,i)
         outco(3,23,i)=inco(3,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(1,i)
         outco(2,24,i)=-inco(1,i)+inco(2,i)
         outco(3,24,i)=+inco(3,i)+0.5_DP
      !*****************************************
      !Cubic 195-230
      CASE (195) !P23
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)

      CASE (196) !F23
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
      
      CASE (197) !I23
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)

      CASE (198) !P2(1)3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP

      CASE (199) !I2(1)3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP

      CASE (200) !Pm-3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)

      CASE(201) !Pn-3
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)+0.5_DP
         outco(2,13,i)=-inco(2,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(1,i)+0.5_DP
         outco(2,15,i)=-inco(2,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.5_DP
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(3,i)+0.5_DP
         outco(2,17,i)=-inco(1,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(3,i)+0.5_DP
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(3,i)+0.5_DP
         outco(2,20,i)=-inco(1,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)+0.5_DP
         outco(2,21,i)=-inco(3,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(2,i)+0.5_DP
         outco(2,22,i)=-inco(3,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.5_DP
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)+0.5_DP
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)+0.5_DP
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)+0.5_DP
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)
         END IF

      CASE (202) !Fm-3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)

      CASE (203) !Fd-3
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)+0.25_DP
         outco(2,13,i)=-inco(2,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=+inco(1,i)+0.25_DP
         outco(2,14,i)=+inco(2,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(1,i)+0.25_DP
         outco(2,15,i)=-inco(2,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.25_DP
         outco(2,16,i)=+inco(2,i)+0.25_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=-inco(3,i)+0.25_DP
         outco(2,17,i)=-inco(1,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(3,i)+0.25_DP
         outco(2,18,i)=+inco(1,i)+0.25_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.25_DP
         outco(2,19,i)=+inco(1,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(3,i)+0.25_DP
         outco(2,20,i)=-inco(1,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.25_DP
         !S=21
         outco(1,21,i)=-inco(2,i)+0.25_DP
         outco(2,21,i)=-inco(3,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(2,i)+0.25_DP
         outco(2,22,i)=-inco(3,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.25_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.25_DP
         outco(2,23,i)=inco(3,i)+0.25_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.25_DP
         outco(2,24,i)=+inco(3,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
         END IF
         
         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.75_DP
         outco(2,2,i)=-inco(2,i)+0.75_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.75_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.75_DP
         outco(3,4,i)=-inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.75_DP
         outco(3,6,i)=-inco(2,i)+0.75_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.75_DP
         outco(2,7,i)=-inco(1,i)+0.75_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.75_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.75_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.75_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.75_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.75_DP
         outco(3,11,i)=-inco(1,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=-inco(3,i)+0.75_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.25_DP
         outco(2,14,i)=+inco(2,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)+0.25_DP
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)+0.25_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)+0.25_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=+inco(3,i)+0.25_DP
         outco(2,19,i)=+inco(1,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)+0.25_DP
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)+0.25_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)+0.25_DP
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)+0.25_DP
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)+0.25_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=+inco(2,i)+0.25_DP
         outco(2,24,i)=+inco(3,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)
         END IF

      CASE (204) !Im-3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(1,i)
         outco(2,16,i)=+inco(2,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)
         outco(2,18,i)=+inco(1,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(2,i)
         outco(2,23,i)=inco(3,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)

      CASE (205) !Pa-3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.5_DP
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)+0.5_DP
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.5_DP
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)+0.5_DP

      CASE (206) !Ia-3
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=-inco(1,i)
         outco(2,13,i)=-inco(2,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=+inco(1,i)+0.5_DP
         outco(2,14,i)=+inco(2,i)
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(1,i)
         outco(2,15,i)=-inco(2,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(1,i)+0.5_DP
         outco(2,16,i)=+inco(2,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=-inco(3,i)
         outco(2,17,i)=-inco(1,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(3,i)+0.5_DP
         outco(2,18,i)=+inco(1,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=+inco(3,i)+0.5_DP
         outco(2,19,i)=+inco(1,i)
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(3,i)
         outco(2,20,i)=-inco(1,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=-inco(2,i)
         outco(2,21,i)=-inco(3,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(2,i)
         outco(2,22,i)=-inco(3,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(2,i)+0.5_DP
         outco(2,23,i)=inco(3,i)+0.5_DP
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=+inco(2,i)+0.5_DP
         outco(2,24,i)=+inco(3,i)
         outco(3,24,i)=-inco(1,i)+0.5_DP

      CASE (207) !P432
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)

      CASE (208) !P4(2)32
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=+inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=+inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP

      CASE (209) !F432
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)

      CASE (210) !F4(1)32
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP

      CASE (211) !I432
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=+inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=+inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)

      CASE (212) !P4(3)32
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.25_DP
         outco(2,13,i)=inco(1,i)+0.75_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.25_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.25_DP
         outco(3,16,i)=+inco(3,i)+0.75_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.25_DP
         outco(2,17,i)=+inco(3,i)+0.75_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.25_DP
         outco(3,18,i)=+inco(2,i)+0.75_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.75_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.25_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.25_DP
         outco(2,21,i)=inco(2,i)+0.75_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.75_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.25_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.25_DP
         outco(3,23,i)=inco(1,i)+0.75_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP

      CASE (213) !P4(1)32
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.75_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=-inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.25_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.25_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.75_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=-inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.25_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.75_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=-inco(1,i)+0.75_DP

      CASE (214) !I4(1)32
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.75_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=-inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.25_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=+inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.25_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=+inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.75_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=-inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.25_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.75_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=-inco(1,i)+0.75_DP

      CASE (215) !P-43m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=-inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=-inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=-inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=-inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=inco(1,i)

      CASE (216) !F-43m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=-inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=-inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=-inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=-inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=inco(1,i)

      CASE (217) !I-43m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=-inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=-inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=-inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=-inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=-inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=-inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=inco(1,i)

      CASE (218) !P-43n
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=-inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=-inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=-inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=-inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=-inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=inco(1,i)+0.5_DP

      CASE (219) !F-43c
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=-inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=-inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=-inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=-inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=-inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=-inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=inco(1,i)+0.5_DP

      CASE (220) !I-43d
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.25_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.75_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=-inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=-inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.25_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=-inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.75_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=-inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.25_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.75_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=-inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=-inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=inco(1,i)+0.75_DP

      CASE (221) !Pm-3m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)

      CASE (222) !Pn-3n
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)+0.5_DP
         outco(2,25,i)=-inco(2,i)+0.5_DP
         outco(3,25,i)=-inco(3,i)+0.5_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.5_DP
         outco(2,29,i)=-inco(1,i)+0.5_DP
         outco(3,29,i)=-inco(2,i)+0.5_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.5_DP
         outco(2,33,i)=-inco(3,i)+0.5_DP
         outco(3,33,i)=-inco(1,i)+0.5_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF

      CASE (223) !Pm-3n
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP

      CASE (224) !Pn-3m
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)+0.5_DP
         outco(2,25,i)=-inco(2,i)+0.5_DP
         outco(3,25,i)=-inco(3,i)+0.5_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.5_DP
         outco(2,29,i)=-inco(1,i)+0.5_DP
         outco(3,29,i)=-inco(2,i)+0.5_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.5_DP
         outco(2,33,i)=-inco(3,i)+0.5_DP
         outco(3,33,i)=-inco(1,i)+0.5_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)+0.5_DP
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)+0.5_DP
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)+0.5_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)+0.5_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)+0.5_DP
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)+0.5_DP
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.5_DP
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)+0.5_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)+0.5_DP
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF

      CASE (225) !Fm-3m
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)

      CASE (226) !Fm-3c
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)+0.5_DP
         outco(2,13,i)=inco(1,i)+0.5_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.5_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.5_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.5_DP
         outco(3,16,i)=inco(3,i)+0.5_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.5_DP
         outco(2,17,i)=+inco(3,i)+0.5_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.5_DP
         outco(3,18,i)=inco(2,i)+0.5_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.5_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.5_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.5_DP
         outco(2,21,i)=inco(2,i)+0.5_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.5_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.5_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.5_DP
         outco(3,23,i)=inco(1,i)+0.5_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP

      CASE (227) !Fd-3m
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
         !S=25
         outco(1,25,i)=-inco(1,i)+0.25_DP
         outco(2,25,i)=-inco(2,i)+0.25_DP
         outco(3,25,i)=-inco(3,i)+0.25_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.25_DP
         outco(2,26,i)=inco(2,i)+0.75_DP
         outco(3,26,i)=-inco(3,i)+0.75_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.75_DP
         outco(2,27,i)=-inco(2,i)+0.75_DP
         outco(3,27,i)=inco(3,i)+0.25_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.75_DP
         outco(2,28,i)=inco(2,i)+0.25_DP
         outco(3,28,i)=inco(3,i)+0.75_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.25_DP
         outco(2,29,i)=-inco(1,i)+0.25_DP
         outco(3,29,i)=-inco(2,i)+0.25_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.75_DP
         outco(2,30,i)=inco(1,i)+0.25_DP
         outco(3,30,i)=inco(2,i)+0.75_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.25_DP
         outco(2,31,i)=inco(1,i)+0.75_DP
         outco(3,31,i)=-inco(2,i)+0.75_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.75_DP
         outco(2,32,i)=-inco(1,i)+0.75_DP
         outco(3,32,i)=inco(2,i)+0.25_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.25_DP
         outco(2,33,i)=-inco(3,i)+0.25_DP
         outco(3,33,i)=-inco(1,i)+0.25_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.75_DP
         outco(2,34,i)=-inco(3,i)+0.75_DP
         outco(3,34,i)=inco(1,i)+0.25_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.75_DP
         outco(2,35,i)=inco(3,i)+0.25_DP
         outco(3,35,i)=inco(1,i)+0.75_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.25_DP
         outco(2,36,i)=inco(3,i)+0.75_DP
         outco(3,36,i)=-inco(1,i)+0.75_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.5_DP
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.5_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.5_DP
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)+0.5_DP
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.5_DP
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.5_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.5_DP
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.5_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.5_DP
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.75_DP
         outco(2,2,i)=-inco(2,i)+0.25_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.25_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.75_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.75_DP
         outco(3,4,i)=-inco(3,i)+0.25_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.75_DP
         outco(3,6,i)=-inco(2,i)+0.25_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.75_DP
         outco(2,7,i)=-inco(1,i)+0.25_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.25_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.75_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.25_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.75_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.75_DP
         outco(3,11,i)=-inco(1,i)+0.25_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.75_DP
         outco(2,12,i)=-inco(3,i)+0.25_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.5_DP
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.5_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.5_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.5_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.5_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.5_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.5_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.5_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.5_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.25_DP
         outco(2,26,i)=inco(2,i)+0.75_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.75_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.25_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.25_DP
         outco(3,28,i)=inco(3,i)+0.75_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.25_DP
         outco(3,30,i)=inco(2,i)+0.75_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.25_DP
         outco(2,31,i)=inco(1,i)+0.75_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.75_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.25_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.75_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.25_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.25_DP
         outco(3,35,i)=inco(1,i)+0.75_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.25_DP
         outco(2,36,i)=inco(3,i)+0.75_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.25_DP
         outco(2,37,i)=-inco(1,i)+0.75_DP
         outco(3,37,i)=inco(3,i)+0.5_DP
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)+0.75_DP
         outco(2,39,i)=inco(1,i)+0.5_DP
         outco(3,39,i)=-inco(3,i)+0.25_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.5_DP
         outco(2,40,i)=-inco(1,i)+0.25_DP
         outco(3,40,i)=-inco(3,i)+0.75_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.25_DP
         outco(2,41,i)=-inco(3,i)+0.75_DP
         outco(3,41,i)=+inco(2,i)+0.5_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.5_DP
         outco(2,42,i)=-inco(3,i)+0.25_DP
         outco(3,42,i)=-inco(2,i)+0.75_DP
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)+0.75_DP
         outco(2,44,i)=+inco(3,i)+0.5_DP
         outco(3,44,i)=-inco(2,i)+0.25_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.25_DP
         outco(2,45,i)=-inco(2,i)+0.75_DP
         outco(3,45,i)=+inco(1,i)+0.5_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.75_DP
         outco(2,46,i)=inco(2,i)+0.5_DP
         outco(3,46,i)=-inco(1,i)+0.25_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.5_DP
         outco(2,47,i)=-inco(2,i)+0.25_DP
         outco(3,47,i)=-inco(1,i)+0.75_DP
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)
         END IF

      CASE (228) !Fd-3c
         IF (unique=='1') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)+0.5_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.5_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)+0.5_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)+0.5_DP
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)+0.5_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.5_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.5_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)+0.5_DP
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)+0.5_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.75_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.25_DP
         outco(2,14,i)=-inco(1,i)+0.25_DP
         outco(3,14,i)=-inco(3,i)+0.25_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.75_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.75_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.75_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.75_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.25_DP
         outco(2,19,i)=-inco(3,i)+0.25_DP
         outco(3,19,i)=-inco(2,i)+0.25_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.75_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.75_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.75_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.75_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.25_DP
         outco(2,24,i)=-inco(2,i)+0.25_DP
         outco(3,24,i)=-inco(1,i)+0.25_DP
         !S=25
         outco(1,25,i)=-inco(1,i)+0.75_DP
         outco(2,25,i)=-inco(2,i)+0.75_DP
         outco(3,25,i)=-inco(3,i)+0.75_DP
         !S=26
         outco(1,26,i)=inco(1,i)+0.75_DP
         outco(2,26,i)=inco(2,i)+0.25_DP
         outco(3,26,i)=-inco(3,i)+0.25_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.25_DP
         outco(2,27,i)=-inco(2,i)+0.25_DP
         outco(3,27,i)=inco(3,i)+0.75_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.25_DP
         outco(2,28,i)=inco(2,i)+0.75_DP
         outco(3,28,i)=inco(3,i)+0.25_DP
         !S=29
         outco(1,29,i)=-inco(3,i)+0.75_DP
         outco(2,29,i)=-inco(1,i)+0.75_DP
         outco(3,29,i)=-inco(2,i)+0.75_DP
         !S=30
         outco(1,30,i)=-inco(3,i)+0.25_DP
         outco(2,30,i)=inco(1,i)+0.75_DP
         outco(3,30,i)=inco(2,i)+0.25_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.75_DP
         outco(2,31,i)=inco(1,i)+0.25_DP
         outco(3,31,i)=-inco(2,i)+0.25_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.25_DP
         outco(2,32,i)=-inco(1,i)+0.25_DP
         outco(3,32,i)=inco(2,i)+0.75_DP
         !S=33
         outco(1,33,i)=-inco(2,i)+0.75_DP
         outco(2,33,i)=-inco(3,i)+0.75_DP
         outco(3,33,i)=-inco(1,i)+0.75_DP
         !S=34
         outco(1,34,i)=inco(2,i)+0.25_DP
         outco(2,34,i)=-inco(3,i)+0.25_DP
         outco(3,34,i)=inco(1,i)+0.75_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.25_DP
         outco(2,35,i)=inco(3,i)+0.75_DP
         outco(3,35,i)=inco(1,i)+0.25_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.75_DP
         outco(2,36,i)=inco(3,i)+0.25_DP
         outco(3,36,i)=-inco(1,i)+0.25_DP
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)+0.5_DP
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.5_DP
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)+0.5_DP
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)+0.5_DP
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)+0.5_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.5_DP
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)+0.5_DP
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)+0.5_DP
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)+0.5_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF

         IF (unique=='2') THEN
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.25_DP
         outco(2,2,i)=-inco(2,i)+0.75_DP
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)+0.75_DP
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.25_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.25_DP
         outco(3,4,i)=-inco(3,i)+0.75_DP
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.25_DP
         outco(3,6,i)=-inco(2,i)+0.75_DP
         !S=7
         outco(1,7,i)=-inco(3,i)+0.25_DP
         outco(2,7,i)=-inco(1,i)+0.75_DP
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)+0.75_DP
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.25_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)+0.75_DP
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.25_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.25_DP
         outco(3,11,i)=-inco(1,i)+0.75_DP
         !S=12
         outco(1,12,i)=-inco(2,i)+0.25_DP
         outco(2,12,i)=-inco(3,i)+0.75_DP
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)+0.5_DP
         outco(2,14,i)=-inco(1,i)+0.5_DP
         outco(3,14,i)=-inco(3,i)+0.5_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.5_DP
         outco(2,19,i)=-inco(3,i)+0.5_DP
         outco(3,19,i)=-inco(2,i)+0.5_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.5_DP
         outco(2,24,i)=-inco(2,i)+0.5_DP
         outco(3,24,i)=-inco(1,i)+0.5_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.75_DP
         outco(2,26,i)=inco(2,i)+0.25_DP
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)+0.25_DP
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.75_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.75_DP
         outco(3,28,i)=inco(3,i)+0.25_DP
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.75_DP
         outco(3,30,i)=inco(2,i)+0.25_DP
         !S=31
         outco(1,31,i)=inco(3,i)+0.75_DP
         outco(2,31,i)=inco(1,i)+0.25_DP
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)+0.25_DP
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.75_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)+0.25_DP
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.75_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.75_DP
         outco(3,35,i)=inco(1,i)+0.25_DP
         !S=36
         outco(1,36,i)=inco(2,i)+0.75_DP
         outco(2,36,i)=inco(3,i)+0.25_DP
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.25_DP
         outco(2,37,i)=-inco(1,i)+0.75_DP
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)+0.5_DP
         outco(2,38,i)=inco(1,i)+0.5_DP
         outco(3,38,i)=inco(3,i)+0.5_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.75_DP
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)+0.25_DP
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)+0.25_DP
         outco(3,40,i)=-inco(3,i)+0.75_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.25_DP
         outco(2,41,i)=-inco(3,i)+0.75_DP
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)+0.25_DP
         outco(3,42,i)=-inco(2,i)+0.75_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.5_DP
         outco(2,43,i)=inco(3,i)+0.5_DP
         outco(3,43,i)=inco(2,i)+0.5_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.75_DP
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)+0.25_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.25_DP
         outco(2,45,i)=-inco(2,i)+0.75_DP
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)+0.75_DP
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)+0.25_DP
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)+0.25_DP
         outco(3,47,i)=-inco(1,i)+0.75_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.5_DP
         outco(2,48,i)=inco(2,i)+0.5_DP
         outco(3,48,i)=inco(1,i)+0.5_DP
         END IF

      CASE (229)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)
         outco(3,3,i)=-inco(3,i)
         !S=4
         outco(1,4,i)=inco(1,i)
         outco(2,4,i)=-inco(2,i)
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)
         outco(2,6,i)=-inco(1,i)
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)
         outco(3,8,i)=-inco(2,i)
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)
         outco(3,10,i)=-inco(1,i)
         !S=11
         outco(1,11,i)=inco(2,i)
         outco(2,11,i)=-inco(3,i)
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)
         !S=13
         outco(1,13,i)=inco(2,i)
         outco(2,13,i)=inco(1,i)
         outco(3,13,i)=-inco(3,i)
         !S=14
         outco(1,14,i)=-inco(2,i)
         outco(2,14,i)=-inco(1,i)
         outco(3,14,i)=-inco(3,i)
         !S=15
         outco(1,15,i)=inco(2,i)
         outco(2,15,i)=-inco(1,i)
         outco(3,15,i)=inco(3,i)
         !S=16
         outco(1,16,i)=-inco(2,i)
         outco(2,16,i)=+inco(1,i)
         outco(3,16,i)=inco(3,i)
         !S=17
         outco(1,17,i)=+inco(1,i)
         outco(2,17,i)=+inco(3,i)
         outco(3,17,i)=-inco(2,i)
         !S=18
         outco(1,18,i)=-inco(1,i)
         outco(2,18,i)=+inco(3,i)
         outco(3,18,i)=inco(2,i)
         !S=19
         outco(1,19,i)=-inco(1,i)
         outco(2,19,i)=-inco(3,i)
         outco(3,19,i)=-inco(2,i)
         !S=20
         outco(1,20,i)=inco(1,i)
         outco(2,20,i)=-inco(3,i)
         outco(3,20,i)=inco(2,i)
         !S=21
         outco(1,21,i)=inco(3,i)
         outco(2,21,i)=inco(2,i)
         outco(3,21,i)=-inco(1,i)
         !S=22
         outco(1,22,i)=inco(3,i)
         outco(2,22,i)=-inco(2,i)
         outco(3,22,i)=inco(1,i)
         !S=23
         outco(1,23,i)=-inco(3,i)
         outco(2,23,i)=inco(2,i)
         outco(3,23,i)=inco(1,i)
         !S=24
         outco(1,24,i)=-inco(3,i)
         outco(2,24,i)=-inco(2,i)
         outco(3,24,i)=-inco(1,i)
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)
         outco(3,27,i)=inco(3,i)
         !S=28
         outco(1,28,i)=-inco(1,i)
         outco(2,28,i)=inco(2,i)
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)
         outco(2,30,i)=inco(1,i)
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)
         outco(3,32,i)=inco(2,i)
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)
         outco(3,34,i)=inco(1,i)
         !S=35
         outco(1,35,i)=-inco(2,i)
         outco(2,35,i)=inco(3,i)
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)
         !S=37
         outco(1,37,i)=-inco(2,i)
         outco(2,37,i)=-inco(1,i)
         outco(3,37,i)=inco(3,i)
         !S=38
         outco(1,38,i)=inco(2,i)
         outco(2,38,i)=inco(1,i)
         outco(3,38,i)=inco(3,i)
         !S=39
         outco(1,39,i)=-inco(2,i)
         outco(2,39,i)=inco(1,i)
         outco(3,39,i)=-inco(3,i)
         !S=40
         outco(1,40,i)=inco(2,i)
         outco(2,40,i)=-inco(1,i)
         outco(3,40,i)=-inco(3,i)
         !S=41
         outco(1,41,i)=-inco(1,i)
         outco(2,41,i)=-inco(3,i)
         outco(3,41,i)=+inco(2,i)
         !S=42
         outco(1,42,i)=inco(1,i)
         outco(2,42,i)=-inco(3,i)
         outco(3,42,i)=-inco(2,i)
         !S=43
         outco(1,43,i)=inco(1,i)
         outco(2,43,i)=inco(3,i)
         outco(3,43,i)=inco(2,i)
         !S=44
         outco(1,44,i)=-inco(1,i)
         outco(2,44,i)=+inco(3,i)
         outco(3,44,i)=-inco(2,i)
         !S=45
         outco(1,45,i)=-inco(3,i)
         outco(2,45,i)=-inco(2,i)
         outco(3,45,i)=+inco(1,i)
         !S=46
         outco(1,46,i)=-inco(3,i)
         outco(2,46,i)=inco(2,i)
         outco(3,46,i)=-inco(1,i)
         !S=47
         outco(1,47,i)=inco(3,i)
         outco(2,47,i)=-inco(2,i)
         outco(3,47,i)=-inco(1,i)
         !S=48
         outco(1,48,i)=inco(3,i)
         outco(2,48,i)=inco(2,i)
         outco(3,48,i)=inco(1,i)

      CASE (230)
         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
         !S=2
         outco(1,2,i)=-inco(1,i)+0.5_DP
         outco(2,2,i)=-inco(2,i)
         outco(3,2,i)=inco(3,i)+0.5_DP
         !S=3
         outco(1,3,i)=-inco(1,i)
         outco(2,3,i)=inco(2,i)+0.5_DP
         outco(3,3,i)=-inco(3,i)+0.5_DP
         !S=4
         outco(1,4,i)=inco(1,i)+0.5_DP
         outco(2,4,i)=-inco(2,i)+0.5_DP
         outco(3,4,i)=-inco(3,i)
         !S=5
         outco(1,5,i)=inco(3,i)
         outco(2,5,i)=inco(1,i)
         outco(3,5,i)=inco(2,i)
         !S=6
         outco(1,6,i)=inco(3,i)+0.5_DP
         outco(2,6,i)=-inco(1,i)+0.5_DP
         outco(3,6,i)=-inco(2,i)
         !S=7
         outco(1,7,i)=-inco(3,i)+0.5_DP
         outco(2,7,i)=-inco(1,i)
         outco(3,7,i)=inco(2,i)+0.5_DP
         !S=8
         outco(1,8,i)=-inco(3,i)
         outco(2,8,i)=inco(1,i)+0.5_DP
         outco(3,8,i)=-inco(2,i)+0.5_DP
         !S=9
         outco(1,9,i)=inco(2,i)
         outco(2,9,i)=inco(3,i)
         outco(3,9,i)=inco(1,i)
         !S=10
         outco(1,10,i)=-inco(2,i)
         outco(2,10,i)=inco(3,i)+0.5_DP
         outco(3,10,i)=-inco(1,i)+0.5_DP
         !S=11
         outco(1,11,i)=inco(2,i)+0.5_DP
         outco(2,11,i)=-inco(3,i)+0.5_DP
         outco(3,11,i)=-inco(1,i)
         !S=12
         outco(1,12,i)=-inco(2,i)+0.5_DP
         outco(2,12,i)=-inco(3,i)
         outco(3,12,i)=inco(1,i)+0.5_DP
         !S=13
         outco(1,13,i)=inco(2,i)+0.75_DP
         outco(2,13,i)=inco(1,i)+0.25_DP
         outco(3,13,i)=-inco(3,i)+0.25_DP
         !S=14
         outco(1,14,i)=-inco(2,i)+0.75_DP
         outco(2,14,i)=-inco(1,i)+0.75_DP
         outco(3,14,i)=-inco(3,i)+0.75_DP
         !S=15
         outco(1,15,i)=inco(2,i)+0.25_DP
         outco(2,15,i)=-inco(1,i)+0.25_DP
         outco(3,15,i)=inco(3,i)+0.75_DP
         !S=16
         outco(1,16,i)=-inco(2,i)+0.25_DP
         outco(2,16,i)=+inco(1,i)+0.75_DP
         outco(3,16,i)=inco(3,i)+0.25_DP
         !S=17
         outco(1,17,i)=+inco(1,i)+0.75_DP
         outco(2,17,i)=+inco(3,i)+0.25_DP
         outco(3,17,i)=-inco(2,i)+0.25_DP
         !S=18
         outco(1,18,i)=-inco(1,i)+0.25_DP
         outco(2,18,i)=+inco(3,i)+0.75_DP
         outco(3,18,i)=inco(2,i)+0.25_DP
         !S=19
         outco(1,19,i)=-inco(1,i)+0.75_DP
         outco(2,19,i)=-inco(3,i)+0.75_DP
         outco(3,19,i)=-inco(2,i)+0.75_DP
         !S=20
         outco(1,20,i)=inco(1,i)+0.25_DP
         outco(2,20,i)=-inco(3,i)+0.25_DP
         outco(3,20,i)=inco(2,i)+0.75_DP
         !S=21
         outco(1,21,i)=inco(3,i)+0.75_DP
         outco(2,21,i)=inco(2,i)+0.25_DP
         outco(3,21,i)=-inco(1,i)+0.25_DP
         !S=22
         outco(1,22,i)=inco(3,i)+0.25_DP
         outco(2,22,i)=-inco(2,i)+0.25_DP
         outco(3,22,i)=inco(1,i)+0.75_DP
         !S=23
         outco(1,23,i)=-inco(3,i)+0.25_DP
         outco(2,23,i)=inco(2,i)+0.75_DP
         outco(3,23,i)=inco(1,i)+0.25_DP
         !S=24
         outco(1,24,i)=-inco(3,i)+0.75_DP
         outco(2,24,i)=-inco(2,i)+0.75_DP
         outco(3,24,i)=-inco(1,i)+0.75_DP
         !S=25
         outco(1,25,i)=-inco(1,i)
         outco(2,25,i)=-inco(2,i)
         outco(3,25,i)=-inco(3,i)
         !S=26
         outco(1,26,i)=inco(1,i)+0.5_DP
         outco(2,26,i)=inco(2,i)
         outco(3,26,i)=-inco(3,i)+0.5_DP
         !S=27
         outco(1,27,i)=inco(1,i)
         outco(2,27,i)=-inco(2,i)+0.5_DP
         outco(3,27,i)=inco(3,i)+0.5_DP
         !S=28
         outco(1,28,i)=-inco(1,i)+0.5_DP
         outco(2,28,i)=inco(2,i)+0.5_DP
         outco(3,28,i)=inco(3,i)
         !S=29
         outco(1,29,i)=-inco(3,i)
         outco(2,29,i)=-inco(1,i)
         outco(3,29,i)=-inco(2,i)
         !S=30
         outco(1,30,i)=-inco(3,i)+0.5_DP
         outco(2,30,i)=inco(1,i)+0.5_DP
         outco(3,30,i)=inco(2,i)
         !S=31
         outco(1,31,i)=inco(3,i)+0.5_DP
         outco(2,31,i)=inco(1,i)
         outco(3,31,i)=-inco(2,i)+0.5_DP
         !S=32
         outco(1,32,i)=inco(3,i)
         outco(2,32,i)=-inco(1,i)+0.5_DP
         outco(3,32,i)=inco(2,i)+0.5_DP
         !S=33
         outco(1,33,i)=-inco(2,i)
         outco(2,33,i)=-inco(3,i)
         outco(3,33,i)=-inco(1,i)
         !S=34
         outco(1,34,i)=inco(2,i)
         outco(2,34,i)=-inco(3,i)+0.5_DP
         outco(3,34,i)=inco(1,i)+0.5_DP
         !S=35
         outco(1,35,i)=-inco(2,i)+0.5_DP
         outco(2,35,i)=inco(3,i)+0.5_DP
         outco(3,35,i)=inco(1,i)
         !S=36
         outco(1,36,i)=inco(2,i)+0.5_DP
         outco(2,36,i)=inco(3,i)
         outco(3,36,i)=-inco(1,i)+0.5_DP
         !S=37
         outco(1,37,i)=-inco(2,i)+0.25_DP
         outco(2,37,i)=-inco(1,i)+0.75_DP
         outco(3,37,i)=inco(3,i)+0.75_DP
         !S=38
         outco(1,38,i)=inco(2,i)+0.25_DP
         outco(2,38,i)=inco(1,i)+0.25_DP
         outco(3,38,i)=inco(3,i)+0.25_DP
         !S=39
         outco(1,39,i)=-inco(2,i)+0.75_DP
         outco(2,39,i)=inco(1,i)+0.75_DP
         outco(3,39,i)=-inco(3,i)+0.25_DP
         !S=40
         outco(1,40,i)=inco(2,i)+0.75_DP
         outco(2,40,i)=-inco(1,i)+0.25_DP
         outco(3,40,i)=-inco(3,i)+0.75_DP
         !S=41
         outco(1,41,i)=-inco(1,i)+0.25_DP
         outco(2,41,i)=-inco(3,i)+0.75_DP
         outco(3,41,i)=+inco(2,i)+0.75_DP
         !S=42
         outco(1,42,i)=inco(1,i)+0.75_DP
         outco(2,42,i)=-inco(3,i)+0.25_DP
         outco(3,42,i)=-inco(2,i)+0.75_DP
         !S=43
         outco(1,43,i)=inco(1,i)+0.25_DP
         outco(2,43,i)=inco(3,i)+0.25_DP
         outco(3,43,i)=inco(2,i)+0.25_DP
         !S=44
         outco(1,44,i)=-inco(1,i)+0.75_DP
         outco(2,44,i)=+inco(3,i)+0.75_DP
         outco(3,44,i)=-inco(2,i)+0.25_DP
         !S=45
         outco(1,45,i)=-inco(3,i)+0.25_DP
         outco(2,45,i)=-inco(2,i)+0.75_DP
         outco(3,45,i)=+inco(1,i)+0.75_DP
         !S=46
         outco(1,46,i)=-inco(3,i)+0.75_DP
         outco(2,46,i)=inco(2,i)+0.75_DP
         outco(3,46,i)=-inco(1,i)+0.25_DP
         !S=47
         outco(1,47,i)=inco(3,i)+0.75_DP
         outco(2,47,i)=-inco(2,i)+0.25_DP
         outco(3,47,i)=-inco(1,i)+0.75_DP
         !S=48
         outco(1,48,i)=inco(3,i)+0.25_DP
         outco(2,48,i)=inco(2,i)+0.25_DP
         outco(3,48,i)=inco(1,i)+0.25_DP
      END SELECT simmetria   
      RETURN
   END SUBROUTINE find_equivalent_tau
END MODULE
