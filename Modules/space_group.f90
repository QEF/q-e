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

  REAL(DP), PARAMETER :: unterz=(1.0_DP)/(3.0_DP)
  REAL(DP), PARAMETER :: duterz=(2.0_DP)/(3.0_DP)
  REAL(DP), PARAMETER :: unsest=(1.0_DP)/(6.0_DP)
  REAL(DP), PARAMETER :: cisest=(5.0_DP)/(6.0_DP)

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

      INTEGER :: k,j
      simmetria: SELECT CASE (space_group_number)
      !*****************************************
      !Triclinic 1-2
      CASE (1)
             CALL find_equiv_1  ( i, inco, outco )
      CASE (2)
             CALL find_equiv_2  ( i, inco, outco )
      CASE (3)
             CALL find_equiv_3  ( i, inco, unique, outco )
      CASE (4)
             CALL find_equiv_4  ( i, inco, unique, outco )
      CASE (5)
             CALL find_equiv_5  ( i, inco, unique, outco )
      CASE (6)
             CALL find_equiv_6  ( i, inco, unique, outco )
      CASE (7)
             CALL find_equiv_7  ( i, inco, unique, outco )
      CASE (8)
             CALL find_equiv_8  ( i, inco, unique, outco )
      CASE (9)
             CALL find_equiv_9  ( i, inco, unique, outco )
      CASE (10)
             CALL find_equiv_10 ( i, inco, unique, outco )
      CASE (11)
             CALL find_equiv_11 ( i, inco, unique, outco )
      CASE (12)
             CALL find_equiv_12 ( i, inco, unique, outco )
      CASE (13)
             CALL find_equiv_13 ( i, inco, unique, outco )
      CASE (14)
             CALL find_equiv_14 ( i, inco, unique, outco )
      CASE (15)
             CALL find_equiv_15 ( i, inco, unique, outco )
      CASE (16) !P222
             CALL find_equiv_16 ( i, inco, outco )
      CASE (17) !P222(1)
             CALL find_equiv_17 ( i, inco, outco )
      CASE (18) !P2(1)2(1)2
             CALL find_equiv_18 ( i, inco, outco )
      CASE (19) !P2(1)2(1)2(1)
             CALL find_equiv_19 ( i, inco, outco )
      CASE (20) !C222(1)
             CALL find_equiv_20 ( i, inco, outco )
      CASE (21) !C222
             CALL find_equiv_21 ( i, inco, outco )
      CASE (22) !F222
             CALL find_equiv_22 ( i, inco, outco )
      CASE (23) !I222
             CALL find_equiv_23 ( i, inco, outco )
      CASE (24) !I2(1)2(1)2(1)
             CALL find_equiv_24 ( i, inco, outco )
      CASE (25) !Pmm2
             CALL find_equiv_25 ( i, inco, outco )
      CASE (26) !Pmc2(1)
             CALL find_equiv_26 ( i, inco, outco )
      CASE (27) !Pcc2
             CALL find_equiv_27 ( i, inco, outco )
      CASE (28) !Pma2
             CALL find_equiv_28 ( i, inco, outco )
      CASE (29) !Pca2(1)
             CALL find_equiv_29 ( i, inco, outco )
      CASE (30) !Pnc2
             CALL find_equiv_30 ( i, inco, outco )
      CASE (31) !Pmn2(1)
             CALL find_equiv_31 ( i, inco, outco )
      CASE (32) !Pba2
             CALL find_equiv_32 ( i, inco, outco )
      CASE (33) !Pna2(1)
             CALL find_equiv_33 ( i, inco, outco )
      CASE (34) !Pnn2
             CALL find_equiv_34 ( i, inco, outco )
      CASE (35) !Cmm2
             CALL find_equiv_35 ( i, inco, outco )
      CASE (36) !Cmc2(1)
             CALL find_equiv_36 ( i, inco, outco )
      CASE (37) !Ccc2
             CALL find_equiv_37 ( i, inco, outco )
      CASE (38) !Amm2
             CALL find_equiv_38 ( i, inco, outco )
      CASE (39) !Abm2
             CALL find_equiv_39 ( i, inco, outco )
      CASE (40) !Ama2
             CALL find_equiv_40 ( i, inco, outco )
      CASE (41) !Aba2
             CALL find_equiv_41 ( i, inco, outco )
      CASE (42) !Fmm2
             CALL find_equiv_42 ( i, inco, outco )
      CASE (43) !Fdd2
             CALL find_equiv_43 ( i, inco, outco )
      CASE (44) !Imm2
             CALL find_equiv_44 ( i, inco, outco )
      CASE (45) !Iba2
             CALL find_equiv_45 ( i, inco, outco )
      CASE (46) !Ima2
             CALL find_equiv_46 ( i, inco, outco )
      CASE (47) !Pmmm
             CALL find_equiv_47 ( i, inco, outco )
      CASE (48) !Pnnn
             CALL find_equiv_48 ( i, inco, unique, outco )
      CASE (49) !Pccm
             CALL find_equiv_49 ( i, inco, outco )
      CASE (50) !Pban
             CALL find_equiv_50 ( i, inco, unique, outco )
      CASE (51) !Pmma
             CALL find_equiv_51 ( i, inco, outco )
      CASE (52) !Pnna
             CALL find_equiv_52 ( i, inco, outco )
      CASE (53) !Pmna
             CALL find_equiv_53 ( i, inco, outco )
      CASE (54) !Pcca
             CALL find_equiv_54 ( i, inco, outco )
      CASE (55) !Pbam
             CALL find_equiv_55 ( i, inco, outco )
      CASE (56) !Pccn
             CALL find_equiv_56 ( i, inco, outco )
      CASE (57) !Pbcm
             CALL find_equiv_57 ( i, inco, outco )
      CASE (58) !Pnnm
             CALL find_equiv_58 ( i, inco, outco )
      CASE (59) !Pmmn
             CALL find_equiv_59 ( i, inco, unique, outco )
      CASE (60) !Pbcn
             CALL find_equiv_60 ( i, inco, outco )
      CASE (61) !Pbca
             CALL find_equiv_61 ( i, inco, outco )
      CASE (62) !Pnma
             CALL find_equiv_62 ( i, inco, outco )
      CASE (63) !Cmcm
             CALL find_equiv_63 ( i, inco, outco )
      CASE (64) !Cmca
             CALL find_equiv_64 ( i, inco, outco )
      CASE (65) !Cmmm
             CALL find_equiv_65 ( i, inco, outco )
      CASE (66) !Cccm
             CALL find_equiv_66 ( i, inco, outco )
      CASE (67) !Cmma
             CALL find_equiv_67 ( i, inco, outco )
      CASE (68) !Ccca
             CALL find_equiv_68 ( i, inco, unique, outco )
      CASE (69) !Fmmm
             CALL find_equiv_69 ( i, inco, outco )
      CASE (70) !Fddd
             CALL find_equiv_70 ( i, inco, unique, outco )
      CASE (71) !Immm
             CALL find_equiv_71 ( i, inco, outco )
      CASE (72) !Ibam
             CALL find_equiv_72 ( i, inco, outco )
      CASE (73) !Ibca
             CALL find_equiv_73 ( i, inco, outco )
      CASE (74) !Imma
             CALL find_equiv_74 ( i, inco, outco )
      CASE (75) !P4
             CALL find_equiv_75 ( i, inco, outco )
      CASE (76) !P4(1)
             CALL find_equiv_76 ( i, inco, outco )
      CASE (77) !P4(2)
             CALL find_equiv_77 ( i, inco, outco )
      CASE (78) !P4(3)
             CALL find_equiv_78 ( i, inco, outco )
      CASE (79) !I4
             CALL find_equiv_79 ( i, inco, outco )
      CASE (80) !I4(1)
             CALL find_equiv_80 ( i, inco, outco )
      CASE (81) !P-4
             CALL find_equiv_81 ( i, inco, outco )
      CASE (82) !I-4
             CALL find_equiv_82 ( i, inco, outco )
      CASE (83) !P4/m
             CALL find_equiv_83 ( i, inco, outco )
      CASE (84) !P(2)/m
             CALL find_equiv_84 ( i, inco, outco )
      CASE (85) !P4/n
             CALL find_equiv_85 ( i, inco, unique, outco )
      CASE (86) !P4(2)/n
             CALL find_equiv_86 ( i, inco, unique, outco )
      CASE (87) !I4/m
             CALL find_equiv_87 ( i, inco, outco )
      CASE (88) !I4(1)/a
             CALL find_equiv_88 ( i, inco, unique, outco )
      CASE (89) !P422
             CALL find_equiv_89 ( i, inco, outco )
      CASE (90) !P42(1)2
             CALL find_equiv_90 ( i, inco, outco )
      CASE (91) !P4(1)22
             CALL find_equiv_91 ( i, inco, outco )
      CASE (92) !P4(1)2(1)2
             CALL find_equiv_92 ( i, inco, outco )
      CASE (93) !P4(2)22
             CALL find_equiv_93 ( i, inco, outco )
      CASE (94) !P4(2)2(1)2
             CALL find_equiv_94 ( i, inco, outco )
      CASE (95) !P4(3)22
             CALL find_equiv_95 ( i, inco, outco )
      CASE (96) !P4(3)2(1)2
             CALL find_equiv_96 ( i, inco, outco )
      CASE (97) !I422
             CALL find_equiv_97 ( i, inco, outco )
      CASE (98) !I4(1)22
             CALL find_equiv_98 ( i, inco, outco )
      CASE (99) !P4mm
             CALL find_equiv_99 ( i, inco, outco )
      CASE (100) !P4bm
             CALL find_equiv_100( i, inco, outco )
      CASE (101) !P4(2)cm
             CALL find_equiv_101( i, inco, outco )
      CASE (102) !P4(2)nm
             CALL find_equiv_102( i, inco, outco )
      CASE (103) !P4cc
             CALL find_equiv_103( i, inco, outco )
      CASE (104) !P4nc
             CALL find_equiv_104( i, inco, outco )
      CASE (105) !P4(2)mc
             CALL find_equiv_105( i, inco, outco )
      CASE (106) !P4(2)bc
             CALL find_equiv_106( i, inco, outco )
      CASE (107) !I4mm
             CALL find_equiv_107( i, inco, outco )
      CASE (108) !I4cm
             CALL find_equiv_108( i, inco, outco )
      CASE (109) !I4(1)md
             CALL find_equiv_109( i, inco, outco )
      CASE (110) !I4(1)cd
             CALL find_equiv_110( i, inco, outco )
      CASE (111) !P-42m
             CALL find_equiv_111( i, inco, outco )
      CASE (112) !P-42c
             CALL find_equiv_112( i, inco, outco )
      CASE (113) !P-42(1)m
             CALL find_equiv_113( i, inco, outco )
      CASE (114) !P-42(1)c
             CALL find_equiv_114( i, inco, outco )
      CASE (115) !P-4m2
             CALL find_equiv_115( i, inco, outco )
      CASE (116) !P-4c2
             CALL find_equiv_116( i, inco, outco )
      CASE (117) !P-4b2
             CALL find_equiv_117( i, inco, outco )
      CASE (118) !P-4n2
             CALL find_equiv_118( i, inco, outco )
      CASE (119) !I-4m2
             CALL find_equiv_119( i, inco, outco )
      CASE (120) !I-4c2
             CALL find_equiv_120( i, inco, outco )
      CASE (121) !I-42m
             CALL find_equiv_121( i, inco, outco )
      CASE (122) !I-42d
             CALL find_equiv_122( i, inco, outco )
      CASE (123) !P4/mmm
             CALL find_equiv_123( i, inco, outco )
      CASE (124) !P4/mcc
             CALL find_equiv_124( i, inco, outco )
      CASE (125) !P4/nbm
             CALL find_equiv_125( i, inco, unique, outco )
      CASE (126) !P4/nnc
             CALL find_equiv_126( i, inco, unique, outco )
      CASE (127) !P4/mbm
             CALL find_equiv_127( i, inco, outco )
      CASE (128) !P4/mnc
             CALL find_equiv_128( i, inco, outco )
      CASE (129)
             CALL find_equiv_129( i, inco, unique, outco )
      CASE (130) !P4/ncc
             CALL find_equiv_130( i, inco, unique, outco )
      CASE (131) !P4(2)/mmc
             CALL find_equiv_131( i, inco, outco )
      CASE (132) !P4(2)mcm
             CALL find_equiv_132( i, inco, outco )
      CASE (133) !P4(2)/nbc
             CALL find_equiv_133( i, inco, unique, outco )
      CASE (134) !P4(2)/nnm
             CALL find_equiv_134( i, inco, unique, outco )
      CASE (135) !P4(2)/mbc
             CALL find_equiv_135( i, inco, outco )
      CASE (136) !P4(2)mnm
             CALL find_equiv_136( i, inco, outco )
      CASE (137) !P4(2)/nmc
             CALL find_equiv_137( i, inco, unique, outco )
      CASE (138) !P4(2)/ncm
             CALL find_equiv_138( i, inco, unique, outco )
      CASE (139) !I4/mmm
             CALL find_equiv_139( i, inco, outco )
      CASE (140) !I4/mcm
             CALL find_equiv_140( i, inco, outco )
      CASE (141) !I4(1)amd
             CALL find_equiv_141( i, inco, unique, outco )
      CASE (142) !I4(1)/acd
             CALL find_equiv_142( i, inco, unique, outco )
      CASE (143) !P3
             CALL find_equiv_143( i, inco, outco )
      CASE (144) !P3(1)
             CALL find_equiv_144( i, inco, outco )
      CASE (145) !P3(2)
             CALL find_equiv_145( i, inco, outco )
      CASE (146) !R3
             CALL find_equiv_146( i, inco, unique, outco )
      CASE (147) !P-3
             CALL find_equiv_147( i, inco, outco )
      CASE (148) !R-3
             CALL find_equiv_148( i, inco, unique, outco )
      CASE (149) !P312
             CALL find_equiv_149( i, inco, outco )
      CASE (150) !P321
             CALL find_equiv_150( i, inco, outco )
      CASE (151) !P3(1)12
             CALL find_equiv_151( i, inco, outco )
      CASE (152) !P3(1)21
             CALL find_equiv_152( i, inco, outco )
      CASE (153) !P3(2)12
             CALL find_equiv_153( i, inco, outco )
      CASE (154) !P3(2)21
             CALL find_equiv_154( i, inco, outco )
      CASE (155) !R32
             CALL find_equiv_155( i, inco, unique, outco )
      CASE (156) !P3m1
             CALL find_equiv_156( i, inco, outco )
      CASE (157) !P31m
             CALL find_equiv_157( i, inco, outco )
      CASE (158) !P3c1
             CALL find_equiv_158( i, inco, outco )
      CASE (159) !P31c
             CALL find_equiv_159( i, inco, outco )
      CASE (160) !R3m
             CALL find_equiv_160( i, inco, unique, outco )
      CASE (161) !R3c
             CALL find_equiv_161( i, inco, unique, outco )
      CASE (162) !P-31m
             CALL find_equiv_162( i, inco, outco )
      CASE (163) !P-31c
             CALL find_equiv_163( i, inco, outco )
      CASE (164) !P-3m1
             CALL find_equiv_164( i, inco, outco )
      CASE (165) !P-3c1
             CALL find_equiv_165( i, inco, outco )
      CASE (166) !R-3m
             CALL find_equiv_166( i, inco, unique, outco )
      CASE (167) !R-3c
             CALL find_equiv_167( i, inco, unique, outco )
      CASE (168) !P6
             CALL find_equiv_168( i, inco, outco )
      CASE (169) !P6(1)
             CALL find_equiv_169( i, inco, outco )
      CASE (170) !P6(5)
             CALL find_equiv_170( i, inco, outco )
      CASE (171) !P6(2)
             CALL find_equiv_171( i, inco, outco )
      CASE (172) !P6(4)
             CALL find_equiv_172( i, inco, outco )
      CASE (173) !P6(3)
             CALL find_equiv_173( i, inco, outco )
      CASE (174) !P-6
             CALL find_equiv_174( i, inco, outco )
      CASE (175) !P6/m
             CALL find_equiv_175( i, inco, outco )
      CASE (176) !P6(3)/m
             CALL find_equiv_176( i, inco, outco )
      CASE (177) !P622
             CALL find_equiv_177( i, inco, outco )
      CASE (178) !P(1)22
             CALL find_equiv_178( i, inco, outco )
      CASE (179) !P6(5)22
             CALL find_equiv_179( i, inco, outco )
      CASE (180) !P6(2)22
             CALL find_equiv_180( i, inco, outco )
      CASE (181) !P6(4)22
             CALL find_equiv_181( i, inco, outco )
      CASE (182) !6(3)22
             CALL find_equiv_182( i, inco, outco )
      CASE (183) !P6mm
             CALL find_equiv_183( i, inco, outco )
      CASE (184) !P6cc
             CALL find_equiv_184( i, inco, outco )
      CASE (185) !P6(3)cm
             CALL find_equiv_185( i, inco, outco )
      CASE (186) !P(3)mc
             CALL find_equiv_186( i, inco, outco )
      CASE (187) !P-6m2
             CALL find_equiv_187( i, inco, outco )
      CASE (188) !P-6c2
             CALL find_equiv_188( i, inco, outco )
      CASE (189) !P-62m
             CALL find_equiv_189( i, inco, outco )
      CASE (190) !P-62c
             CALL find_equiv_190( i, inco, outco )
      CASE (191) !P6/mmm
             CALL find_equiv_191( i, inco, outco )
      CASE (192) !P6/mmc
             CALL find_equiv_192( i, inco, outco )
      CASE (193) !P6(3)/mcm
             CALL find_equiv_193( i, inco, outco )
      CASE (194)
             CALL find_equiv_194( i, inco, outco )
      CASE (195) !P23
             CALL find_equiv_195( i, inco, outco )
      CASE (196) !F23
             CALL find_equiv_196( i, inco, outco )
      CASE (197) !I23
             CALL find_equiv_197( i, inco, outco )
      CASE (198) !P2(1)3
             CALL find_equiv_198( i, inco, outco )
      CASE (199) !I2(1)3
             CALL find_equiv_199( i, inco, outco )
      CASE (200) !Pm-3
             CALL find_equiv_200( i, inco, outco )
      CASE(201) !Pn-3
             CALL find_equiv_201( i, inco, unique, outco )
      CASE (202) !Fm-3
             CALL find_equiv_202( i, inco, outco )
      CASE (203) !Fd-3
             CALL find_equiv_203( i, inco, unique, outco )
      CASE (204) !Im-3
             CALL find_equiv_204( i, inco, outco )
      CASE (205) !Pa-3
             CALL find_equiv_205( i, inco, outco )
      CASE (206) !Ia-3
             CALL find_equiv_206( i, inco, outco )
      CASE (207) !P432
             CALL find_equiv_207( i, inco, outco )
      CASE (208) !P4(2)32
             CALL find_equiv_208( i, inco, outco )
      CASE (209) !F432
             CALL find_equiv_209( i, inco, outco )
      CASE (210) !F4(1)32
             CALL find_equiv_210( i, inco, outco )
      CASE (211) !I432
             CALL find_equiv_211( i, inco, outco )
      CASE (212) !P4(3)32
             CALL find_equiv_212( i, inco, outco )
      CASE (213) !P4(1)32
             CALL find_equiv_213( i, inco, outco )
      CASE (214) !I4(1)32
             CALL find_equiv_214( i, inco, outco )
      CASE (215) !P-43m
             CALL find_equiv_215( i, inco, outco )
      CASE (216) !F-43m
             CALL find_equiv_216( i, inco, outco )
      CASE (217) !I-43m
             CALL find_equiv_217( i, inco, outco )
      CASE (218) !P-43n
             CALL find_equiv_218( i, inco, outco )
      CASE (219) !F-43c
             CALL find_equiv_219( i, inco, outco )
      CASE (220) !I-43d
             CALL find_equiv_220( i, inco, outco )
      CASE (221) !Pm-3m
             CALL find_equiv_221( i, inco, outco )
      CASE (222) !Pn-3n
             CALL find_equiv_222( i, inco, unique, outco )
      CASE (223) !Pm-3n
             CALL find_equiv_223( i, inco, outco )
      CASE (224) !Pn-3m
             CALL find_equiv_224( i, inco, unique, outco )
      CASE (225) !Fm-3m
             CALL find_equiv_225( i, inco, outco )
      CASE (226) !Fm-3c
             CALL find_equiv_226( i, inco, outco )
      CASE (227) !Fd-3m
             CALL find_equiv_227( i, inco, unique, outco )
      CASE (228) !Fd-3c
             CALL find_equiv_228( i, inco, unique, outco )
      CASE (229)
             CALL find_equiv_229( i, inco, outco )
      CASE (230)
             CALL find_equiv_230( i, inco, outco )
      END SELECT simmetria

    END SUBROUTINE find_equivalent_tau

SUBROUTINE find_equiv_1  ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         END DO
END SUBROUTINE find_equiv_1  

SUBROUTINE find_equiv_2  ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

         DO k=1,3
         outco(k,1,i)=inco(k,i)
         outco(k,2,i)=-inco(k,i)
         END DO
      !*****************************************
      !Monoclinic 3-15
END SUBROUTINE find_equiv_2  

SUBROUTINE find_equiv_3  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_3  

SUBROUTINE find_equiv_4  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_4  

SUBROUTINE find_equiv_5  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_5  

SUBROUTINE find_equiv_6  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_6  

SUBROUTINE find_equiv_7  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_7  

SUBROUTINE find_equiv_8  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_8  

SUBROUTINE find_equiv_9  ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_9  

SUBROUTINE find_equiv_10 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_10 

SUBROUTINE find_equiv_11 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_11 

SUBROUTINE find_equiv_12 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_12 

SUBROUTINE find_equiv_13 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_13 

SUBROUTINE find_equiv_14 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_14 

SUBROUTINE find_equiv_15 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_15 

SUBROUTINE find_equiv_16 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_16 

SUBROUTINE find_equiv_17 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_17 

SUBROUTINE find_equiv_18 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_18 

SUBROUTINE find_equiv_19 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_19 

SUBROUTINE find_equiv_20 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_20 

SUBROUTINE find_equiv_21 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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


END SUBROUTINE find_equiv_21 

SUBROUTINE find_equiv_22 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_22 

SUBROUTINE find_equiv_23 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_23 

SUBROUTINE find_equiv_24 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_24 

SUBROUTINE find_equiv_25 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_25 

SUBROUTINE find_equiv_26 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_26 

SUBROUTINE find_equiv_27 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_27 

SUBROUTINE find_equiv_28 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_28 

SUBROUTINE find_equiv_29 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_29 

SUBROUTINE find_equiv_30 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_30 

SUBROUTINE find_equiv_31 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_31 

SUBROUTINE find_equiv_32 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_32 

SUBROUTINE find_equiv_33 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_33 

SUBROUTINE find_equiv_34 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_34 

SUBROUTINE find_equiv_35 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_35 

SUBROUTINE find_equiv_36 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_36 

SUBROUTINE find_equiv_37 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_37 

SUBROUTINE find_equiv_38 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_38 

SUBROUTINE find_equiv_39 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_39 

SUBROUTINE find_equiv_40 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_40 

SUBROUTINE find_equiv_41 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_41 

SUBROUTINE find_equiv_42 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_42 

SUBROUTINE find_equiv_43 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_43 

SUBROUTINE find_equiv_44 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_44 

SUBROUTINE find_equiv_45 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_45 

SUBROUTINE find_equiv_46 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_46 

SUBROUTINE find_equiv_47 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_47 

SUBROUTINE find_equiv_48 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_48 

SUBROUTINE find_equiv_49 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_49 

SUBROUTINE find_equiv_50 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_50 

SUBROUTINE find_equiv_51 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_51 

SUBROUTINE find_equiv_52 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_52 

SUBROUTINE find_equiv_53 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_53 

SUBROUTINE find_equiv_54 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_54 

SUBROUTINE find_equiv_55 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_55 

SUBROUTINE find_equiv_56 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_56 

SUBROUTINE find_equiv_57 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_57 

SUBROUTINE find_equiv_58 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_58 

SUBROUTINE find_equiv_59 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_59 

SUBROUTINE find_equiv_60 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_60 

SUBROUTINE find_equiv_61 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_61 

SUBROUTINE find_equiv_62 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_62 

SUBROUTINE find_equiv_63 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_63 

SUBROUTINE find_equiv_64 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_64 

SUBROUTINE find_equiv_65 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_65 

SUBROUTINE find_equiv_66 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_66 

SUBROUTINE find_equiv_67 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_67 

SUBROUTINE find_equiv_68 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_68 

SUBROUTINE find_equiv_69 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_69 

SUBROUTINE find_equiv_70 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_70 

SUBROUTINE find_equiv_71 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_71 

SUBROUTINE find_equiv_72 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_72 

SUBROUTINE find_equiv_73 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_73 

SUBROUTINE find_equiv_74 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_74 

SUBROUTINE find_equiv_75 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_75 

SUBROUTINE find_equiv_76 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_76 

SUBROUTINE find_equiv_77 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_77 

SUBROUTINE find_equiv_78 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_78 

SUBROUTINE find_equiv_79 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_79 

SUBROUTINE find_equiv_80 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_80 

SUBROUTINE find_equiv_81 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_81 

SUBROUTINE find_equiv_82 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_82 

SUBROUTINE find_equiv_83 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_83 

SUBROUTINE find_equiv_84 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_84 

SUBROUTINE find_equiv_85 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k


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

END SUBROUTINE find_equiv_85 

SUBROUTINE find_equiv_86 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_86 

SUBROUTINE find_equiv_87 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_87 

SUBROUTINE find_equiv_88 ( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_88 

SUBROUTINE find_equiv_89 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_89 

SUBROUTINE find_equiv_90 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_90 

SUBROUTINE find_equiv_91 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_91 

SUBROUTINE find_equiv_92 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_92 

SUBROUTINE find_equiv_93 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_93 

SUBROUTINE find_equiv_94 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_94 

SUBROUTINE find_equiv_95 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_95 

SUBROUTINE find_equiv_96 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_96 

SUBROUTINE find_equiv_97 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_97 

SUBROUTINE find_equiv_98 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_98 

SUBROUTINE find_equiv_99 ( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_99 

SUBROUTINE find_equiv_100( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_100

SUBROUTINE find_equiv_101( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_101

SUBROUTINE find_equiv_102( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_102

SUBROUTINE find_equiv_103( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_103

SUBROUTINE find_equiv_104( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_104

SUBROUTINE find_equiv_105( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_105

SUBROUTINE find_equiv_106( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_106

SUBROUTINE find_equiv_107( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_107

SUBROUTINE find_equiv_108( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_108

SUBROUTINE find_equiv_109( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_109

SUBROUTINE find_equiv_110( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_110

SUBROUTINE find_equiv_111( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_111

SUBROUTINE find_equiv_112( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_112

SUBROUTINE find_equiv_113( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_113

SUBROUTINE find_equiv_114( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_114

SUBROUTINE find_equiv_115( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_115

SUBROUTINE find_equiv_116( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_116

SUBROUTINE find_equiv_117( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_117

SUBROUTINE find_equiv_118( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_118

SUBROUTINE find_equiv_119( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_119

SUBROUTINE find_equiv_120( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_120

SUBROUTINE find_equiv_121( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_121

SUBROUTINE find_equiv_122( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_122

SUBROUTINE find_equiv_123( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_123

SUBROUTINE find_equiv_124( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_124

SUBROUTINE find_equiv_125( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_125

SUBROUTINE find_equiv_126( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_126

SUBROUTINE find_equiv_127( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_127

SUBROUTINE find_equiv_128( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_128

SUBROUTINE find_equiv_129( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_129

SUBROUTINE find_equiv_130( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_130

SUBROUTINE find_equiv_131( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_131

SUBROUTINE find_equiv_132( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_132

SUBROUTINE find_equiv_133( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_133

SUBROUTINE find_equiv_134( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_134

SUBROUTINE find_equiv_135( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_135

SUBROUTINE find_equiv_136( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_136

SUBROUTINE find_equiv_137( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_137

SUBROUTINE find_equiv_138( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_138

SUBROUTINE find_equiv_139( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_139

SUBROUTINE find_equiv_140( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_140

SUBROUTINE find_equiv_141( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_141

SUBROUTINE find_equiv_142( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_142

SUBROUTINE find_equiv_143( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_143

SUBROUTINE find_equiv_144( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_144

SUBROUTINE find_equiv_145( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_145

SUBROUTINE find_equiv_146( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_146

SUBROUTINE find_equiv_147( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_147

SUBROUTINE find_equiv_148( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_148

SUBROUTINE find_equiv_149( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_149

SUBROUTINE find_equiv_150( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_150

SUBROUTINE find_equiv_151( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_151

SUBROUTINE find_equiv_152( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_152

SUBROUTINE find_equiv_153( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_153

SUBROUTINE find_equiv_154( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_154

SUBROUTINE find_equiv_155( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_155

SUBROUTINE find_equiv_156( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_156

SUBROUTINE find_equiv_157( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_157

SUBROUTINE find_equiv_158( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_158

SUBROUTINE find_equiv_159( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_159

SUBROUTINE find_equiv_160( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_160

SUBROUTINE find_equiv_161( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_161

SUBROUTINE find_equiv_162( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_162

SUBROUTINE find_equiv_163( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_163

SUBROUTINE find_equiv_164( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_164

SUBROUTINE find_equiv_165( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_165

SUBROUTINE find_equiv_166( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_166

SUBROUTINE find_equiv_167( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_167

SUBROUTINE find_equiv_168( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_168

SUBROUTINE find_equiv_169( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_169

SUBROUTINE find_equiv_170( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_170

SUBROUTINE find_equiv_171( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_171

SUBROUTINE find_equiv_172( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_172

SUBROUTINE find_equiv_173( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_173

SUBROUTINE find_equiv_174( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_174

SUBROUTINE find_equiv_175( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_175

SUBROUTINE find_equiv_176( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_176

SUBROUTINE find_equiv_177( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_177

SUBROUTINE find_equiv_178( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_178

SUBROUTINE find_equiv_179( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_179

SUBROUTINE find_equiv_180( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_180

SUBROUTINE find_equiv_181( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_181

SUBROUTINE find_equiv_182( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_182

SUBROUTINE find_equiv_183( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_183

SUBROUTINE find_equiv_184( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_184

SUBROUTINE find_equiv_185( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_185

SUBROUTINE find_equiv_186( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_186

SUBROUTINE find_equiv_187( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_187

SUBROUTINE find_equiv_188( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_188

SUBROUTINE find_equiv_189( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_189

SUBROUTINE find_equiv_190( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_190

SUBROUTINE find_equiv_191( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_191

SUBROUTINE find_equiv_192( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_192

SUBROUTINE find_equiv_193( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_193

SUBROUTINE find_equiv_194( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
END SUBROUTINE find_equiv_194

SUBROUTINE find_equiv_195( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_195

SUBROUTINE find_equiv_196( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_196

SUBROUTINE find_equiv_197( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_197

SUBROUTINE find_equiv_198( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_198

SUBROUTINE find_equiv_199( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_199

SUBROUTINE find_equiv_200( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_200

SUBROUTINE find_equiv_201( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_201

SUBROUTINE find_equiv_202( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_202

SUBROUTINE find_equiv_203( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_203

SUBROUTINE find_equiv_204( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_204

SUBROUTINE find_equiv_205( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_205

SUBROUTINE find_equiv_206( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_206

SUBROUTINE find_equiv_207( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_207

SUBROUTINE find_equiv_208( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_208

SUBROUTINE find_equiv_209( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_209

SUBROUTINE find_equiv_210( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_210

SUBROUTINE find_equiv_211( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_211

SUBROUTINE find_equiv_212( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_212

SUBROUTINE find_equiv_213( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_213

SUBROUTINE find_equiv_214( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_214

SUBROUTINE find_equiv_215( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_215

SUBROUTINE find_equiv_216( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_216

SUBROUTINE find_equiv_217( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_217

SUBROUTINE find_equiv_218( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_218

SUBROUTINE find_equiv_219( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_219

SUBROUTINE find_equiv_220( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_220

SUBROUTINE find_equiv_221( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_221

SUBROUTINE find_equiv_222( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_222

SUBROUTINE find_equiv_223( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_223

SUBROUTINE find_equiv_224( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_224

SUBROUTINE find_equiv_225( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_225

SUBROUTINE find_equiv_226( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_226

SUBROUTINE find_equiv_227( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_227

SUBROUTINE find_equiv_228( i, inco, unique, outco )
   CHARACTER(LEN=1), INTENT(in) :: unique
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_228

SUBROUTINE find_equiv_229( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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

END SUBROUTINE find_equiv_229

SUBROUTINE find_equiv_230( i, inco, outco )
   INTEGER,  INTENT(in)  :: i
   REAL(dp), INTENT(in)  :: inco (:,:)
   REAL(dp), INTENT(out) :: outco(:,:,:)
   INTEGER :: k

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
       END SUBROUTINE FIND_EQUIV_230

     END MODULE space_group
