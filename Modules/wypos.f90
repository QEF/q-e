!
! Copyright (C) 2014 Federico Zadra
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE wy_pos
USE kinds,  ONLY : DP
IMPLICIT NONE

SAVE
PRIVATE

PUBLIC wypos

CONTAINS
   SUBROUTINE wypos(tau,wp,inp,space_group_number,uniqueb,&
                                      rhombohedral,origin_choice)

   !-----------------------------------------------------------
   ! Convert atomic positions given in Wyckoff convention:
   ! multiplicity-letter + parameter(s), to crystal positions
   !   wp = Wyckoff label (e.g. 8c)
   !   inp(3) = parameter(s) (if needed)
   !-----------------------------------------------------------

      REAL(DP), DIMENSION(3), INTENT(OUT) :: tau
      REAL(DP), INTENT(IN) :: inp(3)
      CHARACTER(LEN=*), INTENT (IN) :: wp

      INTEGER, INTENT(IN) :: space_group_number
      LOGICAL, INTENT(IN) :: uniqueb, rhombohedral
      INTEGER, INTENT(IN) :: origin_choice

      tau=1.d5

      SELECT CASE (space_group_number)
         CASE (2) !P-1
             CALL wypos_2  ( wp, tau )
         CASE (3) !P2
             CALL wypos_3  ( wp, inp, uniqueb, tau )
         CASE (5) !C2
             CALL wypos_5  ( wp, inp, uniqueb, tau )
         CASE (6) !Pm
             CALL wypos_6  ( wp, inp, uniqueb, tau )
         CASE (8) !Cm
             CALL wypos_8  ( wp, inp, uniqueb, tau )
         CASE (10) !P2/m
             CALL wypos_10 ( wp, inp, uniqueb, tau )
         CASE (11) !P2(1)/m
             CALL wypos_11 ( wp, inp, uniqueb, tau )
         CASE (12) !C2/m
             CALL wypos_12 ( wp, inp, uniqueb, tau )
         CASE (13) !P2/c
             CALL wypos_13 ( wp, inp, uniqueb, tau )
         CASE (14) !-P2(1)/c
             CALL wypos_14 ( wp, inp, uniqueb, tau )
         CASE (15) !C2/c
             CALL wypos_15 ( wp, inp, uniqueb, tau )
         CASE (16) !P222
             CALL wypos_16 ( wp, inp, tau )
         CASE (17) !P222(1)
             CALL wypos_17 ( wp, inp, tau )
         CASE (18) !P2(1)2(1)2
             CALL wypos_18 ( wp, inp, tau )
         CASE (20) !C222(1)
             CALL wypos_20 ( wp, inp, tau )
         CASE (21) !C222
             CALL wypos_21 ( wp, inp, tau )
         CASE (22) !F222
             CALL wypos_22 ( wp, inp, tau )
         CASE (23) !I222
             CALL wypos_23 ( wp, inp, tau )
         CASE (24) !I2(1)2(1)2(1)
             CALL wypos_24 ( wp, inp, tau )
         CASE (25) !Pmm2
             CALL wypos_25 ( wp, inp, tau )
         CASE (26) !Pmc2(1)
             CALL wypos_26 ( wp, inp, tau )
         CASE (27) !Pcc2
             CALL wypos_27 ( wp, inp, tau )
         CASE (28) !Pma2
             CALL wypos_28 ( wp, inp, tau )
         CASE (30) !Pca2(1)
             CALL wypos_30 ( wp, inp, tau )
         CASE (31) !Pmn2(1)
             CALL wypos_31 ( wp, inp, tau )
         CASE (32) !Pba2
             CALL wypos_32 ( wp, inp, tau )
         CASE (34) !Pnn2
             CALL wypos_34 ( wp, inp, tau )
         CASE (35) !Cmm2
             CALL wypos_35 ( wp, inp, tau )
         CASE (36) !Cmc2(1)
             CALL wypos_36 ( wp, inp, tau )
         CASE (37) !Ccc2
             CALL wypos_37 ( wp, inp, tau )
         CASE (38) !Amm2
             CALL wypos_38 ( wp, inp, tau )
         CASE (39) !Aem2
             CALL wypos_39 ( wp, inp, tau )
         CASE (40) !Ama2
             CALL wypos_40 ( wp, inp, tau )
         CASE (41) !Aea2
             CALL wypos_41 ( wp, inp, tau )
         CASE (42) !Fmm2
             CALL wypos_42 ( wp, inp, tau )
         CASE (43) !Fdd2
             CALL wypos_43 ( wp, inp, tau )
         CASE (44) !Imm2
             CALL wypos_44 ( wp, inp, tau )
         CASE (45) !Iba2
             CALL wypos_45 ( wp, inp, tau )
         CASE (46) !Ima2
             CALL wypos_46 ( wp, inp, tau )
         CASE (47) !Pmmm
             CALL wypos_47 ( wp, inp, tau )
         CASE (48) !Pnnn
             CALL wypos_48 ( wp, inp, origin_choice, tau )
         CASE (49) !Pccm
             CALL wypos_49 ( wp, inp, tau )
         CASE (50) !Pban
             CALL wypos_50 ( wp, inp, origin_choice, tau )
         CASE (51) !Pmma
             CALL wypos_51 ( wp, inp, tau )
         CASE (52) !Pnna
             CALL wypos_52 ( wp, inp, tau )
         CASE (53) !Pmna
             CALL wypos_53 ( wp, inp, tau )
         CASE (54) !Pcca
             CALL wypos_54 ( wp, inp, tau )
         CASE (55) !Pbam
             CALL wypos_55 ( wp, inp, tau )
         CASE (56) !Pccn
             CALL wypos_56 ( wp, inp, tau )
         CASE (57) !Pbcm
             CALL wypos_57 ( wp, inp, tau )
         CASE (58) !Pnnm
             CALL wypos_58 ( wp, inp, tau )
         CASE (59) !Pmmn
             CALL wypos_59 ( wp, inp, origin_choice, tau )
         CASE (60) !Pbcn
             CALL wypos_60 ( wp, inp, tau )
         CASE (61) !Pbca
             CALL wypos_61 ( wp, inp, tau )
         CASE (62) !Pnma
             CALL wypos_62 ( wp, inp, tau )
         CASE (63) !Cmcm
             CALL wypos_63 ( wp, inp, tau )
         CASE (64) !Cmce
             CALL wypos_64 ( wp, inp, tau )
         CASE (65) !Cmmm
             CALL wypos_65 ( wp, inp, tau )
         CASE (66) !Cccm
             CALL wypos_66 ( wp, inp, tau )
         CASE (67) !Cmma
             CALL wypos_67 ( wp, inp, tau )
         CASE (68) !Ccce
             CALL wypos_68 ( wp, inp, origin_choice, tau )
         CASE (69) !Fmmm
             CALL wypos_69 ( wp, inp, tau )
         CASE (70) !Fddd
             CALL wypos_70 ( wp, inp, origin_choice, tau )
         CASE (71) !Immm
             CALL wypos_71 ( wp, inp, tau )
         CASE (72) !Ibam
             CALL wypos_72 ( wp, inp, tau )
         CASE (73) !Ibca
             CALL wypos_73 ( wp, inp, tau )
         CASE (74) !Imma
             CALL wypos_74 ( wp, inp, tau )
         CASE (75) !P4
             CALL wypos_75 ( wp, inp, tau )
         CASE (77) !P4(2)
             CALL wypos_77 ( wp, inp, tau )
         CASE (79) !I4(2)
             CALL wypos_79 ( wp, inp, tau )
         CASE (80) !I4(1)
             CALL wypos_80 ( wp, inp, tau )
         CASE (81) !P-4
             CALL wypos_81 ( wp, inp, tau )
         CASE (82) !I-4
             CALL wypos_82 ( wp, inp, tau )
         CASE (83) !P4/m
             CALL wypos_83 ( wp, inp, tau )
         CASE (84)
             CALL wypos_84 ( wp, inp, tau )
         CASE (85)
             CALL wypos_85 ( wp, inp, origin_choice, tau )
         CASE (86)
             CALL wypos_86 ( wp, inp, origin_choice, tau )
         CASE (87) !I4/m
             CALL wypos_87 ( wp, inp, tau )
         CASE (88) !I4(1)/a
             CALL wypos_88 ( wp, inp, origin_choice, tau )
         CASE (89) !P422
             CALL wypos_89 ( wp, inp, tau )
         CASE (90) !P42(1)2
             CALL wypos_90 ( wp, inp, tau )
         CASE (91) !P4(1)22
             CALL wypos_91 ( wp, inp, tau )
         CASE (92) !P4(1)2(1)2
             CALL wypos_92 ( wp, inp, tau )
         CASE (93) !P4(2)22
             CALL wypos_93 ( wp, inp, tau )
         CASE (94) !P4(2)2(1)2
             CALL wypos_94 ( wp, inp, tau )
         CASE (95) !P4(3)22
             CALL wypos_95 ( wp, inp, tau )
         CASE (96) !P4(2)2(1)2
             CALL wypos_96 ( wp, inp, tau )
         CASE (97) !I422
             CALL wypos_97 ( wp, inp, tau )
         CASE (98) !I4(1)22
             CALL wypos_98 ( wp, inp, tau )
         CASE (99) !P4mm
             CALL wypos_99 ( wp, inp, tau )
         CASE (100) !P4bm
             CALL wypos_100( wp, inp, tau )
         CASE (101) !P4(2)cm
             CALL wypos_101( wp, inp, tau )
         CASE (102) !P4(2)nm
             CALL wypos_102( wp, inp, tau )
         CASE (103) !P4cc
             CALL wypos_103( wp, inp, tau )
         CASE (104) !P4nc
             CALL wypos_104( wp, inp, tau )
         CASE (105) !P4(2)mc
             CALL wypos_105( wp, inp, tau )
         CASE (106) !P4(2)bc
             CALL wypos_106( wp, inp, tau )
         CASE (107) !I4mm
             CALL wypos_107( wp, inp, tau )
         CASE (108) !I4cm
             CALL wypos_108( wp, inp, tau )
         CASE (109) !I4(1)md
             CALL wypos_109( wp, inp, tau )
         CASE (110) !I4(1)cd
             CALL wypos_110( wp, inp, tau )
         CASE (111) !P-42m
             CALL wypos_111( wp, inp, tau )
         CASE (112) !P-42c
             CALL wypos_112( wp, inp, tau )
         CASE (113) !P-42(1)m
             CALL wypos_113( wp, inp, tau )
         CASE (114) !P-42(1)c
             CALL wypos_114( wp, inp, tau )
         CASE (115) !P-4m2
             CALL wypos_115( wp, inp, tau )
         CASE (116) !P4c2
             CALL wypos_116( wp, inp, tau )
         CASE (117) !P-4b2
             CALL wypos_117( wp, inp, tau )
         CASE (118) !P-4n2
             CALL wypos_118( wp, inp, tau )
         CASE (119) !I-4m2
             CALL wypos_119( wp, inp, tau )
         CASE (120) !I-4c2
             CALL wypos_120( wp, inp, tau )
         CASE (121) !I-42m
             CALL wypos_121( wp, inp, tau )
         CASE (122) !I-42d
             CALL wypos_122( wp, inp, tau )
         CASE (123) !P4/mmm
             CALL wypos_123( wp, inp, tau )
         CASE (124) !P4/mmc
             CALL wypos_124( wp, inp, tau )
         CASE (125) !P/nbm
             CALL wypos_125( wp, inp, origin_choice, tau )
         CASE (126)
             CALL wypos_126( wp, inp, origin_choice, tau )
         CASE (127) !P4/mbm
             CALL wypos_127( wp, inp, tau )
         CASE (128) !P4/mnc
             CALL wypos_128( wp, inp, tau )
         CASE (129) !P4/nmm
             CALL wypos_129( wp, inp, origin_choice, tau )
         CASE (130) !P4/ncc
             CALL wypos_130( wp, inp, origin_choice, tau )
         CASE (131) !P4(2)/mmc
             CALL wypos_131( wp, inp, tau )
         CASE (132) !P4(2)mcm
             CALL wypos_132( wp, inp, tau )
         CASE (133) !P4(2)/nbc
             CALL wypos_133( wp, inp, origin_choice, tau )
         CASE (134) !P4(2)/nnm
             CALL wypos_134( wp, inp, origin_choice, tau )
         CASE (135) !P3(2)/mbc
             CALL wypos_135( wp, inp, tau )
         CASE (136) !P4(2)/mnm
             CALL wypos_136( wp, inp, tau )
         CASE (137) !P4(2)/nmc
             CALL wypos_137( wp, inp, origin_choice, tau )
         CASE (138) !P4(2)/ncm
             CALL wypos_138( wp, inp, origin_choice, tau )
         CASE (139) !I4/mmm
             CALL wypos_139( wp, inp, tau )
         CASE (140) !I4/mcm
             CALL wypos_140( wp, inp, tau )
         CASE (141) !I4(1)/amd
             CALL wypos_141( wp, inp, origin_choice, tau )
         CASE (142) !I4(1)/acd
             CALL wypos_142( wp, inp, origin_choice, tau )
         CASE(143) !P3
             CALL wypos_143( wp, inp, tau )
         CASE (146) !R3
             CALL wypos_146( wp, inp, rhombohedral, tau )
         CASE(147) !P-3
             CALL wypos_147( wp, inp, tau )
         CASE (148) !R-3
             CALL wypos_148( wp, inp, rhombohedral, tau )
         CASE (149) !P312
             CALL wypos_149( wp, inp, tau )
         CASE (150) !P321
             CALL wypos_150( wp, inp, tau )
         CASE (151) !P3(1)12
             CALL wypos_151( wp, inp, tau )
         CASE (152) !P3(1)21
             CALL wypos_152( wp, inp, tau )
         CASE (153) !P3(2)12
             CALL wypos_153( wp, inp, tau )
         CASE (154) !3(2)21
             CALL wypos_154( wp, inp, tau )
         CASE (155) !R32
             CALL wypos_155( wp, inp, rhombohedral, tau )
         CASE (156) !P-3m1
             CALL wypos_156( wp, inp, tau )
         CASE (157) !P31m
             CALL wypos_157( wp, inp, tau )
         CASE (158) !P3c1
             CALL wypos_158( wp, inp, tau )
         CASE (159) !P31c
             CALL wypos_159( wp, inp, tau )
         CASE (160) !R3m
             CALL wypos_160( wp, inp, rhombohedral, tau )
         CASE (161) !R3c
             CALL wypos_161( wp, inp, rhombohedral, tau )
         CASE (162) !P-31m
             CALL wypos_162( wp, inp, tau )
         CASE (163) !P-31c
             CALL wypos_163( wp, inp, tau )
         CASE (164) !P-3m1
             CALL wypos_164( wp, inp, tau )
         CASE (165) !P-3c1
             CALL wypos_165( wp, inp, tau )
         CASE (166) !R-3m
             CALL wypos_166( wp, inp, rhombohedral, tau )
         CASE (167) !R-3c
             CALL wypos_167( wp, inp, rhombohedral, tau )
         CASE (168) !P6
             CALL wypos_168( wp, inp, tau )
         CASE (171) !P6/m
             CALL wypos_171( wp, inp, tau )
         CASE (172) !P6(4)
             CALL wypos_172( wp, inp, tau )
         CASE (173) !P6(3)
             CALL wypos_173( wp, inp, tau )
         CASE (174) !P-6
             CALL wypos_174( wp, inp, tau )
         CASE (175) !P6/m
             CALL wypos_175( wp, inp, tau )
         CASE (176) !P6(3)/m
             CALL wypos_176( wp, inp, tau )
         CASE (177) !P622
             CALL wypos_177( wp, inp, tau )
         CASE (178) !P6(1)22
             CALL wypos_178( wp, inp, tau )
         CASE (179) !P6(5)22
             CALL wypos_179( wp, inp, tau )
         CASE (180) !P6(2)22
             CALL wypos_180( wp, inp, tau )
         CASE (181) !P6(4)22
             CALL wypos_181( wp, inp, tau )
         CASE (182) !P6(3)22
             CALL wypos_182( wp, inp, tau )
         CASE (183) !P6mm
             CALL wypos_183( wp, inp, tau )
         CASE (184) !P6cc
             CALL wypos_184( wp, inp, tau )
         CASE (185) !P6(3)cm
             CALL wypos_185( wp, inp, tau )
         CASE (186) !P6(3)mc
             CALL wypos_186( wp, inp, tau )
         CASE (187) !P-6m2
             CALL wypos_187( wp, inp, tau )
         CASE (188) !P-6c2
             CALL wypos_188( wp, inp, tau )
         CASE (189) !P-62m
             CALL wypos_189( wp, inp, tau )
         CASE (190) !P-62c
             CALL wypos_190( wp, inp, tau )
         CASE (191) !P6/mmm
             CALL wypos_191( wp, inp, tau )
         CASE (192) !P6/mcc
             CALL wypos_192( wp, inp, tau )
         CASE (193) !P6(3)/mcm
             CALL wypos_193( wp, inp, tau )
         CASE (194) !P6(3)mmc
             CALL wypos_194( wp, inp, tau )
         CASE (195) !P23
             CALL wypos_195( wp, inp, tau )
         CASE (196) !F23
             CALL wypos_196( wp, inp, tau )
         CASE (197) !I23
             CALL wypos_197( wp, inp, tau )
         CASE (198) !P2(1)3
             CALL wypos_198( wp, inp, tau )
         CASE (199) !I2(1)3
             CALL wypos_199( wp, inp, tau )
         CASE (200) !Pm-3
             CALL wypos_200( wp, inp, tau )
         CASE (201) !Pn-3
             CALL wypos_201( wp, inp, origin_choice, tau )
         CASE (202) !Fm-3
             CALL wypos_202( wp, inp, tau )
         CASE (203) !Fd-3
             CALL wypos_203( wp, inp, origin_choice, tau )
         CASE (204) ! Im-3
             CALL wypos_204( wp, inp, tau )
         CASE (205) !Pa-3
             CALL wypos_205( wp, inp, tau )
         CASE (206) !Ia-3
             CALL wypos_206( wp, inp, tau )
         CASE (207) !P432
             CALL wypos_207( wp, inp, tau )
         CASE (208) !P4(2)32
             CALL wypos_208( wp, inp, tau )
         CASE (209) !F432
             CALL wypos_209( wp, inp, tau )
         CASE (210) !F4(1)32
             CALL wypos_210( wp, inp, tau )
         CASE (211) !I432
             CALL wypos_211( wp, inp, tau )
         CASE (212) !P4(3)32
             CALL wypos_212( wp, inp, tau )
         CASE (213) !P4(1)32
             CALL wypos_213( wp, inp, tau )
         CASE (214) !I4(I)32
             CALL wypos_214( wp, inp, tau )
         CASE (215) !P-43m
             CALL wypos_215( wp, inp, tau )
         CASE (216) !F-43m
             CALL wypos_216( wp, inp, tau )
         CASE (217) !I-43m
             CALL wypos_217( wp, inp, tau )
         CASE (218) !P-43n
             CALL wypos_218( wp, inp, tau )
         CASE (219) !F-43c
             CALL wypos_219( wp, inp, tau )
         CASE (220) !I-43d
             CALL wypos_220( wp, inp, tau )
         CASE (221) !Pm-3m
             CALL wypos_221( wp, inp, tau )
         CASE (222) !Pn-3n
             CALL wypos_222( wp, inp, origin_choice, tau )
         CASE (223) !Pm-3n
             CALL wypos_223( wp, inp, tau )
         CASE (224) !Pn-3m
             CALL wypos_224( wp, inp, origin_choice, tau )
         CASE (225) !Fm-3m
             CALL wypos_225( wp, inp, tau )
         CASE (226) !Fm-3c
             CALL wypos_226( wp, inp, tau )
         CASE (227) !Fd-3m
             CALL wypos_227( wp, inp, origin_choice, tau )
         CASE (228) !Fd-3c
             CALL wypos_228( wp, inp, origin_choice, tau )
         CASE (229) !Im-3m
             CALL wypos_229( wp, inp, tau )
         CASE (230) !Ia-3d
             CALL wypos_230( wp, inp, tau )
          CASE DEFAULT
            CALL errore('wypos','group not recognized',1)
          END SELECT

         IF (tau(1)==1.d5.OR.tau(2)==1.d5.OR.tau(3)==1.d5) THEN
            IF (inp(1)==1.d5.OR.inp(2)==1.d5.OR.inp(3)==1.d5) THEN
               CALL errore('wypos','wyckoff position not found',1)
            ELSE
               CALL infomsg('wypos','wyckoff position not found, assuming x y z')
               tau(:)=inp(:)
            END IF
         END IF

      END SUBROUTINE wypos

SUBROUTINE wypos_2  ( wp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ENDIF

          END SUBROUTINE wypos_2

SUBROUTINE wypos_3  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)

   IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

          END SUBROUTINE wypos_3

SUBROUTINE wypos_5  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

          END SUBROUTINE wypos_5

SUBROUTINE wypos_6  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.5_DP
               ENDIF
            ENDIF

          END SUBROUTINE wypos_6

SUBROUTINE wypos_8  ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ENDIF
            ENDIF

          END SUBROUTINE wypos_8

SUBROUTINE wypos_10 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2k') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2l') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2m') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='2n') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='1g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2i') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2j') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2k') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2l') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2n') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.5_DP
               ENDIF
            ENDIF

          END SUBROUTINE wypos_10

SUBROUTINE wypos_11 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

          END SUBROUTINE wypos_11

SUBROUTINE wypos_12 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (uniqueb) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN

               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(2)
                  tau(3)=0.0_DP
               ENDIF
            ENDIF

          END SUBROUTINE wypos_12

SUBROUTINE wypos_13 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN

               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2f') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN

               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF

            ENDIF

          END SUBROUTINE wypos_13

SUBROUTINE wypos_14 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN

               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN

               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ENDIF

            ENDIF

          END SUBROUTINE wypos_14

SUBROUTINE wypos_15 ( wp, inp, uniqueb, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: uniqueb
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (uniqueb) THEN

               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25
               ENDIF

            ELSEIF (.NOT.uniqueb) THEN

               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF

            ENDIF

          END SUBROUTINE wypos_15

SUBROUTINE wypos_16 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2o') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2q') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2r') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2s') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2t') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF


          END SUBROUTINE wypos_16

SUBROUTINE wypos_17 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF

          END SUBROUTINE wypos_17

SUBROUTINE wypos_18 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_18

SUBROUTINE wypos_20 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF

          END SUBROUTINE wypos_20

SUBROUTINE wypos_21 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)

            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_21

SUBROUTINE wypos_22 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF

          END SUBROUTINE wypos_22

SUBROUTINE wypos_23 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_23

SUBROUTINE wypos_24 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_24

SUBROUTINE wypos_25 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

          END SUBROUTINE wypos_25

SUBROUTINE wypos_26 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

          END SUBROUTINE wypos_26

SUBROUTINE wypos_27 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_27

SUBROUTINE wypos_28 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

          END SUBROUTINE wypos_28

SUBROUTINE wypos_30 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_30

SUBROUTINE wypos_31 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF
          END SUBROUTINE wypos_31

SUBROUTINE wypos_32 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_32

SUBROUTINE wypos_34 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

          END SUBROUTINE wypos_34

SUBROUTINE wypos_35 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

          END SUBROUTINE wypos_35

SUBROUTINE wypos_36 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

          END SUBROUTINE wypos_36

SUBROUTINE wypos_37 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_37 

SUBROUTINE wypos_38 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_38 

SUBROUTINE wypos_39 ( wp, inp, tau )
  CHARACTER(LEN=*), INTENT(in)  :: wp
  REAL(dp), INTENT(in) :: inp(3)
  REAL(dp), INTENT(out) :: tau (3)
  
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_39 

SUBROUTINE wypos_40 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_40 

SUBROUTINE wypos_41 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_41 

SUBROUTINE wypos_42 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_42 

SUBROUTINE wypos_43 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_43 

SUBROUTINE wypos_44 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_44 

SUBROUTINE wypos_45 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_45 

SUBROUTINE wypos_46 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_46 

SUBROUTINE wypos_47 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2o') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2q') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2r') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2s') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2t') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4u') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4v') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4w') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4x') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4y') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4z') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_47 

SUBROUTINE wypos_48 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.75_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_48 

SUBROUTINE wypos_49 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4p') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_49 

SUBROUTINE wypos_50 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4j') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4k') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4l') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_50 

SUBROUTINE wypos_51 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_51 

SUBROUTINE wypos_52 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_52 

SUBROUTINE wypos_53 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_53 

SUBROUTINE wypos_54 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_54 

SUBROUTINE wypos_55 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_55 

SUBROUTINE wypos_56 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.75_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_56 

SUBROUTINE wypos_57 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_57 

SUBROUTINE wypos_58 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF



END SUBROUTINE wypos_58 

SUBROUTINE wypos_59 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_59 

SUBROUTINE wypos_60 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ENDIF


END SUBROUTINE wypos_60 

SUBROUTINE wypos_61 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF


END SUBROUTINE wypos_61 

SUBROUTINE wypos_62 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_62 

SUBROUTINE wypos_63 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_63 

SUBROUTINE wypos_64 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_64 

SUBROUTINE wypos_65 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF


END SUBROUTINE wypos_65 

SUBROUTINE wypos_66 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.75_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_66 

SUBROUTINE wypos_67 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_67 

SUBROUTINE wypos_68 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_68 

SUBROUTINE wypos_69 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='16j') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16m') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16o') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_69 

SUBROUTINE wypos_70 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF


            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_70 

SUBROUTINE wypos_71 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_71 

SUBROUTINE wypos_72 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_72 

SUBROUTINE wypos_73 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ENDIF


END SUBROUTINE wypos_73 

SUBROUTINE wypos_74 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_74 

SUBROUTINE wypos_75 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_75 

SUBROUTINE wypos_77 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_77 

SUBROUTINE wypos_79 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_79 

SUBROUTINE wypos_80 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_80 

SUBROUTINE wypos_81 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_81 

SUBROUTINE wypos_82 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_82 

SUBROUTINE wypos_83 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_83 

SUBROUTINE wypos_84 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_84 

SUBROUTINE wypos_85 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_85 

SUBROUTINE wypos_86 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_86 

SUBROUTINE wypos_87 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF


END SUBROUTINE wypos_87 

SUBROUTINE wypos_88 ( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_88 

SUBROUTINE wypos_89 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_89 

SUBROUTINE wypos_90 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_90 

SUBROUTINE wypos_91 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.375_DP
            ENDIF

END SUBROUTINE wypos_91 

SUBROUTINE wypos_92 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_92 

SUBROUTINE wypos_93 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.75_DP
            ENDIF

END SUBROUTINE wypos_93 

SUBROUTINE wypos_94 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_94 

SUBROUTINE wypos_95 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.625_DP
            ENDIF

END SUBROUTINE wypos_95 

SUBROUTINE wypos_96 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_96 

SUBROUTINE wypos_97 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_97 

SUBROUTINE wypos_98 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=-inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.125_DP
            ENDIF

END SUBROUTINE wypos_98 

SUBROUTINE wypos_99 ( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_99 

SUBROUTINE wypos_100( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_100

SUBROUTINE wypos_101( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_101

SUBROUTINE wypos_102( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_102

SUBROUTINE wypos_103( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_103

SUBROUTINE wypos_104( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_104

SUBROUTINE wypos_105( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_105

SUBROUTINE wypos_106( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_106

SUBROUTINE wypos_107( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_107

SUBROUTINE wypos_108( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_108

SUBROUTINE wypos_109( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_109

SUBROUTINE wypos_110( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_110

SUBROUTINE wypos_111( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_111

SUBROUTINE wypos_112( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF


END SUBROUTINE wypos_112

SUBROUTINE wypos_113( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_113

SUBROUTINE wypos_114( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF


END SUBROUTINE wypos_114

SUBROUTINE wypos_115( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_115

SUBROUTINE wypos_116( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF


END SUBROUTINE wypos_116

SUBROUTINE wypos_117( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.5_DP
            ENDIF


END SUBROUTINE wypos_117

SUBROUTINE wypos_118( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_118

SUBROUTINE wypos_119( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_119

SUBROUTINE wypos_120( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_120

SUBROUTINE wypos_121( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_121

SUBROUTINE wypos_122( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8d') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.125_DP
            ENDIF


END SUBROUTINE wypos_122

SUBROUTINE wypos_123( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)

            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4n') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4o') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8r') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8s') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8t') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_123

SUBROUTINE wypos_124( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF


END SUBROUTINE wypos_124

SUBROUTINE wypos_125( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2d') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4h') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_125

SUBROUTINE wypos_126( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)

   IF (origin_choice==1) THEN

               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_126

SUBROUTINE wypos_127( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_127

SUBROUTINE wypos_128( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_128

SUBROUTINE wypos_129( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF


            ELSEIF (origin_choice==2) THEN


               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_129

SUBROUTINE wypos_130( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_130

SUBROUTINE wypos_131( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4k') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8p') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='8q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF


END SUBROUTINE wypos_131

SUBROUTINE wypos_132( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8k') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8l') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8m') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8n') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8o') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_132

SUBROUTINE wypos_133( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.0_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.00_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_133

SUBROUTINE wypos_134( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4g') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8j') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8l') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8m') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF

            ENDIF

END SUBROUTINE wypos_134

SUBROUTINE wypos_135( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_135

SUBROUTINE wypos_136( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_136

SUBROUTINE wypos_137( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_137

SUBROUTINE wypos_138( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.5_DP
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='4d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4e') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8f') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='8g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8i') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_138

SUBROUTINE wypos_139( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16m') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='16n') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_139

SUBROUTINE wypos_140( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8f') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16l') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)+0.5_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_140

SUBROUTINE wypos_141( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='4a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.25_DP
                  tau(3)=0.875_DP
               ELSEIF (TRIM(wp)=='16h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ENDIF

END SUBROUTINE wypos_141

SUBROUTINE wypos_142( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=0.25_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.25_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='16e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)+0.25_DP
                  tau(3)=0.125_DP
               ENDIF

            ENDIF

END SUBROUTINE wypos_142

SUBROUTINE wypos_143( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_143

SUBROUTINE wypos_146( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
            ELSE !If HEXAGONAL
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_146

SUBROUTINE wypos_147( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF


END SUBROUTINE wypos_147

SUBROUTINE wypos_148( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_148

SUBROUTINE wypos_149( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_149

SUBROUTINE wypos_150( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_150

SUBROUTINE wypos_151( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=1.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=5.0_DP/6.0_DP
            ENDIF

END SUBROUTINE wypos_151

SUBROUTINE wypos_152( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=1.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=5.0_DP/6.0_DP
            ENDIF

END SUBROUTINE wypos_152

SUBROUTINE wypos_153( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=2.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=1.0_DP/6.0_DP
            ENDIF

END SUBROUTINE wypos_153

SUBROUTINE wypos_154( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=2.0_DP/3.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=1.0_DP/6.0_DP
            ENDIF

END SUBROUTINE wypos_154

SUBROUTINE wypos_155( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_155

SUBROUTINE wypos_156( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_156

SUBROUTINE wypos_157( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_157

SUBROUTINE wypos_158( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_158

SUBROUTINE wypos_159( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_159

SUBROUTINE wypos_160( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9b') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_160

SUBROUTINE wypos_161( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_161

SUBROUTINE wypos_162( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_162

SUBROUTINE wypos_163( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_163

SUBROUTINE wypos_164( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_164

SUBROUTINE wypos_165( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_165

SUBROUTINE wypos_166( wp, inp, rhombohedral, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (rhombohedral) THEN
               IF (TRIM(wp)=='1a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='1b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='2c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='3d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3e') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6f') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6g') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6h') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='3a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='3b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='9d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='9e') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='18h') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_166

SUBROUTINE wypos_167( wp, inp, rhombohedral,  tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   LOGICAL, INTENT(in) :: rhombohedral
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (rhombohedral) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='2b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6e') THEN
                  tau(1)=inp(1)
                  tau(2)=-inp(1)+0.5_DP
                  tau(3)=0.25_DP
               ENDIF

            ELSEIF (.NOT.rhombohedral) THEN
               IF (TRIM(wp)=='6a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='18d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='18e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_167

SUBROUTINE wypos_168( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_168

SUBROUTINE wypos_171( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_171

SUBROUTINE wypos_172( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_172

SUBROUTINE wypos_173( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_173

SUBROUTINE wypos_174( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_174

SUBROUTINE wypos_175( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_175

SUBROUTINE wypos_176( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_176

SUBROUTINE wypos_177( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_177

SUBROUTINE wypos_178( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='6a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_178

SUBROUTINE wypos_179( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
            IF (TRIM(wp)=='6a') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.75_DP
            ENDIF

END SUBROUTINE wypos_179

SUBROUTINE wypos_180( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_180

SUBROUTINE wypos_181( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='3a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_181

SUBROUTINE wypos_182( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_182

SUBROUTINE wypos_183( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_183

SUBROUTINE wypos_184( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_184

SUBROUTINE wypos_185( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_185

SUBROUTINE wypos_186( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_186

SUBROUTINE wypos_187( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='1e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='2i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3k') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6n') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_187

SUBROUTINE wypos_188( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2f') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4g') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4i') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=-inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_188

SUBROUTINE wypos_189( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_189

SUBROUTINE wypos_190( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=2.0_DP/3.0_DP
               tau(2)=1.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ENDIF


END SUBROUTINE wypos_190

SUBROUTINE wypos_191( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='2e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='3f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='3g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6l') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6m') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12n') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12o') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12p') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12q') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_191

SUBROUTINE wypos_192( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12l') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_192

SUBROUTINE wypos_193( wp, inp, tau )
   
   REAL(dp), INTENT(in) :: inp(3)
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8h') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_193

SUBROUTINE wypos_194( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='2b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2c') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='2d') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='4f') THEN
               tau(1)=1.0_DP/3.0_DP
               tau(2)=2.0_DP/3.0_DP
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=inp(2)
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=inp(1)
               tau(2)=2.0_DP*inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_194

SUBROUTINE wypos_195( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6i') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_195

SUBROUTINE wypos_196( wp, inp, tau )
  CHARACTER(LEN=*), INTENT(in)  :: wp
  REAL(dp), INTENT(in) :: inp(3)
  REAL(dp), INTENT(out) :: tau (3)
  
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_196

SUBROUTINE wypos_197( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ENDIF

END SUBROUTINE wypos_197

SUBROUTINE wypos_198( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_198

SUBROUTINE wypos_199( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12b') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_199

SUBROUTINE wypos_200( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_200

SUBROUTINE wypos_201( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.5_DP
                  tau(3)=0.0_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_201

SUBROUTINE wypos_202( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_202

SUBROUTINE wypos_203( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ENDIF
            ENDIF

END SUBROUTINE wypos_203

SUBROUTINE wypos_204( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_204

SUBROUTINE wypos_205( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF

END SUBROUTINE wypos_205

SUBROUTINE wypos_206( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_206

SUBROUTINE wypos_207( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ENDIF


END SUBROUTINE wypos_207

SUBROUTINE wypos_208( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12k') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='12l') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.5_DP
            ENDIF

END SUBROUTINE wypos_208

SUBROUTINE wypos_209( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_209

SUBROUTINE wypos_210( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='16d') THEN
               tau(1)=0.625_DP
               tau(2)=0.625_DP
               tau(3)=0.625_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF

END SUBROUTINE wypos_210

SUBROUTINE wypos_211( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ENDIF

END SUBROUTINE wypos_211

SUBROUTINE wypos_212( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.625_DP
               tau(2)=0.625_DP
               tau(3)=0.625_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF

END SUBROUTINE wypos_212

SUBROUTINE wypos_213( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.375_DP
               tau(2)=0.375_DP
               tau(3)=0.375_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.875_DP
               tau(2)=0.875_DP
               tau(3)=0.875_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.25_DP
            ENDIF

END SUBROUTINE wypos_213

SUBROUTINE wypos_214( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.875_DP
               tau(2)=0.875_DP
               tau(3)=0.875_DP
            ELSEIF (TRIM(wp)=='12c') THEN
               tau(1)=0.125_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.625_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.25_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF


END SUBROUTINE wypos_214

SUBROUTINE wypos_215( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_215

SUBROUTINE wypos_216( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='4c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='4d') THEN
               tau(1)=0.75_DP
               tau(2)=0.75_DP
               tau(3)=0.75_DP
            ELSEIF (TRIM(wp)=='16e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_216

SUBROUTINE wypos_217( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_217

SUBROUTINE wypos_218( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12g') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ENDIF

END SUBROUTINE wypos_218

SUBROUTINE wypos_219( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_219

SUBROUTINE wypos_220( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='12a') THEN
               tau(1)=0.375_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12b') THEN
               tau(1)=0.875_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='16c') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ENDIF

END SUBROUTINE wypos_220

SUBROUTINE wypos_221( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='1a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='1b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3c') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='3d') THEN
               tau(1)=0.5_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6f') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='12j') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='24l') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='24m') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF


END SUBROUTINE wypos_221

SUBROUTINE wypos_222( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='12e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='24g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='6b') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='8c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='12d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='12e') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='16f') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='24g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.75_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_222

SUBROUTINE wypos_223( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='6d') THEN
               tau(1)=0.25_DP
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='8e') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='12g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12h') THEN
               tau(1)=inp(1)
               tau(2)=0.5_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16i') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24j') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='24k') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_223

SUBROUTINE wypos_224( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
           IF (origin_choice==1) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.75_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='24i') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24j') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='2a') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='4b') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='4c') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='6d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.75_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='8e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='12f') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='12g') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='24h') THEN
                  tau(1)=inp(1)
                  tau(2)=0.25_DP
                  tau(3)=0.75_DP
               ELSEIF (TRIM(wp)=='24i') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=inp(1)+0.5_DP
               ELSEIF (TRIM(wp)=='24j') THEN
                  tau(1)=0.5_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ELSEIF (TRIM(wp)=='24k') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ENDIF
            ENDIF

END SUBROUTINE wypos_224

SUBROUTINE wypos_225( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='4a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='4b') THEN
               tau(1)=0.5_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='32f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=0.5_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='96k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_225

SUBROUTINE wypos_226( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp 
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='8a') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='8b') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.0_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='64g') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96h') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='96i') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_226

SUBROUTINE wypos_227( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.625_DP
                  tau(2)=0.625_DP
                  tau(3)=0.625_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='96h') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.25_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='8a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='8b') THEN
                  tau(1)=0.375_DP
                  tau(2)=0.375_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='16c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='16d') THEN
                  tau(1)=0.5_DP
                  tau(2)=0.5_DP
                  tau(3)=0.5_DP
               ELSEIF (TRIM(wp)=='32e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='48f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(2)
               ELSEIF (TRIM(wp)=='96h') THEN
                  tau(1)=0.0_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_227

SUBROUTINE wypos_228( wp, inp, origin_choice, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   INTEGER, INTENT(in) :: origin_choice
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (origin_choice==1) THEN
               IF (TRIM(wp)=='16a') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='32b') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='32c') THEN
                  tau(1)=0.375_DP
                  tau(2)=0.375_DP
                  tau(3)=0.375_DP
               ELSEIF (TRIM(wp)=='48d') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='64e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='96f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=0.125_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)+0.25_DP
               ENDIF

            ELSEIF (origin_choice==2) THEN
               IF (TRIM(wp)=='16a') THEN
                  tau(1)=0.125_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='32b') THEN
                  tau(1)=0.25_DP
                  tau(2)=0.25_DP
                  tau(3)=0.25_DP
               ELSEIF (TRIM(wp)=='32c') THEN
                  tau(1)=0.0_DP
                  tau(2)=0.0_DP
                  tau(3)=0.0_DP
               ELSEIF (TRIM(wp)=='48d') THEN
                  tau(1)=0.875_DP
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='64e') THEN
                  tau(1)=inp(1)
                  tau(2)=inp(1)
                  tau(3)=inp(1)
               ELSEIF (TRIM(wp)=='96f') THEN
                  tau(1)=inp(1)
                  tau(2)=0.125_DP
                  tau(3)=0.125_DP
               ELSEIF (TRIM(wp)=='96g') THEN
                  tau(1)=0.25_DP
                  tau(2)=inp(1)
                  tau(3)=-inp(1)
               ENDIF
            ENDIF

END SUBROUTINE wypos_228

SUBROUTINE wypos_229( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='2a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='6b') THEN
               tau(1)=0.0_DP
               tau(2)=0.5_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='8c') THEN
               tau(1)=0.25_DP
               tau(2)=0.25_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='12d') THEN
               tau(1)=0.25_DP
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='12e') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16f') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='24g') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.5_DP
            ELSEIF (TRIM(wp)=='24h') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48i') THEN
               tau(1)=0.25_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.5_DP
            ELSEIF (TRIM(wp)=='48j') THEN
               tau(1)=0.0_DP
               tau(2)=inp(1)
               tau(3)=inp(2)
            ELSEIF (TRIM(wp)=='48k') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(2)
            ENDIF

END SUBROUTINE wypos_229

SUBROUTINE wypos_230( wp, inp, tau )
   CHARACTER(LEN=*), INTENT(in)  :: wp
   REAL(dp), INTENT(in) :: inp(3)
   REAL(dp), INTENT(out) :: tau (3)
   
            IF (TRIM(wp)=='16a') THEN
               tau(1)=0.0_DP
               tau(2)=0.0_DP
               tau(3)=0.0_DP
            ELSEIF (TRIM(wp)=='16b') THEN
               tau(1)=0.125_DP
               tau(2)=0.125_DP
               tau(3)=0.125_DP
            ELSEIF (TRIM(wp)=='24c') THEN
               tau(1)=0.125_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='24d') THEN
               tau(1)=0.375_DP
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='32e') THEN
               tau(1)=inp(1)
               tau(2)=inp(1)
               tau(3)=inp(1)
            ELSEIF (TRIM(wp)=='48f') THEN
               tau(1)=inp(1)
               tau(2)=0.0_DP
               tau(3)=0.25_DP
            ELSEIF (TRIM(wp)=='48g') THEN
               tau(1)=0.125_DP
               tau(2)=inp(1)
               tau(3)=-inp(1)+0.25_DP
            ENDIF
          END SUBROUTINE wypos_230
END MODULE
