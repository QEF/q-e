!
! Copyright (C) 2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------
MODULE qe_dft_refs
  !----------------------------------------------------------------------------------
  !! List of references and brief descriptions of available DFTs in QE internal library.
  ! If a functional is never called in QE, please start the %wrn (warning) message
  ! with 'never called' so that the testing program skips it.
  !
  USE qe_dft_list, ONLY: nxc, ncc, ngcx, ngcc, nmeta, n_dft
  !
  SAVE
  !
  TYPE dft_refwrn
     CHARACTER(LEN=150) :: ref
     CHARACTER(LEN=100) :: wrn
  END TYPE dft_refwrn
  !
  ! -- single DFT terms (family-type)
  TYPE(dft_refwrn) :: dft_LDAx(0:nxc),  dft_LDAc(0:ncc),  &
                      dft_GGAx(0:ngcx), dft_GGAc(0:ngcc), &
                      dft_MGGA(0:nmeta)
  ! -- total DFTs
  CHARACTER(LEN=100) :: dft_full_descr(n_dft)
  !
  ! ---- LDA exchange ----
  !
  ! NOX
  DATA dft_LDAx(0)%ref  / 'No LDA exchange.' /
  DATA dft_LDAx(0)%wrn  / 'none' /
  ! SLA
  DATA dft_LDAx(1)%ref  / '[Slater exchange - alpha=2/3]' /
  DATA dft_LDAx(1)%wrn  / 'none' /
  ! SL1
  DATA dft_LDAx(2)%ref  / '[Slater exchange - alpha=1.0]' /
  DATA dft_LDAx(2)%wrn  / 'none' /
  ! RXC
  DATA dft_LDAx(3)%ref  / '[Relativistic Slater]' /
  DATA dft_LDAx(3)%wrn  / 'none' /
  ! OEP
  DATA dft_LDAx(4)%ref  / '[Optimized Effective Potential]' /
  DATA dft_LDAx(4)%wrn  / 'none' /
  ! HF
  DATA dft_LDAx(5)%ref  / '[Hartree-Fock]' /
  DATA dft_LDAx(5)%wrn  / 'none' /
  ! PB0X  (Slater*0.75+HF*0.25) for PBE0 and vdW-DF-cx0 and vdW-DF2-0 etc
  DATA dft_LDAx(6)%ref  / 'J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)' /
  DATA dft_LDAx(6)%wrn  / 'for PBE0 and vdW-DF-cx0 and vdW-DF2-0 etc.' /
  ! B3LP  (Slater*0.80+HF*0.20)
  DATA dft_LDAx(7)%ref  / 'P.J.Stephens, F.J.Devlin, C.F.Chabalowski, M.J.Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  DATA dft_LDAx(7)%wrn  / 'none' /
  ! KZK   Finite-size corrections
  DATA dft_LDAx(8)%ref  / 'H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)' /
  DATA dft_LDAx(8)%wrn  / 'none' /
  ! xxxx [X3LYP_LDA]
  DATA dft_LDAx(9)%ref  / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  DATA dft_LDAx(9)%wrn  / 'none' /
  ! xxxx [KLI]
  DATA dft_LDAx(10)%ref / 'KLI approximation for exx - no ref. available' /
  DATA dft_LDAx(10)%wrn / 'Currently not implemented' /
  !
  !
  !  ---- LDA correlation ----
  !
  ! NOC
  DATA dft_LDAc(0)%ref  / 'No LDA correlation.' /
  DATA dft_LDAc(0)%wrn  / 'none' /
  ! PZ
  DATA dft_LDAc(1)%ref  / 'J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)' /
  DATA dft_LDAc(1)%wrn  / 'none' /
  ! VWN
  DATA dft_LDAc(2)%ref  / 'S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)' /
  DATA dft_LDAc(2)%wrn  / 'none' /
  ! LYP
  DATA dft_LDAc(3)%ref  / 'C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)' /
  DATA dft_LDAc(3)%wrn  / 'none' /
  ! PW
  DATA dft_LDAc(4)%ref  / 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)' /
  DATA dft_LDAc(4)%wrn  / 'none' /
  ! WIG
  DATA dft_LDAc(5)%ref  / 'E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)' /
  DATA dft_LDAc(5)%wrn  / 'none' /
  ! HL
  DATA dft_LDAc(6)%ref  / 'L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)' /
  DATA dft_LDAc(6)%wrn  / 'none' /
  ! OBZ
  DATA dft_LDAc(7)%ref  / 'G.Ortiz and P.Ballone, PRB 50, 1391 (1994)' /
  DATA dft_LDAc(7)%wrn  / 'none' /
  ! OBW
  DATA dft_LDAc(8)%ref  / 'G.Ortiz and P.Ballone, PRB 50, 1391 (1994)' /
  DATA dft_LDAc(8)%wrn  / 'none' /
  ! GL
  DATA dft_LDAc(9)%ref  / 'O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)' /
  DATA dft_LDAc(9)%wrn  / 'none' /
  ! KZK
  DATA dft_LDAc(10)%ref / 'H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)' /
  DATA dft_LDAc(10)%wrn / 'none' /
  ! xxxx [vwn1_rpa]
  DATA dft_LDAc(11)%ref / 'vwn1_rpa' /
  DATA dft_LDAc(11)%wrn / 'none' /
  ! B3LP
  DATA dft_LDAc(12)%ref / 'P.J.Stephens, F.J.Devlin, C.F.Chabalowski, M.J.Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  DATA dft_LDAc(12)%wrn / 'none' /
  ! xxxx [B3LYP-V1R]
  DATA dft_LDAc(13)%ref / 'B3LYP-V1R' /
  DATA dft_LDAc(13)%wrn / 'none' /
  ! xxxx [X3LYP]
  DATA dft_LDAc(14)%ref / 'X3LYP' /
  DATA dft_LDAc(14)%wrn / 'none' /
  !
  !
  !  ---- GGA exchange ----
  !
  ! NOGX
  DATA dft_GGAx(0)%ref  / 'No GGA exchange.' /
  DATA dft_GGAx(0)%wrn  / 'none' /
  ! B88
  DATA dft_GGAx(1)%ref  / 'A.D.Becke, PRA 38, 3098 (1988)' /
  DATA dft_GGAx(1)%wrn  / 'none' /
  ! GGX
  DATA dft_GGAx(2)%ref  / 'J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)' /
  DATA dft_GGAx(2)%wrn  / 'none' /
  ! PBX
  DATA dft_GGAx(3)%ref  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  DATA dft_GGAx(3)%wrn  / 'none' /
  ! REVX
  DATA dft_GGAx(4)%ref  / 'Zhang and Yang, PRL 80, 890 (1998)' /
  DATA dft_GGAx(4)%wrn  / 'none' /
  ! HCTH
  DATA dft_GGAx(5)%ref  / 'Handy et al, JCP 109, 6264 (1998)' /
  DATA dft_GGAx(5)%wrn  / 'none' /
  ! OPTX
  DATA dft_GGAx(6)%ref  / 'Handy et al, JCP 116, 5411 (2002)' /
  DATA dft_GGAx(6)%wrn  / 'OPTX untested! please test' /
  ! void
  DATA dft_GGAx(7)%ref  / 'void' /
  DATA dft_GGAx(7)%wrn  / 'no functional with this index available' /
  ! PB0X
  DATA dft_GGAx(8)%ref  / 'J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)' /
  DATA dft_GGAx(8)%wrn  / 'none' /
  ! B3LP
  DATA dft_GGAx(9)%ref  / 'P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  DATA dft_GGAx(9)%wrn  / 'none' /
  ! PSX
  DATA dft_GGAx(10)%ref / 'J.P. Perdew et al., PRL 100, 136406 (2008)' /
  DATA dft_GGAx(10)%wrn / 'none' /
  ! WCX
  DATA dft_GGAx(11)%ref / 'Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)' /
  DATA dft_GGAx(11)%wrn / 'none' /
  ! HSE
  DATA dft_GGAx(12)%ref / 'Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003), &
                           &Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006)' /
  DATA dft_GGAx(12)%wrn / 'none' /
  ! RW86
  DATA dft_GGAx(13)%ref / 'Eamonn D. Murray et al, J. Chem. Theory Comput. 5, 2754 (2009)' /
  DATA dft_GGAx(13)%wrn / 'none' /
  ! PBE
  DATA dft_GGAx(14)%ref / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  DATA dft_GGAx(14)%wrn / 'none' /
  ! xxxx
  DATA dft_GGAx(15)%ref / 'void' /
  DATA dft_GGAx(15)%wrn / 'no functional available with this ID' /
  ! C09X
  DATA dft_GGAx(16)%ref / 'V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)' /
  DATA dft_GGAx(16)%wrn / 'none' /
  ! SOX
  DATA dft_GGAx(17)%ref / 'Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)' /
  DATA dft_GGAx(17)%wrn / 'none' /
  ! xxxx
  DATA dft_GGAx(18)%ref / 'void' /
  DATA dft_GGAx(18)%wrn / 'no functional available with this ID' /
  ! Q2DX
  DATA dft_GGAx(19)%ref / 'L. Chiodo et al., PRL 108, 126402 (2012)' /
  DATA dft_GGAx(19)%wrn / 'none' /
  ! GAUP
  DATA dft_GGAx(20)%ref / 'J.-W. Song, K. Yamashita, K. Hirao, JCP 135, 071103 (2011)' /
  DATA dft_GGAx(20)%wrn / 'none' /
  ! PW86
  DATA dft_GGAx(21)%ref / 'J.P.Perdew, PRB 33, 8800 (1986)' /
  DATA dft_GGAx(21)%wrn / 'none' /
  ! B86B
  DATA dft_GGAx(22)%ref / 'A.D.Becke, J.Chem.Phys. 85, 7184 (1986)' /
  DATA dft_GGAx(22)%wrn / 'none' /
  ! OBK8
  DATA dft_GGAx(23)%ref / 'Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)' /
  DATA dft_GGAx(23)%wrn / 'none' /
  ! OB86
  DATA dft_GGAx(24)%ref / 'Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)' /
  DATA dft_GGAx(24)%wrn / 'none' /
  ! EVX
  DATA dft_GGAx(25)%ref / 'Engel-Vosko, Phys. Rev. B 47, 13164 (1993)' /
  DATA dft_GGAx(25)%wrn / 'none' /
  ! B86R
  DATA dft_GGAx(26)%ref / 'I. Hamada, Phys. Rev. B 89, 121103(R) (2014)' /
  DATA dft_GGAx(26)%wrn / 'none' /
  ! CX13
  DATA dft_GGAx(27)%ref / 'K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)' /
  DATA dft_GGAx(27)%wrn / 'none' /
  ! X3LP
  DATA dft_GGAx(28)%ref / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  DATA dft_GGAx(28)%wrn / 'none' /
  ! CX0
  DATA dft_GGAx(29)%ref / 'K. Berland, Y. Jiao, J.-H. Lee, T. Rangel, J. B. Neaton &
                           &and P. Hyldgaard, J. Chem. Phys. 146, 234106 (2017)' /
  DATA dft_GGAx(29)%wrn / 'none' /
  ! rPW86+HF/4 (rw86-0) (for DF2-0)
  DATA dft_GGAx(30)%ref / 'K. Berland, Y. Jiao, J.-H. Lee, T. Rangel, J. B.  Neaton &
                           &and P. Hyldgaard, J. Chem. Phys. 146, 234106 (2017)' /
  DATA dft_GGAx(30)%wrn / 'none' /
  ! CX0P  vdW-DF-cx+HF/5 (cx13-0p)
  DATA dft_GGAx(31)%ref / 'Y. Jiao, E. Schröder and P. Hyldgaard, &
                           &J. Chem. Phys. 148, 194115 (2018)' /
  DATA dft_GGAx(31)%wrn / 'none' /
  ! AHCX (part of vdW-DF-ahcx)
  DATA dft_GGAx(32)%ref / 'V. Shukla, Y. Jiao, C.M. Frostenson and Per Hyldgaard &
                           &J. Phys.:Condens. Matter 34, 025902 (2022)' /
  DATA dft_GGAx(32)%wrn / 'none' /
  ! AHF2  (part of vdW-DF2-AH)
  DATA dft_GGAx(33)%ref / 'V. Shukla, Y. Jiao, C.M. Frostenson and Per Hyldgaard &
                           &J. Phys.:Condens. Matter 34, 025902 (2022)' /
  DATA dft_GGAx(33)%wrn / 'none' /
  ! AHPB  (part of PBE-AH)
  DATA dft_GGAx(34)%ref / 'J Chem. Phys. 128, 194105 (2008) + &
                           &J. Phys.:Condens. Matter 34, 025902 (2022); Compare HJS08-PBE' /
  DATA dft_GGAx(34)%wrn / 'none' /
  ! AHPS (part of PBESOL-AH)
  DATA dft_GGAx(35)%ref / 'J Chem. Phys. 128, 194105 (2008) + &
                           &J. Phys.:Condens. Matter 34, 025902 (2022); Compare HJS08-PBESOL' /
  DATA dft_GGAx(35)%wrn / 'none' /
  ! CX14  (reserved P.H.)
  DATA dft_GGAx(36)%ref / 'Reserved, no ref. available' /
  DATA dft_GGAx(36)%wrn / 'never called in QE' /
  ! CX15  (reserved P.H.)
  DATA dft_GGAx(37)%ref / 'Reserved, no ref. available' /
  DATA dft_GGAx(37)%wrn / 'never called in QE' /
  ! BR0 
  DATA dft_GGAx(38)%ref / 'vdW-DF2-b86r+HF/4: Phys. Rev. B 89, 121103(R) (2014)+&
                           &J. Chem. Phys. 146, 234106 (2017), PRX 12, 041003 (2022)' /
  DATA dft_GGAx(38)%wrn / 'none' /
  ! CX16  (reserved P.H.)
  DATA dft_GGAx(39)%ref / 'Reserved, no ref. available' /
  DATA dft_GGAx(39)%wrn / 'none' /
  ! C090
  DATA dft_GGAx(40)%ref / 'vdW-DF-c09+HF/4: Phys. Rev. B 81, 161104(R) (2010)+&
                           &J. Chem. Phys. 146, 234106 (2017)' /
  DATA dft_GGAx(40)%wrn / 'none' /
  ! B86X
  DATA dft_GGAx(41)%ref / '[B86B exchange * 0.75]' /
  DATA dft_GGAx(41)%wrn / 'none' /
  ! B88X
  DATA dft_GGAx(42)%ref / '[Becke88 exchange * 0.50]' /
  DATA dft_GGAx(42)%wrn / 'none' /
  ! BEEX
  DATA dft_GGAx(43)%ref / 'BEE exchange' /
  DATA dft_GGAx(43)%wrn / 'none' /
  ! HHNX
  DATA dft_GGAx(44)%ref / 'Hammer-Hansen-Norskov' /
  DATA dft_GGAx(44)%wrn / 'none' /
  ! W31X  vdW-DF3-opt1 exchange
  DATA dft_GGAx(45)%ref / 'D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)' /
  DATA dft_GGAx(45)%wrn / 'none' /
  ! W32X  vdW-DF3-opt2 exchange
  DATA dft_GGAx(46)%ref / 'D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)' /
  DATA dft_GGAx(46)%wrn / 'none' /
  ! AHRR  (reserved P.H., testing)
  DATA dft_GGAx(47)%ref / 'V. Shukla et al, PRX 12, 041003 (2022)' /
  ! EHPB  (reserved P.H.)
  DATA dft_GGAx(48)%ref / 'Reserved, No ref. available' /
  DATA dft_GGAx(48)%wrn / 'none' /
  ! HJPB  (Short-ranged PBE exchange by HJS08 parameters, cross-reference)
  DATA dft_GGAx(49)%ref / 'Short-ranged pbe exchange as set by Henderson et al, &
                           &J. Chem. Phys. 128, 194105 (2008); Compare PBE-AH' /
  DATA dft_GGAx(49)%wrn / 'none' /
  ! HJPS  (Short-ranged PBEsol exchange by HJS08 param, cross-reference)
  DATA dft_GGAx(50)%ref / 'Short-ranged pbesol exchange as set by Henderson et al, &
                           & J. Chem. Phys. 128, 194105 (2008); Compare PBESOL-AH' /
  DATA dft_GGAx(50)%wrn / 'none' /
  !
  !
  ! ---- GGA correlation ----
  ! NOGC
  DATA dft_GGAc(0)%ref  / 'No GGA correlation - default' /
  DATA dft_GGAc(0)%wrn  / 'none' /
  ! P86   Perdew86
  DATA dft_GGAc(1)%ref  / 'J.P.Perdew, PRB 33, 8822 (1986)' /
  DATA dft_GGAc(1)%wrn  / 'none' /
  ! GGC   Perdew-Wang 91 corr.
  DATA dft_GGAc(2)%ref  / 'J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)' /
  DATA dft_GGAc(2)%wrn  / 'none' /
  ! BLYP  Lee-Yang-Parr
  DATA dft_GGAc(3)%ref  / 'C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)' /
  DATA dft_GGAc(3)%wrn  / 'none' /
  ! PBC   Perdew-Burke-Ernzenhof corr.
  DATA dft_GGAc(4)%ref  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  DATA dft_GGAc(4)%wrn  / 'none' /
  ! HCTH  Cambridge corr, Handy et al.
  DATA dft_GGAc(5)%ref  / 'Handy et al, JCP 109, 6264 (1998)' /
  DATA dft_GGAc(5)%wrn  / 'none' /
  ! xxxx
  DATA dft_GGAc(6)%ref  / 'void' /
  DATA dft_GGAc(6)%wrn  / 'no GGAc functional available with this ID' /
  ! B3LP  b3lyp (Lee-Yang-Parr*0.81)
  DATA dft_GGAc(7)%ref  / 'P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  DATA dft_GGAc(7)%wrn  / 'none' /
  ! PSC   PBEsol corr
  DATA dft_GGAc(8)%ref  / 'J.P. Perdew et al., PRL 100, 136406 (2008)' /
  DATA dft_GGAc(8)%wrn  / 'none' /
  ! PBE   same as PBX, back-compatibility
  DATA dft_GGAc(9)%ref  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  DATA dft_GGAc(9)%wrn  / 'none' /
  ! void
  DATA dft_GGAc(10)%ref / 'void' /
  DATA dft_GGAc(10)%wrn / 'no GGAc functional available with this ID' /
  ! void
  DATA dft_GGAc(11)%ref / 'void' /
  DATA dft_GGAc(11)%wrn / 'no GGAc functional available with this ID' /
  ! Q2DC  correlation grad corr.
  DATA dft_GGAc(12)%ref / 'L. Chiodo et al., PRL 108, 126402 (2012)' /
  DATA dft_GGAc(12)%wrn / 'none' /
  ! X3LC  (Lee-Yang-Parr*0.871)
  DATA dft_GGAc(13)%ref / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  DATA dft_GGAc(13)%wrn / 'none' /
  ! BEEC  beef correlation
  DATA dft_GGAc(14)%ref / 'BEEF correlation' /
  DATA dft_GGAc(14)%wrn / 'none' /
  !
  !
  ! ---- MGGA (exchange+correlation) ----
  !
  ! NONE
  DATA dft_MGGA(0)%ref  / 'No mGGA exchange.' /
  DATA dft_MGGA(0)%wrn  / 'none' /
  ! TPSS
  DATA dft_MGGA(1)%ref  / 'J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, PRL 91, 146401 (2003)' /
  DATA dft_MGGA(1)%wrn  / 'none' /
  ! M06L
  DATA dft_MGGA(2)%ref  / 'Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)' /
  DATA dft_MGGA(2)%wrn  / 'none' /
  ! TB09
  DATA dft_MGGA(3)%ref  / 'F. Tran and P. Blaha, Phys.Rev.Lett. 102, 226401 (2009)' /
  DATA dft_MGGA(3)%wrn  / 'needs Libxc, provides only potential' /
  ! void
  DATA dft_MGGA(4)%ref  / 'void' /
  DATA dft_MGGA(4)%wrn  / 'no MGGA functional available with this ID' /
  ! SCAN
  DATA dft_MGGA(5)%ref  / 'J Sun, A Ruzsinszky and J Perdew, PRL 115, 36402 (2015)' /
  DATA dft_MGGA(5)%wrn  / 'needs Libxc' /
  ! SCA0
  DATA dft_MGGA(6)%ref  / 'K Hui and J-D. Chai, JCP 144, 44114 (2016)' /
  DATA dft_MGGA(6)%wrn  / 'needs Libxc' /
  ! R2SCAN
  DATA dft_MGGA(7)%ref  / 'J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, &
                          &and J. Sun, JPCL 11, 8208 (2020)' /
  DATA dft_MGGA(7)%wrn  / 'needs Libxc' /
  ! RSCAN
  DATA dft_MGGA(8)%ref  / 'A. P. Bartók, J. R. Yates., JCP 150, 161101 (2019)' /
  DATA dft_MGGA(8)%wrn  / 'needs Libxc' /
  !
  !
  ! ---- Full DFTs ----
  !
  ! PZ / LDA
  DATA dft_full_descr(1)  / 'Perdew-Zunger LDA' /
  ! PW
  DATA dft_full_descr(2)  / 'LDA with PW correlation' /
  ! VWN-RPA
  DATA dft_full_descr(3)  / 'VWN LDA using vwn1-rpa parametrization' /
  ! OEP
  DATA dft_full_descr(4)  / 'Optimized Effective Potential. No GC part, no corr. by default' /
  ! KLI
  DATA dft_full_descr(5)  / 'KLI - currently not implemented' /
  ! HF
  DATA dft_full_descr(6)  / 'HF no GC part (nor LDA...) and no correlation by default' /
  ! PBE
  DATA dft_full_descr(7)  / 'Perdew-Burke-Ernzerhof GGA' /
  ! B88
  DATA dft_full_descr(8)  / 'Becke88 (beta=0.0042)' /
  ! BP
  DATA dft_full_descr(9)  / 'Becke-Perdew grad.corr.' /
  ! PW91
  DATA dft_full_descr(10) / 'PW91 (aka GGA)' /
  ! REVPBE
  DATA dft_full_descr(11) / 'revPBE (Zhang-Yang)' /
  ! PBESOL
  DATA dft_full_descr(12) / 'PBEsol' /
  ! BLYP
  DATA dft_full_descr(13) / 'Becke-Lee-Yang-Parr LDA+GGA' /
  ! OPTBK88
  DATA dft_full_descr(14) / 'optB88' /
  ! OPTB86B
  DATA dft_full_descr(15) / 'optB86' /
  ! PBC
  DATA dft_full_descr(16) / 'PBC = PW + PBC' /
  ! HCTH
  DATA dft_full_descr(17) / 'HCTH/120' /
  ! OLYP
  DATA dft_full_descr(18) / 'OLYP = OPTX + LYP' /
  ! WC
  DATA dft_full_descr(19) / 'Wu-Cohen' /
  ! PW86PBE
  DATA dft_full_descr(20) / 'PW86 exchange + PBE correlation' /
  ! B86BPBE
  DATA dft_full_descr(21) / 'B86b exchange + PBE correlation' /
  ! PBEQ2D
  DATA dft_full_descr(22) / 'PBEQ2D' /
  ! SOGGA
  DATA dft_full_descr(23) / 'SOGGA' /
  ! EV93
  DATA dft_full_descr(24) / 'Engel-Vosko' /
  ! RPBE
  DATA dft_full_descr(25) / 'RPBE' /
  ! PBE0
  DATA dft_full_descr(26) / 'PBE0 in: Perdew, Ernzerhof, Burke, JCP 105, 9982 (1996)' /
  ! B86BPBEX
  DATA dft_full_descr(27) / 'B86bPBE hybrid' /
  ! BHAHLYP
  DATA dft_full_descr(28) / 'Becke half-and-half LYP' /
  ! HSE
  DATA dft_full_descr(29) / 'Heyd-Scuseria-Ernzerhof (HSE 06, see references)' /
  ! GAUP / GAUPBE
  DATA dft_full_descr(30) / 'Gau-PBE (also gaup)' /
  ! B3LYP
  DATA dft_full_descr(31) / 'B3LYP' /
  ! B3LYP-V1R
  DATA dft_full_descr(32) / 'B3LYP-VWN1-RPA' /
  ! X3LYP
  DATA dft_full_descr(33) / 'X3LYP' /
  ! TPSS
  DATA dft_full_descr(34) / 'TPSS Meta-GGA' /
  ! M06L
  DATA dft_full_descr(35) / 'M06L Meta-GGA' /
  ! TB09
  DATA dft_full_descr(36) / 'TB09 Meta-GGA - needs Libxc' /
  ! SCAN
  DATA dft_full_descr(37) / 'SCAN Meta-GGA - needs Libxc.' /
  ! SCAN0
  DATA dft_full_descr(38) / 'SCAN0 Meta-GGA - needs Libxc.' /
  ! R2SCAN
  DATA dft_full_descr(39) / 'R2SCAN Meta-GGA - needs Libxc.' /
  ! PBE-AH
  DATA dft_full_descr(40) / 'PBE-AH: HJS implementation (PBE params).' /
  ! PBESOL-AH
  DATA dft_full_descr(41) / 'PBESOL-AH: HJS implementation (PBEsol params).' /
  ! RSCAN
  DATA dft_full_descr(42) / 'RSCAN Meta-GGA - needs Libxc.' /
  !
END MODULE qe_dft_refs

