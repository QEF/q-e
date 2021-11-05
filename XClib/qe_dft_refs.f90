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
  !
  USE qe_dft_list, ONLY: nxc, ncc, ngcx, ngcc, nmeta, n_dft
  !
  SAVE
  !
  ! -- single DFT terms (family-type)
  CHARACTER(LEN=150) :: dft_LDAx_ref(0:nxc),  dft_LDAc_ref(0:ncc),  &
                        dft_GGAx_ref(0:ngcx), dft_GGAc_ref(0:ngcc), &
                        dft_MGGA_ref(0:nmeta)
  ! -- total DFTs
  CHARACTER(LEN=100) :: dft_full_descr(n_dft)
  !
  ! ---- LDA exchange ----
  !
  ! NOX
  DATA dft_LDAx_ref(0)  / 'No LDA exchange.' /
  ! SLA
  DATA dft_LDAx_ref(1)  / '[Slater exchange - alpha=2/3]' /
  ! SL1
  DATA dft_LDAx_ref(2)  / '[Slater exchange - alpha=1.0]' /
  ! RXC
  DATA dft_LDAx_ref(3)  / '[Relativistic Slater]' /
  ! OEP
  DATA dft_LDAx_ref(4)  / '[Optimized Effective Potential]' /
  ! HF
  DATA dft_LDAx_ref(5)  / '[Hartree-Fock]' /
  ! PB0X  (Slater*0.75+HF*0.25) for PBE0 and vdW-DF-cx0 and vdW-DF2-0 etc
  DATA dft_LDAx_ref(6)  / 'J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)' /
  ! B3LP  (Slater*0.80+HF*0.20)
  DATA dft_LDAx_ref(7)  / 'P.J.Stephens, F.J.Devlin, C.F.Chabalowski, M.J.Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  ! KZK   Finite-size corrections
  DATA dft_LDAx_ref(8)  / 'H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)' /
  ! xxxx [X3LYP_LDA]
  DATA dft_LDAx_ref(9)  / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  ! xxxx [KLI]
  DATA dft_LDAx_ref(10) / 'KLI aproximation for exx - currently not implemented' /
  !
  !
  !  ---- LDA correlation ----
  !
  ! NOC
  DATA dft_LDAc_ref(0)  / 'No LDA correlation.' /
  ! PZ
  DATA dft_LDAc_ref(1)  / 'J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)' /
  ! VWN
  DATA dft_LDAc_ref(2)  / 'S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)' /
  ! LYP
  DATA dft_LDAc_ref(3)  / 'C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)' /
  ! PW
  DATA dft_LDAc_ref(4)  / 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)' /
  ! WIG
  DATA dft_LDAc_ref(5)  / 'E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)' /
  ! HL
  DATA dft_LDAc_ref(6)  / 'L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)' /
  ! OBZ
  DATA dft_LDAc_ref(7)  / 'G.Ortiz and P.Ballone, PRB 50, 1391 (1994)' /
  ! OBW
  DATA dft_LDAc_ref(8)  / 'G.Ortiz and P.Ballone, PRB 50, 1391 (1994)' /
  ! GL
  DATA dft_LDAc_ref(9)  / 'O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)' /
  ! KZK
  DATA dft_LDAc_ref(10) / 'H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)' /
  ! xxxx [vwn1_rpa]
  DATA dft_LDAc_ref(11) / 'vwn1_rpa' /
  ! B3LP
  DATA dft_LDAc_ref(12) / 'P.J.Stephens, F.J.Devlin, C.F.Chabalowski, M.J.Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  ! xxxx [B3LYP-V1R]
  DATA dft_LDAc_ref(13) / 'B3LYP-V1R' /
  ! xxxx [X3LYP]
  DATA dft_LDAc_ref(14) / 'X3LYP' /
  !
  !
  !  ---- GGA exchange ----
  !
  ! NOGX
  DATA dft_GGAx_ref(0)  / 'No GGA exchange.' /
  ! B88
  DATA dft_GGAx_ref(1)  / 'A.D.Becke, PRA 38, 3098 (1988)' /
  ! GGX
  DATA dft_GGAx_ref(2)  / 'J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)' /
  ! PBX
  DATA dft_GGAx_ref(3)  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  ! REVX
  DATA dft_GGAx_ref(4)  / 'Zhang and Yang, PRL 80, 890 (1998)' /
  ! HCTH
  DATA dft_GGAx_ref(5)  / 'Handy et al, JCP 109, 6264 (1998)' /
  ! OPTX
  DATA dft_GGAx_ref(6)  / 'Handy et al, JCP 116, 5411 (2002)' /
  ! void
  DATA dft_GGAx_ref(7)  / 'void' /
  ! PB0X
  DATA dft_GGAx_ref(8)  / 'J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)' /
  ! B3LP
  DATA dft_GGAx_ref(9)  / 'P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  ! PSX
  DATA dft_GGAx_ref(10) / 'J.P. Perdew et al., PRL 100, 136406 (2008)' /
  ! WCX
  DATA dft_GGAx_ref(11) / 'Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)' /
  ! HSE
  DATA dft_GGAx_ref(12) / 'Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003), &
                           &Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006)' /
  ! RW86
  DATA dft_GGAx_ref(13) / 'Eamonn D. Murray et al, J. Chem. Theory Comput. 5, 2754 (2009)' /
  ! PBE
  DATA dft_GGAx_ref(14) / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  ! xxxx
  DATA dft_GGAx_ref(15) / 'void' / 
  ! C09X
  DATA dft_GGAx_ref(16) / 'V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)' /
  ! SOX
  DATA dft_GGAx_ref(17) / 'Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)' /
  ! xxxx
  DATA dft_GGAx_ref(18) / 'void' /
  ! Q2DX
  DATA dft_GGAx_ref(19) / 'L. Chiodo et al., PRL 108, 126402 (2012)' /
  ! GAUP
  DATA dft_GGAx_ref(20) / 'J.-W. Song, K. Yamashita, K. Hirao, JCP 135, 071103 (2011)' /
  ! PW86
  DATA dft_GGAx_ref(21) / 'J.P.Perdew, PRB 33, 8800 (1986)' /
  ! B86B
  DATA dft_GGAx_ref(22) / 'A.D.Becke, J.Chem.Phys. 85, 7184 (1986)' /
  ! OBK8
  DATA dft_GGAx_ref(23) / 'Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)' /
  ! OB86
  DATA dft_GGAx_ref(24) / 'Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)' /
  ! EVX
  DATA dft_GGAx_ref(25) / 'Engel-Vosko, Phys. Rev. B 47, 13164 (1993)' /
  ! B86R
  DATA dft_GGAx_ref(26) / 'I. Hamada, Phys. Rev. B 89, 121103(R) (2014)' /
  ! CX13
  DATA dft_GGAx_ref(27) / 'K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)' /
  ! X3LP
  DATA dft_GGAx_ref(28) / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  ! CX0
  DATA dft_GGAx_ref(29) / 'K. Berland, Y. Jiao, J.-H. Lee, T. Rangel, J. B. Neaton &
                           &and P. Hyldgaard, J. Chem. Phys. 146, 234106 (2017)' /
  ! R860
  DATA dft_GGAx_ref(30) / 'rPW86+HF/4 (rw86-0) (for DF0) - no ref. available' /
  ! CX0P  vdW-DF-cx+HF/5 (cx13-0p)
  DATA dft_GGAx_ref(31) / 'Y. Jiao, E. Schröder and P. Hyldgaard, &
                           &J. Chem. Phys. 148, 194115 (2018)' /
  ! AHCX  (reserved PH)
  DATA dft_GGAx_ref(32) / 'vdW-DF-cx based not yet in use' /
  ! AHF2  (reserved PH)
  DATA dft_GGAx_ref(33) / 'vdW-DF2 based not yet in use' /
  ! AHPB  (reserved PH)
  DATA dft_GGAx_ref(34) / 'PBE based not yet in use' /
  ! AHPS
  DATA dft_GGAx_ref(35) / 'PBE-sol based not yet in use' /
  ! CX14  (reserved PH)
  DATA dft_GGAx_ref(36) / 'no ref. available' /
  ! CX15  (reserved PH)
  DATA dft_GGAx_ref(37) / 'no ref. available' /
  ! BR0
  DATA dft_GGAx_ref(38) / 'vdW-DF2-b86r+HF/4 (b86r-0) - no ref. available' /
  ! CX16  (reserved PH)
  DATA dft_GGAx_ref(39) / 'no ref. available' /
  ! C090
  DATA dft_GGAx_ref(40) / 'vdW-DF-c09+HF/4 (c09-0) - no ref. available' /
  ! B86X
  DATA dft_GGAx_ref(41) / '[B86B exchange * 0.75]' /
  ! B88X
  DATA dft_GGAx_ref(42) / '[Becke88 exchange * 0.50]' /
  ! BEEX
  DATA dft_GGAx_ref(43) / 'BEE exchange' /
  ! HHNX
  DATA dft_GGAx_ref(44) / 'Hammer-Hansen-Norskov' /
  ! W31X  vdW-DF3-opt1 exchange
  DATA dft_GGAx_ref(45) / 'D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)' /
  ! W32X  vdW-DF3-opt2 exchange
  DATA dft_GGAx_ref(46) / 'D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)' /
  !
  !
  ! ---- GGA correlation ----
  ! NOGC
  DATA dft_GGAc_ref(0)  / 'No GGA correlation - default' /
  ! P86   Perdew86
  DATA dft_GGAc_ref(1)  / 'J.P.Perdew, PRB 33, 8822 (1986)' /
  ! GGC   Perdew-Wang 91 corr.
  DATA dft_GGAc_ref(2)  / 'J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)' /
  ! BLYP  Lee-Yang-Parr
  DATA dft_GGAc_ref(3)  / 'C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)' /
  ! PBC   Perdew-Burke-Ernzenhof corr.
  DATA dft_GGAc_ref(4)  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  ! HCTH  Cambridge corr, Handy et al.
  DATA dft_GGAc_ref(5)  / 'Handy et al, JCP 109, 6264 (1998)' /
  ! xxxx
  DATA dft_GGAc_ref(6)  / 'void' /
  ! B3LP  b3lyp (Lee-Yang-Parr*0.81)
  DATA dft_GGAc_ref(7)  / 'P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  ! PSC   PBEsol corr
  DATA dft_GGAc_ref(8)  / 'J.P. Perdew et al., PRL 100, 136406 (2008)' /
  ! PBE   same as PBX, back-compatibility
  DATA dft_GGAc_ref(9)  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  ! void
  DATA dft_GGAc_ref(10) / 'void' /
  ! void
  DATA dft_GGAc_ref(11) / 'void' /
  ! Q2DC  correlation grad corr.
  DATA dft_GGAc_ref(12) / 'L. Chiodo et al., PRL 108, 126402 (2012)' /
  ! X3LC  (Lee-Yang-Parr*0.871)
  DATA dft_GGAc_ref(13) / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  ! BEEC  beef correlation
  DATA dft_GGAc_ref(14) / 'BEEF correlation' /
  !
  !
  ! ---- MGGA (exchange+correlation) ----
  !
  ! NONE
  DATA dft_MGGA_ref(0)  / 'No mGGA exchange.' /
  ! TPSS
  DATA dft_MGGA_ref(1)  / 'J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, PRL 91, 146401 (2003)' /
  ! M06L
  DATA dft_MGGA_ref(2)  / 'Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)' /
  ! TB09
  DATA dft_MGGA_ref(3)  / 'F. Tran and P. Blaha, Phys.Rev.Lett. 102, 226401 (2009) - Libxc needed' /
  ! void
  DATA dft_MGGA_ref(4)  / 'void' /
  ! SCAN
  DATA dft_MGGA_ref(5)  / 'J Sun, A Ruzsinszky and J Perdew, PRL 115, 36402 (2015) - Libxc needed' /
  ! SCA0
  DATA dft_MGGA_ref(6)  / 'K Hui and J-D. Chai, JCP 144, 44114 (2016)' /
  ! R2SCAN
  DATA dft_MGGA_ref(7)  / 'J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, &
                          &and J. Sun, JPCL 11, 8208 (2020) - Libxc needed' /
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
  !
END MODULE qe_dft_refs

