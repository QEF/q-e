!
! Copyright (C) 2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------------
MODULE qe_dft_list
  !-----------------------------------------------------------------------------------
  !! Contains the list of functionals available in QE, both single terms (family+type)
  !! and combinations.
  !
  USE dft_setting_params,  ONLY: notset
  
  SAVE
  !
  ! -- single DFT terms (family-type)
  INTEGER, PARAMETER :: nxc=8, ncc=10, ngcx=43, ngcc=10, nmeta=6
  CHARACTER(LEN=4)  :: dft_LDAx_name(0:nxc),  dft_LDAc_name(0:ncc),  &
                       dft_GGAx_name(0:ngcx), dft_GGAc_name(0:ngcc), &
                       dft_MGGA_name(0:nmeta) 
  CHARACTER(LEN=150) :: dft_LDAx_ref(0:nxc),  dft_LDAc_ref(0:ncc),  &
                       dft_GGAx_ref(0:ngcx), dft_GGAc_ref(0:ngcc), &
                       dft_MGGA_ref(0:nmeta)
  !
  ! -- total DFTs
  INTEGER, PARAMETER :: n_dft = 41
  CHARACTER(LEN=10) :: dft_name(n_dft)
  CHARACTER(LEN=10) :: dft_name2(n_dft)
  CHARACTER(LEN=80) :: dft_descr(n_dft)
  INTEGER, DIMENSION(n_dft,6) :: dft_IDs
  !
  !
  ! ---- LDA exchange ----
  !
  !
  DATA dft_LDAx_name(0)  / 'NOX' /
  DATA dft_LDAx_ref(0)   / 'No LDA exchange.' /
  !
  DATA dft_LDAx_name(1)  / 'SLA' /
  DATA dft_LDAx_ref(1)   / '[Slater exchange - alpha=2/3]' /
  !
  DATA dft_LDAx_name(2)  / 'SL1' /
  DATA dft_LDAx_ref(2)   / '[Slater exchange - alpha=1.0]' /
  !
  DATA dft_LDAx_name(3)  / 'RXC' /
  DATA dft_LDAx_ref(3)   / '[Relativistic Slater]' /
  !
  DATA dft_LDAx_name(4)  / 'OEP' /
  DATA dft_LDAx_ref(4)   / '[Optimized Effective Potential]' /
  !
  DATA dft_LDAx_name(5)  / 'HF' /
  DATA dft_LDAx_ref(5)   / '[Hartree-Fock]' /
  !
  DATA dft_LDAx_name(6)  / 'PB0X' /
  DATA dft_LDAx_ref(6)   / 'J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)' /
  !
  DATA dft_LDAx_name(7)  / 'B3LP' /
  DATA dft_LDAx_ref(7)   / 'P.J.Stephens, F.J.Devlin, C.F.Chabalowski, M.J.Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  !
  DATA dft_LDAx_name(8)  / 'KZK' /
  DATA dft_LDAx_ref(8)   / 'H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)' /
  !
  !
  !  ---- LDA correlation ----
  !
  !
  DATA dft_LDAc_name(0)  / 'NOX' /
  DATA dft_LDAc_ref(0)   / 'No LDA correlation.' /
  !
  DATA dft_LDAc_name(1)  / 'PZ' /
  DATA dft_LDAc_ref(1)   / 'J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)' /
  !
  DATA dft_LDAc_name(2)  / 'VWN' /
  DATA dft_LDAc_ref(2)   / 'S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)' /
  !
  DATA dft_LDAc_name(3)  / 'LYP' /
  DATA dft_LDAc_ref(3)   / 'C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)' /
  !
  DATA dft_LDAc_name(4)  / 'PW' /
  DATA dft_LDAc_ref(4)   / 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)' /
  !
  DATA dft_LDAc_name(5)  / 'WIG' /
  DATA dft_LDAc_ref(5)   / 'E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)' /
  !
  DATA dft_LDAc_name(6)  / 'HL' /
  DATA dft_LDAc_ref(6)   / 'L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)' /
  !
  DATA dft_LDAc_name(7)  / 'OBZ' /
  DATA dft_LDAc_ref(7)   / 'G.Ortiz and P.Ballone, PRB 50, 1391 (1994)' /
  !
  DATA dft_LDAc_name(8)  / 'OBW' /
  DATA dft_LDAc_ref(8)   / 'G.Ortiz and P.Ballone, PRB 50, 1391 (1994)' /
  !
  DATA dft_LDAc_name(9)  / 'GL' /
  DATA dft_LDAc_ref(9)   / 'O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)' /
  !
  DATA dft_LDAc_name(10) / 'KZK' /
  DATA dft_LDAc_ref(10)  / 'H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)' /
  !
  !
  !  ---- GGA exchange ----
  !
  !
  DATA dft_GGAx_name(0)  / 'NOGX' /
  DATA dft_GGAx_ref(0)   / 'No GGA exchange.' /
  !
  DATA dft_GGAx_name(1)  / 'B88' /
  DATA dft_GGAx_ref(1)   / 'A.D.Becke, PRA 38, 3098 (1988)' /
  !
  DATA dft_GGAx_name(2)  / 'GGX' /
  DATA dft_GGAx_ref(2)   / 'J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)' /
  !
  DATA dft_GGAx_name(3)  / 'PBX' /
  DATA dft_GGAx_ref(3)   / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  !
  DATA dft_GGAx_name(4)  / 'RPB' /
  DATA dft_GGAx_ref(4)   / 'Zhang and Yang, PRL 80, 890 (1998)' /
  !
  DATA dft_GGAx_name(5)  / 'HCTH' /
  DATA dft_GGAx_ref(5)   / 'Handy et al, JCP 109, 6264 (1998)' /
  !
  DATA dft_GGAx_name(6)  / 'OPTX' /
  DATA dft_GGAx_ref(6)   / 'Handy et al, JCP 116, 5411 (2002)' /
  !
  DATA dft_GGAx_name(7)  / 'PB0X' /
  DATA dft_GGAx_ref(7)   / 'J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)' /
  !
  DATA dft_GGAx_name(8)  / 'B3LP' /
  DATA dft_GGAx_ref(8)   / 'P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  !
  DATA dft_GGAx_name(9)  / 'PSX' /
  DATA dft_GGAx_ref(9)   / 'J.P. Perdew et al., PRL 100, 136406 (2008)' /
  !
  DATA dft_GGAx_name(10) / 'WCX' /
  DATA dft_GGAx_ref(10)  / 'Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)' /
  !
  DATA dft_GGAx_name(11) / 'HSE' /
  DATA dft_GGAx_ref(11)  / 'Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003), &
                           &Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006)' /
  !
  DATA dft_GGAx_name(12) / 'RW86' /
  DATA dft_GGAx_ref(12)  / 'Eamonn D. Murray et al, J. Chem. Theory Comput. 5, 2754 (2009)' /
  !
  DATA dft_GGAx_name(13) / 'PBE' /
  DATA dft_GGAx_ref(13)  / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  !
  DATA dft_GGAx_name(14) / 'C09X' /
  DATA dft_GGAx_ref(14)  / 'V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)' /
  !
  DATA dft_GGAx_name(15) / 'SOX' /
  DATA dft_GGAx_ref(15)  / 'Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)' /
  !
  DATA dft_GGAx_name(16) / 'Q2DX' /
  DATA dft_GGAx_ref(16)  / 'L. Chiodo et al., PRL 108, 126402 (2012)' /
  !
  DATA dft_GGAx_name(17) / 'GAUP' /
  DATA dft_GGAx_ref(17)  / 'J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)' /
  !
  DATA dft_GGAx_name(18)  / 'PW86' /
  DATA dft_GGAx_ref(18) / 'J.P.Perdew, PRB 33, 8800 (1986)' /
  !
  DATA dft_GGAx_name(19) / 'B86B' /
  DATA dft_GGAx_ref(19)  / 'A.D.Becke, J.Chem.Phys. 85, 7184 (1986)' /
  !
  DATA dft_GGAx_name(20) / 'OBK8' /
  DATA dft_GGAx_ref(20)  / 'Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)' /
  !
  DATA dft_GGAx_name(21) / 'OB86' /
  DATA dft_GGAx_ref(21)  / 'Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)' /
  !
  DATA dft_GGAx_name(22) / 'EVX' /
  DATA dft_GGAx_ref(22)  / 'Engel-Vosko, Phys. Rev. B 47, 13164 (1993)' /
  !
  DATA dft_GGAx_name(23) / 'B86R' /
  DATA dft_GGAx_ref(23)  / 'I. Hamada, Phys. Rev. B 89, 121103(R) (2014)' /
  !
  DATA dft_GGAx_name(24) / 'CX13' /
  DATA dft_GGAx_ref(24)  / 'K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)' /
  !
  DATA dft_GGAx_name(25) / 'X3LP' /
  DATA dft_GGAx_ref(25)  / 'X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)' /
  !
  DATA dft_GGAx_name(26) / 'CX0' /
  DATA dft_GGAx_ref(26)  / 'K. Berland, Y. Jiao, J.-H. Lee, T. Rangel, J. B. Neaton &
                           &and P. Hyldgaard, J. Chem. Phys. 146, 234106 (2017)' /
  !
  DATA dft_GGAx_name(27) / 'R860' /
  DATA dft_GGAx_ref(27)  / '' /
  !
  DATA dft_GGAx_name(28) / 'CX0P' /
  DATA dft_GGAx_ref(28)  / 'Y. Jiao, E. Schröder and P. Hyldgaard, &
                           &J. Chem. Phys. 148, 194115 (2018)' /
  !
  DATA dft_GGAx_name(29) / 'AHCX' /
  DATA dft_GGAx_ref(29)  / '' /
  !
  DATA dft_GGAx_name(30) / 'AHF2' /
  DATA dft_GGAx_ref(30)  / '' /
  !
  DATA dft_GGAx_name(31) / 'AHPB' /
  DATA dft_GGAx_ref(31)  / '' /
  !
  DATA dft_GGAx_name(32) / 'AHPS' /
  DATA dft_GGAx_ref(32)  / '' /
  !
  DATA dft_GGAx_name(33) / 'CX14' /
  DATA dft_GGAx_ref(33)  / '' /
  !
  DATA dft_GGAx_name(34) / 'CX15' /
  DATA dft_GGAx_ref(34)  / '' /
  !
  DATA dft_GGAx_name(35) / 'BR0' /
  DATA dft_GGAx_ref(35)  / '' /
  !
  DATA dft_GGAx_name(36) / 'CX16' /
  DATA dft_GGAx_ref(36)  / '' /
  !
  DATA dft_GGAx_name(37) / 'C090' /
  DATA dft_GGAx_ref(37)  / '' /
  !
  DATA dft_GGAx_name(38) / 'B86X' /
  DATA dft_GGAx_ref(38)  / '[B86B exchange * 0.75]' /
  !
  DATA dft_GGAx_name(39) / 'B88X' /
  DATA dft_GGAx_ref(39)  / '[Becke88 exchange * 0.50]' /
  !
  DATA dft_GGAx_name(40) / 'BEEX' /
  DATA dft_GGAx_ref(40)  / '' /
  !
  DATA dft_GGAx_name(41) / 'HHNX' /
  DATA dft_GGAx_ref(41)  / 'Hammer-Hansen-Norskov' /
  !
  DATA dft_GGAx_name(42) / 'W31X' /
  DATA dft_GGAx_ref(42)  / 'D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)' /
  !
  DATA dft_GGAx_name(43) / 'W32X' /
  DATA dft_GGAx_ref(43)  / 'D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)' /
  !
  !
  ! ---- GGA correlation ----
  !
  !
  DATA dft_GGAc_name(0)  / 'NOGC' /
  DATA dft_GGAc_ref(0)   / 'No GGA correlation.' /
  !
  DATA dft_GGAc_name(1)  / 'P86' /
  DATA dft_GGAc_ref(1)   / 'J.P.Perdew, PRB 33, 8822 (1986)' /
  !
  DATA dft_GGAc_name(2)  / 'GGC' /
  DATA dft_GGAc_ref(2)   / 'J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)' /
  !
  DATA dft_GGAc_name(3)  / 'BLYP' /
  DATA dft_GGAc_ref(3)   / 'C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)' /
  !
  DATA dft_GGAc_name(4)  / 'PBC' /
  DATA dft_GGAc_ref(4)   / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  !
  DATA dft_GGAc_name(5)  / 'HCTH' /
  DATA dft_GGAc_ref(5)   / 'Handy et al, JCP 109, 6264 (1998)' /
  !
  DATA dft_GGAc_name(6)  / 'B3LP' /
  DATA dft_GGAc_ref(6)   / 'P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch, &
                           &J.Phys.Chem 98, 11623 (1994)' /
  !
  DATA dft_GGAc_name(7)  / 'PSC' /
  DATA dft_GGAc_ref(7)   / 'J.P. Perdew et al., PRL 100, 136406 (2008)' /
  !
  DATA dft_GGAc_name(8)  / 'PBE' /
  DATA dft_GGAc_ref(8)   / 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' /
  !
  DATA dft_GGAc_name(9)  / 'Q2DC' /
  DATA dft_GGAc_ref(9)   / 'L. Chiodo et al., PRL 108, 126402 (2012)' /
  !
  DATA dft_GGAc_name(10) / 'BEEC' /
  DATA dft_GGAc_ref(10)  / '' /
  !
  !
  ! ---- MGGA (exchange+correlation) ----
  !
  !
  DATA dft_MGGA_name(0)  / 'NONE' /
  DATA dft_MGGA_ref(0)   / 'No mGGA exchange.' /
  !
  DATA dft_MGGA_name(1)  / 'TPSS' /
  DATA dft_MGGA_ref(1)   / 'J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, PRL 91, 146401 (2003)' /
  !
  DATA dft_MGGA_name(2)  / 'M06L' /
  DATA dft_MGGA_ref(2)   / 'Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)' /
  !
  DATA dft_MGGA_name(3)  / 'TB09' /
  DATA dft_MGGA_ref(3)   / 'F Tran and P Blaha, Phys.Rev.Lett. 102, 226401 (2009)' /
  !
  DATA dft_MGGA_name(4)  / 'META' /
  DATA dft_MGGA_ref(4)   / '' /
  !
  DATA dft_MGGA_name(5)  / 'SCAN' /
  DATA dft_MGGA_ref(5)   / 'J Sun, A Ruzsinszky and J Perdew, PRL 115, 36402 (2015)' /
  !
  DATA dft_MGGA_name(6)  / 'SCA0' /
  DATA dft_MGGA_ref(6)   / 'K Hui and J-D. Chai, JCP 144, 44114 (2016)' /
  !
  !
  ! ---- Full DFTs ----
  !
  !
  DATA dft_name(1)     / 'PZ' /
  DATA dft_name2(1)    / 'LDA' /
  DATA dft_IDs(1,1:6)  / 1,1,0,0,0,0 /
  DATA dft_descr(1)    / 'Perdew-Zunger LDA' /
  !
  DATA dft_name(2)     / 'PW' /
  DATA dft_name2(2)    / 'none' /
  DATA dft_IDs(2,1:6)  / 1,4,0,0,0,0 /
  DATA dft_descr(2)    / 'LDA with PW correlation' /
  !
  DATA dft_name(3)     / 'VWN-RPA' /
  DATA dft_name2(3)    / 'none' /
  DATA dft_IDs(3,1:6)  / 1,11,0,0,0,0 /
  DATA dft_descr(3)    / 'VWN LDA using vwn1-rpa parametrization' /
  !
  DATA dft_name(4)     / 'OEP' /
  DATA dft_name2(4)    / 'none' /
  DATA dft_IDs(4,1:6)  / 4,0,0,0,0,0 /
  DATA dft_descr(4)    / 'Optimized Effective Potential. No GC part, no corr. by default' /
  !
  DATA dft_name(5)     / 'KLI' /
  DATA dft_name2(5)    / 'none' /
  DATA dft_IDs(5,1:6)  / 10,0,0,0,0,0 /
  DATA dft_descr(5)    / 'KLI - unavailable' /
  !
  DATA dft_name(6)     / 'HF' /
  DATA dft_name2(6)    / 'none' /
  DATA dft_IDs(6,1:6)  / 5,0,0,0,0,0 /
  DATA dft_descr(6)    / 'HF no GC part (nor LDA...) and no correlation by default' /
  !
  DATA dft_name(7)     / 'PBE' /
  DATA dft_name2(7)    / 'none' /
  DATA dft_IDs(7,1:6)  / 1,4,3,4,0,0 /
  DATA dft_descr(7)    / 'Perdew-Burke-Ernzerhof GGA' /
  !
  DATA dft_name(8)     / 'B88' /
  DATA dft_name2(8)    / 'none' /
  DATA dft_IDs(8,1:6)  / 1,1,1,0,0,0 /
  DATA dft_descr(8)    / 'B88' /
  !
  DATA dft_name(9)     / 'BP' /
  DATA dft_name2(9)    / 'none' /
  DATA dft_IDs(9,1:6)  / 1,1,1,1,0,0 /
  DATA dft_descr(9)    / 'Becke-Perdew grad.corr.' /
  !----
  DATA dft_name(10)    / 'PW91' /
  DATA dft_name2(10)   / 'none' /
  DATA dft_IDs(10,1:6) / 1,4,2,2,0,0 /
  DATA dft_descr(10)   / 'PW91 (aka GGA)' /
  !
  DATA dft_name(11)    / 'REVPBE' /
  DATA dft_name2(11)   / 'none' /
  DATA dft_IDs(11,1:6) / 1,4,4,4,0,0 /
  DATA dft_descr(11)   / 'revPBE (Zhang-Yang)' /
  !
  DATA dft_name(12)    / 'PBESOL' /
  DATA dft_name2(12)   / 'none' /
  DATA dft_IDs(12,1:6) / 1,4,9,7,0,0 /
  DATA dft_descr(12)   / 'PBEsol' /
  !
  DATA dft_name(13)    / 'BLYP' /
  DATA dft_name2(13)   / 'none' /
  DATA dft_IDs(13,1:6) / 1,3,1,3,0,0 /
  DATA dft_descr(13)   / 'Becke-Lee-Yang-Parr LDA+GGA' /
  !
  DATA dft_name(14)    / 'OPTBK88' /
  DATA dft_name2(14)   / 'none' /
  DATA dft_IDs(14,1:6) / 1,4,20,1,0,0 /
  DATA dft_descr(14)   / 'optB88' /
  !
  DATA dft_name(15)    / 'OPTB86B' /
  DATA dft_name2(15)   / 'none' /
  DATA dft_IDs(15,1:6) / 1,4,21,1,0,0 /
  DATA dft_descr(15)   / 'optB86' /
  !
  DATA dft_name(16)    / 'PBC' /
  DATA dft_name2(16)   / 'none' /
  DATA dft_IDs(16,1:6) / 1,4,0,4,0,0 /
  DATA dft_descr(16)   / 'PBC = PW + PBC' /
  !
  DATA dft_name(17)    / 'HCTH' /
  DATA dft_name2(17)   / 'none' /
  DATA dft_IDs(17,1:6) / 0,0,5,5,0,0 /
  DATA dft_descr(17)   / 'HCTH/120' /
  !
  DATA dft_name(18)    / 'OLYP' /
  DATA dft_name2(18)   / 'none' /
  DATA dft_IDs(18,1:6) / 0,3,6,3,0,0 /
  DATA dft_descr(18)   / 'OLYP = OPTX + LYP' /
  !
  DATA dft_name(19)    / 'WC' /
  DATA dft_name2(19)   / 'none' /
  DATA dft_IDs(19,1:6) / 1,4,10,4,0,0 /
  DATA dft_descr(19)   / 'Wu-Cohen' /
  !
  DATA dft_name(20)    / 'PW86PBE' /
  DATA dft_name2(20)   / 'none' /
  DATA dft_IDs(20,1:6) / 1,4,19,4,0,0 /
  DATA dft_descr(20)   / 'PW86 exchange + PBE correlation' /
  !
  DATA dft_name(21)    / 'B86BPBE' /
  DATA dft_name2(21)   / 'none' /
  DATA dft_IDs(21,1:6) / 1,4,22,4,0,0 /
  DATA dft_descr(21)   / 'B86b exchange + PBE correlation' /
  !
  DATA dft_name(22)    / 'PBEQ2D' /
  DATA dft_name2(22)   / 'Q2D' /
  DATA dft_IDs(22,1:6) / 1,4,16,9,0,0 /
  DATA dft_descr(22)   / 'PBEQ2D' /
  !
  DATA dft_name(23)    / 'SOGGA' /
  DATA dft_name2(23)   / 'none' /
  DATA dft_IDs(23,1:6) / 1,4,15,4,0,0 /
  DATA dft_descr(23)   / 'SOGGA' /
  !
  DATA dft_name(24)    / 'EV93' /
  DATA dft_name2(24)   / 'none' /
  DATA dft_IDs(24,1:6) / 1,4,22,0,0,0 /
  DATA dft_descr(24)   / 'Engel-Vosko' /
  !
  DATA dft_name(25)    / 'RPBE' /
  DATA dft_name2(25)   / 'none' /
  DATA dft_IDs(25,1:6) / 1,4,41,4,0,0 /
  DATA dft_descr(25)   / 'RPBE' /
  !
  DATA dft_name(26)    / 'PBE0' /
  DATA dft_name2(26)   / 'none' /
  DATA dft_IDs(26,1:6) / 6,4,7,4,0,0 /
  DATA dft_descr(26)   / 'PBE0 in: Perdew, Ernzerhof, Burke, JCP 105, 9982 (1996)' /
  !
  DATA dft_name(27)    / 'B86BPBEX' /
  DATA dft_name2(27)   / 'none' /
  DATA dft_IDs(27,1:6) / 6,4,38,4,0,0 /
  DATA dft_descr(27)   / 'B86bPBE hybrid' /
  !
  DATA dft_name(28)    / 'BHAHLYP' /
  DATA dft_name2(28)   / 'BHANDHLYP' /
  DATA dft_IDs(28,1:6) / 6,4,39,3,0,0 /
  DATA dft_descr(28)   / 'Becke half-and-half LYP' /
  !
  DATA dft_name(29)    / 'HSE' /
  DATA dft_name2(29)   / 'none' /
  DATA dft_IDs(29,1:6) / 1,4,11,4,0,0 /
  DATA dft_descr(29)   / 'Heyd-Scuseria-Ernzerhof (HSE 06, see references)' /
  ! NOTE ABOUT HSE: there are two slight deviations with respect to the HSE06
  ! functional as it is in Gaussian code (that is considered as the reference
  ! in the chemistry community):
  ! - The range separation in Gaussian is precisely 0.11 bohr^-1,
  !   instead of 0.106 bohr^-1 in this implementation
  ! - The gradient scaling relation is a bit more complicated
  !   [ see: TM Henderson, AF Izmaylov, G Scalmani, and GE Scuseria,
  !          J. Chem. Phys. 131, 044108 (2009) ]
  ! These two modifications accounts only for a 1e-5 Ha difference for a
  ! single He atom. Info by Fabien Bruneval.
  !----
  DATA dft_name(30)    / 'GAUP' /
  DATA dft_name2(30)   / 'GAUPBE' /
  DATA dft_IDs(30,1:6) / 1,4,17,4,0,0 /
  DATA dft_descr(30)   / 'Gau-PBE (also gaup)' /
  !
  DATA dft_name(31)    / 'B3LYP' /
  DATA dft_name2(31)   / 'none' /
  DATA dft_IDs(31,1:6) / 7,12,8,6,0,0 /
  DATA dft_descr(31)   / 'B3LYP' /
  !
  DATA dft_name(32)    / 'B3LYP-V1R' /
  DATA dft_name2(32)   / 'none' /
  DATA dft_IDs(32,1:6) / 7,13,8,6,0,0 /
  DATA dft_descr(32)   / 'B3LYP-VWN1-RPA' /
  !
  DATA dft_name(33)    / 'X3LYP' /
  DATA dft_name2(33)   / 'none' /
  DATA dft_IDs(33,1:6) / 9,14,25,10,0,0 /
  DATA dft_descr(33)   / 'X3LYP' /
  !
  DATA dft_name(34)    / 'TPSS' /
  DATA dft_name2(34)   / 'none' /
  DATA dft_IDs(34,1:6) / 1,4,0,0,1,0 /
  DATA dft_descr(34)   / 'TPSS Meta-GGA' /
  !
  DATA dft_name(35)    / 'TPSS-only' /
  DATA dft_name2(35)   / 'none' /
  DATA dft_IDs(35,1:6) / 0,0,0,0,1,0 /
  DATA dft_descr(35)   / 'TPSS Meta-GGA' /
  !
  DATA dft_name(36)    / 'M06L' /
  DATA dft_name2(36)   / 'none' /
  DATA dft_IDs(36,1:6) / 0,0,0,0,2,0 /
  DATA dft_descr(36)   / 'M06L Meta-GGA' /
  !
  DATA dft_name(37)    / 'TB09' /
  DATA dft_name2(37)   / 'none' /
  DATA dft_IDs(37,1:6) / 0,0,0,0,3,0 /
  DATA dft_descr(37)   / 'TB09 Meta-GGA - needs Libxc' /
  !
  DATA dft_name(38)    / 'SCAN' /
  DATA dft_name2(38)   / 'none' /
  DATA dft_IDs(38,1:6) / 0,0,0,0,5,0 /
  DATA dft_descr(38)   / 'SCAN Meta-GGA - needs Libxc.' /
  !
  DATA dft_name(39)    / 'SCAN0' /
  DATA dft_name2(39)   / 'none' /
  DATA dft_IDs(39,1:6) / 0,0,0,0,6,0 /
  DATA dft_descr(39)   / 'SCAN Meta-GGA - needs Libxc.' /
  !
  DATA dft_name(40)    / 'PZ+META' /
  DATA dft_name2(40)   / 'LDA+META' /
  DATA dft_IDs(40,1:6) / 1,1,0,0,4,0 /
  DATA dft_descr(40)   / 'PZ/LDA + null meta-GGA' /
  !
  DATA dft_name(41)    / 'PBE+META' /
  DATA dft_name2(41)   / 'none' /
  DATA dft_IDs(41,1:6) / 1,4,3,4,4,0 /
  DATA dft_descr(41)   / 'PBE + null meta-GGA' /
  !
  !
CONTAINS
  !
  !------------------------------------------------------------------
  SUBROUTINE get_IDs_from_shortname( name, IDs )
    !---------------------------------------------------------------
    !! Get ID numbers of each family-kind term from the DFT shortname.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(INOUT) :: IDs(6)
    INTEGER :: i
    !
    IDs = notset
    !
    DO i = 1, n_dft
      IF (name==dft_name(i) .OR. name==dft_name2(i)) THEN
        IDs(:) = dft_IDs(i,:)
        EXIT
      ENDIF  
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE
  !
  !------------------------------------------------------------------
  SUBROUTINE get_shortname_from_IDs( IDs, name )
    !---------------------------------------------------------------
    !! Get the DFT shortname from ID numbers of each family-kind term.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: IDs(6)
    CHARACTER(LEN=*), INTENT(INOUT) :: name
    INTEGER :: i
    !
    DO i = 1, n_dft
      IF (ALL(IDs(:)==dft_IDs(i,:))) THEN
        name = dft_name(i)
        EXIT
      ENDIF  
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE
  !
  !
END MODULE

