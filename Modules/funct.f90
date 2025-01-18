!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
MODULE funct
  !-------------------------------------------------------------------
  !! This module contains data defining the DFT functional in use
  !! and a number of functions and subroutines to manage them.
  !! All the data and routines related to LDA, GGA and MGGA 
  !! functionals have been moved into the library XClib.  
  !! Here the combinations with nonlocal functionals are still
  !! managed.
  !
  ! Data are PRIVATE and are accessed and set only by function calls.
  !
  !
  USE io_global,      ONLY: stdout, ionode
  USE kinds,          ONLY: DP
  USE beef_interface, ONLY: beef_set_type
  USE xc_lib
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! subroutines/functions managing dft name and indices
  PUBLIC :: set_dft_from_indices, set_dft_from_name
  PUBLIC :: enforce_input_dft, write_dft_name
  PUBLIC :: get_dft_name
  PUBLIC :: get_dft_short, get_dft_long
  PUBLIC :: get_nonlocc_name
  PUBLIC :: get_inlc
  PUBLIC :: dft_is_nonlocc
  ! driver subroutine computing XC non local
  PUBLIC  :: nlc
  ! XC non local index
  PRIVATE :: inlc
  !
  CHARACTER(LEN=37) :: dft = 'not set'
  !
  ! ------------------------------------------------------------------------
  ! "dft" is the exchange-correlation functional label, as set by the user,
  ! using either set_dft_from_name or set_dft_from_indices. It can contain
  ! either the short names or a series of the keywords listed below.
  ! All operations on names are case-insensitive.
  !
  !           short name       complete name       Short description
  !              "pz"    = "sla+pz"            = Perdew-Zunger LDA
  !              "bp"    = "b88+p86"           = Becke-Perdew grad.corr.
  !              "pw91"  = "sla+pw+ggx+ggc"    = PW91 (aka GGA)
  !              "blyp"  = "sla+b88+lyp+blyp"  = BLYP
  !              "pbe"   = "sla+pw+pbx+pbc"    = PBE
  !              "revpbe"= "sla+pw+revx+pbc"   = revPBE (Zhang-Yang)
  !              "rpbe"  = "sla+pw+hhnx+pbc"   = RPBE (Hammer-Hansen-Norskov)
  !              "pw86pbe" = "sla+pw+pw86+pbc" = PW86 exchange + PBE correlation
  !              "b86bpbe" = "sla+pw+b86b+pbc" = B86b exchange + PBE correlation
  !              "pbesol"= "sla+pw+psx+psc"    = PBEsol
  !              "q2d"   = "sla+pw+q2dx+q2dc"  = PBEQ2D
  !              "hcth"  = "nox+noc+hcth+hcth" = HCTH/120
  !              "olyp"  = "nox+lyp+optx+blyp" = OLYP
  !              "wc"    = "sla+pw+wcx+pbc"    = Wu-Cohen
  !              "sogga" = "sla+pw+sox+pbec"   = SOGGA
  !              "optbk88"="sla+pw+obk8+p86"   = optB88
  !              "optb86b"="sla+pw+ob86+p86"   = optB86
  !              "ev93"  = "sla+pw+evx+nogc"   = Engel-Vosko
  !              "tpss"  = "sla+pw+tpss+tpss"  = TPSS Meta-GGA
  !              "m06l"  = "nox+noc+m6lx+m6lc" = M06L Meta-GGA
  !              "tb09"  = "sla+pw+tb09+tb09"  = TB09 Meta-GGA
  !              "pbe0"  = "pb0x+pw+pb0x+pbc"  = PBE0
  !              "b86bx" = "pb0x+pw+b86x+pbc"  = B86bPBE hybrid
  !              "bhahlyp"="pb0x+pw+b88x+blyp" = Becke half-and-half LYP
  !              "hse"   = "sla+pw+hse+pbc"    = Heyd-Scuseria-Ernzerhof (HSE 06, see note below)
  !              "b3lyp"                        = B3LYP
  !              "b3lyp-v1r"                    = B3LYP-VWN1-RPA
  !              "x3lyp"                        = X3LYP
  !              "vwn-rpa" = VWN LDA using vwn1-rpa parametrization
  !              "gaupbe"= "sla+pw+gaup+pbc"   = Gau-PBE (also "gaup")
  !
  !              "vdw-df"       ="sla+pw+revx+vdw1"      = vdW-DF1
  !              "vdw-df2"      ="sla+pw+rw86+vdw2"      = vdW-DF2
  !              "vdw-df-c09"   ="sla+pw+c09x+vdw1"      = vdW-DF-C09
  !              "vdw-df2-c09"  ="sla+pw+c09x+vdw2"      = vdW-DF2-C09
  !              "vdw-df-obk8"  ="sla+pw+obk8+vdw1"      = vdW-DF-obk8 (optB88-vdW)
  !              "vdw-df-ob86"  ="sla+pw+ob86+vdw1"      = vdW-DF-ob86 (optB86b-vdW)
  !              "vdw-df2-b86r" ="sla+pw+b86r+vdw2"      = vdW-DF2-B86R (rev-vdw-df2)
  !              "vdw-df-cx"    ="sla+pw+cx13+vdW1"      = vdW-DF-cx
  !              "vdw-df-cx0"   ="sla+pw+cx13+vdW1+HF/4" = vdW-DF-cx-0
  !              "vdw-df2-0"    ="sla+pw+rw86+vdw2+HF/4" = vdW-DF2-0
  !              "vdw-df2-br0"  ="sla+pw+b86r+vdW2+HF/4" = vdW-DF2-b86r-0
  !              "vdw-df-c090"  ="sla+pw+c09x+vdw1+HF/4" = vdW-DF-C09-0
  !              "vdw-df3-opt1" ="sla+pw+w31x+w31c"      = vdW-DF3-opt1
  !              "vdw-df3-opt2" ="sla+pw+w32x+w32c"      = vdW-DF3-opt2
  !              "vdw-df-C6"    ="sla+pw+b86r+wc6"       = vdW-DF-C6
  !              "rvv10" = "sla+pw+rw86+pbc+vv10"        = rVV10
  !
  ! Any nonconflicting combination of the following keywords is acceptable:
  !
  ! Exchange:    "nox"    none                           iexch=0
  !              "sla"    Slater (alpha=2/3)             iexch=1 (default)
  !              "sl1"    Slater (alpha=1.0)             iexch=2
  !              "rxc"    Relativistic Slater            iexch=3
  !              "oep"    Optimized Effective Potential  iexch=4
  !              "hf"     Hartree-Fock                   iexch=5
  !              "pb0x"   (Slater*0.75+HF*0.25)          iexch=6 for PBE0 and vdW-DF-cx0 and vdW-DF2-0 etc
  !              "b3lp"   B3LYP(Slater*0.80+HF*0.20)     iexch=7
  !              "kzk"    Finite-size corrections        iexch=8
  !              "x3lp"   X3LYP(Slater*0.782+HF*0.218)   iexch=9
  !              "kli"    KLI aproximation for exx       iexch=10
  !
  ! Correlation: "noc"    none                           icorr=0
  !              "pz"     Perdew-Zunger                  icorr=1 (default)
  !              "vwn"    Vosko-Wilk-Nusair              icorr=2
  !              "lyp"    Lee-Yang-Parr                  icorr=3
  !              "pw"     Perdew-Wang                    icorr=4
  !              "wig"    Wigner                         icorr=5
  !              "hl"     Hedin-Lunqvist                 icorr=6
  !              "obz"    Ortiz-Ballone form for PZ      icorr=7
  !              "obw"    Ortiz-Ballone form for PW      icorr=8
  !              "gl"     Gunnarson-Lunqvist             icorr=9
  !              "kzk"    Finite-size corrections        icorr=10
  !              "vwn-rpa" Vosko-Wilk-Nusair, alt param  icorr=11
  !              "b3lp"   B3LYP (0.19*vwn+0.81*lyp)      icorr=12
  !              "b3lpv1r"  B3LYP-VWN-1-RPA
  !                         (0.19*vwn_rpa+0.81*lyp)      icorr=13
  !              "x3lp"   X3LYP (0.129*vwn_rpa+0.871*lyp)icorr=14
  !
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
  !              "revx"   revised PBE by Zhang-Yang      igcx =4
  !              "hcth"   Cambridge exch, Handy et al    igcx =5
  !              "optx"   Handy's exchange functional    igcx =6
  !              "pb0x"   PBE0 (PBE exchange*0.75)       igcx =8
  !              "b3lp"   B3LYP (Becke88*0.72)           igcx =9
  !              "psx"    PBEsol exchange                igcx =10
  !              "wcx"    Wu-Cohen                       igcx =11
  !              "hse"    HSE screened exchange          igcx =12
  !              "rw86"   revised PW86                   igcx =13
  !              "pbe"    same as PBX, back-comp.        igcx =14
  !              "c09x"   Cooper 09                      igcx =16
  !              "sox"    sogga                          igcx =17
  !              "q2dx"   Q2D exchange grad corr         igcx =19
  !              "gaup"   Gau-PBE hybrid exchange        igcx =20
  !              "pw86"   Perdew-Wang (1986) exchange    igcx =21
  !              "b86b"   Becke (1986) exchange          igcx =22
  !              "obk8"   optB88  exchange               igcx =23
  !              "ob86"   optB86b exchange               igcx =24
  !              "evx"    Engel-Vosko exchange           igcx =25
  !              "b86r"   revised Becke (b86b)           igcx =26
  !              "cx13"   consistent exchange            igcx =27
  !              "x3lp"   X3LYP (Becke88*0.542 +
  !                              Perdew-Wang91*0.167)    igcx =28
  !              "cx0"    vdW-DF-cx+HF/4 (cx13-0)        igcx =29 
  !              "r860"   rPW86+HF/4 (rw86-0)            igcx =30 (for DF0)
  !              "cx0p"   vdW-DF-cx+HF/5 (cx13-0p)       igcx =31 
  !              "ahcx"   vdW-DF-cx based analytic hole  igcx =32 ! Launched vdW-DF-ahcx - PH
  !              "ahf2"   vdW-DF2 based analytic hole    igcx =33 ! Defined vdw-DF2-AH at 0.20 - PH
  !              "ahpb"   PBE based analytic hole        igcx =34 ! PBE-AH (rHJS-PBE) at 0.20 - PH
  !              "ahps"   PBE-sol based analytic hole    igcx =35 ! PBESOL-AH (rHJS-PBEsol) at 0.20 - PH
  !              "cx14"   Exporations                    igcx =36 reserved PH
  !              "cx15"   Exporations                    igcx =37 reserved PH
  !              "br0"    vdW-DF2-b86r+HF/4 (b86r-0)     igcx =38 ! Tested (2022)
  !              "cx16"   Exporations                    igcx =39 reserved PH
  !              "c090"   vdW-DF-c09+HF/4 (c09-0)        igcx =40 
  !              "b86x"   B86b exchange * 0.75           igcx =41
  !              "b88x"   B88 exchange * 0.50            igcx =42
  !              "beex"   BEE exchange                   igcx =43 
  !              "hhnx"   Hammer-Hansen-Norskov          igcx =44
  !              "w31x"   vdW-DF3-opt1 exchange          igcx =45
  !              "w32x"   vdW-DF3-opt2 exchange          igcx =46
  !              "ahbr"   vdW-DF2-ahbr exchange          igcx =47 ! for vdW-DF2-ahbr 
  !              "ehpb"   HSE variant                    igcx =48 ! Reserved PH
  !              "hjpb"   HJS-type PBE cross check       igcx =49 ! Reserved PH
  !              "hjps"   HJS-type PBEsol crosscheck     igcx =50 ! Reserved PH
  !
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Lee-Yang-Parr                  igcc =3
  !              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
  !              "hcth"   Cambridge corr, Handy et al    igcc =5
  !              "b3lp"   B3LYP (Lee-Yang-Parr*0.81)     igcc =7
  !              "psc"    PBEsol corr                    igcc =8
  !              "pbe"    same as PBX, back-comp.        igcc =9
  !              "q2dc"   Q2D correlation grad corr      igcc =12
  !              "x3lp"   X3LYP (Lee-Yang-Parr*0.871)    igcc =13
  !              "beec"   BEE correlation                igcc =14
  !
  ! Meta-GGA functionals
  !              "tpss"   TPSS Meta-GGA                  imeta=1
  !              "m6lx"   M06L Meta-GGA                  imeta=2
  !              "tb09"   TB09 Meta-GGA                  imeta=3
  !              "+meta"  activate MGGA even without MGGA-XC   imeta=4
  !              "scan"   SCAN Meta-GGA                  imeta=5
  !              "sca0"   SCAN0  Meta-GGA                imeta=6
  !              "r2scan" R2SCAN Meta-GGA                imeta=7
  !
  ! van der Waals functionals (nonlocal term only)
  !              "nonlc"  none                           inlc =0 (default)
  !--------------inlc = 1 to inlc = 25 reserved for vdW-DF--------------
  !              "vdw1"   vdW-DF1                        inlc =1
  !              "vdw2"   vdW-DF2                        inlc =2
  !              "w31c"   vdW-DF3-opt1                   inlc =3
  !              "w32c"   vdW-DF3-opt2                   inlc =4
  !              "wc6"    vdW-DF-C6                      inlc =5
  !---------------------------------------------------------------------
  !              "vv10"   rVV10                          inlc =26
  !
  ! Meta-GGA with van der Waals
  !              "rvv10-scan" rVV10 (with b=15.7) and scan inlc=26 (PRX 6, 041005 (2016))
  !
  ! Note: as a rule, all keywords should be unique, and should be different
  ! from the short name, but there are a few exceptions.
  !
  ! References:
  !              pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981)
  !              vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
  !              vwn1-rpa S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
  !              wig     E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)
  !              hl      L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)
  !              gl      O.Gunnarsson and B.I.Lundqvist, PRB 13, 4274 (1976)
  !              pw      J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)
  !              obpz    G.Ortiz and P.Ballone, PRB 50, 1391 (1994)
  !              obpw    as above
  !              b88     A.D.Becke, PRA 38, 3098 (1988)
  !              p86     J.P.Perdew, PRB 33, 8822 (1986)
  !              pw86    J.P.Perdew, PRB 33, 8800 (1986)
  !              b86b    A.D.Becke, J.Chem.Phys. 85, 7184 (1986)
  !              ob86    Klimes, Bowler, Michaelides, PRB 83, 195131 (2011)
  !              b86r    I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
  !              w31x    D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)
  !              w32x    D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)
  !              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  !              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
  !              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
  !              hcth    Handy et al, JCP 109, 6264 (1998)
  !              olyp    Handy et al, JCP 116, 5411 (2002)
  !              revPBE  Zhang and Yang, PRL 80, 890 (1998)
  !              pbesol  J.P. Perdew et al., PRL 100, 136406 (2008)
  !              q2d     L. Chiodo et al., PRL 108, 126402 (2012)
  !              rw86    Eamonn D. Murray et al, J. Chem. Theory Comput. 5, 2754 (2009)
  !              wc      Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
  !              kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
  !              pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
  !              hse     Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
  !                      Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006).
  !              b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch
  !                      J.Phys.Chem 98, 11623 (1994)
  !              x3lyp   X. Xu, W.A Goddard III, PNAS 101, 2673 (2004)
  !              vdW-DF       M. Dion et al., PRL 92, 246401 (2004)
  !                           T. Thonhauser et al., PRL 115, 136402 (2015)
  !              vdW-DF2      Lee et al., Phys. Rev. B 82, 081101 (2010)
  !              rev-vdW-DF2  I. Hamada, Phys. Rev. B 89, 121103(R) (2014)
  !              vdW-DF-cx    K. Berland and P. Hyldgaard, PRB 89, 035412 (2014)
  !              vdW-DF-cx0   K. Berland, Y. Jiao, J.-H. Lee, T. Rangel, J. B. Neaton and P. Hyldgaard,
  !                           J. Chem. Phys. 146, 234106 (2017)
  !              vdW-DF-cx0p  Y. Jiao, E. Schröder and P. Hyldgaard, 
  !                           J. Chem. Phys. 148, 194115 (2018)
  !              vdW-DF-obk8  Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)
  !              vdW-DF-ob86  Klimes et al, Phys. Rev. B, 83, 195131 (2011)
  !              vdW-DF3-opt1 D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)
  !              vdW-DF3-opt2 D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)
  !              vdW-DF-C6    K. Berland, D. Chakraborty, and T. Thonhauser, PRB 99, 195418 (2019)
  !              c09x    V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
  !              tpss    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria,
  !                      PRL 91, 146401 (2003)
  !              tb09    F Tran and P Blaha, Phys.Rev.Lett. 102, 226401 (2009)
  !              scan    J Sun, A Ruzsinszky and J Perdew, PRL 115, 36402 (2015)
  !              scan0   K Hui and J-D. Chai, JCP 144, 44114 (2016)
  !              r2scan  J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew,
  !                      and J. Sun, JPCL 11, 8208 (2020)
  !              sogga   Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)
  !              m06l    Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)
  !              gau-pbe J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)
  !              rVV10   R. Sabatini et al. Phys. Rev. B 87, 041108(R) (2013)
  !              ev93     Engel-Vosko, Phys. Rev. B 47, 13164 (1993)
  !              vdW-DF-ahcx V. Shukla, Y. Jiao, iC. M. Frostenson, and P. Hyldgaard, JPCM 34, 025902 (2022)
  !              vdW-DF2-ah V. Shukla, Y. Jiao, iC. M. Frostenson, and P. Hyldgaard, JPCM 34, 025902 (2022)
  !              PBE-ah V. Shukla, Y. Jiao, iC. M. Frostenson, and P. Hyldgaard, JPCM 34, 025902 (2022)
  !              PBESOL-ah V. Shukla, Y. Jiao, iC. M. Frostenson, and P. Hyldgaard, JPCM 34, 025902 (2022)
  !
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
  !
  ! NOTE FOR LIBXC USERS: to use libxc functionals you must enforce them from input (use
  ! 'input_dft' in &system) and write their IDs in the input string. The only notation
  ! now allowed (v7.0) for input DFTs containing Libxc terms is:
  ! XC-000i-000i-000i-000i-000i-000i
  ! where you put the functional IDs instead of the zeros and an 'L' instead of
  ! 'i' if the functional is from Libxc. The order is the usual one:  
  ! LDAexch - LDAcorr - GGAexch - GGAcorr - MGGAexch - MGGAcorr  
  ! however QE will automatically adjust it if needed. You can skip zero tails (e.g.
  ! you don't need GGA/MGGA slots if the dft is LDA only and so on.
  ! You can use combinations of qe and libxc functionals, when they are compatible.
  ! You can also add vdW terms after it, for example, sla+pw+rw86+vdw2 is:
  ! input_dft='XC-001i-004i-013i-vdw2'.
  ! For more details see the user_guide (in 'Doc' folder).
  !
  INTEGER, PARAMETER :: notset = -1
  !
  INTEGER :: inlc = notset
  !
  LOGICAL :: isnonlocc = .FALSE.
  !
  LOGICAL :: discard_input_dft = .FALSE.
  !
  INTEGER  :: beeftype = -1
  INTEGER  :: beefvdw = 0
  !
  INTEGER, PARAMETER :: ncnl = 26
  CHARACTER(LEN=4) :: nonlocc
  DIMENSION :: nonlocc(0:ncnl)
  !
  DATA nonlocc/ 'NONE', 'VDW1', 'VDW2', 'W31C', 'W32C', 'WC6', 20*'NONE', 'VV10' /
  !
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !! It sets the dft functional IDs and parameters from the input name. It 
    !! directly calls the XClib routines and functions to set the LDA, GGA,
    !! MGGA terms (but not the non-local one).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    !
    ! ... local variables
    !
    INTEGER :: len, l, i
    CHARACTER(len=150) :: dftout, dftout_loc
    LOGICAL :: dft_defined
    LOGICAL :: check_libxc
    !
    CHARACTER(LEN=1), EXTERNAL :: capital
    CHARACTER(LEN=4) :: lda_exch, lda_corr, gga_exch, gga_corr
    !
    INTEGER :: save_inlc, lnt, ln_nlc
    INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac
    !
    ! Exit if set to discard further input dft
    !
    IF ( discard_input_dft ) RETURN
    !
    ! save current status of XC indices
    !
    dft_defined = .FALSE.
    !
    save_inlc  = inlc
    !
    ! convert to uppercase
    !
    len = LEN_TRIM(dft_)
    dftout = ' '
    !
    DO l = 1, len
       dftout(l:l) = capital( dft_(l:l) )
    ENDDO
    !
    !
    ! ----------------------------------------------
    ! NOW WE CHECK ALL THE SHORT NAMES
    ! Note: comparison is done via exact matching
    ! ----------------------------------------------
    !
    SELECT CASE( TRIM(dftout) )
    ! special case : case BEEF (default: BEEF-vdW-DF2)
    CASE('BEEF', 'BEEF-VDW')
       IF (LEN_TRIM(dftout) == 4) THEN
          beeftype = 0
       ELSE
          SELECT CASE(TRIM(dftout(5:)))
             CASE('-VDW')
                beeftype = 0
             CASE DEFAULT
                READ(dftout(5:), '(i1)', IOSTAT=i) beeftype
                IF (i /= 0) CALL errore('set_dft_from_name', &
                                      & 'unknown BEEF type', 1)
          END SELECT
       ENDIF
       IF (.NOT. beef_set_type(beeftype, ionode)) &
       & CALL errore('set_dft_from_name', 'unknown BEEF type number', 1)
       SELECT CASE(beeftype)
          CASE(0)
             ! turn on vdW-DF2 type interactions for BEEF-vdW
             beefvdw = 2
       END SELECT
       dft_defined = xclib_set_dft_IDs(1,4,43,14,0,0)
       inlc = beefvdw
    ! Special case vdW-DF
    CASE( 'VDW-DF' )
       dft_defined = xclib_set_dft_IDs(1,4,4,0,0,0)
       inlc = 1
    ! Special case vdW-DF2
    CASE( 'VDW-DF2' )
       dft_defined = xclib_set_dft_IDs(1,4,13,0,0,0)
       inlc = 2
    ! Special case vdW-DF3-opt1
    CASE( 'VDW-DF3-OPT1' )
       dft_defined = xclib_set_dft_IDs(1,4,45,0,0,0)
       inlc = 3
    ! Special case vdW-DF3-opt2
    CASE( 'VDW-DF3-OPT2' )
       dft_defined = xclib_set_dft_IDs(1,4,46,0,0,0)
       inlc = 4
    ! Special case vdW-DF-C6
    CASE( 'VDW-DF-C6' )
       dft_defined = xclib_set_dft_IDs(1,4,26,0,0,0)
       inlc = 5
    ! Special case vdW-DF with C09 exchange
    CASE( 'VDW-DF-C09' )
       dft_defined = xclib_set_dft_IDs(1,4,16,0,0,0)
       inlc = 1
    ! Special case vdW-DF2 with C09 exchange
    CASE( 'VDW-DF2-C09' )
       dft_defined = xclib_set_dft_IDs(1,4,16,0,0,0)
       inlc = 2
    ! Special case vdW-DF-obk8, or vdW-DF + optB88
    CASE( 'VDW-DF-OBK8' )
       dft_defined = xclib_set_dft_IDs(1,4,23,0,0,0)
       inlc = 1
    ! Special case vdW-DF-ob86, or vdW-DF + optB86
    CASE( 'VDW-DF-OB86' )
       dft_defined = xclib_set_dft_IDs(1,4,24,0,0,0)
       inlc = 1
    ! Special case vdW-DF2 with B86R
    CASE( 'VDW-DF2-B86R' )
       dft_defined = xclib_set_dft_IDs(1,4,26,0,0,0)
       inlc = 2
    ! Special case vdW-DF-CX
    CASE( 'VDW-DF-CX' )
       dft_defined = xclib_set_dft_IDs(1,4,27,0,0,0)
       inlc = 1
    ! Special case vdW-DF-CX0
    CASE( 'VDW-DF-CX0' )
       dft_defined = xclib_set_dft_IDs(6,4,29,0,0,0)
       inlc = 1
    ! Special case vdW-DF-CX0P
    CASE( 'VDW-DF-CX0P' )
       dft_defined = xclib_set_dft_IDs(6,4,31,0,0,0)
       inlc = 1
    ! Special case vdW-DF-AHCX
    CASE( 'VDW-DF-AHCX' )
       dft_defined = xclib_set_dft_IDs(1,4,32,0,0,0)
       inlc = 1
    ! Special case vdW-DF2-0
    CASE( 'VDW-DF2-0' )
       dft_defined = xclib_set_dft_IDs(6,4,30,0,0,0)
       inlc = 2
    ! Special case vdW-DF2-AH
    CASE( 'VDW-DF2-AH' )
       dft_defined = xclib_set_dft_IDs(1,4,33,0,0,0)
       inlc = 2
    ! Special case vdW-DF2-BR0
    CASE( 'VDW-DF2-BR0' )
       dft_defined = xclib_set_dft_IDs(6,4,38,0,0,0)
       inlc = 2
    ! Special case vdW-DF2-AHBR
    CASE( 'VDW-DF2-AHBR' )
       dft_defined = xclib_set_dft_IDs(1,4,47,0,0,0)
       inlc = 2
    ! Special case vdW-DF-C090
    CASE( 'VDW-DF-C090' )
       dft_defined = xclib_set_dft_IDs(6,4,40,0,0,0)
       inlc = 1
    ! Special case rVV10
    CASE( 'RVV10' )
       dft_defined = xclib_set_dft_IDs(1,4,13,4,0,0)
       inlc = 26
    ! Special case rVV10+scan
    CASE( 'RVV10-SCAN' )
       dft_defined = xclib_set_dft_IDs(0,0,0,0,5,0)
       inlc = 26
    !
    CASE( 'REV-VDW-DF2' )
       CALL errore( 'set_dft_from_name', 'obsolete XC label, use VDW-DF2-B86R', 1 )
    !
    CASE( 'VDW-DF3' )
       CALL errore( 'set_dft_from_name', 'obsolete XC label, use VDW-DF-OBK8', 1 )
    !
    CASE( 'VDW-DF4', 'OPTB86B-VDW' )
       CALL errore( 'set_dft_from_name', 'obsolete XC label, use VDW-DF-OB86', 1 )
    ! Special case vdW-DF-X
    CASE( 'VDW-DF-X' )
       CALL errore( 'set_dft_from_name', 'functional not yet implemented', 1 )
    ! Special case vdW-DF-Y
    CASE( 'VDW-DF-Y' )
       CALL errore( 'set_dft_from_name', 'functional not yet implemented', 1 )
    ! Special case vdW-DF-Z
    CASE( 'VDW-DF-Z' )
       CALL errore( 'set_dft_from_name', 'functional not yet implemented', 1 )
    ! Case for old RRKJ format, containing indices instead of label
    CASE DEFAULT
       !
       IF ('INDEX:' == dftout(1:6)) THEN
          READ( dftout(7:18), '(6i2)') iexch, icorr, igcx, igcc, inlc, imeta
          dft_defined = xclib_set_dft_IDs(iexch, icorr, igcx, igcc, imeta, 0)
          CALL xclib_get_name('LDA','EXCH', lda_exch)
          CALL xclib_get_name('LDA','CORR', lda_corr)
          CALL xclib_get_name('GGA','EXCH', gga_exch)
          CALL xclib_get_name('GGA','CORR', gga_corr)
          !
          dftout = TRIM(lda_exch) //'-'// &
                   TRIM(lda_corr) //'-'// &
                   TRIM(gga_exch) //'-'// &
                   TRIM(gga_corr) //'-'// nonlocc(inlc)
       ELSE
          !
          dftout_loc = ''
          inlc = matching( dftout, ncnl, nonlocc )
          IF ( inlc/=0 .AND. dftout(1:3) == 'XC-' ) THEN
            lnt = LEN_TRIM(dftout)
            ln_nlc = LEN_TRIM(nonlocc(inlc))
            dftout_loc(1:lnt-ln_nlc) = dftout(1:lnt-ln_nlc)
          ELSE
            dftout_loc = dftout
          ENDIF
          CALL xclib_set_dft_from_name( TRIM(dftout_loc) )
          dft_defined = .TRUE.
          !
       ENDIF
       !
    END SELECT
    !
    !----------------------------------------------------------------
    ! Last check
    ! No more defaults, the code exits if the dft is not defined
    !----------------------------------------------------------------
    !
    iexch = xclib_get_id('LDA','EXCH')
    icorr = xclib_get_id('LDA','CORR')
    igcx  = xclib_get_id('GGA','EXCH')
    igcc  = xclib_get_id('GGA','CORR')
    imeta = xclib_get_id('MGGA','EXCH')
    imetac = xclib_get_id('MGGA','CORR')
    !
    IF (igcx == 6 .AND. .NOT.xclib_dft_is_libxc('GGA','EXCH') ) &
                CALL infomsg( 'set_dft_from_name', 'OPTX untested! please test' )
    !
    ! check for unrecognized labels
    !
    IF ( iexch<=0 .AND. icorr<=0 .AND. igcx<=0 .AND. igcc<=0 .AND. imeta<=0 .AND. imetac<=0) THEN
       IF ( inlc <= 0 .AND. TRIM(dftout) /= 'NOX-NOC') THEN
          CALL errore( 'set_dft_from_name', TRIM(dftout)//': unrecognized dft', 1 )
       ELSE
          ! if inlc is the only nonzero index the label is likely wrong
          CALL errore( 'set_dft_from_name', TRIM(dftout)//': strange dft, please check', inlc )
       ENDIF
    ENDIF
    !
    ! Fill variables and exit
    !
    dft = dftout
    !
    dft_defined = .TRUE.
    !
    isnonlocc = (inlc > 0)
    !
    CALL xclib_set_auxiliary_flags( isnonlocc )
    !
    ! check non-local term has not been previously set differently
    !
    IF (save_inlc /= notset  .AND. save_inlc /= inlc)   THEN
       WRITE (stdout,*) inlc, save_inlc
       CALL errore( 'set_dft_from_name', ' conflicting values for inlc', 1 )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE set_dft_from_name
  !
  !
  !-----------------------------------------------------------------
  INTEGER FUNCTION matching( dft, n, name )
    !-----------------------------------------------------------------
    !! Looks for matches between the names of each single term of the 
    !! xc-functional and the input dft string.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN):: n
    CHARACTER(LEN=*), INTENT(IN):: name(0:n)
    CHARACTER(LEN=*), INTENT(IN):: dft
    INTEGER :: i
    LOGICAL, EXTERNAL :: matches
    !
    matching = notset
    !
    DO i = n, 0, -1
       IF ( matches(name(i), TRIM(dft)) ) THEN
          !
          IF ( matching == notset ) THEN
           !WRITE(*, '("matches",i2,2X,A,2X,A)') i, name(i), TRIM(dft)
           matching = i
          ELSE
             WRITE(*, '(2(2X,i2,2X,A))') i, TRIM(name(i)), &
                                  matching, TRIM(name(matching))
             CALL errore( 'set_dft', 'two conflicting matching values', 1 )
          ENDIF
       ENDIF
    ENDDO
    !
    IF (matching == notset) matching = 0
    !
  END FUNCTION matching
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE enforce_input_dft( dft_, nomsg )
    !---------------------------------------------------------------------
    !! Translates a string containing the exchange-correlation name
    !! into internal indices and force any subsequent call to 
    !! \(\textrm{set_dft_from_name}\) to return without changing them.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: dft_
    LOGICAL, INTENT(IN), OPTIONAL :: nomsg
    !
    CALL set_dft_from_name( dft_ )
    IF (dft == 'not set') CALL errore( 'enforce_input_dft', 'cannot fix unset dft', 1 )
    discard_input_dft = .TRUE.
    !
    IF ( PRESENT(nomsg) ) RETURN
    !
    WRITE(stdout,'(/,5x,a)') "IMPORTANT: XC functional enforced from input :"
    CALL write_dft_name
    WRITE(stdout,'(5x,a)') "Any further DFT definition will be discarded"
    WRITE(stdout,'(5x,a/)') "Please, verify this is what you really want"
    !
    RETURN
    !
  END SUBROUTINE enforce_input_dft
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION get_inlc()
    !! Get dft index for non-local term.
    INTEGER :: get_inlc
    get_inlc = inlc
    RETURN
  END FUNCTION get_inlc
  !-----------------------------------------------------------------------
  FUNCTION get_nonlocc_name()
    !! Get dft name for non-local term.
    CHARACTER(10) get_nonlocc_name
    get_nonlocc_name = TRIM(nonlocc(inlc))
    RETURN
  END FUNCTION get_nonlocc_name
  !-----------------------------------------------------------------------
  FUNCTION dft_is_nonlocc()
    !! TRUE if dft is non-local.
    LOGICAL :: dft_is_nonlocc
    dft_is_nonlocc = isnonlocc
    RETURN
  END FUNCTION dft_is_nonlocc
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION get_dft_name()
    !! Get the string with the full dft name.
    CHARACTER(LEN=37) :: get_dft_name
    get_dft_name = dft
    RETURN
  END FUNCTION get_dft_name
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE set_dft_from_indices( iexch_, icorr_, igcx_, igcc_, imeta_, inlc_ )
     !--------------------------------------------------------------------
     !! Set dft functional from the IDs of each term - OBSOLESCENT:
     !! for compatibility with old PPs only, metaGGA not accounted for
     !
     IMPLICIT NONE
     !
     INTEGER :: iexch_, icorr_, igcx_, igcc_, imeta_, inlc_
     INTEGER :: iexch, icorr, igcx, igcc, imeta
     CHARACTER(LEN=4) :: lda_exch, lda_corr, gga_exch, gga_corr
     LOGICAL :: dft_defined
     !
     IF ( discard_input_dft ) RETURN
     !
     iexch  = xclib_get_id( 'LDA', 'EXCH' )
     icorr  = xclib_get_id( 'LDA', 'CORR' )
     igcx   = xclib_get_id( 'GGA', 'EXCH' )
     igcc   = xclib_get_id( 'GGA', 'CORR' )
     imeta  = xclib_get_id( 'MGGA','EXCH' )
     !
     IF (iexch == notset) iexch = iexch_
     IF (iexch /= iexch_) THEN
        write (stdout,*) iexch, iexch_
        CALL errore( 'set_dft', ' conflicting values for iexch', 1 )
     ENDIF
     IF (icorr == notset) icorr = icorr_
     IF (icorr /= icorr_) THEN
        write (stdout,*) icorr, icorr_
        CALL errore( 'set_dft', ' conflicting values for icorr', 1 )
     ENDIF
     IF (igcx  == notset) igcx = igcx_
     IF (igcx /= igcx_) THEN
        write (stdout,*) igcx, igcx_
        CALL errore( 'set_dft', ' conflicting values for igcx', 1 )
     ENDIF
     IF (igcc  == notset) igcc = igcc_
     IF (igcc /= igcc_) THEN
        write (stdout,*) igcc, igcc_
        CALL errore( 'set_dft', ' conflicting values for igcc', 1 )
     ENDIF
     IF (imeta  == notset) imeta = imeta_
     IF (imeta /= imeta_) THEN
        write (stdout,*) imeta, imeta_
        CALL errore( 'set_dft', ' conflicting values for imeta', 1 )
     ENDIF     
     IF (imeta /= 0 ) CALL errore( 'set_dft', ' META-GGA not allowed', 1 )
     !
     IF (inlc  == notset) inlc = inlc_
     IF (inlc /= inlc_) THEN
        write (stdout,*) inlc, inlc_
        CALL errore( 'set_dft', ' conflicting values for inlc', 1 )
     ENDIF
     CALL xclib_get_name('LDA','EXCH', lda_exch)
     CALL xclib_get_name('LDA','CORR', lda_corr)
     CALL xclib_get_name('GGA','EXCH', gga_exch)
     CALL xclib_get_name('GGA','CORR', gga_corr)
     !
     dft = TRIM(lda_exch) //'-'// &
           TRIM(lda_corr) //'-'// &
           TRIM(gga_exch) //'-'// &
           TRIM(gga_corr) //'-'// nonlocc(inlc)
     !
     dft_defined = xclib_set_dft_IDs(iexch,icorr,igcx,igcc,imeta,0)
     !
     ! WRITE( stdout,'(a)') dft
     isnonlocc = (inlc > 0)
     CALL xclib_set_auxiliary_flags( isnonlocc )
     !
     RETURN
  END SUBROUTINE set_dft_from_indices
  !
  !
  !-------------------------------------------------------------------------
  FUNCTION get_dft_short()
    !---------------------------------------------------------------------
    !! It gets a short version (if exists) of the name of the dft in use.  
    !! If there is no non-local term directly calls the xclib analogous 
    !! routine (\(\texttt{xclib_get_dft_short}\)).
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=37) :: get_dft_short
    CHARACTER(LEN=37) :: shortname
    INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac
    !
    shortname = 'no shortname'
    !
    IF (inlc == 0) THEN
      shortname = xclib_get_dft_short()
    ELSE
      !
      iexch  = xclib_get_id( 'LDA', 'EXCH' )
      icorr  = xclib_get_id( 'LDA', 'CORR' )
      igcx   = xclib_get_id( 'GGA', 'EXCH' )
      igcc   = xclib_get_id( 'GGA', 'CORR' )
      !
      ! ... inlc==1
      !
      IF (iexch==1 .AND. icorr==4 .AND. igcx==4 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==27 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-CX'
      ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==29 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-CX0'
      ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==31 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-CX0P'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==32 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-AHCX'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==16 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-C09'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==24 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-OB86'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==23 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-OBK8'
      ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==40 .AND. igcc==0 .AND. inlc==1) THEN
         shortname = 'VDW-DF-C090'
      !
      ! ... inlc==2
      !
      ELSEIF (iexch==1 .AND. icorr==4  .AND. igcx==43 .AND. igcc==14 .AND. inlc==2) THEN
         shortname = 'BEEF'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==13 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==16 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2-C09'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==26 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2-B86R'
      ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==30 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2-0'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==33 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2-AH'
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==47 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2-AHBR'
      ELSEIF (iexch==6 .AND. icorr==4 .AND. igcx==38 .AND. igcc==0 .AND. inlc==2) THEN
         shortname = 'VDW-DF2-BR0'
      !
      ! ... inlc==3
      !
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==45 .AND. igcc==0 .AND. inlc==3) THEN
        shortname = 'VDW-DF3-OPT1'
      !
      ! ... inlc==4
      !
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==46 .AND. igcc==0 .AND. inlc==4) THEN
        shortname = 'VDW-DF3-OPT2'
      !
      ! ... inlc==5
      !
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==26 .AND. igcc==0 .AND. inlc==5) THEN
        shortname = 'VDW-DF-C6'
      !
      ! ... inlc==26
      !
      ELSEIF (iexch==1 .AND. icorr==4 .AND. igcx==13 .AND. igcc==4 .AND. inlc==26) THEN
        shortname = 'RVV10'
      !
      ! ... all the other combinations
      !
      ELSE
        shortname = xclib_get_dft_short()
        shortname = TRIM(shortname)//'-'//TRIM(nonlocc(inlc))
      ENDIF
      !
    ENDIF
    !
    get_dft_short = shortname
    !
  END FUNCTION get_dft_short
  !
  !
  !---------------------------------------------------------------------
  FUNCTION get_dft_long()
    !---------------------------------------------------------------------
    !! Returns a string containing the name of each term of the dft functional.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=37) :: get_dft_long
    CHARACTER(LEN=37) :: longname
    !
    !WRITE(longname,'(4a5)') exc(iexch), corr(icorr), gradx(igcx), gradc(igcc)
    !
    longname = xclib_get_dft_long()
    !
    IF ( inlc > 0 ) longname = longname(1:20)//TRIM(nonlocc(inlc))
    !
    get_dft_long = longname
    !
  END FUNCTION get_dft_long
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_dft_name
    !-----------------------------------------------------------------------
    !! Print on output the name of each term of the dft functional.
    !
    IMPLICIT NONE
    !
    INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac
    !
    WRITE( stdout, '(5X,"Exchange-correlation= ",A)') TRIM( dft )
    iexch  = xclib_get_id( 'LDA', 'EXCH' )
    icorr  = xclib_get_id( 'LDA', 'CORR' )
    igcx   = xclib_get_id( 'GGA', 'EXCH' )
    igcc   = xclib_get_id( 'GGA', 'CORR' )
    imeta  = xclib_get_id( 'MGGA','EXCH' )
    imetac = xclib_get_id( 'MGGA','CORR' )
    !
    WRITE( stdout, '(27X,"(",I4,3I4,3I4,")")' ) iexch, icorr, igcx, igcc, inlc, &
                                                imeta, imetac
    IF ( xclib_get_exx_fraction() > 0.0_dp ) WRITE( stdout, &
         '(5X,"EXX-fraction              =",F12.2)') xclib_get_exx_fraction()
    RETURN
  END SUBROUTINE write_dft_name
  !
  !
  !-----------------------------------------------------------------------
  !------- NONLOCAL CORRECTIONS DRIVER ----------------------------------
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE nlc( rho_valence, rho_core, nspin, enl, vnl, v )
    !-----------------------------------------------------------------------
    !! Non-local contribution to the correlation energy.
    !
    !     input      :  rho_valence, rho_core
    !     definition :  E_nl = \int E_nl(rho',grho',rho'',grho'',|r'-r''|) dr
    !     output     :  enl = E^nl_c
    !                   vnl = D(E^nl_c)/D(rho)
    !                   v   = non-local contribution to the potential
    !
    !
    USE vdW_DF, ONLY: xc_vdW_DF, xc_vdW_DF_spin, inlc_ => inlc
    USE rVV10,  ONLY: xc_rVV10
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: rho_valence(:,:), rho_core(:)
    INTEGER,  INTENT(IN)    :: nspin
    REAL(DP), INTENT(INOUT) :: v(:,:)
    REAL(DP), INTENT(INOUT) :: enl, vnl
    !
    IF ( inlc > 0 .AND. inlc < 26 ) THEN
      !
      inlc_ = inlc
      IF ( nspin == 1 ) THEN
         CALL xc_vdW_DF      (rho_valence, rho_core, enl, vnl, v)
      ELSE IF ( nspin == 2 ) THEN
         CALL xc_vdW_DF_spin (rho_valence, rho_core, enl, vnl, v)
      ELSE
         CALL errore ('nlc', 'vdW-DF not available for noncollinear spin case',1)
      END If
      !
    ELSE IF ( inlc == 26 ) THEN
      !
      IF ( xclib_get_id('MGGA','EXCH') == 0 ) THEN
        CALL xc_rVV10 (rho_valence(:,1), rho_core, nspin, enl, vnl, v)
      ELSE
        CALL xc_rVV10 (rho_valence(:,1), rho_core, nspin, enl, vnl, v, 15.7_dp)
      END IF
      !
    ELSE
      !
      CALL errore ('nlc', 'inlc choice for E^nl_c not implemented',1)
      !
    END IF
    !
    RETURN
    !
  END SUBROUTINE nlc
  !
  !
END MODULE funct
