!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
module funct
!-------------------------------------------------------------------
! This module contains data defining the DFT functional in use 
! and a number of functions and subroutines to manage them.
! Data are PRIVATE and are accessed and set only by function calls.
! Basic drivers to compute XC quantities are also included.
!  
!  setting routines:   set_dft_from_name (previously which_dft)
!                      set_dft_from_indices
!                      enforce_input_dft
!                      start_exx
!                      stop_exx
!                      set_finite_size_volume
!  retrieve functions: get_dft_name, get_dft_short, get_dft_long
!                      get_iexch
!                      get_icorr
!                      get_igcx
!                      get_igcc
!                      get_exx_fraction
!                      write_dft_name
!  logical functions:  dft_is_gradient
!                      dft_is_meta
!                      dft_is_hybrid
!                      dft_is_nonlocc
!                      exx_is_active
!                      dft_has_finite_size_correction
!
!  XC computation drivers: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
!  derivatives of XC computation drivers: dmxc, dmxc_spin, dmxc_nc, dgcxc,
!                                         dgcxc_spin
!
  USE io_global, ONLY: stdout
  USE kinds,     ONLY: DP
  IMPLICIT NONE
  PRIVATE
  SAVE
  ! subroutines/functions managing dft name and indices
  PUBLIC  :: set_dft_from_indices, set_dft_from_name
  PUBLIC  :: enforce_input_dft, write_dft_name
  PUBLIC  :: get_dft_name, get_dft_short, get_dft_long,&
             get_nonlocc_name
  PUBLIC  :: get_iexch, get_icorr, get_igcx, get_igcc, get_meta, get_inlc
  PUBLIC  :: dft_is_gradient, dft_is_meta, dft_is_hybrid, dft_is_nonlocc, igcc_is_lyp

  ! additional subroutines/functions for hybrid functionals
  PUBLIC  :: start_exx, stop_exx, get_exx_fraction, exx_is_active
  PUBLIC  :: set_exx_fraction
  PUBLIC  :: set_screening_parameter, get_screening_parameter
  PUBLIC  :: set_gau_parameter, get_gau_parameter

  ! additional subroutines/functions for finite size corrections
  PUBLIC  :: dft_has_finite_size_correction, set_finite_size_volume
  ! rpa specific
  PUBLIC  :: init_dft_exxrpa, enforce_dft_exxrpa

  ! driver subroutines computing XC
  PUBLIC  :: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
  PUBLIC  :: tau_xc , tau_xc_spin, dmxc, dmxc_spin, dmxc_nc
  PUBLIC  :: dgcxc, dgcxc_spin
  PUBLIC  :: d3gcxc       
  PUBLIC  :: nlc
  ! vector XC driver
  PUBLIC  :: evxc_t_vec, gcx_spin_vec
  !
  ! PRIVATE variables defining the DFT functional
  !
  PRIVATE :: dft, iexch, icorr, igcx, igcc, imeta, inlc
  PRIVATE :: discard_input_dft
  PRIVATE :: isgradient, ismeta, ishybrid
  PRIVATE :: exx_fraction, exx_started
  PRIVATE :: has_finite_size_correction, &
             finite_size_cell_volume,  finite_size_cell_volume_set 
  !
  character (len=25) :: dft = 'not set'
  !
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
  !              "revpbe"= "sla+pw+rpb+pbc"    = revPBE (Zhang-Yang)
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
  !              "hse"   = "sla+pw+hse+pbc"    = Heyd-Scuseria-Ernzerhof (HSE 06, see note below)
  !              "b3lyp" = "b3lp+b3lp+b3lp+b3lp"= B3LYP
  !              "b3lypv1r"    = "b3lp+b3lpv1r+b3lp+b3lp"= B3LYP-VWN1-RPA
  !              "x3lyp" = "x3lp+x3lp+x3lp+x3lp"= X3LYP
  !              "vwn-rpa"     = "sla+vwn-rpa" = VWN LDA using vwn1-rpa parametriz
  !              "gaupbe"= "sla+pw+gaup+pbc"   = Gau-PBE (also "gaup")
  !              "vdw-df"       ="sla+pw+rpb +vdw1"   = vdW-DF1
  !              "vdw-df2"      ="sla+pw+rw86+vdw2"   = vdW-DF2
  !              "vdw-df-x"     ="sla+pw+????+vdwx"   = vdW-DF-x, reserved Thonhauser, not implemented
  !              "vdw-df-y"     ="sla+pw+????+vdwy"   = vdW-DF-y, reserved Thonhauser, not implemented
  !              "vdw-df-z"     ="sla+pw+????+vdwz"   = vdW-DF-z, reserved Thonhauser, not implemented
  !              "vdw-df-c09"   ="sla+pw+c09x+vdw1"   = vdW-DF-C09
  !              "vdw-df2-c09"  ="sla+pw+c09x+vdw2"   = vdW-DF2-C09
  !              "vdw-df-cx"    ="sla+pw+cx13+vdW1"   = vdW-DF-cx
  !              "vdw-df-obk8"  ="sla+pw+obk8+vdw1"   = vdW-DF-obk8 (optB88-vdW)
  !              "vdw-df-ob86"  ="sla+pw+ob86+vdw1"   = vdW-DF-ob86 (optB86b-vdW)
  !              "vdw-df2-b86r" ="sla+pw+b86r+vdw2"   = vdW-DF2-B86R (rev-vdw-df2)
  !              "rvv10" = "sla+pw+rw86+pbc+vv10"     = rVV10
  !
  ! Any nonconflicting combination of the following keywords is acceptable:
  !
  ! Exchange:    "nox"    none                           iexch=0
  !              "sla"    Slater (alpha=2/3)             iexch=1 (default)
  !              "sl1"    Slater (alpha=1.0)             iexch=2
  !              "rxc"    Relativistic Slater            iexch=3
  !              "oep"    Optimized Effective Potential  iexch=4
  !              "hf"     Hartree-Fock                   iexch=5
  !              "pb0x"   PBE0 (Slater*0.75+HF*0.25)     iexch=6
  !              "b3lp"   B3LYP(Slater*0.80+HF*0.20)     iexch=7
  !              "kzk"    Finite-size corrections        iexch=8
  !              "x3lp"   X3LYP(Slater*0.782+HF*0.218)   iexch=9
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
  !              "rpb"    revised PBE by Zhang-Yang      igcx =4
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
  !
  ! Meta-GGA functionals
  !              "tpss"   TPSS Meta-GGA                  imeta=1
  !              "m6lx"   M06L Meta-GGA                  imeta=2
  !              "tb09"   TB09 Meta-GGA                  imeta=3
  !              "meta"   turn on MGGA                   imeta=4
  !
  ! Van der Waals functionals (nonlocal term only)
  !              "nonlc"  none                           inlc =0 (default)
  !              "vdw1"   vdW-DF1                        inlc =1
  !              "vdw2"   vdW-DF2                        inlc =2
  !              "vv10"   rVV10                          inlc =3
  !              "vdwx"   vdW-DF-x                       inlc =4, reserved Thonhauser, not implemented
  !              "vdwy"   vdW-DF-y                       inlc =5, reserved Thonhauser, not implemented
  !              "vdwz"   vdW-DF-z                       inlc =6, reserved Thonhauser, not implemented
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
  !              vdW-DF-obk8  Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)
  !              vdW-DF-ob86  Klimes et al, Phys. Rev. B, 83, 195131 (2011)
  !              c09x    V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
  !              tpss    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, 
  !                      PRL 91, 146401 (2003)
  !              tb09    F Tran and P Blaha, Phys.Rev.Lett. 102, 226401 (2009) 
  !              sogga   Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)
  !              m06l    Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)
  !              gau-pbe J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)
  !              rVV10   R. Sabatini et al. Phys. Rev. B 87, 041108(R) (2013)
  !              ev93     Engel-Vosko, Phys. Rev. B 47, 13164 (1993)
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
  ! single He atom. Info by Fabien Bruneval
  !
  integer, parameter:: notset = -1
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !    imeta: type of meta-GGA
  integer :: iexch = notset
  integer :: icorr = notset
  integer :: igcx  = notset
  integer :: igcc  = notset
  integer :: imeta = notset
  integer :: inlc  = notset
  !
  real(DP):: exx_fraction = 0.0_DP
  real(DP):: screening_parameter = 0.0_DP
  real(DP):: gau_parameter = 0.0_DP
  logical :: islda       = .false.
  logical :: isgradient  = .false.
  logical :: ismeta      = .false.
  logical :: ishybrid    = .false.
  logical :: isnonlocc   = .false.
  logical :: exx_started = .false.
  logical :: has_finite_size_correction = .false.
  logical :: finite_size_cell_volume_set = .false.
  real(DP):: finite_size_cell_volume = notset
  logical :: discard_input_dft = .false.
  !
  integer, parameter:: nxc=8, ncc=10, ngcx=27, ngcc=12, nmeta=4, ncnl=6
  character (len=4) :: exc, corr, gradx, gradc, meta, nonlocc
  dimension :: exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc), &
               meta(0:nmeta), nonlocc (0:ncnl)

  data exc  / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP', 'KZK' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' , 'KZK' /

  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX',&
               'xxxx', 'PB0X', 'B3LP','PSX', 'WCX', 'HSE', 'RW86', 'PBE', &
               'xxxx', 'C09X', 'SOX', 'xxxx', 'Q2DX', 'GAUP', 'PW86', 'B86B', &
               'OBK8', 'OB86', 'EVX', 'B86R', 'CX13' / 

  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'NONE',&
               'B3LP', 'PSC', 'PBE', 'xxxx', 'xxxx', 'Q2DC' / 

  data meta  / 'NONE', 'TPSS', 'M06L', 'TB09', 'META' / 

  data nonlocc/'NONE', 'VDW1', 'VDW2', 'VV10', 'VDWX', 'VDWY', 'VDWZ' / 

CONTAINS
  !-----------------------------------------------------------------------
  subroutine set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices iexch, icorr, igcx, igcc
    !
    implicit none
    character(len=*), intent(in) :: dft_
    integer :: len, l, i
    character (len=50):: dftout
    logical :: dft_defined = .false.
    character (len=1), external :: capital
    integer ::  save_iexch, save_icorr, save_igcx, save_igcc, save_meta, save_inlc
    !
    ! Exit if discard_input_dft
    !
    if ( discard_input_dft ) return
    !
    ! save current status of XC indices
    !
    save_iexch = iexch
    save_icorr = icorr
    save_igcx  = igcx
    save_igcc  = igcc
    save_meta  = imeta
    save_inlc  = inlc
    !
    ! convert to uppercase
    !
    len = len_trim(dft_)
    dftout = ' '
    do l = 1, len
       dftout (l:l) = capital (dft_(l:l) )
    enddo
    !
    ! ----------------------------------------------
    ! FIRST WE CHECK ALL THE SHORT NAMES
    ! Note: comparison is done via exact matching
    ! ----------------------------------------------
    !
    ! special cases : PZ  (LDA is equivalent to PZ)
    IF (('PZ' .EQ. TRIM(dftout) ).OR.('LDA' .EQ. TRIM(dftout) )) THEN
       dft_defined = set_dft_values(1,1,0,0,0,0)

    ! special cases : VWN-RPA
    else IF ('VWN-RPA' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,11,0,0,0,0)

    ! special cases : OEP no GC part (nor LDA...) and no correlation by default
    else IF ('OEP' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(4,0,0,0,0,0)

    ! special cases : HF no GC part (nor LDA...) and no correlation by default
    else IF ('HF' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(5,0,0,0,0,0)
       
    else if ('PBE' .EQ. TRIM(dftout) ) then
    ! special case : PBE
       dft_defined = set_dft_values(1,4,3,4,0,0)
       
    ! special case : BP = B88 + P86
    else if ('BP'.EQ. TRIM(dftout) ) then
       dft_defined = set_dft_values(1,1,1,1,0,0)
                     
    ! special case : PW91 = GGX + GGC
    else if ('PW91'.EQ. TRIM(dftout) ) then
       dft_defined = set_dft_values(1,4,2,2,0,0)

    elseif ( 'REVPBE' .EQ. TRIM(dftout) ) then
    ! special case : revPBE
       dft_defined = set_dft_values(1,4,4,4,0,0)

    else if ('PBESOL'.EQ. TRIM(dftout) ) then
    ! special case : PBEsol
       dft_defined = set_dft_values(1,4,10,8,0,0)

    ! special cases : BLYP (note, BLYP=>B88)
    else IF ('BLYP' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,3,1,3,0,0)
       
    else if ('OPTBK88' .EQ. TRIM(dftout)) then
    ! Special case optB88
       dft_defined = set_dft_values(1,4,23,1,0,0)
       
    else if ('OPTB86B' .EQ. TRIM(dftout)) then
    ! Special case optB86b
       dft_defined = set_dft_values(1,4,24,1,0,0)
       
    else if ('PBC'.EQ. TRIM(dftout) ) then
    ! special case : PBC  = PW + PBC 
       dft_defined = set_dft_values(1,4,0,4,0,0)
       
    ! special case : HCTH
    else if ('HCTH'.EQ. TRIM(dftout)) then
       dft_defined = set_dft_values(0,0,5,5,0,0)
       call errore('set_dft_from_name','HCTH yields suspicious results',1)
              
    ! special case : OLYP = OPTX + LYP
    else if ('OLYP'.EQ. TRIM(dftout)) then
       dft_defined = set_dft_values(0,3,6,3,0,0)
       call errore('set_dft_from_name','OLYP yields suspicious results',1)
       
    else if ('WC' .EQ. TRIM(dftout) ) then
    ! special case : Wu-Cohen
       dft_defined = set_dft_values(1,4,11,4,0,0)
       
    elseif ('PW86PBE' .EQ. TRIM(dftout) ) then
       ! special case : PW86PBE 
       dft_defined = set_dft_values(1,4,21,4,0,0)

    elseif ('B86BPBE' .EQ. TRIM(dftout) ) then
       ! special case : B86BPBE
       dft_defined = set_dft_values(1,4,22,4,0,0)

    else if ('PBEQ2D' .EQ. TRIM(dftout) .OR. 'Q2D'.EQ. TRIM(dftout) ) then
    ! special case : PBEQ2D
       dft_defined = set_dft_values(1,4,19,12,0,0)

    ! special case : SOGGA = SOX + PBEc       
    else if ('SOGGA' .EQ. TRIM(dftout) ) then
       dft_defined = set_dft_values(1,4,17,4,0,0)
    
    ! special case : Engel-Vosko
    else if ( 'EV93' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,4,25,0,0,0)

    else if ('RPBE' .EQ. TRIM(dftout) ) then
    ! special case : RPBE
         call errore('set_dft_from_name', &
     &   'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)',1)
     
    else if ('PBE0'.EQ. TRIM(dftout) ) then
    ! special case : PBE0
       dft_defined = set_dft_values(6,4,8,4,0,0)
       
   else if ('HSE' .EQ. TRIM( dftout) ) then
    ! special case : HSE
       dft_defined = set_dft_values(1,4,12,4,0,0)

   else if ( 'GAUP' .EQ. TRIM(dftout) .OR. 'GAUPBE' .EQ. TRIM(dftout) ) then
    ! special case : GAUPBE
       dft_defined = set_dft_values(1,4,20,4,0,0)
       
    else if ('VDW-DF' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF
       dft_defined = set_dft_values(1,4,4,0,1,0)

    else if ('VDW-DF-X' .EQ. TRIM(dftout) ) then
       call errore('set_dft_from_name','functional not yet implemented',1)

    else if ('VDW-DF-Y' .EQ. TRIM(dftout) ) then
       call errore('set_dft_from_name','functional not yet implemented',1)

    else if ('VDW-DF-Z' .EQ. TRIM(dftout) ) then
       call errore('set_dft_from_name','functional not yet implemented',1)

    else if ('VDW-DF-CX' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF-CX
       dft_defined = set_dft_values(1,4,27,0,1,0)

    else if ('VDW-DF-C09'  .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF with C09 exchange
       dft_defined = set_dft_values(1,4,16,0,1,0)
       
    else if ('VDW-DF-OBK8' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF-obk8, or vdW-DF + optB88
       dft_defined = set_dft_values(1,4,23,0,1,0)

    else if ('VDW-DF3' .EQ. TRIM(dftout) ) then
       call errore('set_dft_from_name','obsolete XC label, use VDW-DF-OBK8',1)

    else if ('VDW-DF-OB86' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF-ob86, or vdW-DF + optB86
       dft_defined = set_dft_values(1,4,24,0,1,0)

    else if ('VDW-DF4'.EQ.TRIM(dftout) .OR. 'OPTB86B-VDW'.EQ.TRIM(dftout) ) then
       call errore('set_dft_from_name','obsolete XC label, use VDW-DF-OB86',1)

    else if ('VDW-DF2-C09' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF2 with C09 exchange
       dft_defined = set_dft_values(1,4,16,0,2,0)

   else if ('VDW-DF2' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF2
       dft_defined = set_dft_values(1,4,13,0,2,0)

    else if ('VDW-DF2-B86R' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF2 with B86R
       dft_defined = set_dft_values(1,4,26,0,2,0)
    else if ('REV-VDW-DF2' .EQ. TRIM(dftout) ) then
       call errore('set_dft_from_name','obsolete XC label, use VDW-DF2-B86R',1)

    else if ('RVV10' .EQ. TRIM(dftout) ) then
    ! Special case rVV10
       dft_defined = set_dft_values(1,4,13,4,3,0)
       
    else if ('B3LYP'.EQ. TRIM(dftout) ) then
    ! special case : B3LYP hybrid
       dft_defined = set_dft_values(7,12,9,7,0,0)
     
    else if ('B3LYP-V1R'.EQ. TRIM(dftout) ) then
    ! special case : B3LYP-VWN-1-RPA hybrid
       dft_defined = set_dft_values(7,13,9,7,0,0)
     
    else if ('X3LYP'.EQ. TRIM(dftout) ) then
    ! special case : X3LYP hybrid
       dft_defined = set_dft_values(9,14,28,13,0,0)
     
    ! special case : TPSS meta-GGA Exc
    else IF ('TPSS'.EQ. TRIM(dftout ) ) THEN
       dft_defined = set_dft_values(1,4,7,6,0,1)

    ! special case : M06L Meta GGA
    else if ( 'M06L' .EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(0,0,0,0,0,2)

    ! special case : TB09 meta-GGA Exc
    else IF ('TB09'.EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(0,0,0,0,0,3)

    ! special case : PZ/LDA + null meta-GGA
    else IF (('PZ+META'.EQ. TRIM(dftout)) .or. ('LDA+META'.EQ. TRIM(dftout)) ) THEN
       dft_defined = set_dft_values(1,1,0,0,0,4)

    ! special case : PBE + null meta-GGA
    else IF ('PBE+META'.EQ. TRIM(dftout) ) THEN
       dft_defined = set_dft_values(1,4,3,4,0,4)

    END IF

    !
    ! ----------------------------------------------------------------
    ! If the DFT was not yet defined, check every part of the string
    ! ----------------------------------------------------------------
    !
    if (.not. dft_defined) then
    
      iexch = matching (dftout,  nxc, exc) 
      icorr = matching (dftout,  ncc, corr) 
      igcx  = matching (dftout, ngcx,gradx) 
      igcc  = matching (dftout, ngcc,gradc) 
      imeta = matching (dftout,nmeta, meta) 
      inlc  = matching (dftout, ncnl, nonlocc) 
 
    endif

    ! ----------------------------------------------------------------
    ! Last check
    ! No more defaults, the code exits if the dft is not defined
    ! ----------------------------------------------------------------

    ! Back compatibility - TO BE REMOVED

    if (igcx == 14) igcx = 3 ! PBE -> PBX
    if (igcc ==  9) igcc = 4 ! PBE -> PBC

    if (igcx == 6) &
         call infomsg('set_dft_from_name','OPTX untested! please test')

    ! check for unrecognized labels

    if ( iexch<=0.and.icorr<=0.and.igcx<=0.and.igcc<= 0.and.imeta<=0 ) then
        if ( inlc <= 0 .and. trim(dftout) /= 'NOX-NOC') then
           call errore('set_dft_from_name',trim(dftout)//': unrecognized dft',1)
        else
           ! if inlc is the only nonzero index the label is likely wrong
           call errore('set_dft_from_name',trim(dftout)//': strange dft, please check',inlc)
        endif
    endif
    !
    ! Fill variables and exit
    !
    dft = dftout

    !dft_longname = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
    !     &//gradc (igcc) //'-'// nonlocc(inlc)

    call set_auxiliary_flags
    !
    ! check dft has not been previously set differently 
    !
    if (save_iexch .ne. notset .and. save_iexch .ne. iexch) then
       write (stdout,*) iexch, save_iexch
       call errore('set_dft_from_name',' conflicting values for iexch',1)
    end if
    if (save_icorr .ne. notset .and. save_icorr .ne. icorr) then
       write (stdout,*) icorr, save_icorr
       call errore('set_dft_from_name',' conflicting values for icorr',1)
    end if
    if (save_igcx .ne. notset .and. save_igcx .ne. igcx) then
       write (stdout,*) igcx, save_igcx
       call errore('set_dft_from_name',' conflicting values for igcx',1)
    end if
    if (save_igcc .ne. notset .and. save_igcc .ne. igcc) then
       write (stdout,*) igcc, save_igcc
       call errore('set_dft_from_name',' conflicting values for igcc',1)
    end if
    if (save_meta .ne. notset .and. save_meta .ne. imeta) then
       write (stdout,*) inlc, save_meta
       call errore('set_dft_from_name',' conflicting values for imeta',1)
    end if
    if (save_inlc .ne. notset .and. save_inlc .ne. inlc) then
       write (stdout,*) inlc, save_inlc
       call errore('set_dft_from_name',' conflicting values for inlc',1)
    end if

    return
  end subroutine set_dft_from_name
  !
  integer function matching(dft, n, name)
  implicit none
  integer, intent(in):: n
  character(len=*), intent(in):: name(0:n)
  character(len=*), intent(in):: dft
  logical, external :: matches
  integer :: i
 
  matching = notset
  do i = n, 0, -1
     if (matches (name (i), trim(dft)) ) then
        if ( matching == notset ) then
           ! write(*, '("matches",i2,2X,A,2X,A)') i, name(i), trim(dft)
           matching = i 
        else
           write(*, '(2(2X,i2,2X,A))') i, trim(name(i)), &
                                   matching, trim(name(matching))
           call errore ('set_dft', 'two conflicting matching values', 1)
        end if
     endif
  end do
  if (matching == notset) matching = 0
  !
  end function matching
  !
  !-----------------------------------------------------------------------
  subroutine set_auxiliary_flags
    !-----------------------------------------------------------------------
    ! set logical flags describing the complexity of the xc functional
    ! define the fraction of exact exchange used by hybrid fuctionals
    !
    isnonlocc = (inlc > 0)
    ismeta    = (imeta > 0)
    isgradient= (igcx > 0) .or. (igcc > 0) .or. ismeta .or. isnonlocc
    islda     = (iexch> 0) .and. (icorr > 0) .and. .not. isgradient
    ! PBE0
    IF ( iexch==6 .or. igcx ==8 ) exx_fraction = 0.25_DP
    ! HSE
    IF ( igcx ==12 ) THEN
       exx_fraction = 0.25_DP
       screening_parameter = 0.106_DP
    END IF
    ! gau-pbe
    IF ( igcx ==20 ) THEN
       exx_fraction = 0.24_DP
       gau_parameter = 0.150_DP
    END IF
    ! HF or OEP
    IF ( iexch==4 .or. iexch==5 ) exx_fraction = 1.0_DP
    ! B3LYP or B3LYP-VWN-1-RPA
    IF ( iexch == 7 ) exx_fraction = 0.2_DP
    ! X3LYP
    IF ( iexch == 9 ) exx_fraction = 0.218_DP
    !
    ishybrid = ( exx_fraction /= 0.0_DP )

    has_finite_size_correction = ( iexch==8 .or. icorr==10)

    return
  end subroutine set_auxiliary_flags
  !
  !-----------------------------------------------------------------------
  logical function set_dft_values (i1,i2,i3,i4,i5,i6)
    !-----------------------------------------------------------------------
    !
    implicit none
    integer :: i1,i2,i3,i4,i5,i6
   
    iexch=i1
    icorr=i2
    igcx =i3
    igcc =i4
    inlc =i5
    imeta=i6
    set_dft_values = .true.

    return

  end function set_dft_values

  !-----------------------------------------------------------------------
  subroutine enforce_input_dft (dft_, nomsg)
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices and force any subsequent call to set_dft_from_name
    ! to return without changing them
    !
    implicit none
    character(len=*), intent(in) :: dft_
    logical, intent(in), optional :: nomsg

     call set_dft_from_name (dft_)
     if (dft == 'not set') call errore('enforce_input_dft','cannot fix unset dft',1)
     discard_input_dft = .true.

     if ( present (nomsg) ) return

     write (stdout,'(/,5x,a)') "IMPORTANT: XC functional enforced from input :"
     call write_dft_name
     write (stdout,'(5x,a)') "Any further DFT definition will be discarded"
     write (stdout,'(5x,a/)') "Please, verify this is what you really want"

     return
  end subroutine enforce_input_dft
 
  !-----------------------------------------------------------------------
  subroutine enforce_dft_exxrpa ( )
    !
    implicit none
    !
    !character(len=*), intent(in) :: dft_
    !logical, intent(in), optional :: nomsg

    iexch = 0; icorr = 0; igcx = 0; igcc = 0
    exx_fraction = 1.0_DP
    ishybrid = ( exx_fraction /= 0.0_DP )

    write (stdout,'(/,5x,a)') "XC functional enforced to be EXXRPA"
    call write_dft_name
    write (stdout,'(5x,a)') "!!! Any further DFT definition will be discarded"
    write (stdout,'(5x,a/)') "!!! Please, verify this is what you really want !"

    return
  end subroutine enforce_dft_exxrpa

  !-----------------------------------------------------------------------
  subroutine init_dft_exxrpa ( )
    !
    implicit none
    !
    exx_fraction = 1.0_DP
    ishybrid = ( exx_fraction /= 0.0_DP )

    write (stdout,'(/,5x,a)') "Only exx_fraction is set to 1.d0"
    write (stdout,'(5x,a)') "XC functional still not changed"
    call write_dft_name

    return
  end subroutine init_dft_exxrpa

  !-----------------------------------------------------------------------
  subroutine start_exx 
     if (.not. ishybrid) &
        call errore('start_exx','dft is not hybrid, wrong call',1)
     exx_started = .true.
  end subroutine start_exx
  !-----------------------------------------------------------------------
  subroutine stop_exx 
     if (.not. ishybrid) &
        call errore('stop_exx','dft is not hybrid, wrong call',1)
     exx_started = .false.
  end subroutine stop_exx
  !-----------------------------------------------------------------------
  function exx_is_active ()
     logical exx_is_active
     exx_is_active = exx_started
  end function exx_is_active
  !-----------------------------------------------------------------------
  subroutine set_exx_fraction (exxf_)
     implicit none
     real(DP):: exxf_
     exx_fraction = exxf_
     write (stdout,'(5x,a,f6.2)') 'EXX fraction changed: ',exx_fraction
  end subroutine set_exx_fraction
  !---------------------------------------------------------------------
  subroutine set_screening_parameter (scrparm_)
     implicit none
     real(DP):: scrparm_
     screening_parameter = scrparm_
     write (stdout,'(5x,a,f12.7)') 'EXX Screening parameter changed: ', &
          & screening_parameter
  end subroutine set_screening_parameter
  !----------------------------------------------------------------------
  function get_screening_parameter ()
     real(DP):: get_screening_parameter
     get_screening_parameter = screening_parameter
     return
  end function get_screening_parameter
  !---------------------------------------------------------------------
  subroutine set_gau_parameter (gauparm_)
     implicit none
     real(DP):: gauparm_
     gau_parameter = gauparm_
     write (stdout,'(5x,a,f12.7)') 'EXX Gau parameter changed: ', &
          & gau_parameter
  end subroutine set_gau_parameter
  !----------------------------------------------------------------------
  function get_gau_parameter ()
     real(DP):: get_gau_parameter
     get_gau_parameter = gau_parameter
     return
  end function get_gau_parameter
  !-----------------------------------------------------------------------
  function get_iexch ()
     integer get_iexch
     get_iexch = iexch
     return
  end function get_iexch
  !-----------------------------------------------------------------------
  function get_icorr ()
     integer get_icorr
     get_icorr = icorr
     return
  end function get_icorr
  !-----------------------------------------------------------------------
  function get_igcx ()
     integer get_igcx
     get_igcx = igcx
     return
  end function get_igcx
  !-----------------------------------------------------------------------
  function get_igcc ()
     integer get_igcc
     get_igcc = igcc
     return
  end function get_igcc
  !-----------------------------------------------------------------------
  function get_meta ()
     integer get_meta
     get_meta = imeta
     return
  end function get_meta
  !-----------------------------------------------------------------------
  function get_inlc ()
     integer get_inlc
     get_inlc = inlc
     return
  end function get_inlc
  !-----------------------------------------------------------------------
  function get_nonlocc_name ()
     character(10) get_nonlocc_name
     get_nonlocc_name = TRIM(nonlocc(inlc))
     return
  end function get_nonlocc_name
 !-----------------------------------------------------------------------
  function dft_is_nonlocc ()
    logical :: dft_is_nonlocc
    dft_is_nonlocc = isnonlocc
    return
  end function dft_is_nonlocc
  !-----------------------------------------------------------------------
  function get_exx_fraction ()
     real(DP):: get_exx_fraction
     get_exx_fraction = exx_fraction
     return
  end function get_exx_fraction
  !-----------------------------------------------------------------------
  function get_dft_name ()
     character (len=25) :: get_dft_name
     get_dft_name = dft
     return
  end function get_dft_name
  !-----------------------------------------------------------------------
  function dft_is_gradient ()
     logical :: dft_is_gradient
     dft_is_gradient = isgradient
     return
  end function dft_is_gradient
  !-----------------------------------------------------------------------
  function dft_is_meta ()
     logical :: dft_is_meta
     dft_is_meta = ismeta
     return
  end function dft_is_meta
  !-----------------------------------------------------------------------
  function dft_is_hybrid ()
     logical :: dft_is_hybrid
     dft_is_hybrid = ishybrid
     return
  end function dft_is_hybrid
  !-----------------------------------------------------------------------
  function igcc_is_lyp ()
     logical :: igcc_is_lyp
     igcc_is_lyp = (get_igcc() == 3 .or. get_igcc() == 7)
     return
  end function igcc_is_lyp
  !-----------------------------------------------------------------------
  function dft_has_finite_size_correction ()
     logical :: dft_has_finite_size_correction
     dft_has_finite_size_correction = has_finite_size_correction
     return
  end function  dft_has_finite_size_correction
  !-----------------------------------------------------------------------
  subroutine set_finite_size_volume(volume)
     real, intent (IN) :: volume
     if (.not. has_finite_size_correction) &
         call errore('set_finite_size_volume', &
                      'dft w/o finite_size_correction, wrong call',1)
     if (volume <= 0.d0) &
         call errore('set_finite_size_volume', &
                     'volume is not positive, check omega and/or nk1,nk2,nk3',1)
     finite_size_cell_volume = volume
     finite_size_cell_volume_set = .TRUE.
  end subroutine set_finite_size_volume
  !-----------------------------------------------------------------------
  
  !-----------------------------------------------------------------------
  subroutine set_dft_from_indices(iexch_,icorr_,igcx_,igcc_, inlc_)
     integer :: iexch_, icorr_, igcx_, igcc_, inlc_
     if ( discard_input_dft ) return
     if (iexch == notset) iexch = iexch_
     if (iexch /= iexch_) then
        write (stdout,*) iexch, iexch_
        call errore('set_dft',' conflicting values for iexch',1)
     end if
     if (icorr == notset) icorr = icorr_
     if (icorr /= icorr_) then
        write (stdout,*) icorr, icorr_
        call errore('set_dft',' conflicting values for icorr',1)
     end if
     if (igcx  == notset) igcx = igcx_
     if (igcx /= igcx_) then
        write (stdout,*) igcx, igcx_
        call errore('set_dft',' conflicting values for igcx',1)
     end if
     if (igcc  == notset) igcc = igcc_
     if (igcc /= igcc_) then
        write (stdout,*) igcc, igcc_
        call errore('set_dft',' conflicting values for igcc',1)
     end if
     if (inlc  == notset) inlc = inlc_
     if (inlc /= inlc_) then
        write (stdout,*) inlc, inlc_
        call errore('set_dft',' conflicting values for inlc',1)
     end if
     dft = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
           &//gradc (igcc)//'-'//nonlocc (inlc)
     ! WRITE( stdout,'(a)') dft
     call set_auxiliary_flags
     return
  end subroutine set_dft_from_indices
  !------------------------------------------------------------------------------- 
  function get_dft_short ( )
  !---------------------------------------------------------------------
  !
  character (len=12) :: get_dft_short
  character (len=12) :: shortname
  !
  shortname = 'no shortname'
  if (iexch==1.and.igcx==0.and.igcc==0) then
     shortname = corr(icorr)
  else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
     shortname = 'BLYP'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) then
     shortname = 'B88'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) then
     shortname = 'BP'
  else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
     shortname = 'PW91'
  else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
     shortname = 'PBE'
  else if (iexch==6.and.icorr==4.and.igcx==8.and.igcc==4) then
     shortname = 'PBE0'
  else if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==4) then
     shortname = 'revPBE'
  else if (iexch==1.and.icorr==4.and.igcx==10.and.igcc==8) then
     shortname = 'PBESOL'
  else if (iexch==1.and.icorr==4.and.igcx==19.and.igcc==12) then
     shortname = 'Q2D'
  else if (iexch==1.and.icorr==4.and.igcx==12.and.igcc==4) then
     shortname = 'HSE'
  else if (iexch==1.and.icorr==4.and.igcx==20.and.igcc==4) then
     shortname = 'GAUPBE'
  else if (iexch==1.and.icorr==4.and.igcx==11.and.igcc==4) then
     shortname = 'WC'
  else if (iexch==7.and.icorr==12.and.igcx==9.and. igcc==7) then
     shortname = 'B3LYP'
  else if (iexch==7.and.icorr==13.and.igcx==9.and. igcc==7) then
     shortname = 'B3LYP-V1R'
  else if (iexch==9.and.icorr==14.and.igcx==28.and. igcc==13) then
     shortname = 'X3LYP'
  else if (iexch==0.and.icorr==3.and.igcx==6.and.igcc==3) then
     shortname = 'OLYP'
  else if (iexch==1.and.icorr==4.and.igcx==17.and.igcc==4) then
     shortname = 'SOGGA'
  else if (iexch==1.and.icorr==4.and.igcx==25.and.igcc==0) then
     shortname = 'EV93'
  end if

  if (imeta == 1 ) then
     shortname = 'TPSS'
  else if (imeta == 2) then
     shortname = 'M06L'
  else if (imeta == 3) then
     shortname = 'TB09'
  else if (imeta == 4) then
     shortname = 'META'
  end if

  if ( inlc==1 ) then
     if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==0) then
        shortname = 'VDW-DF'
     else if (iexch==1.and.icorr==4.and.igcx==27.and.igcc==0) then
        shortname = 'VDW-DF-CX'
     else if (iexch==1.and.icorr==4.and.igcx==16.and.igcc==0) then
        shortname = 'VDW-DF-C09'
     else if (iexch==1.and.icorr==4.and.igcx==24.and.igcc==0) then
        shortname = 'VDW-DF-OB86'
     else if (iexch==1.and.icorr==4.and.igcx==23.and.igcc==0) then
        shortname = 'VDW-DF-OBK8'
     end if
  else if ( inlc==2 ) then
     if (iexch==1.and.icorr==4.and.igcx==13.and.igcc==0) then
        shortname = 'VDW-DF2'
     else if (iexch==1.and.icorr==4.and.igcx==16.and.igcc==0) then
        shortname = 'VDW-DF2-C09'
     else if (iexch==1.and.icorr==4.and.igcx==26.and.igcc==0) then
        shortname = 'VDW-DF2-B86R'
     end if
  else if ( inlc==3) then
     shortname = 'RVV10'
  end if
  !
  get_dft_short = shortname
  !
  end function get_dft_short
  !
  !---------------------------------------------------------------------
  function get_dft_long( )
  !---------------------------------------------------------------------
  !
  character (len=25) :: get_dft_long
  character (len=25) :: longname
  !
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)
  if ( imeta > 0 ) then
     longname = longname(1:20)//trim(meta(imeta))
  else if ( inlc > 0 ) then
     longname = longname(1:20)//trim(nonlocc(inlc))
  end if
  get_dft_long = longname

end function get_dft_long

!-----------------------------------------------------------------------
subroutine write_dft_name
!-----------------------------------------------------------------------
   WRITE( stdout, '(5X,"Exchange-correlation      = ",A, &
        &  " (",I2,3I3,2I2,")")') TRIM( dft ), iexch,icorr,igcx,igcc,inlc,imeta
   IF ( get_exx_fraction() > 0.0_dp ) WRITE( stdout, &
        '(5X,"EXX-fraction              =",F12.2)') get_exx_fraction()
   return
end subroutine write_dft_name

!
!-----------------------------------------------------------------------
!-------  LDA DRIVERS --------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine xc (rho, ex, ec, vx, vc)
  !-----------------------------------------------------------------------
  !     lda exchange and correlation functionals - Hartree a.u.
  !
  !     exchange   :  Slater, relativistic Slater
  !     correlation:  Ceperley-Alder (Perdew-Zunger parameters)
  !                   Vosko-Wilk-Nusair
  !                   Lee-Yang-Parr
  !                   Perdew-Wang
  !                   Wigner
  !                   Hedin-Lundqvist
  !                   Ortiz-Ballone (Perdew-Zunger formula)
  !                   Ortiz-Ballone (Perdew-Wang formula)
  !                   Gunnarsson-Lundqvist
  !
  !     input : rho=rho(r)
  !     definitions: E_x = \int E_x(rho) dr, E_x(rho) = rho\epsilon_c(rho)
  !                  same for correlation
  !     output: ex = \epsilon_x(rho) ( NOT E_x(rho) )
  !             vx = dE_x(rho)/drho  ( NOT d\epsilon_x(rho)/drho )
  !             ec, vc as above for correlation
  !
  implicit none

  real(DP) :: rho, ec, vc, ex, vx
  real(DP) :: ec__, vc__
  !
  real(DP), parameter :: small = 1.E-10_DP,  third = 1.0_DP / 3.0_DP, &
       pi34 = 0.6203504908994_DP  ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0_DP
     vc = 0.0_DP
     ex = 0.0_DP
     vx = 0.0_DP
     return
  else
     rs = pi34 / rho**third
     ! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  endif
  !..exchange
  if (iexch == 1) THEN             !  'sla'
     call slater (rs, ex, vx)
  ELSEIF (iexch == 2) THEN         !  'sl1'
     call slater1(rs, ex, vx)
  ELSEIF (iexch == 3) THEN         !  'rxc'
     CALL slater_rxc(rs, ex, vx)
  ELSEIF ((iexch == 4).or.(iexch==5)) THEN  ! 'oep','hf'
     IF (exx_started) then
        ex = 0.0_DP
        vx = 0.0_DP
     else
        call slater (rs, ex, vx)
     endif
  ELSEIF (iexch == 6) THEN         !  'pb0x'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = (1.0_DP - exx_fraction) * ex 
        vx = (1.0_DP - exx_fraction) * vx 
     end if
  ELSEIF (iexch == 7) THEN         !  'B3LYP'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = 0.8_DP * ex 
        vx = 0.8_DP * vx 
     end if
  ELSEIF (iexch == 8) THEN         !  'sla+kzk'
     if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
          'finite size corrected exchange used w/o initialization',1)
     call slaterKZK (rs, ex, vx, finite_size_cell_volume)
     !
  ELSEIF (iexch == 9) THEN         !  'X3LYP'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = 0.782_DP * ex 
        vx = 0.782_DP * vx 
     end if
  else
     ex = 0.0_DP
     vx = 0.0_DP
  endif
  !..correlation
  if (icorr == 1) then
     call pz (rs, 1, ec, vc)
  elseif (icorr == 2) then
     call vwn (rs, ec, vc)
  elseif (icorr == 3) then
     call lyp (rs, ec, vc)
  elseif (icorr == 4) then
     call pw (rs, 1, ec, vc)
  elseif (icorr == 5) then
     call wigner (rs, ec, vc)
  elseif (icorr == 6) then
     call hl (rs, ec, vc)
  elseif (icorr == 7) then
     call pz (rs, 2, ec, vc)
  elseif (icorr == 8) then
     call pw (rs, 2, ec, vc)
  elseif (icorr == 9) then
     call gl (rs, ec, vc)
  elseif (icorr ==10) then
     if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
          'finite size corrected correlation used w/o initialization',1)
     call pzKZK (rs, ec, vc, finite_size_cell_volume)
  elseif (icorr ==11) then
     call vwn1_rpa (rs, ec, vc)
  elseif (icorr ==12) then  ! 'B3LYP'
     call vwn (rs, ec, vc)
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     call lyp( rs, ec__, vc__ )
     ec = ec + 0.81_DP * ec__
     vc = vc + 0.81_DP * vc__
  elseif (icorr ==13) then  ! 'B3LYP-V1R'
     call vwn1_rpa (rs, ec, vc)
     ec = 0.19_DP * ec
     vc = 0.19_DP * vc
     call lyp( rs, ec__, vc__ )
     ec = ec + 0.81_DP * ec__
     vc = vc + 0.81_DP * vc__
  elseif (icorr ==14) then  ! 'X3LYP'
     call vwn1_rpa (rs, ec, vc)
     ec = 0.129_DP * ec
     vc = 0.129_DP * vc
     call lyp( rs, ec__, vc__ )
     ec = ec + 0.871_DP * ec__
     vc = vc + 0.871_DP * vc__
  else
     ec = 0.0_DP
     vc = 0.0_DP
  endif
  !
  return
end subroutine xc
!!!!!!!!!!!!!!SPIN
!-----------------------------------------------------------------------
subroutine xc_spin (rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  implicit none

  real(DP) :: rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw
  real(DP) :: ec__, vcup__, vcdw__
  !
  real(DP), parameter :: small= 1.E-10_DP, third = 1.0_DP/3.0_DP, &
       pi34= 0.6203504908994_DP ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0_DP
     vcup = 0.0_DP
     vcdw = 0.0_DP
     ex = 0.0_DP
     vxup = 0.0_DP
     vxdw = 0.0_DP
     return
  else
     rs = pi34 / rho**third
  endif
  !..exchange
  IF (iexch == 1) THEN      ! 'sla'
     call slater_spin (rho, zeta, ex, vxup, vxdw)
  ELSEIF (iexch == 2) THEN  ! 'sl1'
     call slater1_spin (rho, zeta, ex, vxup, vxdw)
  ELSEIF (iexch == 3) THEN  ! 'rxc'
     call slater_rxc_spin ( rho, zeta, ex, vxup, vxdw )
  ELSEIF ((iexch == 4).or.(iexch==5)) THEN  ! 'oep','hf'
     IF (exx_started) then
        ex   = 0.0_DP
        vxup = 0.0_DP 
        vxdw = 0.0_DP 
     else
        call slater_spin (rho, zeta, ex, vxup, vxdw)
     endif
  ELSEIF (iexch == 6) THEN  ! 'pb0x'
     call slater_spin (rho, zeta, ex, vxup, vxdw)
     if (exx_started) then
        ex   = (1.0_DP - exx_fraction) * ex
        vxup = (1.0_DP - exx_fraction) * vxup 
        vxdw = (1.0_DP - exx_fraction) * vxdw 
     end if
  ELSEIF (iexch == 7) THEN  ! 'B3LYP'
     call slater_spin (rho, zeta, ex, vxup, vxdw)
     if (exx_started) then
        ex   = 0.8_DP * ex
        vxup = 0.8_DP * vxup 
        vxdw = 0.8_DP * vxdw 
     end if
  ELSE
     ex = 0.0_DP
     vxup = 0.0_DP
     vxdw = 0.0_DP
  ENDIF
  !..correlation
  if (icorr == 0) then
     ec = 0.0_DP
     vcup = 0.0_DP
     vcdw = 0.0_DP
  elseif (icorr == 1) then
     call pz_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 2) then
     call vwn_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 3) then
     call lsd_lyp (rho, zeta, ec, vcup, vcdw) ! from CP/FPMD (more_functionals)
  elseif (icorr == 4) then
     call pw_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 12) then ! 'B3LYP'
     call vwn_spin (rs, zeta, ec, vcup, vcdw)
     ec = 0.19_DP * ec
     vcup = 0.19_DP * vcup
     vcdw = 0.19_DP * vcdw
     call lsd_lyp (rho, zeta, ec__, vcup__, vcdw__) ! from CP/FPMD (more_functionals)
     ec = ec + 0.81_DP * ec__
     vcup = vcup + 0.81_DP * vcup__
     vcdw = vcdw + 0.81_DP * vcdw__
  elseif (icorr == 13) then   ! 'B3LYP-V1R'
     call vwn1_rpa_spin (rs, zeta, ec, vcup, vcdw)
     ec = 0.19_DP * ec
     vcup = 0.19_DP * vcup
     vcdw = 0.19_DP * vcdw
     call lsd_lyp (rho, zeta, ec__, vcup__, vcdw__) ! from CP/FPMD (more_functionals)
     ec = ec + 0.81_DP * ec__
     vcup = vcup + 0.81_DP * vcup__
     vcdw = vcdw + 0.81_DP * vcdw__
  else
     call errore ('lsda_functional (xc_spin)', 'not implemented', icorr)
  endif
  !
  return
end subroutine xc_spin
!
!-----------------------------------------------------------------------
subroutine xc_spin_vec (rho, zeta, length, evx, evc)
  !-----------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  implicit none

  integer, intent(in)   :: length
  real(DP), intent(in)  :: rho(length), zeta(length)
  real(DP), intent(out) :: evx(length,3), evc(length,3)
  !
  real(DP), parameter :: small= 1.E-10_DP, third = 1.0_DP/3.0_DP, &
       pi34= 0.6203504908994_DP ! pi34=(3/4pi)^(1/3)
  !
  integer  :: i
  logical  :: comp_energy_loc
  real(DP) :: rs(length)
  !
  !..exchange
  select case (iexch)
  case(1)            ! 'sla'
     call slater_spin_vec (rho, zeta, evx, length)
  case(2)            ! 'sl1'
     do i=1,length
        call slater1_spin (rho(i), zeta(i), evx(i,3), evx(i,1), evx(i,2))
     end do
  case(3)            ! 'rxc'
     do i=1,length
        call slater_rxc_spin (rho(i), zeta(i), evx(i,3), evx(i,1), evx(i,2))
     end do
  case(4,5)          ! 'oep','hf'
     if (exx_started) then
        evx = 0.0_DP
     else
        call slater_spin_vec (rho, zeta, evx, length)
     endif
  case(6)            ! 'pb0x'
     call slater_spin_vec (rho, zeta, evx, length)
     if (exx_started) then
        evx = (1.0_DP - exx_fraction) * evx
     end if
  case(7)            ! 'B3LYP'
     call slater_spin_vec (rho, zeta, evx, length)
     if (exx_started) then
        evx = 0.8_DP * evx
     end if
  case default
     evx = 0.0_DP
  end select

  !..correlation
  where (rho.gt.small)
     rs = pi34 / rho**third
  elsewhere
     rs = 1.0_DP ! just a sane default, results are discarded anyway
  end where

  select case(icorr)
  case (0)
     evc = 0.0_DP
  case (1)
     do i=1,length
        call pz_spin (rs(i), zeta(i), evc(i,3), evc(i,1), evc(i,2))
     end do
  case (2)
     do i=1,length
        call vwn_spin (rs(i), zeta(i), evc(i,3), evc(i,1), evc(i,2))
     end do
  case(3)
     do i=1,length
        call lsd_lyp (rho(i), zeta(i), evc(i,3), evc(i,1), evc(i,2)) ! from CP/FPMD (more_functionals)
     end do
  case(4)
     call pw_spin_vec (rs, zeta, evc, length)
  case default
     call errore ('lsda_functional (xc_spin_vec)', 'not implemented', icorr)
  end select
  !
  where (rho.le.small)
     evx(:,1) = 0.0_DP
     evc(:,1) = 0.0_DP

     evx(:,2) = 0.0_DP
     evc(:,2) = 0.0_DP

     evx(:,3) = 0.0_DP
     evc(:,3) = 0.0_DP
  end where
  !
end subroutine xc_spin_vec
!
!-----------------------------------------------------------------------
!------- GRADIENT CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine gcxc (rho, grho, sx, sc, v1x, v2x, v1c, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     See comments at the beginning of module for implemented cases
  !
  !     input:  rho, grho=|\nabla rho|^2
  !     definition:  E_x = \int E_x(rho,grho) dr
  !     output: sx = E_x(rho,grho)
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             sc, v1c, v2c as above for correlation
  !
  implicit none

  real(DP) :: rho, grho, sx, sc, v1x, v2x, v1c, v2c
  real(DP) :: sx__,v1x__, v2x__
  real(DP) :: sxsr, v1xsr, v2xsr
  real(DP), parameter:: small = 1.E-10_DP

  ! exchange
  if (rho <= small) then
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  elseif (igcx == 1) then
     call becke88 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 2) then
     call ggax (rho, grho, sx, v1x, v2x)
  elseif (igcx == 3) then
     call pbex (rho, grho, 1, sx, v1x, v2x)
  elseif (igcx == 4) then
     call pbex (rho, grho, 2, sx, v1x, v2x)
  elseif (igcx == 5 .and. igcc == 5) then
     call hcth(rho, grho, sx, v1x, v2x)
  elseif (igcx == 6) then
     call optx (rho, grho, sx, v1x, v2x)
  ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
  ! routine: needs kinetic energy density in addition to rho and grad rho
  elseif (igcx == 8) then ! 'PBE0'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if (exx_started) then
        sx  = (1.0_DP - exx_fraction) * sx
        v1x = (1.0_DP - exx_fraction) * v1x
        v2x = (1.0_DP - exx_fraction) * v2x
     end if
  elseif (igcx == 9) then ! 'B3LYP'
     call becke88 (rho, grho, sx, v1x, v2x)
     if (exx_started) then
        sx  = 0.72_DP * sx
        v1x = 0.72_DP * v1x
        v2x = 0.72_DP * v2x
     end if
  elseif (igcx ==10) then ! 'pbesol'
     call pbex (rho, grho, 3, sx, v1x, v2x)
  elseif (igcx ==11) then ! 'wc'
     call wcx (rho, grho, sx, v1x, v2x)
  elseif (igcx ==12) then ! 'pbexsr'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if(exx_started) then
       call pbexsr (rho, grho, sxsr, v1xsr, v2xsr, screening_parameter)
       sx = sx - exx_fraction * sxsr
       v1x = v1x - exx_fraction * v1xsr
       v2x = v2x - exx_fraction * v2xsr
     endif 
  elseif (igcx ==13) then ! 'rPW86'
     call rPW86 (rho, grho, sx, v1x, v2x)
  elseif (igcx ==16) then ! 'C09x'
     call c09x (rho, grho, sx, v1x, v2x)
  elseif (igcx ==17) then ! 'sogga'
     call sogga(rho, grho, sx, v1x, v2x)
  elseif (igcx ==19) then ! 'pbeq2d'
     call pbex (rho, grho, 4, sx, v1x, v2x)
  elseif (igcx ==20) then ! 'gau-pbe'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if(exx_started) then
       call pbexgau (rho, grho, sxsr, v1xsr, v2xsr, gau_parameter)
       sx = sx - exx_fraction * sxsr
       v1x = v1x - exx_fraction * v1xsr
       v2x = v2x - exx_fraction * v2xsr
     endif
  elseif (igcx == 21) then ! 'pw86'
     call pw86 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 22) then ! 'b86b'
     call becke86b (rho, grho, sx, v1x, v2x)
     ! call b86b (rho, grho, 1, sx, v1x, v2x)
  elseif (igcx == 23) then ! 'optB88'
     call pbex (rho, grho, 5, sx, v1x, v2x)
  elseif (igcx == 24) then ! 'optB86b'
     call pbex (rho, grho, 6, sx, v1x, v2x)
     ! call b86b (rho, grho, 2, sx, v1x, v2x)
  elseif (igcx == 25) then ! 'ev93'
     call pbex (rho, grho, 7, sx, v1x, v2x)
  elseif (igcx == 26) then ! 'b86r'
     call b86b (rho, grho, 3, sx, v1x, v2x)
  elseif (igcx == 27) then ! 'cx13'
     call cx13 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 28) then ! 'X3LYP'
     call becke88 (rho, grho, sx, v1x, v2x)
     call pbex (rho, grho, 1, sx__, v1x__, v2x__)
     if (exx_started) then
        sx  = real(0.765*0.709,DP) * sx
        v1x = real(0.765*0.709,DP) * v1x
        v2x = real(0.765*0.709,DP) * v2x
        sx  = sx  + real(0.235*0.709) * sx__
        v1x = v1x + real(0.235*0.709) * v1x__
        v2x = v2x + real(0.235*0.709) * v2x__
     end if
  else
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  endif
  ! correlation
  if (rho.le.small) then
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  elseif (igcc == 1) then
     call perdew86 (rho, grho, sc, v1c, v2c)
  elseif (igcc == 2) then
     call ggac (rho, grho, sc, v1c, v2c)
  elseif (igcc == 3) then
     call glyp (rho, grho, sc, v1c, v2c)
  elseif (igcc == 4) then
     call pbec (rho, grho, 1, sc, v1c, v2c)
  ! igcc == 5 (HCTH) is calculated together with case igcx=5
  ! igcc == 6 (meta-GGA) is treated in a different routine
  elseif (igcc == 7) then !'B3LYP'
     call glyp (rho, grho, sc, v1c, v2c)
     if (exx_started) then
        sc  = 0.81_DP * sc
        v1c = 0.81_DP * v1c
        v2c = 0.81_DP * v2c
     end if
  elseif (igcc == 8) then ! 'PBEsol'
     call pbec (rho, grho, 2, sc, v1c, v2c)
  ! igcc == 9 set to 5, back-compatibility
  ! igcc ==10 set to 6, back-compatibility
  ! igcc ==11 M06L calculated in another routine
  else if (igcc == 12) then ! 'Q2D'
     call pbec (rho, grho, 3, sc, v1c, v2c)
  elseif (igcc == 13) then !'X3LYP'
     call glyp (rho, grho, sc, v1c, v2c)
     if (exx_started) then
        sc  = 0.871_DP * sc
        v1c = 0.871_DP * v1c
        v2c = 0.871_DP * v2c
     end if
  else
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  endif
  !
  return
end subroutine gcxc
!
!!!!!!!!!!!!!!SPIN
!-----------------------------------------------------------------------
subroutine gcx_spin (rhoup, rhodw, grhoup2, grhodw2, &
                     sx, v1xup, v1xdw, v2xup, v2xdw)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !
  implicit none
  !
  !     dummy arguments
  !
  real(DP) :: rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, v1xdw, &
       v2xup, v2xdw
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  real(DP) :: sxsr, v1xupsr, v2xupsr, v1xdwsr, v2xdwsr
  real(DP), parameter :: small = 1.E-10_DP
  real(DP) :: rho, sxup, sxdw
  integer :: iflag
  !
  !
  ! exchange
  rho = rhoup + rhodw
  if (rho <= small .or. igcx == 0) then
     sx = 0.0_DP
     v1xup = 0.0_DP
     v2xup = 0.0_DP
     v1xdw = 0.0_DP
     v2xdw = 0.0_DP
  elseif (igcx == 1) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = sxup + sxdw
  elseif (igcx == 2) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call ggax (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call ggax (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
  elseif (igcx == 3 .or. igcx == 4 .or. igcx == 8 .or. igcx ==10 .or. &
          igcx ==12 .or. igcx ==20 .or. igcx ==23 .or. igcx ==24 .or. igcx == 25) then
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8: PBE0, igcx=10: PBEsol
     ! igcx=12: HSE,  igcx=20: gau-pbe, igcx=23: obk8, igcx=24: ob86, igcx=25: ev93
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
     elseif (igcx == 23) then
        iflag = 5
     elseif (igcx == 24) then
        iflag = 6
     elseif (igcx == 25) then
        iflag = 7
     else
        iflag = 1
     endif
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pbex (2.0_DP * rhoup, 4.0_DP * grhoup2, iflag, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pbex (2.0_DP * rhodw, 4.0_DP * grhodw2, iflag, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
     if (igcx == 8 .and. exx_started ) then
       sx = (1.0_DP - exx_fraction) * sx
       v1xup = (1.0_DP - exx_fraction) * v1xup
       v1xdw = (1.0_DP - exx_fraction) * v1xdw
       v2xup = (1.0_DP - exx_fraction) * v2xup
       v2xdw = (1.0_DP - exx_fraction) * v2xdw
     end if
     if (igcx == 12 .and. exx_started ) then

        call pbexsr_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
                         v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                         screening_parameter)
        sx  = sx - exx_fraction*sxsr
        v1xup = v1xup - exx_fraction*v1xupsr
        v2xup = v2xup - exx_fraction*v2xupsr
        v1xdw = v1xdw - exx_fraction*v1xdwsr
        v2xdw = v2xdw - exx_fraction*v2xdwsr
     end if

     if (igcx == 20 .and. exx_started ) then
        ! gau-pbe
        call pbexgau_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
                         v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                         gau_parameter)
        sx  = sx - exx_fraction*sxsr
        v1xup = v1xup - exx_fraction*v1xupsr
        v2xup = v2xup - exx_fraction*v2xupsr
        v1xdw = v1xdw - exx_fraction*v1xdwsr
        v2xdw = v2xdw - exx_fraction*v2xdwsr
     end if

  elseif (igcx == 9) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = sxup + sxdw

     if (exx_started ) then
       sx = 0.72_DP * sx
       v1xup = 0.72_DP * v1xup
       v1xdw = 0.72_DP * v1xdw
       v2xup = 0.72_DP * v2xup
       v2xdw = 0.72_DP * v2xdw
     end if

  elseif (igcx == 11) then ! 'Wu-Cohen'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call wcx (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call wcx (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 13) then ! 'revised PW86 for vdw-df2'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call rPW86 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call rPW86 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 16) then ! 'c09x for vdw-df-c09.'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call c09x (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call c09x (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 21) then ! 'PW86'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pw86 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pw86 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 22) then ! 'B86B'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke86b (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke86b (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

   elseif (igcx == 26) then ! 'B86R for rev-vdW-DF2'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call b86b (2.0_DP * rhoup, 4.0_DP * grhoup2, 3, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call b86b (2.0_DP * rhodw, 4.0_DP * grhodw2, 3, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  elseif (igcx == 27) then ! 'cx13 for vdw-df-cx'
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call cx13 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call cx13 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw


  ! case igcx == 5 (HCTH) and 6 (OPTX) not implemented
  ! case igcx == 7 (meta-GGA) must be treated in a separate call to another
  ! routine: needs kinetic energy density in addition to rho and grad rho

  else
     call errore ('gcx_spin', 'not implemented', igcx)
  endif
  !
  return
end subroutine gcx_spin
!
!-----------------------------------------------------------------------
subroutine gcx_spin_vec(rhoup, rhodw, grhoup2, grhodw2, &
     sx, v1xup, v1xdw, v2xup, v2xdw, length)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !
  implicit none
  !
  !     dummy arguments
  !
  integer, intent(in) :: length
  real(DP),intent(in) :: rhoup(length), rhodw(length)
  real(DP),intent(in) :: grhoup2(length), grhodw2(length)
  real(DP),intent(out) :: sx(length)
  real(DP),intent(out) :: v1xup(length), v1xdw(length)
  real(DP),intent(out) :: v2xup(length), v2xdw(length)
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  real(DP), parameter :: small = 1.E-10_DP
  real(DP) :: rho(length), sxup(length), sxdw(length)
  integer :: iflag
  integer :: i
  ! only used for HSE (igcx == 12):
  real(DP) :: sxsr, v1xupsr, v2xupsr, v1xdwsr, v2xdwsr
  !
  !
  ! exchange
  rho = rhoup + rhodw
  select case(igcx)
  case(0)
     sx = 0.0_DP
     v1xup = 0.0_DP
     v2xup = 0.0_DP
     v1xdw = 0.0_DP
     v2xdw = 0.0_DP
  case(1)
     do i=1,length
        if (rhoup(i) > small .and. sqrt (abs (grhoup2(i)) ) > small) then
           call becke88_spin (rhoup(i), grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt (abs (grhodw2(i)) ) > small) then
           call becke88_spin (rhodw(i), grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = sxup + sxdw
  case(2)
     do i=1,length
        if (rhoup(i) > small .and. sqrt (abs (grhoup2(i)) ) > small) then
           call ggax (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt (abs (grhodw2(i)) ) > small) then
           call ggax (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
  case(3,4,8,10,12,25)
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8 PBE0, igcx=10: PBEsol, 
     ! igcx=25: EV93
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
     elseif (igcx == 25) then
        iflag = 7
     else
        iflag = 1
     endif

     call pbex_vec (2.0_DP * rhoup, 4.0_DP * grhoup2, iflag, sxup, v1xup, v2xup, length, small)
     call pbex_vec (2.0_DP * rhodw, 4.0_DP * grhodw2, iflag, sxdw, v1xdw, v2xdw, length, small)
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
     if (igcx == 8 .and. exx_started ) then
        sx = (1.0_DP - exx_fraction) * sx
        v1xup = (1.0_DP - exx_fraction) * v1xup
        v1xdw = (1.0_DP - exx_fraction) * v1xdw
        v2xup = (1.0_DP - exx_fraction) * v2xup
        v2xdw = (1.0_DP - exx_fraction) * v2xdw
     end if
     if (igcx == 12 .and. exx_started ) then
        ! in this case the subroutine is not really "vector"
        DO i = 1, length
          call pbexsr_lsd (rhoup(i), rhodw(i), grhoup2(i), grhodw2(i), sxsr,  &
                           v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                           screening_parameter)
        sx(i)  = sx(i) - exx_fraction*sxsr
        v1xup(i) = v1xup(i) - exx_fraction*v1xupsr
        v2xup(i) = v2xup(i) - exx_fraction*v2xupsr
        v1xdw(i) = v1xdw(i) - exx_fraction*v1xdwsr
        v2xdw(i) = v2xdw(i) - exx_fraction*v2xdwsr
        ENDDO
     end if
  case(9)
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i)) ) > small) then
           call becke88_spin (rhoup(i), grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call becke88_spin (rhodw(i), grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = sxup + sxdw

     if (exx_started ) then
        sx = 0.72_DP * sx
        v1xup = 0.72_DP * v1xup
        v1xdw = 0.72_DP * v1xdw
        v2xup = 0.72_DP * v2xup
        v2xdw = 0.72_DP * v2xdw
     end if

  case(11) ! 'Wu-Cohen'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call wcx (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call wcx (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(13) ! 'rPW86 for vdw-df2'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call rPW86 (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call rPW86 (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(16) ! 'c09x for vdw-df-c09'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call c09x (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call c09x (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(21) ! 'pw86'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call pw86 (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call pw86 (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(22) ! 'b86b'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call becke86b (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call becke86b (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(26) ! 'B86R for rev-vdW-DF2'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call b86b (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), 3, sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call b86b (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), 3, sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case(27) ! 'cx13 for vdw-df-cx'
     do i=1,length
        if (rhoup(i) > small .and. sqrt(abs(grhoup2(i))) > small) then
           call cx13 (2.0_DP * rhoup(i), 4.0_DP * grhoup2(i), sxup(i), v1xup(i), v2xup(i))
        else
           sxup(i) = 0.0_DP
           v1xup(i) = 0.0_DP
           v2xup(i) = 0.0_DP
        endif
        if (rhodw(i) > small .and. sqrt(abs(grhodw2(i))) > small) then
           call cx13 (2.0_DP * rhodw(i), 4.0_DP * grhodw2(i), sxdw(i), v1xdw(i), v2xdw(i))
        else
           sxdw(i) = 0.0_DP
           v1xdw(i) = 0.0_DP
           v2xdw(i) = 0.0_DP
        endif
     end do
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw

  case default
     call errore ('gcx_spin_vec', 'not implemented', igcx)
  end select
  !
  if (igcx.ne.0) then
     where (rho.le.small)
        sx = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     end where
  end if
  !
end subroutine gcx_spin_vec
!
!-----------------------------------------------------------------------
subroutine gcc_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for correlations - Hartree a.u.
  !     Implemented:  Perdew86, GGA (PW91), PBE
  !
  implicit none
  !
  !     dummy arguments
  !
  real(DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  ! the total charge
  ! the magnetization
  ! the gradient of the charge squared
  ! exchange and correlation energies
  ! derivatives of correlation wr. rho
  ! derivatives of correlation wr. grho

  real(DP), parameter :: small = 1.E-10_DP, epsr=1.E-6_DP
  !
  if ( abs(zeta) > 1.0_DP ) then
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     v2c = 0.0_DP
     return
  else
     !
     ! ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
     zeta = SIGN( MIN( ABS( zeta ), ( 1.0_DP - epsr ) ) , zeta )
  endif

  if (igcc == 0 .or. rho <= small .or. sqrt(abs(grho)) <= small) then
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     v2c = 0.0_DP
  elseif (igcc == 1) then
     call perdew86_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 2) then
     call ggac_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 4) then
     call pbec_spin (rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 8) then
     call pbec_spin (rho, zeta, grho, 2, sc, v1cup, v1cdw, v2c)
  else
     call errore ('lsda_functionals (gcc_spin)', 'not implemented', igcc)
  endif
  !
  return
end subroutine gcc_spin
!
!   ==================================================================
    SUBROUTINE gcc_spin_more( RHOA, RHOB, GRHOAA, GRHOBB, GRHOAB, &
                              SC, V1CA, V1CB, V2CA, V2CB, V2CAB )
!   ==--------------------------------------------------------------==
!   ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
!   ==                                                              ==
!   ==  EXCHANGE  :  BECKE88                                        ==
!   ==               GGAX                                           ==
!   ==  CORRELATION : PERDEW86                                      ==
!   ==                LEE, YANG & PARR                              ==
!   ==                GGAC                                          ==
!   ==--------------------------------------------------------------==

      IMPLICIT NONE
      REAL(DP) :: RHOA,RHOB,GRHOAA,GRHOBB,GRHOAB
      REAL(DP) :: SC,V1CA,V2CA,V1CB,V2CB,V2CAB

      ! ... Gradient Correction for correlation

      REAL(DP) :: SMALL, RHO
      PARAMETER(SMALL=1.E-20_DP)

      SC=0.0_DP
      V1CA=0.0_DP
      V2CA=0.0_DP
      V1CB=0.0_DP
      V2CB=0.0_DP
      V2CAB=0.0_DP
      IF( igcc == 3 .or. igcc == 7) THEN
        RHO=RHOA+RHOB
        IF(RHO.GT.SMALL) then
             CALL LSD_GLYP(RHOA,RHOB,GRHOAA,GRHOAB,GRHOBB,SC,&
                  V1CA,V2CA,V1CB,V2CB,V2CAB)
             if (igcc == 7 .and. exx_started) then
                SC = 0.81d0*SC
                V1CA = 0.81d0*V1CA
                V2CA = 0.81d0*V2CA
                V1CB = 0.81d0*V1CB
                V2CB = 0.81d0*V2CB
                V2CAB = 0.81d0*V2CAB
             endif
         endif
      ELSE
        CALL errore( " gcc_spin_more ", " gradiet correction not implemented ", 1 )
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE gcc_spin_more
!
!
!-----------------------------------------------------------------------
!------- NONLOCAL CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
subroutine nlc (rho_valence, rho_core, nspin, enl, vnl, v)
  !-----------------------------------------------------------------------
  !     non local correction for the correlation
  !
  !     input:  rho_valence, rho_core
  !     definition:  E_nl = \int E_nl(rho',grho',rho'',grho'',|r'-r''|) dr
  !     output: enl = E_nl
  !             vnl= D(E_x)/D(rho)
  !             v  = Correction to the potential
  !

  USE vdW_DF, ONLY: xc_vdW_DF, xc_vdW_DF_spin, vdw_type
  USE rVV10,  ONLY: xc_rVV10
 
  implicit none
  
  REAL(DP), INTENT(IN) :: rho_valence(:,:), rho_core(:)
  INTEGER, INTENT(IN)  :: nspin
  REAL(DP), INTENT(INOUT) :: v(:,:)
  REAL(DP), INTENT(INOUT) :: enl, vnl

  if ( inlc == 1 .or. inlc == 2 .or. inlc == 4 .or. inlc == 5 .or. inlc == 6 ) then

     vdw_type = inlc
     if( nspin == 1 ) then
        call xc_vdW_DF      (rho_valence, rho_core, enl, vnl, v)
     else if( nspin == 2 ) then
        call xc_vdW_DF_spin (rho_valence, rho_core, enl, vnl, v)
     else
        call errore ('nlc','vdW-DF not available for noncollinear spin case',1)
     end if

  elseif (inlc == 3) then

      call xc_rVV10 (rho_valence, rho_core, nspin, enl, vnl, v)
  
  else
     enl = 0.0_DP
     vnl = 0.0_DP
     v = 0.0_DP
  endif
  !
  return
end subroutine nlc

!
!-----------------------------------------------------------------------
!------- META CORRECTIONS DRIVERS ----------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine tau_xc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     See comments at the beginning of module for implemented cases
  !
  !     input:  rho, grho=|\nabla rho|^2
  !
  !     definition:  E_x = \int e_x(rho,grho) dr
  !
  !     output: sx = e_x(rho,grho) = grad corr
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             v3x= D(E_x)/D(tau)
  !
  !             sc, v1c, v2c as above for correlation
  !
  implicit none

  real(DP) :: rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c
  
  !_________________________________________________________________________
  
  if     (imeta == 1) then
     call tpsscxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 2) then
     call   m06lxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 3) then
     call  tb09cxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
  elseif (imeta == 4) then
     ! do nothing
  else
     call errore('v_xc_meta','wrong igcx and/or igcc',1)
  end if
  
  return
  
end subroutine tau_xc

!
!
!-----------------------------------------------------------------------
subroutine tau_xc_spin (rhoup, rhodw, grhoup, grhodw, tauup, taudw, ex, ec,   &
           &            v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw, v1cup, v1cdw,&
           &            v2cup, v2cdw, v3cup, v3cdw)

!-----------------------------------------------------------------------
  !
  !
  
  implicit none

  real(dp), intent(in)                :: rhoup, rhodw, tauup, taudw
  real(dp), dimension (3), intent(in) :: grhoup, grhodw
  
  real(dp), intent(out)               :: ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
                                      &  v1cup, v1cdw, v3cup, v3cdw
  real(dp), dimension(3), intent(out) :: v2cup, v2cdw
  
  !
  !  Local variables
  !
  integer                 :: ipol
  real(dp)                :: rh, zeta, atau, grhoup2, grhodw2
  real(dp), parameter     :: epsr=1.0d-08, zero=0._dp
  !
  !_____________________________

  grhoup2 = zero
  grhodw2 = zero
  
  v2cup         = zero
  v2cdw         = zero

  do ipol=1,3
     grhoup2 = grhoup2 + grhoup(ipol)**2
     grhodw2 = grhodw2 + grhodw(ipol)**2
  end do

  
  if (imeta == 1) then

     call tpsscx_spin(rhoup, rhodw, grhoup2, grhodw2, tauup,   &
              &  taudw, ex, v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
  
     rh   =  rhoup + rhodw
        
     zeta = (rhoup - rhodw) / rh
     atau =  tauup + taudw    ! KE-density in Hartree

     call tpsscc_spin(rh,zeta,grhoup,grhodw, atau,ec,              &
     &                v1cup,v1cdw,v2cup,v2cdw,v3cup, v3cdw) 
  
  
  elseif (imeta == 2) then
  
     call   m06lxc_spin (rhoup, rhodw, grhoup2, grhodw2, tauup, taudw,      &
            &            ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
            &            v1cup, v1cdw, v2cup(1), v2cdw(1), v3cup, v3cdw)
     
  else
  
     call errore('v_xc_meta','wrong igcx and/or igcc',1)
     
  end if
  
end subroutine tau_xc_spin                
                

!-----------------------------------------------------------------------
!------- DRIVERS FOR DERIVATIVES OF XC POTENTIAL -----------------------
!-----------------------------------------------------------------------
!
      !-----------------------------------------------------------------------
      function dmxc (rho)
        !-----------------------------------------------------------------------
        !
        !  derivative of the xc potential with respect to the local density
        !
        !
        implicit none
        !
        real(DP), intent(in) :: rho
        ! input: the charge density ( positive )
        real(DP) :: dmxc
        ! output: the derivative of the xc potential
        !
        ! local variables
        !
        real(DP) :: dr, vxp, vcp, vxm, vcm, vx, ex, ec, rs
        real(DP), external :: dpz
        integer :: iflg
        !
        real(DP), parameter :: small = 1.E-30_DP, e2 = 2.0_DP, &
             pi34 = 0.75_DP / 3.141592653589793_DP, third = 1.0_DP /3.0_DP
        !
        dmxc = 0.0_DP
        if (rho < small) then
           return
        endif
        !
        !    first case: analytical derivatives available
        !
        if (get_iexch() == 1 .and. get_icorr() == 1) then
           rs = (pi34 / rho) **third
           !..exchange
           call slater (rs, ex, vx)
           dmxc = vx / (3.0_DP * rho)
           !..correlation
           iflg = 2
           if (rs < 1.0_DP) iflg = 1
           dmxc = dmxc + dpz (rs, iflg)
        else
           !
           !     second case: numerical derivatives
           !
           dr = min (1.E-6_DP, 1.E-4_DP * rho)
           call xc (rho + dr, ex, ec, vxp, vcp)
           call xc (rho - dr, ex, ec, vxm, vcm)
           dmxc = (vxp + vcp - vxm - vcm) / (2.0_DP * dr)
        endif
        !
        ! bring to rydberg units
        !
        dmxc = e2 * dmxc
        return
        !
      end function dmxc
      !
      !-----------------------------------------------------------------------
      subroutine dmxc_spin (rhoup, rhodw, dmuxc_uu, dmuxc_ud, dmuxc_du, &
           dmuxc_dd)
      !-----------------------------------------------------------------------
        !  derivative of the xc potential with respect to the local density
        !  spin-polarized case
        !
        implicit none
        !
        real(DP), intent(in) :: rhoup, rhodw
        ! input: spin-up and spin-down charge density
        real(DP), intent(out) :: dmuxc_uu, dmuxc_ud, dmuxc_du, dmuxc_dd
        ! output: up-up, up-down, down-up, down-down derivatives of the
        ! XC functional
        !
        ! local variables
        !
        real(DP) :: rhotot, rs, zeta, fz, fz1, fz2, ex, vx, ecu, ecp, vcu, &
             vcp, dmcu, dmcp, aa, bb, cc, dr, dz, ec, vxupm, vxdwm, vcupm, &
             vcdwm, rho, vxupp, vxdwp, vcupp, vcdwp, zeta_eff
        real(DP), external :: dpz, dpz_polarized
        integer :: iflg
        !
        real(DP), parameter :: small = 1.E-30_DP, e2 = 2.0_DP, &
             pi34 = 0.75_DP / 3.141592653589793_DP, third = 1.0_DP/3.0_DP, &
             p43 = 4.0_DP / 3.0_DP, p49 = 4.0_DP / 9.0_DP, m23 = -2.0_DP / 3.0_DP
        !
        dmuxc_uu = 0.0_DP
        dmuxc_du = 0.0_DP
        dmuxc_ud = 0.0_DP
        dmuxc_dd = 0.0_DP
        !
        rhotot = rhoup + rhodw
        if (rhotot <= small) return
        zeta = (rhoup - rhodw) / rhotot
        
        if (abs (zeta) > 1.0_DP) return
        if (get_iexch() == 1 .and. get_icorr() == 1) then
           !
           !    first case: analytical derivative available
           !
           !..exchange
           rs = (pi34 / (2.0_DP * rhoup) ) **third
           call slater (rs, ex, vx)
           dmuxc_uu = vx / (3.0_DP * rhoup)
           rs = (pi34 / (2.0_DP * rhodw) ) **third
           call slater (rs, ex, vx)
           dmuxc_dd = vx / (3.0_DP * rhodw)
           !..correlation
           rs = (pi34 / rhotot) **third
           iflg = 2
           if (rs < 1.0_DP) iflg = 1
           dmcu = dpz (rs, iflg)
           dmcp = dpz_polarized (rs, iflg)
           call pz (rs, 1, ecu, vcu)
           call pz_polarized (rs, ecp, vcp)
           fz = ( (1.0_DP + zeta) **p43 + (1.0_DP - zeta) **p43 - 2.0_DP) &
                / (2.0_DP**p43 - 2.0_DP)
           fz1 = p43 * ( (1.0_DP + zeta) **third- (1.0_DP - zeta) **third) &
                / (2.0_DP**p43 - 2.0_DP)
           fz2 = p49 * ( (1.0_DP + zeta) **m23 + (1.0_DP - zeta) **m23) &
                / (2.0_DP**p43 - 2.0_DP)
           aa = dmcu + fz * (dmcp - dmcu)
           bb = 2.0_DP * fz1 * (vcp - vcu - (ecp - ecu) ) / rhotot
           cc = fz2 * (ecp - ecu) / rhotot
           dmuxc_uu = dmuxc_uu + aa + (1.0_DP - zeta) * bb + (1.0_DP - zeta)**2 * cc
           dmuxc_du = dmuxc_du + aa + ( - zeta) * bb + (zeta**2 - 1.0_DP) * cc
           dmuxc_ud = dmuxc_du
           dmuxc_dd = dmuxc_dd+aa - (1.0_DP + zeta) * bb + (1.0_DP + zeta)**2 * cc
           
        else
           
           rho = rhoup + rhodw
           dr = min (1.E-6_DP, 1.E-4_DP * rho)
           call xc_spin (rho - dr, zeta, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
           call xc_spin (rho + dr, zeta, ex, ec, vxupp, vxdwp, vcupp, vcdwp)
           dmuxc_uu = (vxupp + vcupp - vxupm - vcupm) / (2.0_DP * dr)
           dmuxc_ud = dmuxc_uu
           dmuxc_dd = (vxdwp + vcdwp - vxdwm - vcdwm) / (2.0_DP * dr)
           dmuxc_du = dmuxc_dd
           ! dz = min (1.d-6, 1.d-4 * abs (zeta) )
           dz = 1.E-6_DP
!
!          If zeta is too close to +-1, the derivative is computed at a slightly
!          smaller zeta
!
           zeta_eff = SIGN( MIN( ABS( zeta ), ( 1.0_DP - 2.0_DP*dz ) ) , zeta )

           call xc_spin (rho, zeta_eff - dz, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
           call xc_spin (rho, zeta_eff + dz, ex, ec, vxupp, vxdwp, vcupp, vcdwp)
           dmuxc_uu = dmuxc_uu + (vxupp + vcupp - vxupm - vcupm) * &
                (1.0_DP - zeta) / rho / (2.0_DP * dz)
           dmuxc_ud = dmuxc_ud- (vxupp + vcupp - vxupm - vcupm) * &
                (1.0_DP + zeta) / rho / (2.0_DP * dz)
           dmuxc_du = dmuxc_du + (vxdwp + vcdwp - vxdwm - vcdwm) * &
                (1.0_DP - zeta) / rho / (2.0_DP * dz)
           dmuxc_dd = dmuxc_dd- (vxdwp + vcdwp - vxdwm - vcdwm) * &
                (1.0_DP + zeta) / rho / (2.0_DP * dz)
        endif
        !
        ! bring to rydberg units
        !
        dmuxc_uu = e2 * dmuxc_uu
        dmuxc_du = e2 * dmuxc_du
        dmuxc_ud = e2 * dmuxc_ud
        dmuxc_dd = e2 * dmuxc_dd
        !
        return
        
      end subroutine dmxc_spin

      !-----------------------------------------------------------------------
      subroutine dmxc_nc (rho, mx, my, mz, dmuxc)
      !-----------------------------------------------------------------------
        !  derivative of the xc potential with respect to the local density
        !  and magnetization
        !  non colinear case
        !
        implicit none
        !
        real(DP), intent(in) :: rho, mx, my, mz
        ! input: charge density and magnetization
        real(DP), intent(out) :: dmuxc(4,4)
        ! output: derivative of XC functional
        !
        ! local variables
        !
        REAL(DP) :: zeta, ex, ec, dr, dz, vxupm, vxdwm, vcupm, &
              vcdwm, vxupp, vxdwp, vcupp, vcdwp, vxup, vxdw, vcup, vcdw
        REAL(DP) :: amag, vs, dvxc_rho, dvxc_mx, dvxc_my, dvxc_mz,  &
                    dbx_rho, dbx_mx, dbx_my, dbx_mz, dby_rho, dby_mx, &
                    dby_my, dby_mz, dbz_rho, dbz_mx, dbz_my, dbz_mz, zeta_eff
        REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP
        !
        !
        dmuxc = 0.0_DP
        !
        IF (rho <= small) RETURN
        amag = sqrt(mx**2+my**2+mz**2) 
        zeta = amag / rho
        
        IF (abs (zeta) > 1.0_DP) RETURN
        CALL xc_spin (rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw)
        vs=0.5_DP*(vxup+vcup-vxdw-vcdw)
   
        dr = min (1.E-6_DP, 1.E-4_DP * rho)
        CALL xc_spin (rho - dr, zeta, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
        CALL xc_spin (rho + dr, zeta, ex, ec, vxupp, vxdwp, vcupp, vcdwp)
        dvxc_rho = ((vxupp + vcupp - vxupm - vcupm)+     &
                    (vxdwp + vcdwp - vxdwm - vcdwm)) / (4.0_DP * dr)
        IF (amag > 1.E-10_DP) THEN
           dbx_rho  = ((vxupp + vcupp - vxupm - vcupm)-     &
                       (vxdwp + vcdwp - vxdwm - vcdwm))* mx / (4.0_DP*dr*amag)
           dby_rho  = ((vxupp + vcupp - vxupm - vcupm)-     &
                       (vxdwp + vcdwp - vxdwm - vcdwm))* my / (4.0_DP*dr*amag)
           dbz_rho  = ((vxupp + vcupp - vxupm - vcupm)-     &
                       (vxdwp + vcdwp - vxdwm - vcdwm))* mz / (4.0_DP*dr*amag)
!           dz = min (1.d-6, 1.d-4 * abs (zeta) )
           dz = 1.0E-6_DP
!
!          If zeta is too close to +-1, the derivative is computed at a slightly
!          smaller zeta
!
           zeta_eff = SIGN( MIN( ABS( zeta ), ( 1.0_DP - 2.0_DP*dz ) ) , zeta )

           CALL xc_spin (rho, zeta_eff - dz, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
           CALL xc_spin (rho, zeta_eff + dz, ex, ec, vxupp, vxdwp, vcupp, vcdwp)

!  The variables are rho and m, so zeta depends on rho
!
           dvxc_rho=dvxc_rho- ((vxupp + vcupp - vxupm - vcupm)+     &
                         (vxdwp + vcdwp - vxdwm - vcdwm))*zeta/rho/(4.0_DP * dz)
           dbx_rho  = dbx_rho-((vxupp + vcupp - vxupm - vcupm)-     &
                    (vxdwp + vcdwp - vxdwm - vcdwm))*mx*zeta/rho/(4.0_DP*dz*amag)
           dby_rho  = dby_rho-((vxupp + vcupp - vxupm - vcupm)-     &
                    (vxdwp + vcdwp - vxdwm - vcdwm))*my*zeta/rho/(4.0_DP*dz*amag)
           dbz_rho  = dbz_rho-((vxupp + vcupp - vxupm - vcupm)-     &
                    (vxdwp + vcdwp - vxdwm - vcdwm))*mz*zeta/rho/(4.0_DP*dz*amag)
!
! here the derivatives with respect to m
!
           dvxc_mx = ((vxupp + vcupp - vxupm - vcupm) + &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*mx/rho/(4.0_DP*dz*amag)
           dvxc_my = ((vxupp + vcupp - vxupm - vcupm) + &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*my/rho/(4.0_DP*dz*amag)
           dvxc_mz = ((vxupp + vcupp - vxupm - vcupm) + &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*mz/rho/(4.0_DP*dz*amag)
           dbx_mx  = (((vxupp + vcupp - vxupm - vcupm) -                 &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*mx**2*amag/rho/       &
                      (4.0_DP*dz) + vs*(my**2+mz**2))/amag**3
           dbx_my  = (((vxupp + vcupp - vxupm - vcupm) -                 &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*mx*my*amag/rho/       &
                      (4.0_DP*dz) - vs*(mx*my))/amag**3
           dbx_mz  = (((vxupp + vcupp - vxupm - vcupm) -                 &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*mx*mz*amag/rho/       &
                      (4.0_DP*dz) - vs*(mx*mz))/amag**3
           dby_mx  = dbx_my
           dby_my  = (((vxupp + vcupp - vxupm - vcupm) -                 &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*my**2*amag/rho/       &
                      (4.0_DP*dz) + vs*(mx**2+mz**2))/amag**3
           dby_mz  = (((vxupp + vcupp - vxupm - vcupm) -                 &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*my*mz*amag/rho/  &
                      (4.0_DP*dz) - vs*(my*mz))/amag**3
           dbz_mx  = dbx_mz
           dbz_my  = dby_mz
           dbz_mz  = (((vxupp + vcupp - vxupm - vcupm) -                 &
                      (vxdwp + vcdwp - vxdwm - vcdwm))*mz**2*amag/rho/       &
                      (4.0_DP*dz) + vs*(mx**2+my**2))/amag**3
           dmuxc(1,1)=dvxc_rho
           dmuxc(1,2)=dvxc_mx 
           dmuxc(1,3)=dvxc_my 
           dmuxc(1,4)=dvxc_mz 
           dmuxc(2,1)=dbx_rho
           dmuxc(2,2)=dbx_mx 
           dmuxc(2,3)=dbx_my 
           dmuxc(2,4)=dbx_mz 
           dmuxc(3,1)=dby_rho
           dmuxc(3,2)=dby_mx 
           dmuxc(3,3)=dby_my 
           dmuxc(3,4)=dby_mz 
           dmuxc(4,1)=dbz_rho
           dmuxc(4,2)=dbz_mx
           dmuxc(4,3)=dbz_my
           dmuxc(4,4)=dbz_mz
        ELSE
           dmuxc(1,1)=dvxc_rho
        ENDIF
        !
        ! bring to rydberg units
        !
        dmuxc = e2 * dmuxc
        !
        RETURN
        
      end subroutine dmxc_nc
      !
      !-----------------------------------------------------------------------
      subroutine dgcxc (r, s2, vrrx, vsrx, vssx, vrrc, vsrc, vssc)
        !-----------------------------------------------------------------------
        USE kinds, only : DP
        implicit none
        real(DP) :: r, s2, vrrx, vsrx, vssx, vrrc, vsrc, vssc
        real(DP) :: dr, s, ds
        
        real(DP) :: sx, sc, v1xp, v2xp, v1cp, v2cp, v1xm, v2xm, v1cm, &
             v2cm
        s = sqrt (s2)
        dr = min (1.d-4, 1.d-2 * r)
        
        ds = min (1.d-4, 1.d-2 * s)
        call gcxc (r + dr, s2, sx, sc, v1xp, v2xp, v1cp, v2cp)
        
        call gcxc (r - dr, s2, sx, sc, v1xm, v2xm, v1cm, v2cm)
        vrrx = 0.5d0 * (v1xp - v1xm) / dr
        
        vrrc = 0.5d0 * (v1cp - v1cm) / dr
        vsrx = 0.25d0 * (v2xp - v2xm) / dr
        
        vsrc = 0.25d0 * (v2cp - v2cm) / dr
        call gcxc (r, (s + ds) **2, sx, sc, v1xp, v2xp, v1cp, v2cp)
        
        call gcxc (r, (s - ds) **2, sx, sc, v1xm, v2xm, v1cm, v2cm)
        vsrx = vsrx + 0.25d0 * (v1xp - v1xm) / ds / s
        
        vsrc = vsrc + 0.25d0 * (v1cp - v1cm) / ds / s
        vssx = 0.5d0 * (v2xp - v2xm) / ds / s
        
        vssc = 0.5d0 * (v2cp - v2cm) / ds / s
        return
      end subroutine dgcxc
      !
      !-----------------------------------------------------------------------
      subroutine dgcxc_spin (rup, rdw, gup, gdw, vrrxup, vrrxdw, vrsxup, &
           vrsxdw, vssxup, vssxdw, vrrcup, vrrcdw, vrscup, vrscdw, vssc, &
           vrzcup, vrzcdw)
        !-----------------------------------------------------------------------
        !
        !    This routine computes the derivative of the exchange and correlatio
        !    potentials with respect to the density, the gradient and zeta
        !
        USE kinds, only : DP
        implicit none
        real(DP), intent(in) :: rup, rdw, gup (3), gdw (3)
        ! input: the charges and the gradient
        real(DP), intent(out):: vrrxup, vrrxdw, vrsxup, vrsxdw, vssxup, &
             vssxdw, vrrcup, vrrcdw, vrscup, vrscdw, vssc, vrzcup, vrzcdw
        ! output: derivatives of the exchange and of the correlation
        !
        !    local variables
        !
        real(DP) :: r, zeta, sup2, sdw2, s2, s, sup, sdw, dr, dzeta, ds, &
             drup, drdw, dsup, dsdw, sx, sc, v1xupp, v1xdwp, v2xupp, v2xdwp, &
             v1xupm, v1xdwm, v2xupm, v2xdwm, v1cupp, v1cdwp, v2cp, v1cupm, &
             v1cdwm, v2cm
        ! charge densities and square gradients
        ! delta charge densities and gra
        ! delta gradients
        ! energies
        ! exchange potentials
        ! exchange potentials
        ! coorelation potentials
        ! coorelation potentials
        real(DP), parameter :: eps = 1.d-6
        !
        r = rup + rdw
        if (r.gt.eps) then
           zeta = (rup - rdw) / r
        else
           zeta = 2.d0
        endif
        sup2 = gup (1) **2 + gup (2) **2 + gup (3) **2
        sdw2 = gdw (1) **2 + gdw (2) **2 + gdw (3) **2

        s2 = (gup (1) + gdw (1) ) **2 + (gup (2) + gdw (2) ) **2 + &
             (gup (3) + gdw (3) ) **2
        sup = sqrt (sup2)
        sdw = sqrt (sdw2)
        s = sqrt (s2)
        !
        !     up part of exchange
        !
        
        if (rup.gt.eps.and.sup.gt.eps) then
           drup = min (1.d-4, 1.d-2 * rup)
           dsup = min (1.d-4, 1.d-2 * sdw)
           !
           !    derivatives of exchange: up part
           !
           call gcx_spin (rup + drup, rdw, sup2, sdw2, sx, v1xupp, v1xdwp, &
                v2xupp, v2xdwp)
           
           call gcx_spin (rup - drup, rdw, sup2, sdw2, sx, v1xupm, v1xdwm, &
                v2xupm, v2xdwm)
           vrrxup = 0.5d0 * (v1xupp - v1xupm) / drup
           vrsxup = 0.25d0 * (v2xupp - v2xupm) / drup

           call gcx_spin (rup, rdw, (sup + dsup) **2, sdw2, sx, v1xupp, &
                v1xdwp, v2xupp, v2xdwp)
 
           call gcx_spin (rup, rdw, (sup - dsup) **2, sdw2, sx, v1xupm, &
                v1xdwm, v2xupm, v2xdwm)
           vrsxup = vrsxup + 0.25d0 * (v1xupp - v1xupm) / dsup / sup
           vssxup = 0.5d0 * (v2xupp - v2xupm) / dsup / sup
        else
           vrrxup = 0.d0
           vrsxup = 0.d0
           vssxup = 0.d0
        endif
        
        if (rdw.gt.eps.and.sdw.gt.eps) then
           drdw = min (1.d-4, 1.d-2 * rdw)
           dsdw = min (1.d-4, 1.d-2 * sdw)
           !
           !    derivatives of exchange: down part
           !
           call gcx_spin (rup, rdw + drdw, sup2, sdw2, sx, v1xupp, v1xdwp, &
                v2xupp, v2xdwp)

           call gcx_spin (rup, rdw - drdw, sup2, sdw2, sx, v1xupm, v1xdwm, &
                v2xupm, v2xdwm)
           vrrxdw = 0.5d0 * (v1xdwp - v1xdwm) / drdw

           vrsxdw = 0.25d0 * (v2xdwp - v2xdwm) / drdw
           call gcx_spin (rup, rdw, sup2, (sdw + dsdw) **2, sx, v1xupp, &
                v1xdwp, v2xupp, v2xdwp)

           call gcx_spin (rup, rdw, sup2, (sdw - dsdw) **2, sx, v1xupm, &
                v1xdwm, v2xupm, v2xdwm)
           vrsxdw = vrsxdw + 0.25d0 * (v1xdwp - v1xdwm) / dsdw / sdw
           vssxdw = 0.5d0 * (v2xdwp - v2xdwm) / dsdw / sdw
        else
           vrrxdw = 0.d0
           vrsxdw = 0.d0
           vssxdw = 0.d0
        endif
        !
        !     derivatives of correlation
        !
        
        if (r.gt.eps.and.abs (zeta) .le.1.d0.and.s.gt.eps) then
           
           dr = min (1.d-4, 1.d-2 * r)
           call gcc_spin (r + dr, zeta, s2, sc, v1cupp, v1cdwp, v2cp)
           
           call gcc_spin (r - dr, zeta, s2, sc, v1cupm, v1cdwm, v2cm)
           vrrcup = 0.5d0 * (v1cupp - v1cupm) / dr
           
           vrrcdw = 0.5d0 * (v1cdwp - v1cdwm) / dr
           
           ds = min (1.d-4, 1.d-2 * s)
           call gcc_spin (r, zeta, (s + ds) **2, sc, v1cupp, v1cdwp, v2cp)

           call gcc_spin (r, zeta, (s - ds) **2, sc, v1cupm, v1cdwm, v2cm)
           vrscup = 0.5d0 * (v1cupp - v1cupm) / ds / s
           vrscdw = 0.5d0 * (v1cdwp - v1cdwm) / ds / s

           vssc = 0.5d0 * (v2cp - v2cm) / ds / s
!           dzeta = min (1.d-4, 1.d-2 * abs (zeta) )

           dzeta = 1.d-6
!
!   If zeta is too close to +-1 the derivative is evaluated at a slightly 
!   smaller  value
!
           zeta = SIGN( MIN( ABS( zeta ), ( 1.0_DP - 2.0_DP*dzeta ) ) , zeta )

           call gcc_spin (r, zeta + dzeta, s2, sc, v1cupp, v1cdwp, v2cp)
           
           call gcc_spin (r, zeta - dzeta, s2, sc, v1cupm, v1cdwm, v2cm)
           vrzcup = 0.5d0 * (v1cupp - v1cupm) / dzeta
           vrzcdw = 0.5d0 * (v1cdwp - v1cdwm) / dzeta
        else
           vrrcup = 0.d0
           vrrcdw = 0.d0
           vrscup = 0.d0
           vrscdw = 0.d0
           vssc = 0.d0
           vrzcup = 0.d0
           vrzcdw = 0.d0

        endif
        return
      end subroutine dgcxc_spin

    !-----------------------------------------------------------------------
    subroutine d3gcxc (r, s2, vrrrx, vsrrx, vssrx, vsssx, &
         vrrrc, vsrrc, vssrc, vsssc )
    !-----------------------------------------------------------------------
    !
    !    wat20101006: Calculates all derviatives of the exchange (x) and 
    !                 correlation (c) potential in third order.  
    !                 of the Exc.
    !
    !    input:       r = rho, s2=|\nabla rho|^2
    !    definition:  E_xc = \int ( f_x(r,s2) + f_c(r,s2) ) dr
    !    output:      vrrrx = d^3(f_x)/d(r)^3
    !                 vsrrx = d^3(f_x)/d(|\nabla r|)d(r)^2 / |\nabla r|
    !                 vssrx = d/d(|\nabla r|) [ &
    !                           d^2(f_x)/d(|\nabla r|)d(r) / |\nabla r| ] &
    !                                                           / |\nabla r|
    !                 vsssx = d/d(|\nabla r|) [ &
    !                           d/d(|\nabla r|) [ & 
    !                           d(f_x)/d(|\nabla r|) / |\nabla r| ] &
    !                                                   / |\nabla r| ] &
    !                                                   / |\nabla r|
    !                 same for (c)
    !
    USE kinds, ONLY : DP
    IMPLICIT NONE
    REAL(DP) :: r, s2, vrrrx, vsrrx, vssrx, vsssx, &
                vrrrc, vsrrc, vssrc, vsssc
    REAL(DP) :: dr, s, ds
    !
    REAL(DP) :: vrrx_rp, vsrx_rp, vssx_rp, vrrc_rp, vsrc_rp, vssc_rp, &
                vrrx_rm, vsrx_rm, vssx_rm, vrrc_rm, vsrc_rm, vssc_rm, &
                vrrx_sp, vsrx_sp, vssx_sp, vrrc_sp, vsrc_sp, vssc_sp, &
                vrrx_sm, vsrx_sm, vssx_sm, vrrc_sm, vsrc_sm, vssc_sm
    !
    s = sqrt (s2)
    dr = min (1.d-4, 1.d-2 * r)
    ds = min (1.d-4, 1.d-2 * s)
    !  
    call dgcxc (r+dr, s2, vrrx_rp, vsrx_rp, vssx_rp, vrrc_rp, vsrc_rp, vssc_rp)
    call dgcxc (r-dr, s2, vrrx_rm, vsrx_rm, vssx_rm, vrrc_rm, vsrc_rm, vssc_rm)
    !  
    call dgcxc (r, (s+ds)**2, vrrx_sp, vsrx_sp, vssx_sp, vrrc_sp, vsrc_sp, vssc_sp)
    call dgcxc (r, (s-ds)**2, vrrx_sm, vsrx_sm, vssx_sm, vrrc_sm, vsrc_sm, vssc_sm)
    !  
    vrrrx = 0.5d0 * (vrrx_rp - vrrx_rm) / dr
    vsrrx = 0.25d0 * (vsrx_rp - vsrx_rm) / dr &
                  + 0.25d0 * (vrrx_sp - vrrx_sm) / ds / s
    vssrx = 0.25d0 * (vssx_rp - vssx_rm) / dr &
                  + 0.25d0 * (vsrx_sp - vsrx_sm) / ds / s
    vsssx = 0.5d0 * (vssx_sp - vssx_sm) / ds / s
    !
    vrrrc = 0.5d0 * (vrrc_rp - vrrc_rm) / dr
    vsrrc = 0.25d0 * (vsrc_rp - vsrc_rm) / dr &
                  + 0.25d0 * (vrrc_sp - vrrc_sm) / ds / s
    vssrc = 0.25d0 * (vssc_rp - vssc_rm) / dr &
                  + 0.25d0 * (vsrc_sp - vsrc_sm) / ds / s
    vsssc = 0.5d0 * (vssc_sp - vssc_sm) / ds / s
    !
    return
    !
  end subroutine d3gcxc
!
!-----------------------------------------------------------------------
!------- VECTOR AND GENERAL XC DRIVERS -------------------------------
!-----------------------------------------------------------------------
!

subroutine evxc_t_vec(rho,rhoc,lsd,length,vxc,exc)
  !---------------------------------------------------------------
  !
  !  this function returns the XC potential in LDA or LSDA approximation
  !
  integer,  intent(in)  :: lsd, length
  real(DP), intent(in)  :: rho(length,2), rhoc(length)
  real(DP), intent(out), optional :: vxc(length,2)
  real(DP), intent(out), optional :: exc(length)
  !
  real(DP) :: arho
  real(DP) :: arhoV(length), zetaV(length)
  real(DP) :: evx(length,3), evc(length,3)
  real(DP) :: ex, ec, vx, vc
  !
  integer :: i
  real(DP), parameter :: e2 = 2.0_dp, eps = 1.e-30_dp

  if (lsd.eq.0) then
     !
     !     LDA case
     !
     do i=1,length
        arho = abs(rho(i,1)+rhoc(i))
        if (arho.gt.eps) then
           call xc(arho,ex,ec,vx,vc)
        else
           ex = 0.0_dp
           ec = 0.0_dp
           vx = 0.0_dp
           vc = 0.0_dp
        end if
        if (present(vxc)) vxc(i,1) = e2*(vx+vc)
        if (present(exc)) exc(i) = e2*(ex+ec)
     end do
  else
     !
     !     LSDA case
     !
     arhoV = abs(rho(:,1)+rho(:,2)+rhoc(:))
     where (arhoV.gt.eps)
        zetaV = (rho(:,1)-rho(:,2)) / arhoV
     elsewhere
        zetaV = 0.0_DP ! just a sane default, results are discarded anyway
     end where
     ! zeta has to stay between -1 and 1, but can get a little
     ! out of bound during the first iterations.
     zetaV = min( 1.0_DP, zetaV)
     zetaV = max(-1.0_DP, zetaV)
     call xc_spin_vec(arhoV, zetaV, length, evx, evc)
     if (present(vxc)) then
        vxc(:,1) = e2*(evx(:,1) + evc(:,1))
        vxc(:,2) = e2*(evx(:,2) + evc(:,2))
     end if
     if (present(exc)) exc = e2*(evx(:,3)+evc(:,3))
  end if

end subroutine evxc_t_vec

end module funct

