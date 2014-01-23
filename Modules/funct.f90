!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
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
!  retrieve functions: get_dft_name
!                      get_iexch
!                      get_icorr
!                      get_igcx
!                      get_igcc
!                      get_exx_fraction
!                      dft_name
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
  PUBLIC  :: enforce_input_dft, write_dft_name, dft_name
  PUBLIC  :: init_dft_exxrpa, enforce_dft_exxrpa
  PUBLIC  :: get_dft_name, get_iexch, get_icorr, get_igcx, get_igcc, get_inlc
  PUBLIC  :: dft_is_gradient, dft_is_meta, dft_is_hybrid, dft_is_nonlocc

  ! additional subroutines/functions for hybrid functionals
  PUBLIC  :: start_exx, stop_exx, get_exx_fraction, exx_is_active
  PUBLIC  :: set_exx_fraction
  PUBLIC  :: set_screening_parameter, get_screening_parameter
  PUBLIC  :: set_gau_parameter, get_gau_parameter

  ! additional subroutines/functions for finite size corrections
  PUBLIC  :: dft_has_finite_size_correction, set_finite_size_volume
  ! driver subroutines computing XC
  PUBLIC  :: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
  PUBLIC  :: tau_xc , tau_xc_spin, dmxc, dmxc_spin, dmxc_nc
  PUBLIC  :: dgcxc, dgcxc_spin
  PUBLIC  :: nlc
  ! general XC driver
  PUBLIC  :: vxc_t, exc_t
  ! vector XC driver
  PUBLIC  :: evxc_t_vec, gcx_spin_vec
  !
  ! PRIVATE variables defining the DFT functional
  !
  PRIVATE :: dft, dft_shortname, iexch, icorr, igcx, igcc, inlc
  PRIVATE :: discard_input_dft
  PRIVATE :: isgradient, ismeta, ishybrid
  PRIVATE :: exx_fraction, exx_started
  PRIVATE :: has_finite_size_correction, &
             finite_size_cell_volume,  finite_size_cell_volume_set 
  !
  character (len=25) :: dft = 'not set'
  character (len=6)  :: dft_shortname = ' '
  !
  ! dft is the exchange-correlation functional, described by
  ! one of the following keywords ("dft_shortname"):
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
  !              "tpss"  = "sla+pw+tpss+tpss"  = TPSS Meta-GGA
  !              "m06l"  = "nox+noc+m6lx+m6lc" = M06L Meta-GGA
  !              "pbe0"  = "pb0x+pw+pb0x+pbc"  = PBE0
  !              "hse"   = "sla+pw+hse+pbc"    = Heyd-Scuseria-Ernzerhof 
  !                                              (HSE 06, see note below)
  !              "b3lyp" = "b3lp+vwn+b3lp+b3lp"= B3LYP
  !              "vdw-df"= "sla+pw+rpb+vdw1"   = vdW-DF
  !              "vdw-df2"="sla+pw+rw86+vdw2"  = vdW-DF2
  !              "vdw-df-c09"="sla+pw+c09x+vdw1"
  !              "vdw-df2-c09"="sla+pw+c09x+vdw2"
  !              "vdw-df3"="sla+pw+obk8+vdw1"  = vdW-DF3
  !              "vdw-df4"="sla+pw+ob86+vdw1"  = vdW-DF4
  !              "optbk88"="sla+pw+obk8"       = optB88
  ! or by any nonconflicting combination of the following keywords
  ! (case-insensitive):
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
  !
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
  !              "rpb"    revised PBE by Zhang-Yang      igcx =4
  !              "hcth"   Cambridge exch, Handy et al    igcx =5
  !              "tpss"   TPSS meta-gga                  igcx =7
  !              "optx"   Handy's exchange functional    igcx =6
  !              "pb0x"   PBE0 (PBE exchange*0.75)       igcx =8
  !              "b3lp"   B3LYP (Becke88*0.72)           igcx =9
  !              "psx"    PBEsol exchange                igcx =10
  !              "wcx"    Wu-Cohen                       igcx =11
  !              "hse"    HSE screened exchange          igcx =12
  !              "rw86"   revised PW86                   igcx =13
  !              "pbe"    same as PBX, back-comp.        igcx =14
  !              "meta"   same as TPSS, back-comp.       igcx =15
  !              "c09x"   Cooper 09                      igcx =16
  !              "sox"    sogga                          igcx =17
  !              "m6lx"   M06L exchange Meta-GGA         igcx =18
  !              "q2dx"   Q2D exchange grad corr         igcx =19
  !              "gaup"   Gau-PBE hybrid exchange        igcx =20
  !              "pw86"   Perdew-Wang (1986) exchange    igcx =21
  !              "b86b"   Becke (1986) exchange          igcx =22
  !              "obk8"   optB88  exchange               igcx =23
  !              "ob86"   optB86b exchange               igcx =24
  !
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Lee-Yang-Parr                  igcc =3
  !              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
  !              "hcth"   Cambridge corr, Handy et al    igcc =5
  !              "tpss"   TPSS meta-gga                  igcc =6
  !              "b3lp"   B3LYP (Lee-Yang-Parr*0.81)     igcc =7
  !              "psc"    PBEsol corr                    igcc =8
  !              "pbe"    same as PBX, back-comp.        igcc =9
  !              "meta"   same as TPSS, back-comp.       igcc =10
  !              "m6lc"   M06L corr  Meta-GGA            igcc =11
  !              "q2dc"   Q2D correlation grad corr      igcc =12
  !
  ! Van der Waals functionals (nonlocal term only)
  !             "nonlc"   none                           inlc =0 (default)
  !              "vdw1"   vdW-DF1                        inlc =1
  !              "vdw2"   vdW-DF2                        inlc =2
  !              "vv10"   rVV10                          inlc =3  
  ! References:
  !              pz      J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) 
  !              vwn     S.H.Vosko, L.Wilk, M.Nusair, Can.J.Phys. 58,1200(1980)
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
  !              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  !              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
  !              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
  !              hcth    Handy et al, JCP 109, 6264 (1998)
  !              olyp    Handy et al, JCP 116, 5411 (2002)
  !              revPBE  Zhang and Yang, PRL 80, 890 (1998)
  !              pbesol  J.P. Perdew et al., PRL 100, 136406 (2008)
  !              q2d     L. Chiodo et al., PRL 108, 126402 (2012)
  !              rw86    E. Amonn D. Murray et al, J. Chem. Theory comp. 5, 2754 (2009) 
  !              wc      Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
  !              kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
  !              pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
  !              hse     Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
  !                      Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006).
  !              b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch
  !                      J.Phys.Chem 98, 11623 (1994)
  !              vdW-DF  M. Dion et al., PRL 92, 246401 (2004)
  !                      T. Thonhauser et al., PRB 76, 125112 (2007)
  !              vdw-DF2 Lee et al., Phys. Rev. B 82, 081101 (2010)
  !              vdw-DF3  Klimes et al, J. Phys. Cond. Matter, 22, 022201 (2010)
  !              vdw-DF4  Klimes et al, Phys. Rev. B, 83, 195131 (2011)
  !              c09x    V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
  !              tpss    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, 
  !                      PRL 91, 146401 (2003)
  !              sogga   Y. Zhao and D. G. Truhlar, JCP 128, 184109 (2008)
  !              m06l    Y. Zhao and D. G. Truhlar, JCP 125, 194101 (2006)
  !              gau-pbe J.-W. Song, K. Yamashita, K. Hirao JCP 135, 071103 (2011)
  !              rVV10   R. Sabatini et al. Phys. Rev. B 87, 041108(R) (2013)
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
  integer :: iexch = notset
  integer :: icorr = notset
  integer :: igcx  = notset
  integer :: igcc  = notset
  integer :: inlc  = notset
  real(DP):: exx_fraction = 0.0_DP
  real(DP):: screening_parameter = 0.0_DP
  real(DP):: gau_parameter = 0.0_DP
  logical :: isgradient  = .false.
  logical :: ismeta      = .false.
  logical :: ishybrid    = .false.
  logical :: exx_started = .false.
  logical :: has_finite_size_correction = .false.
  logical :: finite_size_cell_volume_set = .false.
  real(DP):: finite_size_cell_volume = notset
  
  logical :: isnonlocc       = .false.

  logical :: discard_input_dft = .false.
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !    inlc:  type of non local correction on correlation
  !
  !    ismeta: .TRUE. if gradient correction is of meta-gga type
  !    ishybrid: .TRUE. if the xc functional is an HF+DFT hybrid like
  !              PBE0, B3LYP, HSE or HF itself
  !
  ! see comments above and routine "set_dft_from_name" below 
  !
  ! data
  integer :: nxc, ncc, ngcx, ngcc, ncnl
  parameter (nxc = 8, ncc =10, ngcx =24, ngcc = 12, ncnl=3)
  character (len=4) :: exc, corr
  character (len=4) :: gradx, gradc, nonlocc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc), nonlocc (0: ncnl)

  data exc / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP', 'KZK' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' , 'KZK' /

  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX',&
               'TPSS', 'PB0X', 'B3LP','PSX', 'WCX', 'HSE', 'RW86', 'PBE', &
               'META', 'C09X', 'SOX', 'M6LX', 'Q2DX', 'GAUP', 'PW86', 'B86B', &
               'OBK8','OB86' / 

  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'TPSS',&
               'B3LP', 'PSC', 'PBE', 'META', 'M6LC', 'Q2DC' / 

  data nonlocc / '    ', 'VDW1', 'VDW2', 'VV10' / 

CONTAINS
  !-----------------------------------------------------------------------
  subroutine set_dft_from_name( dft_ )
    !-----------------------------------------------------------------------
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices iexch, icorr, igcx, igcc
    !
    implicit none
    ! input
    character(len=*)               :: dft_
    ! local
    integer :: len, l, i
    character (len=50):: dftout
    logical :: dft_defined = .false.
    logical, external :: matches
    character (len=1), external :: capital
    integer ::  save_iexch, save_icorr, save_igcx, save_igcc, save_inlc
    
    !
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
    ! FIRST WE CHECK ALL THE SPECIAL NAMES
    ! Note: comparison is now done via exact matching
    !       not using function "matches"
    ! ----------------------------------------------
    !

    if ( 'REVPBE' .EQ. TRIM(dftout) ) then
    ! special case : revPBE
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 4)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc, 0)
       dft_defined = .true.
       
    elseif ('PW86PBE' .EQ. TRIM(dftout) ) then
       ! special case : PW86PBE 
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 21)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc, 0)
       dft_defined = .true.

    elseif ('B86BPBE' .EQ. TRIM(dftout) ) then
       ! special case : B86BPBE
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 22)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc, 0)
       dft_defined = .true.

    else if ('RPBE' .EQ. TRIM(dftout)) then
    ! special case : RPBE
         call errore('set_dft_from_name', &
     &   'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)',1)
     
    else if ('PBE0'.EQ. TRIM(dftout) ) then
    ! special case : PBE0
       call set_dft_value (iexch,6)
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 8)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.
       
   else if ('HSE' .EQ. TRIM( dftout) ) then
    ! special case : HSE
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 12)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.

   else if (matches ('GAUP', dftout) ) then
    ! special case : GAUPBE
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 20)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) !Default
       dft_defined = .true.
       
    else if ('PBESOL'.EQ. TRIM(dftout) ) then
    ! special case : PBEsol
       call set_dft_value (iexch,1) !Default    
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,10)
       call set_dft_value (igcc, 8)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.

    else if ('RVV10' .EQ. TRIM(dftout) ) then
    ! Special case rVV10
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 13)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc, 3)
       dft_defined = .true.

    else if ('PBEQ2D' .EQ. TRIM(dftout) .OR. 'Q2D'.EQ. TRIM(dftout) ) then
    ! special case : PBEQ2D
       call set_dft_value (iexch,1) !Default    
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,19)
       call set_dft_value (igcc,12)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.

    else if ('VDW-DF4' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF4, or optB86b+vdW
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 24)
       call set_dft_value (igcc, 0)
       call set_dft_value (inlc, 1)       
       dft_defined = .true.
       
    else if ('VDW-DF3' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF3, or optB88+vdW
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 23)
       call set_dft_value (igcc, 0)
       call set_dft_value (inlc, 1)       
       dft_defined = .true.
       
    else if ('OPTBK88' .EQ. TRIM(dftout)) then
    ! Special case optB88 (without vdW)
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 23)
       call set_dft_value (igcc, 1)
       call set_dft_value (inlc, 0)       
       dft_defined = .true.
       
    else if ('OPTB86B' .EQ. TRIM(dftout)) then
    ! Special case optB86b (without vdW)
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 24)
       call set_dft_value (igcc, 1)
       call set_dft_value (inlc, 0)       
       dft_defined = .true.
       
    else if ('VDW-DF2-C09' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF2 with C09 exchange
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 16)
       call set_dft_value (igcc, 0)
       call set_dft_value (inlc, 2)
       dft_defined = .true.

    else if ('VDW-DF-C09'  .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF with C09 exchange
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 16)
       call set_dft_value (igcc, 0)
       call set_dft_value (inlc, 1)
       dft_defined = .true.

   else if ('VDW-DF2' .EQ. TRIM(dftout) ) then
    ! Special case vdW-DF2
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 13)
       call set_dft_value (igcc, 0)
       call set_dft_value (inlc, 2)
       dft_defined = .true.

       
    else if ('VDW-DF' .EQ. TRIM(dftout)) then
    ! Special case vdW-DF
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
       call set_dft_value (igcx, 4)
       call set_dft_value (igcc, 0)
       call set_dft_value (inlc, 1)       
       dft_defined = .true.

    else if ('PBE' .EQ. TRIM(dftout) ) then
    ! special case : PBE
       call set_dft_value (iexch,1) !Default    
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 3)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.
       
    else if ('WC' .EQ. TRIM(dftout) ) then
    ! special case : Wu-Cohen
       call set_dft_value (iexch,1) !Default       
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,11)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.
       
    else if ('B3LYP'.EQ. TRIM(dftout) ) then
    ! special case : B3LYP hybrid
       call set_dft_value (iexch,7)
       call set_dft_value (icorr,2)
       call set_dft_value (igcx, 9)
       call set_dft_value (igcc, 7)
       call set_dft_value (inlc,0) !Default       
       dft_defined = .true.
       
    else if ('PBC'.EQ. TRIM(dftout) ) then
    ! special case : PBC  = PW + PBC 
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,0) !Default       
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.
       
    ! special case : BP = B88 + P86
    else if ('BP'.EQ. TRIM(dftout) ) then
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,1) !Default
       call set_dft_value (igcx, 1)
       call set_dft_value (igcc, 1)
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.
                     
    ! special case : PW91 = GGX + GGC
    else if ('PW91'.EQ. TRIM(dftout) ) then
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 2)
       call set_dft_value (igcc, 2)
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.
       
    ! special case : HCTH
    else if ('HCTH'.EQ. TRIM(dftout)) then
       call set_dft_value(iexch,0) ! contained in hcth
       call set_dft_value(icorr,0) ! contained in hcth
       call set_dft_value (igcx,5)
       call set_dft_value (igcc,5)
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.
              
    ! special case : OLYP = OPTX + LYP
    else if ('OLYP'.EQ. TRIM(dftout)) then
       call set_dft_value(iexch,0) ! contained in optx
       call set_dft_value(icorr,3)
       call set_dft_value(igcx, 6)
       call set_dft_value(igcc, 3)
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.
       
    ! special case : TPSS meta-GGA Exc
    else IF ('TPSS'.EQ. TRIM(dftout ) ) THEN
       CALL set_dft_value( iexch, 1 )
       CALL set_dft_value( icorr, 4 )
       CALL set_dft_value( igcx,  7 )
       CALL set_dft_value( igcc,  6 )
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.
              
    ! special cases : OEP no GC part (nor LDA...) and no correlation by default
    else IF ('OEP' .EQ. TRIM(dftout) ) THEN
       call set_dft_value (iexch,4) 
       call set_dft_value (icorr, 0)
       CALL set_dft_value( igcx,  0 )
       call set_dft_value (igcc, 0) !Default       
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.

    ! special cases : HF no GC part (nor LDA...) and no correlation by default
    else IF ('HF' .EQ. TRIM(dftout) ) THEN
       call set_dft_value (iexch,5) 
       call set_dft_value (icorr, 0)
       CALL set_dft_value( igcx,  0 )
       call set_dft_value (igcc, 0) !Default       
       call set_dft_value (inlc,0) !Default    
       dft_defined = .true.

    ! special cases : BLYP (note, BLYP=>B88)
    else IF ('BLYP' .EQ. TRIM(dftout) ) THEN
       call set_dft_value (iexch,1) !Default
       call set_dft_value (icorr,3)
       CALL set_dft_value( igcx, 1 )
       call set_dft_value (igcc, 3)
       call set_dft_value (inlc, 0) !Default    
       dft_defined = .true.

    ! special cases : PZ  (LDA is equivalent to PZ)
    else IF (('PZ' .EQ. TRIM(dftout) ).OR.('LDA' .EQ. TRIM(dftout) )) THEN
       call set_dft_value (iexch,1) 
       call set_dft_value (icorr, 1) 
       CALL set_dft_value( igcx,  0)
       call set_dft_value (igcc, 0)      
       call set_dft_value (inlc,0)    
       dft_defined = .true.

    ! special case : SOGGA = SOX + PBEc       
    else if (matches ('SOGGA', dftout) ) then
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,17)
       call set_dft_value (igcc, 4)
       call set_dft_value (inlc,0) ! Default    
       dft_defined = .true.
     
    ! special case : M06L Meta GGA
    else if ( matches( 'M06L', dftout ) ) THEN
       !
       CALL set_dft_value( iexch, 0 ) ! contained in m6lx
       CALL set_dft_value( icorr, 0 ) ! contained in m6lc
       CALL set_dft_value( igcx,  18 )
       CALL set_dft_value( igcc,  11)
       call set_dft_value (inlc,0) ! Default  
       dft_defined = .true.
    
    END IF

    !
    ! ----------------------------------------------------------------
    ! If the DFT was not yet defined, check every part of the string
    ! ----------------------------------------------------------------
    !
    if (.not. dft_defined) then
    
      ! write(*,"(A,A)") "Setting by parts: ", TRIM(dftout)      

      !  exchange
      iexch = notset
      do i = 0, nxc
         if (matches (exc (i), dftout) ) call set_dft_value (iexch, i)
      enddo
      if (iexch .eq. notset) call set_dft_value (iexch,0)

      !  correlation
      icorr = notset
      do i = 0, ncc
         if (matches (corr (i), dftout) ) call set_dft_value (icorr, i)
      enddo
      if (icorr .eq. notset) call set_dft_value (icorr,0)

      !  gradient correction, exchange
      igcx = notset
      do i = 0, ngcx
         if (matches (gradx (i), dftout) ) call set_dft_value (igcx, i)
      enddo
      if (igcx .eq. notset) call set_dft_value (igcx,0)
    
      !  gradient correction, correlation
      igcc = notset
      do i = 0, ngcc
         if (matches (gradc (i), dftout) ) call set_dft_value (igcc, i)
      enddo
      if (igcc .eq. notset) call set_dft_value (igcc,0)
    
      !  non-local correlation
      !     THE LOOP IS REVERSED TO HANDLE THE VDW2 CASE BEFORE THE VDW
      inlc = notset
      do i = ncnl ,1, -1
         if (matches (nonlocc (i), dftout) ) call set_dft_value (inlc, i)
      enddo
      if (inlc .eq. notset) call set_dft_value (inlc,0)
        
    endif

    ! ----------------------------------------------------------------
    ! Last check
    ! No more defaults, the code exit if the dft is not defined
    ! ----------------------------------------------------------------

    ! Back compatibility - TO BE REMOVED

    if (igcx == 13 .and. iexch /= 1) &
         call errore('set_dft_from_name','rPW86 no longer contains Slater exchange, add it explicitly',-igcx)
 
    if (igcx == 14) igcx = 3 ! PBE -> PBX
    if (igcc == 9) igcc = 4  ! PBE -> PBC

    if (igcx == 15) igcx = 7 ! TPSS -> META
    if (igcc == 10) igcc = 6 ! TPSS -> META

    if (igcx == 6) &
         call errore('set_dft_from_name','OPTX untested! please test',-igcx)
         
    if (iexch <=0 .and. &
       icorr <=0 .and. &
       igcx <= 0 .and. &
       igcc <= 0 .and. &
       inlc <= 0) &
           call errore('set_dft_from_name','No dft definition was found',0)

    !
    ! Fill variables and exit
    !
    dft = dftout

    dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
         &//gradc (igcc) //'-'// nonlocc(inlc)


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
    if (save_inlc .ne. notset .and. save_inlc .ne. inlc) then
       write (stdout,*) inlc, save_inlc
       call errore('set_dft_from_name',' conflicting values for inlc',1)
    end if

    return
  end subroutine set_dft_from_name
  !
  !-----------------------------------------------------------------------
  subroutine set_auxiliary_flags
    !-----------------------------------------------------------------------
    ! set logical flags describing the complexity of the xc functional
    ! define the fraction of exact exchange used by hybrid fuctionals
    !
    logical, external :: matches

    !! Reversed as before VDW
    isgradient =  ( (igcx > 0) .or. ( igcc > 0) )

    isnonlocc = (inlc > 0)

    ismeta     =  (igcx == 7) .or. (igcx == 18)

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
    !B3LYP
    IF ( matches( 'B3LP',dft ) .OR. matches( 'B3LYP',dft ) ) &
                                  exx_fraction = 0.2_DP
    ishybrid = ( exx_fraction /= 0.0_DP )

    has_finite_size_correction = ( iexch==8 .or. icorr==10)

    return
  end subroutine set_auxiliary_flags
  !
  !-----------------------------------------------------------------------
  subroutine set_dft_value (m, i)
    !-----------------------------------------------------------------------
    !
    implicit none
    integer :: m, i
    ! local

    if ( m /= notset .and. m /= i) then
         write(*, '(A,2I4)') "parameters", m, i
         call errore ('set_dft_value', 'two conflicting matching values', 1)
    end if
    m = i
    return

  end subroutine set_dft_value

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
  !-----------------------------------------------------------------------
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
  function get_inlc ()
     integer get_inlc
     get_inlc = inlc
     return
  end function get_inlc
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
  !---------------------------------------------------------------------
  subroutine dft_name(iexch_, icorr_, igcx_, igcc_, inlc_, longname_, shortname_)
  !---------------------------------------------------------------------
  ! convert the four indices iexch, icorr, igcx, igcc
  ! into user-readable strings
  !
  implicit none
  integer iexch_, icorr_, igcx_, igcc_, inlc_
  character (len=6) :: shortname_
  character (len=25):: longname_
  !
  if (iexch_==1.and.igcx_==0.and.igcc_==0) then
     shortname_ = corr(icorr_)
  else if (iexch_==1.and.icorr_==3.and.igcx_==1.and.igcc_==3) then
     shortname_ = 'BLYP'
  else if (iexch_==1.and.icorr_==1.and.igcx_==1.and.igcc_==0) then
     shortname_ = 'B88'
  else if (iexch_==1.and.icorr_==1.and.igcx_==1.and.igcc_==1) then
     shortname_ = 'BP'
  else if (iexch_==1.and.icorr_==4.and.igcx_==2.and.igcc_==2) then
     shortname_ = 'PW91'
  else if (iexch_==1.and.icorr_==4.and.igcx_==3.and.igcc_==4) then
     shortname_ = 'PBE'
  else if (iexch_==6.and.icorr_==4.and.igcx_==8.and.igcc_==4) then
     shortname_ = 'PBE0'
  else if (iexch_==1.and.icorr_==4.and.igcx_==4.and.igcc_==4) then
     shortname_ = 'revPBE'
  else if (iexch_==1.and.icorr_==4.and.igcx_==10.and.igcc_==8) then
     shortname_ = 'PBESOL'
  else if (iexch_==1.and.icorr_==4.and.igcx_==19.and.igcc_==12) then
     shortname_ = 'Q2D'
  else if (iexch_==1.and.icorr_==4.and.igcx_==12.and.igcc_==4) then
     shortname_ = 'HSE'
  else if (iexch_==1.and.icorr_==4.and.igcx_==20.and.igcc_==4) then
     shortname_ = 'GAUPBE'
  else if (iexch_==1.and.icorr_==4.and.igcx_==11.and.igcc_==4) then
     shortname_ = 'WC'
  else if (iexch_==7.and.(icorr_==10.or.icorr_==2).and.igcx_==9.and. &
           igcc_==7) then
     shortname_ = 'B3LYP'
  else if (iexch_==0.and.icorr_==3.and.igcx_==6.and.igcc_==3) then
     shortname_ = 'OLYP'
  else if (iexch_==1.and.icorr_==4.and.igcx_==4.and.igcc_==0.and.inlc_==1) then
     shortname_ = 'VDW-DF'
  else if (iexch_==1.and.icorr_==4.and.igcx_==13.and.igcc_==0.and.inlc_==2) then
     shortname_ = 'VDW-DF2'
  else if (iexch_==1.and.icorr_==4.and.igcx_==16.and.igcc_==0.and.inlc_==1) then
     shortname_ = 'VDW-DF-C09'
  else if (iexch_==1.and.icorr_==4.and.igcx_==16.and.igcc_==0.and.inlc_==2) then
     shortname_ = 'VDW-DF2-C09'
  else if (iexch_==1.and.icorr_==4.and.igcx_==13.and.igcc_==4.and.inlc_==3) then
     shortname_ = 'RVV10'
  else if (iexch_==1.and.icorr_==4.and.igcx_==24.and.igcc_==0.and.inlc_==1) then
     shortname_ = 'VDW-DF4'
  else if (iexch_==1.and.icorr_==4.and.igcx_==23.and.igcc_==0.and.inlc_==1) then
     shortname_ = 'VDW-DF3'
  else if (iexch_==0.and.icorr_==0.and.igcx_==18.and.igcc_==11) then
     shortname_ = 'M06L'
  else if (iexch_==1.and.icorr_==4.and.igcx_==17.and.igcc_==4) then
     shortname_ = 'SOGGA'
  else
     shortname_ = ' '
  end if
  write(longname_,'(5a5)') exc(iexch_),corr(icorr_),gradx(igcx_),gradc(igcc_),nonlocc(inlc_)
  
  return
end subroutine dft_name

subroutine write_dft_name
!-----------------------------------------------------------------------
   WRITE( stdout, '(5X,"Exchange-correlation      = ",A, &
        &  " (",I2,3I3,I2,")")') TRIM( dft ), iexch, icorr, igcx, igcc, inlc
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
  ELSEIF (iexch == 7) THEN         !  'b3lyp'
     CALL slater(rs, ex, vx)
     if (exx_started) then
        ex = 0.8_DP * ex 
        vx = 0.8_DP * vx 
     end if
  ELSEIF (iexch == 8) THEN         !  'sla+kzk'
     if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
          'finite size corrected exchange used w/o initialization',1)
     call slaterKZK (rs, ex, vx, finite_size_cell_volume)
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
  ELSEIF (iexch == 7) THEN  ! 'b3lyp'
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
  case(7)            ! 'b3lyp'
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
  elseif (igcx == 8) then ! 'pbe0'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if (exx_started) then
        sx  = (1.0_DP - exx_fraction) * sx
        v1x = (1.0_DP - exx_fraction) * v1x
        v2x = (1.0_DP - exx_fraction) * v2x
     end if
  elseif (igcx == 9) then ! 'b3lyp'
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
  elseif (igcx == 23) then ! 'optB88'
     call pbex (rho, grho, 5, sx, v1x, v2x)
  elseif (igcx == 24) then ! 'optB86b'
     call pbex (rho, grho, 6, sx, v1x, v2x)
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
  elseif (igcx == 3 .or. igcx == 4 .or. igcx == 8 .or. &
          igcx == 10 .or. igcx == 12 .or. igcx == 20) then
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8: PBE0, igcx=10: PBEsol
     ! igcx=12: HSE,  igcx=20: gau-pbe
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
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
  case(3,4,8,10,12)
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8 PBE0, igcx=10: PBEsol
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
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

  USE vdW_DF, ONLY: xc_vdW_DF, vdw_type
  USE rVV10,  ONLY: xc_rVV10
 
  implicit none
  
  REAL(DP), INTENT(IN) :: rho_valence(:,:), rho_core(:)
  INTEGER, INTENT(IN)  :: nspin
  REAL(DP), INTENT(INOUT) :: v(:,:)
  REAL(DP), INTENT(INOUT) :: enl, vnl

  if (inlc == 1 .or. inlc == 2) then
     
     vdw_type = inlc
     call xc_vdW_DF(rho_valence, rho_core, nspin, enl, vnl, v)

  elseif (inlc == 3) then

      call xc_rVV10(rho_valence, rho_core, nspin, enl, vnl, v)
  
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
  
  if     (igcx == 7  .and. igcc == 6) then
  
     call tpsscxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
     
  elseif (igcx == 18 .and. igcc == 11) then
  
     call   m06lxc (rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)
     
  else
  
     call errore('v_xc_meta','wrong igcx and/or igcc',1)
     
  end if
  
  return
  
end subroutine tau_xc

!
!
!-----------------------------------------------------------------------
subroutine tau_xc_spin (rhoup, rhodw, grhoup, grhodw, tauup, taudw, ex, ec,       &
           &            v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw, v1cup, v1cdw,   &
           &            v2cup, v2cdw, v2cup_vec, v2cdw_vec, v3cup, v3cdw)

!-----------------------------------------------------------------------
  !
  !
  
  implicit none

  real(dp), intent(in)                :: rhoup, rhodw, tauup, taudw
  real(dp), dimension (3), intent(in) :: grhoup, grhodw
  
  real(dp), intent(out)               :: ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
                                      &  v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw
  real(dp), dimension(3), intent(out) :: v2cup_vec, v2cdw_vec
  
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
  v2cup_vec (:) = zero
  v2cdw_vec (:) = zero
  
  
  do ipol=1,3
     grhoup2 = grhoup2 + grhoup(ipol)**2
     grhodw2 = grhodw2 + grhodw(ipol)**2
  end do

  
  if (igcx == 7 .and. igcc == 6) then

     call tpsscx_spin(rhoup, rhodw, grhoup2, grhodw2, tauup,   &
              &  taudw, ex, v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
  
     rh   =  rhoup + rhodw
        
     zeta = (rhoup - rhodw) / rh
     atau =  tauup + taudw    ! KE-density in Hartree

     call tpsscc_spin(rh,zeta,grhoup,grhodw, atau,ec,              &
     &                v1cup,v1cdw,v2cup_vec,v2cdw_vec,v3cup, v3cdw) 
  
  
  elseif (igcx == 18 .and. igcc == 11) then
  
     call   m06lxc_spin (rhoup, rhodw, grhoup2, grhodw2, tauup, taudw,      &
            &            ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
            &            v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw)
     
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

!
!-----------------------------------------------------------------------
!------- VECTOR AND GENERAL XC DRIVERS -------------------------------
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------
subroutine vxc_t(rho,rhoc,lsd,vxc)
  !---------------------------------------------------------------
  !
  !  this function returns the XC potential in LDA or LSDA approximation
  !
  use io_global, only : stdout
  use kinds, only : DP
  implicit none
  integer:: lsd
  real(DP):: vxc(2), rho(2),rhoc,arho,zeta
  real(DP):: vx(2), vc(2), ex, ec
  !
  real(DP), parameter :: e2=2.0_dp, eps=1.e-30_dp

  vxc(1)=0.0_dp
  if (lsd.eq.1) vxc(2)=0.0_dp

  if (lsd.eq.0) then
     !
     !     LDA case
     !
     arho=abs(rho(1)+rhoc)
     if (arho.gt.eps) then      
        call xc(arho,ex,ec,vx(1),vc(1))
        vxc(1)=e2*(vx(1)+vc(1))
     endif
  else
     !
     !     LSDA case
     !
     arho = abs(rho(1)+rho(2)+rhoc)
     if (arho.gt.eps) then      
        zeta = (rho(1)-rho(2)) / arho
        ! zeta has to stay between -1 and 1, but can get a little
        ! out the bound during the first iterations.
        if (abs(zeta).gt.1.0_dp) zeta = sign(1._dp, zeta)
        call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
        vxc(1) = e2*(vx(1)+vc(1))
        vxc(2) = e2*(vx(2)+vc(2))
     endif
  endif

  return
end subroutine vxc_t


!---------------------------------------------------------------
function exc_t(rho,rhoc,lsd)
  !---------------------------------------------------------------
  !
  use kinds, only : DP
  implicit none
  integer:: lsd
  real(DP) :: exc_t, rho(2),arho,rhot, zeta,rhoc
  real(DP) :: ex, ec, vx(2), vc(2)

  real(DP),parameter:: e2 =2.0_DP

  exc_t=0.0_DP

  if(lsd == 0) then
     !
     !     LDA case
     !
     rhot = rho(1) + rhoc
     arho = abs(rhot)
     if (arho.gt.1.e-30_DP) then      
        call xc(arho,ex,ec,vx(1),vc(1))
        exc_t=e2*(ex+ec)
     endif
  else
     !
     !     LSDA case
     !
     rhot = rho(1)+rho(2)+rhoc
     arho = abs(rhot)
     if (arho.gt.1.e-30_DP) then      
        zeta = (rho(1)-rho(2)) / arho
        ! In atomic this cannot happen, but in PAW zeta can become
        ! a little larger than 1, or smaller than -1:
        if( abs(zeta) > 1._dp) zeta = sign(1._dp, zeta)
        call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
        exc_t=e2*(ex+ec)
     endif
  endif

  return
end function exc_t

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
