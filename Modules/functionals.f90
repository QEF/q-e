!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
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
!  retrive functions:  get_dft_name
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
  PUBLIC  :: get_dft_name, get_iexch, get_icorr, get_igcx, get_igcc
  PUBLIC  :: dft_is_gradient, dft_is_meta, dft_is_hybrid
  ! additional subroutines/functions for hybrid functionals
  PUBLIC  :: start_exx, stop_exx, get_exx_fraction, exx_is_active
  PUBLIC  :: set_exx_fraction
  PUBLIC  :: set_screening_parameter, get_screening_parameter
  ! additional subroutines/functions for finite size corrections
  PUBLIC  :: dft_has_finite_size_correction, set_finite_size_volume
  ! driver subroutines computing XC
  PUBLIC  :: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
  PUBLIC  :: dmxc, dmxc_spin, dmxc_nc
  PUBLIC  :: dgcxc, dgcxc_spin
  ! general XC driver
  PUBLIC  :: vxc_t, exc_t
  ! vector XC driver
  PUBLIC  :: evxc_t_vec, gcx_spin_vec
  !
  ! PRIVATE variables defining the DFT functional
  !
  PRIVATE :: dft, dft_shortname, iexch, icorr, igcx, igcc
  PRIVATE :: discard_input_dft
  PRIVATE :: isgradient, ismeta, ishybrid
  PRIVATE :: exx_fraction, exx_started
  PRIVATE :: has_finite_size_correction, &
             finite_size_cell_volume,  finite_size_cell_volume_set 
  !
  character (len=20) :: dft = 'not set'
  character (len=6)  :: dft_shortname = ' '
  !
  ! dft is the exchange-correlation functional, described by
  ! any nonconflicting combination of the following keywords
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
  !              "b3lp"   B3LYP (same as "vwn")          icorr=10
  !              "kzk"    Finite-size corrections        icorr=11
  !
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
  !              "rpb"    revised PBE by Zhang-Yang      igcx =4
  !              "hcth"   Cambridge exch, Handy et al    igcx =5
  !              "optx"   Handy's exchange functional    igcx =6
  !              "meta"   TPSS meta-gga                  igcx =7
  !              "pb0x"   PBE0 (PBE exchange*0.75)       igcx =8
  !              "b3lp"   B3LYP (Becke88*0.72)           igcx =9
  !              "psx"    PBEsol exchange                igcx =10
  !              "wcx"    Wu-Cohen                       igcx =11
  !
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Lee-Yang-Parr                  igcc =3
  !              "pbc"    Perdew-Burke-Ernzenhof corr    igcc =4
  !              "hcth"   Cambridge corr, Handy et al    igcc =5
  !              "meta"   TPSS meta-gga                  igcc =6
  !              "b3lp"   B3LYP (Lee-Yang-Parr*0.81)     igcc =7
  !              "psc"    PBEsol corr                    igcc =8
  !
  ! Special cases (dft_shortname):
  !              "bp"    = "b88+p86"           = Becke-Perdew grad.corr.
  !              "pw91"  = "pw +ggx+ggc"       = PW91 (aka GGA)
  !              "blyp"  = "sla+b88+lyp+blyp"  = BLYP
  !              "pbe"   = "sla+pw+pbx+pbc"    = PBE
  !              "revpbe"="sla+pw+rpb+pbc"     = revPBE (Zhang-Yang)
  !              "pbesol"="sla+pw+psx+psc"     = PBEsol
  !              "hcth"  = "nox+noc+hcth+hcth" = HCTH/120
  !              "olyp"  = "nox+lyp+optx+blyp" = OLYP
  !              "tpss"  = "sla+pw+meta+meta"  = TPSS Meta-GGA
  !              "wc"    = "sla+pw+wcx+pbc"    = Wu-Cohen
  !              "pbe0"  = "pb0x+pw+pb0x+pbc"  = PBE0
  !              "hse"   = "moh "              = HSE
  !              "b3lyp" = "b3lp+vwn+b3lp+b3lp"= B3LYP
  !
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
  !              pbe     J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  !              pw91    J.P.Perdew and Y. Wang, PRB 46, 6671 (1992)
  !              blyp    C.Lee, W.Yang, R.G.Parr, PRB 37, 785 (1988)
  !              hcth    Handy et al, JCP 109, 6264 (1998)
  !              olyp    Handy et al, JCP 116, 5411 (2002)
  !              revPBE  Zhang and Yang, PRL 80, 890 (1998)
  !              meta    J.Tao, J.P.Perdew, V.N.Staroverov, G.E. Scuseria, 
  !                      PRL 91, 146401 (2003)
  !              kzk     H.Kwee, S. Zhang, H. Krakauer, PRL 100, 126404 (2008)
  !              pbe0    J.P.Perdew, M. Ernzerhof, K.Burke, JCP 105, 9982 (1996)
  !              hse     Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 118, 8207 (2003)
  !                      Heyd, Scuseria, Ernzerhof, J. Chem. Phys. 124, 219906 (2006).
  !              b3lyp   P.J. Stephens,F.J. Devlin,C.F. Chabalowski,M.J. Frisch
  !                      J.Phys.Chem 98, 11623 (1994)
  !              pbesol  J.P. Perdew et al., PRL 100, 136406 (2008)
  !              wc      Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
  !
  integer, parameter:: notset = -1
  !
  integer :: iexch = notset
  integer :: icorr = notset
  integer :: igcx  = notset
  integer :: igcc  = notset
  real(DP):: exx_fraction = 0.0_DP
  real(DP):: screening_parameter = 0.0_DP
  logical :: isgradient  = .false.
  logical :: ismeta      = .false.
  logical :: ishybrid    = .false.
  logical :: exx_started = .false.
  logical :: has_finite_size_correction = .false.
  logical :: finite_size_cell_volume_set = .false.
  real(DP):: finite_size_cell_volume = notset

  logical :: discard_input_dft = .false.
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlation
  !
  !    ismeta: .TRUE. if gradient correction is of meta-gga type
  !    ishybrid: .TRUE. if the xc functional is an HF+DFT hybrid like
  !              PBE0, B3LYP, HSE or HF itself
  !
  ! see comments above and routine "set_dft_from_name" below 
  !
  ! data
  integer :: nxc, ncc, ngcx, ngcc
  parameter (nxc = 8, ncc =11, ngcx =12, ngcc = 8)
  character (len=4) :: exc, corr
  character (len=4) :: gradx, gradc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)

  data exc / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP', 'KZK' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' , 'B3LP', 'KZK' /
  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX',&
               'META', 'PB0X', 'B3LP','PSX', 'WCX', 'HSE'  /
  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'META',&
                'B3LP', 'PSC' /

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
    logical, external :: matches
    character (len=1), external :: capital
    !
    !
    ! if 
    !
    if ( discard_input_dft ) return
    !
    ! convert to uppercase
    len = len_trim(dft_)
    dftout = ' '
    do l = 1, len
       dftout (l:l) = capital (dft_(l:l) )
    enddo

    !  exchange
    iexch = notset
    do i = 0, nxc
       if (matches (exc (i), dftout) ) call set_dft_value (iexch, i)
    enddo

    !  correlation
    icorr = notset
    do i = 0, ncc
       if (matches (corr (i), dftout) ) call set_dft_value (icorr, i)
    enddo

    !  gradient correction, exchange
    igcx = notset
    do i = 0, ngcx
       if (matches (gradx (i), dftout) ) call set_dft_value (igcx, i)
    enddo

    !  gradient correction, correlation
    igcc = notset
    do i = 0, ngcc
       if (matches (gradc (i), dftout) ) call set_dft_value (igcc, i)
    enddo

    ! special case : BLYP => B88 for gradient correction on exchange
    ! warning: keyword BLYP is used for both the XC functional "BLYP"
    !          and for Lee-Yang-Parr gradient correction to correlation
    !          in the former case, iexch and igcx shouldn't have been set
    if (matches('BLYP', dftout) .and. (iexch == notset .and. igcx == notset)) &
         call set_dft_value (igcx, 1)
    ! special case : various variants of PBE 
    ! As routine matches returns .true. when the first string is contained 
    ! in the second one, all tests on functionals containing PBE as substring
    ! must preceed the test on PBE itself.
    if (matches ('REVPBE', dftout) ) then
    ! special case : revPBE
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 4)
       call set_dft_value (igcc, 4)
    else if (matches('RPBE',dftout)) then
    ! special case : RPBE
         call errore('set_dft_from_name', &
     &   'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)',1)
    else if (matches ('PBE0', dftout) ) then
    ! special case : PBE0
       call set_dft_value (iexch,6)
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 8)
       call set_dft_value (igcc, 4)
   else if (matches ('HSE', dftout) ) then
    ! special case : HSE
!       call set_dft_value (iexch,6)
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 12)
       call set_dft_value (igcc, 4)
    else if (matches ('PBESOL', dftout) ) then
    ! special case : PBEsol
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,10)
       call set_dft_value (igcc, 8)
    else if (matches ('PBE', dftout) ) then
    ! special case : PBE
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 3)
       call set_dft_value (igcc, 4)
    else if (matches ('WC', dftout) ) then
    ! special case : Wu-Cohen
       call set_dft_value (icorr,4)
       call set_dft_value (igcx,11)
       call set_dft_value (igcc, 4)
    else if (matches ('B3LYP', dftout) ) then
    ! special case : B3LYP hybrid
       call set_dft_value (iexch,7)
       !!! cannot use set_dft_value due to conflict with blyp
       icorr = 2
       call set_dft_value (igcx, 9)
       !!! as above
       igcc = 7
    endif

    if (matches ('PBC', dftout) ) then
    ! special case : PBC  = PW + PBC 
       call set_dft_value (icorr,4)
       call set_dft_value (igcc, 4)
    endif

    ! special case : BP = B88 + P86
    if (matches ('BP', dftout) ) then
       call set_dft_value (igcx, 1)
       call set_dft_value (igcc, 1)
    endif

    ! special case : PW91 = GGX + GGC
    if (matches ('PW91', dftout) ) then
       call set_dft_value (igcx, 2)
       call set_dft_value (igcc, 2)
    endif

    ! special case : HCTH already contains LDA exchange and correlation

    if (matches('HCTH',dftout)) then
       call set_dft_value(iexch,0)
       call set_dft_value(icorr,0)
    end if

    ! special case : OPTX already contains LDA exchange
     
    if (matches('OPTX',dftout)) then
       call set_dft_value(iexch,0)
    end if

    ! special case : OLYP = OPTX + LYP

    if (matches('OLYP',dftout)) then
       call set_dft_value(iexch,0)
       call set_dft_value(icorr,3)
       call set_dft_value(igcx,6)
       call set_dft_value(igcc,3)
    end if

    !
    ! ... special case : TPSS meta-GGA Exc
    !
    IF ( matches( 'TPSS', dftout ) ) THEN
       !
       CALL set_dft_value( iexch, 1 )
       CALL set_dft_value( icorr, 4 )
       CALL set_dft_value( igcx,  7 )
       CALL set_dft_value( igcc,  6 )
       !
    END IF
    !
    ! ... special cases : OEP and HF need not GC part (nor LDA...)
    !                     and include no correlation by default
    !
    IF ( matches( 'OEP', dftout ) .OR. matches( 'HF', dftout )) THEN
       !
       CALL set_dft_value( igcx,  0 )
       if (icorr == notset) call set_dft_value (icorr, 0)
       !
    END IF


    if (igcx == 6) &
         call errore('set_dft_from_name','OPTX untested! please test',-igcx)
    ! Default value: Slater exchange
    if (iexch == notset) call set_dft_value (iexch, 1)

    ! Default value: Perdew-Zunger correlation
    if (icorr == notset) call set_dft_value (icorr, 1)

    ! Default value: no gradient correction on exchange
    if (igcx == notset) call set_dft_value (igcx, 0)

    ! Default value: no gradient correction on correlation
    if (igcc == notset) call set_dft_value (igcc, 0)

    dft = dftout

    dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
         &//gradc (igcc)
    ! WRITE( stdout,'(a)') dftout

    call set_auxiliary_flags

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

    isgradient =  (igcx > 0) .or. (igcc > 0) 
    ismeta     =  (igcx == 7)

    ! PBE0
    IF ( iexch==6 .or. igcx ==8 ) exx_fraction = 0.25_DP
    ! HSE
    IF ( igcx ==12 ) THEN
       exx_fraction = 0.25_DP
       screening_parameter = 0.106_DP
    END IF
    ! HF or OEP
    IF ( iexch==4 .or. iexch==5 ) exx_fraction = 1.0_DP
    !B3LYP
    IF ( matches( 'B3LP',dft ) .OR. matches( 'B3LYP',dft ) ) &
                                  exx_fraction = 0.2_DP
    ishybrid = ( exx_fraction /= 0.0_DP )

    has_finite_size_correction = ( iexch==8 .or. icorr==11)

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

    if ( m /= notset .and. m /= i) &
         call errore ('set_dft_value', 'two conflicting matching values', 1)
    m = i
    return

  end subroutine set_dft_value

  !-----------------------------------------------------------------------
  subroutine enforce_input_dft (dft_)
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices and force any subsequent call to set_dft_from_name
    ! to return without changing them
    !
    implicit none
    ! input
    character(len=*) :: dft_
    ! data

     call set_dft_from_name (dft_)
     if (dft == 'not set') call errore('enforce_input_dft','cannot fix unset dft',1)
     discard_input_dft = .true.

     write (stdout,'(/,5x,a)') "!!! XC functional enforced from input :"
     call write_dft_name
     write (stdout,'(5x,a)') "!!! Any further DFT definition will be discarded"
     write (stdout,'(5x,a)') "!!! Please, verify this is what you really want !"

     return
  end subroutine enforce_input_dft
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
     write (stdout,'(5x,a,f)') 'EXX fraction changed: ',exx_fraction
  end subroutine set_exx_fraction
  !---------------------------------------------------------------------
  subroutine set_screening_parameter (scrparm_)
     implicit none
     real(DP):: scrparm_
     screening_parameter = scrparm_
     write (stdout,'(5x,F12.7)') 'Screening parameter changed: ', &
          & screening_parameter
  end subroutine set_screening_parameter
  !----------------------------------------------------------------------
  function get_screening_parameter ()
     real(DP):: get_screening_parameter
     get_screening_parameter = screening_parameter
     return
  end function get_screening_parameter
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
  function get_exx_fraction ()
     real(DP):: get_exx_fraction
     get_exx_fraction = exx_fraction
     return
  end function get_exx_fraction
  !-----------------------------------------------------------------------
  function get_dft_name ()
     character (len=20) :: get_dft_name
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
  subroutine set_dft_from_indices(iexch_,icorr_,igcx_,igcc_)
     integer :: iexch_, icorr_, igcx_, igcc_
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
     dft = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
           &//gradc (igcc)
     ! WRITE( stdout,'(a)') dft
     call set_auxiliary_flags
     return
  end subroutine set_dft_from_indices
  !---------------------------------------------------------------------
  subroutine dft_name(iexch_, icorr_, igcx_, igcc_, longname_, shortname_)
  !---------------------------------------------------------------------
  ! convert the four indices iexch, icorr, igcx, igcc
  ! into user-readable strings
  !
  implicit none
  integer iexch_, icorr_, igcx_, igcc_
  character (len=6) :: shortname_
  character (len=20):: longname_
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
  else if (iexch_==1.and.icorr_==4.and.igcx_==12.and.igcc_==4) then
     shortname_ = 'HSE'
  else if (iexch_==1.and.icorr_==4.and.igcx_==11.and.igcc_==4) then
     shortname_ = 'WC'
  else if (iexch_==7.and.(icorr_==10.or.icorr_==2).and.igcx_==9.and. &
           igcc_==7) then
     shortname_ = 'B3LYP'
  else if (iexch_==0.and.icorr_==3.and.igcx_==6.and.igcc_==3) then
     shortname_ = 'OLYP'
  else
     shortname_ = ' '
  end if
  write(longname_,'(4a5)') exc(iexch_),corr(icorr_),gradx(igcx_),gradc(igcc_)
  
  return
end subroutine dft_name

subroutine write_dft_name
!-----------------------------------------------------------------------
   WRITE( stdout, '(5X,"Exchange-correlation      = ",A, &
        &  " (",4I1,")")') TRIM( dft ), iexch, icorr, igcx, igcc
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
        ex = 0.75_DP * ex 
        vx = 0.75_DP * vx 
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
  elseif (icorr ==10) then ! b3lyp
     call vwn (rs, ec, vc)
  elseif (icorr ==11) then
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
        ex   = 0.75_DP * ex
        vxup = 0.75_DP * vxup 
        vxdw = 0.75_DP * vxdw 
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
     call errore ('lsda_functional', 'not implemented', icorr)
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
        evx = 0.75_DP * evx
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
     call errore ('lsda_functional', 'not implemented', icorr)
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
  elseif (igcx == 8) then ! 'pbe0'
     call pbex (rho, grho, 1, sx, v1x, v2x)
     if (exx_started) then
        sx  = 0.75_DP * sx
        v1x = 0.75_DP * v1x
        v2x = 0.75_DP * v2x
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
  elseif (igcc == 7) then !'B3LYP'
     call glyp (rho, grho, sc, v1c, v2c)
     if (exx_started) then
        sc  = 0.81_DP * sc
        v1c = 0.81_DP * v1c
        v2c = 0.81_DP * v2c
     end if
  elseif (igcc == 8) then ! 'PBEsol'
     call pbec (rho, grho, 2, sc, v1c, v2c)
  else
     ! note that if igcc == 5 the hcth functional is called above
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
          igcx == 10 .or. igcx == 12) then
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8 PBE0, igcx=10: PBEsol
     ! igcx=12: HSE
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
       sx = 0.75_DP * sx
       v1xup = 0.75_DP * v1xup
       v1xdw = 0.75_DP * v1xdw
       v2xup = 0.75_DP * v2xup
       v2xdw = 0.75_DP * v2xdw
     end if
     if (igcx == 12 .and. exx_started ) then

        call pbexsr_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
                         v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
                         screening_parameter)
!        write(*,*) sxsr,v1xsr,v2xsr
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
  case(3,4,8,10)
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
        sx = 0.75_DP * sx
        v1xup = 0.75_DP * v1xup
        v1xdw = 0.75_DP * v1xdw
        v2xup = 0.75_DP * v2xup
        v2xdw = 0.75_DP * v2xdw
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

  case default
     call errore ('gcx_spin', 'not implemented', igcx)
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
     call errore ('lsda_functionals', 'not implemented', igcc)
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
      IF( igcc == 3 ) THEN
        RHO=RHOA+RHOB
        IF(RHO.GT.SMALL)  CALL LSD_GLYP(RHOA,RHOB,GRHOAA,GRHOAB,GRHOBB,SC,&
                  V1CA,V2CA,V1CB,V2CB,V2CAB)
      ELSE
        CALL errore( " gcc_spin_more ", " gradiet correction not implemented ", 1 )
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE gcc_spin_more
!
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
             vcdwm, rho, vxupp, vxdwp, vcupp, vcdwp
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
           call xc_spin (rho, zeta - dz, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
           call xc_spin (rho, zeta + dz, ex, ec, vxupp, vxdwp, vcupp, vcdwp)
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
                    dby_my, dby_mz, dbz_rho, dbz_mx, dbz_my, dbz_mz
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
           CALL xc_spin (rho, zeta - dz, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
           CALL xc_spin (rho, zeta + dz, ex, ec, vxupp, vxdwp, vcupp, vcdwp)

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
           dzeta = min (1.d-4, 1.d-2 * abs (zeta) )

           if (dzeta.lt.1.d-7) dzeta = 1.d-7
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
