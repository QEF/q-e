!
! Copyright (C) 2004 PWSCF group
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
!  retrive functions:  get_dft_name
!                      get_iexch
!                      get_icorr
!                      get_igcx
!                      get_igcc
!                      dft_name
!                      write_dft_name
!  retrivable logical: ishybrid
!                      ismeta
!                      isgradient
!
!  XC computation drivers: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
!
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  PRIVATE
  SAVE
  ! subroutines/functions managing dft name and indices
  PUBLIC  :: set_dft_from_indices, set_dft_from_name
  PUBLIC  :: enforce_input_dft, write_dft_name, dft_name
  PUBLIC  :: get_dft_name, get_iexch, get_icorr, get_igcx, get_igcc
  ! two (still) public variables
  PUBLIC  :: ismeta, ishybrid, isgradient
  ! driver subroutines computing XC
  PUBLIC  :: xc, xc_spin, gcxc, gcx_spin, gcc_spin, gcc_spin_more
  !
  ! PRIVATE variables defining the DFT functional
  !
  PRIVATE :: dft, dft_shortname, iexch, icorr, igcx, igcc
  PRIVATE :: discard_input_dft
  !
  character (len=20) :: dft = 'not set'
  character (len=4)  :: dft_shortname = ' '
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
  !
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbx"    Perdew-Burke-Ernzenhof exch    igcx =3
  !              "rpb"    revised PBE by Zhang-Yang      igcx =4
  !              "hcth"   Cambridge exch, Handy et al    igcx =5
  !              "optx"   Handy's exchange functional    igcx =6
  !
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Lee-Yang-Parr                  igcc =3
  !              "pbx"    Perdew-Burke-Ernzenhof corr    igcc =4
  !              "hcth"   Cambridge corr, Handy et al    igcc =5
  !
  ! Special cases (dft_shortnames):
  !              "bp"   = "b88+p86"         = Becke-Perdew grad.corr.
  !              "pw91" = "pw +ggx+ggc"     = PW91 (aka GGA)
  !              "blyp" = "sla+b88+lyp+blyp"= BLYP
  !              "pbe"  = "sla+pw+pbx+pbc"  = PBE
  !              "revpbe"="sla+pw+rpb+pbc"  = revPBE (Zhang-Yang)
  !              "hcth" = "nox+noc+hcth+hcth"=HCTH/120
  !              "olyp" = "nox+lyp+optx+blyp" !!! UNTESTED !!!
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
  !              oep

  integer, parameter:: notset = -1
  !
  integer :: iexch = notset
  integer :: icorr = notset
  integer :: igcx  = notset
  integer :: igcc  = notset
  logical :: ismeta = .false.
  logical :: ishybrid = .false.
  logical :: isgradient = .false.

  logical :: discard_input_dft = .false.
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlations
  !
  !    ismeta: .TRUE. if gradient correction is of meta-gga type
  !    ishybrid: .TRUE. if the xc finctional is an HF+DFT hybrid like
  !              PBE0 or B3LYP of HF itself
  !
  ! see comments above and routine "set_dft_from_name" below 
  !
  !
  ! data
  integer :: nxc, ncc, ngcx, ngcc
  parameter (nxc = 5, ncc = 9, ngcx = 7, ngcc = 6)
  character (len=4) :: exc, corr
  character (len=4) :: gradx, gradc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)

  data exc / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' /
  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX', 'META' /
  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'META'/

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
    if (matches ('BLYP', dftout) ) call set_dft_value (igcx, 1)

    ! special case : revPBE
    if (matches ('REVPBE', dftout) ) then
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 4)
       call set_dft_value (igcc, 4)
    else if (matches('RPBE',dftout)) then
         call errore('set_dft_from_name', &
     &   'RPBE (Hammer-Hansen-Norskov) not implemented (revPBE is)',1)
   else if (matches ('PBE', dftout) ) then
    ! special case : PBE
       call set_dft_value (icorr,4)
       call set_dft_value (igcx, 3)
       call set_dft_value (igcc, 4)
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
    ! ... special case : OEP is exact exchange no GC part
    !
    IF ( matches( 'OEP', dftout ) ) THEN
       !
       CALL set_dft_value( igcx,  0 )
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

    dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
         &//gradc (igcc)
    ! WRITE( stdout,'(a)') dftout

    dft = dft_

    isgradient =  (igcx > 0) .or. (igcc > 0)
    ismeta     =  (igcx == 7) .or. (igcx == 6 )

    return
  end subroutine set_dft_from_name
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

     write (stdout,'(a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write (stdout,'(a)') " XC functional is enforced to be :"
     write (stdout,'(a)') dft_
     write (stdout,'(a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

     call set_dft_from_name (dft_)
     if (dft == 'not set') call errore('enforce_input_dft','cannot fix unset dft',1)
     discard_input_dft = .true.
     return
  end subroutine enforce_input_dft

  function get_iexch ()
     integer get_iexch
     get_iexch = iexch
     return
  end function get_iexch
  function get_icorr ()
     integer get_icorr
     get_icorr = icorr
     return
  end function get_icorr
  function get_igcx ()
     integer get_igcx
     get_igcx = igcx
     return
  end function get_igcx 
  function get_igcc ()
     integer get_igcc
     get_igcc = igcc
     return
  end function get_igcc 

  function get_dft_name ()
     character (len=20) :: get_dft_name
     get_dft_name = dft
     return
  end function get_dft_name

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
     isgradient =  (igcx > 0) .or. (igcc > 0)
     ismeta     =  (igcx == 7) .or. (igcx == 6 )
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
  character (len=4) :: shortname_
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
  USE kinds
  implicit none

  real(DP) :: rho, ec, vc, ex, vx
  !
  real(DP), parameter :: small = 1.d-10,  third = 1.d0 / 3.d0, &
       pi34 = 0.6203504908994d0  ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0d0
     vc = 0.0d0
     ex = 0.0d0
     vx = 0.0d0
     return
  else
     rs = pi34 / rho**third
     ! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  endif
  !..exchange
  if (iexch == 1) then
     call slater (rs, ex, vx)
  ELSEIF (iexch == 2) THEN
     call slater1(rs, ex, vx)
  ELSEIF (iexch == 3) THEN
     CALL slater_rxc(rs, ex, vx)
  else
     ex = 0.0d0
     vx = 0.0d0
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
  else
     ec = 0.0d0
     vc = 0.0d0
  endif
  !
  return
end subroutine xc
!
!-----------------------------------------------------------------------
subroutine gcxc (rho, grho, sx, sc, v1x, v2x, v1c, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     exchange  :  Becke88
  !                  GGA (Generalized Gradient Approximation), PW91
  !                  PBE
  !                  revPBE
  !     correlation: Perdew86
  !                  GGA (PW91)
  !                  Lee-Yang-Parr
  !                  PBE
  !
  !     input:  rho, grho=|\nabla rho|^2
  !     definition:  E_x = \int E_x(rho,grho) dr
  !     output: sx = E_x(rho,grho)
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             sc, v1c, v2c as above for correlation
  !
  USE kinds
  implicit none

  real(DP) :: rho, grho, sx, sc, v1x, v2x, v1c, v2c
  real(DP), parameter:: small = 1.d-10

  ! exchange
  if (rho <= small) then
     sx = 0.0d0
     v1x = 0.0d0
     v2x = 0.0d0
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
  else
     sx = 0.0d0
     v1x = 0.0d0
     v2x = 0.0d0
  endif
  ! correlation
  if (rho.le.small) then
     sc = 0.0d0
     v1c = 0.0d0
     v2c = 0.0d0
  elseif (igcc == 1) then
     call perdew86 (rho, grho, sc, v1c, v2c)
  elseif (igcc == 2) then
     call ggac (rho, grho, sc, v1c, v2c)
  elseif (igcc == 3) then
     call glyp (rho, grho, sc, v1c, v2c)
  elseif (igcc == 4) then
     call pbec (rho, grho, sc, v1c, v2c)
  else
     ! note that if igcc == 5 the hcth functional is called above
     sc = 0.0d0
     v1c = 0.0d0
     v2c = 0.0d0
  endif
  !
  return
end subroutine gcxc
!
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
  USE kinds
  implicit none

  real(DP) :: rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw
  !
  real(DP), parameter :: small= 1.d-10, third = 1.d0/3.d0, &
       pi34= 0.6203504908994d0 ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0d0
     vcup = 0.0d0
     vcdw = 0.0d0
     ex = 0.0d0
     vxup = 0.0d0
     vxdw = 0.0d0
     return
  else
     rs = pi34 / rho**third
  endif
  !..exchange
  if (iexch == 1) then
     call slater_spin (rho, zeta, ex, vxup, vxdw)
  elseif (iexch == 2) then
     call slater1_spin (rho, zeta, ex, vxup, vxdw)
  ELSEIF (iexch == 3) THEN
     call slater_rxc_spin ( rho, zeta, ex, vxup, vxdw )
  else
     ex = 0.0d0
     vxup = 0.0d0
     vxdw = 0.0d0
  endif
  !..correlation
  if (icorr == 0) then
     ec = 0.0d0
     vcup = 0.0d0
     vcdw = 0.0d0
  elseif (icorr == 1) then
     call pz_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 3) then
     call lsd_lyp (rho, zeta, ec, vcup, vcdw)  !  from CP/FPMD (more_functionals)
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
subroutine gcx_spin (rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, &
     v1xdw, v2xup, v2xdw)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !     Implemented:  Becke88, GGA (PW91), PBE, revPBE
  !
  USE kinds
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
  real(DP), parameter :: small = 1.d-10
  real(DP) :: rho, sxup, sxdw
  integer :: iflag
  !
  !
  ! exchange
  rho = rhoup + rhodw
  if (rho <= small .or. igcx == 0) then
     sx = 0.0d0
     v1xup = 0.0d0
     v2xup = 0.0d0
     v1xdw = 0.0d0
     v2xdw = 0.0d0
  elseif (igcx == 1) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.d0
        v1xup = 0.d0
        v2xup = 0.d0
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.d0
        v1xdw = 0.d0
        v2xdw = 0.d0
     endif
     sx = sxup + sxdw
  elseif (igcx == 2) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call ggax (2.d0 * rhoup, 4.d0 * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.d0
        v1xup = 0.d0
        v2xup = 0.d0
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call ggax (2.d0 * rhodw, 4.d0 * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.d0
        v1xdw = 0.d0
        v2xdw = 0.d0
     endif
     sx = 0.5d0 * (sxup + sxdw)
     v2xup = 2.d0 * v2xup
     v2xdw = 2.d0 * v2xdw
  elseif (igcx == 3 .or. igcx == 4) then
     ! igcx=3: PBE  igcx=4: revised PBE
     if (igcx == 3) then
        iflag = 1
     else
        iflag = 2
     endif
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pbex (2.d0 * rhoup, 4.d0 * grhoup2, iflag, sxup, v1xup, v2xup)
     else
        sxup = 0.d0
        v1xup = 0.d0
        v2xup = 0.d0
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pbex (2.d0 * rhodw, 4.d0 * grhodw2, iflag, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.d0
        v1xdw = 0.d0
        v2xdw = 0.d0
     endif
     sx = 0.5d0 * (sxup + sxdw)
     v2xup = 2.d0 * v2xup
     v2xdw = 2.d0 * v2xdw
  else
     call errore ('gcx_spin', 'not implemented', igcx)
  endif
  !
  return
end subroutine gcx_spin
!
!-----------------------------------------------------------------------
subroutine gcc_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for correlations - Hartree a.u.
  !     Implemented:  Perdew86, GGA (PW91), PBE
  !
  USE kinds
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

  real(DP), parameter :: small = 1.d-10, epsr=1.d-6
  !
  !
  if ( abs(zeta) > 1.d0 ) then
     sc = 0.0d0
     v1cup = 0.0d0
     v1cdw = 0.0d0
     v2c = 0.0d0
     return
  else
     !
     ! ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
     zeta = SIGN( MIN( ABS( zeta ), ( 1.D0 - epsr ) ) , zeta )
  endif

  if (rho <= small .or. sqrt(abs(grho)) <= small) then
     sc = 0.0d0
     v1cup = 0.0d0
     v1cdw = 0.0d0
     v2c = 0.0d0
  elseif (igcc == 1) then
     call perdew86_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 2) then
     call ggac_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 3 .or. igcc > 4) then
     call errore ('lsda_functionals', 'not implemented', igcc)
  elseif (igcc == 4) then
        call pbec_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  else
     sc = 0.0d0
     v1cup = 0.0d0
     v1cdw = 0.0d0
     v2c = 0.0d0
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

      USE kinds, ONLY: DP
      IMPLICIT NONE
      REAL(DP) :: RHOA,RHOB,GRHOAA,GRHOBB,GRHOAB
      REAL(DP) :: SC,V1CA,V2CA,V1CB,V2CB,V2CAB

      ! ... Gradient Correction for correlation

      REAL(DP) :: SMALL, RHO
      PARAMETER(SMALL=1.D-20)

      SC=0.0D0
      V1CA=0.0D0
      V2CA=0.0D0
      V1CB=0.0D0
      V2CB=0.0D0
      V2CAB=0.0D0
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


end module funct
