!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module funct
  IMPLICIT NONE
  PRIVATE
  SAVE
  PUBLIC :: dft, iexch, icorr, igcx, igcc, which_dft
  !
  character (len=20) :: dft
  !
  ! dft is the exchange-correlation functional, described by
  ! any nonconflicting combination of the following keywords
  ! (case-insensitive):
  !
  ! Exchange:    "nox"    none                           iexch=0
  !              "sla"    Slater (alpha=2/3)             iexch=1 (default)
  !              "rxc"    Relativistic Slater (?)        iexch=2
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
  ! Gradient Correction on Exchange:
  !              "nogx"   none                           igcx =0 (default)
  !              "b88"    Becke88 (beta=0.0042)          igcx =1
  !              "ggx"    Perdew-Wang 91                 igcx =2
  !              "pbe"    Perdew-Burke-Ernzenhof         igcx =3
  ! Gradient Correction on Correlation:
  !              "nogc"   none                           igcc =0 (default)
  !              "p86"    Perdew86                       igcc =1
  !              "ggc"    Perdew-Wang 91 corr.           igcc =2
  !              "blyp"   Becke88 + Lee-Yang-Parr        igcc =3
  !              "pbe"    Perdew-Burke-Ernzenhof corr    igcc =4
  !
  ! Special cases: "bp" = "b88+p86"        = Becke-Perdew grad.corr
  !              "pw91" = "pw +ggx+ggc"    = PW91 (aka GGA)
  !              "pbe"  = "sla+pw+ggx+ggc" = PBE
  !
  integer :: iexch, icorr, igcx, igcc
  !
  ! internal indices for exchange-correlation
  !    iexch: type of exchange
  !    icorr: type of correlation
  !    igcx:  type of gradient correction on exchange
  !    igcc:  type of gradient correction on correlations
  ! see comments above and routine "which_dft" below 
  !

CONTAINS
  !-----------------------------------------------------------------------
  subroutine which_dft (dft_)
    !-----------------------------------------------------------------------
    !
    ! translates a string containing the exchange-correlation name
    ! into internal indices iexch, icorr, igcx, igcc
    !
    USE parser, ONLY: matches, capital
    !
    implicit none
    ! input
    character (len=*) :: dft_
    ! data
    integer :: nxc, ncc, ngcx, ngcc
    parameter (nxc = 2, ncc = 9, ngcx = 3, ngcc = 4)
    character (len=3) :: exc, corr
    character (len=4) :: gradx, gradc
    dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)
    ! local
    integer :: len, l, i
    integer, parameter:: notset = -1
    character (len=50):: dftout * 50
    data exc / 'NOX', 'SLA', 'RXC' /
    data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
         'OBW', 'GL' /
    data gradx / 'NOGX', 'B88', 'GGX', 'PBE' /
    data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBE' /
    !
    ! convert to uppercase
    len = len_trim(dft)
    dftout = ' '
    do l = 1, len
       dftout (l:l) = capital (dft (l:l) )
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

    ! special case : PBE
    if (matches ('PBE', dftout) ) then
       call set_dft_value (iexch, 1)
       call set_dft_value (icorr, 4)
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
    !      WRITE( stdout,'(a)') dftout
    return
  end subroutine which_dft
  !
  !-----------------------------------------------------------------------
  subroutine set_dft_value (m, i)
    !-----------------------------------------------------------------------
    !
    implicit none
    integer :: m, i
    ! local
    integer, parameter :: notset = - 1

    if ( m /= notset .and. m /= i) &
         call errore ('set_dft_value', 'two conflicting matching values', 1)
    m = i
    return

  end subroutine set_dft_value
end module funct
