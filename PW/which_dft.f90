!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine which_dft (dft, iexch, icorr, igcx, igcc)
  !-----------------------------------------------------------------------
  !
  use parameters
  use parser
  implicit none
  ! input
  character (len=*) :: dft
  ! output
  integer :: iexch, icorr, igcx, igcc
  ! data
  integer :: nxc, ncc, ngcx, ngcc
  parameter (nxc = 1, ncc = 9, ngcx = 3, ngcc = 4)
  character (len=3) :: exc, corr
  character (len=4) :: gradx, gradc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)
  ! local
  integer :: len, l, i, notset
  character (len=50):: dftout * 50
  data notset / - 1 /
  data exc / 'NOX', 'SLA' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
       'OBW', 'GL' /
  data gradx / 'NOGX', 'B88', 'GGX', 'PBE' /


  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBE' /
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


  if (iexch.eq.notset) call set_dft_value (iexch, 1)
  ! Default value: Perdew-Zunger correlation

  if (icorr.eq.notset) call set_dft_value (icorr, 1)
  ! Default value: no gradient correction on exchange

  if (igcx.eq.notset) call set_dft_value (igcx, 0)
  ! Default value: no gradient correction on correlation

  if (igcc.eq.notset) call set_dft_value (igcc, 0)
  !

  dftout = exc (iexch) //'-'//corr (icorr) //'-'//gradx (igcx) //'-' &
       &//gradc (igcc)
  !cc      write (6,'(a)') dftout
  return
end subroutine which_dft
!
!-----------------------------------------------------------------------
subroutine set_dft_value (m, i)
  !-----------------------------------------------------------------------
  !
  implicit none
  ! input / output
  integer :: m, i
  ! local
  integer :: notset

  parameter (notset = - 1)

  if (m.ne.notset.and.m.ne.i) call errore ('decifra', 'two conflictin &
       &g matching values', 1)

  m = i
  return

end subroutine set_dft_value
!-----------------------------------------------------------------------
logical function matches (string1, string2)
  !-----------------------------------------------------------------------
  !
  implicit none
  character (len=*) :: string1, string2
  integer :: len1, len2, l


  len1 = len_trim(string1)
  len2 = len_trim(string2)
  do l = 1, len2 - len1 + 1
     if (string1 (1:len1) .eq.string2 (l:l + len1 - 1) ) then
        matches = .true.
        return
     endif

  enddo

  matches = .false.
  return
end function matches
!
!-----------------------------------------------------------------------
function capital (character)
  !-----------------------------------------------------------------------
  !
  !   converts character to capital if lowercase
  !   copy character to output in all other cases
  !
  implicit none
  character (len=1) :: capital, character
  !
  character(len=26) :: minuscole='abcdefghijklmnopqrstuvwxyz', &
                       maiuscole='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer :: i
  !
  do i=1,26
     if (character.eq.minuscole(i:i)) then
        capital=maiuscole(i:i)
        return
     end if
  end do
  capital = character
  !
  return
end function capital
