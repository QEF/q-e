!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_fft_dim
  !-----------------------------------------------------------------------
  !     This routine computes the dimensions of the minimum FFT grid
  !     compatible with the input cut-off
  !
  !     NB: The values of nr1, nr2, nr3 are computed only if they are not
  !     given as input parameters. Input values are kept otherwise.
  !
  use pwcom
  use fft_scalar, only: allowed
  implicit none

  integer, parameter :: nmax = 5000
  ! an unreasonably big number for a FFT grid
  !
  ! the values of nr1, nr2, nr3 are computed only if they are not given
  ! as input parameters
  !
  if (nr1 == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     nr1 = int (2 * sqrt (gcutm) * sqrt (at (1, 1) **2 + at (2, 1) &
          **2 + at (3, 1) **2) ) + 1
10   continue
     if (nr1 > nmax) &
          call errore ('set_fft_dim', 'nr1 is unreasonably large', nr1)
     if (allowed (nr1) ) goto 15
     nr1 = nr1 + 1
     goto 10
  else
     if (.not.allowed (nr1) ) call errore ('set_fft_dim', &
          'input nr1 value not allowed', 1)
  endif
15 continue
  !
  !    here we compute nr1s if it is not in input
  !
  if (nr1s == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     if (doublegrid) then
        nr1s = 2 * int (sqrt (gcutms) * sqrt (at (1, 1) **2 + at (2, 1)**2 &
             + at (3, 1) **2) ) + 1
     else
        nr1s = nr1
     endif
11   continue
     if (nr1s > nmax) &
          call errore ('set_fft_dim', 'nr1s is unreasonably large', nr1s)
     if (allowed (nr1s) ) goto 16
     nr1s = nr1s + 1
     goto 11
  endif
16 continue
  !
  if (nr2 == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     nr2 = int (2 * sqrt (gcutm) * sqrt (at (1, 2) **2 + at (2, 2) &
          **2 + at (3, 2) **2) ) + 1
20   continue
     if (nr2 > nmax) &
          call errore ('set_fft_dim', 'nr2 is unreasonably large', nr2)
     if (allowed (nr2) ) goto 25
     nr2 = nr2 + 1
     goto 20
  else
     if (.not.allowed (nr2) ) call errore ('set_fft_dim', &
          'input nr2 value not allowed', 2)
  endif
25 continue
  !
  !    here we compute nr2s if it is not in input
  !
  if (nr2s == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     if (doublegrid) then
        nr2s = 2 * int (sqrt (gcutms) * sqrt (at (1, 2) **2 + at (2, 2) **2 &
             + at (3, 2) **2) ) + 1
     else
        nr2s = nr2
     endif
21   continue
     if (nr2s > nmax) &
          call errore ('set_fft_dim', 'nr2s is unreasonably large', nr2s)
     if (allowed (nr2s) ) goto 26
     nr2s = nr2s + 1
     goto 21
  endif
26 continue
  !
  if (nr3 == 0) then
     !
     ! estimate nr3 and check if it is an allowed value for FFT
     !
     nr3 = int (2 * sqrt (gcutm) * sqrt (at (1, 3) **2 + at (2, 3) &
          **2 + at (3, 3) **2) ) + 1
30   continue
     if (nr3 > nmax) &
          call errore ('set_fft_dim', 'nr3 is unreasonably large', nr3)
     if (allowed (nr3) ) goto 35
     nr3 = nr3 + 1
     goto 30
  else
     if (.not.allowed (nr3) ) call errore ('set_fft_dim', &
          'input nr3 value not allowed', 3)
  endif
35 continue
  !
  !    here we compute nr3s if it is not in input
  !
  if (nr3s == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     if (doublegrid) then
        nr3s = 2 * int (sqrt (gcutms) * sqrt (at (1, 3) **2 + at (2, &
             3) **2 + at (3, 3) **2) ) + 1
     else
        nr3s = nr3
     endif
31   continue
     if (nr3s > nmax) &
          call errore ('set_fft_dim', 'nr3s is unreasonably large', nr3s)
     if (allowed (nr3s) ) goto 36
     nr3s = nr3s + 1
     goto 31
  endif
36 continue
  !
  if (nr1s > nr1 .or. nr2s > nr2 .or. nr3s > nr3) then
     write (6, * ) nr1s, nr2s, nr3s, nr1, nr2, nr3
     call errore ('set_fft_dim', 'smooth grid larger than big grid', 1)
  endif
  !
  if (.not.doublegrid .and. (nr1s.ne.nr1.or.nr2s.ne.nr2.or.nr3s.ne.nr3) ) &
       call errore ('set_fft_dim', 'smooth and big grids must be the same', 2)
  return
end subroutine set_fft_dim

