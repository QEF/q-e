!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------


function allowed (nr)
  !-----------------------------------------------------------------------
#include "machine.h"
  ! find if the fft dimension is a good one
  ! a "bad one" is either not implemented (as on IBM with ESSL)
  ! or implemented but with awful performances (most other cases)
  use parameters
  implicit none
  integer :: nr

  logical :: allowed
  integer :: pwr (5)
  integer :: mr, i, fac, p, maxpwr
  integer :: factors (5)


  data factors / 2, 3, 5, 7, 11 /
  ! find the factors of the fft dimension
  mr = nr
  do i = 1, 5
     pwr (i) = 0
  enddo
  do i = 1, 5
     fac = factors (i)
     maxpwr = nint (log (float (mr) ) / log (float (fac) ) ) + 1
     do p = 1, maxpwr
        if (mr.eq.1) goto 10
        if (mod (mr, fac) .eq.0) then
           mr = mr / fac
           pwr (i) = pwr (i) + 1
        endif
     enddo

  enddo
10 if (nr.ne.mr * 2**pwr (1) * 3**pwr (2) * 5**pwr (3) * 7**pwr (4) &
       * 11**pwr (5) ) call errore ('allowed', 'what ?!?', 1)
  if (mr.ne.1) then
     ! fft dimension contains factors > 11 : no good in any case
     allowed = .false.
  else
#ifdef __CERNFFT
     ! this is for the generic (cernlib) case
     allowed = .true.
     ! specific (machine- and library-dependent cases
#else
#ifdef __AIX
     ! IBM machines with essl libraries

     allowed = pwr (1) .ge.1.and.pwr (2) .le.2.and.pwr (3) &
          .le.1.and.pwr (4) .le.1.and.pwr (5) .le.1.and. ( (pwr (2) &
          .eq.0.and.pwr (3) + pwr (4) + pwr (5) .le.2) .or. (pwr (2) &
          .ne.0.and.pwr (3) + pwr (4) + pwr (5) .le.1) )
#else
     ! fftw and all other cases
     allowed = pwr (4) .eq.0.and.pwr (5) .eq.0
#endif
#endif
  endif

  return
end function allowed

