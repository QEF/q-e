!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function scnds ()
  !
  ! Returns elapsed CPU time in seconds since its first call
  ! uses standard f90 call (with several exceptions)
  !
  use parameters
  implicit none
  real(kind=DP) :: scnds
  !
#if defined(PGI) || defined(DEC)
  real(kind=4) :: etime, tarry(2)
  ! etime:  system function, returns the CPU time in sec.
  ! PGI compiler has no intrinsic f90 cpu_time
#endif
  real(kind=DP) :: t0, t1
  ! t0 contains the time of the first call
  ! t1 contains the present time
  logical :: first=.true.
  save first, t0
  !
#if defined(PGI) || defined (DEC)
  t1 = etime( tarry )
#elif defined(HITACHI)
  call CLOCK(IT,2)
  t1 = it
# else
  call cpu_time(t1)
#endif
  !
  if (first) then
     t0 = t1
     scnds = 0.d0
     first = .false.
  else
     scnds = t1 - t0
  endif
  return
end function scnds

