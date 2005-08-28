!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function spinor(l,j,m,spin)
! This function calculates the numerical coefficient of a spinor
! with orbital angular momentum l, total angular momentum j, 
! projection along z of the total angular momentum m+-1/2. Spin selects
! the up (spin=1) or down (spin=2) coefficient.

use kinds
implicit none

real(DP) :: spinor    
integer :: l, &            ! orbital angular momentum
           m, &            ! projection of the total angular momentum+-1/2
           spin            ! 1 or 2 select the component

real(DP) :: j         ! total angular momentum
real(DP) :: denom     ! denominator

if (spin.ne.1.and.spin.ne.2) call errore('spinor','spin direction unknown',1)
if (m.lt.-l-1.or.m.gt.l) call errore('spinor','m not allowed',1)

denom=1.d0/(2.d0*l+1.d0)
if (abs(j-l-0.5d0).lt.1.d-8) then
   if (spin.eq.1) spinor= sqrt((l+m+1.d0)*denom)
   if (spin.eq.2) spinor= sqrt((l-m)*denom)
elseif (abs(j-l+0.5d0).lt.1.d-8) then
   if (m.lt.-l+1) then
      spinor=0.d0
   else
      if (spin.eq.1) spinor= sqrt((l-m+1.d0)*denom)
      if (spin.eq.2) spinor= -sqrt((l+m)*denom)
   endif
else
   call errore('spinor','j and l not compatible',1)
endif

return
end function spinor
