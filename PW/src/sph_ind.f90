!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function sph_ind(l,j,m,spin)
! This function calculates the m index of the spherical harmonic
! in a spinor with orbital angular momentum l, total angular 
! momentum j, projection along z of the total angular momentum m+-1/2. 
! Spin selects the up (spin=1) or down (spin=2) coefficient.
!
use kinds
implicit none

integer :: sph_ind
integer :: l, &            ! orbital angular momentum
           m, &            ! projection of the total angular momentum+-1/2
           spin            ! 1 or 2 select the component

real(DP) :: j         ! total angular momentum

if (spin.ne.1.and.spin.ne.2) call errore('sph_ind','spin direction unknown',1)
if (m.lt.-l-1.or.m.gt.l) call errore('sph_ind','m not allowed',1)

if (abs(j-l-0.5d0).lt.1.d-8) then
   if (spin.eq.1) sph_ind= m
   if (spin.eq.2) sph_ind= m+1
elseif (abs(j-l+0.5d0).lt.1.d-8) then
   if (m.lt.-l+1) then
      sph_ind=0
   else
      if (spin.eq.1) sph_ind= m-1
      if (spin.eq.2) sph_ind= m
   endif
else
   write(6,*) l, j
   call errore('sph_ind','l and j not compatible',1)
endif
if (sph_ind.lt.-l.or.sph_ind.gt.l) sph_ind=0

return
end function sph_ind
