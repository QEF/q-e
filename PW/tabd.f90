!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine tabd (nt, occd)
!-----------------------------------------------------------------------
!
! This routine is a table (far from being complete) for the total number
! of d electrons in transition metals (PPs usually are built on non
! physical configurations)
!
use pwcom
implicit none
real(kind=DP) :: occd
                         ! output: the total number of d electrons

integer :: nt
if (psd (nt) .eq.'Mn') then
   occd = 5.d0
elseif (psd (nt) .eq.'Fe') then
   occd = 6.d0
elseif (psd (nt) .eq.'Co') then
   occd = 9.d0
elseif (psd (nt) .eq.'Ni') then
   occd = 8.d0
elseif (psd (nt) .eq.'Cu') then
   occd = 10.d0
else
   occd = 0.d0
   call error ('tabd', 'pseudopotential not yet inserted', - 1)

endif
return
end subroutine tabd

