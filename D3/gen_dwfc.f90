!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine gen_dwfc (isw_sl)
  !-----------------------------------------------------------------------
  !
  !    Calculates and writes  | d/du(0) psi(k+q) >
  !
  !    Several cases are possible:
  !  isw_sl = 1   : it calculates | d/du(q) psi_k >
  !  isw_sl = 2   : it calculates | d/du(0) psi_k+q >
  !  isw_sl = 3,4 : it calculates | d/du(0) psi_k >
  !
#include "machine.h"
  use pwcom
  use phcom
  use d3com
#ifdef __PARA
  use para
#endif
  implicit none
  integer isw_sl, nirr_x, irr, irr1, imode0
  ! switch
  ! the number of irreducible representation
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  integer, pointer ::  npert_x (:)
  ! the number of perturbations per IR

  if (isw_sl.eq.1) then
     nirr_x = nirr
     npert_x => npert
  else
     nirr_x = nirrg0
     npert_x => npertg0
  endif
  !
  !    For each irreducible representation we compute the change
  !    of the wavefunctions
  !
  do irr = 1, nirr_x
     imode0 = 0
     do irr1 = 1, irr - 1
        imode0 = imode0 + npert_x (irr1)
     enddo
     if (npert_x (irr) .eq.1) then
        write (6, '(//,5x,"Representation #", i3, &
             &                        " mode # ",i3)') irr, imode0 + 1
     else
        write (6, '(//,5x,"Representation #", i3, &
             &                 " modes # ",3i3)') irr,  (imode0 + irr1, irr1 = &
             & 1, npert_x (irr) )

     endif
     call solve_linter_d3 (irr, imode0, npert_x (irr), isw_sl)
  enddo
  !
  ! Writes FermiEnergy shift on a file
  !
#ifdef __PARA
  if (me.ne.1.or.mypool.ne.1) goto 210
#endif
  if (isw_sl.eq.3.and.degauss.ne.0.d0) then
     rewind (unit = iuef)
     write (iuef) ef_sh
  endif
#ifdef __PARA
210 continue
#endif
  !
  ! closes and opens some units --useful in case of interrupted run--
  !

  call close_open (isw_sl)
  return
end subroutine gen_dwfc
