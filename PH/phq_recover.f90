!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine phq_recover
  !-----------------------------------------------------------------------
  !
  !    This subroutine tests if a recover file exists for the phonon run
  !    and writes the appropriate messages if the run starts from
  !    a given iteration of a given irreducible representation
  !
#include "f_defs.h"
  !
  USE ions_base,     ONLY : nat
  USE io_global,     ONLY : stdout
  use pwcom
  USE kinds,         ONLY : DP
  use phcom
  USE control_flags, ONLY : modenum
  
  implicit none

  integer :: irr, na
  ! counter on representations
  ! counter on atoms
  logical :: exst

  iunrec = 99
  call seqopn (iunrec, 'recover', 'unformatted', exst)
  irr0 = 0
  zstarue0 (:,:) = (0.d0, 0.d0)
  recover = recover .AND. exst
  if (recover) then
     !
     ! irr: state of the calculation
     ! irr > 0 irrep up to irr done
     ! irr = 0 nothing done
     ! irr =-10 to -19 Raman
     ! irr =-20 Electric Field
     !
     read (iunrec) irr0
     !
     ! partially calculated results
     !
     if (okvan.and.irr0>0) read (iunrec) int1, int2
     read (iunrec) dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, &
          zstarue0
     !
     if (irr0 > 0) then
        read (iunrec) done_irr, comp_irr, ifat
        nat_todo = 0
        do na = 1, nat
           if (ifat (na) == 1) then
              nat_todo = nat_todo + 1
              atomo (nat_todo) = na
           endif
        enddo
        all_comp = ( nat_todo == nat )
     end if

     if (irr0 == -20) then
        WRITE( stdout, '(/,4x," Restart in Electric Field calculation")')
     elseif (irr0 > -20 .AND. irr0 <= -10) then
        WRITE( stdout, '(/,4x," Restart in Raman calculation")') 
     elseif (irr0 > 0 .AND. irr0 <= nirr) then
        WRITE( stdout, '(/,4x," Restart in Phonon calculation")')
     else
        call errore ('phq_recover', 'wrong restart file', 1)
     endif
  else
     close (unit = iunrec, status = 'delete')
  endif

  return
end subroutine phq_recover
