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
#include "machine.h"
  !
  USE ions_base,     ONLY : nat
  USE io_global,     ONLY : stdout
  use pwcom
  USE kinds,         ONLY : DP
  use phcom
  USE control_flags, ONLY : iswitch
  
  implicit none

  integer :: ifat0 (nat), comp_irr0 (3 * nat), irr, na
  ! dummy variable to read
  ! dummy variable to read
  ! counter on representations
  ! counter on atoms

  iunrec = 99
  call seqopn (iunrec, 'recover', 'unformatted', recover)
  irr0 = 0
  iter0 = 0
  zstarue0 (:,:) = (0.d0, 0.d0)
  if (recover) then
     if (okvan) read (iunrec) int1, int2

     read (iunrec) dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, &
          zstarue0

     read (iunrec) irr0, iter0, convt, done_irr, comp_irr0, ifat0
     if (iter0 == 0 .and. iswitch == - 4) then
        call errore ('recover', 'recover not possible', 1)
        dyn = dyn00
        do irr = 1, nirr
           done_irr (irr) = 0
        enddo
     endif
     if (iter0 /= 0) then
        !
        !    If iter0.ne.0 we must recover also the irreducible representations
        !    to do and atoms
        !
        do irr = 1, nirr
           comp_irr (irr) = comp_irr0 (irr)
        enddo
        nat_todo = 0
        do na = 1, nat
           ifat (na) = ifat0 (na)
           if (ifat (na) == 1) then
              nat_todo = nat_todo + 1
              atomo (nat_todo) = na
           endif
        enddo
        all_comp = ( nat_todo == nat )
     endif

     if (irr0 == - 1) then
        WRITE( stdout, '(/,4x," Reading only int1 and int2" )')
        close (unit = iunrec, status = 'keep')
     elseif (irr0 == - 2) then
        WRITE( stdout, '(/,4x," Restart from Iteration #",i5, &
             &                   " of Elect. Field")') iter0 + 1
     elseif (irr0 > 0) then
        if (iter0 /= 0) then
           WRITE( stdout, '(/,4x," Restart from Iteration #",i5, &
                &             " of Representation #",i5)') iter0 + 1, irr0
        else
           if (irr0 == nirr) then
              WRITE( stdout, '(/,4x," From recover: Dynamical ", &
                   &       " Matrix calculation ")')
           else
              do irr = 1, nirr
                 if ( (comp_irr (irr) == 1) .and. (done_irr (irr) == 0) ) then
                    WRITE( stdout, '(/,4x," Restart from first iteration", &
                         &      " of Representation # ",i5)') irr
                    goto 1000
                 endif
              enddo
1000          continue
           endif
           close (unit = iunrec, status = 'keep')
        endif
     else
        call errore ('phq_recover', 'wrong irr0', 1)
     endif
  else
     close (unit = iunrec, status = 'delete')
  endif

  return
end subroutine phq_recover
