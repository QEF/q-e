!
! Copyright (C) 2001-2006 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine phq_recover
  !-----------------------------------------------------------------------
  !
  !    This subroutine tests if a restart file exists for the phonon run,
  !    reads data in the header of the file, writes the appropriate message
  !
  !    The restart file is unit "iunrec" (unformatted). Contents:
  ! irr
  !    integer, state of the calculation
  !    irr > 0 irrep up to irr done
  !    irr =-10 to -19 Raman
  !    irr =-20 Electric Field
  ! if (irr>0 .and. okvan) int1, int2
  !    arrays used in phonon calculations with US potentials
  !    (calculated in dvanqq, used in many places)
  ! dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, zstarue0
  !    arrays containing partial results: dyn
  !    or calculated in dynmat0 (not called after restart): dyn00
  !    or calculated in dielec or in zstar_eu: epsilon, zstar*
  !    (not called after a restart if irr>0)
  ! if (irr>0) done_irr, comp_irr, ifat
  !    info on calculated irreps - overrides initialization in phq_setup
  !
  !    phq_readin reads up to here. The following data are read by
  !    routines solve_e, solve_e2, solve_linter:
  ! iter, convt, dr2
  !    info on status of linear-response calculation for a given irrep
  ! dvscfin
  !    self-consistent potential for current iteration and irrep
  ! if (okvan) int3
  !    arrays used with US potentials, calculated by newdq (depends on
  !    self-consistency) and used in adddvscf, called by solve_* 
  !
  !    If a valid restart file is found:
  !    - dynmat0 is not called in any case
  !    - if irr = -20 the electric field calculation (solve_e) restarts from
  !                   the saved value of the potential
  !    - if -10 < irr < -20 solve_e does nothing, the Raman calculation
  !                   (solve_e2), restarts from the saved value of the pot.
  !    - if irr > 0   the entire electric field and Raman section is not
  !                   called, the phonon calculation restarts from irrep irr
  !                   and from the saved value of the potential
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat
  USE io_global,     ONLY : stdout
  USE uspp,          ONLY : okvan
  USE ramanm,        ONLY : lraman, elop, ramtns, eloptns
  USE phcom
  !
  implicit none
  !
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
     read (iunrec) irr0
     !
     ! partially calculated results
     !
     if (okvan.and.irr0>0) read (iunrec) int1, int2
     read (iunrec) dyn, dyn00
     read (iunrec) epsilon, zstareu, zstarue, zstareu0, zstarue0
     IF (irr0>0 .and. lraman) read (iunrec) ramtns
     IF (irr0>0 .and. elop)   read (iunrec) eloptns
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
