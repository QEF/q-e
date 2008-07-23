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
  !    This subroutine tests if a xml restart file and an unformatted
  !    data file exist for the phonon run, reads data in the xml file, 
  !    writes the appropriate message
  !
  !    The unformatted file is unit "iunrec". The xml file is in the
  !    directory prefix.phsave. The xml file contains
  ! where_rec  a string with information of the point where the calculation
  !            stopped
  ! rec_code
  !    integer, state of the calculation
  !    rec_code > 0 phonon (solve_linter 10 or phqscf 20)
  !    rec_code =-10 to -19 Raman
  !    rec_code =-20 Electric Field
  ! dyn, epsilon, zstareu, zstarue0
  !    arrays containing partial results: dyn
  !    or calculated in dielec or in zstar_eu: epsilon, zstar*
  !    (not called after a restart if irr>0)
  ! if (rec_code>0) done_irr, comp_irr
  !    info on calculated irreps - overrides initialization in phq_setup
  !
  !    phq_recover reads up to here. The following data are in the 
  !    unformatted file and are read by
  !    routines solve_e, solve_e2, solve_linter:
  ! iter, dr2
  !    info on status of linear-response calculation for a given irrep.
  !    IMPORTANT: at convergence, iter is reset to 0 so that if the code
  !    restarts it will redo the calculation of the given irrep from the
  !    beginning. The reason for this apparent absurdity is that after
  !    convergence is achieved, files containing information needed for
  !    restarting may be lost (files opened by mix_pot for instance)
  !    or overwritten at the subsequent interation (files containing
  !    dvpsi and dpsi). While not efficient in some specific case, this
  !    is the only safe way to restart without trouble.
  ! dvscfin
  !    self-consistent potential for current iteration and irrep
  ! if (okvan) int1, int2, int3
  !    arrays used with US potentials : int1 and int2 calculated in dvanqq, 
  !    int3 calculatec in newdq (depends upon self-consistency)
  !
  !    If a valid restart file is found:
  !    - dynmat0 is not called in any case
  !    - if rec_code = -20 the electric field calculation (solve_e) 
  !                   restarts from the saved value of the potential
  !    - if -10 < rec_code < -20 solve_e does nothing, the Raman calculation
  !                   (solve_e2), restarts from the saved value of the pot.
  !    - if rec_code > 0   the entire electric field and Raman section is not
  !                   called, the phonon calculation restarts from irrep irr
  !                   and from the saved value of the potential
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE phcom
  USE ph_restart,    ONLY : ph_readfile

  !
  implicit none
  !
  integer :: irr, ierr
  ! counter on representations
  ! error code
  logical :: exst, recover_file
  character(len=256) :: filename

  ierr=0
  IF (recover) THEN 
     CALL ph_readfile('data_u',ierr)
     IF (ierr==0) CALL ph_readfile('data',ierr)
     IF (where_rec=='solve_e...') THEN
        WRITE( stdout, '(/,4x," Restart in Electric Field calculation")')
     ELSEIF (where_rec=='solve_e2..') then
        WRITE( stdout, '(/,4x," Restart in Raman calculation")') 
     ELSEIF (where_rec=='solve_lint'.OR.where_rec=='done_drhod') then
        WRITE( stdout, '(/,4x," Restart in Phonon calculation")')
     ELSE
        call errore ('phq_recover', 'wrong restart data file', -1)
        ierr=1
     ENDIF
  ENDIF
!
!  open the recover file
!
  iunrec = 99
  call seqopn (iunrec, 'recover', 'unformatted', exst)
  recover_file = recover .AND. exst .AND. ierr==0
  if (.not.recover_file) close (unit = iunrec, status = 'delete')

  recover=recover_file

  RETURN
END SUBROUTINE phq_recover
