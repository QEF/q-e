!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE punch( what )
  !----------------------------------------------------------------------------
  !
  ! ... This routine is called at the end of the run to save to a file
  ! ... the information needed for further processing (phonon etc.)
  !
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, nkstot
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho
  USE control_flags,        ONLY : reduce_io, lscf, lbands
  USE wvfct,                ONLY : et, wg, nbnd
  USE wavefunctions_module, ONLY : evc, evc_nc
  USE io_files,             ONLY : prefix, iunpun, iunwfc, nwordwfc
  USE noncollin_module,     ONLY : noncolin
  USE restart_module,       ONLY : writefile_new
  USE mp_global,            ONLY : kunit
  USE pw_restart,           ONLY : pw_writefile
  USE a2F,                  ONLY : la2F
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: what
  INTEGER          :: ik, i, ibnd, kunittmp
  LOGICAL          :: exst
  !
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing output data file ",A)' ) &
      TRIM( prefix ) // '.save'
  !
  kunittmp = 1
  !
  ! ... if the wavefunction has not been written on file, do it now
  !
  IF ( noncolin ) THEN
     !
     IF ( nks == 1 .AND. reduce_io ) &
        CALL davcio( evc_nc, nwordwfc, iunwfc, 1, +1 )
     !
  ELSE
     !
     IF ( nks == 1 .AND. reduce_io ) &
        CALL davcio( evc, nwordwfc, iunwfc, 1, +1 )
     !
  ENDIF
  !
  ! ... The following instruction is used  when more k-points are needed
  ! ... for finite-q phonon calculations (on fine q-grid) then those needed
  ! ... for self-consistency. In such a case, a self-consistent calculation
  ! ... with few k-points is followed by a non-self-consistent one with added
  ! ... k-points, whose weight is set to zero.
  !
  IF ( .NOT. lscf .AND. .NOT. lbands ) CALL sum_band()
  !
  ! ... Write: general variables (including dimensions of the arrays),
  ! ... atomic positions, forces, k-points, eigenvalues
  !
#if defined (__PARA)
  !
  ! ... xk, wk, isk, et, wg are distributed across pools
  ! ... the first node has a complete copy of xk, wk, isk,
  ! ... while eigenvalues et and weights wg must be
  ! ... explicitely collected to the first node
  !
  CALL poolrecover( et, nbnd, nkstot, nks )
  CALL poolrecover( wg, nbnd, nkstot, nks )
  !
  ! ... In parallel execution, only the first node writes this file
  !
  kunittmp = kunit
  !
#endif
  !
  ! ... Write the charge density on a separate file
  !
  CALL io_pot( + 1, 'rho', rho, nspin )
  !
  iunpun = 4
  !
#if defined (__OLDPUNCH)
  !
  CALL writefile_new( TRIM( what ), iunpun, et, wg, kunittmp )
  !
#else
  !
  CALL pw_writefile( TRIM( what ) )
  !
#endif
  !
  IF ( la2F ) CALL enfdos()
  !
  RETURN
  !
END SUBROUTINE punch
