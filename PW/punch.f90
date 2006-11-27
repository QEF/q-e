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
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : prefix, iunpun, iunwfc, nwordwfc
  USE noncollin_module,     ONLY : noncolin
  USE pw_restart,           ONLY : pw_writefile
  USE a2F,                  ONLY : la2F, a2Fsave
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: what
  LOGICAL          :: exst
  !
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing output data file ",A)' ) &
      TRIM( prefix ) // '.save'
  !
  ! ... if the wavefunction has not been written on file, do it now
  !
  IF ( nks == 1 .AND. reduce_io ) &
     CALL davcio( evc, nwordwfc, iunwfc, 1, +1 )
  !
  ! ... The following instruction is used when phonons are calculated
  ! ... one q-wavevector at the time. Weights and occupations have to be
  ! ... recomputed. Note that if you want to use more k-points for the
  ! ... phonon calculation then those needed for self-consistency, you
  ! ... can, by performing a scf with less k-points, followed by a non-scf
  ! ... one with additional k-points, whose weight on input is set to zero
  !
  IF ( .NOT. lscf .AND. .NOT. lbands ) CALL sum_band ( )
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
#endif
  !
  iunpun = 4
  !
  CALL pw_writefile( TRIM( what ) )
  !
  IF ( la2F ) CALL a2Fsave()
  !
  RETURN
  !
END SUBROUTINE punch
