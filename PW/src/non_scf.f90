!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!-----------------------------------------------------------------------
  SUBROUTINE non_scf (ik_)
  !-----------------------------------------------------------------------
  !
  ! ... diagonalization of the KS hamiltonian in the non-scf case
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield, lberry
  USE control_flags,        ONLY : lbands, io_level
  USE ener,                 ONLY : ef
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunwfc, nwordwfc, iunefield
  USE buffers,              ONLY : save_buffer
  USE klist,                ONLY : xk, wk, nks, nkstot
  USE lsda_mod,             ONLY : lsda, nspin
  USE wvfct,                ONLY : nbnd, et, npwx
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: ik_
  !
  ! ... local variables
  !
  INTEGER :: iter = 1, i, ik
  REAL(DP) :: dr2 = 0.d0
  REAL(DP), EXTERNAL :: get_clock
  !
  !
  CALL start_clock( 'electrons' )
  !
  WRITE( stdout, 9002 )
  !
  CALL flush_unit( stdout )
  !
  IF ( lelfield) THEN
     !
     CALL c_bands_efield ( iter, ik_, dr2 )
     !
  ELSE
     !
     CALL c_bands_nscf ( ik_ )
     !
  END IF
  !
  ! ... xk, wk, isk, et, wg are distributed across pools;
  ! ... the first node has a complete copy of xk, wk, isk,
  ! ... while eigenvalues et and weights wg must be
  ! ... explicitly collected to the first node
  ! ... this is done here for et, in weights () for wg
  !
  CALL poolrecover( et, nbnd, nkstot, nks )
  !
  ! ... calculate weights of Kohn-Sham orbitals 
  ! ... may be needed in further calculations such as phonon
  !
  IF ( .NOT. lbands ) CALL weights  ( )
  !
  ! ... Note that if you want to use more k-points for the phonon
  ! ... calculation then those needed for self-consistency, you can,
  ! ... by performing a scf with less k-points, followed by a non-scf
  ! ... one with additional k-points, whose weight on input is set to zero
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  WRITE( stdout, 9102 )
  !
  ! ... write band eigenvalues
  !
  CALL print_ks_energies ( ) 
  !
  ! ... save converged wfc if they have not been written previously
  !
  IF ( nks == 1 .AND. (io_level < 2) .AND. (io_level > -1) ) &
        CALL save_buffer ( evc, nwordwfc, iunwfc, nks )
  !
  ! ... do a Berry phase polarization calculation if required
  !
  IF ( lberry ) CALL c_phase()
  !
  CALL stop_clock( 'electrons' )
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9102 FORMAT(/'     End of band structure calculation' )
  !
END SUBROUTINE non_scf
