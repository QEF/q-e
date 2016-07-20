!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!-----------------------------------------------------------------------
  SUBROUTINE non_scf ( )
  !-----------------------------------------------------------------------
  !
  ! ... diagonalization of the KS hamiltonian in the non-scf case
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield, lberry, lorbm
  USE check_stop,           ONLY : stopped_by_user
  USE control_flags,        ONLY : io_level, conv_elec, lbands
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
  ! ... local variables
  !
  INTEGER :: iter, i
  REAL(DP), EXTERNAL :: get_clock
  !
  !
  CALL start_clock( 'electrons' )
  iter = 1
  !
  WRITE( stdout, 9002 )
  FLUSH( stdout )
  !
  IF ( lelfield) THEN
     !
     CALL c_bands_efield ( iter )
     !
  ELSE
     !
     CALL c_bands_nscf ( )
     !
  END IF
  !
  ! ... check if calculation was stopped in c_bands
  !
  IF ( stopped_by_user ) THEN
     conv_elec=.FALSE.
     RETURN
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
  ! ... calculate weights of Kohn-Sham orbitals (only weights, not Ef,
  ! ... for a "bands" calculation where Ef is read from data file)
  ! ... may be needed in further calculations such as phonon
  !
  IF ( lbands ) THEN
     CALL weights_only  ( )
  ELSE
     CALL weights  ( )
  END IF
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
  ! ... write band eigenvalues (conv_elec is used in print_ks_energies)
  !
  conv_elec = .true.
  CALL print_ks_energies ( ) 
  !
  ! ... save converged wfc if they have not been written previously
  ! ... FIXME: it shouldn't be necessary to do this here
  !
  IF ( nks == 1 .AND. (io_level < 2) .AND. (io_level > -1) ) &
        CALL save_buffer ( evc, nwordwfc, iunwfc, nks )
  !
  ! ... do a Berry phase polarization calculation if required
  !
  IF ( lberry ) CALL c_phase()
  !
  ! ... do an orbital magnetization (Kubo terms) calculation
  !
  IF ( lorbm ) CALL orbm_kubo()
  !
  CALL stop_clock( 'electrons' )
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9102 FORMAT(/'     End of band structure calculation' )
  !
END SUBROUTINE non_scf

