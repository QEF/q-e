!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!-----------------------------------------------------------------------
SUBROUTINE non_scf( )
  !-----------------------------------------------------------------------
  !! Diagonalization of the KS hamiltonian in the non-scf case.
  !
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : lelfield, lberry, lorbm
  USE check_stop,           ONLY : stopped_by_user
  USE control_flags,        ONLY : io_level, conv_elec, lbands, ethr, use_gpu
  USE ener,                 ONLY : ef, ef_up, ef_dw, eband
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : save_buffer
  USE klist,                ONLY : xk, wk, nks, nkstot, two_fermi_energies
  USE lsda_mod,             ONLY : lsda, nspin
  USE wvfct,                ONLY : nbnd, et, npwx
  USE wavefunctions,        ONLY : evc
  USE add_dmft_occ,         ONLY : dmft
  !
  USE exx,                  ONLY : exxinit, aceinit, use_ace
  USE scf,                  ONLY : rho, rho_core, rhog_core, v, vltot, vrs, kedtau
  USE ener,                 ONLY : ehart, etxc, vtxc, epaw
  USE ldaU,                 ONLY : eth
  USE extfield,             ONLY : etotefield
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_variables,        ONLY : okpaw, ddd_paw
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : doublegrid
  USE ions_base,            ONLY : nat
  USE xc_lib,               ONLY : stop_exx, xclib_dft_is
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: iter, i, dr2 = 0.0_dp
  REAL(dp):: ef_scf, ef_scf_up, ef_scf_dw
  REAL(DP), EXTERNAL :: e_band, get_clock
  REAL(DP) :: charge
  REAL(DP) :: etot_cmp_paw(nat,2,2)
  !
  !
  CALL start_clock( 'electrons' )
  iter = 1
  !
  WRITE( stdout, 9002 )
  FLUSH( stdout )
  !
  IF ( lelfield ) THEN
     !
     CALL c_bands_efield( iter )
     !
  ELSE
     !
     CALL c_bands_nscf()
     !
  ENDIF
  !
  ! ... check if calculation was stopped in c_bands
  !
  IF ( stopped_by_user ) THEN
     conv_elec=.FALSE.
     RETURN
  ENDIF
  !
  ! ... xk, wk, isk, et, wg are distributed across pools;
  ! ... the first node has a complete copy of xk, wk, isk,
  ! ... while eigenvalues et and weights wg must be
  ! ... explicitly collected to the first node
  ! ... this is done here for et, in weights () for wg
  !
  CALL poolrecover( et, nbnd, nkstot, nks )
  !
  ! ... the new density is computed here. For PAW:
  ! ... sum_band computes new becsum (stored in uspp modules)
  ! ... and a subtly different copy in rho%bec (scf module)
  !
  IF(xclib_dft_is('hybrid')) THEN 
    CALL sum_band()
  END IF 
  !
  ! ... calculate weights of Kohn-Sham orbitals (only weights, not Ef,
  ! ... for a "bands" calculation where Ef is read from data file)
  ! ... may be needed in further calculations such as phonon
  !
  ! save Fermi energy read from scf calculation
  ef_scf = ef
  ef_scf_up = ef_up
  ef_scf_dw = ef_dw
  IF ( lbands ) THEN
     CALL weights_only( )
  ELSE
     CALL weights( )
     eband =  e_band( )
  ENDIF
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
  ! ... if Ef is re-computed: print original Ef as well
  !
  conv_elec = .TRUE.
  CALL print_ks_energies_nonscf ( ef_scf, ef_scf_up, ef_scf_dw ) 
  !
  ! ... save converged wfc if they have not been written previously
  ! ... FIXME: it shouldn't be necessary to do this here
  !
  IF ( nks == 1 .AND. (io_level < 2) .AND. (io_level > -1) ) &
        CALL save_buffer( evc, nwordwfc, iunwfc, nks )
  !
  ! ... do a Berry phase polarization calculation if required
  !
  IF ( lberry ) CALL c_phase()
  !
  ! ... do an orbital magnetization (Kubo terms) calculation
  !
  IF ( lorbm ) CALL orbm_kubo()
  !
  ! ... for DMFT write everything to file to restart next scf step from here
  !
  IF ( dmft ) THEN
     CALL save_in_electrons( iter-1, dr2, ethr, et )
     RETURN
  ENDIF
  !
  ! ... for exact exchange case update the ACE projector with the actual number of bands 
  !
  IF(xclib_dft_is('hybrid')) THEN 
     !CALL save_buffer( evc, nwordwfc, iunwfc, nks )
     ! I want exx_is_active to be false inside exxinit to allow all relevant initializations
     CALL stop_exx() 
     CALL exxinit(.false., nbnd)
     IF (use_ace) CALL aceinit ( .false. )
     CALL v_of_rho( rho, rho_core, rhog_core, &
         ehart, etxc, vtxc, eth, etotefield, charge, v)
     IF (okpaw) CALL PAW_potential(rho%bec, ddd_paw, epaw,etot_cmp_paw)
     CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &
                   nspin, doublegrid )
     !
     WRITE(stdout,'(5x,"Calculation (EXX) restarted with the new ACE potential")' ) 
     !
     conv_elec = .FALSE.
     CALL c_bands_nscf()
     !
     IF ( stopped_by_user ) THEN
        conv_elec=.FALSE.
        RETURN
     ENDIF
     CALL poolrecover( et, nbnd, nkstot, nks )
     ef_scf = ef
     ef_scf_up = ef_up
     ef_scf_dw = ef_dw
     IF ( lbands ) THEN
        CALL weights_only( )
     ELSE
        CALL weights( )
     ENDIF
     WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
     WRITE( stdout, 9102 )
     conv_elec = .TRUE.
     CALL print_ks_energies_nonscf ( ef_scf, ef_scf_up, ef_scf_dw ) 
     IF ( nks == 1 .AND. (io_level < 2) .AND. (io_level > -1) ) &
           CALL save_buffer( evc, nwordwfc, iunwfc, nks )
     IF ( lberry ) CALL c_phase()
     IF ( lorbm ) CALL orbm_kubo()
  END IF
  !
  CALL stop_clock( 'electrons' )
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9102 FORMAT(/'     End of band structure calculation' )
  !
END SUBROUTINE non_scf
