!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE deallocate_modules_var()
  !----------------------------------------------------------------------------
  !
  USE uspp,       ONLY : beta, qq
  USE qgb_mod,    ONLY : qgb
  USE core,       ONLY : rhocb
  USE cdvan,      ONLY : dbeta
  USE dqgb_mod,   ONLY : dqgb
  !
  USE core,                 ONLY : deallocate_core
  USE cvan,                 ONLY : deallocate_cvan
  USE uspp,                 ONLY : deallocate_uspp
  USE electrons_base,       ONLY : deallocate_elct
  USE efield_module,        ONLY : deallocate_efield
  USE ensemble_dft,         ONLY : deallocate_ensemble_dft
  USE cg_module,            ONLY : deallocate_cg
  USE reciprocal_vectors,   ONLY : deallocate_recvecs
  USE recvecs_indexes,      ONLY : deallocate_recvecs_indexes
  USE local_pseudo,         ONLY : deallocate_local_pseudo
  USE qgb_mod,              ONLY : deallocate_qgb_mod
  USE dqgb_mod,             ONLY : deallocate_dqgb_mod
  USE betax,                ONLY : deallocate_betax
  USE wavefunctions_module, ONLY : deallocate_wavefunctions
  USE wannier_module,       ONLY : deallocate_wannier
  USE fft_types,            ONLY : fft_dlay_descriptor, fft_dlay_deallocate
  USE fft_types,            ONLY : fft_box_deallocate
  USE fft_base,             ONLY : dfftp, dffts, dfftb
  USE stick_base,           ONLY : sticks_deallocate
  USE electrons_module,     ONLY : deallocate_electrons
  USE ions_base,            ONLY : deallocate_ions_base
  USE sic_module,           ONLY : deallocate_sic
  USE polarization,         ONLY : deallocate_polarization
  USE cp_main_variables,    ONLY : deallocate_mainvar
  USE cdvan,                ONLY : deallocate_cdvan
  USE pseudopotential,      ONLY : deallocate_pseudopotential
  USE ions_nose,            ONLY : ions_nose_deallocate
  USE metagga,              ONLY : deallocate_metagga
  USE ions_positions,       ONLY : deallocate_ions_positions
  USE kohn_sham_states,     ONLY : ks_states_closeup
  USE ldau,                 ONLY : deallocate_lda_plus_u
  USE step_penalty,         ONLY : deallocate_step_pen

  !
  IMPLICIT NONE
  !
  !
  IF ( ALLOCATED( beta ) )     DEALLOCATE( beta )
  IF ( ALLOCATED( qgb ) )      DEALLOCATE( qgb )
  IF ( ALLOCATED( qq ) )       DEALLOCATE( qq )
  IF ( ALLOCATED( rhocb ) )    DEALLOCATE( rhocb )
  IF ( ALLOCATED( dqgb ) )     DEALLOCATE( dqgb )
  IF ( ALLOCATED( dbeta ) )    DEALLOCATE( dbeta )
  !
  CALL deallocate_mainvar()
  CALL deallocate_ions_positions()
  CALL deallocate_cvan()
  CALL deallocate_efield( )
  CALL deallocate_ensemble_dft()
  CALL deallocate_cg( )
  CALL deallocate_core()
  CALL deallocate_uspp()
  CALL deallocate_recvecs()
  CALL deallocate_recvecs_indexes()
  CALL deallocate_local_pseudo()
  CALL deallocate_qgb_mod()
  CALL deallocate_dqgb_mod()
  CALL deallocate_cdvan()
  CALL deallocate_betax()
  !
  CALL fft_dlay_deallocate( dfftp )
  CALL fft_dlay_deallocate( dffts )
  CALL fft_box_deallocate( dfftb )
  CALL sticks_deallocate()
  !
  CALL deallocate_ions_base()
  !
  CALL deallocate_wavefunctions()
  CALL deallocate_wannier()
  !
  CALL deallocate_elct()
  CALL deallocate_electrons()
  CALL deallocate_polarization()
  CALL deallocate_pseudopotential()
  !
  CALL deallocate_sic()
  CALL deallocate_metagga()
  CALL ions_nose_deallocate()
  CALL ks_states_closeup()
  !
  CALL deallocate_lda_plus_u()
  CALL deallocate_step_pen()
  !
  RETURN
  !
END SUBROUTINE deallocate_modules_var
