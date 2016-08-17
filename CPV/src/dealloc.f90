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
  USE uspp,       ONLY : beta, dbeta, qq
  USE core,       ONLY : rhocb
  !
  USE core,                 ONLY : deallocate_core
  USE uspp,                 ONLY : deallocate_uspp
  USE electrons_base,       ONLY : deallocate_elct
  USE efield_module,        ONLY : deallocate_efield
  USE ensemble_dft,         ONLY : deallocate_ensemble_dft
  USE cg_module,            ONLY : deallocate_cg
  USE gvect,                ONLY : deallocate_gvect
  USE gvecs,                ONLY : deallocate_gvecs
  USE gvecw,                ONLY : deallocate_gvecw
  USE smallbox_gvec,        ONLY : deallocate_smallbox_gvec
  USE local_pseudo,         ONLY : deallocate_local_pseudo
  USE qgb_mod,              ONLY : deallocate_qgb_mod
  USE betax,                ONLY : deallocate_betax
  USE wavefunctions_module, ONLY : deallocate_wavefunctions
  USE wannier_module,       ONLY : deallocate_wannier
  USE fft_types,            ONLY : fft_type_descriptor, fft_type_deallocate
  USE fft_smallbox_type,    ONLY : fft_box_deallocate
  USE fft_base,             ONLY : dfftp, dffts, dfftb
  USE electrons_module,     ONLY : deallocate_electrons
  USE ions_base,            ONLY : deallocate_ions_base
  ! USE polarization,         ONLY : deallocate_polarization ! obsolescent
  USE cp_main_variables,    ONLY : deallocate_mainvar
  USE pseudopotential,      ONLY : deallocate_pseudopotential
  USE ions_nose,            ONLY : ions_nose_deallocate
  USE metagga,              ONLY : deallocate_metagga
  USE ions_positions,       ONLY : deallocate_ions_positions
  USE kohn_sham_states,     ONLY : ks_states_closeup
  USE ldaU_cp,              ONLY : deallocate_lda_plus_u
  USE step_penalty,         ONLY : deallocate_step_pen
  USE fft_base,             ONLY : pstickdealloc

  !
  IMPLICIT NONE
  !
  !
  IF ( ALLOCATED( beta ) )     DEALLOCATE( beta )
  IF ( ALLOCATED( qq ) )       DEALLOCATE( qq )
  IF ( ALLOCATED( rhocb ) )    DEALLOCATE( rhocb )
  IF ( ALLOCATED( dbeta ) )    DEALLOCATE( dbeta )
  !
  CALL deallocate_mainvar()
  CALL deallocate_ions_positions()
  CALL deallocate_efield( )
  CALL deallocate_ensemble_dft()
  CALL deallocate_cg( )
  CALL deallocate_core()
  CALL deallocate_uspp()
  CALL deallocate_gvect()
  CALL deallocate_gvecs()
  CALL deallocate_gvecw()
  CALL deallocate_smallbox_gvec( )
  CALL deallocate_local_pseudo()
  CALL deallocate_qgb_mod()
  CALL deallocate_betax()
  !
  CALL fft_type_deallocate( dfftp )
  CALL fft_type_deallocate( dffts )
  CALL fft_box_deallocate( dfftb )
  CALL pstickdealloc( )
  !
  CALL deallocate_ions_base()
  !
  CALL deallocate_wavefunctions()
  CALL deallocate_wannier()
  !
  CALL deallocate_elct()
  CALL deallocate_electrons()
  ! CALL deallocate_polarization() ! obsolescent
  CALL deallocate_pseudopotential()
  !
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
