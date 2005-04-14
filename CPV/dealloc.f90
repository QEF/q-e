!
! Copyright (C) 2002-2004 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!
!=======================================================================
!
   SUBROUTINE deallocate_modules_var()
     
      use control_flags, ONLY: lneb

      use core, only: deallocate_core
      use cvan, only: deallocate_cvan
      use uspp, only: deallocate_uspp
      use electrons_base, only: deallocate_elct
      use efield_module, only: deallocate_efield
      use ensemble_dft, only: deallocate_ensemble_dft
      use cg_module, only: deallocate_cg
      use reciprocal_vectors, only: deallocate_recvecs
      use recvecs_indexes,    only: deallocate_recvecs_indexes
      use pseu, only: deallocate_pseu
      use qgb_mod, only: deallocate_qgb_mod
      use dqgb_mod, only: deallocate_dqgb_mod
      use qradb_mod, only: deallocate_qradb_mod
      use dqrad_mod, only: deallocate_dqrad_mod
      use betax, only: deallocate_betax
      use wavefunctions_module, only: deallocate_wavefunctions
      use wannier_module, only: deallocate_wannier
      USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_deallocate
      USE fft_base, ONLY: dfftp, dffts
      use stick_base, only: sticks_deallocate
      USE electrons_module, ONLY: deallocate_electrons
      USE diis, ONLY: deallocate_diis
      USE charge_mix, ONLY: deallocate_charge_mix
      USE chi2, ONLY: deallocate_chi2
      USE ions_base, ONLY: deallocate_ions_base
      USE sic_module, ONLY: deallocate_sic
      USE polarization, ONLY: deallocate_polarization
      USE turbo, ONLY: deallocate_turbo
      USE ions_module, ONLY: deallocate_ions
      USE cp_main_variables, ONLY: deallocate_mainvar
      USE derho, ONLY: deallocate_derho
      USE dpseu, ONLY: deallocate_dpseu
      USE cdvan, ONLY: deallocate_cdvan
      USE pseudopotential, ONLY: deallocate_pseudopotential

      IMPLICIT NONE

      CALL deallocate_mainvar()
      CALL deallocate_cvan()
      CALL deallocate_efield( )
      CALL deallocate_ensemble_dft()
      CALL deallocate_cg( )
      CALL deallocate_elct()
      CALL deallocate_core()
      CALL deallocate_uspp()
      CALL deallocate_recvecs()
      CALL deallocate_recvecs_indexes()
      CALL deallocate_pseu()
      CALL deallocate_qgb_mod()
      CALL deallocate_qradb_mod()
      CALL deallocate_derho()
      CALL deallocate_dqgb_mod()
      CALL deallocate_dpseu()
      CALL deallocate_cdvan()
      CALL deallocate_dqrad_mod()
      CALL deallocate_betax()

      call fft_dlay_deallocate( dfftp )
      call fft_dlay_deallocate( dffts )
      call sticks_deallocate()

      CALL deallocate_wavefunctions()
      CALL deallocate_wannier()

      CALL deallocate_electrons()
      CALL deallocate_polarization()
      CALL deallocate_pseudopotential()
      CALL deallocate_turbo()
      CALL deallocate_diis()

      CALL deallocate_charge_mix()
      CALL deallocate_chi2()
      IF( .NOT. lneb ) THEN
        CALL deallocate_ions_base
      END IF
      CALL deallocate_ions( )
      CALL deallocate_sic( )


     RETURN
   END SUBROUTINE
