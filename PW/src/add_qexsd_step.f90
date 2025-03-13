!
! Copyright (C) 2016-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------
SUBROUTINE add_qexsd_step(i_step)
  !-----------------------------------------------------------------
  !! This routine just calls the routine \(\texttt{qexsd_step_addstep}\)
  !! which adds a new xml element to to the list of steps run by PW.
  !! In this way the \(\texttt{addstep}\) routine in the 
  !! \(\texttt{qexsd_step_addstep}\) routine does not depend on global
  !! variables.  
  !! P. Delugas, April 2016.
  !
  !------------------------------------------------------------------------
  !       START_GLOBAL_VARIABLES ( INTENT (IN) ) 
  !--------------------------------------------------------------------------
  USE kinds,        ONLY: DP
  USE constants,    ONLY: e2
  USE ions_base,    ONLY: tau, nat, nsp, atm, ityp
  USE cell_base,    ONLY: alat, at
  USE ener,         ONLY: etot, eband, ehart, etxc, vtxc, ewld, demet, ef 
  USE klist,        ONLY: degauss, tot_charge
  USE force_mod,    ONLY: force, sigma
  USE control_flags,ONLY: nstep, n_scf_steps, scf_error, conv_elec
  USE fcp_module,   ONLY: fcp_mu, lfcp
  USE extfield,     ONLY: gate, etotgatefield, tefield, etotefield   
  USE control_flags,ONLY: max_xml_steps 
  !-----------------------------------------------------------------------------
  !   END_GLOBAL_VARIABLES
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !   SUBROUTINES FROM MODULES
  !-----------------------------------------------------------------------------
  USE qexsd_module, ONLY: qexsd_step_addstep
  !-----------------------------------------------------------------------------
  IMPLICIT NONE
  ! 
  !-----------------------------------------------------------------------------
  !                     START_INPUT_VARIABLES
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) ::   i_step
  !-----------------------------------------------------------------------------
  !                    END_INPUT_VARIABLES
  !-----------------------------------------------------------------------------
  !
  REAL(DP), TARGET :: potstat_contr_tgt, fcp_force_tgt, fcp_tot_charge_tgt,&
                      demet_tgt, degauss_tgt, gatefield_en_tgt, efield_corr_tgt

  REAL(DP),POINTER :: potstat_contr_ptr, fcp_force_ptr, fcp_tot_charge_ptr,&
                      demet_ptr, degauss_ptr, gatefield_en_ptr, efield_corr_ptr 
  !
  INTEGER :: stride = 1, max_xml_steps_  
  !
  IF ( max_xml_steps > 0 ) THEN 
     stride = nstep/max_xml_steps 
     IF (nstep > stride*max_xml_steps) stride = stride+1 
     max_xml_steps_ = max_xml_steps+2
  ELSE 
     max_xml_steps_ = nstep 
  END IF
  IF (.NOT. ( i_step == 1 .OR. MOD(i_step-1, stride) == 0 .OR. i_step == nstep)) RETURN
  !
  NULLIFY(potstat_contr_ptr, fcp_force_ptr, fcp_tot_charge_ptr, demet_ptr, &
       degauss_ptr, gatefield_en_ptr, efield_corr_ptr)
  !
  IF ( degauss > 0.0d0 ) THEN 
     degauss_tgt = degauss/e2
     demet_tgt = demet/e2
     degauss_ptr => degauss_tgt
     demet_ptr    => demet_tgt
  END IF
  !
  IF ( lfcp ) THEN
     potstat_contr_tgt = fcp_mu * tot_charge / e2
     potstat_contr_ptr => potstat_contr_tgt
     !
     fcp_force_tgt = (fcp_mu - ef)/e2
     fcp_force_ptr  => fcp_force_tgt
     !
     fcp_tot_charge_tgt = tot_charge
     fcp_tot_charge_ptr =>  fcp_tot_charge_tgt
     !
  END IF
  !
  IF ( gate ) THEN 
     gatefield_en_tgt = etotgatefield/e2
     gatefield_en_ptr => gatefield_en_tgt
  END IF
  !
  IF (tefield) THEN 
     efield_corr_tgt = etotefield/e2 
     efield_corr_ptr => efield_corr_tgt
  END IF
  !
  CALL qexsd_step_addstep ( i_step, max_xml_steps_, nsp, atm, ityp, nat, &
       alat*tau, alat, alat*at(:,1), alat*at(:,2), alat*at(:,3), etot/e2,&
       eband/e2, ehart/e2, vtxc/e2, etxc/e2, ewld/e2, degauss_ptr, demet_ptr,&
       force/e2, sigma/e2, conv_elec, n_scf_steps, scf_error, &
       FCP_FORCE  = fcp_force_ptr , FCP_TOT_CHARGE = fcp_tot_charge_ptr,&
       GATEFIELD_EN = gatefield_en_ptr) 
  !
END SUBROUTINE  add_qexsd_step
