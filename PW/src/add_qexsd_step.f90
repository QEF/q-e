!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This routine just calls the routint qexsd_step_addstep which adds a new xml
! element to to the list of step run by pw. In this way the addstep routine in
! the qexsd_step_addstep routine does not depend on global variables. 
!  P.  Delugas April 2016
 
!----------------------------------------------------------------
SUBROUTINE add_qexsd_step(i_step)
#if defined (__OLDXML)
CONTINUE 
#else
!-----------------------------------------------------------------
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
USE control_flags,ONLY: nstep, n_scf_steps, scf_error
USE fcp_variables,ONLY: fcp_mu, lfcpopt, lfcpdyn 
USE extfield,     ONLY: gate, etotgatefield, tefield, etotefield   
!-----------------------------------------------------------------------------
!   END_GLOBAL_VARIABLES
!----------------------------------------------------------------------------- 
!
!-------------------------------------------------------------------------------
!   SUBROUTINES FROM MODULES
!-------------------------------------------------------------------------------
USE qexsd_module, ONLY: qexsd_step_addstep
!-------------------------------------------------------------------------------
!------------------------------------------------------------------------------
IMPLICIT NONE
! 
!-------------------------------------------------------------------------------
!                     START_INPUT_VARIABLES
!-------------------------------------------------------------------------------
INTEGER,INTENT(IN)        ::   i_step
!-------------------------------------------------------------------------------
!                    END_INPUT_VARIABLES
!-------------------------------------------------------------------------------- 
!
REAL(DP),TARGET             :: potstat_contr_tgt, fcp_force_tgt, fcp_tot_charge_tgt,&
                               demet_tgt, degauss_tgt, gatefield_en_tgt, efield_corr_tgt

REAL(DP),POINTER            :: potstat_contr_ptr, fcp_force_ptr, fcp_tot_charge_ptr,&
                                   demet_ptr, degauss_ptr, gatefield_en_ptr, efield_corr_ptr 
!
NULLIFY(potstat_contr_ptr, fcp_force_ptr, fcp_tot_charge_ptr, demet_ptr, degauss_ptr, &
        gatefield_en_ptr, efield_corr_ptr)
!
IF ( degauss > 0.0d0 ) THEN 
   degauss_tgt = degauss/e2
   demet_tgt = demet/e2
   degauss_ptr => degauss_tgt
   demet_ptr    => demet_tgt
END IF    
IF ( lfcpopt .OR. lfcpdyn ) THEN  
   potstat_contr_tgt = ef * tot_charge / e2
   potstat_contr_ptr => potstat_contr_tgt
   !FIXME ( again shouldn't we use Hartree units for this ? )
   fcp_force_tgt = fcp_mu - ef
   fcp_force_ptr  => fcp_force_tgt
   !
   fcp_tot_charge_tgt = tot_charge
   fcp_tot_charge_ptr =>  fcp_tot_charge_tgt
   ! 
END IF
IF ( gate ) THEN 
   gatefield_en_tgt = etotgatefield/e2
   gatefield_en_ptr => gatefield_en_tgt
END IF
IF (tefield) THEN 
   efield_corr_tgt = etotefield/e2 
   efield_corr_ptr => efield_corr_tgt
END IF
CALL qexsd_step_addstep ( i_step, nstep, nsp, atm, ityp, nat, alat*tau, alat, alat*at(:,1),   &
                          alat*at(:,2), alat*at(:,3), etot/e2, eband/e2, ehart/e2, vtxc/e2, etxc/e2, &
                          ewld/e2, degauss_ptr, demet_ptr, force/e2, sigma/e2, n_scf_steps, scf_error, &
                          FCP_FORCE  = fcp_force_ptr , FCP_TOT_CHARGE = fcp_tot_charge_ptr,&
                          GATEFIELD_EN = gatefield_en_ptr) 
#endif    
!
END SUBROUTINE  add_qexsd_step
