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
REAL(DP),ALLOCATABLE            :: potstat_contr_(:), fcp_force_(:), fcp_tot_charge_(:),&
                                   demet_(:), degauss_(:), gatefield_en_(:), efield_corr_(:) 
!
IF ( degauss > 0.0d0 ) THEN 
   ALLOCATE ( degauss_(1), demet_(1)) 
   degauss_ = degauss/e2
   demet_   = demet/e2
END IF    
IF ( lfcpopt .OR. lfcpdyn ) THEN 
   ALLOCATE ( potstat_contr_(1), fcp_force_(1), fcp_tot_charge_(1)) 
   potstat_contr_ = ef * tot_charge / e2 
   !FIXME ( again shouldn't we use Hartree units for this ? ) 
   fcp_force_  = fcp_mu - ef 
   !
   fcp_tot_charge_ = tot_charge 
END IF
IF ( gate ) THEN 
   ALLOCATE (gatefield_en_(1)) 
   gatefield_en_ = etotgatefield/e2 
END IF
IF (tefield) THEN 
   ALLOCATE ( efield_corr_(1))
   efield_corr_(1) = etotefield/e2 
END IF
CALL qexsd_step_addstep ( i_step, nstep, nsp, atm, ityp, nat, alat*tau, alat, alat*at(:,1),   &
                          alat*at(:,2), alat*at(:,3), etot, eband/e2, ehart/e2, vtxc/e2, etxc/e2, &
                          ewld/e2, degauss_, demet_, force/e2, sigma/e2, n_scf_steps, scf_error, &
                          EFIELDCORR = efield_corr_, POTSTAT_CONTR = potstat_contr_,            &
                          FCP_FORCE  = fcp_force_ , FCP_TOT_CHARGE = fcp_tot_charge_,            &
                          GATEFIELD_EN = gatefield_en_) 
#endif    
!
IF (ALLOCATED(potstat_contr_)) DEALLOCATE(potstat_contr_)
IF (ALLOCATED(fcp_force_))  DEALLOCATE(fcp_force_)
IF (ALLOCATED(fcp_tot_charge_)) DEALLOCATE ( fcp_tot_charge_)
IF (ALLOCATED(demet_)) DEALLOCATE (demet_)
IF (ALLOCATED(degauss_)) DEALLOCATE(degauss_) 
IF (ALLOCATED(gatefield_en_)) DEALLOCATE(gatefield_en_)
IF (ALLOCATED(efield_corr_)) DEALLOCATE(efield_corr_) 
END SUBROUTINE  add_qexsd_step
