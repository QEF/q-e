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
#if defined (__XSD)
!-----------------------------------------------------------------
!
!------------------------------------------------------------------------
!       START_GLOBAL_VARIABLES ( INTENT (IN) ) 
!--------------------------------------------------------------------------
USE ions_base,    ONLY: tau, nat, nsp, atm, ityp
USE cell_base,    ONLY: alat, at
USE ener,         ONLY: etot, eband, ehart, etxc, vtxc, ewld, demet, ef 
USE klist,        ONLY: degauss, tot_charge
USE force_mod,    ONLY: force, sigma
USE control_flags,ONLY: nstep, n_scf_steps, scf_error
USE fcp_variables,ONLY: fcp_mu, lfcpopt, lfcpdyn 
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
IF ( lfcpopt .OR. lfcpdyn ) THEN 
   CALL qexsd_step_addstep ( i_step, nstep, nsp, atm, ityp, nat, tau, alat, at(:,1),   &
                          at(:,2), at(:,3), etot, eband, ehart, vtxc, etxc, ewld, degauss, demet, force, sigma,&
                          n_scf_steps, scf_error, &
                          POTSTAT_CONTR = (ef * tot_charge),  FCP_FORCE  = (fcp_mu-ef) , FCP_TOT_CHARGE = tot_charge)
ELSE 
   CALL qexsd_step_addstep ( i_step, nstep, nsp, atm, ityp, nat, tau, alat, at(:,1), at(:,2), at(:,3), etot, eband, &
                             ehart, vtxc, etxc, ewld, degauss, demet, force, sigma, n_scf_steps, scf_error)
END IF 
#else
CONTINUE 
#endif    
END SUBROUTINE  add_qexsd_step
