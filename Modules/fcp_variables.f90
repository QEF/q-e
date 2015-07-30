!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
MODULE fcp_variables
!------------------------------------------------------------------------------!
   !
   USE kinds,      ONLY : DP
   !
   IMPLICIT NONE
   SAVE
   !
   LOGICAL  :: lfcpopt              = .FALSE.
   LOGICAL  :: lfcpdyn              = .FALSE.
   REAL(DP) :: fcp_mu               = 0.0_DP
   REAL(DP) :: fcp_mass             = 10000.0_DP
   REAL(DP) :: fcp_temperature      = 0.0_DP
   REAL(DP) :: fcp_relax_step       = 0.5_DP
   REAL(DP) :: fcp_relax_crit       = 0.001_DP
   REAL(DP) :: fcp_tot_charge_first = 0.0_DP
   REAL(DP) :: fcp_tot_charge_last  = 0.0_DP

END MODULE fcp_variables
