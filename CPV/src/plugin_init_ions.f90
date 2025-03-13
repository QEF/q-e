!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_init_ions( tau )
!----------------------------------------------------------------------------
!
USE kinds,            ONLY : DP
USE fft_base,         ONLY : dfftp
USE ions_base,        ONLY : nat, nsp, na
USE plugin_flags
!
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: tau(3,nat)
!
!
END SUBROUTINE plugin_init_ions
