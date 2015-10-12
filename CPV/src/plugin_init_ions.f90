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
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN) :: tau(3,nat)
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_init_ions
