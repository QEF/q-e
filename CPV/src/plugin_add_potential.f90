!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE plugin_add_potential( v )
!----------------------------------------------------------------------------
! This routine is used to add the plugin potentials to the total electronic potential
!
USE io_global,        ONLY : stdout, ionode
USE kinds,            ONLY : DP
!
USE fft_base,         ONLY : dfftp
USE electrons_base,   ONLY : nspin
!
USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin)
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***

END SUBROUTINE plugin_add_potential
