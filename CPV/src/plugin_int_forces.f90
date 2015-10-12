!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_int_forces(force)
!----------------------------------------------------------------------------
!
USE kinds,             ONLY : DP
USE io_global,         ONLY : stdout
USE ions_base,         ONLY : nat, ityp
USE cell_base,         ONLY : omega
USE fft_base,          ONLY : dfftp
USE fft_interfaces,    ONLY : fwfft
USE electrons_base,    ONLY : nspin
USE gvect,             ONLY : ngm, nl, eigts1, eigts2, eigts3 
!
USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
REAL(DP), INTENT(INOUT) :: force(3,nat)
!
! aux is used to store a possible additional density
! now defined in real space
!
COMPLEX(DP), ALLOCATABLE :: auxg(:), auxr(:)
!
INTEGER  :: ipol, na
! counter on polarization
! counter on atoms
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_int_forces
