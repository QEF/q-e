!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_energy(rhoin,plugin_etot)
!----------------------------------------------------------------------------
! This routine is used to calculate plugin energy related quantities
! that needs to be solved inside the scf cycle
!
USE io_global,        ONLY : stdout, ionode
USE kinds,            ONLY : DP
USE io_files,         ONLY : tmp_dir
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
real(DP), intent(in) :: rhoin(dfftp%nnr,nspin)
real(DP), intent(inout) :: plugin_etot
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_energy
