!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_scf_potential(rhoin,conv_elec,dr2)
!----------------------------------------------------------------------------
! This routine is used to calculate plugin energy related quantities
! that needs to be solved inside the scf cycle
!
USE io_global,        ONLY : stdout, ionode
USE kinds,            ONLY : DP
USE io_files,         ONLY : tmp_dir
!
USE fft_base,         ONLY : dfftp
USE lsda_mod,         ONLY : nspin
USE scf,              ONLY : scf_type, vltot
!
USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
type(scf_type), intent(in) :: rhoin
LOGICAL, intent(in) :: conv_elec
real(DP), intent(in) :: dr2
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_scf_potential
