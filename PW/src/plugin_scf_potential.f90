!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_scf_potential(rhoin,conv_elec,dr2,vltot)
!----------------------------------------------------------------------------
!
! This routine is used to calculate plugin contributions to Kohn-Sham
! potential
!
USE io_global,        ONLY : stdout, ionode
USE kinds,            ONLY : DP
USE fft_base,         ONLY : dfftp
USE lsda_mod,         ONLY : nspin
USE scf,              ONLY : scf_type
USE plugin_flags
!
!
IMPLICIT NONE
!
TYPE(scf_type), INTENT(IN) :: rhoin
LOGICAL, INTENT(IN) :: conv_elec
REAL(DP), INTENT(IN) :: dr2
REAL(DP), INTENT(INOUT) :: vltot(dfftp%nnr)
!
!
END SUBROUTINE plugin_scf_potential
