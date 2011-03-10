!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_initialization()
  !----------------------------------------------------------------------------
  !
  !
  USE io_global,        ONLY : stdout, ionode
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : outdir
  !
  USE plugin_flags
  !
  USE ions_base,        ONLY : amass, ityp, nat
  !
  USE time_step,  ONLY : delt
  !
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: mass(:),charge(:)
  INTEGER  :: na
  !
  IF ( use_plumed .and. ionode ) THEN

    ALLOCATE(mass(nat),charge(nat))

    DO na = 1, nat
    !
    mass(na)   = amass( ityp(na) )
    charge(na) = amass( ityp(na) )
    !
    END DO

    CALL init_metadyn(nat,delt,mass,charge,1,1.0D0,trim(outdir)//"plumed.dat"//char(0));

    DEALLOCATE(mass,charge)
  ENDIF
  !
END SUBROUTINE plugin_initialization
