!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_forces()
  !----------------------------------------------------------------------------
  !
  !  
  USE mp_global,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_bcast
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : outdir
  !
  USE plugin_flags
  !
  USE cell_base,        ONLY : alat, at
  USE ions_base,        ONLY : tau, nat
  USE force_mod,        ONLY : force
  USE control_flags,    ONLY : istep
  !
  IMPLICIT NONE
  !
  REAL(DP) :: at_meta(3,3)
  REAL(DP), ALLOCATABLE :: tau_meta(:,:)
  !
  IF(use_plumed) then
    IF(ionode)THEN
      at_meta=alat*at;
      allocate(tau_meta(3,nat))
      tau_meta=alat*tau
      call meta_force_calculation(at_meta,istep,tau_meta(1,1),0,0,force(1,1),0,0,0)
      deallocate(tau_meta)
    ENDIF
    CALL mp_bcast(force, ionode_id, intra_image_comm)
  ENDIF
  !
END SUBROUTINE plugin_forces
