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
  USE cell_base,        ONLY : at, h, alat
  USE cell_base,        ONLY : s_to_r, r_to_s
  USE ions_base,        ONLY : nat, na, nsp
  USE ions_positions,   ONLY : taus, fion
  USE cp_main_variables, ONLY : nfi
  !
  IMPLICIT NONE
  !
  REAL(DP) :: at_meta(3,3)
  REAL(DP), ALLOCATABLE :: tau_meta(:,:)
  INTEGER :: istep
  !
  IF(use_plumed) then
     IF ( ionode ) THEN
       at_meta(:,:)=alat*at(:,:)

       allocate(tau_meta(3,nat))
       tau_meta(:,:) = 0.0D0

       call s_to_r(taus,tau_meta,na,nsp,h)

       istep = nfi
       !
       ! convert to Rydberg as BOLTZMAN constant has been defined in Ry 
       ! for quantum ESPRESSO
       !
       fion(:,:) = fion(:,:) * 2.0D0
       !
       call meta_force_calculation(at_meta,istep,tau_meta(1,1),0,0,fion(1,1),0,0,0)
       !
       ! convert to Hartree (cp internal units)
       !
       fion(:,:) = fion(:,:) / 2.0D0
       !
       deallocate(tau_meta)
     END IF

     CALL mp_bcast( fion, ionode_id, intra_image_comm )
  ENDIF
  !
END SUBROUTINE plugin_forces
