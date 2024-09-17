!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
MODULE lr_nc_mag
   !
   USE kinds, ONLY : DP
   !
   IMPLICIT NONE
   PUBLIC
   SAVE
   !
   COMPLEX (DP), ALLOCATABLE :: deeq_nc_save(:, :, :, :, :)
   !! deeq_nc_save(:,:,:,:,1) stores the original state of deeq_nc.
   !! deeq_nc_save(:,:,:,:,2) stores deeq_nc computed with flipped magnetic field.
   COMPLEX (DP), ALLOCATABLE :: int1_nc_save(:, :, :, :, :, :)
   !! Same as deeq_nc_save, for int1_nc
   COMPLEX (DP), ALLOCATABLE :: int3_nc_save(:, :, :, :, :, :)
   !! Same as deeq_nc_save, for int3_nc
   !
   !
   CONTAINS
   !
!------------------------------------------------------------------------------
SUBROUTINE lr_apply_time_reversal(first_iter, ind, dvscfins)
!------------------------------------------------------------------------------
   !! Apply time reversal for DFPT in noncollinear magnetic systems. Used in sternheimer_kernel.
   !! If ind == 2, flip to the time-reversed state. Called during initialization.
   !! If ind == 1, revert to the original state. Called during finalization.
   !---------------------------------------------------------------------------
   USE kinds,             ONLY : DP
   USE scf,               ONLY : vrs
   USE uspp,              ONLY : okvan, deeq_nc
   USE noncollin_module,  ONLY : noncolin, domag
   USE lrus,              ONLY : int3_nc
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! True if first iteration. Skip some calculation if true.
   INTEGER, INTENT(IN) :: ind
   !! If 1, flip to the time-reversed state. If 1, revert back to the original state.
   COMPLEX(DP), POINTER, INTENT(INOUT) :: dvscfins(:, :, :)
   !! change of the scf potential (smooth part only, dffts)
   !
   IF (.NOT. (ind == 1 .OR. ind == 2)) CALL errore('lr_apply_time_reversal', &
      'ind must be 1 or 2', 1)
   !
   IF (.NOT. (noncolin.AND.domag)) CALL errore('lr_apply_time_reversal', &
      'This routine is only for noncollinear magnetic systems', 1)
   !
   ! Flip the sign of the magnetic field
   !
   IF (.NOT. first_iter) THEN
      dvscfins(:, 2:4, :) = -dvscfins(:, 2:4, :)
      IF (okvan) int3_nc(:,:,:,:,:) = int3_nc_save(:,:,:,:,:,ind)
   ENDIF
   !
   !$acc kernels
   vrs(:, 2:4) = -vrs(:, 2:4)
   !$acc end kernels
   !
   IF (okvan) THEN
      deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,ind)
      !$acc update device(deeq_nc)
   ENDIF
   !
END SUBROUTINE lr_apply_time_reversal
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
END MODULE lr_nc_mag
!------------------------------------------------------------------------------
