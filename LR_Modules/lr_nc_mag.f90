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
   COMPLEX (DP), ALLOCATABLE :: int1_nc_save(:, :, :, :, :, :)
   COMPLEX (DP), ALLOCATABLE :: int3_save(:, :, :, :, :, :)
   !
   CONTAINS
   !
!------------------------------------------------------------------------------
SUBROUTINE lr_apply_time_reversal(tr, first_iter, dvscfins)
!------------------------------------------------------------------------------
   !! Apply time reversal for DFPT in noncollinear magnetic systems
   !! If tr == .true.,  flip sign of the magnetic field
   !! If tr == .false., revert back to the original state.
   !---------------------------------------------------------------------------
   USE kinds,                ONLY : DP
   USE scf,                  ONLY : vrs
   USE uspp,                 ONLY : okvan, deeq_nc
   USE lrus,                 ONLY : int3_nc
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: tr
   !! If .true., flip magnetic field. If .false., revert back.
   LOGICAL, INTENT(IN) :: first_iter
   !! True if first iteration. Skip some calculation if true.
   COMPLEX(DP), POINTER, INTENT(INOUT) :: dvscfins(:, :, :)
   !! change of the scf potential (smooth part only, dffts)
   !
   ! Flip the sign of the magnetic field
   !
   IF (.NOT. first_iter) THEN
      dvscfins(:, 2:4, :) = -dvscfins(:, 2:4, :)
   ENDIF
   !
   !$acc kernels
   vrs(:, 2:4) = -vrs(:, 2:4)
   !$acc end kernels
   !
   IF (tr) THEN
      !
      ! Set to the time-reversed state
      !
      IF (okvan) THEN
         IF (.NOT. first_iter) THEN
            int3_nc(:,:,:,:,:) = int3_save(:,:,:,:,:,2)
         ENDIF
         !
         deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,2)
         !$acc update device(deeq_nc)
      ENDIF
      !
   ELSE ! tr == .FALSE.
      !
      ! Set back to the original state
      !
      IF (okvan) THEN
         IF (.NOT. first_iter) THEN
            int3_nc(:,:,:,:,:) = int3_save(:,:,:,:,:,1)
         ENDIF
         !
         deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,1)
         !$acc update device(deeq_nc)
      ENDIF
      !
   ENDIF ! tr
   !
END SUBROUTINE lr_apply_time_reversal
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
END MODULE lr_nc_mag
!------------------------------------------------------------------------------
