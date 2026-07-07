!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
MODULE dfpt_type
   !
   USE kinds, ONLY : DP
   !
   IMPLICIT NONE
   !
   ! Data that describes linear response quantities.
   !
   TYPE :: dfpt_data_type
      INTEGER :: npert
      !! Number of perturbations
      COMPLEX(DP), ALLOCATABLE :: drhos(:, :, :)
      !! Change of the charge density (smooth part only, dffts)
      !! If norm conserving PP and doublegrid = .FALSE., identical to drhop.
      !! If USPP or PAW, does not include the augmentation charges. (drhop does.)
      COMPLEX(DP), ALLOCATABLE :: drhop(:, :, :)
      !! Change of the charge density (smooth and hard parts, dfftp)
      COMPLEX(DP), ALLOCATABLE :: dvscfs(:, :, :)
      !! Change of the scf potential (smooth part only, dffts)
      !! If doublegrid = .FALSE., identical to dvscfp, even for USPP or PAW.
      COMPLEX(DP), ALLOCATABLE :: dvscfp(:, :, :)
      !! Change of the scf potential (smooth and hard parts, dfftp)
      COMPLEX(DP), ALLOCATABLE :: dbecsum(:, :, :, :)
      !! Change of becsum. Only used for USPP and PAW.
      !
      ! Variables for metallic systems at q=0
      !
      COMPLEX(DP), ALLOCATABLE :: def(:)
      !! Change of the Fermi energy for each perturbation.
      COMPLEX(DP), ALLOCATABLE :: dn0(:)
      !! Bare charge perturbation for metallic systems at q=0. It is the Delta n_ext(q=0)
      !! term in Eq.(79) of Baroni et al, Rev. Mod. Phys. 73, 515 (2001).
      !! In phonon calculations, this term is zero. It is nonzero when calculating the
      !! screened Coulomb interaction.
      !
      ! Variables allocated and used only for phonon perturbations
      !
      COMPLEX(DP), ALLOCATABLE :: drhoc(:, :)
      !! Change of core charge due to nonlinear core correction.
      !! Size (dfftp%nnr, npert)
      COMPLEX(DP), ALLOCATABLE :: drhop_pulay(:, :, :)
      !! Pulay correction to drhop due to augmentation charge. Used for USPP or PAW.
      !! Size (dfftp%nnr, nspin_mag, npert)
      COMPLEX(DP), ALLOCATABLE :: dbecsum_pulay(:, :, :, :)
      !! Pulay correction to dbecsum due to augmentation charge. Used for PAW only.
      !! Size ((nhm * (nhm + 1))/2, nat, nspin_mag, npert)
   END TYPE dfpt_data_type
   !
   ! Data that describes local density of states at the Fermi level.
   !
   TYPE :: dfpt_ldos_type
      !! Local density of states (LDOS) at Fermi energy
      REAL(DP) :: dos_ef
      !! Total (integrated) density of states at Ef
      COMPLEX(DP), ALLOCATABLE :: ldoss(:, :)
      !! Local DOS at Ef (smooth part only, dffts)
      !! Contains only wavefunction contributions, no augmentation
      COMPLEX(DP), ALLOCATABLE :: ldos(:, :)
      !! Local DOS at Ef (smooth and hard parts, dfftp)
      !! Includes augmentation contributions from becsum_dos
      REAL(DP), ALLOCATABLE :: becsum_dos(:, :, :)
      !! Augmentation occupancy contribution to LDOS at Ef
      !! Size: ((nhm * (nhm + 1))/2, nat, nspin_mag)
      !! Only used for USPP and PAW
      !! For USPP, only used inside localdos, and then deallocated.
      !! For PAW, stored and used in dfpt_kernel as well.
      !! NOTE: Real-valued (unlike dbecsum in dfpt_data_type which is complex)
   END TYPE dfpt_ldos_type
   !
   CONTAINS
   !
   !---------------------------------------------------------------------------
   SUBROUTINE dfpt_dvscfp_to_dvscfs(dfpt_data)
      USE fft_base,              ONLY : dffts, dfftp
      USE fft_interfaces,        ONLY : fft_interpolate
      USE gvecs,                 ONLY : doublegrid
      USE noncollin_module,      ONLY : nspin_mag
      !
      IMPLICIT NONE
      !
      TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
      !
      INTEGER :: ipert, is
      !
      IF (doublegrid) THEN
         DO ipert = 1, dfpt_data%npert
            DO is = 1, nspin_mag
               CALL fft_interpolate(dfftp, dfpt_data%dvscfp(:, is, ipert), &
                                    dffts, dfpt_data%dvscfs(:, is, ipert))
            ENDDO
         ENDDO
      ELSE
         ! If doublegrid is false, dvscfs and dvscfp are the same.
         CALL zcopy(dfftp%nnr*nspin_mag*dfpt_data%npert, dfpt_data%dvscfp, 1, dfpt_data%dvscfs, 1)
      ENDIF
      !
   END SUBROUTINE dfpt_dvscfp_to_dvscfs
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   SUBROUTINE allocate_dfpt_data(dfpt_data, npert)
      !
      USE fft_base,             ONLY : dfftp, dffts
      USE noncollin_module,     ONLY : nspin_mag
      USE uspp_param,           ONLY : nhm
      USE ions_base,            ONLY : nat
      USE klist,                ONLY : ltetra, lgauss
      USE uspp,                 ONLY : okvan
      USE control_lr,           ONLY : lgamma
      !
      IMPLICIT NONE
      !
      TYPE(dfpt_data_type), INTENT(OUT) :: dfpt_data
      INTEGER, INTENT(IN) :: npert
      !
      LOGICAL :: lmetq0
      !! true if xq=(0,0,0) in a metal
      !
      dfpt_data%npert = npert
      !
      ALLOCATE(dfpt_data%drhos(dffts%nnr, nspin_mag, npert))
      ALLOCATE(dfpt_data%drhop(dfftp%nnr, nspin_mag, npert))
      ALLOCATE(dfpt_data%dvscfs(dffts%nnr, nspin_mag, npert))
      ALLOCATE(dfpt_data%dvscfp(dfftp%nnr, nspin_mag, npert))
      dfpt_data%drhos = (0.d0, 0.d0)
      dfpt_data%drhop = (0.d0, 0.d0)
      dfpt_data%dvscfs = (0.d0, 0.d0)
      dfpt_data%dvscfp = (0.d0, 0.d0)
      !
      lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
      IF (lmetq0) THEN
         ALLOCATE(dfpt_data%def(npert))
         ALLOCATE(dfpt_data%dn0(npert))
         dfpt_data%def = (0.d0, 0.d0)
         dfpt_data%dn0 = (0.d0, 0.d0)
      ENDIF
      !
      IF (okvan) THEN
         ALLOCATE(dfpt_data%dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag, npert))
         dfpt_data%dbecsum = (0.d0, 0.d0)
      ELSE
         ! If okvan is false, dbecsum is not used.
         ! TODO: Do not allocate it (need to change the code that uses it).
         ALLOCATE(dfpt_data%dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag, npert))
      ENDIF
      !
   END SUBROUTINE allocate_dfpt_data
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   SUBROUTINE allocate_dfpt_ldos(ldos_data)
      !! Allocate arrays for LDOS calculation
      !
      USE fft_base,             ONLY : dfftp, dffts
      USE noncollin_module,     ONLY : nspin_mag
      USE uspp_param,           ONLY : nhm
      USE ions_base,            ONLY : nat
      !
      IMPLICIT NONE
      !
      TYPE(dfpt_ldos_type), INTENT(OUT) :: ldos_data
      !
      ! Initialize scalar
      ldos_data%dos_ef = 0.d0
      !
      ! Allocate grid quantities
      ALLOCATE(ldos_data%ldoss(dffts%nnr, nspin_mag))
      ALLOCATE(ldos_data%ldos(dfftp%nnr, nspin_mag))
      ldos_data%ldoss = (0.d0, 0.d0)
      ldos_data%ldos = (0.d0, 0.d0)
      !
      ! Always allocate becsum_dos for consistency
      ALLOCATE(ldos_data%becsum_dos((nhm * (nhm + 1))/2, nat, nspin_mag))
      ldos_data%becsum_dos = 0.d0
      !
   END SUBROUTINE allocate_dfpt_ldos
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   SUBROUTINE deallocate_dfpt_data(dfpt_data)
      !
      TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
      !
      IF (ALLOCATED(dfpt_data%drhos)) DEALLOCATE(dfpt_data%drhos)
      IF (ALLOCATED(dfpt_data%drhop)) DEALLOCATE(dfpt_data%drhop)
      IF (ALLOCATED(dfpt_data%dvscfp)) DEALLOCATE(dfpt_data%dvscfp)
      IF (ALLOCATED(dfpt_data%dvscfs)) DEALLOCATE(dfpt_data%dvscfs)
      IF (ALLOCATED(dfpt_data%dbecsum)) DEALLOCATE(dfpt_data%dbecsum)
      IF (ALLOCATED(dfpt_data%def)) DEALLOCATE(dfpt_data%def)
      IF (ALLOCATED(dfpt_data%dn0)) DEALLOCATE(dfpt_data%dn0)
      IF (ALLOCATED(dfpt_data%drhoc)) DEALLOCATE(dfpt_data%drhoc)
      IF (ALLOCATED(dfpt_data%drhop_pulay)) DEALLOCATE(dfpt_data%drhop_pulay)
      IF (ALLOCATED(dfpt_data%dbecsum_pulay)) DEALLOCATE(dfpt_data%dbecsum_pulay)
      !
   END SUBROUTINE deallocate_dfpt_data
   !---------------------------------------------------------------------------
   !
   !---------------------------------------------------------------------------
   SUBROUTINE deallocate_dfpt_ldos(ldos_data)
      !! Deallocate LDOS data structure
      !
      IMPLICIT NONE
      !
      TYPE(dfpt_ldos_type), INTENT(INOUT) :: ldos_data
      !
      IF (ALLOCATED(ldos_data%ldoss)) DEALLOCATE(ldos_data%ldoss)
      IF (ALLOCATED(ldos_data%ldos)) DEALLOCATE(ldos_data%ldos)
      IF (ALLOCATED(ldos_data%becsum_dos)) DEALLOCATE(ldos_data%becsum_dos)
      !
   END SUBROUTINE deallocate_dfpt_ldos
   !---------------------------------------------------------------------------
END MODULE dfpt_type
