!
! Copyright (C) 2003-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_opt_routines
   !---------------------------------------------------------------------------
   !
   ! ... This module contains all subroutines and functions needed for
   ! ... the optimisation of the reaction path (NEB and SMD calculations)
   !
   ! ... Written by Carlo Sbraccia ( 2003-2006 )
   !
   ! ... MDIIS algorithm is implemented by Satomichi Nishihara ( 2016 )
   !
   USE kinds,          ONLY : DP
   USE constants,      ONLY : eps8, eps16, e2, rytoev, fpi
   USE ring_variables, ONLY : pos
   USE io_global,      ONLY : meta_ionode, meta_ionode_id
   USE mp,             ONLY : mp_bcast
   USE mp_world,       ONLY : world_comm
   USE fcp_variables,  ONLY : fcp_mu, fcp_relax, fcp_relax_step, &
                              fcp_mdiis_size, fcp_mdiis_step, &
                              fcp_tot_charge_first, fcp_tot_charge_last
   USE mdiis,          ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis
   USE pimd_variables, ONLY : nbeadMD
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   REAL(DP), ALLOCATABLE :: fcp_neb_nelec(:)
   REAL(DP), ALLOCATABLE :: fcp_neb_ef(:)
   !
   ! ... variables for line-minimisation
   REAL(DP), ALLOCATABLE :: force0(:)
   REAL(DP), ALLOCATABLE :: nelec0(:)
   LOGICAL,  ALLOCATABLE :: firstcall(:)
   !
   ! ... variables for MDIIS
   LOGICAL          :: init_mdiis
   TYPE(mdiis_type) :: mdiist
   !
   PUBLIC :: fcp_neb_nelec, fcp_neb_ef, &
      fcp_opt_allocation, fcp_opt_deallocation
   !
CONTAINS
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_opt_allocation()
      !----------------------------------------------------------------------
      !
      USE ions_base,      ONLY : nat, ityp, zv
      USE klist,          ONLY : tot_charge, nelec
      USE io_files,       ONLY : prefix, tmp_dir
      USE ring_variables, ONLY : restart
      USE ener,           ONLY : ef
      !
      IMPLICIT NONE
      !
      REAL(DP) :: ionic_charge, nelec_, first, last
      INTEGER  :: n, i, ierr
      CHARACTER (LEN=256)   :: tmp_dir_saved
      !
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ALLOCATE( fcp_neb_nelec( nbeadMD ) )
      ALLOCATE( fcp_neb_ef   ( nbeadMD ) )
      !
      IF ( TRIM(fcp_relax) == 'lm' ) THEN
         !
         ALLOCATE( force0    ( nbeadMD ) )
         ALLOCATE( nelec0    ( nbeadMD ) )
         ALLOCATE( firstcall ( nbeadMD ) )
         !
         force0    (:) = 0.0_DP
         nelec0    (:) = 0.0_DP
         firstcall (:) = .TRUE.
         !
      ELSE IF ( TRIM(fcp_relax) == 'mdiis' ) THEN
         !
         init_mdiis = .TRUE.
         CALL allocate_mdiis(mdiist, fcp_mdiis_size, nbeadMD, fcp_mdiis_step, 1)
         !
      END IF
      !
      IF ( restart ) THEN
         !
         tmp_dir_saved = tmp_dir
         !
         DO i = 1, nbeadMD
            !
            tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                 TRIM( int_to_char( i ) ) // "/"
            !
            CALL errore('fcp_opt_routines','XSD implementation pending',1)
            !
            fcp_neb_nelec(i) = nelec
            fcp_neb_ef   (i) = ef * e2 ! factor e2: hartree -> Ry.
            !
         END DO
         !
         tmp_dir = tmp_dir_saved
         !
      ELSE
         !
         ionic_charge = SUM(zv(ityp(1:nat)))
         nelec_ = ionic_charge - tot_charge
         !
         n = nbeadMD
         first = fcp_tot_charge_first
         last  = fcp_tot_charge_last
         DO i = 1, n
            fcp_neb_nelec(i) = nelec_ - (first * (n - i) + last * (i - 1) ) / (n - 1)
         END DO
         !
         fcp_neb_ef(:) = 0.0_DP
         !
      END IF
      !
   END SUBROUTINE fcp_opt_allocation
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_opt_deallocation()
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      IF ( ALLOCATED( fcp_neb_nelec ) ) DEALLOCATE( fcp_neb_nelec )
      IF ( ALLOCATED( fcp_neb_ef    ) ) DEALLOCATE( fcp_neb_ef    )
      IF ( ALLOCATED( force0        ) ) DEALLOCATE( force0        )
      IF ( ALLOCATED( nelec0        ) ) DEALLOCATE( nelec0        )
      IF ( ALLOCATED( firstcall     ) ) DEALLOCATE( firstcall     )
      !
      IF ( init_mdiis ) THEN
         !
         CALL deallocate_mdiis(mdiist)
         !
      END IF
      !
   END SUBROUTINE fcp_opt_deallocation
   !
   !
END MODULE fcp_opt_routines
