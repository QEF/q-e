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
   USE kinds,          ONLY : DP
   USE constants,      ONLY : eps8, eps16, e2, rytoev, fpi
   USE path_variables, ONLY : ds, pos, grad
   USE io_global,      ONLY : meta_ionode, meta_ionode_id
   USE mp,             ONLY : mp_bcast
   USE mp_world,       ONLY : world_comm
   USE fcp_variables,  ONLY : fcp_mu, fcp_relax_step, fcp_tot_charge_first, &
                              fcp_tot_charge_last
   USE path_variables, ONLY : num_of_images
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   REAL(DP), ALLOCATABLE :: fcp_neb_nelec(:)
   REAL(DP), ALLOCATABLE :: fcp_neb_ef(:)
   !
   REAL(DP), ALLOCATABLE :: force0(:)
   REAL(DP), ALLOCATABLE :: nelec0(:)
   LOGICAL,  ALLOCATABLE :: firstcall(:)
   !
   PUBLIC :: fcp_neb_nelec, fcp_neb_ef, fcp_line_minimisation, &
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
      USE pw_restart,     ONLY : pw_readfile
      USE io_files,       ONLY : prefix, tmp_dir
      USE path_variables, ONLY : restart
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
      ALLOCATE( fcp_neb_nelec( num_of_images ) )
      ALLOCATE( fcp_neb_ef   ( num_of_images ) )
      ALLOCATE( force0       ( num_of_images ) )
      ALLOCATE( nelec0       ( num_of_images ) )
      ALLOCATE( firstcall    ( num_of_images ) )
      !

      IF ( restart ) THEN
         !
         tmp_dir_saved = tmp_dir
         !
         DO i = 1, num_of_images
            !
            tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                 TRIM( int_to_char( i ) ) // "/"
            !
            CALL pw_readfile('fcpopt', ierr)
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
         n = num_of_images
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
      force0    (:) = 0.0_DP
      nelec0    (:) = 0.0_DP
      firstcall (:) = .TRUE.
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
   END SUBROUTINE fcp_opt_deallocation
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_line_minimisation()
      !----------------------------------------------------------------------
      !
      USE ions_base, ONLY : nat, ityp, zv
      USE cell_base, ONLY : at, alat
      !
      USE path_variables,       ONLY : frozen
      !
      IMPLICIT NONE
      !
      INTEGER  :: image
      REAL(DP) :: force, ef, nelec, n_tmp, max_tot_charge, capacitance, ionic_charge
      !
      IF ( meta_ionode ) THEN
         !
         DO image = 1, num_of_images
            !
            IF ( frozen(image) ) CYCLE
            !
            ef    = fcp_neb_ef(image)
            nelec = fcp_neb_nelec(image)
            !
            force = fcp_mu - ef
            !
            ! ... assumption: capacitance with vacuum gives the upper bound of 
            ! ... tot_charge difference.
            !
            capacitance = (at(1,1) * at(2,2) - at(2,1) * at(1,2)) * alat**2 &
                          / (alat * at(3,3) / 2._DP ) / fpi
            max_tot_charge = abs( capacitance * force / e2 )
            IF ( firstcall(image) .OR. ABS( force0(image) - force ) < 1.0D-20 ) THEN
               firstcall(image) = .FALSE.
               nelec0(image) = nelec
               force0(image) = force
               nelec = nelec + fcp_relax_step * force
            ELSE
               n_tmp = nelec
               nelec = (nelec * force0(image) - nelec0(image) * force ) &
                    / ( force0(image) - force )
               nelec0(image) = n_tmp
               force0(image) = force
            END IF
            ionic_charge = SUM(zv(ityp(1:nat)))
            !write( *,'(/,5X,"Upper bound for tot_charge:",F12.6)') &
            !       max_tot_charge
            !write( *,'(5X,"Original:",F12.6,"Expected:",F12.6)') &
            !       ionic_charge - nelec0(image), ionic_charge - nelec
            if( nelec-nelec0(image) < -max_tot_charge ) &
                nelec= nelec0(image) - max_tot_charge
            if( nelec-nelec0(image) >  max_tot_charge ) &
                nelec= nelec0(image) + max_tot_charge
            !write( *,'(5X,"Next tot_charge:",F12.6)') ionic_charge - nelec
            !
            fcp_neb_nelec(image) = nelec
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( fcp_neb_nelec, meta_ionode_id, world_comm )
      !
      RETURN
      !
   END SUBROUTINE fcp_line_minimisation
   !
END MODULE fcp_opt_routines
