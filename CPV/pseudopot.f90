!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Tue Nov  2 08:03:11 MET 1999
!  ----------------------------------------------

#include "f_defs.h"

MODULE pseudopotential

  ! ...   declare modules

  USE kinds,   ONLY: DP
  USE splines, ONLY: spline_data

  IMPLICIT NONE
  SAVE

  !  declare module-scope variables

  INTEGER :: nsanl   ! number of atoms of the non local species

  TYPE (spline_data), ALLOCATABLE ::  vps_sp(:)
  TYPE (spline_data), ALLOCATABLE :: dvps_sp(:)
  !
  TYPE (spline_data), ALLOCATABLE ::  wnl_sp(:,:)
  TYPE (spline_data), ALLOCATABLE :: wnla_sp(:,:)
  !
  TYPE (spline_data), ALLOCATABLE :: rhoc1_sp(:)
  TYPE (spline_data), ALLOCATABLE :: rhocp_sp(:)
  !
  REAL(DP), ALLOCATABLE :: xgtab(:)

  LOGICAL               :: tpstab = .TRUE.

!  ----------------------------------------------

CONTAINS

!  ----------------------------------------------

      SUBROUTINE deallocate_pseudopotential

          USE splines,      ONLY: kill_spline
          USE local_pseudo, ONLY: deallocate_local_pseudo
          USE uspp,         ONLY: dvan

          INTEGER :: i, j

          CALL deallocate_local_pseudo()
          !
          IF( ALLOCATED( dvan ) ) DEALLOCATE( dvan )
          IF( ALLOCATED( xgtab ) )  DEALLOCATE( xgtab )
          !
          IF( ALLOCATED( vps_sp  ) ) THEN
            DO i = 1, size(vps_sp)
              CALL kill_spline(vps_sp(i),'a')
            END DO
            DEALLOCATE(vps_sp)
          END IF
          !
          IF( ALLOCATED(dvps_sp) ) THEN
            DO i = 1, size(dvps_sp)
              CALL kill_spline(dvps_sp(i),'a')
            END DO
            DEALLOCATE(dvps_sp)
          END IF
          !
          IF( ALLOCATED(rhoc1_sp) ) THEN
            DO i = 1, size(rhoc1_sp)
              CALL kill_spline(rhoc1_sp(i),'a')
            END DO
            DEALLOCATE(rhoc1_sp)
          END IF
          !
          IF( ALLOCATED(rhocp_sp) ) THEN
            DO i = 1, size(rhocp_sp)
              CALL kill_spline(rhocp_sp(i),'a')
            END DO
            DEALLOCATE(rhocp_sp)
          END IF
          !
          IF( ALLOCATED(wnl_sp) ) THEN
            DO i = 1, size(wnl_sp,2)
              DO j = 1, size(wnl_sp,1)
                CALL kill_spline(wnl_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnl_sp)
          END IF
          !
          IF( ALLOCATED(wnla_sp) ) THEN
            DO i = 1, size(wnla_sp,2)
              DO j = 1, size(wnla_sp,1)
                CALL kill_spline(wnla_sp(j,i),'a')
              END DO
            END DO
            DEALLOCATE(wnla_sp)
          END IF

          RETURN
        END SUBROUTINE deallocate_pseudopotential


END MODULE pseudopotential

!  ----------------------------------------------
