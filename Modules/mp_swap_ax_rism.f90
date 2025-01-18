!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE mp_swap_ax_rism(mp_site, mp_task, lsite, rsite, ltask, rtask, isign)
  !---------------------------------------------------------------------------
  !
  ! ... swap parallel axies of RISM's data (for real numbers)
  ! ... if isign > 0, swap rtask -> rsite.
  ! ... if isign < 0, swap rsite -> rtask.
  !
  USE kinds,   ONLY : DP
  USE mp,      ONLY : mp_max
  USE mp_rism, ONLY : mp_rism_site, mp_rism_task
  USE parallel_include
  !
  IMPLICIT NONE
  !
  TYPE(mp_rism_site), INTENT(IN)    :: mp_site
  TYPE(mp_rism_task), INTENT(IN)    :: mp_task
  INTEGER,            INTENT(IN)    :: lsite
  REAL(DP),           INTENT(INOUT) :: rsite(lsite,1:*)
  INTEGER,            INTENT(IN)    :: ltask
  REAL(DP),           INTENT(INOUT) :: rtask(ltask,1:*)
  INTEGER,            INTENT(IN)    :: isign
  !
  REAL(DP), ALLOCATABLE :: rwork(:)
  !
  ALLOCATE(rwork(mp_task%nvec))
  !
  IF (isign > 0) THEN
    CALL rtask_to_rsite()
  ELSE IF (isign < 0) THEN
    CALL rsite_to_rtask()
  ELSE !IF (isign == 0) THEN
    ! NOP
  END IF
  !
  DEALLOCATE(rwork)
  !
CONTAINS
  !
#if defined (__MPI)
  SUBROUTINE rtask_to_rsite()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: iroot
    INTEGER :: ierr
    !
    DO isite = 1, mp_site%nsite
      ! ... get root in a task group
      IF (mp_site%isite_start <= isite .AND. isite <= mp_site%isite_end) THEN
        iroot = mp_task%me_task + 1
      ELSE
        iroot = 0
      END IF
      !
      CALL mp_max(iroot, mp_task%itask_comm)
      iroot = iroot - 1
      !
      ! ... gather data
      IF (iroot >= 0) THEN
        CALL MPI_GATHERV(rtask(1, isite), mp_task%ilen_vecs(mp_task%me_task + 1), &
           & MPI_DOUBLE_PRECISION, rwork, mp_task%ilen_vecs, mp_task%idis_vecs, &
           & MPI_DOUBLE_PRECISION, iroot, mp_task%itask_comm, ierr)
        !
        IF (ierr /= MPI_SUCCESS) THEN
          CALL errore('mp_swap_ax_rism', 'error at MPI_GATHERV', 1)
        END IF
        !
      ELSE
        rwork(1:mp_task%nvec) = 0.0_DP
      END IF
      !
      ! ... copy gathered data to rsite
      IF (mp_site%isite_start <= isite .AND. isite <= mp_site%isite_end) THEN
        rsite(1:mp_task%nvec, isite - mp_site%isite_start + 1) = rwork(1:mp_task%nvec)
      END IF
      !
    END DO
    !
  END SUBROUTINE rtask_to_rsite
  !
  SUBROUTINE rsite_to_rtask()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: iroot
    INTEGER :: ierr
    !
    DO isite = 1, mp_site%nsite
      ! ... get root in a task group
      IF (mp_site%isite_start <= isite .AND. isite <= mp_site%isite_end) THEN
        iroot = mp_task%me_task + 1
        rwork(1:mp_task%nvec) = rsite(1:mp_task%nvec, isite - mp_site%isite_start + 1)
      ELSE
        iroot = 0
      END IF
      !
      CALL mp_max(iroot, mp_task%itask_comm)
      iroot = iroot - 1
      !
      ! ... scatter data
      IF (iroot >= 0) THEN
        CALL MPI_SCATTERV(rwork, mp_task%ilen_vecs, mp_task%idis_vecs, &
           & MPI_DOUBLE_PRECISION, rtask(1, isite), mp_task%ilen_vecs(mp_task%me_task + 1), &
           & MPI_DOUBLE_PRECISION, iroot, mp_task%itask_comm, ierr)
        !
        IF (ierr /= MPI_SUCCESS) THEN
          CALL errore('mp_swap_ax_rism', 'error at MPI_SCATTERV', 1)
        END IF
        !
      ELSE
        rtask(1:mp_task%nvec, isite) = 0.0_DP
      END IF
      !
    END DO
    !
  END SUBROUTINE rsite_to_rtask
  !
#else
  SUBROUTINE rtask_to_rsite()
    IMPLICIT NONE
    !
    rsite(1:mp_task%nvec, 1:mp_site%nsite) = rtask(1:mp_task%nvec, 1:mp_site%nsite)
  END SUBROUTINE rtask_to_rsite
  !
  SUBROUTINE rsite_to_rtask()
    IMPLICIT NONE
    !
    rtask(1:mp_task%nvec, 1:mp_site%nsite) = rsite(1:mp_task%nvec, 1:mp_site%nsite)
  END SUBROUTINE rsite_to_rtask
#endif
  !
END SUBROUTINE mp_swap_ax_rism
