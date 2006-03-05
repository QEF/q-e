!-----------------------------------------
! VARIABLES FOR TASKGROUPS
! C. Bekas, October 2005
!-----------------------------------------
!
!Variable description
!--------------------
!MAXGRP:	Maximum number of task-groups
!--------------------------------------------

MODULE GROUPS_MODULE

   USE kinds, ONLY: DP
   USE parameters, ONLY: MAXCPU, MAXGRP

   IMPLICIT NONE
   SAVE

   INTEGER, DIMENSION(:), ALLOCATABLE :: ALL_Z_STICKS
   INTEGER, DIMENSION(:), ALLOCATABLE  :: NOLIST, NPLIST, PGROUP
   INTEGER, DIMENSION(:), ALLOCATABLE :: tmp_nsw, tmp_npp, tmp_planes, tmp_revs, recvs, tmp_ismap
   INTEGER, DIMENSION(:), ALLOCATABLE :: ngw_vec !GLOBAL VECTOR OF ALL NGW VECTORS
   COMPLEX(DP), DIMENSION(:,:),  ALLOCATABLE :: tg_betae
   COMPLEX(DP), DIMENSION(:),    ALLOCATABLE :: tg_c2, tg_c3
   REAL(DP),  DIMENSION(:,:), ALLOCATABLE :: tg_rhos
   REAL(DP),  DIMENSION(:),   ALLOCATABLE :: tg_ggp
   INTEGER,  DIMENSION(:),   ALLOCATABLE  ::  nnrsx_vec !INCREASE THIS TO THE MAXIMUM NUMBER OF PROCS IF NEEDED
   REAL(DP), DIMENSION(:,:), ALLOCATABLE  ::  tmp_rhos, local_rhos
   INTEGER  :: recv_cnt(MAXGRP), recv_displ(MAXGRP), send_cnt(MAXGRP), send_displ(MAXGRP)
   INTEGER  :: SZ, RANK, CLOCK1, CLOCK2, CLOCK3, CLOCK4
   INTEGER  :: sticks_index, eig_offset, strd       
   REAL(DP) :: tm_tg, tm_rhoofr

CONTAINS

SUBROUTINE DEALLOCATE_GROUPS

   IMPLICIT NONE

   !DEALLOCATE GROUPS RELATED ARRAYS

   IF (ALLOCATED(ALL_Z_STICKS)) DEALLOCATE(ALL_Z_STICKS)
   IF (ALLOCATED(NOLIST))       DEALLOCATE(NOLIST)
   IF (ALLOCATED(NPLIST))       DEALLOCATE(NPLIST)
   IF (ALLOCATED(PGROUP))       DEALLOCATE(PGROUP)
   IF (ALLOCATED(tmp_nsw))      DEALLOCATE(tmp_nsw)
   IF (ALLOCATED(tmp_npp))      DEALLOCATE(tmp_npp)
   IF (ALLOCATED(tmp_planes))   DEALLOCATE(tmp_planes)
   IF (ALLOCATED(tmp_revs))     DEALLOCATE(tmp_revs)
   IF (ALLOCATED(recvs))        DEALLOCATE(recvs)
   IF (ALLOCATED(tmp_ismap))    DEALLOCATE(tmp_ismap)
   IF (ALLOCATED(ngw_vec))      DEALLOCATE(ngw_vec)
   IF (ALLOCATED(tg_betae))     DEALLOCATE(tg_betae)
   IF (ALLOCATED(tg_c2))        DEALLOCATE(tg_c2)
   IF (ALLOCATED(tg_c3))        DEALLOCATE(tg_c3)
   IF (ALLOCATED(tg_rhos))      DEALLOCATE(tg_rhos)
   IF (ALLOCATED(tg_ggp))       DEALLOCATE(tg_ggp)
   IF (ALLOCATED(nnrsx_vec))    DEALLOCATE(nnrsx_vec)
   IF (ALLOCATED(tmp_rhos))     DEALLOCATE(tmp_rhos)
   IF (ALLOCATED(local_rhos))   DEALLOCATE(local_rhos)

END SUBROUTINE DEALLOCATE_GROUPS

END MODULE GROUPS_MODULE
