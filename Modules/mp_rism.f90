!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE mp_rism
  !--------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_size, mp_rank, mp_barrier, mp_comm_split, mp_sum, mp_comm_free
  USE parallel_include
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... Site groups, for parallelization over sites
  TYPE mp_rism_site
    ! ... MPI's setting
    INTEGER :: nsitg      = 1  ! number of site groups
    INTEGER :: nproc_sitg = 1  ! number of processies within a site group
    INTEGER :: me_sitg    = 0  ! index of the process within a site group
    INTEGER :: root_sitg  = 0  ! index of the root process within a site group
    INTEGER :: my_sitg_id = 0  ! index of my site group
#if defined(__MPI)
    INTEGER :: inter_sitg_comm = MPI_COMM_NULL  ! inter site group communicator
    INTEGER :: intra_sitg_comm = MPI_COMM_NULL  ! intra site group communicator
#else
    INTEGER :: inter_sitg_comm = 0
    INTEGER :: intra_sitg_comm = 0
#endif
    LOGICAL :: inter_sitg_keep = .FALSE.        ! keep inter_sitg_comm by oneself or not
    LOGICAL :: intra_sitg_keep = .FALSE.        ! keep intra_sitg_comm by oneself or not
    ! ... site splitting
    INTEGER :: nsite       = 0  ! total number of sites
    INTEGER :: isite_start = 0  ! starting site index
    INTEGER :: isite_end   = 0  ! ending site index
  END TYPE mp_rism_site
  !
  ! ... A task group, for parallelization over R- and G-spaces
  TYPE mp_rism_task
    ! ... MPI's setting
    INTEGER :: nproc_task = 1  ! number of processies within a task group
    INTEGER :: me_task    = 0  ! index of the process within a task group
    INTEGER :: root_task  = 0  ! index of the root process within a task group
#if defined(__MPI)
    INTEGER :: itask_comm = MPI_COMM_NULL  ! task group communicator
#else
    INTEGER :: itask_comm = 0
#endif
    LOGICAL :: itask_keep = .FALSE.        ! keep itask_comm by oneself or not
    ! ... vector splitting
    INTEGER          :: nvec       = 0  ! total number of vector
    INTEGER          :: ivec_start = 0  ! starting vector index
    INTEGER          :: ivec_end   = 0  ! ending vector index
    INTEGER, POINTER :: ilen_vecs(:)    ! lengths of vectors in all processies
    INTEGER, POINTER :: idis_vecs(:)    ! displacement of vectors in all processies
  END TYPE mp_rism_task
  !
  ! ... public components
  PUBLIC :: mp_rism_site
  PUBLIC :: mp_rism_task
  PUBLIC :: mp_start_rism_task_and_site
  PUBLIC :: mp_start_rism_task_on_site
  PUBLIC :: mp_end_rism
  PUBLIC :: mp_set_index_rism_site
  PUBLIC :: mp_set_index_rism_task
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE mp_start_rism_task_and_site(mp_risms, mp_rismt, parent_comm)
    !--------------------------------------------------------------------------
    !
    ! ... Copy "parent_comm" to task and site groups
    ! ... Requires: parent_comm, typically world_comm = group of all processors
    !
    IMPLICIT NONE
    !
    TYPE(mp_rism_site), INTENT(INOUT) :: mp_risms
    TYPE(mp_rism_task), INTENT(INOUT) :: mp_rismt
    INTEGER,            INTENT(IN)    :: parent_comm
    !
    INTEGER :: parent_nproc
    INTEGER :: parent_mype
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size(parent_comm)
    parent_mype  = mp_rank(parent_comm)
    !
    ! ... create mp_rism_site
    mp_risms%nsitg           = parent_nproc
    mp_risms%nproc_sitg      = 1
    mp_risms%me_sitg         = 0
    mp_risms%root_sitg       = 0
    mp_risms%my_sitg_id      = parent_mype
    mp_risms%inter_sitg_comm = parent_comm
    mp_risms%intra_sitg_comm = MPI_COMM_NULL
    mp_risms%inter_sitg_keep = .FALSE.
    mp_risms%intra_sitg_keep = .FALSE.
    !
    ! ... create mp_rism_task
    mp_rismt%nproc_task      = parent_nproc
    mp_rismt%me_task         = parent_mype
    mp_rismt%root_task       = 0
    mp_rismt%itask_comm      = parent_comm
    mp_rismt%itask_keep      = .FALSE.
    !
#else
    !
    ! ... create mp_rism_site
    mp_risms%nsitg           = 1
    mp_risms%nproc_sitg      = 1
    mp_risms%me_sitg         = 0
    mp_risms%root_sitg       = 0
    mp_risms%my_sitg_id      = 0
    mp_risms%inter_sitg_comm = 0 !MPI_COMM_NULL
    mp_risms%intra_sitg_comm = 0 !MPI_COMM_NULL
    mp_risms%inter_sitg_keep = .FALSE.
    mp_risms%intra_sitg_keep = .FALSE.
    !
    ! ... create mp_rism_task
    mp_rismt%nproc_task      = 1
    mp_rismt%me_task         = 0
    mp_rismt%root_task       = 0
    mp_rismt%itask_comm      = 0 !MPI_COMM_NULL
    mp_rismt%itask_keep      = .FALSE.
    !
#endif
  END SUBROUTINE mp_start_rism_task_and_site
  !
  !--------------------------------------------------------------------------
  SUBROUTINE mp_start_rism_task_on_site(mp_risms, mp_rismt, itask_comm, parent_comm)
    !--------------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into site groups, which include a task group
    ! ... Requires: itask_comm,  task group (read from command line)
    ! ...           parent_comm, typically world_comm = group of all processors
    !
    IMPLICIT NONE
    !
    TYPE(mp_rism_site), INTENT(INOUT) :: mp_risms
    TYPE(mp_rism_task), INTENT(INOUT) :: mp_rismt
    INTEGER,            INTENT(IN)    :: itask_comm
    INTEGER,            INTENT(IN)    :: parent_comm
    !
    INTEGER :: itask_nproc
    INTEGER :: itask_mype
    INTEGER :: parent_nproc
    INTEGER :: parent_mype
    !
#if defined (__MPI)
    !
    itask_nproc  = mp_size(itask_comm)
    itask_mype   = mp_rank(itask_comm)
    parent_nproc = mp_size(parent_comm)
    parent_mype  = mp_rank(parent_comm)
    !
    ! ... create mp_rism_site
    IF (itask_nproc < 1 .OR. itask_nproc > parent_nproc) THEN
      CALL errore('mp_start_rism_task_on_site', 'invalid number of tasks, out of range', 1)
    END IF
    !
    IF (MOD(parent_nproc, itask_nproc) /= 0) THEN
      CALL errore('mp_start_rism_task_on_site', &
                & 'invalid number of tasks, parent_nproc /= nproc_task * nsite', 1)
    END IF
    !
    mp_risms%nsitg           = parent_nproc / itask_nproc
    mp_risms%nproc_sitg      = itask_nproc
    mp_risms%me_sitg         = MOD(parent_mype, itask_nproc)
    mp_risms%root_sitg       = 0
    mp_risms%my_sitg_id      = parent_mype / itask_nproc
    mp_risms%intra_sitg_comm = itask_comm
    mp_risms%intra_sitg_keep = .FALSE.
    !
    CALL mp_barrier(parent_comm)
    CALL mp_comm_split(parent_comm, mp_risms%me_sitg, parent_mype, mp_risms%inter_sitg_comm)
    mp_risms%inter_sitg_keep = .TRUE.
    !
    ! ... create mp_rism_task
    mp_rismt%nproc_task      = itask_nproc
    mp_rismt%me_task         = itask_mype
    mp_rismt%root_task       = 0
    mp_rismt%itask_comm      = itask_comm
    mp_rismt%itask_keep      = .FALSE.
    !
#else
    !
    ! ... create mp_rism_site
    mp_risms%nsitg           = 1
    mp_risms%nproc_sitg      = 1
    mp_risms%me_sitg         = 0
    mp_risms%root_sitg       = 0
    mp_risms%my_sitg_id      = 0
    mp_risms%inter_sitg_comm = 0 !MPI_COMM_NULL
    mp_risms%intra_sitg_comm = 0 !MPI_COMM_NULL
    mp_risms%inter_sitg_keep = .FALSE.
    mp_risms%intra_sitg_keep = .FALSE.
    !
    ! ... create mp_rism_task
    mp_rismt%nproc_task      = 1
    mp_rismt%me_task         = 0
    mp_rismt%root_task       = 0
    mp_rismt%itask_comm      = 0 !MPI_COMM_NULL
    mp_rismt%itask_keep      = .FALSE.
    !
#endif
  END SUBROUTINE mp_start_rism_task_on_site
  !
  !--------------------------------------------------------------------------
  SUBROUTINE mp_end_rism(mp_risms, mp_rismt)
    !--------------------------------------------------------------------------
    !
    ! ... Release communicator, and deallocate memories.
    !
    IMPLICIT NONE
    !
    TYPE(mp_rism_site), INTENT(INOUT) :: mp_risms
    TYPE(mp_rism_task), INTENT(INOUT) :: mp_rismt
    !
#if defined (__MPI)
    !
    ! ... delte mp_rism_site
    IF (mp_risms%inter_sitg_keep .AND. mp_risms%inter_sitg_comm /= MPI_COMM_NULL) THEN
      CALL mp_comm_free(mp_risms%inter_sitg_comm)
    END IF
    !
    IF (mp_risms%intra_sitg_keep .AND. mp_risms%intra_sitg_comm /= MPI_COMM_NULL) THEN
      CALL mp_comm_free(mp_risms%intra_sitg_comm)
    END IF
    !
    ! ... delte mp_rism_task
    IF (mp_rismt%itask_keep .AND. mp_rismt%itask_comm /= MPI_COMM_NULL) THEN
      CALL mp_comm_free(mp_rismt%itask_comm)
    END IF
    !
#endif
    !
    ! ... deallocate arrays
    IF (ASSOCIATED(mp_rismt%ilen_vecs)) THEN
      DEALLOCATE(mp_rismt%ilen_vecs)
    END IF
    !
    IF (ASSOCIATED(mp_rismt%idis_vecs)) THEN
      DEALLOCATE(mp_rismt%idis_vecs)
    END IF
    !
  END SUBROUTINE mp_end_rism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE mp_set_index_rism_site(mp_risms, nsite)
    !--------------------------------------------------------------------------
    !
    ! ... create and keep indexes of site-parallel.
    !
    IMPLICIT NONE
    TYPE(mp_rism_site), INTENT(INOUT) :: mp_risms
    INTEGER,            INTENT(IN)    :: nsite
    !
    INTEGER :: npe
    INTEGER :: myrank
    INTEGER :: rest
    INTEGER :: k
    !
    mp_risms%nsite = nsite
    !
    myrank = mp_risms%my_sitg_id
    npe    = mp_risms%nsitg
    rest   = MOD(nsite, npe)
    k      = INT(nsite / npe)
    !
    ! ... set isite_start, isite_end
    IF (k >= 0) THEN
      IF (rest > myrank) THEN
        mp_risms%isite_start = (myrank    ) * k + (myrank + 1)
        mp_risms%isite_end  =  (myrank + 1) * k + (myrank + 1)
      ELSE
        mp_risms%isite_start = (myrank    ) * k + rest + 1
        mp_risms%isite_end   = (myrank + 1) * k + rest
      END IF
    ELSE
      CALL errore(' mp_set_index_rism_site ', ' too small nsite ', 1)
    END IF
    !
  END SUBROUTINE mp_set_index_rism_site
  !
  !--------------------------------------------------------------------------
  SUBROUTINE mp_set_index_rism_task(mp_rismt, nvec)
    !--------------------------------------------------------------------------
    !
    ! ... create and keep indexes of task-parallel.
    !
    IMPLICIT NONE
    TYPE(mp_rism_task), INTENT(INOUT) :: mp_rismt
    INTEGER,            INTENT(IN)    :: nvec
    !
    INTEGER :: npe
    INTEGER :: myrank
    INTEGER :: rest
    INTEGER :: k
    !
    mp_rismt%nvec = nvec
    !
    myrank = mp_rismt%me_task
    npe    = mp_rismt%nproc_task
    rest   = MOD(nvec, npe)
    k      = INT(nvec / npe)
    !
    IF (k < 1) THEN
      CALL errore('mp_set_index_rism_task', 'too much processies npe > nvec', 1)
    END IF
    !
    ! ... set ivec_start, ivec_end
    IF (k >= 1) THEN
      IF (rest > myrank) THEN
        mp_rismt%ivec_start = (myrank    ) * k + (myrank + 1)
        mp_rismt%ivec_end  =  (myrank + 1) * k + (myrank + 1)
      ELSE
        mp_rismt%ivec_start = (myrank    ) * k + rest + 1
        mp_rismt%ivec_end   = (myrank + 1) * k + rest
      END IF
    ELSE
      CALL errore(' mp_set_index_rism_task ', ' too small nvec ', 1)
    END IF
    !
    ! ... set ilen_vecs
    ALLOCATE(mp_rismt%ilen_vecs(npe))
    mp_rismt%ilen_vecs(:) = 0
    mp_rismt%ilen_vecs(myrank + 1) = mp_rismt%ivec_end - mp_rismt%ivec_start + 1
    CALL mp_sum(mp_rismt%ilen_vecs, mp_rismt%itask_comm)
    !
    ! ... set idis_vecs
    ALLOCATE(mp_rismt%idis_vecs(npe))
    mp_rismt%idis_vecs(:) = 0
    mp_rismt%idis_vecs(myrank + 1) = mp_rismt%ivec_start - 1
    CALL mp_sum(mp_rismt%idis_vecs, mp_rismt%itask_comm)
    !
  END SUBROUTINE mp_set_index_rism_task
  !
END MODULE mp_rism

