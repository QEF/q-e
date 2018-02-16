#if defined(__CUDA)
PROGRAM test_mp_bcast_cm_gpu
!
! Simple program to check the functionalities of test_mp_bcast_i1.
!
    USE cudafor
#if defined(__MPI)
    USE MPI
#endif
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    INTEGER, PARAMETER :: datasize = 10
    !
    ! Stuff for Ranum data
    integer, allocatable :: seed(:)
    integer :: i, n
    REAL :: rnd(datasize,datasize)
    !
    ! test variable
    COMPLEX(8), DEVICE :: cm_d(datasize,datasize)
    COMPLEX(8) :: cm_h(datasize,datasize)
    COMPLEX(8) :: aux_h(datasize,datasize)
    
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    cm_h = mpime + 1
    cm_d = cm_h
    CALL mp_bcast(cm_d, root, world_comm)
    cm_h = cm_d
    !
    CALL test%assert_equal(ALL(cm_h .eq. 1) , .true. )
    !
    cm_h = mpime
    cm_d = cm_h
    CALL mp_bcast(cm_d, nproc-1, world_comm)
    cm_h = cm_d
    !
    CALL test%assert_equal(ALL(cm_h .eq. nproc-1) , .true. )
    !
    ! Test against CPU implementation
    CALL random_seed(size = n)
    ALLOCATE(seed(n))
    CALL random_seed(get=seed)
    WRITE (*, *) "Random seed: ", seed
    DEALLOCATE(seed)
    !
    DO i = 0, nproc-1
      CALL RANDOM_NUMBER(rnd)
      cm_h = DCMPLX ( 10.0 * rnd )
      cm_d = cm_h
      CALL mp_bcast(cm_d, i , world_comm)
      CALL mp_bcast(cm_h, i , world_comm)
      aux_h = cm_d
      CALL test%assert_equal(SUM(cm_h) , SUM(aux_h) )
    END DO
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_cm_gpu
#else
PROGRAM test_mp_bcast_cm_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_cm_gpu
#endif
