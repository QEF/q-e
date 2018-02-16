#if defined(__CUDA)
PROGRAM test_mp_bcast_iv_gpu
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
    REAL :: rnd(datasize)
    !
    ! test variable
    INTEGER, DEVICE :: iv_d(datasize)
    INTEGER :: iv_h(datasize)
    INTEGER :: aux_h(datasize)
    
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    iv_h = mpime + 1
    iv_d = iv_h
    CALL mp_bcast(iv_d, root, world_comm)
    iv_h = iv_d
    !
    CALL test%assert_equal(ALL(iv_h .eq. 1) , .true. )
    !
    iv_h = mpime
    iv_d = iv_h
    CALL mp_bcast(iv_d, nproc-1, world_comm)
    iv_h = iv_d
    !
    CALL test%assert_equal(ALL(iv_h .eq. nproc-1) , .true. )
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
      iv_h = INT ( 10.0 * rnd )
      iv_d = iv_h
      CALL mp_bcast(iv_d, i , world_comm)
      CALL mp_bcast(iv_h, i , world_comm)
      aux_h = iv_d
      CALL test%assert_equal(SUM(iv_h) , SUM(aux_h) )
    END DO
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_iv_gpu
#else
PROGRAM test_mp_bcast_iv_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_iv_gpu
#endif
