#if defined(__CUDA)
PROGRAM test_mp_bcast_im_gpu
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
    ! test variable
    INTEGER, DEVICE :: im_d(10, 10)
    INTEGER :: im_h(10, 10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    im_h(:, :) = mpime
    im_d(:, :) = im_h(:, :)
    CALL mp_bcast(im_d, root, world_comm)
    im_h(:, :) = im_d(:, :)
    !
    CALL test%assert_equal(ALL(im_h .eq. 0) , .true. )
    !
    im_h(:, :) = mpime
    im_d(:, :) = im_h(:, :)
    CALL mp_bcast(im_d, nproc-1, world_comm)
    im_h(:, :) = im_d(:, :)
    !
    CALL test%assert_equal(ALL(im_h .eq. nproc-1) , .true. )
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_im_gpu
#else
PROGRAM test_mp_bcast_im_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_im_gpu
#endif
