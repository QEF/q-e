#if defined(__CUDA)
program test_diaghg_gpu

    USE laxlib_parallel_include
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, nproc, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    !
    me_bgrp = mpime; root_bgrp=root; intra_bgrp_comm=world_comm
    !
    CALL real_1(test)
    !
    CALL complex_1(test)
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
  CONTAINS
  !
  SUBROUTINE real_1(test)
    USE LAXlib
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    ! variables on device
    real(DP), device :: h_d(2,2)
    real(DP), device :: s_d(2,2)
    real(DP), device :: e_d(2)
    real(DP), device :: v_d(2,2)
    ! variables on host
    real(DP) :: h_h(2,2)
    real(DP) :: s_h(2,2)
    real(DP) :: e_h(2)
    real(DP) :: v_h(2,2)
    
    h_h = 0.d0
    h_h(1,1) = 1.d0
    h_h(2,2) = 1.d0
    h_d = h_h
    
    s_h = 0.d0
    s_h(1,1) = 1.d0
    s_h(2,2) = 1.d0
    s_d = s_h
    !
    v_d = 0.d0
    e_d = 0.d0
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    v_h = v_d
    e_h = e_d
    h_h = h_d
    s_h = s_d
    !
    CALL test%assert_close( e_h, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v_h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h_h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(s_h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    !
    v_d = 0.d0
    e_d = 0.d0
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
    v_h = v_d
    e_h = e_d
    h_h = h_d
    s_h = s_d
    !
    CALL test%assert_close( e_h, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v_h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h_h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(s_h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    !
  END SUBROUTINE real_1
  !
  SUBROUTINE complex_1(test)
    USE LAXlib
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    ! variables on device
    complex(DP), device :: h_d(2,2)
    complex(DP), device :: s_d(2,2)
    complex(DP), device :: v_d(2,2)
    real(DP),    device :: e_d(2)
    ! variables on host
    complex(DP) :: h_h(2,2)
    complex(DP) :: s_h(2,2)
    complex(DP) :: v_h(2,2)
    real(DP)    :: e_h(2)


    complex(DP) :: s_save(2,2)
    complex(DP) :: h_save(2,2)


    !
    h_h = 0.d0
    h_h(1,1) = (1.d0,  0.d0)
    h_h(1,2) = (0.d0, -2.d0)
    h_h(2,1) = ( 0.d0, 2.d0)
    h_h(2,2) = ( 5.d0, 0.d0)
    s_h = 0.d0
    s_h(1,1) = (1.d0, 0.d0)
    s_h(2,2) = (1.d0, 0.d0)
    !
    ! save for later comparison
    h_save = h_h
    s_save = s_h
    !
    ! Update device
    h_d = h_h
    s_d = s_h
    !
    v_h = (0.d0, 0.d0)
    e_h = 0.d0
    v_d = v_h; e_d = e_h
    !
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    v_h = v_d
    e_h = e_d
    h_h = h_d
    s_h = s_d
    !
    ! 0.1715728752538099, 5.82842712474619 
    CALL test%assert_close( e_h, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( v_h(:,1), [( 0.d0, -0.9238795325112867d0), (-0.3826834323650898d0, 0.d0)] )
    CALL test%assert_close( v_h(:,2), [( 0.d0, -0.3826834323650898d0), ( 0.9238795325112867d0, 0.d0)] )
    CALL test%assert_close( RESHAPE(h_h, [4]), RESHAPE(h_save, [4]))
    CALL test%assert_close( RESHAPE(s_h, [4]), RESHAPE(s_save, [4]))
    !
    v_h = (0.d0, 0.d0)
    e_h = 0.d0
    !
    ! Update device
    h_d = h_h
    s_d = s_h
    v_d = v_h; e_d = e_h
    CALL diaghg(  2, 2, h_d, s_d, 2, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
    v_h = v_d
    e_h = e_d
    h_h = h_d
    s_h = s_d
    !
    CALL test%assert_close( e_h, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( v_h(:,1), [( 0.d0, -0.9238795325112867d0), (-0.3826834323650898d0, 0.d0)] )
    CALL test%assert_close( v_h(:,2), [( 0.d0, -0.3826834323650898d0), ( 0.9238795325112867d0, 0.d0)] )
    CALL test%assert_close( RESHAPE(h_h, [4]), RESHAPE(h_save, [4]))
    CALL test%assert_close( RESHAPE(s_h, [4]), RESHAPE(s_save, [4]))
    !
  END SUBROUTINE complex_1
end program test_diaghg_gpu
#else
program test_diaghg_gpu
end program test_diaghg_gpu
#endif
