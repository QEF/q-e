program test_diaghg
#if defined(__MPI)
    USE MPI
#endif
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, nproc, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    IMPLICIT NONE
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
    USE la_param, ONLY : DP
    implicit none
    !
    TYPE(tester_t) :: test
    ! real variables
    real(DP) :: h(2,2)
    real(DP) :: s(2,2)
    real(DP) :: e(2)
    real(DP) :: v(2,2)
    
    h = 0.d0
    h(1,1) = 1.d0
    h(2,2) = 1.d0
    s = 0.d0
    s(1,1) = 1.d0
    s(2,2) = 1.d0
    !
    v = 0.d0
    e = 0.d0
    CALL diaghg(  2, 2, h, s, 2, e, v, .false. )
    !
    CALL test%assert_close( e, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(s, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    !
    v = 0.d0
    e = 0.d0
    CALL diaghg(  2, 2, h, s, 2, e, v, .true. )
    !
    CALL test%assert_close( e, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(s, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    !
  END SUBROUTINE real_1
  !
  SUBROUTINE complex_1(test)
    USE LAXlib
    USE la_param, ONLY : DP
    implicit none
    !
    TYPE(tester_t) :: test
    ! real variables
    complex(DP) :: h(2,2)
    complex(DP) :: h_save(2,2)
    complex(DP) :: s(2,2)
    complex(DP) :: s_save(2,2)
    real(DP) :: e(2)
    complex(DP) :: v(2,2)
    !
    h = 0.d0
    h(1,1) = (1.d0,  0.d0)
    h(1,2) = (0.d0, -2.d0)
    h(2,1) = ( 0.d0, 2.d0)
    h(2,2) = ( 5.d0, 0.d0)
    s = 0.d0
    s(1,1) = (1.d0, 0.d0)
    s(2,2) = (1.d0, 0.d0)
    !
    h_save = h
    s_save = s
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diaghg(  2, 2, h, s, 2, e, v, .false. )
    ! 0.1715728752538099, 5.82842712474619 
    CALL test%assert_close( e, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( v(:,1), [( 0.d0, -0.9238795325112867d0), (-0.3826834323650898d0, 0.d0)] )
    CALL test%assert_close( v(:,2), [( 0.d0, -0.3826834323650898d0), ( 0.9238795325112867d0, 0.d0)] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    CALL test%assert_close( RESHAPE(s, [4]), RESHAPE(s_save, [4]))
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diaghg(  2, 2, h, s, 2, e, v, .true. )
    !
    CALL test%assert_close( e, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( v(:,1), [( 0.d0, -0.9238795325112867d0), (-0.3826834323650898d0, 0.d0)] )
    CALL test%assert_close( v(:,2), [( 0.d0, -0.3826834323650898d0), ( 0.9238795325112867d0, 0.d0)] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    CALL test%assert_close( RESHAPE(s, [4]), RESHAPE(s_save, [4]))
    !
  END SUBROUTINE complex_1
end program test_diaghg
