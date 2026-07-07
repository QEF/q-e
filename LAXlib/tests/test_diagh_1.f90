program test_diagh

    USE laxlib_parallel_include
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, world_comm
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
    implicit none
    !
    TYPE(tester_t) :: test
    ! real variables
    real(DP) :: h(2,2)
    real(DP) :: h_save(2,2)
    real(DP) :: e(2)
    real(DP) :: v(2,2)

    ! Test 1: Identity matrix
    h = 0.d0
    h(1,1) = 1.d0
    h(2,2) = 1.d0
    h_save = h
    !
    v = 0.d0
    e = 0.d0
    CALL diagh(  2, 2, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    CALL test%assert_close( e, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(v, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]) )
    !
    ! Test 2: Diagonal matrix with different eigenvalues
    h = 0.d0
    h(1,1) = 2.d0
    h(2,2) = 5.d0
    h_save = h
    !
    v = 0.d0
    e = 0.d0
    CALL diagh(  2, 2, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    CALL test%assert_close( e, [2.d0, 5.d0] )
    CALL test%assert_close( RESHAPE(v, [4]), [1.d0, 0.d0, 0.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]) )
    !
    ! Test 3: Off-diagonal elements
    h = 0.d0
    h(1,1) = 1.d0
    h(1,2) = 2.d0
    h(2,1) = 2.d0
    h(2,2) = 1.d0
    h_save = h
    !
    v = 0.d0
    e = 0.d0
    CALL diagh(  2, 2, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    ! Eigenvalues: 1-2 = -1, 1+2 = 3
    CALL test%assert_close( e, [-1.d0, 3.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]) )
    !
  END SUBROUTINE real_1
  !
  SUBROUTINE complex_1(test)
    USE LAXlib
    implicit none
    !
    TYPE(tester_t) :: test
    ! complex variables
    complex(DP) :: h(2,2)
    complex(DP) :: h_save(2,2)
    real(DP) :: e(2)
    complex(DP) :: v(2,2)
    !
    ! Test 1: Identity matrix
    h = (0.d0, 0.d0)
    h(1,1) = (1.d0, 0.d0)
    h(2,2) = (1.d0, 0.d0)
    h_save = h
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diagh(  2, 2, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    CALL test%assert_close( e, [1.d0, 1.d0] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    !
    ! Test 2: Hermitian matrix with complex off-diagonal
    h = (0.d0, 0.d0)
    h(1,1) = (1.d0,  0.d0)
    h(1,2) = (0.d0, -2.d0)
    h(2,1) = (0.d0,  2.d0)
    h(2,2) = (5.d0,  0.d0)
    h_save = h
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diagh(  2, 2, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    ! Eigenvalues: (1+5 +/- sqrt((1-5)^2 + 4*4))/2 = (6 +/- sqrt(16+16))/2 = (6 +/- sqrt(32))/2
    ! = (6 +/- 5.65685)/2 = 0.1715728752538099, 5.82842712474619
    CALL test%assert_close( e, [0.1715728752538099d0,  5.82842712474619d0] )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    !
    ! Test 3: Subset of eigenvalues (m < n)
    h = (0.d0, 0.d0)
    h(1,1) = (1.d0, 0.d0)
    h(2,2) = (5.d0, 0.d0)
    h_save = h
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    ! Request only 1 eigenvalue
    CALL diagh(  2, 1, h, e, v(:,1:1), me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    CALL test%assert_close( e(1), 1.d0 )
    CALL test%assert_close( RESHAPE(h, [4]), RESHAPE(h_save, [4]))
    !
  END SUBROUTINE complex_1
end program test_diagh
