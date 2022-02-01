program test_diaghg_3

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
    test%tolerance64=1.d-8
    !
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
    USE mp_world, ONLY : mpime
    USE LAXlib
    USE test_io
    implicit none
    !
    TYPE(tester_t) :: test
    ! 
    integer :: ldh, n, m
    real(DP), allocatable :: h(:,:)
    real(DP), allocatable :: h_save(:,:)
    real(DP), allocatable :: s(:,:)
    real(DP), allocatable :: s_save(:,:)
    real(DP), allocatable    :: e(:)
    real(DP), allocatable :: v(:,:)
    real(DP), allocatable    :: e_save(:)
    real(DP), allocatable :: v_save(:,:)
    !
    character(len=20)        :: inputs(2)
    integer :: i, info
    !
    inputs = ["ZnOG1.bin", "ZnOG2.bin"]
    !
    DO i=1, SIZE(inputs)
        !
        CALL read_problem(inputs(i), ldh, n, m, h, s, e, v, info)
        !
        IF (info /= 0) THEN
            IF (mpime == 0) print *, "Test with ", inputs(i), " skipped. Input not found."
            CYCLE
        END IF
        !
        ALLOCATE(h_save, SOURCE=h)
        ALLOCATE(s_save, SOURCE=s)
        ALLOCATE(e_save, SOURCE=e)
        ALLOCATE(v_save, SOURCE=v)
        !
        h_save = h
        s_save = s
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg(  n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
        !
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        h = h_save
        s = s_save
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
        !
        CALL test%assert_close( e(1:m), e_save(1:m))
        DEALLOCATE(h,s,e,v,h_save,s_save,e_save,v_save)
    END DO
    !
  END SUBROUTINE real_1
  !
  SUBROUTINE complex_1(test)
    USE mp_world, ONLY : mpime
    USE LAXlib
    USE test_io
    implicit none
    !
    TYPE(tester_t) :: test
    ! 
    integer :: ldh, n, m
    complex(DP), allocatable :: h(:,:)
    complex(DP), allocatable :: h_save(:,:)
    complex(DP), allocatable :: s(:,:)
    complex(DP), allocatable :: s_save(:,:)
    real(DP), allocatable    :: e(:)
    complex(DP), allocatable :: v(:,:)
    real(DP), allocatable    :: e_save(:)
    complex(DP), allocatable :: v_save(:,:)
    !
    character(len=20)        :: inputs(4)
    integer :: i, info
    !
    inputs = ["ZnOK1.bin ", &
              "ZnOK2.bin ", &
              "SiGeK1.bin", &
              "SiGeK2.bin"]
    !
    DO i=1, SIZE(inputs)
        !
        CALL read_problem(inputs(i), ldh, n, m, h, s, e, v, info)
        !
        IF (info /= 0) THEN
            IF (mpime == 0) print *, "Test with ", inputs(i), " skipped. Input not found."
            CYCLE
        END IF
        !
        ALLOCATE(h_save, SOURCE=h)
        ALLOCATE(s_save, SOURCE=s)
        ALLOCATE(e_save, SOURCE=e)
        ALLOCATE(v_save, SOURCE=v)
        !
        h_save = h
        s_save = s
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg(  n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
        !
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        h = h_save
        s = s_save
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
        !
        CALL test%assert_close( e(1:m), e_save(1:m))
        DEALLOCATE(h,s,e,v,h_save,s_save,e_save,v_save)
    END DO
    !
  END SUBROUTINE complex_1
  
end program test_diaghg_3
