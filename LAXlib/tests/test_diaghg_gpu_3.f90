#ifdef __CUDA
program test_diaghg_gpu_3
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
    USE la_param, ONLY : DP
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
    real(DP), allocatable :: e(:)
    real(DP), allocatable :: v(:,:)
    real(DP), allocatable :: e_save(:)
    real(DP), allocatable :: v_save(:,:)
    real(DP), allocatable :: e_dc(:)      ! stores results from divide and conquer
    real(DP), allocatable :: v_dc(:,:)    !
    ! device copies
    real(DP), allocatable, device :: h_d(:,:)
    real(DP), allocatable, device :: s_d(:,:)
    real(DP), allocatable, device :: e_d(:)
    real(DP), allocatable, device :: v_d(:,:)
    !
    !
    character(len=20) :: inputs(2)
    integer :: i, j, info
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
        CALL diaghg(  n, m, h, s, ldh, e, v, .false. )
        !
        test%tolerance64=1.d-6 ! <- check this
        DO j = 1, m
            CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        test%tolerance64=1.d-8 ! <- check this
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg( n, m, h, s, ldh, e, v, .true. )
        !
        ALLOCATE(v_dc, SOURCE=h_save)
        ALLOCATE(e_dc(n))
        s = s_save
        CALL solve_with_dsygvd(n, v_dc, s, ldh, e_dc)
        s = s_save
        !
        DO j = 1, m
            CALL test%assert_close( v(1:n, j), v_dc(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m))
        !
        ! GPU data & subroutines
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        s = s_save
        h = h_save
        ALLOCATE(e_d, SOURCE=e); ALLOCATE(v_d, SOURCE=v)
        ALLOCATE(h_d, SOURCE=h); ALLOCATE(s_d, SOURCE=s)
        !
        CALL diaghg(  n, m, h_d, s_d, ldh, e_d, v_d, .false. )
        !
        v(1:n, 1:m) = v_d(1:n, 1:m)
        e(1:m)      = e_d(1:m)
        DO j = 1, m
            CALL test%assert_close( v(1:n, j), v_dc(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        v_d = (0.d0, 0.d0)
        e_d = 0.d0
        s_d = s_save
        h_d = h_save
        !
        ! Start from data on the GPU and diagonalize on the CPU
        CALL diaghg( n, m, h_d, s_d, ldh, e_d, v_d, .true. )
        !
        v(1:n, 1:m) = v_d(1:n, 1:m)
        e(1:m)      = e_d(1:m)
        !
        test%tolerance64=1.d-6 ! <- check this
        DO j = 1, m
            CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        test%tolerance64=1.d-8 ! <- check this
        !
        CALL test%assert_close( e(1:m), e_save(1:m))
        !
        DEALLOCATE(h_d, s_d, e_d, v_d)
        DEALLOCATE(h,s,e,v,h_save,s_save,e_save,v_save, v_dc, e_dc)
    END DO
    !
  END SUBROUTINE real_1
  !
  SUBROUTINE complex_1(test)
    USE mp_world, ONLY : mpime
    USE LAXlib
    USE la_param, ONLY : DP
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
    real(DP), allocatable    :: e_dc(:)      ! stores results from divide and conquer
    complex(DP), allocatable :: v_dc(:,:)    !
    
    ! device copies
    complex(DP), allocatable, device :: h_d(:,:)
    complex(DP), allocatable, device :: s_d(:,:)
    real(DP), allocatable, device    :: e_d(:)
    complex(DP), allocatable, device :: v_d(:,:)
    !
    !
    character(len=20)        :: inputs(4)
    integer :: i, j, info
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
        ALLOCATE(e_dc(n))
        ALLOCATE(v_dc, SOURCE=h)
        !
        h_save = h
        s_save = s
        !
        ! == Check CPU interface without and with offloading ==
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg(  n, m, h, s, ldh, e, v, .false. )
        !
        DO j = 1, m
            CALL test%assert_close(  v(1:n, j), v_save(1:n, j) )
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        h = h_save; s = s_save;
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL diaghg( n, m, h, s, ldh, e, v, .true. )
        !
        ! N.B.: GPU eigensolver uses a different algorithm: zhegvd
        s = s_save; e_dc = 0.d0
        CALL solve_with_zhegvd(n, v_dc, s, ldh, e_dc)
        !
        DO j = 1, m
            CALL test%assert_close( v(1:n, j), v_dc(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m))
        !
        ! GPU data & subroutines
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        s = s_save
        h = h_save
        ALLOCATE(e_d, SOURCE=e); ALLOCATE(v_d, SOURCE=v)
        ALLOCATE(h_d, SOURCE=h); ALLOCATE(s_d, SOURCE=s)
        !
        CALL diaghg(  n, m, h_d, s_d, ldh, e_d, v_d, .false. )
        !
        v(1:n, 1:m) = v_d(1:n, 1:m)
        e(1:m)      = e_d(1:m)
        DO j = 1, m
            !CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        v_d = (0.d0, 0.d0)
        e_d = 0.d0
        s_d = s_save
        h_d = h_save
        CALL diaghg( n, m, h_d, s_d, ldh, e_d, v_d, .true. )
        !
        v(1:n, 1:m) = v_d(1:n, 1:m)
        e(1:m)      = e_d(1:m)
        !
        DO j = 1, m
            !CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m))
        !
        DEALLOCATE(h_d, s_d, e_d, v_d)
        DEALLOCATE(h,s,e,v,h_save,s_save,e_save,v_save)
        DEALLOCATE(e_dc, v_dc)
    END DO
    !
  END SUBROUTINE complex_1
  !
end program test_diaghg_gpu_3
#else
program test_diaghg_gpu_3
end program test_diaghg_gpu_3
#endif
