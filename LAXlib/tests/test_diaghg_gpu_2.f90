#if defined(__CUDA)
program test_diaghg_gpu_2
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
  SUBROUTINE complex_1(test)
    USE LAXlib
    USE la_param, ONLY : DP
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    ! 
    integer, parameter :: m_size=1024
    integer :: i
    complex(DP) :: h(m_size,m_size)
    complex(DP), DEVICE :: h_d(m_size,m_size)
    complex(DP) :: h_save(m_size,m_size)
    complex(DP) :: s(m_size,m_size)
    complex(DP), DEVICE :: s_d(m_size,m_size)
    complex(DP) :: s_save(m_size,m_size)
    real(DP) :: e(m_size), e_(m_size)
    real(DP), DEVICE :: e_d(m_size)
    complex(DP) :: v(m_size,m_size)
    complex(DP), DEVICE :: v_d(m_size,m_size)
    real(DP) :: e_save(m_size)
    complex(DP) :: v_save(m_size,m_size)
    !
    CALL hermitian(m_size, h)
    CALL hermitian(m_size, s)
    !
    h_save = h
    s_save = s
    h_d = h     ! copy H and S to device
    s_d = s     ! <----------|
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    v_d = (0.d0, 0.d0)
    e_d = 0.d0
    !
    ! Compare same algorithm starting from data on device ...
    CALL diaghg(  m_size, m_size-1, h_d, s_d, m_size, e_d, v_d )
    e_save = e_d
    v_save = v_d
    ! ... and on the host
    CALL diaghg(  m_size, m_size-1, h, s, m_size, e, v, .true. )
    !
    CALL test%assert_close( RESHAPE(h, [m_size*m_size]), RESHAPE(h_save, [m_size*m_size]))
    CALL test%assert_close( RESHAPE(s, [m_size*m_size]), RESHAPE(s_save, [m_size*m_size]))
    test%tolerance32=1.d-5
    test%tolerance64=1.d-10
    DO i=1, m_size-1
        CALL test%assert_close( v(1:m_size,i), v_save(1:m_size,i))
        CALL test%assert_close( e(i), e_save(i) )
    END DO
    !
    h = h_save
    s = s_save
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diaghg( m_size, m_size-1, h, s, m_size, e, v )
    h = h_save; s = s_save
    CALL solve_with_zhegvd(m_size, h, s, m_size, e_)
    !
    test%tolerance32=1.d-5
    test%tolerance64=1.d-10
    DO i=1, m_size-1
        CALL test%assert_close( v_save(1:m_size,i), h(1:m_size,i) )
        CALL test%assert_close( e(i), e_save(i) )
        CALL test%assert_close( e_(i), e_save(i) )
    END DO
    !
  END SUBROUTINE complex_1
  !
  SUBROUTINE hermitian(mSize, M)
    USE la_param, ONLY : DP
    IMPLICIT NONE
    integer, intent(in) :: msize
    complex(dp), intent(out) :: M(:,:)
    !       
    real(dp), allocatable :: rnd(:)
    complex(dp), allocatable :: tmp(:,:)

    INTEGER :: h, k, j
    !
    ALLOCATE(rnd(mSize*(mSize+1)))
    CALL RANDOM_NUMBER(rnd)
    rnd = 1.d0*rnd - 5.d-1
    !
    M = (0.d0, 0.d0)
    j = 1
    DO k=1,mSize
      DO h=1,mSize
        IF(h>k) THEN
          M(h,k) = CMPLX(rnd(j), rnd(j+1))
          M(k,h) = CONJG(M(h,k))
          j=j+2;
        ELSE IF(k == h) THEN
          M(k,h) = CMPLX(mSize, 0.d0, kind=DP)
        END IF
      END DO
    END DO
    !
  END SUBROUTINE hermitian
  !
  SUBROUTINE solve_with_zhegvd(n, v, s, ldh, e)
    USE la_param, ONLY : DP
    IMPLICIT NONE
    !
    complex(DP) :: v(ldh,n)
    complex(DP) :: s(ldh,n)
    real(DP) :: e(n)
    INTEGER                  :: n
    !
    INTEGER                  :: lwork, lrwork, liwork, info, ldh
    !
    REAL(DP)                 :: abstol
    INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
    REAL(DP),    ALLOCATABLE :: rwork(:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space

    !
    ALLOCATE(work(1), rwork(1), iwork(1))
    CALL ZHEGVD( 1, 'V', 'U', n, v, ldh, &
                s, ldh, e, work, -1, rwork, -1, iwork, -1, info )
    
    IF (info /= 0) print *, "Workspace not computed!"
    
    lwork = work(1)
    lrwork = rwork(1)
    liwork = iwork(1)
    
    DEALLOCATE(work, rwork, iwork)
    ALLOCATE(work(lwork), rwork(lrwork), iwork(liwork))
    
    CALL ZHEGVD( 1, 'V', 'U', n, v, ldh, &
                s, ldh, e, work, lwork, rwork, lrwork, iwork, liwork, info )
    
    DEALLOCATE(work, rwork, iwork)
    !
  END SUBROUTINE solve_with_zhegvd
  
end program test_diaghg_gpu_2
#else
program test_diaghg_gpu_2
end program test_diaghg_gpu_2
#endif
