program test_diaghg_2

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
    !
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
    implicit none
    !
    TYPE(tester_t) :: test
    ! 
    integer, parameter :: m_size=1024
    complex(DP) :: h(m_size,m_size)
    complex(DP) :: h_save(m_size,m_size)
    complex(DP) :: s(m_size,m_size)
    complex(DP) :: s_save(m_size,m_size)
    real(DP)    :: e(m_size)
    complex(DP) :: v(m_size,m_size)
    real(DP)    :: e_save(m_size)
    complex(DP) :: v_save(m_size,m_size)
    integer :: j
    !
    CALL hermitian(m_size, h)
    CALL hermitian(m_size, s)
    !
    h_save = h
    s_save = s
    !
    v = (0.d0, 0.d0)
    e = 0.d0

    CALL diaghg(  m_size, m_size, h, s, m_size, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .false. )
    ! 
    DO j = 1, m_size
       CALL test%assert_close( h(1:m_size, j), h_save(1:m_size, j))
       CALL test%assert_close( s(1:m_size, j), s_save(1:m_size, j))
    END DO
    !
    e_save = e
    v_save = v
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diaghg(  m_size, m_size, h, s, m_size, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, .true. )
    !
    DO j = 1, m_size
       CALL test%assert_close( h(1:m_size, j), h_save(1:m_size, j))
       CALL test%assert_close( s(1:m_size, j), s_save(1:m_size, j))
    END DO

    test%tolerance32=1.e-5
    test%tolerance64=1.d-14
    CALL test%assert_close( e, e_save)
    !

  END SUBROUTINE complex_1
  !
  SUBROUTINE hermitian(mSize, M)
    IMPLICIT NONE
    integer, intent(in) :: mSize
    complex(dp), intent(out) :: M(mSize,mSize)
    !       
    real(dp), allocatable :: rnd(:)
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
  
end program test_diaghg_2
