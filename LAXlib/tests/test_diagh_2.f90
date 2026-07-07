program test_diagh_2

    USE laxlib_parallel_include
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    USE test_helpers,  ONLY : hermitian, symmetric
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
    CALL real_1(test)
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
    real(DP)    :: e(m_size)
    complex(DP) :: v(m_size,m_size)
    real(DP)    :: e_save(m_size)
    complex(DP) :: v_save(m_size,m_size)
    integer :: j
    !
    CALL hermitian(m_size, h)
    !
    h_save = h
    !
    v = (0.d0, 0.d0)
    e = 0.d0
    !
    CALL diagh(  m_size, m_size, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, m_size
       CALL test%assert_close( h(1:m_size, j), h_save(1:m_size, j))
    END DO
    !
    e_save = e
    v_save = v
    !
    ! Test that calling again gives the same results
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diagh(  m_size, m_size, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, m_size
       CALL test%assert_close( h(1:m_size, j), h_save(1:m_size, j))
    END DO
    !
    test%tolerance32=1.e-5
    test%tolerance64=1.d-14
    CALL test%assert_close( e, e_save)
    !
    ! Test subset of eigenvalues
    v = (0.d0, 0.d0)
    e = 0.d0
    CALL diagh(  m_size, m_size/2, h, e, v(:,1:m_size/2), me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    test%tolerance32=1.e-5
    ! Use a looser tolerance for subset eigenvalues since ZHEEVD and ZHEEVX
    ! use different algorithms and may produce slightly different results
    test%tolerance64=1.d-12
    CALL test%assert_close( e(1:m_size/2), e_save(1:m_size/2))
    !
  END SUBROUTINE complex_1
  !
  SUBROUTINE real_1(test)
    USE LAXlib
    implicit none
    !
    TYPE(tester_t) :: test
    !
    integer, parameter :: m_size=1024
    real(DP) :: h(m_size,m_size)
    real(DP) :: h_save(m_size,m_size)
    real(DP) :: e(m_size)
    real(DP) :: v(m_size,m_size)
    real(DP) :: e_save(m_size)
    real(DP) :: v_save(m_size,m_size)
    integer :: j
    !
    CALL symmetric(m_size, h)
    !
    h_save = h
    !
    v = 0.d0
    e = 0.d0
    !
    CALL diagh(  m_size, m_size, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, m_size
       CALL test%assert_close( h(1:m_size, j), h_save(1:m_size, j))
    END DO
    !
    e_save = e
    v_save = v
    !
    ! Test that calling again gives the same results
    v = 0.d0
    e = 0.d0
    CALL diagh(  m_size, m_size, h, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    DO j = 1, m_size
       CALL test%assert_close( h(1:m_size, j), h_save(1:m_size, j))
    END DO
    !
    test%tolerance32=1.e-5
    test%tolerance64=1.d-14
    CALL test%assert_close( e, e_save)
    !
  END SUBROUTINE real_1
  !
end program test_diagh_2
