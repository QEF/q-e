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
!
SUBROUTINE solve_with_dsygvd(n, v, s, ldh, e)
  USE la_param, ONLY : DP
  IMPLICIT NONE
  !
  REAL(DP) :: v(ldh,n)
  REAL(DP) :: s(ldh,n)
  real(DP) :: e(n)
  INTEGER                  :: n
  !
  INTEGER                  :: lwork, liwork, info, ldh
  !
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: work(:)
  ! various work space

  !
  ALLOCATE(work(1), iwork(1))
  CALL dsygvd( 1, 'V', 'U', n, v, ldh, &
              s, ldh, e, work, -1, iwork, -1, info )
  
  IF (info /= 0) print *, "Workspace not computed!"
  
  lwork = work(1)
  liwork = iwork(1)
  
  DEALLOCATE(work, iwork)
  ALLOCATE(work(lwork), iwork(liwork))
  !
  CALL dsygvd( 1, 'V', 'U', n, v, ldh, &
              s, ldh, e, work, lwork, iwork, liwork, info )
  !
  DEALLOCATE(work, iwork)
  !
END SUBROUTINE solve_with_dsygvd
!
! TODO: add check for eigenvalue probelm
