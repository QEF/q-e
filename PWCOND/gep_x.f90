!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine gep_x(n, amt, bmt, eigen, veigen)
!
! It solves GEP: A X = lambda B X using LAPACK routines
!
  USE kinds, only : DP
  USE cond, only : delgep
  implicit none

  integer :: i, n, info, lwork
  complex(DP) :: trywork
  real(DP), allocatable :: rwork(:)
  complex(DP) ::                                   &
                 amt(n,n),  &  ! A
                 bmt(n,n),  &  ! B
                 eigen(n),  &  ! lambda
                 veigen(n,n)   ! X
  complex(DP), allocatable :: alpha(:), beta(:), work(:)

  allocate( alpha( n ) )
  allocate( beta( n ) )
  allocate( rwork( 8*n ) )
!
!  for some reason the lapack routine does not work if the diagonal elements
!  of the matrices are exactly zero. We tested these routines on
!  pc_ifc, ibmsp and origin.
!
  do i=1,n
     amt(i,i)=amt(i,i)+delgep
     bmt(i,i)=bmt(i,i)+delgep
  enddo

  call ZGGEV('N', 'V', n, amt, n, bmt, n, alpha, beta, veigen, n, veigen, &
              n, trywork, -1, rwork, info)

  lwork=abs(trywork)
  allocate( work( lwork ) )

  call ZGGEV('N', 'V', n, amt, n, bmt, n, alpha, beta, veigen, n, veigen, &
              n, work, lwork, rwork, info)

  if(info.ne.0) call errore('gep_x','error on zggev',info)
!
! lambda = alpha / beta
!
  do i=1, n
    eigen(i)=alpha(i)/beta(i)
!    write(6,'(i5, 2f40.20)') i,  DBLE(eigen(i)), AIMAG(eigen(i))
  enddo

  deallocate(work)
  deallocate(rwork)
  deallocate(beta)
  deallocate(alpha)

  return
end subroutine gep_x
