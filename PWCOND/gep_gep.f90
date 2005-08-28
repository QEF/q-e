!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine gep_gep(n, amt, bmt, eigen, veigen)
!
! It solves GEP: A X = lambda B X using routines in GEP.f
!
#include "f_defs.h"
  USE kinds, only : DP
  use cond, only: delgep
implicit none  
  integer :: i, n
  integer, allocatable :: ipiv(:) 
  complex(DP) ::  &
        amt(n,n),      & ! A  
        bmt(n,n),      & ! B  
        eigen(n),      & ! lambda 
        veigen(n,n)      ! X
  complex(DP), allocatable ::                      & 
        alpha(:), beta(:)

  do i=1, n
    amt(i,i) = amt(i,i)+delgep
    bmt(i,i) = bmt(i,i)+delgep
  enddo

  allocate( alpha( n ) )
  allocate( beta( n ) )
  allocate( ipiv( n ) )
  call LZHES(n, amt, n, bmt, n, veigen, n, .true.)
  call LZIT(n, amt, n, bmt, n, veigen, n, .true.,       &
                            ipiv, alpha, beta) 
!
! lambda = alpha / beta 
!
  do i=1, n
    eigen(i)=alpha(i)/beta(i)
  enddo
  deallocate(alpha)
  deallocate(beta)
  deallocate(ipiv) 


  return
end subroutine gep_gep 
