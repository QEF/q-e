!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cgramg1 (lda, nvecx, n, start, finish, psi, spsi, hpsi)  
  !-----------------------------------------------------------------------
  !
  !      This routine orthogonalizes several vectors with the method of
  !      Gram-Schmidt and imposing that  <psi_i|S|psi_j> = delta_ij.
  !      It receives on input the psi and the spsi.
  !      It updates also the Hamiltonian so that it contains the new hpsi.
  !
#include "machine.h"
  use parameters
  use allocate 
  implicit none  
  !
  !     first the dummy variables
  !
  integer :: lda, n, nvecx, start, finish  
  ! input: leading dimension of the vectors
  ! input: physical dimension
  ! input: dimension of psi
  ! input: first vector to orthogonalize
  ! input: last vector to orthogonalize
  complex(kind=DP) :: psi (lda, nvecx), spsi (lda, nvecx), hpsi (lda, &
       nvecx)
  ! inp/out: the vectors to be orthogonalized
  !
  !    two parameters
  !
  integer :: ierrx  
  ! maximum number of errors
  real(kind=DP) :: eps  
  ! a small number
  parameter (ierrx = 3, eps = 1.0d-6)  
  !
  !    here the local variables
  !
  integer :: vec, vecp, ierr  
  ! counter on vectors
  ! counter on vectors
  ! counter on errors
  complex(kind=DP), pointer :: ps (:)
  complex(kind=DP) ::  ZDOTC  
  ! the scalar products
  ! function which computes scalar products

  real(kind=DP) :: norm, DDOT  
  ! the norm of a vector
  ! function computing the dot product of two v

  call start_clock ('cgramg1')  

  call mallocate(ps, finish)  
  do vec = start, finish  
     ierr = 0  
1    continue  
     do vecp = 1, vec  
        ps (vecp) = ZDOTC (n, psi (1, vecp), 1, spsi (1, vec), 1)  
     enddo
#ifdef PARA
     call reduce (2 * vec, ps)  
#endif
     do vecp = 1, vec - 1  
        call ZAXPY (n, - ps (vecp), psi (1, vecp), 1, psi (1, vec), 1)
        call ZAXPY (n, - ps (vecp), spsi (1, vecp), 1, spsi (1, vec), 1)
        call ZAXPY (n, - ps (vecp), hpsi (1, vecp), 1, hpsi (1, vec), 1)
     enddo
     norm = ps (vec) - DDOT (2 * (vec - 1), ps, 1, ps, 1)  
     !        print*,norm,vec
     !        norm = DDOT( 2*n, psi(1,vec), 1, spsi(1,vec), 1 )
#ifdef PARA
     !        call reduce (1,norm)
#endif
     if (norm.lt.0.d0) then  
        print * , 'norma = ', norm, vec  
        call error ('cgramg1', ' negative norm in S ', 1)  
     endif
     norm = 1.d0 / sqrt (norm)  
     call DSCAL (2 * n, norm, psi (1, vec), 1)  
     call DSCAL (2 * n, norm, spsi (1, vec), 1)  
     call DSCAL (2 * n, norm, hpsi (1, vec), 1)  
     if (1.d0 / norm.lt.eps) then  
        ierr = ierr + 1  
        if (ierr.le.ierrx) goto 1  
        call error ('cgramg1', ' absurd correction vector', vec)  
     endif

  enddo
  call mfree(ps)  
  call stop_clock ('cgramg1')  
  return  
end subroutine cgramg1

