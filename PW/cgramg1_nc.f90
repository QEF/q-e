!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cgramg1_nc (lda, nvecx, n, start, finish, psi, spsi, hpsi, npol)
  !-----------------------------------------------------------------------
  !
  !      This routine orthogonalizes several vectors with the method of
  !      Gram-Schmidt and imposing that  <psi_i|S|psi_j> = delta_ij.
  !      It receives on input the psi and the spsi.
  !      It updates also the Hamiltonian so that it contains the new hpsi.
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  !     first the dummy variables
  !
  integer :: lda, n, nvecx, start, finish, npol
  ! input: leading dimension of the vectors
  ! input: physical dimension
  ! input: dimension of psi
  ! input: first vector to orthogonalize
  ! input: last vector to orthogonalize
  ! input: number of coordonates of wfc
  complex(DP) :: psi(lda,npol,nvecx), spsi(lda,npol,nvecx), hpsi(lda, &
       npol,nvecx)
  ! inp/out: the vectors to be orthogonalized
  !
  !    two parameters
  !
  integer :: ierrx
  ! maximum number of errors
  real(DP) :: eps
  ! a small number
  parameter (ierrx = 3, eps = 1.0d-6)
  !
  !    here the local variables
  !
  integer :: vec, vecp, ierr
  ! counter on vectors
  ! counter on vectors
  ! counter on errors
  complex(DP), allocatable :: ps (:)
  complex(DP) ::  ZDOTC
  ! the scalar products
  ! function which computes scalar products

  real(DP) :: norm, DDOT
  ! the norm of a vector
  ! function computing the dot product of two v

  call start_clock ('cgramg1')

  allocate (ps( finish))    
  do vec = start, finish
     ierr = 0
1    continue
     do vecp = 1, vec
        IF (npol ==1) THEN
           ps (vecp) = ZDOTC (n, psi (1, 1, vecp), 1, spsi (1, 1, vec), 1)
        ELSE
           ps (vecp) = ZDOTC (lda*npol, psi (1,1,vecp), 1, spsi (1,1,vec), 1)
        ENDIF
     enddo
#ifdef __PARA
     call reduce (2 * vec, ps)
#endif
     do vecp = 1, vec - 1
        IF (npol ==1) THEN
           call ZAXPY (n, - ps (vecp), psi (1, 1, vecp), 1, psi (1, 1, vec), 1)
           call ZAXPY (n, - ps (vecp), spsi (1,1,vecp), 1, spsi (1, 1, vec), 1)
           call ZAXPY (n, - ps (vecp), hpsi (1,1,vecp), 1, hpsi (1, 1, vec), 1)
        ELSE
           call ZAXPY (lda*npol, - ps(vecp),psi(1,1,vecp),1,psi(1,1,vec),1)
           call ZAXPY (lda*npol, - ps(vecp),spsi(1,1,vecp),1,spsi(1,1,vec),1)
           call ZAXPY (lda*npol, - ps(vecp),hpsi(1,1,vecp),1,hpsi(1,1,vec),1)
        ENDIF
     enddo
     norm = ps (vec) - DDOT (2 * (vec - 1), ps, 1, ps, 1)
     !        print*,norm,vec
     !        norm = DDOT( 2*n, psi(1,vec), 1, spsi(1,vec), 1 )
#ifdef __PARA
     !        call reduce (1,norm)
#endif
     if (norm.lt.0.d0) then
        print * , 'norma = ', norm, vec
        call errore ('cgramg1_nc', ' negative norm in S ', 1)
     endif
     norm = 1.d0 / sqrt (norm)
     IF (npol == 1) THEN
        call DSCAL (2 * n, norm, psi (1, 1, vec), 1)
        call DSCAL (2 * n, norm, spsi (1, 1, vec), 1)
        call DSCAL (2 * n, norm, hpsi (1, 1, vec), 1)
     ELSE
        call DSCAL (2 * lda * npol, norm, psi (1, 1, vec), 1)
        call DSCAL (2 * lda * npol, norm, spsi (1, 1, vec), 1)
        call DSCAL (2 * lda * npol, norm, hpsi (1, 1, vec), 1)
     ENDIF
     if (1.d0 / norm.lt.eps) then
        ierr = ierr + 1
        if (ierr.le.ierrx) goto 1
        call errore ('cgramg1_nc', ' absurd correction vector', vec)
     endif

  enddo
  deallocate(ps)
  call stop_clock ('cgramg1')
  return
end subroutine cgramg1_nc

