!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"
!
!@@@@@
      real(kind=8) function ssum(n,a,nstride)
!
!  wrapper routine for cray scilib function ssum
!
      implicit none
      integer n, nstride
      real(kind=8) a(*)
      integer i
!
      ssum=0.d0
      do i=1,n,nstride
         ssum=ssum+a(i)
      end do
!
      return
      end
!
!@@@@@
      subroutine mxma (a,na,iad,b,nb,ibd,c,nc,icd,nar,nac,nbc)
!
!  wrapper routine for cray scilib matrix-matrix multiplication
!  routine mxma: c=a*b . Uses blas routine dgemm
!  na, nb, nc  = spacing between column elements of a, b ,c resp.
!  iad,ibd,icd = spacing between    row elements of a, b ,c resp.
!  nar=number of rows of a and c
!  nac=number of columns of a, number of rows of b
!  nbc=number of columns of b and c
!
      use io_global, only: stdout

      implicit none
      integer na, iad, nb, ibd, nc, icd, nar, nac, nbc
      real(8) a(iad,nac), b(ibd,nbc), c(icd,nbc)
      character(len=1) mode1, mode2
      integer lda, ldb
!
!  fortran equivalent (a,b,c are one-dimensional arrays)
!
!      real(8) a(iad*nac), b(ibd*nbc), c(icd*nbc)
!      integer i,j,k
!
!      do j=1,nbc
!         do i=1,nar
!            c((i-1)*nc+(j-1)*icd+1)=0.d0
!            do k=1,nac
!               c((i-1)*nc+(j-1)*icd+1) = c((i-1)*nc+(j-1)*icd+1)        &
!     &                                 + a((i-1)*na+(k-1)*iad+1)        &
!     &                                 * b((k-1)*nb+(j-1)*ibd+1)
!            end do
!         end do
!      end do
!
      if ( na.ne.1.and.iad.ne.1 .or.                                    &
     &     nb.ne.1.and.ibd.ne.1 .or. nc.ne.1 ) then
         WRITE( stdout,'(''MXMA : na,nb,nc,iad,ibd,icd,nar,nac,nbc =''/      &
     &                9i8)') na,nb,nc,iad,ibd,icd,nar,nac,nbc
         WRITE( stdout,'(''MXMA : not implemented'')')
         stop
      end if
!
      if (na.eq.1) then
         mode1='N'
         lda=iad
      else if (na.ne.1.and.iad.eq.1) then
         mode1='T'
         lda=na
      end if
!
      if (nb.eq.1) then
         mode2='N'
         ldb=ibd
      else if (nb.ne.1.and.ibd.eq.1) then
         mode2='T'
         ldb=nb
      end if
!
! call to BLAS3 routine GEMM
!
      call DGEMM                                                         &
     &     (mode1,mode2,nar,nbc,nac,1.d0,a,lda,b,ldb,0.d0,c,icd)
!
      return
      end subroutine mxma
