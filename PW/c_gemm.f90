!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef __BENCHLIB

subroutine c_gemm (forma, formb, l, n, m, alpha, a, lda, b, ldb, &
     beta, c, ldc)
  ! driver for CGEMM routine
  ! it conjugate matrixes before calling the SGEMM subroutine in the T/N o
  ! or T T mode in order to exploit cache and stream buffer mechanics
  ! the calling parameter are the same as in the actual sgemm routines
  ! to use it, just replace ZGEMM with c_gemm in the original code...
  implicit none
  character (len=1) :: forma, formb
  real (8) :: alpha, beta
  complex (8) :: a, b, c
  complex (8) , allocatable::auxa (:, :)
  complex (8) , allocatable::auxb (:, :)
  complex (8) , allocatable::auxc (:, :)
  integer :: i, j
  integer :: lda, ldb, ldc, l, m, n

  dimension a (lda, * ), b (ldb, * ), c (ldc, * )
  integer :: irtc
  integer :: mm, nn, kk
  integer :: s2, s3, itype

  integer :: info (2)
  mm = l
  nn = n
  kk = m
  ! defines the algorithm to transpose matrix

  itype = 2
  if (forma.eq.'c'.or.forma.eq.'C') then
     if (formb.eq.'c'.or.formb.eq.'C') then
        allocate (auxc (mm, nn) )
        s2 = irtc ()
        call cgemm ('n', 'n', mm, nn, kk, alpha, b, ldb, a, lda, &
             beta, auxc, mm)
        s3 = irtc ()
        call c_transpose (auxc, mm, c, ldc, nn, mm, itype, info (2) )
        deallocate (auxc)
     else
        allocate (auxa (mm, kk) )
        call c_transpose (a, lda, auxa, mm, mm, kk, itype, info (2) )
        s2 = irtc ()
        call cgemm ('n', 'n', mm, nn, kk, alpha, auxa, mm, b, ldb, &
             beta, c, ldc)
        s3 = irtc ()
        deallocate (auxa)
     endif
  elseif (formb.eq.'c'.or.formb.eq.'C') then
     allocate (auxb (kk, nn) )
     call c_transpose (b, ldb, auxb, kk, kk, nn, itype, info (2) )
     s2 = irtc ()
     call cgemm ('n', 'n', mm, nn, kk, alpha, a, lda, auxb, kk, &
          beta, c, ldc)
     s3 = irtc ()
     deallocate (auxb)
  else
     s2 = irtc ()
     call cgemm (forma, formb, mm, nn, kk, alpha, a, lda, b, ldb, &
          beta, c, ldc)
     s3 = irtc ()
     info (2) = 0

  endif
  return

end subroutine c_gemm

subroutine c_transpose (a, lda, b, ldb, n, m, itype, info)

  implicit none

  include "mpp/shmem.fh"
  integer :: lda, ldb, n, m, itype, info


  complex (8) :: a (lda, * ), b (ldb, * )
  !     driver for matrix trasposition
  !
  !     a  (input ) m x n matrix to be transposed
  !     b  (output) n x m transposed matrix
  !
  !     itype (input )  type of transposition method
  !     info  (output)  routine execution time ( microseconds )
  integer :: irtc, my_pe, lputp
  integer :: i, j, jn, mype

  integer :: s2, s3

  if (itype.eq.0) then
     s2 = irtc ()
     do j = 1, n, 480
        jn = min (n - j + 1, 480)
        do i = 1, m
           call lgetv (a (i, j), lda, jn)
           call lputv (b (j, i), 1, jn)
        enddo
     enddo
10   if (lputp () .ne.0) goto 10

     s3 = irtc ()


  elseif (itype.eq.1) then
     ! traspose with shmem
     mype = my_pe ()
     s2 = irtc ()
     do j = 1, m
        call shmem_iget (b (1, j), a (j, 1), 1, lda, n, mype)
     enddo

     s3 = irtc ()

  elseif (itype.eq.2) then
     s2 = irtc ()
     do j = 1, m
        do i = 1, n
           b (i, j) = CONJG(a (j, i) )
        enddo
     enddo

     s3 = irtc ()

  else
     WRITE( stdout, * ) '*** c_transpose : parameter itype out of range'

     stop

  endif

  info = int (dble (s3 - s2) * 3.333d-3)
  return

end subroutine c_transpose
#else
subroutine cgemmdummy
  return
end subroutine cgemmdummy
#endif
