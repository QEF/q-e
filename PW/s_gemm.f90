!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef T3D_BENCHLIB


SUBROUTINE S_gemm (FORMA, FORMB, L, N, M, alpha, A, LDA, B, LDB, &
     beta, C, LDC)
  ! driver for SGEMM routine
  ! it transposes matrixes before calling the SGEMM subroutine in the T/N
  ! or T T mode in order to exploit cache and stream buffer mechanics
  ! written by C.Cavazzoni and S.Cozzini 18/12/97
  ! the calling parameter are the same as in the actual sgemm routines
  ! to use it, just replace SGEMM with S_gemm in the original code...
  IMPLICIT NONE  
  CHARACTER (LEN=1) :: FORMA, FORMB  
  real (8) :: alpha, beta  
  REAL (8) :: A, B, C  
  real (8) , allocatable::auxa (:, :)  
  real (8) , allocatable::auxb (:, :)  
  real (8) , allocatable::auxc (:, :)  
  integer :: i, j  
  INTEGER :: LDA, LDB, LDC, L, M, N  

  DIMENSION A (LDA, * ), B (LDB, * ), C (LDC, * )  
  INTEGER :: IRTC  
  INTEGER :: MM, NN, KK  
  INTEGER :: S2, S3, itype  

  INTEGER :: INFO (2)  
  MM = L  
  NN = N  
  KK = M  
  ! defines the algorithm to transpose matrix

  itype = 0  
  if (forma.eq.'t'.or.forma.eq.'T') then  
     if (formb.eq.'t'.or.formb.eq.'T') then  
        allocate (auxc (MM, NN) )  
        S2 = IRTC ()  
        CALL SGEMM ('N', 'N', MM, NN, KK, alpha, B, LDB, A, LDA, &
             beta, AUXC, MM)
        S3 = IRTC ()  
        call transpose (AUXC, MM, C, LDC, NN, MM, itype, INFO (2) )  
        deallocate (auxc)  
     else  
        allocate (auxa (MM, KK) )  
        call transpose (A, LDA, AUXA, MM, MM, KK, itype, INFO (2) )  
        S2 = IRTC ()  
        CALL SGEMM ('N', 'N', MM, NN, KK, alpha, auxa, MM, B, LDB, &
             beta, C, LDC)
        S3 = IRTC ()  
        deallocate (auxa)  
     endif
  elseif (formb.eq.'t'.or.formb.eq.'T') then  
     allocate (auxb (KK, NN) )  
     call transpose (B, LDB, AUXB, KK, KK, NN, itype, INFO (2) )  
     S2 = IRTC ()  
     CALL SGEMM ('N', 'N', MM, NN, KK, alpha, A, LDA, auxb, KK, &
          beta, C, LDC)
     S3 = IRTC ()  
     deallocate (auxb)  
  else  
     S2 = IRTC ()  
     CALL SGEMM (forma, formb, MM, NN, KK, alpha, A, LDA, B, LDB, &
          beta, C, LDC)
     S3 = IRTC ()  
     INFO (2) = 0  

  endif
  RETURN  




END SUBROUTINE S_gemm

subroutine transpose (A, LDA, B, LDB, N, M, ITYPE, INFO)  

  IMPLICIT NONE  

  INCLUDE "mpp/shmem.fh"  
  integer :: lda, ldb, n, m, itype, info  



  real (8) :: A (LDA, * ), B (LDB, * )  
  !     DRIVER FOR MATRIX TRASPOSITIN
  !
  !     A  (INPUT ) MxN MATRIX TO BE TRANSPOSED
  !     B  (OUTPUT) NxM TRANSPOSED MATRIX
  !
  !     ITYPE (INPUT )  TYPE OF TRANSPOSITION METHOD
  !     INFO  (OUTPUT)  ROUTINE EXECUTION TIME ( MICROSECONDS )
  integer :: irtc, my_pe, lputp  
  integer :: i, j, jn, mype  

  integer :: s2, s3  

  if (itype.eq.0) then  
     S2 = IRTC ()  
     do j = 1, N, 480  
        jn = min (n - j + 1, 480)  
        do i = 1, M  
           call lgetv (A (I, J), LDA, JN)  
           call lputv (B (J, I), 1, JN)  
        enddo
     enddo
10   IF (LPUTP () .NE.0) GoTo 10  

     S3 = IRTC ()  

  elseif (itype.eq.1) then  
     ! traspose with shmem
     MYPE = MY_PE ()  
     S2 = IRTC ()  
     DO J = 1, M  
        CALL SHMEM_IGET (B (1, J), A (J, 1), 1, LDA, N, MYPE)  
     ENDDO

     S3 = IRTC ()  

  elseif (itype.eq.2) then  
     S2 = irtc ()  
     DO J = 1, M  
        !DIR$ CACHE_BYPASS A,B
        DO I = 1, N  
           B (I, J) = A (J, I)  
        ENDDO
     ENDDO
     S3 = IRTC ()  

  else  
     write (6,  * ) ' *** TRANSPOSE : PARAMETER ITYPE ', 'OUT OF RANGE'  

     stop  

  endif

  INFO = INT (DBLE (S3 - S2) * 3.333d-3)  
  return  
end subroutine transpose
#else
subroutine sgemmdummy  
  return  
end subroutine sgemmdummy
#endif
