!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if ! defined __AIX
      SUBROUTINE DGEMUL(A,LDA,FORMA,B,LDB,FORMB,C,LDC,L,M,N)
#else
      SUBROUTINE DGEMUL_WRAPPER(A,LDA,FORMA,B,LDB,FORMB,C,LDC,L,M,N)
#endif
      IMPLICIT NONE
      CHARACTER*1 FORMA,FORMB
      REAL*8 A,B,C
      INTEGER LDA,LDB,LDC,L,M,N
      DIMENSION A(*),B(*),C(*)
      CALL DGEMM(FORMA,FORMB,L,N,M,1.0D0,A,LDA,B,LDB,0.0D0,C,LDC)
      RETURN
      END 
