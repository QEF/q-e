!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! Author: Ivan Carnimeo
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! General-purpose routines for scalar products, printing,
! linear-algebra operators for exact-exchange and localization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matcalc (label, DoE, PrtMat, ninner, n, m, U, V, mat, ee)
  !
  USE kinds,    ONLY : dp
  USE io_global,ONLY : stdout
  USE becmod,   ONLY : calbec
  USE wvfct,    ONLY : current_k, wg
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  CHARACTER(len=*), INTENT(IN) :: label
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: PrtMat, ninner, n, m
  COMPLEX(DP), INTENT(IN) :: U(ninner,n), V(ninner,m)
  REAL(DP), INTENT(OUT) :: mat(n,m), ee
  !
  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock('matcalc')

  string = 'M-'
  mat = 0.0_dp
  CALL calbec(ninner, U, V, mat, m)

  IF( PrtMat .ge.2 ) CALL matprt(string//label,n,m,mat)

  IF(DoE) THEN
     IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
     string = 'E-'
     ee = 0.0_dp
     DO i = 1,n
        ee = ee + wg(i,current_k)*mat(i,i)
     ENDDO
     IF ( PrtMat .ge. 1 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE matcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matcalc_k (label, DoE, PrtMat, ik, ninner, n, m, U, V, mat, ee)
  !
  USE kinds,                ONLY : dp
  USE io_global,ONLY : stdout
  USE wvfct,                ONLY : wg, npwx
  USE becmod,               ONLY : calbec
  USE noncollin_module,     ONLY : noncolin, npol
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: PrtMat, ik, ninner, n, m
  COMPLEX(dp), INTENT(IN) :: U(ninner,n), V(ninner,m)
  COMPLEX(dp), INTENT(OUT):: mat(n,m)
  REAL(DP), INTENT(OUT) :: ee
  CHARACTER(len=*), INTENT(IN) :: label

  INTEGER :: i
  CHARACTER(len=2) :: string

  CALL start_clock('matcalc')

  string = 'M-'
  mat = (0.0_dp, 0.0_dp)
  IF(noncolin) THEN
    noncolin = .false.
    CALL calbec(ninner, U, V, mat, m)
    noncolin = .true.
  ELSE
    CALL calbec(ninner, U, V, mat, m)
  ENDIF

  IF(DoE) THEN
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    IF( PrtMat > 1 ) CALL matprt_k(string//label,n,m,mat)
    string = 'E-'
    ee = 0.0_dp
    DO i = 1,n
      ee = ee + wg(i,ik)*DBLE(mat(i,i))
    ENDDO
    IF ( PrtMat > 0 ) WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE matcalc_k

#if defined(__CUDA)
SUBROUTINE matcalc_k_gpu (DoE,ik, ninner, n, m, U, V, mat, ee)
  !
  USE kinds,                ONLY : dp
  USE io_global,ONLY : stdout
  USE wvfct,                ONLY : wg, npwx
  USE wvfct_gpum,           ONLY : wg_d
  USE becmod_subs_gpum,     ONLY : calbec_gpu
  USE noncollin_module,     ONLY : noncolin, npol
  IMPLICIT NONE
  !
  ! compute the (n,n) matrix representation <U|V>
  ! and energy from V (m,n) and U(m,n)
  !
  LOGICAL, INTENT(IN) :: DoE
  INTEGER, INTENT(IN) :: ik, ninner, n, m
  COMPLEX(dp), INTENT(IN), DEVICE :: U(ninner,n), V(ninner,m)
  COMPLEX(dp), INTENT(OUT), DEVICE:: mat(n,m)
  REAL(DP), INTENT(OUT) :: ee

  INTEGER :: i

  CALL start_clock('matcalc')

  mat = (0.0_dp, 0.0_dp)
  IF(noncolin) THEN
    noncolin = .false.
    CALL calbec_gpu(ninner, U, V, mat, m)
    noncolin = .true.
  ELSE
    CALL calbec_gpu(ninner, U, V, mat, m)
  ENDIF

  IF(DoE) THEN
    allocate(wg_d,source=wg)
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    ee = 0.0_dp
    !$cuf kernel do (1) 
    DO i = 1,n
      ee = ee + wg_d(i,ik)*DBLE(mat(i,i))
    ENDDO
    deallocate(wg_d)
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE matcalc_k_gpu
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MatPrt(label,n,m,A)
  USE kinds, ONLY : dp
  USE io_global,ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n,m
  REAL(dp), INTENT(IN):: A(n,m)
  INTEGER :: i
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(stdout,'(A)') label
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f16.10)'
  DO i = 1,n
    WRITE(stdout,frmt) A(i,:)
  ENDDO

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matprt_k(label,n,m,A)
  USE kinds, ONLY : dp
  USE io_global,ONLY : stdout
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n,m
  COMPLEX(dp), INTENT(IN):: A(n,m)
  INTEGER :: i
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(stdout,'(A)') label//'(real)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
     WRITE(stdout,frmt) dreal(A(i,:))
  ENDDO

  WRITE(stdout,'(A)') label//'(imag)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
     WRITE(stdout,frmt) aimag(A(i,:))
  ENDDO
END SUBROUTINE matprt_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MatInv(MShape,n,A)
  USE kinds, ONLY : dp
  IMPLICIT NONE
!
! given a real square matrix A, returns the inverse in A
! in the same shape as the input matrix, as indicated by MShape
!    MShape = 'L' ... A is Lower Triangular (allocated in square shape)
!             'U' ... A is Upper Triangular (allocated in square shape)
!             'G' ... A is a general matrix
!
  INTEGER, INTENT(IN) :: n
  REAL(dp), INTENT(INOUT):: A(n,n)
  CHARACTER (LEN=1) :: MShape
  INTEGER :: INFO, LWORK
  INTEGER, ALLOCATABLE :: IPIV(:)
  REAL(DP), ALLOCATABLE :: WORK(:)

  IF(MShape.eq.'L'.or.MShape.eq.'U') then 
    INFO = -1
    CALL DTRTRI( MShape, 'N', n, A, n, INFO )
    CALL errinfo('DTRTRI','inversion failed in MatInv.',INFO)
  ELSEIF(MShape.eq.'G') then 
    LWORK = 3*n
    ALLOCATE( IPIV(n), WORK(LWORK) ) 
    INFO = -1
    CALL DGETRF( n, n, A, n, IPIV, INFO )
    CALL errinfo('DGETRF','LU decomposition failed in MatInv.',INFO)
    INFO = -1
    CALL DGETRI( n, A, n, IPIV, WORK, LWORK, INFO )
    CALL errinfo('DGETRI','inversion failed in MatInv.',INFO)
    DEALLOCATE( IPIV, WORK ) 
  ELSE
    call errore('MatInv', 'Wrong MShape.', 1) 
  END IF 

END SUBROUTINE MatInv 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MatChol(n,A)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(dp), INTENT(INOUT):: A(n,n)
!
! given a (real, positive definite) matrix A, returns the Cholesky factor in A
! (only the Lower Triangular part of the input matrix is considered)
!
  integer :: INFO

  INFO = -1
  CALL DPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('DPOTRF','Cholesky failed in MatChol.',INFO)

END SUBROUTINE MatChol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE invchol(n,A)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(dp), INTENT(IN):: A(n,n)
!
! given a matrix A, returns the inverse of the Cholesky decomposition of A
! for real matrices
!
  integer :: INFO

  INFO = -1
  CALL DPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('DPOTRF','Cholesky failed in invchol.',INFO)
  INFO = -1
  CALL DTRTRI( 'L', 'N', n, A, n, INFO )
  CALL errinfo('DTRTRI','inversion failed in invchol.',INFO)

END SUBROUTINE invchol

SUBROUTINE invchol_k(n,A)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX(dp), INTENT(IN):: A(n,n)
  !
  ! given a matrix A, returns the inverse of the Cholesky decomposition of A
  ! for cholesky matrices
  !
  INTEGER :: INFO

  INFO = -1
  CALL ZPOTRF( 'L', n, A, n, INFO )
  CALL errinfo('ZPOTRF','Cholesky failed in invchol.',INFO)
  INFO = -1
  CALL ZTRTRI( 'L', 'N', n, A, n, INFO )
  CALL errinfo('ZTRTRI','inversion failed in invchol.',INFO)

END SUBROUTINE invchol_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE errinfo(routine,message,INFO)
  IMPLICIT NONE
  INTEGER :: INFO
  CHARACTER(len=*) :: routine,message

  IF(INFO/=0) THEN
     WRITE(*,*) routine,' exited with INFO= ',INFO
     CALL errore(routine,message,1)
  ENDIF

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MatSymm( MShape, How, Mat, n )
USE kinds, ONLY : dp
!
! Symmetrize the (square) matrix Mat
!       How    = 'U'  ... copying the upper block into the lower block
!                'L'  ... copying the lower block into the upper block
!                'S'  ... averaging 
!
!       MShape = 'U'  ... return the Upper Triangular (Zeros in Lower)
!                'L'  ... return the Lower Triangular (Zeros in Upper)
!                'S'  ... return the Square symmetric matrix
!
IMPLICIT NONE
  INTEGER :: n, i, j 
  REAL(DP) ::  Mat(n,n) 
  REAL(DP), ALLOCATABLE :: MatT(:,:)
  REAL(DP), PARAMETER :: Zero=0.0d0, Two=2.0d0
  CHARACTER(LEN=1) :: How, MShape

  ALLOCATE( MatT(n,n) )
 
! Properly fill the lower triangular of MatT
  MatT = Zero 
  IF(How.eq.'L') then ! use lower
    do i = 1, n
      MatT(i,i) = Mat(i,i)
      do j = i+1, n 
        MatT(j,i) = Mat(j,i)
      end do        
    end do        
  ELSE IF( How.eq.'U' ) then ! use upper
    do i = 1, n
      MatT(i,i) = Mat(i,i)
      do j = i+1, n
        MatT(j,i) = Mat(i,j)
      end do        
    end do        
  ELSE IF( How.eq.'S' ) then ! use average 
    do i = 1, n
      MatT(i,i) = Mat(i,i)
      do j = i+1, n
        MatT(j,i) = (Mat(i,j) + Mat(j,i))  / Two
      end do        
    end do        
  ELSE
    Call errore('MatSymm','Wrong How in MatSymm.',1)
  END IF 

! Properly copy the results in Mat
  Mat = Zero 
  IF(MShape.eq.'L') then ! return lower 
    Mat=MatT
  ELSE IF(MShape.eq.'U') then ! return upper 
    do i = 1, n
      Mat(i,i) = MatT(i,i)
      do j = i+1, n
        Mat(i,j) = MatT(j,i)   
      end do        
    end do        
  ELSE IF(MShape.eq.'S') then ! return square
    Mat=MatT
    do i = 1, n
      do j = i+1, n
        Mat(i,j) = MatT(j,i)   
      end do        
    end do        
  ELSE
    Call errore('MatSymm','Wrong MShape in MatSymm.',1)
  END IF 

  DEALLOCATE( MatT )

END SUBROUTINE MatSymm  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MatTrp( Mat, n )
USE kinds, ONLY : dp
!
! Transpose the (square) matrix Mat
!
IMPLICIT NONE
  INTEGER :: n, i, j 
  REAL(DP) ::  IJth, JIth
  REAL(DP), INTENT(INOUT) ::  Mat(n,n)

  do i = 1, n
    do j = i+1, n 
      IJth = Mat(i,j)
      JIth = Mat(j,i)
      Mat(i,j) = JIth 
      Mat(j,i) = IJth 
    end do        
  end do        

END SUBROUTINE MatTrp  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
