!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE basic_algebra_routines
  !----------------------------------------------------------------------------
  !
  USE parameters,  ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE OPERATOR( .dot. )
     !
     MODULE PROCEDURE internal_dot_product
     !
  END INTERFACE
  !  
  INTERFACE OPERATOR( * )
     !
     MODULE PROCEDURE matrix_times_vector, vector_times_matrix
     !
  END INTERFACE
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     FUNCTION internal_dot_product( vector1, vector2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector1(:), vector2(:)
       REAL (KIND=DP)             :: internal_dot_product
       REAL (KIND=DP)             :: DDOT
       !
       EXTERNAL   DDOT
       !
       !
       internal_dot_product = DDOT( SIZE( vector1 ), vector1, 1, vector2, 1 )
       !
     END FUNCTION internal_dot_product
     !
     !     
     !----------------------------------------------------------------------- 
     FUNCTION norm( vector )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector(:)
       REAL (KIND=DP)             :: norm
       REAL (KIND=DP)             :: DNRM2
       !
       EXTERNAL   DNRM2
       !
       !
       norm = DNRM2( SIZE( vector ), vector, 1 )
       !
     END FUNCTION norm
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION matrix_times_vector( matrix , vector )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector(:)
       REAL (KIND=DP), INTENT(IN) :: matrix(:,:)
       REAL (KIND=DP)             :: matrix_times_vector(SIZE( vector ))
       INTEGER                    :: dim
       !
       EXTERNAL   DGEMV
       !
       !
       dim = SIZE( vector )
       !
       CALL DGEMV( 'N', dim, dim, 1.D0, matrix, dim, vector, 1, &
                   0.D0, matrix_times_vector, 1 )
       !
     END FUNCTION  matrix_times_vector
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION vector_times_matrix( vector , matrix )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector(:)
       REAL (KIND=DP), INTENT(IN) :: matrix(:,:)
       REAL (KIND=DP)             :: vector_times_matrix(SIZE( vector ))
       INTEGER                    :: i, dim
       !
       !
       dim = SIZE( vector )
       !
       DO i = 1, dim
          !
          vector_times_matrix(i) = vector(:) .dot. matrix(:,i)
          !
       END DO
       !
     END FUNCTION vector_times_matrix
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION matrix( vector1 , vector2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector1(:), vector2(:)
       REAL (KIND=DP)             :: matrix(SIZE( vector1 ),SIZE( vector2 ))
       INTEGER                    :: dim1, dim2
       !
       EXTERNAL   DGER
       !
       !
       dim1 = SIZE( vector1 )
       dim2 = SIZE( vector2 )
       !
       CALL DGER( dim1, dim2, 1.D0, &
                  vector1, 1, vector2, 1, matrix, MAX( dim1, dim2 )  )
       !
     END FUNCTION matrix
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION identity( dim )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: dim
       REAL(KIND=DP)       :: identity(dim,dim)
       INTEGER             :: i
       !
       !
       identity = 0.D0
       !
       DO i = 1, dim
          !
          identity(i,i) = 1.D0
          !
       END DO
       !
     END FUNCTION identity
     !    
END MODULE basic_algebra_routines
