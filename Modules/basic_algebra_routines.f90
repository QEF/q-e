!
! Copyright (C) 2003-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE basic_algebra_routines
  !----------------------------------------------------------------------------
  !! This module contains a limited number of functions and operators
  !! for vectorial algebra. Wherever possible the appropriate BLAS routine
  !! (always the double precision version) is used.
  !
  !! List of public methods :
  !
  !! * \(x.\text{dot}.y\): dot product between vectors (\(\langle x|y\rangle\));
  !! * \(x.\text{ext}.y\): external (vector) product between vectors 
  !!                       (\(\langle x|y\rangle\));
  !! * \(\text{norm}(x)\): norm of a vector ( \(\sqrt{\langle x|x \rangle}\) );
  !! * \(A.\text{times}.x\): matrix-vector multiplication ( \(A|x\rangle\) );
  !! * \(x.\text{times}.A\): vector-matrix multiplication ( \(\langle x|A\) );
  !! * \(\text{matrix}(x,y)\): vector-vector multiplication ( |x\langle\rangle y| );
  !! * \(\text{identity}(N)\): identity matrix of rank N.
  !
  !! Written by Carlo Sbraccia (16/12/2003)
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE OPERATOR( .dot. )
     !
     MODULE PROCEDURE dot_product_
     !
  END INTERFACE
  !
  INTERFACE OPERATOR( .ext. )
     !
     MODULE PROCEDURE external_product_
     !
  END INTERFACE
  !  
  INTERFACE OPERATOR( .times. )
     !
     MODULE PROCEDURE matrix_times_vector, vector_times_matrix
     !
  END INTERFACE
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     FUNCTION dot_product_( vec1, vec2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vec1(:), vec2(:)
       REAL(DP)             :: dot_product_
       !
       REAL(DP) :: ddot
       EXTERNAL    ddot
       !
       dot_product_ = ddot( SIZE( vec1 ), vec1, 1, vec2, 1 )
       !
     END FUNCTION dot_product_
     !
     !-----------------------------------------------------------------------
     FUNCTION external_product_( vec1, vec2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vec1(:), vec2(:)
       REAL(DP)             :: external_product_(SIZE( vec1 ))
       !
       !
       external_product_(1) = + vec1(2)*vec2(3) - vec1(3)*vec2(2)
       external_product_(2) = - vec1(1)*vec2(3) + vec1(3)*vec2(1)
       external_product_(3) = + vec1(1)*vec2(2) - vec1(2)*vec2(1)
       !
     END FUNCTION external_product_
     !     
     !----------------------------------------------------------------------- 
     FUNCTION norm( vec )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vec(:)
       REAL(DP)             :: norm
       !
       REAL(DP) :: dnrm2
       EXTERNAL    dnrm2   
       !
       norm = dnrm2( SIZE( vec ), vec, 1 )
       !
     END FUNCTION norm
     !
     !-----------------------------------------------------------------------
     FUNCTION matrix_times_vector( mat, vec )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vec(:)
       REAL(DP), INTENT(IN) :: mat(:,:)
       REAL(DP)             :: matrix_times_vector(SIZE( vec ))
! gfortran hack
       REAL(DP)             :: aux(SIZE( vec ))
       INTEGER              :: dim1
       !
       dim1 = SIZE( vec )
       !
       CALL DGEMV( 'N', dim1, dim1, 1.0_DP, mat, dim1, vec, 1, 0.0_DP, &
                   aux, 1 ) 
! gfortran hack 
!                  matrix_times_vector, 1 )
       matrix_times_vector = aux
       !
     END FUNCTION  matrix_times_vector
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION vector_times_matrix( vec, mat )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vec(:)
       REAL(DP), INTENT(IN) :: mat(:,:)
       REAL(DP)             :: vector_times_matrix(SIZE( vec ))
! gfortran hack
       REAL(DP)             :: aux(SIZE( vec ))
       INTEGER              :: dim1
       !
       dim1 = SIZE( vec )
       !
       CALL DGEMV( 'T', dim1, dim1, 1.0_DP, mat, dim1, vec, 1, 0.0_DP, &
                   aux, 1) 
! gfortran hack 
!                  vector_times_matrix, 1 )
       vector_times_matrix = aux
       !
     END FUNCTION vector_times_matrix
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION matrix( vec1, vec2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vec1(:), vec2(:)
       REAL(DP)             :: matrix(SIZE( vec1 ),SIZE( vec2 ))
       INTEGER              :: dim1, dim2
       !
       dim1 = SIZE( vec1 )
       dim2 = SIZE( vec2 )
       !
       matrix = 0.0_DP
       CALL DGER( dim1, dim2, 1.0_DP, vec1, 1, vec2, 1, matrix, dim1 )
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
       REAL(DP)            :: identity(dim,dim)
       INTEGER             :: i
       !
       identity = 0.0_DP
       !
       FORALL( i = 1:dim ) identity(i,i) = 1.0_DP
       !
     END FUNCTION identity
     !    
END MODULE basic_algebra_routines
