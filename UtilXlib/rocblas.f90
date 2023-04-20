! Copyright (C) 2022 Advanced Micro Devices, Inc. All Rights Reserved.  
!  
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!  
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!  
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#if defined(__ROCBLAS)

#ifdef rocblas_ILP64
#define rocblas_int c_int64_t
#else
#define rocblas_int c_int32_t
#endif
    
MODULE rocblas

    USE ISO_C_BINDING
    USE util_param, ONLY : DP
    IMPLICIT NONE

    TYPE(C_PTR) :: handle
    INTEGER, PARAMETER :: ROCBLAS_STATUS_SUCCESS = 0
    LOGICAL :: isactive = .FALSE. 

    INTERFACE    

    FUNCTION rocblas_create_handle(handle) BIND(C,NAME="rocblas_create_handle")
    USE ISO_C_BINDING
    IMPLICIT NONE
    TYPE(C_PTR) :: handle
    INTEGER :: rocblas_create_handle
    END FUNCTION rocblas_create_handle

    FUNCTION rocblas_destroy_handle(handle) BIND(C,NAME="rocblas_destroy_handle")
    use iso_c_binding
    IMPLICIT NONE
    TYPE(C_PTR), VALUE :: handle
    INTEGER :: rocblas_destroy_handle
    END FUNCTION rocblas_destroy_handle

    END INTERFACE

    CONTAINS

        FUNCTION rocblas_get_operation(op)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: op
            INTEGER :: rocblas_get_operation
            SELECT CASE(OP)
            CASE ('N')
                rocblas_get_operation = 111
            CASE ('T')
                rocblas_get_operation = 112
            CASE ('C')
                rocblas_get_operation = 113
            END SELECT
        END FUNCTION

        SUBROUTINE rocblas_print_status(stat)
            IMPLICIT NONE
            INTEGER :: STAT

            SELECT CASE (STAT)
            CASE (0)
               PRINT *, "rocBLAS: success"
            CASE (1)
               PRINT *, "rocBLAS: invalid_handle"
            CASE (2)
               PRINT *, "rocBLAS: not implemented" 
            CASE (3)
               PRINT *, "rocBLAS: invalid pointer" 
            CASE (4)
               PRINT *, "rocBLAS: invalid size" 
            CASE (5)
               PRINT *, "rocBLAS: memory error" 
            CASE (6)
               PRINT *, "rocBLAS: internal error" 
            CASE (7)
               PRINT *, "rocBLAS: perf degraded" 
            CASE (8)
               PRINT *, "rocBLAS: query mismatch" 
            CASE (9)
               PRINT *, "rocBLAS: size increased" 
            CASE (10)
               PRINT *, "rocBLAS: size unchanged" 
            CASE (11)
               PRINT *, "rocBLAS: invalid value" 
            CASE (12)
               PRINT *, "rocBLAS: continue" 
            CASE (13)
               PRINT *, "rocBLAS: check numerics fail" 
            END SELECT
        END SUBROUTINE

        SUBROUTINE rocblas_init()
            IMPLICIT NONE
            INTEGER :: stat

            stat = rocblas_create_handle(handle)

            IF (isactive) THEN
                print *,"WARNING:: rocBLAS has already been initialized."
                RETURN
            END IF

            IF (stat /= ROCBLAS_STATUS_SUCCESS) THEN
                print *, "ERROR:: rocblas_init: "
                call rocblas_print_status(stat)
                STOP
            END IF

            isactive = .TRUE.

        END SUBROUTINE

        SUBROUTINE rocblas_destroy()
            IMPLICIT NONE
            INTEGER :: stat

            stat = rocblas_destroy_handle(handle)

            IF (stat /= ROCBLAS_STATUS_SUCCESS) THEN
                PRINT *, "ERROR:: rocblas_destroy: "
                CALL rocblas_print_status(stat)
            END IF

            isactive = .FALSE.
        END SUBROUTINE

        SUBROUTINE rocblas_check(stat, msg)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: stat
            CHARACTER(len=*), INTENT(IN) :: msg

            IF (stat /= ROCBLAS_STATUS_SUCCESS) THEN
                write (0, *) "rocBLAS Error: ", msg
                CALL rocblas_print_status(stat)
            END IF; 

        END SUBROUTINE


END MODULE rocblas

MODULE rocblas_utils
    USE rocblas, only : handle, rocblas_check, ROCBLAS_STATUS_SUCCESS, rocblas_get_operation, &
                        rocblas_print_status
    USE util_param, only : DP
    USE iso_c_binding
    IMPLICIT NONE

    INTERFACE    
        FUNCTION rocblas_dgemm_(handle, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
                               BIND(C, NAME="rocblas_dgemm") 
                USE ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR), VALUE :: handle
                INTEGER(C_INT), VALUE :: transA
                INTEGER(C_INT), VALUE :: transB
                INTEGER(rocblas_int), VALUE :: m, n, k
                REAL(c_double) :: alpha
                TYPE(C_PTR), VALUE :: A
                INTEGER(rocblas_int), VALUE :: lda
                TYPE(C_PTR), VALUE :: B
                INTEGER(rocblas_int), VALUE :: ldb
                REAL(c_double) :: beta
                TYPE(C_PTR), VALUE :: C
                INTEGER(rocblas_int), VALUE :: ldc
                INTEGER :: rocblas_dgemm_
        END FUNCTION rocblas_dgemm_

        FUNCTION rocblas_zgemm_(handle, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
                               BIND(C, NAME="rocblas_zgemm") 
                USE ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR), VALUE :: handle
                INTEGER(C_INT), VALUE :: transA
                INTEGER(C_INT), VALUE :: transB
                INTEGER(rocblas_int), VALUE :: m, n, k
                COMPLEX(c_double_complex) :: alpha
                TYPE(C_PTR), VALUE :: A
                INTEGER(rocblas_int), VALUE :: lda
                TYPE(C_PTR), VALUE :: B
                INTEGER(rocblas_int), VALUE :: ldb
                COMPLEX(c_double_complex) :: beta
                TYPE(C_PTR), VALUE :: C
                INTEGER(rocblas_int), VALUE :: ldc
                INTEGER :: rocblas_zgemm_
        END FUNCTION rocblas_zgemm_
    END INTERFACE

    INTERFACE rocblas_dgemm
        MODULE PROCEDURE rocblas_ddgemm, rocblas_dzgemm1, rocblas_dzgemm2, rocblas_dzgemm3    
    END INTERFACE
 
    INTERFACE 
        FUNCTION rocblas_dger_(handle, m, n, alpha, x, incx, y, incy, A, lda) &
                               BIND(C, NAME="rocblas_dger") 
                USE ISO_C_BINDING
                IMPLICIT NONE
                TYPE(C_PTR), VALUE :: handle
                INTEGER(rocblas_int), VALUE :: m, n
                REAL(c_double) :: alpha
                TYPE(C_PTR), VALUE :: x
                INTEGER(rocblas_int), VALUE :: incx
                TYPE(C_PTR), VALUE :: y
                INTEGER(rocblas_int), VALUE :: incy
                TYPE(C_PTR), VALUE :: A
                INTEGER(rocblas_int), VALUE :: lda
                INTEGER :: rocblas_dger_
        END FUNCTION rocblas_dger_

    END INTERFACE

    INTERFACE rocblas_dger
        MODULE PROCEDURE rocblas_dger1, rocblas_dzger
    END INTERFACE  

    CONTAINS

        SUBROUTINE rocblas_zgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: transA, transB
            INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
            COMPLEX(DP), INTENT(IN) :: alpha, beta
            COMPLEX(DP), INTENT(IN), TARGET :: A(lda,*)
            COMPLEX(DP), INTENT(IN), TARGET :: B(ldb,*)
            COMPLEX(DP), INTENT(INOUT), TARGET :: C(ldc,*)
            INTEGER :: stat
            INTEGER(c_int) :: itransA, itransB
            INTEGER(rocblas_int) :: rm, rn, rk, rlda, rldb, rldc

            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rk = int(k, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))
            rldb = int(ldb, kind(rocblas_int))
            rldc = int(ldc, kind(rocblas_int))
            itransA = rocblas_get_operation(transA)
            itransB = rocblas_get_operation(transB)

            !$omp target data use_device_ptr(A, B, C)
            stat = rocblas_zgemm_(handle, itransA, itransB,    &
                                 rm, rn, rk,                   &
                                 alpha, c_loc(A), rlda,        &
                                 c_loc(B), rldb, beta,         &
                                 c_loc(C), rldc)
            !$omp end target data

            CALL rocblas_check(stat, "zgemm")

        END SUBROUTINE

        SUBROUTINE rocblas_ddgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: transA, transB
            INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
            REAL(DP), INTENT(IN) :: alpha, beta
            REAL(DP), INTENT(IN), TARGET :: A(lda,*)
            REAL(DP), INTENT(IN), TARGET :: B(ldb,*)
            REAL(DP), INTENT(INOUT), TARGET :: C(ldc,*)
            INTEGER :: stat
            INTEGER(c_int) :: itransA, itransB
            INTEGER :: rm, rn, rk, rlda, rldb, rldc

            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rk = int(k, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))
            rldb = int(ldb, kind(rocblas_int))
            rldc = int(ldc, kind(rocblas_int))
            itransA = rocblas_get_operation(transA)
            itransB = rocblas_get_operation(transB)

            !$omp target data use_device_ptr(A, B, C)
            stat = rocblas_dgemm_(handle, itransA, itransB,    &
                                 rm, rn, rk,                   &
                                 alpha, c_loc(A), rlda,        &
                                 c_loc(B), rldb, beta,         &
                                 c_loc(C), rldc);
            !$omp end target data
            
            CALL rocblas_check(stat, "ddgemm")

        END SUBROUTINE

        SUBROUTINE rocblas_dzgemm(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: transA, transB
            INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
            REAL(DP), INTENT(IN) :: alpha, beta
            COMPLEX(DP), INTENT(IN), TARGET :: A(lda,*)
            COMPLEX(DP), INTENT(IN), TARGET :: B(ldb,*)
            COMPLEX(DP), INTENT(INOUT), TARGET :: C(ldc,*)
            INTEGER :: stat
            INTEGER(c_int) :: itransA, itransB
            INTEGER :: rm, rn, rk, rlda, rldb, rldc

            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rk = int(k, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))
            rldb = int(ldb, kind(rocblas_int))
            rldc = int(ldc, kind(rocblas_int))
            itransA = rocblas_get_operation(transA)
            itransB = rocblas_get_operation(transB)

            !$omp target data use_device_ptr(A, B, C)
            stat = rocblas_dgemm_(handle, itransA, itransB,  &
                                 rm, rn, rk,                 &
                                 alpha, c_loc(A), rlda,      &
                                 c_loc(B), rldb, beta,       &
                                 c_loc(C), rldc);            
            !$omp end target data

            CALL rocblas_check(stat, "dzgemm")

        END SUBROUTINE

        SUBROUTINE rocblas_dzgemm1(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: transA, transB
            INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
            REAL(DP), INTENT(IN) :: alpha, beta
            REAL(DP), INTENT(IN), TARGET :: A(lda,*)
            COMPLEX(DP), INTENT(IN), TARGET :: B(ldb,*)
            COMPLEX(DP), INTENT(INOUT), TARGET :: C(ldc,*)
            INTEGER :: stat
            INTEGER(c_int) :: itransA, itransB
            INTEGER :: rm, rn, rk, rlda, rldb, rldc

            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rk = int(k, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))
            rldb = int(ldb, kind(rocblas_int))
            rldc = int(ldc, kind(rocblas_int))
            itransA = rocblas_get_operation(transA)
            itransB = rocblas_get_operation(transB)

            !$omp target data use_device_ptr(A, B, C)
            stat = rocblas_dgemm_(handle, itransA, itransB,  &
                                 rm, rn, rk,                 &
                                 alpha, c_loc(A), rlda,      &
                                 c_loc(B), rldb, beta,       &
                                 c_loc(C), rldc);
            !$omp end target data

            CALL rocblas_check(stat, "dzgemm1")

        END SUBROUTINE

        SUBROUTINE rocblas_dzgemm2(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: transA, transB
            INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
            REAL(DP), INTENT(IN) :: alpha, beta
            COMPLEX(DP), INTENT(IN), TARGET :: A(lda,*)
            REAL(DP), INTENT(IN), TARGET :: B(ldb,*)
            COMPLEX(DP), INTENT(INOUT), TARGET :: C(ldc,*)
            INTEGER :: stat
            INTEGER(c_int) :: itransA, itransB
            INTEGER :: rm, rn, rk, rlda, rldb, rldc

            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rk = int(k, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))
            rldb = int(ldb, kind(rocblas_int))
            rldc = int(ldc, kind(rocblas_int))
            itransA = rocblas_get_operation(transA)
            itransB = rocblas_get_operation(transB)

            !$omp target data use_device_ptr(A, B, C)
            stat = rocblas_dgemm_(handle, itransA, itransB,  &
                                 rm, rn, rk,                 &
                                 alpha, c_loc(A), rlda,      &
                                 c_loc(B), rldb, beta,       &
                                 c_loc(C), rldc);
            !$omp end target data

            CALL rocblas_check(stat, "dzgemm2")

        END SUBROUTINE

        SUBROUTINE rocblas_dzgemm3(transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            IMPLICIT NONE
            CHARACTER, INTENT(IN) :: transA, transB
            INTEGER, INTENT(IN) :: m, n, k, lda, ldb, ldc
            REAL(DP), INTENT(IN) :: alpha, beta
            COMPLEX(DP), INTENT(IN), TARGET :: A(lda,*)
            COMPLEX(DP), INTENT(IN), TARGET :: B(ldb,*)
            REAL(DP), INTENT(INOUT), TARGET :: C(ldc,*)
            INTEGER :: stat
            INTEGER(c_int) :: itransA, itransB
            INTEGER :: rm, rn, rk, rlda, rldb, rldc

            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rk = int(k, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))
            rldb = int(ldb, kind(rocblas_int))
            rldc = int(ldc, kind(rocblas_int))
            itransA = rocblas_get_operation(transA)
            itransB = rocblas_get_operation(transB)

            !$omp target data use_device_ptr(A, B, C)
            stat = rocblas_dgemm_(handle, itransA, itransB,  &
                                 rm, rn, rk,                 &
                                 alpha, c_loc(A), rlda,      &
                                 c_loc(B), rldb, beta,       &
                                 c_loc(C), rldc);
            !$omp end target data

            CALL rocblas_check(stat, "dzgemm3")

        END SUBROUTINE

        SUBROUTINE rocblas_dger1(m, n, alpha, x, incx, y, incy, A, lda) 
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: m, n, incx, incy, lda
            REAL(DP), INTENT(IN) :: alpha 
            REAL(DP), INTENT(IN), TARGET :: x(m), y(n)
            REAL(DP), INTENT(INOUT), TARGET :: A(lda,*)
            INTEGER :: rm, rn, rincx, rincy, rlda
            INTEGER :: stat
            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rincx = int(incx, kind(rocblas_int))
            rincy = int(incx, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))

            !$omp target data use_device_ptr(A, x, y)
            stat = rocblas_dger_(handle, rm, rn, alpha, c_loc(x), rincx, c_loc(y), rincy, &
                                 c_loc(A), rlda)
            !$omp end target data
            CALL rocblas_check(stat, "DGER1")

        END SUBROUTINE

        SUBROUTINE rocblas_dzger(m, n, alpha, x, incx, y, incy, A, lda) 
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: m, n, incx, incy, lda
            REAL(DP), INTENT(IN) :: alpha 
            COMPLEX(DP), INTENT(IN) :: x(m), y(n)
            REAL(DP), INTENT(INOUT) :: A(lda,*)
            INTEGER :: rm, rn, rincx, rincy, rlda
            INTEGER :: stat
            rm = int(m, kind(rocblas_int))
            rn = int(n, kind(rocblas_int))
            rincx = int(incx, kind(rocblas_int))
            rincy = int(incx, kind(rocblas_int))
            rlda = int(lda, kind(rocblas_int))

            !$omp target data use_device_ptr(A, x, y)
            stat = rocblas_dger_(handle, rm, rn, alpha, c_loc(x), rincx, c_loc(y), rincy, &
                                 c_loc(A), rlda)
            !$omp end target data
            CALL rocblas_check(stat, "DZGER")

        END SUBROUTINE

END MODULE rocblas_utils
#endif
