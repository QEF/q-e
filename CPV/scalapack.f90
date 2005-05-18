!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      MODULE scalapack

        USE kinds

        IMPLICIT NONE
        SAVE
 
        PRIVATE

        INTERFACE redistribute
          MODULE PROCEDURE redistribute_rmatrix, redistribute_cmatrix, &
            redistribute_imatrix
        END INTERFACE
        INTERFACE pmatmul
          MODULE PROCEDURE real_parallel_matmul
        END INTERFACE
        INTERFACE ptranspose
          MODULE PROCEDURE real_parallel_mtran
        END INTERFACE
        INTERFACE pdiagonalize
          MODULE PROCEDURE real_parallel_mdiag
        END INTERFACE

        PUBLIC :: pmatmul, ptranspose, pdiagonalize

      CONTAINS

        SUBROUTINE print_rmatrix(a,comment)
          USE parallel_types, ONLY: real_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor
          TYPE (real_parallel_matrix) :: a
          CHARACTER(len=*), INTENT(IN) :: comment
          REAL (dbl), ALLOCATABLE :: work(:)
          INTEGER :: desc_a(9)
          CALL pblas_descriptor(desc_a, a%desc)
          ALLOCATE(work(a%desc%nxblk))
#if defined __SCALAPACK
          CALL PDLAPRNT( a%desc%nx, a%desc%ny, A%m(1,1), 1, 1, desc_a, 0, 0, &
            comment, 6, WORK(1) )
#endif
          DEALLOCATE(work)
          RETURN
        END SUBROUTINE print_rmatrix


        SUBROUTINE redistribute_rmatrix(a,b)
          USE parallel_types, ONLY: real_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor
          TYPE (real_parallel_matrix) :: a, b
          INTEGER :: desc_a(9), desc_b(9)
          CALL pblas_descriptor(desc_a, a%desc)
          CALL pblas_descriptor(desc_b, b%desc)
#if defined __SCALAPACK
          CALL PDGEMR2D(a%desc%nx, a%desc%ny, a%m(1,1), 1, 1, desc_a, &
            b%m(1,1), 1, 1, desc_b, b%desc%grid%context)
#endif
          RETURN
        END SUBROUTINE redistribute_rmatrix

        SUBROUTINE redistribute_cmatrix(a,b)
          USE parallel_types, ONLY: complex_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor
          TYPE (complex_parallel_matrix) :: a, b
          INTEGER :: desc_a(9), desc_b(9)
          CALL pblas_descriptor(desc_a, a%desc)
          CALL pblas_descriptor(desc_b, b%desc)
#if defined __SCALAPACK
          CALL PZGEMR2D(a%desc%nx, a%desc%ny, a%m(1,1), 1, 1, desc_a, &
            b%m(1,1), 1, 1, desc_b, b%desc%grid%context)
#endif
          RETURN
        END SUBROUTINE redistribute_cmatrix

        SUBROUTINE redistribute_imatrix(a,b)
          USE parallel_types, ONLY: integer_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor
          TYPE (integer_parallel_matrix) :: a, b
          INTEGER :: desc_a(9), desc_b(9)
          CALL pblas_descriptor(desc_a, a%desc)
          CALL pblas_descriptor(desc_b, b%desc)
#if defined __SCALAPACK
          CALL PIGEMR2D(a%desc%nx, a%desc%ny, a%m(1,1), 1, 1, desc_a, &
            b%m(1,1), 1, 1, desc_b, b%desc%grid%context)
#endif
          RETURN
        END SUBROUTINE redistribute_imatrix

        SUBROUTINE real_parallel_matmul(a,b,c,transa,transb,eq_dist, &
          alphax,betax,iax,jax,ibx,jbx,icx,jcx)
          USE parallel_types, ONLY: real_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor
          TYPE (real_parallel_matrix) :: a, b, c
          CHARACTER(LEN=*), INTENT(IN) :: transa, transb
          LOGICAL, INTENT(IN), OPTIONAL :: eq_dist
          REAL(dbl), INTENT(IN), OPTIONAL :: alphax, betax
          INTEGER, INTENT(IN), OPTIONAL :: iax, jax, ibx, jbx, icx, jcx
          INTEGER :: desc_a(9), desc_b(9), desc_c(9)
          INTEGER :: ia, ja, ib, jb, ic, jc, k
          REAL(dbl) :: alpha, beta
          ia = 1; ja = 1; ib = 1; jb = 1; ic = 1; jc = 1
          alpha = 1.0_dbl; beta = 0.0_dbl
          IF(PRESENT(iax)) ia = iax
          IF(PRESENT(jax)) ja = jax
          IF(PRESENT(ibx)) ib = ibx
          IF(PRESENT(jbx)) jb = jbx
          IF(PRESENT(icx)) ic = icx
          IF(PRESENT(jcx)) jc = jcx
          IF(PRESENT(alphax)) alpha = alphax
          IF(PRESENT(betax)) beta = betax
          IF((transa(1:1) .EQ. 'T') .OR. (transa(1:1) .EQ. 'C') .OR. &
             (transa(1:1) .EQ. 't') .OR. (transa(1:1) .EQ. 'c') ) THEN
            k = a%desc%nx
          ELSE
            k = a%desc%ny
          END IF
          IF(PRESENT(eq_dist) .AND. eq_dist ) THEN
            CALL pblas_descriptor(desc_a, a%desc)
#if defined __SCALAPACK
            CALL pdgemm(transa, transb, c%desc%nx, c%desc%ny, k, &
              alpha, a%m(1,1), ia, ja, DESC_A, b%m(1,1), ib, jb, DESC_A, &
              beta, c%m(1,1), ic, jc, DESC_A)
#endif
          ELSE
            CALL pblas_descriptor(desc_a, a%desc)
            CALL pblas_descriptor(desc_b, b%desc)
            CALL pblas_descriptor(desc_c, b%desc)
#if defined __SCALAPACK
            CALL pdgemm(transa, transb, c%desc%nx, c%desc%ny, k, &
              alpha, a%m(1,1), ia, ja, DESC_A, b%m(1,1), ib, jb, DESC_B, &
              beta, c%m(1,1), ic, jc, DESC_C)
#endif
          END IF
          RETURN
        END SUBROUTINE real_parallel_matmul


        SUBROUTINE real_parallel_mtran(a,c,alphax,betax,iax,jax,icx,jcx)
          USE parallel_types, ONLY: real_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor

          TYPE (real_parallel_matrix) :: a, c
          REAL(dbl), INTENT(IN), OPTIONAL :: alphax, betax
          INTEGER, INTENT(IN), OPTIONAL :: iax, jax, icx, jcx
          INTEGER :: desc_a(9), desc_c(9)
          INTEGER :: ia, ja, ic, jc
          REAL(dbl) :: alpha, beta
          ia = 1; ja = 1; ic = 1; jc = 1
          alpha = 1.0_dbl; beta = 0.0_dbl
          IF(PRESENT(iax)) ia = iax
          IF(PRESENT(jax)) ja = jax
          IF(PRESENT(icx)) ic = icx
          IF(PRESENT(jcx)) jc = jcx
          IF(PRESENT(alphax)) alpha = alphax
          IF(PRESENT(betax)) beta = betax
          CALL pblas_descriptor(desc_a, a%desc)
          CALL pblas_descriptor(desc_c, c%desc)
#if defined __SCALAPACK
          CALL pdtran(c%desc%nx, c%desc%ny, &
              alpha, a%m(1,1), ia, ja, DESC_A, &
              beta, c%m(1,1), ic, jc, DESC_C)
#endif
          RETURN
        END SUBROUTINE real_parallel_mtran

        SUBROUTINE real_parallel_mdiag(uplo, a, w, z, iax, jax, izx, jzx)
          USE parallel_types, ONLY: real_parallel_matrix
          USE descriptors_module, ONLY: pblas_descriptor
          USE io_global, ONLY: stdout
          CHARACTER(LEN=*), INTENT(IN) :: uplo
          TYPE (real_parallel_matrix) :: a
          REAL(dbl), INTENT(IN) :: w(:)
          TYPE (real_parallel_matrix), OPTIONAL :: z
          INTEGER, INTENT(IN), OPTIONAL :: iax, jax, izx, jzx
          INTEGER :: desc_a(9), desc_z(9)
          INTEGER :: ia, ja, iz, jz, idum
          INTEGER :: lwork, info, mw, nz, liwork, iwork_size
          INTEGER, ALLOCATABLE :: ifail(:), iclustr(:)
          INTEGER, ALLOCATABLE :: iwork(:)
          REAL(dbl), ALLOCATABLE :: rwork(:)
          REAL(dbl) :: zdum, rdum, rwork_size, dum
          REAL(dbl), ALLOCATABLE :: gap(:) 
          CHARACTER(LEN=1) :: jobz

          ia = 1; ja = 1; iz = 1; jz = 1
          IF( PRESENT(z) ) THEN
            jobz = 'V'
            CALL pblas_descriptor(desc_a, a%desc)
            CALL pblas_descriptor(desc_z, z%desc)
          ELSE
            jobz = 'N'
            CALL pblas_descriptor(desc_a, a%desc)
            desc_z = desc_a
          END IF
          IF(PRESENT(iax)) ia = iax
          IF(PRESENT(jax)) ja = jax
          IF(PRESENT(izx)) iz = izx
          IF(PRESENT(jzx)) jz = jzx

          lwork      =    -1; liwork     = -1; info = 0
          mw         =     0; nz         =  0; idum = 0; dum = 0.0d0;
          rwork_size = 0.0d0; iwork_size =  0

          ALLOCATE(ifail(a%desc%nx))
          ALLOCATE(iclustr(2 * a%desc%grid%npx * a%desc%grid%npy ) )
          ALLOCATE(gap( a%desc%grid%npx * a%desc%grid%npy ) )

#if defined __SCALAPACK
!          CALL PDSYEV( jobz, uplo, a%desc%nx, a%m(1,1), ia, ja, desc_a, &
!            w(1), zdum, iz, jz, desc_z, rwork_size, lwork, info )

          CALL PDSYEVX(JOBZ, 'A', UPLO, a%desc%nx, a%m(1,1), ia, ja, &
            desc_a, DUM, DUM, IDUM, IDUM, 0.0d0, MW, NZ, w(1), 0.0d0, &
            zdum, iz, jz, desc_z, rwork_size, lwork, iwork_size, &
            liwork, IFAIL, ICLUSTR, GAP, INFO )

          IF(INFO.NE.0) THEN
            WRITE( stdout,101) INFO, IFAIL
 101        FORMAT (' *** ERROR IN ROUTINE real_parallel_mdiag: 1 ***  ', 2I5)
            STOP
          END IF

#endif

          lwork = 2 * MAX(1,INT(rwork_size + 1))
          ALLOCATE( rwork(lwork) )
          liwork = 2 * MAX(1,iwork_size)
          ALLOCATE( iwork(liwork) )

          IF(PRESENT (z)) THEN
#if defined __SCALAPACK
!            CALL PDSYEV( jobz, uplo, a%desc%nx, a%m(1,1), ia, ja, desc_a, &
!              w(1), z%m(1,1), iz, jz, desc_z, rwork(1), lwork, info )

            CALL PDSYEVX(JOBZ, 'A', UPLO, a%desc%nx, a%m(1,1), ia, ja, &
              desc_a, DUM, DUM, IDUM, IDUM, 0.0d0, MW, NZ, w(1), 0.0d0, &
              z%m(1,1), iz, jz, desc_z, rwork(1), lwork, iwork(1), &
              liwork, IFAIL, ICLUSTR, GAP, INFO )

#endif
          ELSE
#if defined __SCALAPACK
!            CALL PDSYEV( jobz, uplo, a%desc%nx, a%m(1,1), ia, ja, desc_a, &
!              w(1), zdum, iz, jz, desc_z, rwork(1), lwork, info )

            CALL PDSYEVX(JOBZ, 'A', UPLO, a%desc%nx, a%m(1,1), ia, ja, &
              desc_a, DUM, DUM, IDUM, IDUM, 0.0d0, MW, NZ, w(1), 0.0d0, &
              zdum, iz, jz, desc_z, rwork(1), lwork, iwork(1), &
              liwork, IFAIL, ICLUSTR, GAP, INFO )

#endif
          END IF
          IF(INFO.NE.0) THEN
            WRITE( stdout,102) INFO, IFAIL
 102        FORMAT (' *** ERROR IN ROUTINE real_parallel_mdiag: 2 ***  ', 2I5)
            STOP
          END IF
          DEALLOCATE(rwork)
          DEALLOCATE(iwork)
          DEALLOCATE(ifail, iclustr, gap)
          RETURN
        END SUBROUTINE real_parallel_mdiag

      END MODULE scalapack
