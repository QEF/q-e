MODULE additional_cusolver_subs
   USE kinds,     ONLY : DP
   USE mp_bands,  ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
   USE mp,        ONLY : mp_bcast
   USE io_global, ONLY : stdout
#if defined(__CUDA)
   USE cudafor
   USE cusolverdn
#endif

   PRIVATE
   PUBLIC init_cusolver_handle, rdiagh_gpu, cdiagh_gpu
#if defined(__CUDA)
   PUBLIC cuSolverHandle, cuSolverInitialized
#endif

#if defined(__CUDA)
   TYPE(cusolverDnHandle), SAVE :: cuSolverHandle
   LOGICAL, SAVE :: cuSolverInitialized = .FALSE.
#endif

   CONTAINS
      SUBROUTINE init_cusolver_handle
         IMPLICIT NONE
         INTEGER :: info
#if defined(__CUDA)
         IF (cuSolverInitialized) RETURN

         info = cusolverDnCreate(cuSolverHandle)
         IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL errore( 'nao', 'cusolverDnCreate',  ABS(info))
         cuSolverInitialized = .TRUE.
#endif
      END SUBROUTINE init_cusolver_handle

      SUBROUTINE rdiagh_gpu(n, h_d, ldh, e_d, v_d)
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: n, ldh
         REAL(DP), INTENT(INOUT) :: h_d(ldh,n), e_d(n), v_d(ldh,n)
#if defined(__CUDA)

         INTEGER :: info, devInfo, lwork
         INTEGER, DEVICE :: devInfo_d
         REAL(DP), ALLOCATABLE :: work_d(:)
         ATTRIBUTES(DEVICE) :: h_d, e_d, v_d, work_d

         IF (me_bgrp == root_bgrp) THEN
            CALL init_cusolver_handle

            v_d(1:n,1:n) = h_d(1:n,1:n)
            info = cusolverDnDsyevd_bufferSize(&
               cuSolverHandle,  CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,&
               n, v_d, ldh,  e_d, lwork)
            IF (info /= CUSOLVER_STATUS_SUCCESS) THEN
               CALL errore("rdiagh_gpu", "cusolverDnDsyevd_bufferSize fail", ABS(info))
            END IF

            ALLOCATE(work_d(lwork))

            info = cusolverDnDsyevd(&
               cuSolverHandle,  CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,&
               n, v_d, ldh, e_d, work_d, lwork, devInfo_d)
            IF (info /= CUSOLVER_STATUS_SUCCESS) THEN
               WRITE(stdout, "('ERROR: cusolverDnDsyevd info: ', I0)") info
               CALL errore("rdiagh_gpu", "cusolverDnDsyevd fail", ABS(info))
            END IF
            devInfo = devInfo_d
            IF (devInfo /= 0) THEN
               WRITE(stdout, "('ERROR: cusolverDnDsyevd devInfo: ', I0)") devInfo
               CALL errore("rdiagh_gpu", "cusolverDnDsyevd fail", ABS(devInfo))
            END IF

            DEALLOCATE(work_d)
         END IF
         CALL mp_bcast(e_d, root_bgrp, intra_bgrp_comm)
         CALL mp_bcast(v_d, root_bgrp, intra_bgrp_comm)
#else
         CALL errore("rdiagh_gpu", "called gpu function without gpu", 1)
#endif
      END SUBROUTINE rdiagh_gpu

      SUBROUTINE cdiagh_gpu(n, h_d, ldh, e_d, v_d)
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: n, ldh
         REAL(DP),    INTENT(INOUT) :: e_d(n)
         COMPLEX(DP), INTENT(INOUT) :: h_d(ldh,n), v_d(ldh,n)
#if defined(__CUDA)

         INTEGER :: info, devInfo, lwork
         INTEGER, DEVICE :: devInfo_d
         COMPLEX(DP), ALLOCATABLE :: work_d(:)
         ATTRIBUTES(DEVICE) :: h_d, e_d, v_d, work_d

         IF (me_bgrp == root_bgrp) THEN
            CALL init_cusolver_handle

            v_d(1:n,1:n) = h_d(1:n,1:n)
            info = cusolverDnZheevd_bufferSize(&
               cuSolverHandle,  CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,&
               n, v_d, ldh,  e_d, lwork)
            IF (info /= CUSOLVER_STATUS_SUCCESS) THEN
               CALL errore("cdiagh_gpu", "cusolverDnZheevd_bufferSize fail", ABS(info))
            END IF

            ALLOCATE(work_d(lwork))

            info = cusolverDnZheevd(&
               cuSolverHandle,  CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER,&
               n, v_d, ldh, e_d, work_d, lwork, devInfo_d)
            IF (info /= CUSOLVER_STATUS_SUCCESS) THEN
               WRITE(stdout, "('ERROR: cusolverDnDsyevd info: ', I0)") info
               CALL errore("cdiagh_gpu", "cusolverDnZheevd fail", ABS(info))
            END IF
            devInfo = devInfo_d
            IF (devInfo /= 0) THEN
               WRITE(stdout, "('ERROR: cusolverDnDsyevd devInfo: ', I0)") devInfo
               CALL errore("cdiagh_gpu", "cusolverDnZheevd fail", ABS(devInfo))
            END IF

            DEALLOCATE(work_d)
         END IF
         CALL mp_bcast(e_d, root_bgrp, intra_bgrp_comm)
         CALL mp_bcast(v_d, root_bgrp, intra_bgrp_comm)
#else
         CALL errore("cdiagh_gpu", "called gpu function without gpu", 1)
#endif
      END SUBROUTINE cdiagh_gpu
END MODULE additional_cusolver_subs
