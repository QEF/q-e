#if !defined(__CUDA)
#define cublasZgemv zgemv
#define cublasZcopy zcopy
#endif

MODULE oscdft_wavefunction_subs_gpu
#if defined (__OSCDFT)
   USE kinds, ONLY : DP
   USE oscdft_context, ONLY : oscdft_context_type
   USE oscdft_indices, ONLY : oscdft_indices_type,&
                              oscdft_constr_indices_type,&
                              oscdft_orbital_indices_type
   USE oscdft_input,   ONLY : oscdft_input_type

   PRIVATE
   PUBLIC oscdft_ortho_swfc_gpu, oscdft_get_overlap_gpu, oscdft_lowdin_ortho_overlap_gpu

   INTERFACE oscdft_get_overlap_gpu
      MODULE PROCEDURE get_overlap_gpu_ctx, get_overlap_gpu_internal
   END INTERFACE oscdft_get_overlap_gpu
   INTERFACE oscdft_lowdin_ortho_overlap_gpu
      MODULE PROCEDURE lowdin_ortho_overlap_gpu_ctx, lowdin_ortho_overlap_gpu_internal
   END INTERFACE oscdft_lowdin_ortho_overlap_gpu

   CONTAINS
      SUBROUTINE oscdft_ortho_swfc_gpu(npwx, npw, m, wfc, swfc, normalize_only)
         USE kinds,            ONLY : DP
         USE mp_bands,         ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
         USE mp,               ONLY : mp_sum
         USE gvect,            ONLY : gstart
         USE control_flags,    ONLY : gamma_only
         USE additional_cusolver_subs, ONLY : cdiagh_gpu
#if defined(__CUDA)
         USE cublas
#endif

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: npwx, npw, m
         COMPLEX(DP), INTENT(INOUT) :: wfc (npwx,m)
         COMPLEX(DP), INTENT(INOUT) :: swfc(npwx,m)
         LOGICAL, INTENT(IN)        :: normalize_only

         COMPLEX(DP) :: temp
         COMPLEX(DP), ALLOCATABLE :: work(:,:), overlap(:,:)
         REAL(DP),    ALLOCATABLE :: e(:), overlap_gam(:,:)
         INTEGER :: i, j, k


         ALLOCATE(overlap(m,m))
         ALLOCATE(work(m,m))
         ALLOCATE(e(m))
         !$acc data create(overlap, work, e) present(wfc(:,:), swfc(:,:))

         !$acc kernels
         overlap(:,:) = (0.D0,0.D0)
         work(:,:) = (0.D0,0.D0)
         !$acc end kernels

         IF (gamma_only) THEN
            ALLOCATE(overlap_gam(m,m))
            !$acc data create(overlap_gam)

            !$acc host_data use_device(wfc(:,:), swfc(:,:), overlap_gam)
            CALL MYDGEMM('T', 'N', m, m, 2*npw, 2.D0, wfc, 2*npwx, swfc, 2*npwx, 0.D0, overlap_gam, m)
            IF (gstart == 2) THEN
               CALL MYDGER(m, m, -1.D0, wfc, 2*npwx, swfc, 2*npwx, overlap_gam, m)
            END IF
            CALL mp_sum(overlap_gam, intra_bgrp_comm)
            !$acc end host_data

            !$acc kernels
            overlap = CMPLX(overlap_gam, 0.D0, kind=DP)
            !$acc end kernels

            !$acc end data
            DEALLOCATE(overlap_gam)
         ELSE
            !$acc host_data use_device(wfc(:,:), swfc(:,:), overlap)
            CALL MYZGEMM(&
               'C', 'N', m, m, npw, (1.D0,0.D0), wfc, npwx,&
               swfc, npwx, (0.D0,0.D0), overlap, m)
            CALL mp_sum(overlap, intra_bgrp_comm)
            !$acc end host_data
         END IF

         IF (normalize_only) THEN
            !$acc parallel loop collapse(2) private(i, j)
            DO i=1,m
               DO j=1,m
                  IF (i .NE. j) THEN
                     overlap(i,j) = CMPLX(0.D0,0.D0, kind=DP)
                  END IF
               END DO
            END DO
         END IF

         !$acc host_data use_device(overlap, e, work)
         CALL cdiagh_gpu(m, overlap, m, e, work)
         !$acc end host_data

         ! calc O^(-1/2 T)
         !$acc parallel loop collapse(2) private(i, j, temp, k)
         DO i=1,m
            DO j=1,m
               IF (j < i) CYCLE
               temp = (0.D0, 0.D0)
               DO k=1,m
                  temp = temp + work(j,k) * (1.D0/SQRT(CMPLX(e(k),0.D0,kind=DP))) * CONJG(work(i,k))
               END DO
               overlap(i,j) = temp
               IF (j /= i) overlap(j,i) = CONJG(temp)
            END DO
         END DO

         ! lowdin ortho
         ! swfc = swfc * O^(-1/2 T)
         ! wfc  = wfc  * O^(-1/2 T)
         DO i=1,npw
            !$acc kernels
            work(:,1:2) = (0.D0,0.D0)
            !$acc end kernels

            !$acc host_data use_device(overlap, swfc(:,:), wfc(:,:), work)
            CALL cublasZgemv('n', m, m, (1.D0, 0.D0), overlap,&
                             m, swfc(i,1), npwx, (0.D0, 0.D0), work(1,1), 1)
            CALL cublasZgemv('n', m, m, (1.D0, 0.D0), overlap,&
                             m,  wfc(i,1), npwx, (0.D0, 0.D0), work(1,2), 1)
            CALL cublasZcopy(m, work(1,1), 1, swfc(i,1), npwx)
            CALL cublasZcopy(m, work(1,2), 1,  wfc(i,1), npwx)
            !$acc end host_data
         END DO
         !$acc end data
         DEALLOCATE(overlap, work, e)
      END SUBROUTINE oscdft_ortho_swfc_gpu

      SUBROUTINE get_overlap_gpu_ctx(ctx, ik, oatwfcS, wfcatom, swfcatom)
         USE symm_base, ONLY : nsym
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         INTEGER,                   INTENT(IN)    :: ik, oatwfcS(nsym,ctx%idx%orbs%norbs)
         COMPLEX(DP),               INTENT(IN)    :: wfcatom(:,:), swfcatom(:,:)

         !$acc data present(wfcatom, swfcatom)
         CALL get_overlap_gpu_internal(ctx%inp, ctx%idx, ik, oatwfcS, wfcatom, swfcatom)
         !$acc end data
      END SUBROUTINE get_overlap_gpu_ctx

      SUBROUTINE get_overlap_gpu_internal(inp, idx, ik, oatwfcS, wfcatom, swfcatom)
         USE control_flags, ONLY : gamma_only
         USE mp,            ONLY : mp_sum
         USE mp_bands,      ONLY : intra_bgrp_comm
         USE symm_base,     ONLY : nsym
         USE gvect,         ONLY : gstart
         USE klist,         ONLY : ngk
         USE wvfct,         ONLY : npwx
         IMPLICIT NONE

         TYPE(oscdft_input_type),           INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), TARGET, INTENT(INOUT) :: idx

         INTEGER,                           INTENT(IN)    :: ik, oatwfcS(nsym,idx%orbs%norbs)
         COMPLEX(DP),                       INTENT(IN)    :: wfcatom(:,:), swfcatom(:,:)

         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs


         INTEGER :: row, row_orb, row_m, row_off,&
                    col, col_orb, col_m, col_off,&
                    ioscdft, nsdim, isym, npw

         REAL(DP),    POINTER :: overlap_gam(:,:,:,:)
         COMPLEX(DP), POINTER :: overlap_k  (:,:,:,:)

         !$acc data present(wfcatom, swfcatom)
         orbs   => idx%orbs
         npw = ngk(ik)
         IF (gamma_only) THEN
            overlap_gam => idx%overlap_gam(:,:,:,:,ik)
            !$acc enter data create(overlap_gam(:,:,:,:))
            !$acc kernels
            overlap_gam = 0.D0
            !$acc end kernels
         ELSE
            overlap_k => idx%overlap_k(:,:,:,:,ik)
            !$acc enter data create(overlap_k(:,:,:,:))
            !$acc kernels
            overlap_k = (0.D0, 0.D0)
            !$acc end kernels
         END IF

         DO ioscdft=1,inp%noscdft
            DO isym=1,nsym
               row = 1
               DO row_orb=orbs%iorb_start(ioscdft),orbs%iorb_end(ioscdft)
                  row_m = 2 * orbs%l(row_orb) + 1
                  row_off = oatwfcS(isym,row_orb) + 1

                  col = 1
                  DO col_orb=orbs%iorb_start(ioscdft),orbs%iorb_end(ioscdft)
                     col_m = 2 * orbs%l(col_orb) + 1
                     col_off = oatwfcS(isym,col_orb) + 1

                     IF (gamma_only) THEN
                        !$acc host_data use_device(wfcatom, swfcatom, overlap_gam(:,:,:,:))
                        CALL MYDGEMM('T', 'N', row_m, col_m, 2*npwx, 2.D0,&
                           wfcatom(:,row_off), 2*npwx,&
                           swfcatom(:,col_off), 2*npwx,&
                           0.D0, overlap_gam(row,col,isym,ioscdft),&
                           idx%max_ns_dim)
                        IF (gstart == 2) THEN
                           CALL MYDGER(row_m, col_m, -1.D0,&
                              wfcatom(:,row_off), 2*npwx,&
                              swfcatom(:,col_off), 2*npwx,&
                              overlap_gam(row,col,isym,ioscdft),&
                              idx%max_ns_dim)
                        END IF
                        !$acc end host_data
                     ELSE
                        !$acc host_data use_device(wfcatom, swfcatom, overlap_k(:,:,:,:))
                        CALL MYZGEMM('C', 'N', row_m, col_m, npw, (1.D0,0.D0),&
                           wfcatom(:,row_off), npwx,&
                           swfcatom(:,col_off), npwx,&
                           (0.D0,0.D0), overlap_k(row,col,isym,ioscdft),&
                           idx%max_ns_dim)
                        !$acc end host_data
                     END IF

                     col = col + col_m
                  END DO

                  row = row + row_m
               END DO
            END DO
         END DO
         IF (gamma_only) THEN
            !$acc exit data copyout(overlap_gam(:,:,:,:))
            CALL mp_sum(idx%overlap_gam(:,:,:,:,ik), intra_bgrp_comm)
         ELSE
            !$acc exit data copyout(overlap_k(:,:,:,:))
            CALL mp_sum(idx%overlap_k  (:,:,:,:,ik), intra_bgrp_comm)
         END IF
         !$acc end data
      END SUBROUTINE get_overlap_gpu_internal

      SUBROUTINE lowdin_ortho_overlap_gpu_ctx(ctx)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx

         CALL lowdin_ortho_overlap_gpu_internal(ctx%inp, ctx%idx)
      END SUBROUTINE lowdin_ortho_overlap_gpu_ctx

      SUBROUTINE lowdin_ortho_overlap_gpu_internal(inp, idx)
         USE control_flags, ONLY : gamma_only
         USE klist,         ONLY : nks
         USE symm_base,     ONLY : nsym
         USE mp_bands,      ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
         USE additional_cusolver_subs, ONLY : cdiagh_gpu
         IMPLICIT NONE
         ! similar to oscdft_ortho_swfc_gpu

         TYPE(oscdft_input_type),   INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), INTENT(INOUT) :: idx

         COMPLEX(DP), ALLOCATABLE :: mat(:,:), eigv(:,:)
         REAL(DP),    ALLOCATABLE :: e(:)

         INTEGER     :: max_ns_dim, ik, ioscdft, isym, nsdim, i, j, k
         COMPLEX(DP) :: temp

         max_ns_dim = idx%max_ns_dim
         IF (max_ns_dim <= 0) THEN
            RETURN
         END IF
         ALLOCATE(mat(max_ns_dim,max_ns_dim), eigv(max_ns_dim,max_ns_dim), e(max_ns_dim))
         !$acc data create(mat, eigv, e)

         DO ik=1,nks
            DO ioscdft=1,inp%noscdft
               DO isym=1,nsym
                  nsdim = idx%ns_dim(ioscdft)

                  IF (gamma_only) THEN
                     mat(1:nsdim,1:nsdim) = CMPLX(idx%overlap_gam(1:nsdim,1:nsdim,isym,ioscdft,ik), 0.D0, kind=DP)
                  ELSE
                     mat(1:nsdim,1:nsdim) = idx%overlap_k(1:nsdim,1:nsdim,isym,ioscdft,ik)
                  END If
                  !$acc update device(mat)

                  !$acc host_data use_device(mat, e, eigv)
                  CALL cdiagh_gpu(nsdim, mat, max_ns_dim, e, eigv)
                  !$acc end host_data

                  !$acc parallel loop collapse(2) private(i, j, temp)
                  DO i=1,nsdim
                     DO j=1,nsdim
                        IF ( j < i ) CYCLE
                        temp = (0.D0,0.D0)
                        DO k=1,nsdim
                           temp = temp + eigv(j,k) * (1.D0/SQRT(CMPLX(e(k),0.D0,kind=DP))) * CONJG(eigv(i,k))
                        END DO
                        mat(j,i) = temp
                        IF (i.NE.j) mat(i,j) = CONJG(temp)
                     END DO
                  END DO

                  !$acc update host(mat)
                  idx%coeffs(1:nsdim,1:nsdim,isym,ioscdft,ik) = mat(1:nsdim,1:nsdim)
               END DO
            END DO
         END DO

         !$acc end data
         DEALLOCATE(mat, eigv, e)
      END SUBROUTINE lowdin_ortho_overlap_gpu_internal
#endif
END MODULE oscdft_wavefunction_subs_gpu
