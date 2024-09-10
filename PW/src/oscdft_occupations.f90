MODULE oscdft_occupations
#if defined (__OSCDFT)
   USE kinds,     ONLY : DP
   USE io_global, ONLY : ionode, ionode_id, stdout
   USE mp,        ONLY : mp_bcast
   USE mp_images, ONLY : intra_image_comm
   USE oscdft_enums
   USE oscdft_context,           ONLY : oscdft_context_type, oscdft_ns_type
   USE oscdft_wavefunction,      ONLY : oscdft_wavefunction_type
   USE oscdft_wavefunction_subs, ONLY : oscdft_get_buffer
   USE oscdft_input,             ONLY : oscdft_input_type
   USE oscdft_indices,           ONLY : oscdft_indices_type, oscdft_orbital_indices_type

   IMPLICIT NONE

   PRIVATE
   PUBLIC oscdft_new_ns, oscdft_get_occupation_numbers

   INTERFACE oscdft_new_ns
      MODULE PROCEDURE new_ns, new_ns_ctx, new_ns_ctx_debug
   END INTERFACE oscdft_new_ns
   INTERFACE oscdft_get_occupation_numbers
      MODULE PROCEDURE get_occ, get_occ_ctx
   END INTERFACE oscdft_get_occupation_numbers

   CONTAINS
      SUBROUTINE new_ns(inp, idx, wfcS, nst, wfc_evc)
         USE control_flags, ONLY : use_gpu
         IMPLICIT NONE
         TYPE(oscdft_input_type),                  INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), TARGET,        INTENT(INOUT) :: idx
         TYPE(oscdft_wavefunction_type),           INTENT(INOUT) :: wfcS
         TYPE(oscdft_ns_type),                     INTENT(INOUT) :: nst
         TYPE(oscdft_wavefunction_type), OPTIONAL, INTENT(INOUT) :: wfc_evc

         ! passing optional variables is ok 12.5.2.12.4
         IF (use_gpu) THEN
            CALL new_ns_normal_gpu(inp, idx, wfcS, nst, wfc_evc)
         ELSE
            CALL new_ns_normal(inp, idx, wfcS, nst, wfc_evc)
         END IF
      END SUBROUTINE new_ns

      SUBROUTINE new_ns_ctx(ctx)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx

         CALL new_ns(ctx%inp, ctx%idx, ctx%wfcS, ctx%nst)
      END SUBROUTINE new_ns_ctx

      SUBROUTINE new_ns_ctx_debug(ctx, debug_info)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         CHARACTER(LEN=*),          INTENT(IN)    :: debug_info

         IF (ctx%inp%debug_print) WRITE(stdout, 100) debug_info
         CALL new_ns(ctx%inp, ctx%idx, ctx%wfcS, ctx%nst)

         100 FORMAT("OSCDFT DEBUG: oscdft_new_ns from ", A)
      END SUBROUTINE new_ns_ctx_debug

      SUBROUTINE new_ns_normal(inp, idx, wfcS, nst, wfc_evc)
         USE klist,               ONLY : nks, ngk
         USE wvfct,               ONLY : nbnd, wg
         USE io_files,            ONLY : nwordwfc, iunwfc
         USE buffers,             ONLY : get_buffer
         USE becmod,              ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
         USE lsda_mod,            ONLY : isk, nspin
         USE wavefunctions,       ONLY : evc
         USE control_flags,       ONLY : gamma_only
         USE mp_pools,            ONLY : inter_pool_comm
         USE mp,                  ONLY : mp_sum
         USE symm_base,           ONLY : nsym
         USE oscdft_wavefunction, ONLY : check_bec_type_unallocated

         IMPLICIT NONE
         TYPE(oscdft_input_type),                  INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), TARGET,        INTENT(INOUT) :: idx
         TYPE(oscdft_wavefunction_type),           INTENT(INOUT) :: wfcS
         TYPE(oscdft_ns_type),                     INTENT(INOUT) :: nst
         TYPE(oscdft_wavefunction_type), OPTIONAL, INTENT(INOUT) :: wfc_evc
         TYPE(bec_type)                                          :: proj
         REAL(DP)                                                :: temp
         REAL(DP),    ALLOCATABLE                                :: nr(:,:,:,:)
         COMPLEX(DP), ALLOCATABLE                                :: nrtemp(:,:), nrtemp2(:,:)
         INTEGER                                                 :: ik, ioscdft, curr_dim, col, row,&
                                                                    ibnd, npw, maxl, isym,&
                                                                    col_off, row_off
         TYPE(oscdft_orbital_indices_type), POINTER              :: orbs

         orbs => idx%orbs

         ! There should be a better way for nr storage
         CALL start_clock("oscdft_ns")
         CALL check_bec_type_unallocated(proj)
         CALL allocate_bec_type(wfcS%n, nbnd, proj)
         maxl = idx%max_ns_dim
         ALLOCATE(nr(maxl,maxl,nsym,inp%noscdft))
         IF (inp%orthogonalize_ns) THEN
            ALLOCATE(nrtemp(maxl,maxl), nrtemp2(maxl,maxl))
         END IF

         IF (.NOT.ALLOCATED(evc)) THEN
            CALL errore("oscdft_new_ns", "evc not allocated", 1)
         END IF
         ! IF (inp%debug_print) WRITE(stdout, 300) "nbnd", nbnd

         nr(:,:,:,:) = 0.D0
         DO ik=1,nks
            npw = ngk(ik)
            IF (nks > 1) THEN
               IF (PRESENT(wfc_evc)) THEN
                  CALL get_buffer(evc, wfc_evc%nword, wfc_evc%iun, ik)
               ELSE
                  CALL get_buffer(evc, nwordwfc, iunwfc, ik)
               END IF
               CALL oscdft_get_buffer(wfcS, ik)
            END IF
            CALL calbec(npw, wfcS%wfc, evc, proj)
            DO ioscdft=1,inp%noscdft
               IF (inp%spin_index(ioscdft) /= isk(ik)) CYCLE
               curr_dim = idx%ns_dim(ioscdft)
               DO isym=1,nsym
                  IF (inp%orthogonalize_ns) THEN
                     nrtemp = (0.D0,0.D0)
                  END IF
                  DO col=1,curr_dim
                     col_off = wfcS%get_offset(idx%orbs, col, ioscdft, isym)
                     DO row=col,curr_dim
                        row_off = wfcS%get_offset(idx%orbs, row, ioscdft, isym)
                        temp = 0.D0
                        IF (gamma_only) THEN
!$omp parallel do reduction(+:temp)
                           DO ibnd=1,nbnd
                              temp = temp + proj%r(col_off,ibnd)*proj%r(row_off,ibnd)*wg(ibnd,ik)
                           END DO
!$omp end parallel do
                        ELSE
!$omp parallel do reduction(+:temp)
                           DO ibnd=1,nbnd
                              temp = temp + DBLE(proj%k(col_off,ibnd)*CONJG(proj%k(row_off,ibnd)))*wg(ibnd,ik)
                           END DO
!$omp end parallel do
                        END IF
                        IF (inp%orthogonalize_ns) THEN
                           nrtemp(row,col) = CMPLX(temp, 0.D0, kind=DP)
                           nrtemp(col,row) = CMPLX(temp, 0.D0, kind=DP)
                        ELSE
                           nr(row,col,isym,ioscdft) = nr(row,col,isym,ioscdft) + temp
                        END IF
                     END DO
                     ! IF (gamma_only) THEN
                     !    WRITE(stdout, 400) ik, col, (SQRT(wg(ibnd,ik)) * proj%r(col_off,ibnd), ibnd=1,nbnd)
                     ! ELSE
                     !    WRITE(stdout, 401) ik, col, (SQRT(wg(ibnd,ik)) * proj%k(col_off,ibnd), ibnd=1,nbnd)
                     ! END IF
                  END DO

                  IF (inp%orthogonalize_ns) THEN
                     CALL ZGEMM('N', 'N', curr_dim, curr_dim, curr_dim,&
                        (1.D0,0.D0), nrtemp, maxl,&
                        idx%coeffs(:,:,isym,ioscdft,ik), maxl,&
                        (0.D0,0.D0), nrtemp2, maxl)
                     CALL ZGEMM('C', 'N', curr_dim, curr_dim, curr_dim,&
                        (1.D0,0.D0), idx%coeffs(:,:,isym,ioscdft,ik), maxl,&
                        nrtemp2, maxl, (0.D0,0.D0), nrtemp, maxl)
                     nr(1:curr_dim,1:curr_dim,isym,ioscdft) = nr(1:curr_dim,1:curr_dim,isym,ioscdft) +&
                        DBLE(nrtemp(1:curr_dim,1:curr_dim))
                  END IF
               END DO
            END DO
         END DO
         CALL deallocate_bec_type(proj)
         CALL mp_sum(nr, inter_pool_comm)

         IF (inp%orthogonalize_ns) THEN
            DEALLOCATE(nrtemp, nrtemp2)
         END IF

         400 FORMAT("OSCDFT DEBUG: n2(", I3, ",", I2, "): ", *(F6.3, :, " "))
         401 FORMAT("OSCDFT DEBUG: n2(", I3, ",", I2, "): ", *(SS, F6.3, SP, F6.3, "i", SS, :, " "))


         CALL oscdft_fill_nr_upper(inp, idx, nr)
         CALL oscdft_symmetrize_ns(inp, idx, nst%ns, nr)
         DEALLOCATE(nr)

         CALL stop_clock("oscdft_ns")
         100 FORMAT("OSCDFT DEBUG: ", A, ": ", ES14.7)
         300 FORMAT("OSCDFT DEBUG: ", A, ": ", I5)
         102 FORMAT("OSCDFT DEBUG: nst%nr(", I5, ", ", I5, ")")
         103 FORMAT("OSCDFT DEBUG: nst%ns(", I5, ")")
         200 FORMAT(*(F10.7, " "))
      END SUBROUTINE new_ns_normal

      SUBROUTINE oscdft_fill_nr_upper(inp, idx, nr)
         USE symm_base,           ONLY : nsym
         IMPLICIT NONE
         TYPE(oscdft_input_type),                  INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), TARGET,        INTENT(INOUT) :: idx
         REAL(DP),                                 INTENT(INOUT) :: nr(:,:,:,:)

         INTEGER :: ioscdft, curr_dim, isym, col, row

         DO ioscdft=1,inp%noscdft
            curr_dim = idx%ns_dim(ioscdft)
            DO isym=1,nsym
               DO col=1,curr_dim
                  DO row=col+1,curr_dim
                     nr(col,row,isym,ioscdft) = nr(row,col,isym,ioscdft)
                  END DO
               END DO
            END DO
         END DO
         ! DO ioscdft=1,inp%noscdft
         !    curr_dim = idx%ns_dim(ioscdft)
         !    DO isym=1,nsym
         !       WRITE(stdout, 102) isym, ioscdft
         !       DO row=1,curr_dim
         !          WRITE(stdout, 200) nr(row,1:curr_dim,isym,ioscdft)
         !       END DO
         !    END DO
         ! END DO

         200 FORMAT(*(F10.7, " "))
      END SUBROUTINE oscdft_fill_nr_upper

      SUBROUTINE oscdft_symmetrize_ns(inp, idx, ns, nr)
         USE symm_base,           ONLY : nsym, d1, d2, d3
         IMPLICIT NONE

         TYPE(oscdft_input_type),                  INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), TARGET,        INTENT(INOUT) :: idx
         REAL(DP),                                 INTENT(IN)    :: nr(:,:,:,:)
         REAL(DP),                                 INTENT(OUT)   :: ns(:,:,:)
         TYPE(oscdft_orbital_indices_type), POINTER              :: orbs

         INTEGER  :: ioscdft, curr_dim, isym, col, row,&
                     col_orb, col_off, col_l, col_base,&
                     row_orb, row_off, row_l, row_base,&
                     m1, m0
         REAL(DP) :: psum, left_multiplier, right_multiplier

         orbs => idx%orbs

         ns(:,:,:) = 0.D0
         DO ioscdft=1,inp%noscdft
            curr_dim = idx%ns_dim(ioscdft)
            DO col=1,curr_dim
               col_orb  = orbs%ins2iorb(col,ioscdft)
               col_off  = orbs%ins2ioff(col,ioscdft)
               col_l    = orbs%l(col_orb)
               col_base = col - col_off
               DO row=1,curr_dim
                  row_orb  = orbs%ins2iorb(row,ioscdft)
                  row_off  = orbs%ins2ioff(row,ioscdft)
                  row_l    = orbs%l(row_orb)
                  row_base = row - row_off
                  DO isym=1,nsym
                     DO m1=1,2*col_l+1
                        DO m0=1,2*row_l+1
                           IF (row_l == 0) THEN
                              left_multiplier = 1.D0
                           ELSE IF (row_l == 1) THEN
                              left_multiplier = d1(m0,row_off,isym)
                           ELSE IF (row_l == 2) THEN
                              left_multiplier = d2(m0,row_off,isym)
                           ELSE IF (row_l == 3) THEN
                              left_multiplier = d3(m0,row_off,isym)
                           ELSE
                              CALL errore("oscdft_new_ns", "angular momentum not implemented",&
                                          ABS(row_l))
                           END IF
                           IF (col_l == 0) THEN
                              right_multiplier = 1.D0
                           ELSE IF (col_l == 1) THEN
                              right_multiplier = d1(m1,col_off,isym)
                           ELSE IF (col_l == 2) THEN
                              right_multiplier = d2(m1,col_off,isym)
                           ELSE IF (col_l == 3) THEN
                              right_multiplier = d3(m1,col_off,isym)
                           ELSE
                              CALL errore("oscdft_new_ns", "angular momentum not implemented",&
                                          ABS(col_l))
                           END IF
                           ns(row,col,ioscdft) = ns(row,col,ioscdft)+&
                              left_multiplier * nr(row_base+m0,col_base+m1,isym,ioscdft) *&
                              right_multiplier / nsym
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! DO ioscdft=1,inp%noscdft
         !    curr_dim = idx%ns_dim(ioscdft)
         !    WRITE(stdout, 103) ioscdft
         !    DO row=1,curr_dim
         !       WRITE(stdout, 200) ns(row,1:curr_dim,ioscdft)
         !    END DO
         ! END DO
         DO ioscdft=1,inp%noscdft
            curr_dim = idx%ns_dim(ioscdft)
            DO col=1,curr_dim
               DO row=col+1,curr_dim
                  psum = ABS(ns(row,col,ioscdft) - ns(col,row,ioscdft))
                  IF (psum > 1.D-10) THEN
                     WRITE(stdout, *) ioscdft, row, col
                     WRITE(stdout, *) ns(row,col,ioscdft)
                     WRITE(stdout, *) ns(col,row,ioscdft)
                     CALL errore("oscdft_new_ns", "non hermitian matrix", 1)
                  ELSE
                     ns(row,col,ioscdft) = 0.5D0 * (ns(row,col,ioscdft)+&
                                                        ns(col,row,ioscdft))
                     ns(col,row,ioscdft) = ns(row,col,ioscdft)
                  END IF
               END DO
            END DO
         END DO
         ! DO ioscdft=1,inp%noscdft
         !    curr_dim = idx%ns_dim(ioscdft)
         !    WRITE(stdout, 103) ioscdft
         !    DO row=1,curr_dim
         !       WRITE(stdout, 200) ns(row,1:curr_dim,ioscdft)
         !    END DO
         ! END DO

         100 FORMAT("OSCDFT DEBUG: ", A, ": ", ES14.7)
         300 FORMAT("OSCDFT DEBUG: ", A, ": ", I5)
         102 FORMAT("OSCDFT DEBUG: nst%nr(", I5, ", ", I5, ")")
         103 FORMAT("OSCDFT DEBUG: nst%ns(", I5, ")")
         200 FORMAT(*(F10.7, " "))
      END SUBROUTINE oscdft_symmetrize_ns

      SUBROUTINE new_ns_normal_gpu(inp, idx, wfcS, nst, wfc_evc)
         USE klist,               ONLY : nks, ngk
         USE wvfct,               ONLY : nbnd, wg
         USE io_files,            ONLY : nwordwfc, iunwfc
         USE buffers,             ONLY : get_buffer
         USE lsda_mod,            ONLY : isk, nspin
         USE wavefunctions,       ONLY : evc
         USE control_flags,       ONLY : gamma_only, offload_type
         USE mp_pools,            ONLY : inter_pool_comm
         USE mp,                  ONLY : mp_sum
         USE symm_base,           ONLY : nsym
         USE oscdft_wavefunction, ONLY : check_bec_type_unallocated
         USE becmod,              ONLY : bec_type, calbec,&
                                         allocate_bec_type_acc,&
                                         deallocate_bec_type_acc

         IMPLICIT NONE
         TYPE(oscdft_input_type),                  INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), TARGET,        INTENT(INOUT) :: idx
         TYPE(oscdft_wavefunction_type), TARGET,   INTENT(INOUT) :: wfcS
         TYPE(oscdft_ns_type),                     INTENT(INOUT) :: nst
         TYPE(oscdft_wavefunction_type), OPTIONAL, INTENT(INOUT) :: wfc_evc
         TYPE(bec_type), TARGET                                  :: proj
         REAL(DP)                                                :: temp
         REAL(DP),    ALLOCATABLE                                :: nr(:,:,:,:)
         COMPLEX(DP), ALLOCATABLE                                :: nrtemp(:,:), nrtemp2(:,:)
         INTEGER                                                 :: ik, ioscdft, curr_dim, col, row,&
                                                                    ibnd, npw, maxl, isym,&
                                                                    col_off, row_off
         TYPE(oscdft_orbital_indices_type), POINTER              :: orbs
         COMPLEX(DP), POINTER :: wfcS_wfc(:,:)
         REAL(DP),    POINTER :: proj_r(:,:)
         COMPLEX(DP), POINTER :: proj_k(:,:)

         orbs => idx%orbs

         CALL start_clock_gpu("oscdft_ns")
         CALL check_bec_type_unallocated(proj)
         CALL allocate_bec_type_acc(wfcS%n, nbnd, proj)

         IF (gamma_only) THEN
            proj_r => proj%r
         ELSE
            proj_k => proj%k
         END IF

         maxl = idx%max_ns_dim
         ALLOCATE(nr(maxl,maxl,nsym,inp%noscdft))
         IF (inp%orthogonalize_ns) THEN
            ALLOCATE(nrtemp(maxl,maxl), nrtemp2(maxl,maxl))
         END IF

         IF (.NOT.ALLOCATED(evc)) THEN
            CALL errore("oscdft_new_ns", "evc not allocated", 1)
         END IF

         wfcS_wfc => wfcS%wfc
         !$acc data create(wfcS_wfc(:,:)) present_or_create(evc) present_or_copyin(wg)

         nr(:,:,:,:) = 0.D0
         DO ik=1,nks
            npw = ngk(ik)
            IF (nks > 1) THEN
               IF (PRESENT(wfc_evc)) THEN
                  CALL get_buffer(evc, wfc_evc%nword, wfc_evc%iun, ik)
               ELSE
                  CALL get_buffer(evc, nwordwfc, iunwfc, ik)
               END IF
               CALL oscdft_get_buffer(wfcS, ik)
            END IF
            !$acc update device(evc)
            !$acc update device(wfcS_wfc(:,:))

            CALL calbec(offload_type, npw, wfcS_wfc, evc, proj)
            DO ioscdft=1,inp%noscdft
               IF (inp%spin_index(ioscdft) /= isk(ik)) CYCLE
               curr_dim = idx%ns_dim(ioscdft)
               DO isym=1,nsym
                  IF (inp%orthogonalize_ns) THEN
                     nrtemp = (0.D0,0.D0)
                  END IF
                  DO col=1,curr_dim
                     col_off = wfcS%get_offset(idx%orbs, col, ioscdft, isym)
                     DO row=col,curr_dim
                        row_off = wfcS%get_offset(idx%orbs, row, ioscdft, isym)
                        temp = 0.D0
                        IF (gamma_only) THEN
                           !$acc parallel loop reduction(+:temp) present(wg, proj_r)
                           DO ibnd=1,nbnd
                              temp = temp + proj_r(col_off,ibnd)*&
                                            proj_r(row_off,ibnd)*&
                                            wg(ibnd,ik)
                           END DO
                        ELSE
                           !$acc parallel loop reduction(+:temp) present(wg, proj_k)
                           DO ibnd=1,nbnd
                              temp = temp + DBLE(proj_k(col_off,ibnd)*&
                                                 CONJG(proj_k(row_off,ibnd)))*&
                                            wg(ibnd,ik)
                           END DO
                        END IF
                        IF (inp%orthogonalize_ns) THEN
                           nrtemp(row,col) = CMPLX(temp, 0.D0, kind=DP)
                           nrtemp(col,row) = CMPLX(temp, 0.D0, kind=DP)
                        ELSE
                           nr(row,col,isym,ioscdft) = nr(row,col,isym,ioscdft) + temp
                        END IF
                     END DO
                  END DO
                  IF (inp%orthogonalize_ns) THEN
                     CALL ZGEMM('N', 'N', curr_dim, curr_dim, curr_dim,&
                        (1.D0,0.D0), nrtemp, maxl,&
                        idx%coeffs(:,:,isym,ioscdft,ik), maxl,&
                        (0.D0,0.D0), nrtemp2, maxl)
                     CALL ZGEMM('C', 'N', curr_dim, curr_dim, curr_dim,&
                        (1.D0,0.D0), idx%coeffs(:,:,isym,ioscdft,ik), maxl,&
                        nrtemp2, maxl, (0.D0,0.D0), nrtemp, maxl)
                     nr(1:curr_dim,1:curr_dim,isym,ioscdft) = nr(1:curr_dim,1:curr_dim,isym,ioscdft) +&
                        DBLE(nrtemp(1:curr_dim,1:curr_dim))
                  END IF
               END DO
            END DO
         END DO
         !$acc end data
         CALL deallocate_bec_type_acc(proj)
         CALL mp_sum(nr, inter_pool_comm)

         IF (inp%orthogonalize_ns) THEN
            DEALLOCATE(nrtemp, nrtemp2)
         END IF

         CALL oscdft_fill_nr_upper(inp, idx, nr)
         CALL oscdft_symmetrize_ns(inp, idx, nst%ns, nr)
         DEALLOCATE(nr)
         CALL stop_clock_gpu("oscdft_ns")
      END SUBROUTINE new_ns_normal_gpu

      FUNCTION get_occ_sum(numbers, index_sum, ioscdft) RESULT(res)
         IMPLICIT NONE

         REAL(DP) :: res
         REAL(DP), INTENT(IN) :: numbers(:)
         INTEGER,  INTENT(IN) :: index_sum(:), ioscdft
         INTEGER              :: osum, isum

         res = numbers(ioscdft)
         DO isum=2,index_sum(1)
            osum = index_sum(isum+1)
            res = res + numbers(osum)
         END DO
         RETURN
      END FUNCTION get_occ_sum

      SUBROUTINE get_occ(inp, indx, nst, skip_print_occ)
         IMPLICIT NONE

         TYPE(oscdft_input_type),   INTENT(INOUT)         :: inp
         TYPE(oscdft_indices_type), INTENT(INOUT), TARGET :: indx
         TYPE(oscdft_ns_type),      INTENT(INOUT)         :: nst
         LOGICAL,                   INTENT(IN)            :: skip_print_occ
         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs

         LOGICAL  :: add_to_max
         INTEGER  :: iconstr, ioscdft, curr_dim, lda, info, i, idx, row,&
                     j,angular_momentum, max_idx, jconstr, joscdft,&
                     permute_idx(indx%nconstr), isum, osum, oidx,&
                     iorb, jorb, ina, jna, i_n, j_n, i_l, j_l
         REAL(DP) :: eigval(indx%max_ns_dim),&
                     eigvect(indx%max_ns_dim,indx%max_ns_dim),&
                     temp(indx%max_ns_dim,indx%max_ns_dim),&
                     temp2(indx%max_ns_dim,indx%max_ns_dim),&
                     work(15 * indx%max_ns_dim),&
                     sum_value, max_val, overlap, occup_sum_val,&
                     occup_sum_val_test, alpha, beta
         CHARACTER(LEN=4) :: spin_label(1:2)=(/&
            'UP  ', 'DOWN' /)
         CHARACTER(LEN=3) :: constr_label(0:8)=(/&
            'F  ', 'T  ', 'LE ', 'GE ', 'LE2', 'GE2', 'LE3', 'GE3', 'D  '/)

         orbs => indx%orbs

         lda = indx%max_ns_dim
         iconstr = 0
         DO ioscdft=1,inp%noscdft
            IF (inp%constraint_applied(ioscdft) /= 0) iconstr = iconstr + 1
            curr_dim = indx%ns_dim(ioscdft)
            eigvect(1:curr_dim,1:curr_dim) = nst%ns(1:curr_dim,1:curr_dim,ioscdft)
            CALL DSYEV("V", "U",&
                       curr_dim, eigvect, lda,&
                       eigval, work, 15*lda, info)
            IF (info /= 0) THEN
               WRITE(stdout, *) 'OSCDFT ERROR: dsyev info: ', info
               CALL errore("get_occupation_number", "dsyev fail", abs(info))
            END IF

            CALL mp_bcast(eigval(:),    ionode_id, intra_image_comm)
            CALL mp_bcast(eigvect(:,:), ionode_id, intra_image_comm)

            sum_value = SUM(eigval(1:curr_dim))
            IF (inp%print_occup(ioscdft).AND..NOT.skip_print_occ) THEN
               IF (ioscdft == 1) WRITE(stdout, *) ""
               WRITE(stdout, 100) ioscdft,&
                                  TRIM(constr_label(inp%constraint_applied(ioscdft))),&
                                  TRIM(spin_label(inp%spin_index(ioscdft))),&
                                  TRIM(inp%orbital_desc(ioscdft)),&
                                  sum_value
               WRITE(stdout, 101)
               WRITE(stdout, 200) eigval(1:curr_dim)
               IF (inp%print_occup_eigvects) THEN
                  WRITE(stdout, 102)
                  DO row=1,curr_dim
                     WRITE(stdout, 200) eigvect(row,1:curr_dim)
                  END DO
               END IF
               IF (inp%print_occup_matrix) THEN
                  WRITE(stdout, 103)
                  DO row=1,curr_dim
                     WRITE(stdout, 200) nst%ns(row,1:curr_dim,ioscdft)
                  END DO
               END IF
               ! WRITE(stdout, 104)
               ! WRITE(stdout, *) ""
            END IF
            nst%eigvals(:,ioscdft) = eigval(:)
            nst%eigvects(:,:,ioscdft) = eigvect(:,:)
            oidx = inp%occup_index(ioscdft)
            IF (oidx == OCCUP_TRACE) THEN
               nst%numbers(ioscdft) = sum_value
            ELSE IF (oidx == OCCUP_SUM) THEN
               nst%numbers(ioscdft) = eigval(inp%occup_index_sum(2,ioscdft))
            ELSE IF ((oidx > 0).AND.(oidx <= curr_dim)) THEN
               nst%numbers(ioscdft) = eigval(oidx)
            END IF
         END DO

         IF (inp%swapping_technique == OSCDFT_PERMUTE .AND. &
             ANY(inp%occup_index(indx%iconstr2ioscdft) > 0) .AND. &
             nst%eigvects_set) THEN
            WRITE(stdout, 501)
         END IF

         iconstr = 0
         alpha = 0.D0
         beta  = 1.D0
         DO ioscdft=1,inp%noscdft
            IF (inp%constraint_applied(ioscdft) /= 0) THEN
               iconstr = iconstr + 1
               curr_dim = indx%ns_dim(ioscdft)
               lda = indx%max_ns_dim
               eigval(:) = nst%eigvals(:,ioscdft)
               eigvect(:,:) = nst%eigvects(:,:,ioscdft)

               oidx = inp%occup_index(ioscdft)
               SELECT CASE (inp%swapping_technique)
                  CASE (OSCDFT_NONE)
                     nst%occup_eigvects(:,:,iconstr) = alpha*nst%occup_eigvects(:,:,iconstr)+&
                                                       beta*eigvect(:,:)
                     nst%occup_eigvals(:,iconstr) = eigval(:)
                     IF (oidx == OCCUP_TRACE) THEN
                        nst%occup_numbers(iconstr) = nst%numbers(ioscdft)
                     ELSE IF (oidx == OCCUP_SUM) THEN
                        occup_sum_val = get_occ_sum(nst%numbers, inp%occup_index_sum(:,ioscdft), ioscdft)
                        nst%occup_numbers(iconstr) = occup_sum_val
                     ELSE
                        nst%occup_numbers(iconstr) = eigval(inp%occup_index(ioscdft))
                     END IF
                  CASE (OSCDFT_PERMUTE)
                     IF (oidx == OCCUP_TRACE) THEN
                        nst%occup_eigvects(:,:,iconstr) = alpha*nst%occup_eigvects(:,:,iconstr)+beta*eigvect(:,:)
                        nst%occup_eigvals(:,iconstr) = eigval(:)
                        nst%occup_numbers(iconstr) = nst%numbers(ioscdft)
                     ELSE IF (oidx == OCCUP_SUM) THEN
                        occup_sum_val = get_occ_sum(nst%numbers, inp%occup_index_sum(:,ioscdft), ioscdft)
                        nst%occup_numbers(iconstr) = occup_sum_val
                        nst%occup_eigvals(:,iconstr) = eigval(:)
                        nst%occup_eigvects(:,:,iconstr) = alpha*nst%occup_eigvects(:,:,iconstr)+beta*eigvect(:,:)
                     ELSE
                        IF (nst%eigvects_set) THEN
                           max_val = -1
                           max_idx = 0
                           i = inp%occup_index(ioscdft)
                           DO idx=1,curr_dim
                              overlap = ABS(SUM(nst%occup_eigvects(1:curr_dim,i,iconstr)*&
                                                eigvect(1:curr_dim,idx)))
                              iorb = orbs%ins2iorb(idx,ioscdft)
                              ina  = orbs%iat(iorb)
                              i_n  = orbs%n(iorb)
                              i_l  = orbs%l(iorb)
                              IF (overlap > max_val) THEN
                                 add_to_max = .TRUE.
                                 DO jconstr=1,iconstr-1
                                    joscdft = indx%iconstr2ioscdft(jconstr)
                                    jorb = orbs%ins2iorb(permute_idx(jconstr),joscdft)
                                    jna  = orbs%iat(jorb)
                                    j_n  = orbs%n(jorb)
                                    j_l  = orbs%l(jorb)
                                    IF (ina == jna .AND. i_n == j_n .AND. i_l == j_l .AND.&
                                        permute_idx(jconstr) == idx) THEN
                                       add_to_max = .FALSE.
                                    END IF
                                 END DO
                                 IF (add_to_max) THEN
                                    max_idx = idx
                                    max_val = overlap
                                 END IF
                              END IF
                           END DO
                           idx = max_idx
                           CALL mp_bcast(idx, ionode_id, intra_image_comm)
                           WRITE(stdout, 500) ioscdft, idx, max_val

                           DO j=1,curr_dim
                              IF (j == i) THEN
                                 nst%occup_eigvects(:,j,iconstr) = alpha*nst%occup_eigvects(:,j,iconstr)+&
                                                                   beta*eigvect(:,idx)
                              ELSE IF (j == idx) THEN
                                 nst%occup_eigvects(:,j,iconstr) = alpha*nst%occup_eigvects(:,j,iconstr)+&
                                                                   beta*eigvect(:,i)
                              ELSE
                                 nst%occup_eigvects(:,j,iconstr) = alpha*nst%occup_eigvects(:,j,iconstr)+&
                                                                   beta*eigvect(:,j)
                              END IF
                           END DO
                           nst%occup_eigvals(:,iconstr) = eigval(:)
                           nst%occup_eigvals(i,iconstr) = eigval(idx)
                           nst%occup_eigvals(idx,iconstr) = eigval(i)
                           nst%occup_numbers(iconstr) = eigval(idx)
                           permute_idx(iconstr) = idx
                        ELSE
                           nst%occup_eigvects(1:curr_dim,1:curr_dim,iconstr) =&
                              alpha*nst%occup_eigvects(1:curr_dim,1:curr_dim,iconstr)+&
                              beta*eigvect(1:curr_dim,1:curr_dim)
                           nst%occup_eigvals(1:curr_dim,iconstr) = eigval(1:curr_dim)
                           nst%occup_numbers(iconstr) = eigval(inp%occup_index(ioscdft))
                        END IF
                     END IF
               END SELECT
            END IF
         END DO
         nst%eigvects_set = .true.
         100 FORMAT("OSCDFT: occupation #", I3, ": ", A, " ", A4, " ", A, "; tr[ns]: ", F8.4)
         101 FORMAT("OSCDFT: eigenvalue")
         102 FORMAT("OSCDFT: eigenvector (column)")
         103 FORMAT("OSCDFT: matrix (before diag.)")
         200 FORMAT("OSCDFT:   ", *(F7.4, :, " "))
         ! 101 FORMAT("OSCDFT: +-------eigenvalue tr[ns]:", F8.4, "-------")
         ! 102 FORMAT("OSCDFT: +--------------eigenvector---------------")
         ! 103 FORMAT("OSCDFT: +-----------------matrix-----------------")
         ! 104 FORMAT("OSCDFT: +----------------------------------------")
         ! 200 FORMAT("OSCDFT: | ", *(F7.4, :, " "))
         500 FORMAT("OSCDFT: occupation #", I3, ": following occup #", I0, " with overlap of ", F7.4)
         501 FORMAT("OSCDFT: PERMUTE")
      END SUBROUTINE get_occ

      SUBROUTINE get_occ_ctx(ctx, skip_print_occ)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         LOGICAL,                   INTENT(IN)    :: skip_print_occ

         CALL get_occ(ctx%inp, ctx%idx, ctx%nst, skip_print_occ)
      END SUBROUTINE get_occ_ctx
#endif
END MODULE oscdft_occupations
