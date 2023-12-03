MODULE oscdft_wavefunction_subs
#if defined (__OSCDFT)
   USE kinds, ONLY : DP
   USE oscdft_wavefunction, ONLY : oscdft_wavefunction_type
   USE io_global, ONLY : stdout
   USE oscdft_context, ONLY : oscdft_context_type
   USE oscdft_indices, ONLY : oscdft_indices_type,&
                              oscdft_constr_indices_type,&
                              oscdft_orbital_indices_type
   USE oscdft_input,   ONLY : oscdft_input_type
   USE oscdft_enums

   PRIVATE
   PUBLIC oscdft_init_wavefunctions,&
          oscdft_debug_print_wavefunctions,&
          oscdft_fill_wavefunctions,&
          oscdft_destroy_wavefunctions,&
          oscdft_ortho_swfc,&
          oscdft_get_overlap,&
          oscdft_write_overlap,&
          oscdft_poolrecover_overlap,&
          oscdft_lowdin_ortho_overlap,&
          oscdft_get_buffer

   INTERFACE oscdft_init_wavefunctions
      MODULE PROCEDURE init_ctx, init_internal
   END INTERFACE oscdft_init_wavefunctions
   INTERFACE oscdft_get_overlap
      MODULE PROCEDURE get_overlap_ctx, get_overlap_internal
   END INTERFACE oscdft_get_overlap
   INTERFACE oscdft_write_overlap
      MODULE PROCEDURE write_overlap_ctx, write_overlap_internal
   END INTERFACE oscdft_write_overlap
   INTERFACE oscdft_lowdin_ortho_overlap
      MODULE PROCEDURE lowdin_ortho_overlap_ctx, lowdin_ortho_overlap_internal
   END INTERFACE oscdft_lowdin_ortho_overlap

   CONTAINS
      SUBROUTINE init_ctx(ctx, wfc, iun, extension, constr_index, use_sym, have_S)
         IMPLICIT NONE
         TYPE(oscdft_context_type),      INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_wavefunction_type), INTENT(INOUT)         :: wfc

         INTEGER,          INTENT(IN) :: iun
         CHARACTER(LEN=*), INTENT(IN) :: extension
         LOGICAL,          INTENT(IN) :: constr_index, use_sym, have_S

         CALL init_internal(ctx%inp, ctx%idx, wfc, iun, extension, constr_index, use_sym, have_S)
      END SUBROUTINE init_ctx

      SUBROUTINE init_internal(inp, idx, wfc, iun, extension, constr_index, use_sym, have_S)
         USE uspp_param,       ONLY : upf
         USE ions_base,        ONLY : ityp, nat
         USE symm_base,        ONLY : nsym
         USE control_flags,    ONLY : io_level
         USE buffers,          ONLY : open_buffer
         USE basis,            ONLY : natomwfc
         USE noncollin_module, ONLY : npol
         USE wvfct,            ONLY : npwx

         IMPLICIT NONE
         TYPE(oscdft_input_type),        INTENT(INOUT)         :: inp
         TYPE(oscdft_indices_type),      INTENT(INOUT), TARGET :: idx
         TYPE(oscdft_wavefunction_type), INTENT(INOUT)         :: wfc

         INTEGER,          INTENT(IN) :: iun
         CHARACTER(LEN=*), INTENT(IN) :: extension
         LOGICAL,          INTENT(IN) :: constr_index, use_sym, have_S

         TYPE(oscdft_constr_indices_type),  POINTER :: constr
         TYPE(oscdft_orbital_indices_type), POINTER :: orbs

         INTEGER :: counter, m, na, nt, iwfc, n, l, isym, iorb,&
                    ioscdft, iconstr, iconstr_orbital, nnsym
         LOGICAL :: adv, test, exst

         constr => idx%constr
         orbs   => idx%orbs

         IF (wfc%initialized) RETURN

         wfc%iun = iun
         wfc%constr_index = constr_index
         wfc%use_sym = use_sym
         wfc%have_S = have_S

         wfc%n = 0
         m = MERGE(constr%norbs, orbs%norbs, constr_index)
         nnsym = MERGE(nsym, 1, use_sym)
         ALLOCATE(wfc%offset(nnsym,m), wfc%source(nnsym,m))
         wfc%offset = 0
         wfc%source = 0

         counter = 0
         DO na=1,nat
            nt = ityp(na)
            DO iwfc=1,upf(nt)%nwfc
               IF (upf(nt)%oc(iwfc) < 0.D0) CYCLE
               n = idx%nchi(iwfc,nt)
               l = upf(nt)%lchi(iwfc)
               adv = .false.
               DO iorb=1,orbs%norbs
                  ioscdft = orbs%iorb2ioscdft(iorb)
                  test = (orbs%n(iorb) == n).AND.(orbs%l(iorb) == l)
                  IF (constr_index) THEN
                     IF (inp%constraint_applied(ioscdft) == CONSTR_FALSE) CYCLE
                     iconstr = idx%ioscdft2iconstr(ioscdft)
                     iconstr_orbital = constr%iorb2icorb(iorb)
                     IF (use_sym) THEN
                        DO isym=1,nsym
                           IF (test .AND. (orbs%iat_sym(isym,iorb) == na)) THEN
                              adv = .true.
                              wfc%source(isym,iconstr_orbital) = counter
                              wfc%offset(isym,iconstr_orbital) = wfc%n
                           END IF
                        END DO
                     ELSE
                        IF (test .AND. (orbs%iat(iorb) == na)) THEN
                           adv = .true.
                           wfc%source(1,iconstr_orbital) = counter
                           wfc%offset(1,iconstr_orbital) = wfc%n
                        END IF
                     END IF
                  ELSE
                     IF (use_sym) THEN
                        DO isym=1,nsym
                           IF (test .AND. (orbs%iat_sym(isym,iorb) == na)) THEN
                              adv = .true.
                              wfc%source(isym,iorb) = counter
                              wfc%offset(isym,iorb) = wfc%n
                           END IF
                        END DO
                     ELSE
                        IF (test .AND. (orbs%iat(iorb) == na)) THEN
                           adv = .true.
                           wfc%source(1,iorb) = counter
                           wfc%offset(1,iorb) = wfc%n
                        END IF
                     END IF
                  END IF
               END DO
               counter = counter + 2 * l + 1
               IF (adv) wfc%n = wfc%n + 2 * l + 1
            END DO
         END DO
         IF (counter /= natomwfc) THEN
            CALL errore("oscdft_init_wavefunctions", "internal error: counter /= natomwfc", counter)
         END IF

         IF (wfc%iun /= 0) THEN
            wfc%nword = npwx * npol * wfc%n
            IF (wfc%n > 0) THEN
               CALL open_buffer(wfc%iun, extension, wfc%nword, io_level, exst)
               ALLOCATE(wfc%wfc(npwx*npol, wfc%n))
               wfc%wfc = (0.D0, 0.D0)
            END IF
         END IF
         wfc%initialized = .true.
      END SUBROUTINE init_internal

      SUBROUTINE oscdft_debug_print_wavefunctions(wfc, wfc_name, wfc_desc)
         USE symm_base, ONLY : nsym
         IMPLICIT NONE
         TYPE(oscdft_wavefunction_type), INTENT(IN) :: wfc
         CHARACTER(LEN=*),               INTENT(IN) :: wfc_name
         CHARACTER(LEN=*), OPTIONAL,     INTENT(IN) :: wfc_desc

         INTEGER :: isym

         IF (PRESENT(wfc_desc)) THEN
            WRITE(stdout, 100) wfc_name, wfc%n, TRIM(wfc_desc)
         ELSE
            WRITE(stdout, 100) wfc_name, wfc%n
         END IF
         IF (wfc%use_sym) THEN
            DO isym=1,nsym
               WRITE(stdout, 200) isym, wfc%offset(isym,:)
            END DO
            DO isym=1,nsym
               WRITE(stdout, 201) isym, wfc%source(isym,:)
            END DO
         ELSE
            WRITE(stdout, 300) wfc%offset(1,:)
            WRITE(stdout, 301) wfc%source(1,:)
         END IF
         WRITE(stdout, *) ""
         100 FORMAT("OSCDFT DEBUG: ", A, ": ", I0, " atomic wfcs", :, "; desc: ", A)
         200 FORMAT("OSCDFT DEBUG: |-offset(", I2, "): ", *(I0, :, " "))
         201 FORMAT("OSCDFT DEBUG: |-source(", I2, "): ", *(I0, :, " "))
         300 FORMAT("OSCDFT DEBUG: |-offset: ", *(I0, :, " "))
         301 FORMAT("OSCDFT DEBUG: |-source: ", *(I0, :, " "))
      END SUBROUTINE oscdft_debug_print_wavefunctions

      SUBROUTINE oscdft_fill_wavefunctions(idx, wfc, ik, wfcatom, swfcatom)
         USE symm_base,        ONLY : nsym
         USE buffers,          ONLY : save_buffer
         USE noncollin_module, ONLY : npol
         USE wvfct,            ONLY : npwx
         USE basis,            ONLY : natomwfc
         USE klist,            ONLY : nks

         IMPLICIT NONE
         TYPE(oscdft_indices_type),      INTENT(INOUT), TARGET :: idx
         TYPE(oscdft_wavefunction_type), INTENT(INOUT)         :: wfc
         INTEGER,                        INTENT(IN)            :: ik
         COMPLEX(DP),                    INTENT(IN)            ::  wfcatom(npwx*npol,natomwfc)
         COMPLEX(DP),                    INTENT(IN)            :: swfcatom(npwx*npol,natomwfc)

         TYPE(oscdft_constr_indices_type),  POINTER :: constr
         TYPE(oscdft_orbital_indices_type), POINTER :: orbs

         INTEGER :: nnsym, isym, iorb, iconstr_orb, m, m1, m2, n1, n2

         constr => idx%constr
         orbs   => idx%orbs

         IF (wfc%iun == 0) RETURN
         IF (wfc%n <= 0) RETURN

         nnsym = MERGE(nsym, 1, wfc%use_sym)
         wfc%wfc = (0.D0, 0.D0)
         IF (wfc%constr_index) THEN
            DO iconstr_orb=1,constr%norbs
               iorb = constr%icorb2iorb(iconstr_orb)
               m = 2 * orbs%l(iorb) + 1
               DO isym=1,nnsym
                  m1 = wfc%offset(isym,iconstr_orb) + 1
                  m2 = wfc%offset(isym,iconstr_orb) + m
                  n1 = wfc%source(isym,iconstr_orb) + 1
                  n2 = wfc%source(isym,iconstr_orb) + m
                  IF (wfc%have_S) THEN
                     wfc%wfc(:,m1:m2) = swfcatom(:,n1:n2)
                  ELSE
                     wfc%wfc(:,m1:m2) =  wfcatom(:,n1:n2)
                  END IF
               END DO
            END DO
         ELSE
            DO iorb=1,orbs%norbs
               m = 2 * orbs%l(iorb) + 1
               DO isym=1,nnsym
                  m1 = wfc%offset(isym,iorb) + 1
                  m2 = wfc%offset(isym,iorb) + m
                  n1 = wfc%source(isym,iorb) + 1
                  n2 = wfc%source(isym,iorb) + m
                  IF (wfc%have_S) THEN
                     wfc%wfc(:,m1:m2) = swfcatom(:,n1:n2)
                  ELSE
                     wfc%wfc(:,m1:m2) =  wfcatom(:,n1:n2)
                  END IF
               END DO
            END DO
         END IF
         IF (nks > 1) THEN
            CALL save_buffer(wfc%wfc, wfc%nword, wfc%iun, ik)
         END IF
      END SUBROUTINE oscdft_fill_wavefunctions

      SUBROUTINE oscdft_destroy_wavefunctions(wfc, status)
         USE buffers, ONLY : close_buffer
         IMPLICIT NONE

         TYPE(oscdft_wavefunction_type), INTENT(INOUT) :: wfc
         CHARACTER(LEN=*),               INTENT(IN)    :: status

         IF (wfc%iun /= 0) THEN
            CALL close_buffer(wfc%iun, status)
         END IF
         IF (ALLOCATED(wfc%wfc)) DEALLOCATE(wfc%wfc)
         IF (ALLOCATED(wfc%offset)) DEALLOCATE(wfc%offset)
         IF (ALLOCATED(wfc%source)) DEALLOCATE(wfc%source)
      END SUBROUTINE oscdft_destroy_wavefunctions

      SUBROUTINE oscdft_ortho_swfc(npwx, npw, m, wfc, swfc, normalize_only)
         USE kinds,            ONLY : DP
         USE mp_bands,         ONLY : intra_bgrp_comm
         USE mp,               ONLY : mp_sum
         USE gvect,            ONLY : gstart
         USE control_flags,    ONLY : gamma_only, use_gpu
         USE mp_bands,         ONLY : intra_bgrp_comm, me_bgrp, root_bgrp

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: npwx, npw, m
         COMPLEX(DP), INTENT(INOUT) :: wfc (npwx,m)
         COMPLEX(DP), INTENT(INOUT) :: swfc(npwx,m)
         LOGICAL, INTENT(IN)        :: normalize_only

         COMPLEX(DP) :: temp
         COMPLEX(DP), ALLOCATABLE ::  work (:,:), overlap (:,:), s(:,:)
         REAL(DP) , ALLOCATABLE  :: e (:), overlap_gam(:,:)
         !$acc declare device_resident(work, overlap, e, overlap_gam, s)
         INTEGER :: i, j, k

         !$acc data present(wfc, swfc)
         ALLOCATE (overlap(m,m))
         ALLOCATE (work   (m,m))
         ALLOCATE (e      (m))
         ALLOCATE (s      (m, m))

         !$acc kernels
         overlap(:,:) = (0.d0,0.d0)
         work(:,:)    = (0.d0,0.d0)
         !$acc end kernels

         IF (gamma_only) THEN
            ALLOCATE(overlap_gam(m,m))
            !$acc host_data use_device(wfc, swfc, overlap_gam)
            CALL MYDGEMM('T', 'N', m, m, 2*npw, 2.D0, wfc, 2*npwx, swfc,&
               2*npwx, 0.D0, overlap_gam, m)
            IF (gstart == 2) THEN
               CALL MYDGER(m, m, -1.D0, wfc, 2*npwx, swfc, 2*npwx, overlap_gam, m)
            END IF
            CALL mp_sum(overlap_gam, intra_bgrp_comm)
            !$acc end host_data

            !$acc kernels
            overlap = CMPLX(overlap_gam, 0.D0, kind=DP)
            !$acc end kernels

            DEALLOCATE(overlap_gam)
         ELSE
            !$acc host_data use_device(wfc, swfc, overlap)
            CALL MYZGEMM('C', 'N', m, m, npw, (1.D0, 0.D0), wfc, npwx,&
               swfc, npwx, (0.D0, 0.D0), overlap, m)
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

         ! calc O^(-1/2 T)
         IF (use_gpu) THEN
            !$acc kernels
            s(:,:) = CMPLX(0.d0,0.d0, kind=dp)
            DO i = 1, m
               s(i,i) = CMPLX(1.d0,0.d0, kind=dp)
            ENDDO
            !$acc end kernels

            !$acc host_data use_device(overlap, e, work, s)
            CALL laxlib_cdiaghg_gpu(m, m, overlap, s, m, e, work, me_bgrp,&
                                    root_bgrp, intra_bgrp_comm)
            !$acc end host_data
         ELSE
            CALL cdiagh (m, overlap, m, e, work)
         END IF

         ! calc O^(-1/2 T)
         !$acc parallel loop collapse(2) private(i, j, temp, k)
         DO i=1,m
            DO j=1,m
               IF (j < i) CYCLE
               temp = (0.d0, 0.d0)
               DO k=1,m
                  temp = temp + work (j, k) * (1.d0/SQRT(CMPLX(e(k),0.D0,kind=DP))) * CONJG (work (i, k) )
               END DO
               overlap (i, j) = temp
               IF (j /= i) overlap (j, i) = CONJG (temp)
            END DO
         END DO

         DEALLOCATE(work)

         ALLOCATE(work(npwx,m))

         ! lowdin ortho
         ! swfc = swfc * O^(-1/2 T)
         ! wfc  = wfc  * O^(-1/2 T)

         !$acc host_data use_device(overlap, swfc, work)
         CALL MYZGEMM('N', 'T', npwx, m, m, (1.D0, 0.D0), swfc, npwx,&
            overlap, m, (0.D0, 0.D0), work, npwx)
         !$acc end host_data
         !$acc parallel loop collapse(2)
         DO j=1,m
            DO i=1,npwx
               swfc(i,j) = work(i,j)
            END DO
         END DO

         !$acc host_data use_device(overlap, wfc, work)
         CALL MYZGEMM('N', 'T', npwx, m, m, (1.D0, 0.D0), wfc, npwx,&
            overlap, m, (0.D0, 0.D0), work, npwx)
         !$acc end host_data
         !$acc parallel loop collapse(2)
         DO j=1,m
            DO i=1,npwx
               wfc(i,j) = work(i,j)
            END DO
         END DO

         DEALLOCATE(overlap, work, e, s)
         !$acc end data
      END SUBROUTINE oscdft_ortho_swfc

      SUBROUTINE get_overlap_ctx(ctx, ik, oatwfcS, wfcatom, swfcatom)
         USE symm_base, ONLY : nsym
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         INTEGER,                   INTENT(IN)    :: ik, oatwfcS(nsym,ctx%idx%orbs%norbs)
         COMPLEX(DP),               INTENT(IN)    :: wfcatom(:,:), swfcatom(:,:)

         CALL get_overlap_internal(ctx%inp, ctx%idx, ik, oatwfcS, wfcatom, swfcatom)
      END SUBROUTINE get_overlap_ctx

      SUBROUTINE get_overlap_internal(inp, idx, ik, oatwfcS, wfcatom, swfcatom)
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
                    ioscdft, isym, npw

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
                           0.D0, idx%overlap_gam(row,col,isym,ioscdft,ik),&
                           idx%max_ns_dim)
                        IF (gstart == 2) THEN
                           CALL MYDGER(row_m, col_m, -1.D0,&
                              wfcatom(:,row_off), 2*npwx,&
                              swfcatom(:,col_off), 2*npwx,&
                              idx%overlap_gam(row,col,isym,ioscdft,ik),&
                              idx%max_ns_dim)
                        END IF
                        !$acc end host_data
                     ELSE
                        !$acc host_data use_device(wfcatom, swfcatom, overlap_k(:,:,:,:))
                        CALL MYZGEMM('C', 'N', row_m, col_m, npw, (1.D0,0.D0),&
                           wfcatom(:,row_off), npwx,&
                           swfcatom(:,col_off), npwx,&
                           (0.D0,0.D0), idx%overlap_k(row,col,isym,ioscdft,ik),&
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
      END SUBROUTINE get_overlap_internal

      SUBROUTINE write_overlap_ctx(ctx)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx

         CALL write_overlap_internal(ctx%inp, ctx%idx)
      END SUBROUTINE write_overlap_ctx

      SUBROUTINE write_overlap_internal(inp, idx)
         USE control_flags, ONLY : gamma_only
         USE klist,         ONLY : nkstot
         USE symm_base,     ONLY : nsym
         IMPLICIT NONE

         TYPE(oscdft_input_type),   INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), INTENT(INOUT) :: idx
         INTEGER :: ik, ioscdft, isym, nsdim, row

         DO ik=1,nkstot
            DO ioscdft=1,inp%noscdft
               DO isym=1,nsym
                  nsdim = idx%ns_dim(ioscdft)
                  DO row=1,nsdim
                     IF (gamma_only) THEN
                        WRITE(stdout, 300) row,isym,ioscdft,ik, idx%overlap_gam(row,1:nsdim,isym,ioscdft,ik)
                     ELSE
                        WRITE(stdout, 301) row,isym,ioscdft,ik, idx%overlap_k  (row,1:nsdim,isym,ioscdft,ik)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         DO ik=1,nkstot
            DO ioscdft=1,inp%noscdft
               DO isym=1,nsym
                  nsdim = idx%ns_dim(ioscdft)
                  DO row=1,nsdim
                     WRITE(stdout, 302) row,isym,ioscdft,ik, idx%coeffs(row,1:nsdim,isym,ioscdft,ik)
                  END DO
               END DO
            END DO
         END DO
         300 FORMAT("OSCDFT DEBUG: overlap(",I2, ",:,", I2, ",", I3,",", I3, "): ", *(F6.3, :, " "))
         301 FORMAT("OSCDFT DEBUG: overlap(",I2, ",:,", I2, ",", I3,",", I3, "): ", *(SS, F6.3, SP, F6.3, "i", SS, :, " "))
         302 FORMAT("OSCDFT DEBUG: coeffs(",I2, ",:,", I2, ",", I3,",", I3, "): ", *(SS, F6.3, SP, F6.3, "i", SS, :, " "))
      END SUBROUTINE write_overlap_internal

      SUBROUTINE oscdft_poolrecover_overlap(idx)
         USE control_flags, ONLY : gamma_only
         USE klist,         ONLY : nks, nkstot
         USE mp_pools,      ONLY : npool

         TYPE(oscdft_indices_type), INTENT(INOUT) :: idx
         INTEGER                                  :: length

         IF (npool <= 1) RETURN
         IF (gamma_only) THEN
            length = SIZE(idx%overlap_gam(:,:,:,:,1))
            CALL poolrecover(idx%overlap_gam, length, nkstot, nks)
         ELSE
            length = SIZE(idx%overlap_k  (:,:,:,:,1))
            CALL poolrecover(idx%overlap_k, 2*length, nkstot, nks)
         END IF
      END SUBROUTINE oscdft_poolrecover_overlap

      SUBROUTINE lowdin_ortho_overlap_ctx(ctx)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx

         CALL lowdin_ortho_overlap_internal(ctx%inp, ctx%idx)
      END SUBROUTINE lowdin_ortho_overlap_ctx

      SUBROUTINE lowdin_ortho_overlap_internal(inp, idx)
         USE control_flags, ONLY : gamma_only
         USE klist,         ONLY : nks
         USE symm_base,     ONLY : nsym
         IMPLICIT NONE

         TYPE(oscdft_input_type),   INTENT(INOUT) :: inp
         TYPE(oscdft_indices_type), INTENT(INOUT) :: idx

         COMPLEX(DP), ALLOCATABLE :: mat(:,:), eigv(:,:)
         REAL(DP),    ALLOCATABLE :: e(:)

         INTEGER     :: max_ns_dim, ik, ioscdft, isym, nsdim, i, j, k
         COMPLEX(DP) :: temp

         max_ns_dim = idx%max_ns_dim
         IF (max_ns_dim <= 0) RETURN

         ALLOCATE(mat(max_ns_dim,max_ns_dim), eigv(max_ns_dim,max_ns_dim), e(max_ns_dim))

         DO ik=1,nks
            DO ioscdft=1,inp%noscdft
               DO isym=1,nsym
                  nsdim = idx%ns_dim(ioscdft)

                  IF (gamma_only) THEN
                     mat(1:nsdim,1:nsdim) = CMPLX(idx%overlap_gam(1:nsdim,1:nsdim,isym,ioscdft,ik), 0.D0, kind=DP)
                  ELSE
                     mat(1:nsdim,1:nsdim) = idx%overlap_k(1:nsdim,1:nsdim,isym,ioscdft,ik)
                  END IF

                  CALL cdiagh(nsdim, mat, max_ns_dim, e, eigv)

                  ! idx%coeffs = O^(-1/2), where O = <atomic wfc|S|atomic wfc> (from idx%overlap_k/gam)

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

                  idx%coeffs(1:nsdim,1:nsdim,isym,ioscdft,ik) = mat(1:nsdim,1:nsdim)
               END DO
            END DO
         END DO

         DEALLOCATE(mat, eigv, e)
      END SUBROUTINE lowdin_ortho_overlap_internal

      SUBROUTINE oscdft_get_buffer(wfc, ik)
         USE buffers, ONLY : get_buffer
         IMPLICIT NONE
         TYPE(oscdft_wavefunction_type), INTENT(INOUT) :: wfc
         INTEGER, INTENT(IN)                           :: ik

         IF (wfc%nword == 0) RETURN
         CALL get_buffer(wfc%wfc, wfc%nword, wfc%iun, ik)
      END SUBROUTINE oscdft_get_buffer
#endif
END MODULE oscdft_wavefunction_subs
