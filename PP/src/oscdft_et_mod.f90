MODULE oscdft_et_mod
#if defined (__OSCDFT)
   USE kinds,                    ONLY : DP
   USE io_global,                ONLY : ionode, ionode_id, stdout
   USE oscdft_wavefunction,      ONLY : oscdft_wavefunction_type
   USE oscdft_input,             ONLY : oscdft_input_type
   USE oscdft_indices,           ONLY : oscdft_indices_type,&
                                        oscdft_orbital_indices_type,&
                                        oscdft_constr_indices_type
   USE oscdft_context,           ONLY : oscdft_ns_type
   USE oscdft_pp_mod,            ONLY : oscdft_copy_wavefunctions
   USE oscdft_wavefunction_subs, ONLY : oscdft_get_buffer
   USE oscdft_enums

   PRIVATE
   PUBLIC oscdft_file_type,&
          oscdft_wavefunction_type,&
          oscdft_read_file, get_et,&
          oscdft_close_file,&
          print_debug, oscdft_et_print_clock,&
          oscdft_nelup_neldw_from_input

   LOGICAL, SAVE :: print_debug = .FALSE.

   TYPE oscdft_file_type
      CHARACTER(LEN=256)             :: prefix, dir
      TYPE(oscdft_wavefunction_type) :: wfc, swfc, wfcO, wfcS
      TYPE(oscdft_input_type)        :: inp
      TYPE(oscdft_indices_type)      :: idx
      TYPE(oscdft_ns_type)           :: nst
      REAL(DP)                       :: energy
      REAL(DP), ALLOCATABLE          :: multipliers(:)
      REAL(DP), ALLOCATABLE          :: wg(:,:), wk(:)
   END TYPE oscdft_file_type

   CONTAINS
      SUBROUTINE oscdft_close_file(fil)
         USE oscdft_wavefunction_subs, ONLY : oscdft_destroy_wavefunctions
         IMPLICIT NONE

         TYPE(oscdft_file_type), INTENT(INOUT) :: fil

         CALL oscdft_destroy_wavefunctions(fil%wfc , "DELETE")
         CALL oscdft_destroy_wavefunctions(fil%swfc, "DELETE")
         CALL oscdft_destroy_wavefunctions(fil%wfcO, "DELETE")
         CALL oscdft_destroy_wavefunctions(fil%wfcS, "DELETE")
      END SUBROUTINE oscdft_close_file

      SUBROUTINE oscdft_get_swfc(swfc, wfc, extension_swfc)
         USE buffers,          ONLY : open_buffer, get_buffer, save_buffer
         USE wvfct,            ONLY : nbnd, npwx
         USE noncollin_module, ONLY : npol
         USE klist,            ONLY : nks, xk, ngk, igk_k
         USE control_flags,    ONLY : io_level
         USE uspp,             ONLY : nkb, vkb
         USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type,&
                                      becp, calbec, bec_type
         USE uspp_init,        ONLY : init_us_2

         IMPLICIT NONE
         TYPE(oscdft_wavefunction_type), INTENT(INOUT) :: swfc, wfc
         CHARACTER(LEN=*), INTENT(IN)                  :: extension_swfc
         LOGICAL                                       :: exst
         INTEGER                                       :: ik, npw

         swfc%n = nbnd
         swfc%nword = nbnd * npwx * npol
         CALL open_buffer(swfc%iun, extension_swfc, swfc%nword, io_level, exst)
         ALLOCATE(swfc%wfc(npwx * npol, nbnd))
         CALL allocate_bec_type(nkb, nbnd, becp)
         DO ik=1,nks
            npw = ngk(ik)
            IF (nks > 1) THEN
               CALL get_buffer(wfc%wfc, wfc%nword, wfc%iun, ik)
            END IF
            CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
            CALL calbec(npw, vkb, wfc%wfc, becp)
            CALL s_psi(npwx, npw, nbnd, wfc%wfc, swfc%wfc)
            IF (nks > 1) THEN
               CALL save_buffer(swfc%wfc, swfc%nword, swfc%iun, ik)
            END IF
         END DO
         CALL deallocate_bec_type(becp)
      END SUBROUTINE oscdft_get_swfc

      SUBROUTINE oscdft_read_oscdft_save(f)
         USE klist,     ONLY : nelup, neldw, nkstot
         USE mp_world,  ONLY : world_comm
         USE mp,        ONLY : mp_bcast
         USE wvfct,     ONLY : wg,nbnd
         USE lsda_mod,  ONLY : isk
         IMPLICIT NONE

         TYPE(oscdft_file_type), INTENT(INOUT) :: f
         INTEGER, EXTERNAL                     :: find_free_unit
         INTEGER                               :: u, ik, ibnd
         REAL(DP)                              :: nelupDP, neldwDP
         LOGICAL                               :: exst
         CHARACTER(LEN=256)                    :: filename

         ALLOCATE(f%multipliers(f%idx%nconstr))
         IF (ionode) THEN
            filename = TRIM(f%dir)//TRIM(f%prefix)//".save/oscdft_save"
            u = find_free_unit()
            INQUIRE(file=TRIM(filename), exist=exst)
            IF (.NOT.exst) CALL errore("read_oscdft", "missing "//TRIM(filename)//" file", 1)
            OPEN(unit=u, FILE=TRIM(filename), status="old")

            READ(u, 101) nelupDP, neldwDP
            READ(u, 101) f%energy
            IF (f%idx%nconstr > 0) THEN
               READ(u, 101) f%multipliers(1:f%idx%nconstr)
            END IF

            CLOSE(u)
         END IF

         nelup = 0.D0
         neldw = 0.D0
         DO ik=1,nkstot
            IF (isk(ik).EQ.1) THEN
               nelup = nelup + SUM(wg(1:nbnd,ik))
            ELSE
               neldw = neldw + SUM(wg(1:nbnd,ik))
            END IF
         END DO

         CALL mp_bcast(nelup,         ionode_id, world_comm)
         CALL mp_bcast(neldw,         ionode_id, world_comm)
         CALL mp_bcast(f%energy,      ionode_id, world_comm)
         IF (f%idx%nconstr > 0) THEN
            CALL mp_bcast(f%multipliers, ionode_id, world_comm)
         END IF

         IF (ionode) THEN
            WRITE(stdout, 102) "nel", nelup, neldw
            WRITE(stdout, 102) "total_energy", f%energy
            IF (f%idx%nconstr > 0) THEN
               WRITE(stdout, 102) "multipliers", f%multipliers(1:f%idx%nconstr)
            END IF
         ENDIF

         101 FORMAT(*(ES14.7, " "))
         102 FORMAT("OSCDFT: ", A, ": ", *(ES14.7, " "))
      END SUBROUTINE oscdft_read_oscdft_save

      SUBROUTINE oscdft_read_file(initial, final, ierr)
         USE io_files,             ONLY : nwordwfc, prefix, tmp_dir, wfc_dir
         USE scf,                  ONLY : rho
         USE gvect,                ONLY : g, ngm
         USE noncollin_module,     ONLY : npol
         USE control_flags,        ONLY : io_level
         USE wvfct,                ONLY : nbnd, npwx
         USE input_parameters,     ONLY : nat
         USE lsda_mod,             ONLY : nspin
         USE oscdft_context,       ONLY : oscdft_init_indices
         USE ions_base,            ONLY : ntyp => nsp
         USE symm_base,            ONLY : d1, d2, d3
         USe oscdft_input,         ONLY : oscdft_read_input

         IMPLICIT NONE

         TYPE(oscdft_file_type), INTENT(INOUT) :: initial, final
         INTEGER, INTENT(OUT)                  :: ierr
         LOGICAL                               :: exst, needwf
         INTEGER                               :: ik, i, nt
         CHARACTER(LEN=256)                    :: input_dft

         CALL start_clock("read_file")
         initial%wfc%iun = 41
         initial%swfc%iun = 42
         final%wfc%iun = 43
         final%swfc%iun = 44

         ierr = 0
         IF (ionode) THEN
            WRITE(stdout, 100) "Reading data from directory:", TRIM(initial%dir)
         END IF
         tmp_dir = initial%dir
         prefix = initial%prefix
         wfc_dir = initial%dir

         needwf = .true.
         CALL read_file_new(needwf) ! calls init_us_1, init_at_1
         IF (.NOT.needwf) THEN
            CALL errore("oscdft_read_file", "wavefunction not collected", 1)
         END IF

         nwordwfc = nbnd * npwx * npol
         ! io_level = 1
         CALL oscdft_copy_wavefunctions('wfc_initial', initial%wfc)
         CALL oscdft_get_swfc(initial%swfc, initial%wfc, 'swfc_initial')

         !
         tmp_dir = final%dir
         prefix  = final%prefix
         wfc_dir = final%dir

         CALL oscdft_copy_wavefunctions('wfc_final',   final%wfc)
         CALL oscdft_get_swfc(final%swfc,  final%wfc, 'swfc_final')

         IF (ionode) THEN
            WRITE(stdout, 100) "Reading data from directory:", TRIM(final%dir)
         END IF
         ! CALL read_final_dir2(final)

         CALL oscdft_read_input(initial%inp, TRIM(initial%dir)//'/'//TRIM(initial%prefix)//'.save/oscdft.in')
         CALL oscdft_read_input(final%inp, TRIM(final%dir)//'/'//TRIM(final%prefix)//'.save/oscdft.in')

         CALL oscdft_init_indices(initial%idx, initial%inp)
         CALL oscdft_init_indices(final%idx,   final%inp)

         CALL d_matrix(d1, d2, d3)

         WRITE(stdout, 200)
         CALL oscdft_read_oscdft_save(initial)
         WRITE(stdout, 201)
         CALL oscdft_read_oscdft_save(final)

         CALL oscdft_init_wfcO(initial, 45, "wfcO_initial", 46, "wfcS_initial")
         CALL oscdft_init_wfcO(final, 47, "wfcO_final", 48, "wfcS_final")

         CALL stop_clock("read_file")
         100 FORMAT("OSCDFT: ", A, A)
         200 FORMAT("OSCDFT: reading initial info");
         201 FORMAT("OSCDFT: reading final info");
      END SUBROUTINE oscdft_read_file

      SUBROUTINE oscdft_init_wfcO(f, iunO, extensionO, iunS, extensionS)
         USE basis,                    ONLY : natomwfc, swfcatom
         USE control_flags,            ONLY : gamma_only, io_level
         USE noncollin_module,         ONLY : npol
         USE wvfct,                    ONLY : npwx
         USE buffers,                  ONLY : open_buffer, get_buffer, save_buffer
         USE ions_base,                ONLY : nat, ityp
         USE uspp_param,               ONLY : upf
         USE uspp_init,                ONLY : init_us_2
         USE mp,                       ONLY : mp_barrier
         USE mp_pools,                 ONLY : inter_pool_comm, npool, my_pool_id, me_pool
         USE klist,                    ONLY : nks, xk, ngk, igk_k
         USE uspp,                     ONLY : nkb, vkb
         USE becmod,                   ONLY : allocate_bec_type, deallocate_bec_type,&
                                              becp, calbec, bec_type
         USE symm_base,                ONLY : nsym
         USE oscdft_wavefunction_subs, ONLY : oscdft_ortho_swfc, oscdft_get_overlap,&
                                              oscdft_write_overlap, oscdft_poolrecover_overlap,&
                                              oscdft_lowdin_ortho_overlap, oscdft_init_wavefunctions,&
                                              oscdft_fill_wavefunctions, oscdft_debug_print_wavefunctions
         IMPLICIT NONE

         TYPE(oscdft_file_type), TARGET, INTENT(INOUT) :: f
         TYPE(oscdft_indices_type),         POINTER    :: idx
         TYPE(oscdft_constr_indices_type),  POINTER    :: constr
         TYPE(oscdft_orbital_indices_type), POINTER    :: orbs
         TYPE(oscdft_input_type),           POINTER    :: inp
         INTEGER, INTENT(IN)      :: iunO, iunS
         CHARACTER(LEN=*)         :: extensionO, extensionS
         INTEGER                  :: npw, ik
         LOGICAL                  :: exst, in_wfcO, in_wfcS, test
         COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)

         inp  => f%inp
         idx  => f%idx
         constr => idx%constr
         orbs   => idx%orbs

         CALL start_clock("oscdft_wfcO")

         CALL oscdft_init_wavefunctions(inp, idx, f%wfcO, iunO, extensionO, .true., .false., .true.)
         CALL oscdft_init_wavefunctions(inp, idx, f%wfcS, iunS, extensionS, .false., .true., .true.)

         IF (print_debug) THEN
            CALL oscdft_debug_print_wavefunctions(f%wfcO, TRIM(extensionO), "S*atwfc to calc Hab")
            CALL oscdft_debug_print_wavefunctions(f%wfcS, TRIM(extensionS), "S*atwfc to calc ns")
         END IF

         ALLOCATE(wfcatom(npwx*npol,natomwfc),swfcatom(npwx*npol,natomwfc))
         CALL allocate_bec_type(nkb, natomwfc, becp)
         DO ik=1,nks
            CALL atomic_wfc(ik, wfcatom)
            npw = ngk(ik)
            CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
            CALL calbec(npw, vkb, wfcatom, becp)
            CALL s_psi(npwx, npw, natomwfc, wfcatom, swfcatom)

            IF (inp%orthogonalize_swfc) THEN
               CALL oscdft_ortho_swfc(npwx, npw, natomwfc, wfcatom, swfcatom, .false.)
            ELSE IF (inp%normalize_swfc) THEN
               CALL oscdft_ortho_swfc(npwx, npw, natomwfc, wfcatom, swfcatom, .true.)
            END IF

            CALL oscdft_fill_wavefunctions(idx, f%wfcO, ik, wfcatom, swfcatom)
            CALL oscdft_fill_wavefunctions(idx, f%wfcS, ik, wfcatom, swfcatom)

            IF (inp%orthogonalize_ns) THEN
               CALL oscdft_get_overlap(inp, idx, ik, f%wfcS%source, wfcatom, swfcatom)
            END IF
         END DO
         CALL deallocate_bec_type(becp)
         DEALLOCATE(wfcatom, swfcatom)

         IF (inp%orthogonalize_ns) THEN
            CALL oscdft_poolrecover_overlap(idx)
            CALL oscdft_lowdin_ortho_overlap(inp, idx)
            ! CALL oscdft_write_overlap(inp, idx)
         END IF
         CALL stop_clock("oscdft_wfcO")
         100 FORMAT("OSCDFT DEBUG: ", A, ": ", I6)
         102 FORMAT("OSCDFT DEBUG: ", A, ": ", *(I6, :, " "))
      END SUBROUTINE oscdft_init_wfcO

      SUBROUTINE get_dot_products(left, right, lSr, lWr)
         USE becmod,            ONLY : allocate_bec_type,&
                                       deallocate_bec_type,&
                                       bec_type, calbec
         USE control_flags,     ONLY : gamma_only
         USE mp_pools,          ONLY : inter_pool_comm
         USE mp_bands,          ONLY : intra_bgrp_comm
         USE wvfct,             ONLY : nbnd, npwx, wg
         USE klist,             ONLY : nelup, neldw, nks, wk,&
                                       igk_k, xk, ngk
         USE lsda_mod,          ONLY : isk
         USE mp,                ONLY : mp_sum
         USE mp_pools,          ONLY : inter_pool_comm
         USE noncollin_module,  ONLY : nspin_lsda
         IMPLICIT NONE

         TYPE(oscdft_file_type), INTENT(INOUT)  :: left, right
         COMPLEX(DP), INTENT(OUT)               :: lSr(nks), lWr(right%idx%nconstr, nks)
         TYPE(bec_type)                         :: projS, projL, projR
         COMPLEX(DP), ALLOCATABLE               :: work(:), Smatrix(:,:),&
                                                   Wmatrix(:,:),&
                                                   vl(:,:), vr(:,:)
         COMPLEX(DP)                            :: curr_output, eigval(nbnd),&
                                                   vltemp, vrtemp, rwork(2*nbnd),&
                                                   multiplier
         REAL(DP)                               :: occ_const
         INTEGER                                :: info, i, j, ik, h, k,&
                                                   ispin, m, worksize, npw,&
                                                   ipiv(nbnd), col, row, curr_dim,&
                                                   iconstr, ioscdft, oidx, osum,&
                                                   h_off, k_off

         lSr(:) = (0.D0, 0.D0)
         lWr(:,:) = (0.D0, 0.D0)
         CALL allocate_bec_type(nbnd, nbnd, projS)
         CALL allocate_bec_type(right%wfcO%n, nbnd, projL)
         CALL allocate_bec_type(right%wfcO%n, nbnd, projR)
         ALLOCATE(work(16*nbnd),&
                  Smatrix(nbnd,nbnd),&
                  Wmatrix(nbnd,nbnd),&
                  vl(nbnd,nbnd),&
                  vr(nbnd,nbnd))
         DO ik=1,nks
            ispin = isk(ik)
            m = NINT(MERGE(nelup, neldw, ispin.EQ.1))
            npw = ngk(ik)
            IF (nks > 1) THEN
               CALL oscdft_get_buffer(left%wfc, ik)
               CALL oscdft_get_buffer(right%swfc, ik)
               CALL oscdft_get_buffer(right%wfc, ik)
               CALL oscdft_get_buffer(right%wfcO, ik)
            END IF
            CALL calbec(npw, left%wfc%wfc, right%swfc%wfc, projS)
            CALL calbec(npw, right%wfcO%wfc, left%wfc%wfc, projL)
            CALL calbec(npw, right%wfcO%wfc, right%wfc%wfc, projR)

            IF (gamma_only) THEN
               Smatrix(1:m,1:m) = CMPLX(projS%r(1:m,1:m), 0.D0, kind=DP)
            ELSE
               Smatrix(1:m,1:m) = projS%k(1:m,1:m)
            END IF
            CALL ZGEEV('V', 'V', m, Smatrix, nbnd, eigval,&
               vl, nbnd, vr, nbnd, work, 16 * nbnd, rwork, info)
            curr_output = (1.D0, 0.D0)
            DO i=1,m
               curr_output = curr_output * eigval(i)
            END DO
            lSr(ik) = curr_output

            IF (gamma_only) THEN
               Smatrix(1:m,1:m) = CMPLX(projS%r(1:m,1:m), 0.D0, kind=DP)
            ELSE
               Smatrix(1:m,1:m) = projS%k(1:m,1:m)
            END IF
            CALL ZGETRF(m, m, Smatrix, nbnd, ipiv, info)
            CALL ZGETRI(m, Smatrix, nbnd, ipiv, work, 16 * nbnd, info)

            Smatrix(1:m,1:m) = Smatrix(1:m,1:m) * curr_output
            ! Smatrix TRANSPOSED is now the minor (determinant with i,j row col removed) multiplied by (-1)**(i+j)

            DO iconstr=1,right%idx%nconstr
               Wmatrix(:,:) = (0.D0, 0.D0)
               ioscdft = right%idx%iconstr2ioscdft(iconstr)
               IF (ispin.EQ.right%inp%spin_index(ioscdft)) THEN
                  curr_dim = right%idx%ns_dim(ioscdft)
                  oidx = right%inp%occup_index(ioscdft)
                  DO h=1,curr_dim
                     h_off = right%wfcO%get_offset(right%idx%constr, h, iconstr, -1)
                     DO k=1,curr_dim
                        k_off = right%wfcO%get_offset(right%idx%constr, k, iconstr, -1)
                        IF (oidx == OCCUP_TRACE) THEN
                           occ_const = MERGE(1.D0, 0.D0, h == k)
                        ELSE IF (oidx == OCCUP_SUM) THEN
                           osum = right%inp%occup_index_sum(2,ioscdft)
                           occ_const = right%nst%occup_eigvects(h,osum,iconstr)*&
                              right%nst%occup_eigvects(k,osum,iconstr)
                        ELSE
                           occ_const = right%nst%occup_eigvects(h,oidx,iconstr)*&
                              right%nst%occup_eigvects(k,oidx,iconstr)
                        END IF
                        multiplier = CMPLX(occ_const * right%multipliers(iconstr), 0.D0, kind=DP)

                        DO col=1,m
                           DO row=1,m
                              IF (gamma_only) THEN
                                 Wmatrix(row,col) = Wmatrix(row,col)+&
                                    multiplier*wg(col,ik)*&
                                    CMPLX(projL%r(h_off,row)*projR%r(k_off,col), 0.D0, kind=DP)
                              ELSE
                                 Wmatrix(row,col) = Wmatrix(row,col)+&
                                    multiplier*wg(col,ik)*&
                                    CONJG(projL%k(h_off,row))*projR%k(k_off,col)
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END IF
               DO i=1,m
                  DO j=1,m
                     lWr(iconstr,ik) = lWr(iconstr,ik) + Wmatrix(j,i) * Smatrix(i,j)
                  END DO
               END DO
            END DO
         END DO
         CALL deallocate_bec_type(projS)
         CALL deallocate_bec_type(projL)
         CALL deallocate_bec_type(projR)
         DEALLOCATE(work, Smatrix, Wmatrix, vl, vr)
         101 FORMAT("eigval(", I3, "): ", 2ES14.7, "; ABS: ", ES14.7)
         103 FORMAT("output for k: ", I2, " is ", 2ES14.7)
         104 FORMAT("i: ", I3, "; ipiv: ", I3)

         200 FORMAT("ik: ", I0, "; nocc: ", I0, "; nbnd: ", I0)
         !201 FORMAT("<psi_l(", I0, ")|psi_r(", I0, "): ", F10.7)
         !202 FORMAT("<psi_l(", I0, ")|psi_r(", I0, "): ", F10.7,"+", F10.7, "i")
         201 FORMAT(F10.7, ": <psi_l(", I5, ")|psi_r(", I5, ")>_", I0)
         202 FORMAT(F10.7, "+", F10.7, "i", ": <psi_l(", I5, ")|psi_r(", I5, ")>_", I0)
         203 FORMAT(F10.7, "+", F10.7, "i", ": <psi_l(", I5, ")|W|psi_r(", I5, ")>_", I0)
         301 FORMAT(F10.7, ": (", F10.7, ",", F10.7, "): eigval (", I5, ")_", I0)
         302 FORMAT("EIGVAL(", I5, ")_", I0, ": ", F10.7, ": (", F10.7, ",", F10.7, ")")
         303 FORMAT(F10.7, ": (", F10.7, ",", F10.7, "): eigvect (", I5, ",", I5, ")_", I0)
      END SUBROUTINE get_dot_products

      SUBROUTINE getHaa(f,Haa)
         IMPLICIT NONE
         TYPE(oscdft_file_type), INTENT(INOUT) :: f
         COMPLEX(DP), INTENT(OUT)              :: Haa
         REAL(DP)                              :: toAdd
         INTEGER                               :: iconstr, ioscdft, oidx, osum
         Haa = CMPLX(f%energy, 0.D0, kind=DP)
         DO iconstr=1,f%idx%nconstr
            ioscdft = f%idx%iconstr2ioscdft(iconstr)
            oidx = f%inp%occup_index(ioscdft)
            IF (oidx == OCCUP_SUM) THEN
               osum = f%inp%occup_index_sum(2,ioscdft)
               toAdd = f%multipliers(iconstr) * f%nst%occup_eigvals(osum,iconstr)
            ELSE
               toAdd = f%multipliers(iconstr) * f%inp%target_occup(ioscdft)
            END IF
            ! toAdd = f%multipliers(iconstr) * f%inp%target_occup(ioscdft)
            ! IF (f%inp%occup_index(ioscdft) == OCCUP_SUM) THEN
            !    toAdd = toAdd / f%inp%occup_index_sum(1,ioscdft)
            ! END IF
            Haa = Haa + CMPLX(toAdd, 0.D0, kind=DP)
         END DO
      END SUBROUTINE getHaa

      SUBROUTINE getWaa(f,Waa)
         IMPLICIT NONE
         TYPE(oscdft_file_type), INTENT(INOUT) :: f
         COMPLEX(DP), INTENT(OUT)              :: Waa
         REAL(DP)                              :: toAdd
         INTEGER                               :: iconstr, ioscdft, oidx, osum
         Waa = (0.D0, 0.D0)
         DO iconstr=1,f%idx%nconstr
            ioscdft = f%idx%iconstr2ioscdft(iconstr)
            oidx = f%inp%occup_index(ioscdft)
            IF (oidx == OCCUP_SUM) THEN
               osum = f%inp%occup_index_sum(2,ioscdft)
               toAdd = f%multipliers(iconstr) * f%nst%occup_eigvals(osum,iconstr)
            ELSE
               toAdd = f%multipliers(iconstr) * f%inp%target_occup(ioscdft)
            END IF
            ! toAdd = f%multipliers(iconstr) * f%inp%target_occup(ioscdft)
            ! IF (f%inp%occup_index(ioscdft) == OCCUP_SUM) THEN
            !    toAdd = toAdd / f%inp%occup_index_sum(1,ioscdft)
            ! END IF
            Waa = Waa + CMPLX(toAdd, 0.D0, kind=DP)
         END DO
      END SUBROUTINE getWaa

      SUBROUTINE oscdft_write_occupations
         USE wvfct,                ONLY : nbnd, wg
         USE klist,                ONLY : nelup, neldw, nkstot, xk
         USE lsda_mod,             ONLY : lsda
         IMPLICIT NONE
         INTEGER :: ibnd, jbnd, ik

         IF (.NOT.ionode) RETURN
         WRITE(stdout, 100) nelup, neldw, nbnd

         DO ik=1,nkstot
            IF (lsda) THEN
               IF (ik  == 1) WRITE(stdout, 101)
               IF (ik  == (1 + nkstot/2)) WRITE(stdout, 102) 
            END IF
            WRITE(stdout, 103) ik, xk(1:3,ik)
            DO ibnd=1,nbnd,8
               jbnd = MIN(nbnd,ibnd+7)
               WRITE(stdout, 104) ibnd, jbnd, wg(ibnd:jbnd,ik)
            END DO
         END DO

         100 FORMAT("nelup: ", F0.2, "; neldw: ", F0.2, "; nbnd: ", I0)
         101 FORMAT("------ SPIN UP ----------")
         102 FORMAT("------ SPIN DOWN --------")
         103 FORMAT("k(", I0, ") = ", 3F7.4)
         104 FORMAT("wg(", I5, ":", I5, "): ", 8(F7.4, " "))
      END SUBROUTINE oscdft_write_occupations

      SUBROUTINE get_et(initial, final)
         USE noncollin_module,   ONLY : nspin_lsda
         USE klist,              ONLY : nelup, neldw, wk, nkstot, nks, xk
         USE oscdft_occupations, ONLY : oscdft_new_ns, oscdft_get_occupation_numbers
         USE oscdft_context,     ONLY : oscdft_alloc_nst
         USE mp_pools,           ONLY : me_pool, npool, inter_pool_comm, my_pool_id
         USE mp,                 ONLY : mp_barrier
         USE lsda_mod,           ONLY : isk
         USE wvfct,              ONLY : nbnd

         IMPLICIT NONE
         TYPE(oscdft_file_type), INTENT(INOUT) :: initial, final
         COMPLEX(DP)                          :: fdoti(nkstot),&
                                                 idotf(nkstot),&
                                                 idoti(nkstot),&
                                                 fdotf(nkstot),&
                                                 Sfi, Sif, Sii, Sff,&
                                                 Wii, Wff,&
                                                 Wif, Wfi,&
                                                 Hii, Hff,&
                                                 Hfi, Hif,&
                                                 Hab, Sab,&
                                                 H(2,2), S(2,2),&
                                                 W(2,2), HC(2,2),&
                                                 CHC(2,2), eigvect(2,2),&
                                                 eigval(2), swork(2,2), work(192), rwork(4)
         COMPLEX(DP), ALLOCATABLE             :: iWi(:,:), iWf(:,:), fWi(:,:), fWf(:,:)
         INTEGER                              :: info, i, ispin, ik, n, ipool, iconstr, iik

         CALL start_clock("get_et")
         IF (print_debug) THEN
            WRITE(stdout, 102) "nelup", nelup
            WRITE(stdout, 102) "neldw", neldw
            WRITE(stdout, 101) "nkstot", nkstot
         END IF

         CALL oscdft_alloc_nst(initial%nst, initial%idx%max_ns_dim,&
                               initial%idx%nconstr, initial%inp%noscdft)
         CALL oscdft_alloc_nst(final%nst, final%idx%max_ns_dim,&
                               final%idx%nconstr, final%inp%noscdft)

         CALL oscdft_write_occupations

         WRITE(stdout, 600)
         CALL oscdft_new_ns(initial%inp, initial%idx,&
            initial%wfcS, initial%nst, initial%wfc)
         CALL oscdft_get_occupation_numbers(initial%inp, initial%idx, initial%nst, .false.)

         WRITE(stdout, 601)
         CALL oscdft_new_ns(final%inp, final%idx,&
            final%wfcS, final%nst, final%wfc)
         CALL oscdft_get_occupation_numbers(final%inp, final%idx, final%nst, .false.)
         WRITE(stdout, *) ""

         ALLOCATE(iWi(initial%idx%nconstr, nkstot),&
                  iWf(  final%idx%nconstr, nkstot),&
                  fWi(initial%idx%nconstr, nkstot),&
                  fWf(  final%idx%nconstr, nkstot))

         CALL get_dot_products(final,  initial,  fdoti, fWi)
         CALL get_dot_products(initial, final,   idotf, iWf)
         CALL get_dot_products(initial, initial, idoti, iWi)
         CALL get_dot_products(final,   final,   fdotf, fWf)

         IF (print_debug) THEN
            DO ipool=0,(npool-1)
               IF ((ipool.EQ.my_pool_id).AND.(me_pool.EQ.0)) THEN
                  WRITE(*, 201) "ipool", ipool
                  WRITE(*, 205) "isk",   isk(1:nks)
                  WRITE(*, 204) "fdoti", fdoti(1:nks)
                  WRITE(*, 204) "idotf", idotf(1:nks)
                  WRITE(*, 204) "idoti", idoti(1:nks)
                  WRITE(*, 204) "fdotf", fdotf(1:nks)
                  WRITE(*, 204) "fWi",   fWi(:,1:nks)
                  WRITE(*, 204) "iWf",   iWf(:,1:nks)
                  WRITE(*, 204) "iWi",   iWi(:,1:nks)
                  WRITE(*, 204) "fWf",   fWf(:,1:nks)
               ENDIF
               CALL SLEEP(1)
               CALL mp_barrier(inter_pool_comm)
            ENDDO
         END IF

         CALL poolrecover(fdoti, 2, nkstot, nks)
         CALL poolrecover(idotf, 2, nkstot, nks)
         CALL poolrecover(idoti, 2, nkstot, nks)
         CALL poolrecover(fdotf, 2, nkstot, nks)
         CALL poolrecover(fWi, 2 * initial%idx%nconstr, nkstot, nks)
         CALL poolrecover(iWf, 2 *   final%idx%nconstr, nkstot, nks)
         CALL poolrecover(iWi, 2 * initial%idx%nconstr, nkstot, nks)
         CALL poolrecover(fWf, 2 *   final%idx%nconstr, nkstot, nks)
         CALL poolrecover(wk, 1, nkstot, nks)

         WRITE(stdout, 104) "fdoti", fdoti(1:nkstot)
         WRITE(stdout, 104) "idotf", idotf(1:nkstot)
         WRITE(stdout, 104) "idoti", idoti(1:nkstot)
         WRITE(stdout, 104) "fdotf", fdotf(1:nkstot)
         DO ik=1,nkstot
            WRITE(stdout, 400) ik
            WRITE(stdout, 104) "fdoti", fdoti(ik)
            WRITE(stdout, 104) "idotf", idotf(ik)
            WRITE(stdout, 104) "idoti", idoti(ik)
            WRITE(stdout, 104) "fdotf", fdotf(ik)
            WRITE(stdout, 104) "fWi", fWi(:,ik)
            WRITE(stdout, 104) "iWf", iWf(:,ik)
            WRITE(stdout, 104) "iWi", iWi(:,ik)
            WRITE(stdout, 104) "fWf", fWf(:,ik)
         END DO

         IF (print_debug) THEN
            DO ik=1,nkstot
               WRITE(stdout, 500) ik, isk(ik), ik, xk(:,ik), ik, wk(ik)
            END DO
         END IF

         IF (nspin_lsda.EQ.1) THEN
            DO ik=1,nkstot
               fdoti(ik) = wk(ik) * fdoti(ik)**2
               idotf(ik) = wk(ik) * idotf(ik)**2
               idoti(ik) = wk(ik) * idoti(ik)**2
               fdotf(ik) = wk(ik) * fdotf(ik)**2
            ENDDO
            Sif = SUM(idotf)
            Sfi = SUM(fdoti)
            Sii = SUM(idoti)
            Sff = SUM(fdotf)

            Wif = SUM(iWf)
            Wfi = SUM(fWi)
            Wii = SUM(iWi)
            Wff = SUM(fWf)
         ELSE
            n = nkstot/2
            Sif = (0.D0, 0.D0)
            Sfi = (0.D0, 0.D0)
            Sii = (0.D0, 0.D0)
            Sff = (0.D0, 0.D0)
            DO ik=1,n
               Sif = Sif + wk(ik) * idotf(ik) * idotf(ik+n)
               Sfi = Sfi + wk(ik) * fdoti(ik) * fdoti(ik+n)
               Sii = Sii + wk(ik) * idoti(ik) * idoti(ik+n)
               Sff = Sff + wk(ik) * fdotf(ik) * fdotf(ik+n)
            ENDDO

            Wif = (0.D0, 0.0)
            Wfi = (0.D0, 0.0)
            Wii = (0.D0, 0.0)
            Wff = (0.D0, 0.0)
            DO ik=1,nkstot
               iik = ik+n
               IF (ik > n) THEN
                  iik = ik - n
               ENDIF
               DO iconstr=1,initial%idx%nconstr
                  Wii = Wii + iWi(iconstr,ik) * idoti(iik)
                  Wfi = Wfi + fWi(iconstr,ik) * fdoti(iik)
               ENDDO
               DO iconstr=1,final%idx%nconstr
                  Wif = Wif + iWf(iconstr,ik) * idotf(iik)
                  Wff = Wff + fWf(iconstr,ik) * fdotf(iik)
               ENDDO
            ENDDO
         ENDIF

         WRITE(stdout, *) ""

         ! IF (nspin_lsda.EQ.2) THEN
         !    Sif = idotf(1) * idotf(2)
         !    Sfi = fdoti(1) * fdoti(2)
         ! ELSE
         !    Sif = idotf(1) * idotf(1)
         !    Sfi = fdoti(1) * fdoti(1)
         ! ENDIF
         ! Wif = DBLE(iWf)
         ! Wfi = DBLE(fWi)

         Sab = 0.5D0 * (CONJG(Sfi) + Sif)
         WRITE(stdout, 104) "Sif", Sif
         WRITE(stdout, 104) "Sfi", Sfi
         WRITE(stdout, 104) "Sab", Sab

         S(1,1) = 1.D0
         S(2,2) = 1.D0
         S(1,2) = Sab
         S(2,1) = CONJG(Sab)

         CALL getHaa(initial, Hii)
         CALL getHaa(final, Hff)

         ! Hfi = initial%energy * Sfi
         ! Hif = final%energy * Sif
         ! Hab = DBLE(Hfi + Hif) * 0.5D0
         Hfi = Hii * Sfi - Wfi
         Hif = Hff * Sif - Wif
         Hab = 0.5D0 * (CONJG(Hfi) + Hif)
         WRITE(stdout, 104) "Hfi", Hfi
         WRITE(stdout, 104) "Hif", Hif
         WRITE(stdout, 104) "Hab", Hab

         ! H(1,1) = Hii - Sab * Hab
         ! H(2,2) = Hff - Sab * Hab
         ! H(1,2) = Hab - Sab * Hff
         ! H(2,1) = Hab - Sab * Hii
         ! H = H / (1 - Sab * Sab)
         H(1,1) = initial%energy
         H(2,2) = final%energy
         H(1,2) = Hab
         H(2,1) = CONJG(Hab)

         WRITE(stdout, 104) "Wii", Wii
         WRITE(stdout, 104) "Wff", Wff
         CALL getWaa(initial, Wii)
         CALL getWaa(final, Wff)

         WRITE(stdout, 104) "Wii (check with iWi)", Wii
         WRITE(stdout, 104) "Wff (check with fWf)", Wff

         W(1,1) = Wii
         W(2,2) = Wff
         W(1,2) = DBLE(Wif + Wfi) * 0.5D0
         W(1,2) = 0.5D0 * (CONJG(Wfi) + Wif)
         W(2,1) = CONJG(W(1,2))

         eigvect(:,:) = W(:,:)
         swork(:,:) = S(:,:)
         CALL ZHEGV(1,'V','U',2,eigvect,2,swork,2,eigval,work,192,rwork,info)
         WRITE(stdout, 101) "info", info

         CALL ZGEMM('N','N',2,2,2,(1.D0,0.D0),&
            H,2,eigvect,2,&
            (0.D0,0.D0),HC,2)
         CALL ZGEMM('C','N',2,2,2,(1.D0,0.D0),&
            eigvect,2,HC,2,&
            (0.D0,0.D0),CHC,2)
         DO i=1,2
            WRITE(stdout, 104) "H", H(i,:)
         ENDDO
         DO i=1,2
            WRITE(stdout, 104) "W", W(i,:)
         ENDDO
         DO i=1,2
            WRITE(stdout, 104) "S", S(i,:)
         ENDDO
         DO i=1,2
            WRITE(stdout, 104) "C", eigvect(i,:)
         ENDDO
         DO i=1,2
            WRITE(stdout, 104) "CHC", CHC(i,:)
         ENDDO

         WRITE(stdout, 102) "H_II", initial%energy
         WRITE(stdout, 102) "H_FF", final%energy
         WRITE(stdout, 104) "S_IF", Sab
         WRITE(stdout, 104) "H_IF", Hab

         WRITE(stdout, 102) "ABS(Hab)", ABS(CHC(1,2))

         DEALLOCATE(iWi, iWf, fWi, fWf)

         CALL stop_clock("get_et")
         100 FORMAT("OSCDFT: ", A, ": ", 2(ES14.7, " "), "; ", 2(ES14.7, " "))
         101 FORMAT("OSCDFT: ", A, ": ", I5)
         102 FORMAT("OSCDFT: ", A, ": ", ES14.7)
         103 FORMAT("OSCDFT: ", A, ": ", 2(ES14.7, " "))
         104 FORMAT("OSCDFT: ", A, ": ", *(ES14.7, " + ", ES14.7, " i, ", :))
         105 FORMAT("OSCDFT: ", A, ": ", *(I5, " ", :))

         200 FORMAT("OSCDFT DEBUG: ", A, ": ", 2(ES14.7, " "), "; ", 2(ES14.7, " "))
         201 FORMAT("OSCDFT DEBUG: ", A, ": ", I5)
         202 FORMAT("OSCDFT DEBUG: ", A, ": ", ES14.7)
         203 FORMAT("OSCDFT DEBUG: ", A, ": ", 2(ES14.7, " "))
         204 FORMAT("OSCDFT DEBUG: ", A, ": ", *(ES14.7, " + ", ES14.7, " i, ", :))
         205 FORMAT("OSCDFT DEBUG: ", A, ": ", *(I5, " ", :))

         400 FORMAT("OSCDFT: ik: ", I3)
         500 FORMAT("OSCDFT DEBUG: isk(", I3, "): ", I1,&
            "; xk(:,", I3, "): ", 3(ES14.7, " "),&
            "; wk(", I3, "): ", ES14.7)
         600 FORMAT("OSCDFT: Initial State");
         601 FORMAT("OSCDFT: Final State");
      END SUBROUTINE get_et

      SUBROUTINE oscdft_nelup_neldw_from_input(nelup_inp, neldw_inp)
         USE klist,     ONLY : nelup, neldw
         IMPLICIT NONE

         REAL(DP), INTENT(IN) :: nelup_inp, neldw_inp

         IF (nelup_inp /= 0.D0) THEN
            WRITE(stdout, 100) nelup_inp
            nelup = nelup_inp
         END IF
         IF (neldw_inp /= 0.D0) THEN
            WRITE(stdout, 101) neldw_inp
            neldw = neldw_inp
         END IF

         100 FORMAT("Replacing nelup from input: ", F0.4)
         101 FORMAT("Replacing neldw from input: ", F0.4)
      END SUBROUTINE oscdft_nelup_neldw_from_input

      SUBROUTINE oscdft_et_print_clock
         IMPLICIT NONE

         WRITE(stdout, '(/,5X,"OSCDFT_ET routines")')
                          !oscdft_xxxxx 12 chars max
         CALL print_clock("read_file")
         CALL print_clock("oscdft_wfcO")
         CALL print_clock("oscdft_ns")
         CALL print_clock("get_et")
         CALL print_clock("oscdft_rho")
      END SUBROUTINE oscdft_et_print_clock
#endif
END MODULE oscdft_et_mod
