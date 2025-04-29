MODULE oscdft_forces_subs
#if defined (__OSCDFT)
   USE kinds,                    ONLY : DP
   USE io_global,                ONLY : stdout
   USE oscdft_wavefunction,      ONLY : oscdft_wavefunction_type
   USE oscdft_wavefunction_subs, ONLY : oscdft_get_buffer
   USE oscdft_context,           ONLY : oscdft_context_type,&
                                        oscdft_ns_type
   USE oscdft_input,             ONLY : oscdft_input_type
   USE oscdft_indices,           ONLY : oscdft_indices_type,&
                                        oscdft_constr_indices_type,&
                                        oscdft_orbital_indices_type
   USE oscdft_forces
   USE oscdft_enums
   PRIVATE
   PUBLIC oscdft_apply_forces, oscdft_print_forces

   CONTAINS
      SUBROUTINE oscdft_apply_forces(ctx)
         USE force_mod, ONLY : force
         USE ions_base, ONLY : nat, ityp

         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (.NOT.ctx%inp%skip_forces .AND. ctx%idx%nconstr > 0) THEN
            IF (.NOT. ALLOCATED(ctx%force_oscdft)) ALLOCATE(ctx%force_oscdft(3,nat))
            ctx%force_oscdft(:,:) = 0.D0
            CALL oscdft_get_forces(ctx, ctx%force_oscdft)

            force(:,:) = force(:,:) + ctx%force_oscdft(:,:)
         END IF
      END SUBROUTINE oscdft_apply_forces

      SUBROUTINE oscdft_print_forces(ctx)
         USE ions_base, ONLY : nat, ityp
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp
         INTEGER                                  :: na

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (.NOT.ctx%inp%skip_forces .AND. ctx%idx%nconstr > 0) THEN
            WRITE(stdout, '(/,5x, "OS-CDFT contribution to forces:")')
            DO na=1,nat
               WRITE(stdout, 101) na, ityp(na), ctx%force_oscdft(:,na)
            END DO
         END IF
         101 FORMAT(5X,'atom ',I4,' type ',I2,'   oscdft force = ',3F14.8)
      END SUBROUTINE oscdft_print_forces

      SUBROUTINE oscdft_get_forces(ctx, force_oscdft)
         USE ions_base,          ONLY : nat, ityp
         USE lsda_mod,           ONLY : isk, lsda
         USE klist,              ONLY : ngk, nks, igk_k, xk
         USE becmod,             ONLY : bec_type, becp, calbec,&
                                        allocate_bec_type_acc, &
                                        deallocate_bec_type_acc
         USE uspp_param,         ONLY : upf, nh
         USE wavefunctions,      ONLY : evc
         USE control_flags,      ONLY : gamma_only, use_gpu, offload_type
         USE io_files,           ONLY : iunwfc, nwordwfc
         USE wvfct,              ONLY : nbnd, npwx, wg
         USE uspp,               ONLY : nkb, vkb, ofsbeta, qq_at
         USE oscdft_wfcO,        ONLY : orthoOwfc
         USE oscdft_occupations, ONLY : oscdft_new_ns,&
                                        oscdft_get_occupation_numbers
         USE buffers,            ONLY : get_buffer
         USE mp,                 ONLY : mp_sum
         USE mp_pools,           ONLY : inter_pool_comm
         USE mp_images,          ONLY : intra_image_comm
         USE symm_base,          ONLY : nsym
         USE uspp_init,          ONLY : init_us_2

         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         REAL(DP),                  INTENT(OUT)           :: force_oscdft(3,nat)
         TYPE(oscdft_input_type),        POINTER          :: inp
         TYPE(oscdft_indices_type),      POINTER          :: idx
         TYPE(oscdft_ns_type),           POINTER          :: nst
         TYPE(oscdft_forces_type),       POINTER          :: forces
         TYPE(oscdft_wavefunction_type), POINTER          :: wfcF

         INTEGER, EXTERNAL        :: set_Hubbard_l
         INTEGER                  :: ispin, npw, beta_indices_all(nat), beta_indices(nat),&
                                     curr_dim, h, k, ibnd, iconstr, ih, jh,&
                                     ijkb0, ik, ioscdft, ipos, na, nt, ierr, oidx, i,&
                                     isum, osum, angular_momentum, isym
         TYPE(bec_type)           :: proj
         REAL(DP)                 :: multiplier
         COMPLEX(DP), ALLOCATABLE :: spsi(:,:)
         REAL(DP),    ALLOCATABLE :: force_oscdft_sym(:,:,:),&
                                     dns(:,:,:,:)
         REAL(DP)                 :: temp_gam
         COMPLEX(DP)              :: temp_k
         COMPLEX(DP), POINTER     :: wfcF_wfc(:,:)

         inp    => ctx%inp
         idx    => ctx%idx
         nst    => ctx%nst
         forces => ctx%forces
         wfcF   => forces%wfcO

         CALL start_clock("oscdft_force")

         ! recalculate atomic wavefunctions
         CALL orthoOwfc(ctx)

         CALL allocate_bec_type_acc(nkb, nbnd, becp)
         CALL allocate_bec_type_acc(wfcF%n, nbnd, proj)

         ! allocate derivatives
         IF (gamma_only) THEN
            CALL oscdft_derivatives_gamma_alloc(forces%deriv_gamma,&
                                                wfcF%n, nbnd, nkb)
         ELSE
            CALL oscdft_derivatives_k_alloc(forces%deriv_k,&
                                            wfcF%n, nbnd, nkb)
         ENDIF

         ! get ns and occupations
         CALL oscdft_new_ns(ctx, "oscdft_get_forces")
         CALL oscdft_get_occupation_numbers(ctx, .true.)

         ALLOCATE(spsi(npwx,nbnd))
         ALLOCATE(force_oscdft_sym(3,nsym,nat),&
                  dns(2*idx%max_ns_dim+1,2*idx%max_ns_dim+1,nsym,idx%nconstr))
         wfcF_wfc => wfcF%wfc
         !$acc data create(spsi, wfcF_wfc(:,:)) present_or_create(evc) present_or_copyin(wg)

         force_oscdft_sym(:,:,:) = 0.D0
         ispin = 1
         DO ik=1,nks
            IF (lsda) ispin = isk(ik)
            npw = ngk(ik)

            IF (nks > 1) THEN
               CALL get_buffer(evc, nwordwfc, iunwfc, ik)
               CALL oscdft_get_buffer(wfcF, ik)
            ENDIF
            !$acc update device(wfcF_wfc(:,:), evc)

            CALL init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
            ! proj = <wfcF|S|psi>
            CALL calbec(offload_type, npw, vkb, evc, becp)
            IF (use_gpu) THEN
               CALL s_psi_acc(npwx, npw, nbnd, evc, spsi)
            ELSE
               CALL s_psi(npwx, npw, nbnd, evc, spsi)
            END IF
            CALL calbec(offload_type, npw, wfcF_wfc, spsi, proj)

            ! calculate <beta|psi>
            IF (gamma_only) THEN
               CALL calc_betapsi_gamma(forces%deriv_gamma%betapsi, becp%r)
            ELSE
               CALL calc_betapsi_k(forces%deriv_k%betapsi, becp%k)
            END IF

            DO ipos=1,3 ! x y z
               IF (gamma_only) THEN
                  CALL oscdft_get_derivatives_gamma(ctx, spsi, ik, ispin, ipos)
               ELSE
                  CALL oscdft_get_derivatives_k(ctx, spsi, ik, ispin, ipos)
               ENDIF
               DO na=1,nat
                  IF (gamma_only) THEN
                     CALL oscdft_dndtau_gamma(ctx, forces%deriv_gamma, proj%r, na, ipos, ispin, ik, dns)
                  ELSE
                     CALL oscdft_dndtau_k(ctx, forces%deriv_k, proj%k, na, ipos, ispin, ik, dns)
                  ENDIF
                  DO iconstr=1,idx%nconstr
                     ioscdft = idx%iconstr2ioscdft(iconstr)
                     oidx = inp%occup_index(ioscdft)
                     curr_dim = idx%ns_dim(ioscdft)

                     ! WRITE(stdout, 100) ik, ipos, na, iconstr
                     ! DO h=1,curr_dim
                     !    WRITE(stdout, 700) dns(h,1:curr_dim,iconstr)
                     ! ENDDO
                     DO isym=1,nsym
                        DO h=1,curr_dim
                           DO k=1,curr_dim
                              IF (oidx.EQ.OCCUP_TRACE) THEN
                                 multiplier = ctx%multipliers(iconstr)*&
                                              MERGE(1.D0, 0.D0, h.EQ.k)
                              ELSE IF (oidx.EQ.OCCUP_SUM) THEN
                                 osum = inp%occup_index_sum(2,ioscdft)
                                 multiplier = ctx%multipliers(iconstr)*&
                                              nst%occup_eigvects(h,osum,ioscdft)*&
                                              nst%occup_eigvects(k,osum,ioscdft)
                                 ! DO isum=2,inp%occup_index_sum(1,ioscdft)
                                 !    osum = inp%occup_index_sum(2,ioscdft)
                                 !    multiplier = multiplier+&
                                 !                 ctx%multipliers(iconstr)*&
                                 !                 nst%eigvects(h,osum,ioscdft)*&
                                 !                 nst%eigvects(k,osum,ioscdft)
                                 ! ENDDO
                              ELSE
                                 multiplier = ctx%multipliers(iconstr)*&
                                              nst%occup_eigvects(h,oidx,iconstr)*&
                                              nst%occup_eigvects(k,oidx,iconstr)
                              END IF
                              force_oscdft_sym(ipos,isym,na) =&
                                 force_oscdft_sym(ipos,isym,na) - multiplier * dns(h,k,isym,iconstr)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         IF (gamma_only) THEN
            CALL oscdft_derivatives_gamma_dealloc(forces%deriv_gamma)
         ELSE
            CALL oscdft_derivatives_k_dealloc(forces%deriv_k)
         ENDIF

         !$acc end data
         DEALLOCATE(spsi, dns)

         CALL deallocate_bec_type_acc(becp)
         CALL deallocate_bec_type_acc(proj)

         CALL mp_sum(force_oscdft_sym, inter_pool_comm)

         CALL oscdft_symmetrize(nat,force_oscdft, force_oscdft_sym)

         DEALLOCATE(force_oscdft_sym)

         CALL stop_clock("oscdft_force")

         100 FORMAT("OSCDFT DEBUG: dns; ik: ", I2, ", ipos: ", I1, ", na: ", I2, ", iconstr: ", I2)
         200 FORMAT("OSCDFT DEBUG: h: ", I1, "; k: ", I1, "; multipliers: ", ES14.7,&
                    "; h eigvect: ", ES14.7, "; k eigvect: ", ES14.7,&
                    "; multiplier: ", ES14.7, "; curr_force: ", ES14.7)
         300 FORMAT("OSCDFT DEBUG: force loop: ik: ", I5, "; ipos: ", I5)
         700 FORMAT(*(ES14.7, " "))
      END SUBROUTINE oscdft_get_forces

      SUBROUTINE calc_betapsi_gamma(betapsi, becp_r)
         USE wvfct,              ONLY : nbnd
         USE ions_base,          ONLY : nat, ityp
         USE uspp,               ONLY : ofsbeta, qq_at
         USE uspp_param,         ONLY : nh
         USE becmod,             ONLY : becp
         IMPLICIT NONE

         REAL(DP), INTENT(OUT) :: betapsi(:,:)
         REAL(DP), INTENT(IN)  :: becp_r(:,:)
         INTEGER               :: ibnd, na, nt, ijkb0, ih, jh
         REAL(DP)              :: temp



         !$acc data present(betapsi, qq_at, becp_r)

         !$acc kernels
         betapsi(:,:) = 0.D0
         !$acc end kernels


#if defined(__CUDA)
!$acc parallel loop private(na, ijkb0, nt, ih, temp, jh)
#else
!$omp parallel do private(na, ijkb0, nt, ih, temp, jh)
#endif
               DO ibnd=1,nbnd
                  DO na=1,nat
                     nt = ityp(na)
                     ijkb0 = ofsbeta(na)
                     DO ih=1,nh(nt)
                        temp = 0.D0
                        DO jh=1,nh(nt)
                           temp = temp + qq_at(ih,jh,na) * becp_r(ijkb0+jh,ibnd)
                        END DO
                        betapsi(ijkb0+ih,ibnd) = betapsi(ijkb0+ih,ibnd) + temp
                     END DO
                  END DO
               END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif
         !$acc end data
      END SUBROUTINE calc_betapsi_gamma

      SUBROUTINE calc_betapsi_k(betapsi, becp_k)
         USE wvfct,              ONLY : nbnd
         USE ions_base,          ONLY : nat, ityp
         USE uspp,               ONLY : ofsbeta, qq_at
         USE uspp_param,         ONLY : nh
         USE becmod,             ONLY : becp
         IMPLICIT NONE

         COMPLEX(DP), INTENT(OUT) :: betapsi(:,:)
         COMPLEX(DP), INTENT(IN)  :: becp_k(:,:)
         INTEGER                  :: ibnd, na, nt, ijkb0, ih, jh
         COMPLEX(DP)              :: temp



         !$acc data present(betapsi, qq_at, becp_k)

         !$acc kernels
         betapsi(:,:) = 0.D0
         !$acc end kernels

#if defined(__CUDA)
!$acc parallel loop private(na, ijkb0, nt, ih, temp, jh)
#else
!$omp parallel do private(na, ijkb0, nt, ih, temp, jh)
#endif
               DO ibnd=1,nbnd
                  DO na=1,nat
                     nt = ityp(na)
                     ijkb0 = ofsbeta(na)
                     DO ih=1,nh(nt)
                        temp = (0.D0, 0.D0)
                        DO jh=1,nh(nt)
                           temp = temp + qq_at(ih,jh,na) * becp_k(ijkb0+jh,ibnd)
                        END DO
                        betapsi(ijkb0+ih,ibnd) = betapsi(ijkb0+ih,ibnd) + temp
                     END DO
                  END DO
               END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif
         !$acc end data
      END SUBROUTINE calc_betapsi_k

      SUBROUTINE oscdft_symmetrize(nat, dst, src)
         USE symm_base, ONLY : nsym, s, irt
         USE cell_base, ONLY : at, bg
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: nat
         REAL(DP), INTENT(OUT) :: dst(3,nat)
         REAL(DP), INTENT(IN)  :: src(3,nsym,nat)

         REAL(DP), ALLOCATABLE :: tmp(:,:,:)
         INTEGER               :: na, isym, nb

         IF (nsym == 1) THEN
            dst(1:3,1:nat) = src(1:3,1,1:nat)
         ELSE
            ALLOCATE(tmp(3,nsym,nat))
            DO na=1,nat
               DO isym=1,nsym
                  tmp(:,isym,na) = src(1,isym,na)*at(1,:) + &
                                   src(2,isym,na)*at(2,:) + &
                                   src(3,isym,na)*at(3,:)
               END DO
            END DO

            dst(:,:) = 0.D0
            DO na=1,nat
               DO isym=1,nsym
                  nb = irt(isym,na)
                  dst(:,na) = dst(:,na)+&
                     s(:,1,isym) * tmp(1,isym,nb)+&
                     s(:,2,isym) * tmp(2,isym,nb)+&
                     s(:,3,isym) * tmp(3,isym,nb)
               END DO
            END DO
            tmp(:,1,:) = dst(:,:) / DBLE(nsym)

            DO na=1,nat
               dst(:,na) = tmp(1,1,na)*bg(:,1)+&
                           tmp(2,1,na)*bg(:,2)+&
                           tmp(3,1,na)*bg(:,3)
            END DO

            DEALLOCATE(tmp)
         ENDIF
      END SUBROUTINE oscdft_symmetrize

      SUBROUTINE oscdft_get_derivatives_gamma(ctx, spsi, ik, ispin, ipos)
         !TODO: should use BLAS
         USE becmod,        ONLY : calbec
         USE wvfct,         ONLY : npwx, nbnd
         USE klist,         ONLY : ngk, igk_k
         USE wavefunctions, ONLY : evc
         USE ions_base,     ONLY : nat, ityp
         USE uspp,          ONLY : nkb, vkb, ofsbeta, qq_at
         USE cell_base,     ONLY : tpiba
         USE uspp_param,    ONLY : nh
         USE gvect,         ONLY : g
         USE symm_base,     ONLY : nsym
         USE mp_bands,      ONLY : intra_bgrp_comm
         USE control_flags, ONLY : offload_type

         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),             POINTER     :: inp
         TYPE(oscdft_indices_type),           POINTER     :: idx
         TYPE(oscdft_constr_indices_type),    POINTER     :: constr
         TYPE(oscdft_wavefunction_type),      POINTER     :: wfcF
         TYPE(oscdft_derivatives_gamma_type), POINTER     :: deriv
         COMPLEX(DP), INTENT(IN)                          :: spsi(:,:)
         INTEGER, INTENT(IN)                              :: ik, ipos, ispin
         REAL(DP), ALLOCATABLE    :: dbetapsi_temp(:,:) ! size(nkb,nbnd)
         COMPLEX(DP), ALLOCATABLE :: dwfc(:,:),& ! size(npwx,nwfcO)
                                     dbeta(:,:) ! size(npwx,nkb)
         INTEGER                  :: npw, ibnd, iconstr, ih, jh, ijkb0, ikb,&
                                     ioscdft, ipw, ldim, m1, na, nt, isym, m1_off
         REAL(DP)                 :: temp

         ! openacc does not like deriv%...
         REAL(DP), POINTER        :: dwfatpsi(:,:),&
                                     wfatbeta(:,:),&
                                     wfatdbeta(:,:),&
                                     dbetapsi(:,:)
         COMPLEX(DP), POINTER :: wfcF_wfc(:,:)

         inp    => ctx%inp
         idx    => ctx%idx
         constr => idx%constr
         wfcF   => ctx%forces%wfcO
         deriv  => ctx%forces%deriv_gamma
         wfcF_wfc => ctx%forces%wfcO%wfc

         npw = ngk(ik)

         dwfatpsi  => deriv%dwfatpsi
         wfatbeta  => deriv%wfatbeta
         wfatdbeta => deriv%wfatdbeta
         dbetapsi  => deriv%dbetapsi

         !$acc data present(dwfatpsi, wfatbeta, wfatdbeta, dbetapsi, wfcF_wfc(:,:))

         ALLOCATE(dwfc(npwx,wfcF%n))
         !$acc data create(dwfc)

         !$acc kernels
         dwfatpsi(:,:) = 0.D0
         wfatbeta(:,:) = 0.D0
         wfatdbeta(:,:) = 0.D0
         dbetapsi(:,:) = 0.D0
         dwfc = (0.D0, 0.D0)
         !$acc end kernels

         ! only set this if na is constrained
         ! others is (0.D0, 0.D0)
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            ldim = idx%ns_dim(ioscdft)
            IF (inp%spin_index(ioscdft).EQ.ispin) THEN
               DO isym=1,nsym
                  DO m1=1,ldim
                     m1_off = wfcF%get_offset(constr, m1, iconstr, isym)
#if defined(__CUDA)
                     !$acc parallel loop present(g, igk_k)
#else
!$omp parallel do
#endif
                     DO ipw=1,npw
                        dwfc(ipw,m1_off) =&
                           (0.D0,-1.D0)*tpiba*g(ipos,igk_k(ipw,ik)) *&
                           wfcF_wfc(ipw,m1_off)
                     END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif
                  END DO
               END DO
            END IF
         END DO

         CALL calbec(offload_type, npw, dwfc, spsi, dwfatpsi)
         !$acc end data
         DEALLOCATE(dwfc)

         ALLOCATE(dbeta(npwx,nkb),dbetapsi_temp(nkb,nbnd))
         !$acc data create(dbeta, dbetapsi_temp)
         !$acc kernels
         dbeta = (0.D0,0.D0)
         dbetapsi_temp = (0.D0, 0.D0)
         !$acc end kernels

         CALL calbec(offload_type, npw, wfcF_wfc, vkb, wfatbeta)

#if defined(__CUDA)
         !$acc parallel loop collapse(2) present(g, igk_k, vkb, dbeta)
#endif
         DO ikb=1,nkb
#if !defined(__CUDA)
!$omp parallel do
#endif
            DO ipw=1,npw
               dbeta(ipw,ikb) = (0.D0,-1.D0) * tpiba * g(ipos,igk_k(ipw,ik)) * vkb(ipw,ikb)
            END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif
         END DO

         CALL calbec(offload_type, npw, dbeta, evc, dbetapsi_temp)
         CALL calbec(offload_type, npw, wfcF_wfc, dbeta, wfatdbeta)

#if defined(__CUDA)
!$acc parallel loop private(na, nt, ijkb0, ih, jh, temp)
#else
!$omp parallel do private(na, nt, ijkb0, ih, jh, temp)
#endif
         DO ibnd=1,nbnd
            DO na=1,nat
               nt = ityp(na)
               ijkb0 = ofsbeta(na)
               DO ih=1,nh(nt)
                  temp = 0.D0
                  DO jh=1,nh(nt)
                     temp = temp + qq_at(ih,jh,na) * dbetapsi_temp(ijkb0+jh,ibnd)
                  END DO
                  dbetapsi(ijkb0+ih,ibnd) = dbetapsi(ijkb0+ih,ibnd) + temp
               END DO
            END DO
         END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif

         !$acc end data
         DEALLOCATE(dbeta,dbetapsi_temp)
         !$acc end data
      END SUBROUTINE oscdft_get_derivatives_gamma

      SUBROUTINE oscdft_get_derivatives_k(ctx, spsi, ik, ispin, ipos)
         !TODO: should use BLAS
         USE becmod,        ONLY : calbec
         USE wvfct,         ONLY : nbnd, npwx
         USE klist,         ONLY : ngk, igk_k
         USE wavefunctions, ONLY : evc
         USE ions_base,     ONLY : nat, ityp
         USE uspp,          ONLY : nkb
         USE cell_base,     ONLY : tpiba
         USE uspp,          ONLY : vkb, ofsbeta, qq_at
         USE gvect,         ONLY : g
         USE uspp_param,    ONLY : nh
         USE symm_base,     ONLY : nsym
         USE control_flags, ONLY : offload_type

         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),         POINTER         :: inp
         TYPE(oscdft_indices_type),       POINTER         :: idx
         TYPE(oscdft_wavefunction_type),  POINTER         :: wfcF
         TYPE(oscdft_constr_indices_type), POINTER       :: constr
         TYPE(oscdft_derivatives_k_type), POINTER         :: deriv
         COMPLEX(DP), INTENT(IN)  :: spsi(:,:)
         INTEGER,     INTENT(IN)  :: ik, ipos, ispin
         COMPLEX(DP), ALLOCATABLE :: dbetapsi_temp(:,:) ! size(nkb,nbnd)
         COMPLEX(DP), ALLOCATABLE :: dwfc(:,:),& ! size(npwx,nwfcO)
                                     dbeta(:,:) ! size(npwx,nkb)
         INTEGER                  :: npw, ibnd, iconstr, ioscdft, na, ldim, ih, jh, nt,&
                                     ijkb0, ikb, ipw, m1, isym, m1_off
         COMPLEX(DP)              :: temp

         ! openacc does not like deriv%...
         COMPLEX(DP), POINTER     :: dwfatpsi(:,:),&
                                     wfatbeta(:,:),&
                                     wfatdbeta(:,:),&
                                     dbetapsi(:,:)

         COMPLEX(DP), POINTER :: wfcF_wfc(:,:)

         inp    => ctx%inp
         idx    => ctx%idx
         wfcF   => ctx%forces%wfcO
         deriv  => ctx%forces%deriv_k
         constr => idx%constr
         wfcF_wfc => ctx%forces%wfcO%wfc

         npw = ngk(ik)

         dwfatpsi  => deriv%dwfatpsi
         wfatbeta  => deriv%wfatbeta
         wfatdbeta => deriv%wfatdbeta
         dbetapsi  => deriv%dbetapsi

         !$acc data present(dwfatpsi, wfatbeta, wfatdbeta, dbetapsi, wfcF_wfc(:,:))

         ALLOCATE(dwfc(npwx,wfcF%n))
         !$acc data create(dwfc)

         !$acc kernels
         dwfatpsi(:,:) =  (0.D0, 0.D0)
         wfatbeta(:,:) =  (0.D0, 0.D0)
         wfatdbeta(:,:) = (0.D0, 0.D0)
         dbetapsi(:,:) =  (0.D0, 0.D0)
         dwfc = (0.D0, 0.D0)
         !$acc end kernels

         ! only set this if na is constrained
         ! others is (0.D0, 0.D0)
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            ldim = idx%ns_dim(ioscdft)
            IF (inp%spin_index(ioscdft).EQ.ispin) THEN
               DO isym=1,nsym
                  DO m1=1,ldim
                     m1_off = wfcF%get_offset(constr, m1, iconstr, isym)
#if defined(__CUDA)
!$acc parallel loop present(g, igk_k)
#else
!$omp parallel do
#endif
                     DO ipw=1,npw
                        dwfc(ipw,m1_off) =&
                           (0.D0,-1.D0)*tpiba*g(ipos,igk_k(ipw,ik)) * wfcF_wfc(ipw,m1_off)
                     END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif
                  END DO
               END DO
            END IF
         END DO

         CALL calbec(offload_type, npw, dwfc, spsi, dwfatpsi)
         !$acc end data
         DEALLOCATE(dwfc)

         ALLOCATE(dbeta(npwx,nkb),dbetapsi_temp(nkb,nbnd))
         !$acc data create(dbeta, dbetapsi_temp)
         !$acc kernels
         dbeta = (0.D0, 0.D0)
         dbetapsi_temp = (0.D0, 0.D0)
         !$acc end kernels

         CALL calbec(offload_type, npw, wfcF_wfc, vkb, wfatbeta)

#if defined(__CUDA)
         !$acc parallel loop collapse(2) present(g, igk_k, vkb, dbeta)
#endif
         DO ikb=1,nkb
#if !defined(__CUDA)
!$omp parallel do
#endif
            DO ipw=1,npw
               dbeta(ipw,ikb) = (0.D0,-1.D0) * tpiba * g(ipos,igk_k(ipw,ik)) * vkb(ipw,ikb)
            ENDDO
#if !defined(__CUDA)
!$omp end parallel do
#endif
         ENDDO

         CALL calbec(offload_type, npw, dbeta, evc, dbetapsi_temp)
         CALL calbec(offload_type, npw, wfcF_wfc, dbeta, wfatdbeta)

#if defined(__CUDA)
!$acc parallel loop private(na, nt, ijkb0, ih, jh, temp)
#else
!$omp parallel do private(na, nt, ijkb0, ih, jh, temp)
#endif
         DO ibnd=1,nbnd
            DO na=1,nat
               nt = ityp(na)
               ijkb0 = ofsbeta(na)
               DO ih=1,nh(nt)
                  temp = (0.D0, 0.D0)
                  DO jh=1,nh(nt)
                     temp = temp + qq_at(ih,jh,na) * dbetapsi_temp(ijkb0+jh,ibnd)
                  END DO
                  dbetapsi(ijkb0+ih,ibnd) = dbetapsi(ijkb0+ih,ibnd) + temp
               END DO
            END DO
         END DO
#if !defined(__CUDA)
!$omp end parallel do
#endif

         !$acc end data
         DEALLOCATE(dbeta,dbetapsi_temp)
         !$acc end data
      END SUBROUTINE oscdft_get_derivatives_k

      SUBROUTINE oscdft_dndtau_gamma(ctx, deriv, proj, na, ipos, ispin, ik, dns)
         USE wvfct,      ONLY : nbnd, npwx, wg
         USE uspp,       ONLY : nkb, ofsbeta
         USE ions_base,  ONLY : ityp
         USE uspp_param, ONLY : nh
         USE symm_base,  ONLY : nsym
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET   :: ctx
         TYPE(oscdft_input_type),           POINTER         :: inp
         TYPE(oscdft_indices_type),         POINTER         :: idx
         TYPE(oscdft_orbital_indices_type), POINTER         :: orbs
         TYPE(oscdft_constr_indices_type),  POINTER         :: constr
         TYPE(oscdft_wavefunction_type),    POINTER         :: wfcF
         TYPE(oscdft_derivatives_gamma_type), INTENT(INOUT), TARGET :: deriv
         REAL(DP), INTENT(IN)  :: proj(:,:)
         INTEGER,  INTENT(IN)  :: na, ipos, ispin, ik
         REAL(DP), INTENT(OUT) :: dns(:,:,:,:)
         INTEGER               :: curr_dim, iconstr, ioscdft, ibnd, m1, m2,&
                                  ijkb0, nt, isym, m1_off, m2_off, m1_na, m2_na
         REAL(DP)              :: temp, temp2

         REAL(DP), POINTER        :: dwfatpsi(:,:),&
                                     wfatbeta(:,:),&
                                     wfatdbeta(:,:),&
                                     dbetapsi(:,:),&
                                     betapsi(:,:),&
                                     wfatdpsi(:,:)

         inp    => ctx%inp
         idx    => ctx%idx
         orbs   => idx%orbs
         wfcF   => ctx%forces%wfcO
         constr => idx%constr

         dwfatpsi  => deriv%dwfatpsi
         wfatbeta  => deriv%wfatbeta
         wfatdbeta => deriv%wfatdbeta
         dbetapsi  => deriv%dbetapsi
         betapsi   => deriv%betapsi
         wfatdpsi  => deriv%wfatdpsi

         !$acc data present(dwfatpsi, wfatbeta, wfatdbeta,&
         !$acc&             dbetapsi, betapsi, wfatdpsi,&
         !$acc&             proj) present_or_copyin(wg)

         ijkb0 = ofsbeta(na)
         nt = ityp(na)
         !$acc host_data use_device(wfatdbeta,&
         !$acc&                     betapsi,&
         !$acc&                     wfatdpsi,&
         !$acc&                     wfatbeta,&
         !$acc&                     dbetapsi)
         CALL MYDGEMM('N','N',wfcF%n,nbnd,nh(nt),1.D0,&
                      wfatdbeta(1,ijkb0+1),wfcF%n,&
                      betapsi(ijkb0+1,1),nkb,&
                      0.D0,wfatdpsi,wfcF%n)
         CALL MYDGEMM('N','N',wfcF%n,nbnd,nh(nt),1.D0,&
                      wfatbeta(1,ijkb0+1),wfcF%n,&
                      dbetapsi(ijkb0+1,1),nkb,&
                      1.D0,wfatdpsi,wfcF%n)
         !$acc end host_data

         dns(:,:,:,:) = 0.D0
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            curr_dim = idx%ns_dim(ioscdft)
            IF (inp%spin_index(ioscdft).EQ.ispin) THEN
               DO isym=1,nsym
                  DO m1=1,curr_dim
                     m1_off = wfcF%get_offset(constr, m1, iconstr, isym)
                     m1_na  = orbs%iat_sym(isym,constr%ins2iorb(m1,iconstr))
                     DO m2=m1,curr_dim
                        m2_off = wfcF%get_offset(constr, m2, iconstr, isym)
                        m2_na  = orbs%iat_sym(isym,constr%ins2iorb(m2,iconstr))
                        temp = 0.D0
#if defined(__CUDA)
                        !$acc parallel loop reduction(+:temp) private(temp2)
#else
                        !$omp parallel do reduction(+:temp) private(temp2)
#endif
                        DO ibnd=1,nbnd
                           temp2 = (proj(m1_off,ibnd)*wfatdpsi(m2_off,ibnd)+&
                              wfatdpsi(m1_off,ibnd)*proj(m2_off,ibnd))
                           IF (m1_na.EQ.na) THEN
                              temp2 = temp2 + dwfatpsi(m1_off,ibnd)*proj(m2_off,ibnd)
                           END IF
                           IF (m2_na.EQ.na) THEN
                              temp2 = temp2 + proj(m1_off,ibnd)*dwfatpsi(m2_off,ibnd)
                           END IF
                           temp = temp + wg(ibnd,ik)*temp2
                        END DO
#if !defined(__CUDA)
                        !$omp end parallel do
#endif
                        dns(m1,m2,isym,iconstr) = temp
                        IF (m1.NE.m2) THEN
                           dns(m2,m1,isym,iconstr) = dns(m1,m2,isym,iconstr)
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         END DO
         !$acc end data
      END SUBROUTINE oscdft_dndtau_gamma

      SUBROUTINE oscdft_dndtau_k(ctx, deriv, proj, na, ipos, ispin, ik, dns)
         USE wvfct,      ONLY : nbnd, npwx, wg
         USE uspp,       ONLY : nkb, ofsbeta
         USE ions_base,  ONLY : ityp
         USE uspp_param, ONLY : nh
         USE symm_base,  ONLY : nsym
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),           POINTER       :: inp
         TYPE(oscdft_indices_type),         POINTER       :: idx
         TYPE(oscdft_wavefunction_type),    POINTER       :: wfcF
         TYPE(oscdft_orbital_indices_type), POINTER       :: orbs
         TYPE(oscdft_constr_indices_type),  POINTER       :: constr
         TYPE(oscdft_derivatives_k_type), INTENT(INOUT), TARGET :: deriv
         COMPLEX(DP), INTENT(IN)  :: proj(:,:)
         INTEGER,     INTENT(IN)  :: na, ipos, ispin, ik
         REAL(DP),    INTENT(OUT) :: dns(:,:,:,:)
         INTEGER                  :: curr_dim, iconstr, ioscdft, ibnd, m1, m2,&
                                     ijkb0, nt, isym, m1_off, m2_off, m1_na, m2_na
         REAL(DP)                 :: temp, temp2
         COMPLEX(DP), POINTER     :: dwfatpsi(:,:),&
                                     wfatbeta(:,:),&
                                     wfatdbeta(:,:),&
                                     dbetapsi(:,:),&
                                     betapsi(:,:),&
                                     wfatdpsi(:,:)

         inp    => ctx%inp
         idx    => ctx%idx
         orbs   => idx%orbs
         wfcF   => ctx%forces%wfcO
         constr => idx%constr

         dwfatpsi  => deriv%dwfatpsi
         wfatbeta  => deriv%wfatbeta
         wfatdbeta => deriv%wfatdbeta
         dbetapsi  => deriv%dbetapsi
         betapsi   => deriv%betapsi
         wfatdpsi  => deriv%wfatdpsi
         !$acc data present(dwfatpsi, wfatbeta, wfatdbeta,&
         !$acc&             dbetapsi, betapsi, wfatdpsi,&
         !$acc&             proj) present_or_copyin(wg)

         ijkb0 = ofsbeta(na)
         nt = ityp(na)
         !$acc host_data use_device(wfatdbeta,&
         !$acc&                     betapsi,&
         !$acc&                     wfatdpsi,&
         !$acc&                     wfatbeta,&
         !$acc&                     dbetapsi)
         CALL MYZGEMM('N','N',wfcF%n,nbnd,nh(nt),(1.D0,0.D0),&
                      wfatdbeta(1,ijkb0+1),wfcF%n,&
                      betapsi(ijkb0+1,1),nkb,&
                      (0.D0,0.D0),wfatdpsi,wfcF%n)
         CALL MYZGEMM('N','N',wfcF%n,nbnd,nh(nt),(1.D0,0.D0),&
                      wfatbeta(1,ijkb0+1),wfcF%n,&
                      dbetapsi(ijkb0+1,1),nkb,&
                      (1.D0,0.D0),wfatdpsi,wfcF%n)
         !$acc end host_data

         dns(:,:,:,:) = 0.D0
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            curr_dim = idx%ns_dim(ioscdft)
            IF (inp%spin_index(ioscdft).EQ.ispin) THEN
               DO isym=1,nsym
                  DO m1=1,curr_dim
                     m1_off = wfcF%get_offset(constr, m1, iconstr, isym)
                     m1_na  = orbs%iat_sym(isym,constr%ins2iorb(m1,iconstr))
                     DO m2=m1,curr_dim
                        m2_off = wfcF%get_offset(constr, m2, iconstr, isym)
                        m2_na  = orbs%iat_sym(isym,constr%ins2iorb(m2,iconstr))
                        temp = (0.D0, 0.D0)
#if defined(__CUDA)
                        !$acc parallel loop reduction(+:temp) private(temp2)
#else
                        !$omp parallel do reduction(+:temp) private(temp2)
#endif
                        DO ibnd=1,nbnd
                           temp2 = DBLE(proj(m1_off,ibnd)*CONJG(wfatdpsi(m2_off,ibnd))+&
                              wfatdpsi(m1_off,ibnd)*CONJG(proj(m2_off,ibnd)))
                           IF (m1_na.EQ.na) THEN
                              temp2 = temp2 + DBLE(dwfatpsi(m1_off,ibnd)*CONJG(proj(m2_off,ibnd)))
                           END IF
                           IF (m2_na.EQ.na) THEN
                              temp2 = temp2 + DBLE(proj(m1_off,ibnd)*CONJG(dwfatpsi(m2_off,ibnd)))
                           END IF
                           temp = temp + wg(ibnd,ik)*temp2
                        END DO
#if !defined(__CUDA)
                        !$omp end parallel do
#endif
                        dns(m1,m2,isym,iconstr) = temp
                        IF (m1.NE.m2) THEN
                           dns(m2,m1,isym,iconstr) = dns(m1,m2,isym,iconstr)
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         END DO
         !$acc end data
      END SUBROUTINE oscdft_dndtau_k
#endif
END MODULE oscdft_forces_subs
