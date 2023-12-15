!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine localdos (ldos, ldoss, becsum1, dos_ef)
  !-----------------------------------------------------------------------
  !
  !    This routine compute the local and total density of state at Ef
  !
  !    Note: this routine use psic as auxiliary variable. it should alread
  !          be defined
  !
  !    NB: this routine works only with gamma
  !
  !
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : omega
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE ener,             ONLY : ef
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : invfft, fft_interpolate
  USE buffers,          ONLY : get_buffer
  USE gvecs,            ONLY : doublegrid
  USE klist,            ONLY : xk, wk, ngk, igk_k, degauss, ngauss, ltetra
  USE lsda_mod,         ONLY : nspin, lsda, current_spin, isk
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wvfct,            ONLY : nbnd, npwx, et
  USE becmod,           ONLY : calbec, bec_type, allocate_bec_type_acc, deallocate_bec_type_acc
  USE wavefunctions,    ONLY : evc, psic, psic_nc
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc_d
#endif
  USE uspp,             ONLY : okvan, nkb, vkb
  USE uspp_param,       ONLY : upf, nh, nhm
  USE qpoint,           ONLY : nksq, ikks
  USE control_lr,       ONLY : nbnd_occ
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE mp_pools,         ONLY : inter_pool_comm
  USE mp,               ONLY : mp_sum
  USE dfpt_tetra_mod,   ONLY : dfpt_tetra_delta
  USE uspp_init,        ONLY : init_us_2
  USE control_flags,    ONLY : offload_type

  implicit none

  complex(DP) :: ldos (dfftp%nnr, nspin_mag), ldoss (dffts%nnr, nspin_mag)
  ! output: the local density of states at Ef
  ! output: the local density of states at Ef without augmentation
  REAL(DP) :: becsum1 ((nhm * (nhm + 1))/2, nat, nspin_mag)
  ! output: the local becsum at ef
  real(DP) :: dos_ef, check
  ! output: the density of states at Ef
  !
  !    local variables for Ultrasoft PP's
  !
  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, nt
  ! counters
  complex(DP), allocatable :: becsum1_nc(:,:,:,:)
  TYPE(bec_type) :: becp
  !
  ! local variables
  !
  real(DP) :: weight, w1, wdelta
  ! weights
  real(DP), external :: w0gauss
  !
  integer :: npw, ik, is, ig, ibnd, j, is1, is2, v_siz
  ! counters
  integer :: ios
  ! status flag for i/o
  !
  !  initialize ldos and dos_ef
  !
  ! For device buffer
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nl_d(:)
  !
  nl_d  => dffts%nl_d
  evc_d = evc
#else
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
#endif
  v_siz = dffts%nnr

  call start_clock ('localdos')
  IF (noncolin) THEN
     allocate (becsum1_nc( (nhm * (nhm + 1)) / 2, nat, npol, npol))
     becsum1_nc=(0.d0,0.d0)
  ENDIF

  call allocate_bec_type_acc (nkb, nbnd, becp)

  becsum1 (:,:,:) = 0.d0
  ldos (:,:) = (0d0, 0.0d0)
  ldoss(:,:) = (0d0, 0.0d0)
  dos_ef = 0.d0
  !
  !  loop over kpoints
  !
  !$acc data create(psic, psic_nc) copy(ldoss) 
  do ik = 1, nksq
     if (lsda) current_spin = isk (ikks(ik))
     npw = ngk(ikks(ik))
     weight = wk (ikks(ik))
     !
     ! unperturbed wfs in reciprocal space read from unit iuwfc
     !
     if (nksq > 1) then
             call get_buffer (evc, lrwfc, iuwfc, ikks(ik))
#if defined(__CUDA)
             evc_d = evc
#endif
     endif
     call init_us_2 (npw, igk_k(1,ikks(ik)), xk (1, ikks(ik)), vkb, .true.)
     !
     !$acc data copyin(evc) present(vkb, becp)
     call calbec ( offload_type, npw, vkb, evc, becp)
     !$acc end data
     !
     do ibnd = 1, nbnd_occ (ikks(ik))
        !
        if(ltetra) then
           wdelta = dfpt_tetra_delta(ibnd,ikks(ik))
        else
           wdelta = w0gauss ( (ef-et(ibnd,ikks(ik))) / degauss, ngauss) / degauss
        end if
        !
        w1 = weight * wdelta / omega
        !
        ! unperturbed wf from reciprocal to real space
        !
        IF (noncolin) THEN
           !$acc kernels
           psic_nc = (0.d0, 0.d0)
           !$acc end kernels
           !$acc parallel loop present(igk_k, psic_nc)
           do ig = 1, npw
#if defined(__CUDA)
              psic_nc (nl_d (igk_k(ig,ikks(ik))), 1 ) = evc_d (ig, ibnd)
              psic_nc (nl_d (igk_k(ig,ikks(ik))), 2 ) = evc_d (ig+npwx, ibnd)
#else
              psic_nc (nl_d (igk_k(ig,ikks(ik))), 1 ) = evc (ig, ibnd)
              psic_nc (nl_d (igk_k(ig,ikks(ik))), 2 ) = evc (ig+npwx, ibnd)
#endif
           enddo
           !$acc end parallel loop
           !$acc host_data use_device(psic_nc)
           CALL invfft ('Rho', psic_nc(:,1), dffts)
           CALL invfft ('Rho', psic_nc(:,2), dffts)
           !$acc end host_data
           !$acc parallel loop present(psic_nc, ldoss)
           do j = 1, v_siz
              ldoss (j, 1) = ldoss (j, 1) + &
                    w1 * ( DBLE(psic_nc(j,1))**2+AIMAG(psic_nc(j,1))**2 + &
                           DBLE(psic_nc(j,2))**2+AIMAG(psic_nc(j,2))**2)
           enddo
           !$acc end parallel loop
           IF (nspin_mag==4) THEN
              !$acc parallel loop present(psic_nc, ldoss)
              DO j = 1, v_siz
              !
                 ldoss(j,2) = ldoss(j,2) + w1*2.0_DP* &
                             (DBLE(psic_nc(j,1))* DBLE(psic_nc(j,2)) + &
                             AIMAG(psic_nc(j,1))*AIMAG(psic_nc(j,2)))

                 ldoss(j,3) = ldoss(j,3) + w1*2.0_DP* &
                             (DBLE(psic_nc(j,1))*AIMAG(psic_nc(j,2)) - &
                              DBLE(psic_nc(j,2))*AIMAG(psic_nc(j,1)))

                 ldoss(j,4) = ldoss(j,4) + w1* &
                             (DBLE(psic_nc(j,1))**2+AIMAG(psic_nc(j,1))**2 &
                             -DBLE(psic_nc(j,2))**2-AIMAG(psic_nc(j,2))**2)
              !
              END DO
              !$acc end parallel loop
           END IF
        ELSE
           !$acc kernels
           psic (:) = (0.d0, 0.d0)
           !$acc end kernels
           !$acc parallel loop present(psic)
           do ig = 1, npw
#if defined(__CUDA)
              psic (nl_d (igk_k(ig,ikks(ik)) ) ) = evc_d (ig, ibnd)
#else
              psic (nl_d (igk_k(ig,ikks(ik)) ) ) = evc (ig, ibnd)
#endif
           enddo
           !$acc end parallel loop
           !$acc host_data use_device(psic)
           CALL invfft ('Rho', psic, dffts)
           !$acc end host_data
           !$acc parallel loop present(ldoss, psic)
           do j = 1, v_siz
              ldoss (j, current_spin) = ldoss (j, current_spin) + &
                    w1 * ( DBLE ( psic (j) ) **2 + AIMAG (psic (j) ) **2)
           enddo
           !$acc end parallel loop
        END IF
        !
        !    If we have a US pseudopotential we compute here the becsum term
        !
        if(noncolin) then
        !$acc update self(becp%nc)
        else
        !$acc update self(becp%k)
        endif
        !
        w1 = weight * wdelta
        ijkb0 = 0
        do nt = 1, ntyp
           if (upf(nt)%tvanp ) then
              do na = 1, nat
                 if (ityp (na) == nt) then
                    ijh = 1
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       IF (noncolin) THEN
                          DO is1=1,npol
                             DO is2=1,npol
                                becsum1_nc (ijh, na, is1, is2) = &
                                becsum1_nc (ijh, na, is1, is2) + w1 * &
                                 (CONJG(becp%nc(ikb,is1,ibnd))* &
                                        becp%nc(ikb,is2,ibnd))
                             END DO
                          END DO
                       ELSE
                          becsum1 (ijh, na, current_spin) = &
                            becsum1 (ijh, na, current_spin) + w1 * &
                             DBLE (CONJG(becp%k(ikb,ibnd))*becp%k(ikb,ibnd) )
                       ENDIF
                       ijh = ijh + 1
                       do jh = ih + 1, nh (nt)
                          jkb = ijkb0 + jh
                          IF (noncolin) THEN
                             DO is1=1,npol
                                DO is2=1,npol
                                   becsum1_nc(ijh,na,is1,is2) = &
                                      becsum1_nc(ijh,na,is1,is2) + w1* &
                                      (CONJG(becp%nc(ikb,is1,ibnd))* &
                                             becp%nc(jkb,is2,ibnd) )
                                END DO
                             END DO
                          ELSE
                             becsum1 (ijh, na, current_spin) = &
                               becsum1 (ijh, na, current_spin) + w1 * 2.d0 * &
                                DBLE(CONJG(becp%k(ikb,ibnd))*becp%k(jkb,ibnd) )
                          END IF
                          ijh = ijh + 1
                       enddo
                    enddo
                    ijkb0 = ijkb0 + nh (nt)
                 endif
              enddo
           else
              do na = 1, nat
                 if (ityp (na) == nt) ijkb0 = ijkb0 + nh (nt)
              enddo
           endif
        enddo
        dos_ef = dos_ef + weight * wdelta
     enddo

  enddo
  !$acc end data
  if (doublegrid) then
     do is = 1, nspin_mag
        call fft_interpolate (dffts, ldoss (:, is), dfftp, ldos (:, is))
     enddo
  else
     ldos (:,:) = ldoss (:,:)
  endif

  IF (noncolin.and.okvan) THEN
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==nt) THEN
                 IF (upf(nt)%has_so) THEN
                    CALL transform_becsum_so(becsum1_nc,becsum1,na)
                 ELSE
                    CALL transform_becsum_nc(becsum1_nc,becsum1,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF

  call addusldos (ldos, becsum1)
  !
  !    Collects partial sums on k-points from all pools
  !
  call mp_sum ( ldoss, inter_pool_comm )
  call mp_sum ( ldos, inter_pool_comm )
  call mp_sum ( dos_ef, inter_pool_comm )
  call mp_sum ( becsum1, inter_pool_comm )
  !check
  !      check =0.d0
  !      do is=1,nspin_mag
  !         call fwfft('Rho',ldos(:,is),dfftp)
  !         check = check + omega* DBLE(ldos(nl(1),is))
  !         call invfft('Rho',ldos(:,is),dfftp)
  !      end do
  !      WRITE( stdout,*) ' check ', check, dos_ef
  !check
  !
  IF (noncolin) deallocate(becsum1_nc)
  call deallocate_bec_type_acc(becp)

  call stop_clock ('localdos')
  return
end subroutine localdos
